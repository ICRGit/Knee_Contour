$!
$!         _/_/_/  _/_/_/    _/
$!          _/    _/    _/  _/           Image Processing Language
$!         _/    _/_/_/    _/
$!        _/    _/        _/             (c)  Nicolas Vilayphiou
$!     _/_/_/  _/        _/_/_/_/             Scanco Medical
$!
$!  IPL Batch Scanco
$!
$!  This script is based on IPLV6_UPAT_E5E6_XT2.COM (Version 6.2.2012)
$!
$!  Steven Boyd, May 30 2016
$!
$! Alice Escher, Andres Kroker, and Sarah Manske, August 15, 2016
$! Adapted to work on the knee:
$! increased threshold to 550 mg/HA
$! decreased the slicewise percentage to 1% from 50%
$!
$! Aug 25 - cleaned up
$!
$! Oct 14 - used erosion instead of peel_iter to ensure minimal cortical thickness
$!          since peel_iter only works in 2D. This perviously led to missing cortex
$!          at the weight bearing regions since cortex there is mostly in Z-plane
$!
$! March 24, 2017. Andres Kroker. improve memory management
$! -------------------------------------------------------------------------------
$! Define Parameters
$! -------------------------------------------------------------------------------
$! Step 1
$ IPL_SIGMA0	:=	2
$ IPL_SUPPORT0	:=	3
$ IPL_LOWER0	:=	550			! high threshold to get mostly cortical
$ IPL_PEEL0		:=  4			! min cortical thickness; TODO: check if 6 is better than 8
$ IPL_MISC1_0	:=	300000		! min number, tibia at 200000 and radius at 800
$ IPL_MISC1_1	:=  50			! close distance, 30 for radius, 50 for tibia
$!==========================================================================
$!
$! Script Breakdown - from original, August 15 version is only Step 1
$!
$! STEP 1: Cortical and Trabecular contours definition
$! STEP 2: Segmentation
$! STEP 3: ROI Creation (not yet implemented)
$!
$!==========================================================================
$!
$! STEP 1: Creates CORT_MASK.GOBJ & TRAB_MASK.GOBJ, and their related AIMs
$!
$!==========================================================================
$!
$ IPL_BATCH

! this is the greyscale aim from the ISQ cropped with a rough gobj around the single bone
/read
  -name                      aim
  -filename                  "IPL_AIM

! Periosteal gobj created in the periosteal contouring script
/gobj_to_aim
  -gobj_filename             "IPL_GOBJ0
  -output                    all
  -peel_iter                 0

!! Gross Gaussian Segmentation to eliminate thin (mostly Trabeculae) structures
!! Threshold is high to start with cortical structure as a basis
/gobj_maskaimpeel_ow
  -input_output              aim
  -gobj_filename             "IPL_GOBJ0
  -peel_iter                 0

/bounding_box_cut
  -input                     aim
  -output                    aim_bbc
  -z_only                    false
  -border                    10 10 20		! changed from 0to stop dilation errors

! CHANGE: delete
/delete aim

/seg_gauss
  -input                     aim_bbc
  -output                    cort
  -sigma                     "IPL_SIGMA0		!originally 2
  -support                   "IPL_SUPPORT0		!originally 3
  -lower_in_perm_aut_al      "IPL_LOWER0		!originally 500
  -upper_in_perm_aut_al      3000
  -value_in_range            127
  -unit                      2

! CHANGE: delete
/delete aim_bbc

! This needs to be written as a temp file to be used as a gobj_maskaim peel ow file later
! to ensure minimum cortical thickness. Peel iter does not work since it's only applied in 2D
/erosion
  -input                     all
  -output                    all_eroded
  -erode_distance            "IPL_PEEL0
  -use_previous_margin       false
  -metric                    11

/togobj_from_aim
  -input                     all_eroded
  -gobj_filename             "IPL_FNAME7
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1

! CHANGE: delete
/delete all_eroded

!! Start working on the trabecular mask
/subtract_aims
  -input1                    all		!peri mask
  -input2                    cort
  -output                    trab

!! Substracting a minimium cortical thickness (i.e. 6 vox)
/gobj_maskaimpeel_ow
  -input_output              trab
  -gobj_filename             "IPL_FNAME7
  -peel_iter                 0        		!orig 6, now using temp file that's already eroded by 6

/set trab 0 127

/rename trab cort

/cl_ow_rank_extract
  -input_output              cort
  -first_rank                1
  -last_rank                 1
  -connect_boundary          false
  -value_in_range            127

!! Removing the background
/gobj_maskaimpeel_ow
  -input_output              cort
  -gobj_filename             "IPL_GOBJ0
  -peel_iter                 0
!! End of Cortex_mask_off alternative sequence  !!

/subtract_aims
  -input1                    all
  -input2                    cort
  -output                    trab

/cl_ow_rank_extract
  -input_output              trab
  -first_rank                1
  -last_rank                 1
  -connect_boundary          false
  -value_in_range            127

!! 1st Erosion-CL-Dilation loop on Trab mask to remove Tb speckles not connected to it.
!! Erosion distance is small
/erosion
  -input                     trab
  -output                    trab_ero
  -erode_distance            3			!orig 3; run 4 this was at 4
  -use_previous_margin       false
  -metric                    11

/cl_ow_rank_extract
  -input_output              trab_ero
  -first_rank                1
  -last_rank                 1
  -connect_boundary          false
  -value_in_range            127

/dilation
  -input                     trab_ero
  -output                    trab_dil
  -dilate_distance           6			!orig at 3, same as erosion
  -continuous_at_boundary    0 0 0
  -use_previous_margin       false
  -metric                    11

!!==========================================================================
!!
!! Large close/open sequence to smooth the contour
!!
!!==========================================================================
/close
  -input                     trab_dil
  -output                    trab_close
  -close_distance            15
  -continuous_at_boundary    0 0 0
  -metric                    11

/delete trab_dil
/delete trab_dil_big
/delete trab_ero
/delete cort

!! GOBJ and BOUND used to bring back the white box size to its original size
!! (size is being modified by the dilation processes)
/gobj_maskaimpeel_ow
  -input_output              trab_close
  -gobj_filename             "IPL_GOBJ0
  -peel_iter                 "IPL_PEEL0

/bounding_box_cut
  -input                     trab_close
  -output                    trab_bbc
  -z_only                    false
  -border                    20 20 20	!changed from 0

/delete trab_close

! uncomment for debugging
!/write  trab_bbc    		"ipl_fname_s2aim

! uncomment for debugging
!/togobj_from_aim
!  -input                     trab_bbc
!  -gobj_filename             "ipl_fname_s2
!  -min_elements              0
!  -max_elements              0
!  -curvature_smooth          1

!! At this stage, the Tb mask (trab_bbc) is close to be correct,
!! but sharp corners may have been deleted by the open process
!! From here we put those corners back.
/open
  -input                     trab_bbc
  -output                    trab_open
  -open_distance             15
  -metric                    11

/set trab_open 127 0

/set trab_bbc 127 0

/subtract_aims
  -input1                    trab_bbc
  -input2                    trab_open
  -output                    corners

! CHANGE: delete
/delete trab_bbc

!! 2nd Erosion-CL-Dilation loop, for removing small speckles and keep only sharp corner lost by the
!! large open process
/erosion
  -input                     corners
  -output                    corn_ero
  -erode_distance            3
  -use_previous_margin       false
  -metric                    11


! min_number = ¨ipl_misc1_0
! ipl_misc1_0 = 800 for the Radius
! ipl_misc1_0 = 200000 for the Tibia
/cl_nr_extract
  -input                     corn_ero
  -output                    corn_cl
  -min_number                "IPL_MISC1_0
  -max_number                0
  -value_in_range            127

/delete corn_ero

/dilation
  -input                     corn_cl
  -output                    corn_major
  -dilate_distance           3
  -continuous_at_boundary    0 0 0
  -use_previous_margin       false
  -metric                    11

! min_number = ¨ipl_misc1_0
! ipl_misc1_0 = 800 for the Radius
! ipl_misc1_0 = 200000 for the Tibia
/cl_nr_extract
  -input                     corn_major
  -output                    corn_cl
  -min_number                "IPL_MISC1_0
  -max_number                0
  -value_in_range            127

/delete corn_major

! CL to handle case where dilation of null AIM give a full AIM
/cl_nr_extract
  -input                     corn_cl
  -output                    corners
  -min_number                1
  -max_number                500000
  -value_in_range            127

/delete corn_cl

!! Adding the corners to the almost finished trabecular mask
/add_aims
  -input1                    corners
  -input2                    trab_open
  -output                    trab

! CHANGE: delete
/delete trab_open
! CHANGE: delete
/delete corners

! ipl_misc1_1 = 30 for the Radius
! ipl_misc1_1 = 50 for the Tibia
/close
  -input                     trab
  -output                    trab_close
  -close_distance            "IPL_MISC1_1
  -continuous_at_boundary    0 0 0
  -metric                    11

!!==========================================================================
!!
!! Heuristics: Ensuring a minimum cortical thickness, slicewise CL clean up
!! and fixing up the box size of the contour object
!!
!!==========================================================================

/gobj_maskaimpeel_ow
  -input_output              trab_close
  -gobj_filename             "IPL_FNAME7
  -peel_iter                 0            !orig 6, now using temp file that's already eroded by 6

/cl_slicewise_extractow
  -input_output              trab_close
  -lo_vol_fract_in_perc      1				!modified to 10 from 50 due to 2nd condyle
  -up_vol_fract_in_perc      100
  -value_in_range            127

/subtract_aims
  -input1                    all
  -input2                    trab_close
  -output                    cort

/delete trab_close

/cl_slicewise_extractow
  -input_output              cort
  -lo_vol_fract_in_perc      1 			!modified to 1 from 50 due to 2nd condyle
  -up_vol_fract_in_perc      100
  -value_in_range            127

/bounding_box_cut
  -input                     cort
  -output                    cort_final
  -z_only                    false
  -border                    0 0 0
!! End of Heuristics

! CHANGE: delete
/delete cort

!!==========================================================================
!!
!! Heuristics: TRAB_MASK.GOBJ = GOBJ - CORT_MASK.GOBJ   !!!
!!
!!==========================================================================
/subtract_aims
  -input1                    all
  -input2                    cort_final
  -output                    trab_final

!! Creating the final contours
! Cortex will be created in step 3 after manual correction of trab.
!/togobj_from_aim
!  -input                     cort_final
!  -gobj_filename             "IPL_FNAME0
!  -min_elements              0
!  -max_elements              0
!  -curvature_smooth          1

/togobj_from_aim
  -input                     trab_final
  -gobj_filename             "IPL_FNAME1
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1
..
$
$ DELETE 'IPL_FNAME7'; ! Delete temp gobj
$ EXIT
