
SUBROUTINE soil_color_refl(L,soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb)
! ======================================================================
! Guess the soil color (reflectance) based on the land cover types
! 
! Created by Yongjiu Dai, 03/2014
! ======================================================================
use precision

IMPLICIT NONE
      integer, intent(in) :: L  ! land cover types (GLCC USGS/MODIS IGBP)
      real(r8), intent(out) :: soil_s_v_alb ! albedo of visible of the saturated soil 
      real(r8), intent(out) :: soil_d_v_alb ! albedo of visible of the dry soil 
      real(r8), intent(out) :: soil_s_n_alb ! albedo of near infrared of the saturated soil
      real(r8), intent(out) :: soil_d_n_alb ! albedo of near infrared of the dry soil

      integer :: isc             ! soil color

      real(r8) soil_s_v_refl(20) ! Saturated visible soil reflectance
      real(r8) soil_d_v_refl(20) ! Dry visible soil reflectance
      real(r8) soil_s_n_refl(20) ! Saturated near infrared soil reflectance
      real(r8) soil_d_n_refl(20) ! Dry near infrared soil reflectance

! ----------------------------------------------------------------------
! The soil color and reflectance is from the work: 
! Peter J. Lawrence and Thomas N. Chase, 2007: 
! Representing a MODIS consistent land surface in the Community Land Model (CLM 3.0): 
! Part 1 generating MODIS consistent land surface parameters

      soil_s_v_refl = (/ 0.26, 0.24, 0.22, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, &
                         0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04 /)

      soil_d_v_refl = (/ 0.37, 0.35, 0.33, 0.31, 0.30, 0.29, 0.28, 0.27, 0.26, 0.25, &
                         0.24, 0.23, 0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15 /)

      soil_s_n_refl = (/ 0.52, 0.48, 0.44, 0.40, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, &
                         0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08 /)

      soil_d_n_refl = (/ 0.63, 0.59, 0.55, 0.51, 0.49, 0.47, 0.45, 0.43, 0.41, 0.39, &
                         0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19 /)

! Guessed soil color of the Lawrence's classification

      if(L.eq. 0) isc = 1  ! 0  Ocean (not used)
      if(L.eq. 1) isc = 16 ! 1  Urban and Built-Up Land
      if(L.eq. 2) isc = 3  ! 2  Dryland Cropland and Pasture 
      if(L.eq. 3) isc = 9  ! 3  Irrigated Cropland and Pasture
      if(L.eq. 4) isc = 10 ! 4  Mixed Dryland/Irrigated Cropland and Pasture
      if(L.eq. 5) isc = 4  ! 5  Cropland/Grassland Mosaic
      if(L.eq. 6) isc = 6  ! 6  Cropland/Woodland Mosaic
      if(L.eq. 7) isc = 2  ! 7  Grassland
      if(L.eq. 8) isc = 8  ! 8  Shrubland
      if(L.eq. 9) isc = 7  ! 9  Mixed Shrubland/Grassland
      if(L.eq.10) isc = 5  !10  Savanna
      if(L.eq.11) isc = 19 !11  Deciduous Broadleaf Forest 
      if(L.eq.12) isc = 20 !12  Deciduous Needleleaf Forest 
      if(L.eq.13) isc = 18 !13  Evergreen Broadleaf Forest  
      if(L.eq.14) isc = 17 !14  Evergreen Needleleaf Forest 
      if(L.eq.15) isc = 16 !15  Mixed Forest
      if(L.eq.16) isc = 1  !16  Water Bodies (not used)
      if(L.eq.17) isc = 15 !17  Herbaceous Wetland 
      if(L.eq.18) isc = 14 !18  Wooded Wetland
      if(L.eq.19) isc = 1  !19  Barren or Sparsely Vegetated
      if(L.eq.20) isc = 12 !20  Herbaceous Tundra
      if(L.eq.21) isc = 12 !21  Wooded Tundra
      if(L.eq.22) isc = 13 !22  Mixed Tundra
      if(L.eq.23) isc = 11 !23  Bare Ground Tundra  
      if(L.eq.24) isc = 1  !24  Snow or Ice (not used)


      soil_s_v_alb = soil_s_v_refl(isc)
      soil_d_v_alb = soil_d_v_refl(isc)
      soil_s_n_alb = soil_s_n_refl(isc)
      soil_d_n_alb = soil_d_n_refl(isc)

END SUBROUTINE soil_color_refl
! ----------------------------------------------------------------------
! EOP
