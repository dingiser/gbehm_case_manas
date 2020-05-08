SUBROUTINE format_clmforcing(forc_t        , &
                             forc_q        , &
                             forc_psrf     , &
                             forc_solarin  , &
                             forc_prec     , &
                             forc_pco2m    , &
                             forc_po2m     , &
                             forc_prc      , &
                             forc_prl      , &
                             forc_pbot     , &
                             forc_sols     , &
                             forc_soll     , &
                             forc_solsd    , &
                             forc_solld    , &
                             forc_hgt_u    , &
                             forc_hgt_t    , &
                             forc_hgt_q    , &
                             forc_rhoair     )

! ====================================================================================
! Modified by LIHONGYI, 2016-5-16
! ------------------------------------------
! Format special forcing data needed in CoLM module. Quality control should be 
! finished before this routine.
! The following data is additionally needed:
!
! forc_pco2m : partial pressure of CO2 at observational height [pa]
! forc_po2m  : partial pressure of O2 at observational height [pa]
! forc_prc   : convective precipitation [mm/s]
! forc_prl   : large scale precipitation [mm/s]
! forc_pbot  : atmosphere pressure at the bottom of the atmos.model level [pa]
! forc_sols  : atm vis direct beam solar rad onto srf [W/m2]
! forc_soll  : atm nir direct beam solar rad onto srf [W/m2]
! forc_solsd : atm vis diffuse solar rad onto srf [W/m2]
! forc_solld : atm nir diffuse solar rad onto srf [W/m2]
! forc_hgt_u : observational height of wind [m]
! forc_hgt_t : observational height of temperature [m]
! forc_hgt_q : observational height of humidity [m]
! forc_rhoair: density air [kg/m3]
! ====================================================================================

  use precision
  use PhysicalConstants, only: rgas
  ! use MOD_TimeInvariants
  ! use timemanager

  IMPLICIT NONE
  real(r8), INTENT(in) :: &
      forc_t       , &! temperature at agcm reference height [kelvin]
      forc_q       , &! specific humidity at agcm reference height [kg/kg]
      forc_psrf    , &! atmosphere pressure at the surface [pa]
      forc_solarin , &! solar rad onto srf [W/m2]
      forc_prec       ! precipitation [mm/s]

! Forcing
! ----------------------
  real(r8), INTENT(out) :: &
      forc_pco2m   , &! partial pressure of CO2 at observational height [pa]
      forc_po2m    , &! partial pressure of O2 at observational height [pa]
      forc_prc     , &! convective precipitation [mm/s]
      forc_prl     , &! large scale precipitation [mm/s]
      forc_pbot    , &! atmosphere pressure at the bottom of the atmos.model level [pa]
      forc_sols    , &! atm vis direct beam solar rad onto srf [W/m2]
      forc_soll    , &! atm nir direct beam solar rad onto srf [W/m2]
      forc_solsd   , &! atm vis diffuse solar rad onto srf [W/m2]
      forc_solld   , &! atm nir diffuse solar rad onto srf [W/m2]
      forc_hgt_u   , &! observational height of wind [m]
      forc_hgt_t   , &! observational height of temperature [m]
      forc_hgt_q   , &! observational height of humidity [m]
      forc_rhoair     ! density air [kg/m3]

! local variables
!----------------------
      ! real(r8) :: coszen  ! cosine of solar zenith angle
      ! real(r8) :: calday    ! Julian cal day (1.xx to 365.xx)
      ! real(r8) :: sunang, cloud, difrat, vnrat
      real(r8) :: hsolar, ratio_rvrf 
      !real(r8) :: orb_coszen

!------------------------------------------------------------

      forc_pco2m  = forc_psrf *398.03e-06 
      forc_po2m   = forc_psrf *0.209       
      forc_pbot   = forc_psrf 
      forc_prl    = forc_prec * 2/3. ! ?
      forc_prc    = forc_prec * 1/3. ! ?
      forc_hgt_u  = 10.0 !m
      forc_hgt_t  = 2.0
      forc_hgt_q  = 2.0      

! #if(defined USE_QIAN_DATA)
         !---------------------------------------------------------------
         ! NOTE: codes from CLM4.5-CESM1.2.0
         ! relationship between incoming NIR or VIS radiation and ratio of
         ! direct to diffuse radiation calculated based on one year's worth of
         ! hourly CAM output from CAM version cam3_5_55
         !---------------------------------------------------------------

      hsolar = forc_solarin*0.5_R8
    ! NIR (dir, diff)
      ratio_rvrf = min(0.99_R8,max(0.29548_R8 + 0.00504_R8*hsolar  &
         -1.4957e-05_R8*hsolar**2 + 1.4881e-08_R8*hsolar**3,0.01_R8))
      forc_soll  = ratio_rvrf*hsolar
      forc_solld = (1._R8 - ratio_rvrf)*hsolar
 
    ! VIS (dir, diff)
      ratio_rvrf = min(0.99_R8,max(0.17639_R8 + 0.00380_R8*hsolar  &
         -9.0039e-06_R8*hsolar**2 + 8.1351e-09_R8*hsolar**3,0.01_R8))
      forc_sols  = ratio_rvrf*hsolar
      forc_solsd = (1._R8 - ratio_rvrf)*hsolar
       
      forc_rhoair= (forc_pbot &
                    - 0.378*forc_q*forc_pbot/(0.622+0.378*forc_q))&
                    / (rgas*forc_t)  
      if (forc_rhoair >1.29 .or. forc_rhoair <0.1) then
        forc_rhoair = 1.29
      endif
         
END SUBROUTINE format_clmforcing