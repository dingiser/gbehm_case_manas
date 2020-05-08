SUBROUTINE regular_forcing(idate,lon_points,lat_points, &
                           forc_xy_t        , &
                           forc_xy_q        , &
                           forc_xy_psrf     , &
                           forc_xy_solarin  , &
                           forc_xy_frl      , &
                           forc_xy_prec     , &
                           forc_xy_us       , &
                           forc_xy_vs       , &
                           forc_pco2m       , &
                           forc_po2m        , &
                           forc_us          , &
                           forc_vs          , &
                           forc_t           , &
                           forc_q           , &
                           forc_rh          , &
                           forc_prc         , &
                           forc_prl         , &
                           forc_psrf        , &
                           forc_pbot        , &
                           forc_sols        , &
                           forc_soll        , &
                           forc_solsd       , &
                           forc_solld       , &
                           forc_frl         , &
                           forc_hgt_u       , &
                           forc_hgt_t       , &
                           forc_hgt_q       , &
                           forc_rhoair      )

  use precision
  use PhysicalConstants, only: rgas, grav
  ! use MOD_TimeInvariants
  use timemanager

  IMPLICIT NONE
  integer,  INTENT(in) :: idate(3)
  integer,  INTENT(in) :: lon_points
  integer,  INTENT(in) :: lat_points
  ! integer,  INTENT(in) :: lons(lon_points) 
  
  real(r8), INTENT(in) :: &
      forc_xy_t       (lon_points,lat_points), &! temperature at agcm reference height [kelvin]
      forc_xy_q       (lon_points,lat_points), &! specific humidity at agcm reference height [kg/kg]
      forc_xy_psrf    (lon_points,lat_points), &! atmosphere pressure at the surface [pa]
      forc_xy_solarin (lon_points,lat_points), &! solar rad onto srf [W/m2]
      forc_xy_frl     (lon_points,lat_points), &! atmospheric infrared (longwave) radiation [W/m2]
      forc_xy_prec    (lon_points,lat_points), &! precipitation [mm/s]
      forc_xy_us      (lon_points,lat_points), &! wind speed in eastward direction [m/s]
      forc_xy_vs      (lon_points,lat_points)   ! wind speed in northward direction [m/s]

! Forcing
! ----------------------
  real(r8), INTENT(out) :: &
      forc_pco2m (lon_points,lat_points)  , &! partial pressure of CO2 at observational height [pa]
      forc_po2m  (lon_points,lat_points)  , &! partial pressure of O2 at observational height [pa]
      forc_us    (lon_points,lat_points)  , &! wind speed in eastward direction [m/s]
      forc_vs    (lon_points,lat_points)  , &! wind speed in northward direction [m/s]
      forc_t     (lon_points,lat_points)  , &! temperature at agcm reference height [kelvin]
      forc_q     (lon_points,lat_points)  , &! specific humidity at agcm reference height [kg/kg]
      forc_rh    (lon_points,lat_points)  , &! relative humidity[0-1]
      forc_prc   (lon_points,lat_points)  , &! convective precipitation [mm/s]
      forc_prl   (lon_points,lat_points)  , &! large scale precipitation [mm/s]
      forc_psrf  (lon_points,lat_points)  , &! atmosphere pressure at the surface [pa]
      forc_pbot  (lon_points,lat_points)  , &! atmosphere pressure at the bottom of the atmos.model level [pa]
      forc_sols  (lon_points,lat_points)  , &! atm vis direct beam solar rad onto srf [W/m2]
      forc_soll  (lon_points,lat_points)  , &! atm nir direct beam solar rad onto srf [W/m2]
      forc_solsd (lon_points,lat_points)  , &! atm vis diffuse solar rad onto srf [W/m2]
      forc_solld (lon_points,lat_points)  , &! atm nir diffuse solar rad onto srf [W/m2]
      forc_frl   (lon_points,lat_points)  , &! atmospheric infrared (longwave) radiation [W/m2]
      forc_hgt_u (lon_points,lat_points)  , &! observational height of wind [m]
      forc_hgt_t (lon_points,lat_points)  , &! observational height of temperature [m]
      forc_hgt_q (lon_points,lat_points)  , &! observational height of humidity [m]
      forc_rhoair(lon_points,lat_points)     ! density air [kg/m3]
! ----------------F2PY definition -----------------

! local variables
!----------------------
      integer  :: i, j
      ! real(r8) :: coszen  ! cosine of solar zenith angle
      real(r8) :: calday    ! Julian cal day (1.xx to 365.xx)
      real(r8) :: sunang, cloud, difrat, vnrat
      real(r8) :: a, hsolar, ratio_rvrf 

      !real(r8) :: orb_coszen

!------------------------------------------------------------

      forc_pco2m (:,:) = forc_xy_psrf (:,:)*398.03e-06 
      forc_po2m  (:,:) = forc_xy_psrf (:,:)*0.209       
      forc_us    (:,:) = forc_xy_us   (:,:) 
      forc_frl   (:,:) = forc_xy_frl  (:,:)      
      forc_vs    (:,:) = forc_xy_vs   (:,:)
      forc_t     (:,:) = forc_xy_t    (:,:)
      forc_q     (:,:) = forc_xy_q    (:,:)
      forc_psrf  (:,:) = forc_xy_psrf (:,:)
      forc_pbot  (:,:) = forc_xy_psrf (:,:)
      forc_prl   (:,:) = forc_xy_prec * 2/3. ! ?
      forc_prc   (:,:) = forc_xy_prec * 1/3. ! ?
      forc_hgt_u (:,:) = 10.0 !m
      forc_hgt_t (:,:) = 2.0
      forc_hgt_q (:,:) = 2.0
      
      !calday = calendarday(idate, lons(1)) 

! #if(defined USE_QIAN_DATA)
         !---------------------------------------------------------------
         ! NOTE: codes from CLM4.5-CESM1.2.0
         ! relationship between incoming NIR or VIS radiation and ratio of
         ! direct to diffuse radiation calculated based on one year's worth of
         ! hourly CAM output from CAM version cam3_5_55
         !---------------------------------------------------------------

         do j = 1, lat_points
            do i = 1, lon_points
               hsolar = forc_xy_solarin(i,j)*0.5_R8
             ! NIR (dir, diff)
               ratio_rvrf = min(0.99_R8,max(0.29548_R8 + 0.00504_R8*hsolar  &
                  -1.4957e-05_R8*hsolar**2 + 1.4881e-08_R8*hsolar**3,0.01_R8))
               forc_soll (i,j) = ratio_rvrf*hsolar
               forc_solld(i,j) = (1._R8 - ratio_rvrf)*hsolar

             ! VIS (dir, diff)
               ratio_rvrf = min(0.99_R8,max(0.17639_R8 + 0.00380_R8*hsolar  &
                  -9.0039e-06_R8*hsolar**2 + 8.1351e-09_R8*hsolar**3,0.01_R8))
               forc_sols (i,j) = ratio_rvrf*hsolar
               forc_solsd(i,j) = (1._R8 - ratio_rvrf)*hsolar

         if(forc_t(i,j) <= 0.) forc_t(i,j) = 270.      
! The standard measuring conditions for temperature are two meters above the ground
! Scientists have measured the most frigid temperature ever 
! recorded on the continent's eastern highlands: about (180K) colder than dry ice.
         if(forc_t(i,j) < 180.) forc_t(i,j) = 180.
! the highest air temp was found in Kuwait 326 K, Sulaibya 2012-07-31; 
! Pakistan, Sindh 2010-05-26; Iraq, Nasiriyah 2011-08-03
         if(forc_t(i,j) > 326.) forc_t(i,j) = 326.
         
         forc_rhoair(i,j) = (forc_pbot(i,j) &
                         - 0.378*forc_q(i,j)*forc_pbot(i,j)/(0.622+0.378*forc_q(i,j)))&
                         / (rgas*forc_t(i,j))  
         if (forc_rhoair(i,j) >1.29 .or. forc_rhoair(i,j) <0.1) then
          ! print*,'error forc_rhoair=',forc_rhoair(i,j),forc_pbot(i,j),forc_q(i,j),forc_t(i,j) 
          forc_rhoair(i,j) = 1.29
         endif
! some assumption to avoid the NAN value.
         if ( isnan(forc_pco2m (i,j)))  forc_pco2m (i,j) = 101325.*398.03e-06 ! [pa]
         if ( isnan(forc_po2m  (i,j)))  forc_po2m  (i,j) = 101325.*0.209      ! [pa]         
         if ( isnan(forc_us    (i,j)).OR. abs(forc_us(i,j))>100.)  forc_us(i,j) = 2.! wind speed in eastward direction [m/s]         
         if ( isnan(forc_vs    (i,j)).OR. abs(forc_vs(i,j))>100.)  forc_vs(i,j) = 2.! wind speed in northward direction [m/s]         
         if ( isnan(forc_t     (i,j)))  forc_t     (i,j) = 273.! temperature at agcm reference height [kelvin]         
         if ( isnan(forc_q     (i,j)) .or. (forc_q     (i,j))<=0.)  forc_q     (i,j) = 0.001! specific humidity at agcm reference height [kg/kg]         
         if ( isnan(forc_prc   (i,j)) .OR. forc_prc(i,j)<0. .OR. forc_prc(i,j)>100.)  forc_prc   (i,j) = 0.! convective precipitation [mm/s]         
         if ( isnan(forc_prl   (i,j)) .OR. forc_prl(i,j)<0. .OR. forc_prl(i,j)>100.)  forc_prl   (i,j) = 0.! large scale precipitation [mm/s]         
         if ( isnan(forc_psrf  (i,j)) .OR. forc_psrf  (i,j)<300. .OR. forc_psrf  (i,j)>120000.)  forc_psrf  (i,j) = 101325.![pa]         
         if ( isnan(forc_pbot  (i,j)) .OR. forc_pbot  (i,j)<300. .OR. forc_pbot  (i,j)>120000.)  forc_pbot  (i,j) = 101325.! [pa]         
         if ( isnan(forc_sols  (i,j)))  forc_sols  (i,j) = 100.! atm vis direct beam solar rad onto srf [W/m2]         
         if ( isnan(forc_soll  (i,j)))  forc_soll  (i,j) = 100.! atm nir direct beam solar rad onto srf [W/m2]         
         if ( isnan(forc_solsd (i,j)))  forc_solsd (i,j) = 50.! atm vis diffuse solar rad onto srf [W/m2]         
         if ( isnan(forc_solld (i,j)))  forc_solld (i,j) = 50.! atm nir diffuse solar rad onto srf [W/m2]         
         if ( isnan(forc_frl   (i,j)))  forc_frl   (i,j) = 100.! atmospheric infrared (longwave) radiation [W/m2]         
         if ( isnan(forc_hgt_u (i,j)))  forc_hgt_u (i,j) = 10.! observational height of wind [m]         
         if ( isnan(forc_hgt_t (i,j)))  forc_hgt_t (i,j) = 2.! observational height of temperature [m]         
         if ( isnan(forc_hgt_q (i,j)))  forc_hgt_q (i,j) = 2.! observational height of humidity [m]         
         if ( isnan(forc_rhoair(i,j)))  forc_rhoair(i,j) = 1.29! density air [kg/m3]         

         call rhtoe(forc_q(i,j),forc_t(i,j),forc_psrf(i,j),forc_rh(i,j))               
         enddo
      enddo
         
END SUBROUTINE regular_forcing

  subroutine rhtoe(qe,airtk,pres,rh)
    use precision
    IMPLICIT NONE
    real(r8), intent(in) :: airtk,&!K
                            pres   !Pa
    real(r8), intent(inout) :: qe  !kg/kg
                            
    real(r8), intent(out):: rh     !relative humidity[0-1]
    real(r8) es,airtc
    airtc = airtk-273.15
    es = 611.*exp(17.27*airtc/(237.3+airtc))
    rh = qe*pres/(0.622*es)
    if (rh>1.) then
        ! print*,'running error in calculating RH',airtc,rh,pres,qe,es
        rh = 0.95
        qe = 0.622*es*rh/pres
    endif
    
    if (rh<0.) then
        ! print*,'running error in calculating RH',airtc,rh,pres,qe,es
        rh = 0.05
        qe = 0.622*es*rh/pres
    endif    
  end subroutine rhtoe