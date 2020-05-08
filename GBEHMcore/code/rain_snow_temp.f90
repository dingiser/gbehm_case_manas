
 SUBROUTINE rain_snow_temp (forc_t,forc_q,forc_rh,forc_psrf,forc_prc,forc_prl,tcrit,elev,&
                            prc_rain,prc_snow,prl_rain,prl_snow,t_precip,bifall)

!=======================================================================
! Determine precipitation type and precipitation water temperature
!
!---Original author: Yongjiu Dai, 09/1999; 08/31/2002, 04/2014
!---Major modified : Hongyi Li, 11/2016
!                    use Ding's (2014) scheme depending on elevation 
!                    and relative humidity. Wet-bulb temperature (not
!                    air temperature) are used as critical temperature 
!                    to determine precipitation type. 
!=======================================================================
!
  use precision
  use PhysicalConstants, only : tfrz
 
  IMPLICIT NONE
 
! ------------------------ Dummy Argument ------------------------------

  real(r8), INTENT(in) :: forc_t     ! temperature at agcm reference height [kelvin]
  real(r8), INTENT(in) :: forc_q     ! specific humidity at agcm reference height [kg/kg]
  real(r8), INTENT(in) :: forc_rh    ! relative humidity at agcm reference height [0-1]
  real(r8), INTENT(in) :: forc_psrf  ! atmosphere pressure at the surface [pa]
  real(r8), INTENT(in) :: forc_prc   ! convective precipitation [mm/s]
  real(r8), INTENT(in) :: forc_prl   ! large scale precipitation [mm/s]

  real(r8), INTENT(in) :: tcrit      ! critical temp. to determine rain or snow
  real(r8), INTENT(in) :: elev       ! elevation [m]

  real(r8), INTENT(out) :: prc_rain  ! convective rainfall [kg/(m2 s)]
  real(r8), INTENT(out) :: prc_snow  ! convective snowfall [kg/(m2 s)]
  real(r8), INTENT(out) :: prl_rain  ! large scale rainfall [kg/(m2 s)]
  real(r8), INTENT(out) :: prl_snow  ! large scale snowfall [kg/(m2 s)]
  real(r8), INTENT(out) :: t_precip  ! snowfall/rainfall temperature [kelvin]
  real(r8), INTENT(out) :: bifall    ! bulk density of newly fallen dry snow [kg/m3]
  real(r8) :: flfall                 ! fraction of liquid water within falling precip.
  real(r8) :: t_wetbulb              ! wet-bulb temperature
  real(r8) :: tc,tmin,tmax,deltat,deltas,c1,c2,c3,c4,c5
  integer, parameter :: use_colm_division_method = 0

  call wetbulb(forc_t,forc_psrf,forc_q,t_wetbulb)

  if (use_colm_division_method == 0) then 
    !-----------------------------------------------------------------------
    !--The following codes is from the Dai's COLM scheme.
    !   the upper limit of air temperature is set for snowfall, this cut-off 
    !   was selected based on Fig. 1, Plate 3-1, of Snow Hydrology (1956).
    !   the percentage of liquid water by mass, which is arbitrarily set to 
    !   vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.

      if(forc_t>tfrz+tcrit)then
        prc_rain = forc_prc       ! convective rainfall (mm/s)
        prl_rain = forc_prl       ! large scale rainfall (mm/s)
        prc_snow = 0.             ! convective snowfall (mm/s)
        prl_snow = 0.             ! large scale snowfall (mm/s)
        flfall = 1.               ! fraction of liquid water within falling precip.
        bifall = 169.             ! (not used)
      else
        if(forc_t<=tfrz)then
          flfall = 0.
        else if(forc_t>tfrz .and. forc_t<=tfrz+tcrit)then
          flfall = max(0.,-54.632+0.2*forc_t)
          flfall = min(0.4,flfall)
        else
          flfall = 0.4
        endif 
        prc_rain = forc_prc*flfall
        prl_rain = forc_prl*flfall
        prc_snow = forc_prc*(1.-flfall)                 
        prl_snow = forc_prl*(1.-flfall) 
      endif
  else
    !-----------------------------------------------------------------------
    !--Coded by HY., 2016-11-14 PM.
    !  Use Ding BH's(2014) scheme.
    !  Ding, Baohong, et al. Journal of Hydrology 513 (2014): 154-163.
    !  --"The dependence of precipitation types on surface elevation and meteorological conditions and its parameterization." 

      c1 = -9.614; c2 = 16.06; c3 = 0.0885; c4 = -0.1042; c5 = -5.87
      tc = min(2.5,c1*forc_rh**2+c2*forc_rh+c3*(elev*0.001)**2+c4*(elev*0.001)+c5)
      if (forc_rh<=0.78) then 
        if(t_wetbulb>tfrz+tc)then
          prc_rain = forc_prc       ! convective rainfall (mm/s)
          prl_rain = forc_prl       ! large scale rainfall (mm/s)
          prc_snow = 0.             ! convective snowfall (mm/s)
          prl_snow = 0.             ! large scale snowfall (mm/s)
          flfall = 1.               ! fraction of liquid water within falling precip.
        else
          prc_rain = 0.             ! convective rainfall (mm/s)
          prl_rain = 0.             ! large scale rainfall (mm/s)
          prc_snow = forc_prc       ! convective snowfall (mm/s)
          prl_snow = forc_prl       ! large scale snowfall (mm/s)
          flfall = 0.               ! fraction of liquid water within falling precip.
        endif
      else ! forc_rh>0.78
        deltas = max(2.374-1.634*forc_rh,0.74)
        deltat = 0.215-0.099*forc_rh + 1.018*forc_rh*forc_rh

        tmin = min(0.,tc - deltas*log(exp(deltat/deltas)-2.*exp(deltat/deltas)))
        tmax = min(4.,2.*tc - tmin)

        if (tmax<tmin) then
          print*,'debug critical precipitation temperature: tmax,tmin', tmax, tmin
        endif 

        if(t_wetbulb>tfrz+tmax)then
          prc_rain = forc_prc       ! convective rainfall (mm/s)
          prl_rain = forc_prl       ! large scale rainfall (mm/s)
          prc_snow = 0.             ! convective snowfall (mm/s)
          prl_snow = 0.             ! large scale snowfall (mm/s)
          flfall = 1.               ! fraction of liquid water within falling precip.
        else if(t_wetbulb<tfrz+tmin) then
          prc_rain = 0.             ! convective rainfall (mm/s)
          prl_rain = 0.             ! large scale rainfall (mm/s)
          prc_snow = forc_prc       ! convective snowfall (mm/s)
          prl_snow = forc_prl       ! large scale snowfall (mm/s)
          flfall = 0.               ! fraction of liquid water within falling precip.
        else 
          flfall = 0.4
          prc_rain = forc_prc*flfall
          prl_rain = forc_prl*flfall
          prc_snow = forc_prc*(1.-flfall)                 
          prl_snow = forc_prl*(1.-flfall) 
        endif 
      endif
  endif

    ! -------------------------------------------------------------
    ! calculate bulk density of newly fallen dry snow [kg/m3]
    ! -------------------------------------------------------------     
    ! use Alta relationship, Anderson(1976); LaChapelle(1961), 
    ! U.S.Department of Agriculture Forest Service, Project F, 
    ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

    if(forc_t>tfrz+2.)then
      bifall = 169.
    else if(forc_t>tfrz-15.)then
      bifall = 50.+1.7*(forc_t-tfrz+15.)**1.5
    else
      bifall = 50.
    endif

    ! -------------------------------------------------------------
    ! calculate temperature of rainfall or snowfall
    ! -------------------------------------------------------------
    ! this scheme is from the Original colm, I dont know it is right or not

    t_precip = t_wetbulb
    if (forc_t > 275.65) then
       if (t_precip < tfrz) t_precip = tfrz
    else
       t_precip = min(tfrz,t_precip)
       if(flfall > 0.0)then
         t_precip = tfrz - sqrt((1.0/flfall)-1.0)/100.0
       endif
    endif

 END SUBROUTINE rain_snow_temp
