!#include <define.h>

SUBROUTINE CLMMAIN (icol,irow,deltim, oro, &
           nl,nl_soil, nl_lake, maxsnl, &
           lon, dlon, dlat, ivt, itypwat, &
           lakedepth, dz_lake, & ! new lake scheme         
           soil_s_v_alb, soil_d_v_alb, soil_s_n_alb, soil_d_n_alb,  &! soil information and lake depth
           wsat,wrsd,watern,alpha,slope,length,elev,Ds,Dr,Dg,kground, &!GBHM parameters
           porsl,        psi0,         bsw,          hksati,        &
           wfld,         adjfac,&
           csol,         dksatu,       dkdry,        rootfr,        &         
           z0m,          displa,       sqrtdi,                      &! vegetation information
           effcon,       vmax25,       slti,         hlti,          &
           shti,         hhti,         trda,         trdm,          &
           trop,         gradm,        binter,       extkn,         &
           chil,         ref,          tran,                        &         
           forc_pco2m,   forc_po2m,    forc_us,      forc_vs,       &! atmospheric forcing
           forc_t,       forc_q,       forc_prc,     forc_prl,      &
           forc_rain,    forc_snow,    forc_psrf,    forc_pbot,     &
           forc_sols,    forc_soll,    forc_solsd,   forc_solld,    &
           forc_frl,     forc_hgt_u,   forc_hgt_t,   forc_hgt_q,    &
           forc_rhoair,  forc_rh,                                   &         
           idate,                                                   &
           z_soisno,     dz_soisno,    t_soisno,     wliq_soisno,   &
           wice_soisno,  sst,          ss_f,         smf_soisno,    &
           smf_ice,&
           t_grnd,       tlsun,        tlsha,        ldew,          &
           sag,          scv,          snowdp,       fveg,          &
           fsno,         sigf,         green,        lai,           &
           sai,          coszen,       albg,         albv,          &
           alb,          ssun,         ssha,         thermk,        &
           extkb,        extkd,                                     &
           zwt,          wa,                                        &
           t_lake,       lake_icefrac,                              &         
           laisun,       laisha,                                    &! additional diagnostic variables for output
           Drw,          qground,      qsmelt,       qsubl,         &
           qfros,        snoweva,&
           rstfac,       h2osoi,       wat,          qlat,          &         
           taux,         tauy,         fsena,        fevpa,         &! FLUXES
           lfevpa,       fsenl,        fevpl,        etr,           &
           fseng,        fevpg,        olrg,         fgrnd,         &
           trad,         tref,         qref,         rsur,          &
           rnof,         qintr,        qinfl,        qdrip,         &
           rst,          assim,        respc,        sabvsun,       &
           sabvsha,      sabg,         sr,           solvd,         &
           solvi,        solnd,        solni,        srvd,          &
           srvi,         srnd,         srni,         solvdln,       &
           solviln,      solndln,      solniln,      srvdln,        &
           srviln,       srndln,       srniln,       qcharge,       &
           xerr,         zerr,                                      &         
           zlnd,         zsno,         csoilc,       dewmx,         &! TUNABLE modle constants
           wtfact,       capr,         cnfac,        ssi,           &
           wimp,         pondmx,       smpmax,       smpmin,        &
           trsmx0,       tcrit,&
           rssca)!,                                     &           
           ! emis,         z0ma,         zol,          rib)!,           &! additional variables required by coupling with WRF model
           ! ustar,        qstar,        tstar,                       &
           ! fm,           fh,           fq)

!=======================================================================
!
! Main subroutine, advance time information
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002, 12/2012, 02/2014
!
!    FLOW DIAGRAM FOR CLMMAIN
!
!    CLMMAIN ===> netsolar                 |> all surface
!                 rain_snow_temp           !> all surface
!
!                 LEAF_interception        |]  
!                 newsnow                  |] itypwat = 0 (soil ground)
!                 THERMAL                  |]         = 1 (urban & built-up)
!                 WATER                    |]         = 2 (wetland)
!                 snowcompaction           |]         = 3 (land ice)
!                 snowlayerscombine        |]         = 4 (lake)
!                 snowlayersdivide         |]
!                 snowage                  |]  
!
!                 newsnow_lake             |]
!                 laketem                  |] lake scheme
!                 snowwater_lake           |]
!
!                 SOCEAN                   |> ocean and sea ice
!                
!                 orb_coszen_func          |> all surface
!                 EcoModel (LAI_empirical) |> land 
!                 snowfraction             |> land
!                 albland                  |> land 
!                 albocean                 |> ocean & sea ice
!
!=======================================================================

  use precision
  use PhysicalConstants, only : tfrz, denh2o, denice
  use SOIL_SNOW_hydrology
  use SNOW_Layers_CombineDivide
  use GLACIER
  use LAKE
  use SIMPLE_OCEAN
  use ALBEDO
  use MOD_THERMAL
  use timemanager
 
  implicit none
 
! ------------------------ Dummy Argument ------------------------------
  real(r8),INTENT(in) :: deltim  !seconds in a time step [second]
  ! logical, INTENT(in) :: doalb   !true if time for surface albedo calculation
  ! logical, INTENT(in) :: dolai   !true if time for leaf area index calculation
  ! logical, INTENT(in) :: dosst   !true to update sst/ice/snow before calculation

  integer, INTENT(in) :: & 
        nl         , &! total layers from soil surface to acquifer bottom
        nl_soil    , &! number of soil layers
        nl_lake    , &! number of lake layers
        maxsnl        ! maximum number of snow layers
        ! ipatch        ! index i ; changed by Hongyi 2015-9-5
  integer, INTENT(in) :: icol,irow
  real(r8), INTENT(in) :: &
        lon        , &
        dlon       , &! logitude in radians
        dlat          ! latitude in radians

  integer, INTENT(in) :: & 
        ivt        , &! land cover type of USGS classification or others
        itypwat       ! land water type (0=soil, 1=urban and built-up, 
                      ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)
! Parameters
! ----------------------
  real(r8), INTENT(in) :: &
        lakedepth       , &! lake depth (m)
        dz_lake(nl_lake)   ! lake layer thickness (m)

  real(r8), INTENT(in) :: &
        rssca              ! remotely sensed snow cover area (1: snow; 0: snow free; 2: no avialable)

  real(r8), INTENT(in) :: &        ! soil physical parameters and lake info
        soil_s_v_alb   , &! albedo of visible of the saturated soil
        soil_d_v_alb   , &! albedo of visible of the dry soil
        soil_s_n_alb   , &! albedo of near infrared of the saturated soil
        soil_d_n_alb   , &! albedo of near infrared of the dry soil
        porsl(nl) , &! fraction of soil that is voids [-]
        psi0(nl)  , &! minimum soil suction [mm]
        bsw(nl)   , &! clapp and hornbereger "b" parameter [-]
        wfld(nl)   , &! water field capacity
        hksati(nl), &! hydraulic conductivity at saturation [mm h2o/s]
        csol(nl)  , &! heat capacity of soil solids [J/(m3 K)]
        dksatu(nl), &! thermal conductivity of saturated soil [W/m-K]
        dkdry(nl) , &! thermal conductivity for dry soil  [J/(K s m)]
        rootfr(nl_soil)   ! fraction of roots in each soil layer

  real(r8), INTENT(in) :: &        ! vegetation static, dynamic, derived parameters
        z0m        , &! aerodynamic roughness length [m]
        displa     , &! displacement height [m]
        sqrtdi     , &! inverse sqrt of leaf dimension [m**-0.5]
        effcon     , &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25     , &! maximum carboxylation rate at 25 C at canopy top
        slti       , &! slope of low temperature inhibition function      [s3] 
        hlti       , &! 1/2 point of low temperature inhibition function  [s4]
        shti       , &! slope of high temperature inhibition function     [s1]
        hhti       , &! 1/2 point of high temperature inhibition function [s2]
        trda       , &! temperature coefficient in gs-a model             [s5]
        trdm       , &! temperature coefficient in gs-a model             [s6]
        trop       , &! temperature coefficient in gs-a model          
        gradm      , &! conductance-photosynthesis slope parameter
        binter     , &! conductance-photosynthesis intercep
        extkn      , &! coefficient of leaf nitrogen allocation
        chil       , &! leaf angle distribution factor
        ref(2,2)   , &! leaf reflectance (iw=iband, il=life and dead)
        tran(2,2)     ! leaf transmittance (iw=iband, il=life and dead)

  real(r8), INTENT(in) :: &        ! tunable parameters
        zlnd       , &!roughness length for soil [m]
        zsno       , &!roughness length for snow [m]
        csoilc     , &!drag coefficient for soil under canopy [-]
        dewmx      , &!maximum dew
        wtfact     , &!fraction of model area with high water table
        capr       , &!tuning factor to turn first layer T into surface T
        cnfac      , &!Crank Nicholson factor between 0 and 1
        ssi        , &!irreducible water saturation of snow
        wimp       , &!water impremeable if porosity less than wimp
        pondmx     , &!ponding depth (mm)
        smpmax     , &!wilting point potential in mm
        smpmin     , &!restriction for min of soil poten.  (mm)
        trsmx0     , &!max transpiration for moist soil+100% veg.  [mm/s]
        tcrit      , &!critical temp. to determine rain or snow
        adjfac        !an adjustment parameter to solve the calibration of surface evaporation

! Forcing
! ----------------------
  real(r8), INTENT(in) :: &
        forc_pco2m , &! partial pressure of CO2 at observational height [pa]
        forc_po2m  , &! partial pressure of O2 at observational height [pa]
        forc_us    , &! wind speed in eastward direction [m/s]
        forc_vs    , &! wind speed in northward direction [m/s]
        forc_t     , &! temperature at agcm reference height [kelvin]
        forc_q     , &! specific humidity at agcm reference height [kg/kg]
        forc_rh    , &! relative humidity at agcm reference height [0-1]
        forc_prc   , &! convective precipitation [mm/s]
        forc_prl   , &! large scale precipitation [mm/s]
        forc_psrf  , &! atmosphere pressure at the surface [pa]
        forc_pbot  , &! atmosphere pressure at the bottom of the atmos. model level [pa]
        forc_sols  , &! atm vis direct beam solar rad onto srf [W/m2]
        forc_soll  , &! atm nir direct beam solar rad onto srf [W/m2]
        forc_solsd , &! atm vis diffuse solar rad onto srf [W/m2]
        forc_solld , &! atm nir diffuse solar rad onto srf [W/m2]
        forc_frl   , &! atmospheric infrared (longwave) radiation [W/m2]
        forc_hgt_u , &! observational height of wind [m]
        forc_hgt_t , &! observational height of temperature [m]
        forc_hgt_q , &! observational height of humidity [m]
        forc_rhoair   ! density air [kg/m3]

! GBHM parameters
! ----------------------------------------------------------------------        
  real(r8), INTENT(in) :: wsat,wrsd,watern,alpha !GBHM parameter 
  real(r8), intent(in) :: & ! added by HONGYILI, 2015.12.13
              ss_f,&
              Ds, &       ! depth of topsoil(m)
              Dr, &       ! depth of river (m)
              Dg, &       ! depth of unconfined acquifer (m)
              length, &   ! length of hillslope (m)
              slope,  &   ! slope of hillslope (m)
              elev,   &   ! elevation (m)
              kground     ! GW hydraulic conductivity (m/s)  
! Variables required for restart run
! ----------------------------------------------------------------------
  integer, INTENT(in) :: &
        idate(3)      ! next time-step /year/julian day/second in a day/

  real(r8), INTENT(inout) :: oro  ! ocean(0)/seaice(2)/ flag
  real(r8), INTENT(inout) :: &
        z_soisno(maxsnl+1:nl)   , &! layer depth (m)
        dz_soisno(maxsnl+1:nl)  , &! layer thickness (m)
        t_soisno(maxsnl+1:nl)   , &! soil + snow layer temperature [K]
        wliq_soisno(maxsnl+1:nl), &! liquid water (kg/m2)
        wice_soisno(maxsnl+1:nl), &! ice lens (kg/m2)
        smf_soisno (maxsnl:nl_soil+1), &! snowmelt fraction in snow-soil (0-1).Added by HY,2016-11-26
        smf_ice    (1:nl_soil),&! snowmelt fraction in soil ice  (0-1).Added by HY,2016-11-30
        t_lake(nl_lake)       ,&! lake temperature (kelvin)
        lake_icefrac(nl_lake) ,&! lake mass fraction of lake layer that is frozen
        t_grnd     , &! ground surface temperature [k]
        tlsun      , &! sunlit leaf temperature [K]
        tlsha      , &! shaded leaf temperature [K]
        ldew       , &! depth of water on foliage [kg/m2/s]
        sag        , &! non dimensional snow age [-]
        scv        , &! snow mass (kg/m2)
        snowdp     , &! snow depth (m)
        sst        , &! surface water storage (mm)
        zwt        , &! the depth to water table [m]
        wa         , &! water storage in aquifer [mm]
        fveg       , &! fraction of vegetation cover
        fsno       , &! fractional snow cover
        sigf       , &! fraction of veg cover, excluding snow-covered veg [-]
        green      , &! greenness
        lai        , &! leaf area index
        sai        , &! stem area index
        coszen     , &! cosine of solar zenith angle
        albg(2,2)  , &! albedo, ground [-]
        albv(2,2)  , &! albedo, vegetation [-]
        alb(2,2)   , &! averaged albedo [-]
        ssun(2,2)  , &! sunlit canopy absorption for solar radiation
        ssha(2,2)  , &! shaded canopy absorption for solar radiation
        thermk     , &! canopy gap fraction for tir radiation
        extkb      , &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd         ! diffuse and scattered diffuse PAR extinction coefficient
  real(r8), INTENT(inout) :: & !Added for GBHM
        Drw           ! water depth in the river of the flow interval (m)
    
     
! additional diagnostic variables for output
  real(r8), INTENT(out) :: &
        laisun     , &! sunlit leaf area index
        laisha     , &! shaded leaf area index
        rstfac     , &! factor of soil water stress 
        wat        , &! total water storage
        h2osoi(nl)! volumetric soil water in layers [m3/m3]

! Fluxes
! ----------------------------------------------------------------------
  real(r8), INTENT(out) :: &
        taux       , &! wind stress: E-W [kg/m/s**2]
        tauy       , &! wind stress: N-S [kg/m/s**2]
        fsena      , &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa      , &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa     , &! latent heat flux from canopy height to atmosphere [W/2]
        fsenl      , &! ensible heat from leaves [W/m2]
        fevpl      , &! evaporation+transpiration from leaves [mm/s]
        etr        , &! transpiration rate [mm/s]
        fseng      , &! sensible heat flux from ground [W/m2]
        fevpg      , &! evaporation heat flux from ground [mm/s]
        olrg       , &! outgoing long-wave radiation from ground+canopy
        fgrnd      , &! ground heat flux [W/m2]
        xerr       , &! water balance error at current time-step [mm/s]
        zerr       , &! energy balnce errore at current time-step [W/m2]
        tref       , &! 2 m height air temperature [K]
        qref       , &! 2 m height air specific humidity
        trad       , &! radiative temperature [K]
        rsur       , &! surface runoff (mm h2o/s)
        rnof       , &! total runoff (mm h2o/s)
        qlat       , &! lateral runoff(positive = out of soil column) (mm H2O /s)
        qground    , &! ground runoff(positive = out of soil column) (mm H2O /s)
        qintr      , &! interception (mm h2o/s)
        qinfl      , &! inflitration (mm h2o/s)
        qdrip      , &! throughfall (mm h2o/s)
        qcharge    , &! groundwater recharge [mm/s]
        qsmelt     , &! total snowmelt flux (mm H2O /s)
        rst        , &! canopy stomatal resistance 
        assim      , &! canopy assimilation
        respc      , &! canopy respiration
        sabvsun    , &! solar absorbed by sunlit vegetation [W/m2]
        sabvsha    , &! solar absorbed by shaded vegetation [W/m2]
        sabg       , &! solar absorbed by ground  [W/m2]
        sr         , &! total reflected solar radiation (W/m2)
        solvd      , &! incident direct beam vis solar radiation (W/m2)
        solvi      , &! incident diffuse beam vis solar radiation (W/m2)
        solnd      , &! incident direct beam nir solar radiation (W/m2)
        solni      , &! incident diffuse beam nir solar radiation (W/m2)
        srvd       , &! reflected direct beam vis solar radiation (W/m2)
        srvi       , &! reflected diffuse beam vis solar radiation (W/m2)
        srnd       , &! reflected direct beam nir solar radiation (W/m2)
        srni       , &! reflected diffuse beam nir solar radiation (W/m2)
        solvdln    , &! incident direct beam vis solar radiation at local noon(W/m2)
        solviln    , &! incident diffuse beam vis solar radiation at local noon(W/m2)
        solndln    , &! incident direct beam nir solar radiation at local noon(W/m2)
        solniln    , &! incident diffuse beam nir solar radiation at local noon(W/m2)
        srvdln     , &! reflected direct beam vis solar radiation at local noon(W/m2)
        srviln     , &! reflected diffuse beam vis solar radiation at local noon(W/m2)
        srndln     , &! reflected direct beam nir solar radiation at local noon(W/m2)
        srniln     , &! reflected diffuse beam nir solar radiation at local noon(W/m2)
        qsubl      , &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros      , &! surface dew added to snow pack (mm h2o /s) [+]
        forc_rain  , &! rain [mm/s]
        forc_snow  , &! snow [mm/s]
        snoweva       ! total snow evaporsublimation flux (mm H2O /s)

        
   real(r8) :: &     !, INTENT(out)
        emis       , &! averaged bulk surface emissivity
        z0ma       , &! effective roughness [m]
        zol        , &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib        , &! bulk Richardson number in surface layer
        ustar      , &! u* in similarity theory [m/s]
        qstar      , &! q* in similarity theory [kg/kg]
        tstar      , &! t* in similarity theory [K]
        fm         , &! integral of profile function for momentum
        fh         , &! integral of profile function for heat
        fq            ! integral of profile function for moisture

! ----------------------- Local  Variables -----------------------------
   real(r8), parameter :: spval = -1.e36_r8  ! a special value for missing and filling use
   real(r8) :: &
        calday     , &! Julian cal day (1.xx to 365.xx)
        endwb      , &! water mass at the end of time step
        errore     , &! energy balnce errore (Wm-2)
        errorw     , &! water balnce errore (mm)
        fiold(maxsnl+1:nl), &! fraction of ice relative to the total water
        w_old      , &! liquid water mass of the column at the previous time step (mm)
        orb_coszen , &! cosine of the solar zenith angle
        sabvg      , &! solar absorbed by ground + vegetation [W/m2]
        parsun     , &! PAR by sunlit leaves [W/m2]
        parsha     , &! PAR by shaded leaves [W/m2]
        qseva      , &! ground surface evaporation rate (mm h2o/s)
        qsdew      , &! ground surface dew formation (mm h2o /s) [+]
        rootr(1:nl_soil)  , &! root resistance of a layer, all layers add to 1.0
        scvold     , &! snow cover for previous time step [mm]
        sm         , &! rate of snowmelt [kg/(m2 s)]
        ssw        , &! water volumetric content of soil surface layer [m3/m3]
        tssub(7)   , &! surface/sub-surface temperatures [K]
        tssea      , &! sea surface temperature [K]
        totwb      , &! water mass at the begining of time step
        wt         , &! fraction of vegetation buried (covered) by snow [-]
        zi_soisno(maxsnl:nl) ! interface level below a "z" level (m)

   real(r8) :: &
        prc_rain   , &! convective rainfall [kg/(m2 s)]
        prc_snow   , &! convective snowfall [kg/(m2 s)]
        prl_rain   , &! large scale rainfall [kg/(m2 s)]
        prl_snow   , &! large scale snowfall [kg/(m2 s)]
        t_precip   , &! snowfall/rainfall temperature [kelvin]
        bifall     , &! bulk density of newly fallen dry snow [kg/m3]
        pg_rain    , &! rainfall onto ground including canopy runoff [kg/(m2 s)]
        pg_snow       ! snowfall onto ground including canopy runoff [kg/(m2 s)]

  integer snl      , &! number of snow layers
        imelt(maxsnl+1:nl), &! flag for: melting=1, freezing=2, Nothing happended=0
        lb         , &! lower bound of arrays
        j             ! do looping index

  real(r8) :: a, aa 
  ! Some debug variables, by Hongyi
  real(r8) :: debug_prcsnow,debug_prcrain,debug_prlsnow,debug_prlrain
!======================================================================
!  [1] Solar absorbed by vegetation and ground
!      and precipitation information (rain/snow fall and precip temperature
!======================================================================

      call netsolar (idate,dlon,deltim,&
                     itypwat,sigf,albg,albv,alb,ssun,ssha,&
                     forc_sols,forc_soll,forc_solsd,forc_solld,&
                     parsun,parsha,sabvsun,sabvsha,sabg,sabvg,sr,&
                     solvd,solvi,solnd,solni,srvd,srvi,srnd,srni,&
                     solvdln,solviln,solndln,solniln,srvdln,srviln,srndln,srniln)

      call rain_snow_temp (forc_t,forc_q,forc_rh,forc_psrf,forc_prc,forc_prl,tcrit,elev,&
                           prc_rain,prc_snow,prl_rain,prl_snow,t_precip,bifall)

      ! use remotely sensed snow data to modify the precipitation devision. HYLi, 2016-8-11
      ! if (rssca<2.) then 
      !   if (rssca==1.) then
      !     prc_snow = prc_snow + prc_rain
      !     prl_snow = prl_snow + prl_rain
      !     prc_rain = 0.
      !     prl_rain = 0.
      !     t_precip = min(t_precip,270.)
        ! elseif (rssca==0..and. forc_t<(tfrz+tcrit)) then
        !   prc_snow = prc_snow + prc_rain
        !   prl_snow = prl_snow + prl_rain
        !   prc_rain = 0.
        !   prl_rain = 0.  
        !   t_precip = min(t_precip,270.)
      !   endif  
      ! endif     

      forc_rain = prc_rain + prl_rain
      forc_snow = prc_snow + prl_snow

!======================================================================

                       !         / SOIL GROUND        (itypwat = 0)
if(itypwat<=2)then     ! <=== is - URBAN and BUILT-UP (itypwat = 1)
                       !         \ WETLAND            (itypwat = 2)

!======================================================================
                          !initial set
      scvold = scv        !snow mass at previous time step

      snl = 0
      do j=maxsnl+1,0
         if(wliq_soisno(j)+wice_soisno(j)>0.) snl=snl-1
      enddo

      zi_soisno(0)=0.
      if(snl<0)then
      do j = -1, snl, -1
         zi_soisno(j)=zi_soisno(j+1)-dz_soisno(j+1)
      enddo
      endif
      do j = 1,nl
         zi_soisno(j)=zi_soisno(j-1)+dz_soisno(j)
      enddo

      totwb = ldew + scv + sum(wice_soisno(1:nl_soil)+wliq_soisno(1:nl_soil)) + wa + sst
      fiold(:) = 0.0
      if (snl <0 ) then
         fiold(snl+1:0)=wice_soisno(snl+1:0)/(wliq_soisno(snl+1:0)+wice_soisno(snl+1:0))
      endif

!----------------------------------------------------------------------
! [2] Canopy interception and precipitation onto ground surface
!----------------------------------------------------------------------

      call LEAF_interception (deltim,dewmx,chil,sigf,lai,sai,tlsun,&
                              prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,pg_rain,pg_snow,qintr)

! Debug, 2016-11-4. Ignoring leaf interception, debug code:16110401
! debug_prcsnow = 0.
! debug_prcrain = 0.
! debug_prlsnow = 0.
! debug_prlrain = 0.      
! call LEAF_interception (deltim,dewmx,chil,sigf,lai,sai,tlsun,&
!                               debug_prcrain,debug_prcsnow,debug_prlrain,debug_prlsnow,&
!                               ldew,pg_rain,pg_snow,qintr)      
! pg_rain = prc_rain+prl_rain+pg_rain
! pg_snow = prc_snow+prl_snow+pg_snow
! End debug
      qdrip = pg_rain + pg_snow

!----------------------------------------------------------------------
! [3] Initilize new snow nodes for snowfall / sleet
!----------------------------------------------------------------------
! print*, 'debug albedo 1' , scv  ,pg_snow,snowdp,forc_t,forc_prc+forc_prl

      call newsnow (itypwat,maxsnl,deltim,t_grnd,pg_rain,pg_snow,bifall,&
                    t_precip,zi_soisno(:0),z_soisno(:0),dz_soisno(:0),t_soisno(:0),&
                    wliq_soisno(:0),wice_soisno(:0),fiold(:0),snl,sag,scv,snowdp,fsno)

!----------------------------------------------------------------------
! [4] Energy and Water balance 
!----------------------------------------------------------------------
      lb  = snl + 1           ! lower bound of array 
      CALL THERMAL (itypwat ,lb      ,nl      ,nl_soil        ,deltim          ,&
           trsmx0           ,zlnd             ,zsno           ,csoilc          ,&
           dewmx            ,capr             ,cnfac          ,csol            ,&
           wrsd             ,adjfac           ,smf_soisno(-1:nl_soil) ,smf_ice  ,&
           porsl            ,psi0             ,bsw            ,dkdry           ,&
           dksatu           ,lai              ,laisun         ,laisha          ,&
           sai              ,z0m              ,displa         ,sqrtdi          ,&
           rootfr           ,rstfac           ,effcon         ,vmax25          ,&
           slti             ,hlti             ,shti           ,hhti            ,&
           trda             ,trdm             ,trop           ,gradm           ,&
           binter           ,extkn            ,forc_hgt_u     ,forc_hgt_t      ,&
           forc_hgt_q       ,forc_us          ,forc_vs        ,forc_t          ,&
           forc_q           ,forc_rhoair      ,forc_psrf      ,forc_pco2m      ,&
           forc_po2m        ,coszen           ,parsun         ,parsha          ,&
           sabvsun          ,sabvsha          ,sabg           ,forc_frl        ,&
           extkb            ,extkd            ,thermk         ,fsno            ,&
           sigf             ,dz_soisno(lb:)   ,z_soisno(lb:)  ,zi_soisno(lb-1:),&
           tlsun            ,tlsha            ,t_soisno(lb:)  ,wice_soisno(lb:),&
           wliq_soisno(lb:) ,ldew             ,scv            ,snowdp          ,&
           imelt(lb:)       ,taux             ,tauy           ,fsena           ,&
           fevpa            ,lfevpa           ,fsenl          ,fevpl           ,&
           etr              ,fseng            ,fevpg          ,olrg            ,&
           fgrnd            ,rootr(1:nl_soil) ,qseva          ,qsdew           ,&
           qsubl            ,qfros            ,sm             ,tref            ,&
           qref             ,trad             ,rst            ,assim           ,&
           respc            ,errore           ,emis           ,z0ma            ,&
           zol              ,rib              ,ustar          ,qstar           ,&
           tstar            ,fm               ,fh             ,fq              )

      CALL WATER (ivt,itypwat,lb,nl_soil,deltim  ,&    !MODIFIRD, 20151213
             z_soisno(lb:nl_soil)    ,dz_soisno(lb:nl_soil)   ,zi_soisno(lb-1:nl_soil)    ,&
             sst         ,ss_f        ,wfld(1:nl_soil),&
             wsat        ,wrsd        ,watern      ,alpha   ,slope  ,&
             length      ,Ds          ,Dr          ,Dg      ,kground,&
             bsw(1:nl_soil),porsl(1:nl_soil),psi0(1:nl_soil),hksati(1:nl_soil),rootr(1:nl_soil)  ,&
             t_soisno(lb:nl_soil)    ,wliq_soisno(lb:nl_soil) ,wice_soisno(lb:nl_soil) ,pg_rain ,sm,&
             etr         ,qseva       ,qsdew       ,qsubl   ,qfros  ,&
             rsur        ,rnof        ,qinfl       ,wtfact  ,pondmx ,&
             ssi         ,wimp        ,smpmin      ,zwt     ,wa     ,&
             Drw         ,qcharge     ,qlat        ,qground ,qsmelt ,&
             smf_soisno(-1:),snoweva)  
      if(snl<0)then
         ! Compaction rate for snow 
         ! Natural compaction and metamorphosis. The compaction rate
         ! is recalculated for every new timestep
         lb  = snl + 1   ! lower bound of array 
         call snowcompaction (lb,deltim,&
                         imelt(lb:0),fiold(lb:0),t_soisno(lb:0),&
                         wliq_soisno(lb:0),wice_soisno(lb:0),dz_soisno(lb:0))

         ! Combine thin snow elements
         lb = maxsnl + 1
         call snowlayerscombine (lb,snl,&
                         z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:1),&
                         wliq_soisno(lb:1),wice_soisno(lb:1),t_soisno(lb:1),scv,snowdp)

         ! Divide thick snow elements
         if(snl<0) &
         call snowlayersdivide (lb,snl,&
                         z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
                         wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0))
      endif
      ! Set zero to the empty node
      if (snl > maxsnl) then
         wice_soisno(maxsnl+1:snl) = 0.
         wliq_soisno(maxsnl+1:snl) = 0.
         t_soisno   (maxsnl+1:snl) = 0.
         z_soisno   (maxsnl+1:snl) = 0.
         dz_soisno  (maxsnl+1:snl) = 0.
      endif

      lb = snl + 1
      t_grnd = t_soisno(lb)

      ! ----------------------------------------
      ! energy balance
      ! ----------------------------------------
      zerr=errore
!#if(defined CLMDEBUG)
      if(abs(errore)>.5 .and. qground>0.)then
         write(6,*) 'Warning: energy balance violation ',errore,ivt,icol,irow,z0m,snowdp
      endif
!#endif

      ! ----------------------------------------
      ! water balance
      ! ----------------------------------------
      endwb=sum(wice_soisno(1:nl_soil)+wliq_soisno(1:nl_soil))+ldew+scv + wa+sst
      errorw=(endwb-totwb)-(forc_prc+forc_prl-fevpa-rnof)*deltim 
      if(itypwat==2) errorw=0.    ! wetland
      xerr=errorw/deltim

!#if(defined CLMDEBUG)
      if(abs(errorw)>1.e-3 .and.ivt /= 14) then
         write(6,*) 'Warning: water balance violation(BUG in CLMAIN 558)',icol,irow
         write(6,*) 'variables:', icol,irow
         write(6,*) 'errorw:'   ,errorw  
         write(6,*) 'sigf:'   ,sigf 
         write(6,*) 'z0m:'   ,z0m
         write(6,*) 'snowdp:'   ,snowdp 
         write(6,*) 'ivt:'      ,ivt, etr,rootfr ,sum(rootfr)  
         write(6,*) 'forc_prec:',(forc_prc+forc_prl)*deltim
         write(6,*) 'forcing:',forc_t,forc_q,forc_rh,forc_psrf
         write(6,*) 'runoff:' ,rsur,qcharge,qlat,qground,qsmelt
         write(6,*) 'wa:'       ,wa    
         write(6,*) 'ldew:'     ,ldew     
         write(6,*) 'fevpa:'    ,fevpa *deltim ,fevpl*deltim ,fevpg*deltim
         write(6,*) 'rnof:'     ,rnof  *deltim  
         write(6,*) 'scv:'      ,scv     
         write(6,*) 'scvold:'   ,scvold  
         write(6,*) 'endwb:'    ,endwb   
         write(6,*) 'totwb:'    ,totwb   
         write(6,*) 'wice_soisno(0:2):',wice_soisno(0:2)
         write(6,*) 'wliq_soisno(0:2):',wliq_soisno(0:2)
         !stop
      endif
!#endif

!======================================================================

else if(itypwat == 3)then   ! <=== is LAND ICE (glacier/ice sheet) (itypwat = 3) 

!======================================================================
                            !initial set
      scvold = scv          !snow mass at previous time step

      snl = 0
      do j=maxsnl+1,0
         if(wliq_soisno(j)+wice_soisno(j)>0.) snl=snl-1
      enddo

      zi_soisno(0)=0.
      if(snl<0)then
      do j = -1, snl, -1
         zi_soisno(j)=zi_soisno(j+1)-dz_soisno(j+1)
      enddo
      endif
      do j = 1,nl_soil
         zi_soisno(j)=zi_soisno(j-1)+dz_soisno(j)
      enddo

      totwb = scv + sum(wice_soisno(1:)+wliq_soisno(1:))+sst+wa
      fiold(:) = 0.0
      if (snl <0 ) then
         fiold(snl+1:0)=wice_soisno(snl+1:0)/(wliq_soisno(snl+1:0)+wice_soisno(snl+1:0))
      endif

      pg_rain = prc_rain + prl_rain
      pg_snow = prc_snow + prl_snow

      !----------------------------------------------------------------
      ! Initilize new snow nodes for snowfall / sleet
      !----------------------------------------------------------------

      call newsnow (itypwat,maxsnl,deltim,t_grnd,pg_rain,pg_snow,bifall,&
                    t_precip,zi_soisno(:0),z_soisno(:0),dz_soisno(:0),t_soisno(:0),&
                    wliq_soisno(:0),wice_soisno(:0),fiold(:0),snl,sag,scv,snowdp,fsno)

      !----------------------------------------------------------------
      ! Energy and Water balance 
      !----------------------------------------------------------------
      lb  = snl + 1            ! lower bound of array 

      CALL GLACIER_TEMP (lb    ,nl_soil     ,deltim     ,&
                   zlnd        ,zsno        ,capr       ,cnfac       ,&
                   forc_hgt_u  ,forc_hgt_t  ,forc_hgt_q ,forc_us     ,&
                   forc_vs     ,forc_t      ,forc_q     ,forc_rhoair ,&
                   forc_psrf   ,coszen      ,sabg       ,forc_frl    ,&
                   fsno,dz_soisno(lb:),z_soisno(lb:),zi_soisno(lb-1:),&
                   t_soisno(lb:),wice_soisno(lb:),wliq_soisno(lb:)   ,&
                   scv         ,snowdp      ,imelt(lb:) ,taux        ,&
                   tauy        ,fsena       ,fevpa      ,lfevpa      ,&
                   fseng       ,fevpg       ,olrg       ,fgrnd       ,&
                   qseva       ,qsdew       ,qsubl      ,qfros       ,&
                   sm          ,tref        ,qref       ,trad        ,&
                   errore      ,emis        ,z0ma       ,zol         ,&
                   rib         ,ustar       ,qstar      ,tstar       ,&
                   fm          ,fh          ,fq)


      CALL GLACIER_WATER (nl_soil,maxsnl,deltim                      ,&
                    forc_t,&
                   z_soisno    ,dz_soisno   ,zi_soisno  ,t_soisno    ,&
                   wliq_soisno ,wice_soisno ,pg_rain    ,pg_snow     ,&
                   sm          ,scv         ,snowdp     ,imelt       ,&
                   fiold       ,snl         ,qseva      ,qsdew       ,&
                   qsubl       ,qfros       ,rsur       ,rnof        ,&
                   ssi         ,wimp        ,wa )


      lb = snl + 1
      t_grnd = t_soisno(lb)
      qsmelt = rnof ! Assume runoff from glacier is all snowmelt, 2016-11-29. by LIHY

      ! ----------------------------------------
      ! energy and water balance check
      ! ----------------------------------------
      zerr=errore

      endwb=scv+sum(wice_soisno(1:)+wliq_soisno(1:))+sst+wa
      errorw=(endwb-totwb)-(pg_rain+pg_snow-fevpa-rnof)*deltim
      xerr=errorw/deltim

!======================================================================

else if(itypwat == 4)then   ! <=== is LAND WATER BODIES (lake, reservior and river) (itypwat = 4) 

!======================================================================

      snl = 0
      do j=maxsnl+1,0
         if (wliq_soisno(j)+wice_soisno(j) > 0.) then
            snl=snl-1
         endif
      enddo

      zi_soisno(0) = 0.
      if (snl <0 ) then
         do j = -1, snl, -1
            zi_soisno(j)=zi_soisno(j+1)-dz_soisno(j+1)
         enddo
      endif

      do j = 1,nl_soil
         zi_soisno(j)=zi_soisno(j-1)+dz_soisno(j)
      enddo

      scvold = scv          ! snow mass at previous time step
      fiold(:) = 0.0
      if (snl < 0) then
         fiold(snl+1:0)=wice_soisno(snl+1:0)/(wliq_soisno(snl+1:0)+wice_soisno(snl+1:0))
      endif

      w_old = sum(wliq_soisno(snl+1:))

      pg_rain = prc_rain + prl_rain
      pg_snow = prc_snow + prl_snow

      CALL newsnow_lake ( &
           maxsnl       , nl_lake      , deltim          , dz_lake , &
           pg_rain      , pg_snow      , t_precip        , bifall  , &
           t_lake       , zi_soisno(:0), z_soisno(:0)    , &         
           dz_soisno(:0), t_soisno(:0) , wliq_soisno(:0) , wice_soisno(:0) ,&
           fiold(:0)    , snl          , sag             , scv             ,&
           snowdp       , lake_icefrac )

 
      CALL laketem ( &
           itypwat      , maxsnl       , nl_soil         , nl_lake         ,&
           dlat         , deltim       , forc_hgt_u      , forc_hgt_t      ,&
           forc_hgt_q   , forc_us      , forc_vs         , forc_t          ,&
           forc_q       , forc_rhoair  , forc_psrf       , forc_sols       ,&
           forc_soll    , forc_solsd   , forc_solld      , sabg            ,&
           forc_frl     , dz_soisno    , z_soisno        , zi_soisno       ,&
           dz_lake      , lakedepth    , csol            , porsl           ,&
           dkdry        , dksatu       , &
           t_grnd       , scv          , snowdp          , t_soisno        ,&
           wliq_soisno  , wice_soisno  , imelt           , t_lake          ,&
           lake_icefrac , &
           taux         , tauy         , fsena                             ,&
           fevpa        , lfevpa       , fseng           , fevpg           ,&
           qseva        , qsubl        , qsdew           , qfros           ,&
           olrg         , fgrnd        , tref            , qref            ,&
           trad         , emis         , z0ma            , zol             ,&
           rib          , ustar        , qstar           , tstar           ,&
           fm           , fh           , fq              , sm )

      CALL snowwater_lake ( &
           maxsnl       , nl_soil      , nl_lake         , deltim          ,&
           ssi          , wimp         , porsl           , pg_rain         ,&
           pg_snow      , dz_lake      , imelt(:0)       , fiold(:0)       ,&
           qseva        , qsubl        , qsdew           , qfros           ,&
           z_soisno     , dz_soisno   , zi_soisno        , t_soisno        ,&
           wice_soisno  , wliq_soisno , t_lake           , lake_icefrac    ,&
           fseng        , fgrnd       , snl              , scv             ,&
           snowdp       , sm )

      ! We assume the land water bodies have zero extra liquid water capacity 
      ! (i.e.,constant capacity), all excess liquid water are put into the runoff,
      ! this unreasonable assumption should be updated in the future version
      a = (sum(wliq_soisno(snl+1:))-w_old)/deltim
      aa = qseva-(qsubl-qsdew)
      rsur = max(0., pg_rain - aa - a)
      rnof = rsur

      ! Set zero to the empty node
      if (snl > maxsnl) then
         wice_soisno(maxsnl+1:snl) = 0.
         wliq_soisno(maxsnl+1:snl) = 0.
         t_soisno   (maxsnl+1:snl) = 0.
         z_soisno   (maxsnl+1:snl) = 0.
         dz_soisno  (maxsnl+1:snl) = 0.
      endif


!======================================================================

else                     ! <=== is OCEAN (itypwat >= 99) 

!======================================================================
! simple ocean-sea ice model

    tssea = t_grnd
    tssub (1:7) = t_soisno (1:7) 
    CALL SOCEAN (deltim,oro,forc_hgt_u,forc_hgt_t,forc_hgt_q,&
                 forc_us,forc_vs,forc_t,forc_t,forc_rhoair,forc_psrf,&
                 sabg,forc_frl,tssea,tssub(1:7),scv,&
                 taux,tauy,fsena,fevpa,lfevpa,fseng,fevpg,tref,qref,&
                 z0ma,zol,rib,ustar,qstar,tstar,fm,fh,fq,emis,olrg)
                 
               ! null data for sea component
                 z_soisno   (:) = 0.0
                 dz_soisno  (:) = 0.0
                 t_soisno   (:) = 0.0
                 t_soisno (1:7) = tssub(1:7)
                 wliq_soisno(:) = 0.0
                 wice_soisno(:) = 0.0
                 t_grnd  = tssea
                 snowdp  = scv/1000.*20. 

                 trad    = tssea
                 fgrnd   = 0.0
                 rsur    = 0.0
                 rnof    = 0.0

!======================================================================

endif


!======================================================================
! Preparation for the next time step
! 1) time-varying parameters for vegatation
! 2) fraction of snow cover 
! 3) solar zenith angle and
! 4) albedos 
!======================================================================

    ! cosine of solar zenith angle 
    ! call calendarday_func(idate, lon,calday)
    calday = float(idate(2)) + float(idate(3))/86400. ! modified by Hongyi Li, based on gmt time not local time
    call orb_coszen_func(calday,dlon,dlat,coszen)

    if(itypwat <= 5)then   ! LAND
!#if(defined DYN_PHENOLOGY)
       ! need to update lai and sai, fveg, green, they are done once in a day only
       ! if(dolai)then
          call LAI_empirical(ivt,nl_soil,rootfr,t_soisno(1:),lai,sai,fveg,green)
       ! endif
!#endif

       ! fraction of snow cover.
       call snowfraction (fveg,z0m,zlnd,scv,snowdp,wt,sigf,fsno)

       ! update the snow age 
       if(snl==0) sag=0.
       call snowage (deltim,t_grnd,scv,scvold,sag)

       ! water volumetric content of soil surface layer [m3/m3]
       ssw = min(1.,1.e-3*wliq_soisno(1)/dz_soisno(1))
       if(itypwat >=3) ssw = 1.0
       ! albedos 
       ! we supposed call it every time-step, because 
       ! other vegeation related parameters are needed to create
       !*if(doalb)then
            call albland (itypwat,&
                 soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb,&
                 chil,ref,tran,fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,t_grnd,&
                 alb,albg,albv,ssun,ssha,thermk,extkb,extkd)
       !*endif
! print*,'debug albg,albv,alb',albg,albv,alb       
    else                   ! OCEAN
       sag = 0.0
       !*if(doalb)then
            call albocean (oro,scv,coszen,alb)
       !*endif
    endif

    ! zero-filling set for glacier/ice-sheet/land water bodies/ocean components
    if(itypwat > 2)then
       lai = 0.0
       sai = 0.0
       laisun = 0.0
       laisha = 0.0
       green = 0.0
       fveg = 0.0
       sigf = 0.0

       albg(:,:) = alb(:,:)
       albv(:,:) = 0.0
       ssun(:,:) = 0.0
       ssha(:,:) = 0.0
       thermk = 0.0
       extkb = 0.0
       extkd = 0.0

       tlsun   = forc_t
       tlsha   = forc_t
       ldew    = 0.0
       fsenl   = 0.0
       fevpl   = 0.0
       etr     = 0.0
       assim   = 0.0
       respc   = 0.0

       zerr=0.
       xerr=0.

       qinfl = 0.
       qdrip = forc_rain + forc_snow
       qintr = 0.
       h2osoi = 0. 
       rstfac = 0.
       zwt   = 0.
       wa    = 0.1*Dg*1000.
       qcharge = 0.
    endif
! ---------------------------    
! update swe information
    ! snl = 0
    ! do j=maxsnl+1,0
       ! if(wliq_soisno(j)+wice_soisno(j)>0.) snl=snl-1
    ! enddo
    ! if (snl<0) then
        ! lb = snl+1
        ! scv = sum(wliq_soisno(lb:0) + wice_soisno(lb:0))
    ! else
        ! scv = 0.
    ! endif
    
    h2osoi = wliq_soisno(1:)/(dz_soisno(1:)*denh2o) + wice_soisno(1:)/(dz_soisno(1:)*denice)
    wat = sum(wice_soisno(1:)+wliq_soisno(1:))+ldew+scv + wa
!----------------------------------------------------------------------



END SUBROUTINE CLMMAIN