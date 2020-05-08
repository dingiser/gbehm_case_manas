  PROGRAM CLMGBHM
! ======================================================================
! Reference: 
!     [1] Dai et al., 2003: The Common Land Model (CoLM). 
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. J. Climate, 17: 2281-2299.
!     [3] Dai et al., 2014: The Terrestrial Modeling System (TMS).
!
!     Created by Yongjiu Dai, Februay 2004
!     Revised by Yongjiu Dai and Hua Yuan, April 2014
!     Revised by Hongyi Li, September 2015
! ======================================================================

      use precision
      use PhysicalConstants
      !use timemanager
      use omp_lib

      IMPLICIT NONE
!if(defined USGS_CLASSIFICATION)
      integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
!endif
!if(defined IGBP_CLASSIFICATION)
 !     integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
!endif
      integer, parameter :: nl_soil  = 10   ! number of soil layers
      integer, parameter :: nl_lake  = 10   ! number of lake layers
      integer, parameter :: maxsnl   = -5   ! max number of snow layers
      integer, parameter :: lon_points = 100 ! number of lon points
      integer, parameter :: lat_points = 100 ! number of lat points
      real(r8),parameter :: deltim = 3600  ! seconds in a time-step
  ! ------------------------ Dummy Argument ------------------------------

  integer   :: idate(3) ! model calendar for next time step (year, julian day, seconds)


  ! logical,  INTENT(in) :: dolai    ! true if time for time-varying vegetation paramter
  ! logical,  INTENT(in) :: doalb    ! true if time for surface albedo calculation
  ! logical,  INTENT(in) :: dosst    ! true if time for update sst/ice/snow

  real(r8)  :: &
        lons(lon_points)                  , &! longitude in degree
        dlon(lon_points,lat_points)       , &! logitude in radians
        dlat(lon_points,lat_points)       , &   ! latitude in radians
        latixy(lon_points,lat_points)     , &
        longxy(lon_points,lat_points) 
  integer :: & 
        ivt(lon_points,lat_points)        , &! land cover type of USGS classification or others
        itypwat(lon_points,lat_points)    ! land water type (0=soil, 1=urban and built-up, 
                                       ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)
  real(r8)  :: &
        lakedepth(lon_points,lat_points)       , &! lake depth (m)
        dz_lake(nl_lake,lon_points,lat_points), &! lake layer thickness (m)
        ! soil physical parameters and lake info
        soil_s_v_alb(lon_points,lat_points)   , &! albedo of visible of the saturated soil
        soil_d_v_alb(lon_points,lat_points)   , &! albedo of visible of the dry soil
        soil_s_n_alb(lon_points,lat_points)   , &! albedo of near infrared of the saturated soil
        soil_d_n_alb(lon_points,lat_points)   , &! albedo of near infrared of the dry soil
        porsl(nl_soil,lon_points,lat_points) , &! fraction of soil that is voids [-]
        psi0(nl_soil,lon_points,lat_points)  , &! minimum soil suction [mm]
        bsw(nl_soil,lon_points,lat_points)   , &! clapp and hornbereger "b" parameter [-]
        hksati(nl_soil,lon_points,lat_points), &! hydraulic conductivity at saturation [mm h2o/s]
        csol(nl_soil,lon_points,lat_points)  , &! heat capacity of soil solids [J/(m3 K)]
        dksatu(nl_soil,lon_points,lat_points), &! thermal conductivity of saturated soil [W/m-K]
        dkdry(nl_soil,lon_points,lat_points) , &! thermal conductivity for dry soil  [J/(K s m)]
        rootfr(nl_soil,lon_points,lat_points), &! fraction of roots in each soil layer

        gravel_grid(nl_soil,lon_points,lat_points) , &
        sand_grid  (nl_soil,lon_points,lat_points) , &
        clay_grid  (nl_soil,lon_points,lat_points) , &
        SOC_grid   (nl_soil,lon_points,lat_points) , &
        BD_grid    (nl_soil,lon_points,lat_points) , &
        soil_t_grid(nl_soil,lon_points,lat_points) , &! soil layer temperature (K)
        soil_w_grid(nl_soil,lon_points,lat_points) , &! soil layer wetness (-)
        snow_d_grid(lon_points,lat_points)         , &! snow depth (m)        
        ! vegetation static, dynamic, derived parameters
        fveg (lon_points,lat_points)        , &! fraction of vegetation cover
        green(lon_points,lat_points)        , &! greenness
        lai  (lon_points,lat_points)        , &! leaf area index
        sai  (lon_points,lat_points)        , &! stem area index   
        
        z0m   (lon_points,lat_points)       , &! aerodynamic roughness length [m]
        displa(lon_points,lat_points)       , &! displacement height [m]
        sqrtdi(lon_points,lat_points)       , &! inverse sqrt of leaf dimension [m**-0.5]
        effcon(lon_points,lat_points)       , &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25(lon_points,lat_points)       , &! maximum carboxylation rate at 25 C at canopy top
        slti  (lon_points,lat_points)       , &! slope of low temperature inhibition function      [s3] 
        hlti  (lon_points,lat_points)       , &! 1/2 point of low temperature inhibition function  [s4]
        shti  (lon_points,lat_points)       , &! slope of high temperature inhibition function     [s1]
        hhti  (lon_points,lat_points)       , &! 1/2 point of high temperature inhibition function [s2]
        trda  (lon_points,lat_points)       , &! temperature coefficient in gs-a model             [s5]
        trdm  (lon_points,lat_points)       , &! temperature coefficient in gs-a model             [s6]
        trop  (lon_points,lat_points)       , &! temperature coefficient in gs-a model          
        gradm (lon_points,lat_points)       , &! conductance-photosynthesis slope parameter
        binter(lon_points,lat_points)       , &! conductance-photosynthesis intercep
        extkn (lon_points,lat_points)       , &! coefficient of leaf nitrogen allocation
        chil  (lon_points,lat_points)       , &! leaf angle distribution factor
        ref   (2,2,lon_points,lat_points)   , &! leaf reflectance (iw=iband, il=life and dead)
        tran  (2,2,lon_points,lat_points)   , &! leaf transmittance (iw=iband, il=life and dead)   

        ! tunable parameters
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


        z_soisno(maxsnl+1:nl_soil,lon_points,lat_points)   , &! layer depth (m)
        dz_soisno(maxsnl+1:nl_soil,lon_points,lat_points)  , &! layer thickness (m)
        t_soisno(maxsnl+1:nl_soil,lon_points,lat_points)   , &! soil + snow layer temperature [K]
        wliq_soisno(maxsnl+1:nl_soil,lon_points,lat_points), &! liquid water (kg/m2)
        wice_soisno(maxsnl+1:nl_soil,lon_points,lat_points), &! ice lens (kg/m2)

        t_lake(nl_lake,lon_points,lat_points)       ,&! lake temperature (kelvin)
        lake_icefrac(nl_lake,lon_points,lat_points) ,&! lake mass fraction of lake layer that is frozen

        h2osoi(nl_soil,lon_points,lat_points), &
        t_grnd(lon_points,lat_points)     , &! ground surface temperature [k]
        tlsun(lon_points,lat_points)      , &! sunlit leaf temperature [K]
        tlsha(lon_points,lat_points)      , &! shaded leaf temperature [K]
        ldew(lon_points,lat_points)       , &! depth of water on foliage [kg/m2/s]
        sag(lon_points,lat_points)        , &! non dimensional snow age [-]
        scv(lon_points,lat_points)        , &! snow mass (kg/m2)
        snowdp(lon_points,lat_points)     , &! snow depth (m)
        zwt(lon_points,lat_points)        , &! the depth to water table [m]
        wa(lon_points,lat_points)         , &! water storage in aquifer [mm]

        fsno(lon_points,lat_points)       , &! fractional snow cover
        sigf(lon_points,lat_points)       , &! fraction of veg cover, excluding snow-covered veg [-]

 
        coszen(lon_points,lat_points)     , &! cosine of solar zenith angle
        albg(2,2,lon_points,lat_points)   , &! albedo, ground [-]
        albv(2,2,lon_points,lat_points)   , &! albedo, vegetation [-]
        alb(2,2,lon_points,lat_points)    , &! averaged albedo [-]
        ssun(2,2,lon_points,lat_points)   , &! sunlit canopy absorption for solar radiation
        ssha(2,2,lon_points,lat_points)   , &! shaded canopy absorption for solar radiation
        thermk(lon_points,lat_points)     , &! canopy gap fraction for tir radiation
        extkb(lon_points,lat_points)      , &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd(lon_points,lat_points)      , &! diffuse and scattered diffuse PAR extinction coefficient
                    
        ! Additional variables required by reginal model (WRF & RSM) 
        trad(lon_points,lat_points),                   &! radiative temperature of surface [K]
        tref(lon_points,lat_points),                   &! 2 m height air temperature [kelvin]
        qref(lon_points,lat_points),                   &! 2 m height air specific humidity
        rst(lon_points,lat_points),                    &! canopy stomatal resistance (s/m)
        emis(lon_points,lat_points),                   &! averaged bulk surface emissivity
        z0ma(lon_points,lat_points),                   &! effective roughness [m]
        zol(lon_points,lat_points),                    &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib(lon_points,lat_points),                    &! bulk Richardson number in surface layer
        ustar(lon_points,lat_points),                  &! u* in similarity theory [m/s]
        qstar(lon_points,lat_points),                  &! q* in similarity theory [kg/kg]
        tstar(lon_points,lat_points),                  &! t* in similarity theory [K]
        fm(lon_points,lat_points),                     &! integral of profile function for momentum
        fh(lon_points,lat_points),                     &! integral of profile function for heat
        fq(lon_points,lat_points)                       ! integral of profile function for moisture

  real(r8) :: &
        taux(lon_points,lat_points)       , &! wind stress: E-W [kg/m/s**2]
        tauy(lon_points,lat_points)       , &! wind stress: N-S [kg/m/s**2]
        fsena(lon_points,lat_points)      , &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa(lon_points,lat_points)      , &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa(lon_points,lat_points)     , &! latent heat flux from canopy height to atmosphere [W/2]
        fsenl(lon_points,lat_points)      , &! ensible heat from leaves [W/m2]
        fevpl(lon_points,lat_points)      , &! evaporation+transpiration from leaves [mm/s]
        etr(lon_points,lat_points)        , &! transpiration rate [mm/s]
        fseng(lon_points,lat_points)      , &! sensible heat flux from ground [W/m2]
        fevpg(lon_points,lat_points)      , &! evaporation heat flux from ground [mm/s]
        olrg(lon_points,lat_points)       , &! outgoing long-wave radiation from ground+canopy
        fgrnd(lon_points,lat_points)      , &! ground heat flux [W/m2]
        xerr(lon_points,lat_points)       , &! water balance error at current time-step [mm/s]
        zerr(lon_points,lat_points)       , &! energy balnce errore at current time-step [W/m2]

        rsur(lon_points,lat_points)       , &! surface runoff (mm h2o/s)
        rnof(lon_points,lat_points)       , &! total runoff (mm h2o/s)
        qintr(lon_points,lat_points)      , &! interception (mm h2o/s)
        qinfl(lon_points,lat_points)      , &! inflitration (mm h2o/s)
        qdrip(lon_points,lat_points)      , &! throughfall (mm h2o/s)
        qcharge(lon_points,lat_points)    , &! groundwater recharge [mm/s]
       
        assim(lon_points,lat_points)      , &! canopy assimilation
        respc(lon_points,lat_points)      , &! canopy respiration

        sabvsun(lon_points,lat_points)    , &! solar absorbed by sunlit vegetation [W/m2]
        sabvsha(lon_points,lat_points)    , &! solar absorbed by shaded vegetation [W/m2]
        sabg(lon_points,lat_points)       , &! solar absorbed by ground  [W/m2]
        sr(lon_points,lat_points)         , &! total reflected solar radiation (W/m2)
        solvd(lon_points,lat_points)      , &! incident direct beam vis solar radiation (W/m2)
        solvi(lon_points,lat_points)      , &! incident diffuse beam vis solar radiation (W/m2)
        solnd(lon_points,lat_points)      , &! incident direct beam nir solar radiation (W/m2)
        solni(lon_points,lat_points)      , &! incident diffuse beam nir solar radiation (W/m2)
        srvd(lon_points,lat_points)       , &! reflected direct beam vis solar radiation (W/m2)
        srvi(lon_points,lat_points)       , &! reflected diffuse beam vis solar radiation (W/m2)
        srnd(lon_points,lat_points)       , &! reflected direct beam nir solar radiation (W/m2)
        srni(lon_points,lat_points)       , &! reflected diffuse beam nir solar radiation (W/m2)
        solvdln(lon_points,lat_points)    , &! incident direct beam vis solar radiation at local noon(W/m2)
        solviln(lon_points,lat_points)    , &! incident diffuse beam vis solar radiation at local noon(W/m2)
        solndln(lon_points,lat_points)    , &! incident direct beam nir solar radiation at local noon(W/m2)
        solniln(lon_points,lat_points)    , &! incident diffuse beam nir solar radiation at local noon(W/m2)
        srvdln(lon_points,lat_points)     , &! reflected direct beam vis solar radiation at local noon(W/m2)
        srviln(lon_points,lat_points)     , &! reflected diffuse beam vis solar radiation at local noon(W/m2)
        srndln(lon_points,lat_points)     , &! reflected direct beam nir solar radiation at local noon(W/m2)
        srniln(lon_points,lat_points)      ! reflected diffuse beam nir solar radiation at local noon(W/m2)
        
! Forcing
! ----------------------
  real(r8) :: &
      forc_xy_t       (lon_points,lat_points), &! temperature at agcm reference height [kelvin]
      forc_xy_q       (lon_points,lat_points), &! specific humidity at agcm reference height [kg/kg]
      forc_xy_psrf    (lon_points,lat_points), &! atmosphere pressure at the surface [pa]
      forc_xy_solarin (lon_points,lat_points), &! solar rad onto srf [W/m2]
      forc_xy_frl     (lon_points,lat_points), &! atmospheric infrared (longwave) radiation [W/m2]
      forc_xy_prec    (lon_points,lat_points), &! precipitation [mm/s]
      forc_xy_us      (lon_points,lat_points), &! wind speed in eastward direction [m/s]
      forc_xy_vs      (lon_points,lat_points), &! wind speed in northward direction [m/s]
              laisun(lon_points,lat_points)     , &! sunlit leaf area index
        laisha(lon_points,lat_points)     , &! shaded leaf area index
        rstfac(lon_points,lat_points)     , &! factor of soil water stress 
        wat(lon_points,lat_points)        , &! total water storage
        forc_rain(lon_points,lat_points)  , &! rain [mm/s]
        forc_snow(lon_points,lat_points)     ! snow [mm/s]        
! ----------------local variables ---------------------------------

      ! character(LEN=256) :: site  ! site name

      integer :: edate(3)         ! calendar (year, julian day, seconds)
      integer :: pdate(3)         ! calendar (year, julian day, seconds)
      logical :: solarin_all_band ! downward solar in broad band
      logical :: greenwich        ! greenwich time
                                  !
      ! logical :: doalb            ! true => start up the surface albedo calculation
      ! logical :: dolai            ! true => start up the time-varying vegetation paramter
      ! logical :: dosst            ! true => update sst/ice/snow
      logical :: lwrite           ! true: write output  file frequency
      logical :: rwrite           ! true: write restart file frequency
                                  !
      integer :: istep            ! looping step
      integer :: nac              ! number of accumulation
      real(r8) :: oro(lon_points,lat_points)           ! ocean(0)/seaice(2)/ flag
      integer :: Julian_1day_p, Julian_1day 
      integer :: Julian_8day_p, Julian_8day 
      integer :: s_year, s_julian, s_seconds
      integer :: e_year, e_julian, e_seconds
      integer :: p_year, p_julian, p_seconds
      integer :: i, j


      s_year    =  2000
      s_julian  =  1
      s_seconds =  0

      e_year    =  2001
      e_julian  =  1
      e_seconds =  0
      
      idate(1) = s_year; idate(2) = s_julian; idate(3) = s_seconds
      edate(1) = e_year; edate(2) = e_julian; edate(3) = e_seconds
      ! pdate(1) = p_year; pdate(2) = p_julian; pdate(3) = p_seconds
      
      !call adj2end(edate)
      ! call adj2end(pdate)
      

      ! ptstamp = pdate

      oro = 1.
      ! dolai=.true.
      ! doalb=.true.
      ! dosst=.true.
      
           forc_xy_t(:,:)  = 270.
           forc_xy_q(:,:)  =1.0
           forc_xy_psrf(:,:) =680.
           forc_xy_solarin(:,:) =300.
           forc_xy_frl(:,:) =100.
           forc_xy_prec(:,:) =0.1
           forc_xy_us(:,:) =3.
           forc_xy_vs(:,:) =3.
           
           ivt(:,:) = 1
           latixy(:,:) = 30.
           longxy(:,:) = 50.
           snow_d_grid(:,:) = 1.
           lakedepth(:,:) = 0.
           soil_t_grid(:,:,:) = 270.
           soil_w_grid(:,:,:) = 0.5
           gravel_grid(:,:,:) = 20.
           sand_grid(:,:,:) = 10.
           clay_grid(:,:,:) = 30.
           SOC_grid(:,:,:)=20.
           BD_grid(:,:,:)=20.
! Initial 
! =====================================
      call initialize (  &
         ! INPUTS
           idate,lon_points,lat_points,                     &
           nl_soil,maxsnl,nl_lake,ivt, latixy,longxy, &
           snow_d_grid,lakedepth,&
           gravel_grid,sand_grid,clay_grid,SOC_grid,BD_grid,&
           soil_t_grid,soil_w_grid,&
         ! OUTPUTS                   
           lons, dlon, dlat, itypwat, &
           dz_lake, & ! new lake scheme       

         ! soil information and lake depth
           soil_s_v_alb, soil_d_v_alb, soil_s_n_alb, soil_d_n_alb,  &
           porsl,   psi0,    bsw,     hksati,   &
           csol,    dksatu,  dkdry,   rootfr,   &

         ! vegetation information
           z0m,          displa,       sqrtdi,                      &
           effcon,       vmax25,       slti,         hlti,          &
           shti,         hhti,         trda,         trdm,          &
           trop,         gradm,        binter,       extkn,         &
           chil,         ref,          tran,                        &

         ! land surface variables required for restart
           z_soisno,     dz_soisno,    t_soisno,     wliq_soisno,   &
           wice_soisno,  &

           t_grnd,       tlsun,        tlsha,        ldew,          &
           sag,          scv,          snowdp,       fveg,          &
           fsno,         sigf,         green,        lai,           &
           sai,          coszen,       albg,         albv,          &
           alb,          ssun,         ssha,         thermk,        &
           extkb,        extkd,                                     &

           zwt,          wa,                                        &
           t_lake,       lake_icefrac,                              &

         ! FLUXES
           qref, rst,trad,tref,&
           ! taux,         tauy,         fsena,        fevpa,         &
           ! lfevpa,       fsenl,        fevpl,        etr,           &
           ! fseng,        fevpg,        olrg,         fgrnd,         &
           ! trad,         tref,         qref,         rsur,          &
           ! rnof,         qintr,        qinfl,        qdrip,         &
           ! rst,          assim,        respc,        sabvsun,       &
           ! sabvsha,      sabg,         sr,           solvd,         &
           ! solvi,        solnd,        solni,        srvd,          &
           ! srvi,         srnd,         srni,         solvdln,       &
           ! solviln,      solndln,      solniln,      srvdln,        &
           ! srviln,       srndln,       srniln,       qcharge,       &
           ! xerr,         zerr,                                      &

         ! TUNABLE modle constants
           zlnd,         zsno,         csoilc,       dewmx,         &
           wtfact,       capr,         cnfac,        ssi,           &
           wimp,         pondmx,       smpmax,       smpmin,        &
           trsmx0,       tcrit,                                     & 

         ! additional variables required by coupling with WRF model 
           emis,         z0ma,         zol,          rib,           &
           ustar,        qstar,        tstar,                       &
           fm,           fh,           fq )

     
! ======================================================================
! begin time stepping loop
! ======================================================================

      istep = 1

      TIMELOOP : DO while (istep < 100)
print*, 'TIMELOOP = ', istep
       ! Calendar for NEXT time step
       ! ----------------------------------------------------------------------
        ! CALL TICKTIME (deltim,idate)

       ! Call clm driver
       ! ----------------------------------------------------------------------
         CALL CLMDRIVER (nl_soil,maxsnl,nl_lake,lon_points,lat_points,deltim,&
           oro,&
           lons, dlon, dlat, ivt, itypwat, &
           lakedepth, dz_lake, & ! new lake scheme

         ! soil information and lake depth
           soil_s_v_alb, soil_d_v_alb, soil_s_n_alb, soil_d_n_alb,  &
           porsl,        psi0,         bsw,          hksati,        &
           csol,         dksatu,       dkdry,        rootfr,        &

         ! vegetation information
           z0m,          displa,       sqrtdi,                      &
           effcon,       vmax25,       slti,         hlti,          &
           shti,         hhti,         trda,         trdm,          &
           trop,         gradm,        binter,       extkn,         &
           chil,         ref,          tran,                        &

         ! atmospheric forcing
           forc_xy_t        , &
           forc_xy_q        , &
           forc_xy_psrf     , &
           forc_xy_solarin  , &
           forc_xy_frl      , &
           forc_xy_prec     , &
           forc_xy_us       , &
           forc_xy_vs       , &
         ! land surface variables required for restart
           idate,                                                   &
           z_soisno,     dz_soisno,    t_soisno,     wliq_soisno,   &
           wice_soisno,  &

           t_grnd,       tlsun,        tlsha,        ldew,          &
           sag,          scv,          snowdp,       fveg,          &
           fsno,         sigf,         green,        lai,           &
           sai,          coszen,       albg,         albv,          &
           alb,          ssun,         ssha,         thermk,        &
           extkb,        extkd,                                     &

           zwt,          wa,                                        &
           t_lake,       lake_icefrac,                              &

         ! additional diagnostic variables for output
           laisun,       laisha,                                    &
           rstfac,       h2osoi,       wat,                         &

         ! FLUXES
           forc_rain,    forc_snow,                                 &
           taux,         tauy,         fsena,        fevpa,         &
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

         ! TUNABLE modle constants
           zlnd,         zsno,         csoilc,       dewmx,         &
           wtfact,       capr,         cnfac,        ssi,           &
           wimp,         pondmx,       smpmax,       smpmin,        &
           trsmx0,       tcrit,                                     & 

         ! additional variables required by coupling with WRF model 
           emis,         z0ma,         zol,          rib,           &
           ustar,        qstar,        tstar,                       &
           fm,           fh,           fq )


       ! Get leaf area index
       ! ----------------------------------------------------------------------
! #if(!defined DYN_PHENOLOGY)
       !READ in Leaf area index and stem area index
       ! Update every 8 days (time interval of the MODIS LAI data) 
       ! ----------------------------------------------------------------------
         ! Julian_8day = int(calendarday(idate)-1)/8*8 + 1
         ! if(Julian_8day /= Julian_8day_p)then
            ! CALL LAI_readin (lon_points,lat_points,&
                             ! Julian_8day,numpatch,dir_model_landdata)
         ! endif
! #else
       !Update once a day
         ! dolai = .false.
         ! Julian_1day = int(calendarday(idate)-1)/1*1 + 1
         ! if(Julian_1day /= Julian_1day_p)then
            ! dolai = .true.
         ! endif
! #endif

       ! Mapping subgrid patch [numpatch] vector of subgrid points to 
       !     -> [lon_points]x[lat_points] grid average
       ! ----------------------------------------------------------------------


         ! do j = 1, lat_points
            ! do i = 1, lon_points
               ! if(a_rnof(i,j) < 1.e-10) a_rnof(i,j) = 0.
            ! enddo
         ! enddo

! #if(defined CaMa_Flood)
!       Simulating the hydrodynamics in continental-scale rivers
!       ----------------------------------------------------------------------
         ! r2roffin(:,:) = a_rnof(:,:)*86400.  ! total runoff [mm/s] -> [mm/day]
         ! CALL CaMaMAIN(iyyyy,imm,idd,istep,r2roffin)
! #endif

       ! Logical idenfication for writing output and restart file
         ! lwrite = .false.
         ! rwrite = .false.
         ! CALL lpwrite(idate,deltim,lwrite,rwrite)

       ! Write out the model variables for restart run and the histroy file
       ! ----------------------------------------------------------------------
         ! if ( lwrite ) then

            ! if ( .NOT. (itstamp<=ptstamp) ) then
               ! CALL flxwrite (idate,nac,nac_ln,lon_points,lat_points,nl_soil,maxsnl,nl_lake,dir_output,site)
            ! endif

          ! Setting for next time step 
          ! ----------------------------------------------------------------------
            ! call FLUSH_2D_Fluxes
            ! nac = 0; nac_ln(:,:) = 0
         ! endif

         ! if ( rwrite ) then
            ! if ( .NOT. (itstamp<=ptstamp) ) then
               ! CALL WRITE_TimeVariables (idate,dir_restart_hist,site)
            ! endif
         ! endif

         istep = istep + 1

      END DO TIMELOOP
 
!#ifdef usempi
     ! call mpi_finalize(ierr)
!#endif
      write(6,*) 'CLM Execution Completed'

  END PROGRAM CLMGBHM
! ----------------------------------------------------------------------
! EOP
