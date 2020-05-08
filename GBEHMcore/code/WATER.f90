SUBROUTINE WATER (mask,run_pbsm,nl_soil,maxsnl,nl_lake,lon_points,lat_points,deltim,oro,&
           lons, dlon, dlat, ivt, itypwat, lakedepth, dz_lake, & 
           soil_s_v_alb, soil_d_v_alb, soil_s_n_alb, soil_d_n_alb,  &
           wsat,wrsd,watern,alpha,slope,length,Ds,Dr,Dg,kground,    &
           porsl,        psi0,         bsw,          hksati,        &
           wfld,      csol,         dksatu,       dkdry,        &
           rootfr,        &         
           z0m,          displa,       sqrtdi,                      &
           effcon,       vmax25,       slti,         hlti,          &
           shti,         hhti,         trda,         trdm,          &
           trop,         gradm,        binter,       extkn,         &
           chil,         ref,          tran,                        & 
           forc_pco2m,   forc_po2m,    forc_us,      forc_vs,       &
           forc_t,       forc_q,       forc_rh,      forc_prc,      &
           forc_prl,     forc_psrf,    forc_pbot,    forc_sols,     &
           forc_soll,    forc_solsd,   forc_solld,   forc_frl,      &
           forc_hgt_u,   forc_hgt_t,   forc_hgt_q,   forc_rhoair ,  &       
           idate,                                                   &
           z_soisno,     dz_soisno,    t_soisno,     wliq_soisno,   &
           wice_soisno,  sst, &
           t_grnd,       tlsun,        tlsha,        ldew,          &
           sag,          scv,          snowdp,       fveg,          &
           fsno,         sigf,         green,        lai,           &
           sai,          coszen,       albg,         albv,          &
           alb,          ssun,         ssha,         thermk,        &
           extkb,        extkd,                                     &
           zwt,          wa,                                        &
           t_lake,       lake_icefrac,                              &         
           laisun,       laisha,                                    &
           region_qsub,  net_blow,     blowing_add_region,          &
           rstfac,       h2osoi,       wat,          qlat,          &         
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
           zlnd,         zsno,         csoilc,       dewmx,         &
           wtfact,       capr,         cnfac,        ssi,           &
           wimp,         pondmx,       smpmax,       smpmin,        &
           trsmx0,       tcrit,        adjfac,                      & 
           Drw_gbhm,     Drw,          qground,      qsmelt,        &
           run_gbhm,     subbasin,     area,         dx,            &
           Dr_gbhm,      wriver,       s0,           mnroughness,   &
           nflow,        psubbasin,    nbasinup,     nsub,          &
           ngrid,        grid_row,     grid_col,     start,         &
           year,         month,        day,          hour,          &
           nc,           startyear,    endyear,                     & 
           dayinmonth,   startmonth,   endmonth,     startday,      &
           endday,       idc,          ihc,          start_sub,     &
           end_sub,      qr1,          qlin1,        Qd,            &
           Qm,           q1,           qr2,          q2,          Qh,&
           ssf) 
!=======================================================================
!
! CLM MODEL DRIVER
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002, 03/2014
! Revised : Hongyi Li, 09/20/2015
! variables behind 'run_gbhm' belong to GBHM module
!=======================================================================

 use precision
 use PhysicalConstants, only : tfrz, rgas, vonkar
 use omp_lib
 use gbhm_para, only: np,nx

 IMPLICIT NONE

  ! ------------------------ Dummy Argument ------------------------------
  integer,  INTENT(in) :: nl_soil  ! number of soil layers
  integer,  INTENT(in) :: nl_lake  ! number of lake layers
  integer,  INTENT(in) :: maxsnl   ! max number of snow layers
  integer,  INTENT(in) :: lon_points  ! number of lon points
  integer,  INTENT(in) :: lat_points  ! number of lat points
  integer,  INTENT(in) :: idate(3) ! model calendar for next time step (year, julian day, seconds)
  integer,  INTENT(in) :: run_pbsm ! 0=dont run PBSM;1=run PBSM
  real(r8), INTENT(in) :: deltim   ! seconds in a time-step

  real(r8), INTENT(in) :: &
        lons(lon_points),                   &! longitude in degree
        dlon(lon_points,lat_points)       , &! logitude in radians
        dlat(lon_points,lat_points)          ! latitude in radians

  integer, INTENT(in) :: & 
        mask(lon_points,lat_points)       , &
        ivt(lon_points,lat_points)        , &! land cover type of USGS classification or others
        itypwat(lon_points,lat_points)       ! land water type (0=soil, 1=urban and built-up, 
                                       ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)
                                       
! Parameters
! ----------------------
  real(r8), INTENT(in) :: &
        lakedepth(lon_points,lat_points)      , &! lake depth (m)
        dz_lake(nl_lake,lon_points,lat_points), &! lake layer thickness (m)
        soil_s_v_alb(lon_points,lat_points)   , &! albedo of visible of the saturated soil
        soil_d_v_alb(lon_points,lat_points)   , &! albedo of visible of the dry soil
        soil_s_n_alb(lon_points,lat_points)   , &! albedo of near infrared of the saturated soil
        soil_d_n_alb(lon_points,lat_points)   , &! albedo of near infrared of the dry soil
        porsl(nl_soil,lon_points,lat_points)  , &! fraction of soil that is voids [-]
        psi0(nl_soil,lon_points,lat_points)   , &! minimum soil suction [mm]
        wfld(nl_soil,lon_points,lat_points)   , &! water field capicity
        bsw(nl_soil,lon_points,lat_points)    , &! clapp and hornbereger "b" parameter [-]
        hksati(nl_soil,lon_points,lat_points) , &! hydraulic conductivity at saturation [mm h2o/s]
        csol(nl_soil,lon_points,lat_points)   , &! heat capacity of soil solids [J/(m3 K)]
        dksatu(nl_soil,lon_points,lat_points) , &! thermal conductivity of saturated soil [W/m-K]
        dkdry(nl_soil,lon_points,lat_points)  , &! thermal conductivity for dry soil  [J/(K s m)]
        rootfr(nl_soil,lon_points,lat_points) , &! fraction of roots in each soil layer
        z0m(lon_points,lat_points)        , &! aerodynamic roughness length [m]
        displa(lon_points,lat_points)     , &! displacement height [m]
        sqrtdi(lon_points,lat_points)     , &! inverse sqrt of leaf dimension [m**-0.5]
        effcon(lon_points,lat_points)     , &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25(lon_points,lat_points)     , &! maximum carboxylation rate at 25 C at canopy top
        slti(lon_points,lat_points)       , &! slope of low temperature inhibition function      [s3] 
        hlti(lon_points,lat_points)       , &! 1/2 point of low temperature inhibition function  [s4]
        shti(lon_points,lat_points)       , &! slope of high temperature inhibition function     [s1]
        hhti(lon_points,lat_points)       , &! 1/2 point of high temperature inhibition function [s2]
        trda(lon_points,lat_points)       , &! temperature coefficient in gs-a model             [s5]
        trdm(lon_points,lat_points)       , &! temperature coefficient in gs-a model             [s6]
        trop(lon_points,lat_points)       , &! temperature coefficient in gs-a model          
        gradm(lon_points,lat_points)      , &! conductance-photosynthesis slope parameter
        binter(lon_points,lat_points)     , &! conductance-photosynthesis intercep
        extkn(lon_points,lat_points)      , &! coefficient of leaf nitrogen allocation
        chil(lon_points,lat_points)       , &! leaf angle distribution factor
        ref(2,2,lon_points,lat_points)    , &! leaf reflectance (iw=iband, il=life and dead)
        tran(2,2,lon_points,lat_points)      ! leaf transmittance (iw=iband, il=life and dead)

  real(r8), intent(in) :: & ! added by HONGYILI, 2015.12.13
        wsat   (lon_points,lat_points)    , &
        wrsd   (lon_points,lat_points)    , &
        watern (lon_points,lat_points)    , &
        alpha  (lon_points,lat_points)    , &
        Ds     (lon_points,lat_points)    , &   ! depth of topsoil(m)
        Dr     (lon_points,lat_points)    , &   ! depth of river (m)
        Dg     (lon_points,lat_points)    , &   ! depth of unconfined acquifer (m)
        length (lon_points,lat_points)    , &   ! average hillslope length (m)
        slope  (lon_points,lat_points)    , &   ! slope of hillslope (m)
        kground(lon_points,lat_points)          ! GW hydraulic conductivity (m/s)      
        
  real(r8), INTENT(in) :: &! tunable parameters  
        ssf(lon_points,lat_points),&      
        zlnd       , &!roughness length for soil [m]
        zsno       , &!roughness length for snow [m]
        csoilc     , &!drag coefficient for soil under canopy [-]
        adjfac     , &!an adjustment parameter to solve the calibration of surface evaporation
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
        tcrit         !critical temp. to determine rain or snow

! Forcing
! ----------------------
  ! real(r8), INTENT(in) :: &
  !     forc_xy_t       (lon_points,lat_points), &! temperature at agcm reference height [kelvin]
  !     forc_xy_q       (lon_points,lat_points), &! specific humidity at agcm reference height [kg/kg]
  !     forc_xy_psrf    (lon_points,lat_points), &! atmosphere pressure at the surface [pa]
  !     forc_xy_solarin (lon_points,lat_points), &! solar rad onto srf [W/m2]
  !     forc_xy_frl     (lon_points,lat_points), &! atmospheric infrared (longwave) radiation [W/m2]
  !     forc_xy_prec    (lon_points,lat_points), &! precipitation [mm/s]
  !     forc_xy_us      (lon_points,lat_points), &! wind speed in eastward direction [m/s]
  !     forc_xy_vs      (lon_points,lat_points)   ! wind speed in northward direction [m/s]
  real(r8), INTENT(in) :: &
        forc_pco2m(lon_points,lat_points) , &! partial pressure of CO2 at observational height [pa]
        forc_po2m(lon_points,lat_points)  , &! partial pressure of O2 at observational height [pa]
        forc_us(lon_points,lat_points)    , &! wind speed in eastward direction [m/s]
        forc_vs(lon_points,lat_points)    , &! wind speed in northward direction [m/s]
        forc_t(lon_points,lat_points)     , &! temperature at agcm reference height [kelvin]
        forc_q(lon_points,lat_points)     , &! specific humidity at agcm reference height [kg/kg]
        forc_rh(lon_points,lat_points)    , &! relative humidity [0-1]
        forc_prc(lon_points,lat_points)   , &! convective precipitation [mm/s]
        forc_prl(lon_points,lat_points)   , &! large scale precipitation [mm/s]
        forc_psrf(lon_points,lat_points)  , &! atmosphere pressure at the surface [pa]
        forc_pbot(lon_points,lat_points)  , &! atmosphere pressure at the bottom of the atmos. model level [pa]
        forc_sols(lon_points,lat_points)  , &! atm vis direct beam solar rad onto srf [W/m2]
        forc_soll(lon_points,lat_points)  , &! atm nir direct beam solar rad onto srf [W/m2]
        forc_solsd(lon_points,lat_points) , &! atm vis diffuse solar rad onto srf [W/m2]
        forc_solld(lon_points,lat_points) , &! atm nir diffuse solar rad onto srf [W/m2]
        forc_frl(lon_points,lat_points)   , &! atmospheric infrared (longwave) radiation [W/m2]
        forc_hgt_u(lon_points,lat_points) , &! observational height of wind [m]
        forc_hgt_t(lon_points,lat_points) , &! observational height of temperature [m]
        forc_hgt_q(lon_points,lat_points) , &! observational height of humidity [m]
        forc_rhoair(lon_points,lat_points)   ! density air [kg/m3]

! Variables required for restart run
! ----------------------------------------------------------------------

  real(r8), INTENT(inout) :: oro(lon_points,lat_points)  ! ocean(0)/seaice(2)/ flag
  real(r8), INTENT(inout) :: &
        z_soisno(maxsnl+1:nl_soil,lon_points,lat_points)   , &! layer depth (m)
        dz_soisno(maxsnl+1:nl_soil,lon_points,lat_points)  , &! layer thickness (m)
        t_soisno(maxsnl+1:nl_soil,lon_points,lat_points)   , &! soil + snow layer temperature [K]
        wliq_soisno(maxsnl+1:nl_soil,lon_points,lat_points), &! liquid water (kg/m2)
        wice_soisno(maxsnl+1:nl_soil,lon_points,lat_points), &! ice lens (kg/m2)
        smf_soisno (maxsnl+1:nl_soil,lon_points,lat_points), &! snowmelt fraction in snow-soil (0-1).Added by HY,2016-11-26
        t_lake(nl_lake,lon_points,lat_points)       ,&! lake temperature (kelvin)
        lake_icefrac(nl_lake,lon_points,lat_points) ,&! lake mass fraction of lake layer that is frozen
        t_grnd(lon_points,lat_points)     , &! ground surface temperature [k]
        tlsun(lon_points,lat_points)      , &! sunlit leaf temperature [K]
        tlsha(lon_points,lat_points)      , &! shaded leaf temperature [K]
        ldew(lon_points,lat_points)       , &! depth of water on foliage [kg/m2/s]
        sag(lon_points,lat_points)        , &! non dimensional snow age [-]
        scv(lon_points,lat_points)        , &! snow mass (kg/m2)
        snowdp(lon_points,lat_points)     , &! snow depth (m)
        sst(lon_points,lat_points)        , &! surface water storage (mm)
        zwt(lon_points,lat_points)        , &! the depth to water table [m]
        wa(lon_points,lat_points)         , &! water storage in aquifer [mm]
        fveg(lon_points,lat_points)       , &! fraction of vegetation cover
        fsno(lon_points,lat_points)       , &! fractional snow cover
        sigf(lon_points,lat_points)       , &! fraction of veg cover, excluding snow-covered veg [-]
        green(lon_points,lat_points)      , &! greenness
        lai(lon_points,lat_points)        , &! leaf area index
        sai(lon_points,lat_points)        , &! stem area index
        coszen(lon_points,lat_points)     , &! cosine of solar zenith angle
        albg(2,2,lon_points,lat_points)  , &! albedo, ground [-]
        albv(2,2,lon_points,lat_points)  , &! albedo, vegetation [-]
        alb(2,2,lon_points,lat_points)   , &! averaged albedo [-]
        ssun(2,2,lon_points,lat_points)  , &! sunlit canopy absorption for solar radiation
        ssha(2,2,lon_points,lat_points)  , &! shaded canopy absorption for solar radiation
        thermk(lon_points,lat_points)     , &! canopy gap fraction for tir radiation
        extkb(lon_points,lat_points)      , &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd(lon_points,lat_points)         ! diffuse and scattered diffuse PAR extinction coefficient
        
  real(r8), INTENT(inout) :: & !Added for GBHM
        Drw_gbhm(nc,nx)           ! water depth in the river of the flow interval (m)        
     
! additional diagnostic variables for output
  real(r8), INTENT(out) :: &
        laisun(lon_points,lat_points)     , &! sunlit leaf area index
        laisha(lon_points,lat_points)     , &! shaded leaf area index
        rstfac(lon_points,lat_points)     , &! factor of soil water stress 
        wat(lon_points,lat_points)        , &! total water storage
        h2osoi(nl_soil,lon_points,lat_points)! volumetric soil water in layers [m3/m3]

! Fluxes
! ----------------------------------------------------------------------
  real(r8), INTENT(out) :: &
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
        tref(lon_points,lat_points)       , &! 2 m height air temperature [K]
        qref(lon_points,lat_points)       , &! 2 m height air specific humidity
        trad(lon_points,lat_points)       , &! radiative temperature [K]
        rsur(lon_points,lat_points)       , &! surface runoff (mm h2o/s)
        rnof(lon_points,lat_points)       , &! total runoff (mm h2o/s)
        qlat(lon_points,lat_points)       , &! lateral runoff (mm h2o/s)
        qintr(lon_points,lat_points)      , &! interception (mm h2o/s)
        qinfl(lon_points,lat_points)      , &! inflitration (mm h2o/s)
        qdrip(lon_points,lat_points)      , &! throughfall (mm h2o/s)
        qcharge(lon_points,lat_points)    , &! groundwater recharge [mm/s]
        rst(lon_points,lat_points)        , &! canopy stomatal resistance 
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
        srniln(lon_points,lat_points)     , &! reflected diffuse beam nir solar radiation at local noon(W/m2)
        forc_rain(lon_points,lat_points)  , &! rain [mm/s]
        forc_snow(lon_points,lat_points)  !, &! snow [mm/s]
        ! emis(lon_points,lat_points)       , &! averaged bulk surface emissivity
        ! z0ma(lon_points,lat_points)       , &! effective roughness [m]
        ! zol(lon_points,lat_points)        , &! dimensionless height (z/L) used in Monin-Obukhov theory
        ! rib(lon_points,lat_points)        !, &! bulk Richardson number in surface layer
        ! ustar(lon_points,lat_points)      , &! u* in similarity theory [m/s]
        ! qstar(lon_points,lat_points)      , &! q* in similarity theory [kg/kg]
        ! tstar(lon_points,lat_points)      , &! t* in similarity theory [K]
        ! fm(lon_points,lat_points)         , &! integral of profile function for momentum
        ! fh(lon_points,lat_points)         , &! integral of profile function for heat
        ! fq(lon_points,lat_points)            ! integral of profile function for moisture
    
  real(r8), intent(out) :: &                !Added for MODEL development
        qground    (lon_points,lat_points), &
        qsmelt     (lon_points,lat_points), &! total snowmelt flux (mm H2O /s)
        region_qsub(lon_points,lat_points), &
        net_blow   (lon_points,lat_points), &
        blowing_add_region! 
  ! GBHM module needed
  integer, INTENT(in) :: nc, startyear, endyear       
  character*6 , INTENT(in) ::  subbasin(nc)      ! name of sub-basin     
  real(r8), INTENT(in) :: &  
        area        (lon_points,lat_points) , &   ! area of the local grid (m2)
        dx          (nc,nx)     , &   ! flow interval length (m)
        Dr_gbhm     (nc,nx)     , &   ! river depth in the flow interval (m)
        wriver      (nc,nx)     , &   ! width of river (m)
        s0          (nc,nx)     , &   ! slope of river bed (ND)      
        mnroughness (nc,nx)           ! Manning's roughness

  integer, INTENT(in) :: & 
        run_gbhm           , &        ! run:1; Dont run:0
        nflow    (nc)      , &        ! total number of flow-intervals in a sub-basin
        psubbasin(nc)      , & 
        nbasinup (nc,4)    , &
        nsub               , &        ! total number of sub-basins
        ngrid   (nc,nx)    , &        ! number of grids in this flow-interval
        grid_row(nc,nx,np) , &        ! row of grids in this flow-interval
        grid_col(nc,nx,np) , &        ! column of grids in this flow-interval
        year               , &        ! calendar year   (19XX)
        month              , &        ! month in a year (1-12)
        day                , &        ! day in a month  (1-30)
        hour               , &        ! hour in a day   (1-24)
        dayinmonth(12)     , &        ! days of a month
        startmonth         , &        ! start month in the hydro-year
        endmonth           , &        ! end   month in the hydro-year
        startday           , &        ! start day in the first month of the hydro-year
        endday             , &        ! end   day in the last  month of the hydro-year
        idc                , &        ! contineous day in a year (1-366)
        ihc                , &        ! contineous hour in a year (1-366*24)
        start_sub, end_sub            ! for partial simulation of a basin

  integer, INTENT(inout) :: start     ! 0 - faulse,    1 - true

  real(r8), INTENT(out)   ::  q2   (nc,nx)       ! discharge of current time step (m^3/s)
  real(r8), INTENT(out)   ::  Qh   (nc,8800)     ! hourly mean discharge (m3/s)
  real(r8), INTENT(inout) ::  qr1  (nc,nx)       ! discharge of reservoir flowout last time step(m^3/s)
  real(r8), INTENT(inout) ::  qlin1(nc,nx)       ! lateral inflow of last time step (m^3/m/s)
  real(r8), INTENT(inout) ::  Qd   (nc,366)      ! daily average discharge (m^3/s)
  real(r8), INTENT(inout) ::  Qm   (nc,12)       ! monthly mean discharge (m3/s)
  real(r8), INTENT(inout) ::  q1   (nc,nx)       ! discharge of last time step (m^3/s)
  real(r8), INTENT(inout) ::  qr2  (nc,nx)       ! discharge of reservoir flowout current time step(m^3/s)
  real(r8), INTENT(out)   ::  Drw  (lon_points,lat_points)
!====================F2PY===========================
!f2py intent(in) run_pbsm,nl_soil,nl_lake,maxsnl,lon_points,lat_points,mask
!f2py intent(in) nc,startyear,endyear   
!f2py intent(in) idate,deltim,lons,dlon,dlat,ivt,itypwat,lakedepth ,dz_lake
!f2py intent(in) soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb 
!f2py intent(in) wfld,wsat,wrsd,watern,alpha,slope,length,Ds,Dr,Dg,kground  
!f2py intent(in) porsl,psi0,bsw,hksati,csol,dksatu,dkdry,rootfr,z0m    
!f2py intent(in) displa,sqrtdi,effcon,vmax25,slti,hlti,shti,hhti,trda
!f2py intent(in) trdm,trop,gradm,binter,extkn,chil,ref,tran,zlnd,zsno,adjfac             
!f2py intent(in) csoilc ,dewmx,wtfact ,capr,cnfac,ssi,wimp,pondmx           
!f2py intent(in) smpmax,smpmin,trsmx0,tcrit,forc_xy_q,forc_xy_t      
!f2py intent(in) forc_pco2m,forc_po2m,forc_us,forc_vs,forc_t,forc_q,forc_rh,forc_prc,forc_prl
!f2py intent(in) forc_psrf,forc_pbot,forc_sols,forc_soll,forc_solsd,forc_solld,forc_frl
!f2py intent(in) forc_hgt_u,forc_hgt_t,forc_hgt_q,forc_rhoair
! ----------------------------------------------------------------------------
!f2py intent(in,out)  z_soisno,dz_soisno,t_soisno,wliq_soisno,wice_soisno,sst
!f2py intent(in,out)  t_lake,lake_icefrac,oro,t_grnd ,tlsun,tlsha,ldew,sag,scv     
!f2py intent(in,out)  snowdp,zwt,wa,Drw_gbhm,fveg,fsno,sigf,green,lai,sai,coszen,albg,albv,alb
!f2py intent(in,out)  ssun,ssha,thermk,extkb,extkd
!f2py intent(out)     h2osoi,laisun,laisha ,rstfac,wat 
!f2py intent(out)     taux,tauy,fsena,fevpa,lfevpa,fsenl,fevpl,etr,fseng,fevpg,olrg,fgrnd  
!f2py intent(out)     xerr,zerr,tref,qref,trad,rsur,rnof,qlat,qground,qsmelt 
!f2py intent(out)     qintr,qinfl,qdrip,qcharge,rst,assim,respc,sabvsun,sabvsha,sabg,sr       
!f2py intent(out)     solvd ,solvi,solnd,solni,srvd,srvi,srnd ,srni,solvdln,solviln,solndln  
!f2py intent(out)     solniln ,srvdln ,srviln,srndln,srniln,forc_rain,forc_snow 
!f2py intent(out)     region_qsub, net_blow,blowing_add_region,Drw
!--------------------GBHM module-------------------------
!f2py intent(in) subbasin,nflow,psubbasin,nbasinup,nsub,ngrid,grid_row,grid_col,run_gbhm
!f2py intent(in) year,month,day,hour,dayinmonth,startmonth,endmonth,startday,endday
!f2py intent(in) idc,ihc,start_sub,end_sub,area,dx,Dr_gbhm,wriver,mnroughness,s0,ssf
!f2py intent(in,out) start,q1,qr1,qr2,qlin1,Qd,Qm
!f2py intent(out) q2,Qh
!--------------------------------------------------------
!f2py depend(lon_points) lons   
!f2py depend(lon_points,lat_points) region_qsub, net_blow,dlon,dlat,ivt,mask,itypwat,lakedepth 
!f2py depend(lon_points,lat_points) soil_s_v_alb,soil_d_v_alb, soil_s_n_alb,soil_d_n_alb 
!f2py depend(lon_points,lat_points) wsat,wrsd,watern,alpha,slope,length,Ds,Dr,Dg,kground,Drw
!f2py depend(lon_points,lat_points) z0m,displa,sqrtdi,effcon,vmax25,slti,hlti,shti,hhti,trda,trdm,trop
!f2py depend(lon_points,lat_points) gradm,binter,extkn,chil,ref,tran
!f2py depend(lon_points,lat_points) forc_pco2m,forc_po2m,forc_us,forc_vs,forc_t,forc_q,forc_rh,forc_prc,forc_prl
!f2py depend(lon_points,lat_points) forc_psrf,forc_pbot,forc_sols,forc_soll,forc_solsd,forc_solld,forc_frl
!f2py depend(lon_points,lat_points) forc_hgt_u,forc_hgt_t,forc_hgt_q,forc_rhoair,oro
!f2py depend(lon_points,lat_points) t_grnd,tlsun,tlsha,ldew,sag,scv,snowdp,zwt,wa,sst,ssf 
!f2py depend(lon_points,lat_points) fveg,fsno,sigf,green,lai,sai,coszen 
!f2py depend(lon_points,lat_points) albg,albv,alb,ssun,ssha,thermk
!f2py depend(lon_points,lat_points) extkb,extkd,laisun,laisha,rstfac  
!f2py depend(lon_points,lat_points) wat,taux,tauy,fsena,fevpa,lfevpa ,fsenl,fevpl   
!f2py depend(lon_points,lat_points) etr,fseng,fevpg,olrg,fgrnd,xerr,zerr,tref,qref
!f2py depend(lon_points,lat_points) trad,rsur,rnof,qlat,qground,qsmelt,qintr,qinfl,qdrip,qcharge,rst,assim,respc   
!f2py depend(lon_points,lat_points) sabvsun,sabvsha,sabg,sr,solvd,solvi,solnd,solni  
!f2py depend(lon_points,lat_points) srvd,srvi,srnd,srni,solvdln,solviln,solndln,solniln,srvdln,srviln,srndln,srniln
!f2py depend(lon_points,lat_points) forc_rain,forc_snow
!f2py depend(nl_soil,lon_points,lat_points)        wfld,h2osoi,porsl, psi0,bsw,hksati,csol,dksatu,dkdry,rootfr
!f2py depend(maxsnl,nl_soil,lon_points,lat_points) z_soisno,dz_soisno,t_soisno,wliq_soisno,wice_soisno
!f2py depend(nl_lake,lon_points,lat_points)        dz_lake,t_lake,lake_icefrac      
!-------------------GBHM module--------------------------
!f2py depend(lon_points,lat_points) area
!f2py depend(nc) subbasin,nflow,psubbasin 
!f2py depend(nc,4) nbasinup
!f2py depend(nc,8800) Qh
!f2py depend(nc,366)  Qd
!f2py depend(nc,12)   Qm 
!f2py depend(nc,nx)   dx,Dr_gbhm,Drw_gbhm,wriver,s0,mnroughness,ngrid,q2,qr1,qlin1,q1,qr2 
!f2py depend(nc,nx,np) grid_row,grid_col   
!===================END F2PY========================        

! -------------- Local varaibles -------------------
  integer :: i,j
  real(r8):: temppass , &
        !For blowing snow module
        region_grainsize(lon_points,lat_points), &
        region_fetch(lon_points,lat_points)
  ! GBHM needed      
  real(r8)  runoff    (lon_points,lat_points) , &   ! surface runoff
      runoff_inter(lon_points,lat_points) , &   ! Subsurface flow In each grid
      runoff_g    (lon_points,lat_points) , &   ! exchanges between groundwater and river In each grid
      slope_gbhm  (lon_points,lat_points) , &
      length_gbhm (lon_points,lat_points) , &
      dl_sst      (lon_points,lat_points) , &
      dl_Drw      (lon_points,lat_points)
  integer       isub,iflow,ig             ! sub-basin number      

! -----------------------------------------------------------------------------------------
! Run blowing snow module
! -----------------------------------------------------------------------------------------
      if (run_pbsm == 1) then
        region_grainsize(:,:) = 0.0001 !m. just assumed.
        region_fetch    (:,:) = 300.   !m
        call blowing_snow(lon_points,lat_points,nl_soil,maxsnl,deltim,forc_hgt_u,& ! input
            forc_us,forc_vs,forc_t,forc_rh,region_grainsize,region_fetch,&         ! input
            t_soisno,wliq_soisno,wice_soisno,dz_soisno,z_soisno,scv, snowdp,&      ! in-out 
            region_qsub,net_blow,blowing_add_region)
      endif 

      if (run_gbhm == 1) then
        dl_Drw(:,:) =0.0
        do isub = 1, nsub
          do iflow = 1, nflow(isub)
            do ig = 1,ngrid(isub,iflow)      
              dl_Drw(grid_row(isub,iflow,ig),grid_col(isub,iflow,ig)) = Drw_gbhm(isub,iflow) 
              Drw(:,:) = dl_Drw(lon_points:1:-1,:)
            enddo
          enddo
        enddo
      else
        Drw(:,:) =0.0
      endif 
! -----------------------------------------------------------------------------------------
! CLM model at one dimension
! -----------------------------------------------------------------------------------------  
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
      DO j = 1, lat_points                                        
        DO i = 1, lon_points 
         if (mask(i,j)==0) then
            rsur(i,j) = 0.
            qlat(i,j) = 0.        
            qground(i,j) = 0.
         else
            CALL CLMMAIN (i,j,deltim,oro(i,j), nl_soil, nl_lake, maxsnl,&
          ! The following inputs are defined in common module variables
            lons(i), dlon(i,j), dlat(i,j), ivt(i,j), itypwat(i,j), &
            lakedepth(i,j), dz_lake(1:,i,j), &                     ! new lake scheme
    
          ! SOIL INFORMATION AND LAKE DEPTH
            soil_s_v_alb(i,j), soil_d_v_alb(i,j), soil_s_n_alb(i,j), soil_d_n_alb(i,j), &
            wsat(i,j),wrsd(i,j),watern(i,j),alpha(i,j),slope(i,j),length(i,j),  &!GBHM parameters
            Ds(i,j),           Dr(i,j),           Dg(i,j),           kground(i,j),&
            porsl(1:,i,j),     psi0(1:,i,j),      bsw(1:,i,j),       hksati(1:,i,j),    &
            wfld(1:,i,j),      adjfac,&
            csol(1:,i,j),      dksatu(1:,i,j),    dkdry(1:,i,j),     rootfr(1:,i,j),    &
    
          ! VEGETATION INFORMATION
            z0m(i,j),          displa(i,j),       sqrtdi(i,j),                          &
            effcon(i,j),       vmax25(i,j),       slti(i,j),         hlti(i,j),         &
            shti(i,j),         hhti(i,j),         trda(i,j),         trdm(i,j),         &
            trop(i,j),         gradm(i,j),        binter(i,j),       extkn(i,j),        &
            chil(i,j),         ref(1:,1:,i,j),    tran(1:,1:,i,j),                      &
    
          ! ATMOSPHERIC FORCING
            forc_pco2m(i,j),   forc_po2m(i,j),    forc_us(i,j),      forc_vs(i,j),      &
            forc_t(i,j),       forc_q(i,j),       forc_prc(i,j),     forc_prl(i,j),     &
            forc_rain(i,j),    forc_snow(i,j),    forc_psrf(i,j),    forc_pbot(i,j),    &
            forc_sols(i,j),    forc_soll(i,j),    forc_solsd(i,j),   forc_solld(i,j),   &
            forc_frl(i,j),     forc_hgt_u(i,j),   forc_hgt_t(i,j),   forc_hgt_q(i,j),   &
            forc_rhoair(i,j),                                                           &
    
          ! LAND SURFACE VARIABLES REQUIRED FOR RESTART
            idate,                                                                      &
            z_soisno(maxsnl+1:,i,j),            dz_soisno(maxsnl+1:,i,j),               &
            t_soisno(maxsnl+1:,i,j),            wliq_soisno(maxsnl+1:,i,j),             &
            wice_soisno(maxsnl+1:,i,j),                                                 &
            sst(i,j),          ssf(i,j),                                                &
            t_grnd(i,j),       tlsun(i,j),        tlsha(i,j),        ldew(i,j),         &
            sag(i,j),          scv(i,j),          snowdp(i,j),       fveg(i,j),         &
            fsno(i,j),         sigf(i,j),         green(i,j),        lai(i,j),          &
            sai(i,j),          coszen(i,j),       albg(1:,1:,i,j),   albv(1:,1:,i,j),   &
            alb(1:,1:,i,j),    ssun(1:,1:,i,j),   ssha(1:,1:,i,j),   thermk(i,j),       &
            extkb(i,j),        extkd(i,j),        &
           
            zwt(i,j),          wa(i,j),                                             &
            t_lake(1:,i,j),    lake_icefrac(1:,i,j),                                & ! new lake scheme
    
          ! additional diagnostic variables for output
            laisun(i,j),       laisha(i,j),                                         &
            Drw(i,j),          qground(i,j),      qsmelt(i,j),                          &
            rstfac(i,j),       h2osoi(1:,i,j),    wat(i,j),          qlat(i,j),         &
    
          ! FLUXES
            taux(i,j),         tauy(i,j),         fsena(i,j),        fevpa(i,j),        &
            lfevpa(i,j),       fsenl(i,j),        fevpl(i,j),        etr(i,j),          &
            fseng(i,j),        fevpg(i,j),        olrg(i,j),         fgrnd(i,j),        &
            trad(i,j),         tref(i,j),         qref(i,j),         rsur(i,j),         &
            rnof(i,j),         qintr(i,j),        qinfl(i,j),        qdrip(i,j),        &
            rst(i,j),          assim(i,j),        respc(i,j),        sabvsun(i,j),      &
            sabvsha(i,j),      sabg(i,j),         sr(i,j),           solvd(i,j),        &
            solvi(i,j),        solnd(i,j),        solni(i,j),        srvd(i,j),         &
            srvi(i,j),         srnd(i,j),         srni(i,j),         solvdln(i,j),      &
            solviln(i,j),      solndln(i,j),      solniln(i,j),      srvdln(i,j),       &
            srviln(i,j),       srndln(i,j),       srniln(i,j),       qcharge(i,j),      &
            xerr(i,j),         zerr(i,j),                                           &
    
          ! TUNABLE modle constants
            zlnd,            zsno,            csoilc,          dewmx,           &
            wtfact,          capr,            cnfac,           ssi,             &
            wimp,            pondmx,          smpmax,          smpmin,          &
            trsmx0,          tcrit)!,                                             &
    
          ! additional variables required by coupling with WRF model
            ! emis(i,j),         z0ma(i,j),         zol(i,j),          rib(i,j))!,        &
            ! ustar(i,j),        qstar(i,j),        tstar(i,j),                         &
            ! fm(i,j),           fh(i,j),           fq(i,j))
            ! z_sno (maxsnl+1:0,i,j) = z_soisno (maxsnl+1:0,i,j)
            ! dz_sno(maxsnl+1:0,i,j) = dz_soisno(maxsnl+1:0,i,j)
          endif 
        ENDDO
      ENDDO
!$OMP END PARALLEL DO                          

! -----------------------------------------------------------------------------------------
! Run GBHM module, 2016-3-16
! Common variables with CLM module: Drw, sst
!   note: the array structure of Drw and sst is different in CLM from GBHM      
! INPUT:
! OUTPUT:      
! -----------------------------------------------------------------------------------------

      if (run_gbhm == 1) then
        ! from mm/s to m/s, to meet the unit used in GBHM
        runoff      (:,:) = rsur   (lon_points:1:-1,:)*0.001 
        runoff_inter(:,:) = qlat   (lon_points:1:-1,:)*0.001
        runoff_g    (:,:) = qground(lon_points:1:-1,:)*0.001        
        dl_sst      (:,:) = sst    (lon_points:1:-1,:)
        slope_gbhm  (:,:) = slope  (lon_points:1:-1,:)
        length_gbhm (:,:) = length (lon_points:1:-1,:)

        do isub = 1, nsub
          if (start == 1) then 
            call hillslope_route_model(lon_points,lat_points,nc,subbasin,isub,nflow,&
                psubbasin,nbasinup,nsub,ngrid,grid_row,grid_col,start, &
                year,month,day,hour,dayinmonth,startyear,endyear,startmonth,endmonth,&
                startday,endday, deltim, &
                idc,ihc,start_sub,end_sub,runoff_inter,runoff_g,runoff,slope_gbhm,length_gbhm,&
                area,dx,Dr_gbhm,s0,wriver,mnroughness,&
                q1,q2,qr1,qr2,qlin1,Qh,Qd,Qm,Drw_gbhm,dl_sst)
            start = 0
          endif  
          if (start == 0) then 
            call hillslope_route_model(lon_points,lat_points,nc,subbasin,isub,nflow,&
                psubbasin,nbasinup,nsub,ngrid,grid_row,grid_col,start, &
                year,month,day,hour,dayinmonth,startyear,endyear,startmonth,endmonth,&
                startday,endday, deltim, &
                idc,ihc,start_sub,end_sub,runoff_inter,runoff_g,runoff,slope_gbhm,length_gbhm,&
                area,dx,Dr_gbhm,s0,wriver,mnroughness,&
                q1,q2,qr1,qr2,qlin1,Qh,Qd,Qm,Drw_gbhm,dl_sst)
          endif
        enddo  
      endif

END SUBROUTINE WATER
! --------- EOP --------