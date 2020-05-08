module initial
!-----------------------------------------------------------------------
  use precision
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: initialize
  public :: IniTimeVar
  public :: snow_ini
  public :: polint
!-----------------------------------------------------------------------

  CONTAINS

SUBROUTINE initialize (  &
           idate,lon_points,lat_points, nl,                 &
           nl_soil,maxsnl,nl_lake,mxy_patch, latixy,longxy, &
           snow_d_grid,lakedepth,&
           gravel_grid,sand_grid,clay_grid,SOC_grid,BD_grid,&
           soil_t_grid,soil_w_grid,&
           lons, dlon, dlat, itypwat, &
           dz_lake, &
           soil_s_v_alb, soil_d_v_alb, soil_s_n_alb, soil_d_n_alb,  &
           porsl_grid,   psi0_grid,    bsw_grid,     hksati_grid,   &
           csol_grid,    dksatu_grid,  dkdry_grid,   rootfr_grid,   &
           z0m,          displa,       sqrtdi,                      &
           effcon,       vmax25,       slti,         hlti,          &
           shti,         hhti,         trda,         trdm,          &
           trop,         gradm,        binter,       extkn,         &
           chil,         ref,          tran,                        &
           z_soisno,     dz_soisno,    t_soisno,     wliq_soisno,   &
           wice_soisno,  smf_soisno,   smf_ice,&
           t_grnd,       tlsun,        tlsha,        ldew,          &
           sag,          scv,          snowdp,       fveg,          &
           fsno,         sigf,         green,        lai,           &
           sai,          coszen,       albg,         albv,          &
           alb,          ssun,         ssha,         thermk,        &
           extkb,        extkd,                                     &
           zwt,          wa,                                        &
           t_lake,       lake_icefrac,                              &
           qref,         rst,          trad,         tref,          &
           zlnd,         zsno,         csoilc,       dewmx,         &
           wtfact,       capr,         cnfac,        ssi,           &
           wimp,         pondmx,       smpmax,       smpmin,        &
           trsmx0,       tcrit,                                     &
           emis,         z0ma,         zol,          rib,           &
           ustar,        qstar,        tstar,                       &
           fm,           fh,           fq )


! ======================================================================
! initialization routine for land surface model.
!
! Created and revised by Yongjiu Dai, 09/15/1999,08/30/2002,03/2014
! Major Revised by Hongyi  Li ,       09/2015
! ======================================================================
   use precision
   use PhysicalConstants
   use timemanager
   IMPLICIT NONE

! ----------------------------------------------------------------------
   integer, INTENT(in)    :: lon_points ! number of longitude points on model grid
   integer, INTENT(in)    :: lat_points ! number of latitude points on model grid
   integer, INTENT(in)    :: nl         ! number of total ground layers,30
   integer, INTENT(in)    :: nl_soil    ! number of soil layers
   integer, INTENT(in)    :: maxsnl     ! max number of snow layers
   integer, INTENT(in)    :: nl_lake    ! number of land water bodies' layers
   integer, INTENT(in)    :: idate(3)   ! year, julian day, seconds of the starting time
   integer, INTENT(in)    :: mxy_patch(lon_points,lat_points)          !index for vegetation type [-]
   real(r8),INTENT(in)    :: latixy   (lon_points,lat_points)          ! latitude in degree
   real(r8),INTENT(in)    :: longxy   (lon_points,lat_points)          ! longitude in degree
   real(r8),INTENT(in)    :: gravel_grid(nl,lon_points,lat_points)
   real(r8),INTENT(in)    :: sand_grid  (nl,lon_points,lat_points)
   real(r8),INTENT(in)    :: clay_grid  (nl,lon_points,lat_points)
   real(r8),INTENT(in)    :: SOC_grid   (nl,lon_points,lat_points)
   real(r8),INTENT(in)    :: BD_grid    (nl,lon_points,lat_points)
   real(r8),INTENT(in)    :: soil_t_grid(nl,lon_points,lat_points)! soil layer temperature (K)
   real(r8),INTENT(in)    :: soil_w_grid(nl,lon_points,lat_points)! soil layer wetness (-)
   real(r8),INTENT(in)    :: lakedepth  (lon_points,lat_points)        ! lake depth (m)
   real(r8),INTENT(in)    :: snow_d_grid(lon_points,lat_points)        ! snow depth (m)
   real(r8), INTENT(out) :: &
        lons   (lon_points),                       &! longitude in degree
        dlat   (lon_points,lat_points),            &! latitude in radians
        dlon   (lon_points,lat_points),            &! longitude in radians
        dz_lake(nl_lake,lon_points,lat_points)      ! lake layer thickness (m)
   integer, INTENT(out) ::  itypwat(lon_points,lat_points)
   real(r8), INTENT(out) :: &! soil physical parameters and lake info
        soil_s_v_alb(lon_points,lat_points)   , &! albedo of visible of the saturated soil
        soil_d_v_alb(lon_points,lat_points)   , &! albedo of visible of the dry soil
        soil_s_n_alb(lon_points,lat_points)   , &! albedo of near infrared of the saturated soil
        soil_d_n_alb(lon_points,lat_points)   , &! albedo of near infrared of the dry soil
        porsl_grid  (nl,lon_points,lat_points)   , &! fraction of soil that is voids [-]
        psi0_grid   (nl,lon_points,lat_points)   , &! minimum soil suction [mm]
        bsw_grid    (nl,lon_points,lat_points)   , &! clapp and hornbereger "b" parameter [-]
        hksati_grid (nl,lon_points,lat_points)   , &! hydraulic conductivity at saturation[mm h2o/s]
        csol_grid   (nl,lon_points,lat_points)   , &! heat capacity of soil solids [J/(m3 K)]
        dksatu_grid (nl,lon_points,lat_points)   , &! thermal conductivity of saturated soil [W/m-K]
        dkdry_grid  (nl,lon_points,lat_points)   , &! thermal conductivity for dry soil  [J/(K s m)]
        rootfr_grid (nl_soil,lon_points,lat_points)      ! fraction of roots in each soil layer
   real(r8), INTENT(out) :: &        ! vegetation static, dynamic, derived parameters
        fveg  (lon_points,lat_points)       , &! fraction of vegetation cover
        green (lon_points,lat_points)       , &! greenness
        lai   (lon_points,lat_points)       , &! leaf area index
        sai   (lon_points,lat_points)       , &! stem area index
        z0m   (lon_points,lat_points)       , &! aerodynamic roughness length [m]
        displa(lon_points,lat_points)       , &! displacement height [m]
        sqrtdi(lon_points,lat_points)       , &! inverse sqrt of leaf dimension [m**-0.5]
        effcon(lon_points,lat_points)       , &! quantum efficiency of RuBP regeneration(mol CO2/mol )
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
        tran  (2,2,lon_points,lat_points)      ! leaf transmittance (iw=iband, il=life and dead)
   real(r8), INTENT(out) :: &        ! tunable parameters
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
        tcrit         !critical temp. to determine rain or snow
   real(r8), INTENT(out) :: &
        z_soisno(maxsnl+1:nl,lon_points,lat_points)   , &! layer depth (m)
        dz_soisno(maxsnl+1:nl,lon_points,lat_points)  , &! layer thickness (m)
        t_soisno(maxsnl+1:nl,lon_points,lat_points)   , &! soil + snow layer temperature [K]
        wliq_soisno(maxsnl+1:nl,lon_points,lat_points), &! liquid water (kg/m2)
        wice_soisno(maxsnl+1:nl,lon_points,lat_points), &! ice lens (kg/m2)
        smf_soisno (maxsnl:nl_soil+1,lon_points,lat_points), &! snowmelt fraction in snow-soil (0-1).Added by HY,2016-11-26
        smf_ice (1:nl_soil,lon_points,lat_points), &! snowmelt fraction in soilice (0-1).Added by HY,2016-11-26
        t_lake(nl_lake,lon_points,lat_points)       ,&! lake temperature (kelvin)
        lake_icefrac(nl_lake,lon_points,lat_points) ,&! lake mass fraction of lake layer that is frozen
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
        extkd(lon_points,lat_points)         ! diffuse and scattered diffuse PAR extinction coefficient
   real(r8), INTENT(out) :: &         ! Additional variables required by reginal model (WRF & RSM)
        trad(lon_points,lat_points),                   &! radiative temperature of surface [K]
        tref(lon_points,lat_points),                   &! 2 m height air temperature [kelvin]
        qref(lon_points,lat_points),                   &! 2 m height air specific humidity
        rst(lon_points,lat_points),                    &! canopy stomatal resistance (s/m)
        emis(lon_points,lat_points),                   &! averaged bulk surface emissivity
        z0ma(lon_points,lat_points),                   &! effective roughness [m]
        zol(lon_points,lat_points),                    &! dimensionless height (z/L)
        rib(lon_points,lat_points),                    &! bulk Richardson number in surface layer
        ustar(lon_points,lat_points),                  &! u* in similarity theory [m/s]
        qstar(lon_points,lat_points),                  &! q* in similarity theory [kg/kg]
        tstar(lon_points,lat_points),                  &! t* in similarity theory [K]
        fm(lon_points,lat_points),                     &! integral of profile function for momentum
        fh(lon_points,lat_points),                     &! integral of profile function for heat
        fq(lon_points,lat_points)                       ! integral of profile function for moisture
!====================F2PY===========================
!f2py intent(in)   lon_points,lat_points,nl,nl_soil,maxsnl,nl_lake,idate,mxy_patch,latixy,longxy,gravel_grid
!f2py intent(in)   sand_grid,clay_grid,SOC_grid,BD_grid,soil_t_grid
!f2py intent(in)   soil_w_grid, lakedepth,snow_d_grid
!f2py intent(out)  lons ,dlat ,dlon ,dz_lake,itypwat
!f2py intent(out)  soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb
!f2py intent(out)  porsl_grid,psi0_grid,bsw_grid ,hksati_grid,csol_grid,dksatu_grid,dkdry_grid,rootfr_grid
!f2py intent(out)  fveg,green ,lai ,sai
!f2py intent(out)  z0m ,displa,sqrtdi,effcon,vmax25
!f2py intent(out)  slti,hlti,shti,hhti,trda,trdm,trop,gradm,binter,extkn,chil,ref,tran,zlnd,zsno,csoilc
!f2py intent(out)  dewmx,wtfact,capr,cnfac,ssi,wimp,pondmx,smpmax,smpmin,trsmx0,tcrit
!f2py intent(out)  z_soisno,dz_soisno,t_soisno,wliq_soisno,wice_soisno,smf_soisno,smf_ice
!f2py intent(out)  t_lake,lake_icefrac,t_grnd,tlsun,tlsha,ldew
!f2py intent(out)  sag,scv,snowdp,zwt,wa,fsno,sigf,coszen
!f2py intent(out)  albg,albv,alb,ssun,ssha,thermk,extkb,extkd,trad,tref
!f2py intent(out)  qref,rst,emis,z0ma,zol,rib,ustar,qstar,tstar,fm,fh,fq
!--------------------------------------------------------
!f2py depend(lon_points,lat_points) mxy_patch,latixy,longxy,lakedepth,snow_d_grid,dlat,dlon,itypwat
!f2py depend(lon_points,lat_points) soil_s_v_alb,soil_d_v_alb,soil_s_n_alb, soil_d_n_alb
!f2py depend(lon_points,lat_points) fveg,green,lai,sai
!f2py depend(lon_points,lat_points) z0m,displa,sqrtdi,effcon,vmax25,slti,hlti,shti,hhti,trda,trdm,trop
!f2py depend(lon_points,lat_points) gradm,binter,extkn,chil,ref,tran,t_grnd,tlsun,tlsha,ldew
!f2py depend(lon_points,lat_points) sag,scv,snowdp,zwt,wa,fsno,sigf
!f2py depend(lon_points,lat_points) coszen,albg,albv,alb,ssun,ssha,thermk,extkb,extkd,trad,tref,qref
!f2py depend(lon_points,lat_points) rst,emis,z0ma,zol,rib,ustar,qstar,tstar,fm,fh,fq
!f2py depend(nl_soil,lon_points,lat_points) rootfr_grid,smf_ice ,smf_soisno
!f2py depend(nl,lon_points,lat_points) gravel_grid,sand_grid,clay_grid,SOC_grid,BD_grid
!f2py depend(nl,lon_points,lat_points) soil_t_grid,soil_w_grid
!f2py depend(nl,lon_points,lat_points) porsl_grid,psi0_grid,bsw_grid,hksati_grid,csol_grid,dksatu_grid,dkdry_grid
!f2py depend(maxsnl,nl,lon_points,lat_points) z_soisno,  dz_soisno,t_soisno,wliq_soisno,wice_soisno
!f2py depend(nl_lake,lon_points,lat_points) t_lake,lake_icefrac,dz_lake
!===================END F2PY========================
! ------------------------ local variables -----------------------------
! surface classification and soil information
  real(r8)  :: zsoi(1:nl+2)      ! soil layer depth [m]
  real(r8)  :: dzsoi(1:nl)       ! soil node thickness [m]
  real(r8)  :: zsoih(0:nl+1)     ! interface level below a zsoi level [m]
  real(r8)  :: latdeg  (lat_points)   ! latitude in degree
  real(r8)  :: londeg  (lon_points)   ! longitude in degree
  real(r8)  :: calday                 ! Julian cal day (1.xx to 365.xx)
  real(r8)  :: pi                     ! pie
  real(r8)  :: deltash                ! tmpororay variable for soil depth calculation
  integer   :: i,j,k,l,m,nsl          ! indices
  integer   :: custom_soil            ! define customed soil layer division or default exponential increase used in CoLM

      custom_soil = 0  !use customed soil depth(5-10cm/layer)
! ----------------------------------
! 1. Initialize TUNABLE constants
! ----------------------------------
      zlnd   = 0.01    !Roughness length for soil [m]
      zsno   = 0.0024  !Roughness length for snow [m], default:0.0024
      csoilc = 0.004   !Drag coefficient for soil under canopy [-]
      dewmx  = 0.1     !maximum dew
      wtfact = 0.38    !Maximum saturated fraction (global mean; see Niu et al., 2005)
      capr   = 0.34    !Tuning factor to turn first layer T into surface T
      cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
      ssi    = 1.   !Irreducible water saturation of snow,0.033
      wimp   = 0.05    !Water impremeable if porosity less than wimp
      pondmx = 10.0    !Ponding depth (mm)
      smpmax = -1.5e5  !Wilting point potential in mm
      smpmin = -1.e8   !Restriction for min of soil poten. (mm)
      trsmx0 = 2.e-4   !Max transpiration for moist soil+100% veg. [mm/s]
      tcrit  = 0.      !critical temp. to determine rain or snow
! ---------------------------------------------------------------
!2. INITIALIZE TIME INVARIANT VARIABLES
! ---------------------------------------------------------------
! 2.1 Define the soil and lake layers's thickness
! Modified by HY, to meet the layer rule in CLM.
! [0-45,45-91,91-166,166-289,289-493,493-829,829-1383] (\mm)
! ...............................................................

      do nsl = 1, nl+2
         zsoi(nsl) = 0.025*(exp(0.5*(nsl-0.5))-1.)  ! node depths
      end do

      zsoih(0) = 0.
      do nsl = 1, nl+1
         zsoih(nsl) = 0.5*(zsoi(nsl)+zsoi(nsl+1))   ! interface depths
      enddo
      zsoih(0:8) = (/0.,0.0451,0.0906,0.1655,0.289,0.493,0.829,1.383,2.5/)
if(custom_soil==1) then
      ! same as CLM in the first 3 soil layers
      zsoih(0:4) = (/0.,0.0175,0.0451,0.0906,0.1655/)
      zsoih(5:23) = (/0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,&
                     1.6,1.7,1.8,1.9,2.0/)
      zsoih(24:34) = (/2.2,2.6,3.0,3.5,4.5,6.0,8.0,11.0,15.0,20.0,30.0/)
      zsoih(35:38) = (/50.0,80.0,120.0,150./)
endif

      ! dzsoi(1) = zsoih(1)
      do nsl = 1, nl
         dzsoi(nsl) = zsoih(nsl)-zsoih(nsl-1) ! thickness b/n two interfaces
      end do
      ! zsoi(1)  =  0.5*dzsoi(1)
      do nsl = 1, nl
         zsoi(nsl) =0.5*(zsoih(nsl)+zsoih(nsl-1))  ! node depths
      end do

! ...............................................................
! 2.2 Plant time-invariant variables (based on the look-up tables)
! ...............................................................
      do i = 1, lon_points
        do j = 1,lat_points
          CALL IniTimeConst(nl_soil,zsoih,mxy_patch(i,j),itypwat(i,j),&
               z0m(i,j),displa(i,j),sqrtdi(i,j),effcon(i,j),vmax25(i,j),slti(i,j),hlti(i,j),&
               shti(i,j),hhti(i,j),trda(i,j),trdm(i,j),trop(i,j),gradm(i,j),binter(i,j),extkn(i,j),&
               chil(i,j),ref(1:,1:,i,j),tran(1:,1:,i,j),rootfr_grid(1:,i,j))
         ! soil layer division at grid scale
         z_soisno (1:nl ,i, j) = zsoi (1:nl)
         dz_soisno(1:nl ,i, j) = dzsoi(1:nl)
        enddo
      enddo

      ! get grid latitudes and longitudes
      latdeg = latixy(1,:)
      londeg = longxy(:,1)
      lons   = londeg

      ! convert latitudes and longitudes from degress to radians
      pi = 4.*atan(1.)
      dlat(:,:) = latixy(:,:)*pi/180.
      dlon(:,:) = longxy(:,:)*pi/180.

! ................................
! 2.3 cosine of solar zenith angle
! ................................
      call calendarday_func(idate, londeg(1),calday)
      do i = 1, lon_points
        do j = 1,lat_points
          call orb_coszen_func(calday,dlon(i,j),dlat(i,j),coszen(i,j))
        enddo
      enddo

! ----------------------------------------------------------------------
! 3. INITIALIZE TIME-VARYING VARIABLES
! ----------------------------------------------------------------------
! ...................
! 3.1 LEAF area index
! ...................

                ! porsl (nsl,npatch) =    soil_theta_s_l  (m,i,j)               ! cm/cm
                ! psi0  (nsl,npatch) =    soil_psi_s_l    (m,i,j) * 10.         ! cm -> mm
                ! bsw   (nsl,npatch) = 1./soil_lambda_l(m,i,j)                  ! dimensionless
                ! hksati(nsl,npatch) =    soil_k_s_l      (m,i,j) * 10./86400.  ! cm/day -> mm/s
                ! csol  (nsl,npatch) =    soil_csol_l     (m,i,j)               ! J/(m2 K)
                ! dksatu(nsl,npatch) =    soil_tksatu_l   (m,i,j)               ! W/(m K)
                ! dkdry (nsl,npatch) =    soil_tkdry_l    (m,i,j)               ! W/(m K)
    ! CREAT fraction of vegetation cover, greenness, leaf area index, stem index
      lai(:,:)=0.0; sai(:,:)=0.0; green(:,:)=0.0; fveg(:,:)=0.0
      do i = 1, lon_points
        do j = 1,lat_points
         do l = 1, nl
            t_soisno(l,i,j) = soil_t_grid(l,i,j)
            CALL soil_hydraulic_parameters( sand_grid(l,i,j),clay_grid(l,i,j),&
                          SOC_grid(l,i,j),BD_grid(l,i,j),zsoi(l), &
                          porsl_grid(l,i,j),psi0_grid(l,i,j),bsw_grid(l,i,j),hksati_grid(l,i,j))
            CALL soil_thermal_parameter(gravel_grid(l,i,j),sand_grid(l,i,j),&
                          clay_grid(l,i,j),SOC_grid(l,i,j),BD_grid(l,i,j),porsl_grid(l,i,j),zsoi(l),&
                                   csol_grid(l,i,j),dksatu_grid(l,i,j),dkdry_grid(l,i,j))
         enddo
       ! Call Ecological Model()
         if(mxy_patch(i,j)>0 .and. mxy_patch(i,j)<25)then
            CALL lai_empirical(mxy_patch(i,j),nl_soil,rootfr_grid(1:nl_soil,i,j),&
                               t_soisno(1:nl_soil,i,j),lai(i,j),sai(i,j),fveg(i,j),green(i,j))
         elseif (mxy_patch(i,j)>=25)then
            CALL lai_empirical(7,nl_soil,rootfr_grid(1:nl_soil,i,j),&
                               t_soisno(1:nl_soil,i,j),lai(i,j),sai(i,j),fveg(i,j),green(i,j))
         endif
        enddo
      enddo


      bsw_grid    = 1./bsw_grid               ! dimensionless
      psi0_grid   = psi0_grid * 10.           ! cm -> mm
      hksati_grid = hksati_grid * 10./86400.  ! cm/day -> mm/s
      do i = 1, lon_points
        do j = 1,lat_points
          if(mxy_patch(i,j)>0 .and. mxy_patch(i,j)<25)then
            CALL soil_color_refl(mxy_patch(i,j) &
              ,soil_s_v_alb(i,j),soil_d_v_alb(i,j),soil_s_n_alb(i,j),soil_d_n_alb(i,j))
          elseif (mxy_patch(i,j)>=25)then
            CALL soil_color_refl(7 &
              ,soil_s_v_alb(i,j),soil_d_v_alb(i,j),soil_s_n_alb(i,j),soil_d_n_alb(i,j))
          endif
          CALL iniTimeVar(nl,maxsnl,itypwat(i,j)&
              ,porsl_grid(1:,i,j)&
              ,soil_s_v_alb(i,j),soil_d_v_alb(i,j),soil_s_n_alb(i,j),soil_d_n_alb(i,j)&
              ,nl,zsoi,soil_t_grid(1:,i,j),soil_w_grid(1:,i,j),snow_d_grid(i,j)&
              ,z0m(i,j),zlnd,chil(i,j),ref(1:,1:,i,j),tran(1:,1:,i,j)&
              ,z_soisno(maxsnl+1:,i,j),dz_soisno(maxsnl+1:,i,j)&
              ,t_soisno(maxsnl+1:,i,j),wliq_soisno(maxsnl+1:,i,j),wice_soisno(maxsnl+1:,i,j)&
              ,zwt(i,j),wa(i,j)&
              ,t_grnd(i,j),tlsun(i,j),tlsha(i,j),ldew(i,j),sag(i,j),scv(i,j)&
              ,snowdp(i,j),fveg(i,j),fsno(i,j),sigf(i,j),green(i,j),lai(i,j),sai(i,j),coszen(i,j)&
              ,albg(1:,1:,i,j),albv(1:,1:,i,j),alb(1:,1:,i,j),ssun(1:,1:,i,j),ssha(1:,1:,i,j)&
              ,thermk(i,j),extkb(i,j),extkd(i,j)&
              ,trad(i,j),tref(i,j),qref(i,j),rst(i,j),emis(i,j),z0ma(i,j),zol(i,j),rib(i,j)&
              ,ustar(i,j),qstar(i,j),tstar(i,j),fm(i,j),fh(i,j),fq(i,j))
        enddo
      enddo


    ! ------------------------------------------
    ! PLEASE
    ! PLEASE UPDATE
    ! PLEASE UPDATE when have the observed lake status
      t_lake      (:,:,:) = 285.
      lake_icefrac(:,:,:) = 0.
      dz_lake     (:,:,:) = 1.
      smf_soisno  (:,:,:) = 0.
      smf_ice     (:,:,:) = 0.
    ! ------------------------------------------

END SUBROUTINE initialize

SUBROUTINE IniTimeVar(nl_soil,maxsnl,itypwat&
                     ,porsl,soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb&
                     ,nl_soil_ini,soil_z,soil_t,soil_w,snow_d&
                     ,z0m,zlnd,chil,ref,tran,z_soisno,dz_soisno&
                     ,t_soisno,wliq_soisno,wice_soisno,zwt,wa&
                     ,t_grnd,tlsun,tlsha,ldew,sag,scv&
                     ,snowdp,fveg,fsno,sigf,green,lai,sai,coszen&
                     ,albg,albv,alb,ssun,ssha,thermk,extkb,extkd&
                     ,trad,tref,qref,rst,emis,z0ma,zol,rib&
                     ,ustar,qstar,tstar,fm,fh,fq)

!=======================================================================
! Created by Yongjiu Dai, 09/15/1999
! Revised by Yongjiu Dai, 08/30/2002
!                         03/2014
! Revised by Hongyi Li,   09/2015
!=======================================================================

  use precision
  use PhysicalConstants, only : tfrz
  use ALBEDO

  IMPLICIT NONE

  integer, INTENT(in) ::        &!
        nl_soil,                &! soil layer number
        maxsnl,                 &! maximum snow layer number
        itypwat                  ! index for land cover type [-]

  real(r8), INTENT(in) ::       &!
        fveg,                   &! fraction of vegetation cover
        green,                  &! leaf greenness
        lai,                    &! leaf area index
        sai,                    &! stem area index
        coszen,                 &! cosine of solar zenith angle
        soil_s_v_alb,           &! albedo of visible of the saturated soil
        soil_d_v_alb,           &! albedo of visible of the dry soil
        soil_s_n_alb,           &! albedo of near infrared of the saturated soil
        soil_d_n_alb,           &! albedo of near infrared of the dry soil
        z0m,                    &! aerodynamic roughness length [m]
        zlnd,                   &! aerodynamic roughness length over soil surface [m]
        chil,                   &! leaf angle distribution factor
        ref (2,2),              &! leaf reflectance (iw=iband, il=life and dead)
        tran(2,2),              &! leaf transmittance (iw=iband, il=life and dead)
        porsl(1:nl_soil)         ! porosity of soil


  integer, INTENT(in) :: nl_soil_ini
  real(r8), INTENT(in) ::       &!
        soil_z(nl_soil_ini),    &! soil layer depth for initial (m)
        soil_t(nl_soil_ini),    &! soil temperature from initial file (K)
        soil_w(nl_soil_ini),    &! soil wetness from initial file (-)
        snow_d                   ! snow depth (m)

  real(r8), INTENT(inout) ::    &!
        z_soisno (maxsnl+1:nl_soil),   &! node depth [m]
        dz_soisno(maxsnl+1:nl_soil)     ! layer thickness [m]

  real(r8), INTENT(out) ::      &!
        t_soisno (maxsnl+1:nl_soil), &! soil temperature [K]
        wliq_soisno(maxsnl+1:nl_soil), &! liquid water in layers [kg/m2]
        wice_soisno(maxsnl+1:nl_soil), &! ice lens in layers [kg/m2]
        t_grnd,                 &! ground surface temperature [K]
        tlsun,                  &! sunlit leaf temperature [K]
        tlsha,                  &! shaded leaf temperature [K]
        ldew,                   &! depth of water on foliage [mm]
        sag,                    &! non dimensional snow age [-]
        scv,                    &! snow cover, water equivalent [mm]
        snowdp,                 &! snow depth [meter]
        fsno,                   &! fraction of snow cover on ground
        sigf,                   &! fraction of veg cover, excluding snow-covered veg [-]
        albg(2,2),              &! albedo, ground [-]
        albv(2,2),              &! albedo, vegetation [-]
        alb (2,2),              &! averaged albedo [-]
        ssun(2,2),              &! sunlit canopy absorption for solar radiation
        ssha(2,2),              &! shaded canopy absorption for solar radiation
        thermk,                 &! canopy gap fraction for tir radiation
        extkb,                  &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,                  &! diffuse and scattered diffuse PAR extinction coefficient
        wa,                     &! water storage in aquifer [mm]
        zwt,                    &! the depth to water table [m]
        trad,                   &! radiative temperature of surface [K]
        tref,                   &! 2 m height air temperature [kelvin]
        qref,                   &! 2 m height air specific humidity
        rst,                    &! canopy stomatal resistance (s/m)
        emis,                   &! averaged bulk surface emissivity
        z0ma,                   &! effective roughness [m]
        zol,                    &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,                    &! bulk Richardson number in surface layer
        ustar,                  &! u* in similarity theory [m/s]
        qstar,                  &! q* in similarity theory [kg/kg]
        tstar,                  &! t* in similarity theory [K]
        fm,                     &! integral of profile function for momentum
        fh,                     &! integral of profile function for heat
        fq                       ! integral of profile function for moisture

        integer j, snl
        real(r8) wet(nl_soil), wt, ssw, oro, rhosno_ini, a
!-----------------------------------------------------------------------

  if(itypwat <= 5)then ! land grid
     rhosno_ini = 250.

     do j = 1, nl_soil
        call polint(soil_z,soil_t,nl_soil_ini,z_soisno(j),t_soisno(j))
        call polint(soil_z,soil_w,nl_soil_ini,z_soisno(j),wet(j))
        a           = min(soil_t(1),soil_t(2),soil_t(3))-5.
        t_soisno(j) = max(t_soisno(j), a)
        a           = max(soil_t(1),soil_t(2),soil_t(3))+5.
        t_soisno(j) = min(t_soisno(j), a)
        a      = min(soil_w(1),soil_w(2),soil_w(3))
        wet(j) = max(wet(j), a, 0.1)
        a      = max(soil_w(1),soil_w(2),soil_w(3))
        wet(j) = min(wet(j), a, 0.5)
        wet(j) = max(wet(j), a, 0.1)

        if (t_soisno(j) .gt. tfrz) then
           wliq_soisno(j) = wet(j)*dz_soisno(j)*1000.
!          wliq_soisno(j) = porsl(j)*wet(j)*dz_soisno(j)*1000.
           wice_soisno(j) = 0.
        else
           wliq_soisno(j) = 0.07*dz_soisno(j)*1000.
           ! wliq_soisno(j) = 0.5*wet(j)*dz_soisno(j)*1000.
           wice_soisno(j) = wet(j)*dz_soisno(j)*1000.
!          wliq_soisno(j) = porsl(j)*wet(j)*dz_soisno(j)*1000.
        endif
! print*,'debug3',wliq_soisno
! print*,'debug4', wice_soisno
     enddo

     snowdp = snow_d
     sag    = 0.
     scv    = snowdp*rhosno_ini

     call snowfraction (fveg,z0m,zlnd,scv,snowdp,wt,sigf,fsno)
     call snow_ini (itypwat,maxsnl,snowdp,snl,z_soisno,dz_soisno)
     if(snl.lt.0)then
        do j = snl+1, 0
           t_soisno(j) = min(tfrz-1., t_soisno(1))
           wliq_soisno(j) = 0.
           wice_soisno(j) = dz_soisno(j)*rhosno_ini         ! m * kg m-3 = kg m-2
        enddo
     endif
! ---------- added 20151219
! water table depth (initially at 1.0 m below the model bottom; wa when zwt
!                    is below the model bottom zi(nl_soil)
     wa  = 4800.                             ! assuming aquifer capacity is 5000 mm
     zwt = (25. + z_soisno(nl_soil))+dz_soisno(nl_soil)/2. - wa/1000./0.2 ! to result in zwt = zi(nl_soil) + 1.0 m
! ---------------------------------------------


     if(snl>maxsnl)then
        t_soisno (maxsnl+1:snl) = -999.
        wice_soisno(maxsnl+1:snl) = 0.
        wliq_soisno(maxsnl+1:snl) = 0.
        z_soisno   (maxsnl+1:snl) = 0.
        dz_soisno  (maxsnl+1:snl) = 0.
     endif
     ldew  = 0.
     tlsun = t_soisno(1)
     tlsha = t_soisno(1)
     t_grnd = t_soisno(1)

! surface albedo
     ssw = min(1.,1.e-3*wliq_soisno(1)/dz_soisno(1))
     call albland (itypwat,soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb,&
                   chil,ref,tran,fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,t_grnd,&
                   alb,albg,albv,ssun,ssha,thermk,extkb,extkd)

  else                 ! ocean grid
     t_soisno(:) = 300.
     wice_soisno(:) = 0.
     wliq_soisno(:) = 1000.
     z_soisno (maxsnl+1:0) = 0.
     dz_soisno(maxsnl+1:0) = 0.
     sigf   = 0.
     fsno   = 0.
     ldew   = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.
     tlsun  = 300.
     tlsha  = 300.
     t_grnd = 300.

     oro = 0
     call albocean (oro,scv,coszen,alb)
     albg(:,:) = alb(:,:)
     albv(:,:) = 0.0
     ssun(:,:) = 0.0
     ssha(:,:) = 0.0
     thermk = 0.0
     extkb = 0.0
     extkd = 0.0
  endif

! Additional variables required by reginal model (WRF & RSM)
! totally arbitrarily assigned here
  trad  = t_grnd
  tref  = t_grnd
  qref  = 0.3
  rst   = 1.e36
  emis  = 1.0
  z0ma  = 0.01
  zol   = -1.0
  rib   = -0.1
  ustar = 0.25
  qstar = 0.001
  tstar = -1.5
  fm    = alog(30.)
  fh    = alog(30.)
  fq    = alog(30.)

END SUBROUTINE IniTimeVar


SUBROUTINE snow_ini(itypwat,maxsnl,snowdp,snl,z_soisno,dz_soisno)

! Snow spatial discretization initially

  use precision
  implicit none

  integer, intent(in) :: maxsnl  ! maximum of snow layers
  integer, intent(in) :: itypwat ! index for land cover type [-]
  real(r8), intent(in) :: snowdp ! snow depth [m]
  real(r8), intent(out) :: z_soisno (maxsnl+1:0) ! node depth [m]
  real(r8), intent(out) :: dz_soisno(maxsnl+1:0) ! layer thickness [m]
  integer, intent(out) :: snl ! number of snow layer
  real(r8) zi
  integer i
!-----------------------------------------------------------------------

  dz_soisno(:0) = 0.
  z_soisno(:0) = 0.
  snl = 0
  if(itypwat.le.3)then ! non water bodies

     if(snowdp.lt.0.01)then
        snl = 0
     else
        if(snowdp>=0.01 .and. snowdp<=0.03)then
           snl = -1
           dz_soisno(0)  = snowdp
        else if(snowdp>0.03 .and. snowdp<=0.04)then
           snl = -2
           dz_soisno(-1) = snowdp/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.04 .and. snowdp<=0.07)then
           snl = -2
           dz_soisno(-1) = 0.02
           dz_soisno( 0) = snowdp - dz_soisno(-1)
        else if(snowdp>0.07 .and. snowdp<=0.12)then
           snl = -3
           dz_soisno(-2) = 0.02
           dz_soisno(-1) = (snowdp - 0.02)/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.12 .and. snowdp<=0.18)then
           snl = -3
           dz_soisno(-2) = 0.02
           dz_soisno(-1) = 0.05
           dz_soisno( 0) = snowdp - dz_soisno(-2) - dz_soisno(-1)
        else if(snowdp>0.18 .and. snowdp<=0.29)then
           snl = -4
           dz_soisno(-3) = 0.02
           dz_soisno(-2) = 0.05
           dz_soisno(-1) = (snowdp - dz_soisno(-3) - dz_soisno(-2))/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.29 .and. snowdp<=0.41)then
           snl = -4
           dz_soisno(-3) = 0.02
           dz_soisno(-2) = 0.05
           dz_soisno(-1) = 0.11
           dz_soisno( 0) = snowdp - dz_soisno(-3) - dz_soisno(-2) - dz_soisno(-1)
        else if(snowdp>0.41 .and. snowdp<=0.64)then
           snl = -5
           dz_soisno(-4) = 0.02
           dz_soisno(-3) = 0.05
           dz_soisno(-2) = 0.11
           dz_soisno(-1) = (snowdp - dz_soisno(-4) - dz_soisno(-3) - dz_soisno(-2))/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.64)then
           snl = -5
           dz_soisno(-4) = 0.02
           dz_soisno(-3) = 0.05
           dz_soisno(-2) = 0.11
           dz_soisno(-1) = 0.23
           dz_soisno( 0) = snowdp - dz_soisno(-4) - dz_soisno(-3) - dz_soisno(-2) - dz_soisno(-1)
        endif

        zi = 0.
        do i = 0, snl+1, -1
           z_soisno(i) = zi - dz_soisno(i)/2.
           ! In origin CLM: zi = -zi-dz_soisno(i)
           zi = zi-dz_soisno(i) ! Corrected by HY, change '-zi' to 'zi',2016-4-6
        enddo
     endif

  endif

END SUBROUTINE snow_ini
!-----------------------------------------------------------------------
! EOP


SUBROUTINE polint(xa,ya,n,x,y)

! Given arrays xa and ya, each of length n, and gi
! value y, and an error estimate dy. If P (x) is the p
! P (xa(i)) = ya(i), i = 1, . . . , n, then the returned value
! (from: "Numerical Recipes")

  use precision
  implicit none
  integer n,NMAX
  real(r8) dy,x,y,xa(n),ya(n)
  parameter (NMAX=10)      !Largest anticipated val
  integer i,m,ns
  real(r8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(x-xa(1))

  do i=1,n       !Here we find the index ns of the closest table entry,
     dift=abs(x-xa(i))
     if(dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)  !and initialize the tableau of c's and d's.
     d(i)=ya(i)
  enddo

  y=ya(ns)       !This is the initial approximation to y.
  ns=ns-1

  do m=1,n-1  !For each column of the tableau,
     do i=1,n-m   !we loop over the current c's and d's and update them.
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) print*, 'failure in polint'  !two input xa's are identical.
        den=w/den
        d(i)=hp*den                                !here the c's and d's are updated.
        c(i)=ho*den
     enddo
     if(2*ns.lt.n-m)then  !After each column in the tableau is completed, we decide
        dy=c(ns+1)        !which correction, c or d, we want to add to our accumulating
     else                 !value of y, i.e., which path to take through
        dy=d(ns)          !the tableau-forking up or down. We do this in such a
        ns=ns-1           !way as to take the most "straight line" route through the
     endif                !tableau to its apex, updating ns accordingly to keep track
     y=y+dy               !of where we are. This route keeps the partial approximations
  enddo                   !centered (insofar as possible) on the target x. T he
                          !last dy added is thus the error indication.
END SUBROUTINE polint

end module initial
!-----------------------------------------------------------------------
! EOP

