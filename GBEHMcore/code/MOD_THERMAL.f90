MODULE MOD_THERMAL

!-----------------------------------------------------------------------
 use precision
 IMPLICIT NONE
 SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: THERMAL
  public :: groundfluxes
  public :: groundtem


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------
 SUBROUTINE THERMAL (itypwat     ,lb    ,nl   ,nl_soil     ,deltim     ,&
                     trsmx0      ,zlnd        ,zsno        ,csoilc     ,&
                     dewmx       ,capr        ,cnfac       ,csol       ,&
                     wrsd        ,adjfac      ,smf_soisno  ,smf_ice    ,&
                     porsl       ,psi0        ,bsw         ,dkdry      ,&
                     dksatu      ,lai         ,laisun      ,laisha     ,&
                     sai         ,z0m         ,displa      ,sqrtdi     ,&
                     rootfr      ,rstfac      ,effcon      ,vmax25     ,&
                     slti        ,hlti        ,shti        ,hhti       ,&
                     trda        ,trdm        ,trop        ,gradm      ,&
                     binter      ,extkn       ,forc_hgt_u  ,forc_hgt_t ,&
                     forc_hgt_q  ,forc_us     ,forc_vs     ,forc_t     ,&
                     forc_q      ,forc_rhoair ,forc_psrf   ,forc_pco2m ,&
                     forc_po2m   ,coszen      ,parsun      ,parsha     ,&
                     sabvsun     ,sabvsha     ,sabg        ,frl        ,&
                     extkb       ,extkd       ,thermk      ,fsno       ,&
                     sigf        ,dz_soisno   ,z_soisno    ,zi_soisno  ,&
                     tlsun       ,tlsha       ,t_soisno    ,wice_soisno,&
                     wliq_soisno ,ldew        ,scv         ,snowdp     ,&
                     imelt       ,taux        ,tauy        ,fsena      ,&
                     fevpa       ,lfevpa      ,fsenl       ,fevpl      ,&
                     etr         ,fseng       ,fevpg       ,olrg       ,&
                     fgrnd       ,rootr       ,qseva       ,qsdew      ,&
                     qsubl       ,qfros       ,sm          ,tref       ,&
                     qref        ,trad        ,rst         ,assim      ,&
                     respc       ,errore      ,emis        ,z0ma       ,&
                     zol         ,rib         ,ustar       ,qstar      ,&
                     tstar       ,fm          ,fh          ,fq          )

!=======================================================================
! this is the main subroutine to execute the calculation 
! of thermal processes and surface fluxes
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
! 
! FLOW DIAGRAM FOR THERMAL.F90
! 
! THERMAL ===> qsadv
!              groundfluxes
!              eroot                      |dewfraction
!              leaftemone |               |qsadv
!              leaftemtwo |  ---------->  |moninobukini
!                                         |moninobuk
!                                         |ASSIM_STOMATA_conductance
!
!              groundTem     ---------->   meltf
!                               
!=======================================================================

  use precision
  use PhysicalConstants, only : denh2o,roverg,hvap,hsub,rgas,cpair,&
                                stefnc,denice,tfrz,vonkar,grav 
  use FRICTION_VELOCITY
  use LEAF_temperature

  ! IMPLICIT NONE
 
!---------------------Argument------------------------------------------

  integer, INTENT(in) :: &
        lb,          &! lower bound of array 
        nl,          &! total ground layers from soil surface to aquifer bottom
        nl_soil,     &! upper bound of array
        itypwat       ! land water type (0=soil, 1=urban or built-up, 2=wetland,
                      !                  3=glacier/ice sheet, 4=land water bodies)
  real(r8), INTENT(in) :: &
        deltim,      &! model time step [second]
        trsmx0,      &! max transpiration for moist soil+100% veg.  [mm/s]
        zlnd,        &! roughness length for soil [m]
        zsno,        &! roughness length for snow [m]
        csoilc,      &! drag coefficient for soil under canopy [-]
        dewmx,       &! maximum dew
        capr,        &! tuning factor to turn first layer T into surface T
        cnfac         ! Crank Nicholson factor between 0 and 1
        ! soil physical parameters
  real(r8), INTENT(in) :: &        
        csol(1:nl), &! heat capacity of soil solids [J/(m3 K)]
        porsl(1:nl),&! soil porosity [-]
        psi0(1:nl), &! soil water suction, negative potential [m]
        bsw(1:nl),  &! clapp and hornbereger "b" parameter [-]
        dkdry(1:nl),&! thermal conductivity of dry soil [W/m-K]
        dksatu(1:nl), &! thermal conductivity of saturated soil [W/m-K]
        wrsd,&
        adjfac  ! an adjustment parameter to solve the calibration of surface evaporation

        ! vegetation parameters
    real(r8), INTENT(in) :: &      
        lai,         &! adjusted leaf area index for seasonal variation [-]
        sai,         &! stem area index  [-]
        z0m,         &! roughness length, momentum [m]
        displa,      &! displacement height [m]
        sqrtdi,      &! inverse sqrt of leaf dimension [m**-0.5]
        rootfr(1:nl_soil),&! root fraction        
        effcon,      &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25,      &! maximum carboxylation rate at 25 C at canopy top
        slti,        &! slope of low temperature inhibition function      [s3]
        hlti,        &! 1/2 point of low temperature inhibition function  [s4]
        shti,        &! slope of high temperature inhibition function     [s1]
        hhti,        &! 1/2 point of high temperature inhibition function [s2]
        trda,        &! temperature coefficient in gs-a model             [s5]
        trdm,        &! temperature coefficient in gs-a model             [s6]
        trop,        &! temperature coefficient in gs-a model          
        gradm,       &! conductance-photosynthesis slope parameter
        binter,      &! conductance-photosynthesis intercept
        extkn         ! coefficient of leaf nitrogen allocation

        ! atmospherical variables and observational height
    real(r8), INTENT(in) :: &         
        forc_hgt_u,  &! observational height of wind [m]
        forc_hgt_t,  &! observational height of temperature [m]
        forc_hgt_q,  &! observational height of humidity [m]
        forc_us,     &! wind component in eastward direction [m/s]
        forc_vs,     &! wind component in northward direction [m/s]
        forc_t,      &! temperature at agcm reference height [kelvin]
        forc_q,      &! specific humidity at agcm reference height [kg/kg]
        forc_rhoair, &! density air [kg/m3]
        forc_psrf,   &! atmosphere pressure at the surface [pa]
        forc_pco2m,  &! CO2 concentration in atmos. (pascals)
        forc_po2m     ! O2 concentration in atmos. (pascals)

        ! radiative fluxes
    real(r8), INTENT(in) :: &         
        coszen,      &! cosine of the solar zenith angle
        parsun,      &! photosynthetic active radiation by sunlit leaves (W m-2)
        parsha,      &! photosynthetic active radiation by shaded leaves (W m-2)
        sabvsun,     &! solar radiation absorbed by vegetation [W/m2]
        sabvsha,     &! solar radiation absorbed by vegetation [W/m2]
        sabg,        &! solar radiation absorbed by ground [W/m2]
        frl,         &! atmospheric infrared (longwave) radiation [W/m2]
        extkb,       &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,       &! diffuse and scattered diffuse PAR extinction coefficient
        thermk        ! canopy gap fraction for tir radiation

        ! state variable (1)
    real(r8), INTENT(in) :: &         
        fsno,        &! fraction of ground covered by snow
        sigf,        &! fraction of veg cover, excluding snow-covered veg [-]
        dz_soisno(lb:nl),  &! layer thickiness [m]
        z_soisno (lb:nl),  &! node depth [m]
        zi_soisno(lb-1:nl)  ! interface depth [m]

        ! state variables (2)
  real(r8), INTENT(inout) :: &
        tlsun,       &! sunlit leaf temperature [K]
        tlsha,       &! shaded leaf temperature [K]
        t_soisno(lb:nl),   &! soil temperature [K]
        wice_soisno(lb:nl),&! ice lens [kg/m2]
        wliq_soisno(lb:nl),&! liqui water [kg/m2]
        smf_soisno (-1:nl_soil),&! snowmelt fraction in snow-soil (0-1).Added by HY,2016-11-26
        smf_ice    (1:nl_soil),&! snowmelt fraction in soil ice  (0-1).Added by HY,2016-11-30
        ldew,        &! depth of water on foliage [kg/(m2 s)] 
        scv,         &! snow cover, water equivalent [mm, kg/m2]
        snowdp        ! snow depth [m]

  integer, INTENT(out) :: & 
       imelt(lb:nl) ! flag for melting or freezing [-]
  
  real(r8), INTENT(out) :: &
       laisun,       &! sunlit leaf area index
       laisha,       &! shaded leaf area index
       rstfac         ! factor of soil water stress 
 
        ! Output fluxes
  real(r8), INTENT(out) :: &
        taux,        &! wind stress: E-W [kg/m/s**2]
        tauy,        &! wind stress: N-S [kg/m/s**2]
        fsena,       &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa,       &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa,      &! latent heat flux from canopy height to atmosphere [W/m2]
        fsenl,       &! ensible heat from leaves [W/m2]
        fevpl,       &! evaporation+transpiration from leaves [mm/s]
        etr,         &! transpiration rate [mm/s]
        fseng,       &! sensible heat flux from ground [W/m2]
        fevpg,       &! evaporation heat flux from ground [mm/s]
        olrg,        &! outgoing long-wave radiation from ground+canopy
        fgrnd,       &! ground heat flux [W/m2]
        rootr(1:nl_soil),&! root resistance of a layer, all layers add to 1
        qseva,       &! ground surface evaporation rate (mm h2o/s)
        qsdew,       &! ground surface dew formation (mm h2o /s) [+]
        qsubl,       &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros,       &! surface dew added to snow pack (mm h2o /s) [+]
        sm,          &! rate of snowmelt [kg/(m2 s)]
        tref,        &! 2 m height air temperature [kelvin]
        qref,        &! 2 m height air specific humidity
        trad,        &! radiative temperature [K]
        rst,         &! stomatal resistance (s m-1)
        assim,       &! assimilation
        respc         ! respiration

       ! additional variables required by coupling with WRF or RSM model
    real(r8), INTENT(out) :: &        
        emis,        &! averaged bulk surface emissivity
        z0ma,        &! effective roughness [m]
        zol,         &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,         &! bulk Richardson number in surface layer
        ustar,       &! u* in similarity theory [m/s]
        qstar,       &! q* in similarity theory [kg/kg]
        tstar,       &! t* in similarity theory [K]
        fm,          &! integral of profile function for momentum
        fh,          &! integral of profile function for heat
        fq            ! integral of profile function for moisture

  ! integer, INTENT(out) :: & 
        ! ipatch        ! patch index

!---------------------Local Variables-----------------------------------

  integer i,j,ii

  real(r8) :: &
       cgrnd,        &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
       cgrndl,       &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
       cgrnds,       &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
       degdT,        &! d(eg)/dT
       dqgdT,        &! d(qg)/dT
       dlrad,        &! downward longwave radiation blow the canopy [W/m2]
       eg,           &! water vapor pressure at temperature T [pa]
       egsmax,       &! max. evaporation which soil can provide at one time step
       egidif,       &! the excess of evaporation over "egsmax"
       emg,          &! ground emissivity (0.97 for snow, glaciers and water surface; 0.96 for soil and wetland)
       errore,       &! energy balnce error [w/m2]
       etrc,         &! maximum possible transpiration rate [mm/s]
       fac,          &! soil wetness of surface layer
       fact(lb:nl), &! used in computing tridiagonal matrix
       fsun,         &! fraction of sunlit canopy
       hr,           &! relative humidity
       htvp,         &! latent heat of vapor of water (or sublimation) [j/kg]
       olru,         &! olrg excluding dwonwelling reflection [W/m2]
       olrb,         &! olrg assuming blackbody emission [W/m2]
       psit,         &! negative potential of soil
       par,          &! PAR absorbed by canopy [W/m2]
       qg,           &! ground specific humidity [kg/kg]
       qsatg,        &! saturated humidity [kg/kg]
       qsatgdT,      &! d(qsatg)/dT
       qred,         &! soil surface relative humidity
       sabv,         &! solar absorbed by canopy [W/m2]
       thm,          &! intermediate variable (forc_t+0.0098*forc_hgt_t)
       th,           &! potential temperature (kelvin)
       thv,          &! virtual potential temperature (kelvin)
       tl,           &! leaf temperature
       t_grnd,       &! ground surface temperature [K]
       t_soisno_bef(lb:nl), &! soil/snow temperature before update
       tinc,         &! temperature difference of two time step
       ur,           &! wind speed at reference height [m/s]
       ulrad,        &! upward longwave radiation above the canopy [W/m2]
       wice0(lb:nl),&! ice mass from previous time-step
       wliq0(lb:nl),&! liquid mass from previous time-step
       liqdiff(1:nl),&! 
       icediff(1:nl),&! 
       wx,           &! patitial volume of ice and water of surface layer
       xmf            ! total latent heat of phase change of ground water

  real(r8) :: z0ma_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g
  real(r8) :: fm10m,fm_g,fh_g,fq_g,fh2m,fq2m,um,obu
  real(r8), parameter :: watmin = 0.01  ! Added by HY, 2016-4-25
                                        ! Limit irreduciable wrapping liquid water 
                                        ! a tunable constant
                                        ! there soil moisture is limited.
                                        ! snow moisture is not considered to be limited.
  real(r8) :: tmp_liq_smf,tmp_ice_smf
  real(r8) :: debugtmp1,debugtmp2,sumsmf
!=======================================================================
! [1] Initial set and propositional variables
!=======================================================================

      ! fluxes 
      taux   = 0.;  tauy   = 0.    
      fsena  = 0.;  fevpa  = 0.  
      lfevpa = 0.;  fsenl  = 0.    
      fevpl  = 0.;  etr    = 0.  
      fseng  = 0.;  fevpg  = 0.    
      dlrad  = 0.;  ulrad  = 0. 
      cgrnds = 0.;  cgrndl = 0.    
      cgrnd  = 0.;  tref   = 0. 
      qref   = 0.;  rst    = 2.0e4
      assim  = 0.;  respc  = 0. 

      emis   = 0.;  z0ma   = 0.
      zol    = 0.;  rib    = 0.
      ustar  = 0.;  qstar  = 0.
      tstar  = 0.;  rootr  = 0.

      ! temperature and water mass from previous time step
      t_grnd = t_soisno(lb)
      t_soisno_bef(lb:) = t_soisno(lb:)
      wice0(lb:) = wice_soisno(lb:)
      wliq0(lb:) = wliq_soisno(lb:)

      ! emissivity
      emg = 0.96
      if(scv>0. .OR. itypwat==3) emg = 0.97

      ! latent heat, assumed that the sublimation occured only as wliq_soisno=0
      htvp = hvap
      if(wliq_soisno(lb)<=0. .and. wice_soisno(lb)>0.) htvp = hsub

      ! potential temperatur at the reference height
      thm = forc_t + 0.0098*forc_hgt_t              ! intermediate variable equivalent to
                                                    ! forc_t*(pgcm/forc_psrf)**(rgas/cpair)
      th = forc_t*(100000./forc_psrf)**(rgas/cpair) ! potential T
      thv = th*(1.+0.61*forc_q)             ! virtual potential T
      ur = max(0.1,sqrt(forc_us*forc_us+forc_vs*forc_vs))   ! limit set to 0.1

!=======================================================================
! [2] specific humidity and its derivative at ground surface
!=======================================================================

      qred = 1.
      call qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)

      if(itypwat<=1)then            ! soil ground
         wx   = (wliq_soisno(1)/denh2o + wice_soisno(1)/denice)/dz_soisno(1)
         if(porsl(1)<1.e-6)then     ! bed rock
            fac  = 0.001
         else 
            fac  = min(1.,adjfac*wx/porsl(1))
            fac  = max( fac, 0.001 )
         endif

         psit = psi0(1) * fac ** (- bsw(1) )   ! psit = max(smpmin, psit)
         ! psit = 0.001*psi0(1) * fac ** (- bsw(1) )   ! Debug by HY, 2016-4-5; *0.001, because psi0 is in unit mm.
         psit = max( -1.e8, psit )
         hr   = exp(psit/roverg/t_grnd)
         qred = (1.-fsno)*hr + fsno
      endif

      qg = qred*qsatg  
      dqgdT = qred*qsatgdT

      if(qsatg > forc_q .and. forc_q > qred*qsatg)then
        qg = forc_q; dqgdT = 0.
      endif

!=======================================================================
! [3] Compute sensible and latent fluxes and their derivatives with respect 
!     to ground temperature using ground temperatures from previous time step.
!=======================================================================

      if(sigf <= 0.999) then
         call groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q, &
                            forc_us,forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, &
                            ur,thm,th,thv,t_grnd,qg,dqgdT,htvp, &
                            fsno,sigf,cgrnd,cgrndl,cgrnds, &
                            taux,tauy,fsena,fevpa,fseng,fevpg,tref,qref, &
              z0ma_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g,fm_g,fh_g,fq_g)
      end if
 
!=======================================================================
! [4] Canopy temperature, fluxes from the canopy
!=======================================================================

      par = parsun + parsha
      sabv = sabvsun + sabvsha

      if(sigf >= 0.001) then

         ! soil water strees factor on stomatal resistance
         call eroot (nl_soil,trsmx0,porsl,bsw,psi0,rootfr,&
                     dz_soisno,t_soisno,wliq_soisno,rootr,etrc,rstfac)

         ! fraction of sunlit and shaded leaves of canopy
         fsun = ( 1. - exp(-min(extkb*lai,40.))) / max( min(extkb*lai,40.), 1.e-6 )

         if(coszen<=0.0 .OR. sabv<1.) fsun = 0.
         
         laisun = lai*fsun
         laisha = lai*(1-fsun)

         if(fsun.le.0.1)then

            fsun = 0.
            tl = tlsha

            call leaftemone (deltim     ,csoilc  ,dewmx     ,htvp        ,&
                 lai        ,sai        ,displa  ,sqrtdi    ,z0m         ,&
                 effcon     ,vmax25     ,slti    ,hlti      ,shti        ,&
                 hhti       ,trda       ,trdm    ,trop      ,gradm       ,&
                 binter     ,extkn      ,extkb   ,extkd     ,forc_hgt_u  ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us ,forc_vs   ,thm         ,&
                 th         ,thv        ,forc_q  ,forc_psrf ,forc_rhoair ,&
                 par        ,sabv       ,frl     ,thermk    ,rstfac      ,&
                 forc_po2m  ,forc_pco2m ,sigf    ,etrc      ,t_grnd      ,&
                 qg         ,dqgdT      ,emg     ,tl        ,ldew        ,&
                 taux       ,tauy       ,fseng   ,fevpg     ,cgrnd       ,&
                 cgrndl     ,cgrnds     ,tref    ,qref      ,rst         ,&
                 assim      ,respc      ,fsenl   ,fevpl     ,etr         ,&
                 dlrad      ,ulrad      ,z0ma    ,zol       ,rib         ,&
                 ustar      ,qstar      ,tstar   ,fm        ,fh          ,&
                 fq)

                 tlsun = tl
                 tlsha = tl

         else

            call leaftemtwo (deltim     ,csoilc  ,dewmx     ,htvp        ,&
                 lai        ,sai        ,displa  ,sqrtdi    ,z0m         ,&
                 effcon     ,vmax25     ,slti    ,hlti      ,shti        ,&
                 hhti       ,trda       ,trdm    ,trop      ,gradm       ,&
                 binter     ,extkn      ,extkb   ,extkd     ,forc_hgt_u  ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us ,forc_vs   ,thm         ,&
                 th         ,thv        ,forc_q  ,forc_psrf ,forc_rhoair ,&
                 parsun     ,parsha     ,sabvsun ,sabvsha   ,frl         ,&
                 fsun       ,thermk     ,rstfac  ,forc_po2m ,forc_pco2m  ,&
                 sigf       ,etrc       ,t_grnd  ,qg        ,dqgdT       ,&
                 emg        ,tlsun      ,tlsha   ,ldew      ,taux        ,&
                 tauy       ,fseng      ,fevpg   ,cgrnd     ,cgrndl      ,&
                 cgrnds     ,tref       ,qref    ,rst       ,assim       ,&
                 respc      ,fsenl      ,fevpl   ,etr       ,dlrad       ,&
                 ulrad      ,z0ma       ,zol     ,rib       ,ustar       ,&
                 qstar      ,tstar      ,fm      ,fh        ,fq)
         endif

      endif

    ! equate canopy temperature to air over bareland.
    ! required as sigf=0 carried over to next time step
      if(sigf < 0.001)then
         tlsun = forc_t
         tlsha = forc_t
         laisun  = 0.
         laisha  = 0.
         ldew = 0.
      endif
!=======================================================================
! [5] Gound temperature
!=======================================================================

      call groundtem (itypwat,lb,nl,deltim,&
                      capr,cnfac,csol,porsl,dkdry,dksatu,&
                      sigf,dz_soisno,z_soisno,zi_soisno,&
                      t_soisno,wice_soisno,wliq_soisno,scv,snowdp,&
                      frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg,&
                      imelt,sm,xmf,fact,psi0,bsw,wrsd)
if(lb<=0)then
  if(sum(wliq_soisno(lb:0))>0.)then
    smf_soisno(-1) = (smf_soisno(-1)*sum(wliq0(lb:0))+(sum(wliq_soisno(lb:0))-sum(wliq0(lb:0))))&
                    /sum(wliq_soisno(lb:0))
    smf_soisno(-1) = min(1.,smf_soisno(-1))
  else
    smf_soisno(-1) = 0.
  endif    
endif             
      ! To calculate the smf_soisno due to soil ice melt. By HY.,2016-11-28
      !--Begin update
      liqdiff = wliq_soisno(1:) - wliq0(1:)
      icediff = wice_soisno(1:) - wice0(1:)
! if (abs(sum(icediff(1:)+liqdiff(1:)))>1E-10) then
!   print*,'wrong diff'      
!   ! pause
! endif
! debugtmp1 = sum(smf_soisno(1:)*wliq0(1:)+smf_ice(1:)*wice0(1:))
      do j =1,nl_soil
        ! if (abs(icediff(j)+liqdiff(j))>1E-10) pause
        sumsmf = smf_soisno(j)*wliq_soisno(j)+smf_ice(j)*wice_soisno(j)
        if (liqdiff(j)>0. .and. wliq_soisno(j)>1E-8)then
          ! soil ice melting
          tmp_liq_smf   = smf_soisno(j)
          smf_soisno(j) = min(1.,(smf_soisno(j)*wliq0(j)+smf_ice(j)*liqdiff(j))/wliq_soisno(j))
          smf_soisno(j) = max(min(tmp_liq_smf,smf_ice(j)),smf_soisno(j))
        elseif (liqdiff(j)<0. .and. wice_soisno(j)>1E-8) then
          tmp_ice_smf   = smf_ice(j)
          smf_ice(j)    = min(1.,(smf_ice(j)*wice0(j)-smf_soisno(j)*liqdiff(j))/wice_soisno(j))
          smf_ice(j)    = max(min(tmp_ice_smf,smf_soisno(j)),smf_ice(j))
        endif
      enddo
! debugtmp2 = sum(smf_soisno(1:)*wliq_soisno(1:)+smf_ice(1:)*wice_soisno(1:))
! if(abs(debugtmp2-debugtmp1)>1E-8) print*,'MOD_THERMAL debug',debugtmp2,debugtmp1,debugtmp2-debugtmp1
      !--end update
  
!=======================================================================
! [6] Correct fluxes to present soil temperature
!=======================================================================

      t_grnd = t_soisno(lb)
      tinc = t_soisno(lb) - t_soisno_bef(lb)
! print*,'test1',fevpg *3600     
      fseng = fseng + tinc*cgrnds 
      fevpg = fevpg + tinc*cgrndl      
! print*,'test2',fevpg *3600
! calculation of evaporative potential; flux in kg m-2 s-1.  
! egidif holds the excess energy if all water is evaporated
! during the timestep.  this energy is later added to the sensible heat flux.

      egsmax = (wice_soisno(lb)+wliq_soisno(lb)) / deltim
      ! Keep the moisture of top soil layer, Added by HY, 2016-4-25
      ! if (lb>0) then
      !   egsmax = egsmax-watmin*porsl(1)*dz_soisno(1)*1000. / deltim
      !   egsmax = max(egsmax,0.)
      ! endif
      
      egidif = max( 0., fevpg - egsmax )
      fevpg  = min ( fevpg, egsmax )
      fseng  = fseng + htvp*egidif

! total fluxes to atmosphere
      fsena = fsenl + fseng
      fevpa = fevpl + fevpg
      lfevpa= hvap*fevpl + htvp*fevpg   ! W/m^2 (accouting for sublimation)
      
      qseva = 0.
      qsubl = 0.
      qfros = 0.
      qsdew = 0.

      if(fevpg >= 0.)then
! not allow for sublimation in melting (melting ==> evap. ==> sublimation)
         qseva = min(wliq_soisno(lb)/deltim, fevpg)
         qsubl = fevpg - qseva
      else
         if(t_grnd < tfrz)then
            qfros = abs(fevpg)
         else
            qsdew = abs(fevpg)
         endif
      endif

! ground heat flux
      fgrnd = sabg + dlrad + (1.-sigf)*emg*frl &
            - emg*stefnc*t_soisno_bef(lb)**3*(t_soisno_bef(lb) + 4.*tinc) &
            - (fseng+fevpg*htvp)

! outgoing long-wave radiation from canopy + ground
      olrg = ulrad &
           + (1.-sigf)*(1.-emg)*frl &
           + (1.-sigf)*emg*stefnc * t_soisno_bef(lb)**4 &
! for conservation we put the increase of ground longwave to outgoing
           + 4.*emg*stefnc*t_soisno_bef(lb)**3*tinc

! averaged bulk surface emissivity 
      olrb = stefnc*t_soisno_bef(lb)**3*((1.-sigf)*t_soisno_bef(lb) + 4.*tinc)
      olru = ulrad + emg*olrb
      olrb = ulrad + olrb
      emis = olru / olrb

! radiative temperature
      trad = (olrg/stefnc)**0.25

! additonal variables required by WRF and RSM model
      if(sigf < 0.001)then
         ustar = ustar_g
         tstar = tstar_g
         qstar = qstar_g
         rib   = rib_g
         zol   = zol_g
         z0ma  = z0ma_g
         fm    = fm_g
         fh    = fh_g
         fq    = fq_g

      else if(sigf <= 0.7)then
         z0ma  = sigf*z0ma  + (1.-sigf)*z0ma_g

       ! assumed um ~= ur here
         um = ur
         ustar =   sqrt(max(1.e-6,sqrt(taux*taux+tauy*tauy))/forc_rhoair)
         tstar = - fsena/(cpair*ustar*forc_rhoair)
         qstar = - fevpa/(ustar*forc_rhoair)
         zol = (forc_hgt_u-displa)*vonkar*grav*(tstar+0.61*th*qstar)/(ustar**2*thv)
         if(zol .ge. 0.)then   !stable
            zol = min(2.,max(zol,1.e-6))
         else                  !unstable
            zol = max(-100.,min(zol,-1.e-6))
         endif

         obu = (forc_hgt_u-displa)/zol
         call moninobuk(forc_hgt_u,forc_hgt_t,forc_hgt_q,displa,z0ma,z0ma,z0ma,obu,um,&
              ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

      else
         ustar = ustar
         tstar = tstar
         qstar = qstar
         rib   = rib
         zol   = zol
         z0ma  = z0ma
         fm    = fm
         fh    = fh
         fq    = fq

      endif

!=======================================================================
! [7] energy balance error
!=======================================================================

      errore = sabv + sabg + frl - olrg - fsena - lfevpa - xmf
      do j = lb, nl
         errore = errore - (t_soisno(j)-t_soisno_bef(j))/fact(j)
      enddo
      
     if(abs(errore)>.5)then 
     
      ! open(6,file = 'log.txt', access='append', status='old')
             ! write(9, '(3i6,f16.3)') year,im,id, Qd(isub,j)
      write(6,*) 'THERMAL.F90 : energy  balance violation'
      write(6,*) errore,sabv,sabg,frl,olrg,fsenl,fseng,hvap*fevpl,htvp*fevpg,xmf,t_soisno(0), &
                 t_soisno(1),t_soisno(5),lb,wice_soisno(lb),wliq_soisno(lb),ulrad,t_soisno_bef(lb),snowdp
      ! close(6)
      do ii = lb, nl
        if (t_soisno(ii)< 180. .or. t_soisno(ii)>360.) then
          print*,'bad surface temperature',t_soisno(ii),'layer',ii
          print*,'tref=', tref
          print*,'forcing=', forc_us, forc_vs, forc_t, forc_q, forc_rhoair,forc_psrf,&
                  forc_pco2m,forc_po2m,sabvsun,sabvsha,sabg, frl   
          if (forc_t<290. .and. forc_t>255.) then
            t_soisno(ii) = forc_t
          else
            t_soisno(ii) = max(255.,t_soisno_bef(ii))
            t_soisno(ii) = min(290.,t_soisno(ii))
          endif  
        endif
      end do     
     endif

100  format(10(f15.3))

 END SUBROUTINE THERMAL

  subroutine groundfluxes (zlnd, zsno, hu, ht, hq, &
                          us, vs, tm, qm, rhoair, psrf, &
                          ur, thm, th, thv, t_grnd, qg, dqgdT, htvp, &
                          fsno, sigf, cgrnd, cgrndl, cgrnds, & 
                          taux, tauy, fsena, fevpa, fseng, fevpg, tref, qref, &
                          z0ma, zol, rib, ustar, qstar, tstar, fm, fh, fq)

!=======================================================================
! this is the main subroutine to execute the calculation of thermal processes
! and surface fluxes
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!=======================================================================

  use precision
  use PhysicalConstants, only : cpair,vonkar,grav
  use FRICTION_VELOCITY
  ! implicit none
 
!----------------------- Dummy argument --------------------------------
  real(r8), INTENT(in) :: &
        zlnd,     &! roughness length for soil [m]
        zsno       ! roughness length for snow [m]

        ! atmospherical variables and observational height
  real(r8), INTENT(in) :: &        
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq,       &! observational height of humidity [m]
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        tm,       &! temperature at agcm reference height [kelvin] [not used]
        qm,       &! specific humidity at agcm reference height [kg/kg]
        rhoair,   &! density air [kg/m3]
        psrf,     &! atmosphere pressure at the surface [pa] [not used]
        fsno,     &! fraction of ground covered by snow
        sigf,     &! fraction of veg cover, excluding snow-covered veg [-]
        ur,       &! wind speed at reference height [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht)
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)
        t_grnd,   &! ground surface temperature [K]
        qg,       &! ground specific humidity [kg/kg]
        dqgdT,    &! d(qg)/dT
        htvp       ! latent heat of vapor of water (or sublimation) [j/kg]

  real(r8), INTENT(out) :: &
        taux,     &! wind stress: E-W [kg/m/s**2]
        tauy,     &! wind stress: N-S [kg/m/s**2]
        fsena,    &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa,    &! evapotranspiration from canopy height to atmosphere [mm/s]
        fseng,    &! sensible heat flux from ground [W/m2]
        fevpg,    &! evaporation heat flux from ground [mm/s]
        cgrnd,    &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,   &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,   &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,     &! 2 m height air temperature [kelvin]
        qref,     &! 2 m height air humidity
        z0ma,     &! effective roughness [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq         ! integral of profile function for moisture

!------------------------ LOCAL VARIABLES ------------------------------
  integer niters, &! maximum number of iterations for surface temperature
       iter,      &! iteration index
       nmozsgn     ! number of times moz changes sign

  real(r8) :: &
       beta,      &! coefficient of conective velocity [-]
       displax,   &! zero-displacement height [m]
       dth,       &! diff of virtual temp. between ref. height and surface
       dqh,       &! diff of humidity between ref. height and surface
       dthv,      &! diff of vir. poten. temp. between ref. height and surface
       obu,       &! monin-obukhov length (m)
       obuold,    &! monin-obukhov length from previous iteration
       ram,       &! aerodynamical resistance [s/m]
       rah,       &! thermal resistance [s/m]
       raw,       &! moisture resistance [s/m]
       raih,      &! temporary variable [kg/m2/s]
       raiw,      &! temporary variable [kg/m2/s]
       fh2m,      &! relation for temperature at 2m
       fq2m,      &! relation for specific humidity at 2m
       fm10m,     &! integral of profile function for momentum at 10m
       thvstar,   &! virtual potential temperature scaling parameter
       um,        &! wind speed including the stablity effect [m/s]
       wc,        &! convective velocity [m/s]
       wc2,       &! wc**2
       zeta,      &! dimensionless height used in Monin-Obukhov theory
       zii,       &! convective boundary height [m]
       zldis,     &! reference height "minus" zero displacement heght [m]
       z0mg,      &! roughness length over ground, momentum [m]
       z0hg,      &! roughness length over ground, sensible heat [m]
       z0qg        ! roughness length over ground, latent heat [m]

!----------------------- Dummy argument --------------------------------
! initial roughness length
      if(fsno > 0.)then
         z0mg = zsno
         z0hg = z0mg
         z0qg = z0mg
      else
         z0mg = zlnd
         z0hg = z0mg
         z0qg = z0mg
      endif

! potential temperatur at the reference height
      beta = 1.      ! -  (in computing W_*)
      zii = 1000.    ! m  (pbl height)
      z0ma = z0mg

!-----------------------------------------------------------------------
!     Compute sensible and latent fluxes and their derivatives with respect 
!     to ground temperature using ground temperatures from previous time step.
!-----------------------------------------------------------------------
! Initialization variables
      nmozsgn = 0
      obuold = 0.

      dth   = thm-t_grnd
      dqh   = qm-qg
      dthv  = dth*(1.+0.61*qm)+0.61*th*dqh
      zldis = hu-0.
      if (zldis<=0.) zldis = 0.1

      call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)
 
! Evaluated stability-dependent variables using moz from prior iteration
      niters=6

      !----------------------------------------------------------------
      ITERATION : do iter = 1, niters         ! begin stability iteration
      !----------------------------------------------------------------
         displax = 0.
         call moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
                        ustar,fh2m,fq2m,fm10m,fm,fh,fq)

         tstar = vonkar/fh*dth
         qstar = vonkar/fq*dqh

         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg

         thvstar=tstar+0.61*th*qstar
         zeta=zldis*vonkar*grav*thvstar/(ustar**2*thv)
         if(zeta >= 0.) then     !stable
           zeta = min(2.,max(zeta,1.e-6))
         else                    !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
         endif
         obu = zldis/zeta

         if(zeta >= 0.)then
           um = max(ur,0.1)
         else
           wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
          wc2 = beta*beta*(wc*wc)
           um = sqrt(ur*ur+wc2)
         endif

         if (obuold*obu < 0.) nmozsgn = nmozsgn+1
         if(nmozsgn >= 4) EXIT

         obuold = obu

      !----------------------------------------------------------------
      enddo ITERATION                         ! end stability iteration
      !----------------------------------------------------------------

! Get derivative of fluxes with repect to ground temperature
      ram    = 1./(ustar*ustar/um)
      rah    = 1./(vonkar/fh*ustar) 
      raw    = 1./(vonkar/fq*ustar) 

      raih   = (1.-sigf)*rhoair*cpair/rah
      raiw   = (1.-sigf)*rhoair/raw          
      cgrnds = raih
      cgrndl = raiw*dqgdT
      cgrnd  = cgrnds + htvp*cgrndl

      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! surface fluxes of momentum, sensible and latent 
! using ground temperatures from previous time step
      taux   = -(1.-sigf)*rhoair*us/ram        
      tauy   = -(1.-sigf)*rhoair*vs/ram
      fseng  = -raih*dth
      fevpg  = -raiw*dqh

      fsena  = fseng
      fevpa  = fevpg

! 2 m height air temperature
      tref   = (1.-sigf)*(thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar))
      qref   = (1.-sigf)*( qm + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar))

 end subroutine groundfluxes

 
 subroutine groundtem (itypwat,lb,nl_soil,deltim, &
                       capr,cnfac,csol,porsl,dkdry,dksatu, &
                       sigf,dz_soisno,z_soisno,zi_soisno,&
                       t_soisno,wice_soisno,wliq_soisno,scv,snowdp, &
                       frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg, &
                       imelt,sm,xmf,fact,psi0,bsw,wrsd)

!=======================================================================
! Snow and soil temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
!   the formulation used in SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction
!   in 10 soil layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two neighbor layers
!   (j, j+1) are derived from an assumption that the flux across the interface
!   is equal to that from the node j to the interface and the flux from the
!   interface to the node j+1. The equation is solved using the Crank-Nicholson
!   method and resulted in a tridiagonal system equation.
!
! Phase change (see meltf.F90)
! 
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!=======================================================================

  use precision
  use PhysicalConstants, only : stefnc
  use soil_thermal_parameters

  implicit none

  integer, INTENT(in) :: lb           !lower bound of array
  integer, INTENT(in) :: nl_soil      !upper bound of array
  integer, INTENT(in) :: itypwat      !land water type (0=soil,1=urban or built-up,2=wetland,
                                      !3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) :: deltim      !seconds in a time step [second]
  real(r8), INTENT(in) :: capr        !tuning factor to turn first layer T into surface T
  real(r8), INTENT(in) :: cnfac       !Crank Nicholson factor between 0 and 1

  real(r8), INTENT(in) :: csol(1:nl_soil)  !heat capacity of soil solids [J/(m3 K)]
  real(r8), INTENT(in) :: porsl(1:nl_soil) !soil porosity [-]
  real(r8), INTENT(in) :: psi0 (1:nl_soil) !soil water suction, negative potential [m]
  real(r8), INTENT(in) :: bsw  (1:nl_soil) !clapp and hornbereger "b" parameter [-]
  real(r8), INTENT(in) :: wrsd
  real(r8), INTENT(in) :: dkdry(1:nl_soil) !thermal conductivity of dry soil [W/m-K]
  real(r8), INTENT(in) :: dksatu(1:nl_soil)!thermal conductivity of saturated soil [W/m-K]

  real(r8), INTENT(in) :: sigf     !fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: dz_soisno(lb:nl_soil)   !layer thickiness [m]
  real(r8), INTENT(in) :: z_soisno (lb:nl_soil)   !node depth [m]
  real(r8), INTENT(in) :: zi_soisno(lb-1:nl_soil) !interface depth [m]

  real(r8), INTENT(in) :: sabg             !total solar radiation absorbed by ground [W/m2]
  real(r8), INTENT(in) :: frl              !atmospheric infrared (longwave) radiation [W/m2]
  real(r8), INTENT(in) :: dlrad            !downward longwave radiation blow the canopy [W/m2]
  real(r8), INTENT(in) :: fseng            !sensible heat flux from ground [W/m2]
  real(r8), INTENT(in) :: fevpg            !evaporation heat flux from ground [mm/s]
  real(r8), INTENT(in) :: cgrnd            !deriv. of soil energy flux wrt to soil temp [w/m2/k]
  real(r8), INTENT(in) :: htvp             !latent heat of vapor of water (or sublimation) [j/kg]
  real(r8), INTENT(in) :: emg              !ground emissivity (0.97 for snow,

  real(r8), INTENT(inout) :: t_soisno (lb:nl_soil)   !soil temperature [K]
  real(r8), INTENT(inout) :: wice_soisno(lb:nl_soil) !ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq_soisno(lb:nl_soil) !liqui water [kg/m2]
  real(r8), INTENT(inout) :: scv      !snow cover, water equivalent [mm, kg/m2]
  real(r8), INTENT(inout) :: snowdp   !snow depth [m]

  real(r8), INTENT(out) :: sm         !rate of snowmelt [kg/(m2 s)]
  real(r8), INTENT(out) :: xmf        !total latent heat of phase change of ground water
  real(r8), INTENT(out) :: fact(lb:nl_soil)  !used in computing tridiagonal matrix
  integer,  INTENT(out) :: imelt(lb:nl_soil) !flag for melting or freezing [-]

!------------------------ local variables ------------------------------
  real(r8) sab_soisno(lb:1)   !solar radiation absorbed by snow layers and the 1th ground [W/m2]

  real(r8) cv(lb:nl_soil)     ! heat capacity [J/(m2 K)]
  real(r8) tk(lb:nl_soil)     ! thermal conductivity [W/(m K)]

  real(r8) at(lb:nl_soil)     !"a" vector for tridiagonal matrix
  real(r8) bt(lb:nl_soil)     !"b" vector for tridiagonal matrix
  real(r8) ct(lb:nl_soil)     !"c" vector for tridiagonal matrix
  real(r8) rt(lb:nl_soil)     !"r" vector for tridiagonal solution

  real(r8) fn  (lb:nl_soil)   ! heat diffusion through the layer interface [W/m2]
  real(r8) fn1 (lb:nl_soil)   ! heat diffusion through the layer interface [W/m2]
  real(r8) dzm                ! used in computing tridiagonal matrix
  real(r8) dzp                ! used in computing tridiagonal matrix

  real(r8) t_soisno_bef(lb:nl_soil) ! soil/snow temperature before update
  real(r8) hs                 ! net energy flux into the surface (w/m2)
  real(r8) dhsdt              ! d(hs)/dT
  real(r8) brr(lb:nl_soil)    ! temporay set

  integer i,j

!=======================================================================
      sab_soisno(:) = 0.
      if (lb<1) then
        call solar_radiation_absorption(lb,nl_soil,dz_soisno,sabg,wice_soisno,wliq_soisno,sab_soisno)
      endif
! heat capacity 
      call hCapacity (itypwat,lb,nl_soil,csol,porsl,wice_soisno,wliq_soisno,scv,dz_soisno,cv)

! thermal conductivity
      if(zi_soisno(0) < 0.)then
         print*,'[groundtem],zi_soisno(0)',zi_soisno(0)
         stop
      endif
      
      call hConductivity (itypwat,lb,nl_soil,&
                          dkdry,dksatu,porsl,dz_soisno,z_soisno,zi_soisno,&
                          t_soisno,wice_soisno,wliq_soisno,tk)

! net ground heat flux into the surface and its temperature derivative

      if (lb<1) then
        hs = sab_soisno(lb) + dlrad &
           + (1.-sigf)*emg*frl - emg*stefnc*t_soisno(lb)**4 &
           - (fseng+fevpg*htvp) 
      else
        hs = sabg + dlrad &
         + (1.-sigf)*emg*frl - emg*stefnc*t_soisno(lb)**4 &
         - (fseng+fevpg*htvp) 
      endif

      dhsdT = - cgrnd - 4.*emg * stefnc * t_soisno(lb)**3
      t_soisno_bef(lb:) = t_soisno(lb:)

      j       = lb
      fact(j) = deltim / cv(j) &
              * dz_soisno(j) / (0.5*(z_soisno(j)-zi_soisno(j-1)+capr*(z_soisno(j+1)-zi_soisno(j-1))))

      do j = lb + 1, nl_soil
         fact(j) = deltim/cv(j)
      enddo

      do j = lb, nl_soil - 1
        fn(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z_soisno(j+1)-z_soisno(j))
      enddo
      fn(nl_soil) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
      j     = lb
      dzp   = z_soisno(j+1)-z_soisno(j)
      at(j) = 0.
      bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
      ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = t_soisno(j) + fact(j)*( hs - dhsdT*t_soisno(j) + cnfac*fn(j) )


      do j = lb + 1, nl_soil - 1
         dzm   = (z_soisno(j)-z_soisno(j-1))
         dzp   = (z_soisno(j+1)-z_soisno(j))
         at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
         bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
         ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
         rt(j) = t_soisno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
         ! if this is a snow layer or the top soil layer,
         ! add absorbed solar flux to factor 'rt'
         if (j <= 1) then
            rt(j) = rt(j) + (fact(j)*sab_soisno(j))
         endif         
      end do

      j     =  nl_soil
      dzm   = (z_soisno(j)-z_soisno(j-1))
      at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
      ct(j) = 0.
      rt(j) = t_soisno(j) - cnfac*fact(j)*fn(j-1)

! solve for t_soisno
      i = size(at)
      call tridia (i ,at ,bt ,ct ,rt ,t_soisno) 

!=======================================================================
! melting or freezing 
!=======================================================================

      do j = lb, nl_soil - 1
         fn1(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z_soisno(j+1)-z_soisno(j))
      enddo
      fn1(nl_soil) = 0.

      j = lb
      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

      do j = lb + 1, nl_soil
         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
      enddo

      ! call meltf (lb,nl_soil,deltim, &
      !             fact(lb:),brr(lb:),hs,dhsdT, &
      !             t_soisno_bef(lb:),t_soisno(lb:),wliq_soisno(lb:),wice_soisno(lb:),imelt(lb:), &
      !             scv,snowdp,sm,xmf)
      call meltf (itypwat,lb,nl_soil,deltim,psi0,porsl,wrsd,bsw,&
                   fact(lb:),brr(lb:),hs,dhsdT, dz_soisno(lb:),&
                   t_soisno_bef(lb:),t_soisno(lb:),wliq_soisno(lb:),wice_soisno(lb:),imelt(lb:), &
                   scv,snowdp,sm,xmf,sab_soisno)
!-----------------------------------------------------------------------

 end subroutine groundtem

 SUBROUTINE solar_radiation_absorption(lb,nl_soil,dz_soisno,sabg,wice_soisno,wliq_soisno,sab_soisno)
  use precision

  implicit none

  integer,  INTENT(in) :: lb           !lower bound of array
  integer,  INTENT(in) :: nl_soil      !upper bound of array
  real(r8), INTENT(in) :: dz_soisno(lb:nl_soil)   !layer thickiness [m]
  real(r8), INTENT(in) :: sabg           !total solar radiation absorbed by ground [W/m2]
  real(r8), INTENT(in) :: wice_soisno(lb:nl_soil) !ice lens [kg/m2]
  real(r8), INTENT(in) :: wliq_soisno(lb:nl_soil) !liqui water [kg/m2]
  real(r8), INTENT(out):: sab_soisno(lb:1) !solar radiation absorbed by snow layers and the 1th ground [W/m2]

  real(r8) snow_density(lb:0)
  real(r8) cevis,cenir
  real(r8) ds,sab
  integer i

  ds = 0.001 !assume 1 mm snow grainsize
  cenir = 1.
  sab_soisno(:)=0.
  sab = 0.

  if (sabg.le.0.) then
    sab_soisno(:)=0.
  else
    do i = lb,0
      snow_density(i) = (wice_soisno(i)+wliq_soisno(i))/max(0.0001,dz_soisno(i))
      snow_density(i) = max(snow_density(i),100.)
      cevis           = 0.003795*snow_density(i)*(ds**(-0.5))
      if (i==lb) then
        sab_soisno(i) = sabg*min(1.,exp(-1.*cevis*dz_soisno(lb)))*exp(-1.*cenir*dz_soisno(lb))
      else
        sab_soisno(i) = (sabg-sab)*min(1.,exp(-1.*cevis*dz_soisno(i)))
      endif
      sab = sab + sab_soisno(i)
      if ((sabg-sab).le.0.1) then
        sab_soisno(i) = sab_soisno(i)+(sabg-sab)
        return
      endif
    end do
    sab_soisno(1)=max(sabg-sab,0.)
  endif

 END SUBROUTINE solar_radiation_absorption

 END MODULE MOD_THERMAL