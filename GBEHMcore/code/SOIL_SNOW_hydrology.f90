!#include <define.h>
MODULE SOIL_SNOW_hydrology

!-----------------------------------------------------------------------
  !use precision
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: WATER
  public :: snowwater
  public :: soilwater

! PRIVATE MEMBER FUNCTIONS:
  public :: surfacerunoff
  public :: subsurfacerunoff

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  subroutine WATER (ivt,itypwat,lb,nl_soil,deltim  ,&
             z_soisno    ,dz_soisno   ,zi_soisno   ,&
             sst         ,ss_f,        wfld,       &
             wsat        ,wrsd        ,watern      ,alpha   ,slope  ,&
             length      ,Ds          ,Dr          ,Dg      ,kground,&
             bsw         ,porsl       ,psi0        ,hksati  ,rootr  ,&
             t_soisno    ,wliq_soisno ,wice_soisno ,pg_rain ,sm     ,&
             etr         ,qseva       ,qsdew       ,qsubl   ,qfros  ,&
             rsur        ,rnof        ,qinfl       ,wtfact  ,pondmx ,&
             ssi         ,wimp        ,smpmin      ,zwt     ,wa     ,&
             Drw         ,qcharge     ,qlat        ,qground ,qsmelt ,&
             smf_soisno  ,snoweva)

!=======================================================================
! this is the main subroutine to execute the calculation of
! hydrological processes
!
! Original author: Yongjiu Dai, /09/1999/, /08/2002/, /04/2014/
! Modified       : Hongyi Li, 2015.8 - 2016.12
!
!---------------------------
! CODE Directory
!---------------------------
!--WATER
!    |----snowwater
!    |----surfacerunoff
!    |----soilwater
!    |----subsurfacerunoff
!
!=======================================================================

  use precision
  use PhysicalConstants, only : denice, denh2o, tfrz
  use topoConst,         only : aniks,surfns,sstmaxs

  implicit none

!-----------------------Argument---------- ------------------------------
  integer, INTENT(in) :: &
        lb               , &! lower bound of array
        nl_soil          , &! upper bound of array
        ivt              , &! land cover type of USGS classification or others
        itypwat             ! land water type (0=soil, 1=urban or built-up, 2=wetland,
                            ! 3=land ice, 4=land water bodies, 99=ocean
  real(r8), INTENT(in) :: &
        deltim           , &! time step (s)
        wtfact           , &! fraction of model area with high water table
        pondmx           , &! ponding depth (mm),Depth of water standing on the surface
        ssi              , &! irreducible water saturation of snow
        wimp             , &! water impremeable if porosity less than wimp
        smpmin           , &! restriction for min of soil poten. (mm)
        z_soisno (lb:nl_soil)   , &! layer depth (m)
        dz_soisno(lb:nl_soil)   , &! layer thickness (m)
        zi_soisno(lb-1:nl_soil) , &! interface level below a "z" level (m)
        bsw(1:nl_soil)   , &! Clapp-Hornberger "B"
        porsl(1:nl_soil) , &! saturated volumetric soil water content(porosity)
        psi0(1:nl_soil)  , &! saturated soil suction (mm) (NEGATIVE)
        hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
        rootr(1:nl_soil) , &! root resistance of a layer, all layers add to 1.0
        wfld(1:nl_soil) , & ! water field capacity
        t_soisno(lb:nl_soil), &! soil/snow skin temperature (K)
        pg_rain          , &! rainfall after removal of interception (mm h2o/s)
        sm               , &! snow melt (mm h2o/s)
        etr              , &! actual transpiration (mm h2o/s)
        qseva            , &! ground surface evaporation rate (mm h2o/s)
        qsdew            , &! ground surface dew formation (mm h2o /s) [+]
        qsubl            , &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros               ! surface dew added to snow pack (mm h2o /s) [+]

  real(r8), INTENT(in) :: wsat,wrsd,watern,alpha !GBHM parameter

  real(r8),intent(in):: & ! added by HONGYILI, 2015.12.13
              ss_f,&
              Ds, &       ! depth of topsoil(m)
              Dr, &       ! depth of river (m)
              Dg, &       ! depth of unconfined acquifer (m)
              length, &   ! length of hillslope (m)
              slope,  &   ! slope of hillslope (m)
              kground,&   ! GW hydraulic conductivity (m/s)
              Drw         ! water depth in the river of the flow interval (m)

  real(r8), INTENT(inout) :: &
        wice_soisno(lb:nl_soil) , &! ice lens (kg/m2)
        wliq_soisno(lb:nl_soil) , &! liquid water (kg/m2)
        smf_soisno (-1:nl_soil+1) , &! snowmelt fraction in snow-soil (0-1).Added by HY,2016-11-26
        zwt              , &! the depth from ground (soil) surface to water table [m]
        wa               , &! water storage in aquifer [mm]
        sst                 ! surface water storage (mm)
                            ! Different from old GBHM,Changed to in grid scale,2015-12-13

  real(r8), INTENT(out) :: &
        rsur             , &! surface runoff (mm h2o/s)
        rnof             , &! total runoff (mm h2o/s)
        qinfl            , &! infiltration rate (mm h2o/s)
        qcharge          , &! groundwater recharge (positive to aquifer) [mm/s]
        qlat             , &! lateral runoff(positive = out of soil column) (mm H2O /s)
        qground          , &! ground runoff(positive = out of soil column) (mm H2O /s)
        qsmelt           , &! total snowmelt flux (mm H2O /s)
        snoweva             ! total snow evaporsublimation flux (mm H2O /s)

!====================================================================
!-----------------------Local Variables------------------------------
!
  integer j                 ! loop counter
  integer debug_hydro

  real(r8) :: &
  eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
       hk(1:nl_soil)     , &! hydraulic conductivity [mm h2o/s]
       hksati_m(1:nl_soil),&! hydraulic conductivity [m h2o/s]
       dwat(1:nl_soil)   , &! change in soil water
       gwat              , &! net water input from top (mm/s)
       rsubst            , &! subsurface runoff (mm h2o/s)
       vol_liq(1:nl_soil), &! partitial volume of liquid water in layer
       vol_ice(1:nl_soil), &! partitial volume of ice lens in layer
       icefrac(1:nl_soil), &! ice fraction (-)
       zmm (1:nl_soil)   , &! layer depth (mm)
       dzmm(1:nl_soil)   , &! layer thickness (mm)
       zimm(0:nl_soil)      ! interface level below a "z" level (mm)
  real(r8) :: err_solver, w_sum,dl_drw
  real(r8) :: anik,surface_n,sstmax,GWcs
  real(r8) :: oldice(lb:nl_soil),oldliq(lb:nl_soil)
  real(r8) :: sml_sst,smf_snow,sst0,tmpsst ! temporal VAR for calculate snowmelt contribution
  real(r8) :: soilsm_old,soilsm_new,err_soilsm    ! for checking the balance of snow melt water in soil
  real(r8) :: evpsm,qinsm
  real(r8) debugtmp1,debugtmp2,debugtmp3,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

  debug_hydro = 1
  if (debug_hydro == 1) then
    oldice = wice_soisno(:)
    oldliq = wliq_soisno(:)
  endif

  ! From GBHM defination
  anik      = aniks  (ivt+1)
  surface_n = surfns (ivt+1)
  sstmax    = sstmaxs(ivt+1)
  dl_drw    = Drw
!=======================================================================
! [1] update the liquid water within snow layer and the water onto soil
!=======================================================================
  ! For water balance check, the sum of water in soil column before the calcultion
  w_sum = sum(wliq_soisno(1:)) + sum(wice_soisno(1:)) + wa + sst

  ! For snowmelt water balance checking
  soilsm_old = sum(wliq_soisno(1:)*smf_soisno(1:nl_soil)) + wa*smf_soisno(nl_soil+1)
  qsmelt  = 0.
  snoweva = qsubl
! if(lb<=0)then
! debugtmp1 = sum(wliq_soisno(lb:0)*smf_soisno(-1))+smf_soisno(0)*sst
! elseif (lb>=1)then
!   debugtmp1=smf_soisno(0)*sst
! endif

  if (lb>=1)then
     gwat = pg_rain + sm - qseva !-etr
     if (sm>0.)then
       smf_soisno(-1) = min(1.,sm/(pg_rain+sm))
     else
      smf_soisno(-1) = 0.
     endif
  else
    ! Only deal with the rain through snowpack, should be improved further.
    if (pg_rain>0.) then
      ! smf_snow = 1. ! assume that liquid water in snow are all from snowmelt, including rain leave in snowpack
      smf_snow = smf_soisno(-1) ! assume that liquid water in snow are all from snowmelt, including rain leave in snowpack
      if ((pg_rain*deltim+sum(wliq_soisno(lb:0))+qsdew*deltim) >0.001) then
        smf_soisno(-1) = (smf_snow*sum(wliq_soisno(lb:0))+qsdew*deltim) &
                        /(pg_rain*deltim+sum(wliq_soisno(lb:0))+qsdew*deltim)
      else
        smf_soisno(-1) = 1.
      endif
      smf_soisno(-1) = min(smf_soisno(-1),1.)
    else
      smf_soisno(-1) = 1.
    endif

    snoweva = snoweva + smf_soisno(-1)*qseva

    call snowwater (lb,deltim,ssi,wimp,&
                     pg_rain,qseva,qsdew,qsubl,qfros,&
                     dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),gwat)
  endif

!=======================================================================
! [2] surface runoff and infiltration
!=======================================================================

  if(itypwat<=1)then   ! soil ground only

      ! porosity of soil, partitial volume of ice and liquid
      do j = 1, nl_soil
         vol_ice(j) = min(porsl(j), wice_soisno(j)/(dz_soisno(j)*denice))
         eff_porosity(j) =  max(0.01, porsl(j)-vol_ice(j))
         vol_liq(j) = min(eff_porosity(j), wliq_soisno(j)/(dz_soisno(j)*denh2o))
         if(porsl(j) < 1.e-6)then
            icefrac(j) = 0.
         else
            icefrac(j) = min(1.,vol_ice(j)/porsl(j))
         endif
      enddo

      ! surface runoff including water table and surface staturated area
      rsur = 0.
      sst0 = sst
      sst  = sst + gwat*deltim ! Modified by HY: add variable 'sst'.
      tmpsst = sst

      ! calcultion of snowmelt proporation in SST. HY, 2016-11-28
      !- Begin-----
      if (gwat>0.)then
        sml_sst = (sst0*smf_soisno(0) + gwat*deltim*smf_soisno(-1))!assume fully mixed in sst.
        if(lb>0)then
          snoweva = snoweva + smf_soisno(-1)*qseva
        endif
      else
        sml_sst = (sst0+gwat*deltim)*smf_soisno(0)
        sml_sst = max(0.,sml_sst)
        if(lb>0)then
          snoweva = snoweva + sm + min(abs(gwat),sst0/deltim)*smf_soisno(0)
        endif
      endif
      !- End ------
      if (sst>0.1E-6)then
        smf_soisno(0) = min(1.,sml_sst/sst)
      elseif (sst<0.1E-6 .and. sst>0.1E-10) then
        smf_soisno(0) = min(1.,(sml_sst*1E10)/(sst*1E10))
      endif


      if (tmpsst > 0.) then
        call gbhm_surfacerunoff (deltim,nl_soil,wtfact,wimp,hksati,dz_soisno(1:),&
                            eff_porosity,icefrac,zwt,sst,sstmax,surface_n,slope,length, &
                            rsur,qinfl)
        qsmelt = qsmelt + rsur*smf_soisno(0)
      else
        rsur = 0.
        qinfl= sst/deltim
        sst  = 0.
        smf_soisno(0) = 0.
      endif

! if(lb<=0)then
! debugtmp2 = sum(wliq_soisno(lb:0)*smf_soisno(-1))+smf_soisno(0)*sst
! elseif (lb>=1)then
!   debugtmp2=smf_soisno(0)*sst
! endif
! debugtmp3 = debugtmp1-qsmelt*deltim-snoweva*deltim-max(qinfl,0.)*smf_soisno(0)*deltim
! if(debugtmp3-debugtmp2>1E-8) print*,'surface', debugtmp1,debugtmp2,debugtmp3
!=======================================================================
! [3] determine the change of soil water
!=======================================================================

      ! convert length units from m to mm
      zmm(1:)  = z_soisno(1:) *1000.
      dzmm(1:) = dz_soisno(1:)*1000.
      zimm(0:) = zi_soisno(0:)*1000.

      ! added by HY, 2016-4-20
      GWcs = 0.1
      hksati_m = hksati(1:nl_soil)*0.001 ! transfer from mm to m
      ! call gbhm_recharge(deltim,nl_soil,zwt,Dg,Ds,wa,GWcs,wsat,wrsd,watern,alpha,wfld,&
      !                 hksati_m,dz_soisno(1:nl_soil),wliq_soisno(1:nl_soil),wice_soisno(1:nl_soil),&
      !                 qcharge)
      ! --above added by HY, 2016-4-20---
! debugtmp1 = sum(smf_soisno(1:nl_soil)*wliq_soisno(1:nl_soil)) +wa*smf_soisno(nl_soil+1)+sst*smf_soisno(0)
      call soilwater(nl_soil,deltim,wimp,wrsd,smpmin,&
                     qinfl,etr,z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
                     t_soisno(1:),vol_liq,vol_ice,icefrac,eff_porosity,&
                     porsl,hksati,bsw,psi0,rootr,&
                     zwt,hk,dwat,qcharge,wliq_soisno(1:))

      qinsm = max(0.,qinfl*smf_soisno(0)*deltim)
      ! call smfsoil(nl_soil,deltim,qinfl,etr,dzmm(1:),&
      !               rootr(1:),dwat,qcharge,wliq_soisno(1:),wa,smf_soisno(0:))

      call smsoil(nl_soil,deltim,qinfl,etr,dzmm(1:),&
                    rootr(1:),dwat,qcharge,wliq_soisno(1:),wa,smf_soisno(0:nl_soil+1),evpsm)

      if (tmpsst<=0.) then
        ! calculate snow evaporation after smf_soisno are calculated.
        snoweva = snoweva + evpsm/deltim + sum(etr*rootr(1:nl_soil)*smf_soisno(1:nl_soil))
      endif
! debugtmp2 = sum(smf_soisno(1:nl_soil)*wliq_soisno(1:nl_soil))+(wa+qcharge*deltim)*smf_soisno(nl_soil+1)+sst*smf_soisno(0)
! tmp1 = sum(smf_soisno(1:nl_soil)*wliq_soisno(1:nl_soil))
! tmp2 = (wa+qcharge*deltim)*smf_soisno(nl_soil+1)
! tmp3 = sst*smf_soisno(0)
! if(debugtmp1-evpsm-sum(rootr(:)*smf_soisno(1:nl_soil)*etr*deltim)-debugtmp2>1.E-8) &
! print*,'debugtmp1,',debugtmp1+qinsm-evpsm-debugtmp2, &
! sum(rootr(:)*smf_soisno(1:nl_soil)*etr*deltim)

!=======================================================================
! [4] subsurface runoff and the corrections
!=======================================================================

      ! Modified by HONGYILI,20151213
      CALL subsurfacerunoff(nl_soil,deltim,pondmx,hk,hksati,anik,&
                              wsat,wrsd,watern,alpha,slope,length,&
                              ss_f,Ds,Dr,Dg,kground,t_soisno(1:),&
                              eff_porosity,icefrac,wfld,smf_soisno(0:nl_soil+1),&
                              dz_soisno(1:),zi_soisno(0:),wice_soisno(1:),wliq_soisno(1:),&
                              porsl,psi0,bsw,zwt,wa,dl_drw,&
                              qcharge,rsubst,qlat,qground,qsmelt,sst)

! debugtmp3 = sum(smf_soisno(1:nl_soil)*wliq_soisno(1:nl_soil))+(wa)*smf_soisno(nl_soil+1) +sst*smf_soisno(0)
! tmp4 = sum(smf_soisno(1:nl_soil)*wliq_soisno(1:nl_soil))
! tmp5 = (wa)*smf_soisno(nl_soil+1)
! tmp6 = sst*smf_soisno(0)
! if(debugtmp2-qsmelt*deltim-debugtmp3>1E-8) then
!   print*, 'error in subsurfacerunoff',debugtmp2-debugtmp3,qsmelt*deltim,qcharge,wa,tmp4-tmp1,tmp5-tmp2,tmp6-tmp3
!   pause
! endif
      ! total runoff (mm/s)
      rnof = rsubst + rsur

      ! Renew the ice and liquid mass due to condensation
      if(lb >= 1)then
         ! make consistent with how evap_grnd removed in infiltration
         wliq_soisno(1) = max(0., wliq_soisno(1) + qsdew * deltim)
         wice_soisno(1) = max(0., wice_soisno(1) + (qfros-qsubl) * deltim)
         ! if((wice_soisno(1) + (qfros-qsubl) * deltim)<0.)  then
          ! print*,'potential mass unbalance'
        ! endif
      end if

      if(lb >= 1)then
         err_solver = (sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa+sst) - w_sum &
                    - (gwat+qsdew+qfros-qsubl-rnof-etr)*deltim!
      else
         err_solver = (sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa+sst) - w_sum &
                    - (gwat-etr-rnof)*deltim!
      endif

! #if(defined CLMDEBUG)
     if(abs(err_solver) > 1.e-3 )then
        write(6,*) 'Warning: water balance violation after all soilwater calculation'
        write(6,*) 'err_solver:' ,err_solver
        write(6,*) 'ivt:'      ,ivt,itypwat,lb
        write(6,*) 'wa:'    ,wa,porsl*dzmm
        write(6,*) 'runoff:', rnof*deltim,rsubst*deltim,rsur*deltim,gwat*deltim,qlat*deltim,qsmelt*deltim,qground*deltim,sm*deltim
        write(6,*) 'evaporation:',etr*deltim,qseva*deltim,qsubl*deltim,qsdew*deltim,qfros*deltim,qcharge*deltim
        write(6,*) 'precipitation:',pg_rain*deltim
        write(6,*) 'w_sum0:'     ,w_sum
        write(6,*) 'w_sum1:'   ,(sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa+sst)
        write(6,*) 'mass loss(-):' ,(gwat+qsdew+qfros-qsubl-rnof-etr)*deltim,(gwat-etr-rnof)*deltim
        write(6,*) 'wice_soisno(0:2):',wice_soisno(-1:2)
        write(6,*) 'wliq_soisno(0:2):',wliq_soisno(-1:2)
        write(6,*) 'oldice:',oldice(-1:2)
        write(6,*) 'oldliq:',oldliq(-1:2)
     endif
! #endif
  ! For snowmelt water balance checking
  ! soilsm_new = sum(wliq_soisno(1:)*smf_soisno(1:nl_soil)) + wa*smf_soisno(nl_soil+1)
  ! if (qinfl>0.) then
  !   err_soilsm = (soilsm_old+qinsm-qsmelt*deltim- &
  !               etr*deltim*sum(rootr(1:)*smf_soisno(1:nl_soil))) - soilsm_new
  ! else
  !   err_soilsm = (soilsm_old-evpsm-qsmelt*deltim- &
  !               etr*deltim*sum(rootr(1:)*smf_soisno(1:nl_soil)/sum(rootr(1:)))) - soilsm_new
  ! endif
  ! if(abs(err_soilsm) > 1.e-3 )then
  !   write(6,*) 'Warning: water balance violation in soil snowmelt calculation'
  !   write(6,*) 'err_soilsm',err_soilsm
  ! endif



!=======================================================================
! [6] assumed hydrological scheme for the wetland and glacier
!=======================================================================

  else
      if(itypwat==2)then        ! WETLAND
         qinfl = 0.
         rsur = max(0.,gwat)
         rsubst = 0.
         rnof = 0.
         do j = 1, nl_soil
            if(t_soisno(j)>tfrz)then
              wice_soisno(j) = 0.0
              wliq_soisno(j) = porsl(j)*dz_soisno(j)*1000.
            endif
         enddo
      endif
      if(itypwat==3)then        ! LAND ICE
         rsur = max(0.0,gwat)
         qinfl = 0.
         rsubst = 0.
         rnof = rsur
         qsmelt = rsur ! Added by HY.,20161225
         wice_soisno(1:nl_soil) = dz_soisno(1:nl_soil)*1000.
         wliq_soisno(1:nl_soil) = 0.0
      endif

      wa  = 0.1*Dg*1000.!4800
      zwt = 0.
      qcharge = 0.

  endif

!-----------------------------------------------------------------------

  end subroutine WATER



  subroutine snowwater (lb,deltim,ssi,wimp, &
                        pg_rain,qseva,qsdew,qsubl,qfros, &
                        dz_soisno,wice_soisno,wliq_soisno,qout_snowb)

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, /09/1999; /04/2014
!
! Water flow wihtin snow is computed by an explicit and non-physical based scheme,
! which permits a part of liquid water over the holding capacity (a tentative value
! is used, i.e., equal to 0.033*porosity) to percolate into the underlying layer,
! except the case of that the porosity of one of the two neighboring layers is
! less than 0.05, the zero flow is assumed. The water flow out of the bottom
! snow pack will participate as the input of the soil water and runoff.
!
!-----------------------------------------------------------------------

  use precision
  use PhysicalConstants, only : denice, denh2o  ! physical constant
  implicit none

!----------------------- dummy argument --------------------------------
  integer, INTENT(in) :: &
        lb          ! lower bound of array

  real(r8), INTENT(in) :: &
        deltim,    &! seconds in a time step (s)
        ssi,       &! irreducible water saturation of snow
        wimp,      &! water impremeable if porosity less than wimp
        dz_soisno(lb:0),  &! layer thickness (m)
        pg_rain,   &! rainfall after removal of interception (mm h2o/s)
        qseva,     &! ground surface evaporation rate (mm h2o/s)
        qsdew,     &! ground surface dew formation (mm h2o /s) [+]
        qsubl,     &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros       ! surface dew added to snow pack (mm h2o /s) [+]

  real(r8), INTENT(inout) :: &
        wice_soisno(lb:0),&! ice lens (kg/m2)
        wliq_soisno(lb:0)  ! liquid water (kg/m2)

  real(r8), INTENT(out) :: &
        qout_snowb  ! rate of water out of snow bottom (mm/s)

!----------------------- local variables --------------------------------
  integer j         ! k do loop/array indices

  real(r8) :: &
       qin,        &! water flow into the elmement (mm/s)
       qout,       &! water flow out of the elmement (mm/s)
       zwice,      &! the sum of ice mass of snow cover (kg/m2)
       wgdif,      &! ice mass after minus sublimation
    vol_liq(lb:0), &! partitial volume of liquid water in layer
    vol_ice(lb:0), &! partitial volume of ice lens in layer
 eff_porosity(lb:0) ! effective porosity = porosity - vol_ice

!=======================================================================
! renew the mass of ice lens (wice_soisno) and liquid (wliq_soisno) in the surface snow layer,
! resulted by sublimation (frost) / evaporation (condense)

      wgdif = wice_soisno(lb) + (qfros - qsubl)*deltim
      wice_soisno(lb) = wgdif
      if(wgdif < 0.)then
         wice_soisno(lb) = 0.
         wliq_soisno(lb) = wliq_soisno(lb) + wgdif
      endif
      wliq_soisno(lb) = wliq_soisno(lb) + (pg_rain + qsdew - qseva)*deltim
      wliq_soisno(lb) = max(0., wliq_soisno(lb))

! Porosity and partitial volume
      do j = lb, 0
         vol_ice(j) = min(1., wice_soisno(j)/(dz_soisno(j)*denice))
         eff_porosity(j) = max(0.01, 1. - vol_ice(j))
         vol_liq(j) = min(eff_porosity(j), wliq_soisno(j)/(dz_soisno(j)*denh2o))
      enddo

! Capillary force within snow could be two or more orders of magnitude
! less than those of gravity, this term may be ignored.
! Here we could keep the garavity term only. The genernal expression
! for water flow is "K * ss**3", however, no effective paramterization
! for "K". Thus, a very simple treatment (not physical based) is introduced:
! when the liquid water of layer exceeds the layer's holding
! capacity, the excess meltwater adds to the underlying neighbor layer.

      qin = 0.
      do j= lb, 0
         wliq_soisno(j) = wliq_soisno(j) + qin

         if(j <= -1)then
         ! no runoff over snow surface, just ponding on surface
           if(eff_porosity(j)<wimp .OR. eff_porosity(j+1)<wimp)then
             qout = 0.
           else
             qout = max(0.,(vol_liq(j)-ssi*eff_porosity(j))*dz_soisno(j))
             qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dz_soisno(j+1))
           endif
         else
           qout = max(0.,(vol_liq(j)-ssi*eff_porosity(j))*dz_soisno(j))
         endif

         qout = qout*1000.
         wliq_soisno(j) = wliq_soisno(j) - qout
         qin = qout

      enddo

      qout_snowb = qout/deltim

  end subroutine snowwater

  subroutine gbhm_surfacerunoff (deltim,nl_soil,wtfact,wimp,hksati,dz_soisno,&
                            eff_porosity,icefrac,zwt,sst,sstmax,surface_n,slp,length, &
                            rsur,qinfl)

!=======================================================================
! the original code was provide by Robert E. Dickinson based on following clues:
! a water table level determination level added including highland and
! lowland levels and fractional area of wetland (water table above the surface.
! Runoff is parametrized from the lowlands in terms of precip incident on
! wet areas and a base flow, where these are estimated using ideas from TOPMODEL.
!
! Author : Yongjiu Dai, 07/29/2002, Guoyue Niu, 06/2012
! Modified by Hongyi Li, 03/2016
!=======================================================================

  use precision
  implicit none

!-----------------------Arguments---------------------------------------

  integer, INTENT(in) :: nl_soil   ! number of soil layers
  real(r8), INTENT(in) :: &
        deltim       , &
        wtfact       , &! fraction of model area with high water table
        wimp         , &! water impremeable if porosity less than wimp
    hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
 dz_soisno(1:nl_soil), &! layer thickness (m)
 eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
   icefrac(1:nl_soil), &! ice fraction (-)
         zwt         , &! the depth from ground (soil) surface to water table [m]
         surface_n   , &
         slp         , &
         sstmax      , &
         length
  real(r8), INTENT(inout) :: sst  ! surface liquid water storage (mm)

  real(r8), INTENT(out) :: rsur , & ! surface runoff (mm h2o/s)
                           qinfl    ! infiltration into surface soil layer (mm h2o/s)

!-----------------------Local Variables---------------------------------

  real(r8) qinmax       ! maximum infiltration capability
  real(r8) fsat         ! fractional area with water table at surface
  real(r8) sstrate      ! sst/deltim
  real(r8) q_hillslope  ! flow at one hillslope,(m3/m)
  real(r8) water_depth  ! (mm)
  real(r8) soil_con_f   ! soil conservation factor (-)
  real(r8) waterhead
  real(r8) detension
  real(r8) rsurtmp,power,sst0

  real(r8), parameter :: fff = 0.5   ! runoff decay factor (m-1)

!-----------------------End Variable List-------------------------------
      sst0 = sst
! ---------------------------------------------------------------------
! 1. Determining infiltration and update sst (partitially from CLM)
! ---------------------------------------------------------------------
      sstrate = sst/deltim
      ! fraction of saturated area
      ! fsat = wtfact*min(1.0,exp(-0.5*fff*zwt))
      ! fsat = 0.
      ! Maximum infiltration capacity
      ! qinmax = minval(10.**(-6.0*icefrac(1:3))*hksati(1:3))
      ! qinmax = minval(hksati(1:3))
      qinmax = max(hksati(1),0.)

      if(eff_porosity(1)<wimp) qinmax = 0.
      ! Surface runoff
      ! rsurtmp = fsat*max(0.0,sstrate) + (1.-fsat)*max(0.,sstrate-qinmax)
      rsurtmp = max(0.,sstrate-qinmax)
      ! rsurtmp = min(sstrate,rsurtmp)!,0.)!Modified by HY, option value 0 is deleted, because I think rsur should be positive.
      ! rsurtmp = max(rsurtmp,0.)
      ! infiltration into surface soil layer
      qinfl = sstrate - rsurtmp
      sst   = sst - qinfl*deltim
      if (sst.gt.-1.0E-10 .and. sst.lt.0.) then
        qinfl = qinfl+sst/deltim
        sst = 0.
      elseif (sst.lt.-1.0E-10 ) then
        print*, 'debug in line 482: negative sst value',sst, qinfl*deltim
        qinfl = qinfl+sst/deltim
        sst = 0.
      endif
! ----------------------------------------------------------------------
! 2. surface routing: steady constant sheet flow (From GBHM)
! ----------------------------------------------------------------------
      soil_con_f  = 1.
      detension   = sstmax*soil_con_f
      detension   = amax1(3.0, detension)
      water_depth = amax1(0.0, (sst-detension) )
      ! water_depth = max(0.0, sst) !Debug, by HY. To testify whether the detension promote evaporation or not.
      q_hillslope = 0.0
      if( water_depth .ge. 0.0 ) then
        sst         = sst-water_depth
        !transfer to surface runoff (m)
        water_depth = 0.001 * water_depth
        waterhead   = slp! + water_depth/length
        power       = 1.6667
        q_hillslope = deltim * sqrt(waterhead)*water_depth**power/surface_n
        if(q_hillslope .le. 0.1E-20) q_hillslope = 0.0
        ! for agricultural fields
        ! if(iland.eq.6) then
            ! q_hillslope = q_hillslope * amax1(0.3, (1.0-LAI/LAImax))
        ! endif
        q_hillslope = amin1(q_hillslope, water_depth*length)
        ! q_hillslope = water_depth*length!Debug by HY,2016-4-4,assuming all water depth flow out
        water_depth = water_depth - q_hillslope/length

        ! update surface storage
        sst = sst+1000.0*water_depth
        water_depth = 0.0
        rsur = (q_hillslope/length/deltim)*1000.
        if (abs(rsur+sst/deltim+qinfl-sst0/deltim).gt.0.00001) &
            print*,'debug: balance error in line 518'
      else
        rsur = 0.0
      endif

  end subroutine gbhm_surfacerunoff

  subroutine surfacerunoff (nl_soil,wtfact,wimp,bsw,porsl,psi0,hksati,&
                            z_soisno,dz_soisno,zi_soisno,&
                            eff_porosity,icefrac,zwt,gwat,rsur,qinfl)

!=======================================================================
! the original code was provide by Robert E. Dickinson based on following clues:
! a water table level determination level added including highland and
! lowland levels and fractional area of wetland (water table above the surface.
! Runoff is parametrized from the lowlands in terms of precip incident on
! wet areas and a base flow, where these are estimated using ideas from TOPMODEL.
!
! Author : Yongjiu Dai, 07/29/2002, Guoyue Niu, 06/2012
! Major modified by Hongyi Li, 2015
!=======================================================================

  use precision
  implicit none

!-----------------------Arguments---------------------------------------

  integer, INTENT(in) :: nl_soil   ! number of soil layers
  real(r8), INTENT(in) :: &
        wtfact,        &! fraction of model area with high water table
        wimp,          &! water impremeable if porosity less than wimp
       bsw(1:nl_soil), &! Clapp-Hornberger "B"
     porsl(1:nl_soil), &! saturated volumetric soil water content(porosity)
      psi0(1:nl_soil), &! saturated soil suction (mm) (NEGATIVE)
    hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
  z_soisno(1:nl_soil), &! layer depth (m)
 dz_soisno(1:nl_soil), &! layer thickness (m)
 zi_soisno(0:nl_soil), &! interface level below a "z" level (m)
 eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
   icefrac(1:nl_soil), &! ice fraction (-)
        gwat,          &! net water input from top
         zwt            ! the depth from ground (soil) surface to water table [m]

  real(r8), INTENT(out) :: rsur, & ! surface runoff (mm h2o/s)
                           qinfl   ! infiltration into surface soil layer (mm h2o/s)

!-----------------------Local Variables---------------------------------

  real(r8) qinmax       ! maximum infiltration capability
  real(r8) fsat         ! fractional area with water table at surface

  real(r8), parameter :: fff = 0.5   ! runoff decay factor (m-1)

!-----------------------End Variable List-------------------------------

!  fraction of saturated area
      fsat = wtfact*min(1.0,exp(-0.5*fff*zwt))
      ! fsat = 0.

! Maximum infiltration capacity
      qinmax = minval(10.**(-6.0*icefrac(1:3))*hksati(1:3))
      ! qinmax = minval(hksati(1:3))
      ! if(eff_porosity(1)<wimp) qinmax = 0.

! Surface runoff
      rsur = fsat*max(0.0,gwat) + (1.-fsat)*max(0.,gwat-qinmax)
      rsur = min(gwat,rsur)!,0.)!Modified by Hongyi Li, option value 0 is deleted.because I think rsur should be positive.
! infiltration into surface soil layer
      qinfl = gwat - rsur
  end subroutine surfacerunoff

  subroutine soilwater(nl_soil,deltim,wimp,wrsd,smpmin,&
                       qinfl,etr,z_soisno,dz_soisno,zi_soisno,&
                       t_soisno,vol_liq,vol_ice,icefrac,eff_porosity,&
                       porsl,hksati,bsw,psi0,rootr,&
                       zwt,hk,dwat,qcharge,wliq_soisno)
!-----------------------------------------------------------------------
! REVISION HISTORY:
! Original author : Yongjiu Dai
! Modified by :     Hongyi Li
!
! 2018-03-02, Hongyi Li:
!   Parameterization in CLM4.5 is coupled with the original CoLM codes.
!   Now the aquifer flow at vertical direction is solved with the soil
!   flow processes. In the old codes, qcharge is not related to the water
!   table. The new codes deal with the aquifer as an important factor in
!   determining QCHARGE. This parameterization scheme is similar with other
!   hydrological model such as GBEHM.
!   The variable names are yet the same with the old CoLM codes, and the
!   whole codes structure is adjusted refer to CLM 4.5 and CoLM.
!-----------------------------------------------------------------------
! INTRODUCTION:
! Soil moisture is predicted from a 10-layer model (as with soil
! temperature), in which the vertical soil moisture transport is governed
! by infiltration, runoff, gradient diffusion, gravity, and root
! extraction through canopy transpiration. The net water applied to the
! surface layer is the snowmelt plus precipitation plus the throughfall
! of canopy dew minus surface runoff and evaporation.
!
! The vertical water flow in an unsaturated porous media is described by
! Darcy's law, and the hydraulic conductivity and the soil negative
! potential vary with soil water content and soil texture based on the work
! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
! integrated over the layer thickness, in which the time rate of change in
! water mass must equal the net flow across the bounding interface, plus the
! rate of internal source or sink. The terms of water flow across the layer
! interfaces are linearly expanded by using first-order Taylor expansion.
! The equations result in a tridiagonal system equation.
!
! Note: length units here are all millimeter
! (in temperature subroutine uses same soil layer
! structure required but lengths are m)
!
! Richards equation:
!
! d wat     d     d psi
! ----- =  -- [ k(----- - 1) ] + S
!   dt     dz       dz
!
! where: wat = volume of water per volume of soil (mm**3/mm**3)
! psi = soil matrix potential (mm)
! dt  = time step (s)
! z   = depth (mm) (positive downward)
! dz  = thickness (mm)
! qin = inflow at top (mm h2o /s)
! qout= outflow at bottom (mm h2o /s)
! s   = source/sink flux (mm h2o /s)
! k   = hydraulic conductivity (mm h2o /s)
!
!                       d qin                  d qin
! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
!                       d wat(j-1)             d wat(j)
!                ==================|=================
!                                  < qin
!
!                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
!
!                                  > qout
!                ==================|=================
!                        d qout               d qout
! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!                        d wat(j)             d wat(j+1)
!
!
! Solution: linearize k and psi about d wat and use tridiagonal
! system of equations to solve for d wat,
! where for layer j
!
!
! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
!-----------------------------------------------------------------------
    use precision
    use PhysicalConstants , only : grav,hfus,tfrz,denh2o,denice

    IMPLICIT NONE

    integer, INTENT(in) :: nl_soil  ! number of soil layers
    real(r8), INTENT(in) :: deltim  ! land model time step (sec)
    real(r8), INTENT(in) :: wimp    ! water impremeable if porosity less than wimp
    real(r8), INTENT(in) :: wrsd
    real(r8), INTENT(in) :: smpmin  ! restriction for min of soil potential (mm)

    real(r8), INTENT(in) :: qinfl   ! infiltration (mm H2O /s)
    real(r8), INTENT(in) :: etr     ! vegetation transpiration (mm H2O/s) (+ = to atm)

    real(r8), INTENT(in) :: z_soisno (1:nl_soil) ! layer depth (m)
    real(r8), INTENT(in) :: dz_soisno(1:nl_soil) ! layer thickness (m)
    real(r8), INTENT(in) :: zi_soisno(0:nl_soil) ! interface level below a "z" level (m)

    real(r8), INTENT(in) :: t_soisno (1:nl_soil) ! soil temperature (Kelvin)
    real(r8), INTENT(in) :: vol_liq  (1:nl_soil) ! liquid volumetric water content
    real(r8), INTENT(in) :: vol_ice  (1:nl_soil) ! ice volumetric water content
    real(r8), INTENT(in) :: icefrac  (1:nl_soil)
    real(r8), INTENT(in) :: eff_porosity(1:nl_soil) ! effective porosity = porosity - vol_ice

    real(r8), INTENT(in) :: porsl  (1:nl_soil) ! volumetric soil water at saturation (porosity)
    real(r8), INTENT(in) :: hksati (1:nl_soil) ! hydraulic conductivity at saturation (mm H2O /s)
    real(r8), INTENT(in) :: bsw    (1:nl_soil) ! Clapp and Hornberger "b"
    real(r8), INTENT(in) :: psi0   (1:nl_soil) ! minimum soil suction (mm) [-]
    real(r8), INTENT(in) :: rootr  (1:nl_soil) ! effective fraction of roots in each soil layer
    real(r8), INTENT(in) :: zwt                ! the depth from ground (soil) surface to water table [m]

    real(r8), intent(out) :: hk(1:nl_soil)     ! hydraulic conductivity [mm h2o/s]
    real(r8), intent(out) :: dwat(1:nl_soil)   ! change of soil water [m3/m3]
    real(r8), INTENT(out) :: qcharge         ! aquifer recharge rate (positive to aquifer) (mm/s)
    real(r8), INTENT(inout) :: wliq_soisno(1:nl_soil)
    !
    ! local arguments
    !
    integer  :: j                 ! do loop indices
    real(r8) :: amx(1:nl_soil+1)    ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(1:nl_soil+1)    ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(1:nl_soil+1)    ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(1:nl_soil+1)    ! "r" forcing term of tridiagonal matrix
    real(r8) :: zmm(1:nl_soil+1)  ! layer depth [mm]
    real(r8) :: dzmm(1:nl_soil+1) ! layer thickness [mm]
    real(r8) :: zimm(0:nl_soil)   ! layer interface depth [mm]
    real(r8) :: zwtmm             ! water table depth [mm]
    real(r8) :: zq(1:nl_soil+1)     ! equilibrium matric potential for each layer [mm]
    real(r8) :: vol_eq(1:nl_soil+1) ! equilibrium volumetric water content
    ! real(r8) :: den(1:nl_soil)    ! used in calculating qin, qout
    real(r8) :: den    ! used in calculating qin, qout
    real(r8) :: alpha  ! used in calculating qin, qout
    ! real(r8) :: num    ! used in calculating qin, qout
    ! real(r8) :: alpha(1:nl_soil)  ! used in calculating qin, qout
    real(r8) :: qin(1:nl_soil+1)    ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout(1:nl_soil+1)   ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: dqidw0(1:nl_soil+1) ! d(qin)/d(vol_liq(j-1))
    real(r8) :: dqidw1(1:nl_soil+1) ! d(qin)/d(vol_liq(j))
    real(r8) :: dqodw1(1:nl_soil+1) ! d(qout)/d(vol_liq(j))
    real(r8) :: dqodw2(1:nl_soil+1) ! d(qout)/d(vol_liq(j+1))
    real(r8) :: dsmpdw(1:nl_soil+1) ! d(smp)/d(vol_liq)
    real(r8) :: s_node            ! soil wetness
    real(r8) :: s1                ! "s" at interface of layer
    real(r8) :: s2                ! k*s**(2b+2)
    real(r8) :: smp(1:nl_soil)    ! soil matrix potential [mm]
    ! real(r8) :: hk(1:nl_soil)     ! hydraulic conductivity [mm h2o/s]
    real(r8) :: dhkdw(1:nl_soil) ! d(hk)/d(vol_liq)
    ! real(r8) :: dhkdw1(1:nl_soil) ! d(hk)/d(vol_liq(j))
    ! real(r8) :: dhkdw2(1:nl_soil) ! d(hk)/d(vol_liq(j+1))
    real(r8) :: dwat_soi_acq(1:nl_soil+1)
    real(r8) :: imped(1:nl_soil)  !
    real(r8) :: vol_liq_zwt       ! volumetric water content at water table depth
    real(r8) :: errorw            ! mass balance error for this time step
    real(r8) :: sdamp             ! extrapolates soiwat dependence of evaporation

    integer  :: jwt               ! index of the soil layer right above the water table (-)
    integer  :: nl ! the number of soil and acquifer layers

    real(r8), parameter :: e_ice=6.0      !soil ice impedance factor
    real(r8) :: hksat_inter,qmax,qmin,temq(1:nl_soil),smp1
    real(r8) :: tempi,temp0,voleq1,dsmpdw_aquifer
    real(r8) :: wh_zwt,ka,wh
!-----------------------------------------------------------------------
    ! the nl is the temperature layers which is different from wetness layers (soil layer)
    ! sdamp is an ajustable parameter in soil surface
    nl = nl_soil
    sdamp = 0._r8

    !compute jwt index
    ! The layer index of the first unsaturated layer,
    ! i.e., the layer right above the water table
    jwt = nl_soil
    ! allow jwt to equal zero when zwt is in top layer
    do j = 1, nl_soil
       if(zwt <= zi_soisno(j)) then
          jwt = j-1
          exit
       end if
    enddo

    ! transfer depth metric from m to mm, by HYLI
    zmm (1:nl_soil) = z_soisno (1:nl_soil)*1000.
    dzmm(1:nl_soil) = dz_soisno(1:nl_soil)*1000.
    zimm(1:nl_soil) = zi_soisno(1:nl_soil)*1000.
    zimm(0) = 0.0
    zwtmm   = zwt*1000.
    ! aquifer layer, by HYLI
    zmm(nl_soil+1) = 0.5*(1000.*zwt + zmm(nl_soil))
    if(jwt < nl_soil) then
     dzmm(nl_soil+1) = dzmm(nl_soil)
    else
     dzmm(nl_soil+1) = (1000.*zwt - zmm(nl_soil))
    end if


    ! the following codes is modified from CLM 4.5, by HYLI
    ! compute Volumetric Water Content at water table depth (vol_liq_zwt)
    !  (mainly for case when t < tfrz)
    !  this will only be used when zwt is below the soil column
    ! NOTE: 'psi0' in CoLM is equivalent to '-1*sucsat' in CLM 4.5
    !
    ! BEGIN CALC vol_liq_zwt
    vol_liq_zwt = porsl(nl_soil)
    if(t_soisno(jwt+1) < tfrz) then
      vol_liq_zwt = vol_liq(nl_soil)
      do j = nl_soil,nl
         if(zwt <= zi_soisno(j)) then
            smp1 = hfus*(tfrz-t_soisno(j))/(grav*t_soisno(j)) * 1000._r8  !(mm)
            smp1 = max(-psi0(nl_soil),smp1)
            vol_liq_zwt = porsl(nl_soil)*(smp1/(-psi0(nl_soil)))**(-1._r8/bsw(nl_soil))
            ! for temperatures close to tfrz, limit vol_liq_zwt to total water content
            vol_liq_zwt = min(vol_liq_zwt, 0.5*(porsl(nl_soil) + vol_liq(nl_soil)) )
            exit
         endif
      enddo
    endif
    ! END CALC vol_liq_zwt

    ! calculate the equilibrium water content (zq) based on the water table depth
    CALC_zq: do j=1,nl_soil
      if ((zwtmm .le. zimm(j-1))) then
         vol_eq(j) = porsl(j)
      ! use the weighted average from the saturated part (depth > wtd) and the equilibrium solution for the
      ! rest of the layer
      else if ((zwtmm .lt. zimm(j)) .and. (zwtmm .gt. zimm(j-1))) then
         tempi = 1.0_r8
         temp0 = (((psi0(j)-zwtmm+zimm(j-1))/psi0(j)))**(1._r8-1._r8/bsw(j))
         voleq1 = psi0(j)*porsl(j)/(1._r8-1._r8/bsw(j))/(zwtmm-zimm(j-1))*(tempi-temp0)
         vol_eq(j) = (voleq1*(zwtmm-zimm(j-1)) + porsl(j)*(zimm(j)-zwtmm))/(zimm(j)-zimm(j-1))
         vol_eq(j) = min(porsl(j),vol_eq(j))
         vol_eq(j) = max(vol_eq(j),0.0_r8)
      else
         tempi = (((psi0(j)-zwtmm+zimm(j))/psi0(j)))**(1._r8-1._r8/bsw(j))
         temp0 = (((psi0(j)-zwtmm+zimm(j-1))/psi0(j)))**(1._r8-1._r8/bsw(j))
         vol_eq(j) = psi0(j)*porsl(j)/(1._r8-1._r8/bsw(j))/(zimm(j)-zimm(j-1))*(tempi-temp0)
         vol_eq(j) = max(vol_eq(j),0.0_r8)
         vol_eq(j) = min(porsl(j),vol_eq(j))
      endif
      zq(j) = psi0(j)*(max(vol_eq(j)/porsl(j),0.01_r8))**(-bsw(j))
      zq(j) = max(smpmin, zq(j))
    end do CALC_zq

    ! If water table is below soil column calculate zq for the 11th layer
    j = nl_soil
    if(jwt == nl_soil) then
      tempi = 1._r8
      temp0 = (((psi0(j)-zwtmm+zimm(j))/psi0(j)))**(1._r8-1._r8/bsw(j))
      vol_eq(j+1) = psi0(j)*porsl(j)/(1._r8-1._r8/bsw(j))/(zwtmm-zimm(j))*(tempi-temp0)
      vol_eq(j+1) = max(vol_eq(j+1),0.0_r8)
      vol_eq(j+1) = min(porsl(j),vol_eq(j+1))
      zq(j+1) = psi0(j)*(max(vol_eq(j+1)/porsl(j),0.01_r8))**(-bsw(j))
      zq(j+1) = max(smpmin, zq(j+1))
    end if


    CALC_hk_dhkdw: do j = 1, nl_soil
      ! Hydraulic conductivity and soil matric potential and their derivatives
      ! ---[hk, dhkdw]
      s1 = 0.5_r8*(vol_liq(j) + vol_liq(min(nl_soil, j+1))) / &
       (0.5_r8*(porsl(j)+porsl(min(nl_soil, j+1))))
      s1 = min(1._r8, s1)
      s2 = hksati(j)*s1**(2._r8*bsw(j)+2._r8)
      imped(j)=10.**(-e_ice*(0.5*(icefrac(j)+icefrac(min(nl_soil,j+1)))))
      hk(j) = imped(j)*s1*s2
      dhkdw(j) = imped(j)*(2._r8*bsw(j)+3._r8)*s2* &
                (1._r8/(porsl(j)+porsl(min(nl_soil, j+1))))

      ! Compute matric potential and derivative based on liquid water content only
      ! 'I disabled the codes in range (t_soisno(j)>tfrz). HYLI'
      ! ---[smp, dsmpdw]
      ! if(t_soisno(j)>tfrz) then
      if(porsl(j)<1.e-6)then     ! bed rock
          s_node = 0.001
          smp(j) = psi0(j)
          dsmpdw(j) = 0.
      else
          s_node = max(vol_liq(j)/porsl(j),0.01)
          s_node = min(1.0,s_node)
          smp(j) = psi0(j)*s_node**(-bsw(j))
          smp(j) = max(smpmin,smp(j))
          dsmpdw(j) = -bsw(j)*smp(j)/(s_node*porsl(j))
      endif
      ! else
      !    ! when ice is present, the matric potential is only related to temperature
      !    ! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
      !    ! Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm
      !    smp(j) = 1.e3 * 0.3336e6/9.80616*(t_soisno(j)-tfrz)/t_soisno(j)
      !    smp(j) = max(smpmin, smp(j))        ! Limit soil suction
      !    dsmpdw(j) = 0.
      ! endif
    end do CALC_hk_dhkdw

    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node j=1 (top)

    j = 1
    qin(j) = qinfl

    den       = (zmm(j+1)-zmm(j))
    alpha     = ( (smp(j+1)-smp(j)) - (zq(j+1)-zq(j)) )/den
    qout(j)   = -hk(j)*alpha
    dqodw1(j) = -(alpha*dhkdw(j) - hk(j)*dsmpdw(j)/den)
    dqodw2(j) = -(alpha*dhkdw(j) + hk(j)*dsmpdw(j+1)/den)

    amx(j) = 0._r8
    bmx(j) = dzmm(j)*(sdamp+1._r8/deltim) + dqodw1(j)
    cmx(j) = dqodw2(j)
    rmx(j) = qin(j) - qout(j) - etr*rootr(j)

    ! Nodes j=2 to j=nl_soil-1
    do j = 2, nl_soil - 1
      den       = (zmm(j) - zmm(j-1))
      alpha     = ( (smp(j)-smp(j-1)) - (zq(j)-zq(j-1)) )/den
      qin(j)    = -hk(j-1)*alpha
      dqidw0(j) = -(alpha*dhkdw(j-1) - hk(j-1)*dsmpdw(j-1)/den)
      dqidw1(j) = -(alpha*dhkdw(j-1) + hk(j-1)*dsmpdw(j)/den)

      den       = (zmm(j+1)-zmm(j))
      alpha     = ( (smp(j+1)-smp(j)) - (zq(j+1)-zq(j)) )/den
      qout(j)   = -hk(j)*alpha
      dqodw1(j) = -(alpha*dhkdw(j) - hk(j)*dsmpdw(j)/den)
      dqodw2(j) = -(alpha*dhkdw(j) + hk(j)*dsmpdw(j+1)/den)

      amx(j) = -dqidw0(j)
      bmx(j) =  dzmm(j)/deltim - dqidw1(j) + dqodw1(j)
      cmx(j) =  dqodw2(j)
      rmx(j) =  qin(j) - qout(j) - etr*rootr(j)
    end do

    ! Node j=nl_soil (bottom)
    j = nl_soil
    !--- set the inflow related coefficients which are not changed with water table position.
    den       = (zmm(j) - zmm(j-1))
    alpha     = ( (smp(j)-smp(j-1)) - (zq(j)-zq(j-1)) )/den
    qin(j)    = -hk(j-1)*alpha
    dqidw0(j) = -(alpha*dhkdw(j-1) - hk(j-1)*dsmpdw(j-1)/den)
    dqidw1(j) = -(alpha*dhkdw(j-1) + hk(j-1)*dsmpdw(j)/den)

    if(j > jwt) then ! water table is in soil column
      qout(j)   = 0.
      dqodw1(j) = 0.
      ! dqodw2(j) = 0.

      amx(j) = -dqidw0(j)
      bmx(j) =  dzmm(j)/deltim - dqidw1(j) + dqodw1(j)
      cmx(j) =  0._r8
      rmx(j) =  qin(j) - qout(j) - etr*rootr(j)
      ! next set up aquifer layer; hydrologically inactive
      amx(j+1) = 0._r8
      bmx(j+1) = dzmm(j+1)/deltim
      cmx(j+1) = 0._r8
      rmx(j+1) = 0._r8
    else ! water table is below soil column
      ! compute aquifer soil moisture as average of layer 10 and saturation
      s_node = max(0.5*((vol_liq_zwt+vol_liq(j))/porsl(j)), 0.01_r8)
      s_node = min(1.0_r8, s_node)

      ! compute smp for aquifer layer
      smp1 = psi0(j)*s_node**(-bsw(j))
      smp1 = max(smpmin, smp1)
      ! compute dsmpdw for aquifer layer
      dsmpdw_aquifer = -bsw(j)*smp1/(s_node*porsl(j))

      ! first set up bottom layer of soil column

      den       = (zmm(j+1)-zmm(j))
      alpha     = ( (smp1-smp(j)) - (zq(j+1)-zq(j)) )/den
      qout(j)   = -hk(j)*alpha
      dqodw1(j) = -(alpha*dhkdw(j) - hk(j)*dsmpdw(j)/den)
      dqodw2(j) = -(alpha*dhkdw(j) + hk(j)*dsmpdw_aquifer/den)

      amx(j) = -dqidw0(j)
      bmx(j) =  dzmm(j)/deltim - dqidw1(j) + dqodw1(j)
      cmx(j) =  dqodw2(j)
      rmx(j) =  qin(j) - qout(j) - etr*rootr(j)

      ! next set up aquifer layer; alpha unchanged, qin=qout
      qin(j+1)    = qout(j)
      dqidw0(j+1) = dqodw1(j)
      dqidw1(j+1) = dqodw2(j)
      qout(j+1)   =  0._r8  ! zero-flow bottom boundary condition
      dqodw1(j+1) =  0._r8  ! zero-flow bottom boundary condition

      amx(j+1) = -dqidw0(j+1)
      bmx(j+1) =  dzmm(j+1)/deltim - dqidw1(j+1) + dqodw1(j+1)
      cmx(j+1) =  0._r8
      rmx(j+1) =  qin(j+1) - qout(j+1)

    endif

    ! Solve for dwat
    call tridia (nl_soil+1, amx, bmx, cmx, rmx, dwat_soi_acq)
    dwat(1:nl_soil) = dwat_soi_acq(1:nl_soil)

    ! Renew the mass of liquid water
    ! also compute qcharge from dwat in aquifer layer
    ! update the mass of liquid water for soil layers
    do j= 1, nl_soil
      wliq_soisno(j) = wliq_soisno(j)+dwat(j)*dzmm(j)
    enddo
    ! calculate qcharge for case jwt < nl_soil
    if(jwt < nl_soil) then
      wh_zwt = 0._r8   !since wh_zwt = psi0 - zq_zwt, where zq_zwt = psi0
      ! Recharge rate qcharge to groundwater (positive to aquifer)
      s_node = max(vol_liq(jwt+1)/porsl(jwt+1),0.01_r8)
      s1 = min(1._r8, s_node)
      !scs: this is the expression for unsaturated hk
      ka = imped(jwt+1)*hksati(jwt+1) &
           *s1**(2._r8*bsw(jwt+1)+3._r8)
      ! Recharge rate qcharge to groundwater (positive to aquifer)
      smp1 = max(smpmin, smp(max(1,jwt)))
      wh      = smp1 - zq(max(1,jwt))
      !scs: original formulation
      if(jwt == 0) then
         qcharge = -ka * (wh_zwt-wh)  /((zwt+1.e-3)*1000._r8)
      else
         !             qcharge = -ka * (wh_zwt-wh)/((zwt-z(jwt))*1000._r8)
         !scs: 1/2, assuming flux is at zwt interface, saturation deeper than zwt
         qcharge = -ka * (wh_zwt-wh)/((zwt-z_soisno(jwt))*1000._r8*2.0)
      endif
      ! To limit qcharge  (for the first several timesteps)
      qcharge = max(-10.0_r8/deltim,qcharge)
      qcharge = min( 10.0_r8/deltim,qcharge)
    else
      ! if water table is below soil column, compute qcharge from dwat_soi_acq(11)
      qcharge = dwat_soi_acq(nl_soil+1)*dzmm(nl_soil+1)/deltim
    endif

  end subroutine soilwater


  subroutine subsurfacerunoff(nl_soil,deltim,pondmx,hk,hksati,anik,&
                              wsat,wrsd,watern,alpha,slope,length,&
                              ss_f,Ds,Dr,Dg,kground,t_soisno,&
                              eff_porosity,icefrac,wfld,smf_soisno, &
                              dz_soisno,zi_soisno,wice_soisno,wliq_soisno,&
                              porsl,psi0,bsw,zwt,wa,Drw,&
                              qcharge,rsubst,qlat,qground,qsmelt,sst)
! -------------------------------------------------------------------------

    use precision
    use PhysicalConstants, only : tfrz
!
! ARGUMENTS:
    IMPLICIT NONE

    integer,  INTENT(in) :: nl_soil      ! number of soil layers [-]
    real(r8), INTENT(in) :: deltim       ! land model time step (sec)
    real(r8), INTENT(in) :: pondmx       !
    real(r8), INTENT(in) :: anik
    real(r8), INTENT(in) :: hk(1:nl_soil)! hydraulic conductivity [mm h2o/s]
    real(r8), INTENT(in) :: hksati (1:nl_soil) ! hydraulic conductivity at saturation (mm H2O /s)
    real(r8), INTENT(in) :: eff_porosity(1:nl_soil) ! effective porosity = porosity - vol_ice
    real(r8), INTENT(in) :: icefrac(1:nl_soil)      ! ice fraction (-)
    real(r8), INTENT(in) :: wfld(1:nl_soil)  ! water field capacity
    real(r8), INTENT(in) :: dz_soisno  (1:nl_soil)  ! layer depth (m)
    real(r8), INTENT(in) :: zi_soisno  (0:nl_soil)  ! interface level below a "z" level (m)
    real(r8), INTENT(in) :: porsl(1:nl_soil)        ! volumetric soil water at saturation (porosity)
    real(r8), INTENT(in) :: psi0(1:nl_soil)         ! minimum soil suction (mm) [-]
    real(r8), INTENT(in) :: bsw(1:nl_soil)          ! Clapp and Hornberger "b"
    real(r8), INTENT(in) :: qcharge   ! aquifer recharge rate (positive to aquifer) (mm/s)
    real(r8), INTENT(in) :: t_soisno (1:nl_soil) ! soil temperature (Kelvin)

    real(r8), INTENT(in) :: ss_f,wsat,wrsd,watern,alpha !GBHM parameter
    real(r8), INTENT(in) :: & ! added by HONGYILI, 2015.12.13
                Ds, &       ! depth of topsoil(m)
                Dr, &       ! depth of river (m)
                Dg, &       ! depth of unconfined acquifer (m)
                length, &   ! length of hillslope (m)
                slope,  &   ! slope of hillslope (m)
                kground     ! GW hydraulic conductivity (m/s)

    real(r8), INTENT(inout) :: smf_soisno(0:nl_soil+1)
    real(r8), INTENT(inout) :: wice_soisno(1:nl_soil)  ! ice lens (kg/m2)
    real(r8), INTENT(inout) :: wliq_soisno(1:nl_soil)  ! liquid water (kg/m2)
    real(r8), INTENT(inout) :: zwt       ! the depth from ground (soil) surface to water table [m]
    real(r8), INTENT(inout) :: wa        ! water in the unconfined aquifer (mm)
    real(r8), INTENT(inout) :: sst       ! surface liquid water storage (mm)
    real(r8), INTENT(inout) :: qsmelt
    real(r8), INTENT(inout) :: Drw ! water depth in the river of the flow interval (m)
                                   ! Different from old GBHM,Changed to in grid scale,2015-12-13
    real(r8), INTENT(out)   :: rsubst    ! drainage drainage (positive = out of soil column) (mm H2O /s)
    real(r8), INTENT(out)   :: qlat      ! lateral runoff(positive = out of soil column) (mm H2O /s)
    real(r8), INTENT(out)   :: qground   ! ground runoff(positive = out of soil column) (mm H2O /s)

!
! LOCAL ARGUMENTS
!

    integer  :: j                ! indices
    integer  :: jwt              ! index of the soil layer right above the water table (-)
    real(r8) :: xs               ! water needed to bring soil moisture to watmin (mm)
    real(r8) :: dzmm(1:nl_soil)  ! layer thickness (mm)
    real(r8) :: sm_soil(1:nl_soil) ! meltwater in soil layers(mm)
    real(r8) :: xsi              ! excess soil water above saturation at layer i (mm)
    real(r8) :: xsia             ! available pore space at layer i (mm)
    real(r8) :: xs1              ! excess soil water above saturation at layer 1 (mm)
    real(r8) :: ws               ! summation of pore space of layers below water table (mm)
    real(r8) :: s_node           ! soil wetness (-)
    real(r8) :: available_wliq_soisno     ! available soil liquid water in a layer
    real(r8) :: qcharge_tot      !
    real(r8) :: qcharge_layer    !
    real(r8) :: drainage         !
    real(r8) :: drainage_tot     !
    real(r8) :: drainage_layer   !
    real(r8) :: s_y              !
    real(r8) :: rous             ! specific yield [-]

    real(r8) :: wt
    real(r8) :: wtsub
    real(r8) :: dzsum
    real(r8) :: icefracsum
    real(r8) :: fracice_rsub
    real(r8) :: imped
    real(r8) :: tmpporsl
    real(r8) :: wa0,groundmelt

    real(r8), parameter :: watmin = 1E-16 ! Limit irreduciable wrapping liquid water,0.01
                                          ! a tunable constant
    real(r8), parameter :: rsbmx  = 5.0   ! baseflow coefficient [mm/s]
    real(r8), parameter :: timean = 10.5  ! global mean topographic index

    integer,  parameter :: use_clm_drainage = 0 !1-use;0-not use, instead GBHM
    integer :: isat
! -------------------------------------------------------------------------

!   ! Convert layer thicknesses from m to mm

    dzmm(:) = dz_soisno(:)*1000.
    ! update the jwt according to zwt
    call zwt_to_jwt(nl_soil,zwt,zi_soisno,jwt)

    sm_soil(1:nl_soil) = wliq_soisno(1:nl_soil)*smf_soisno(1:nl_soil)
!============================== QCHARGE =========================================
! Water table changes due to qcharge
! use analytical expression for aquifer specific yield

    rous = porsl(nl_soil)*(1.-(1.-1.e3*zwt/psi0(nl_soil))**(-1./bsw(nl_soil)))
    rous = max(rous,0.02)
    ! rous = 0.1


    ! wa = wa + qcharge*deltim
!
!---------------------------------------
    ! water table is below the soil column
    if(jwt == nl_soil) then
      wa  = wa + qcharge*deltim
      wa0 = wa      
      zwt = max(0.,zwt - (qcharge*deltim)/1000./rous)
    else
    ! water table within soil layers 1-9
    ! try to raise water table to account for qcharge

      qcharge_tot = qcharge * deltim

      if(qcharge_tot > 0.) then ! rising water table
          do j = jwt+1, 1,-1
             ! use analytical expression for specific yield

             s_y = porsl(j) * (1.-(1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
             s_y=max(s_y,0.02)
             qcharge_layer = min(qcharge_tot,(s_y*(zwt-zi_soisno(j-1))*1.e3))
             qcharge_layer = max(qcharge_layer,0.)

             zwt = max(0.,zwt - qcharge_layer/s_y/1000.)

             qcharge_tot = qcharge_tot - qcharge_layer
             if (qcharge_tot <= 0.) exit
          enddo
      else ! deepening water table (negative qcharge)

          do j = jwt+1, nl_soil
             ! use analytical expression for specific yield
             s_y = porsl(j) * (1.-(1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
             s_y=max(s_y,0.02)
             qcharge_layer = max(qcharge_tot,-(s_y*(zi_soisno(j) - zwt)*1.e3))
             qcharge_layer = min(qcharge_layer,0.)
             qcharge_tot = qcharge_tot - qcharge_layer

             if (qcharge_tot >= 0.) then
                zwt = max(0.,zwt - qcharge_layer/s_y/1000.)
                exit
             else
                zwt = zi_soisno(j)
             endif
          enddo
          if (qcharge_tot > 0.) zwt = max(0.,zwt - qcharge_tot/1000./rous)
      endif
      ! update the jwt according to zwt
      call zwt_to_jwt(nl_soil,zwt,zi_soisno,jwt)
    endif
    wa0 = wa

    ! update the jwt according to zwt

!============================== QGROUND[GBHM]====================================
    CALL ground_runoff(deltim,nl_soil,Ds,Dr,Dg,kground,length,slope,anik,& !IN
                            dz_soisno,wice_soisno(1:),hksati,wsat,&        !IN
                            porsl,psi0,bsw, &
                            wliq_soisno,zwt,wa,Drw,&                       !IN-OUT
                            qground,isat)                                  !OUT

!=======================================================================
! [GBHM] determine the lateral flow from soil layer.(For GBHM)
! Added by Hongyi,2015.11.14
!=======================================================================
    CALL lateral_runoff(hk,anik,nl_soil,jwt,dz_soisno,eff_porosity,deltim,&
                        ss_f,wfld,wsat,wrsd,watern,alpha,slope,length, &
                        t_soisno, wliq_soisno,qlat,smf_soisno(1:),qsmelt)


! ======================================================================
! [Default] not use the following code in origin CLM, instead from GBHM
! Modified by Hongyi Li, 2015-12-20
! ======================================================================
    if (use_clm_drainage==1) then
    !-- Topographic runoff  --------------------------------------------
        dzsum = 0.
        icefracsum = 0.
        do j = max(jwt,1), nl_soil
        dzsum = dzsum + dzmm(j)
        icefracsum = icefracsum + icefrac(j) * dzmm(j)
        end do
        ! add ice impedance factor to baseflow
        fracice_rsub = max(0.,exp(-3.*(1.-(icefracsum/dzsum))) &
                        -exp(-3.))/(1.0-exp(-3.))
        imped = max(0.,1.-fracice_rsub)
        ! drainage (positive = out of soil column)
        drainage = imped * 5.5e-3 * exp(-2.5*zwt)

    !-- Water table is below the soil column  --------------------------
        if(jwt == nl_soil) then
        wa = wa - drainage * deltim
        zwt = max(0.,zwt + (drainage * deltim)/1000./rous)
        wliq_soisno(nl_soil) = wliq_soisno(nl_soil) + max(0.,(wa-5000.))
        wa = min(wa, 5000.)
        else
    !-- Water table within soil layers 1-9  ----------------------
    !============================== RSUB_TOP =====================
        !-- Now remove water via drainage
        drainage_tot = - drainage * deltim
        do j = jwt+1, nl_soil
            ! use analytical expression for specific yield
            s_y = porsl(j) * ( 1. - (1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
            s_y = max(s_y,0.02)

            drainage_layer = max(drainage_tot, -(s_y*(zi_soisno(j)-zwt)*1.e3))
            drainage_layer = min(drainage_layer,0.)
            wliq_soisno(j) = wliq_soisno(j) + drainage_layer

            drainage_tot = drainage_tot - drainage_layer

            if(drainage_tot >= 0.)then
                zwt = max(0.,zwt - drainage_layer/s_y/1000.)
                exit
            else
                zwt = zi_soisno(j)
            endif
        enddo

    !-- Remove residual drainage  ---
        zwt = max(0.,zwt - drainage_tot/1000./rous)
        wa = wa + drainage_tot

    !-- Recompute jwt  ----
        call zwt_to_jwt(nl_soil,zwt,zi_soisno,jwt)

        end if   ! end of jwt if construct
    else
        drainage = 0.
    !-- Recompute jwt  ----
        call zwt_to_jwt(nl_soil,zwt,zi_soisno,jwt)
    end if

! ---------------------------------------------------------------
    zwt = max(0.0,zwt)
    zwt = min(80.,zwt)

    rsubst = drainage + qlat + qground  !changed by Hongyi, '+ qlat', 20151114
                                        !changed by Hongyi, '+ qground', 20151213

    ! groundmelt = qground*smf_soisno(nl_soil+1) !assume ground melt is homogeneous flow out
    groundmelt = min(qground,wa0*smf_soisno(nl_soil+1)/deltim)
    groundmelt = max(0.,groundmelt)
    qsmelt = qsmelt + groundmelt

    if (wa>0.) then
      smf_soisno(nl_soil+1) = (wa0*smf_soisno(nl_soil+1)-groundmelt*deltim)/wa
      smf_soisno(nl_soil+1) = min(1.,smf_soisno(nl_soil+1))
      smf_soisno(nl_soil+1) = max(0.,smf_soisno(nl_soil+1))
    else
      smf_soisno(nl_soil+1) = 0.
    endif

    ! Correction [1]
    ! NON-physically based corection on wliq_soisno
    ! excessive water above saturation added to the above unsaturated layer like a bucket
    ! if column over saturated, excess water goes to runoff

    do j = nl_soil,2,-1
       tmpporsl = max(0.,(porsl(j)*dzmm(j)-wice_soisno(j)))
       xsi = max(wliq_soisno(j)-tmpporsl,0.)
       if (xsi>0.) then
        wliq_soisno(j)   = min(tmpporsl, wliq_soisno(j))
        sm_soil(j)       = max(0.,(wliq_soisno(j)-xsi)*smf_soisno(j))
        sm_soil(j-1)     = smf_soisno(j-1)*wliq_soisno(j-1) + smf_soisno(j)*xsi
        wliq_soisno(j-1) = wliq_soisno(j-1) + xsi
        smf_soisno(j-1)  = min(1.,sm_soil(j-1)/max(sm_soil(j-1),wliq_soisno(j-1)))
       endif
    end do

    !--Major revised by HY. 2016-11-14 AM
    !
    !----1. Delete pondmx, because 'pondmx' has been addressed in GBHM module.
    !    xs1 = wliq_soisno(1) - pondmx+(porsl(1)*dzmm(1)-wice_soisno(1))
    !
    !----2. Because snow sublimation/condensation are counted to the soil surface mass in
    !    WATER subroutine:
    !    --- wice_soisno(1) = max(0., wice_soisno(1) + (qfros-qsubl) * deltim)
    !    So there are some water balance bugs ocuured in case ice mass is more than PORSL
    !    can hold:
    !    --- (porsl(1)*dzmm(1)-wice_soisno(1))<0.
    !    --- xs1 = wliq_soisno(1) - (porsl(1)*dzmm(1)-wice_soisno(1)) > 0.
    !    In this case, xs1 can make some runoff in qlat and cause balance error without
    !    liquid water in surface soil.
    !    I set a discriminant to limit the false liquid water increase:
    !    --- tmpporsl = max(0.,(porsl(1)*dzmm(1)-wice_soisno(1)))
    ! ------------------------------------------------------------------------------------

    tmpporsl = max(0.,(porsl(1)*dzmm(1)-wice_soisno(1)))
    xs1 = wliq_soisno(1) - tmpporsl
    if(xs1 > 0.)then
       ! wliq_soisno(1) = pondmx+porsl(1)*dzmm(1)-wice_soisno(1)!
       wliq_soisno(1) = tmpporsl
       ! qlat   = qlat   + xs1 / deltim
       ! rsubst = rsubst + xs1 / deltim
       ! qsmelt = qsmelt + xs1*smf_soisno(1)/deltim
       smf_soisno(0) = min(1.,(sst*smf_soisno(0)+xs1*smf_soisno(1))/(sst+xs1))
       sst = sst + xs1
    else
       xs1 = 0.
    endif

    ! Correction [2]
    ! NON-physically based corection on wliq_soisno
    ! Limit wliq_soisno to be greater than or equal to watmin.
    ! Get water needed to bring wliq_soisno equal watmin from lower layer.
    ! If insufficient water in soil layers, get from aquifer water
    ! Rewrite it by HY, 2016-4-26.

    xs = 0.
    if (wliq_soisno(1) < 0.) then
      ! print*,'debug,wliq_soisno(1) < 0.'
      if ( wliq_soisno(2)>watmin-wliq_soisno(1)) then
        wliq_soisno(2) = wliq_soisno(2) + wliq_soisno(1) - watmin
      else
        xs = xs + wliq_soisno(1) - watmin
      endif
      wliq_soisno(1) = watmin
    endif

    do j = 2, nl_soil-1
       if (wliq_soisno(j) < 0.) then
        ! print*,'debug,wliq_soisno(j) < 0.',j
          if ( wliq_soisno(j+1)>watmin-wliq_soisno(j)) then
            wliq_soisno(j+1) = wliq_soisno(j+1) + wliq_soisno(j) - watmin
          elseif (wliq_soisno(j-1)>watmin-wliq_soisno(j)) then
            wliq_soisno(j-1) = wliq_soisno(j-1) + wliq_soisno(j) - watmin
          else
            xs = xs + wliq_soisno(j) - watmin
          endif
          wliq_soisno(j) = watmin
       endif
    enddo

    if (wliq_soisno(nl_soil) < 0.) then
      ! print*,'debug,wliq_soisno(n) < 0.'
      if ( wliq_soisno(nl_soil-1)>watmin-wliq_soisno(nl_soil)) then
        wliq_soisno(nl_soil-1) = wliq_soisno(nl_soil-1) + wliq_soisno(nl_soil) - watmin
      else
        xs = xs + wliq_soisno(nl_soil) - watmin
      endif
      wliq_soisno(nl_soil) = watmin
    endif

    ! Adjust sub-surface runoff
    if (xs<-1.0e-10) then
      ! 'xs' is the water close error. Assume it should be balanced by the neighbour layer.
      available_wliq_soisno = 0.
      do j = 1, nl_soil
        if (wliq_soisno(j)>watmin) then
          available_wliq_soisno = available_wliq_soisno + wliq_soisno(j)
        endif
      enddo
      if (available_wliq_soisno>abs(xs)) then
        do j = 1, nl_soil
          if (wliq_soisno(j)>watmin) then
            wliq_soisno(j) = wliq_soisno(j) - (abs(xs)*wliq_soisno(j)/available_wliq_soisno)
          endif
        enddo
      elseif (qground + xs/deltim >0.) then
        rsubst = rsubst + xs/deltim
        qground   = qground   + xs / deltim
      elseif (qlat + xs/deltim >0.) then
        rsubst = rsubst + xs/deltim
        qlat   = qlat   + xs / deltim
      elseif (wa+xs >0.) then
        wa = wa+xs
        ! print*,'debug wa, negative soil mositure'
      else
        ! print*,'debug wa, negative soil mositure,no solved xs:',xs
      endif
    endif

  end subroutine subsurfacerunoff


  subroutine lateral_runoff(ksoil,anik,nl_soil,isat,deltz,eff_porosity,dt,ss_f,&
                            wfld,wsat,wrsd,watern,alpha,slope,length, t_soisno,&
                            wliq_soisno,qlat,smf_soisno,qsmelt)
    ! Added by Hongyi, for coupling GBHM
    !------------------------------------------------------------------------
    ! Subsurface flow from top saturated zone above the groundwater level
    ! (GBHM)
    ! wsat,wrsd,watern,alpha,slope. These 5 parameters is from local file.
    ! ksoil is the hydraulic conductivity of soil, which is from GBHM scheme
    ! or CLM scheme.
    !------------------------------------------------------------------------
    use precision
    use PhysicalConstants, only : tfrz
    implicit none
    integer,  INTENT(in)    :: nl_soil,isat
    real(r8), INTENT(in)    :: wsat,wrsd,watern,alpha,slope,length,dt
    real(r8), INTENT(in)    :: anik,ss_f
    real(r8), INTENT(in)    :: ksoil(1:nl_soil)       ! hydraulic conductivity of each layer (mm/s)
    real(r8), INTENT(in)    :: deltz(1:nl_soil)       ! m
    real(r8), INTENT(in)    :: wfld(1:nl_soil)        ! water field capacity
    real(r8), INTENT(in)    :: eff_porosity(1:nl_soil)! effective porosity = porosity - vol_ice
    real(r8), INTENT(in)    :: smf_soisno(1:nl_soil+1)
    real(r8), INTENT(in)    :: t_soisno(1:nl_soil)    ! soil temperature (Kelvin)
    real(r8), INTENT(inout) :: wliq_soisno(1:nl_soil) ! liquid water (kg/m2)
    real(r8), INTENT(inout) :: qsmelt
    real(r8), INTENT(out)   :: qlat
    ! real(r8)  wfld                             ! soil moisture at field capacity
    real(r8)  suction, w(1:nl_soil),mksoil(1:nl_soil)
    real(r8)  tmp,tmpw,tmpqsub,f,tempwfld(1:nl_soil)
    integer   i,j

    ! transfer hydraulic conductivity and soil water depth to GBHM format
    !--------------------------------------------------------------------
    mksoil(1:nl_soil) = ksoil(1:nl_soil)/1000. !mm/s -> m/s
    w(1:nl_soil)      = wliq_soisno(1:nl_soil)/1000./deltz(1:nl_soil)!kg/m2 -> %
    !--------------------------------------------------------------------

    ! suction   = -1.02  ! suction, meter water-head
    ! call MoistureFromSuction_V(wfld,wsat,wrsd,watern,alpha,suction)!move this routine out, just keep wfld, and control the value range of wfld.
    qlat= 0.
    tmp = 0.
    j = min(isat, nl_soil)
    ! j = nl_soil
    do i = 1, j
      tmp = tmp + deltz(i)
      if(((w(i)-wfld(i)) .gt. 0.1E-3) .and. t_soisno(i) .gt.tfrz) then
        ! The following codes is from Wanglei's WEB-DHM model.
        f = -log(1.0/anik)/5.0
        tmpqsub = mksoil(i)*slope*dt
        tmpqsub = tmpqsub*anik *exp(-f*tmp)* (deltz(i)+ss_f*length)/length
        if(tmpqsub .lt. 0.1E-20) tmpqsub = 0.0
        tmpw    = (w(i)-wfld(i))*deltz(i)
        tmpqsub = amin1(tmpqsub, tmpw)
        tmpqsub = amax1(tmpqsub, 0.0)
        w(i)    = w(i) - tmpqsub/deltz(i)
        if (w(i).lt.0.) print*,'debug: negative soil moisture in line 1567'
        wliq_soisno(i)=w(i)*deltz(i)*1000.! m -> kg/m2
        qlat = qlat + tmpqsub

        qsmelt = qsmelt+(tmpqsub* 1000./dt)*smf_soisno(i) !Added in 2016-11-28, HYLI

      endif
    enddo
    qlat = qlat * 1000. / dt              ! m -> mm/sec
  end subroutine lateral_runoff

  subroutine smfsoil(nl_soil,deltim,qinfl,etr,dzmm,rootr,dwat,qcharge,wliq_soisno,wa,smf_soisno)
    use precision
    IMPLICIT NONE
    integer,  INTENT(in) :: nl_soil  ! number of soil layers
    real(r8), INTENT(in) :: deltim  ! land model time step (sec)
    real(r8), INTENT(in) :: qinfl   ! infiltration (mm H2O /s)
    real(r8), INTENT(in) :: etr     ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), INTENT(in) :: dzmm(1:nl_soil) ! layer thickness (mm)
    real(r8), INTENT(in) :: rootr    (1:nl_soil) ! effective fraction of roots in each soil layer
    real(r8), intent(in) :: dwat     (1:nl_soil) ! change of soil water [m3/m3]
    real(r8), INTENT(in) :: qcharge              ! aquifer recharge rate (positive to aquifer) (mm/s)
    real(r8), INTENT(in) :: wliq_soisno(1:nl_soil)
    real(r8), INTENT(in) :: wa

    real(r8), INTENT(inout) :: smf_soisno(0:nl_soil+1)

    integer j
    real(r8) qintmp(1:nl_soil+1), qin(1:nl_soil+1)
    real(r8) qinup(1:nl_soil+1),qinsm_down(1:nl_soil+1),qinsm_up(1:nl_soil+1)
    real(r8) oldliq(1:nl_soil)
    real(r8) debugtmp

    qin(:) = 0.
    qin(1) = qinfl * deltim
    oldliq(:) = 0.
    oldliq(1) = wliq_soisno(1) - dwat(1)*dzmm(1) + etr*rootr(1)*deltim
    oldliq(1) = max(0.,oldliq(1))
    ! qin(nl_soil+1) = qcharge*deltim ! if this is equal to qin (-1) ?

    qintmp(:) = qin(:)
    qinsm_down(1) = qin(1)*smf_soisno(0)
    qinup(:) = 0.
    if (qin(1)<0.) then
      qintmp(1) = 0.
      qinsm_down(1) = 0.
      qinsm_up(1)   = -qin(1)*smf_soisno(1)
      qinup = -1.*qin(1)
    endif

    do j= 2, nl_soil+1
      qin(j) = qin(j-1) - dwat(j-1)*dzmm(j-1) - etr*rootr(j-1)*deltim
      if (j<=nl_soil) then
        oldliq(j) = wliq_soisno(j) - dwat(j)*dzmm(j) + etr*rootr(j)*deltim
        oldliq(j) = max(0.,oldliq(j))
      endif
      if (qin(j)<0.) then
        qintmp(j) = 0.
        qinup(j)  = -1.*qin(j)
      else
        qintmp(j) = qin(j)
        qinup(j)  = 0.
      endif
      if (qin(j-1)>0.)then
        if ((oldliq(j-1)+qintmp(j-1))>0.1E-8) then
          smf_soisno(j-1) = (smf_soisno(j-1)*oldliq(j-1) + &
                            qinsm_down(j-1))/(oldliq(j-1)+qintmp(j-1))
          smf_soisno(j-1) = min(smf_soisno(j-1),1.)
          oldliq(j-1) = oldliq(j-1)+qintmp(j-1)
        else
          smf_soisno(j-1) = 0.
          oldliq(j-1) = max(0.,oldliq(j-1)+qintmp(j-1))
        endif
      endif
      qinsm_down(j) = smf_soisno(j-1)*qintmp(j)

      if ((j==nl_soil+1).and.(qin(j)>0.)) then
        if ((wa+qintmp(j))>0.1E-6)then
          smf_soisno(j) = min(1.,(wa*smf_soisno(j)+qinsm_down(j))/(wa+qintmp(j)))
        else
          smf_soisno(j) = min(1.,(wa*smf_soisno(j)+qinsm_down(j))*1.E10/((wa+qintmp(j))*1.E10))
        endif
      endif
    enddo

    do j= nl_soil+1,2, -1
      if(qin(j)<0.)then
        qinsm_up(j) = smf_soisno(j)*qinup(j)
        if ((oldliq(j-1)+qinup(j))>0.1E-8) then
          smf_soisno(j-1) = (smf_soisno(j-1)*oldliq(j-1) &
                            + qinsm_up(j))/(oldliq(j-1)+qinup(j))
          smf_soisno(j-1) = min(smf_soisno(j-1),1.)
        else
          smf_soisno(j-1) = 0.
        endif
      endif
    enddo
  end subroutine smfsoil

  subroutine smsoil(nl_soil,deltim,qinfl,etr,dzmm,rootr,dwat,qcharge,wliq_soisno,wa,smf_soisno,evpsm)
    use precision
    IMPLICIT NONE
    integer,  INTENT(in) :: nl_soil ! number of soil layers
    real(r8), INTENT(in) :: deltim  ! land model time step (sec)
    real(r8), INTENT(in) :: qinfl   ! infiltration (mm H2O /s)
    real(r8), INTENT(in) :: etr     ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), INTENT(in) :: dzmm(1:nl_soil) ! layer thickness (mm)
    real(r8), INTENT(in) :: rootr    (1:nl_soil) ! effective fraction of roots in each soil layer
    real(r8), intent(in) :: dwat     (1:nl_soil) ! change of soil water [m3/m3]
    real(r8), INTENT(in) :: qcharge              ! aquifer recharge rate (positive to aquifer) (mm/s)
    real(r8), INTENT(in) :: wliq_soisno(1:nl_soil)
    real(r8), INTENT(in) :: wa

    real(r8), INTENT(inout) :: smf_soisno(0:nl_soil+1)
    real(r8), INTENT(out) :: evpsm
    real(r8)  sm_soisno(1:nl_soil+1)

    integer  j
    integer  i,ni
    real(r8) qindown(1:nl_soil+1), qin(1:nl_soil+1)
    real(r8) qinup(1:nl_soil+1)
    real(r8) oldliq(1:nl_soil+1)
    real(r8) tmpliq(1:nl_soil+1) ! intermediate liquid layer (mm)
    real(r8) tmpliq_t(1:nl_soil+1) ! intermediate liquid layer (mm)
    real(r8) tmpsmelt(1:nl_soil+1) ! intermediate snowmelt liquid layer (mm)
    real(r8) tmpsmelt_t(1:nl_soil+1) ! intermediate snowmelt liquid layer (mm)
    real(r8) sumliq_j,residual_liq,upward_wliq_tmp(1:nl_soil+1)
    real(r8) wliq_tmp(1:nl_soil+1),upward_wliq(1:nl_soil+1),upward_sm(1:nl_soil+1)
    real(r8) soilsm_old,soilsm_new,err_soilsm,insmelt,outsmelt  ! for checking the balance of snow melt water in soil
    real(r8) adjust_wa,tmpetr
    real(r8) debugtmp1,debugtmp2,debugtmp3,debugtmp4

    tmpetr = etr
    ! tmpetr = max(0.,etr)
    ! if(etr>0. .and. sum(rootr(:))<=0.) tmpetr=0.
    qin(:) = 0.
    qin(1) = qinfl * deltim
    oldliq(:) = 0.
    oldliq(1) = wliq_soisno(1) - dwat(1)*dzmm(1) !+ tmpetr*rootr(1)*deltim
    oldliq(1) = max(0.,oldliq(1))
    oldliq(nl_soil+1) = wa
    adjust_wa = wa
    ! qin(nl_soil+1) = qcharge*deltim ! if this is equal to qin (-1) ?

    qindown(:) = qin(:)
    qinup(:) = 0.
    if (qin(1)<0.) then
      qindown(1) = 0.
      qinup(1) = -1.*qin(1)
    endif

    do j= 2, nl_soil+1
      qin(j) = qin(j-1) - dwat(j-1)*dzmm(j-1) - tmpetr*rootr(j-1)*deltim
      if (j<=nl_soil) then
        oldliq(j) = wliq_soisno(j) - dwat(j)*dzmm(j) !+ tmpetr*rootr(j)*deltim
if (oldliq(j)<0.) print*,'debug oldliq',oldliq
        oldliq(j) = max(0.,oldliq(j))
      endif
      if (qin(j)<0.) then
        qindown(j) = 0.
        qinup(j)  = -1.*qin(j)
      else
        qindown(j) = qin(j)
        qinup(j)  = 0.
      endif

      if(oldliq(j)<=0.) smf_soisno(j) = 0.

      if ((j==nl_soil+1)) then
        ! if(abs(qin(j)-qcharge*deltim)>0.1E-6) then
        !   print*,'error in sm_soisno groundwater module ',qin(j),qcharge*deltim,etr
        ! endif
        if(qinup(j)>wa) then
          ! print*,'error in qinup(j)>wa ',qinup(j),qcharge*deltim,wa,smf_soisno(j),sm_soisno(j)
          adjust_wa = qinup(j)
          oldliq(nl_soil+1) = adjust_wa
          ! smf_soisno(j) = 0.
        endif
      endif

    enddo

if(abs(sum(oldliq(1:nl_soil))-qin(nl_soil+1)+qin(1)-tmpetr*deltim-sum(wliq_soisno(:)))>1e-8) then
  print*,'debug balance',sum(oldliq(1:nl_soil)),qin(nl_soil+1),qin(1),tmpetr*deltim,sum(wliq_soisno(:)),&
          sum(oldliq(:))-wa+qin(1)-tmpetr*deltim-sum(wliq_soisno(:))
  ! pause
endif
    wliq_tmp(:) = oldliq(:)
    sm_soisno(1:nl_soil+1) = smf_soisno(1:nl_soil+1)*oldliq(1:nl_soil+1)
    soilsm_old = sum(sm_soisno(1:nl_soil+1))
    insmelt = qindown(1)*smf_soisno(0)
    do j= 1, nl_soil
      if(j==1)then
        ni = 1
        tmpliq_t(:)   = 0.
        tmpsmelt_t(:) = 0.
        tmpliq(j)     = oldliq(j)
        tmpsmelt(j)   = sm_soisno(j)
        if (qindown(1)>0.) then
          tmpliq(j+1)   = qindown(1)
          tmpsmelt(j+1) = qindown(1)*smf_soisno(0)
          ni = 2
        endif
      else
        tmpliq(:)   = tmpliq_t(:)
        tmpsmelt(:) = tmpsmelt_t(:)
      endif

      if(qindown(j+1)<=0.)then
        wliq_tmp(j)   = sum(tmpliq(1:ni))
        sm_soisno(j)  = sum(tmpsmelt(1:ni))
        tmpliq_t(1)   = oldliq(j+1)
        tmpsmelt_t(1) = sm_soisno(j+1)
        ni = 1
      else
        sumliq_j = 0.
        do i=1,ni
          sumliq_j  = sumliq_j+tmpliq(i)
          if (sumliq_j>=qindown(j+1)) then
            residual_liq    = sumliq_j-qindown(j+1)
            ! update snowmelt proporation in j layer
            tmpliq_t(i+1)   = tmpliq(i)-residual_liq
            tmpliq_t(1)     = oldliq(j+1)
            if(tmpliq(i)>1E-15)then
              tmpsmelt_t(i+1) = (min(1.,(tmpliq_t(i+1)/tmpliq(i))))*tmpsmelt(i)
            else
              tmpsmelt_t(i+1) = 0.
            endif
            ! tmpsmelt_t(i+1) = tmpsmelt(i)-(tmpliq_t(i+1)/tmpliq(i))*tmpsmelt(i)
            tmpsmelt_t(1)   = sm_soisno(j+1)
            if(tmpliq(i)>1E-15)then
              sm_soisno(j)  = (min(1.,residual_liq/tmpliq(i)))*tmpsmelt(i)+sum(tmpsmelt(i+1:ni))
            else
              sm_soisno(j)  = sum(tmpsmelt(i+1:ni))
            endif
            wliq_tmp(j)     = residual_liq + sum(tmpliq(i+1:ni))
            ni = i+1
            exit
          endif
          tmpliq_t(i+1)   = tmpliq(i)
          tmpsmelt_t(i+1) = tmpsmelt(i)
        enddo
      endif

      if ((j==nl_soil).and.(qin(j+1)>0.)) then
        sm_soisno(j+1) = sum(tmpsmelt_t(1:ni))
        wliq_tmp(j+1)  = sum(tmpliq_t(1:ni))
if(sm_soisno(j+1)>wliq_tmp(j+1)) sm_soisno(j+1)=wliq_tmp(j+1)
      endif
    enddo
! debugtmp2 = sum(sm_soisno(:))
! if (abs(debugtmp2-debugtmp1)>0.1E-10) print*,'debug temp 1-2',abs(debugtmp2-debugtmp1)
    !----calculate snowmelt distribution at upward direction
! print*,'debug1',oldliq,wliq_tmp
! debugtmp3 = sum(sm_soisno(:))
    upward_wliq(1:nl_soil) = wliq_tmp(nl_soil:1:-1)
    upward_sm(1:nl_soil)   = sm_soisno(nl_soil:1:-1)
    upward_wliq(nl_soil+1) = 0.
    upward_sm(nl_soil+1)   = 0.
    qinup(1:nl_soil+1)     = qinup(nl_soil+1:1:-1)
    do j= 1, nl_soil
      if(j==1)then
        ni = 1
        tmpliq_t(:)   = 0.
        tmpsmelt_t(:) = 0.
        tmpliq(j)     = upward_wliq(j)
        tmpsmelt(j)   = upward_sm(j)
        if (qinup(1)>0.) then
          tmpliq(j+1)   = qinup(1)
          tmpsmelt(j+1) = qinup(1)*smf_soisno(nl_soil+1)
          sm_soisno(nl_soil+1) = sm_soisno(nl_soil+1)-tmpsmelt(j+1)
          wliq_tmp(nl_soil+1)  = adjust_wa-qinup(1)
          ! if(wa-qinup(1)<0.) print*,'debug waqinup',wa,qinup(1)
          ni = 2
        endif
      else
        tmpliq(:)   = tmpliq_t(:)
        tmpsmelt(:) = tmpsmelt_t(:)
      endif
      if(qinup(j+1)<=0.)then
        upward_wliq_tmp(j) = sum(tmpliq(1:ni))
        upward_sm(j)  = sum(tmpsmelt(1:ni))
        tmpliq_t(1)   = upward_wliq(j+1)
        tmpsmelt_t(1) = upward_sm(j+1)
        ni = 1

      else
        sumliq_j = 0.
        do i=1,ni
          sumliq_j  = sumliq_j+tmpliq(i)
! print*,'sumliq_j',sumliq_j,tmpliq(i),qinup(j+1)
          residual_liq = sumliq_j-qinup(j+1)
          if (abs(residual_liq)<1E-15) then
            residual_liq = 0.
            qinup(j+1) = sumliq_j
          endif
          if (residual_liq>=0._r8) then
            ! update snowmelt proporation in j layer
            tmpliq_t(i+1)   = tmpliq(i)-residual_liq
            tmpliq_t(1)     = upward_wliq(j+1)
            if(tmpliq(i)>1E-15)then
              tmpsmelt_t(i+1) = (min(1.,(tmpliq_t(i+1)/tmpliq(i))))*tmpsmelt(i)
            else
              tmpsmelt_t(i+1) = 0.
            endif
            ! tmpsmelt_t(i+1) = tmpsmelt(i)-(tmpliq_t(i+1)/tmpliq(i))*tmpsmelt(i)
            tmpsmelt_t(1)   = upward_sm(j+1)
            if(tmpliq(i)>1E-15)then
              ! upward_sm(j)    = (residual_liq/tmpliq(i))*tmpsmelt(i)+sum(tmpsmelt(i+1:ni))
              upward_sm(j)  = (min(1.,residual_liq/tmpliq(i)))*tmpsmelt(i)+sum(tmpsmelt(i+1:ni))
            else
              upward_sm(j)  = sum(tmpsmelt(i+1:ni))
            endif
            upward_wliq_tmp(j)= residual_liq + sum(tmpliq(i+1:ni))
! print*,'inerror',j,i,ni,residual_liq,upward_sm(j) ,tmpliq(i),tmpsmelt(i) ,tmpsmelt(i+1:ni)
            ni = i+1
            exit
          endif
          tmpliq_t(i+1)   = tmpliq(i)
          tmpsmelt_t(i+1) = tmpsmelt(i)
        enddo
      endif

      if ((j==nl_soil).and.(qinup(j+1)>0.)) then
        evpsm = sum(tmpsmelt_t(1:ni))
      else
        evpsm =0.
      endif
    enddo
    sm_soisno(1:nl_soil) = upward_sm(nl_soil:1:-1)
    upward_wliq_tmp(1:nl_soil) = upward_wliq_tmp(nl_soil:1:-1)
    upward_wliq_tmp(nl_soil+1) = wliq_tmp(nl_soil+1)
! print*,tmpsmelt_t(1:ni),'tmpsmelt'
! print*,evpsm ,'evpsm'
! print*,'sm_soisno',sm_soisno
! debugtmp4 = sum(sm_soisno(:))+evpsm
! if (abs(debugtmp3-debugtmp4)>0.1E-10) print*,'debug temp 3-4',abs(debugtmp3-debugtmp4),evpsm,sm_soisno
! if (abs(debugtmp3-debugtmp4)>0.1E-10) print*,'debug temp 3-4',upward_sm,tmpsmelt_t,ni
! if (abs(debugtmp3-debugtmp4)>0.1E-10) print*,'debug temp 3-4',i,sumliq_j,qinup(8),residual_liq,tmpliq(1),tmpsmelt_t(1) ,upward_sm(8)

    do i=1,nl_soil+1
      if(sm_soisno(i)>upward_wliq_tmp(i))then
        ! print*,'sm_soisno(i)>upward_wliq_tmp(i)',i,sm_soisno(i),upward_wliq_tmp(i),smf_soisno(i)
      endif
      if (upward_wliq_tmp(i)>0.)then
        smf_soisno(i) = min(sm_soisno(i)/upward_wliq_tmp(i),1.)
      else
        smf_soisno(i) = 0.
      endif
    enddo
    ! soilsm_new = sum(sm_soisno(1:nl_soil+1))
    soilsm_new = sum(smf_soisno(1:nl_soil+1)*upward_wliq_tmp(1:nl_soil+1))
    outsmelt   = evpsm
    err_soilsm = abs(soilsm_old+insmelt-outsmelt-soilsm_new)
    ! if (err_soilsm>1E-8) then
    !   print*,'error in smsoil module',err_soilsm,insmelt,evpsm,qin(1),soilsm_old,soilsm_new
    !   print*,'error in smsoil module','wrong',sm_soisno,upward_wliq_tmp
    ! endif
  end subroutine smsoil

  subroutine ground_runoff(dt,nl_soil,Ds,Dr,Dg,kground,length,slope,anik,&
                            dz_soisno,wice_soisno,hksati,wsat,&
                            porsl,psi0,bsw,&
                            wliq_soisno,zwt,wa,Drw,&
                            qsub,isat)
    !----------------------------------------------------------------------
    ! calculation of exchange rate between groundwater and river: qsub
    ! Yang D. et al., Hydrological Processes, 14, 403-416,2000
    !----------------------------------------------------------------------
    use precision
    use PhysicalConstants, only: denice
    implicit none
    integer,intent(in) :: &
                nl_soil     ! Number of soil layers
    real(r8),intent(in):: &
                dt, &       ! time step
                Ds, &       ! depth of topsoil
                Dr, &       ! depth of river (m)
                Dg, &       ! depth of unconfined acquifer (m)
                length, &   ! length of hillslope (m)
                slope,  &   ! slope of hillslope (m)
                kground,&   ! GW hydraulic conductivity (m/s),from prefile.
                anik,   &   ! Soil anitropic ratio (>=1.0)
                dz_soisno(1:nl_soil),  & ! layer depth (m)
                wice_soisno(1:nl_soil),&
                hksati   (1:nl_soil),  & !In CLM: hksati(1:nl_soil),hydraulic conductivity at saturation (mm h2o/s)
                wsat

    real(r8), INTENT(in) :: porsl(1:nl_soil)        ! volumetric soil water at saturation (porosity)
    real(r8), INTENT(in) :: psi0 (1:nl_soil)        ! minimum soil suction (mm) [-]
    real(r8), INTENT(in) :: bsw  (1:nl_soil)        ! Clapp and Hornberger "b"

    real(r8),intent(inout):: &
                wliq_soisno(1:nl_soil), &
                zwt , &     ! the depth from ground (soil) surface to water table [m]
                wa  , &     ! water in the unconfined aquifer (mm)
                Drw         ! water depth in the river of the flow interval (m)
    real(r8), intent(out):: &
                qsub
    integer, intent(out)::  isat
    real(r8)::  k0   (1:nl_soil),  &! saturated hydraulic conductivity of each layer (m/s)
                deltz(1:nl_soil),  &! depth of each layer (m)
                wice (1:nl_soil),  &! Volumetric ice content, %
                w    (1:nl_soil)
    real(r8)::  tmp,tmpqsub,tmpw,avkg,tmp_D,excess,tmpkg
    real(r8)::  GWcs,&   ! GW storage coefficient
                Dgl ,&   ! depth to GW level (m)
                GWst     ! Groundwater storage (m)
    integer ::  i,nuz,GWkey

    !k0 and hksati is same, but in gbhm the unit of k0 is m/s, while hksati is mm/s in CLM.
    deltz(1:nl_soil) = dz_soisno(1:nl_soil)
    w(1:nl_soil)     = wliq_soisno(1:nl_soil)/1000./deltz(1:nl_soil)!kg/m2 -> %
    wice (1:nl_soil) = wice_soisno(1:nl_soil)/denice/deltz(1:nl_soil)!kg/m2 -> %
    k0(1:nl_soil)    = hksati(1:nl_soil)/1000.
    tmpkg            = kground
    ! tmpkg            = k0(nl_soil)!*0.01
! print*, kground, k0(nl_soil)
    Dgl     = zwt
    ! GWcs    = porsl(nl_soil)*(1.-(1.-1.e3*zwt/psi0(nl_soil))**(-1./bsw(nl_soil)))
    ! GWcs    = max(GWcs,0.02)
    ! GWcs    = min(GWcs,0.4)
    GWcs  = 0.1 ! From GBHM, why?
    GWst    = wa/1000.
    nuz     = nl_soil

    isat    = nuz+1
    qsub    = 0.0
    tmpqsub = 0.0

    if((Dgl-Ds) .ge. 0.0 ) then
    ! groundwater table is bolow the toplayer
        avkg = tmpkg
        call gwriv(Dgl,length,slope,Ds,Dg,Dr,Drw,avkg,tmpqsub)   ! m3/s/m
        qsub = tmpqsub * dt / length                  ! m
! print*,'debug',qsub,GWst
        if(abs(qsub) .lt. 0.1E-20) qsub = 0.0
    else
    ! groundwater table reaches the toplayer
        if(Drw.ge.Dr .and. Dgl.le.0.0) goto 199
        tmp_D = 0.0
        avkg  = 0.0
        do i = nuz, 1, -1
          tmp_D = tmp_D + deltz(i)
          if( (Ds-tmp_D) .le. 0.5) then
              avkg  = avkg + anik*k0(i)*deltz(i)
          else
              avkg  = avkg + k0(i)*deltz(i)
          endif
          if(tmp_D .ge. (Ds-Dgl)) goto 50
        end do
50      if(Dr .gt. Ds) then
            tmp_D = tmp_D + Dr-Ds
            avkg  = avkg + tmpkg*(Dr-Ds)
        endif
        avkg  = avkg / tmp_D
        call gwriv(Dgl,length,slope,Ds,Dg,Dr,Drw,avkg,tmpqsub)   ! m3/s/m
        qsub = tmpqsub * dt / length                  ! m
        if(abs(qsub) .lt. 0.1E-20) qsub = 0.0
    endif

!   Renewing the groundwater table and toplayer soil moisture
    if(qsub .lt. 0.0) then      ! River infiltrates into aquifer
      GWst  = GWst - qsub
      if(GWst .le. Dg*GWcs) then! Ground water lever is below the top soil layer
        Dgl = Dgl + qsub/GWcs
        if(Dgl .lt. 0.0) print *, &
          'wrong in renewing groundwater table ... 1', Dgl
      else                      ! Ground water infilatrate into the soil layer
        tmp    = GWst - Dg*GWcs
        GWst   = Dg*GWcs
        Dgl    = Ds
        w(nuz) = w(nuz) + tmp / deltz(nuz)
        tmpw   = porsl(i)-wice(nuz)
        excess = (w(nuz)-tmpw) * deltz(nuz)
        excess = amax1(excess, 0.0)
        w(nuz) = w(nuz) - excess / deltz(nuz)
        GWkey  = 0
        if(abs(w(nuz)-porsl(i)) .lt. 0.1E-6) then
          Dgl = Ds - deltz(nuz)
          GWkey = 1
        else
          GWkey = 0
        endif
        do i=nuz-1, 1, -1
            w(i)   = w(i) + excess / deltz(i)
            tmpw   = porsl(i)-wice(i)
            excess = (w(i)-tmpw) * deltz(i)
            excess = amax1(excess, 0.0)
            w(i)   = w(i) - excess /deltz(i)
            if(GWkey.eq.1 .and. abs(w(i)-porsl(i)) .lt. 0.1E-5) then
                Dgl = Dgl-deltz(i)
            else
                GWkey = 0
            endif
        end do
      endif
  elseif(qsub .gt. 0.0) then                ! Aquifer flows out
        if( Dgl .ge. Ds-0.0001 ) then
            GWst  = GWst - qsub
            Dgl   = Dgl + qsub/GWcs
            if(Dgl.gt.(Ds+Dg) .or. GWst.lt.0.0) then
                qsub= qsub+GWst
                Dgl = Ds+Dg
                GWst= 0.0
            end if
        else
            tmp_D = 0.0
            isat  = nuz+1
            do i = nuz, 1, -1
                tmp_D = tmp_D + deltz(i)
                if( abs(tmp_D - (Ds-Dgl)) .le. 0.1E-3) goto 100
            end do
            goto 110
100         if(abs(w(i)-porsl(i)) .le. 0.1E-3)   isat = i
            goto 120
110         if(abs(w(i+1)-porsl(i)) .le. 0.1E-3) then
                isat = i+1
            end if
            !print *,"wrong in Dgl & qsub>0",isat,Dgl

120         tmp = 0.0
            do i = isat, nuz
                if(w(i)-porsl(i) .gt. 0.1E-3) &
                print *, "wrong in GW flow", isat, i,Dgl,w(i),porsl(i)
                tmp  = tmp + GWcs*deltz(i)
                w(i) = w(i) - GWcs
                Dgl  = Dgl + deltz(i)
                if(tmp .ge. qsub) goto 150
            end do
150         if(i.le.nuz) w(i) = w(i) + (tmp-qsub)/deltz(i)
            if(i.gt.nuz) then
                GWst  = GWst - (qsub-tmp)
                Dgl   = Dgl + (qsub-tmp)/GWcs
            end if
            if(i.gt.1 .and. w(i-1)-porsl(i) .gt. 0.1E-3) &
            print *,"w(i)>porsl(i)", i-1, w(i-1),porsl(i),isat,nuz
            if(Dgl.gt.Ds .and. abs(GWst-Dg*GWcs).lt.0.1E-3) Dgl=Ds
        endif
        if(Dgl.gt.(Ds+Dg) .or. GWst.lt.0.0) then
            qsub= qsub+GWst
            Dgl = Ds+Dg
            GWst= 0.0
        end if
    endif
    wliq_soisno(1:nl_soil)=w(1:nl_soil)*deltz(1:nl_soil)*1000.! m -> kg/m2
    qsub = qsub * 1000. / dt              ! m -> mm/sec
    wa   = GWst*1000.
    zwt  = Dgl
! print*,'qsub,wa,zwt', qsub,wa,zwt
199 return
  end subroutine ground_runoff


  subroutine MoistureFromSuction_V(w,wsat,wrsd,n,alpha,ps)
  !From GBHM by Dr.Dawen Yang
    use precision
    implicit none
    real(r8) w,wsat,wrsd,alpha,m,n,ps
    real(r8) se,tmpe,tmpps

    tmpps = 100.0*ps ! m->cm
    m     = 1.0-1.0/n
    tmpe  = 1.0+(alpha*abs(tmpps))**n
    se    = (1.0/tmpe)**m
    w     = se*(wsat-wrsd)+wrsd
    if(w .gt. wsat) w = wsat
    if(w .lt. wrsd) w = wrsd
    return
  end subroutine MoistureFromSuction_V

! ----------------------------------------------------------------------
!
!  Created by Yang D.; introduced by Li H. from GBHM,2015-12
!
!  calculation of exchange rate between groundwater and river: qsub
!  Yang D. et al., Hydrological Processes, 14, 403-416, 2000
!  The datum is sitted at the bottom of the unconfined aquifer
!  Varibles:
!  Dtg:      depth to groundwater (m)
!  length:   length of hillslope (m)
!  slope:    slope of hillslope (m)
!  Ds:       depth of top soil (m)
!  Dg:       depth of unconfined groundwater acquifer below topsoil (m)
!  Dr:       depth of river (m)
!  Drw:      depth of river water (m)
!  kg:       hydraulic conductivity (m/sec)
!  conlen:   contact length between the river and aquifer (m)
!  grad:     gradient of water head
!  Q:        discharge exchanged between aquifer and river (m^3/sec)
! ----------------------------------------------------------------------
  subroutine gwriv(Dtg,length,slope,Ds,Dg,Dr,Drw,kg,Q)
    use precision
    implicit none
    real(r8), intent(in) :: &
            Dtg,    &   !depth to groundwater (m)
            length, &   !length of hillslope  (m)
            slope,  &   !slope of hillslope   (m)
            Ds,     &   !depth of top soil    (m)
            Dg,     &   !depth of unconfined groundwater acquifer below topsoil (m)
            Dr,     &   !depth of river       (m)
            kg          !hydraulic conductivity (m/sec)
    real(r8), intent(inout) :: &
            Drw         !depth of river water (m)
    real(r8), intent(out) :: &
              Q           !discharge exchanged between aquifer and river (m^3/sec)
    real(r8) H1,H2,hs1,hs2,grad,conlen,Hrd

    Drw = amax1(Drw, 0.0)
    Drw = amin1(Dr, Drw)

    ! Distance from datum to riverbed (m)
    Hrd = Ds+Dg-Dr
    if(Dtg .lt. Dr-Drw) then
      H1  = 0.5*length*slope + Ds+Dg-Dtg
    else
      ! waterhead of groundwater (m)
      H1  = sqrt(0.5*length*slope) + Ds+Dg-Dtg
    endif
    ! saturated acquifer depth (m)
    hs1=H1-Hrd

    if(Hrd.ge.0.0) then
        ! waterhead of river (m)
        H2 =Hrd+Drw
        ! water depth in river (m)
        hs2=Drw
    else
        hs2=amax1(0.0, Drw)
        H2 =amax1(Hrd+Drw, 0.0)
    endif
    ! gradient of waterhead
    grad   =(H1-H2)/(0.5*length)
    ! contact length between the river and aquifer	(m)
    conlen =0.5*(abs(hs1)+hs2)
    ! discharge per unit width of a hillslope (m3/sec/m)
    Q=kg*grad*conlen
    if(Drw .le. 5E-3 .and. Q.lt.0.0) Q=0.0
    return
  end subroutine gwriv

  subroutine gbhm_recharge(dt,nl_soil,Dgl,Dg,Ds,wa,GWcs,wsat,wrsd,watern,alpha,wfld,&
                      k0,deltz,wliq_soisno,wice_soisno,&
                      rech)

  ! ----------------------------------------------------------------------
  ! calculation of recharge rate to groundwater:rech
  ! ----------------------------------------------------------------------
    use precision
    use PhysicalConstants, only: denice
    implicit none
    integer,  intent(in) :: &
               nl_soil
    real(r8), intent(in) :: &
               dt, &
               k0(1:nl_soil),          &
               deltz(1:nl_soil),       &
               wice_soisno(1:nl_soil), &
               Ds        , &
               Dg        , &!depth of unconfined acquifer (m)
               GWcs      , &
               wsat      , &
               wrsd      , &
               watern    , &
               alpha     , &
               wfld(1:nl_soil)
    real(r8), intent(in) :: &
               Dgl,   &!depth to GW level (m)
               wa
    real(r8), intent(in) :: &
               wliq_soisno(1:nl_soil)  ! liquid water (kg/m2)
    real(r8), intent(out) :: &
               rech
    integer  nuz,i
    real(r8) GWst,tmp,ps1,ps2,tmp2,tmpw,deficit,w_sum0,w_sum1,erro
    real(r8) wliq(1:nl_soil) !Soil moisture,%
    real(r8) wice(1:nl_soil)
    real(r8) hk_soil_end   !the hydraulic conductivity at the end soil layer [mm h2o/s]
    real(r8) max_rech_soil ! the maximum recharge from the end soil layer
    real(r8) max_rech_aqf  ! the maximum recharge filled by the void of aquifer

    rech = 0.0
    nuz  = nl_soil
    GWst = wa/1000. ! from mm to m
    wliq(1:nl_soil) = wliq_soisno(1:nl_soil)*0.001/deltz(1:nl_soil)
    wice(1:nl_soil) = wice_soisno(1:nl_soil)/denice/deltz(1:nl_soil)
    if(Dgl .ge. Ds) then
       !calculate suction from soil moisture by Van Genuchten's equation
       call SuctionFromMoisture_V(wliq(nuz),wsat,wrsd,watern,alpha,ps1)  ! m
       tmp = wfld(nl_soil)*(1.0 - (Dgl-Ds)/Dg)
       tmp = amax1(tmp, wrsd)
       call SuctionFromMoisture_V(tmp, wsat,wrsd,watern,alpha,ps2)    ! m
       call conductivity_V(k0(nuz),wsat,wrsd,watern,wliq(nuz),hk_soil_end)
       tmp  = 0.5*(deltz(nuz) + (Dgl-Ds))
       rech = hk_soil_end*(1.-(ps2-ps1)/tmp)*dt                ! m

       if(rech .ge. 0.0) then
         max_rech_soil = max((wliq(nuz)-wfld(nuz))*deltz(nuz),0.0)
         max_rech_aqf  = max(Dg*GWcs-GWst, 0.0)
         rech = min(rech, max_rech_soil, max_rech_aqf)
       else
         max_rech_aqf  = min(0.0, -1.0*GWst)
         ! rech = max(rech, -1.0*GWst)
         ! tmpw = wsat-wice(nuz)
         max_rech_soil = min(0.0, -1.*(wsat-wice(nuz)-wliq(nuz))*deltz(nuz))
         rech = max(rech, max_rech_aqf, max_rech_soil)
       endif
    endif
       rech = rech*1000./dt !transfer to mm/s
  end subroutine gbhm_recharge

! ----------------------------------------------------------------------
! calculate suction from soil moisture by Van Genuchten's equation
! ----------------------------------------------------------------------
subroutine SuctionFromMoisture_V(w,wsat,wrsd,n,alpha,ps)
    use precision
    implicit none
  real(r8) w,wsat,wrsd,alpha,m,n,se,ps
  real(r8) tmpe, tmpps,tmpw
  m = 1.0-1.0/n
    tmpw = w
  if(tmpw.ge.wsat) then
    ps=0.0
  elseif(tmpw.lt.wsat) then
    if(tmpw.lt.wrsd+0.001) tmpw=wrsd+0.001
    se=(tmpw-wrsd)/(wsat-wrsd)
    tmpe=se**(1.0/m)
    tmpps=-(1.0/tmpe-1.0)**(1.0/n)/alpha
      ps = tmpps/100.0 ! cm->m
  end if
  return
end subroutine SuctionFromMoisture_V

  subroutine tridia (n, a, b, c, r, u)

      use precision
      implicit none
      integer, intent(in) :: n        !length of diagonal element vector
      real(r8), intent(in) :: a(1:n)  !subdiagonal elements
      real(r8), intent(in) :: b(1:n)  !diagonal elements
      real(r8), intent(in) :: c(1:n)  !superdiagonal elements
      real(r8), intent(in) :: r(1:n)  !right hand side
      real(r8), intent(out) :: u(1:n) !solution vector

      integer j
      real(r8) gam(1:n),bet
! -----------------------------------------------------------------

      bet = b(1)
      u(1) = r(1) / bet
      do j = 2, n
            gam(j) = c(j-1) / bet
            bet = b(j) - a(j) * gam(j)
            u(j) = (r(j) - a(j)*u(j-1)) / bet

      end do
      do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
      end do

  end subroutine tridia

  SUBROUTINE TRDIG (N,DL,DM,DU,RS,X,MARK)
    use precision
    implicit none
!  	parameter (N=3)
      INTEGER MARK,N,I
      REAL(r8) DL(1:N),DM(1:N),DU(1:N),RS(1:N),X(1:N)
      real(r8) ROW,D
      MARK = -1
      IF (N .LT. 3) goto 999

!   checking for strong nonsingularity with N=1
      MARK = 0
      ROW = ABS(DM(1)) + ABS(DU(1))
      IF (ROW .EQ. 0.0E0) goto 999
      D = 1.0E0/ROW
      IF (ABS(DM(1))*D .LE. 1.0E-20) goto 999

!   factoring A while checking for strong nonsingularity
      DL(1) = 0.0E0
      DU(N) = 0.0E0
      DU(1) = DU(1)/DM(1)
      DO I=2,N,1
         ROW = ABS(DL(I)) + ABS(DM(I)) + ABS(DU(I))
         IF (ROW .EQ. 0.0E0) goto 999
         D = 1.0E0/ROW
         DM(I) = DM(I) - DL(I) * DU(I-1)
         IF (ABS(DM(I))*D .LE. 1.0E-20) goto 999
         IF (I .LT. N) THEN
            DU(I) = DU(I)/DM(I)
         ENDIF
      END DO
      MARK=1

999   continue

!  If MARK = 1, update the right hand side and solve via backsubstitution
      IF (MARK .EQ. 1) THEN
        RS(1) = RS(1)/DM(1)
        DO I=2,N,1
          RS(I) = (RS(I) - DL(I) * RS(I-1)) / DM(I)
        END DO
!       backsubstitution
        X(N) = RS(N)
        DO I=N-1,1,-1
          X(I) = RS(I) - DU(I) * X(I+1)
        END DO
      ENDIF
!	print *,MARK, X
      RETURN
END SUBROUTINE TRDIG
!
!-----------------------------------------------------------------------
!
!  calculate soil hydraulic conductivity by Van Genuchten's equation
!
!-----------------------------------------------------------------------
!
  subroutine conductivity_V(k0,wsat,wrsd,n,w,k)
  use precision
  implicit none
  real(r8),intent(in) :: k0,wsat,wrsd,n,w
  real(r8),intent(out):: k !m/s
  real(r8) m,tmpw,se

  m=1.0-1.0/n
  tmpw = w
  if(tmpw.gt.wsat) tmpw=wsat

    if(tmpw.le.wrsd) then
    tmpw=wrsd+1.0E-10
  elseif(tmpw.lt.wrsd+0.0001)then
    tmpw=wrsd+0.0001
  end if

  se=(tmpw-wrsd)/(wsat-wrsd)
  k=k0*sqrt(se)*(1.0-(1.0-se**(1.0/m))**m)**2
  if(k.lt.0.0 .or. k.gt.k0) &
        print *,'wrong in calculating conductivity',w,k
  end subroutine conductivity_V


  SUBROUTINE zwt_to_jwt(nl_soil,zwt,zi_soisno,jwt)
  use precision
  implicit none
  integer,  intent(in) :: nl_soil
  real(r8), intent(in) :: &
              zwt        ,&
              zi_soisno(0:nl_soil) ! interface level below a "z" level (m)
  integer, intent(out) :: jwt
  integer j
  ! jwt: The layer index of the first unsaturated layer,
  ! i.e., the layer right above the water table
  jwt = nl_soil
  do j = 1, nl_soil
    if(zwt <= zi_soisno(j)) then
      jwt = j-1
      exit
    end if
  enddo
  end SUBROUTINE zwt_to_jwt

END MODULE SOIL_SNOW_hydrology
! --------- EOP ----------
