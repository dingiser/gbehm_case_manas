
 subroutine meltf (itypwat,lb,nl_soil,deltim,psi0,porsl, wrsd,bsw,&
                   fact,brr,hs,dhsdT, dz_soisno,&
                   t_soisno_bef,t_soisno,wliq_soisno,wice_soisno,imelt, &
                   scv,snowdp,sm,xmf,sab_soisno)

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, /09/1999/, /03/2014/
! Revision History:
! Hongyi Li, add supercool water scheme from CLM,    2016-4-29 
! Hongyi Li, add solar radiation in different layer, 2018-2-7 
!
!-----------------------------------------------------------------------
! calculation of the phase change within snow and soil layers:
! 
! (1) check the conditions which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e., melting),
!     or layer temperature is less than the freezing point
!     and the liquid water mass is not equal to zero (i.e., freezing);
! (2) assess the rate of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point;
! (3) re-adjust the ice and liquid mass, and the layer temperature
!
!-----------------------------------------------------------------------

  use precision
  use PhysicalConstants, only : tfrz, hfus, grav, denh2o, denice
  IMPLICIT NONE

!-----------------------------------------------------------------------

   integer, INTENT(in) :: nl_soil             ! upper bound of array (i.e., soil layers)
   integer, INTENT(in) :: lb                  ! lower bound of array (i.e., snl +1)
   integer, INTENT(in) :: itypwat             ! land water type (0=soil, 1=urban or built-up, 2=wetland,
                      !                        3=glacier/ice sheet, 4=land water bodies)
  real(r8), INTENT(in) :: deltim              ! time step [second]
  real(r8), INTENT(in) :: t_soisno_bef(lb:nl_soil)  ! temperature at previous time step [K]
  real(r8), INTENT(in) :: brr (lb:nl_soil)    ! 
  real(r8), INTENT(in) :: fact(lb:nl_soil)    ! temporary variables
  real(r8), INTENT(in) :: hs                  ! net ground heat flux into the surface
  real(r8), INTENT(in) :: dhsdT               ! temperature derivative of "hs"
  real(r8), INTENT(in) :: psi0(1:nl_soil)     ! minimum soil suction [mm]
  real(r8), INTENT(in) :: porsl(1:nl_soil)    ! saturated volumetric soil water content(porosity)
  real(r8), INTENT(in) :: wrsd
  real(r8), INTENT(in) :: dz_soisno(lb:nl_soil) ! layer thickness (m)
  real(r8), INTENT(in) :: bsw(1:nl_soil)      ! Clapp-Hornberger "B"
  real(r8), INTENT(in) :: sab_soisno(lb:1)    !solar radiation absorbed by snow layers and the 1th ground [W/m2]

  real(r8), INTENT(inout) :: t_soisno (lb:nl_soil) ! temperature at current time step [K]
  real(r8), INTENT(inout) :: wice_soisno(lb:nl_soil) ! ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq_soisno(lb:nl_soil) ! liquid water [kg/m2]
  real(r8), INTENT(inout) :: scv              ! snow mass [kg/m2]
  real(r8), INTENT(inout) :: snowdp           ! snow depth [m]

  real(r8), INTENT(out) :: sm                 ! rate of snowmelt [mm/s, kg/(m2 s)]
  real(r8), INTENT(out) :: xmf                ! total latent heat of phase change
   integer, INTENT(out) :: imelt(lb:nl_soil)  ! flag for melting or freezing [-]

! Local 
  real(r8) :: hm(lb:nl_soil)                  ! energy residual [W/m2]
  real(r8) :: xm(lb:nl_soil)                  ! metling or freezing within a time step [kg/m2]
  real(r8) :: heatr                           ! energy residual or loss after melting or freezing
  real(r8) :: temp1                           ! temporary variables [kg/m2]
  real(r8) :: temp2                           ! temporary variables [kg/m2]

  real(r8), dimension(lb:nl_soil) :: wmass0, wice0, wliq0
  real(r8) :: propor, tinc, we, scvold,supercool(nl_soil)
  real(r8) :: smp   
  integer j

!-----------------------------------------------------------------------

  sm = 0.
  xmf = 0.
  supercool(1:nl_soil) = 0.0 ! added by HY, 2016-4-28
  do j = lb, nl_soil
     imelt(j) = 0
     hm(j) = 0.
     xm(j) = 0.
     wice0(j) = wice_soisno(j)
     wliq0(j) = wliq_soisno(j)
     wmass0(j) = wice_soisno(j) + wliq_soisno(j)
  enddo

  scvold=scv
  we=0.
  if(lb<=0) we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))

  do j = lb, nl_soil
     ! Melting identification
     ! if ice exists above melt point, melt some to liquid.
     if(wice_soisno(j) > 0. .and. t_soisno(j) > tfrz)then
        imelt(j) = 1
        t_soisno(j) = tfrz
     endif

     ! Freezing identification
     ! if liquid exists below melt point, freeze some to ice.
     if(j<1 .and. wliq_soisno(j) > 0. .and. t_soisno(j) < tfrz) then
        imelt(j) = 2
        t_soisno(j) = tfrz
     endif

     ! The following code is from Gaobing or CLM???, Added and modified by HY, 2016-4-28
     ! from Zhao (1997) and Koren (1999)     
!           if (iland /=10) then
     if(itypwat<=4)then   ! soil ground only, need to check it carefully.
        if(j>=1 .and. t_soisno(j) < tfrz) then
           smp = hfus*(tfrz-t_soisno(j))/(grav*t_soisno(j)) * 1000.0  !(mm)
           supercool(j) = porsl(j)*(smp/psi0(j))**(-1.0/bsw(j))
           supercool(j) = supercool(j)*dz_soisno(j)*1000.0       ! (mm)
           supercool(j) = max(supercool(j),wrsd*dz_soisno(j)*1000.0)
           ! supercool(j) = max(supercool(j),0.07*dz_soisno(j)*1000.0)
        endif
     endif
 
     if (j>=1 .and. wliq_soisno(j) > supercool(j) .AND. t_soisno(j) < tfrz) then
        imelt(j) = 2
        t_soisno(j) = tfrz
     endif     
  enddo

! If snow exists, but its thickness less than the critical value (0.01 m)
  if(lb == 1 .and. scv > 0.)then
     if(t_soisno(1) > tfrz)then
        imelt(1) = 1
        t_soisno(1) = tfrz
     endif
  endif

! Calculate the energy surplus and loss for melting and freezing
  do j = lb, nl_soil
     if(imelt(j) > 0)then
        tinc = t_soisno(j)-t_soisno_bef(j)
        if(j > lb)then
           hm(j) = brr(j) - tinc/fact(j) 
           ! snow layer or top soil layer (add solar radiation)
           if (j <= 1) hm(j) = hm(j) + sab_soisno(j)           
        else
           hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j) 
        endif
     endif
  enddo

  do j = lb, nl_soil
     if(imelt(j) == 1 .and. hm(j) < 0.) then
       hm(j) = 0.
       imelt(j) = 0
     endif
! this error was checked carefully, it results from the the computed error
! of "Tridiagonal-Matrix" in subroutine "thermal".
     if(imelt(j) == 2 .and. hm(j) > 0.) then
       hm(j) = 0.
       imelt(j) = 0
     endif
  enddo

! The rate of melting and freezing
  do j = lb, nl_soil
     if(imelt(j) > 0 .and. abs(hm(j)) > .0) then
        xm(j) = hm(j)*deltim/hfus                    ! kg/m2

        ! if snow exists, but its thickness less than the critical value (1 cm)
        ! Note: more work is need on how to tune the snow depth at this case
        if(j == 1 .and. lb == 1 .and. scv > 0. .and. xm(j) > 0.)then
           temp1 = scv                               ! kg/m2
           scv = max(0.,temp1-xm(j))
           propor = scv/temp1
           snowdp = propor * snowdp
           heatr = hm(j) - hfus*(temp1-scv)/deltim   ! W/m2
           if(heatr > 0.) then
              xm(j) = heatr*deltim/hfus              ! kg/m2
              hm(j) = heatr                          ! W/m2
           else
              xm(j) = 0.
              hm(j) = 0.
           endif
           sm = max(0.,(temp1-scv))/deltim           ! kg/(m2 s)
           xmf = hfus*sm
        endif

        heatr = 0.
        if(xm(j) > 0.) then
           wice_soisno(j) = max(0., wice0(j)-xm(j))
           heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
        else if (xm(j) < 0.) then
          if (j <= 0 ) then ! snow layer or other land type except soil land..(or. itypwat>1)
            wice_soisno(j) = min(wmass0(j), wice0(j)-xm(j)) 
            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
          else
            if (wmass0(j) < supercool(j)) then
              wice_soisno(j) = 0.
            else
              wice_soisno(j) = min(wmass0(j) - supercool(j),wice0(j)-xm(j))
              ! wice_soisno(j) = min((wmass0(j)/denh2o-wrsd*dz_soisno(j))*denice, wice_soisno(j)) ! changed by gaobing  
            endif
            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim  
          endif
        endif

        wliq_soisno(j) = max(0.,wmass0(j)-wice_soisno(j))

        if(abs(heatr) > 0.)then
           if(j > lb)then
              t_soisno(j) = t_soisno(j) + fact(j)*heatr
           else
              t_soisno(j) = t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
           endif
           if (j<=0) then ! limit the judge scope in snow layers
              if(wliq_soisno(j)*wice_soisno(j) > 0.) t_soisno(j) = tfrz
           endif 
        endif

        xmf = xmf + hfus * (wice0(j)-wice_soisno(j))/deltim

        if(imelt(j) == 1 .and. j < 1) &
        sm = sm + max(0.,(wice0(j)-wice_soisno(j)))/deltim  

     endif
  enddo

  !scvold=scv
  if(lb<=0) then
  we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))-we
     if(abs(we)>1.e-6) then
        print*, 'meltf err : ', we
     endif
  endif

!------------------------------------------

!   sm = 0.
!   xmf = 0.
!   do j = lb, nl_soil
!      imelt(j) = 0
!      hm(j) = 0.
!      xm(j) = 0.
!      wice0(j) = wice_soisno(j)
!      wliq0(j) = wliq_soisno(j)
!      wmass0(j) = wice_soisno(j) + wliq_soisno(j)
!   enddo

!   scvold=scv
!   we=0.
!   if(lb<=0) we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))

!   do j = lb, nl_soil
!      ! Melting identification
!      ! if ice exists above melt point, melt some to liquid.
!      if(wice_soisno(j) > 0. .and. t_soisno(j) > tfrz)then
!         imelt(j) = 1
!         t_soisno(j) = tfrz
!      endif

!      ! Freezing identification
!      ! if liquid exists below melt point, freeze some to ice.
!      if(wliq_soisno(j) > 0. .and. t_soisno(j) < tfrz) then
!         imelt(j) = 2
!         t_soisno(j) = tfrz
!      endif
!   enddo

! ! If snow exists, but its thickness less than the critical value (0.01 m)
!   if(lb == 1 .and. scv > 0.)then
!      if(t_soisno(1) > tfrz)then
!         imelt(1) = 1
!         t_soisno(1) = tfrz
!      endif
!   endif

! ! Calculate the energy surplus and loss for melting and freezing
!   do j = lb, nl_soil
!      if(imelt(j) > 0)then
!         tinc = t_soisno(j)-t_soisno_bef(j)
!         if(j > lb)then
!            hm(j) = brr(j) - tinc/fact(j) 
!         else
!            hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j) 
!         endif
!      endif
!   enddo

!   do j = lb, nl_soil
!      if(imelt(j) == 1 .and. hm(j) < 0.) then
!        hm(j) = 0.
!        imelt(j) = 0
!      endif
! ! this error was checked carefully, it results from the the computed error
! ! of "Tridiagonal-Matrix" in subroutine "thermal".
!      if(imelt(j) == 2 .and. hm(j) > 0.) then
!        hm(j) = 0.
!        imelt(j) = 0
!      endif
!   enddo

! ! The rate of melting and freezing
!   do j = lb, nl_soil
!      if(imelt(j) > 0 .and. abs(hm(j)) > .0) then
!         xm(j) = hm(j)*deltim/hfus                    ! kg/m2

!         ! if snow exists, but its thickness less than the critical value (1 cm)
!         ! Note: more work is need on how to tune the snow depth at this case
!         if(j == 1 .and. lb == 1 .and. scv > 0. .and. xm(j) > 0.)then
!            temp1 = scv                               ! kg/m2
!            scv = max(0.,temp1-xm(j))
!            propor = scv/temp1
!            snowdp = propor * snowdp
!            heatr = hm(j) - hfus*(temp1-scv)/deltim   ! W/m2
!            if(heatr > 0.) then
!               xm(j) = heatr*deltim/hfus              ! kg/m2
!               hm(j) = heatr                          ! W/m2
!            else
!               xm(j) = 0.
!               hm(j) = 0.
!            endif
!            sm = max(0.,(temp1-scv))/deltim           ! kg/(m2 s)
!            xmf = hfus*sm
!         endif

!         heatr = 0.
!         if(xm(j) > 0.) then
!            wice_soisno(j) = max(0., wice0(j)-xm(j))
!            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
!         else
!            wice_soisno(j) = min(wmass0(j), wice0(j)-xm(j))
!            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim  
!         endif

!         wliq_soisno(j) = max(0.,wmass0(j)-wice_soisno(j))

!         if(abs(heatr) > 0.)then
!            if(j > lb)then
!               t_soisno(j) = t_soisno(j) + fact(j)*heatr
!            else
!               t_soisno(j) = t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
!            endif
!            if(wliq_soisno(j)*wice_soisno(j) > 0.) t_soisno(j) = tfrz
!         endif

!         xmf = xmf + hfus * (wice0(j)-wice_soisno(j))/deltim

!         if(imelt(j) == 1 .and. j < 1) &
!         sm = sm + max(0.,(wice0(j)-wice_soisno(j)))/deltim  

!      endif
!   enddo

!   !scvold=scv
!   if(lb<=0) then
!   we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))-we
!      if(abs(we)>1.e-6) then
!         print*, 'meltf err : ', we
!      endif
!   endif


 end subroutine meltf
