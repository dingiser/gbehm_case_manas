  subroutine blowing_snow(dim_y,dim_x,nl_soil,maxsnl,deltim,forc_hgt_u,& ! input
            region_ws_x,region_ws_y,region_airt,region_rh,region_grainsize,region_fetch,& !input
            t_soisno,wliq_soisno,wice_soisno,dz_soisno,z_soisno,scv, snowdep,& !in-out 
            region_qsub,net_blow,blowing_add_region) !output
    use precision
    use blowing_varpar
    use SNOW_Layers_CombineDivide
    IMPLICIT NONE
     
    integer,  intent(in) :: dim_y,dim_x,nl_soil,maxsnl
    real(r8), intent(in) :: deltim 
    real(r8), intent(in) :: region_ws_x(dim_y,dim_x),region_ws_y(dim_y,dim_x),forc_hgt_u(dim_y,dim_x),&
                           region_fetch(dim_y,dim_x)
    real(r8), intent(in) :: region_airt(dim_y,dim_x)  ! air temperature, K
    real(r8), intent(in) :: region_rh(dim_y,dim_x)
    real(r8), intent(in) :: region_grainsize(dim_y,dim_x)
    real(r8), intent(out):: &
              region_qsub(dim_y,dim_x), &
              net_blow(dim_y,dim_x)   , &  !net accumulation because of blowing snow. 
              blowing_add_region
    real(r8), intent(inout) :: t_soisno(maxsnl+1:nl_soil,dim_y,dim_x)
    real(r8), intent(inout) :: wliq_soisno(maxsnl+1:nl_soil,dim_y,dim_x), wice_soisno(maxsnl+1:nl_soil,dim_y,dim_x)
    real(r8), intent(inout) :: dz_soisno(maxsnl+1:nl_soil,dim_y,dim_x), z_soisno(maxsnl+1:nl_soil,dim_y,dim_x)
    real(r8), intent(inout) :: scv(dim_y,dim_x),snowdep(dim_y,dim_x)
    real(r8) tempmass,ratio     
    real(r8) :: region_qtx(dim_y,dim_x),region_qty(dim_y,dim_x),region_qt(dim_y,dim_x)!temp_swe(dim_y,dim_x),
    real(r8) :: region_qsub_max(dim_y,dim_x),region_qt_salt_x(dim_y,dim_x),region_qt_salt_y(dim_y,dim_x)
    real(r8) :: region_qt_susp_x(dim_y,dim_x),region_qt_susp_y(dim_y,dim_x)
    real(r8) :: grid_qt_max_x(3),grid_qt_max_y(3),grid_qt_x,grid_qt_y,grid_ws_x(3),grid_ws_y(3),&
                grid_swe_x(3),grid_swe_y(3),fetch
    real(r8) :: delta_swe,temp1,temp2,temp3,temp4,tempD,debug1
    real(r8) :: zi_soisno(maxsnl:nl_soil) ! interface level below a "z" level (m)
    integer :: i,j,ii,jj,snl,lb,k

!f2py intent(in) dim_y,dim_x,nl_soil,maxsnl
!f2py intent(in) deltim 
!f2py intent(in) region_ws_x,region_ws_y,forc_hgt_u,region_fetch
!f2py intent(in) region_airt
!f2py intent(in) region_rh
!f2py intent(in) region_grainsize
!f2py intent(in,out) t_soisno,wliq_soisno, wice_soisno
!f2py intent(in,out) dz_soisno, z_soisno
!f2py intent(in,out) scv,snowdep
!f2py intent(out) region_qsub,net_blow,blowing_add_region
!f2py depend(dim_y,dim_x) region_ws_x,region_ws_y,forc_hgt_u,region_fetch,region_airt 
!f2py depend(dim_y,dim_x) region_rh ,region_grainsize,region_qsub,net_blow,scv,snowdep   
!f2py depend(maxsnl,nl_soil,dim_y,dim_x) t_soisno,wliq_soisno, wice_soisno,dz_soisno, z_soisno 
    ! grid_size = 1000. !m. temporal assumption
    
    call region_blowing_max(dim_y,dim_x,region_ws_x,region_ws_y,region_airt,region_rh,region_grainsize,&
            forc_hgt_u,region_qsub_max,region_qt_salt_x,region_qt_susp_x,region_qt_salt_y,region_qt_susp_y)
    blowing_add_region = 0. !区域外吹雪的质量补充
    do j=1,dim_x
        do i=1,dim_y
            if (scv(i,j).le.0.00001) then 
                region_qtx(i,j)  = 0.
                region_qty(i,j)  = 0.
                region_qsub(i,j) = 0.
            else
                !在边界上，如果格网有雪，则假设上风向有雪；如果格网无雪，则假设上风向无雪。
                !可结合遥感数据获取更准确的判断。
                if (i==1 .or. j==1) then
                    region_qt(i,j) = ((region_qt_salt_x(i,j) + region_qt_susp_x(i,j))**2&
                                + (region_qt_salt_y(i,j)+region_qt_susp_y(i,j))**2)**0.5
                    blowing_add_region = blowing_add_region+region_qt(i,j)
                else if (i== dim_y .or. j==dim_x) then
                    region_qt(i,j) = ((region_qt_salt_x(i,j) + region_qt_susp_x(i,j))**2&
                                + (region_qt_salt_y(i,j)+region_qt_susp_y(i,j))**2)**0.5
                    blowing_add_region = blowing_add_region-region_qt(i,j)
                else                !分为8个方向区间求解，会更加符合实际；但这里为了简化计算，仍然采用四方向，可采用降低网格分辨率解决这一问题
                    fetch = region_fetch(i,j)
                    grid_swe_x = (/scv(i,j-1),scv(i,j),scv(i,j+1)/)
                    grid_ws_x = (/region_ws_x(i,j-1),region_ws_x(i,j),region_ws_x(i,j+1)/)
                    grid_qt_max_x = (/region_qt_salt_x(i,j-1),&
                                    region_qt_salt_x(i,j),&
                                    region_qt_salt_x(i,j+1)/) &
                                    +(/region_qt_susp_x(i,j-1),&
                                    region_qt_susp_x(i,j),&
                                    region_qt_susp_x(i,j+1)/) 
                    call grid_blowing1(grid_swe_x,grid_ws_x,grid_qt_max_x,grid_size,fetch, & 
                                grid_qt_x)
                    region_qtx(i,j) = grid_qt_x  
                    
                    grid_swe_y = (/scv(i-1,j),scv(i,j),scv(i+1,j)/)
                    grid_ws_y = (/region_ws_y(i-1,j),region_ws_y(i,j),region_ws_y(i+1,j)/)
                    grid_qt_max_y = (/region_qt_salt_y(i-1,j),&
                                        region_qt_salt_y(i  ,j),&
                                        region_qt_salt_y(i+1,j)/)&
                                    +(/region_qt_susp_y(i-1,j),&
                                        region_qt_susp_y(i  ,j),&
                                        region_qt_susp_y(i+1,j)/) 
                    call grid_blowing1(grid_swe_y,grid_ws_y,grid_qt_max_y,grid_size,fetch, & 
                                        grid_qt_y) 
                    region_qty(i,j) = grid_qt_y                            
                    region_qt(i,j) = (grid_qt_x**2+grid_qt_y**2)**0.5
                end if
                if ((abs(grid_qt_max_x(2))+abs(grid_qt_max_y(2))) <= 0.00000000001) then
                    region_qsub(i,j) = 0.
                else
                    region_qsub(i,j) = region_qsub_max(i,j)
                    ! region_qsub(i,j) = region_qsub_max(i,j)*&
                                        ! (abs(region_qtx(i,j))+abs(region_qty(i,j)))/ &
                                        ! (abs(grid_qt_max_x(2))+abs(grid_qt_max_y(2)))
                endif  
            endif
            ! if the availabe snow mass is less than the computed blowing snow, limit the value
            ! if ((region_qtx(i,j)+region_qty(i,j)+scv(i,j)).lt.0.) then
                ! region_qtx(i,j) = region_qtx(i,j)*scv(i,j)/abs(region_qtx(i,j)+region_qty(i,j))
                ! region_qty(i,j) = region_qty(i,j)*scv(i,j)/abs(region_qtx(i,j)+region_qty(i,j))
            ! endif
       end do
    end do
    
!update SWE distribution

    do j=1,dim_x
        do i=1,dim_y
            net_blow(i,j) = 0.
            
            if (i==1 .or. j==1) then
                delta_swe=0.5*deltim* &
                (region_qtx(i,j)-region_qtx(i,j+1)+region_qty(i,j)-region_qty(i+1,j))&
                -region_qsub(i,j)*deltim               
            else if (i== dim_y .or. j==dim_x) then
                delta_swe=0.5*deltim* &
                (region_qtx(i,j-1)-region_qtx(i,j)+region_qty(i-1,j)-region_qty(i,j))&
                -region_qsub(i,j)*deltim
            else
                delta_swe = 0.5*deltim* &
                (region_qtx(i,j-1)-region_qtx(i,j+1)+region_qty(i-1,j)-region_qty(i+1,j))&
                -region_qsub(i,j)*deltim
            end if

            snl = 0
            do jj=maxsnl+1,0
                if((wliq_soisno(jj,i,j)+wice_soisno(jj,i,j))>0.) then
                   snl=snl-1
                   ! print*,snl,wliq_soisno(jj,i,j)+wice_soisno(jj,i,j)  
                endif
            enddo
            lb=snl+1 
            !update SWE
            if (snl<0) then
                lb = snl+1
                scv(i,j) = sum(wliq_soisno(lb:0,i,j) + wice_soisno(lb:0,i,j))
            else
                scv(i,j) = 0.
            endif
            
            if (abs(delta_swe).lt.0.00001) then
                delta_swe = 0.        
            elseif (delta_swe.gt.0.) then     !accumulation  
                ! -------------------------------------------------------------------------
                ! In snow-free region, the added blowing snow is considered as a new layer
                ! -------------------------------------------------------------------------         
                if (snl.eq.0) then
                    snl = -1
                    lb  = 0
                    dz_soisno  (lb,i,j) = delta_swe/200. !mm to m
                    z_soisno   (lb,i,j) = -0.5*dz_soisno(lb,i,j)
                    wice_soisno(lb,i,j) = delta_swe
                    wliq_soisno(lb,i,j) = 0.
                    t_soisno   (lb,i,j) = min(region_airt(i,j),273.1)                    
                else
                    dz_soisno(lb,i,j) = dz_soisno(lb,i,j)* &
                                    (delta_swe+wice_soisno(lb,i,j)+wliq_soisno(lb,i,j))&
                                    /(wice_soisno(lb,i,j)+wliq_soisno(lb,i,j))
                    if (lb.ge.0) then
                        z_soisno (lb,i,j) = - 0.5*dz_soisno(lb,i,j)
                    else    
                        z_soisno (lb,i,j) = z_soisno(lb+1,i,j) - 0.5*dz_soisno(lb,i,j)
                    endif    
                    wice_soisno(lb,i,j) = wice_soisno(lb,i,j) + delta_swe
              
                endif    
                net_blow(i,j) = delta_swe
            elseif (delta_swe.lt.0. .and. (scv(i,j)+delta_swe).gt.0.) then
                tempmass = 0.
                k = 0
                ii= 0
                do while (tempmass<=abs(delta_swe))
                    k  = k+1
                    ii = snl+k
                    tempmass = wliq_soisno(ii,i,j) + wice_soisno(ii,i,j) + tempmass
                end do
if (ii>0) print*,'wrong ii in line 177',&
ii,snl,lb,tempmass,delta_swe,scv(i,j),wliq_soisno(-4:1,i,j),wice_soisno(-4:1,i,j)
                !noting the unit of deltaZ is 'm'
                ratio = (tempmass+delta_swe)/(wliq_soisno(ii,i,j) + wice_soisno(ii,i,j))
                wliq_soisno(lb:ii-1,i,j) = 0.
                wice_soisno(lb:ii-1,i,j) = 0.
                if (ratio>1..or.ratio<0.) then
print*,'ratio=',ratio                
                    ! In case that leaved snow is too thin or other cases,
                    ! assuming the whole ii snow layer is blowed off, update drifted snow mass.
                    net_blow(i,j) = -1.*tempmass
                    wliq_soisno(ii,i,j) = 0.
                    wice_soisno(ii,i,j) = 0.
                    snl = ii               !since number of snow layer is from 0s 
                else
                    ! In normal case
                    dz_soisno(ii,i,j)   = dz_soisno(ii,i,j)*ratio
                    net_blow(i,j)       = delta_swe
                    wliq_soisno(ii,i,j) = ratio*wliq_soisno(ii,i,j)
                    wice_soisno(ii,i,j) = ratio*wice_soisno(ii,i,j)    
                    ! snowdep(i,j)        = sum(dz_soisno(ii:0,i,j))
                    if (ii.ge.0) then
                        z_soisno(ii,i,j) =  - 0.5*dz_soisno(ii,i,j)
                    else    
                        z_soisno(ii,i,j) = z_soisno(ii+1,i,j) - 0.5*dz_soisno(ii,i,j) &
                                          - 0.5*dz_soisno(ii+1,i,j)
                    endif                       
                    snl = ii-1             !since number of snow layer is from 0s  
                endif
            elseif (delta_swe.lt.0. .and. (scv(i,j)+delta_swe).le.0.) then 
                if (scv(i,j).le.0 .or. lb.ge.1) then 
                    lb = 0
                    scv(i,j) = 0.
                else
                    wliq_soisno(lb:0,i,j) = 0.
                    wice_soisno(lb:0,i,j) = 0. 
                endif
                net_blow   (i,j) = -1.*scv(i,j)
                snl = 0

                if (i==1 .or. j==1 .or. i== dim_y .or. j==dim_x) then

                else
                    temp1 = (delta_swe/deltim)
                    temp2 = 0.5*(abs(region_qtx(i,j-1))+abs(region_qtx(i,j+1))+abs(region_qty(i-1,j))+&
                                 abs(region_qty(i+1,j))) +abs(region_qsub(i,j))
                    temp3 = abs(scv(i,j)+delta_swe)
                    tempD = (temp2-temp1)/2.
                    ! if (scv(i,j).le.0.) then
                    
if (tempD<0.) print*,tempD,temp3,temp2,temp1,scv(i,j),delta_swe,temp3/deltim
                    if (region_qtx(i,j-1)<0. .and. tempD.gt.0.) then
                        temp4 = -temp3*0.5*region_qtx(i,j-1)/tempD
                        region_qtx(i,j-1) = region_qtx(i,j-1)  +  &
                                            temp3*0.5*region_qtx(i,j-1)/tempD/deltim
                        call add_mass(temp4,wliq_soisno(:,i,j-1),wice_soisno(:,i,j-1),dz_soisno(:,i,j-1),&
                                     z_soisno(:,i,j-1),t_soisno(:,i,j-1),scv(i,j-1),region_airt(i,j-1),nl_soil,maxsnl)
                                    
                                     
                        net_blow(i,j-1)   = net_blow(i,j-1)+temp4                        
                        ! scv(i,j-1)        = scv(i,j-1)-temp3*0.5*region_qtx(i,j-1)/tempD
                        ! net_blow(i,j-1)   = net_blow(i,j-1)-temp3*0.5*region_qtx(i,j-1)/tempD
                    endif
                    if (region_qtx(i,j+1)>0..and. tempD.gt.0.) then
                        region_qtx(i,j+1) = region_qtx(i,j+1) - temp3*0.5*region_qtx(i,j+1)/tempD/deltim
                    endif 
    
                    if (region_qty(i-1,j)<0..and. tempD.gt.0.) then
                        temp4 = abs(temp3*0.5*region_qty(i-1,j)/tempD)
                        region_qty(i-1,j) = region_qty(i-1,j)+ temp3*0.5*region_qty(i-1,j)/tempD/deltim
                        call add_mass(temp4,wliq_soisno(:,i-1,j),wice_soisno(:,i-1,j),dz_soisno(:,i-1,j),&
                                     z_soisno(:,i-1,j),t_soisno(:,i-1,j),scv(i-1,j),region_airt(i-1,j),nl_soil,maxsnl)
                        net_blow(i-1,j)   = net_blow(i-1,j)+temp4
                    endif   
    
                    if (region_qty(i+1,j)>0..and. tempD.gt.0.) then
                        region_qty(i+1,j) = region_qty(i+1,j) - temp3*0.5*region_qty(i+1,j)/tempD/deltim
                    endif 
! debug1 = region_qsub(i,j)                    
                    region_qsub(i,j) = region_qsub(i,j)-(temp3/deltim)*region_qsub(i,j)/tempD
                    
! if (region_qsub(i,j).lt.-0.1E-8) print*,'error in line 264',region_qsub(i,j),debug1,temp3/deltim,tempD
                    if (region_qsub(i,j).lt.0.) region_qsub(i,j) = 0.
                endif
            end if
! if (region_qsub(i,j).gt.0.00001) print*,'good',region_qsub(i,j)
        end do
    end do

    do j=1,dim_x
        do i=1,dim_y  
            snl = 0
            do jj=maxsnl+1,0
                if((wliq_soisno(jj,i,j)+wice_soisno(jj,i,j))>0.) then
                   snl=snl-1
                endif
            enddo        
            if(snl<0)then
                ! Combine thin snow elements
                lb = maxsnl + 1
                call snowlayerscombine (lb,snl,&
                                z_soisno(lb:1,i,j),dz_soisno(lb:1,i,j),zi_soisno(lb-1:1),&
                                wliq_soisno(lb:1,i,j),wice_soisno(lb:1,i,j),t_soisno(lb:1,i,j),scv(i,j),snowdep(i,j))            
                ! Divide thick snow elements
                if(snl<0) &
                call snowlayersdivide (lb,snl,&
                                z_soisno(lb:0,i,j),dz_soisno(lb:0,i,j),zi_soisno(lb-1:0),&
                                wliq_soisno(lb:0,i,j),wice_soisno(lb:0,i,j),t_soisno(lb:0,i,j))
            endif
            
            snl = 0
            do jj=maxsnl+1,0
                if((wliq_soisno(jj,i,j)+wice_soisno(jj,i,j))>0.) then
                   snl=snl-1
                endif
            enddo        
            
            ! Set zero to the empty node
            if (snl > maxsnl) then
                snl = min(0,snl)
                wice_soisno(maxsnl+1:snl,i,j) = 0.
                wliq_soisno(maxsnl+1:snl,i,j) = 0.
                t_soisno   (maxsnl+1:snl,i,j) = 0.
                z_soisno   (maxsnl+1:snl,i,j) = 0.
                dz_soisno  (maxsnl+1:snl,i,j) = 0.
            endif
! scv(i,j) = sum(wliq_soisno(maxsnl+1:0,i,j) + wice_soisno(maxsnl+1:0,i,j))
            if (snl<0) then
                lb = snl+1
                scv(i,j) = sum(wliq_soisno(lb:0,i,j) + wice_soisno(lb:0,i,j))
                
                !Examining the snow conditions
                do jj=snl+1,0
                    if(wliq_soisno(jj,i,j)<0. .or. wice_soisno(jj,i,j)<0.) then
                    print*,snl,jj,wliq_soisno(jj,i,j),wice_soisno(jj,i,j),delta_swe  
                    endif
                    if (dz_soisno(jj,i,j).le.0.) &
                        print*,'wrong dz',snl,jj,&
                                dz_soisno(jj,i,j),wliq_soisno(jj,i,j),wice_soisno(jj,i,j)
                    if (z_soisno(jj,i,j)>0.) &
                        print*,'wrong z',z_soisno(jj,i,j),snl,jj,dz_soisno(jj,i,j),&
                                         wliq_soisno(jj,i,j),wice_soisno(jj,i,j)
                    if (0.001*(wice_soisno(jj,i,j)+wliq_soisno(jj,i,j))/dz_soisno(jj,i,j)>0.9) &
                        print*,'wrong density',wliq_soisno(jj,i,j),wice_soisno(jj,i,j),dz_soisno(jj,i,j)
                    if (0.001*(wice_soisno(jj,i,j)+wliq_soisno(jj,i,j))/dz_soisno(jj,i,j)<0.) &
                        print*,'wrong density',wliq_soisno(jj,i,j),wice_soisno(jj,i,j),dz_soisno(jj,i,j)
                enddo                 
            else
                scv(i,j) = 0.
            endif
            ! --------------------------------------------------------
            ! Examining error value of ice or liq value.
            ! --------------------------------------------------------

        end do
    end do
  end subroutine blowing_snow
  
  subroutine region_blowing_max(dim_y,dim_x,region_ws_x,region_ws_y,region_airt,region_rh,region_grainsize,forc_hgt_u,&
            region_qsub,region_qt_salt_x,region_qt_susp_x,region_qt_salt_y,region_qt_susp_y)
    use precision
    IMPLICIT NONE
    
    integer, intent(in) :: dim_y,dim_x
    real(r8), intent(in) :: region_ws_x(dim_y,dim_x),region_ws_y(dim_y,dim_x),forc_hgt_u(dim_y,dim_x)
    real(r8), intent(in) :: region_airt(dim_y,dim_x),region_rh(dim_y,dim_x),region_grainsize(dim_y,dim_x)
    real(r8), intent(out):: region_qsub(dim_y,dim_x),region_qt_salt_x(dim_y,dim_x),region_qt_susp_x(dim_y,dim_x)
    real(r8), intent(out):: region_qt_salt_y(dim_y,dim_x),region_qt_susp_y(dim_y,dim_x)
    real(r8) :: region_qt_salt(dim_y,dim_x),region_qt_susp(dim_y,dim_x)
    real(r8) :: region_ws(dim_y,dim_x)
    integer :: i,j

    do j=1,dim_x
        do i=1,dim_y
            region_ws(i,j)=(region_ws_x(i,j)**2 + region_ws_y(i,j)**2)**0.5
            call blowing_max(abs(region_ws(i,j)),region_airt(i,j),region_rh(i,j),region_grainsize(i,j),&
            forc_hgt_u(i,j),region_qsub(i,j),region_qt_salt(i,j),region_qt_susp(i,j))
            region_qt_salt_x(i,j)=region_qt_salt(i,j)*region_ws_x(i,j)/region_ws(i,j)
            region_qt_salt_y(i,j)=region_qt_salt(i,j)*region_ws_y(i,j)/region_ws(i,j)
            region_qt_susp_x(i,j)=region_qt_susp(i,j)*region_ws_x(i,j)/region_ws(i,j)
            region_qt_susp_y(i,j)=region_qt_susp(i,j)*region_ws_y(i,j)/region_ws(i,j)
        end do
    end do

  end subroutine region_blowing_max

! ---------------------------------------------------------------------------
! Calculating max blowing mass, not the real blowing mass 
! ---------------------------------------------------------------------------
  subroutine blowing_max(windspeed,airt,rh,grainsize,zw,qsub,qt_salt,qt_sups)
    use precision
    use blowing_varpar
    IMPLICIT NONE
    
    real(r8),intent(in) :: windspeed,airt,rh,grainsize,zw
    real(r8),intent(out):: qsub,qt_salt,qt_sups
    real(r8) :: z,zz,zmax,mmz,pi,dens_salt,dens_susp,sigmaz,zo,dens_susp_top,sup_d
    real(r8) :: ustar,ut,utstar,unstar,hstar,sigma2,svdens,diff,lamd,NUSS,M,R
    real(r8) :: deltah,alpha,rm,dmmt,uz
    integer :: i
    
    call blowing_para(windspeed,airt,rh,grainsize,zw,zo,ustar,ut,utstar,unstar,hstar,sigma2,svdens,diff,lamd,NUSS)
    qt_salt = (0.08951*utstar/ustar)*(ustar**2-unstar**2-utstar**2)
    
    zmax  = 5.
    M     = 18.01 !molecular weight of water (kg/kmole)
    R     = 8313. !universcal gas constant(J/kmol K)    
    pi    = 3.14153
    rm    = grainsize/2.
    deltah= 0.1
    qsub  = 0.
    qt_sups = 0.
    do i = 1, int(anint(zmax/deltah))
        z = i*deltah - 0.5*deltah
        alpha = 4.08+12.6*z
        mmz = (4./3.)*(pi*rm**3)*(1.+3./alpha+2./(alpha**2))
        sigmaz = sigma2*(1.02-0.027*log(z))
        dmmt = (2.*pi*rm*sigmaz)    &
                /(Liv*(Liv*M/(R*airt)-1.)/(lamd*airt*NUSS)+1./(diff*svdens*NUSS))
        if (z<hstar) then  !跃移层
            dens_salt = qt_salt/(2.8*utstar*hstar)
            qsub = qsub+(dens_salt/mmz)*dmmt*deltah
        else if (z>=hstar .and. z<= zmax) then !悬浮层
            dens_susp = 0.18*exp(-1.55*(0.05628*ustar)**(-0.544)-z**(-0.544))
            qsub = qsub+(dens_susp/mmz)*dmmt*deltah
            uz = ustar*log(z/zo)/0.4
            qt_sups = qt_sups+uz*dens_susp*deltah
        end if        
    end do

    dens_susp_top = 0.18*exp(-1.55*(0.05628*ustar)**(-0.544)-zmax**(-0.544))
    if (dens_susp_top>0.) then !认为超过5m的吹雪都会在空中被升华掉(pomeroy,1993)
        do i = 1, int(anint((15.-zmax)/0.2))
            zz = zmax+0.2*i-0.1
            sup_d = (0.18*exp(-1.55*(0.05628*ustar)**(-0.544)-zz**(-0.544)))
            if (sup_d > 0.) then 
                qsub = qsub + 0.2*sup_d
            end if
        end do
    end if
    
    if (abs(windspeed)<= 0.5*abs(ut)) then
        qsub   =0.
        qt_salt=0.
        qt_sups=0.
    endif
if (qsub.lt.0.) print*,'error in line 427'

  end subroutine blowing_max

! --------------------------------------------------------------
! Calculating real blowing snow at spatical grids
! --------------------------------------------------------------  
  subroutine grid_blowing(grid_ws,grid_qt_max,grid_size,fetch, & 
                           grid_qt)
    use precision
    IMPLICIT NONE
    
    real(r8), intent(in) :: grid_ws(3),grid_qt_max(3),grid_size
    real(r8), intent(in) :: fetch
    real(r8), intent(out):: grid_qt
    real(r8) :: windspeed,left_ws,right_ws
    real(r8) :: qtmax(3),gsi,spatial_step,qtx,x1,upper_qtmax
    integer :: i
    
    !在一维方向上求解
    windspeed = grid_ws(2) !the central wind speed
    qtmax = abs(grid_qt_max)
    spatial_step = 10. !m
    !首先需要判断上风向网格，并决定上风向吹入量    
    left_ws  = grid_ws(1)
    right_ws = grid_ws(3)
    upper_qtmax = 0.
    if (windspeed>0.) then
        if (left_ws >  0.) then
            upper_qtmax = qtmax(1)            
        else
            upper_qtmax = 0.
        end if 
    end if

    if (windspeed<0.) then
        if (right_ws < 0.) then
            upper_qtmax = qtmax(3)            
        else
            upper_qtmax = 0.
        end if 
    end if

    if (qtmax(2)<0.00001) then 
        grid_qt = 0.
    else if ((1.-upper_qtmax/qtmax(2))<=0.0001) then 
        grid_qt = qtmax(2) * grid_size !kg/m/s
    else
        x1 = -1.*fetch/3.*log(1.-upper_qtmax/qtmax(2))    
        qtx = 0.
        do i = 1,int(anint((grid_size)/spatial_step))
            gsi = x1+spatial_step*(i-0.5)
            qtx = qtx + qtmax(2)*(1.-exp(-3.*gsi/fetch))*spatial_step
        end do 
        grid_qt = qtx
    end if
   
    if (windspeed<0.) then !以向右为正
        grid_qt = -1.*abs(grid_qt)
    end if
    grid_qt = grid_qt/grid_size/grid_size
  end subroutine grid_blowing

! --------------------------------------------------------------
! Calculating real blowing snow at spatical grids
! --------------------------------------------------------------  
  subroutine grid_blowing1(grid_swe,grid_ws,grid_qt_max,grid_size,fetch, & 
                           grid_qt)
    use precision
    IMPLICIT NONE
    
    real(r8), intent(in) :: grid_swe(3),grid_ws(3),grid_qt_max(3),grid_size
    real(r8), intent(in) :: fetch
    real(r8), intent(out):: grid_qt
    real(r8) :: windspeed,left_ws,right_ws
    real(r8) :: qtmax(3),gsi,spatial_step,qtx,x1,upper_qtmax,out_qt
    integer :: i
    !以三个连续的格网为一个单元。out_qt代表下一个单元的输入
    !在一维方向上求解
    windspeed = grid_ws(2) !the central wind speed
    qtmax = abs(grid_qt_max)
    spatial_step = 10. !m
    !首先需要判断上风向网格，并决定上风向吹入量    
    left_ws  = grid_ws(1)
    right_ws = grid_ws(3)
    upper_qtmax = 0.
    if (windspeed>0.) then
        if (left_ws >  0. .and. grid_swe(1)>0.) then
            upper_qtmax = qtmax(1)            
        else
            upper_qtmax = 0.
        end if 
    end if

    if (windspeed<0.) then
        if (right_ws < 0. .and. grid_swe(3)>0.) then
            upper_qtmax = qtmax(3)            
        else
            upper_qtmax = 0.
        end if 
    end if

    if (qtmax(2)<=0.0001) then
        grid_qt = 0.
    elseif ((1.-upper_qtmax/qtmax(2))<=0.0001) then 
    !如果超过最大承载力的话，以当前网格的最大承载力来替代为输雪通量。
        grid_qt = qtmax(2) * grid_size !kg/m/s
    else
    !考虑上一网格的吹雪量，计算实际的当前网格吹雪量
        x1 = -1.*fetch/3.*log(1.-upper_qtmax/qtmax(2))    
        qtx = 0.
        do i = 1,int(anint((grid_size)/spatial_step))
            gsi = x1+spatial_step*(i-0.5)
            qtx = qtx + qtmax(2)*(1.-exp(-3.*gsi/fetch))*spatial_step
        end do 
        grid_qt = qtx
    end if
   
    if (windspeed<0.) then !以向右为正
        grid_qt = -1.*abs(grid_qt)
    end if
    grid_qt = grid_qt/grid_size/grid_size
  end subroutine grid_blowing1
  
  subroutine blowing_para(windspeed,airt,rh, grainsize, zw,&  !input
                          zo,ustar,ut,utstar,unstar,hstar,sigma2,svdens,diff,lamd,NUSS) !output
    use PhysicalConstants, only: Tfrz
    use precision
    IMPLICIT NONE

    !windspeed, 2m height
    real(r8), intent(in)::   &
             windspeed,      & !wind speed (m/s)
             zw,             &
             airt,           & !air temperature (c degree) 
             rh,             & !relative humility(0-1)
             grainsize         !snow grain size
    real(r8), intent(out)::zo,ustar,ut,utstar,unstar,hstar,sigma2,svdens,diff,lamd,NUSS
    real(r8) ::zo0,ustar0,vr,re,vair,csalt,R,M,es 
    
    M = 18.01 !molecular weight of water (kg/kmole)
    R = 8313. !universcal gas constant(J/kmol K)
    vair = 1.88E-5 !空气动力粘度m^2s-1
    csalt = 0.68    
    
    !求得摩阻风速
    zo0 = 0.5E-4
    ustar0 = windspeed*0.4/log(zw/zo0)
    !风吹雪发生情况下，z0有较大变化，需要更新如下
    zo = (0.1203*ustar0**2)/(2.*9.8)
    ustar = windspeed*0.4/log(zw/zo)!摩阻风速
    
    ut = 9.43+0.18*(airt-Tfrz)+0.0033*(airt-Tfrz)**2 !10m临界风速,需要后期加入雪龄等因素，给出一个更好的临界风速判断
    utstar = ut*0.4/log(10./zo)                  !由10m临界风速求得临界摩阻风速
    unstar = ustar*(1.-1./1.638)                 !non-erodible friction velocity
    
    hstar = 0.08436*ustar**1.27 !悬浮层下界及跃移层上界
    
    sigma2=1.-rh !underdaturation at 2m
    es = 611.15*exp(22.452*(airt-273.15)/airt) !sat pressure
    svdens = es*M/(R*airt)!sat density
    diff = 2.06E-5*(airt/273.15)**1.75 
    lamd = 0.00063*airt+0.0673
    
    vr = csalt*ustar + 2.3*utstar
    re = grainsize*vr/vair
    NUSS = 1.79+0.606*re**0.5
  end subroutine blowing_para
  
  subroutine max_blowing_tran(windspeed,airtc,max_trans)
    use precision
    IMPLICIT NONE
    
    real(r8), intent(in)  :: windspeed, airtc
    real(r8), intent(out) :: max_trans
    max_trans=(windspeed**4)*(1710.+1.36*airtc)*10E-9
  end subroutine max_blowing_tran  
  
  subroutine add_mass(delta_swe,wliq_soisno,wice_soisno,dz_soisno,z_soisno,t_soisno,&
                      scv,airt,nl_soil,maxsnl)
    use precision
    IMPLICIT NONE
    integer,  intent(in) :: maxsnl,nl_soil
    real(r8), intent(in) :: delta_swe,airt
    real(r8), intent(inout) :: t_soisno(maxsnl+1:nl_soil)
    real(r8), intent(inout) :: wliq_soisno(maxsnl+1:nl_soil), wice_soisno(maxsnl+1:nl_soil)
    real(r8), intent(inout) :: dz_soisno(maxsnl+1:nl_soil), z_soisno(maxsnl+1:nl_soil)
    real(r8), intent(inout) :: scv
    integer jj, lb,snl
    
    snl = 0
    do jj=maxsnl+1,0
        if((wliq_soisno(jj)+wice_soisno(jj))>0.) then
           snl=snl-1
           ! print*,snl,wliq_soisno(jj,i,j)+wice_soisno(jj,i,j)  
        endif
    enddo
    lb=snl+1 
    
    if (delta_swe.gt.0.) then     !accumulation  
        ! -------------------------------------------------------------------------
        ! In snow-free region, the added blowing snow is considered as a new layer
        ! -------------------------------------------------------------------------         
        if (snl.eq.0) then
            lb  = 0
            dz_soisno  (lb) = delta_swe/200. !mm to m
            z_soisno   (lb) = -0.5*dz_soisno(lb)
            wice_soisno(lb) = delta_swe
            wliq_soisno(lb) = 0.
            t_soisno   (lb) = min(airt,273.0)                    
        else
            dz_soisno(lb) = dz_soisno(lb)* &
                        (delta_swe+wice_soisno(lb)+wliq_soisno(lb))&
                        /(wice_soisno(lb)+wliq_soisno(lb))
            if (lb.ge.0) then
                z_soisno (lb) = - 0.5*dz_soisno(lb)
            else    
                z_soisno (lb) = z_soisno(lb+1) - 0.5*dz_soisno(lb)
            endif    
            wice_soisno(lb) = wice_soisno(lb) + delta_swe
        endif 
        if (snl<0) then
            lb = snl+1
            scv = sum(wliq_soisno(lb:0) + wice_soisno(lb:0))
        else
            scv = 0.
        endif        
if (z_soisno(lb)>0.) print*,'wrong z line 633',z_soisno(lb),snl,lb,dz_soisno(lb),&
                                         wliq_soisno(lb),wice_soisno(lb)              
    else 
        print*,'error in add mass module. negative value of addmass.'
    endif
  end subroutine add_mass                