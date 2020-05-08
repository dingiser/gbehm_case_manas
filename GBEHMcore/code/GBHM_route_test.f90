! subroutine GBHM_route(runoff_inter,runoff_g,runoff,slp,length,area)
    ! implicit none
    ! include 'dims2.inc'

    ! real, INTENT(in) :: &        
      ! runoff_inter(nrow,ncol) , &      !Subsurface flow In each grid
      ! runoff_g    (nrow,ncol) , &      !exchanges between groundwater and river In each grid
      ! slp         (nrow,ncol) , &
      ! length      (nrow,ncol) , &      ! average hillslope length (m)
      ! area        (nrow,ncol)          ! area of the local grid (m**2)
    ! real, INTENT(inout) ::  &
      ! runoff      (nrow,ncol)      


      
    ! character*6 subbasin(nc)   ! name of sub-basin 
    ! integer isub               ! sub-basin number
    ! integer iflow              ! flow-interval number
    ! integer nflow(nc)          ! total number of flow-intervals in a sub-basin
    ! integer psubbasin(nc)      ! subbasin p number 
    ! integer pbasinup(nc,4)     ! upbasin of p number 
    ! integer nbasinup(nc,4)     ! upbasin of n number 
                               
    ! real    dx(nc,nx)          ! flow interval length (m)
    ! real    Dr(nc,nx)          ! river depth in the flow interval (m)
    ! real    s0(nc,nx)          ! slope of river bed (ND)
    ! real    b(nc,nx)           ! width of river (m)
    ! real    roughness(nc,nx)   ! Manning's roughness

    ! integer ngrid(nc,nx)       ! number of grids in this flow-interval
    ! integer grid_row(nc,nx,np) ! row of grids in this flow-interval
    ! integer grid_col(nc,nx,np) ! column of grids in this flow-interval

! c       data and time variables
! c-------------------------------------------------------------------
    ! real    dt             ! time step (second)
    ! integer year, year_tmp ! calendar year   (19XX)
    ! integer month          ! month in a year (1-12)
    ! integer day            ! day in a month  (1-30)
    ! integer hour           ! hour in a day   (1-24)
    ! integer dayinmonth(12) ! days of a month
    ! integer hydroyear      ! year for hydro simulation
    ! integer startmonth     ! start month in the hydro-year
    ! integer endmonth       ! end   month in the hydro-year
    ! integer startday       ! start day in the first month of the hydro-year
    ! integer endday         ! end   day in the last  month of the hydro-year
    ! integer im             !,iy, id, ih ! for year, month, day and hour loop
    ! integer idc            ! contineous day in a year (1-366)
    ! integer ihc            ! contineous hour in a year (1-366*24)

    ! integer start, finish  ! 0 - faulse,    1 - true
    ! integer start_sub, end_sub,nsub  ! for partial simulation of a basin

! c     Other variables
! c -------------------------------------------------------------------     
    ! character*2    ch2, a2           ! a 2-character variable
    ! character*4    ch4               ! a 4-character variable
    ! character*200  atmp
    ! real           tmp, tmp1, tmp2   ! , value, mean  ! temporary variables
    ! integer        i, j, k           ! temporary variables
    ! integer        ir, ic, idd       !, iyy, imm
    ! integer        ia1, ia2, ib1, ib2, ic1, ic2,ir2
    ! integer        itmp ,tmp_num
    ! real           tmp_rain,tmp_tmin,tmp_tmax,tmp_evap
    ! real           ran0    
      
    ! common /date1  / year, month, day, hour, idc,ihc
    ! common /date2  / hydroyear,startmonth,startday,endmonth,endday
    ! common /date3  / dayinmonth
    ! common /river1    /  nsub,subbasin
    ! common /river2    /  nflow, dx, Dr
    ! common /river3    /  s0, b, roughness
    ! common /river4    /  psubbasin,nbasinup
    ! common /simulation/  start, finish, dt, start_sub, end_sub    
    ! common /grid_attrib/ ngrid,grid_row,grid_col

! c     Initialize the model parameters. 
! c     Most of them are global parameters passed through common blocks.
! c-----------------------------------------------------------------------
    ! data dayinmonth /31,28,31,30,31,30,31,31,30,31,30,31/
    ! data startmonth /1 /
    ! data startday   /1 /
    ! data endmonth   /12 /
    ! data endday     /31/ 
! c    To define partial simulation of this basin
! c----------------------------------------------------------------------
    ! start_sub = 1
    ! end_sub   = nc      ! total is 251                          !!!change 1            
    ! dt        = 3600.   ! time step of hydrological simulation,unit: second
      
! c-----------------------------------------------------------------------
! c           0、main code                  
! c-----------------------------------------------------------------------
   
    ! call read_pbasin()
    ! call read_catchmentpara()

    ! print * ,'model initiation'
    ! start  = 1
    ! finish = 0
    ! do isub = start_sub, end_sub
      ! call hillslope_model(isub,runoff_inter,runoff_g,runoff,slp,length,area)
      ! call river_routing(isub) 
    ! end do

    ! do 888 hydroyear = startyear, endyear
      ! start  = 0
      ! finish = 0
      ! year = hydroyear
      ! write(ch4, '(i4.4)') year
      ! idc=0
      ! ihc=0
     
      ! do 777 im = startmonth, endmonth
        ! month = im             
        ! if(mod(year,4).eq.0 .and. month.eq.2) dayinmonth(month)=29
        ! if(mod(year,4).ne.0 .and. month.eq.2) dayinmonth(month)=28
        
        ! do 666 day = 1, dayinmonth(month)
           ! idc = idc + 1          
            ! print *,'running at: ', year, month,day
            ! do 555 hour = 1,24 
               ! ihc = ihc + 1
                    ! do isub = start_sub,end_sub
                        ! call hillslope_model(isub,runoff_inter,runoff_g,runoff,slp,length,area)
                        ! call river_routing(isub)
                    ! enddo
! 555         continue   
! 666        continue
! 777      continue
! 888   continue

! end


subroutine strlen(str, l1,l2)
    character str*200
    integer i,l1,l2,k
    k=0
    do i = 1, 200
      if(k.eq.0 .and. str(i:i).NE.' ') then
        l1=i
        k=1
      elseif(k.eq.1 .and. str(i:i).EQ.' ') then
        l2 = i-1
        return
      endif
    end do
    l2 = i
    return
end
! -----------------------------------------------------------------
!          2. read sub_catchment parameters
! -----------------------------------------------------------------
!      num_flow:  number of flow interval for each sub-catchment
!      ngrid:     number of grids in this flow interval
!      rdx:       length of flow intervals (m)
!      s0:        river-bed slope of the flow interval
!      b:         river width of the flow interval (m)
!      roughness: river roughness (Manning's coefficient) of the flow interval
!      Dr:        river depth of the flow interval (m)
!      grid_row:  row number of grid in this flow interval
!      grid_col:  column number of grid in this flow interval
! ------------------------------------------------------------------
subroutine read_basin(nrows,ncols,&
                      subbasin,nsub,nflow,dx,Dr,s0,b,roughness,ngrid,grid_row,grid_col,&
                      psubbasin,nbasinup,pbasinup, &
                      Dr_grid)
    implicit none
    include 'dims2.inc'
    
    integer,INTENT(in) :: ncols,nrows
    real,INTENT(out) :: &
                dx       (nc,nx) ,&   ! flow interval length (m)
                Dr       (nc,nx) ,&   ! river depth in the flow interval (m)
                s0       (nc,nx) ,&   ! slope of river bed (ND)
                b        (nc,nx) ,&   ! width of river (m)
                roughness(nc,nx) ,&   ! Manning's roughness
                Dr_grid (nrows,ncols)
    integer,INTENT(out) :: &
                nflow(nc)         ,&  ! total number of flow-intervals in a sub-basin    
                ngrid(nc,nx)      ,&  ! number of grids in this flow-interval
                grid_row(nc,nx,np),&  ! row of grids in this flow-interval
                grid_col(nc,nx,np),&  ! column of grids in this flow-interval
                nsub              ,&  ! total number of sub-basins
                psubbasin(nc)     ,&  ! subbasin p number      
                pbasinup(nc,4)    ,&  ! upbasin of p number     
                nbasinup(nc,4)        ! upbasin of n number  
                
    character*6,INTENT(out) ::  subbasin(nc)      ! name of sub-basin 
!f2py intent(in)  nrows,ncols  
!f2py intent(out) subbasin,nsub,nflow,dx,Dr,s0,b,roughness,ngrid,grid_row,grid_col
!f2py intent(out) psubbasin,nbasinup,pbasinup,Dr_grid                
    
    character*200 para_dir          ! directory for parameters
    character *200 infile 
    character *4  ch4
    integer       isub              ! sub-basin number
    integer       i,j,k,ia1,ia2,ib1,ib2,ig
    integer       iflow  ! flow-interval number
    
    para_dir   ="../parameter/"
    call strlen(para_dir,ib1,ib2)
    open(1,file = para_dir(ib1:ib2)//'subbasin.dat', status='old')
    do i = 1,nc
       read(1,*)nsub,psubbasin(i),(pbasinup(i,j),j=1,3)!3?
    enddo
    close(1)

    nbasinup=0
    do i =1,nsub
      write(ch4,'(i4.4)')psubbasin(i)
      subbasin(i)='ws'//ch4
    enddo
    
    do i = 1,nsub
      do k =1,4
        if(pbasinup(i,k) >0) then
          do j =1,nsub
            if(psubbasin(j)== pbasinup(i,k)) then
              nbasinup(i,k)=j              
            endif
          enddo
        endif
      enddo
    enddo

    Dr_grid(:,:) = 0.0
    do isub = 1, nsub
        infile = subbasin(isub) // '_river'
        call strlen(infile, ia1, ia2)
        call strlen(para_dir,ib1,ib2)
        open(3,file=para_dir(ib1:ib2)//infile(ia1:ia2),status='old')
        read(3,*) nflow(isub)
        !--------------------------------------------------------------------------
        ! An example of the '_river' file.
        ! 9(ngrid)   1899.50(dx) 0.134942(s0) 96.25(b) 0.01000(roughness) 15.25(Dr)
        !            15          22          15          23          15          24
        !            15          25          15          26          15          27
        !            15          28          16          27          16          28 
        !--------------------------------------------------------------------------
        do iflow = 1, nflow(isub)
            read(3,*) ngrid    (isub,iflow), dx(isub,iflow),&
                      s0       (isub,iflow), b (isub,iflow),&
                      roughness(isub,iflow), Dr(isub,iflow)    
            
            if( dx(isub,iflow).lt.0.1) then
                dx(isub,iflow) = 3000.0
            endif
            dx(isub,iflow) = dx(isub,iflow)*1.50  !why 1.50?
            
            if(s0(isub,iflow).le.0.1E-5) then   
                s0(isub,iflow)=0.00001
                print *, 'wrong in s0',s0(isub,iflow)  
            endif            
            if(s0(isub,iflow).eq.-9999.0) then 
                s0(isub,iflow)=0.00001 
                print *, 'wrong in s0',s0(isub,iflow)  
            endif            
    
            if(ngrid(isub,iflow) .gt. np) then
                print *,'more than np grids:',ngrid(isub,iflow),isub,iflow,np
            endif
    
            read(3,*) (grid_row(isub,iflow,j),grid_col(isub,iflow,j),&
                    j = 1, ngrid(isub,iflow))
            do ig = 1,ngrid(isub,iflow)      
                Dr_grid(grid_row(isub,iflow,ig),grid_col(isub,iflow,ig)) = Dr(isub,iflow)  
            end do
        end do
        close(3)
    end do
    print * ,'ok2, finish of ws000_river'

end

! -----------------------------------------------------------------
!          3. A physically-based hillslope hydrological model
! -----------------------------------------------------------------
subroutine hillslope_route_model(subbasin,isub,nflow,psubbasin,nbasinup,nsub,ngrid,grid_row,grid_col,start, &
year,month,day,hour,dayinmonth,startmonth,endmonth,startday,endday, dt, &
idc,ihc,start_sub,end_sub,runoff_inter,runoff_g,runoff,slp,length,area,dx,Dr,s0,b,roughness,&
q1,q2,qr1,qr2,qlin1,Qh,Qd,Qm,Drw,sst)
    implicit none
    include 'dims2.inc'
    
    character*6 , INTENT(in) ::  subbasin(nc)      ! name of sub-basin     
    real, INTENT(in) :: &        
      runoff_inter(nrow,ncol) , &   ! Subsurface flow In each grid
      runoff_g    (nrow,ncol) , &   ! exchanges between groundwater and river In each grid
      slp         (nrow,ncol) , &
      length      (nrow,ncol) , &   ! average hillslope length (m)
      area        (nrow,ncol) , &   ! area of the local grid (m2)
      dx          (nc,nx)     , &   ! flow interval length (m)
      Dr          (nc,nx)     , &   ! river depth in the flow interval (m)
      b           (nc,nx)     , &   ! width of river (m)
      roughness   (nc,nx)     , &   ! Manning's roughness
      dt                            ! time step (hour unit in route subroutine?)      
    integer, INTENT(in) :: & 
      isub               , &        ! sub-basin number
      nflow    (nc)      , &        ! total number of flow-intervals in a sub-basin
      psubbasin(nc)      , & 
      nbasinup (nc,4)    , &
      nsub               , &        ! total number of sub-basins
      ngrid   (nc,nx)    , &        ! number of grids in this flow-interval
      grid_row(nc,nx,np) , &        ! row of grids in this flow-interval
      grid_col(nc,nx,np) , &        ! column of grids in this flow-interval
      start,               &        ! 0 - faulse,    1 - true
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
    real, INTENT(in) :: &
      s0          (nc,nx), &        ! slope of river bed (ND)
      runoff      (nrow,ncol) 
      

!      Discharge variables
! --------------------------------------------------------------------------
real, INTENT(out)   ::  q2   (nc,nx)       ! discharge of current time step (m^3/s)
real, INTENT(out)   ::  Qh   (nc,8800)     ! hourly mean discharge (m3/s)
real, INTENT(inout) ::  qr1  (nc,nx)       ! discharge of reservoir flowout last time step(m^3/s)
real, INTENT(inout) ::  qlin1(nc,nx)       ! lateral inflow of last time step (m^3/m/s)
real, INTENT(inout) ::  Qd   (nc,366)      ! daily average discharge (m^3/s)
real, INTENT(inout) ::  Qm   (nc,12)       ! monthly mean discharge (m3/s)
real, INTENT(inout) ::  Drw  (nc,nx)       ! water depth in the river of the flow interval (m)
real, INTENT(inout) ::  sst  (nrow,ncol) 
real, INTENT(inout) ::  q1   (nc,nx)       ! discharge of last time step (m^3/s)
real, INTENT(inout) ::  qr2  (nc,nx)       ! discharge of reservoir flowout current time step(m^3/s)

    ! save q2       ! discharge of current time step (m^3/s)
    ! save qr1      ! discharge of reservoir flowout last time step(m^3/s)
    ! save qlin1    ! lateral inflow of last time step (m^3/m/s)
    ! save qlin2    ! lateral inflow of current time step (m^3/m/s)
    ! save Qh       ! hourly mean discharge (m3/s)
    ! save Qd       ! daily average discharge (m^3/s)
    ! save Qm       ! monthly mean discharge (m3/s)

!f2py intent(in) subbasin,isub,nflow,psubbasin,nbasinup,nsub,ngrid,grid_row,grid_col,start
!f2py intent(in) year,month,day,hour,dayinmonth,startmonth,endmonth,startday,endday, dt
!f2py intent(in) idc,ihc,start_sub,end_sub,nflow,runoff_inter,runoff_g
!f2py intent(in) slp,length,area,dx,Dr,b,roughness    
!f2py intent(in) runoff,s0
!f2py intent(in,out) q1,qr1,qr2,qlin1,Qd,Qm,Drw,sst
!f2py intent(out) q2,Qh
    real    qlin2(nc,nx)! lateral inflow of current time step (m^3/m/s)
    real    qin(nx)     ! lateral inflow of one flow interval(m3/sec/m),unit (m3/s/m, flow into river)
    integer iflow       ! flow-interval number
    real    q_hillslope ! surface runoff of one simulation 
    real    qin_tmp
    real    qin_inter   !Subsurface flow In each grid
    real    qin_g       !exchanges between groundwater and river In each grid
    real    qin_sfc  
    real    q_surf 
    real    water_depth
    real    surface_n
    real    waterhead
    real    power    
    ! real    Drw(nc,nx)         ! water depth in the river of the flow interval (m)
    integer ir,ic,ig,iy, im, id  !,ih ! for year, month, day and hour loop


    integer l1, l2, l3
    integer level
    integer p1(9), p2(9), p3(9)

! c     Other variables
! c----------------------------------------------------------------------

    character*200  para_dir    ! directory for parameters
    character*200  dem_dir    ! directory for dem
    character*200  data_dir    ! directory for input data
    character*200  result1_dir ! directory for storing simulation result
    character*200  result2_dir ! directory for storing simulation result
    character*200  simul_dir   ! directory for storing temporal result

    character*4    ch4         ! a 4-character variable

    real    tmp, value         ! temporary variables
    integer i, j, k, m         !, n      ! temporary variables
    integer tmpj,tmpim,tmpid
    integer ia1, ia2, ib1, ib2, ic1, ic2
    real    criterion, h1, h2, f, df
    real    q1_upper, q2_upper,rs    !, i_sub,  lo_dis
    integer ii1, ii2           ! , ii3
    real    Qtotal,basinarea
    real    b_tmp              ! adjust riverbed-width during heavy-flood
    real    rs_tmp             ! adjust rs during heavy-flood
    integer inicon
! 3.1  model initiation
      inicon=1
      if (start .eq. 1) then
        inicon=0 
        goto 789
      endif
! -------------------------------------
!       3.2  start simulation
! -------------------------------------
      do 999 iflow=1,nflow(isub)
         qin(iflow)=  0.0            ! total lateral inflow (m3/s/m)
         qin_tmp   =  0.0            ! lateral inflow from present flow-interval
        do 666 ig=1,ngrid(isub,iflow)
         ir=grid_row(isub,iflow,ig)
         ic=grid_col(isub,iflow,ig)

         qin_inter =  0.0     !Subsurface flow In each grid
         qin_g     =  0.0     !exchanges between groundwater and river In each grid
         qin_sfc   =  0.0     !surface runoff from a flow interval In each grid
         
         qin_inter = runoff_inter(ir,ic)*length(ir,ic)          !to m/s !change2 (输入COLM中像元径流)
         qin_g     = runoff_g(ir,ic)*length(ir,ic)                !m/s
         
         water_depth = sst(ir,ic)+runoff(ir,ic)
         q_surf    =amax1(0.0,  water_depth)          !m  (surface runoff)
         q_hillslope=0.0                              !m^3/m

         if( q_surf .le. 0.1E-8 ) goto 500                    
           water_depth = q_surf
           ! runoff(ir,ic)=runoff(ir,ic) - water_depth

           surface_n = 0.02      !calibrate 0.1-0.2
           waterhead = slp(ir,ic) 
           power     = 1.6667
           
          q_hillslope = dt * sqrt(waterhead)* &
                          water_depth**power/surface_n     ! m3/m, one hillslope
          if(q_hillslope .le. 0.1E-20) q_hillslope = 0.0

          q_hillslope = amin1(q_hillslope, &
                                water_depth*length(ir,ic))
          water_depth = water_depth - q_hillslope/length(ir,ic) 
          ! runoff(ir,ic) =runoff(ir,ic)+water_depth
          sst(ir,ic)    = water_depth
          water_depth   = 0.0
          qin_sfc       = q_hillslope/dt
qin_sfc =  runoff(ir,ic)* length(ir,ic)  
     
! cc500        qin_tmp = qin_tmp +(qin_sfc
! cc     :           + qin_inter + qin_g)* area(ir,ic)/length(ir,ic)        ! m3/s单宽流量, one flow-interval
500        qin_tmp = qin_tmp +(qin_sfc &
                 + qin_inter + qin_g)* 1000.*1000./length(ir,ic)      !why 1000?

666   continue
!/dx(isub,iflow)???
      qin(iflow) = qin_tmp!/dx(isub,iflow) ! m3/s, total lateral inflow of one flow-interval
999   continue
      ! return
! end


! ------------------------------------------------------------------------------
!          4.*** River Routing Model ***   (Kinematic Wave method)                                  
! ------------------------------------------------------------------------------
! subroutine river_routing(isub)
    ! implicit none
    ! include 'dims2.inc'

! BEGIN initialize model
! ----------------------------------------------------------------     
789 simul_dir  = "../simulation/"
    result1_dir= "../result_river/"
    if(month.eq.startmonth .and.  day.eq.startday .and. hour.eq.1) then
      do i = 1, nsub
        do j= 1, 366
          Qd(i,j)=0.0
        end do
      end do    
      do i = 1, nsub
          do j= 1, 12
              Qm(i,j)=0.0
          end do
      end do
    endif

    if(start .eq. 1) then
      do i = 1, nsub
        do j = 1, nflow(isub)
          qlin1(i,j) = 0.0
        end do
      end do

      if(inicon.eq.1) then
        call strlen(simul_dir, ia1,ia2)
        open(3,file = simul_dir(ia1:ia2)//subbasin(isub)//'I_flow2', status='old')
        read (3,*)(i,q1(isub,j), j=1, nflow(isub))
        close(3)
      else
        q1(isub,1)=0.5
        do iflow = 2, nflow(isub)
          q1(isub,iflow) = q1(isub,iflow-1)+0.4
        end do
      endif
      !calculate initial river water depth
      do iflow = 1, nflow(isub)
         criterion = 0.01
         k = 1
         h1  = q1(isub,iflow) / b(isub,iflow)
 5       tmp = roughness(isub,iflow) * q1(isub,iflow)/ sqrt(s0(isub,iflow))
         f   = b(isub,iflow)*h1 -  (tmp**0.6)*((b(isub,iflow)+2.0*h1)**0.4)
         if(k.gt.1 .and. abs(f) .lt. criterion) goto 8
         df = b(isub,iflow)-0.8*((tmp/(b(isub,iflow)+2.0*h1))**0.6)
         h2 = h1-f/df
         h1 = h2
         if(k.ge.10) goto 8
         k  = k+1
         goto 5
 8       Drw(isub,iflow) = h2        
      end do

      do iflow = 1, nflow(isub)
        qr1(isub,iflow) = q1(isub,iflow)
      end do
    endif 

    if(start .eq. 1) return
    
! END initialize model
! ---------------------------------------------------------------------
    if(month.eq.startmonth.and.day.eq.startday.and.ihc.eq.1) then
        do iflow = 1, nflow(isub)
           qlin1(isub,iflow)= qin(iflow)/dx(isub,iflow)
        end do
    end if

      do iflow = 1, nflow(isub)
        qlin2(isub,iflow) = qin(iflow)/dx(isub,iflow)
      end do

      m = 1!=1???????????
      do 555 j = 1, m                  !time loop
        do 333 iflow = 1, nflow(isub)  !river segment loop
!   define the river network
      q1_upper=0.
      q2_upper=0.
      if(iflow .eq. 1) then  
          if(psubbasin(isub) > 1000 .and. psubbasin(isub) < 2000) then
            q1_upper=0.
            q2_upper=0.
          else
              do ii1 =1,4
                if(nbasinup(isub,ii1) > 0) then
                    ii2=nbasinup(isub,ii1)
                    q1_upper=q1_upper+qr1(ii2,nflow(ii2))
                    q2_upper=q2_upper+qr2(ii2,nflow(ii2)) !where initial qr2?
                endif
            enddo
          endif
      else
        q1_upper=qr1(isub,iflow-1)
        q2_upper=qr2(isub,iflow-1)
      endif                         
!   end of river network definition
      rs=roughness(isub,iflow)

      b_tmp = b(isub,iflow)
      rs_tmp = rs
      if(Drw(isub,iflow) .gt. 0.5*Dr(isub,iflow)) then        
        b_tmp=1.1*b(isub,iflow)
        rs_tmp = 1.5*rs
      endif
      if(Drw(isub,iflow) .gt. 1.0*Dr(isub,iflow)) then
        b_tmp=1.5*b(isub,iflow)
        rs_tmp = 3.0*rs
      endif

      call nkws(dt,dx(isub,iflow),b_tmp &
              ,s0(isub,iflow),rs_tmp, qlin1(isub,iflow),&
              qlin2(isub,iflow),q1_upper,q1(isub,iflow),&
              q2_upper,q2(isub,iflow),Drw(isub,iflow))

      qr2(isub,iflow)=q2(isub,iflow)

      if(isub.eq.1 .and. iflow.eq.1 .and. q2_upper.gt.0.0) then     
        print *,"out nkws", dt,dx(isub,iflow),b(isub,iflow),             &
          s0(isub,iflow),rs,qlin1(isub,iflow),qlin2(isub,iflow),q1_upper &
         ,q1(isub,iflow),q2_upper,q2(isub,iflow),Drw(isub,iflow),idc,ihc
        ! stop
      endif

 333  continue   
 
      do iflow = 1, nflow(isub)
         qlin1(isub,iflow) = qlin2(isub,iflow)
         q1   (isub,iflow) = q2   (isub,iflow)
         qr1  (isub,iflow) = qr2  (isub,iflow)
      end do

 555  continue

! discharge output
! ---------------------------------------------------------------------
    Qh(isub,ihc)  = qr2(isub,nflow(isub))
    Qd(isub,idc)  = Qd(isub,idc) + qr2(isub,nflow(isub))/24.0
    Qm(isub,month)= Qm(isub,month) + qr2(isub,nflow(isub))/ (float(dayinmonth(month))*24.0)
! print*,'qr2=',qr2
    if(month.eq.endmonth .and. day.eq.endday .and. hour.eq.24) then
        call strlen(simul_dir, ia1,ia2)
        open(9,file = simul_dir(ia1:ia2)//subbasin(isub)//'I_flow2', status='replace')
            do iflow = 1, nflow(isub)
                write(9,110) iflow, q2(isub,iflow)
 110         format(2x,i4,f20.3)
            end do
        close(9)

        call strlen(result1_dir,ib1,ib2)
        if(year.eq.startyear ) then
            open(9,file = result1_dir(ib1:ib2)//subbasin(isub)//'.daily', status='unknown')
        else
            open(9,file = result1_dir(ib1:ib2)//subbasin(isub)//'.daily', access='append', status='old')
        endif
        
        !................................................................
        ! The following codes is to make sure a correct order of Qm data,
        ! in case startday or startmonth is not equal to 1
        tmpj = 0
        if (startmonth.gt.1 .or. startday.gt.1) then
            do tmpim = 1,startmonth
                do tmpid = 1, dayinmonth(tmpim)
                    tmpj = tmpj+1
                enddo
            enddo
            tmpj = tmpj-dayinmonth(tmpim)+startday-1
        endif
        !.................................................................
        
        j = 0
        do im = startmonth, endmonth
            do id = 1, dayinmonth(im)
                j = j+1
! print*,'Qd1=',Qd(isub,j)
                write(9, '(3i6,f16.3)') year,im,id, Qd(isub,j+tmpj)!
            end do
        end do
            close(9)
    endif    
        
    return
end subroutine hillslope_route_model

subroutine nkws(dt, dx, b, s0, roughness,qlin0, qlin, Q01,Q02,Q1, Q2, y)
!**********************************************************************
!                                                                     *
!                 ***** Kinematic Wave Model *****                    *
!             (Nonlinear Scheme using Newton's method)                *
!                                                                     *
!**********************************************************************
!
!	Definition of Variales:
!
!            (time)
!              ^ 
!              |
!	          Q1(known)       Q2(unknown)
!
!               Q01(known)      Q02(known)    --> (distance)
!
!	Q01, Q02 - known, discharge of last time step
!	Q1       - known, discharge of current time step
! Q2       - unknown, discharge of current time step
!	qlin0    - known, lateral inflow of last time step
!	qlin     - known, lateral inflow of current time step
!
!	y  -- water depth
!	p  -- wetted perimeter
!	b  -- river width
!	s0 -- river slope
!
!**********************************************************************
    real Q01, Q02, Q1, Q2
    real qlin0, qlin
    real dt
    real dx, b, s0, roughness
    real y, p, beta, criterion

    beta = 0.6
    criterion = 0.0001
    k   = 1
    h1  = Q02/b
    
 15 tmp = roughness*Q02/sqrt(s0)
 
    f   = b*h1 - (tmp**0.6)*((b+2.0*h1)**0.4)
    if(k.gt.1 .and. abs(f) .lt. criterion) goto 18
    df = b - 0.8*((tmp/(b+2.0*h1))**0.6)
    h2 = h1 - f/df
    h1 = h2
    if(k.ge.30) goto 18
    k  = k+1
    goto 15
    
 18 y = h2
 
    p = b+2.0*y  ! for rectangular channel
    alfa = (roughness*p**(2.0/3.0)/sqrt(s0))**0.6
    if((Q02+Q1).le.0.0) then
      cc = 0
    else
      cc = (0.5*(Q02+Q1))**(beta-1.0)
    endif
    aa   = dt*Q1/dx + alfa*beta*Q02*cc + 0.5*(qlin0+qlin)*dt
    bb   = dt/dx + alfa*beta*cc
    qq1  = aa/bb

    if(qq1 .le. 0.1e-5)  qq1=0.1e-5 
   
    ! Using Newton's method to calculate discharge
    ctmp = dt*Q1/dx + alfa*Q02**beta + 0.5*dt*(qlin+qlin0)
    k= 1
    
 20 f= dt*qq1/dx + alfa*qq1**beta-ctmp
 
    if((k.gt.1).and.(abs(f).le.criterion)) goto 30
    df  = dt/dx + alfa*beta*qq1**(beta-1.0)
    qq2 = qq1 - f/df     
!print *,k,ctmp,f,df,alfa ,qq1,qq2

    if(qq2 .le. 0.1e-5) then
      qq2=0.0
      goto 30
    endif
    
    qq1 = qq2
    if(k .ge. 30) goto 30
    k = k+1
    goto 20
    
 30 Q2 = qq2
 
    k   = 1
    h1  = Q2/b
    
 45 tmp = roughness*Q2/sqrt(s0)
 
    f   = b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
    if(k.gt.1 .and. abs(f) .lt. criterion) goto 50
    df = b-0.8*((tmp/(b+2.0*h1))**0.6)
    h2 = h1-f/df
    h1 = h2
    if(k.ge.30) goto 50
    k = k+1
    goto 45
    
 50 y = h2
 
  ! print *,"nkws q=" ,q2,"  drw=",y
    return
end subroutine nkws