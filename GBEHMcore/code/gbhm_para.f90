module gbhm_para
!-----------------------------------------------------------------------
!
! Author: Hongyi Li (lihongyi@lzb.ac.cn)
!
!-----------------------------------------------------------------------
    use precision
    implicit none
    save
! !PUBLIC DATA MEMBERS:
    public
    !_____________________
    !physical constants
    ! real(8), parameter :: Lf  = 333600.       !Latent heat of fusion
    ! character(len=12) :: output_path = 'data\output\'
      ! integer, parameter :: nrow = 157, ncol = 226    ! number of rows and cols of the whole area
      integer, parameter :: nc = 17 ,  &! nc, number of sub-catchments 
                            nx = 300,  &! nx, number of flow interval in a sub-catchment
                            np = 500    ! np, number of grids in a flow-interval      
      ! integer, parameter :: startyear = 2008, endyear = 2008
      ! integer, parameter :: iflai = 1   ! 1 for using observed lai 

 
! !Public member function
!     public :: read_basin
!   contains    
! ! -----------------------------------------------------------------
! !          read sub_catchment parameters
! ! -----------------------------------------------------------------
! !      num_flow:  number of flow interval for each sub-catchment
! !      ngrid:     number of grids in this flow interval
! !      rdx:       length of flow intervals (m)
! !      s0:        river-bed slope of the flow interval
! !      b:         river width of the flow interval (m)
! !      roughness: river roughness (Manning's coefficient) of the flow interval
! !      Dr:        river depth of the flow interval (m)
! !      grid_row:  row number of grid in this flow interval
! !      grid_col:  column number of grid in this flow interval
! ! ------------------------------------------------------------------
! subroutine read_basin(nrows,ncols,&
!                       subbasin,nsub,nflow,dx,Dr,s0,b,roughness,ngrid,grid_row,grid_col,&
!                       psubbasin,nbasinup,pbasinup, &
!                       Dr_grid)
!     implicit none
!     ! include 'dims2.inc'
    
!     integer,INTENT(in) :: ncols,nrows
!     real,INTENT(out) :: &
!                 dx       (nc,nx) ,&   ! flow interval length (m)
!                 Dr       (nc,nx) ,&   ! river depth in the flow interval (m)
!                 s0       (nc,nx) ,&   ! slope of river bed (ND)
!                 b        (nc,nx) ,&   ! width of river (m)
!                 roughness(nc,nx) ,&   ! Manning's roughness
!                 Dr_grid (nrows,ncols)
!     integer,INTENT(out) :: &
!                 nflow(nc)         ,&  ! total number of flow-intervals in a sub-basin    
!                 ngrid(nc,nx)      ,&  ! number of grids in this flow-interval
!                 grid_row(nc,nx,np),&  ! row of grids in this flow-interval
!                 grid_col(nc,nx,np),&  ! column of grids in this flow-interval
!                 nsub              ,&  ! total number of sub-basins
!                 psubbasin(nc)     ,&  ! subbasin p number      
!                 pbasinup(nc,4)    ,&  ! upbasin of p number     
!                 nbasinup(nc,4)        ! upbasin of n number  
                
!     character*6,INTENT(out) ::  subbasin(nc)      ! name of sub-basin 
! !f2py intent(in)  nrows,ncols  
! !f2py intent(out) subbasin,nsub,nflow,dx,Dr,s0,b,roughness,ngrid,grid_row,grid_col
! !f2py intent(out) psubbasin,nbasinup,pbasinup,Dr_grid                
    
!     character*200 para_dir          ! directory for parameters
!     character *200 infile 
!     character *4  ch4
!     integer       isub              ! sub-basin number
!     integer       i,j,k,ia1,ia2,ib1,ib2,ig
!     integer       iflow  ! flow-interval number
    
!     para_dir   ="../parameter/"
!     call strlen(para_dir,ib1,ib2)
!     open(1,file = para_dir(ib1:ib2)//'subbasin.dat', status='old')
!     do i = 1,nc
!        read(1,*)nsub,psubbasin(i),(pbasinup(i,j),j=1,3)!3?
!     enddo
!     close(1)

!     nbasinup=0
!     do i =1,nsub
!       write(ch4,'(i4.4)')psubbasin(i)
!       subbasin(i)='ws'//ch4
!     enddo
    
!     do i = 1,nsub
!       do k =1,4
!         if(pbasinup(i,k) >0) then
!           do j =1,nsub
!             if(psubbasin(j)== pbasinup(i,k)) then
!               nbasinup(i,k)=j              
!             endif
!           enddo
!         endif
!       enddo
!     enddo

!     Dr_grid(:,:) = 0.0
!     do isub = 1, nsub
!         infile = subbasin(isub) // '_river'
!         call strlen(infile, ia1, ia2)
!         call strlen(para_dir,ib1,ib2)
!         open(3,file=para_dir(ib1:ib2)//infile(ia1:ia2),status='old')
!         read(3,*) nflow(isub)
!         !--------------------------------------------------------------------------
!         ! An example of the '_river' file.
!         ! 9(ngrid)   1899.50(dx) 0.134942(s0) 96.25(b) 0.01000(roughness) 15.25(Dr)
!         !            15          22          15          23          15          24
!         !            15          25          15          26          15          27
!         !            15          28          16          27          16          28 
!         !--------------------------------------------------------------------------
!         do iflow = 1, nflow(isub)
!             read(3,*) ngrid    (isub,iflow), dx(isub,iflow),&
!                       s0       (isub,iflow), b (isub,iflow),&
!                       roughness(isub,iflow), Dr(isub,iflow)    
            
!             if( dx(isub,iflow).lt.0.1) then
!                 dx(isub,iflow) = 3000.0
!             endif
!             dx(isub,iflow) = dx(isub,iflow)*1.50  !why 1.50?
            
!             if(s0(isub,iflow).le.0.1E-5) then   
!                 s0(isub,iflow)=0.00001
!                 print *, 'wrong in s0',s0(isub,iflow)  
!             endif            
!             if(s0(isub,iflow).eq.-9999.0) then 
!                 s0(isub,iflow)=0.00001 
!                 print *, 'wrong in s0',s0(isub,iflow)  
!             endif            
    
!             if(ngrid(isub,iflow) .gt. np) then
!                 print *,'more than np grids:',ngrid(isub,iflow),isub,iflow,np
!             endif
    
!             read(3,*) (grid_row(isub,iflow,j),grid_col(isub,iflow,j),&
!                     j = 1, ngrid(isub,iflow))
!             do ig = 1,ngrid(isub,iflow)      
!                 Dr_grid(grid_row(isub,iflow,ig),grid_col(isub,iflow,ig)) = Dr(isub,iflow)  
!             end do
!         end do
!         close(3)
!     end do
!     print * ,'ok2, finish of ws000_river'
         
! end subroutine read_basin   


end module gbhm_para
