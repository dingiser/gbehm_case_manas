
! --------------------------------------------------------
! MODULE NANE:
!     time manager model
!
! PURPOSE :
!     To provide some basic operations for time stamp
!
! Initial author: Hua Yuan, /04/2014/
! Revised :       Hongyi LI, 10/2015
! --------------------------------------------------------

MODULE timemanager

   use precision
   implicit none

   logical, save :: isgreenwich

CONTAINS

   SUBROUTINE adj2begin(idate)

      implicit none
      integer, intent(inout) :: idate(3)

      if (idate(3) == 86400) then
         idate(3) = 0
         idate(2) = idate(2) + 1
         if (MOD(idate(1),4)==0 .AND. idate(2)==367) then
            idate(1) = idate(1) + 1; idate(2) = 1
         end if
         if ( MOD(idate(1),4)/=0 .AND. idate(2)==366) then
            idate(1) = idate(1) + 1; idate(2) = 1
         end if
      end if

   END SUBROUTINE adj2begin

   SUBROUTINE adj2end(idate)

      implicit none
      integer, intent(inout) :: idate(3)

      if (idate(3) == 0) then
         idate(3) = 86400
         idate(2) = idate(2) - 1
         if (idate(2) == 0) then
            idate(1) = idate(1) - 1
            if ( MOD(idate(1),4)==0) then
               idate(2) = 366
            else
               idate(2) = 365
            end if
         end if
      end if

   END SUBROUTINE adj2end

   SUBROUTINE localtime2gmt(idate, longi)

      implicit none
      integer, intent(inout) :: idate(3)
      real(r8),intent(in)    :: longi

      integer  maxday
      real(r8) tdiff

      tdiff = longi/15.*3600.
      idate(3) = idate(3) - int(tdiff)

      if (idate(3) < 0) then

         idate(3) = 86400 + idate(3)
         idate(2) = idate(2) - 1

         if (idate(2) < 1) then
            idate(1) = idate(1) - 1
            if ( MOD(idate(1),4)==0 ) then
               idate(2) = 366
            else
               idate(2) = 365
            endif
         endif
      endif

      if (idate(3) > 86400) then

         idate(3) = idate(3) - 86400
         idate(2) = idate(2) + 1

         if ( MOD(idate(1),4)==0 ) then
            maxday = 366
         else
            maxday = 365
         endif

         if(idate(2) > maxday) then
            idate(1) = idate(1) + 1
            idate(2) = 1
         endif
      endif

   END SUBROUTINE localtime2gmt

   SUBROUTINE ticktime(deltim, idate)

      implicit none

      real(r8),INTENT(in)    :: deltim
      integer, INTENT(inout) :: idate(3)
      integer maxday

      idate(3) = idate(3) + nint(deltim)
      if (idate(3) > 86400) then

         idate(3) = idate(3) - 86400
         idate(2) = idate(2) + 1

         if ( MOD(idate(1),4)==0) then
            maxday = 366
         else
            maxday = 365
         endif

         if(idate(2) > maxday) then
            idate(1) = idate(1) + 1
            idate(2) = 1
         endif
      endif

   END SUBROUTINE ticktime

   SUBROUTINE calendarday_func(date, longi,calendarday)

      implicit none
      integer, intent(in) :: date(3)
      real(r8),optional   :: longi
      real(r8),intent(out):: calendarday
      integer idate(3)
      real(r8) longitude

      idate(:) = date(:)

      longitude = longi
      ! transfer local time to GMT time.
      call localtime2gmt(idate, longitude)

      calendarday = float(idate(2)) + float(idate(3))/86400.

   END SUBROUTINE calendarday_func

   ! real(r8) FUNCTION calendarday_stamp(stamp, long)

      ! implicit none
      ! type(timestamp), intent(in) :: stamp
      ! real(r8),        optional   :: long

      ! integer idate(3)
      ! real(r8) longitude

      ! idate(1) = stamp%year
      ! idate(2) = stamp%day
      ! idate(3) = stamp%sec

      ! if (.NOT. present(long)) then
         ! longitude = long
      ! else
         ! longitude = long
      ! end if

      ! if ( .not. isgreenwich ) then
         ! call localtime2gmt(idate, longitude)
      ! end if

      ! calendarday_stamp = float(idate(2)) + float(idate(3))/86400.
      ! return

   ! END FUNCTION calendarday_stamp

END MODULE timemanager
