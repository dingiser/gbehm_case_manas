module precision
!-------------------------------------------------------------------------------
! Purpose:
!	Define the precision to use for floating point and integer operations
!	throughout the model.
!-------------------------------------------------------------------------------
  integer, parameter :: r4 = selected_real_kind(5)
  integer, parameter :: r8 = selected_real_kind(8,70)!8
  integer, parameter :: i8 = selected_int_kind(13)
end module precision
