!> Deliveres the double precision accuracy for the project.

module accuracy
  
  implicit none

  !> Kind for double precision.
  integer, parameter :: dp = selected_real_kind(12,99)

end module accuracy
