MODULE precision

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: sp
   PUBLIC :: dp
   PUBLIC :: qp

   !INTEGER, PARAMETER :: sp = KIND(1.0)
   !INTEGER, PARAMETER :: dp = KIND(1.d0)
   ! Taken from http://fortranwiki.org/fortran/show/Real+precision
   INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)    ! Single precision
   INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)  ! Double precision
   INTEGER, PARAMETER :: qp = selected_real_kind(33, 4931) ! Quadruple precision

END MODULE precision
