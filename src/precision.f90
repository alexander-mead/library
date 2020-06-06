MODULE precision

   !USE, INTRINSIC :: iso_fortran_env

   IMPLICIT NONE

   PRIVATE

   ! Integer precision
   PUBLIC :: int4
   PUBLIC :: int8

   ! Real precision
   PUBLIC :: sp
   PUBLIC :: dp
   PUBLIC :: qp

   ! Integer precision
   INTEGER, PARAMETER :: int4 = selected_int_kind(4)
   INTEGER, PARAMETER :: int8 = selected_int_kind(8)

   ! Hmmm....
   !INTEGER, PARAMETER :: sp = KIND(1.0)
   !INTEGER, PARAMETER :: dp = KIND(1.d0)

   ! Taken from http://fortranwiki.org/fortran/show/Real+precision
   INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)    ! Single precision
   INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)  ! Double precision
   INTEGER, PARAMETER :: qp = selected_real_kind(33, 4931) ! Quadruple precision

   ! Fortran 2008 http://fortranwiki.org/fortran/show/Real+precision
   ! Needs to USE, INTRINSIC :: iso_fortran_env which conflicts with default double precision
   !INTEGER, PARAMETER :: sp = REAL32
   !INTEGER, PARAMETER :: dp = REAL64
   !INTEGER, PARAMETER :: qp = REAL128
  
END MODULE precision
