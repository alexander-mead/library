PROGRAM fitting_demo

   ! TODO: Most of these are actually minimization tests and should be moved to minimization_demo.f90

   USE minimization
   USE calculus
   USE fitting

   IMPLICIT NONE
   INTEGER :: itest

   ! Initial white space
   WRITE(*,*)

   ! Choose demo
   WRITE (*, *) 'Choose demo'
   WRITE (*, *) '1 - Demo fit constant '
   READ (*, *) itest
   WRITE (*, *)

   IF (itest == 1) THEN

      CALL demo_fit_constant()

   ELSE 

      STOP 'FITTING_DEMO: Error, itest not specified correctly'

   END IF

CONTAINS

   SUBROUTINE demo_fit_constant()

      IMPLICIT NONE
      REAL, ALLOCATABLE :: a(:), w(:)
      REAL :: fit
      INTEGER :: i, j
      INTEGER, PARAMETER :: ndemo=2 ! Number of demos
      INTEGER, PARAMETER :: n=5 ! Array size to be fitted with a constant

      DO j=1,ndemo

         ALLOCATE(a(n),w(n))
         a(1)=1.
         a(2)=0.8
         a(3)=1.2
         a(4)=1.5
         a(5)=0.5
         w=1.
         IF(j==2) w(5)=0.

         WRITE(*,*) 'Demo: ', j
         WRITE(*,*) 'Array and weights'
         DO i=1,n
            WRITE(*,*) i, a(i), w(i)
         END DO
         WRITE(*,*)

         fit=fit_constant(a,w,n)
         WRITE(*,*) 'Fitted constant:', fit
         WRITE(*,*)

         DEALLOCATE(a,w)

      END DO

   END SUBROUTINE demo_fit_constant

END PROGRAM fitting_demo
