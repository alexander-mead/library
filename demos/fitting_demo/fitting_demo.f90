PROGRAM fitting_demo

   ! TODO: Most of these are actually minimization tests and should be moved to minimization_demo.f90

   USE minimization
   USE calculus
   USE fitting

   IMPLICIT NONE
   REAL :: df, min, min2(2)
   REAL :: k(100), d2(100)!, d2mod(100)
   INTEGER :: i, ngrid, ngrid3(3)
   INTEGER :: itest
   REAL :: x1, x2, x3, f1, f2, f3
   REAL :: x, xmin, hi, lo, ref
   REAL :: xmin3(3), lo3(3), hi3(3), ref3(3)
   REAL :: a, b, c

   !Accuracy parameter
   REAL, PARAMETER :: acc = 1e-5

   !Initial guess
   REAL, PARAMETER :: x0 = 3.

   !One-dimensional gradient fit

   ! Initial white space
   WRITE(*,*)

   ! Choose demo
   WRITE (*, *) 'Choose demo'
   !WRITE (*, *) '1 - Demo gradient minimization 1D'
   !WRITE (*, *) '2 - Demo gradient minimization 2D'
   !WRITE (*, *) '3 - Demo grid minimization 1D'
   !WRITE (*, *) '4 - Demo grid minimization 2D'
   !WRITE (*, *) '5 - Demo grid minimization 3D'
   !WRITE (*, *) '6 - Demo quadratic minimization 1D'
   !WRITE (*, *) '7 - Demo adaptive minimization 1D'
   !WRITE (*, *) '8 - Demo adaptive minimization 3D'
   WRITE (*, *) '9 - Demo fit constant '
   READ (*, *) itest
   WRITE (*, *)

   IF (itest == 1) THEN

      ! WRITE (*, *) 'Initial guess:', x0

      ! df = derivative(quadratic, x0, acc)
      ! WRITE (*, *) 'Derivative at initial guess', df

      ! min = gradient_minimization_1D(quadratic, x0, acc)

      ! WRITE (*, *) 'Gradient minimum:', min
      ! !WRITE(*,*) 'Error:', min/quadratic_extrema()
      ! WRITE (*, *)

   ELSE IF (itest == 2) THEN

      ! min2 = gradient_minimization_2D(test, 10., 3., acc)
      ! WRITE (*, *) min2(1)
      ! WRITE (*, *) min2(2)

   ELSE IF (itest == 3) THEN

      ! DO i = 1, 100
      !    k(i) = float(i)
      !    d2(i) = k(i)**3.2
      ! END DO

      ! CALL grid_minimization_1D(0., 5., 501, k, d2)

   ELSE IF (itest == 4) THEN

!!$     DO i=1,100
!!$        k(i)=float(i)
!!$        d2(i)=k(i)**3.+k(i)**0.5
!!$     END DO
!!$
!!$     CALL fit2(0.,5.,501,0.,5.,501,k,d2)

   ELSE IF (itest == 5) THEN

!!$     DO i=1,100
!!$        k(i)=float(i)
!!$        d(i)=k(i)**2.+sin(3.*k(i))+5.*k(i)
!!$     END DO
!!$
!!$     CALL fit(0.,10.,101,0.,10.,101,0.,10.,101,k,d)

   ELSE IF (itest == 6) THEN

      ! x1 = 1.
      ! x2 = 2.
      ! x3 = 3.

      ! f1 = quadratic(x1)
      ! f2 = quadratic(x2)
      ! f3 = quadratic(x3)

      ! x = quadratic_minimization_1D(x1, f1, x2, f2, x3, f3)

      ! WRITE (*, *) 'Minimum:', x
      ! WRITE (*, *)

   ELSE IF (itest == 7) THEN

      ! a = 1.
      ! b = -4.
      ! c = 2.

      ! WRITE (*, *)
      ! WRITE (*, *) 'Adaptive grid search algorithm'
      ! WRITE (*, *)

      ! CALL adaptive_minimization_1D(xmin, quadratic, -10., 10., 10., 10, 5)

      ! WRITE (*, *) 'Best fit minimum at:', xmin
      ! WRITE (*, *) 'True miniumum at:', quadratic_extrema()
      ! WRITE (*, *)

   ELSE IF (itest == 8) THEN

      ! a = 3.
      ! b = 4.
      ! c = 5.

      ! xmin = 1.
      ! lo = -10.
      ! hi = 10.
      ! ref = 5.
      ! ngrid = 10

      ! WRITE (*, *)
      ! WRITE (*, *) 'Adaptive grid search algorithm'
      ! WRITE (*, *)

      ! CALL adaptive_minimization_3D(xmin3, quadratic3, lo3, hi3, ref3, ngrid3, 5)

      ! WRITE (*, *) 'Best fit minimum at:', xmin
      ! WRITE (*, *) 'True miniumum at:', a, b, c
      ! WRITE (*, *)

   ELSE IF (itest == 9) THEN

      CALL demo_fit_constant()

   ELSE 

      STOP 'FITTING_DEMO: Error, itest not specified correctly'

   END IF

CONTAINS

   ! REAL FUNCTION quadratic(x)

   !    IMPLICIT NONE
   !    REAL, INTENT(IN) :: x

   !    quadratic = a*x**2+b*x+c

   ! END FUNCTION quadratic

   ! REAL FUNCTION quadratic3(x)

   !    IMPLICIT NONE
   !    REAL, INTENT(IN) :: x(3)

   !    quadratic3 = (x(1)-a)**2.+(x(2)-b)**2.+(x(3)-c)**2.

   ! END FUNCTION quadratic3

   ! REAL FUNCTION quadratic_extrema()

   !    IMPLICIT NONE

   !    quadratic_extrema = -b/(2.*a)

   ! END FUNCTION quadratic_extrema

   ! FUNCTION test(x, y)

   !    IMPLICIT NONE
   !    REAL :: test
   !    REAL, INTENT(IN) :: x, y
   !    REAL :: A, B, C, D, E, F

   !    A = 1.
   !    B = 1.
   !    C = 0.
   !    D = 0.
   !    E = 0.
   !    F = 2.

   !    test = A*(x**2.)+B*(y**2.)+C*(x*y)+D*x+E*y+F

   ! END FUNCTION test

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
