PROGRAM minimization_demo

   USE minimization
   USE basic_operations

   IMPLICIT NONE 
   INTEGER :: imode  
   LOGICAL, PARAMETER :: verbose_examples = .TRUE.

   CALL read_command_argument(1, imode, '', 0)

   IF (imode == 0) THEN
      WRITE(*, *)
      WRITE(*, *) 'Chooose example'
      WRITE(*, *) ' 1 - Nelder-Mead 1D'
      WRITE(*, *) ' 2 - Nelder-Mead 2D'
      WRITE(*, *) ' 3 - Gradient minimization 1D'
      WRITE(*, *) ' 4 - Gradient minimization 2D'
      WRITE(*, *) ' 5 - Grid minimization 1D'
      WRITE(*, *) ' 6 - Quadratic minimization 1D'
      WRITE(*, *) ' 7 - Gradient minimization 1D'
      WRITE(*, *) ' 8 - Adaptive minimization 1D'
      WRITE(*, *) ' 9 - Adaptive minimization 3D'
      WRITE(*, *) '10 - Nelder-Mead multiple 2D'
      READ(*, *) imode
      WRITE(*, *)
   END IF

   IF (imode == 1) THEN     
      CALL Nelder_Mead_example_1D(verbose_examples)
   ELSE IF (imode == 2) THEN
      CALL Nelder_Mead_example_2D(verbose_examples)
   ELSE IF (imode == 3) THEN
      CALL gradient_minimization_example_1D()
   ELSE IF (imode == 4) THEN
      CALL gradient_minimization_example_2D()
   ELSE IF (imode == 5) THEN
      CALL grid_minimization_example_1D()
   ELSE IF (imode == 6) THEN
      CALL quadratic_minimization_example_1D()
   ELSE IF (imode == 7) THEN
      CALL gradient_minimization_example_1D()
   ELSE IF (imode == 8) THEN
      CALL adaptive_minimization_example_1D()
   ELSE IF (imode == 9) THEN
      CALL adaptive_minimization_example_3D()
   ELSE IF (imode == 10) THEN
      CALL Nelder_Mead_multiple_example_2D(verbose_examples)
   ELSE
      STOP 'MINIMIZATION_DEMO: Error, imode not specified correctly'
   END IF
   
   CONTAINS

   SUBROUTINE Nelder_Mead_example_1D(verbose)

      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: x(:), dx(:)
      REAL :: fom
      INTEGER :: n
      REAL, PARAMETER :: tol = 1e-8

      n=1
      ALLOCATE(x(n), dx(n))

      ! Initial guess
      x(1) = 1.
      dx(1) = 1.
      WRITE(*, *) 'Initial guess'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *)

      CALL Nelder_Mead(x, dx, n, fom, func_1D, tol, verbose)
      WRITE(*, *) 'Minimization found at'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *)

   END SUBROUTINE Nelder_Mead_example_1D

   REAL FUNCTION func_1D(x, n)

      REAL, INTENT(IN) :: x(n)
      INTEGER, INTENT(IN) :: n
      REAL, PARAMETER :: a = -3.

      func_1D = (x(1)-a)**4

   END FUNCTION func_1D

   SUBROUTINE Nelder_Mead_example_2D(verbose)

      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: x(:), dx(:)
      REAL :: fom
      INTEGER :: n
      REAL, PARAMETER :: tol = 1e-8

      n=2
      ALLOCATE(x(n), dx(n))

      ! Initial guess
      x(1) = 1.
      x(2) = 1.
      dx(1) = 1.
      dx(2) = 1.
      WRITE(*, *) 'Initial guess'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *) 'y:', x(2)
      WRITE(*, *)

      CALL Nelder_Mead(x, dx, n, fom, func_2D, tol, verbose)
      WRITE(*, *) 'Minimization found at'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *) 'y:', x(2)
      WRITE(*, *)

   END SUBROUTINE Nelder_Mead_example_2D

   SUBROUTINE Nelder_Mead_multiple_example_2D(verbose)

      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: x(:), xmin(:), xmax(:), dx(:)
      REAL :: fom
      INTEGER :: n
      REAL, PARAMETER :: tol = 1e-8
      INTEGER, PARAMETER :: m = 10

      n=2
      ALLOCATE(x(n), xmin(n), xmax(n), dx(n))

      ! Initial guess
      x = 0.
      xmin = 0.
      xmax = 6.
      dx(1) = 1.
      dx(2) = 1.
      WRITE(*, *) 'Initial guess'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *) 'y:', x(2)
      WRITE(*, *)

      CALL Nelder_Mead_multiple(x, xmin, xmax, dx, n, fom, m, func_2D, tol, verbose)
      WRITE(*, *) 'Minimization found at'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *) 'y:', x(2)
      WRITE(*, *)

   END SUBROUTINE Nelder_Mead_multiple_example_2D

   REAL FUNCTION func_2D(x, n)

      USE special_functions
      REAL, INTENT(IN) :: x(n)
      INTEGER, INTENT(IN) :: n
      REAL, PARAMETER :: a = -5.
      REAL, PARAMETER :: b = 6.

      !func_2D = (x(1)-a)**2+(x(2)-b)**2 ! Quadratic
      !func_2D = Rosenbrock(x(1), x(2))
      func_2D = Himmelblau(x(1), x(2))

   END FUNCTION func_2D

   SUBROUTINE gradient_minimization_example_1D()

      USE calculus
      REAL :: min, df
      REAL, PARAMETER :: x0 = 3.
      REAL, PARAMETER :: acc = 1e-5

      STOP 'NOT TESTED'

      WRITE (*, *) 'Initial guess:', x0

      df = derivative(quadratic, x0, acc)
      WRITE (*, *) 'Derivative at initial guess', df

      min = gradient_minimization_1D(quadratic, x0, acc)

      WRITE (*, *) 'Gradient minimum:', min
      !WRITE(*,*) 'Error:', min/quadratic_extrema()
      WRITE (*, *)

   END SUBROUTINE gradient_minimization_example_1D

   REAL FUNCTION quadratic(x)

      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, PARAMETER :: a=1.
      REAL, PARAMETER :: b=1.
      REAL, PARAMETER :: c=1.

      quadratic = a*x**2+b*x+c

   END FUNCTION quadratic

   SUBROUTINE gradient_minimization_example_2D()

      REAL :: min(2)
      REAL, PARAMETER :: acc = 1e-5

      STOP 'NOT TESTED'

      min = gradient_minimization_2D(test, 10., 3., acc)
      WRITE (*, *) min(1)
      WRITE (*, *) min(2)

   END SUBROUTINE gradient_minimization_example_2D

   REAL FUNCTION test(x, y)

      REAL, INTENT(IN) :: x, y
      REAL :: A, B, C, D, E, F

      A = 1.
      B = 1.
      C = 0.
      D = 0.
      E = 0.
      F = 2.

      test = A*(x**2.)+B*(y**2.)+C*(x*y)+D*x+E*y+F

   END FUNCTION test

   SUBROUTINE grid_minimization_example_1D()

      INTEGER :: i
      REAL :: k(100), d2(100)

      STOP 'NOT TESTED'

      DO i = 1, 100
         k(i) = float(i)
         d2(i) = k(i)**3.2
      END DO

      CALL grid_minimization_1D(0., 5., 501, k, d2)

   END SUBROUTINE grid_minimization_example_1D

   SUBROUTINE quadratic_minimization_example_1D()

      REAL :: x1, x2, x3
      REAL :: f1, f2, f3
      REAL :: x

      STOP 'NOT TESTED'

      x1 = 1.
      x2 = 2.
      x3 = 3.

      f1 = quadratic(x1)
      f2 = quadratic(x2)
      f3 = quadratic(x3)

      x = quadratic_minimization_1D(x1, f1, x2, f2, x3, f3)

      WRITE (*, *) 'Minimum:', x
      WRITE (*, *)

   END SUBROUTINE quadratic_minimization_example_1D

   SUBROUTINE adaptive_minimization_example_1D()

      REAL :: a, b, c
      REAL :: xmin

      STOP 'NOT TESTED'

      a = 1.
      b = -4.
      c = 2.

      WRITE (*, *)
      WRITE (*, *) 'Adaptive grid search algorithm'
      WRITE (*, *)

      CALL adaptive_minimization_1D(xmin, quadratic, -10., 10., 10., 10, 5)

      WRITE (*, *) 'Best fit minimum at:', xmin
      WRITE (*, *) 'True miniumum at:', quadratic_extrema(a, b)
      WRITE (*, *)

   END SUBROUTINE adaptive_minimization_example_1D

   REAL FUNCTION quadratic_extrema(a, b)

      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b

      quadratic_extrema = -b/(2.*a)

   END FUNCTION quadratic_extrema

   SUBROUTINE adaptive_minimization_example_3D()

      REAL :: a, b, c
      REAL :: xmin, lo, hi, ref
      REAL :: xmin3(3), lo3(3), hi3(3), ref3(3)
      INTEGER :: ngrid, ngrid3(3)

      STOP 'NOT TESTED'

      a = 3.
      b = 4.
      c = 5.

      xmin = 1.
      lo = -10.
      hi = 10.
      ref = 5.
      ngrid = 10

      WRITE (*, *)
      WRITE (*, *) 'Adaptive grid search algorithm'
      WRITE (*, *)

      CALL adaptive_minimization_3D(xmin3, quadratic3, lo3, hi3, ref3, ngrid3, 5)

      WRITE (*, *) 'Best fit minimum at:', xmin
      WRITE (*, *) 'True miniumum at:', a, b, c
      WRITE (*, *)

   END SUBROUTINE adaptive_minimization_example_3D

   REAL FUNCTION quadratic3(x)

      IMPLICIT NONE
      REAL, INTENT(IN) :: x(3)
      REAL, PARAMETER :: a = 1.
      REAL, PARAMETER :: b = 1.
      REAL, PARAMETER :: c = 1.

      quadratic3 = (x(1)-a)**2.+(x(2)-b)**2.+(x(3)-c)**2.

   END FUNCTION quadratic3

END PROGRAM minimization_demo