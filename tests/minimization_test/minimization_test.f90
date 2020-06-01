PROGRAM minimization_test

   USE basic_operations
   USE minimization

   IMPLICIT NONE

   INTEGER :: imode
   LOGICAL, PARAMETER :: verbose_tests = .FALSE.
   INTEGER, PARAMETER :: ntest = 2
   LOGICAL :: ifail

   WRITE (*, *)

   DO imode = 1, ntest

      IF (imode == 1) THEN     
         CALL test_minimization_1D(verbose_tests, ifail)
      ELSE IF (imode == 2) THEN
         CALL test_minimization_2D(verbose_tests, ifail)
      ELSE
         STOP 'MINIMIZATION_DEMO: Error, imode not specified correctly'
      END IF

      IF (ifail .EQV. .TRUE.) THEN
         STOP 'MINIMIZATION_TEST: Failed'
      END IF

   END DO
   
   CONTAINS

   SUBROUTINE test_minimization_1D(verbose, fail)

      LOGICAL, INTENT(IN) :: verbose
      LOGICAL, INTENT(OUT) :: fail
      REAL, ALLOCATABLE :: x(:), dx(:)
      REAL :: fom
      INTEGER :: n
      REAL, PARAMETER :: tol = 1e-8
      REAL, PARAMETER :: ttol = tol
      REAL, PARAMETER :: xmin = -3.

      fail = .FALSE.

      n = 1
      ALLOCATE(x(n), dx(n))

      ! Initial guess
      x(1) = 1.
      dx(1) = 1.
      IF (verbose) THEN
         WRITE(*, *) 'TEST_MINIMIZATION_1D: Initial guess'
         WRITE(*, *) 'TEST_MINIMIZATION_1D: x:', x(1)
         WRITE(*, *)
      END IF

      CALL Nelder_Mead(x, dx, fom, func_1D, tol, verbose)
      IF (verbose) THEN
         WRITE(*, *) 'TEST_MINIMIZATION_1D: Minimization found at'
         WRITE(*, *) 'TEST_MINIMIZATION_1D: x:', x(1)
         WRITE(*, *) 'TEST_MINIMIZATION_1D: Should be: ', xmin
         WRITE(*, *)
      END IF

      IF (.NOT. requal(x(1), xmin, ttol)) fail = .TRUE.

      IF (fail) THEN
         WRITE (*, *) 'TEST_MINIMIZATION_1D: Fail'
         STOP
      ELSE
         WRITE (*, *) 'TEST_MINIMIZATION_1D: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_minimization_1D

   REAL FUNCTION func_1D(x)

      REAL, INTENT(IN) :: x(:)
      REAL, PARAMETER :: a = -3.0

      func_1D = (x(1)-a)**4

   END FUNCTION func_1D

   SUBROUTINE test_minimization_2D(verbose, fail)

      LOGICAL, INTENT(IN) :: verbose
      LOGICAL, INTENT(OUT) :: fail
      REAL, ALLOCATABLE :: x(:), dx(:)
      REAL :: fom
      INTEGER :: n
      REAL, PARAMETER :: tol = 1e-8
      REAL, PARAMETER :: ttol = 1e-4
      REAL, PARAMETER :: xmin = -5.
      REAL, PARAMETER :: ymin = 6.

      fail = .FALSE.

      n=2
      ALLOCATE(x(n), dx(n))

      ! Initial guess
      x(1) = 1.
      x(2) = 1.
      dx(1) = 1.
      dx(2) = 1.
      IF (verbose) THEN
         WRITE(*, *) 'TEST_MINIMIZATION_2D: Initial guess'
         WRITE(*, *) 'TEST_MINIMIZATION_2D: x:', x(1)
         WRITE(*, *) 'TEST_MINIMIZATION_2D: y:', x(2)
         WRITE(*, *)
      END IF

      CALL Nelder_Mead(x, dx, fom, func_2D, tol, verbose)
      IF (verbose) THEN
         WRITE(*, *) 'TEST_MINIMIZATION_2D: Minimization found at'
         WRITE(*, *) 'TEST_MINIMIZATION_2D: x:', x(1)
         WRITE(*, *) 'TEST_MINIMIZATION_2D: y:', x(2)
         WRITE(*, *)
      END IF

      IF ((.NOT. requal(x(1), xmin, ttol)) .OR. (.NOT. requal(x(2), ymin, ttol))) fail = .TRUE.

      IF (fail) THEN
         WRITE (*, *) 'TEST_MINIMIZATION_2D: Fail'
         STOP
      ELSE
         WRITE (*, *) 'TEST_MINIMIZATION_2D: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_minimization_2D

   REAL FUNCTION func_2D(x)

      USE special_functions
      REAL, INTENT(IN) :: x(:)
      REAL, PARAMETER :: a = -5.
      REAL, PARAMETER :: b = 6.0

      func_2D = (x(1)-a)**2+(x(2)-b)**2 ! Quadratic
      !func_2D = Rosenbrock(x(1), x(2))
      !func_2D = Himmelblau(x(1), x(2))

   END FUNCTION func_2D

END PROGRAM minimization_test