PROGRAM minimization_demo

   USE minimization

   IMPLICIT NONE 
   INTEGER :: imode  
   LOGICAL, PARAMETER :: verbose_examples = .TRUE.

   WRITE(*, *)
   WRITE(*, *) 'Chooose example'
   WRITE(*, *) '1 - 1D example'
   WRITE(*, *) '2 - 2D example'
   READ(*, *) imode
   WRITE(*, *)

   IF(imode == 1) THEN     
      CALL example_1D(verbose_examples)
   ELSE IF(imode == 2) THEN
      CALL example_2D(verbose_examples)
   ELSE
      STOP 'MINIMIZATION_DEMO: Error, imode not specified correctly'
   END IF
   
   CONTAINS

   SUBROUTINE example_1D(verbose)

      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: x(:), dx(:)
      INTEGER :: n
      REAL, PARAMETER :: tol = 1e-5

      n=1
      ALLOCATE(x(n), dx(n))

      ! Initial guess
      x(1) = 1.
      dx(1) = 1.
      WRITE(*, *) 'Initial guess'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *)

      CALL Nelder_Mead(x, dx, func_1D, n, tol, verbose)
      WRITE(*, *) 'Minimization found at'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *)

   END SUBROUTINE

   REAL FUNCTION func_1D(x, n)

      REAL, INTENT(IN) :: x(n)
      INTEGER, INTENT(IN) :: n
      REAL, PARAMETER :: a = -3.

      func_1D = (x(1)-a)**4

   END FUNCTION func_1D

   SUBROUTINE example_2D(verbose)

      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: x(:), dx(:)
      INTEGER :: n
      REAL, PARAMETER :: tol = 1e-5

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

      CALL Nelder_Mead(x, dx, func_2D, n, tol, verbose)
      WRITE(*, *) 'Minimization found at'
      WRITE(*, *) 'x:', x(1)
      WRITE(*, *) 'y:', x(2)
      WRITE(*, *)

   END SUBROUTINE

   REAL FUNCTION func_2D(x, n)

      REAL, INTENT(IN) :: x(n)
      INTEGER, INTENT(IN) :: n
      REAL, PARAMETER :: a = -5.
      REAL, PARAMETER :: b = 6.

      func_2D = (x(1)-a)**2+(x(2)-b)**2

   END FUNCTION func_2D

END PROGRAM