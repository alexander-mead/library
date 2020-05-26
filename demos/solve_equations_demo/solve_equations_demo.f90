PROGRAM solve_equations_test

   USE basic_operations
   USE array_operations
   USE solve_equations

   IMPLICIT NONE
   INTEGER :: imode

   WRITE(*, *)

   CALL read_command_argument(1, imode, '', -1)

   IF (imode == -1) THEN
      WRITE(*, *) 'imode = 1: Numerical solve'
      WRITE(*, *) 'imode = 2: Quadratic solve'
      READ(*, *) imode
      WRITE(* ,*)
   END IF

   IF (imode == 1) THEN
      CALL numerical_solve_demo()
   ELSE IF (imode == 2) THEN
      CALL quadratic_solve_demo()
   ELSE
      STOP 'Error. imode specified incorreclty'
   END IF

CONTAINS

   SUBROUTINE quadratic_solve_demo()

      IMPLICIT NONE
      REAL :: x(2)
      REAL, PARAMETER :: a = 1.
      REAL, PARAMETER :: b = -5.
      REAL, PARAMETER :: c = 6.

      WRITE(*, *) 'Solving quadratic: a*x^2 + b*x + c = 0'
      WRITE(*, *) 'a =', a
      WRITE(*, *) 'b =', b
      WRITE(*, *) 'c =', c
      
      x = solve_quadratic(a, b, c)

      WRITE(*, *) 'Solutions:', x(1), x(2)
      WRITE(*, *)

   END SUBROUTINE quadratic_solve_demo

   SUBROUTINE numerical_solve_demo()

      IMPLICIT NONE
      REAL, ALLOCATABLE :: x(:), LHS(:), RHS(:)
      REAL :: sol
      INTEGER :: i

      REAL, PARAMETER :: xmin = 0.
      REAL, PARAMETER :: xmax = 2.
      INTEGER, PARAMETER :: n = 128

      REAL, PARAMETER :: alpha = 1.14619322062058

      WRITE (*, *) 'Testing numerical routines to solve equations'
      WRITE (*, *)

      WRITE (*, *) 'Solving equation x = ln(2+x)'
      WRITE (*, *) 'Wolfram alpha solution: x =', alpha
      WRITE (*, *)

      !Fill array of x value
      CALL fill_array(xmin, xmax, x, n)

      !This example program solves the equation log(2+x) = x
      !This equation has solution x=1.14619 (Wolfram Alpha) but no analytic solution
      !Fill L(x)=log(2+x) and R(x)=x
      ALLOCATE (LHS(n), RHS(n))
      DO i = 1, n
         LHS(i) = log(2.+x(i))
         RHS(i) = x(i)
      END DO

      WRITE (*, *) 'Solving using bisect solve'
      sol = bisect_solve(x, RHS-LHS, n, 1e-6)
      WRITE (*, *) 'Solution:', sol
      WRITE (*, *) 'Accuracy:', -1.+sol/alpha
      WRITE (*, *)

      WRITE (*, *) 'Solving using find solve'
      sol = find_solve(x, RHS-LHS, n)
      WRITE (*, *) 'Solution:', sol
      WRITE (*, *) 'Accuracy:', -1.+sol/alpha
      WRITE (*, *)

   END SUBROUTINE numerical_solve_demo

END PROGRAM solve_equations_test
