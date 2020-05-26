MODULE solve_equations

   ! TODO: Write a wrapper for solve routines

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: solve_quadratic
   PUBLIC :: find_solve
   PUBLIC :: bisect_solve

CONTAINS

   FUNCTION solve_quadratic(a, b, c)

      IMPLICIT NONE
      REAL :: solve_quadratic(2)
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, INTENT(IN) :: c
      REAL :: disc, root

      disc = b**2-4.*a*c
      IF(disc > 0.) THEN
         root = sqrt(disc)
      ELSE
         STOP 'SOLVE_QUADRATIC: Error, your quadratic has imaginary solutions'
      END IF

      solve_quadratic(1) = (root-b)/(2.*a)
      solve_quadratic(2) = (-root-b)/(2.*a)

   END FUNCTION solve_quadratic

   REAL FUNCTION find_solve(xtab, ytab, n)

      ! Solves y(x)=0. for x
      USE interpolate
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: xtab(n)
      REAL, INTENT(IN) :: ytab(n)
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2
      REAL, PARAMETER :: a = 0. ! y(x)=a; a=0 here

      find_solve = find(a, ytab, xtab, n, iorder, ifind, imeth)

   END FUNCTION find_solve

   REAL FUNCTION bisect_solve(xtab, ytab, n, acc)

      ! Solves y(x)=0 for x, f(x) should be monotonic and cross f=0. once only
      USE basic_operations
      USE interpolate
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: xtab(n)
      REAL, INTENT(IN) :: ytab(n)
      REAL, INTENT(IN) :: acc
      REAL :: x1, x2, y1, y2, x, y
      INTEGER :: i
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      ! Initial values taken from top and bottom of table
      x1 = xtab(1)
      x2 = xtab(n)
      y1 = ytab(1)
      y2 = ytab(n)

      ! Now iterate until desired accuracy is reached
      i = 0
      DO
         i = i+1
         x = 0.5*(x1+x2)
         y = find(x, xtab, ytab, n, iorder, ifind, imeth)
         IF (abs(y) < acc) THEN
            EXIT
         ELSE IF (positive(y1) .EQV. positive(y)) THEN
            x1 = x
            y1 = y
         ELSE
            x2 = x
            y2 = y
         END IF
      END DO

      bisect_solve = x

   END FUNCTION bisect_solve

END MODULE solve_equations
