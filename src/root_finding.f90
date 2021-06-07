MODULE root_finding

   ! TODO: Write a wrapper for solve routines
   USE basic_operations
   USE table_integer
   USE interpolate

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: solve_quadratic
   PUBLIC :: solve_gradient
   PUBLIC :: solve_find
   PUBLIC :: solve_bisect_tables
   PUBLIC :: solve_bisect_function

   INTEGER, PARAMETER :: iorder_bisect = 3
   INTEGER, PARAMETER :: ifind_bisect = ifind_split
   INTEGER, PARAMETER :: iinterp_bisect = iinterp_Lagrange

CONTAINS

   FUNCTION solve_quadratic(a, b, c)

      ! Get the two real solutions of a quadratic ax^2 + bx + c = 0
      ! TODO: Solve for imaginary solutions too
      REAL :: solve_quadratic(2)
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, INTENT(IN) :: c
      REAL :: disc, root

      disc = b**2-4.*a*c
      IF(disc >= 0.) THEN
         root = sqrt(disc)
      ELSE
         STOP 'SOLVE_QUADRATIC: Error, your quadratic has imaginary solutions'
      END IF

      solve_quadratic(1) = (root-b)/(2.*a)
      solve_quadratic(2) = (-root-b)/(2.*a)

   END FUNCTION solve_quadratic

   REAL FUNCTION solve_gradient(x1in, x2in, f, acc_opt)

      ! Solves for x0 with f(x0) = 0; initial guesses x1 and x2
      ! e.g., f(x) = x^2-4 would find either x = -2 or x = 2 depending on the initial guess
      USE special_functions
      REAL, INTENT(IN) :: x1in, x2in
      REAL, EXTERNAL :: f
      REAL, OPTIONAL, INTENT(IN) :: acc_opt
      REAL :: x1, x2, f1, f2, a0, a1, xnew, acc
      REAL, PARAMETER :: acc_def = 1e-6
      INTERFACE
         FUNCTION f(x)
            REAL, INTENT(IN) :: x
         END FUNCTION f
      END INTERFACE

      acc = default_or_optional(acc_def, acc_opt)
      x1 = x1in; x2 = x2in
      f1 = f(x1); f2 = f(x2)
      DO
         IF ((min(abs(f1), abs(f2))) <= acc) EXIT
         CALL fix_polynomial(a1, a0, [x1, x2], [f1, f2])
         xnew = -a0/a1 ! x intercept
         IF (abs(f1) < abs(f2)) THEN
            x2 = xnew; f2 = f(x2) ! Guess 1 was better, so replace 2...
         ELSE
            x1 = xnew; f1 = f(x1) ! ...otherwise guess 2 was better, so replace 1
         END IF
      END DO

      IF (abs(f1) <= acc) THEN
         solve_gradient = x1
      ELSE
         solve_gradient = x2
      END IF

   END FUNCTION solve_gradient

   REAL FUNCTION solve_find(xtab, ytab)

      ! Solves y(x) = 0 for x
      ! Can also do x = exp(solve_find(log(xtab), ytab) if better
      ! e.g., consider y(x) = x-4; solve_find(x, y) will return 4 (i.e., value when y=0) when x, y tabulated
      ! e.g., consider y(x) = tan(x)-x; solve_find(x, y) will return 4.493 (or 0) iff this value is bracketed by tabulations
      REAL, INTENT(IN) :: xtab(:)
      REAL, INTENT(IN) :: ytab(:)
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2
      REAL, PARAMETER :: a = 0. ! y(x)=a; a=0 here

      solve_find = find(a, ytab, xtab, size(xtab), iorder, ifind, imeth)

   END FUNCTION solve_find

   REAL FUNCTION solve_bisect_tables(xtab, ytab, acc)

      ! Solves y(x)=0 for x, f(x) should be monotonic and cross f=0 once only
      USE interpolate
      REAL, INTENT(IN) :: xtab(:)
      REAL, INTENT(IN) :: ytab(:)
      REAL, INTENT(IN) :: acc
      REAL :: x1, x2, y1, y2, x, y
      INTEGER :: i, n
      INTEGER, PARAMETER :: iorder = iorder_bisect
      INTEGER, PARAMETER :: ifind = ifind_bisect
      INTEGER, PARAMETER :: iinterp = iinterp_bisect

      ! Check array sizes
      n = size(xtab)
      IF (n /= size(ytab)) STOP 'BISECT_SOLVE: Error, x and y arrays must be the same size'

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
         y = find(x, xtab, ytab, n, iorder, ifind, iinterp)
         IF (abs(y) < acc) THEN
            EXIT
         ELSE IF (positive(y1) .AND. positive(y)) THEN
            x1 = x
            y1 = y
         ELSE
            x2 = x
            y2 = y
         END IF
      END DO

      solve_bisect_tables = x

   END FUNCTION solve_bisect_tables

   REAL FUNCTION solve_bisect_function(x1_ini, x2_ini, f, acc)

      ! Solve f(x0)=0 for x0 using a bisect method, f(x) should be monotonic and cross 0 only once
      REAL, INTENT(IN) :: x1_ini ! Lower limit for initial guess for range of root
      REAL, INTENT(IN) :: x2_ini ! Upper limit for initial guess for range of root
      REAL, EXTERNAL :: f        ! Function of which to find root
      REAL, INTENT(IN) :: acc    ! Accuracy in f(x) ~ 0
      REAL :: x1, x2, xm
      REAL :: f1, f2, fm

      INTERFACE
         FUNCTION f(x)
            REAL, INTENT(IN) :: x
         END FUNCTION f
      END INTERFACE

      ! Calcualte the values of the function at the initial guess points
      x1 = x1_ini; x2 = x2_ini
      f1 = f(x1); f2 = f(x2)
      IF (positive(f1) .AND. positive(f2)) STOP 'FIND_ROOT: Error, initial guesses must bracket the root'

      ! Loop until convergence achieved
      DO

         ! Calculate the mid-point between x1 and x2 and the value of the function here
         xm = (x1+x2)/2.
         fm = f(xm)

         ! Check for convergence, or update the range between which the root is to be found
         IF (abs(fm) < acc) THEN
            ! Convergence has been met
            EXIT
         ELSE IF (positive(fm) .AND. positive(f1)) THEN
            ! If f(x1) and f(xm) have the same sign then the root is between xm and x2, so update x1, f1 ...
            x1 = xm; f1 = fm
         ELSE
            ! ...otherwise update x2 and f2.
            x2 = xm; f2 = fm
         END IF

      END DO

      solve_bisect_function = xm

   END FUNCTION solve_bisect_function

END MODULE root_finding
