MODULE minimization

   USE sorting
   USE statistics
   IMPLICIT NONE

   PRIVATE

   PUBLIC :: Nelder_Mead

   INTEGER, PARAMETER :: isort_Nelder_Mead = isort_bubble
   REAL, PARAMETER :: alpha_Nelder_Mead = 1.  ! Reflection coefficient (alpha > 0; standard alpha = 1)
   REAL, PARAMETER :: gamma_Nelder_Mead = 2.  ! Expansion coefficient (gamma > 1; standard gamma = 2)
   REAL, PARAMETER :: rho_Nelder_Mead = 0.5   ! Contraction coefficient (0 < rho < 0.5; standard rho = 0.5)
   REAL, PARAMETER :: sigma_Nelder_Mead = 0.5 ! Shrink coefficient (standard sigma = 0.5)

   CONTAINS

   SUBROUTINE Nelder_Mead(x, dx, f, n, tol, verbose)

      ! Nelder-Mead simplex for fiding minima of a function
      ! Coded up using https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
      IMPLICIT NONE
      REAL, INTENT(INOUT) :: x(n)
      REAL, INTENT(IN) :: dx(n)
      REAL, EXTERNAL :: f
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: tol
      LOGICAL, INTENT(IN) :: verbose
      REAL :: xx(n+1, n), ff(n+1)
      REAL :: xo(n), xr(n), xe(n), xc(n)
      REAL :: fo, fr, fe, fc
      INTEGER :: i, j, ii
      REAL, PARAMETER :: alpha = alpha_Nelder_Mead  ! Reflection coefficient (alpha > 0; standard alpha = 1)
      REAL, PARAMETER :: gamma = gamma_Nelder_Mead  ! Expansion coefficient (gamma > 1; standard gamma = 2)
      REAL, PARAMETER :: rho = rho_Nelder_Mead      ! Contraction coefficient (0 < rho < 0.5; standard rho = 0.5)
      REAL, PARAMETER :: sigma = sigma_Nelder_Mead  ! Shrink coefficient (standard sigma = 0.5)

      INTERFACE
         FUNCTION f(x, n)
            REAL, INTENT(IN) :: x(n)
            INTEGER, INTENT(IN) :: n
         END FUNCTION f
      END INTERFACE

      ! Set initial test points
      DO i = 1, n+1
         DO j = 1, n
            xx(i, j) = x(j)
            IF(i == j+1) xx(i, j) = xx(i, j)+dx(j)
         END DO
      END DO

      ! Evaluate function at initial points
      DO i = 1, n+1
         ff(i) = f(xx(i, :), n)
      END DO

      ii = 0
      DO

         ii = ii+1

         ! Sort the points from best to worst
         CALL Nelder_Mead_sort(xx, ff, n)

         IF(verbose) WRITE(*, *) ii, ff(1), (xx(i, 1), i = 1, n)

         ! Decide on convergence
         IF(Nelder_Mead_termination(ff, n, tol)) THEN
            DO i = 1, n
               x(i) = xx(1, i)
            END DO
            EXIT
         END IF

         ! Calculate centroid of 1...n (not n+1; the worst point)
         xo = Nelder_Mead_centroid(xx, n)

         ! Calculate the reflected point of n+1 about the centroid
         xr = xo+alpha*(xo-xx(n+1, :))
         fr = f(xr, n)

         IF (ff(1) <= fr .AND. fr < ff(n)) THEN
            ! Keep the reflected point if it is not the best but better than the second worst
            xx(n+1, :) = xr
            ff(n+1) = fr
            CYCLE
         ELSE IF (fr < ff(1)) THEN
            ! Calculate the expanded point if the reflected point is the best so far
            xe = xo+gamma*(xr-xo)
            fe = f(xe, n)
            IF(fe < fr) THEN
               ! Keep the expansed point if it is better than the reflected ...
               xx(n+1, :) = xe
               ff(n+1) = fe
            ELSE
               ! ... otherwise keep the reflected.
               xx(n+1, :) = xr
               ff(n+1) = fr
            END IF
            CYCLE
         ELSE
            ! Here it is certain that the reflected point is worse than the second worst point
            ! Calculate the contracted point
            xc = xo+rho*(xx(n+1, :)-xo)
            fc = f(xc, n)
            IF(fc < ff(n+1)) THEN
               ! Keep the contracted point if it is better than the worst point
               xx(n+1, :) = xc
               ff(n+1) = fc
               CYCLE
            ELSE
               ! The contracted point is worst, and we must shrink
               ! Calculate the shrinkage for all except the best point
               DO i = 2, n+1
                  xx(i, :) = xx(1, :)+sigma*(xx(i, :)-xx(1, :))
                  ff(i) = f(xx(i, :), n)
               END DO
            END IF
         END IF

         ! ! Shrinkage
         ! DO i = 2, n+1
         !    xx(i, :) = xx(1, :)+sigma*(xx(i, :)-xx(1, :))
         !    ff(i) = f(xx(i, :), n)
         ! END DO

      END DO

      ! Report the minimization point
      x = xx(1, :)

   END SUBROUTINE Nelder_Mead

   FUNCTION Nelder_Mead_centroid(x, n)

      ! Calculate the centroid of all points except n+1
      IMPLICIT NONE
      REAL :: Nelder_Mead_centroid(n)
      REAL, INTENT(IN) :: x(n+1, n)
      INTEGER, INTENT(IN) :: n
      INTEGER :: i

      DO i = 1, n
         Nelder_Mead_centroid(i) = mean(x(:,i), n)
      END DO

   END FUNCTION Nelder_Mead_centroid

   SUBROUTINE Nelder_Mead_sort(x, f, n)

      ! Sort the points into order from best to worst
      IMPLICIT NONE
      REAL, INTENT(INOUT) :: x(n+1, n)
      REAL, INTENT(INOUT) :: f(n+1)
      INTEGER, INTENT(IN) :: n
      INTEGER :: i, j(n+1)
      INTEGER, PARAMETER :: isort = isort_Nelder_Mead

      CALL index(f, j, n+1, isort)
      CALL reindex(f, j, n+1)
      DO i = 1, n
         CALL reindex(x(:,i), j, n+1)
      END DO

   END SUBROUTINE Nelder_Mead_sort

   LOGICAL FUNCTION Nelder_Mead_termination(f, n, tol)

      ! Determine if the minimization has converged
      IMPLICIT NONE
      REAL, INTENT(IN) :: f(n+1)
      INTEGER, INTENT(IN) :: n
      !REAL, PARAMETER :: tolerance = tolerance_Nelder_Mead
      REAL, INTENT(IN) :: tol
      REAL :: sigma

      ! Calculate the standard deviation of all points
      sigma = standard_deviation(f, n+1)

      ! Decide on termination
      IF(sigma <= tol) THEN
         Nelder_Mead_termination = .TRUE.
      ELSE
         Nelder_Mead_termination = .FALSE.
      END IF

   END FUNCTION Nelder_Mead_termination

END MODULE minimization