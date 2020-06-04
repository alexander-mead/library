MODULE minimization

   USE sorting
   USE statistics

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: Nelder_Mead
   PUBLIC :: Nelder_Mead_multiple
   PUBLIC :: Nelder_Mead_centroid
   PUBLIC :: Nelder_Mead_sort
   PUBLIC :: Nelder_Mead_termination
   PUBLIC :: adaptive_minimization_1D
   PUBLIC :: adaptive_minimization_3D
   PUBLIC :: grid_minimization_1D
   PUBLIC :: gradient_minimization_1D
   PUBLIC :: gradient_minimization_2D
   PUBLIC :: quadratic_extremum_1D
   PUBLIC :: find_array_maximum

   INTEGER, PARAMETER :: isort_Nelder_Mead = isort_bubble
   REAL, PARAMETER :: alpha_Nelder_Mead = 1.  ! Reflection coefficient (alpha > 0; standard alpha = 1)
   REAL, PARAMETER :: gamma_Nelder_Mead = 2.  ! Expansion coefficient (gamma > 1; standard gamma = 2)
   REAL, PARAMETER :: rho_Nelder_Mead = 0.5   ! Contraction coefficient (0 < rho < 0.5; standard rho = 0.5)
   REAL, PARAMETER :: sigma_Nelder_Mead = 0.5 ! Shrink coefficient (standard sigma = 0.5)

   CONTAINS

   SUBROUTINE Nelder_Mead_multiple(x, xmin, xmax, dx, fom, m, f, tol, verbose)

      USE random_numbers
      REAL, INTENT(OUT) :: x(:)
      REAL, INTENT(IN) :: xmin(:)
      REAL, INTENT(IN) :: xmax(:)
      REAL, INTENT(IN) :: dx(:)
      REAL, INTENT(OUT) :: fom
      INTEGER, INTENT(IN) :: m
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: tol
      LOGICAL, INTENT(IN) :: verbose
      REAL :: fy
      REAL, ALLOCATABLE :: y(:)
      INTEGER :: i, j, n

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin(:)           
         END FUNCTION f
      END INTERFACE

      n = size(x)
      IF (n /= size(xmin) .OR. n /= size(xmax) .OR. n /= size(dx)) THEN
         STOP 'NELDER_MEAD_MULTIPLE: Error; x, xmin, xmax, dx should all be the same size'
      END IF
      ALLOCATE(y(n))

      fom = huge(fom)
      DO i = 1, m

         DO j = 1, n
            y(j) = random_uniform(xmin(j), xmax(j))
         END DO

         CALL Nelder_Mead(y, dx, fy, f, tol, verbose)

         IF (fy < fom) THEN
            x = y
            fom = fy
         END IF

      END DO

      WRITE (*, *) 'NELDER_MEAD_MULTIPLE: Best fit:', fom, (x(j), j = 1, n)
      WRITE (*, *)


   END SUBROUTINE Nelder_Mead_multiple

   SUBROUTINE Nelder_Mead(x, dx, fom, f, tol, verbose)

      ! Nelder-Mead simplex for fiding minima of a function
      ! Coded up using https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
      ! x should be x(n)
      ! dx should be dx(n)
      REAL, INTENT(INOUT) :: x(:)
      REAL, INTENT(IN) :: dx(:)
      REAL, INTENT(OUT) :: fom
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: tol
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: xx(:, :), ff(:)
      REAL, ALLOCATABLE :: xo(:), xr(:), xe(:), xc(:)
      REAL :: fr, fe, fc
      INTEGER :: i, j, ii, n
      CHARACTER(len=256) :: operation, nstring
      REAL, PARAMETER :: alpha = alpha_Nelder_Mead  ! Reflection coefficient (alpha > 0; standard alpha = 1)
      REAL, PARAMETER :: gamma = gamma_Nelder_Mead  ! Expansion coefficient (gamma > 1; standard gamma = 2)
      REAL, PARAMETER :: rho = rho_Nelder_Mead      ! Contraction coefficient (0 < rho < 0.5; standard rho = 0.5)
      REAL, PARAMETER :: sigma = sigma_Nelder_Mead  ! Shrink coefficient (standard sigma = 0.5)

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin(:)           
         END FUNCTION f
      END INTERFACE

      ! Check array sizes
      n = size(x)
      IF (n /= size(dx)) STOP 'NELDER_MEAD: Error, x and dx must be the same size'
      ALLOCATE(xx(n+1, n), ff(n+1))
      ALLOCATE(xo(n), xr(n), xe(n), xc(n))

      ! Set initial test points
      DO i = 1, n+1
         DO j = 1, n
            xx(i, j) = x(j)
            IF(i == j+1) xx(i, j) = xx(i, j)+dx(j)
         END DO
      END DO

      ! Evaluate function at initial points
      DO i = 1, n+1
         ff(i) = f(xx(i, :))
      END DO

      ! Start the minimization
      ii = 0
      operation = 'Starting'
      WRITE(nstring, *) n+1
      DO

         ii = ii+1

         ! Sort the points from best to worst
         CALL Nelder_Mead_sort(xx, ff)
         
         IF (verbose) WRITE(*, '(A16,I10,'//trim(nstring)//'F15.7)') TRIM(operation), ii, ff(1), (xx(1, i), i = 1, n)

         ! Decide on convergence
         IF (Nelder_Mead_termination(ff, tol)) THEN
            DO i = 1, n
               x(i) = xx(1, i)
            END DO
            EXIT
         END IF

         ! Calculate centroid of 1...n (not n+1; the worst point) through which to reflect the worst point
         xo = Nelder_Mead_centroid(xx)

         ! Calculate the reflected point of n+1 about the centroid
         xr = xo+alpha*(xo-xx(n+1, :))
         fr = f(xr)

         IF (fr < ff(1)) THEN
            ! if the reflected point is the best so far then calculate the expanded point 
            xe = xo+gamma*(xr-xo)
            fe = f(xe)
            IF(fe < fr) THEN
               ! Keep the expansed point if it is better than the reflected ...
               operation = 'Expanding'
               xx(n+1, :) = xe
               ff(n+1) = fe
            ELSE
               ! ... otherwise keep the reflected.
               operation = 'Reflecting'
               xx(n+1, :) = xr
               ff(n+1) = fr
            END IF
            CYCLE
         ELSE IF (ff(1) <= fr .AND. fr < ff(n)) THEN
            ! If the reflected point is not the best, but better than the second worst then keep the reflected point
            operation = 'Reflecting'
            xx(n+1, :) = xr
            ff(n+1) = fr
            CYCLE
         ELSE
            ! Here it is certain that the reflected point is at least as worse than the second worst point
            ! Calculate the contracted point
            xc = xo+rho*(xx(n+1, :)-xo)
            fc = f(xc)
            IF(fc < ff(n+1)) THEN
               ! Keep the contracted point if it is better than the worst point
               operation = 'Contracting'
               xx(n+1, :) = xc
               ff(n+1) = fc
               CYCLE
            ELSE
               ! The contracted point is the worst, and we must shrink
               ! Calculate the shrinkage for all except the best point
               operation = 'Shrinking'
               DO i = 2, n+1
                  xx(i, :) = xx(1, :)+sigma*(xx(i, :)-xx(1, :))
                  ff(i) = f(xx(i, :))
               END DO
            END IF
         END IF

      END DO

      ! Report the minimization point
      x = xx(1, :)
      fom = ff(1)

      IF (verbose) THEN
         !WRITE(*, '(A16,I10,'//trim(nstring)//'F15.7)') 'Best', ii, fom, (x(1, i), i = 1, n)
         WRITE(*, *) 'NELDER_MEAD: Done'
         WRITE(*, *)
      END IF

   END SUBROUTINE Nelder_Mead

   FUNCTION Nelder_Mead_centroid(x) result(centroid)

      ! Calculate the centroid of all points except the worst point, which is n+1
      ! x Should be x(n+1, n)
      ! Output should be centroid(n)
      REAL, INTENT(IN) :: x(:, :)
      REAL :: centroid(size(x, 2))
      REAL, ALLOCATABLE :: xx(:, :)
      INTEGER :: i, n

      ! Check array sizes
      n = size(x, 2)
      IF(n+1 /= size(x, 1)) STOP 'NELDER_MEAD_CENTROID: Error array x is wrong'

      ! This ugly stuff is needed because the final element is not counted in the mean
      ALLOCATE(xx(n, n))
      DO i = 1, n
         xx(i, :) = x(i, :)
      END DO

      ! Calculate the mean
      DO i = 1, n
         centroid(i) = mean(xx(:, i))
      END DO

   END FUNCTION Nelder_Mead_centroid

   SUBROUTINE Nelder_Mead_sort(x, f)

      ! Sort the points into order from best to worst
      ! x should be x(n+1, n), f should be f(n+1)
      REAL, INTENT(INOUT) :: x(:, :)
      REAL, INTENT(INOUT) :: f(:)
      INTEGER :: i, n
      INTEGER, ALLOCATABLE :: j(:)
      INTEGER, PARAMETER :: isort = isort_Nelder_Mead

      ! Check array sizes
      n = size(x, 2)
      IF (n+1 /= size(x, 1) .OR. n+1 /= size(f)) THEN
         STOP 'NELDER_MEAD_SORT: Error, arrays are the wrong size'
      END IF
      ALLOCATE(j(n+1))

      ! Carry out the reindexing
      CALL index(f, j, isort)
      CALL reindex(f, j)
      DO i = 1, n
         CALL reindex(x(:,i), j)
      END DO

   END SUBROUTINE Nelder_Mead_sort

   LOGICAL FUNCTION Nelder_Mead_termination(f, tol)

      ! Determine if the minimization has converged
      ! f should be f(n+1)
      REAL, INTENT(IN) :: f(:)
      REAL, INTENT(IN) :: tol
      REAL :: sigma

      ! Calculate the standard deviation of all n+1 points
      sigma = standard_deviation(f)

      ! Decide on termination
      IF(sigma <= tol) THEN
         Nelder_Mead_termination = .TRUE.
      ELSE
         Nelder_Mead_termination = .FALSE.
      END IF

   END FUNCTION Nelder_Mead_termination

   SUBROUTINE adaptive_minimization_1D(xmin, func, min, max, ref, ngrid, nref)

      USE basic_operations
      REAL, INTENT(OUT) :: xmin
      REAL, EXTERNAL :: func
      REAL, INTENT(IN) :: min
      REAL, INTENT(IN) :: max
      REAL, INTENT(IN) :: ref
      INTEGER, INTENT(IN) :: ngrid
      INTEGER, INTENT(IN) :: nref
      REAL :: fom, fommin, range, x
      REAL :: x1, x2
      INTEGER :: i, j, imin

      INTERFACE
         FUNCTION func(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION func
      END INTERFACE

      !xmin is output
      !min/max give you the initial bounds for the search
      !ref is the refinement level for each separate grid search
      !n is the number of points per grid search
      !m is the number of refinement levels

      !Needs to be set to prevent warning
      imin = 1

      !Needs to be set to prevent warnings
      fommin = 0.

      x1 = min
      x2 = max

      DO j = 1, nref

         DO i = 1, ngrid

            x = progression(x1, x2, i, ngrid)
            fom = func(x)

            IF (i == 1 .OR. fom < fommin) THEN
               fommin = fom
               xmin = x
               imin = i
            END IF

         END DO

         IF (imin == 1 .OR. imin == ngrid) STOP 'Error: best fit at edge of range'

         IF (j .NE. nref) THEN
            range = x2-x1
            x1 = xmin-range/ref
            x2 = xmin+range/ref
         END IF

      END DO

   END SUBROUTINE adaptive_minimization_1D

   SUBROUTINE adaptive_minimization_3D(xmin, func, min, max, ref, ngrid, nref)

      USE basic_operations
      REAL, INTENT(OUT) :: xmin(3)
      REAL, INTENT(IN) :: min(3), max(3), ref(3)
      INTEGER, INTENT(IN) :: ngrid(3), nref
      REAL, EXTERNAL :: func
      REAL :: fom, fommin, range(3), x(3)
      REAL :: x1(3), x2(3)
      INTEGER :: i, j, k, l, imin, jmin, kmin

      !xmin is output
      !min/max give you the initial bounds for the search
      !ref is the refinement level for each separate grid search
      !n is the number of points per grid search
      !m is the number of refinement levels

      !Need to be set to prevent warnings
      imin = 1
      jmin = 1
      kmin = 1

      !Needs to be set to prevent warnings
      fommin = 0.

      x1 = min
      x2 = max

      DO l = 1, nref

         DO i = 1, ngrid(1)

            !x(1)=x1(1)+(x2(1)-x1(1))*float(i-1)/float(ngrid(1)-1)
            x(1) = progression(x1(1), x2(1), i, ngrid(1))

            DO j = 1, ngrid(2)

               !x(2)=x1(2)+(x2(2)-x1(2))*float(j-1)/float(ngrid(2)-1)
               x(2) = progression(x1(2), x2(2), j, ngrid(2))

               DO k = 1, ngrid(3)

                  !x(3)=x1(3)+(x2(3)-x1(3))*float(k-1)/float(ngrid(3)-1)
                  x(3) = progression(x1(3), x2(3), k, ngrid(3))

                  fom = func(x)

                  IF ((i == 1 .AND. j == 1 .AND. k == 1) .OR. fom < fommin) THEN
                     fommin = fom
                     xmin = x
                     imin = i
                     jmin = j
                     kmin = k
                  END IF

               END DO

            END DO

         END DO

         IF (imin == 1 .OR. imin == ngrid(1)) STOP 'Error: best fit x(1) at edge of range'
         IF (jmin == 1 .OR. jmin == ngrid(2)) STOP 'Error: best fit x(2) at edge of range'
         IF (kmin == 1 .OR. kmin == ngrid(3)) STOP 'Error: best fit x(3) at edge of range'

         IF (l .NE. nref) THEN
            range = x2-x1
            x1 = xmin-range/ref
            x2 = xmin+range/ref
         END IF

      END DO

   END SUBROUTINE adaptive_minimization_3D

   SUBROUTINE grid_minimization_1D(Amin, Amax, Asteps, x, y)

      !This fits a one parameter model
      USE basic_operations
      REAL, INTENT(IN) :: Amin, Amax
      REAL, INTENT(IN) :: x(:), y(:)
      INTEGER, INTENT(IN) :: Asteps
      REAL :: ymod, fom, fombest, A, Abest
      INTEGER :: i, j

      fombest = 10000000.

      DO j = 1, Asteps

         !A=Amin+(Amax-Amin)*(float(j-1)/(Asteps-1))
         A = progression(Amin, Amax, j, Asteps)

         fom = 0.

         DO i = 1, size(x)

            !This is the model that you wish to fit
            ymod = x(i)**A

            fom = fom+(ymod-y(i))**2.

         END DO

         IF (fom < fombest) THEN
            fombest = fom
            Abest = A
         END IF

      END DO

      WRITE (*, *) 'Best fit is:', Abest

   END SUBROUTINE grid_minimization_1D

!!$  SUBROUTINE grid_fit_2(Amin,Amax,Asteps,Bmin,Bmax,Bsteps,x,y)
!!$
!!$    IMPLICIT NONE
!!$
!!$    REAL, INTENT(IN) :: Amin, Amax, Bmin, Bmax
!!$    REAL, INTENT(IN) :: x(:), y(:)
!!$    INTEGER, INTENT(IN) :: Asteps, Bsteps
!!$    REAL :: fom, fombest, A, Abest, B, Bbest, ymod
!!$    INTEGER :: i, j, k
!!$
!!$    !This fits a two parameter model
!!$
!!$    fombest=10000000.
!!$
!!$    DO j=1,Asteps
!!$
!!$       A=Amin+(Amax-Amin)*(float(j-1)/(float(Asteps-1)))
!!$
!!$       DO k=1,Bsteps
!!$
!!$          B=Bmin+(Bmax-Bmin)*(float(k-1)/(float(Bsteps-1)))
!!$
!!$          fom=0.
!!$          DO i=1,size(x)
!!$
!!$             !This is the model that you wish to fit
!!$             ymod=model(A,B,x(i))
!!$
!!$             !This looks at the raw difference between models!
!!$             fom=fom+(-1.+ymod/y(i))**2.
!!$
!!$             !         WRITE(*,*) y(i), A, B, ymod, fom
!!$
!!$          END DO
!!$
!!$          IF(fom<fombest) THEN
!!$             fombest=fom
!!$             Abest=A
!!$             Bbest=B
!!$          END IF
!!$
!!$       END DO
!!$    END DO
!!$
!!$    WRITE(*,*) 'Best fit A is:', Abest
!!$    WRITE(*,*) 'Best fit B is:', Bbest
!!$
!!$    OPEN(8,file='results.dat')
!!$    DO i=1,size(x)
!!$       WRITE(8,*) x(i), y(i), model(Abest,Bbest,x(i))
!!$    END DO
!!$    CLOSE(8)
!!$
!!$  END SUBROUTINE grid_fit_2

!!$  SUBROUTINE grid_fit_3(Amin,Amax,Asteps,Bmin,Bmax,Bsteps,Cmin,Cmax,Csteps,x,y)
!!$
!!$    IMPLICIT NONE
!!$
!!$    REAL, INTENT(IN) :: Amin, Amax, Bmin, Bmax, Cmin, Cmax
!!$    REAL, INTENT(IN) :: x(:), y(:)
!!$    INTEGER, INTENT(IN) :: Asteps, Bsteps, Csteps
!!$    REAL :: fom, fombest, A, Abest, B, Bbest, C, Cbest, ymod
!!$    INTEGER :: i, j, k, l
!!$
!!$    !This fits a two parameter model
!!$
!!$    fombest=10000000.
!!$
!!$    DO j=1,Asteps
!!$
!!$       IF(Asteps==1) THEN
!!$          A=Amin
!!$       ELSE
!!$          A=Amin+(Amax-Amin)*(float(j-1)/(Asteps-1))
!!$       END IF
!!$
!!$       DO k=1,Bsteps
!!$
!!$          IF(Bsteps==1) THEN
!!$             B=Bmin
!!$          ELSE
!!$             B=Bmin+(Bmax-Bmin)*(float(k-1)/(Bsteps-1))
!!$          END IF
!!$
!!$          DO l=1,Csteps
!!$
!!$             IF(Csteps==1) THEN
!!$                C=Cmin
!!$             ELSE
!!$                C=Cmin+(Cmax-Cmin)*(float(l-1)/(Csteps-1))
!!$             END IF
!!$
!!$             fom=0.
!!$             DO i=1,size(x)
!!$
!!$                !This is the model that you wish to fit
!!$                ymod=model(A,B,C,x(i))
!!$
!!$                !This looks at the raw difference between models!
!!$                fom=fom+(-1.+ymod/y(i))**2.
!!$
!!$                !         WRITE(*,*) y(i), A, B, ymod, fom
!!$
!!$             END DO
!!$
!!$             !         WRITE(10,*) A, B, C, fom
!!$
!!$             IF(fom<fombest) THEN
!!$                fombest=fom
!!$                Abest=A
!!$                Bbest=B
!!$                Cbest=C
!!$             END IF
!!$
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    WRITE(*,*) 'Best fit A is:', Abest
!!$    WRITE(*,*) 'Best fit B is:', Bbest
!!$    WRITE(*,*) 'Best fit C is:', Cbest
!!$
!!$    OPEN(8,file='results.dat')
!!$    DO i=1,size(x)
!!$       WRITE(8,*) x(i), y(i), model(Abest,Bbest,Cbest,x(i))
!!$    END DO
!!$    CLOSE(8)
!!$
!!$  END SUBROUTINE grid_fit_3

   REAL FUNCTION gradient_minimization_1D(f, x, acc)

      ! Finds the minimum of the function 'f' starting at value x0 with accuracy acc!
      USE calculus
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: acc
      REAL :: xold, xnew, df
      INTEGER :: i

      INTEGER, PARAMETER :: n = 1000
      REAL, PARAMETER :: fac = 1e3

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION f
      END INTERFACE

      STOP 'GRADIENT_FIT_1: This does not work!'

      xold = x

      DO i = 1, n

         !       df=derivative(f,xold,acc/100.)
         !       xnew=xold-df*acc*100.

         WRITE (*, *) i, xold, xnew

         df = derivative(f, xold, acc/fac)
         xnew = xold-df*acc*fac

         WRITE (*, *) i, xold, xnew

         IF (i > 1 .AND. abs(xnew/xold-1.) < acc) THEN
            gradient_minimization_1D = xnew
            EXIT
         ELSE
            xold = xnew
         END IF

      END DO

   END FUNCTION gradient_minimization_1D

   FUNCTION gradient_minimization_2D(f, x, y, acc)

      ! Finds the minimum of the function 'f' starting at value x0 with accuracy acc
      USE calculus
      REAL :: gradient_minimization_2D(2)
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: acc
      REAL :: xold, xnew, yold, ynew, dfx, dfy
      INTEGER :: i

      INTERFACE
         FUNCTION f(xin, yin)
            REAL, INTENT(IN) :: xin
            REAL, INTENT(IN) :: yin
         END FUNCTION f
      END INTERFACE

      STOP 'GRADIENT_FIT_2: This does not work!'

      xold = x
      yold = y

      DO i = 1, 1000

         dfx = derivative(f, xold, yold, acc/100., dim=1)
         dfy = derivative(f, xold, yold, acc/100., dim=2)

         xnew = xold-dfx*acc*100.
         ynew = yold-dfy*acc*100.

         IF (i > 1 .AND. abs(xnew/xold-1.) < acc .AND. abs(ynew/yold-1.) < acc) THEN
            gradient_minimization_2D(1) = xnew
            gradient_minimization_2D(2) = ynew
            EXIT
         ELSE
            xold = xnew
            yold = ynew
         END IF

      END DO

   END FUNCTION gradient_minimization_2D

   REAL FUNCTION quadratic_extremum_1D(x1, y1, x2, y2, x3, y3)

      ! This calculates the extrema of a function under the assumption that it is quadratic
      ! Takes 3 points to form a quadratic and then read off minimum
      USE special_functions
      REAL, INTENT(IN) :: x1, x2, x3, y1, y2, y3
      REAL :: a2, a1, a0

      CALL fix_polynomial(a2, a1, a0, [x1, x2, x3], [y1, y2, y3])
      quadratic_extremum_1D = -a1/(2.*a2)

   END FUNCTION quadratic_extremum_1D

   REAL FUNCTION find_array_maximum(x, y)

      ! From an array y(x) finds the x location of the maximum treating y(x) as a continuous function
      USE special_functions
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: y(:)
      REAL :: x1, x2, x3, y1, y2, y3
      INTEGER :: imax(1), i, n

      n = size(x)
      IF(n /= size(y)) STOP 'MAXIMUM: Error, x and y must be the same size'

      ! Need this to stop a compile-time warning
      find_array_maximum = 0.

      ! Integer maximum location
      imax = maxloc(y)
      i = imax(1)

      IF (i == 1 .OR. i == n) THEN

         STOP 'MAXIMUM: Error, maximum array value is at one end of the array'

      ELSE

         ! Get the x positions
         x1 = x(i-1)
         x2 = x(i)
         x3 = x(i+1)

         ! Get the y values
         y1 = y(i-1)
         y2 = y(i)
         y3 = y(i+1)

         find_array_maximum = quadratic_extremum_1D(x1, y1, x2, y2, x3, y3)

      END IF

   END FUNCTION find_array_maximum

END MODULE minimization
