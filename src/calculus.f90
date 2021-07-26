MODULE calculus

   ! TODO: Make jmin -> jmax parameters in header
   ! TODO: Number of steps is 2^(j-1) which may exceed integer range for large j

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: derivative
   PUBLIC :: integrate
   PUBLIC :: integrate_basic
   PUBLIC :: integrate_log
   PUBLIC :: integrate_cubic
   PUBLIC :: integrate_jac
   PUBLIC :: integrate_monte_carlo

   ! Differentiation
   INTEGER, PARAMETER :: min_deriv_iter = 5
   INTEGER, PARAMETER :: max_deriv_iter = 100

   ! Integration
   INTEGER, PARAMETER :: jmin_integrate = 5  ! 2^jmin minimum steps
   INTEGER, PARAMETER :: jmax_integrate = 20 ! 2^jmax maximum steps
   INTEGER, PARAMETER :: ninit_integrate = 5

   INTERFACE derivative
      MODULE PROCEDURE derivative_1D
      MODULE PROCEDURE derivative_2D
   END INTERFACE derivative

CONTAINS

   REAL FUNCTION derivative_1D(f, x, acc)

      ! Calculates the derivative of a function 'f' at the point x to accuracy acc!
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: acc
      REAL :: dnew, dold, dx
      INTEGER :: i
      INTEGER, PARAMETER :: imin = min_deriv_iter ! Minimum number of iterations
      INTEGER, PARAMETER :: n = max_deriv_iter ! Maximum number of iterations

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION f
      END INTERFACE

      dold = 0.
      dx = 1. ! Is this a good choice?

      DO i = 1, n

         dnew = (f(x+dx/2.)-f(x-dx/2.))/dx ! New, using equal sided derivative

         IF (i >= imin .AND. abs(dnew/dold-1.) < acc) THEN
            EXIT
         ELSE IF (i == n) THEN
            STOP 'DERIVATIVE: Error, maximum number of iterations exceeded'
         ELSE
            dold = dnew
            dx = dx/2.
         END IF

      END DO

      derivative_1D = dnew

   END FUNCTION derivative_1D

   REAL FUNCTION derivative_2D(f, x, y, acc, dim)

      ! Calculates the derivative of a function 'f' at the point x, y along dimension dim to accuracy acc!
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: dim
      REAL :: dnew, dold, dx, dy
      INTEGER :: i
      INTEGER, PARAMETER :: n = max_deriv_iter

      INTERFACE
         FUNCTION f(xin, yin)
            REAL, INTENT(IN) :: xin
            REAL, INTENT(IN) :: yin
         END FUNCTION f
      END INTERFACE

      dold = 0.
      dx = 4.
      dy = 4.

      DO i = 1, n

         IF (dim == 1) THEN

            dnew = (f(x+dx, y)-f(x, y))/dx

            IF (i > 1 .AND. abs(dnew/dold-1.) < acc) THEN
               EXIT
            ELSE IF (i == n) THEN
               STOP 'DERIVATIVE_X: Error, maximum number of iterations exceeded'
            ELSE
               dold = dnew
               dx = dx/2.
            END IF

         ELSE IF (dim == 2) THEN

            dnew = (f(x, y+dy)-f(x, y))/dy

            IF (i > 1 .AND. abs(dnew/dold-1.) < acc) THEN
               EXIT
            ELSE IF (i == n) THEN
               STOP 'DERIVATIVE_Y: Error, maximum number of iterations exceeded'
            ELSE
               dold = dnew
               dy = dy/2.
            END IF

         ELSE

            STOP 'DERIVATIVE_2D: Error, dim specified incorrectly'

         END IF

      END DO

      derivative_2D = dnew

   END FUNCTION derivative_2D

   REAL FUNCTION derivative_y(f, x, y, acc)

      ! Calculates the derivative of a function 'f' at the point x to accuracy acc!
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: acc
      REAL :: dnew, dold, dy
      INTEGER :: i
      INTEGER, PARAMETER :: n = max_deriv_iter

      INTERFACE
         FUNCTION f(xin, yin)
            REAL, INTENT(IN) :: xin
            REAL, INTENT(IN) :: yin
         END FUNCTION f
      END INTERFACE

      dy = 4.
      dold = 0.

      DO i = 1, n

         dnew = (f(x, y+dy)-f(x, y))/dy

         IF (i > 1 .AND. abs(dnew/dold-1.) < acc) THEN
            !derivative_y=dnew
            EXIT
         ELSE IF (i == n) THEN
            STOP 'DERIVATIVE_Y: Error, maximum number of iterations exceeded'
         ELSE
            dold = dnew
            dy = dy/2.
         END IF

      END DO

      derivative_y = dnew

   END FUNCTION derivative_y

   REAL FUNCTION integrate_basic(a, b, f, n, iorder)

      ! Integrates between a and b with n points; not adaptive so no error control
      USE basic_operations
      REAL, INTENT(IN) :: a ! Integration limits
      REAL, INTENT(IN) :: b ! Integration limits
      REAL, EXTERNAL :: f
      INTEGER, INTENT(IN) :: n ! Number of points
      INTEGER, INTENT(IN) :: iorder ! Order for integration
      INTEGER :: i
      REAL :: x, dx, weight
      REAL :: sum

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         integrate_basic = 0.

      ELSE

         ! Set the sum variable
         sum = 0.d0

         DO i = 1, n

            !x=a+(b-a)*real(i-1)/real(n-1)
            x = progression(a, b, i, n)

            IF (iorder == 1) THEN
               ! Composite trapezium weights
               IF (i == 1 .OR. i == n) THEN
                  weight = 0.5
               ELSE
                  weight = 1.
               END IF
            ELSE IF (iorder == 2) THEN
               ! Composite extended formula weights
               IF (i == 1 .OR. i == n) THEN
                  weight = 0.4166666666
               ELSE IF (i == 2 .OR. i == n-1) THEN
                  weight = 1.0833333333
               ELSE
                  weight = 1.
               END IF
            ELSE IF (iorder == 3) THEN
               ! Composite Simpson weights
               IF (i == 1 .OR. i == n) THEN
                  weight = 0.375
               ELSE IF (i == 2 .OR. i == n-1) THEN
                  weight = 1.1666666666
               ELSE IF (i == 3 .OR. i == n-2) THEN
                  weight = 0.9583333333
               ELSE
                  weight = 1.
               END IF
            ELSE
               STOP 'INTEGERATE_BASIC: Error, order specified incorrectly'
            END IF

            sum = sum+weight*f(x)

         END DO

         dx = (b-a)/real(n-1)
         integrate_basic = real(sum)*dx

         !WRITE(*,*) 'INTEGRATE_BASIC: Order:', iorder
         !WRITE(*,*) 'INTEGRATE_BASIC: Nint:', n

      END IF

   END FUNCTION integrate_basic

   REAL FUNCTION integrate(a, b, f, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      USE basic_operations
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      LOGICAL :: pass
      REAL :: sum_n, sum_2n, sum_new, sum_old
      INTEGER, PARAMETER :: jmin = jmin_integrate
      INTEGER, PARAMETER :: jmax = jmax_integrate

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.d0
         sum_n = 0.d0
         sum_old = 0.d0
         sum_new = 0.d0

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a)
               f2 = f(b)
               sum_2n = 0.5*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = progression(a, b, i, n)
                  fx = f(x)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF(abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO

         integrate = real(sum_new)

      END IF

   END FUNCTION integrate

   REAL FUNCTION integrate_log(a, b, f, acc, iorder, ilog)

      ! Integrates between a and b until desired accuracy is reached
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ilog
      INTEGER :: i, j, n
      REAL :: x, weight, dx, lima, limb
      LOGICAL :: pass
      REAL :: sum_old, sum_new
      INTEGER, PARAMETER :: jmin = jmin_integrate
      INTEGER, PARAMETER :: jmax = jmax_integrate
      INTEGER, PARAMETER :: ninit = ninit_integrate

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         integrate_log = 0.d0

      ELSE

         ! Set the sum variables
         sum_old = 0.d0
         sum_new = 0.d0

         IF (ilog == 1) THEN
            lima = log(a)
            limb = log(b)
         ELSE
            lima = a
            limb = b
         END IF

         DO j = 1, jmax

            n = ninit*(2**(j-1))

            DO i = 1, n

               x = lima+(limb-lima)*real(i-1)/real(n-1)

               IF (ilog == 1) THEN
                  x = exp(x)
               END IF

               IF (iorder == 1) THEN
                  ! Composite trapezium weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.5
                  ELSE
                     weight = 1.
                  END IF
               ELSE IF (iorder == 2) THEN
                  ! Composite extended formula weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.4166666666
                  ELSE IF (i == 2 .OR. i == n-1) THEN
                     weight = 1.0833333333
                  ELSE
                     weight = 1.
                  END IF
               ELSE IF (iorder == 3) THEN
                  ! Composite Simpson weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.375
                  ELSE IF (i == 2 .OR. i == n-1) THEN
                     weight = 1.1666666666
                  ELSE IF (i == 3 .OR. i == n-2) THEN
                     weight = 0.9583333333
                  ELSE
                     weight = 1.
                  END IF
               ELSE
                  STOP 'INTEGERATE_LOG: Error, order specified incorrectly'
               END IF

               IF (ilog == 0) THEN
                  sum_new = sum_new+weight*f(x)
               ELSE IF (ilog == 1) THEN
                  sum_new = sum_new+weight*f(x)*x
               ELSE
                  STOP 'INTEGRATE_LOG: Error, ilog specified incorrectly'
               END IF

            END DO

            dx = (limb-lima)/real(n-1)
            sum_new = sum_new*dx

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF(abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_LOG: Integration timed out'       
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_new = 0.
            END IF

         END DO
         
         integrate_log = real(sum_new)

      END IF

   END FUNCTION integrate_log

   REAL FUNCTION integrate_cubic(a, b, f, acc)

      ! Integrates between a and b until desired accuracy is reached!
      ! Fits a cubic between successive 4 points
      ! Only useful if points are not eqaully spaced, thus this routine is probably redundant
      USE special_functions
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: acc
      INTEGER :: i, j, nint, nsec
      REAL :: a3, a2, a1, a0
      REAL :: x1, x2, x3, x4
      REAL :: y1, y2, y3, y4
      REAL :: sum_old, sum_new
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = jmin_integrate
      INTEGER, PARAMETER :: jmax = jmax_integrate
      INTEGER, PARAMETER :: ni = ninit_integrate  ! CARE: This was previously set to 1

      INTERFACE
         FUNCTION f(x)
            REAL, INTENT(IN) :: x
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         integrate_cubic = 0.

      ELSE

         ! Set the sum variables
         sum_old = 0.d0
         sum_new = 0.d0

         DO j = 1, jmax

            ! This is the number of cubic sections (each of which has four function evaluations)
            nsec = ni*2**(j-1)

            ! Number of function evaluation points so as to be able to fit a cubic (4,7,10 ...)
            nint = 3*nsec+1

            DO i = 1, nsec

               IF (i == 1) THEN

                  x1 = a+(b-a)*float(3*(i-1)+1-1)/float(nint-1)
                  y1 = f(x1)

               ELSE

                  x1 = x4
                  y1 = y4

               END IF

               x2 = a+(b-a)*float(3*(i-1)+2-1)/float(nint-1)
               x3 = a+(b-a)*float(3*(i-1)+3-1)/float(nint-1)
               x4 = a+(b-a)*float(3*(i-1)+4-1)/float(nint-1)

               y2 = f(x2)
               y3 = f(x3)
               y4 = f(x4)

               CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [y1, y2, y3, y4])

               ! Add the (analytical) intergal of a cubic between points x1 and x4 to the total
               sum_new = sum_new+(a3/4.)*(x4**4.-x1**4.)+(a2/3.)*(x4**3.-x1**3.)+(a1/2.)*(x4**2.-x1**2.)+a0*(x4-x1)

            END DO

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF(abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_CUBIC: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_new = 0.
            END IF

         END DO

         integrate_cubic = real(sum_new)

      END IF

   END FUNCTION integrate_cubic

   REAL FUNCTION integrate_jac(a, b, f, acc, iorder, g, gi, dg)

      ! Integrates between a and b until desired accuracy is reached
      ! Uses a Jacobian to speed up the integration
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      REAL, EXTERNAL :: g
      REAL, EXTERNAL :: gi
      REAL, EXTERNAL :: dg
      INTEGER :: i, j, n
      REAL :: dy, alim, blim
      REAL :: x, y, weight
      REAL :: sum_old, sum_new
      LOGICAL :: pass     

      INTEGER, PARAMETER :: jmin = jmin_integrate
      INTEGER, PARAMETER :: jmax = jmax_integrate
      INTEGER, PARAMETER :: ninit = ninit_integrate

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION f
         FUNCTION g(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION g
         FUNCTION gi(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION gi
         FUNCTION dg(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION dg
      END INTERFACE

      IF (a == b) THEN

         integrate_jac = 0.

      ELSE

         ! Set the sum variables
         sum_old = 0.d0
         sum_new = 0.d0

         alim = g(a)
         blim = g(b)

         DO j = 1, jmax

            n = ninit*(2**(j-1))

            DO i = 1, n

               y = alim+(blim-alim)*real(i-1)/real(n-1)

               IF (iorder == 1) THEN
                  ! Composite trapezium weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.5
                  ELSE
                     weight = 1.
                  END IF
               ELSE IF (iorder == 2) THEN
                  ! Composite extended formula weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.4166666666
                  ELSE IF (i == 2 .OR. i == n-1) THEN
                     weight = 1.0833333333
                  ELSE
                     weight = 1.
                  END IF
               ELSE IF (iorder == 3) THEN
                  ! Composite Simpson weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.375
                  ELSE IF (i == 2 .OR. i == n-1) THEN
                     weight = 1.1666666666
                  ELSE IF (i == 3 .OR. i == n-2) THEN
                     weight = 0.9583333333
                  ELSE
                     weight = 1.
                  END IF
               ELSE
                  STOP 'INTEGERATE_JAC: Error, order specified incorrectly'
               END IF

               x = gi(y)

               sum_new = sum_new+weight*f(x)/dg(x)

            END DO

            dy = (blim-alim)/real(n-1)
            sum_new = sum_new*dy

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF(abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_JAC: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF(pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_new = 0.
            END IF

         END DO
         
         integrate_jac = real(sum_new)

      END IF

   END FUNCTION integrate_jac

   REAL FUNCTION integrate_monte_carlo(a, b, f, n)

      ! Integrates between a and b with n points; not adaptive so no error control
      USE random_numbers
      REAL, INTENT(IN) :: a ! Integration limits
      REAL, INTENT(IN) :: b ! Integration limits
      REAL, EXTERNAL :: f
      INTEGER, INTENT(IN) :: n ! Number of points
      INTEGER :: i
      REAL :: x, dx
      REAL :: sum

      INTERFACE
         FUNCTION f(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         integrate_monte_carlo = 0.

      ELSE

         sum = 0.d0
         DO i = 1, n
            x = random_uniform(a, b)
            sum = sum+f(x)
         END DO

         dx = (b-a)/real(n)

         sum = sum*dx

         integrate_monte_carlo = real(sum)

      END IF

   END FUNCTION integrate_monte_carlo

END MODULE calculus
