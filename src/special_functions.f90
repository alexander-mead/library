MODULE special_functions

   USE precision
   USE constants
   USE basic_operations

   IMPLICIT NONE

   PRIVATE

   ! Integer functions
   PUBLIC :: triangle_number
   PUBLIC :: factorial
   PUBLIC :: get_factorials
   PUBLIC :: Fibonacci
   PUBLIC :: get_Fibonaccis
   PUBLIC :: binomial_coefficient

   ! Real functions
   PUBLIC :: polynomial
   PUBLIC :: fix_polynomial
   PUBLIC :: centred_polynomial
   PUBLIC :: fix_centred_polynomial
   PUBLIC :: Legendre_polynomial
   PUBLIC :: Lagrange_polynomial
   PUBLIC :: Si
   PUBLIC :: Ci
   PUBLIC :: Bessel
   PUBLIC :: spherical_Bessel
   PUBLIC :: sinc
   PUBLIC :: wk_tophat
   PUBLIC :: wk_tophat_deriv
   PUBLIC :: wk_tophat_dderiv
   PUBLIC :: cbrt
   PUBLIC :: Heaviside

   ! Minima testing
   PUBLIC :: Rosenbrock
   PUBLIC :: Himmelblau

   ! Probability distributions
   PUBLIC :: Gaussian_distribution
   PUBLIC :: lognormal_distribution
   PUBLIC :: uniform_distribution
   PUBLIC :: Rayleigh_distribution
   PUBLIC :: exponential_distribution
   PUBLIC :: Lorentzian_distribution
   PUBLIC :: polynomial_distribution

   ! Integer probability distributions
   PUBLIC :: Poisson_distribution
   PUBLIC :: Bernoulli_distribution
   PUBLIC :: geometric_distribution
   PUBLIC :: shifted_geometric_distribution
   PUBLIC :: Binomial_distribution

   ! Complex functions
   PUBLIC :: complex_number
   PUBLIC :: complex_phase

   ! Silly functions
   PUBLIC :: apodise
   PUBLIC :: smooth_apodise
   PUBLIC :: blob
   PUBLIC :: smooth_blob
   PUBLIC :: sigmoid_tanh
   PUBLIC :: sigmoid_log

   ! Taylor expansion below this
   REAL, PARAMETER :: dx_sinc = 1e-3
   REAL, PARAMETER :: dx_tophat = 1e-3
   REAL, PARAMETER :: dx_Bessel = 1e-3

   ! Numerical approximation parameters
   REAL, PARAMETER :: x0_SiCi = 4.
   REAL, PARAMETER :: xbig_Bessel = 1e15

   INTERFACE polynomial
      MODULE PROCEDURE linear_polynomial
      MODULE PROCEDURE quadratic_polynomial
      MODULE PROCEDURE cubic_polynomial
   END INTERFACE polynomial

   INTERFACE centred_polynomial
      MODULE PROCEDURE centred_linear_polynomial
      MODULE PROCEDURE centred_quadratic_polynomial
      MODULE PROCEDURE centred_cubic_polynomial
   END INTERFACE centred_polynomial

   INTERFACE fix_polynomial
      MODULE PROCEDURE fix_linear_polynomial
      MODULE PROCEDURE fix_quadratic_polynomial
      MODULE PROCEDURE fix_cubic_polynomial
   END INTERFACE fix_polynomial

   INTERFACE fix_centred_polynomial
      MODULE PROCEDURE fix_centred_linear_polynomial
      MODULE PROCEDURE fix_centred_quadratic_polynomial
      MODULE PROCEDURE fix_centred_cubic_polynomial
   END INTERFACE fix_centred_polynomial

   INTERFACE cbrt
      MODULE PROCEDURE cbrt_real
      MODULE PROCEDURE cbrt_int
   END INTERFACE cbrt

CONTAINS

   REAL FUNCTION Heaviside(x, Hzero_opt)

      REAL, INTENT(IN) :: x
      REAL, OPTIONAL, INTENT(IN) :: Hzero_opt
      REAL, PARAMETER :: Hzero_def = 0.5

      IF (x == 0.) THEN
         Heaviside = default_or_optional(Hzero_def, Hzero_opt)
      ELSE IF (x < 0.) THEN
         Heaviside = 0.
      ELSE
         Heaviside = 1.
      END IF

   END FUNCTION Heaviside

   COMPLEX FUNCTION complex_number(r, theta)

      ! Complex number r*e^{i theta} in Fortran format
      ! Analogy of inbuilt cmplx(x, y) function
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: theta

      complex_number = r*cos(theta)+(0., 1.)*sin(theta)

   END FUNCTION complex_number

   REAL FUNCTION complex_phase(z)
            
      ! Phase of complex number z in radians
      ! Result lies in the range -pi to pi
      ! Commented-out bit would return result in 0 to 2pi
      COMPLEX, INTENT(IN) :: z

      complex_phase = atan2(aimag(z), real(z))
      !IF (complex_phase < 0.) complex_phase = twopi-complex_phase

   END FUNCTION complex_phase

   ELEMENTAL REAL FUNCTION cbrt_real(x)

      ! Cube root, analogy of sqrt()
      REAL, INTENT(IN) :: x

      cbrt_real = x**(1./3.)

   END FUNCTION cbrt_real

   ELEMENTAL REAL FUNCTION cbrt_int(x)

      ! Cube root, analogy of sqrt()
      INTEGER, INTENT(IN) :: x

      cbrt_int = x**(1./3.)

   END FUNCTION cbrt_int

   ELEMENTAL REAL FUNCTION linear_polynomial(x, a1, a0)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a0

      linear_polynomial = a1*x+a0

   END FUNCTION linear_polynomial

   ELEMENTAL REAL FUNCTION quadratic_polynomial(x, a2, a1, a0)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: a2
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a0

      quadratic_polynomial = a2*x**2+a1*x+a0

   END FUNCTION quadratic_polynomial

   ELEMENTAL REAL FUNCTION cubic_polynomial(x, a3, a2, a1, a0)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: a3
      REAL, INTENT(IN) :: a2
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a0

      cubic_polynomial = a3*x**3+a2*x**2+a1*x+a0

   END FUNCTION cubic_polynomial

   ELEMENTAL REAL FUNCTION centred_linear_polynomial(x, x0, a1, a0)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x0
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a0
      REAL :: xx

      xx = x-x0
      centred_linear_polynomial = linear_polynomial(xx, a1, a0)

   END FUNCTION centred_linear_polynomial

   ELEMENTAL REAL FUNCTION centred_quadratic_polynomial(x, x0, a2, a1, a0)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x0
      REAL, INTENT(IN) :: a2
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a0
      REAL :: xx

      xx = x-x0
      centred_quadratic_polynomial = quadratic_polynomial(xx, a2, a1, a0)

   END FUNCTION centred_quadratic_polynomial

   ELEMENTAL REAL FUNCTION centred_cubic_polynomial(x, x0, a3, a2, a1, a0)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x0
      REAL, INTENT(IN) :: a3
      REAL, INTENT(IN) :: a2
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a0
      REAL :: xx

      xx = x-x0
      centred_cubic_polynomial = cubic_polynomial(xx, a3, a2,  a1, a0)

   END FUNCTION centred_cubic_polynomial

   SUBROUTINE fix_linear_polynomial(a1, a0, x, y)

      ! Given xi, yi i=1,2 fixes a line between these points
      REAL, INTENT(OUT) :: a1
      REAL, INTENT(OUT) :: a0 
      REAL, INTENT(IN) :: x(2)
      REAL, INTENT(IN) :: y(2)

      a1 = (y(2)-y(1))/(x(2)-x(1))
      a0 = y(1)-a1*x(1)

   END SUBROUTINE fix_linear_polynomial

   SUBROUTINE fix_quadratic_polynomial(a2, a1, a0, x, y)

      ! Given xi, yi i=1,2,3 fixes a quadratic between these points
      REAL, INTENT(OUT) :: a2
      REAL, INTENT(OUT) :: a1
      REAL, INTENT(OUT) :: a0
      REAL, INTENT(IN) :: x(3)
      REAL, INTENT(IN) :: y(3)

      a2 = ((y(2)-y(1))/(x(2)-x(1))-(y(3)-y(1))/(x(3)-x(1)))/(x(2)-x(3))
      a1 = (y(2)-y(1))/(x(2)-x(1))-a2*(x(2)+x(1))
      a0 = y(1)-a2*x(1)**2-a1*x(1)

   END SUBROUTINE fix_quadratic_polynomial

   SUBROUTINE fix_cubic_polynomial(a3, a2, a1, a0, x, y)

      ! Given xi, yi i=1,2,3,4 fixes a cubic between these points
      REAL, INTENT(OUT) :: a3
      REAL, INTENT(OUT) :: a2
      REAL, INTENT(OUT) :: a1
      REAL, INTENT(OUT) :: a0
      REAL, INTENT(IN) :: x(4)
      REAL, INTENT(IN) :: y(4)
      REAL :: f1, f2, f3

      f1 = (y(4)-y(1))/((x(4)-x(2))*(x(4)-x(1))*(x(4)-x(3)))
      f2 = (y(3)-y(1))/((x(3)-x(2))*(x(3)-x(1))*(x(4)-x(3)))
      f3 = (y(2)-y(1))/((x(2)-x(1))*(x(4)-x(3)))*(1./(x(4)-x(2))-1./(x(3)-x(2)))

      a3 = f1-f2-f3

      f1 = (y(3)-y(1))/((x(3)-x(2))*(x(3)-x(1)))
      f2 = (y(2)-y(1))/((x(2)-x(1))*(x(3)-x(2)))
      f3 = a3*(x(3)+x(2)+x(1))

      a2 = f1-f2-f3

      f1 = (y(4)-y(1))/(x(4)-x(1))
      f2 = a3*(x(4)**2+x(4)*x(1)+x(1)**2)
      f3 = a2*(x(4)+x(1))

      a1 = f1-f2-f3

      a0 = y(1)-a3*x(1)**3-a2*x(1)**2-a1*x(1)

   END SUBROUTINE fix_cubic_polynomial

   SUBROUTINE fix_centred_linear_polynomial(a1, a0, x0, x, y)

      REAL, INTENT(OUT) :: a1
      REAL, INTENT(OUT) :: a0
      REAL, INTENT(IN) :: x0
      REAL, INTENT(IN) :: x(2)
      REAL, INTENT(IN) :: y(2)

      CALL fix_linear_polynomial(a1, a0, x-x0, y)

   END SUBROUTINE fix_centred_linear_polynomial

   SUBROUTINE fix_centred_quadratic_polynomial(a2, a1, a0, x0, x, y)

      REAL, INTENT(OUT) :: a2
      REAL, INTENT(OUT) :: a1
      REAL, INTENT(OUT) :: a0
      REAL, INTENT(IN) :: x0
      REAL, INTENT(IN) :: x(3)
      REAL, INTENT(IN) :: y(3)

      CALL fix_quadratic_polynomial(a2, a1, a0, x-x0, y)

   END SUBROUTINE fix_centred_quadratic_polynomial

   SUBROUTINE fix_centred_cubic_polynomial(a3, a2, a1, a0, x0, x, y)

      REAL, INTENT(OUT) :: a3
      REAL, INTENT(OUT) :: a2
      REAL, INTENT(OUT) :: a1
      REAL, INTENT(OUT) :: a0
      REAL, INTENT(IN) :: x0
      REAL, INTENT(IN) :: x(4)
      REAL, INTENT(IN) :: y(4)

      CALL fix_cubic_polynomial(a3, a2, a1, a0, x-x0, y)

   END SUBROUTINE fix_centred_cubic_polynomial

   INTEGER FUNCTION triangle_number(n)

      ! Calculates the nth triangle number
      ! T(1) = 1, T(2) = 3, T(3) = 6, T(4) = 10, ..., T(n)=(1/2)*n*(n+1)
      INTEGER, INTENT(IN) :: n

      triangle_number = n*(n+1)/2

   END FUNCTION triangle_number

   SUBROUTINE get_Fibonaccis(F, n)

      ! Provides a sequence of the first n Fibonacci numbers
      ! F(0)=0 is not provided
      ! F(1)=1; F(2)=1; F(3)=2; F(4)=3; ...; F(n)=F(n-1)+F(n-2)
      INTEGER, INTENT(OUT) :: F(:)
      INTEGER, INTENT(IN) :: n   
      INTEGER :: i

      IF (size(F) /= n) STOP 'GET_FIBONACCIS: Error, F should be of size n'

      IF (n <= 0) THEN
         STOP 'GET_FIBONACCIS: Error, this cannot be called for n<=0'
      ELSE
         DO i = 1, n
            IF (i == 1 .OR. i == 2) THEN
               F(i) = 1
            ELSE
               F(i) = F(i-1)+F(i-2)
            END IF
         END DO
      END IF

   END SUBROUTINE get_Fibonaccis

   INTEGER FUNCTION Fibonacci(n)

      ! Returns the nth Fibonacci number
      ! F(0)=0, F(1)=1, F(2)=1, F(3)=2, F(4)=3, ..., F(n)=F(n-1)+F(n-2)
      INTEGER, INTENT(IN) :: n
      INTEGER :: F(n)

      IF (n < 0) THEN
         STOP 'FIBONACCI: Error, Fibonacci numbers undefined for n<0'
      ELSE IF (n == 0) THEN
         Fibonacci = 0
      ELSE
         CALL get_Fibonaccis(F, n)
         Fibonacci = F(n)
      END IF

   END FUNCTION Fibonacci

   SUBROUTINE get_factorials(f, n)

      ! Provides a sequence of factorial numbers up to n
      ! f(0)=1 is not provided
      ! f(1)=1, f(2)=2, f(3)=6, f(4)=24, ..., f(n)=n*f(n-1)
      ! TODO: Should this really be INT8 here?    
      INTEGER(int8), INTENT(OUT) :: f(:)
      INTEGER, INTENT(IN) :: n
      INTEGER :: i

      IF (size(f) /= n) STOP 'GET_FIBONACCIS: Error, F should be of size n'

      IF (n <= 0) THEN
         STOP 'GET_FACTORIALS: Error, this cannot be called for n<=0'
      ELSE
         DO i = 1, n
            IF (i == 1) THEN
               f(i) = 1
            ELSE
               f(i) = i*f(i-1)
            END IF
         END DO
      END IF

   END SUBROUTINE get_factorials

   INTEGER(int8) FUNCTION factorial(n)

      ! Calculates the nth factorial number
      INTEGER, INTENT(IN) :: n
      INTEGER(int8) :: f8(n)

      IF (n < 0) THEN
         STOP 'FACTORIAL: Error, factorials not defined for n<0'
      ELSE IF (n == 0) THEN
         factorial = 1
      ELSE
         CALL get_factorials(f8, n)
         factorial = f8(n)
      END IF

   END FUNCTION factorial

   INTEGER FUNCTION multiply_integers(a, b)

      ! Multiply the integers a through to b (inclusive)
      INTEGER, INTENT(IN) :: a, b
      INTEGER :: i, m

      m = 1
      DO i = a, b
         m = m*i
      END DO
      multiply_integers = m

   END FUNCTION multiply_integers

   INTEGER FUNCTION binomial_coefficient(n, k)

      ! Evalues binomial coefficient (n, k): n-choose-k
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: k

      !choose = factorial(n)/(factorial(k)*factorial(n-k)) ! Inefficient
      binomial_coefficient = multiply_integers(n-k+1, n)/factorial(k)

   END FUNCTION binomial_coefficient

   REAL FUNCTION sigmoid_tanh(x)

      ! A function that smoothly transitions from 0 to 1 around x
      REAL, INTENT(IN) :: x

      sigmoid_tanh = 0.5*(1.+tanh(x))

   END FUNCTION sigmoid_tanh

   REAL FUNCTION sigmoid_log(x, a)

      ! A function that smoothly transitions from 0 to 1 around x
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: a
      REAL :: xp, xm

      xp = x**a
      xm = x**(-a)

      sigmoid_log = 0.5*(1.+(xp-xm)/(xp+xm))

   END FUNCTION sigmoid_log

   REAL FUNCTION Legendre_polynomial(n, x)

      ! Returns the nth order Legendre polynomail: P_n(x)
      REAL, INTENT(IN) :: x
      INTEGER, INTENT(IN) :: n

      IF (n == 0) THEN
         Legendre_polynomial = 1.
      ELSE IF (n == 1) THEN
         Legendre_polynomial = x
      ELSE IF (n == 2) THEN
         Legendre_polynomial = (3.*x**2-1.)/2.
      ELSE IF (n == 3) THEN
         Legendre_polynomial = (5.*x**3-3.*x)/2.
      ELSE IF (n == 4) THEN
         Legendre_polynomial = (35.*x**4-30.*x**2+3.)/8.
      ELSE
         STOP 'LEGENDRE_POLYNOMIAL: polynomial of this order not stored'
      END IF

   END FUNCTION Legendre_polynomial

   REAL FUNCTION Lagrange_polynomial(x, xv, yv)

      ! Computes the result of the nth order Lagrange polynomial at point x, L(x)
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xv(:)
      REAL, INTENT(IN) :: yv(:)
      REAL :: l(size(xv))
      INTEGER :: i, j, n

      n = size(xv)
      IF (n /= size(yv)) STOP 'LAGRANGE_POLYNOMIAL: Error, xv and yv should have the same size'
      n = n-1

      IF (n == 0) THEN
         Lagrange_polynomial = yv(1)
      ELSE IF (n == 1) THEN
         Lagrange_polynomial = Lagrange_polynomial_1(x, xv, yv)
      ELSE IF (n == 2) THEN
         Lagrange_polynomial = Lagrange_polynomial_2(x, xv, yv)
      ELSE IF (n == 3) THEN
         Lagrange_polynomial = Lagrange_polynomial_3(x, xv, yv)
      ELSE

         ! Initialise variables, one for sum and one for multiplication
         Lagrange_polynomial = 0.
         l = 1.

         ! Loops to find the polynomials, one is a sum and one is a multiple
         DO i = 0, n
            DO j = 0, n
               IF (i .NE. j) l(i+1) = l(i+1)*(x-xv(j+1))/(xv(i+1)-xv(j+1))
            END DO
            Lagrange_polynomial = Lagrange_polynomial+l(i+1)*yv(i+1)
         END DO

      END IF

   END FUNCTION Lagrange_polynomial

   REAL FUNCTION Lagrange_polynomial_1(x, xv, yv)

      ! Dedicated function for linear interpolation polynomial
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xv(2)
      REAL, INTENT(IN) :: yv(2)
      REAL :: x01, x02
      REAL :: x12
      REAL :: f1, f2

      x01 = x-xv(1)
      x02 = x-xv(2)

      x12 = xv(1)-xv(2)

      f1 = x02*yv(1)
      f2 = x01*yv(2)

      Lagrange_polynomial_1 = (f1-f2)/x12

   END FUNCTION Lagrange_polynomial_1

   REAL FUNCTION Lagrange_polynomial_2(x, xv, yv)

      ! Dedicated function for quadratic interpolation polynomial
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xv(3)
      REAL, INTENT(IN) :: yv(3)
      REAL :: x01, x02, x03
      REAL :: x12, x13, x23
      REAL :: f1, f2, f3

      x01 = x-xv(1)
      x02 = x-xv(2) 
      x03 = x-xv(3)

      x12 = xv(1)-xv(2)
      x13 = xv(1)-xv(3)
      x23 = xv(2)-xv(3)

      f1 = yv(1)*x02*x03/(x12*x13)
      f2 = x01*yv(2)*x03/(x12*x23)
      f3 = x01*x02*yv(3)/(x13*x23)

      Lagrange_polynomial_2 = f1-f2+f3

   END FUNCTION Lagrange_polynomial_2

   REAL FUNCTION Lagrange_polynomial_3(x, xv, yv)

      ! Dedicated function for cubic interpolation polynomial
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xv(4)
      REAL, INTENT(IN) :: yv(4)
      REAL :: x01, x02, x03, x04
      REAL :: x12, x13, x14, x23, x24, x34
      REAL :: f1, f2, f3, f4

      x01 = x-xv(1)
      x02 = x-xv(2) 
      x03 = x-xv(3)
      x04 = x-xv(4)

      x12 = xv(1)-xv(2)
      x13 = xv(1)-xv(3)
      x14 = xv(1)-xv(4)
      x23 = xv(2)-xv(3)
      x24 = xv(2)-xv(4)
      x34 = xv(3)-xv(4)

      f1 = yv(1)*x02*x03*x04/(x12*x13*x14)
      f2 = x01*yv(2)*x03*x04/(x12*x23*x24)
      f3 = x01*x02*yv(3)*x04/(x13*x23*x34)
      f4 = x01*x02*x03*yv(4)/(x14*x24*x34)

      Lagrange_polynomial_3 = f1-f2+f3-f4

   END FUNCTION Lagrange_polynomial_3

   REAL FUNCTION sinc(x)

      ! sinc function: sin(x)/x
      ! TODO: Is the Taylor expansion here unnecessary?
      REAL, INTENT(IN) :: x
      REAL, PARAMETER :: dx = dx_sinc ! small |x| below which to use Taylor expansion

      IF (abs(x) < dx) THEN
         sinc = 1.-(x**2)/6.+(x**4)/120.
      ELSE
         sinc = sin(x)/x
      END IF

   END FUNCTION sinc

   REAL FUNCTION wk_tophat(x)

      ! The normlaised Fourier Transform of a spherical top-hat
      REAL, INTENT(IN) :: x
      REAL, PARAMETER :: dx = dx_tophat ! Taylor expansion for |x|<dx

      ! Taylor expansion used for low x to avoid cancelation problems
      IF (abs(x) < dx) THEN
         wk_tophat = 1.-x**2/10.
      ELSE
         wk_tophat = (3./x**3)*(sin(x)-x*cos(x))
      END IF

   END FUNCTION wk_tophat

   REAL FUNCTION wk_tophat_deriv(x)

      ! The derivative of a normlaised Fourier Transform of a spherical top-hat
      REAL, INTENT(IN) :: x
      REAL, PARAMETER :: dx = dx_tophat ! Taylor expansion for |x|<dx

      ! Taylor expansion used for low x to avoid cancelation problems
      IF (abs(x) < dx) THEN
         wk_tophat_deriv = -x/5.+x**3/70.
      ELSE
         wk_tophat_deriv = (3./x**4)*((x**2-3.)*sin(x)+3.*x*cos(x))
      END IF

   END FUNCTION wk_tophat_deriv

   REAL FUNCTION wk_tophat_dderiv(x)

      ! The second derivative of a normlaised Fourier Transform of a spherical top-hat
      REAL, INTENT(IN) :: x
      REAL, PARAMETER :: dx = dx_tophat ! Taylor expansion for |x|<dx

      ! Taylor expansion used for low x to avoid cancelation problems
      IF (abs(x) < dx) THEN
         wk_tophat_dderiv = -0.2+3.*x**2/70.
      ELSE
         wk_tophat_dderiv = (3./x**5)*((12.-5.*x**2)*sin(x)+x*(x**2-12.)*cos(x))
      END IF

   END FUNCTION wk_tophat_dderiv

   REAL FUNCTION apodise(x, x1, x2, n)

      ! Apodises a function between x1 and x2
      ! Goes to one smoothly at x1
      ! Goes to zero linearly at x2, so the gradient change is discontinous
      ! n govenrns the severity of the transition
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x1
      REAL, INTENT(IN) :: x2
      REAL, INTENT(IN) :: n

      IF (n <= 0.) STOP 'APODISE: Error, n must be greater than zero'

      IF (x < x1) THEN
         apodise = 1.
      ELSE IF (x > x2) THEN
         apodise = 0.
      ELSE
         apodise = cos((pi/2.)*(x-x1)/(x2-x1))
         apodise = apodise**n
      END IF

   END FUNCTION apodise

   REAL FUNCTION smooth_apodise(x, x1, x2, n)

      ! Apodises a function between x1 and x2
      ! Goes to one smoothly at x1 and zero smoothly at x2
      ! n govenrns the severity of the transition
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x1
      REAL, INTENT(IN) :: x2
      REAL, INTENT(IN) :: n

      IF (n <= 0.) STOP 'APODISE: Error, n must be greater than zero'

      IF (x < x1) THEN
         smooth_apodise = 1.
      ELSE IF (x > x2) THEN
         smooth_apodise = 0.
      ELSE
         smooth_apodise = 0.5*(1.+cos(pi*(x-x1)/(x2-x1)))
         smooth_apodise = smooth_apodise**n
      END IF

   END FUNCTION smooth_apodise

   REAL FUNCTION blob(x, x1, x2, n)

      ! Makes a blob between x1 and x2, with zero elsewhere
      ! Blob goes to zero linearly at x1 and x2, so the gradient change is discontinous
      ! n governs the severity (blobiness) of the blob
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x1
      REAL, INTENT(IN) :: x2
      REAL, INTENT(IN) :: n

      IF (n <= 0.) STOP 'APODISE: Error, n must be greater than zero'

      IF (x < x1) THEN
         blob = 0.
      ELSE IF (x > x2) THEN
         blob = 0.
      ELSE
         blob = sin(pi*(x-x1)/(x2-x1))
         blob = blob**n
      END IF

   END FUNCTION blob

   REAL FUNCTION smooth_blob(x, x1, x2, n)

      ! Makes a blob between x1 and x2, with zero elsewhere
      ! Blob goes to zero smoothly at x1 and x2
      ! n governs the severity (blobiness) of the blob
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x1
      REAL, INTENT(IN) :: x2
      REAL, INTENT(IN) :: n

      IF (n <= 0.) STOP 'APODISE: Error, n must be greater than zero'

      IF (x < x1) THEN
         smooth_blob = 0.
      ELSE IF (x > x2) THEN
         smooth_blob = 0.
      ELSE
         smooth_blob = (1.+cos(twopi*(x-x1)/(x2-x1)))/2.
         smooth_blob = (1.-smooth_blob)**n
      END IF

   END FUNCTION smooth_blob

   REAL FUNCTION Si(x)

      ! Returns the 'sine integral' function: Si(x)=int_0^x sin(t)/t dt
      REAL, INTENT(IN) :: x
      REAL :: x2, y, f, g, si8
      REAL, PARAMETER :: x0 = x0_SiCi ! Transition between two different approximations

      ! Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
      IF (abs(x) <= x0) THEN

         x2 = x*x

         si8 = x*(1.d0+x2*(-4.54393409816329991d-2+x2*(1.15457225751016682d-3 &
               +x2*(-1.41018536821330254d-5+x2*(9.43280809438713025d-8+x2*(-3.53201978997168357d-10 &
               +x2*(7.08240282274875911d-13+x2*(-6.05338212010422477d-16))))))))/ &
               (1.+x2*(1.01162145739225565d-2+x2*(4.99175116169755106d-5+ &
               x2*(1.55654986308745614d-7+x2*(3.28067571055789734d-10+x2*(4.5049097575386581d-13 &
               +x2*(3.21107051193712168d-16)))))))

         Si = real(si8)

      ELSE IF (abs(x) > x0) THEN

         y = 1.d0/(x*x)

         f = (1.d0+y*(7.44437068161936700618d2+y*(1.96396372895146869801d5+ &
               y*(2.37750310125431834034d7+y*(1.43073403821274636888d9+y*(4.33736238870432522765d10 &
               +y*(6.40533830574022022911d11+y*(4.20968180571076940208d12+ &
               y*(1.00795182980368574617d13+y*(4.94816688199951963482d12+ &
               y*(-4.94701168645415959931d11)))))))))))/(x*(1.+y*(7.46437068161927678031d2+ &
               y*(1.97865247031583951450d5+y*(2.41535670165126845144d7+ &
               y*(1.47478952192985464958d9+y*(4.58595115847765779830d10+ &
               y*(7.08501308149515401563d11+y*(5.06084464593475076774d12+ &
               y*(1.43468549171581016479d13+y*(1.11535493509914254097d13)))))))))))

         g = y*(1.d0+y*(8.1359520115168615d2+y*(2.35239181626478200d5+ &
               y*(3.12557570795778731d7+y*(2.06297595146763354d9+y*(6.83052205423625007d10+ &
               y*(1.09049528450362786d12+y*(7.57664583257834349d12+y*(1.81004487464664575d13+ &
               y*(6.43291613143049485d12+y*(-1.36517137670871689d12)))))))))))/ &
               (1.+y*(8.19595201151451564d2+y*(2.40036752835578777d5+y*(3.26026661647090822d7 &
               +y*(2.23355543278099360d9+y*(7.87465017341829930d10+y*(1.39866710696414565d12 &
               +y*(1.17164723371736605d13+y*(4.01839087307656620d13+y*(3.99653257887490811d13))))))))))

         Si = real(pi/2.d0-f*cos(x)-g*sin(x))

      ELSE

         STOP 'SI: Something went very wrong'

      END IF

   END FUNCTION Si

   REAL FUNCTION Ci(x)

      ! Returns the 'cosine integral' function Ci(x): -int_x^inf cos(t)/t dt
      REAL, INTENT(IN) :: x
      REAL :: x2, y, f, g, ci8
      REAL, PARAMETER :: x0 = x0_SiCi ! Transition between two different approximations

      ! Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
      IF (abs(x) <= x0) THEN

         x2 = x*x

         ci8 = em+log(x)+x2*(-0.25d0+x2*(7.51851524438898291d-3+x2*(-1.27528342240267686d-4 &
               +x2*(1.05297363846239184d-6+x2*(-4.68889508144848019d-9+x2*(1.06480802891189243d-11 &
               +x2*(-9.93728488857585407d-15)))))))/(1.+x2*(1.1592605689110735d-2+ &
               x2*(6.72126800814254432d-5+x2*(2.55533277086129636d-7+x2*(6.97071295760958946d-10+ &
               x2*(1.38536352772778619d-12+x2*(1.89106054713059759d-15+x2*(1.39759616731376855d-18))))))))

         Ci = real(ci8)

      ELSE IF (abs(x) > x0) THEN

         y = 1./(x*x)

         f = (1.d0+y*(7.44437068161936700618d2+y*(1.96396372895146869801d5+ &
               y*(2.37750310125431834034d7+y*(1.43073403821274636888d9+y*(4.33736238870432522765d10 &
               +y*(6.40533830574022022911d11+y*(4.20968180571076940208d12+y*(1.00795182980368574617d13 &
               +y*(4.94816688199951963482d12+y*(-4.94701168645415959931d11)))))))))))/ &
               (x*(1.+y*(7.46437068161927678031d2+y*(1.97865247031583951450d5+ &
               y*(2.41535670165126845144d7+y*(1.47478952192985464958d9+ &
               y*(4.58595115847765779830d10+y*(7.08501308149515401563d11+y*(5.06084464593475076774d12 &
               +y*(1.43468549171581016479d13+y*(1.11535493509914254097d13)))))))))))

         g = y*(1.d0+y*(8.1359520115168615d2+y*(2.35239181626478200d5+y*(3.12557570795778731d7 &
               +y*(2.06297595146763354d9+y*(6.83052205423625007d10+ &
               y*(1.09049528450362786d12+y*(7.57664583257834349d12+ &
               y*(1.81004487464664575d13+y*(6.43291613143049485d12+y*(-1.36517137670871689d12))))))))))) &
               /(1.+y*(8.19595201151451564d2+y*(2.40036752835578777d5+ &
               y*(3.26026661647090822d7+y*(2.23355543278099360d9+y*(7.87465017341829930d10 &
               +y*(1.39866710696414565d12+y*(1.17164723371736605d13+y*(4.01839087307656620d13+y*(3.99653257887490811d13))))))))))

         Ci = real(f*sin(x)-g*cos(x))

      ELSE

         STOP 'CI: Something went very wrong'

      END IF

   END FUNCTION Ci

   REAL FUNCTION Bessel(n, x)

      ! Returns the Bessel function of order 'n'
      ! Wraps the Fortran intrinsic functions
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: x  
      REAL, PARAMETER :: xbig = xbig_Bessel ! Set to zero for large values

      IF (x > xbig) THEN

         ! To stop it going mental for very large values of x
         Bessel = 0.

      ELSE

         IF (n < 0) THEN
            WRITE(*,*) 'BESSEL: Order:', n
            STOP 'BESSEL: cannot call for negative n'
         END IF

         IF (n == 0) THEN
            Bessel = Bessel_J0(x)
         ELSE IF (n == 1) THEN
            Bessel = Bessel_J1(x)
         ELSE
            Bessel = Bessel_JN(n, x)
         END IF

      END IF

   END FUNCTION Bessel

   REAL FUNCTION spherical_Bessel(x, n)

      ! Spherical Bessel function of order n: j_n(x)
      ! Unlike the standard Bessel functions, these are expressible in terms of sin(x) and cos(x)
      ! Copied from https://en.wikipedia.org/wiki/Bessel_function
      ! Limits from https://www.wolframalpha.com/input/?i=series+expand+j_0%28x%29
      REAL, INTENT(IN) :: x
      INTEGER, INTENT(IN) :: n
      REAL, PARAMETER :: dx = dx_Bessel

      IF (abs(x) < dx) THEN
         IF (n == 0) THEN
            spherical_Bessel = 1.-x**2/6.
         ELSE IF (n == 1) THEN
            spherical_Bessel = x/3.-x**3/30.
         ELSE IF (n == 2) THEN
            spherical_Bessel = x**2/15.-x**4/210.
         ELSE IF (n == 3) THEN
            spherical_Bessel = x**3/105.-x**5/1890.
         ELSE
            WRITE(*, *) 'SPHERICAL_BESSEL: n:', n
            STOP 'SPHERICAL_BESSEL: Error, value of n not currently supported'
         END IF
      ELSE
         IF (n == 0) THEN
            spherical_Bessel = sin(x)/x
         ELSE IF (n == 1) THEN
            spherical_Bessel = sin(x)/x**2-cos(x)/x
         ELSE IF (n == 2) THEN
            spherical_Bessel = (3./x**2-1.)*sin(x)/x-3.*cos(x)/x**2
         ELSE IF (n == 3) THEN
            spherical_Bessel = (15./x**3-6./x)*sin(x)/x-(15./x**2-1.)*cos(x)/x
         ELSE
            WRITE(*, *) 'SPHERICAL_BESSEL: n:', n
            STOP 'SPHERICAL_BESSEL: Error, value of n not currently supported'
         END IF
      END IF

   END FUNCTION spherical_Bessel

   REAL FUNCTION Gaussian_distribution(x, mu, sigma)

      ! Returns the integral-normalised Gaussian
      REAL, INTENT(IN) :: x     ! [-inf:inf]
      REAL, INTENT(IN) :: mu    ! Mean value
      REAL, INTENT(IN) :: sigma ! Root-variance
      REAL :: f1, f2

      f1 = exp(-((x-mu)**2)/(2.*sigma**2))
      f2 = sigma*sqrt(twopi)

      Gaussian_distribution = f1/f2

   END FUNCTION Gaussian_distribution

   REAL FUNCTION lognormal_distribution(x, mean, sd)

      ! Returns integral-normalised lognormal distribution
      REAL, INTENT(IN) :: x    ! x [0,inf]
      REAL, INTENT(IN) :: mean ! Mean value of x
      REAL, INTENT(IN) :: sd   ! Standard deviation of x
      REAL :: mu, sigma

      IF (mean <= 0.) STOP 'LOGNORMAL_DISTRIBUTION: Error, mean cannot be less than or equal to zero'
      IF (sd < 0.)    STOP 'LOGNORMAL_DISTRIBUTION: Error, standard deviation cannot be less than zero'

      mu = log(mean/sqrt(1.+(sd/mean)**2))
      sigma = sqrt(log(1.+(sd/mean)**2))
      lognormal_distribution = Gaussian_distribution(log(x), mu, sigma)/x

   END FUNCTION lognormal_distribution

   REAL FUNCTION uniform_distribution(x, x1, x2)

      ! Returns integral-normalised one-dimensional top-hat function between x1 and x2
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x1 ! Lower limit
      REAL, INTENT(IN) :: x2 ! Upper limit

      IF (x < x1 .OR. x > x2) THEN
         uniform_distribution = 0.
      ELSE
         uniform_distribution = 1./(x2-x1)
      END IF

   END FUNCTION uniform_distribution

   REAL FUNCTION Rayleigh_distribution(x, sigma)

      ! Returns integral-normalised Rayleigh distribution
      REAL, INTENT(IN) :: x     ! [0:inf]
      REAL, INTENT(IN) :: sigma ! Sigma parameter (*not* root-variance for this distribution)

      IF (sigma <= 0.) STOP 'RAYLEIGH_DISTRIBUTION: Error, sigma cannot be less than or equal to zero'

      IF (x < 0.) THEN
         STOP 'RAYLEIGH_DISTRIBUTION: Error, x cannot be less than zero'
      ELSE
         Rayleigh_distribution = x*exp(-(x**2)/(2.*(sigma**2)))/(sigma**2)
      END IF

   END FUNCTION Rayleigh_distribution

   REAL FUNCTION Poisson_distribution(n, nbar)

      ! Normalised discrete Poisson probability distribution
      INTEGER, INTENT(IN) :: n ! Number of events to evaluate P_n at, n>=0
      REAL, INTENT(IN) :: nbar ! Mean number of events >0

      Poisson_distribution = exp(-nbar)*(nbar**n)/int(factorial(n))

   END FUNCTION Poisson_distribution

   REAL FUNCTION Bernoulli_distribution(k, p)

      ! One trial with probability of success: p
      INTEGER, INTENT(IN) :: k ! 0 or 1 only
      REAL, INTENT(IN) :: p    ! Must be between 0 and 1.

      Bernoulli_distribution = (p**k)*(1.-p)**(1-k)
      !IF (k == 0) THEN
      !   Bernoulli_distribution = 1.-p
      !ELSE IF (k == 1) THEN
      !   Bernoulli_distribution = p
      !ELSE
      !   STOP 'BERNOULLI_DISTRIBUTION: Error, k can only be zero or one'
      !END IF

   END FUNCTION Bernoulli_distribution

   REAL FUNCTION geometric_distribution(k, p)

      ! Probability for number of failures until the first success
      ! Each of k trials is independent and has chance of success k
      INTEGER, INTENT(IN) :: k ! Must be 0 or greater
      REAL, INTENT(IN) :: p    ! Must be between 0 and 1.

      geometric_distribution = p*(1.-p)**k

   END FUNCTION geometric_distribution

   REAL FUNCTION shifted_geometric_distribution(k, p)

      ! Probability for number of trials until the first success
      ! Each of k trials is independent and has chance of success k
      INTEGER, INTENT(IN) :: k ! Must be 1 or greater
      REAL, INTENT(IN) :: p    ! Must be between 0 and 1.

      !shifted_geometric_distribution = p*(1.-p)**(k-1)
      shifted_geometric_distribution = geometric_distribution(k-1, p)

   END FUNCTION shifted_geometric_distribution

   REAL FUNCTION Binomial_distribution(k, n, p)

      INTEGER, INTENT(IN) :: k
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: p

      Binomial_distribution = binomial_coefficient(n, k)*(p**k)*(1.-p)**(n-k)

   END FUNCTION Binomial_distribution

   REAL FUNCTION exponential_distribution(x, mean)

      ! Returns integral-normalised exponential distribution
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: mean

      exponential_distribution = exp(-x/mean)/mean

   END FUNCTION exponential_distribution

   REAL FUNCTION Lorentzian_distribution(x)

      ! Returns integral-normalised Lorentzian distribution
      REAL, INTENT(IN) :: x

      Lorentzian_distribution = 2./(pi*(1.+x**2))

   END FUNCTION Lorentzian_distribution

   REAL FUNCTION polynomial_distribution(x, n)

      ! Returns integral-normalised polynomial distribution
      REAL, INTENT(IN) :: x ! x[0->1]
      REAL, INTENT(IN) :: n ! Polynomail order [n>-1]

      IF (n < -1) STOP 'POLYNOMIAL_DISTRIBUTION: Error, index is less than -1'

      IF (x < 0. .OR. x > 1.) THEN
         STOP 'POLYNOMIAL_DISTRIBUTION: Error, x is outside range 0 to 1'
      ELSE
         polynomial_distribution = (n+1.)*x**n
      END IF

   END FUNCTION polynomial_distribution

   REAL FUNCTION Rosenbrock(x, y)

      ! https://en.wikipedia.org/wiki/Rosenbrock_function
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y

      Rosenbrock = (1.-x)**2+100.*(y-x**2)**2

   END FUNCTION Rosenbrock

   REAL FUNCTION Himmelblau(x, y)
   
      ! https://en.wikipedia.org/wiki/Himmelblau%27s_function
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y

      Himmelblau = (x**2+y-11.)**2+(x+y**2-7.)**2

   END FUNCTION Himmelblau

END MODULE special_functions
