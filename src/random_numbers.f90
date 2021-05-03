MODULE random_numbers

   USE basic_operations

   IMPLICIT NONE

   PRIVATE

   ! Setting routine
   PUBLIC :: random_generator_seed

   ! Integer distributions
   PUBLIC :: random_integer
   PUBLIC :: random_sign
   PUBLIC :: dice
   PUBLIC :: random_Bernoulli
   PUBLIC :: random_twopoint
   PUBLIC :: random_binomial
   PUBLIC :: random_Poisson

   ! Real number distributions
   PUBLIC :: random_unit
   PUBLIC :: random_uniform
   PUBLIC :: random_Rayleigh
   PUBLIC :: random_Lorentzian
   PUBLIC :: random_Gaussian
   PUBLIC :: random_Gaussian_pair
   PUBLIC :: random_lognormal
   PUBLIC :: random_exponential
   PUBLIC :: random_polynomial
   PUBLIC :: random_spherical_theta

   ! Complex number disitrubutions
   PUBLIC :: random_complex_unit

   ! Draw from any real-number distribution
   PUBLIC :: accept_reject

   INTERFACE random_twopoint
      MODULE PROCEDURE random_twopoint_real
      MODULE PROCEDURE random_twopoint_integer
   END INTERFACE random_twopoint

CONTAINS

   ! SUBROUTINE RNG_set(seed, verbose)

   !    ! Seeds the RNG
   !    ! TODO: Retire rand, replace with random_number
   !    IMPLICIT NONE
   !    INTEGER, INTENT(IN) :: seed
   !    LOGICAL, OPTIONAL, INTENT(IN) :: verbose
   !    INTEGER :: int, timearray(3)
   !    REAL(kind=4) :: rand ! Necessary to define for ifort, also the *4 is necessary

   !    IF (present_and_correct(verbose)) THEN
   !       WRITE (*, *) 'RNG_SET: Initialising random number generator'
   !       WRITE (*, *) 'RNG_SET: Seed:', seed
   !    END IF

   !    IF (seed == 0) THEN

   !       ! This fills the time array using the system clock!
   !       ! If called within the same second the numbers will be identical!
   !       CALL itime(timeArray)

   !       ! This then initialises the generator!
   !       int = floor(rand(timeArray(1)+timeArray(2)+timeArray(3)))

   !    ELSE

   !       ! In this case you can keep track of the seed
   !       int = floor(rand(seed))

   !    END IF

   !    IF (present_and_correct(verbose)) THEN
   !       WRITE (*, *) 'RNG_SET: Done'
   !       WRITE (*, *)
   !    END IF

   ! END SUBROUTINE RNG_set

   SUBROUTINE random_generator_seed(seed, verbose)

      ! TODO: Understand seeding properly
      INTEGER, INTENT(IN) :: seed
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: n
      INTEGER, ALLOCATABLE :: seedy(:)

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'RANDOM_GENERATOR_SEED: Initialising random number generator'
         IF(seed .NE. 0) WRITE (*, *) 'RANDOM_GENERATOR_SEED: Seed:', seed
      END IF

      IF (seed == 0) THEN
         CALL random_seed()
      ELSE
         CALL random_seed(size=n)
         ALLOCATE(seedy(n))
         seedy = seed
         CALL random_seed(put=seedy)
      END IF

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'RANDOM_GENERATOR_SEED: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE random_generator_seed

   INTEGER FUNCTION random_integer(i1, i2)

      ! Picks an integer with uniform random probability between i1 and i2 spaced with 1
      INTEGER, INTENT(IN) :: i1 ! Lower bound
      INTEGER, INTENT(IN) :: i2 ! Upper bound
      INTEGER, PARAMETER :: seed = 0

      random_integer = i1-1+ceiling(random_unit()*real(1+i2-i1))
      IF (random_integer == i1-1) random_integer = i1

   END FUNCTION random_integer

   INTEGER FUNCTION dice(ndice, dmin_opt, dmax_opt)

      ! Get a total for rolling ndice
      INTEGER, INTENT(IN) :: ndice               ! Number of dice to roll   
      INTEGER, OPTIONAL, INTENT(IN) :: dmin_opt  ! Minimum value on di
      INTEGER, OPTIONAL, INTENT(IN) :: dmax_opt  ! Maximum value on di (assumes all integers on di between dmin and dmax)    
      INTEGER :: i
      INTEGER :: dmin, dmax
      INTEGER, PARAMETER :: dmin_def = 1
      INTEGER, PARAMETER :: dmax_def = 6

      ! Check number of dice
      IF (ndice < 0) STOP 'DICE: Error, number of rolls must be positive'

      ! Default or optional values for the minimum and maximum on the di
      ! Note that all contiguous integers are taken between dmin and dmax
      dmin = default_or_optional(dmin_def, dmin_opt)
      dmax = default_or_optional(dmax_def, dmax_opt)

      ! Roll the dice and sum the score
      dice = 0
      DO i = 1, ndice
         dice = dice+random_integer(dmin, dmax)
      END DO

   END FUNCTION dice

   REAL FUNCTION random_twopoint_real(a, b, p)

      ! Pick 'a' with probability 'p' and 'b' with probability '1-p'
      REAL, INTENT(IN) :: a, b
      REAL, INTENT(IN) :: p
      REAL :: r

      r = random_unit()
      IF (r <= p) THEN
         random_twopoint_real = a
      ELSE
         random_twopoint_real = b
      END IF

   END FUNCTION random_twopoint_real

   INTEGER FUNCTION random_twopoint_integer(a, b, p)

      ! Pick 'a' with probability '1-p' and 'b' with probability 'p'
      INTEGER, INTENT(IN) :: a, b
      REAL, INTENT(IN) :: p
      REAL :: r

      r = random_unit()
      IF (r <= p) THEN
         random_twopoint_integer = a
      ELSE
         random_twopoint_integer = b
      END IF

   END FUNCTION random_twopoint_integer

   INTEGER FUNCTION random_Bernoulli(p)

      ! Unity with probability p, otherwise zero
      ! PDF: P_k = p^k (1-p)^(1-k) with k = 0, 1
      REAL, INTENT(IN) :: p ! Should be between 0 and 1

      random_Bernoulli=random_twopoint(1, 0, p)

   END FUNCTION random_Bernoulli

   INTEGER FUNCTION random_sign()

      ! Returns either +1 or -1 with equal probability
      random_sign = random_twopoint(-1, 1, 0.5)

   END FUNCTION random_sign

   INTEGER FUNCTION random_binomial(n, p)

      ! Generates a random number from a Binomial distribution with n trials, each with probability of success p
      ! Note that random_Binomial will always be between [0, n]; closer to n with higher p
      INTEGER, INTENT(IN) :: n ! Number of trials
      REAL, INTENT(IN) :: p    ! Probability of success of each trial
      INTEGER :: sum, i

      sum = 0
      DO i = 1, n
         sum = sum+random_Bernoulli(p)
      END DO
      random_binomial = sum

   END FUNCTION random_binomial

   INTEGER FUNCTION random_Poisson(mean)

      ! Generate a random number from a Poisson distribution 
      REAL, INTENT(IN) :: mean
      INTEGER :: k
      REAL :: L, p

      L = exp(-mean)

      k = 0
      p = 1.
      DO
         p = p*random_unit()
         IF (p < L) THEN
            EXIT
         ELSE
            k = k+1
         END IF
      END DO

      random_Poisson = k

   END FUNCTION random_Poisson

   REAL FUNCTION random_unit()

      ! Produces a uniform-random number between 0 and 1

      CALL random_number(random_unit)

   END FUNCTION random_unit

   REAL FUNCTION random_uniform(x1, x2)

      ! Produces a uniform random number between x1 and x2
      REAL, INTENT(IN) :: x1 ! Lower bound
      REAL, INTENT(IN) :: x2 ! Upper bound

      random_uniform = x1+(x2-x1)*random_unit()

   END FUNCTION random_uniform

   REAL FUNCTION random_Rayleigh(sigma)

      ! Produces a Rayleigh-distributed random number
      USE constants
      REAL, INTENT(IN) :: sigma        ! Sigma parameter (*not* standard deviation of the distribution)
      REAL, PARAMETER :: small = 1e-10 ! To avoid ever getting a log(0) call

      ! Problems if small=0. because log(0.) gets called sometimes
      random_Rayleigh = sigma*sqrt(-2.*log(random_uniform(small, one)))

   END FUNCTION random_Rayleigh

   REAL FUNCTION random_Lorentzian()

      ! Produces a Lorentzian-distributed random number
      USE constants

      random_Lorentzian = tan(random_uniform(zero, pi/2.))

   END FUNCTION random_Lorentzian

   FUNCTION random_Gaussian_pair(mean, sigma)

      ! Gets a pair of independent Gaussian random numbers
      ! This uses the Box-Muller method (https://en.wikipedia.org/wiki/Box-Muller_transform)
      USE constants
      REAL :: random_Gaussian_pair(2)
      REAL, INTENT(IN) :: mean  ! Mean of the distribution
      REAL, INTENT(IN) :: sigma ! Root-variance of the distribution
      REAL :: r, theta

      r = random_Rayleigh(sigma)
      theta = random_uniform(zero, twopi)

      ! Both of these numbers are independently Gaussian
      random_Gaussian_pair(1) = r*sin(theta)+mean
      random_Gaussian_pair(2) = r*cos(theta)+mean

   END FUNCTION random_Gaussian_pair

   REAL FUNCTION random_Gaussian(mean, sigma)

      ! Gets a single Gaussian random number
      ! This is wasteful as there is a second, independent Gaussian random number that is thrown away
      REAL, INTENT(IN) :: mean  ! Mean of the distribution
      REAL, INTENT(IN) :: sigma ! Root-variance of the distribution
      REAL :: G(2)

      G = random_Gaussian_pair(mean, sigma)
      random_Gaussian = G(1)

   END FUNCTION random_Gaussian

   REAL FUNCTION random_lognormal(mean, sd)

      ! Gets a single log-normal random number
      ! TODO: Could do random_lognormal_pair
      REAL, INTENT(IN) :: mean ! Mean of the distribution
      REAL, INTENT(IN) :: sd   ! Standard deviation of the distribution
      REAL :: mu, sigma

      mu = log(mean/sqrt(1.+(sd/mean)**2))
      sigma = sqrt(log(1.+(sd/mean)**2))
      random_lognormal = exp(random_Gaussian(mu, sigma))

   END FUNCTION random_lognormal

   REAL FUNCTION random_exponential(mean)

      ! Produces a exponentially-distributed random number
      REAL, INTENT(IN) :: mean         ! Mean of the distribution
      REAL, PARAMETER :: small = 1e-10 ! Needed because problems here if log(0) is ever called

      random_exponential = -mean*log(random_uniform(small, 1.))

   END FUNCTION random_exponential

   REAL FUNCTION random_polynomial(n)

      ! Generate a polynomailly distributed number [x:0->1]
      REAL, INTENT(IN) :: n ! Order for the polynomial [-1:inf]

      !IF (n <= -1) STOP 'RANDOM_POLYNOMIAL: Error, n is less than or equal to -1'
      random_polynomial = random_unit()**(1./(n+1))

   END FUNCTION random_polynomial

   REAL FUNCTION random_spherical_theta()

      ! A random spherical-polar angle such that the solid-angle is uniformally populated
      random_spherical_theta = acos(random_uniform(-1., 1.))

   END FUNCTION random_spherical_theta

   COMPLEX FUNCTION random_complex_unit()

      ! Get a complex phase with theta between 0 and 2pi
      ! Generates a unit amplitude complex number with random phase
      ! TODO: Is 0 actually counted twice because 0 and 2pi are identical?
      USE constants
      REAL :: theta

      theta = random_uniform(0., twopi)
      random_complex_unit = cmplx(cos(theta), sin(theta))

   END FUNCTION random_complex_unit

   REAL FUNCTION accept_reject(func, x1, x2, fmax)

      ! Simple one-dimensional accept-reject algorithm for drawing random numbers from func(x) between x1 and x2
      ! TODO: Increase to n-dimensions
      ! TODO: Include more complicated bounding structure (at the moment it is just a box)
      REAL, EXTERNAL :: func   ! Function to sample from
      REAL, INTENT(IN) :: x1   ! Lower bound for function
      REAL, INTENT(IN) :: x2   ! Upper bound for function
      REAL, INTENT(IN) :: fmax ! Maximum value of the function in the interval x1 to x2
      REAL :: x, y, f

      INTERFACE
         FUNCTION func(xin)
            REAL, INTENT(IN) :: xin
         END FUNCTION func
      END INTERFACE

      ! Try until a value is accepted
      DO

         ! Draw uniform-random x value and function height
         x = random_uniform(x1, x2)
         y = random_uniform(0., fmax)

         ! Evaulate the function at the random x value
         f = func(x)

         ! Decide whether or not to accept
         IF (f > fmax) THEN
            STOP 'ACCEPT_REJECT: Error, your function is not bounded by fmax'
         ELSE IF (y <= f) THEN
            accept_reject = x
            EXIT
         END IF

      END DO

   END FUNCTION accept_reject

END MODULE random_numbers
