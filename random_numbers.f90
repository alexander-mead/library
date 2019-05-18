MODULE random_numbers

  !TODO: Think about using intrinsic 'random_number' instead of 'rand()'
  USE logical_operations
  
  IMPLICIT NONE

  PRIVATE

  ! Setting routine
  PUBLIC :: RNG_set

  ! Integer distributions
  PUBLIC :: random_integer
  PUBLIC :: random_sign
  PUBLIC :: dice

  ! Real number distributions
  PUBLIC :: random_uniform
  PUBLIC :: random_Rayleigh
  PUBLIC :: random_Lorentzian
  PUBLIC :: random_Gaussian
  PUBLIC :: random_lognormal
  PUBLIC :: random_exponential
  PUBLIC :: random_polynomial
  PUBLIC :: random_theta

  ! Draw from any real-number distribution
  PUBLIC :: accept_reject
  
CONTAINS

  SUBROUTINE RNG_set(seed,verbose)

    ! Seeds the RNG
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: seed
    LOGICAL, OPTIONAL, INTENT(IN) :: verbose
    INTEGER :: int, timearray(3)
    REAL*4 :: rand ! Necessary to define for ifort, also the *4 is necessary

    IF(present_and_correct(verbose)) THEN
       WRITE(*,*) 'RNG_SET: Initialising random number generator'
       WRITE(*,*) 'RNG_SET: Seed:', seed
    END IF

    IF(seed==0) THEN
       
       ! This fills the time array using the system clock!
       ! If called within the same second the numbers will be identical!
       CALL itime(timeArray)
       
       ! This then initialises the generator!
       int=floor(rand(timeArray(1)+timeArray(2)+timeArray(3)))
       
    ELSE
       
       ! In this case you can keep track of the seed
       int=floor(rand(seed))
       
    END IF

    IF(present_and_correct(verbose)) THEN
       WRITE(*,*) 'RNG_SET: Done'
       WRITE(*,*)
    END IF
       
  END SUBROUTINE RNG_set

  INTEGER FUNCTION random_integer(i1,i2)

    ! Picks an integer with uniform random probability between i1 and i2 spaced with 1
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1 ! Lower bound
    INTEGER, INTENT(IN) :: i2 ! Upper bound
    REAL*4 :: rand ! Necessary to define for ifort

    random_integer=i1-1+ceiling(rand(0)*real(1+i2-i1))
    IF(random_integer==i1-1) random_integer=i1

  END FUNCTION random_integer

  INTEGER FUNCTION dice(dmin,dmax,ndice)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dmin
    INTEGER, INTENT(IN) :: dmax
    INTEGER, INTENT(IN) :: ndice
    INTEGER :: i

    IF(ndice<0) STOP 'DICE: Error, number of rolls must be positive or zero'
    
    ! Roll the dice and sum the score
    dice=0
    DO i=1,ndice
       dice=dice+random_integer(dmin,dmax)
    END DO
    
  END FUNCTION dice

  INTEGER FUNCTION random_sign()

    ! Returns either +1 or -1 with equal probability
    IMPLICIT NONE

    random_sign=random_integer(0,1)

    IF(random_sign==0) random_sign=-1
    
  END FUNCTION random_sign

  REAL FUNCTION random_uniform(x1,x2)

    ! Produces a uniform random number between x1 and x2
    IMPLICIT NONE
    REAL, INTENT(IN) :: x1 ! Lower bound
    REAL, INTENT(IN) :: x2 ! Upper bound
    REAL*4 :: rand !I think this needs to be defined for ifort

    ! rand is some inbuilt function
    random_uniform=x1+(x2-x1)*(rand(0))

  END FUNCTION random_uniform

  FUNCTION random_Rayleigh(sigma)

    ! Produces a Rayleigh-distributed random number
    IMPLICIT NONE
    REAL :: random_Rayleigh
    REAL, INTENT(IN) :: sigma
    REAL, PARAMETER :: small=1e-10

    ! Problems if small=0. because log(0.) gets called sometimes
    random_Rayleigh=sigma*sqrt(-2.*log(random_uniform(small,1.)))

  END FUNCTION random_Rayleigh

  FUNCTION random_Lorentzian()

    ! Produces a Lorentzian-distributed random number
    USE constants
    IMPLICIT NONE
    REAL :: random_Lorentzian

    random_Lorentzian=tan(random_uniform(0.,pi/2.))

  END FUNCTION random_Lorentzian

  FUNCTION random_Gaussian_pair(mean,sigma)
    
    ! Gets a pair of Gaussian random numbers
    USE constants
    IMPLICIT NONE
    REAL :: random_Gaussian_pair(2)
    REAL, INTENT(IN) :: mean, sigma
    REAL :: r, theta

    r=random_Rayleigh(sigma)
    theta=random_uniform(0.,twopi)

    ! Both of these numbers are Gaussian
    random_Gaussian_pair(1)=r*sin(theta)+mean
    random_Gaussian_pair(2)=r*cos(theta)+mean

  END FUNCTION random_Gaussian_pair

  FUNCTION random_Gaussian(mean,sigma)
    
    ! Gets a single Gaussian random number
    IMPLICIT NONE
    REAL :: random_Gaussian
    REAL, INTENT(IN) :: mean, sigma
    REAL :: G(2)

    ! This is wasteful as there is a second, independent Gaussian random number
    G=random_Gaussian_pair(mean,sigma)
    random_Gaussian=G(1)

  END FUNCTION random_Gaussian

  FUNCTION random_lognormal(mean_x,sigma_lnx)
    
    ! Gets a single Gaussian random number
    ! mean_x: <x>
    ! sigma_lnx: rms of the logarithm of x
    IMPLICIT NONE
    REAL :: random_lognormal
    REAL, INTENT(IN) :: mean_x, sigma_lnx
    REAL :: mu, sigma, G(2)

    sigma=sigma_lnx
    mu=log(mean_x)-0.5*sigma**2

    ! This is wasteful as there is a second, independent Gaussian random number
    G=random_Gaussian_pair(mu,sigma)
    random_lognormal=exp(G(1))

  END FUNCTION random_lognormal

  REAL FUNCTION random_exponential(mean)

    ! Produces a exponentially-distributed random number
    IMPLICIT NONE
    REAL, INTENT(IN) :: mean
    REAL, PARAMETER :: small=1e-10 ! Introducted because there will be problems here if log(0) is ever called
  
    random_exponential=-mean*log(random_uniform(small,1.))

  END FUNCTION random_exponential

  REAL FUNCTION random_polynomial(n)

    ! Generate a polynomailly distributed number [x:0->1]
    IMPLICIT NONE
    REAL, INTENT(IN) :: n

    random_polynomial=(random_uniform(0.,1.))**(1./(n+1))

  END FUNCTION random_polynomial

  REAL FUNCTION random_theta()

    ! A random spherical angle such that the space is equally populated
    IMPLICIT NONE

    random_theta=acos(random_uniform(-1.,1.))

  END FUNCTION random_theta

  REAL FUNCTION accept_reject(func,x1,x2,fmax)

    ! Simple one-dimensional accept-reject algorithm for drawing random numbers from func(x1->x2)
    ! TODO: Increase to n-dimensions
    ! TODO: Include more complicated bounding structure (at the moment it is just a box)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x1, x2
    REAL, INTENT(IN) :: fmax
    REAL :: x, y, f
   
    INTERFACE
       REAL FUNCTION func(x)
         REAL, INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE

    ! Try until you accept an x value
    DO

       ! Draw uniform random x value and function height
       x=random_uniform(x1,x2)
       y=random_uniform(0.,fmax)

       ! Evaulate the function at the random x value
       f=func(x)

       ! Decide whether or not to accept
       IF(f>fmax) THEN
          STOP 'ACCEPT_REJECT: Error, your function is not bounded by fmax'
       ELSE IF(y<=f) THEN
          accept_reject=x
          EXIT
       END IF

    END DO
    
  END FUNCTION accept_reject

END MODULE random_numbers
