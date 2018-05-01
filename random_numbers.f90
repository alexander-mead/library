MODULE random_numbers

  USE constants
  
  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE RNG_set(seed)

    !Seeds the RNG
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: seed
    INTEGER :: int, timearray(3)
    REAL :: rand !Necessary to define for ifort 

    WRITE(*,*) 'RNG_SET: Initialising random number generator'
    WRITE(*,*) 'RNG_SET: seed:', seed

    IF(seed==0) THEN
       !This fills the time array using the system clock!
       !If called within the same second the numbers will be identical!
       CALL itime(timeArray)
       !This then initialises the generator!
       int=FLOOR(rand(timeArray(1)+timeArray(2)+timeArray(3)))
    ELSE
       !In this case you can keep track of the seed
       int=FLOOR(rand(seed))
    END IF
    WRITE(*,*) 'RNG_SET: done'
    WRITE(*,*)

  END SUBROUTINE RNG_set

  FUNCTION random_integer(i1,i2)

    !Picks an integer with uniform random probability between i1 and i2 spaced with 1
    IMPLICIT NONE
    INTEGER :: random_integer
    INTEGER, INTENT(IN) :: i1, i2
    INTEGER :: n
    REAL :: ran
    REAL :: rand !Necessary to define for ifort

    !Range for the number
    n=1+i2-i1

    !Random number between 0 and 1
    !Do I really need both 'ran' and 'rand'?
    ran=rand(0)

    random_integer=i1-1+CEILING(ran*float(n))

  END FUNCTION random_integer

  FUNCTION random_uniform(x1,x2)

    !Produces a uniform random number between x1 and x2
    IMPLICIT NONE
    REAL :: random_uniform
    REAL, INTENT(IN) :: x1,x2
    REAL :: rand !I think this needs to be definted for ifort

    !Rand is some inbuilt function
    random_uniform=x1+(x2-x1)*(rand(0))

  END FUNCTION random_uniform

  FUNCTION random_Rayleigh(sigma)

    !Produces a Rayleigh-distributed random number
    IMPLICIT NONE
    REAL :: random_Rayleigh
    REAL, INTENT(IN) :: sigma
    REAL, PARAMETER :: small=1e-10

    !Problems if small=0. because log(0.) gets called sometimes
    random_Rayleigh=sigma*sqrt(-2.*log(random_uniform(small,1.)))

  END FUNCTION random_Rayleigh

  FUNCTION random_Lorentzian()

    !Produces a Lorentzian-distributed random number
    IMPLICIT NONE
    REAL :: random_Lorentzian

    random_Lorentzian=tan(random_uniform(0.,pi/2.))

  END FUNCTION random_Lorentzian

  FUNCTION random_Gaussian(mean,sigma)

    !This is wasteful as r*sin(theta) is also Gaussian (independantly)!
    !Could be converted to a fuction with two outputs

    IMPLICIT NONE
    REAL :: random_Gaussian
    REAL, INTENT(IN) :: mean, sigma
    REAL :: r, theta

    r=random_Rayleigh(sigma)
    theta=random_uniform(0.,2.*pi)

    !Both of these numbers are Gaussian
    !random_Gaussian=r*sin(theta)+mean
    random_Gaussian=r*cos(theta)+mean

  END FUNCTION random_Gaussian

  FUNCTION random_Poisson(mean)

    !Produces a Poisson-distributed random number
    IMPLICIT NONE
    REAL :: random_Poisson
    REAL, INTENT(IN) :: mean
    REAL, PARAMETER :: small=1e-10

    !small is introducted because there will be problems here if log(0) is ever called
    random_Poisson=-mean*log(random_uniform(small,1.))

  END FUNCTION random_Poisson

  FUNCTION random_polynomial(n)

    IMPLICIT NONE
    REAL :: random_polynomial, n

    random_polynomial=(random_uniform(0.,1.))**(1./(n+1))

  END FUNCTION random_polynomial

END MODULE random_numbers
