MODULE random_numbers

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

  FUNCTION uniform(x1,x2)

    !Produces a uniform random number between x1 and x2
    IMPLICIT NONE
    REAL :: uniform
    REAL, INTENT(IN) :: x1,x2
    REAL :: rand !I think this needs to be definted for ifort

    !Rand is some inbuilt function
    uniform=x1+(x2-x1)*(rand(0))

  END FUNCTION uniform

  FUNCTION Rayleigh(sigma)

    !Produces a Rayleigh-distributed random number
    IMPLICIT NONE
    REAL :: Rayleigh
    REAL, INTENT(IN) :: sigma
    REAL, PARAMETER :: small=1e-10

    !Problems if small=0. because log(0.) gets called sometimes
    Rayleigh=sigma*sqrt(-2.*log(uniform(small,1.)))

  END FUNCTION Rayleigh

  FUNCTION Lorentzian()

    !Produces a Lorentzian-distributed random number
    IMPLICIT NONE
    REAL :: Lorentzian
    REAL, PARAMETER :: pi=3.141592654

    Lorentzian=tan(uniform(0.,pi/2.))

  END FUNCTION Lorentzian

  FUNCTION Gaussian(mean,sigma)

    !This is wasteful as r*sin(theta) is also Gaussian (independantly)!
    !Could be converted to a fuction with two outputs

    IMPLICIT NONE
    REAL :: Gaussian
    REAL, INTENT(IN) :: mean, sigma
    REAL :: r, theta
    REAL, PARAMETER :: pi=3.141592654

    r=Rayleigh(sigma)
    theta=uniform(0.,2.*pi)

    !Both of these numbers are Gaussian
    !   Gaussian=r*sin(theta)+mean
    Gaussian=r*cos(theta)+mean

  END FUNCTION Gaussian

  FUNCTION Poisson(mean)

    !Produces a Poisson-distributed random number
    IMPLICIT NONE
    REAL :: Poisson
    REAL, INTENT(IN) :: mean
    REAL, PARAMETER :: small=1e-10

    !There will be problems here is log(0) is ever called
    poisson=-mean*log(uniform(small,1.))

  END FUNCTION Poisson

  FUNCTION polynomial(n)

    IMPLICIT NONE
    REAL :: polynomial, n

    polynomial=(uniform(0.,1.))**(1./(n+1))

  END FUNCTION polynomial

END MODULE random_numbers
