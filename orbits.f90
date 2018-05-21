MODULE orbits

  USE vectors
  IMPLICIT NONE
  REAL, PARAMETER :: G_orbit=1. !Gravitational constant for orbit problems
  REAL, PARAMETER :: M_Jupiter=0.000954265748 !Jupiter mass in Solar masses

CONTAINS

  REAL FUNCTION Kepler_period(a,M)

    !Calculates the period of an orbit by Kepler's third law
    !a - semi-major axis in AU
    !M - total mass in Solar masses
    IMPLICIT NONE
    REAL, INTENT(IN) :: a, M
    
    Kepler_period=sqrt(a**3./(G_orbit*M))

  END FUNCTION kepler_period

  REAL FUNCTION radius_theta(theta,e)

    IMPLICIT NONE
    REAL, INTENT(IN) :: theta, e

    radius_theta=1./((1.+e*cos(theta))**2)

  END FUNCTION radius_theta

  FUNCTION orbit_position(a,e,theta,delta)

    IMPLICIT NONE
    REAL :: orbit_position(3)
    REAL, INTENT(IN) :: a, e, theta, delta
    REAL :: l

    l=a*(1.-e**2.)

    orbit_position(1)=l*cos(theta)/(1.+e*cos(theta))
    orbit_position(2)=l*sin(theta)/(1.+e*cos(theta))
    orbit_position(3)=0.

    orbit_position=rotate_vector(orbit_position,zhat,delta)

  END FUNCTION orbit_position

  FUNCTION orbit_velocity(a,e,theta,delta,cm)

    IMPLICIT NONE
    REAL :: orbit_velocity(3)
    REAL, INTENT(IN) :: a, e, theta, delta, cm
    REAL :: l, b, w

    l=a*(1.-e**2.)
    b=l/sqrt(1.-e**2.)
    w=sqrt(G_orbit*cm/(a**3))

    orbit_velocity(1)=(-w*a*b/l)*sin(theta)
    orbit_velocity(2)=(w*a*b/l)*(cos(theta)+e)
    orbit_velocity(3)=0.

    orbit_velocity=rotate_vector(orbit_velocity,zhat,delta)

  END FUNCTION orbit_velocity

END MODULE orbits
