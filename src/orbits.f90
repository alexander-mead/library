MODULE orbits

   ! TODO: Shift everything to 3D
   ! TODO: Add inclination and longitude of ascending node
   ! TODO: Add eccentric and mean anomalies

   USE constants
   USE vectors

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: Kepler_period
   PUBLIC :: orbit_radius
   PUBLIC :: orbit_position
   PUBLIC :: orbit_velocity
   PUBLIC :: reduced_energy
   PUBLIC :: reduced_angular_momentum
   PUBLIC :: Laplace_Runge_Lenz
   PUBLIC :: semi_latus_rectum
   PUBLIC :: semi_major_axis
   PUBLIC :: semi_minor_axis
   PUBLIC :: polar_angle
   PUBLIC :: true_anomaly
   PUBLIC :: argument_of_periapsis
   PUBLIC :: eccentricity

   REAL, PARAMETER :: G_orbit = 1.              ! Gravitational constant for sensible units for orbits
   REAL, PARAMETER :: M_Jupiter = 9.54265748e-4 ! Jupiter mass [Msun]
   REAL, PARAMETER :: M_Earth = 3.003e-6        ! Earth mass [Msun]

CONTAINS

   !!! Misc !!!

   REAL FUNCTION Kepler_period(a, M)

      ! Period of an orbit [yr]
      REAL, INTENT(IN) :: a ! Semi-major axis [au]
      REAL, INTENT(IN) :: M ! Central mass [Msun]

      Kepler_period = sqrt(a**3/(G_orbit*M))

   END FUNCTION kepler_period

!!$  REAL FUNCTION radius_theta(theta,e)
!!$
!!$    ! Calculate r(theta) for a Keplerian orbit [au]
!!$    ! BUG: Should this not be l/(1+e*cos(t))
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: theta ! Plane-polar angle; true anomaly [rad]
!!$    REAL, INTENT(IN) :: e     ! Eccentricity
!!$
!!$    radius_theta=1./(1.+e*cos(theta))**2
!!$
!!$  END FUNCTION radius_theta

   !!! Calculate orbit properties from orbital elements !!!

   REAL FUNCTION orbit_radius(a, e, n)

      ! Calculate orbit x,y,z [au]
      REAL, INTENT(IN) :: a ! Semi-major axis [au]
      REAL, INTENT(IN) :: e ! Eccentricity
      REAL, INTENT(IN) :: n ! True anomaly [rad]
      REAL :: l

      l = a*(1.-e**2) ! Semi-latus rectum
      orbit_radius = l/(1.+e*cos(n))

   END FUNCTION orbit_radius

   FUNCTION orbit_position(a, e, n, w) RESULT(position)

      ! Planar orbit posititon x,y,z [au]
      ! TODO: Generalise to 3D
      REAL, INTENT(IN) :: a ! Semi-major axis [au]
      REAL, INTENT(IN) :: e ! Eccentricity
      REAL, INTENT(IN) :: n ! True anomaly [rad]
      REAL, INTENT(IN) :: w ! Argument of periapsis [rad]
      REAL :: position(3)
      REAL :: l, t

      IF (a < 0. .AND. e < 1.) STOP 'ORBIT_POSITION: Error, if a<0 then e>1 must be enforced'
      IF (e > 1. .AND. a > 0.) STOP 'ORBIT_POSITION: Error, if a<0 then e>1 must be enforced'

      t = n+w         ! Polar angle
      l = a*(1.-e**2) ! Semi-latus rectum

      position(1) = cos(t)*orbit_radius(a, e, n)
      position(2) = sin(t)*orbit_radius(a, e, n)
      position(3) = 0.

   END FUNCTION orbit_position

   FUNCTION orbit_velocity(a, e, n, w, M) RESULT(velocity)

      ! Planar orbit velocity vx,vy,vz [2piau/yr]
      ! TODO: Generalise to 3D
      REAL, INTENT(IN) :: a ! Semi-major axis [au]
      REAL, INTENT(IN) :: e ! Eccentricity
      REAL, INTENT(IN) :: n ! True anomaly [rad]
      REAL, INTENT(IN) :: w ! Argument of periapsis [rad]
      REAL, INTENT(IN) :: M ! Central mass [Msun]
      REAL :: velocity(3)
      REAL :: t, l, h

      IF (a < 0. .AND. e < 1.) STOP 'ORBIT_VELOCITY: Error, if a<0 then e>1 must be enforced'
      IF (e > 1. .AND. a > 0.) STOP 'ORBIT_VELOCITY: Error, if a<0 then e>1 must be enforced'

      t = n+w               ! Polar angle
      l = a*(1.-e**2)       ! Semi-latus rectum
      h = sqrt(G_orbit*M*l) ! Reduced angular momentum modulus

      velocity(1) = -sin(t)-e*(sin(t)*cos(n)-sin(n)*cos(t))
      velocity(2) = cos(t)+e*(cos(t)*cos(n)+sin(t)*sin(n))
      velocity(3) = 0.

      velocity = velocity*h/l

   END FUNCTION orbit_velocity

   !!! !!!

   !!! Calculate orbital elements from position, velocity and mass !!!

   REAL FUNCTION reduced_energy(x, v, M)

      ! Conserved reduced orbital energy
      ! This is negative for bound orbits and positive for unbound orbits
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity []
      REAL, INTENT(IN) :: M    ! Central mass [Msun]

      reduced_energy = 0.5*dot_product(v, v)-G_orbit*M/modulus(x)

   END FUNCTION reduced_energy

   FUNCTION reduced_angular_momentum(x, v) RESULT(h)

      ! Conserved reduced angular momentum vector
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity [au]
      REAL :: h(3)

      h = cross_product(x, v)

   END FUNCTION reduced_angular_momentum

   FUNCTION Laplace_Runge_Lenz(x, v, M) RESULT(LRL)

      ! Conserved Laplace-Runge-Lenz vector
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity [2piau/yr]
      REAL, INTENT(IN) :: M    ! Central mass [Msun]
      REAL :: LRL(3)
      REAL :: h(3)

      h = reduced_angular_momentum(x, v)
      LRL = cross_product(v, h)/(G_orbit*M)-x/modulus(x)

   END FUNCTION Laplace_Runge_Lenz

   REAL FUNCTION semi_latus_rectum(x, v, M)

      ! Calculate the semi-latus rectum [au]
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity []
      REAL, INTENT(IN) :: M    ! Central mass [Msun]
      REAL :: h(3)

      h = reduced_angular_momentum(x, v)
      semi_latus_rectum = dot_product(h, h)/(G_orbit*M)

   END FUNCTION semi_latus_rectum

   REAL FUNCTION semi_major_axis(x, v, M)

      ! Semi-major axis [au]
      ! Negative for unbound orbits
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity []
      REAL, INTENT(IN) :: M    ! Central mass [Msun]

      !semi_major_axis=semi_latus_rectum(x,v,M)/(1.-eccentricity(x,v,M)**2)
      semi_major_axis = -G_orbit*M/(2.*reduced_energy(x, v, M))

   END FUNCTION semi_major_axis

   REAL FUNCTION semi_minor_axis(x, v, M)

      ! Semi-minor axis [au]
      ! Assumes that the orbit is bound. Will break for unbound orbits
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity []
      REAL, INTENT(IN) :: M    ! Central mass [Msun]

      semi_minor_axis = semi_latus_rectum(x, v, M)/sqrt(1.-eccentricity(x, v, M)**2)

   END FUNCTION semi_minor_axis

   REAL FUNCTION eccentricity(x, v, M)

      ! Orbit eccentricity
      ! If e<1 the orbit is bound, if e>1 the orbit is unbound
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity [2piau/yr]
      REAL, INTENT(IN) :: M    ! Central mass [Msun]
      REAL :: e(3)

      e = Laplace_Runge_Lenz(x, v, M)
      eccentricity = modulus(e)

   END FUNCTION eccentricity

   REAL FUNCTION polar_angle(x)

      ! Plane-polar orbit angle [rad]
      ! TODO: Assumes orbit is in the x-y plane
      ! Gives result from 0 to 2pi
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]

      polar_angle = atan2(x(2), x(1))
      polar_angle = shift_angle_to_circle(polar_angle)

   END FUNCTION polar_angle

   REAL FUNCTION true_anomaly(x, v, M)

      ! Orbital true anomaly [rad]
      ! TODO: Assumes orbit is in the x-y plane
      ! Gives result from 0 to 2pi
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity [2piau/yr]
      REAL, INTENT(IN) :: M    ! Central mass [Msun]

      true_anomaly = polar_angle(x)-argument_of_periapsis(x, v, M)
      true_anomaly = shift_angle_to_circle(true_anomaly)

   END FUNCTION true_anomaly

   REAL FUNCTION argument_of_periapsis(x, v, M)

      ! Orbit argument of periapsis [rad]
      ! TODO: Assumes orbit is in the x-y plane
      ! Gives result from 0 to 2pi
      REAL, INTENT(IN) :: x(3) ! Orbit position [au]
      REAL, INTENT(IN) :: v(3) ! Orbit velocity [2piau/yr]
      REAL, INTENT(IN) :: M    ! Central mass [Msun]
      REAL :: e(3)

      e = Laplace_Runge_Lenz(x, v, M)
      IF (modulus(e) == 0.) THEN
         argument_of_periapsis = 0. ! Reference direction is then \hat{x}
      ELSE
         argument_of_periapsis = atan2(e(2), e(1))
         argument_of_periapsis = shift_angle_to_circle(argument_of_periapsis)
      END IF

   END FUNCTION argument_of_periapsis

   REAL FUNCTION inclination(x, v)

      REAL, INTENT(IN) :: x(3), v(3)
      REAL :: h(3)

      h = reduced_angular_momentum(x, v)
      inclination = acos(dot_product(h, zhat))
   
   END FUNCTION inclination

   !!! !!!

END MODULE orbits
