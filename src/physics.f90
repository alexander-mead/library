MODULE physics

   USE constants

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: black_body_nu
   PUBLIC :: black_body_lambda
   PUBLIC :: wein_law_nu
   PUBLIC :: wein_law_lambda

CONTAINS

   REAL FUNCTION black_body_nu(nu, T)

      ! The radiance from a blackbody [Wm^-2 Hz^-1 Sr^-1]
      REAL, INTENT(IN) :: nu ! Frequency [Hz]
      REAL, INTENT(IN) :: T ! Black-body temperature [K]
      REAL :: a, x

      a = 2.*h_Planck*nu**3/c_light**2

      x = h_Planck*nu/(kB*T)

      black_body_nu = a/(exp(x)-1.)

   END FUNCTION black_body_nu

   REAL FUNCTION black_body_lambda(lambda, T)

      ! The radiance from a blackbody [Wm^-2 m^-1 Sr^-1]
      REAL, INTENT(IN) :: lambda ! Wavelength [m]
      REAL, INTENT(IN) :: T ! Black-body temperature [K]
      REAL :: nu

      nu = c_light/lambda

      black_body_lambda = black_body_nu(nu, T)*nu/lambda

   END FUNCTION black_body_lambda

   REAL FUNCTION wein_law_nu(T)

      ! The peak emission frequency from a blackbody [Hz]
      REAL, INTENT(IN) :: T ! Black-body temperature [K]
      REAL, PARAMETER :: a = 5.879e10 ! Wein constant [Hz/K]

      wein_law_nu = a*T

   END FUNCTION wein_law_nu

   REAL FUNCTION wein_law_lambda(T)

      ! The peak emission wavelength from a blackbody [m]
      REAL, INTENT(IN) :: T ! Black-body temperature [K]
      REAL, PARAMETER :: b = 2.8977729e-3 ! Wein constant [Km]

      wein_law_lambda = b/T

   END FUNCTION wein_law_lambda

END MODULE physics
