MODULE physics

  USE constants
  IMPLICIT NONE

CONTAINS

  REAL FUNCTION black_body_nu(nu,T)

    IMPLICIT NONE
    REAL, INTENT(IN) :: nu ! Frequency (Hz^-1)
    REAL, INTENT(IN) :: T  ! Black-body temperature [K]
    REAL :: a, x

    a=2.*h_Planck*nu**3/c_light**2

    x=h_Planck*nu/(kB*T)

    black_body_nu=a/(exp(x)-1.)

  END FUNCTION black_body_nu
  
END MODULE physics
