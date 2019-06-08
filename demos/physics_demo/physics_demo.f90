PROGRAM black_body_test

  USE constants
  USE physics
  USE cosmology_functions
  IMPLICIT NONE
  REAL :: nu, z, R, T
  INTEGER :: icos
  TYPE(cosmology) :: cosm
  REAL :: a
  REAL :: I
  LOGICAL, PARAMETER :: verbose=.TRUE.

  ! Initial white space
  WRITE(*,*)

  icos=1
  CALL assign_cosmology(icos,cosm,verbose)
  CALL init_cosmology(cosm)

  z=0.5 ! Source redshift
  T=100. ! Source tempertaure [K]
  R=0.1 ! Source physical size [Mpc]
  nu=353.e9 ! Observing frequency [Hz]

  WRITE(*,*) 'Source properties'
  WRITE(*,*) 'Redshift:', z
  WRITE(*,*) 'Physical size [Mpc/h]:', R
  WRITE(*,*) 'Observing frequency [GHz]:', nu/1e9
  WRITE(*,*)

  a=scale_factor_z(z)

  I=black_body_nu((1.+z)*nu,T)*(1.+z)*R**2/luminosity_distance(a,cosm)**2

  WRITE(*,*) 'Observed properties'
  WRITE(*,*) 'Irradiance [Jy]:', I/Jansky
  WRITE(*,*)
  
END PROGRAM black_body_test
