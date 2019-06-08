PROGRAM cosmology_functions_demo

  USE constants
  USE array_operations
  USE cosmology_functions
  
  IMPLICIT NONE
  INTEGER :: icosmo
  TYPE(cosmology) :: cosm
  REAL :: a, z, k, r, sigv, sigv100
  REAL :: t1, t2, t3, t4, xi
  INTEGER :: i
  CHARACTER(len=256) :: cosmo

  LOGICAL, PARAMETER :: test_background=.TRUE.
  LOGICAL, PARAMETER :: test_power=.TRUE.
  LOGICAL, PARAMETER :: test_correlation=.FALSE.
  LOGICAL, PARAMETER :: test_sigma=.TRUE.
  LOGICAL, PARAMETER :: test_sigmaV=.TRUE.

  REAL, PARAMETER :: amin=1e-5
  REAL, PARAMETER :: amax=1
  INTEGER, PARAMETER :: na=256

  REAL, PARAMETER :: kmin=1e-4
  REAL, PARAMETER :: kmax=1e2
  INTEGER, PARAMETER :: nk=256

  REAL, PARAMETER :: rmin=1e-4
  REAL, PARAMETER :: rmax=1e3
  INTEGER, PARAMETER :: nr=256
  
  LOGICAL, PARAMETER :: verbose=.TRUE.

  !Initial white space
  WRITE(*,*)

  CALL get_command_argument(1,cosmo)
  IF(cosmo=='') THEN
     icosmo=-1
  ELSE
     READ(cosmo,*) icosmo
  END IF

  ! Initial white space
  WRITE(*,*)

  ! Assign the cosmological model
  CALL assign_cosmology(icosmo,cosm,verbose)

  ! Initialise the cosmology
  CALL init_cosmology(cosm)

  ! Write the cosmology to screen
  CALL print_cosmology(cosm)

  ! Test background quantities
  IF(test_background) THEN
     IF(verbose) WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Testing and writing background quantities'
     OPEN(10,file='Omegas.dat')
     OPEN(11,file='Hubble.dat')
     OPEN(12,file='growth.dat')
     OPEN(13,file='wa.dat')
     OPEN(14,file='distance.dat')
     OPEN(15,file='delta_c.dat')
     OPEN(16,file='Delta_v.dat')
     OPEN(17,file='time.dat')
     DO i=1,na
        a=progression_log(amin,amax,i,na)
        z=redshift_a(a)
        WRITE(10,*) a, &
             Omega_m(a,cosm), &
             Omega_c(a,cosm), &
             Omega_b(a,cosm), &
             Omega_r(a,cosm), &
             Omega_g(a,cosm), &
             Omega_nu(a,cosm), &
             Omega_v(a,cosm), &
             Omega_w(a,cosm), &
             Omega(a,cosm)
        WRITE(11,*) a, sqrt(Hubble2(a,cosm))
        WRITE(12,*) a, grow(a,cosm), ungrow(a,cosm), growth_rate(a,cosm), acc_growth(a,cosm), grow_Linder(a,cosm), grow_CPT(a,cosm), growth_rate_Linder(a,cosm)
        WRITE(13,*) a, w_de(a,cosm), w_de_total(a,cosm), w_eff(a,cosm)
        WRITE(14,*) a, comoving_distance(a,cosm), comoving_angular_distance(a,cosm), physical_angular_distance(a,cosm), luminosity_distance(a,cosm)
        WRITE(15,*) a, dc_NakamuraSuto(a,cosm), dc_Mead(a,cosm), dc_spherical(a,cosm)
        WRITE(16,*) a, Dv_BryanNorman(a,cosm), Dv_Mead(a,cosm), Dv_spherical(a,cosm)
        WRITE(17,*) a, cosmic_time(a,cosm), look_back_time(a,cosm)
     END DO
     CLOSE(10)
     CLOSE(11)
     CLOSE(12)
     CLOSE(13)
     CLOSE(14)
     CLOSE(15)
     CLOSE(16)
     CLOSE(17)
     IF(verbose) THEN
        WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Done'
        WRITE(*,*)
     END IF
  END IF

  ! Test power spectrum
  IF(test_power) THEN
     IF(verbose) WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Testing and writing power spectrum'
     OPEN(10,file='power.dat')
     DO i=1,nk
        k=progression_log(kmin,kmax,i,nk)
        WRITE(10,*) k, p_lin(k,1.,cosm), p_lin(k,0.5,cosm)
     END DO
     CLOSE(10)
     IF(verbose) THEN
        WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Done'
        WRITE(*,*)
     END IF
  END IF

  ! Test correlation function
  IF(test_correlation) THEN
     IF(verbose) WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Testing and writing correlation function'
     OPEN(10,file='correlation.dat')
     DO i=1,nk
        r=progression_log(rmin,rmax,i,nr)
        xi=xi_lin(r,1.,cosm)
        WRITE(*,*)  r, xi, 4.*pi*r**2*xi
        WRITE(10,*) r, xi, 4.*pi*r**2*xi
     END DO
     CLOSE(10)
     IF(verbose) THEN
        WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Done'
        WRITE(*,*)
     END IF
  END IF

  ! Test sigma(R)
  IF(test_sigma) THEN
     IF(verbose) WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Testing and writing sigma(R)'
     OPEN(10,file='sigma.dat')
     DO i=1,nr
        r=progression_log(rmin,rmax,i,nr)
        WRITE(10,*) r, sigma(r,1.,cosm), sigma(r,0.5,cosm)
     END DO
     CLOSE(10)
     IF(verbose) THEN
        WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Done'
        WRITE(*,*)
     END IF
  END IF

  ! Test sigmaV(R)
  IF(test_sigmav) THEN
     WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: Testing and writing sigma_v(R)'
     CALL CPU_TIME(t1)
     sigv=sigmaV(0.,1.0,cosm)
     CALL CPU_TIME(t2)
     sigv100=sigmaV(100.,1.0,cosm)
     CALL CPU_TIME(t3)
     WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: sigmaV(R=0 Mpc/h,   a=1.0):', sigv
     WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: sigmaV(R=0 Mpc/h,   a=1.0) time [s]:', t2-t1
     WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: sigmaV(R=100 Mpc/h, a=1.0):', sigv100
     WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: sigmaV(R=100 Mpc/h, a=1.0) time [s]:', t3-t2
     OPEN(10,file='sigmaV.dat')
     DO i=1,nr
        r=progression_log(rmin,rmax,i,nr)
        !WRITE(*,*)  r, sigmaV(r,1.0,cosm), sigmaV(r,0.5,cosm)
        WRITE(10,*) r, sigmaV(r,1.0,cosm), sigmaV(r,0.5,cosm)
     END DO
     CLOSE(10)
     CALL CPU_TIME(t4)
     WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: SigmaV total time [s]:', t4-t1  
     WRITE(*,*) 'COSMOLOGY_FUNCTIONS_DEMO: SigmaV done'
     WRITE(*,*)
  END IF
  
END PROGRAM cosmology_functions_demo