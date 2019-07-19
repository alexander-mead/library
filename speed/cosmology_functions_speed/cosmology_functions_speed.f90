PROGRAM cosmology_speed

   USE array_operations
   USE logical_operations
   USE cosmology_functions
   IMPLICIT NONE
   REAL, ALLOCATABLE :: r(:), sig(:), sigV(:), nef(:)
   REAL :: t1, t2, crap, calc
   INTEGER :: i, itest
   TYPE(cosmology) :: cosm
   INTEGER :: icosmo=1
   REAL, PARAMETER :: a = 1.0
   REAL, PARAMETER :: rmin = 1e-4
   REAL, PARAMETER :: rmax = 1e3
   INTEGER, PARAMETER :: nr = 1024
   LOGICAL, PARAMETER :: verbose_test = .FALSE.

   ! Assign cosmology
   CALL read_command_argument(1, icosmo, 'Please specify cosmological model')
   CALL assign_cosmology(icosmo, cosm, verbose=.TRUE.)
   CALL init_cosmology(cosm)
   

   ! Fill array of R values for test
   CALL fill_array(log(rmin), log(rmax), r, nr)
   r = exp(r)

   ! Ensures that sigma(R) is initialised
   crap = sigma_all(rmin, a, cosm)

   WRITE(*,*) 'Scale factor for test:', a
   WRITE(*,*) 'Minimum R for test [Mpc/h]:', rmin
   WRITE(*,*) 'Maximum R for test [Mpc/h]:', rmax
   WRITE(*,*) 'Number of points in R for test:', nr
   WRITE(*,*)

   ALLOCATE(sig(nr), sigV(nr), nef(nr))

   DO itest = 1,3
      IF(itest==1) WRITE(*,*) 'Sigma test starting'
      IF(itest==2) WRITE(*,*) 'SigmaV test starting'
      IF(itest==3) WRITE(*,*) 'neff test starting'
      CALL CPU_TIME(t1)
      DO i = 1, nr
         IF(itest==1) THEN
            calc = sigma_all(r(i), a, cosm)
            sig(i) = calc
         END IF
         IF(itest==2) THEN
            sigV(i) = sigmaV(r(i), a, cosm)
            sigV(i) = calc
         END IF
         IF(itest==3) THEN
            calc = neff(r(i), a, cosm)
            nef(i) = calc
         END IF
         IF(verbose_test) WRITE(*,*) i, r(i), calc
      END DO
      CALL CPU_TIME(t2)
      IF(itest==1) WRITE(*,*) 'Sigma speed test done'
      IF(itest==2) WRITE(*,*) 'SigmaV speed test done'
      IF(itest==3) WRITE(*,*) 'neff speed test done'
      WRITE(*,*) 'Time for sigma test [s]:', t2-t1 
      WRITE(*,*)
   END DO

END PROGRAM cosmology_speed