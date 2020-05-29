PROGRAM GRF_test

  USE constants
  USE random_numbers
  USE camb_stuff
  USE field_operations
  USE fft

  ! Test the generation of a Gaussian random field
  IMPLICIT NONE
  REAL, ALLOCATABLE :: d(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: dd(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: dk1(:,:,:), dk2(:,:,:)
  REAL, ALLOCATABLE :: k(:), Pk(:), sig(:)
  INTEGER, ALLOCATABLE :: n(:)
  INTEGER :: i, nz, ncamb
  CHARACTER(len=256) :: infile, outfile

  !Parameters
  INTEGER, PARAMETER :: iseed=0 ! Seed for random numbers (seed 1 is weird)
  INTEGER, PARAMETER :: m=128 ! Mesh size for field
  REAL, PARAMETER :: L=1000. ! Length for box
  REAL, PARAMETER :: shot=0. ! Shot noise is zero because there are no particles
  LOGICAL, PARAMETER :: use_average=.FALSE. ! Use average modes or actually make GRF
  REAL, PARAMETER :: kmin=twopi/L ! Minimum wavenumber for power spectrum
  REAL, PARAMETER :: kmax=REAL(m)*pi/L ! Maximum wavenumber for power spectrum
  INTEGER, PARAMETER :: nk=100 ! Number of k-space points for power spectrum
  LOGICAL, PARAMETER :: wasteful=.FALSE. ! Use wasteful FFTW or not

  ! Initial white space
  WRITE(*,*)

  ! Write information to screen
  WRITE(*,*) 'GRF_TEST: Gaussian Random Field test'
  WRITE(*,*) 'GRF_TEST: Random seed:', iseed
  WRITE(*,*) 'GRF_TEST: Mesh size:', m
  WRITE(*,*) 'GRF_TEST: Box size [Mpc/h]:', L
  WRITE(*,*)

  ! Set the random-number generator
  !CALL RNG_set(iseed)
  CALL RNG_seed(iseed)

  ! Read in a CAMB power spectrum
  infile='LCDM_matterpower.txt'
  CALL read_CAMB_Pk(k,Pk,ncamb,infile)
  k=log(k)
  Pk=log(Pk)

  ! Make the Gaussian random field
  ALLOCATE(d(m,m,m))
  CALL make_Gaussian_random_field(d,m,L,k,Pk,ncamb,use_average)

  ! Deallocate the CAMB power arrays
  DEALLOCATE(k,Pk)

  ! Write the field out to disk
  outfile='data/field_slice.dat'
  nz=1
  CALL write_3D_field_projection_ascii(d,m,L,nz,outfile)

  ! Do either the wasteful complex -> complex transform or the better real -> complex
  IF(wasteful) THEN

     !Transfer the field into Fourier space
     ALLOCATE(dk1(m,m,m),dk2(m,m,m))
     dk1=d
     CALL fft3(dk1,dk2,m,m,m,-1)
     dk1=dk2
     DEALLOCATE(dk2)

     ! Compute the power spectrum  
     CALL compute_power_spectrum(dk1,dk1,m,L,kmin,kmax,nk,k,Pk,n,sig)

  ELSE

     ALLOCATE(dd(m,m,m),dk1(m/2+1,m,m))
     dd=d
     CALL fft3(dd,dk1,m,m,m,-1)
     !CALL compute_power_spectrum_real(dk1,dk1,m,L,kmin,kmax,nk,k,Pk,n,sig)
     CALL compute_power_spectrum(dk1,dk1,m,L,kmin,kmax,nk,k,Pk,n,sig)

  END IF

  ! Write out power spectrum to file
  OPEN(7,file='data/power.dat')
  DO i=1,nk
     WRITE(7,*) k(i), Pk(i), shot, n(i), sig(i)
  END DO
  CLOSE(7)
  
END PROGRAM GRF_test
