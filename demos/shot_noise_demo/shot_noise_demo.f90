PROGRAM shot_noise_demo

  ! This test makes a random distribution of particles and then assigns them some weight
  ! It then compares the simple and the fancy (correct) shot noise calculations

  USE array_operations
  USE random_numbers
  USE constants
  USE simulations
  USE fft
  USE statistics
  USE field_operations
  
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:,:), x1(:,:), x2(:,:)
  REAL, ALLOCATABLE :: d(:,:,:), d1(:,:,:), d2(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: dk(:,:,:), dk_out(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: dk1(:,:,:), dk1_out(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: dk2(:,:,:), dk2_out(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: dc(:,:,:)
  REAL, ALLOCATABLE :: k(:), Pk(:), sig(:), w(:), w1(:), w2(:)
  REAL, ALLOCATABLE :: k_11(:), Pk_11(:), sig_11(:)
  REAL, ALLOCATABLE :: k_12(:), Pk_12(:), sig_12(:)
  REAL, ALLOCATABLE :: k_22(:), Pk_22(:), sig_22(:)
  INTEGER, ALLOCATABLE :: nm(:), nm_11(:), nm_12(:), nm_22(:)
  INTEGER :: i, itest
  REAL :: dbar, w1bar, w2bar, w1_total, w2_total
  REAL :: sn1, sn2, snn1, snn2
  REAL :: sn_11, sn_12, sn_22, snn_11, snn_12, snn_22

  INTEGER, PARAMETER :: m=128 ! Mesh size for density field
  INTEGER, PARAMETER :: n=128**3 ! Number of random particles
  REAL, PARAMETER :: L=1. ! Box size
  REAL, PARAMETER :: kmin=twopi/L ! Minimum wavenumber for power spectrum measurement
  REAL, PARAMETER :: kmax=pi*m/L ! Maximum wavenumber for power spectrum measurement
  INTEGER, PARAMETER :: nk=100 ! Number of bins for power spectrum measurement
  INTEGER, PARAMETER :: iseed=1 ! Seed for random-number generator
  INTEGER :: ibin=2 ! Binning scheme
  LOGICAL, PARAMETER :: sharpen_field=.TRUE. ! Sharpen the field or not to account for binning (essential!)
  INTEGER, PARAMETER :: irandom=1 ! Choose random distribution for weights (1 - uniform, 2 - exponential)
  REAL, PARAMETER :: wmin=0. ! Minimum particle uniform weight
  REAL, PARAMETER :: wmax=1. ! Maximum particle uniform weight
  REAL, PARAMETER :: mu=1. ! Mean for exponentially distributed weights
  LOGICAL, PARAMETER :: complex_fft=.TRUE. ! Set the FFT to use (c2c or r2c)
  !INTEGER, PARAMETER :: itest=1 ! Select test (1 - simple vs fancy shot noise, 2 - cross spectra shot noise)
  INTEGER, PARAMETER :: n1=128**3 ! Number of particles in field 1
  INTEGER, PARAMETER :: n2=64**3 ! Number of particles in field 2
  INTEGER, PARAMETER :: n12=32**3 ! Number of particles shared between fields 1 and 2
  LOGICAL, PARAMETER :: all=.TRUE.
  LOGICAL, PARAMETER :: periodic=.TRUE.
  LOGICAL, PARAMETER :: verbose=.TRUE.

  ! Welcome message
  WRITE(*,*)
  WRITE(*,*) 'SHOT_NOISE_DEMO: Running test'
  WRITE(*,*) 'SHOT_NOISE_DEMO: Number of particles:', n
  WRITE(*,*) 'SHOT_NOISE_DEMO: Box size:', L
  IF(irandom==1) THEN
     WRITE(*,*) 'SHOT_NOISE_DEMO: Uniform random weights'
     WRITE(*,*) 'SHOT_NOISE_DEMO: Minimum uniform weight:', wmin
     WRITE(*,*) 'SHOT_NOISE_DEMO: Minimum uniform weight:', wmax
  ELSE IF(irandom==2) THEN
     WRITE(*,*) 'SHOT_NOISE_DEMO: Exponentially randomly distributed weights'
     WRITE(*,*) 'SHOT_NOISE_DEMO: Mean for exponential distribution', mu
  ELSE
     STOP 'SHOT_NOISE_DEMO: Error, irandom assigned incorrectly'
  END IF
  WRITE(*,*) 'SHOT_NOISE_DEMO: Mesh size:', m
  WRITE(*,*) 'SHOT_NOISE_DEMO: Seed for random numbers:', iseed
  WRITE(*,*) 'SHOT_NOISE_DEMO: Field being sharpened to account for binning:', sharpen_field
  WRITE(*,*) 'SHOT_NOISE_DEMO: Binning scheme:', ibin
  WRITE(*,*) 'SHOT_NOISE_DEMO: Minimum wavenumber for power', kmin
  WRITE(*,*) 'SHOT_NOISE_DEMO: Maximum wavenumber for power', kmax
  WRITE(*,*) 'SHOT_NOISE_DEMO: Number of bins for power:', nk
  WRITE(*,*)

  ! Set the random number generator
  CALL RNG_set(iseed)

  ! Loop over tests
  DO itest=1,2

     IF(itest==1) THEN

        ! Make a random distribution of particles
        ALLOCATE(x(3,n))  
        CALL generate_randoms(x,n,L)

        ! Assign weights to particles
        ALLOCATE(w(n))
        w=1.
        DO i=1,n
           IF(irandom==1) THEN
              w(i)=random_uniform(wmin,wmax)
           ELSE IF(irandom==2) THEN
              w(i)=random_exponential(mu)
           END IF
        END DO

        ! Bin the particles and make a density field
        ALLOCATE(d(m,m,m))
        CALL particle_bin(x,n,L,w,d,m,ibin,all,periodic,verbose)
        dbar=SUM_DOUBLE(w,n)/REAL(m)**3
        d=d/dbar

        ! Do the Fourier Transforms
        IF(complex_fft) THEN
           ALLOCATE(dk(m,m,m),dk_out(m,m,m))
           dk=d
           CALL fft3(dk,dk_out,m,m,m,-1)
           dk=dk_out
           DEALLOCATE(dk_out)
        ELSE
           ALLOCATE(dc(m,m,m),dk(m/2+1,m,m))
           dc=d
           CALL fft3(dc,dk,m,m,m,-1)
           DEALLOCATE(dc)
        END IF
        DEALLOCATE(d)

        ! Sharpen the density field to account for the binning
        IF(sharpen_field) THEN
           IF(complex_fft) THEN
              CALL sharpen_k(dk,m,m,ibin)
           ELSE
              CALL sharpen_k(dk,m/2+1,m,ibin)
           END IF
        END IF

        ! Do the power spectrum computation
        !IF(complex_fft) THEN
        CALL compute_power_spectrum(dk,dk,m,L,kmin,kmax,nk,k,Pk,nm,sig)  
        !ELSE
        !   CALL compute_power_spectrum_real(dk,dk,m,L,kmin,kmax,nk,k,Pk,nm,sig)
        !END IF
        DEALLOCATE(dk)

        ! Write the data to file
        OPEN(7,file='power.dat')
        sn1=shot_noise_simple(L,INT8(n))
        sn2=shot_noise_mass(L,w,n)
        DO i=1,nk
           IF(nm(i)==0) CYCLE
           snn1=shot_noise_k(k(i),sn1)
           snn2=shot_noise_k(k(i),sn2)
           WRITE(7,*) k(i), Pk(i), snn1, nm(i), sig(i), snn2
        END DO
        CLOSE(7)

     ELSE IF(itest==2) THEN

        ! Make a random distribution of particles
        ALLOCATE(x1(3,n1),x2(3,n2))  
        CALL generate_randoms(x1,n1,L)
        CALL generate_randoms(x2,n2,L)

        ! Assign some random weights
        ALLOCATE(w1(n1),w2(n2))
        DO i=1,n1
           w1(i)=random_exponential(mu)
        END DO
        DO i=1,n2
           w2(i)=random_uniform(0.,wmax)
        END DO

        w1_total=sum_double(w1,n1)
        w2_total=sum_double(w2,n2)
        w1bar=w1_total/REAL(m)**3 ! Mean value of w1 in a cell
        w2bar=w2_total/REAL(m)**3 ! Mean value of w2 in a cell

        ! Make a some number of the particle positions to be shared between the fields
        DO i=1,n12
           x2(:,i)=x1(:,i)
        END DO

        ! Calculate the constant shot noise terms 
        sn_11=shot_noise(w1/w1_total,w1/w1_total,n1,L)
        sn_12=shot_noise(w1/w1_total,w2/w2_total,n12,L)
        sn_22=shot_noise(w2/w2_total,w2/w2_total,n2,L)

        ! Bin the particles and make density fields
        ALLOCATE(d1(m,m,m),d2(m,m,m))
        CALL particle_bin(x1,n1,L,w1,d1,m,ibin,all,periodic,verbose) ! Bin particles using w1 weights
        CALL particle_bin(x2,n2,L,w2,d2,m,ibin,all,periodic,verbose) ! Bin particles using w2 weights
        d1=d1/w1bar
        d2=d2/w2bar

        ! Do the Fourier Transforms
        ALLOCATE(dk1(m,m,m),dk2(m,m,m),dk1_out(m,m,m),dk2_out(m,m,m))
        dk1=d1
        dk2=d2
        CALL fft3(dk1,dk1_out,m,m,m,-1)
        CALL fft3(dk2,dk2_out,m,m,m,-1)
        dk1=dk1_out
        dk2=dk2_out
        DEALLOCATE(dk1_out,dk2_out,d1,d2)

        ! Sharpen the density field to account for the binning
        CALL sharpen_k(dk1,m,m,ibin)
        CALL sharpen_k(dk2,m,m,ibin)

        ! Do the power spectrum computation
        CALL compute_power_spectrum(dk1,dk1,m,L,kmin,kmax,nk,k_11,Pk_11,nm_11,sig_11)
        CALL compute_power_spectrum(dk1,dk2,m,L,kmin,kmax,nk,k_12,Pk_12,nm_12,sig_12)
        CALL compute_power_spectrum(dk2,dk2,m,L,kmin,kmax,nk,k_22,Pk_22,nm_22,sig_22)
        DEALLOCATE(dk1,dk2)

        ! Write the data to file
        OPEN(7,file='power_11.dat')
        OPEN(8,file='power_12.dat')
        OPEN(9,file='power_22.dat')
        DO i=1,nk
           IF(nm_11(i) .NE. 0) THEN
              snn_11=shot_noise_k(k_11(i),sn_11)
              WRITE(7,*) k_11(i), Pk_11(i), snn_11, nm_11(i), sig_11(i)
           END IF
           IF(nm_12(i) .NE. 0) THEN
              snn_12=shot_noise_k(k_12(i),sn_12)
              WRITE(8,*) k_12(i), Pk_12(i), snn_12, nm_12(i), sig_12(i)
           END IF
           IF(nm_22(i) .NE. 0) THEN
              snn_22=shot_noise_k(k_22(i),sn_22)
              WRITE(9,*) k_22(i), Pk_22(i), snn_22, nm_22(i), sig_22(i)
           END IF
        END DO
        CLOSE(7)
        CLOSE(8)
        CLOSE(9)

     ELSE

        STOP 'SHOT_NOISE_DEMO: Error, itest specified incorrectly'

     END IF

  END DO
  
CONTAINS

END PROGRAM shot_noise_demo
