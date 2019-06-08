PROGRAM power_test

  USE field_operations
  USE random_numbers
  USE constants
  USE fft
  IMPLICIT NONE
  REAL, ALLOCATABLE :: k_array(:), Pk_array(:), sig_array(:)
  DOUBLE PRECISION, ALLOCATABLE :: d(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: dc(:,:,:), dk(:,:,:)
  INTEGER, ALLOCATABLE :: n_array(:)
  INTEGER :: i, j, k

  INTEGER, PARAMETER :: m=32
  REAL, PARAMETER :: L=1.
  INTEGER, PARAMETER :: nk=32
  REAL, PARAMETER :: kmin=twopi/L
  REAL, PARAMETER :: kmax=REAL(m)*pi/L

  WRITE(*,*)
  
  ! Create array of random data
  WRITE(*,*) 'Creating random data'
  ALLOCATE(d(m,m,m))
  DO k=1,m
     DO j=1,m
        DO i=1,m
           d(i,j,k)=random_uniform(0.,1.)
        END DO
     END DO
  END DO
  WRITE(*,*) 'Done'
  WRITE(*,*)

  DO j=1,2

     IF(j==1) THEN

        WRITE(*,*) 'POWER_TEST: Wasteful complex power calculation'
        ALLOCATE(dc(m,m,m))
        dc=d
        ALLOCATE(dk(m,m,m))
        CALL fft3(dc,dk,m,m,m,-1)        
        WRITE(*,*) 'POWER_TEST: Wasteful complex power calculation output'

     ELSE IF(j==2) THEN

        WRITE(*,*) 'POWER_TEST: Real power calculation'
        ALLOCATE(dk(m/2+1,m,m))
        CALL fft3(d,dk,m,m,m,-1)       
        WRITE(*,*) 'POWER_TEST: Real power calculation output'
        
     END IF

     CALL compute_power_spectrum(dk,dk,m,L,kmin,kmax,nk,k_array,Pk_array,n_array,sig_array)
     DEALLOCATE(dk)

     ! Write power to screen
     WRITE(*,*) '==============================================================================='
     WRITE(*,*) '           k                 P(k)   shot-noise               n            sigma'
     WRITE(*,*) '==============================================================================='
     IF(j==1) OPEN(7,file='power_complex.dat')
     IF(j==2) OPEN(7,file='power_real.dat')
     DO i=1,nk
        IF(n_array(i) .NE. 0) THEN
           WRITE(*,*) k_array(i), Pk_array(i), 0., n_array(i), sig_array(i)
           WRITE(7,*) k_array(i), Pk_array(i), 0., n_array(i), sig_array(i)
        END IF
     END DO
     CLOSE(7)
     WRITE(*,*) '==============================================================================='
     WRITE(*,*)

  END DO

END PROGRAM power_test
