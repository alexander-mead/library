PROGRAM fft_test

  USE fft
  USE random_numbers
  USE logical_operations
  
  IMPLICIT NONE
  DOUBLE COMPLEX, ALLOCATABLE :: dc(:,:,:), dk1(:,:,:), dk2(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: d(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: s(:,:), sk(:,:)
  REAL, ALLOCATABLE :: r(:)
  COMPLEX, ALLOCATABLE :: rc3(:), rc3_out(:)
  DOUBLE COMPLEX, ALLOCATABLE :: rcd1(:), rcd2(:), rcd2_out(:)
  DOUBLE PRECISION, ALLOCATABLE :: rd1(:)
  REAL :: kx1, ky1, kx2, ky2, k1, k2, crap
  INTEGER :: i, j, k, m, ii, jj, itest
  CHARACTER(len=256) :: test

  INTEGER, PARAMETER :: iseed=1

  CALL get_command_argument(1,test)
  IF(test=='') THEN
     itest=-1
  ELSE
     READ(test,*) itest
  END IF
  ! Choose test

  IF(itest==-1) THEN
     WRITE(*,*) '1 - Test 3D transforms'
     WRITE(*,*) '2 - Test 1D transforms'
     WRITE(*,*) '3 - Check to see which array elements are the same in a real transform'
     READ(*,*) itest
     WRITE(*,*)
  END IF
 

  ! Seed the random number generator to create random data
  CALL RNG_set(iseed)

  ! 3D test
  IF(itest==1) THEN

     WRITE(*,*) 'FFT_TEST: Testing 3D transforms'
     
     ! Create a random data array
     m=32
     WRITE(*,*) 'FFT_TEST: Mesh size:', m
     WRITE(*,*)
     
     ALLOCATE(d(m,m,m))
     DO k=1,m
        DO j=1,m
           DO i=1,m
              d(i,j,k)=random_uniform(0.,1.)
           END DO
        END DO
     END DO

     WRITE(*,*) 'FFT_TEST: Random input data'
     WRITE(*,*) '==========================='
     DO i=1,m/2
        WRITE(*,fmt='(I10,F20.10)') i, d(i,1,1)
     END DO
     WRITE(*,*) '==========================='     
     WRITE(*,*) 'FFT_TEST: Done'
     WRITE(*,*)

     ! Do the real -> complex transform
     ALLOCATE(dk1(m/2+1,m,m))
     CALL fft3(d,dk1,m,m,m,-1)
     dk1=dk1/REAL(m**3)

     ! Do the complex -> complex transform
     ALLOCATE(dc(m,m,m),dk2(m,m,m))
     dc=d
     CALL fft3(dc,dk2,m,m,m,-1)
     dk2=dk2/REAL(m**3)

     WRITE(*,*) 'FFT_TEST: Checking forward 3D transforms'
     WRITE(*,*) '========================================'
     WRITE(*,*) '        i           r2c: real      r2c: imaginary           c2c: real      c2c: imaginary'
     DO i=1,m/2
        WRITE(*,fmt='(I10,4F20.10)') i, dk1(i,1,1), dk2(i,1,1)
     END DO
     WRITE(*,*) '========================================'
     WRITE(*,*) 'FFT_TEST: Done'
     WRITE(*,*)

     CALL fft3(d,dk1,m,m,m,1)

     CALL fft3(dk2,dc,m,m,m,1)

     ! Write out arrays for comparison
     WRITE(*,*) 'FFT_TEST: Checking reverse transforms'
     WRITE(*,*) '====================================='
     WRITE(*,*) '        i           r2c: real      r2c: imaginary           c2c: real      c2c: imaginary'
     DO i=1,m/2
        WRITE(*,fmt='(I10,4F20.10)') i, d(i,1,1), 0.0, dc(i,1,1)
     END DO
     WRITE(*,*) '====================================='
     WRITE(*,*) 'FFT_TEST: Done'
     WRITE(*,*)

  ELSE IF(itest==2) THEN

     ! Allocate the data array and populate with random numbers
     m=16
     WRITE(*,*) 'FFT_TEST: Testing 1D transforms'
     WRITE(*,*) 'FFT_TEST: Mesh size:', m
     WRITE(*,*)
     
     ALLOCATE(r(m))
     WRITE(*,*) 'FFT_TEST: Random input data'
     WRITE(*,*) '==========================='     
     DO i=1,m
        r(i)=random_uniform(0.,1.)
        WRITE(*,fmt='(I5,F16.8)') i, r(i)
     END DO
     WRITE(*,*) '==========================='
     WRITE(*,*) 'FFT_TEST: Done'
     WRITE(*,*)

     ! Allocate the double and complex array for the real -> complex transform
     ALLOCATE(rd1(m),rcd1(m/2+1))
     rd1=r
     CALL fft1(rd1,rcd1,m,-1)
     rcd1=rcd1/REAL(m)

     ! Allocate complex arrays for the complex -> complex transform
     ALLOCATE(rcd2(m),rcd2_out(m))
     rcd2=r
     CALL fft1(rcd2,rcd2_out,m,-1)
     rcd2=rcd2_out
     DEALLOCATE(rcd2_out)
     rcd2=rcd2/REAL(m)

     ALLOCATE(rc3(m),rc3_out(m))
     rc3=r
     CALL fft1(rc3,rc3_out,m,-1)
     rc3=rc3_out
     DEALLOCATE(rc3_out)
     rc3=rc3/REAL(m)

     ! Write out arrays for comparison
     WRITE(*,*) 'FFT_TEST: Checking forward transforms'
     WRITE(*,*) '====================================='
     !WRITE(*,*) '        i           r2c: real      r2c: imaginary           c2c: real      c2c: imaginary    single: real  single: imaginary'
     DO i=1,m/2+1
        WRITE(*,fmt='(I5,6F16.8)') i, rcd1(i), rcd2(i), rc3(i)
     END DO
     WRITE(*,*) '====================================='
     WRITE(*,*) 'FFT_TEST: Done'
     WRITE(*,*)

     CALL fft1(rd1,rcd1,m,1)

     ALLOCATE(rcd2_out(m))
     CALL fft1(rcd2,rcd2_out,m,1)
     rcd2=rcd2_out
     DEALLOCATE(rcd2_out)

     ALLOCATE(rc3_out(m))
     CALL fft1(rc3,rc3_out,m,1)
     rc3=rc3_out
     DEALLOCATE(rc3_out)

     ! Write out arrays for comparison
     WRITE(*,*) 'FFT_TEST: Checking reverse transforms'
     WRITE(*,*) '====================================='
     !WRITE(*,*) '        i           r2c: real      r2c: imaginary           c2c: real      c2c: imaginary'
     DO i=1,m
        WRITE(*,fmt='(I5,6F16.8)') i, rd1(i), 0.0, rcd2(i), rc3(i)
     END DO
     WRITE(*,*) '====================================='
     WRITE(*,*) 'FFT_TEST: Done'
     WRITE(*,*)

  ELSE IF(itest==3) THEN

     ! Allocate arrays and make random data
     m=6
     IF(odd(m)) STOP 'FFT_TEST: Error, m should be even'
     WRITE(*,*) 'FFT_TEST: Testing 2D transforms for identical elements'
     WRITE(*,*) 'FFT_TEST: Mesh size:', m
     WRITE(*,*)
     
     ALLOCATE(s(m,m),sk(m,m))
     DO j=1,m
        DO i=1,m
           s(i,j)=random_uniform(0.,1.)
        END DO
     END DO

     ! Do the FFT
     CALL fft2(s,sk,m,m,1)

     ! Write complex data to screen
     WRITE(*,*) '============================================'
     WRITE(*,*) 'FFT_TEST: Complex Fourier Transform elements'
     WRITE(*,*) '============================================'
     DO j=1,m
        DO i=1,m
           WRITE(*,*) i, j, sk(i,j)
        END DO
     END DO
     WRITE(*,*) '============================================'
     WRITE(*,*)

     ! Check for array elements that are equal
     WRITE(*,*) '============================='
     WRITE(*,*) 'Array elements that are equal'
     WRITE(*,*) '============================='
     DO j=1,m
        DO i=1,m/2+1
           DO jj=1,m
              DO ii=1,m/2+1
                 IF(i==ii .AND. j==jj) CYCLE
                 IF(sk(i,j)==CONJG(sk(ii,jj))) THEN
                    CALL k_fft(i,j,1,m,kx1,ky1,crap,k1,1.)
                    CALL k_fft(ii,jj,1,m,kx2,ky2,crap,k2,1.)
                    WRITE(*,*) i, j, 'is conjugate to', ii, jj
                    WRITE(*,*) 'Element 1:', sk(i,j)
                    WRITE(*,*) 'Element 2:', sk(ii,jj)
                    WRITE(*,*) 'Vector 1: kx:', kx1, 'ky:', ky1, '|k|:', k1
                    WRITE(*,*) 'Vector 2: ky:', kx2, 'ky:', ky2, '|k|:', k2
                    WRITE(*,*)
                 END IF
              END DO
           END DO
        END DO
     END DO
     WRITE(*,*) '============================='
     WRITE(*,*)

  ELSE

     STOP 'FFT_TEST: Error, itest not specified correctly'
     
  END IF

END PROGRAM fft_test
