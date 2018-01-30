PROGRAM fft_test

  USE fft
  IMPLICIT NONE
  DOUBLE COMPLEX, ALLOCATABLE :: d(:,:,:)
  INTEGER :: i, m

  m=50

  ALLOCATE(d(m,m,m))
  d=1.
  
  CALL fft3(d,d,m,m,m,-1)
  
  DO i=1,10
     WRITE(*,*) i, d(i,1,1)
  END DO

END PROGRAM fft_test
