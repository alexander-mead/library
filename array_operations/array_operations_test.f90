PROGRAM array_operations_test

  USE array_operations
  
  IMPLICIT NONE
  REAL, ALLOCATABLE :: a(:), b(:), c(:,:,:), d(:), bins(:)
  REAL :: min, max
  INTEGER :: na, i, nb, j, k, itest, nc, n, ilog
  INTEGER, ALLOCATABLE :: hist(:)
  
  WRITE(*,*)
  WRITE(*,*) 'Array tools'
  WRITE(*,*) '1 - Test reverse'
  WRITE(*,*) '2 - Test reduce'
  WRITE(*,*) '3 - Test splay'
  WRITE(*,*) '4 - Test binning'
  WRITE(*,*) '5 - Test fill'
  WRITE(*,*) '6 - Test amputate'
  READ(*,*) itest
  WRITE(*,*)

  IF(itest==1) THEN

     WRITE(*,*) 'Dimensions of array to reverse (filled up by i)'
     READ(*,*) na

     ALLOCATE(a(na))

     WRITE(*,*) 'Original'
     DO i=1,na
        a(i)=i
        WRITE(*,*) a(i)
     END DO

     CALL reverse(a,na)

     WRITE(*,*) 'Reversed'
     DO i=1,na
        WRITE(*,*) a(i)
     END DO

  ELSE IF(itest==2) THEN

     WRITE(*,*) 'Original array size'
     READ(*,*) na
     WRITE(*,*) 'Size to be reduced to'
     READ(*,*) nb

     ALLOCATE(a(na),b(nb))

     WRITE(*,*) 'Original'
     DO i=1,na
        a(i)=i
        WRITE(*,*) a(i)
     END DO

     CALL reduce(a,na,b,nb)

     WRITE(*,*) 'Reduced'
     DO i=1,nb
        WRITE(*,*) b(i)
     END DO
     WRITE(*,*)

  ELSE IF(itest==3) THEN

     WRITE(*,*) 'Splay array dimensions (it is a cube)'
     READ(*,*) nc

     ALLOCATE(c(nc,nc,nc))

     WRITE(*,*) 'Original cube'
     DO i=1,nc
        DO j=1,nc
           DO k=1,nc
              c(i,j,k)=i*j*k
              WRITE(*,*) i, j, k, i*j*k
           END DO
        END DO
     END DO

     d=splay(c,nc,nc,nc)

     WRITE(*,*) 'Splayed result'
     DO i=1,nc**3 
        WRITE(*,*) i, d(i)
     END DO
     WRITE(*,*)

  ELSE IF(itest==4) THEN

     WRITE(*,*) 'Dimension of array to bin'
     READ(*,*) na

     ALLOCATE(a(na))

     DO i=1,na
        a(i)=exp(-float(i))
     END DO

     WRITE(*,*) 'Number of bins:'
     READ(*,*) n

     CALL binning(a,MINVAL(a),MAXVAL(a),n,b,d,0)

     DO i=1,n
        WRITE(*,*) b(i), d(i)
     END DO

  ELSE IF(itest==5) THEN

     WRITE(*,*) 'min'
     READ(*,*) min
     WRITE(*,*) 'max'
     READ(*,*) max
     WRITE(*,*) 'number'
     READ(*,*) n
     WRITE(*,*) '0 - Linear fill'
     WRITE(*,*) '1 - Log fill'
     READ(*,*) ilog
     WRITE(*,*)

     IF(ilog==0) THEN
        CALL fill_array(min,max,a,n)
     ELSE IF(ilog==1) THEN
        CALL fill_array(log(min),log(max),a,n)
        a=exp(a)
     ELSE
        STOP
     END IF

     DO i=1,n
        WRITE(*,*) i, a(i)
     END DO
     WRITE(*,*)

  ELSE IF(itest==6) THEN

     CALL fill_array(1.,10.,a,10)

     DO i=1,n
        WRITE(*,*) i, a(i)
     END DO

     CALL amputate(a,10,5)

     DO i=1,n
        WRITE(*,*) i, a(i)
     END DO   

!!$  ELSE IF(itest==6) THEN
!!$
!!$     n=10
!!$     ALLOCATE(a(n))
!!$
!!$     a(1)=1.
!!$     a(2)=3.
!!$     a(3)=2.
!!$     a(4)=5.
!!$     a(5)=1.
!!$     a(6)=2.
!!$     a(7)=4.
!!$     a(8)=2.
!!$     a(9)=1.
!!$     a(10)=2.
!!$
!!$     WRITE(*,*) 'Data:', a
!!$     WRITE(*,*) 'Mean:', mean(a,n)
!!$     WRITE(*,*) 'Variance:', variance(a,n)
!!$     WRITE(*,*)
!!$
!!$  ELSE IF(itest==7) THEN
!!$
!!$     WRITE(*,*)
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Welcome'
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Calculating various one-point statistics'
!!$     WRITE(*,*)
!!$
!!$     !Seed the random number generator
!!$     CALL RNG_set(iseed)
!!$
!!$     !Allocate the 'data' array
!!$     ALLOCATE(A(m))
!!$
!!$     !Fill the 'data' array
!!$     DO i=1,m
!!$        A(i)=Gaussian(x0,sigma)
!!$     END DO
!!$
!!$     !Calculate and write out some useful statistics
!!$     min=MINVAL(A)
!!$     max=MAXVAL(A)
!!$     avg=mean(A,m)
!!$     std=sqrt(variance(A,m))
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Number of points:', m
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Minimum value:', REAL(min)
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Maximum value:', REAL(max)
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Mean value:', REAL(avg)
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Standard deviation:', REAL(std)
!!$     WRITE(*,*)
!!$
!!$     !Create the one-point PDF
!!$     CALL histogram(dmin,dmax,bins,hist,nh,A,m)
!!$
!!$     !Writes out the normalised histogram to a file
!!$     outfile='histogram.dat'
!!$     OPEN(7,file=outfile)
!!$     DO i=1,nh
!!$        x=(bins(i)+bins(i+1))/2.
!!$        dx=bins(i+1)-bins(i)
!!$        WRITE(7,*) x, hist(i)/(REAL(m)*dx)
!!$     END DO
!!$     CLOSE(7)
!!$
!!$     WRITE(*,*) 'HISTOGRAM_TEST: Done'
!!$     WRITE(*,*)

  ELSE

     STOP 'Test not specified correctly'

  END IF

END PROGRAM array_operations_test
