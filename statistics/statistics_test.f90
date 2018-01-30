PROGRAM statistics_test

  USE statistics
  USE table_integer
  USE random_numbers
  USE array_operations

  IMPLICIT NONE
  REAL, ALLOCATABLE :: A(:), bins(:)
  INTEGER, ALLOCATABLE :: hist(:)
  CHARACTER(len=256) :: infile, outfile, mesh, min_value, max_value, number_of_bins
  REAL :: min, max, avg, std, x, dx
  LOGICAL :: lexist
  INTEGER :: i, n
  INTEGER :: itest

  REAL, PARAMETER :: x0=0.
  REAL, PARAMETER :: sigma=1.
  REAL, PARAMETER :: dmin=-5.
  REAL, PARAMETER :: dmax=5.
  INTEGER, PARAMETER :: nh=200
  INTEGER, PARAMETER :: m=1000000
  INTEGER, PARAMETER :: iseed=0

  WRITE(*,*) 'Select test'
  WRITE(*,*) '1 - Histogram'
  WRITE(*,*) '2 - Mean and variance'
  READ(*,*) itest
  WRITE(*,*)

  IF(itest==1) THEN

     WRITE(*,*)
     WRITE(*,*) 'HISTOGRAM_TEST: Welcome'
     WRITE(*,*) 'HISTOGRAM_TEST: Calculating various one-point statistics'
     WRITE(*,*)

     !Seed the random number generator
     CALL RNG_set(iseed)

     !Allocate the 'data' array
     ALLOCATE(A(m))

     !Fill the 'data' array
     DO i=1,m
        A(i)=Gaussian(x0,sigma)
     END DO

     !Calculate and write out some useful statistics
     min=MINVAL(A)
     max=MAXVAL(A)
     avg=mean(A,m)
     std=sqrt(variance(A,m))
     WRITE(*,*) 'HISTOGRAM_TEST: Number of points:', m
     WRITE(*,*) 'HISTOGRAM_TEST: Minimum value:', REAL(min)
     WRITE(*,*) 'HISTOGRAM_TEST: Maximum value:', REAL(max)
     WRITE(*,*) 'HISTOGRAM_TEST: Mean value:', REAL(avg)
     WRITE(*,*) 'HISTOGRAM_TEST: Standard deviation:', REAL(std)
     WRITE(*,*)

     !Create the one-point PDF
     CALL histogram(dmin,dmax,bins,hist,nh,A,m)

     !Writes out the normalised histogram to a file
     outfile='histogram.dat'
     OPEN(7,file=outfile)
     DO i=1,nh
        x=(bins(i)+bins(i+1))/2.
        dx=bins(i+1)-bins(i)
        WRITE(7,*) x, hist(i)/(REAL(m)*dx)
     END DO
     CLOSE(7)

     WRITE(*,*) 'HISTOGRAM_TEST: Done'
     WRITE(*,*)

  ELSE IF(itest==2) THEN

     n=10
     ALLOCATE(a(n))

     a(1)=1.
     a(2)=3.
     a(3)=2.
     a(4)=5.
     a(5)=1.
     a(6)=2.
     a(7)=4.
     a(8)=2.
     a(9)=1.
     a(10)=2.

     WRITE(*,*) 'Data:', a
     WRITE(*,*) 'Mean:', mean(a,n)
     WRITE(*,*) 'Variance:', variance(a,n)
     WRITE(*,*)

  ELSE

     STOP 'STATISTICS_TEST: Error, itest not specified correctly'

  END IF

END PROGRAM statistics_test
