PROGRAM statistics_test

   USE statistics
   USE table_integer
   USE random_numbers
   USE array_operations

   IMPLICIT NONE
   REAL, ALLOCATABLE :: A(:), bins(:)
   INTEGER, ALLOCATABLE :: hist(:)
   CHARACTER(len=256) :: outfile
   REAL :: min, max, avg, std, x, dx, mu, sigma2
   INTEGER :: i
   INTEGER :: itest

   REAL, PARAMETER :: x0 = 0.
   REAL, PARAMETER :: sigma = 1.
   REAL, PARAMETER :: dmin = -5.
   REAL, PARAMETER :: dmax = 5.
   INTEGER, PARAMETER :: nh = 200
   INTEGER, PARAMETER :: m = 1000000
   INTEGER, PARAMETER :: iseed = 0
   INTEGER, PARAMETER :: n = 1000000

   WRITE (*, *)
   WRITE (*, *) 'STATISTICS_TEST: Select test'
   WRITE (*, *) '1 - Histogram'
   WRITE (*, *) '2 - Mean and variance'
   READ (*, *) itest
   WRITE (*, *)

   !Seed the random number generator
   CALL random_generator_seed(iseed)

   IF (itest == 1) THEN

      WRITE (*, *)
      WRITE (*, *) 'STATISTICS_TEST: Welcome'
      WRITE (*, *) 'STATISTICS_TEST: Calculating various one-point statistics'

      !Allocate the 'data' array
      ALLOCATE (A(m))

      !Fill the 'data' array
      WRITE (*, *) 'STATISTICS_TEST: Creating Gaussian random numbers'
      WRITE (*, *) 'STATISTICS_TEST: Mean value:', x0
      WRITE (*, *) 'STATISTICS_TEST: Standard deviation:', sigma
      DO i = 1, m
         A(i) = random_Gaussian(x0, sigma)
      END DO

      !Calculate and write out some useful statistics
      min = MINVAL(A)
      max = MAXVAL(A)
      avg = mean(A)
      std = sqrt(variance(A))
      WRITE (*, *) 'STATISTICS_TEST: Number of points:', m
      WRITE (*, *) 'STATISTICS_TEST: Minimum value:', REAL(min)
      WRITE (*, *) 'STATISTICS_TEST: Maximum value:', REAL(max)
      WRITE (*, *) 'STATISTICS_TEST: Mean value:', REAL(avg)
      WRITE (*, *) 'STATISTICS_TEST: Standard deviation:', REAL(std)
      WRITE (*, *)

      !Create the one-point PDF
      CALL histogram(dmin, dmax, bins, hist, nh, A)

      !Writes out the normalised histogram to a file
      outfile = 'histogram.dat'
      WRITE (*, *) 'STATISTICS_TEST: Writing data'
      OPEN (7, file=outfile)
      DO i = 1, nh
         x = (bins(i)+bins(i+1))/2.
         dx = bins(i+1)-bins(i)
         WRITE (7, *) x, hist(i)/(REAL(m)*dx)
      END DO
      CLOSE (7)

      WRITE (*, *) 'STATISTICS_TEST: Done'
      WRITE (*, *)

   ELSE IF (itest == 2) THEN

      WRITE (*, *) 'STATISTICS_TEST: Generating n uniform random numbers, n:', n
      ALLOCATE (a(n))
      DO i = 1, n
         a(i) = random_uniform(0., 1.)
      END DO

      !WRITE(*,*) 'Data:', a
      mu = mean(a)
      sigma2 = variance(a)
      WRITE (*, *) 'STATISTICS_TEST: Expected mean:', 0.5
      WRITE (*, *) 'STATISTICS_TEST: Expected variance:', 1./12.
      WRITE (*, *) 'STATISTICS_TEST: Expected standard deviation', sqrt(1./12.)
      WRITE (*, *) 'STATISTICS_TEST: Mean:', mu
      WRITE (*, *) 'STATISTICS_TEST: Variance:', sigma2
      WRITE (*, *) 'STATISTICS_TEST: Standard deviation:', sqrt(sigma2)
      WRITE (*, *)

   ELSE

      STOP 'STATISTICS_TEST: Error, itest not specified correctly'

   END IF

END PROGRAM statistics_test
