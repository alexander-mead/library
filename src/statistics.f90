MODULE statistics

   USE basic_operations
   USE array_operations
   USE sorting

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: mean
   PUBLIC :: geometric_mean
   PUBLIC :: harmonic_mean
   PUBLIC :: polynomial_mean
   PUBLIC :: weighted_mean
   PUBLIC :: variance
   PUBLIC :: standard_deviation
   PUBLIC :: root_mean_square
   PUBLIC :: histogram
   PUBLIC :: advanced_histogram
   PUBLIC :: calculate_confidence
   PUBLIC :: parameter_probability
   PUBLIC :: percentile
   PUBLIC :: percentiles

   ! Parameters
   INTEGER, PARAMETER :: percentile_sort = isort_bubble

CONTAINS

   FUNCTION percentiles(x, f) 
  
      ! For input distribution or measuments x(n) calculates the x value above which a fraction 'f' of x(n) lie
      ! Can be used to calculate interquartile ranges and such things
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: f(:)
      REAL :: percentiles(size(f))
      REAL, ALLOCATABLE :: y(:)
      INTEGER :: i1, i2, nx, j
      REAL :: in
      INTEGER, PARAMETER :: isort = percentile_sort

      y = x
      CALL sort(y, isort)
      nx = size(y)

      DO j = 1, size(f)
         IF (.NOT. between(f(j), 0., 1.)) THEN
            ERROR STOP 'PERCENTILE: Error, f must be between 0 and 1'
         ELSE IF (f(j) == 0.) THEN
            percentiles(j) = minval(x)
         ELSE IF (f(j) == 1.) THEN
            percentiles(j) = maxval(x)
         ELSE        
            in = 1.+f(j)*(nx-1)
            i1 = floor(in)
            i2 = ceiling(in)
            percentiles(j) = (y(i1)*(i2-in)+y(i2)*(in-i1)) ! Linear interpolation
         END IF
      END DO

   END FUNCTION percentiles

   REAL FUNCTION percentile(x, f)
  
      ! For input distribution or measuments x(n) calculates the x value above which a fraction 'f' of x(n) lie
      ! Can be used to calculate interquartile ranges and such things
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: f
      REAL, ALLOCATABLE :: y(:)
      INTEGER :: i1, i2, n
      REAL :: in
      INTEGER, PARAMETER :: isort = percentile_sort

      IF (.NOT. between(f, 0., 1.)) THEN
         ERROR STOP 'PERCENTILE: Error, f must be between 0 and 1'
      ELSE IF (f == 0.) THEN
         percentile = minval(x)
      ELSE IF (f == 1.) THEN
         percentile = maxval(x)
      ELSE
         y = x
         CALL sort(y, isort)
         n = size(y)
         in = 1.+f*(n-1)
         i1 = floor(in)
         i2 = ceiling(in)
         percentile = (y(i1)*(i2-in)+y(i2)*(in-i1)) ! Linear interpolation
      END IF

   END FUNCTION percentile

   REAL FUNCTION mean(x)

      REAL, INTENT(IN) :: x(:)

      mean = sum(x)/size(x)

   END FUNCTION mean

   REAL FUNCTION geometric_mean(x)

      REAL, INTENT(IN) :: x(:)

      !geometric_mean = product(x)**(1./size(x)) ! Probably insane to do this
      geometric_mean = exp(mean(log(x)))

   END FUNCTION geometric_mean

   REAL FUNCTION harmonic_mean(x)

      REAL, INTENT(IN) :: x(:)

      !harmonic_mean = 1./mean(1./x)
      harmonic_mean = polynomial_mean(x, -1)

   END FUNCTION harmonic_mean

   REAL FUNCTION polynomial_mean(x, n)

      REAL, INTENT(IN) :: x(:)
      INTEGER, INTENT(IN) :: n

      polynomial_mean = mean(x**n)**n

   END FUNCTION polynomial_mean

   REAL FUNCTION weighted_mean(x, w)

      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: w(:)

      IF (size(x) .NE. size(w)) STOP 'WEIGHTED_MEAN: Error, data and weight arrays must be the same size'
      weighted_mean = sum(w*x)/sum(w)

   END FUNCTION weighted_mean

   REAL FUNCTION variance(x)

      ! Note that depending on the application one might want to use n-1 here
      ! Difference is between sample variance and population varaince
      ! See https://en.wikipedia.org/wiki/Bessel%27s_correction
      ! Sample is simply the variance of your sample 1/N
      ! However, if you are estimating the population variance then you need 1/(N-1)
      REAL, INTENT(IN) :: x(:)

      variance = mean((x-mean(x))**2)
      !variance = mean(x**2)-mean(x)**2 ! Could also use this

   END FUNCTION variance

   REAL FUNCTION standard_deviation(x)

      ! Square root of variance
      REAL, INTENT(IN) :: x(:)

      standard_deviation = sqrt(variance(x))

   END FUNCTION standard_deviation

   REAL FUNCTION root_mean_square(x)

      ! Same as the standard deviation if the mean is zero
      REAL, INTENT(IN) :: x(:)

      root_mean_square = sqrt(mean(x**2))

   END FUNCTION root_mean_square

   SUBROUTINE histogram(xmin, xmax, x, hist, n, data)

      USE table_integer
      REAL, INTENT(IN) :: xmin ! Minimum x value
      REAL, INTENT(IN) :: xmax ! Maximum x value
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)       ! Output array of bin edges, size n+1
      INTEGER, ALLOCATABLE, INTENT(OUT) :: hist(:) ! Output integer array of bin counts, size n
      INTEGER, INTENT(IN) :: n    ! Number of bins
      REAL, INTENT(IN) :: data(:) ! Data to be binned
      INTEGER :: i, j, m

      m = size(data)

      WRITE (*, *) 'HISTOGRAM: Assiging arrays'

      ! Fill the table for the xrange and allocate the histogram array
      CALL fill_array(xmin, xmax, x, n+1)

      ! Set the histogram to zero
      IF (ALLOCATED(hist)) DEALLOCATE (hist)
      ALLOCATE (hist(n))
      hist = 0

      WRITE (*, *) 'HISTOGRAM: Constructing histogram'

      ! Make the histogram from the data
      DO i = 1, m
         IF (data(i) < xmin .OR. data(i) > xmax) THEN
            CYCLE
         ELSE
            ! TODO: Check ifind_linear is correct here
            j = find_table_integer(data(i), x, ifind_linear)
            hist(j) = hist(j)+1
         END IF
      END DO

      WRITE (*, *) 'HISTOGRAM: Fraction of data assigned to histogram:', real(sum(hist))/real(m)
      WRITE (*, *) 'HISTOGRAM: Done'
      WRITE (*, *)

   END SUBROUTINE histogram

   SUBROUTINE advanced_histogram(x_data, y_data, w_data, x_bins, y_bins, w_bins)

      USE table_integer
      REAL, INTENT(IN) :: x_data(:) ! x values for data (can be unordered)
      REAL, INTENT(IN) :: y_data(:) ! y values for data (can be unordered)
      REAL, INTENT(IN) :: w_data(:) ! weight for data
      REAL, INTENT(IN) :: x_bins(:) ! Ordered array of bin edges
      REAL, ALLOCATABLE, INTENT(OUT) :: y_bins(:) ! Output histogram heights
      REAL, ALLOCATABLE, INTENT(OUT) :: w_bins(:) ! Output sum of weights entering histogram
      INTEGER :: i_data, i_bins, n_data, n_bins
      INTEGER, PARAMETER :: ifind = ifind_split

      n_data = size(x_data)
      IF (n_data /= size(y_data)) STOP 'ADVANCED_HISTOGRAM: Error, x and y data should be the same size'
      IF (n_data /= size(w_data)) STOP 'ADVANCED_HISTOGRAM: Error, x and data weights should be the same size'

      n_bins = size(x_bins)-1
      ALLOCATE(y_bins(n_bins), w_bins(n_bins))
      y_bins = 0.
      w_bins = 0.

      DO i_data = 1, n_data
         i_bins = find_table_integer(x_data(i_data), x_bins, ifind)
         IF (i_bins > 0 .AND. i_bins <= n_bins) THEN
            y_bins(i_bins) = y_bins(i_bins)+y_data(i_data)*w_data(i_data)
            w_bins(i_bins) = w_bins(i_bins)+w_data(i_data)
         END IF
      END DO

   END SUBROUTINE advanced_histogram

   SUBROUTINE cumulative_distribution(x, p, c)

      USE calculus_table
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: p(:)
      REAL, INTENT(OUT) :: c(:)
      REAL :: norm
      INTEGER :: i, n
      INTEGER, PARAMETER :: iorder=0 ! Zeroth-order (histogram) integration and plenty of points are best here

      n = size(x)
      IF (n /= size(p) .OR. n /= size(c)) STOP 'CUMULATIVE_DISTRIBUTION: Error, x, p and c should be same size'

      norm = integrate_table(x, p, 1, n, iorder)

      DO i = 1, n
         c(i) = integrate_table(x, p, 1, i, iorder)
      END DO

      c = c/norm

   END SUBROUTINE cumulative_distribution

   SUBROUTINE calculate_confidence(x, p, one_sigma, two_sigma)

      USE interpolate
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: p(:)
      REAL, INTENT(OUT) :: one_sigma(2)
      REAL, INTENT(OUT) :: two_sigma(2)
      REAL :: ans, ci
      REAL, ALLOCATABLE :: c(:)
      INTEGER :: i, n
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      n = size(x)
      IF (n /= size(p)) STOP 'CALCULATE_CONFIDENCE: Error, x and p should be same size'
      ALLOCATE(c(n))

      CALL cumulative_distribution(x, p, c)

      DO i = 1, 4

         IF (i == 1) THEN
            ci = 0.5*(1.-erf(2./sqrt(2.))) ! ~0.025
         ELSE IF (i == 2) THEN
            ci = 1.-ci ! ~0.975
         ELSE IF(i == 3) THEN
            ci = 0.5*(1.-erf(1./sqrt(2.))) ! ~0.16
         ELSE IF(i == 4) THEN
            ci = 1.-ci ! ~0.84
         ELSE
            STOP 'CALCULATE_CONFIDENCE: Error, something went wrong'
         END IF

         ans = find(ci, c, x, n, iorder, ifind, imeth)

         IF(i==1) THEN
            two_sigma(1) = ans
         ELSE IF(i==2) THEN
            two_sigma(2) = ans
         ELSE IF(i==3) THEN
            one_sigma(1) = ans
         ELSE IF(i==4) THEN
            one_sigma(2) = ans
         ELSE
            STOP 'CALCULATE_CONFIDENCE: Error, something went wrong'
         END IF

      END DO

   END SUBROUTINE calculate_confidence

   REAL FUNCTION parameter_probability(p, model, x_data, y_data, sigma_data)

      ! Calculate the probability of parameter p in model fitting the data assuming data uncorrelated
      REAL, INTENT(IN) :: p             ! Parameter value
      REAL, EXTERNAL :: model           ! Model function
      REAL, INTENT(IN) :: x_data(:)     ! x values for data
      REAL, INTENT(IN) :: y_data(:)     ! y values for data
      REAL, INTENT(IN) :: sigma_data(:) ! Error bar for data
      REAL :: y, chi2
      INTEGER :: i, n

      INTERFACE
         FUNCTION model(x, pin)
            REAL, INTENT(IN) :: x
            REAL, INTENT(IN) :: pin
         END FUNCTION model
      END INTERFACE

      n = size(x_data)
      IF (n /= size(y_data) .OR. n /= size(sigma_data)) THEN
         STOP 'PARAMETER_PROBABILITY: Error, x_data, y_data and sigma_data should all be the same size'
      END IF

      chi2=0.
      DO i=1,n
         y = model(x_data(i), p)
         chi2 = chi2+((y-y_data(i))/sigma_data(i))**2
      END DO

      parameter_probability=exp(-chi2/2.)

   END FUNCTION parameter_probability

END MODULE statistics
