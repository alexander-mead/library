MODULE statistics

   USE table_integer
   USE array_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: mean
   PUBLIC :: variance
   PUBLIC :: standard_deviation
   PUBLIC :: histogram
   PUBLIC :: calculate_confidence
   PUBLIC :: parameter_probability

CONTAINS

   REAL FUNCTION mean(x, n)

      IMPLICIT NONE
      REAL, INTENT(IN) :: x(n)
      INTEGER, INTENT(IN) :: n
      DOUBLE PRECISION :: sum
      INTEGER :: i

      sum = 0.d0
      DO i = 1, n
         sum = sum+x(i)
      END DO

      mean = real(sum)/real(n)

   END FUNCTION mean

   REAL FUNCTION variance(x, n)

      ! Note that depending on the application one might want to use n-1 here
      ! Difference is between sample variance and population varaince
      ! See https://en.wikipedia.org/wiki/Bessel%27s_correction
      ! Sample is simply the variance of your sample 1/N
      ! However, if you are estimating the population variance then you need 1/(N-1)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x(n)
      INTEGER, INTENT(IN) :: n
      DOUBLE PRECISION :: sum
      REAL :: avg
      INTEGER :: i

      avg = mean(x, n)

      sum = 0.d0
      DO i = 1, n
         sum = sum+(x(i)-avg)**2
      END DO

      variance = real(sum)/real(n)

   END FUNCTION variance

   REAL FUNCTION standard_deviation(x, n)

      IMPLICIT NONE
      REAL, INTENT(IN) :: x(n)
      INTEGER, INTENT(IN) :: n

      standard_deviation = sqrt(variance(x, n))

   END FUNCTION standard_deviation

   SUBROUTINE histogram(xmin, xmax, x, hist, n, data, m)

      IMPLICIT NONE
      REAL, INTENT(IN) :: xmin     ! Minimum x value
      REAL, INTENT(IN) :: xmax     ! Maximum x value
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)       ! Output array of bin edges, size n+1
      INTEGER, ALLOCATABLE, INTENT(OUT) :: hist(:) ! Output integer array of bin counts, size n
      INTEGER, INTENT(IN) :: n     ! Number of bins
      REAL, INTENT(IN) :: data(m)  ! Data to be binned
      INTEGER, INTENT(IN) :: m     ! Number of points to be binned
      INTEGER :: i, j

      WRITE (*, *) 'HISTOGRAM: Assiging arrays'

      !Fill the table for the xrange and allocate the histogram array
      CALL fill_array(xmin, xmax, x, n+1)

      !Set the histogram to zero
      IF (ALLOCATED(hist)) DEALLOCATE (hist)
      ALLOCATE (hist(n))
      hist = 0

      WRITE (*, *) 'HISTOGRAM: Constructing histogram'

      !Make the histogram from the data
      DO i = 1, m
         IF (data(i) < xmin .OR. data(i) > xmax) THEN
            CYCLE
         ELSE
            ! TODO: Check ifind_linear is correct here
            j = find_table_integer(data(i), x, n, ifind_linear)
            hist(j) = hist(j)+1
         END IF
      END DO

      WRITE (*, *) 'HISTOGRAM: Fraction of data assigned to histogram:', real(sum(hist))/real(m)
      WRITE (*, *) 'HISTOGRAM: Done'
      WRITE (*, *)

   END SUBROUTINE histogram

   SUBROUTINE cumulative_distribution(x, p, c, n)

      USE calculus_table
      IMPLICIT NONE
      REAL, INTENT(IN) :: x(n)
      REAL, INTENT(IN) :: p(n)
      REAL, INTENT(OUT) :: c(n)
      INTEGER, INTENT(IN) :: n
      REAL :: norm
      INTEGER :: i
      INTEGER, PARAMETER :: iorder=0 ! Zeroth-order (histogram) integration and plenty of points are best here

      norm=integrate_table(x, p, n, 1, n, iorder)

      DO i=1,n
         c(i)=integrate_table(x, p, n, 1, i, iorder)
      END DO

      c=c/norm

   END SUBROUTINE cumulative_distribution

   SUBROUTINE calculate_confidence(x, p, n, one_sigma, two_sigma)

      USE interpolate
      IMPLICIT NONE
      REAL, INTENT(IN) :: x(n)
      REAL, INTENT(IN) :: p(n)
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(OUT) :: one_sigma(2)
      REAL, INTENT(OUT) :: two_sigma(2)
      REAL :: c(n), ans, ci
      INTEGER :: i
      INTEGER, PARAMETER :: iorder=3
      INTEGER, PARAMETER :: ifind=3
      INTEGER, PARAMETER :: imeth=2

      CALL cumulative_distribution(x, p, c, n)

      DO i=1,4

         IF(i==1) THEN
            ci = 0.5*(1.-erf(2./sqrt(2.))) ! ~0.025
         ELSE IF(i==2) THEN
            ci = 1.-ci ! ~0.975
         ELSE IF(i==3) THEN
            ci = 0.5*(1.-erf(1./sqrt(2.))) ! ~0.16
         ELSE IF(i==4) THEN
            ci = 1.-ci ! ~0.84
         ELSE
            STOP 'CALCULATE_CONFIDENCE: Error, something went wrong'
         END IF

         ans = find(ci, c, x, n, iorder, ifind, imeth)

         IF(i==1) THEN
            two_sigma(1)=ans
         ELSE IF(i==2) THEN
            two_sigma(2)=ans
         ELSE IF(i==3) THEN
            one_sigma(1) = ans
         ELSE IF(i==4) THEN
            one_sigma(2) = ans
         ELSE
            STOP 'CALCULATE_CONFIDENCE: Error, something went wrong'
         END IF

      END DO

   END SUBROUTINE calculate_confidence

   REAL FUNCTION parameter_probability(p, model, x_data, y_data, sigma_data, n)

      ! Calculate the probability of parameter p in model fitting the data assuming data uncorrelated
      IMPLICIT NONE
      REAL, INTENT(IN) :: p             ! Parameter value
      REAL, EXTERNAL :: model           ! Model function
      REAL, INTENT(IN) :: x_data(n)     ! x values for data
      REAL, INTENT(IN) :: y_data(n)     ! y values for data
      REAL, INTENT(IN) :: sigma_data(n) ! Error bar for data
      INTEGER, INTENT(IN) :: n          ! Number of data points
      REAL :: y, chi2
      INTEGER :: i

      INTERFACE
         FUNCTION model(x, p)
            REAL, INTENT(IN) :: x
            REAL, INTENT(IN) :: p
         END FUNCTION model
      END INTERFACE

      chi2=0.
      DO i=1,n
         y = model(x_data(i), p)
         chi2 = chi2+((y-y_data(i))/sigma_data(i))**2
      END DO

      parameter_probability=exp(-chi2/2.)

   END FUNCTION parameter_probability

END MODULE statistics
