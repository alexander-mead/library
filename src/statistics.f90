MODULE statistics

   USE table_integer
   USE array_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: mean
   PUBLIC :: variance
   PUBLIC :: histogram

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
            j = select_table_integer(data(i), x, n, 1)
            hist(j) = hist(j)+1
         END IF
      END DO

      WRITE (*, *) 'HISTOGRAM: Fraction of data assigned to histogram:', real(sum(hist))/real(m)
      WRITE (*, *) 'HISTOGRAM: Done'
      WRITE (*, *)

   END SUBROUTINE histogram

END MODULE statistics
