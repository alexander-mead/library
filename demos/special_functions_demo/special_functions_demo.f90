PROGRAM special_functions_demo

   USE special_functions

   IMPLICIT NONE
   REAL :: x, f
   REAL :: xmin, xmax
   INTEGER :: i, n, m
   REAL :: x1, x2, pow
   REAL :: mean, sigma
   INTEGER :: test

   WRITE (*, *) 'Choose a special function to test'
   WRITE (*, *) ' 1 - Apodise'
   WRITE (*, *) ' 2 - Blob'
   WRITE (*, *) ' 3 - Smooth apodise'
   WRITE (*, *) ' 4 - Smooth blob'
   WRITE (*, *) ' 5 - Sinc'
   WRITE (*, *) ' 6 - Top hat Fourier transform'
   WRITE (*, *) ' 7 - J0'
   WRITE (*, *) ' 8 - J1'
   WRITE (*, *) ' 9 - J2'
   WRITE (*, *) '10 - Log normal'
   WRITE (*, *) '11 - Gamma'
   WRITE (*, *) '12 - Fibonacci'
   READ (*, *) test
   WRITE (*, *)

   IF (test <= 11) THEN

      xmin = 0.
      xmax = 1e1
      n = 501

      IF (test == 1 .OR. test == 3) THEN
         !Parameters for the apodise functions
         x1 = 4.
         x2 = 6.
         pow = 0.5
      ELSE IF (test == 2 .OR. test == 4) THEN
         !Parameters for the blob function
         x1 = 2.
         x2 = 8.
         pow = 0.1
      ELSE IF (test == 10) THEN
         !Parameters for the log-normal function
         xmin = 1e-3
         xmax = 1e1
         mean = 2.
         sigma = 0.5
      END IF

      OPEN (7, file='results.dat')
      DO i = 1, n
         x = xmin+(xmax-xmin)*float(i-1)/float(n-1)
         IF (test == 1) THEN
            f = apodise(x, x1, x2, pow)
         ELSE IF (test == 2) THEN
            f = blob(x, x1, x2, pow)
         ELSE IF (test == 3) THEN
            f = smooth_apodise(x, x1, x2, pow)
         ELSE IF (test == 4) THEN
            f = smooth_blob(x, x1, x2, pow)
         ELSE IF (test == 5) THEN
            f = sinc(x)
         ELSE IF (test == 6) THEN
            f = wk_tophat(x)
         ELSE IF (test == 7 .OR. test == 8 .OR. test == 9) THEN
            IF (test == 7) m = 0
            IF (test == 8) m = 1
            IF (test == 9) m = 2
            f = Bessel(m, x)
         ELSE IF (test == 10) THEN
            f = Lognormal(x, mean, sigma)
         ELSE IF (test == 11) THEN
            f = Gamma(x)
         ELSE
            STOP 'SPECIAL_FUNCTION_TEST: Error, function not specified correctly'
         END IF
         WRITE (*, *) x, f
         WRITE (7, *) x, f
      END DO
      CLOSE (7)

   ELSE IF (test == 12) THEN

      CALL Fibonacci_demo()

   ELSE

      STOP 'SPECIAL_FUNCTION_TEST: Error, test not specified correctly'

   END IF

CONTAINS

   SUBROUTINE Fibonacci_demo()

      IMPLICIT NONE
      INTEGER :: i
      INTEGER, ALLOCATABLE :: F(:)
      INTEGER, PARAMETER :: n = 10

      ALLOCATE (F(n))
      CALL get_Fibonaccis(F, n)
      DO i = 1, n
         WRITE (*, *) i, F(i), Fibonacci(i)
      END DO
      WRITE (*, *)

   END SUBROUTINE

END PROGRAM special_functions_demo
