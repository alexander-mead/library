PROGRAM special_functions_demo

   USE special_functions
   USE array_operations

   IMPLICIT NONE
   INTEGER :: itest

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
   WRITE (*, *) '13 - Cubic fix polynomials'
   READ (*, *) itest
   WRITE (*, *)

   IF (is_in_array(itest, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])) THEN
      CALL write_function(itest)
   ELSE IF (itest == 12) THEN
      CALL Fibonacci_demo() 
   ELSE IF (itest == 13) THEN
      CALL fix_polynomail_demo()
   ELSE
      STOP 'SPECIAL_FUNCTION_DEMO: Error, test not specified correctly'
   END IF

CONTAINS

   SUBROUTINE write_function(test)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: test
      REAL :: x, f
      REAL :: xmin, xmax
      INTEGER :: i, n, m
      REAL :: x1, x2, pow
      REAL :: mean, sigma

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
            f = Lognormal_distribution(x, mean, sigma)
         ELSE IF (test == 11) THEN
            f = Gamma(x)
         ELSE
            STOP 'SPECIAL_FUNCTION_TEST: Error, function not specified correctly'
         END IF
         WRITE (*, *) x, f
         WRITE (7, *) x, f
      END DO
      CLOSE (7)

   END SUBROUTINE

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

   END SUBROUTINE Fibonacci_demo

   SUBROUTINE fix_polynomail_demo()

      IMPLICIT NONE
      REAL :: x1, x2, x3, x4, y1, y2, y3, y4
      REAL :: x, f, L
      INTEGER :: i, n
      REAL :: a, b, c, d
      REAL :: f1, f2, f3, f4, L1, L2, L3, L4
      INTEGER :: imode

      WRITE(*,*)
      WRITE(*,*) 'Cubic fit tester'
      WRITE(*,*) '=========='
      WRITE(*,*)
      WRITE(*,*) '1 - Line'
      WRITE(*,*) '2 - Slight curve'
      WRITE(*,*) '3 - Curve'
      WRITE(*,*) '4 - Crazy curve'
      WRITE(*,*) '5 - Crazier curve'
      WRITE(*,*) '6 - Craziest curve'
      READ(*,*) imode
      WRITE(*,*)
    
      x1=1.
      x2=2.
      x3=3.
      x4=4.
    
      IF(imode==1) THEN
    
         y1=1.
         y2=2.
         y3=3.
         y4=4.
    
      ELSE IF(imode==2) THEN
    
         y1=1.
         y2=1.8
         y3=3.2
         y4=4.
    
      ELSE IF(imode==3) THEN
    
         y1=1.
         y2=1.5
         y3=3.5
         y4=4.
    
      ELSE IF(imode==4) THEN
    
         y1=1.
         y2=1.
         y3=4.
         y4=4.
    
      ELSE IF(imode==5) THEN
    
         y1=1.
         y2=3.
         y3=2.
         y4=4.
    
      ELSE IF(imode==6) THEN
    
         y1=1.
         y2=1.2
         y3=0.6
         y4=4.
    
      ELSE
    
         STOP 'Curve not specified correctly'
    
      END IF
    
      !CALL fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
      CALL fix_cubic(a, b, c, d, [x1, x2, x3, x4], [y1, y2, y3, y4])
    
      OPEN(7,file='points.dat')
      WRITE(7,*) x1, y1
      WRITE(7,*) x2, y2
      WRITE(7,*) x3, y3
      WRITE(7,*) x4, y4
      CLOSE(7)
    
      n=100
      OPEN(7,file='cubic.dat')
      DO i=1,n
         x=x1+(x4-x1)*float(i-1)/float(n-1)
         f=a*(x**3.)+b*(x**2.)+c*x+d
         L=Lagrange_polynomial(x,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
         WRITE(7,*) x, f, L
      END DO
      CLOSE(7)
    
      f1=a*(x1**3.)+b*(x1**2.)+c*x1+d
      f2=a*(x2**3.)+b*(x2**2.)+c*x2+d
      f3=a*(x3**3.)+b*(x3**2.)+c*x3+d
      f4=a*(x4**3.)+b*(x4**2.)+c*x4+d
    
      L1=Lagrange_polynomial(x1,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
      L2=Lagrange_polynomial(x2,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
      L3=Lagrange_polynomial(x3,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
      L4=Lagrange_polynomial(x4,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
    
      WRITE(*,*) '======================================='
      WRITE(*,*) '      x_i       y_i      y(x)      L(x)'
      WRITE(*,*) '======================================='
      WRITE(*,fmt='(4F10.5)') x1, y1, f1, L1
      WRITE(*,fmt='(4F10.5)') x2, y2, f2, L2
      WRITE(*,fmt='(4F10.5)') x3, y3, f3, L3
      WRITE(*,fmt='(4F10.5)') x4, y4, f4, L4
      WRITE(*,*) '======================================='

   END SUBROUTINE

END PROGRAM special_functions_demo
