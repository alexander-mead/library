MODULE calculus_table

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: derivative_table
   PUBLIC :: integrate_table

   INTERFACE integrate_table
      MODULE PROCEDURE integrate_table_1D
      MODULE PROCEDURE integrate_table_2D
   END INTERFACE integrate_table

CONTAINS

   REAL FUNCTION derivative_table(x, xin, yin, iorder, ifind)

      ! Given two arrays x and y such that y=y(x) this uses interpolation to calculate the derivative y'(x_i) at position x_i
      USE table_integer
      USE special_functions
      USE array_operations
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xin(:)
      REAL, INTENT(IN) :: yin(:)
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      REAL, ALLOCATABLE ::  xtab(:), ytab(:)
      REAL :: a, b, c, d
      REAL :: x1, x2, x3, x4
      REAL :: y1, y2, y3, y4
      INTEGER :: i, n

      ! This version interpolates if the value is off either end of the array!
      ! Care should be chosen to insert x, xtab, ytab as log if this might give better!
      ! Results from the interpolation!

      ! ifind = 1 => find x in xtab by crudely searching
      ! ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
      ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

      ! iorder = 1 => linear interpolation
      ! iorder = 2 => quadratic interpolation
      ! iorder = 3 => cubic interpolation

      n = size(xin)
      IF (n /= size(yin)) STOP 'DERIVATIVE_TABLE: Input x and y arrays should be of the same size'
      ALLOCATE(xtab(n), ytab(n))

      xtab = xin
      ytab = yin

      IF (xtab(1) > xtab(n)) THEN
         ! Reverse the arrays in this case
         CALL reverse_array(xtab)
         CALL reverse_array(ytab)
      END IF

      IF (iorder == 1) THEN

         IF (n < 2) STOP 'DERIVATIVE_TABLE: Not enough points in your table for linear interpolation'

         IF (x <= xtab(2)) THEN

            x2 = xtab(2)
            x1 = xtab(1)

            y2 = ytab(2)
            y1 = ytab(1)

         ELSE IF (x >= xtab(n-1)) THEN

            x2 = xtab(n)
            x1 = xtab(n-1)

            y2 = ytab(n)
            y1 = ytab(n-1)

         ELSE

            i = find_table_integer(x, xtab, ifind)

            x2 = xtab(i+1)
            x1 = xtab(i)

            y2 = ytab(i+1)
            y1 = ytab(i)

         END IF

         CALL fix_polynomial(a, b, [x1, x2], [y1, y2])
         derivative_table = a

      ELSE IF (iorder == 2) THEN

         IF (n < 3) STOP 'DERIVATIVE_TABLE: Not enough points in your table'

         IF (x <= xtab(2) .OR. x >= xtab(n-1)) THEN

            IF (x <= xtab(2)) THEN

               x3 = xtab(3)
               x2 = xtab(2)
               x1 = xtab(1)

               y3 = ytab(3)
               y2 = ytab(2)
               y1 = ytab(1)

            ELSE IF (x >= xtab(n-1)) THEN

               x3 = xtab(n)
               x2 = xtab(n-1)
               x1 = xtab(n-2)

               y3 = ytab(n)
               y2 = ytab(n-1)
               y1 = ytab(n-2)

            ELSE

               STOP 'DERIVATIVE_TABLE: Error, something went wrong'

            END IF

            CALL fix_polynomial(a, b, c, [x1, x2, x3], [y1, y2, y3])

            derivative_table = 2.*a*x+b

         ELSE

            i = find_table_integer(x, xtab, ifind)

            x1 = xtab(i-1)
            x2 = xtab(i)
            x3 = xtab(i+1)
            x4 = xtab(i+2)

            y1 = ytab(i-1)
            y2 = ytab(i)
            y3 = ytab(i+1)
            y4 = ytab(i+2)

            ! In this case take the average of two separate quadratic spline values

            derivative_table = 0.

            CALL fix_polynomial(a, b, c, [x1, x2, x3], [y1, y2, y3])
            derivative_table = derivative_table+(2.*a*x+b)/2.

            CALL fix_polynomial(a, b, c, [x2, x3, x4], [y2, y3, y4])
            derivative_table = derivative_table+(2.*a*x+b)/2.

         END IF

      ELSE IF (iorder == 3) THEN

         IF (n < 4) STOP 'DERIVATIVE_TABLE: Not enough points in your table'

         IF (x <= xtab(3)) THEN

            x4 = xtab(4)
            x3 = xtab(3)
            x2 = xtab(2)
            x1 = xtab(1)

            y4 = ytab(4)
            y3 = ytab(3)
            y2 = ytab(2)
            y1 = ytab(1)

         ELSE IF (x >= xtab(n-2)) THEN

            x4 = xtab(n)
            x3 = xtab(n-1)
            x2 = xtab(n-2)
            x1 = xtab(n-3)

            y4 = ytab(n)
            y3 = ytab(n-1)
            y2 = ytab(n-2)
            y1 = ytab(n-3)

         ELSE

            i = find_table_integer(x, xtab, ifind)

            x1 = xtab(i-1)
            x2 = xtab(i)
            x3 = xtab(i+1)
            x4 = xtab(i+2)

            y1 = ytab(i-1)
            y2 = ytab(i)
            y3 = ytab(i+1)
            y4 = ytab(i+2)

         END IF

         CALL fix_polynomial(a, b, c, d, [x1, x2, x3, x4], [y1, y2, y3, y4])
         derivative_table = 3.*a*(x**2.)+2.*b*x+c

      ELSE

         STOP 'DERIVATIVE_TABLE: Error, order not specified correctly'

      END IF

   END FUNCTION derivative_table

   REAL FUNCTION integrate_table_1D(x, y, n1, n2, iorder)
      
      ! Integrates tables y(x)dx
      USE special_functions
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: y(:)
      INTEGER, INTENT(IN) :: n1
      INTEGER, INTENT(IN) :: n2
      INTEGER, INTENT(IN) :: iorder
      REAL :: a, b, c, d, h
      REAL :: q1, q2, q3, qi, qf
      REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
      REAL :: sum
      INTEGER :: i, i1, i2, i3, i4, n

      n = size(x)
      IF (n /= size(y)) STOP 'INTEGRATE_TABLE_1D: Error, x and y should be the same size'
      IF (n1 < 1) STOP 'INTEGRATE_TABLE_1D: Error, n1 cannot be less than 1'
      IF (n2 > n) STOP 'INTEGRATE_TABLE_1D: Error, n2 cannot be greater than n'
      IF (n2 < n1) STOP 'INTEGRATE_TABLE_1D: Error, n2 cannot be less than n1'

      sum = 0.d0

      ! I think if n1=n2 then the result will just be zero anyway
      !IF(n2<=n1) STOP 'INTEGRATE_TABLE: Error n2 must be greater than n1'

      IF (n1 == n2) THEN

         integrate_table_1D = 0.

      ELSE IF(iorder==0) THEN

         ! Data must be evenly spaced for zeroth order to work
         ! x coordinates must represent bin centres
         ! Note from the range x(1)->x(2) there are only n-1 bins
         ! The 1st bin spills out by h/2 below x(1) and the nth by h/2 above x(n)
         ! This is only useful for summing histograms etc.
         h=(x(n)-x(1))/real(n-1)

         ! Each rectangle contributes y(i)*dx to the integral
         DO i = n1, n2
            sum=sum+y(i)
         END DO
         sum=sum*h

      ELSE IF (iorder == 1) THEN

         ! Sums over all Trapezia (a+b)*h/2
         DO i = n1, n2-1
            a = y(i+1)
            b = y(i)
            h = x(i+1)-x(i)
            sum = sum+(a+b)*h
         END DO
         sum=sum/2.

      ELSE IF (iorder == 2) THEN

         DO i = n1, n2-2

            x1 = x(i)
            x2 = x(i+1)
            x3 = x(i+2)

            y1 = y(i)
            y2 = y(i+1)
            y3 = y(i+2)

            CALL fix_polynomial(a, b, c, [x1, x2, x3], [y1, y2, y3])

            q1 = a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
            q2 = a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
            q3 = a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3

            ! Takes value for first and last sections but averages over sections where you
            ! have two independent estimates of the area
            IF (n == 3) THEN
               sum = sum+q3-q1
            ELSE IF (i == 1) THEN
               sum = sum+(q2-q1)+(q3-q2)/2.d0
            ELSE IF (i == n-2) THEN
               sum = sum+(q2-q1)/2.d0+(q3-q2)
            ELSE
               sum = sum+(q3-q1)/2.
            END IF

         END DO

      ELSE IF (iorder == 3) THEN

         DO i = n1, n2-1

            ! First choose the integers used for defining cubics for each section
            ! First and last are different because the section does not lie in the *middle* of a cubic

            IF (i == 1) THEN

               i1 = 1
               i2 = 2
               i3 = 3
               i4 = 4

            ELSE IF (i == n-1) THEN

               i1 = n-3
               i2 = n-2
               i3 = n-1
               i4 = n

            ELSE

               i1 = i-1
               i2 = i
               i3 = i+1
               i4 = i+2

            END IF

            x1 = x(i1)
            x2 = x(i2)
            x3 = x(i3)
            x4 = x(i4)

            y1 = y(i1)
            y2 = y(i2)
            y3 = y(i3)
            y4 = y(i4)

            CALL fix_polynomial(a, b, c, d, [x1, x2, x3, x4], [y1, y2, y3, y4])

            ! These are the limits of the particular section of integral
            xi = x(i)
            xf = x(i+1)

            qi = a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
            qf = a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf

            sum = sum+qf-qi

         END DO

      ELSE

         STOP 'INTEGRATE_TABLE_1D: Error, order not specified correctly'

      END IF

      integrate_table_1D = sum

   END FUNCTION integrate_table_1D

   REAL FUNCTION integrate_table_2D(x, y, F)

      ! A crude integration scheme for tabulated 2D functions
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: y(:)
      REAL, INTENT(IN) :: F(:, :)
      INTEGER :: ix, iy, nx, ny
      REAL :: sum
      REAL :: dx, dy

      nx = size(x)
      ny = size(y)
      IF (nx /= size(F, 1)) STOP 'INTEGRATE_TABLE_2D: Size of x must be the same size as the first dimension of F'
      IF (ny /= size(F, 2)) STOP 'INTEGRATE_TABLE_2D: Size of y must be the same size as the second dimension of F'

      ! Set the sum variable to zero
      sum = 0.

      ! Loop over function, calculate dx, dy, sum to get integral
      DO iy = 1, ny-1

         dy = y(iy+1)-y(iy)

         DO ix = 1, nx-1

            dx = x(ix+1)-x(ix)

            sum = sum+dx*dy*(F(ix, iy)+F(ix+1, iy)+F(ix, iy+1)+F(ix+1, iy+1))/4.

         END DO

      END DO

      integrate_table_2D = real(sum)

   END FUNCTION integrate_table_2D

END MODULE calculus_table
