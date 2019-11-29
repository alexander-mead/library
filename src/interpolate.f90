MODULE interpolate

   USE fix_polynomial
   USE table_integer
   USE array_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: find
   PUBLIC :: interpolate_array
   PUBLIC :: iinterp_polynomial
   PUBLIC :: iinterp_Lagrange

   INTEGER, PARAMETER :: iinterp_polynomial = 1
   INTEGER, PARAMETER :: iinterp_Lagrange = 2

   INTERFACE find
      MODULE PROCEDURE find_1D
      MODULE PROCEDURE find_2D
      MODULE PROCEDURE find_3D
   END INTERFACE find

CONTAINS

   REAL FUNCTION find_1D(x, xin, yin, n, iorder, ifind, iinterp)

      ! Given two arrays x and y this routine interpolates to find the y_i value at position x_i
      ! Care should be chosen to insert x, xtab, ytab as log if this might give beter results
      ! Extrapolates if the value is off either end of the array   
      ! If the value required is off the table edge the extrapolation is always linear
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xin(n)
      REAL, INTENT(IN) :: yin(n)
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      REAL ::  xtab(n), ytab(n)
      REAL :: a, b, c, d
      REAL :: x1, x2, x3, x4
      REAL :: y1, y2, y3, y4
      REAL :: L1, L2
      REAL :: dx
      INTEGER :: i

      ! iorder = 0 => constant interpolation (histogram data)
      ! iorder = 1 => linear interpolation
      ! iorder = 2 => quadratic interpolation
      ! iorder = 3 => cubic interpolation

      ! ifind = 1 => find x in xtab quickly assuming the table is linearly spaced
      ! ifind = 2 => find x in xtab by crudely searching from x(1) to x(n)
      ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

      ! iinterp = 1 => Uses standard polynomials for interpolation
      ! iinterp = 2 => Uses Lagrange polynomials for interpolation

      xtab = xin
      ytab = yin

      IF (xtab(1) > xtab(n)) THEN
         ! Reverse the arrays in this case
         CALL reverse_array(xtab, n)
         CALL reverse_array(ytab, n)
      END IF

      IF(iorder == 0) THEN

         dx = xtab(2)-xtab(1) ! Assumes bins are equally spaced
         xtab = xtab - dx/2.  ! Assumes that x coordinates give left edge of histogram
         
         IF(x < xtab(1)) THEN
            ! Outside lower boundary
            find_1D = 0.
         ELSE IF(x >= xtab(n)+dx) THEN
            ! Outside upper boundary
            find_1D = 0.
         ELSE
            i = find_table_integer(x, xtab, n, ifind)
            find_1D = ytab(i)
         END IF

      ELSE

         IF (x < xtab(1)) THEN

            ! Do a linear extrapolation beyond the table boundary

            x1 = xtab(1)
            x2 = xtab(2)

            y1 = ytab(1)
            y2 = ytab(2)

            IF (iinterp == iinterp_polynomial) THEN
               CALL fix_line(a, b, x1, y1, x2, y2)
               find_1D = a*x+b
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 1, (/x1, x2/), (/y1, y2/))
            ELSE
               STOP 'FIND_1D: Error, method not specified correctly'
            END IF

         ELSE IF (x > xtab(n)) THEN

            ! Do a linear extrapolation beyond the table boundary

            x1 = xtab(n-1)
            x2 = xtab(n)

            y1 = ytab(n-1)
            y2 = ytab(n)

            IF (iinterp == iinterp_polynomial) THEN
               CALL fix_line(a, b, x1, y1, x2, y2)
               find_1D = a*x+b
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 1, (/x1, x2/), (/y1, y2/))
            ELSE
               STOP 'FIND_1D: Error, method not specified correctly'
            END IF

         ELSE IF (iorder == 1) THEN

            IF (n < 2) STOP 'FIND_1D: Not enough points in your table for linear interpolation'

            IF (x <= xtab(2)) THEN

               x1 = xtab(1)
               x2 = xtab(2)

               y1 = ytab(1)
               y2 = ytab(2)

            ELSE IF (x >= xtab(n-1)) THEN

               x1 = xtab(n-1)
               x2 = xtab(n)

               y1 = ytab(n-1)
               y2 = ytab(n)

            ELSE

               i = find_table_integer(x, xtab, n, ifind)

               x1 = xtab(i)
               x2 = xtab(i+1)

               y1 = ytab(i)
               y2 = ytab(i+1)

            END IF

            IF (iinterp == iinterp_polynomial) THEN
               CALL fix_line(a, b, x1, y1, x2, y2)
               find_1D = a*x+b
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 1, (/x1, x2/), (/y1, y2/))
            ELSE
               STOP 'FIND_1D: Error, method not specified correctly'
            END IF

         ELSE IF (iorder == 2) THEN

            IF (n < 3) STOP 'FIND_1D: Not enough points in your table'

            IF (x <= xtab(2) .OR. x >= xtab(n-1)) THEN

               IF (x <= xtab(2)) THEN

                  x1 = xtab(1)
                  x2 = xtab(2)
                  x3 = xtab(3)

                  y1 = ytab(1)
                  y2 = ytab(2)
                  y3 = ytab(3)

               ELSE IF (x >= xtab(n-1)) THEN

                  x1 = xtab(n-2)
                  x2 = xtab(n-1)
                  x3 = xtab(n)

                  y1 = ytab(n-2)
                  y2 = ytab(n-1)
                  y3 = ytab(n)

               END IF

               IF (iinterp == iinterp_polynomial) THEN
                  CALL fix_quadratic(a, b, c, x1, y1, x2, y2, x3, y3)
                  find_1D = a*(x**2)+b*x+c
               ELSE IF (iinterp == iinterp_Lagrange) THEN
                  find_1D = Lagrange_polynomial(x, 2, (/x1, x2, x3/), (/y1, y2, y3/))
               ELSE
                  STOP 'FIND_1D: Error, method not specified correctly'
               END IF

            ELSE

               i = find_table_integer(x, xtab, n, ifind)

               x1 = xtab(i-1)
               x2 = xtab(i)
               x3 = xtab(i+1)
               x4 = xtab(i+2)

               y1 = ytab(i-1)
               y2 = ytab(i)
               y3 = ytab(i+1)
               y4 = ytab(i+2)

               IF (iinterp == iinterp_polynomial) THEN
                  ! In this case take the average of two separate quadratic spline values
                  CALL fix_quadratic(a, b, c, x1, y1, x2, y2, x3, y3)
                  find_1D = (a*x**2+b*x+c)/2.
                  CALL fix_quadratic(a, b, c, x2, y2, x3, y3, x4, y4)
                  find_1D = find_1D+(a*x**2+b*x+c)/2.
               ELSE IF (iinterp == iinterp_Lagrange) THEN
                  ! In this case take the average of two quadratic Lagrange polynomials
                  L1 = Lagrange_polynomial(x, 2, (/x1, x2, x3/), (/y1, y2, y3/))
                  L2 = Lagrange_polynomial(x, 2, (/x2, x3, x4/), (/y2, y3, y4/))
                  find_1D = (L1+L2)/2.
               ELSE
                  STOP 'FIND_1D: Error, method not specified correctly'
               END IF

            END IF

         ELSE IF (iorder == 3) THEN

            IF (n < 4) STOP 'FIND_1D: Not enough points in your table'

            IF (x <= xtab(3)) THEN

               x1 = xtab(1)
               x2 = xtab(2)
               x3 = xtab(3)
               x4 = xtab(4)

               y1 = ytab(1)
               y2 = ytab(2)
               y3 = ytab(3)
               y4 = ytab(4)

            ELSE IF (x >= xtab(n-2)) THEN

               x1 = xtab(n-3)
               x2 = xtab(n-2)
               x3 = xtab(n-1)
               x4 = xtab(n)

               y1 = ytab(n-3)
               y2 = ytab(n-2)
               y3 = ytab(n-1)
               y4 = ytab(n)

            ELSE

               i = find_table_integer(x, xtab, n, ifind)

               x1 = xtab(i-1)
               x2 = xtab(i)
               x3 = xtab(i+1)
               x4 = xtab(i+2)

               y1 = ytab(i-1)
               y2 = ytab(i)
               y3 = ytab(i+1)
               y4 = ytab(i+2)

            END IF

            IF (iinterp == iinterp_polynomial) THEN
               CALL fix_cubic(a, b, c, d, x1, y1, x2, y2, x3, y3, x4, y4)
               find_1D = a*x**3+b*x**2+c*x+d
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 3, (/x1, x2, x3, x4/), (/y1, y2, y3, y4/))
            ELSE
               STOP 'FIND_1D: Error, method not specified correctly'
            END IF

         ELSE

            STOP 'FIND_1D: Error, interpolation order specified incorrectly'

         END IF

      END IF

   END FUNCTION find_1D

   REAL FUNCTION find_2D(x, xin, y, yin, fin, nx, ny, iorder, ifind, iinterp)

      ! A 2D interpolation routine to find value f(x,y) at position x, y
      ! Care should be chosen to insert x, xtab, ytab as log if this might give better results
      ! Extrapolates if the value is off either end of the array
      ! If the value required is off the table edge the extrapolation is always linear
      ! TODO: Loops over coordinates to avoid repetition?
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xin(nx)
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: yin(ny)
      REAL, INTENT(IN) :: fin(nx, ny)
      INTEGER, INTENT(IN) :: nx
      INTEGER, INTENT(IN) :: ny
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      REAL ::  xtab(nx), ytab(ny), ftab(nx, ny)
      REAL :: a, b, c, d
      REAL :: x1, x2, x3, x4
      REAL :: y1, y2, y3, y4
      REAL :: f11, f12, f13, f14
      REAL :: f21, f22, f23, f24
      REAL :: f31, f32, f33, f34
      REAL :: f41, f42, f43, f44
      REAL :: f10, f20, f30, f40
      REAL :: f01, f02, f03, f04
      INTEGER :: i1, i2, i3, i4
      INTEGER :: j1, j2, j3, j4
      REAL :: findx, findy, V(2,2)
      INTEGER :: i, j, ix, iy, ix1, ix2, iy1, iy2

      ! iorder = 1 => linear interpolation
      ! iorder = 2 => quadratic interpolation
      ! iorder = 3 => cubic interpolation

      ! ifind = 1 => find x in xtab by crudely searching from x(1) to x(n)
      ! ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
      ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

      ! iinterp = 1 => Uses cubic polynomials for interpolation
      ! iinterp = 2 => Uses Lagrange polynomials for interpolation
      IF (iinterp == iinterp_Lagrange) STOP 'FIND_2D: No Lagrange polynomials for you'

      xtab = xin
      ytab = yin
      ftab = fin

      IF (xtab(1) > xtab(nx)) STOP 'FIND_2D: x array in wrong order'
      IF (ytab(1) > ytab(ny)) STOP 'FIND_2D: y array in wrong order'

      IF(iorder == 0) THEN

         ix = find_table_integer(x, xin, nx, ifind)
         IF(ix == 0) ix = 1
         iy = find_table_integer(y, yin, ny, ifind)
         IF(iy == 0) iy = 1

         find_2D = fin(ix, iy)

      ELSE IF(iorder == 1) THEN

         IF (nx < 2) STOP 'FIND_2D: Not enough x points in your array for linear interpolation'
         IF (ny < 2) STOP 'FIND_2D: Not enough y points in your array for linear interpolation'

         !! Get the x,y values !!

         ! Get the integer coordinates in the x direction
         ix = find_table_integer(x, xin, nx, ifind)
         IF(ix==0) THEN
            ix=1
         ELSE IF(ix==nx) THEN
            ix=nx-1
         END IF
         ix1 = ix
         ix2 = ix+1

         ! Get the x values at the corners
         x1 = xin(ix1)
         x2 = xin(ix2)

         ! Get the integer coordinates in the y direction
         iy = find_table_integer(y, yin, ny, ifind)
         IF(iy==0) THEN
            iy=1
         ELSE IF(iy==ny) THEN
            iy=ny-1
         END IF
         iy1 = iy
         iy2 = iy+1

         ! Get the y values at the corners
         y1 = yin(iy1)
         y2 = yin(iy2)

         ! Interpolation is function values at corners weighted by opposite area
         ! TODO: Clever loop here? Maybe not, its nice and transparent without the loop
         V(1, 1) = (x2-x)*(y2-y)*fin(ix1, iy1)
         V(1, 2) = (x2-x)*(y-y1)*fin(ix1, iy2)
         V(2, 1) = (x-x1)*(y2-y)*fin(ix2, iy1)
         V(2, 2) = (x-x1)*(y-y1)*fin(ix2, iy2)

         ! Normalisation
         find_2D = sum(V)/((x2-x1)*(y2-y1))

      ELSE IF (iorder == 2) THEN

         STOP 'FIND_2D: Quadratic 2D interpolation not implemented - also probably pointless'

      ELSE IF (iorder == 3) THEN

         ! No cubic extrapolation implemented if the desired point is outside x AND y array boundary corners
         IF ((x < xtab(1) .OR. x > xtab(nx)) .AND. (y > ytab(ny) .OR. y < ytab(1))) THEN
            WRITE (*, *) 'FIND_2D: array xmin:', xtab(1)
            WRITE (*, *) 'FIND_2D: array xmax:', xtab(nx)
            WRITE (*, *) 'FIND_2D: requested x:', x
            WRITE (*, *) 'FIND_2D: array ymin:', ytab(1)
            WRITE (*, *) 'FIND_2D: array ymax:', ytab(ny)
            WRITE (*, *) 'FIND_2D: requested y:', y
            STOP 'FIND_2D: Desired point is outside x AND y array range'
         END IF

         IF (x < xtab(1) .OR. x > xtab(nx)) THEN

            IF (nx < 2) STOP 'FIND_2D: Not enough x points in your array for linear interpolation'
            IF (ny < 4) STOP 'FIND_2D: Not enough y points in your array for cubic interpolation'

            ! x is off the table edge

            IF (x < xtab(1)) THEN
               i1 = 1
               i2 = 2
            ELSE
               i1 = nx-1
               i2 = nx
            END IF

            x1 = xtab(i1)
            x2 = xtab(i2)

            IF (y <= ytab(4)) THEN
               j = 2
            ELSE IF (y >= ytab(ny-3)) THEN
               j = ny-2
            ELSE
               j = find_table_integer(y, ytab, ny, ifind)
            END IF

            j1 = j-1
            j2 = j
            j3 = j+1
            j4 = j+2

            y1 = ytab(j1)
            y2 = ytab(j2)
            y3 = ytab(j3)
            y4 = ytab(j4)

            f11 = ftab(i1, j1)
            f12 = ftab(i1, j2)
            f13 = ftab(i1, j3)
            f14 = ftab(i1, j4)

            f21 = ftab(i2, j1)
            f22 = ftab(i2, j2)
            f23 = ftab(i2, j3)
            f24 = ftab(i2, j4)

            !! y interpolation

            CALL fix_cubic(a, b, c, d, y1, f11, y2, f12, y3, f13, y4, f14)
            f10 = a*y**3+b*y**2+c*y+d

            CALL fix_cubic(a, b, c, d, y1, f21, y2, f22, y3, f23, y4, f24)
            f20 = a*y**3+b*y**2+c*y+d

            !!

            !! x interpolation

            CALL fix_line(a, b, x1, f10, x2, f20)
            find_2D = a*x+b

            !!

         ELSE IF (y < ytab(1) .OR. y > ytab(ny)) THEN

            ! y is off the table edge

            IF (nx < 4) STOP 'FIND_2D: Not enough x points in your array for cubic interpolation'
            IF (ny < 2) STOP 'FIND_2D: Not enough y points in your array for linear interpolation'

            IF (x <= xtab(4)) THEN
               i = 2
            ELSE IF (x >= xtab(nx-3)) THEN
               i = nx-2
            ELSE
               i = find_table_integer(x, xtab, nx, ifind)
            END IF

            i1 = i-1
            i2 = i
            i3 = i+1
            i4 = i+2

            x1 = xtab(i1)
            x2 = xtab(i2)
            x3 = xtab(i3)
            x4 = xtab(i4)

            IF (y < ytab(1)) THEN
               j1 = 1
               j2 = 2
            ELSE
               j1 = ny-1
               j2 = ny
            END IF

            y1 = ytab(j1)
            y2 = ytab(j2)

            f11 = ftab(i1, j1)
            f21 = ftab(i2, j1)
            f31 = ftab(i3, j1)
            f41 = ftab(i4, j1)

            f12 = ftab(i1, j2)
            f22 = ftab(i2, j2)
            f32 = ftab(i3, j2)
            f42 = ftab(i4, j2)

            ! x interpolation

            CALL fix_cubic(a, b, c, d, x1, f11, x2, f21, x3, f31, x4, f41)
            f01 = a*x**3+b*x**2+c*x+d

            CALL fix_cubic(a, b, c, d, x1, f12, x2, f22, x3, f32, x4, f42)
            f02 = a*x**3+b*x**2+c*x+d

            ! y interpolation

            CALL fix_line(a, b, y1, f01, y2, f02)
            find_2D = a*y+b

         ELSE

            ! Points exists within table boundardies (normal)

            IF (nx < 4) STOP 'FIND_2D: Not enough x points in your array for cubic interpolation'
            IF (ny < 4) STOP 'FIND_2D: Not enough y points in your array for cubic interpolation'

            IF (x <= xtab(4)) THEN
               i = 2
            ELSE IF (x >= xtab(nx-3)) THEN
               i = nx-2
            ELSE
               i = find_table_integer(x, xtab, nx, ifind)
            END IF

            i1 = i-1
            i2 = i
            i3 = i+1
            i4 = i+2

            x1 = xtab(i1)
            x2 = xtab(i2)
            x3 = xtab(i3)
            x4 = xtab(i4)

            IF (y <= ytab(4)) THEN
               j = 2
            ELSE IF (y >= ytab(ny-3)) THEN
               j = ny-2
            ELSE
               j = find_table_integer(y, ytab, ny, ifind)
            END IF

            j1 = j-1
            j2 = j
            j3 = j+1
            j4 = j+2

            y1 = ytab(j1)
            y2 = ytab(j2)
            y3 = ytab(j3)
            y4 = ytab(j4)

            !

            f11 = ftab(i1, j1)
            f12 = ftab(i1, j2)
            f13 = ftab(i1, j3)
            f14 = ftab(i1, j4)

            f21 = ftab(i2, j1)
            f22 = ftab(i2, j2)
            f23 = ftab(i2, j3)
            f24 = ftab(i2, j4)

            f31 = ftab(i3, j1)
            f32 = ftab(i3, j2)
            f33 = ftab(i3, j3)
            f34 = ftab(i3, j4)

            f41 = ftab(i4, j1)
            f42 = ftab(i4, j2)
            f43 = ftab(i4, j3)
            f44 = ftab(i4, j4)

            ! x interpolation

            CALL fix_cubic(a, b, c, d, x1, f11, x2, f21, x3, f31, x4, f41)
            f01 = a*x**3+b*x**2+c*x+d

            CALL fix_cubic(a, b, c, d, x1, f12, x2, f22, x3, f32, x4, f42)
            f02 = a*x**3+b*x**2+c*x+d

            CALL fix_cubic(a, b, c, d, x1, f13, x2, f23, x3, f33, x4, f43)
            f03 = a*x**3+b*x**2+c*x+d

            CALL fix_cubic(a, b, c, d, x1, f14, x2, f24, x3, f34, x4, f44)
            f04 = a*x**3+b*x**2+c*x+d

            CALL fix_cubic(a, b, c, d, y1, f01, y2, f02, y3, f03, y4, f04)
            findy = a*y**3+b*y**2+c*y+d

            ! y interpolation

            CALL fix_cubic(a, b, c, d, y1, f11, y2, f12, y3, f13, y4, f14)
            f10 = a*y**3+b*y**2+c*y+d

            CALL fix_cubic(a, b, c, d, y1, f21, y2, f22, y3, f23, y4, f24)
            f20 = a*y**3+b*y**2+c*y+d

            CALL fix_cubic(a, b, c, d, y1, f31, y2, f32, y3, f33, y4, f34)
            f30 = a*y**3+b*y**2+c*y+d

            CALL fix_cubic(a, b, c, d, y1, f41, y2, f42, y3, f43, y4, f44)
            f40 = a*y**3+b*y**2+c*y+d

            CALL fix_cubic(a, b, c, d, x1, f10, x2, f20, x3, f30, x4, f40)
            findx = a*x**3+b*x**2+c*x+d

            ! Final result is an average over each direction
            find_2D = (findx+findy)/2.

         END IF

      ELSE

         STOP 'FIND_2D: order for interpolation not specified correctly'

      END IF

   END FUNCTION find_2D

   ! REAL FUNCTION find_3D(x, xin, y, yin, z, zin, fin, nx, ny, nz, iorder, ifind, iinterp)

   !    ! A 3D interpolation routine to find value f(x,y,z) given a function evalated on arrays
   !    ! The linear version implemented here is also know as 'trilinear interpolation'
   !    ! TODO: Implement loops over coordinates to avoid repetition
   !    IMPLICIT NONE
   !    REAL, INTENT(IN) :: x
   !    REAL, INTENT(IN) :: xin(nx)
   !    REAL, INTENT(IN) :: y
   !    REAL, INTENT(IN) :: yin(ny)
   !    REAL, INTENT(IN) :: z
   !    REAL, INTENT(IN) :: zin(nz)
   !    REAL, INTENT(IN) :: fin(nx, ny, nz)
   !    INTEGER, INTENT(IN) :: nx
   !    INTEGER, INTENT(IN) :: ny
   !    INTEGER, INTENT(IN) :: nz
   !    INTEGER, INTENT(IN) :: iorder
   !    INTEGER, INTENT(IN) :: ifind
   !    INTEGER, INTENT(IN) :: iinterp
   !    REAL :: x1, x2, y1, y2, z1, z2
   !    REAL :: F(2, 2, 2)
   !    REAL :: V(2, 2, 2)
   !    INTEGER :: ix, ix1, ix2, iy, iy1, iy2, iz, iz1, iz2

   !    ! If the value required is off the table edge then this will halt

   !    ! iorder = 1 => linear interpolation

   !    ! ifind = 1 => find x in xtab by crudely searching from x(1) to x(n)
   !    ! ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
   !    ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

   !    ! iinterp = 1 => Uses cubic polynomials for interpolation

   !    IF (iinterp .NE. 1) STOP 'FIND_3D: No Lagrange polynomials for you, only regular polynomails, iinterp=1'

   !    IF (xin(1) > xin(nx)) STOP 'FIND_3D: x array in wrong order'
   !    IF (yin(1) > yin(ny)) STOP 'FIND_3D: y array in wrong order'
   !    IF (zin(1) > zin(nz)) STOP 'FIND_3D: z array in wrong order'

   !    IF (iorder == 1) THEN

   !       IF (nx < 2) STOP 'FIND_3D: Not enough x points in your array for linear interpolation'
   !       IF (ny < 2) STOP 'FIND_3D: Not enough y points in your array for linear interpolation'
   !       IF (nz < 2) STOP 'FIND_3D: Not enough z points in your array for linear interpolation'

   !       !! Get the x,y,z values !!

   !       ! Get the integer coordinates in the x direction
   !       ix = find_table_integer(x, xin, nx, ifind)
   !       IF(ix==0) THEN
   !          ix=1
   !       ELSE IF(ix==nx) THEN
   !          ix=nx-1
   !       END IF
   !       ix1 = ix
   !       ix2 = ix+1

   !       ! Get the x values
   !       x1 = xin(ix1)
   !       x2 = xin(ix2)

   !       ! Get the integer coordinates in the y direction
   !       iy = find_table_integer(y, yin, ny, ifind)
   !       IF(iy==0) THEN
   !          iy=1
   !       ELSE IF(iy==ny) THEN
   !          iy=ny-1
   !       END IF
   !       iy1 = iy
   !       iy2 = iy+1

   !       ! Get the y values
   !       y1 = yin(iy1)
   !       y2 = yin(iy2)

   !       ! Get the integer coordinates in the z direction
   !       iz = find_table_integer(z, zin, nz, ifind)
   !       IF(iz==0) THEN
   !          iz=1
   !       ELSE IF(iz==nz) THEN
   !          iz=nz-1
   !       END IF
   !       iz1 = iz
   !       iz2 = iz+1

   !       ! Get the z values
   !       z1 = zin(iz1)
   !       z2 = zin(iz2)

   !       !! Function values !!

   !       F(1, 1, 1) = fin(ix1, iy1, iz1)
   !       F(1, 1, 2) = fin(ix1, iy1, iz2)
   !       F(1, 2, 1) = fin(ix1, iy2, iz1)
   !       F(1, 2, 2) = fin(ix1, iy2, iz2)
   !       F(2, 1, 1) = fin(ix2, iy1, iz1)
   !       F(2, 1, 2) = fin(ix2, iy1, iz2)
   !       F(2, 2, 1) = fin(ix2, iy2, iz1)
   !       F(2, 2, 2) = fin(ix2, iy2, iz2)

   !       !! !!

   !       V(1, 1, 1) = (x2-x)*(y2-y)*(z2-z)*F(1, 1, 1)
   !       V(1, 1, 2) = (x2-x)*(y2-y)*(z-z1)*F(1, 1, 2)
   !       V(1, 2, 1) = (x2-x)*(y-y1)*(z2-z)*F(1, 2, 1)
   !       V(1, 2, 2) = (x2-x)*(y-y1)*(z-z1)*F(1, 2, 2)
   !       V(2, 1, 1) = (x-x1)*(y2-y)*(z2-z)*F(2, 1, 1)
   !       V(2, 1, 2) = (x-x1)*(y2-y)*(z-z1)*F(2, 1, 2)
   !       V(2, 2, 1) = (x-x1)*(y-y1)*(z2-z)*F(2, 2, 1)
   !       V(2, 2, 2) = (x-x1)*(y-y1)*(z-z1)*F(2, 2, 2)

   !       find_3D = sum(V)/((x2-x1)*(y2-y1)*(z2-z1))

   !    ELSE

   !       STOP 'FIND_3D: order for interpolation not specified correctly, only linear implemented'

   !    END IF

   ! END FUNCTION find_3D

   REAL FUNCTION find_3D(x, xin, y, yin, z, zin, fin, nx, ny, nz, iorder, ifind, iinterp)

      ! A 3D interpolation routine to find value f(x,y,z) given a function evalated on arrays
      ! The linear version implemented here is also know as 'trilinear interpolation'
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xin(nx)
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: yin(ny)
      REAL, INTENT(IN) :: z
      REAL, INTENT(IN) :: zin(nz)
      REAL, INTENT(IN) :: fin(nx, ny, nz)
      INTEGER, INTENT(IN) :: nx
      INTEGER, INTENT(IN) :: ny
      INTEGER, INTENT(IN) :: nz
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      REAL :: xx(3), x12(3, 2), Dx(3)
      INTEGER :: ix(3), ix12(3, 2), nnx(3)
      INTEGER :: i, j, k, ii, jj, kk, d, di
      REAL :: F(2, 2, 2)
      REAL :: V(2, 2, 2)

      ! If the value required is off the table edge then this will halt

      ! iorder = 1 => linear interpolation

      ! ifind = 1 => find x in xtab by crudely searching from x(1) to x(n)
      ! ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
      ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

      ! iinterp = 1 => Uses cubic polynomials for interpolation

      IF (iinterp .NE. 1) STOP 'FIND_3D: No Lagrange polynomials for you, only regular polynomails, iinterp=1'

      IF (xin(1) > xin(nx)) STOP 'FIND_3D: x array in wrong order'
      IF (yin(1) > yin(ny)) STOP 'FIND_3D: y array in wrong order'
      IF (zin(1) > zin(nz)) STOP 'FIND_3D: z array in wrong order'

      IF(iorder == 0) THEN

         DO d = 1, 3
            IF(d==1) ix(d) = find_table_integer(x, xin, nx, ifind)
            IF(d==2) ix(d) = find_table_integer(y, yin, ny, ifind)
            IF(d==3) ix(d) = find_table_integer(z, zin, nz, ifind)
         END DO

         DO d = 1, 3
            IF(ix(d)==0) ix(d)=1
         END DO

         find_3D = fin(ix(1), ix(2), ix(3))

      ELSE IF (iorder == 1) THEN

         xx(1)=x
         xx(2)=y
         xx(3)=z

         nnx(1)=nx
         nnx(2)=ny
         nnx(3)=nz

         DO d = 1, 3
            IF (nnx(d) < 2) STOP 'FIND_3D: Not enough points in your array for linear interpolation'
         END DO

         !! Get the x,y,z values !!

         ! Get the integer coordinates in each direction
         DO d = 1, 3

            ! Get the integer coordinates of the 'cube' corners that encompass the point in each dimensions
            IF(d == 1) ix(1) = find_table_integer(x, xin, nx, ifind)
            IF(d == 2) ix(2) = find_table_integer(y, yin, ny, ifind)
            IF(d == 3) ix(3) = find_table_integer(z, zin, nz, ifind)

            ! Correct for the case of these integers coordinates being at the edge of the box
            IF(ix(d) == 0) THEN
               ix(d) = 1
            ELSE IF(ix(d) == nnx(d)) THEN
               ix(d) = nnx(d)-1
            END IF
            ix12(d,1) = ix(d)
            ix12(d,2) = ix(d)+1

            ! Get the function value at the corners of the 'cube'
            IF(d == 1) THEN
               x12(d,1) = xin(ix12(d,1))
               x12(d,2) = xin(ix12(d,2))
            ELSE IF(d == 2) THEN
               x12(d,1) = yin(ix12(d,1))
               x12(d,2) = yin(ix12(d,2))
            ELSE IF(d == 3) THEN
               x12(d,1) = zin(ix12(d,1))
               x12(d,2) = zin(ix12(d,2))
            END IF

         END DO 

         ! Do the interpolation
         DO k = 1, 2
            DO j = 1, 2
               DO i = 1, 2

                  ! Evaluate the function at the cube corners
                  F(i, j, k) = fin(ix12(1,i), ix12(2,j), ix12(3,k))

                  ! Calculate the volume of the solid in the opposite corner
                  ii = 3-i ! Swaps 1 <--> 2 (e.g., ii is 2 if i is 1 and visa versa)
                  jj = 3-j ! Swaps 1 <--> 2
                  kk = 3-k ! Swaps 1 <--> 2
                  DO d = 1, 3
                     IF(d==1) di=ii
                     IF(d==2) di=jj
                     IF(d==3) di=kk
                     Dx(d) = (x12(d,di)-xx(d))*(2*di-3) ! Last part gets the sign right (x2-x) or (x-x1), not negatives (x-x2)!!
                  END DO

                  ! Weight the corner function values by their corresponding weights (opposite cubeoid)
                  V(i, j, k) = F(i, j, k)
                  DO d = 1, 3
                     V(i, j, k)=V(i, j, k)*Dx(d)
                  END DO

               END DO
            END DO
         END DO

         ! Finally normalise the result by dividing by the cube volume
         find_3D = sum(V)
         DO d = 1,3
            find_3D = find_3D/(x12(d,2)-x12(d,1))
         END DO

      ELSE

         STOP 'FIND_3D: order for interpolation not specified correctly, only linear implemented'

      END IF

   END FUNCTION find_3D

   SUBROUTINE interpolate_array(x1, y1, n1, x2, y2, n2, iorder, ifind, iinterp)

      ! Interpolates array 'x1-y1' onto new 'x' values x2 and output y2
      ! TODO: This could be more efficient because it currently does 'find integer' every time
      IMPLICIT NONE
      REAL, INTENT(IN) :: x1(n1), y1(n1), x2(n2)
      REAL, INTENT(OUT) :: y2(n2)
      INTEGER, INTENT(IN) :: n1, n2
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      INTEGER :: i
    
      DO i = 1, n2
         y2(i) = find(x2(i), x1, y1, n1, iorder, ifind, iinterp)
      END DO

   END SUBROUTINE interpolate_array

END MODULE interpolate
