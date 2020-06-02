MODULE interpolate

   USE table_integer
   USE array_operations
   USE special_functions

   IMPLICIT NONE

   PRIVATE

   ! Functions
   PUBLIC :: find
   PUBLIC :: interpolate_array
   PUBLIC :: init_interpolator
   PUBLIC :: evaluate_interpolator

   ! Integers for interpolation polynomials
   PUBLIC :: iinterp_polynomial
   PUBLIC :: iinterp_Lagrange

   ! Types
   PUBLIC :: interpolator

   INTEGER, PARAMETER :: iinterp_polynomial = 1
   INTEGER, PARAMETER :: iinterp_Lagrange = 2

   LOGICAL, PARAMETER :: centred_interpolator = .TRUE.

   INTERFACE find
      MODULE PROCEDURE find_1D
      MODULE PROCEDURE find_2D
      MODULE PROCEDURE find_3D
   END INTERFACE find

   INTERFACE init_interpolator
      MODULE PROCEDURE init_interpolator_1D
   END INTERFACE init_interpolator

   INTERFACE evaluate_interpolator
      MODULE PROCEDURE evaluate_interpolator_1D
   END INTERFACE evaluate_interpolator

   TYPE interpolator

      REAL, ALLOCATABLE :: x(:), xmid(:)
      REAL, ALLOCATABLE :: a0(:), a1(:), a2(:), a3(:)
      INTEGER :: iorder, ifind
      INTEGER :: n
      LOGICAL :: logx, logf
      LOGICAL :: extrap

   END TYPE interpolator

CONTAINS

   REAL FUNCTION find_1D(x, xin, yin, n, iorder, ifind, iinterp)

      ! Given two arrays x and y this routine interpolates to find the y_i value at position x_i
      ! Care should be chosen to insert x, xtab, ytab as log if this might give beter results
      ! Extrapolates if the value is off either end of the array   
      ! If the value required is off the table edge the extrapolation is always linear
      REAL, INTENT(IN) :: x
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: xin(n)
      REAL, INTENT(IN) :: yin(n)
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      REAL ::  xtab(n), ytab(n)
      REAL :: a0, a1, a2, a3
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

      ! Reverse the arrays in this case
      IF (xtab(1) > xtab(n)) THEN       
         CALL reverse_array(xtab)
         CALL reverse_array(ytab)
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
            i = find_table_integer(x, xtab, ifind)
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
               CALL fix_linear(a1, a0, [x1, x2], [y1, y2])        
               find_1D = linear_polynomial(x, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 1, [x1, x2], [y1, y2])
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
               CALL fix_linear(a1, a0, [x1, x2], [y1, y2])
               find_1D = linear_polynomial(x, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 1, [x1, x2], [y1, y2])
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

               i = find_table_integer(x, xtab, ifind)

               x1 = xtab(i)
               x2 = xtab(i+1)

               y1 = ytab(i)
               y2 = ytab(i+1)

            END IF

            IF (iinterp == iinterp_polynomial) THEN
               CALL fix_linear(a1, a0, [x1, x2], [y1, y2])  
               find_1D = linear_polynomial(x, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 1, [x1, x2], [y1, y2])
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

               ELSE

                  STOP 'FIND_1D: Error, something went wrong'

               END IF

               IF (iinterp == iinterp_polynomial) THEN
                  CALL fix_quadratic(a2, a1, a0, [x1, x2, x3], [y1, y2, y3])
                  find_1D = quadratic_polynomial(x, a2, a1, a0)
               ELSE IF (iinterp == iinterp_Lagrange) THEN
                  find_1D = Lagrange_polynomial(x, 2, [x1, x2, x3], [y1, y2, y3])
               ELSE
                  STOP 'FIND_1D: Error, method not specified correctly'
               END IF

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

               IF (iinterp == iinterp_polynomial) THEN
                  ! In this case take the average of two separate quadratic spline values
                  find_1D = 0.
                  CALL fix_quadratic(a2, a1, a0, [x1, x2, x3], [y1, y2, y3])     
                  find_1D = find_1D+quadratic_polynomial(x, a2, a1, a0)/2.
                  CALL fix_quadratic(a2, a1, a0, [x2, x3, x4], [y2, y3, y4])
                  find_1D = find_1D+quadratic_polynomial(x, a2, a1, a0)/2.
               ELSE IF (iinterp == iinterp_Lagrange) THEN
                  ! In this case take the average of two quadratic Lagrange polynomials
                  L1 = Lagrange_polynomial(x, 2, [x1, x2, x3], [y1, y2, y3])
                  L2 = Lagrange_polynomial(x, 2, [x2, x3, x4], [y2, y3, y4])
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

            IF (iinterp == iinterp_polynomial) THEN
               CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [y1, y2, y3, y4])
               find_1D = cubic_polynomial(x, a3, a2, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, 3, [x1, x2, x3, x4], [y1, y2, y3, y4])
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
      INTEGER, INTENT(IN) :: nx
      INTEGER, INTENT(IN) :: ny
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xin(nx)
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: yin(ny)
      REAL, INTENT(IN) :: fin(nx, ny)
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      REAL ::  xtab(nx), ytab(ny), ftab(nx, ny)
      REAL :: a3, a2, a1, a0
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

         ix = find_table_integer(x, xin, ifind)
         IF(ix == 0) ix = 1
         iy = find_table_integer(y, yin, ifind)
         IF(iy == 0) iy = 1

         find_2D = fin(ix, iy)

      ELSE IF(iorder == 1) THEN

         IF (nx < 2) STOP 'FIND_2D: Not enough x points in your array for linear interpolation'
         IF (ny < 2) STOP 'FIND_2D: Not enough y points in your array for linear interpolation'

         !! Get the x,y values !!

         ! Get the integer coordinates in the x direction
         ix = find_table_integer(x, xin, ifind)
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
         iy = find_table_integer(y, yin, ifind)
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
               j = find_table_integer(y, ytab, ifind)
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

            CALL fix_cubic(a3, a2, a1, a0, [y1, y2, y3, y4], [f11, f12, f13, f14])
            f10 = cubic_polynomial(y, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [y1, y2, y3, y4], [f21, f22, f23, f24])      
            f20 = cubic_polynomial(y, a3, a2, a1, a0)

            !!

            !! x interpolation

            CALL fix_linear(a1, a0, [x1, x2], [f10, f20])
            find_2D = linear_polynomial(x, a1, a0)

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
               i = find_table_integer(x, xtab, ifind)
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

            CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f11, f21, f31, f41])
            f01 = cubic_polynomial(x, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f12, f22, f32, f42])     
            f02 = cubic_polynomial(x, a3, a2, a1, a0)

            ! y interpolation

            CALL fix_linear(a1, a0, [y1, y2], [f01, f02])
            find_2D = linear_polynomial(y, a1, a0)

         ELSE

            ! Points exists within table boundardies (normal)

            IF (nx < 4) STOP 'FIND_2D: Not enough x points in your array for cubic interpolation'
            IF (ny < 4) STOP 'FIND_2D: Not enough y points in your array for cubic interpolation'

            IF (x <= xtab(4)) THEN
               i = 2
            ELSE IF (x >= xtab(nx-3)) THEN
               i = nx-2
            ELSE
               i = find_table_integer(x, xtab, ifind)
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
               j = find_table_integer(y, ytab, ifind)
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

            CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f11, f21, f31, f41])
            f01 = cubic_polynomial(x, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f12, f22, f32, f42]) 
            f02 = cubic_polynomial(x, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f13, f23, f33, f43])        
            f03 = cubic_polynomial(x, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f14, f24, f34, f44])
            f04 = cubic_polynomial(x, a3, a2, a1, a0)
            
            CALL fix_cubic(a3, a2, a1, a0, [y1, y2, y3, y4], [f01, f02, f03, f04])
            findy = cubic_polynomial(y, a3, a2, a1, a0)
            
            ! y interpolation

            CALL fix_cubic(a3, a2, a1, a0, [y1, y2, y3, y4], [f11, f12, f13, f14])
            f10 = cubic_polynomial(y, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [y1, y2, y3, y4], [f21, f22, f23, f24])
            f20 = cubic_polynomial(y, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [y1, y2, y3, y4], [f31, f32, f33, f34])
            f30 = cubic_polynomial(y, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [y1, y2, y3, y4], [f41, f42, f43, f44])
            f40 = cubic_polynomial(y, a3, a2, a1, a0)

            CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f10, f20, f30, f40])
            findx = cubic_polynomial(x, a3, a2, a1, a0)

            ! Final result is an average over each direction
            find_2D = (findx+findy)/2.

         END IF

      ELSE

         STOP 'FIND_2D: order for interpolation not specified correctly'

      END IF

   END FUNCTION find_2D

   REAL FUNCTION find_3D(x, xin, y, yin, z, zin, fin, nx, ny, nz, iorder, ifind, iinterp)

      ! A 3D interpolation routine to find value f(x,y,z) given a function evalated on arrays
      ! The linear version implemented here is also know as 'trilinear interpolation'
      INTEGER, INTENT(IN) :: nx
      INTEGER, INTENT(IN) :: ny
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xin(nx)
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: yin(ny)
      REAL, INTENT(IN) :: z
      REAL, INTENT(IN) :: zin(nz)
      REAL, INTENT(IN) :: fin(nx, ny, nz)
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
            IF(d==1) ix(d) = find_table_integer(x, xin, ifind)
            IF(d==2) ix(d) = find_table_integer(y, yin, ifind)
            IF(d==3) ix(d) = find_table_integer(z, zin, ifind)
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
            IF(d == 1) ix(1) = find_table_integer(x, xin, ifind)
            IF(d == 2) ix(2) = find_table_integer(y, yin, ifind)
            IF(d == 3) ix(3) = find_table_integer(z, zin, ifind)

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
                     IF(d==1) THEN
                        di=ii
                     ELSE IF(d==2) THEN
                        di=jj
                     ELSE IF(d==3) THEN
                        di=kk
                     ELSE
                        STOP 'FIND_3D: Error, something is very wrong'
                     END IF
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
      INTEGER, INTENT(IN) :: n1
      REAL, INTENT(IN) :: x1(n1)
      REAL, INTENT(IN) :: y1(n1)
      INTEGER, INTENT(IN) :: n2
      REAL, INTENT(IN) :: x2(n2)
      REAL, INTENT(OUT) :: y2(n2)
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      INTEGER :: i
    
      DO i = 1, n2
         y2(i) = find(x2(i), x1, y1, n1, iorder, ifind, iinterp)
      END DO

   END SUBROUTINE interpolate_array

   SUBROUTINE init_interpolator_1D(x, f, iorder, interp, extrap, logx, logf)

      ! Initialise an interpolator
      USE basic_operations
      REAL, INTENT(IN) :: x(:)                  ! Input data x
      REAL, INTENT(IN) :: f(:)                  ! Input data f(x)
      INTEGER, INTENT(IN) :: iorder             ! Order at which to create interpolator
      TYPE(interpolator), INTENT(OUT) :: interp ! Interpolator type
      LOGICAL, OPTIONAL, INTENT(IN) :: extrap   ! Should interpolator extrapolate beyond x?
      LOGICAL, OPTIONAL, INTENT(IN) :: logx     ! Should interpolator take the logarithm of x?
      LOGICAL, OPTIONAL, INTENT(IN) :: logf     ! Should interpolator take the logarithm of y?
      REAL, ALLOCATABLE :: xx(:), ff(:)
      REAL :: a0, a1, a2, a3
      REAL :: f1, f2, f3, f4
      REAL :: x1, x2, x3, x4
      REAL :: xm
      INTEGER :: i, n
      INTEGER :: i1, i2, i3, i4
      INTEGER, PARAMETER :: ifind_default = ifind_split

      ! Check the sizes of the input data arrays
      n = size(x)
      IF (n /= size(f)) STOP 'INIT_INTERPOLATOR: Error, input x and f data should be the same size'
      interp%n = n
      
      ALLOCATE(xx(n), ff(n))

      IF (present_and_correct(extrap)) THEN
         interp%extrap = .TRUE.
      ELSE
         interp%extrap = .FALSE.
      END IF

      ! Sort out logs
      IF(present_and_correct(logx)) THEN
         xx = log(x)
         interp%logx = .TRUE.
      ELSE
         xx = x
         interp%logx = .FALSE.
      END IF
      IF(present_and_correct(logf)) THEN
         ff = log(f)
         interp%logf = .TRUE.
      ELSE
         ff = f
         interp%logf = .FALSE.
      END IF

      ! Allocate arrays for the interpolator x data
      ALLOCATE(interp%x(n))
      interp%x = xx

      IF (centred_interpolator) ALLOCATE(interp%xmid(n-1))

      ! Allocate arrays for the interpolator coefficients
      ALLOCATE(interp%a0(n-1), interp%a1(n-1), interp%a2(n-1), interp%a3(n-1))
      interp%a0 = 0.
      interp%a1 = 0.
      interp%a2 = 0.
      interp%a3 = 0.

      ! Set the interpolation order
      interp%iorder = iorder

      ! If the x data is regular spaced then remember this, otherwise default find
      IF (regular_spacing(x)) THEN
         interp%ifind = ifind_linear
      ELSE
         interp%ifind = ifind_default
      END IF

      ! Default values because a3, a2 will be zero for linear polynomials
      a0 = 0.
      a1 = 0.
      a2 = 0.
      a3 = 0.

      ! Loop over all n-1 sections of the input data
      DO i = 1, n-1

         IF (iorder == 1) THEN

            i1 = i
            i2 = i+1

            x1 = xx(i1)
            x2 = xx(i2)

            f1 = ff(i1)
            f2 = ff(i2)
            
            IF (centred_interpolator) THEN
               xm = (xx(i1)+xx(i2))/2.
               CALL fix_centred_linear(a1, a0, xm, [x1, x2], [f1, f2])
            ELSE
               CALL fix_linear(a1, a0, [x1, x2], [f1, f2])
            END IF

         ELSE IF (iorder == 3) THEN

            ! Deal with the indices for the first and last section and general case
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

            x1 = xx(i1)
            x2 = xx(i2)
            x3 = xx(i3)
            x4 = xx(i4)

            f1 = ff(i1)
            f2 = ff(i2)
            f3 = ff(i3)
            f4 = ff(i4)

            IF (centred_interpolator) THEN
               xm = (xx(i)+xx(i+1))/2.
               CALL fix_centred_cubic(a3, a2, a1, a0, xm, [x1, x2, x3, x4], [f1, f2, f3, f4])
            ELSE
               CALL fix_cubic(a3, a2, a1, a0, [x1, x2, x3, x4], [f1, f2, f3, f4])
            END IF          

         ELSE

            STOP 'INIT_INTERPOLATOR: Error, your order is not supported'

         END IF

         ! Fill the polynomial coefficients in the interpolation type
         IF(centred_interpolator) interp%xmid(i) = xm
         interp%a0(i) = a0
         interp%a1(i) = a1
         interp%a2(i) = a2
         interp%a3(i) = a3

      END DO

   END SUBROUTINE init_interpolator_1D

   REAL FUNCTION evaluate_interpolator_1D(x, interp)

      ! Evaluates the value f(x) from the interpolator
      REAL, INTENT(IN) :: x
      TYPE(interpolator), INTENT(IN) :: interp
      INTEGER :: i, n
      REAL :: xx
      LOGICAL, PARAMETER :: centred = .TRUE.

      ! Change input x to log if necessary
      IF (interp%logx) THEN
         xx = log(x)
      ELSE
         xx = x
      END IF

      n = interp%n

      ! Calculate the index to use
      IF (xx < interp%x(1)) THEN
         IF(interp%extrap) THEN
            i = 1
         ELSE
            STOP 'EVALUATE_INTERPOLATOR_1D: Error, desired x value is below range'
         END IF
      ELSE IF (xx == interp%x(n)) THEN
         i = n-1 ! Cannot use find_table_integer here because this would return 'n', rather than 'n-1'
      ELSE IF (xx > interp%x(n)) THEN
         IF (interp%extrap) THEN
            i = n-1
         ELSE
            STOP 'EVALUATE_INTERPOLATOR_1D: Error, desired x value is above range'
         END IF
      ELSE
         i = find_table_integer(xx, interp%x, interp%ifind)
      END IF

      ! Evaluate the interpolation
      IF (interp%iorder == 1) THEN
         IF (centred_interpolator) THEN
            evaluate_interpolator_1D = centred_linear_polynomial(xx, interp%xmid(i), interp%a1(i), interp%a0(i))
         ELSE
            evaluate_interpolator_1D = linear_polynomial(xx, interp%a1(i), interp%a0(i))
         END IF
      ELSE IF (interp%iorder == 3) THEN
         IF (centred_interpolator) THEN
            evaluate_interpolator_1D = centred_cubic_polynomial(xx, interp%xmid(i), interp%a3(i), interp%a2(i), interp%a1(i), interp%a0(i))
         ELSE
            evaluate_interpolator_1D = cubic_polynomial(xx, interp%a3(i), interp%a2(i), interp%a1(i), interp%a0(i))
         END IF
      ELSE
         STOP 'EVALUATE_INTERPOLATOR_1D: Error, your polynomial order is not supported'
      END IF

      ! Exponentiate result if necessary
      IF (interp%logf) evaluate_interpolator_1D = exp(evaluate_interpolator_1D)

   END FUNCTION evaluate_interpolator_1D

   ! SUBROUTINE rebin_interpolate(x, y, xmin, xmax, n1, n2)

   !    IMPLICIT NONE
   !    REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)
   !    REAL, ALLOCATABLE, INTENT(INOUT) :: y(:)
   !    REAL, INTENT(IN) :: xmin
   !    REAL, INTENT(IN) :: xmax
   !    INTEGER, INTENT(IN) :: n1
   !    INTEGER, INTENT(IN) :: n2

   !    ! Rebin on a log-linear axis
   !    IF (rebin) THEN
   !       kmin = kmin_rebin
   !       kmax = kmax_rebin
   !       nk = nk_rebin
   !       CALL fill_array(log(kmin), log(kmax), k2, nk)
   !       k2 = exp(k2)
   !       ALLOCATE (P2(nk))
   !       CALL interpolate_array(log(k), log(P), n, log(k2), P2, nk, 3, 3, 2)
   !       P2 = exp(P2)
   !       DEALLOCATE (k, P)
   !       ALLOCATE (k(nk), P(nk))
   !       k = k2
   !       P = P2
   !       n = nk
   !       DEALLOCATE (k2, P2)
   !    END IF

   ! END SUBROUTINE rebin_interpolate

END MODULE interpolate
