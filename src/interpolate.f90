MODULE interpolate

   ! TODO: Combine the integer finding that is carried out separately in 'find' and 'evaluate_interpolator'
   ! TODO: Combine the interpolation that is carried out separatly in 'find' and 'evaluate_interpolator'
   ! TODO: Functions to return min/max x, y, ... values from interpolators

   USE basic_operations
   USE array_operations
   USE table_integer
   USE special_functions

   IMPLICIT NONE

   PRIVATE

   ! Functions
   PUBLIC :: find
   PUBLIC :: interpolate_array
   PUBLIC :: rebin_array
   PUBLIC :: init_interpolator
   PUBLIC :: evaluate_interpolator
   PUBLIC :: inverse_interpolator

   ! Switches for interpolation polynomials
   PUBLIC :: iinterp_polynomial
   PUBLIC :: iinterp_Lagrange
   PUBLIC :: iinterp_centred

   ! Switches for extrapolation
   PUBLIC :: iextrap_no
   PUBLIC :: iextrap_std
   PUBLIC :: iextrap_lin
   PUBLIC :: iextrap_zero
   PUBLIC :: iextrap_near

   ! Types
   PUBLIC :: interpolator1D
   PUBLIC :: interpolator2D
   PUBLIC :: interpolator3D

   ! Interpolation methods
   INTEGER, PARAMETER :: iinterp_polynomial = 1
   INTEGER, PARAMETER :: iinterp_Lagrange = 2
   INTEGER, PARAMETER :: iinterp_centred = 3

   ! Extrapolation methods
   INTEGER, PARAMETER :: iextrap_no = 0
   INTEGER, PARAMETER :: iextrap_std = 1
   INTEGER, PARAMETER :: iextrap_lin = 2
   INTEGER, PARAMETER :: iextrap_zero = 3
   INTEGER, PARAMETER :: iextrap_near = 4

   ! Default schemes
   INTEGER, PARAMETER :: ifind_interpolator_default = ifind_split
   INTEGER, PARAMETER :: ifind_inverse_interpolator = ifind_split
   INTEGER, PARAMETER :: iinterp_inverse_interpolator = iinterp_Lagrange

   INTERFACE find
      MODULE PROCEDURE find_1D
      MODULE PROCEDURE find_2D
      MODULE PROCEDURE find_3D
   END INTERFACE find

   INTERFACE init_interpolator
      MODULE PROCEDURE init_interpolator_1D
      MODULE PROCEDURE init_interpolator_2D
      MODULE PROCEDURE init_interpolator_3D
   END INTERFACE init_interpolator

   INTERFACE evaluate_interpolator
      MODULE PROCEDURE evaluate_interpolator_1D
      MODULE PROCEDURE evaluate_interpolator_2D
      MODULE PROCEDURE evaluate_interpolator_3D
   END INTERFACE evaluate_interpolator

   INTERFACE inverse_interpolator
      MODULE PROCEDURE inverse_interpolator_1D
   END INTERFACE inverse_interpolator

   TYPE interpolator1D
      REAL, ALLOCATABLE :: x(:), x0(:)
      REAL, ALLOCATABLE :: f(:)
      REAL, ALLOCATABLE :: a0(:), a1(:), a2(:), a3(:)
      REAL :: xmin, xmax
      INTEGER :: iorder, iextrap
      INTEGER :: ifind
      INTEGER :: n
      LOGICAL :: logx, logf
      LOGICAL :: store
   END TYPE interpolator1D

   TYPE interpolator2D
      REAL, ALLOCATABLE :: x(:), y(:)
      REAL, ALLOCATABLE :: f(:, :)
      REAL, ALLOCATABLE :: x0(:, :), y0(:, :), xl0(:, :), yl0(:, :)
      REAL, ALLOCATABLE :: ax0(:, :), ax1(:, :), ax2(:, :), ax3(:, :)
      REAL, ALLOCATABLE :: ay0(:, :), ay1(:, :), ay2(:, :), ay3(:, :)
      REAL, ALLOCATABLE :: bx0(:, :), bx1(:, :), by0(:, :), by1(:, :)
      REAL :: xmin, xmax, ymin, ymax
      INTEGER :: iorder, iextrap
      INTEGER :: ifindx, ifindy 
      INTEGER :: nx, ny
      LOGICAL :: logx, logy, logf
      LOGICAL :: store
      TYPE(interpolator1D) :: interp1D
   END TYPE interpolator2D

   TYPE interpolator3D
      REAL, ALLOCATABLE :: x(:), y(:), z(:)
      REAL, ALLOCATABLE :: f(:, :, :)
      !REAL, ALLOCATABLE :: x0(:, :), y0(:, :), xl0(:, :), yl0(:, :)
      !REAL, ALLOCATABLE :: ax0(:, :), ax1(:, :), ax2(:, :), ax3(:, :)
      !REAL, ALLOCATABLE :: ay0(:, :), ay1(:, :), ay2(:, :), ay3(:, :)
      !REAL, ALLOCATABLE :: bx0(:, :), bx1(:, :), by0(:, :), by1(:, :)
      REAL :: xmin, xmax, ymin, ymax, zmin, zmax
      INTEGER :: iorder, iextrap
      INTEGER :: ifindx, ifindy, ifindz
      INTEGER :: nx, ny, nz
      LOGICAL :: logx, logy, logz, logf
      LOGICAL :: store
   END TYPE interpolator3D

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
      REAL :: dx, x0
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
      ! iinterp = 3 => Uses centred polynomials for interpolation

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
               CALL fix_polynomial(a1, a0, [x1, x2], [y1, y2])        
               find_1D = polynomial(x, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, [x1, x2], [y1, y2])
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
               CALL fix_polynomial(a1, a0, [x1, x2], [y1, y2])
               find_1D = polynomial(x, a1, a0)
            ELSE IF (iinterp == iinterp_centred) THEN
               x0 = (x1+x2)/2.
               CALL fix_centred_polynomial(a1, a0, x0, [x1, x2], [y1, y2])
               find_1D = centred_polynomial(x, x0, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, [x1, x2], [y1, y2])
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
               CALL fix_polynomial(a1, a0, [x1, x2], [y1, y2])  
               find_1D = polynomial(x, a1, a0)
            ELSE IF (iinterp == iinterp_centred) THEN
               x0 = (x1+x2)/2.
               CALL fix_centred_polynomial(a1, a0, x0, [x1, x2], [y1, y2])
               find_1D = centred_polynomial(x, x0, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, [x1, x2], [y1, y2])
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
                  CALL fix_polynomial(a2, a1, a0, [x1, x2, x3], [y1, y2, y3])
                  find_1D = polynomial(x, a2, a1, a0)
               ELSE IF (iinterp == iinterp_centred) THEN
                  x0 = (x1+x2+x3)/3.
                  CALL fix_centred_polynomial(a2, a1, a0, x0, [x1, x2, x3], [y1, y2, y3])
                  find_1D = centred_polynomial(x, x0, a2, a1, a0)
               ELSE IF (iinterp == iinterp_Lagrange) THEN
                  find_1D = Lagrange_polynomial(x, [x1, x2, x3], [y1, y2, y3])
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
                  CALL fix_polynomial(a2, a1, a0, [x1, x2, x3], [y1, y2, y3])     
                  find_1D = find_1D+polynomial(x, a2, a1, a0)/2.
                  CALL fix_polynomial(a2, a1, a0, [x2, x3, x4], [y2, y3, y4])
                  find_1D = find_1D+polynomial(x, a2, a1, a0)/2.
               ELSE IF (iinterp == iinterp_centred) THEN
                  ! In this case take the average of two separate quadratic spline values
                  find_1D = 0.
                  x0 = (x1+x2+x3)/3.
                  CALL fix_centred_polynomial(a2, a1, a0, x0, [x1, x2, x3], [y1, y2, y3])     
                  find_1D = find_1D+centred_polynomial(x, x0, a2, a1, a0)/2.
                  x0 = (x2+x3+x4)/3.
                  CALL fix_centred_polynomial(a2, a1, a0, x0, [x2, x3, x4], [y2, y3, y4])
                  find_1D = find_1D+centred_polynomial(x, x0, a2, a1, a0)/2.
               ELSE IF (iinterp == iinterp_Lagrange) THEN
                  ! In this case take the average of two quadratic Lagrange polynomials
                  L1 = Lagrange_polynomial(x, [x1, x2, x3], [y1, y2, y3])
                  L2 = Lagrange_polynomial(x, [x2, x3, x4], [y2, y3, y4])
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
               CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [y1, y2, y3, y4])
               find_1D = polynomial(x, a3, a2, a1, a0)
            ELSE IF (iinterp == iinterp_centred) THEN
               x0 = (x1+x2+x3+x4)/4.
               CALL fix_centred_polynomial(a3, a2, a1, a0, x0, [x1, x2, x3, x4], [y1, y2, y3, y4])
               find_1D = centred_polynomial(x, x0, a3, a2, a1, a0)
            ELSE IF (iinterp == iinterp_Lagrange) THEN
               find_1D = Lagrange_polynomial(x, [x1, x2, x3, x4], [y1, y2, y3, y4])
            ELSE
               STOP 'FIND_1D: Error, method not specified correctly'
            END IF

         ELSE

            STOP 'FIND_1D: Error, interpolation order specified incorrectly'

         END IF

      END IF

   END FUNCTION find_1D

   REAL FUNCTION find_2D(x, xin, y, yin, fin, nx, ny, iorder, ifindx, ifindy, iinterp)

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
      INTEGER, INTENT(IN) :: ifindx
      INTEGER, INTENT(IN) :: ifindy
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
      LOGICAL, PARAMETER :: xycubic = .FALSE.

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

      !IF (xtab(1) > xtab(nx)) STOP 'FIND_2D: x array in wrong order'
      IF (xtab(1) > xtab(nx)) THEN
         CALL reverse_array(xtab)
         CALL reverse_array(ftab, 1)
      END IF

      !IF (ytab(1) > ytab(ny)) STOP 'FIND_2D: y array in wrong order'
      IF (ytab(1) > ytab(ny)) THEN
         CALL reverse_array(ytab)
         CALL reverse_array(ftab, 2)
      END IF

      IF(iorder == 0) THEN

         ix = find_table_integer(x, xtab, ifindx)
         IF(ix == 0) ix = 1
         iy = find_table_integer(y, ytab, ifindy)
         IF(iy == 0) iy = 1

         find_2D = fin(ix, iy)

      ELSE IF(iorder == 1) THEN

         IF (nx < 2) STOP 'FIND_2D: Not enough x points in your array for linear interpolation'
         IF (ny < 2) STOP 'FIND_2D: Not enough y points in your array for linear interpolation'

         !! Get the x,y values !!

         ! Get the integer coordinates in the x direction
         ix = find_table_integer(x, xtab, ifindx)
         IF(ix==0) THEN
            ix=1
         ELSE IF(ix==nx) THEN
            ix=nx-1
         END IF
         ix1 = ix
         ix2 = ix+1

         ! Get the x values at the corners
         x1 = xtab(ix1)
         x2 = xtab(ix2)

         ! Get the integer coordinates in the y direction
         iy = find_table_integer(y, ytab, ifindy)
         IF(iy==0) THEN
            iy=1
         ELSE IF(iy==ny) THEN
            iy=ny-1
         END IF
         iy1 = iy
         iy2 = iy+1

         ! Get the y values at the corners
         y1 = ytab(iy1)
         y2 = ytab(iy2)

         ! Interpolation is function values at corners weighted by opposite area
         ! TODO: Clever loop here? Maybe not, its nice and transparent without the loop
         V(1, 1) = (x2-x)*(y2-y)*ftab(ix1, iy1)
         V(1, 2) = (x2-x)*(y-y1)*ftab(ix1, iy2)
         V(2, 1) = (x-x1)*(y2-y)*ftab(ix2, iy1)
         V(2, 2) = (x-x1)*(y-y1)*ftab(ix2, iy2)

         ! Normalisation
         find_2D = sum(V)/((x2-x1)*(y2-y1))

      ELSE IF (iorder == 3) THEN

         ! No cubic extrapolation implemented if the desired point is outside x AND y array boundary corners
         IF ((x < xtab(1) .OR. x > xtab(nx)) .AND. (y > ytab(ny) .OR. y < ytab(1))) THEN
            WRITE (*, *) 'FIND_2D: array xmin:', xtab(1)
            WRITE (*, *) 'FIND_2D: array xmax:', xtab(nx)
            WRITE (*, *) 'FIND_2D: requested x:', x
            WRITE (*, *) 'FIND_2D: array ymin:', ytab(1)
            WRITE (*, *) 'FIND_2D: array ymax:', ytab(ny)
            WRITE (*, *) 'FIND_2D: requested y:', y
            STOP 'FIND_2D: Desired point is outside both x and y array range'
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
               j = find_table_integer(y, ytab, ifindy)
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

            CALL fix_polynomial(a3, a2, a1, a0, [y1, y2, y3, y4], [f11, f12, f13, f14])
            f10 = polynomial(y, a3, a2, a1, a0)

            CALL fix_polynomial(a3, a2, a1, a0, [y1, y2, y3, y4], [f21, f22, f23, f24])      
            f20 = polynomial(y, a3, a2, a1, a0)

            !!

            !! x interpolation

            CALL fix_polynomial(a1, a0, [x1, x2], [f10, f20])
            find_2D = polynomial(x, a1, a0)

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
               i = find_table_integer(x, xtab, ifindx)
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

            CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [f11, f21, f31, f41])
            f01 = polynomial(x, a3, a2, a1, a0)

            CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [f12, f22, f32, f42])     
            f02 = polynomial(x, a3, a2, a1, a0)

            ! y interpolation

            CALL fix_polynomial(a1, a0, [y1, y2], [f01, f02])
            find_2D = polynomial(y, a1, a0)

         ELSE

            ! Points exists within table boundardies (normal)

            IF (nx < 4) STOP 'FIND_2D: Not enough x points in your array for cubic interpolation'
            IF (ny < 4) STOP 'FIND_2D: Not enough y points in your array for cubic interpolation'

            IF (x <= xtab(4)) THEN
               i = 2
            ELSE IF (x >= xtab(nx-3)) THEN
               i = nx-2
            ELSE
               i = find_table_integer(x, xtab, ifindx)
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
               j = find_table_integer(y, ytab, ifindy)
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

            CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [f11, f21, f31, f41])
            f01 = polynomial(x, a3, a2, a1, a0)

            CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [f12, f22, f32, f42]) 
            f02 = polynomial(x, a3, a2, a1, a0)

            CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [f13, f23, f33, f43])        
            f03 = polynomial(x, a3, a2, a1, a0)

            CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [f14, f24, f34, f44])
            f04 = polynomial(x, a3, a2, a1, a0)
            
            CALL fix_polynomial(a3, a2, a1, a0, [y1, y2, y3, y4], [f01, f02, f03, f04])
            findy = polynomial(y, a3, a2, a1, a0)

            IF(xycubic) THEN
            
               ! y interpolation

               CALL fix_polynomial(a3, a2, a1, a0, [y1, y2, y3, y4], [f11, f12, f13, f14])
               f10 = polynomial(y, a3, a2, a1, a0)

               CALL fix_polynomial(a3, a2, a1, a0, [y1, y2, y3, y4], [f21, f22, f23, f24])
               f20 = polynomial(y, a3, a2, a1, a0)

               CALL fix_polynomial(a3, a2, a1, a0, [y1, y2, y3, y4], [f31, f32, f33, f34])
               f30 = polynomial(y, a3, a2, a1, a0)

               CALL fix_polynomial(a3, a2, a1, a0, [y1, y2, y3, y4], [f41, f42, f43, f44])
               f40 = polynomial(y, a3, a2, a1, a0)

               CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [f10, f20, f30, f40])
               findx = polynomial(x, a3, a2, a1, a0)

               ! Final result is an average over each direction

            ELSE

               findx = findy

            END IF

            find_2D = (findx+findy)/2.

         END IF

      ELSE

         WRITE(*, *) 'FIND_2D: Order:', iorder
         STOP 'FIND_2D: Order not supported'

      END IF

   END FUNCTION find_2D

   REAL FUNCTION find_3D(x, xin, y, yin, z, zin, fin, nx, ny, nz, iorder, ifindx, ifindy, ifindz, iinterp)

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
      INTEGER, INTENT(IN) :: ifindx
      INTEGER, INTENT(IN) :: ifindy
      INTEGER, INTENT(IN) :: ifindz
      INTEGER, INTENT(IN) :: iinterp
      REAL :: xx(3), x12(3, 2), Dx(3)
      INTEGER :: ix(3), ix12(3, 2), nnx(3)
      INTEGER :: i, j, k, ii, jj, kk, d, di
      REAL :: F(2, 2, 2)
      REAL :: V(2, 2, 2)

      ! If the value required is off the table edge then this will halt

      ! iorder = 0 => nearest-neighbour interpolation
      ! iorder = 1 => linear interpolation

      ! ifind = 1 => find x in xtab by crudely searching from x(1) to x(n)
      ! ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
      ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

      ! iinterp = 1 => Use polynomials for interpolation

      IF (iinterp .NE. 1) STOP 'FIND_3D: No Lagrange polynomials for you, only regular polynomails, iinterp=1'

      IF (xin(1) > xin(nx)) STOP 'FIND_3D: x array in wrong order'
      IF (yin(1) > yin(ny)) STOP 'FIND_3D: y array in wrong order'
      IF (zin(1) > zin(nz)) STOP 'FIND_3D: z array in wrong order'

      IF(iorder == 0) THEN

         DO d = 1, 3
            IF(d == 1) ix(d) = find_table_integer(x, xin, ifindx)
            IF(d == 2) ix(d) = find_table_integer(y, yin, ifindy)
            IF(d == 3) ix(d) = find_table_integer(z, zin, ifindz)
         END DO

         DO d = 1, 3
            IF(ix(d) == 0) ix(d)=1
         END DO

         find_3D = fin(ix(1), ix(2), ix(3))

      ELSE IF (iorder == 1) THEN

         xx(1) = x
         xx(2) = y
         xx(3) = z

         nnx(1) = nx
         nnx(2) = ny
         nnx(3) = nz

         DO d = 1, 3
            IF (nnx(d) < 2) STOP 'FIND_3D: Not enough points in your array for linear interpolation'
         END DO

         !! Get the x,y,z values !!

         ! Get the integer coordinates in each direction
         DO d = 1, 3

            ! Get the integer coordinates of the 'cube' corners that encompass the point in each dimensions
            IF(d == 1) ix(1) = find_table_integer(x, xin, ifindx)
            IF(d == 2) ix(2) = find_table_integer(y, yin, ifindy)
            IF(d == 3) ix(3) = find_table_integer(z, zin, ifindz)

            ! Correct for the case of these integers coordinates being at the edge of the box
            IF(ix(d) == 0) THEN
               ix(d) = 1
            ELSE IF(ix(d) == nnx(d)) THEN
               ix(d) = nnx(d)-1
            END IF
            ix12(d, 1) = ix(d)
            ix12(d, 2) = ix(d)+1

            ! Get the function value at the corners of the 'cube'
            IF(d == 1) THEN
               x12(d, 1) = xin(ix12(d, 1))
               x12(d, 2) = xin(ix12(d, 2))
            ELSE IF(d == 2) THEN
               x12(d, 1) = yin(ix12(d, 1))
               x12(d, 2) = yin(ix12(d, 2))
            ELSE IF(d == 3) THEN
               x12(d, 1) = zin(ix12(d, 1))
               x12(d, 2) = zin(ix12(d, 2))
            END IF

         END DO 

         ! Do the interpolation
         DO k = 1, 2
            DO j = 1, 2
               DO i = 1, 2

                  ! Evaluate the function at the cube corners
                  F(i, j, k) = fin(ix12(1, i), ix12(2, j), ix12(3, k))

                  ! Calculate the volume of the solid in the opposite corner
                  ii = 3-i ! Swaps 1 <--> 2 (e.g., ii is 2 if i is 1 and visa versa)
                  jj = 3-j ! Swaps 1 <--> 2
                  kk = 3-k ! Swaps 1 <--> 2
                  DO d = 1, 3
                     IF(d == 1) THEN
                        di = ii
                     ELSE IF(d == 2) THEN
                        di = jj
                     ELSE IF(d == 3) THEN
                        di = kk
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

   SUBROUTINE interpolate_array(x_old, y_old, x_new, y_new, iorder, ifind, iinterp, logx, logy)

      ! Interpolates array 'x1-y1' onto new 'x' values x2 and output y2
      ! TODO: This could be more efficient because it currently does 'find integer' every time
      REAL, INTENT(IN) :: x_old(:)
      REAL, INTENT(IN) :: y_old(:)
      REAL, INTENT(IN) :: x_new(:)
      REAL, INTENT(OUT) :: y_new(:)
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      LOGICAL, OPTIONAL, INTENT(IN) :: logx
      LOGICAL, OPTIONAL, INTENT(IN) :: logy
      REAL, ALLOCATABLE :: xx_old(:), yy_old(:)
      REAL :: xx
      INTEGER :: i, n_old, n_new

      n_old = size(x_old)
      IF (n_old /= size(y_old)) STOP 'INTERPOLATE_ARRAY: error, x_old and y_old must be the same size'

      n_new = size(x_new)
      IF (n_new /= size(y_new)) STOP 'INTERPOLATE_ARRAY: error, x_new and y_new must be the same size'

      IF (present_and_correct(logx)) THEN
         xx_old = log(x_old)
      ELSE
         xx_old = x_old
      END IF

      IF (present_and_correct(logy)) THEN
         yy_old = log(y_old)
      ELSE
         yy_old = y_old
      END IF
    
      DO i = 1, n_new
         IF (present_and_correct(logx)) THEN
            xx = log(x_new(i))
         ELSE
            xx = x_new(i)
         END IF
         y_new(i) = find(xx, xx_old, yy_old, n_old, iorder, ifind, iinterp)
      END DO

      IF (present_and_correct(logy)) THEN
         y_new = exp(y_new)
      END IF

   END SUBROUTINE interpolate_array

   SUBROUTINE rebin_array(xmin_new, xmax_new, n_new, x, f, iorder, ifind, iinterp, logx, logf)

      ! TODO: Remove. Interpolate array should be sufficient
      REAL, INTENT(IN) :: xmin_new
      REAL, INTENT(IN) :: xmax_new
      INTEGER, INTENT(IN) :: n_new
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: f(:)
      INTEGER, INTENT(IN) :: iorder
      INTEGER, INTENT(IN) :: ifind
      INTEGER, INTENT(IN) :: iinterp
      LOGICAL, OPTIONAL, INTENT(IN) :: logx
      LOGICAL, OPTIONAL, INTENT(IN) :: logf
      REAL, ALLOCATABLE :: x_old(:), f_old(:)
      INTEGER :: n_old

      n_old = size(x)
      IF (n_old /= size(f, 1)) STOP 'REBIN_ARRAY: Error, x and f must be the same size'

      x_old = x
      f_old = f

      IF (present_and_correct(logx)) THEN
         CALL fill_array_log(xmin_new, xmax_new, x, n_new)
      ELSE
         CALL fill_array(xmin_new, xmax_new, x, n_new)
      END IF
      DEALLOCATE(f)
      ALLOCATE(f(n_new))

      CALL interpolate_array(x_old, f_old, x, f, iorder, ifind, iinterp, logx, logf)

   END SUBROUTINE rebin_array

   SUBROUTINE init_interpolator_1D(x, f, interp, iorder, iextrap, store, logx, logf)

      ! Initialise an interpolator
      ! TODO: Should be loops and arrays rather than x1, x2, x3, x4 etc.
      REAL, INTENT(IN) :: x(:)                  ! Input data x
      REAL, INTENT(IN) :: f(:)                  ! Input data f(x)
      TYPE(interpolator1D), INTENT(OUT) :: interp ! Interpolator type
      INTEGER, INTENT(IN) :: iorder             ! Order at which to create interpolator   
      INTEGER, INTENT(IN) :: iextrap            ! Extrapolation scheme
      LOGICAL, INTENT(IN) :: store              ! Do we store values
      LOGICAL, OPTIONAL, INTENT(IN) :: logx     ! Should interpolator take the logarithm of x?
      LOGICAL, OPTIONAL, INTENT(IN) :: logf     ! Should interpolator take the logarithm of f(x)?
      REAL, ALLOCATABLE :: xx(:), ff(:)
      REAL :: a0, a1, a2, a3
      REAL :: b0, b1, b2
      REAL :: f1, f2, f3, f4
      REAL :: x1, x2, x3, x4
      REAL :: xm
      INTEGER :: i, n, ii, nn
      INTEGER :: i1, i2, i3, i4
      INTEGER, PARAMETER :: ifind_default = ifind_interpolator_default

      ! Should the first and last sections of cubic interpolation be quadratic
      LOGICAL, PARAMETER :: quadratic_end_cubic = .FALSE. 

      ! Deallocate arrays if they are already allocated
      CALL if_allocated_deallocate(interp%x)
      CALL if_allocated_deallocate(interp%x0)
      CALL if_allocated_deallocate(interp%a0)
      CALL if_allocated_deallocate(interp%a1)
      CALL if_allocated_deallocate(interp%a2)
      CALL if_allocated_deallocate(interp%a3)

      ! Check the sizes of the input data arrays
      n = size(x)
      IF (n /= size(f)) STOP 'INIT_INTERPOLATOR_1D: Error, input x and f data should be the same size'
      interp%n = n

      IF (n == 1) THEN
         STOP 'INIT_INTERPOLATOR_1D: Error, array is size 1'
      ELSE IF (n < iorder+1) THEN
         STOP 'INIT_INTERPOLATOR_1D: Error, not enough points for your order of interpolation'
      ELSE

         ! Set the internal variables
         interp%iextrap = iextrap
         interp%iorder = iorder

         ! Fill the internal arrays
         xx = x
         ff = f

         ! Reverse x and f if necessary
         IF (xx(1) > xx(n)) THEN
            CALL reverse_array(xx)
            CALL reverse_array(ff)
         END IF

         ! Set the internal xmin and xmax
         interp%xmin = xx(1)
         interp%xmax = xx(n)

         ! Sort out logs and x arrays
         IF (present_and_correct(logx)) THEN
            xx = log(xx)
            interp%logx = .TRUE.
         ELSE
            interp%logx = .FALSE.
         END IF

         ! Sort out logs and f array
         IF (present_and_correct(logf)) THEN
            ff = log(ff)
            interp%logf = .TRUE.
         ELSE
            interp%logf = .FALSE.
         END IF
         
         ! If the x data is regular spaced then remember this, otherwise default find
         IF (regular_spacing(xx)) THEN
            interp%ifind = ifind_linear
         ELSE
            interp%ifind = ifind_default
         END IF
         interp%store = store

         ! Set the internal 
         interp%x = xx
         interp%f = ff

         IF (interp%store) THEN

            ! Allocate arrays for the interpolator coefficients
            IF (iextrap == iextrap_lin) THEN
               nn = n+1
            ELSE
               nn = n-1
            END IF
            ALLOCATE(interp%a0(nn), interp%a1(nn))
            interp%a0 = 0.
            interp%a1 = 0.
            IF (iorder >= 2) THEN
               ALLOCATE(interp%a2(nn))
               interp%a2 = 0.
            END IF
            IF (iorder >= 3) THEN
               ALLOCATE(interp%a3(nn))
               interp%a3 = 0.
            END IF
            ALLOCATE(interp%x0(nn))

            ! Default values because a3, a2 will be zero for linear polynomials
            a0 = 0.
            a1 = 0.
            a2 = 0.
            a3 = 0.

            ! Default values
            x1 = 0.
            x2 = 0.
            x3 = 0.
            x4 = 0.
            
            ! Default values
            f1 = 0.
            f2 = 0.
            f3 = 0.
            f4 = 0.

            ! Loop over all n-1 sections of the input data
            DO i = 0, n

               IF ((i == 0 .OR. i == n) .AND. iextrap /= iextrap_lin) CYCLE

               IF (iorder == 1 .OR. ((i == 0 .OR. i == n) .AND. iextrap == iextrap_lin)) THEN

                  IF (i == 0) THEN
                     i1 = 1
                     i2 = 2
                  ELSE IF (i == n) THEN
                     i1 = n-1
                     i2 = n
                  ELSE
                     i1 = i
                     i2 = i+1
                  END IF

                  x1 = xx(i1)
                  x2 = xx(i2)

                  f1 = ff(i1)
                  f2 = ff(i2)
                  
                  xm = (x1+x2)/2.
                  CALL fix_centred_polynomial(a1, a0, xm, [x1, x2], [f1, f2])

               ELSE IF (iorder == 2) THEN

                  i1 = i-1
                  i2 = i
                  i3 = i+1
                  i4 = i+2
                  IF (i == 1) THEN
                     i1 = i1+1
                     i2 = i2+1
                     i3 = i3+1
                     i4 = 0
                  ELSE IF (i == n-1) THEN
                     i1 = i1-1
                     i2 = i2-1
                     i3 = i3-1
                     i4 = 0
                  END IF

                  x1 = xx(i1)
                  x2 = xx(i2)
                  x3 = xx(i3)
                  IF (i4 .NE. 0) x4 = xx(i4)

                  f1 = ff(i1)
                  f2 = ff(i2)
                  f3 = ff(i3)
                  IF (i4 .NE. 0) f4 = ff(i4)

                  CALL fix_polynomial(a2, a1, a0, [x1, x2, x3], [f1, f2, f3])
                  IF (i4 .NE. 0) THEN
                     CALL fix_polynomial(b2, b1, b0, [x2, x3, x4], [f2, f3, f4])
                     a2 = (a2+b2)/2.
                     a1 = (a1+b1)/2.
                     a0 = (a0+b0)/2.
                  END IF

               ELSE IF (iorder == 3) THEN

                  ! Deal with the indices for the first and last section and general case
                  i1 = i-1
                  i2 = i
                  i3 = i+1
                  i4 = i+2
                  IF (i == 1) THEN
                     ! Should be 1, 2, 3, 4
                     i1 = i1+1
                     i2 = i2+1 
                     i3 = i3+1 
                     i4 = i4+1 
                  ELSE IF (i == n-1) THEN
                     ! Should be n-3, n-2, n-1, n
                     i1 = i1-1
                     i2 = i2-1
                     i3 = i3-1
                     i4 = i4-1
                  END IF

                  x1 = xx(i1)
                  x2 = xx(i2)
                  x3 = xx(i3)
                  x4 = xx(i4)

                  f1 = ff(i1)
                  f2 = ff(i2)
                  f3 = ff(i3)
                  f4 = ff(i4)

                  xm = (xx(i)+xx(i+1))/2.
                  IF (quadratic_end_cubic .AND. i == 1) THEN
                     CALL fix_centred_polynomial(a2, a1, a0, xm, [x1, x2, x3], [f1, f2, f3]) 
                     a3 = 0.
                  ELSE IF (quadratic_end_cubic .AND. i == n-1) THEN
                     CALL fix_centred_polynomial(a2, a1, a0, xm, [x2, x3, x4], [f2, f3, f4]) 
                     a3 = 0.
                  ELSE
                     CALL fix_centred_polynomial(a3, a2, a1, a0, xm, [x1, x2, x3, x4], [f1, f2, f3, f4])  
                  END IF    

               ELSE

                  STOP 'INIT_INTERPOLATOR: Error, your order is not supported'

               END IF

               ! Fill the polynomial coefficients in the interpolation type
               IF (iextrap == iextrap_lin) THEN
                  ii = i+1
               ELSE
                  ii = i
               END IF
               interp%x0(ii) = xm
               interp%a0(ii) = a0
               interp%a1(ii) = a1
               IF(iorder >= 2 .AND. (i .NE. 0 .AND. i .NE. n)) interp%a2(ii) = a2
               IF(iorder >= 3 .AND. (i .NE. 0 .AND. i .NE. n)) interp%a3(ii) = a3

            END DO

         END IF

      END IF

   END SUBROUTINE init_interpolator_1D

   REAL FUNCTION evaluate_interpolator_1D(x, interp)

      ! Evaluates the value f(x) from the interpolator
      REAL, INTENT(IN) :: x
      TYPE(interpolator1D), INTENT(IN) :: interp
      INTEGER :: i, n
      REAL :: xx

      ! Change input x to log if necessary
      xx = x
      IF (interp%logx) xx = log(xx)

      n = interp%n

      IF (interp%store) THEN

         ! Calculate the index to use
         IF (xx < interp%x(1)) THEN
            IF (interp%iextrap == iextrap_no) THEN
               STOP 'EVALUATE_INTERPOLATOR_1D: Error, desired x value is below range'
            ELSE IF(interp%iextrap == iextrap_std) THEN
               i = 1
            ELSE IF (interp%iextrap == iextrap_lin) THEN
               i = 0
            ELSE IF (interp%iextrap == iextrap_zero) THEN
               evaluate_interpolator_1D = 0.
               RETURN
            ELSE IF (interp%iextrap == iextrap_near) THEN
               evaluate_interpolator_1D = interp%f(1)
               IF (interp%logf) evaluate_interpolator_1D = exp(evaluate_interpolator_1D)
               RETURN
            ELSE
               STOP 'EVALUATE_INTERPOLATOR_1D: Error, extrapolation scheme specified incorrectly'
            END IF
         ELSE IF (xx == interp%x(n)) THEN
            i = n-1 ! Cannot use find_table_integer here because this would return 'n', rather than 'n-1'
         ELSE IF (xx > interp%x(n)) THEN
            IF (interp%iextrap == iextrap_no) THEN
               STOP 'EVALUATE_INTERPOLATOR_1D: Error, desired x value is above range'
            ELSE IF (interp%iextrap == iextrap_std) THEN
               i = n-1
            ELSE IF(interp%iextrap == iextrap_lin) THEN
               i = n
            ELSE IF (interp%iextrap == iextrap_zero) THEN
               evaluate_interpolator_1D = 0.
               RETURN
            ELSE IF (interp%iextrap == iextrap_near) THEN
               evaluate_interpolator_1D = interp%f(n)
               RETURN
            ELSE
               STOP 'EVALUATE_INTERPOLATOR_1D: Error, extrapolation scheme specified incorrectly'
            END IF
         ELSE
            i = find_table_integer(xx, interp%x, interp%ifind)
         END IF

         IF (interp%iextrap == iextrap_lin) i = i+1

         ! Evaluate the interpolation
         IF (interp%iorder == 1) THEN
            evaluate_interpolator_1D = centred_polynomial(xx, interp%x0(i), interp%a1(i), interp%a0(i))
         ELSE IF (interp%iorder == 2) THEN
            evaluate_interpolator_1D = centred_polynomial(xx, interp%x0(i), interp%a2(i), interp%a1(i), interp%a0(i))
         ELSE IF (interp%iorder == 3) THEN
            evaluate_interpolator_1D = centred_polynomial(xx, interp%x0(i), interp%a3(i), interp%a2(i), interp%a1(i), interp%a0(i))
         ELSE
            STOP 'EVALUATE_INTERPOLATOR_1D: Error, your polynomial order is not supported'
         END IF

      ELSE
         evaluate_interpolator_1D = find(xx, interp%x, interp%f, interp%n, interp%iorder, interp%ifind, iinterp_Lagrange)
      END IF

      ! Exponentiate result if necessary
      IF (interp%logf) evaluate_interpolator_1D = exp(evaluate_interpolator_1D)

   END FUNCTION evaluate_interpolator_1D

   REAL FUNCTION inverse_interpolator_1D(f, interp)

      ! If y = f(x) this returns x = f^-1(y), where y = f
      REAL, INTENT(IN) :: f
      TYPE(interpolator1D), INTENT(IN) :: interp
      REAL :: ff, x
      INTEGER :: ifind = ifind_inverse_interpolator
      INTEGER :: iinterp = iinterp_inverse_interpolator

      IF (interp%logf) ff = log(f)
      x = find(ff, interp%f, interp%x, interp%n, interp%iorder, ifind, iinterp)

      IF(interp%logx) THEN
         inverse_interpolator_1D = exp(x)
      ELSE
         inverse_interpolator_1D = x
      END IF

   END FUNCTION inverse_interpolator_1D

   SUBROUTINE init_interpolator_2D(x, y, f, interp, iorder, iextrap, store, logx, logy, logf)

      ! Initialise a 2D interpolator
      ! TODO: Should end sections be quadratic when doing cubic interpolation?
      ! TODO: Linear extrapolation unnecessarily fills entire centre of array with linear interpolation
      ! TODO: Linear extrapolation only needs to do rows ix = 1, nx-1, nx and same for y
      REAL, INTENT(IN) :: x(:)                    ! Input data x
      REAL, INTENT(IN) :: y(:)                    ! Input data y
      REAL, INTENT(IN) :: f(:, :)                 ! Input data f(x, y)
      TYPE(interpolator2D), INTENT(OUT) :: interp ! Interpolator type
      INTEGER, INTENT(IN) :: iorder               ! Order at which to create interpolator     
      INTEGER, INTENT(IN) :: iextrap              ! Should interpolator extrapolate beyond x and y?
      LOGICAL, INTENT(IN) :: store                ! Should we store polynomial coefficients?
      LOGICAL, OPTIONAL, INTENT(IN) :: logx       ! Should interpolator take the logarithm of x?
      LOGICAL, OPTIONAL, INTENT(IN) :: logy       ! Should interpolator take the logarithm of y?
      LOGICAL, OPTIONAL, INTENT(IN) :: logf       ! Should interpolator take the logarithm of f?
      INTEGER :: ix, iy, nx, ny
      INTEGER :: jx(4), jy(4)
      INTEGER :: i
      REAL :: xx(4), yy(4), fx(4), fy(4)
      REAL :: x0, y0
      REAL :: a0, a1, a2, a3
      INTEGER, PARAMETER :: ifind_default = ifind_interpolator_default

      CALL if_allocated_deallocate(interp%x)
      CALL if_allocated_deallocate(interp%y)
      CALL if_allocated_deallocate(interp%f)
      CALL if_allocated_deallocate(interp%ax0)
      CALL if_allocated_deallocate(interp%ax1)
      CALL if_allocated_deallocate(interp%ax2)
      CALL if_allocated_deallocate(interp%ax3)
      CALL if_allocated_deallocate(interp%ay0)
      CALL if_allocated_deallocate(interp%ay1)
      CALL if_allocated_deallocate(interp%ay2)
      CALL if_allocated_deallocate(interp%ay3)

      IF (iorder > 3 .OR. iorder < 0) STOP 'INIT_INTERPOLATOR_2D: Error, this order is not supported'

      ! Check array sizes
      nx = size(x)
      IF(nx /= size(f, 1)) STOP 'INIT_INTERPOLATOR_2D: Error, x should be the same size as first dimension of f'
      interp%nx = nx
      ny = size(y)
      IF(ny /= size(f, 2)) STOP 'INIT_INTERPOLATOR_2D: Error, y should be the same size as second dimension of f'
      interp%ny = ny

      IF (nx == 1 .AND. ny == 1) THEN
         STOP 'INIT_INTERPOLATOR_2D: Error, your array is 1 x 1'
      ELSE IF (nx == 1) THEN
         CALL init_interpolator_1D(y, f(1, :), interp%interp1D, iorder, iextrap, store, logy, logf)
      ELSE IF (ny == 1) THEN
         CALL init_interpolator_1D(x, f(:, 1), interp%interp1D, iorder, iextrap, store, logx, logf)
      ELSE IF (nx < iorder+1 .OR. ny < iorder+1) THEN
         STOP 'INIT_INTERPOLATOR_2D: Error, you do not have enough points for the order of interpolatation'
      ELSE

         ! Sort out x axis      
         IF(present_and_correct(logx)) THEN
            interp%x = log(x)
            interp%logx = .TRUE.
         ELSE
            interp%x = x
            interp%logx = .FALSE.
         END IF
         IF (regular_spacing(interp%x)) THEN
            interp%ifindx = ifind_linear
         ELSE
            interp%ifindx = ifind_default
         END IF

         ! Sort out y axis     
         IF(present_and_correct(logy)) THEN
            interp%y = log(y)
            interp%logy = .TRUE.
         ELSE
            interp%y = y
            interp%logy = .FALSE.
         END IF
         IF (regular_spacing(interp%y)) THEN
            interp%ifindy = ifind_linear
         ELSE
            interp%ifindy = ifind_default
         END IF

         ! Sort out f
         IF(present_and_correct(logf)) THEN
            interp%f = log(f)
            interp%logf = .TRUE.
         ELSE
            interp%f = f
            interp%logf = .FALSE.
         END IF

         ! Reverse arrays if necessary
         IF (interp%x(1) > interp%x(nx)) THEN
            CALL reverse_array(interp%x)
            CALL reverse_array(interp%f, 1)
         END IF
         IF (interp%y(1) > interp%y(ny)) THEN
            CALL reverse_array(interp%y)
            CALL reverse_array(interp%f, 2)
         END IF

         ! Set the interp%xmin, xmax, ymin, ymax values
         interp%xmin = interp%x(1)
         interp%xmax = interp%x(nx)
         interp%ymin = interp%y(1)
         interp%ymax = interp%y(ny)
         IF (interp%logx) THEN
            interp%xmin = exp(interp%xmin)
            interp%xmax = exp(interp%xmax)
         END IF
         IF (interp%logy) THEN
            interp%ymin = exp(interp%ymin)
            interp%ymax = exp(interp%ymax)
         END IF

         ! Set internal variables
         interp%iorder = iorder
         interp%iextrap = iextrap
         interp%store = store

         IF (interp%store) THEN 

            ALLOCATE(interp%x0(nx-1, ny), interp%y0(nx, ny-1))
            ALLOCATE(interp%ax0(nx-1, ny), interp%ax1(nx-1, ny), interp%ax2(nx-1, ny), interp%ax3(nx-1, ny))
            ALLOCATE(interp%ay0(nx, ny-1), interp%ay1(nx, ny-1), interp%ay2(nx, ny-1), interp%ay3(nx, ny-1))
            interp%x0 = 0.
            interp%ax0 = 0.
            interp%ax1 = 0.
            interp%ax2 = 0.
            interp%ax3 = 0.
            interp%y0 = 0.
            interp%ay0 = 0.
            interp%ay1 = 0.
            interp%ay2 = 0.
            interp%ay3 = 0.
            
            IF ((interp%iorder == 1) .OR. (iextrap == iextrap_lin)) THEN

               DO iy = 1, ny

                  ! TODO: This should work, but it does not
                  !IF (iorder > 1 .AND. iy > 2 .AND. iy < ny-1) CYCLE

                  DO ix = 1, nx

                     ! TODO: This should work, but it does not
                     !IF (iorder > 1 .AND. ix > 2 .AND. ix < nx-1) CYCLE

                     IF (ix /= nx) THEN

                        ! Indices run sequentially but must be within 1, nx
                        DO i = 1, 2
                           jx(i) = ix+i-1
                        END DO

                        ! Get x values and function values running along x direction
                        DO i = 1, 2
                           xx(i) = interp%x(jx(i))
                           fx(i) = interp%f(jx(i), iy)
                        END DO

                        ! Calculate the mid point for polynomial
                        x0 = (xx(1)+xx(2))/2.
                        interp%x0(ix, iy) = x0

                        ! Fix polynomial along the x direction
                        CALL fix_centred_polynomial(a1, a0, x0, xx, fx)
                        interp%ax1(ix, iy) = a1
                        interp%ax0(ix, iy) = a0

                     END IF

                     IF (iy /= ny) THEN

                        ! Indices run sequentially but must be within 1, ny
                        DO i = 1, 2
                           jy(i) = iy+i-1
                        END DO

                        ! Get y values and function values running along y direction
                        DO i = 1, 2
                           yy(i) = interp%y(jy(i))
                           fy(i) = interp%f(ix, jy(i))
                        END DO

                        ! Calculate the mid point for polynomial
                        y0 = (yy(1)+yy(2))/2.
                        interp%y0(ix, iy) = y0

                        ! Fix polynomial along the y direction
                        CALL fix_centred_polynomial(a1, a0, y0, yy, fy)
                        interp%ay1(ix, iy) = a1
                        interp%ay0(ix, iy) = a0

                     END IF

                  END DO
               END DO

               IF ((iextrap == iextrap_lin) .AND. (iorder > 1)) THEN

                  interp%bx0 = interp%ax0
                  interp%by0 = interp%ay0
                  interp%bx1 = interp%ax1
                  interp%by1 = interp%ay1
                  interp%xl0 = interp%x0
                  interp%yl0 = interp%y0

               END IF

            END IF

            IF (interp%iorder == 2 .OR. interp%iorder == 3) THEN

               DO iy = 1, ny
                  DO ix = 1, nx

                     IF (ix /= nx) THEN

                        ! Indices run sequentially but must be within 1, nx
                        DO i = 1, 4
                           jx(i) = ix+(i-2)
                        END DO
                        IF (ix == 1)    jx = jx+1
                        IF (ix == nx-1) jx = jx-1

                        ! Get x values and function values running along x direction
                        DO i = 1, 4
                           xx(i) = interp%x(jx(i))
                           fx(i) = interp%f(jx(i), iy)
                        END DO

                        ! Calculate the mid point for polynomial
                        x0 = (xx(2)+xx(3))/2.
                        interp%x0(ix, iy) = x0

                        ! Fix polynomial along the x direction
                        CALL fix_centred_polynomial(a3, a2, a1, a0, x0, xx, fx)
                        interp%ax3(ix, iy) = a3
                        interp%ax2(ix, iy) = a2
                        interp%ax1(ix, iy) = a1
                        interp%ax0(ix, iy) = a0

                     END IF

                     IF (iy /= ny) THEN

                        ! Indices run sequentially but must be within 1, ny
                        DO i = 1, 4
                           jy(i) = iy+(i-2)
                        END DO
                        IF (iy == 1)    jy = jy+1
                        IF (iy == ny-1) jy = jy-1

                        ! Get y values and function values running along y direction
                        DO i = 1, 4
                           yy(i) = interp%y(jy(i))
                           fy(i) = interp%f(ix, jy(i))
                        END DO

                        ! Calculate the mid point for polynomial
                        y0 = (yy(2)+yy(3))/2.
                        interp%y0(ix, iy) = y0

                        ! Fix polynomial along the y direction
                        CALL fix_centred_polynomial(a3, a2, a1, a0, y0, yy, fy)
                        interp%ay3(ix, iy) = a3
                        interp%ay2(ix, iy) = a2
                        interp%ay1(ix, iy) = a1
                        interp%ay0(ix, iy) = a0

                     END IF

                  END DO
               END DO

            END IF
            
         END IF

      END IF

   END SUBROUTINE init_interpolator_2D

   REAL FUNCTION evaluate_interpolator_2D(x, y, interp)

      ! TODO: Quadratics at end sections?
      ! TODO: Should linear interpolation turn of xycubic?
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      TYPE(interpolator2D), INTENT(IN) :: interp
      INTEGER, PARAMETER :: iinterp = iinterp_polynomial ! No Lagrange polynomials in 2D
      INTEGER :: ix, iy, nx, ny, i, j, mx, my
      REAL :: xx, yy, ffx, ffy
      REAL :: a3, a2, a1, a0, x0, y0
      INTEGER, ALLOCATABLE :: jx(:), jy(:)
      REAL, ALLOCATABLE :: fx(:), fy(:), xxx(:), yyy(:)
      LOGICAL :: outx, outy
      LOGICAL, PARAMETER :: xycubic = .FALSE.

      xx = x
      IF (interp%logx) xx = log(xx)
      nx = interp%nx

      yy = y
      IF (interp%logy) yy = log(yy)
      ny = interp%ny

      IF (nx == 1) THEN

         ! Function is actually 1D in y, so evaluate this
         evaluate_interpolator_2D = evaluate_interpolator_1D(y, interp%interp1D)

      ELSE IF (ny == 1) THEN

         ! Function is actually 1D in x, so evaluate this
         evaluate_interpolator_2D = evaluate_interpolator_1D(x, interp%interp1D)

      ELSE IF (interp%store) THEN

         IF (xx < interp%x(1) .OR. xx > interp%x(nx)) THEN
            outx = .TRUE.
         ELSE
            outx = .FALSE.
            END IF
         IF (yy < interp%y(1) .OR. yy > interp%y(ny)) THEN
            outy = .TRUE.
         ELSE
            outy = .FALSE.
         END IF

         IF (outx .OR. outy) THEN
            IF (interp%iextrap == iextrap_no) THEN
               STOP 'EVALUATE_INTERPOLATOR_2D: Error, point is outside interpolator range'
            ELSE IF (interp%iextrap == iextrap_zero) THEN
               evaluate_interpolator_2D = 0.
               RETURN
            ELSE IF (interp%iextrap == iextrap_near) THEN
               ix = find_table_integer(xx, interp%x, interp%ifindx)
               iy = find_table_integer(yy, interp%y, interp%ifindy)
               evaluate_interpolator_2D = interp%f(ix, iy)
               IF(interp%logf) evaluate_interpolator_2D = exp(evaluate_interpolator_2D)
               RETURN
            END IF
         END IF

         IF (interp%iorder == 1 .OR. interp%iorder == 2) THEN
            mx = 2
            my = 2
         ELSE IF (interp%iorder == 3) THEN
            IF (outx) THEN
               mx = 2
            ELSE
               mx = 4
            END IF
            IF (outy) THEN
               my = 2
            ELSE
               my = 4
            END IF
         ELSE
            STOP 'EVALUATE_INTERPOLATOR_2D: Error, something went wrong with iorder'
         END IF
         ALLOCATE(jx(mx), fx(mx), xxx(mx))
         ALLOCATE(jy(my), fy(my), yyy(my))

         ix = find_table_integer(xx, interp%x, interp%ifindx)
         IF (ix == 0)  ix = 1
         IF (ix == nx) ix = nx-1
         IF (mx == 2) THEN
            jx(1) = ix
            jx(2) = ix+1
         ELSE IF (mx == 4) THEN
            DO i = 1, mx
               jx(i) = ix+(i-2)
            END DO
            IF (ix == 1)    jx = jx+1
            IF (ix == nx-1) jx = jx-1  
         ELSE
            STOP 'EVALUATE_INTERPOLATOR_2D: Error, something went wrong with iorder'
         END IF

         iy = find_table_integer(yy, interp%y, interp%ifindy)
         IF (iy == 0)  iy = 1
         IF (iy == ny) iy = ny-1
         IF (my == 2) THEN
            jy(1) = iy
            jy(2) = iy+1
         ELSE IF (my == 4) THEN 
            DO i = 1, my
               jy(i) = iy+(i-2)
            END DO 
            IF (iy == 1)    jy = jy+1
            IF (iy == ny-1) jy = jy-1
         ELSE
            STOP 'EVALUATE_INTERPOLATOR_2D: Error, something went wrong with iorder'
         END IF
         
         DO i = 1, mx
            j = jx(i)
            xxx(i) = interp%x(j)
            IF (interp%iorder .NE. 1 .AND. interp%iextrap == iextrap_lin .AND. outy) THEN
               y0 = interp%yl0(j, iy)
               a3 = 0.
               a2 = 0.
               a1 = interp%by1(j, iy)
               a0 = interp%by0(j, iy)
            ELSE
               y0 = interp%y0(j, iy)
               a3 = interp%ay3(j, iy)
               a2 = interp%ay2(j, iy)
               a1 = interp%ay1(j, iy)
               a0 = interp%ay0(j, iy)
            END IF
            fx(i) = centred_polynomial(yy, y0, a3, a2, a1, a0)
         END DO
         ffy = Lagrange_polynomial(xx, xxx, fx)

         IF (xycubic) THEN
            DO i = 1, my
               j = jy(i)
               yyy(i) = interp%y(j)
               IF (interp%iorder .NE. 1 .AND. interp%iextrap == iextrap_lin .AND. outx) THEN
                  x0 = interp%xl0(ix, j)
                  a3 = 0.
                  a2 = 0.
                  a1 = interp%bx1(ix, j)
                  a0 = interp%bx0(ix, j)
               ELSE
                  x0 = interp%x0(ix, j)
                  a3 = interp%ax3(ix, j)
                  a2 = interp%ax2(ix, j)
                  a1 = interp%ax1(ix, j)
                  a0 = interp%ax0(ix, j)
               END IF
               fy(i) = centred_polynomial(xx, x0, a3, a2, a1, a0)                  
            END DO
            ffx = Lagrange_polynomial(yy, yyy, fy)
         ELSE
            ffx = ffy               
         END IF

         evaluate_interpolator_2D = (ffx+ffy)/2.
         IF (interp%logf) evaluate_interpolator_2D = exp(evaluate_interpolator_2D)

      ELSE

         evaluate_interpolator_2D = find(xx, interp%x, yy, interp%y, interp%f, interp%nx, interp%ny, &
                     interp%iorder, &
                     interp%ifindx, &
                     interp%ifindy, &
                     iinterp) 
         IF (interp%logf) evaluate_interpolator_2D = exp(evaluate_interpolator_2D)

      END IF

   END FUNCTION evaluate_interpolator_2D

   SUBROUTINE init_interpolator_3D(x, y, z, f, interp, iorder, iextrap, store, logx, logy, logz, logf)

      ! Initialise a 3D interpolator
      ! TODO: Checks on array sizes
      ! TODO: Remove unnecessary ALLOCATE statements
      REAL, INTENT(IN) :: x(:)                    ! Input data x
      REAL, INTENT(IN) :: y(:)                    ! Input data y
      REAL, INTENT(IN) :: z(:)                    ! Input data z
      REAL, INTENT(IN) :: f(:, :, :)              ! Input data f(x, y, z)
      TYPE(interpolator3D), INTENT(OUT) :: interp ! Interpolator type
      INTEGER, INTENT(IN) :: iorder               ! Order at which to create interpolator     
      INTEGER, INTENT(IN) :: iextrap              ! Should interpolator extrapolate beyond x and y?
      LOGICAL, INTENT(IN) :: store                ! Should we store polynomial coefficients?
      LOGICAL, OPTIONAL, INTENT(IN) :: logx       ! Should interpolator take the logarithm of x?
      LOGICAL, OPTIONAL, INTENT(IN) :: logy       ! Should interpolator take the logarithm of y?
      LOGICAL, OPTIONAL, INTENT(IN) :: logz       ! Should interpolator take the logarithm of z?
      LOGICAL, OPTIONAL, INTENT(IN) :: logf       ! Should interpolator take the logarithm of f?
      INTEGER :: nx, ny, nz
      INTEGER, PARAMETER :: ifind_default = ifind_interpolator_default

      CALL if_allocated_deallocate(interp%x)
      CALL if_allocated_deallocate(interp%y)
      CALL if_allocated_deallocate(interp%z)
      CALL if_allocated_deallocate(interp%f)

      IF (iorder > 3 .OR. iorder < 0) STOP 'INIT_INTERPOLATOR_2D: Error, this order is not supported'

      ! Sort out x axis
      nx = size(x)
      IF(nx /= size(f, 1)) STOP 'INIT_INTERPOLATOR_2D: Error, x should be the same size as first dimension of f'
      interp%nx = nx
      ALLOCATE(interp%x(nx))
      IF(present_and_correct(logx)) THEN
         interp%x = log(x)
         interp%logx = .TRUE.
      ELSE
         interp%x = x
         interp%logx = .FALSE.
      END IF
      IF (regular_spacing(interp%x)) THEN
         interp%ifindx = ifind_linear
      ELSE
         interp%ifindx = ifind_default
      END IF

      ! Sort out y axis
      ny = size(y)
      IF(ny /= size(f, 2)) STOP 'INIT_INTERPOLATOR_3D: Error, y should be the same size as second dimension of f'
      interp%ny = ny
      ALLOCATE(interp%y(ny))
      IF(present_and_correct(logy)) THEN
         interp%y = log(y)
         interp%logy = .TRUE.
      ELSE
         interp%y = y
         interp%logy = .FALSE.
      END IF
      IF (regular_spacing(interp%y)) THEN
         interp%ifindy = ifind_linear
      ELSE
         interp%ifindy = ifind_default
      END IF

      ! Sort out z axis
      nz = size(z)
      IF(nz /= size(f, 3)) STOP 'INIT_INTERPOLATOR_3D: Error, z should be the same size as third dimension of f'
      interp%nz = nz
      ALLOCATE(interp%z(nz))
      IF(present_and_correct(logz)) THEN
         interp%z = log(z)
         interp%logz = .TRUE.
      ELSE
         interp%z = z
         interp%logz = .FALSE.
      END IF
      IF (regular_spacing(interp%z)) THEN
         interp%ifindz = ifind_linear
      ELSE
         interp%ifindz = ifind_default
      END IF

      ! Sort out f
      ALLOCATE(interp%f(nx, ny, nz))
      IF(present_and_correct(logf)) THEN
         interp%f = log(f)
         interp%logf = .TRUE.
      ELSE
         interp%f = f
         interp%logf = .FALSE.
      END IF

      ! Reverse arrays if necessary
      IF (interp%x(1) > interp%x(nx)) THEN
         CALL reverse_array(interp%x)
         CALL reverse_array(interp%f, 1)
      END IF
      IF (interp%y(1) > interp%y(ny)) THEN
         CALL reverse_array(interp%y)
         CALL reverse_array(interp%f, 2)
      END IF
      IF (interp%z(1) > interp%z(nz)) THEN
         CALL reverse_array(interp%z)
         CALL reverse_array(interp%f, 3)
      END IF

      ! Set the interp%xmin, xmax, ymin, ymax values
      interp%xmin = interp%x(1)
      interp%xmax = interp%x(nx)
      interp%ymin = interp%y(1)
      interp%ymax = interp%y(ny)
      interp%zmin = interp%z(1)
      interp%zmax = interp%z(nz)
      IF (interp%logx) THEN
         interp%xmin = exp(interp%xmin)
         interp%xmax = exp(interp%xmax)
      END IF
      IF (interp%logy) THEN
         interp%ymin = exp(interp%ymin)
         interp%ymax = exp(interp%ymax)
      END IF
      IF (interp%logz) THEN
         interp%zmin = exp(interp%zmin)
         interp%zmax = exp(interp%zmax)
      END IF

      ! Set internal variables
      interp%iorder = iorder
      interp%iextrap = iextrap
      interp%store = store

   END SUBROUTINE init_interpolator_3D

   REAL FUNCTION evaluate_interpolator_3D(x, y, z, interp)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: z
      TYPE(interpolator3D), INTENT(IN) :: interp
      REAL :: xx, yy, zz
      INTEGER :: nx, ny, nz
      INTEGER :: iorder
      INTEGER, PARAMETER :: iinterp = iinterp_polynomial ! No Lagrange polynomials in 3D
      LOGICAL, PARAMETER :: xycubic = .FALSE.

      xx = x
      IF (interp%logx) xx = log(xx)
      nx = interp%nx

      yy = y
      IF (interp%logy) yy = log(yy)
      ny = interp%ny

      zz = z
      IF (interp%logz) zz = log(zz)
      nz = interp%nz

      IF (interp%store) THEN

         ! IF (outside_array_range(xx, interp%x) .OR. outside_array_range(yy, interp%y) .OR. outside_array_range(zz, interp%z)) THEN

         !    IF (interp%iextrap == iextrap_no) THEN
         !       STOP 'EVALUATE_INTERPOLATOR_3D: Error, you have requested a point outside the range of the interpolator'
         !    ELSE IF (interp%iextrap == iextrap_zero) THEN
         !       evaluate_interpolator_3D = 0.
         !    ELSE IF (interp%iextrap == iextrap_near) THEN
         !       STOP 'EVALUATE_INTERPOLATOR_3D: Error, nearest extrapolation not implemented yet'
         !    ELSE IF (interp%iextrap == iextrap_lin) THEN
         !       STOP 'EVALUATE_INTERPOLATOR_3D: Error, linear extrapolation not implemented yet'
         !    ELSE IF (interp%iextrap == iextrap_std) THEN
         !       STOP 'EVALUATE_INTERPOLATOR_3D: Error, standard extrapolation not implemented yet'
         !    ELSE
         !       STOP 'EVALUATE_INTERPOLATOR_3D: Something went wrong with extrapolation methods'
         !    END IF

         ! ELSE

         STOP 'EVALUATE_INTERPOLATOR_3D: Storage of interpolator not implemented yet'

         ! END IF

      ELSE

         IF (is_in_array(interp%iextrap, [iextrap_lin, iextrap_no, iextrap_zero, iextrap_near]) .AND. &
            (outside_array_range(xx, interp%x) .OR. outside_array_range(yy, interp%y) .OR. outside_array_range(zz, interp%z))) THEN

            IF (interp%iextrap == iextrap_no) THEN
               STOP 'EVALUATE_INTERPOLATOR_3D: Error, you have requested a point outside the range of the interpolator'
            ELSE IF (interp%iextrap == iextrap_zero) THEN
               evaluate_interpolator_3D = 0.
               RETURN
            ELSE IF (interp%iextrap == iextrap_lin) THEN
               iorder = 1
            ELSE IF (interp%iextrap == iextrap_near) THEN
               iorder = 0
            ELSE
               STOP 'EVALUATE_INTERPOLATOR_3D: Something went wrong with extrapolation methods'
            END IF

         ELSE

            iorder = interp%iorder

         END IF

            evaluate_interpolator_3D = find(xx, interp%x, yy, interp%y, zz, interp%z, interp%f, &
                        interp%nx, interp%ny, interp%nz, &
                        iorder, &
                        interp%ifindx, &
                        interp%ifindy, &
                        interp%ifindz, &
                        iinterp)

      END IF

      IF (interp%logf) evaluate_interpolator_3D = exp(evaluate_interpolator_3D)

   END FUNCTION evaluate_interpolator_3D

END MODULE interpolate
