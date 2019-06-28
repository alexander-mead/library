PROGRAM interpolate_demo

   USE interpolate
   USE constants
   USE array_operations
   USE field_operations

   IMPLICIT NONE
   REAL :: xmin, xmax, x, ymin, ymax, y, zmin, zmax, z
   REAL, ALLOCATABLE :: xtab(:), ytab(:), ztab(:), f(:, :), g(:, :, :)
   REAL :: f_int1, f_int3, f_true
   REAL :: lin, quad, cube, tru
   REAL :: rlin, rquad, rcube
   INTEGER :: i, j, n, m, nxtab, nytab
   INTEGER :: ix, iy, iz, nx, ny, nz
   INTEGER :: iexample, imeth, itest

   WRITE (*, *)
   WRITE (*, *) 'Routines for interpolating from a table of data'
   WRITE (*, *)

   WRITE (*, *) 'Choose test'
   WRITE (*, *) '1 - 1D'
   WRITE (*, *) '2 - 2D'
   WRITE (*, *) '3 - 3D: Linear'
   WRITE (*, *) '4 - 3D: Sine'
   READ (*, *) itest
   WRITE (*, *)

   IF (itest == 1) THEN

      WRITE (*, *) 'Number of points for tables (low is a good test):'
      READ (*, *) n
      WRITE (*, *)

      WRITE (*, *) 'Method to test'
      WRITE (*, *) '1 - Standard polynomial fitting'
      WRITE (*, *) '2 - Lagrange polynomial fitting'
      READ (*, *) imeth
      WRITE (*, *)

      ALLOCATE (xtab(n), ytab(n))

      WRITE (*, *) '0 - Linear (0 to 3)'
      WRITE (*, *) '1 - Quadratic (0 to 3)'
      WRITE (*, *) '2 - Cubic (0 to 3)'
      WRITE (*, *) '3 - Sin (0 to pi)'
      WRITE (*, *) '4 - Exp (0 to 3)'
      READ (*, *) iexample
      WRITE (*, *)

      IF (iexample == 0) THEN
         xmin = 0.
         xmax = 3.
      ELSE IF (iexample == 1) THEN
         xmin = 0.
         xmax = 3.
      ELSE IF (iexample == 2) THEN
         xmin = 0.
         xmax = 3.
      ELSE IF (iexample == 3) THEN
         xmin = 0.
         xmax = pi
      ELSE IF (iexample == 4) THEN
         xmin = 0.
         xmax = 3.
      ELSE
         STOP 'Error, example not specified correctly'
      END IF

      DO i = 1, n
         xtab(i) = xmin+(xmax-xmin)*float(i-1)/float(n-1)
      END DO

      IF (iexample == 0) ytab = xtab
      IF (iexample == 1) ytab = xtab**2
      IF (iexample == 2) ytab = xtab**3
      IF (iexample == 3) ytab = sin(xtab)
      IF (iexample == 4) ytab = exp(xtab)

      OPEN (7, file='table.dat')
      DO i = 1, n
         WRITE (7, *) xtab(i), ytab(i)
      END DO
      CLOSE (7)

      m = 10*n

      WRITE (*, *) 'Writing tables'
      OPEN (7, file='results.dat')
      OPEN (8, file='ratio.dat')
      DO i = 1, m
         x = xmin+(xmax-xmin)*float(i-1)/float(m-1)
         lin = find(x, xtab, ytab, n, 1, 2, imeth)
         quad = find(x, xtab, ytab, n, 2, 2, imeth)
         cube = find(x, xtab, ytab, n, 3, 2, imeth)
         tru = truth(iexample, x)
         WRITE (7, *) x, lin, quad, cube, tru
         IF (tru == 0.) THEN
            rlin = 1.
            rquad = 1.
            rcube = 1.
         ELSE
            rlin = lin/tru
            rquad = quad/tru
            rcube = cube/tru
         END IF
         WRITE (8, *) x, rlin, rquad, rcube
      END DO
      CLOSE (7)
      CLOSE (8)
      WRITE (*, *) 'Done'
      WRITE (*, *)

      DO

         WRITE (*, *) 'Value to *find* for tests (-1 exits):'
         READ (*, *) x

         IF (x == -1.) EXIT

         WRITE (*, *) 'Linear interpolation:', find(x, xtab, ytab, n, 1, 3, imeth)
         WRITE (*, *) 'Quadratic interpolation:', find(x, xtab, ytab, n, 2, 3, imeth)
         WRITE (*, *) 'Cubic interpolation:', find(x, xtab, ytab, n, 3, 3, imeth)
         WRITE (*, *) 'Truth:', truth(iexample, x)
         WRITE (*, *)

      END DO

      WRITE (*, *)

   ELSE IF (itest == 2) THEN

      ! x range for tabulated function to interpolate
      nxtab = 4
      xmin = 0.
      xmax = 1.
      CALL fill_array(xmin, xmax, xtab, nxtab)

      ! y range for tabulated function to interpolate
      nytab = 4
      ymin = 0.
      ymax = 1.
      CALL fill_array(ymin, ymax, ytab, nytab)

      ! Fill array for tabulated function
      ALLOCATE (f(nxtab, nytab))
      DO i = 1, nxtab
         DO j = 1, nytab
            f(i, j) = func(xtab(i), ytab(j))
            WRITE (*, *) i, j, f(i, j)
         END DO
      END DO

      ! One-off point to test
      x = 0.9
      y = 0.3
      WRITE (*, *) 'x:', x
      WRITE (*, *) 'y:', y
      WRITE (*, *) 'Linear Interpolation:', find(x, xtab, y, ytab, f, nxtab, nytab, 1, 3, 1)
      WRITE (*, *) 'Cubic Interpolation:', find(x, xtab, y, ytab, f, nxtab, nytab, 3, 3, 1)
      WRITE (*, *) 'Truth:', func(x, y)
      WRITE (*, *)

      ! x range for test
      xmin = 0.
      xmax = 1.2
      nx = 100

      ! y range for test
      ymin = 0.
      ymax = 1.
      ny = 100

      ! Test
      OPEN (7, file='results_linear.dat')
      OPEN (8, file='results_cubic.dat')
      DO i = 1, nx
         DO j = 1, ny

            x = cell_position(i, xmax-xmin, nx)+xmin
            y = cell_position(j, ymax-ymin, ny)+ymin

            f_true = func(x, y)
            f_int1 = find(x, xtab, y, ytab, f, nxtab, nytab, 1, 3, 1)
            f_int3 = find(x, xtab, y, ytab, f, nxtab, nytab, 3, 3, 1)

            WRITE (7, *) x, y, f_true, f_int1, f_int1/f_true
            WRITE (8, *) x, y, f_true, f_int3, f_int3/f_true

         END DO
      END DO
      CLOSE (7)
      CLOSE (8)

   ELSE IF (itest == 3 .OR. itest == 4) THEN

      IF (itest == 3) THEN

         xmin = 0.
         xmax = 1.
         nx = 11

         ymin = 0.
         ymax = 1.
         ny = 11

         zmin = 0.
         zmax = 1.
         nz = 11

      ELSE IF (itest == 4) THEN

         xmin = 0.
         xmax = pi/2.
         nx = 11

         ymin = 0.
         ymax = pi/2.
         ny = 11

         zmin = 0.
         zmax = pi/2.
         nz = 11

      ELSE
         STOP 'INTERPOLATE_DEMO: Error, itest not specified correctly'
      END IF

      CALL fill_array(xmin, xmax, xtab, nx)
      CALL fill_array(ymin, ymax, ytab, ny)
      CALL fill_array(zmin, zmax, ztab, nz)

      ALLOCATE (g(nx, ny, nz))
      DO ix = 1, nx
         DO iy = 1, ny
            DO iz = 1, nz
               x = xtab(ix)
               y = ytab(iy)
               z = ztab(iz)
               IF (itest == 3) THEN
                  g(ix, iy, iz) = linear_func_3D(x, y, z)
               ELSE IF (itest == 4) THEN
                  g(ix, iy, iz) = sine_func_3D(x, y, z)
               ELSE
                  STOP 'INTERPOLATE_DEMO: Error, itest not specified correctly'
               END IF
            END DO
         END DO
      END DO

      x = 1.5
      y = 1.5
      z = 1.5
      WRITE(*,*) 'x:', x
      WRITE(*,*) 'y:', y
      WRITE(*,*) 'z:', z
      f_int1 = find(x, xtab, y, ytab, z, ztab, g, nx, ny, nz, iorder=1, ifind=3, imeth=1)
      IF (itest == 3) THEN
         f_true = linear_func_3D(x, y, z)
      ELSE IF (itest == 4) THEN
         f_true = sine_func_3D(x, y, z)
      ELSE
         STOP 'INTERPOLATE_DEMO: Error, itest not specified correctly'
      END IF

      WRITE (*, *) 'Truth:', f_true
      WRITE (*, *) 'Interpolation:', f_int1
      WRITE (*, *) 'Ratio:', f_true/f_int1
      WRITE (*, *)

   ELSE

      STOP 'INTERPOLATE_DEMO: Error, itest not specified correctly'

   END IF

CONTAINS

   FUNCTION truth(iexample, x)

      IMPLICIT NONE
      REAL :: x, truth
      INTEGER :: iexample

      IF (iexample == 0) THEN
         truth = x
      ELSE IF (iexample == 1) THEN
         truth = x**2
      ELSE IF (iexample == 2) THEN
         truth = x**3
      ELSE IF (iexample == 3) THEN
         truth = sin(x)
      ELSE IF (iexample == 4) THEN
         truth = exp(x)
      ELSE
         STOP 'Example not specified correctly'
      END IF

   END FUNCTION truth

   REAL FUNCTION func(x, y)

      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y

      func = sin(x*y)+1.

   END FUNCTION func

   REAL FUNCTION linear_func_3D(x, y, z)

      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: z

      linear_func_3D = x+y+z

   END FUNCTION linear_func_3D

   REAL FUNCTION sine_func_3D(x, y, z)

      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: z

      sine_func_3D = sin(x)*sin(y)*sin(z)

   END FUNCTION sine_func_3D

END PROGRAM interpolate_demo
