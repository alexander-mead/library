PROGRAM interpolate_demo

   USE interpolate
   USE constants
   USE random_numbers
   USE basic_operations
   USE array_operations
   USE field_operations

   IMPLICIT NONE

   INTEGER :: test

   WRITE (*, *)
   WRITE (*, *) 'Routines for interpolating from a table of data'
   WRITE (*, *)

   CALL read_command_argument(1, test, '', def=-1)

   IF (test == -1) THEN
      WRITE (*, *) 'Choose test'
      WRITE (*, *) '1 - Interpolate 1D'
      WRITE (*, *) '2 - Interpolate 2D'
      WRITE (*, *) '3 - Interpolate 3D: Linear'
      WRITE (*, *) '4 - Interpolate 3D: Sine'
      WRITE (*, *) '5 - Interpolator 1D'
      READ (*, *) test
      WRITE (*, *)
   END IF

   IF (test == 1) THEN
      CALL interpolate_demo_1D()
   ELSE IF (test == 2) THEN
      CALL interpolate_demo_2D()
   ELSE IF (test == 3) THEN
      CALL interpolate_demo_3D(test)
   ELSE IF (test == 4) THEN
      CALL interpolate_demo_3D(test)
   ELSE IF (test == 5) THEN
      CALL interpolator_demo_1D
   ELSE
      STOP 'INTERPOLATE_DEMO: Error, demo specified incorrectly'
   END IF

   CONTAINS

   SUBROUTINE interpolator_demo_1D()

      TYPE(interpolator1D) :: interp
      REAL, ALLOCATABLE :: x(:), f(:)
      REAL :: xv, fv
      INTEGER :: ix
      REAL, PARAMETER :: xmin = 0.
      REAL, PARAMETER :: xmax = pi
      INTEGER, PARAMETER :: nx = 9
      INTEGER, PARAMETER :: nnx = 10*nx
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: iextrap = iextrap_standard

      CALL fill_array(xmin, xmax, x, nx)
      f = sin(x)

      CALL init_interpolator(x, f, interp, iorder, iextrap, store=.TRUE.)

      DO ix = 1, nnx
         xv = progression(xmin, xmax, ix, nnx)
         fv = evaluate_interpolator(xv, interp)
         WRITE(*, *) ix, xv, fv, sin(xv)!, fv/sin(xv)
      END DO

   END SUBROUTINE interpolator_demo_1D

   SUBROUTINE interpolate_demo_1D()

      REAL, ALLOCATABLE :: xtab(:), ytab(:)
      REAL :: lin, quad, cube
      REAL :: rlin, rquad, rcube
      REAL :: tru
      REAL :: x, xmin, xmax
      INTEGER :: i, m, n
      INTEGER :: iexample, imeth

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

      !DO i = 1, n
      !   xtab(i) = xmin+(xmax-xmin)*float(i-1)/float(n-1)
      !END DO
      CALL fill_array(xmin, xmax, xtab, n)

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
         !x = xmin+(xmax-xmin)*float(i-1)/float(m-1)
         x = progression(xmin, xmax, i, m)
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

         WRITE (*, *) 'Constant interpolation:', find(x, xtab, ytab, n, 0, 3, imeth)
         WRITE (*, *) 'Linear interpolation:', find(x, xtab, ytab, n, 1, 3, imeth)
         WRITE (*, *) 'Quadratic interpolation:', find(x, xtab, ytab, n, 2, 3, imeth)
         WRITE (*, *) 'Cubic interpolation:', find(x, xtab, ytab, n, 3, 3, imeth)
         WRITE (*, *) 'Truth:', truth(iexample, x)
         WRITE (*, *)

      END DO

      WRITE (*, *)

   END SUBROUTINE interpolate_demo_1D
   
   SUBROUTINE interpolate_demo_2D()

      REAL, ALLOCATABLE :: f(:, :)
      REAL :: f_int0, f_int1, f_int3, f_true
      REAL :: x, xmin, xmax
      REAL :: y, ymin, ymax
      REAL, ALLOCATABLE :: xtab(:), ytab(:)
      INTEGER :: nxtab, nytab, nx, ny
      INTEGER :: i, j

      ! x range for tabulated function to interpolate
      nxtab = 11
      xmin = 0.
      xmax = 1.
      CALL fill_array(xmin, xmax, xtab, nxtab)

      ! y range for tabulated function to interpolate
      nytab = 11
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

      ! Random point to test
      !CALL RNG_set(seed=0)
      CALL random_generator_seed(seed=0)
      x = random_uniform(xmin, xmax)
      y = random_uniform(ymin, ymax)
      WRITE (*, *) 'x:', x
      WRITE (*, *) 'y:', y
      WRITE (*, *) 'Constant interpolation:', find(x, xtab, y, ytab, f, nxtab, nytab, 0, 3, 3, 1)
      WRITE (*, *) 'Linear interpolation:', find(x, xtab, y, ytab, f, nxtab, nytab, 1, 3, 3, 1)
      WRITE (*, *) 'Cubic interpolation:', find(x, xtab, y, ytab, f, nxtab, nytab, 3, 3, 3, 1)
      WRITE (*, *) 'Truth:', func(x, y)
      WRITE (*, *)

      ! x range for test
      xmin = 0.
      xmax = 1.0
      nx = 100

      ! y range for test
      ymin = 0.
      ymax = 1.
      ny = 100

      ! Test
      OPEN (7, file='results_2D.dat')
      DO i = 1, nx
         DO j = 1, ny

            x = cell_position(i, xmax-xmin, nx)+xmin
            y = cell_position(j, ymax-ymin, ny)+ymin

            f_true = func(x, y)
            f_int0 = find(x, xtab, y, ytab, f, nxtab, nytab, 0, 3, 3, 1)
            f_int1 = find(x, xtab, y, ytab, f, nxtab, nytab, 1, 3, 3, 1)
            f_int3 = find(x, xtab, y, ytab, f, nxtab, nytab, 3, 3, 3, 1)

            WRITE (7, *) x, y, f_true, f_int0, f_int1, f_int3

         END DO
      END DO
      CLOSE (7)

   END SUBROUTINE interpolate_demo_2D

   SUBROUTINE interpolate_demo_3D(itest)

      INTEGER, INTENT(IN) :: itest
      REAL, ALLOCATABLE :: g(:, :, :)
      REAL :: f_int1, f_true
      REAL :: x, xmin, xmax
      REAL :: y, ymin, ymax
      REAL :: z, zmin, zmax
      REAL, ALLOCATABLE :: xtab(:), ytab(:), ztab(:)
      INTEGER :: ix, iy, iz
      INTEGER :: nx, ny, nz


      IF (itest == 3) THEN

         xmin = 0.
         xmax = 1.
         nx = 11

         ymin = -0.5
         ymax = 0.5
         ny = 11

         zmin = 0.
         zmax = 2.
         nz = 21

      ELSE IF (itest == 4) THEN

         xmin = 0.
         xmax = 0.5*pi
         nx = 11

         ymin = 0.
         ymax = pi
         ny = 21

         zmin = 0.
         zmax = 1.5*pi
         nz = 31

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
      f_int1 = find(x, xtab, y, ytab, z, ztab, g, nx, ny, nz, iorder=1, ifind=3, iinterp=1)
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

   END SUBROUTINE interpolate_demo_3D

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

      !func = sin(x*y)+1.
      func = 1.+sin(2.*pi*x)+1.+sin(2.*pi*y)

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
