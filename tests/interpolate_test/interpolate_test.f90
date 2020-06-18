PROGRAM interpolate_test

   USE constants
   USE basic_operations
   USE array_operations
   USE special_functions
   USE interpolate

   IMPLICIT NONE
   LOGICAL :: ifail

   WRITE(*, *)

   CALL test_interpolate_1D(ifail)
   CALL test_interpolate_1D_function(ifail)
   CALL test_interpolate_2D(ifail)
   CALL test_interpolator_1D(ifail)
   CALL test_interpolator_2D(ifail)

   CONTAINS

   SUBROUTINE test_interpolator_1D(fail)

      LOGICAL, INTENT(OUT) :: fail
      REAL, ALLOCATABLE :: x(:)
      REAL, ALLOCATABLE :: f(:)
      REAL :: xv, fv, ft
      INTEGER :: i, itest
      TYPE(interpolator1D) :: interp
      INTEGER :: iextrap, iorder
      LOGICAL :: logtest, store
      REAL :: xmin, xmax
      REAL :: xxmin, xxmax
      INTEGER, PARAMETER :: n = 8
      INTEGER, PARAMETER :: m = 128
      REAL, PARAMETER :: tol = 1e-8
      REAL, PARAMETER :: a = 2.
      INTEGER, PARAMETER :: ntest = 9

      fail = .FALSE.

      DO itest = 1, ntest

         IF (itest == 1 .OR. itest == 8 .OR. itest == 9) THEN
            xmin = 3.
            xmax = 0.
            IF (itest == 9) THEN
               xmin = 3.
               xmax = 0.
            END IF
            xxmin = xmin
            xxmax = xmax
            logtest = .FALSE.
            iorder = 3
            iextrap = iextrap_no
            IF(itest == 1) store = .TRUE.
            IF(itest == 8) store = .FALSE.
         ELSE IF (itest == 2) THEN
            xmin = 1e-3
            xmax = 3.
            xxmin = xmin
            xxmax = xmax
            logtest = .TRUE.
            iorder = 3
            iextrap = iextrap_no
            store = .TRUE.
         ELSE IF (itest == 3) THEN
            xmin = 0.
            xmax = 3.
            xxmin = xmin
            xxmax = xmax
            logtest = .FALSE.
            iorder = 3
            iextrap = iextrap_linear
            store = .TRUE.
         ELSE IF (itest == 4) THEN
            xmin = 0.
            xmax = 3.
            xxmin = xmin
            xxmax = xmax
            logtest = .FALSE.
            iorder = 3
            iextrap = iextrap_standard
            store = .TRUE.
         ELSE IF (itest == 5) THEN
            xmin = 0.
            xmax = 3.
            xxmin = xmin
            xxmax = xmax
            logtest = .FALSE.
            iorder = 2
            iextrap = iextrap_no
            store = .TRUE.
         ELSE IF (itest == 6) THEN
            xmin = 0.
            xmax = 3.
            xxmin = xmin
            xxmax = xmax
            logtest = .FALSE.
            iorder = 1
            iextrap = iextrap_no
            store = .TRUE.
         ELSE IF (itest == 7) THEN
            xmin = 0.
            xmax = 3.
            xxmin = -3.
            xxmax = 6.
            logtest = .FALSE.
            iorder = 1
            iextrap = iextrap_linear
            store = .TRUE.
         ELSE
            STOP 'TEST_INTERPOLATOR_1D: Something went wrong'
         END IF

         CALL fill_array(xmin, xmax, x, n)
         CALL safe_allocate(f, n)
         f = a*x**iorder

         CALL init_interpolator(x, f, interp, iorder, iextrap, store=store, logx=logtest, logf=logtest)

         DO i = 1, m
            xv = progression(xxmin, xxmax, i, m)
            fv = evaluate_interpolator(xv, interp)
            ft = a*xv**iorder
            IF(.NOT. requal(fv, ft, tol)) THEN
               fail = .TRUE.
               WRITE(*, *) 'TEST_INTERPOLATOR_1D: Fail'
               WRITE(*, *) 'TEST_INTERPOLATOR_1D: itest:', itest
               WRITE(*, *) 'TEST_INTERPOLATOR_1D: i:', i
               WRITE(*, *) 'TEST_INTERPOLATOR_1D: x:', xv
               WRITE(*, *) 'TEST_INTERPOLATOR_1D: f true:', ft
               WRITE(*, *) 'TEST_INTERPOLATOR_1D: f interp:', fv
               STOP
            END IF
         END DO

      END DO

      IF (fail) THEN
         WRITE(*, *) 'TEST_INTERPOLATOR_1D: Fail'
      ELSE
         WRITE(*, *) 'TEST_INTERPOLATOR_1D: Pass'
         WRITE(*, *)
      END IF

   END SUBROUTINE test_interpolator_1D

   SUBROUTINE test_interpolator_2D(fail)

      LOGICAL, INTENT(OUT) :: fail
      REAL, ALLOCATABLE :: x(:), y(:)
      REAL, ALLOCATABLE :: f(:, :)
      REAL :: xv, yv, fv, ft
      INTEGER :: ix, iy, itest
      TYPE(interpolator2D) :: interp
      INTEGER :: iextrap, iorder
      LOGICAL :: logtest, store
      REAL :: xmin, xmax, ymin, ymax
      REAL :: xmin_test, xmax_test, ymin_test, ymax_test
      INTEGER, PARAMETER :: nx = 8
      INTEGER, PARAMETER :: ny = 10
      INTEGER, PARAMETER :: mx = 64
      INTEGER, PARAMETER :: my = 32
      REAL, PARAMETER :: tol = 1e-8
      REAL, PARAMETER :: a = 2.
      INTEGER, PARAMETER :: ntest = 10

      fail = .FALSE.

      DO itest = 1, ntest

         IF (is_in_array(itest, [1, 3, 4, 5, 6, 7, 8, 9, 10])) THEN
            xmin = 0.
            xmax = 3.
            ymin = 0.
            ymax = 5.
            logtest = .FALSE.
            IF (itest == 3) THEN
               iorder = 2
            ELSE IF (is_in_array(itest, [1, 4, 5, 6, 9])) THEN
               iorder = 3
            ELSE
               iorder = 1
            END IF
            IF (itest == 9 .OR. itest == 10) THEN
               iextrap = iextrap_linear
            ELSE
               iextrap = iextrap_standard
            END IF
            store = .TRUE.
         ELSE IF (itest == 2) THEN
            xmin = 1e-3
            xmax = 3.
            ymin = 1e-3
            ymax = 3.
            logtest = .TRUE.
            iorder = 3
            iextrap = iextrap_standard
            store = .TRUE.
         ELSE
            STOP 'TEST_INTERPOLATOR_2D: Something went wrong'
         END IF

         CALL fill_array(xmin, xmax, x, nx)
         CALL fill_array(ymin, ymax, y, ny)
         IF (itest == 4 .OR. itest == 5) THEN
            CALL reverse_array(x)
         ELSE IF (itest == 4 .OR. itest ==6) THEN
            CALL reverse_array(y)
         END IF
         CALL safe_allocate(f, nx, ny)        
         DO iy = 1, ny
            DO ix = 1, nx
               f(ix, iy) = test_function_2D(x(ix), y(iy), itest)
            END DO
         END DO
         CALL init_interpolator(x, y, f, interp, iorder, iextrap, store=store, logx=logtest, logy=logtest, logf=logtest)

         IF (itest == 8 .OR. itest == 9 .OR. itest == 10) THEN
            xmin_test = xmin-1.
            xmax_test = xmax+1.
            ymin_test = ymin-1.
            ymax_test = ymax+1.
         ELSE
            xmin_test = xmin
            xmax_test = xmax
            ymin_test = ymin
            ymax_test = ymax
         END IF
         DO ix = 1, mx
            DO iy = 1, my
               xv = progression(xmin_test, xmax_test, ix, mx)
               yv = progression(ymin_test, ymax_test, iy, my)
               fv = evaluate_interpolator(xv, yv, interp)
               ft = test_function_2D(xv, yv, itest)
               IF(.NOT. requal(fv, ft, tol)) THEN
                  fail = .TRUE.
                  WRITE(*, *) 'TEST_INTERPOLATOR_2D: Fail'
                  WRITE(*, *) 'TEST_INTERPOLATOR_2D: itest:', itest
                  WRITE(*, *) 'TEST_INTERPOLATOR_2D: ix, iy:', ix, iy
                  WRITE(*, *) 'TEST_INTERPOLATOR_2D: x, y:', xv, yv
                  WRITE(*, *) 'TEST_INTERPOLATOR_2D: f true:', ft
                  WRITE(*, *) 'TEST_INTERPOLATOR_2D: f interp:', fv
                  STOP
               END IF
            END DO
         END DO

      END DO

      IF (fail) THEN
         WRITE(*, *) 'TEST_INTERPOLATOR_2D: Fail'
      ELSE
         WRITE(*, *) 'TEST_INTERPOLATOR_2D: Pass'
         WRITE(*, *)
      END IF

   END SUBROUTINE test_interpolator_2D  

   REAL FUNCTION test_function_2D(x, y, itest) RESULT(f)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      INTEGER, INTENT(IN) :: itest

      IF (itest == 1 .OR. itest == 4 .OR. itest == 5 .OR. itest == 6) THEN
         f = x**2+2.*y**2*3.*x*y
      ELSE IF (itest == 2) THEN
         f = 12.*x*y
      ELSE IF (itest == 3) THEN
         f = x+3.*y+1.
      ELSE IF (itest == 7 .OR. itest == 8 .OR. itest == 9 .OR. itest == 10) THEN
         f = 3.*x+4.*y
      ELSE
         STOP 'TEST_FUNCTION_2D: Error, itest specified incorrectly'
      END IF

   END FUNCTION test_function_2D

   SUBROUTINE test_interpolate_1D(fail)

      LOGICAL, INTENT(OUT) :: fail
      REAL :: xmin, xmax
      REAL :: a3, a2, a1, a0
      REAL :: xv, yv, yt
      REAL, ALLOCATABLE :: x(:), y(:)
      INTEGER :: itest, i, iorder1, iorder2
      INTEGER :: iorder, ifind, iinterp
      INTEGER, PARAMETER :: n = 8
      INTEGER, PARAMETER :: m = 128
      REAL, PARAMETER :: tol = 1e-6

      fail = .FALSE.

      DO itest = 1, 3

         IF (itest == 1) THEN
            xmin = 0.
            xmax = 10.
            a1 = -4.
            a0 = 9.
            iorder1 = 1
            iorder2 = 3
         ELSE IF (itest == 2) THEN
            xmin = -3.
            xmax = 3.
            a2 = 2.
            a1 = -4.
            a0 = 0.1
            iorder1 = 2
            iorder2 = 3
         ELSE IF (itest == 3) THEN
            xmin = -10.
            xmax = 10.
            a3 = 1.
            a2 = -6.
            a1 = -4.
            a0 = 0.5
            iorder1 = 3
            iorder2 = 3
         ELSE
            STOP 'TEST_INTERPOLATE_1D: Error, something went wrong'
         END IF

         CALL fill_array(xmin, xmax, x, n)
         CALL safe_allocate(y, n)

         IF (itest == 1) THEN
            y = polynomial(x, a1, a0)
         ELSE IF (itest == 2) THEN
            y = polynomial(x, a2, a1, a0)
         ELSE IF (itest == 3) THEN
            y = polynomial(x, a3, a2, a1, a0)
         ELSE
            STOP 'TEST_INTERPOLATE_1D: Error, something went wrong'
         END IF

         DO iorder = iorder1, iorder2
            DO ifind = 1, 3
               DO iinterp = 1, 2

                  DO i = 1, m

                     xv = progression(xmin, xmax, i, m)
                     yv = find(xv, x, y, n, iorder, ifind, iinterp)

                     IF(itest == 1) THEN
                        yt = polynomial(xv, a1, a0)
                     ELSE IF (itest == 2) THEN
                        yt = polynomial(xv, a2, a1, a0)
                     ELSE IF (itest == 3) THEN
                        yt = polynomial(xv, a3, a2, a1, a0)
                     ELSE
                        STOP 'TEST_INTERPOLATE_1D: Error, something went wrong'
                     END IF

                     IF(.NOT. requal(yv, yt, tol)) THEN
                        fail = .TRUE.
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: Test failed:', itest
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: Order:', iorder
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: Find:', ifind
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: iinterp:', iinterp
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: i:', i
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: x:', xv
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: Interpolated y:', yv
                        WRITE(*, *) 'TEST_INTERPOLATE_1D: True y:', yt
                        STOP
                     END IF

                  END DO

               END DO
            END DO
         END DO

      END DO

      IF (.NOT. fail) THEN
         WRITE(*, *) 'TEST_INTERPOLATE_1D: Pass'
         WRITE(*, *)
      END IF

   END SUBROUTINE test_interpolate_1D

   SUBROUTINE test_interpolate_1D_function(fail)

      LOGICAL, INTENT(OUT) :: fail
      REAL :: xmin, xmax
      REAL :: a3, a2, a1, a0
      REAL :: xv, yv, yt
      REAL, ALLOCATABLE :: x(:), y(:)
      INTEGER :: itest, i, j
      INTEGER :: iorder, ifind, iinterp
      INTEGER, PARAMETER :: n = 32
      INTEGER, PARAMETER :: m = 128
      REAL, PARAMETER :: tol = 1e-3
      INTEGER, PARAMETER :: ntest = 2

      fail = .FALSE.

      DO itest = 1, ntest

         IF (itest == 1) THEN
            xmin = -1.
            xmax = 2.
            iorder = 3
         ELSE IF (itest == 2) THEN
            xmin = 0.
            xmax = pi
            iorder = 3
         ELSE
            STOP 'TEST_INTERPOLATE_1D_FUNCTION: Error, something went wrong'
         END IF

         CALL fill_array(xmin, xmax, x, n)
         CALL safe_allocate(y, n)

         IF (itest == 1) THEN
            y = interpolate_test_cubic(x)
         ELSE IF (itest == 2) THEN
            y = interpolate_test_sin(x)
         ELSE
            STOP 'TEST_INTERPOLATE_1D_FUNCTION: Error, something went wrong'
         END IF

         DO ifind = 1, 3
            DO iinterp = 1, 2

               DO i = 1, m

                  xv = progression(xmin-0.01, xmax+0.01, i, m)
                  yv = find(xv, x, y, n, iorder, ifind, iinterp)

                  IF (itest == 1) THEN
                     yt = interpolate_test_cubic(xv)
                  ELSE IF (itest == 2) THEN
                     yt = interpolate_test_sin(xv)
                  ELSE
                     STOP 'TEST_INTERPOLATE_1D_FUNCTION: Error, something went wrong'
                  END IF

                  !WRITE(*, *) i, xv, yv/yt

                  IF(.NOT. requal(yv, yt, tol)) THEN
                     fail = .TRUE.
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: Test failed:', itest
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: Order:', iorder
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: Find:', ifind
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: iinterp:', iinterp
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: i:', i
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: x:', xv
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: Interpolated y:', yv
                     WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: True y:', yt
                     STOP
                  END IF

               END DO

            END DO
         END DO

      END DO

      IF (.NOT. fail) THEN
         WRITE(*, *) 'TEST_INTERPOLATE_1D_FUNCTION: Pass'
         WRITE(*, *)
      END IF

   END SUBROUTINE test_interpolate_1D_function

   ELEMENTAL REAL FUNCTION interpolate_test_cubic(x)

   REAL, INTENT(IN) :: x

      IF (x < 0.) THEN
         interpolate_test_cubic = 1.
      ELSE
         interpolate_test_cubic = x**3+1.
      END IF

   END FUNCTION interpolate_test_cubic

   ELEMENTAL REAL FUNCTION interpolate_test_sin(x)

      REAL, INTENT(IN) :: x

      IF (x < 0.) THEN
         interpolate_test_sin = x
      ELSE IF (x > pi) THEN
         interpolate_test_sin = pi-x
      ELSE
         interpolate_test_sin = sin(x)
      END IF

   END FUNCTION interpolate_test_sin

   SUBROUTINE test_interpolate_2D(fail)

      LOGICAL, INTENT(OUT) :: fail
      REAL :: xmin, xmax, ymin, ymax
      REAL :: a3, a2, a1, a0
      REAL :: xv, yv, fv, ft
      REAL, ALLOCATABLE :: x(:), y(:), f(:, :)
      INTEGER :: itest, ix, iy
      INTEGER :: iorder, ifind
      INTEGER, PARAMETER :: ni = 8
      INTEGER, PARAMETER :: nt = 64
      REAL, PARAMETER :: tol = 1e-8
      INTEGER, PARAMETER :: ntest = 1
      INTEGER, PARAMETER :: iinterp = 1

      fail = .FALSE.

      DO itest = 1, ntest

         IF (itest == 1) THEN
            xmin = 0.
            xmax = 4.
            ymin = xmin
            ymax = xmax
            iorder = 3
         ELSE
            STOP 'TEST_INTERPOLATE_2D: Error, something went wrong'
         END IF

         CALL fill_array(xmin, xmax, x, ni)
         CALL fill_array(ymin, ymax, y, ni)
         ALLOCATE(f(ni, ni))

         IF (itest == 1) THEN
            DO ix = 1, ni
               DO iy = 1, ni
                  f(ix, iy) = interpolate_test_surface(x(ix), y(iy))
               END DO
            END DO
         ELSE
            STOP 'TEST_INTERPOLATE_2D: Error, something went wrong'
         END IF


         DO ifind = 1, 3

               DO ix = 1, nt
                  DO iy = 1, nt

                     xv = progression(xmin, xmax, ix, nt)
                     yv = progression(ymin, ymax, iy, nt)
                     fv = find(xv, x, yv, y, f, ni, ni, iorder, ifind, ifind, iinterp)

                     IF (itest == 1) THEN
                        ft = interpolate_test_surface(xv, yv)
                     ELSE
                        STOP 'TEST_INTERPOLATE_2D: Error, something went wrong'
                     END IF

                     IF(.NOT. requal(fv, ft, tol)) THEN
                        fail = .TRUE.
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: Test failed:', itest
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: Order:', iorder
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: Find:', ifind
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: iinterp:', iinterp
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: ix, iy:', ix, iy
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: x, y:', xv, yv
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: Interpolated f:', fv
                        WRITE(*, *) 'TEST_INTERPOLATE_2D: True f:', ft
                        STOP
                     END IF

                  END DO
               END DO

         END DO

      END DO

      IF (.NOT. fail) THEN
         WRITE(*, *) 'TEST_INTERPOLATE_2D: Pass'
         WRITE(*, *)
      END IF

   END SUBROUTINE test_interpolate_2D

   ELEMENTAL REAL FUNCTION interpolate_test_surface(x, y)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y

      interpolate_test_surface = (x-2.)**2+(y-2.)**2+8.
      !interpolate_test_surface = x

   END FUNCTION interpolate_test_surface

END PROGRAM interpolate_test