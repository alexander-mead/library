PROGRAM interpolate_test

   USE basic_operations
   USE array_operations
   USE special_functions
   USE interpolate

   IMPLICIT NONE
   LOGICAL :: ifail

   WRITE(*, *)

   CALL test_interpolate_1D(ifail)

   CONTAINS

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
            y = linear_polynomial(x, a1, a0)
         ELSE IF (itest == 2) THEN
            y = quadratic_polynomial(x, a2, a1, a0)
         ELSE IF (itest == 3) THEN
            y = cubic_polynomial(x, a3, a2, a1, a0)
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
                        yt = linear_polynomial(xv, a1, a0)
                     ELSE IF (itest == 2) THEN
                        yt = quadratic_polynomial(xv, a2, a1, a0)
                     ELSE IF (itest == 3) THEN
                        yt = cubic_polynomial(xv, a3, a2, a1, a0)
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

END PROGRAM interpolate_test