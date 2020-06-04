PROGRAM special_functions_test

   USE basic_operations
   USE special_functions

   IMPLICIT NONE
   LOGICAL :: ifail

   WRITE(*, *)

   CALL test_fix_linear(ifail)
   CALL test_fix_quadratic(ifail)
   CALL test_fix_cubic(ifail)
   CALL test_factorial(ifail)

   contains

   SUBROUTINE test_fix_linear(fail)

      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: fail
      REAL :: x1, y1, x2, y2
      REAL :: a0, a1, b0, b1
      INTEGER :: itest
      INTEGER, PARAMETER :: ntest = 3
      REAL, PARAMETER :: tol = 1e-6

      fail = .FALSE.

      DO itest = 1, ntest

         IF (itest == 1) THEN
            x1 = 0.
            x2 = 1.
            b0 = 3.
            b1 = 2.     
         ELSE IF (itest == 2) THEN
            x1 = 4.
            x2 = 10.
            b0 = -2.
            b1 = 0.  
         ELSE IF (itest == 3) THEN
            x1 = -3.
            x2 = 19.
            b0 = 0.5
            b1 = 2.5     
         ELSE
            STOP 'TEST_FIX_LINEAR: Something went wrong'
         END IF

         y1 = b1*x1+b0
         y2 = b1*x2+b0

         CALL fix_polynomial(a1, a0, [x1, x2], [y1, y2])

         IF ((.NOT. requal(a1, b1, tol)) .OR. (.NOT. requal(a0, b0, tol))) THEN
            WRITE(*, *) 'TEST_FIX_LINEAR: Target a1:', b1
            WRITE(*, *) 'TEST_FIX_LINEAR: Target a0:', b0
            WRITE(*, *) 'TEST_FIX_LINEAR: Fix a1:', a1
            WRITE(*, *) 'TEST_FIX_LINEAR: Fix a0:', a0
            STOP 'TEST_FIX_LINEAR: Fail'
         END IF

      END DO

      WRITE(*, *) 'TEST_FIX_LINEAR: Pass'
      WRITE(*, *)

   END SUBROUTINE test_fix_linear

   SUBROUTINE test_fix_quadratic(fail)

      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: fail
      REAL :: x1, y1, x2, y2, x3, y3
      REAL :: a0, a1, b0, b1, a2, b2
      INTEGER :: itest
      INTEGER, PARAMETER :: ntest = 1
      REAL, PARAMETER :: tol = 1e-6

      fail = .FALSE.

      DO itest = 1, ntest

         IF (itest == 1) THEN
            b2 = 4.
            b1 = -3.
            b0 = 0.
            x1 = 0.
            x2 = -3.
            x3 = 1.
         ELSE
            STOP 'TEST_FIX_QUADRATIC: Something went wrong'
         END IF

         y1 = polynomial(x1, b2, b1, b0)
         y2 = polynomial(x2, b2, b1, b0)
         y3 = polynomial(x3, b2, b1, b0)

         CALL fix_polynomial(a2, a1, a0, [x1, x2, x3], [y1, y2, y3])

         IF ((.NOT. requal(a2, b2, tol)) .OR. (.NOT. requal(a1, b1, tol)) .OR. (.NOT. requal(a0, b0, tol))) THEN
            WRITE(*, *) 'TEST_FIX_QUADRATIC: Target a2:', b2
            WRITE(*, *) 'TEST_FIX_QUADRATIC: Target a1:', b1
            WRITE(*, *) 'TEST_FIX_QUADRATIC: Target a0:', b0
            WRITE(*, *) 'TEST_FIX_QUADRATIC: Fix a2:', a2
            WRITE(*, *) 'TEST_FIX_QUADRATIC: Fix a1:', a1
            WRITE(*, *) 'TEST_FIX_QUADRATIC: Fix a0:', a0
            STOP 'TEST_FIX_QUADRATIC: Fail'
         END IF

      END DO

      WRITE(*, *) 'TEST_FIX_QUADRATIC: Pass'
      WRITE(*, *)

   END SUBROUTINE test_fix_quadratic

   SUBROUTINE test_fix_cubic(fail)

      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: fail
      REAL :: x1, x2, x3, x4
      REAL :: y1, y2, y3, y4
      REAL :: a0, a1, a2, a3
      REAL :: b0, b1, b2, b3
      INTEGER :: itest
      INTEGER, PARAMETER :: ntest = 1
      REAL, PARAMETER :: tol = 1e-6

      fail = .FALSE.

      DO itest = 1, ntest

         IF (itest == 1) THEN
            b3 = -2.
            b2 = 4.
            b1 = -3.
            b0 = 0.
            x1 = 0.
            x2 = -3.
            x3 = 1.
            x4 = 6.
         ELSE
            STOP 'TEST_FIX_CUBIC: Something went wrong'
         END IF

         y1 = polynomial(x1, b3, b2, b1, b0)
         y2 = polynomial(x2, b3, b2, b1, b0)
         y3 = polynomial(x3, b3, b2, b1, b0)
         y4 = polynomial(x4, b3, b2, b1, b0)

         CALL fix_polynomial(a3, a2, a1, a0, [x1, x2, x3, x4], [y1, y2, y3, y4])

         IF ((.NOT. requal(a3, b3, tol)) .OR. &
            (.NOT. requal(a2, b2, tol)) .OR. &
            (.NOT. requal(a1, b1, tol)) .OR. &
            (.NOT. requal(a0, b0, tol))) THEN
            WRITE(*, *) 'TEST_FIX_CUBIC: Target a3:', b3
            WRITE(*, *) 'TEST_FIX_CUBIC: Target a2:', b2
            WRITE(*, *) 'TEST_FIX_CUBIC: Target a1:', b1
            WRITE(*, *) 'TEST_FIX_CUBIC: Target a0:', b0
            WRITE(*, *) 'TEST_FIX_CUBIC: Fix a3:', a3
            WRITE(*, *) 'TEST_FIX_CUBIC: Fix a2:', a2
            WRITE(*, *) 'TEST_FIX_CUBIC: Fix a1:', a1
            WRITE(*, *) 'TEST_FIX_CUBIC: Fix a0:', a0
            STOP 'TEST_FIX_CUBIC: Fail'
         END IF

      END DO

      WRITE(*, *) 'TEST_FIX_CUBIC: Pass'
      WRITE(*, *)

   END SUBROUTINE test_fix_cubic

   SUBROUTINE test_factorial(fail)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      INTEGER :: itest
      INTEGER :: n, ft, fv
      INTEGER, PARAMETER :: ntest = 4

      fail = .FALSE.

      DO itest = 1, ntest

         IF (itest == 1) THEN
            n = 0
            ft = 1
         ELSE IF (itest == 2) THEN
            n = 1
            ft = 1
         ELSE IF (itest == 3) THEN
            n = 4
            ft = 24
         ELSE IF (itest == 4) THEN
            n = 10
            ft = 3628800
         ELSE
            STOP 'TEST_FACTORIAL: Error, test not specified correctly'
         END IF

         fv = factorial(n)

         IF (ft .NE. fv) THEN
            fail = .TRUE.
            STOP 'TEST_FACTORIAL: Fail'
         END IF

      END DO

      WRITE(*, *) 'TEST_FACTORIAL: Pass'
      WRITE(*, *)

   END SUBROUTINE test_factorial

END PROGRAM special_functions_test