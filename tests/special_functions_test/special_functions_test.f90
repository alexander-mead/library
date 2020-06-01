PROGRAM special_functions_test

   USE basic_operations
   USE special_functions

   IMPLICIT NONE
   LOGICAL :: ifail

   WRITE(*, *)

   CALL test_fix_linear(ifail)

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

         CALL fix_linear(a1, a0, x1, y1, x2, y2)

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

END PROGRAM special_functions_test