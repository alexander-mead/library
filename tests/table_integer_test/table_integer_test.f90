PROGRAM table_integer_test

   USE table_integer

   IMPLICIT NONE

   WRITE(*, *)

   CALL test_table_integer()

contains

   SUBROUTINE test_table_integer()

      REAL :: x(5) = [1., 2., 3., 4., 5.] ! Neeeds to be linearly spaced for ifind=1 to work
      REAL :: xfind
      INTEGER :: i, itrue, itest, ifind
      LOGICAL :: fail
      INTEGER, PARAMETER :: nfind = 3
      INTEGER, PARAMETER :: ntest = 5

      fail = .FALSE.

      DO itest = 1, ntest

         IF(itest == 1) THEN
            xfind = 0.
            itrue = 0
         ELSE IF(itest == 2) THEN
            xfind = 105.1
            itrue = 5
         ELSE IF(itest == 3) THEN
            xfind = 1.4
            itrue = 1
         ELSE IF(itest == 4) THEN
            xfind = 1.1
            itrue = 1
         ELSE IF(itest == 5) THEN
            xfind = 4.1
            itrue = 4
         ELSE
            STOP 'TEST_TABLE_INTEGER: Test not recognised'
         END IF

         DO ifind = 1, nfind
            i = find_table_integer(xfind, x, ifind)
            IF (i .NE. itrue) THEN
               fail = .TRUE.
               WRITE(*, *) 'TEST_TABLE_INTEGER: Test:', itest
               WRITE(*, *) 'TEST_TABLE_INTEGER: ifind:', ifind
               WRITE(*, *) 'TEST_TABLE_INTEGER: x:', xfind
               WRITE(*, *) 'TEST_TABLE_INTEGER: Result:', i
               WRITE(*, *) 'TEST_TABLE_INTEGER: Truth:', itrue
               STOP 'TEST_TABLE_INTEGER: Failed'
            END IF
         END DO

      END DO

      WRITE(*, *) 'TEST_TABLE_INTEGER: Pass'
      WRITE(*, *)

   END SUBROUTINE test_table_integer

END PROGRAM table_integer_test