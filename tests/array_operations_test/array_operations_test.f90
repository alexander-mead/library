PROGRAM array_operations_test

   USE array_operations
   USE basic_operations

   IMPLICIT NONE

   LOGICAL :: fail_main
   LOGICAL, PARAMETER :: verbose_main = .FALSE.

   ! Initial white space
   WRITE (*, *)

   ! Tests
   CALL test_reverse_array(fail_main, verbose_main)
   CALL test_swap_arrays(fail_main, verbose_main)
   CALL test_remove_repeated_entries(fail_main, verbose_main)
   CALL test_remove_array_element(fail_main, verbose_main)
   CALL test_repeated_entries(fail_main, verbose_main)
   CALL test_number_of_appearances(fail_main, verbose_main)
   CALL test_array_positions(fail_main, verbose_main)
   CALL test_unique_entries(fail_main, verbose_main)
   CALL test_amputate_array(fail_main, verbose_main)
   CALL test_unique_index(fail_main, verbose_main)
   CALL test_fill_array(fail_main, verbose_main)
   CALL test_splay(fail_main, verbose_main)
   CALL test_reduce_array(fail_main, verbose_main)
   CALL test_regular_spacing(fail_main)
   CALL test_insert_in_array(fail_main)

CONTAINS

   SUBROUTINE test_insert_in_array(fail)

      LOGICAL, INTENT(INOUT) :: fail
      REAL, ALLOCATABLE :: a(:), b(:)
      REAL :: x
      INTEGER :: n

      fail = .FALSE.

      n = 3
      ALLOCATE(a(n), b(n+1))
      a(1) = 1.
      a(2) = 1.2
      a(3) = 2.

      x = 1.1
      CALL insert_in_array(x, 2, a)

      b(1) = 1.
      b(2) = 1.1
      b(3) = 1.2
      b(4) = 2.

      IF(size(a) /= n+1) THEN
         fail = .TRUE.
         STOP 'TEST_INSERT_IN_ARRAY: Error output array is the wrong size'
      END IF

      IF(any(a /= b)) THEN
         fail = .TRUE.
         STOP 'TEST_INSERT_IN_ARRAY: Error output array is wrong'
      END IF

      IF(.NOT. fail) THEN
         WRITE(*, *) 'TEST_INSERT_IN_ARRAY: Pass'
         WRITE(*, *)
      END IF

   END SUBROUTINE test_insert_in_array

   SUBROUTINE test_reverse_array(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:)
      INTEGER :: n = 5

      ALLOCATE (a(n))
      a(1) = 1.
      a(2) = 2.
      a(3) = 3.
      a(4) = 4.
      a(5) = 5.

      IF (verbose) THEN
         WRITE (*, *) 'TEST_REVERSE_ARRAY: Original array:'
         CALL write_array_list(a)
      END IF

      CALL reverse_array(a)

      IF (verbose) THEN
         WRITE (*, *) 'TEST_REVERSE_ARRAY: Reversed array:'
         CALL write_array_list(a)
      END IF

      fail = .FALSE.
      IF (a(1) .NE. 5.) fail = .TRUE.
      IF (a(2) .NE. 4.) fail = .TRUE.
      IF (a(3) .NE. 3.) fail = .TRUE.
      IF (a(4) .NE. 2.) fail = .TRUE.
      IF (a(5) .NE. 1.) fail = .TRUE.

      IF (fail) THEN
         WRITE (*, *) 'TEST_REMOVE_ARRAY_ELEMENT: a:', a
         STOP 'TEST_REVERSE_ARRAY: Failed'
      ELSE
         WRITE (*, *) 'TEST_REVERSE_ARRAY: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_reverse_array

   SUBROUTINE test_swap_arrays(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:), b(:)
      INTEGER :: n = 3

      n = 3
      ALLOCATE (a(n), b(n))
      a(1) = 1.
      a(2) = 2.
      a(3) = 3.

      IF (verbose) THEN
         WRITE (*, *) 'TEST_SWAP_ARRAYS: Array a'
         CALL write_array_list(a)
      END IF

      b(1) = 4.
      b(2) = 5.
      b(3) = 6.

      IF (verbose) THEN
         WRITE (*, *) 'TEST_SWAP_ARRAYS: Array b'
         CALL write_array_list(b)
      END IF

      IF (verbose) WRITE (*, *) 'TEST_SWAP_ARRAYS: Swapping arrays'
      CALL swap_arrays(a, b)
      IF (verbose) THEN
         WRITE (*, *) 'TEST_SWAP_ARRAYS: Arrays swapped'
         WRITE (*, *)
      END IF

      IF (verbose) THEN
         WRITE (*, *) 'TEST_SWAP_ARRAYS: Array a'
         CALL write_array_list(a)
      END IF

      IF (verbose) THEN
         WRITE (*, *) 'TEST_SWAP_ARRAYS: Array b'
         CALL write_array_list(b)
      END IF

      fail = .FALSE.
      IF (a(1) .NE. 4.) fail = .TRUE.
      IF (a(2) .NE. 5.) fail = .TRUE.
      IF (a(3) .NE. 6.) fail = .TRUE.
      IF (b(1) .NE. 1.) fail = .TRUE.
      IF (b(2) .NE. 2.) fail = .TRUE.
      IF (b(3) .NE. 3.) fail = .TRUE.

      IF (fail) THEN
         WRITE (*, *) 'TEST_REMOVE_ARRAY_ELEMENT: a:', a
         WRITE (*, *) 'TEST_REMOVE_ARRAY_ELEMENT: b:', b
         STOP 'TEST_SWAP_ARRAYS: Failed'
      ELSE
         WRITE (*, *) 'TEST_SWAP_ARRAYS: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_swap_arrays

   SUBROUTINE test_remove_repeated_entries(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:)
      INTEGER :: n, m

      n = 9
      ALLOCATE (a(n))
      a(1) = 1.
      a(2) = 1.
      a(3) = 1.
      a(4) = 2.
      a(5) = 3.
      a(6) = 4.
      a(7) = 5.
      a(8) = 6.
      a(9) = 6.

      IF (verbose) CALL write_array_list(a)

      CALL remove_repeated_array_elements(a, m)
      IF (verbose) THEN
         WRITE (*, *) 'TEST_REMOVE_REPEATED_ENTRIES: Removed repeated entries'
         WRITE (*, *)
      END IF

      IF (verbose) CALL write_array_list(a)

      fail = .FALSE.
      IF (m .NE. 6) fail = .TRUE.
      IF (size(a) .NE. 6) fail = .TRUE.
      IF (a(1) .NE. 1.) fail = .TRUE.
      IF (a(2) .NE. 2.) fail = .TRUE.
      IF (a(3) .NE. 3.) fail = .TRUE.
      IF (a(4) .NE. 4.) fail = .TRUE.
      IF (a(5) .NE. 5.) fail = .TRUE.
      IF (a(6) .NE. 6.) fail = .TRUE.

      IF (fail) THEN
         WRITE (*, *) 'TEST_REMOVE_ARRAY_ELEMENT: a:', a
         STOP 'TEST_REMOVE_REPEATED_ENTRIES: Failed'
      ELSE
         WRITE (*, *) 'TEST_REMOVE_REPEATED_ENTRIES: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_remove_repeated_entries

   SUBROUTINE test_remove_array_element(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:)
      INTEGER :: n

      n = 3
      ALLOCATE (a(n))
      a(1) = 2.
      a(2) = 3.
      a(3) = 4.

      IF (verbose) CALL write_array_list(a)

      CALL remove_array_element(2, a)
      n = n-1
      IF (verbose) THEN
         WRITE (*, *) 'TEST_REMOVE_ARRAY_ELEMENT: Removed second element'
         WRITE (*, *)
      END IF

      IF (verbose) CALL write_array_list(a)

      fail = .FALSE.
      IF (size(a) .NE. 2) fail = .TRUE.
      IF (a(1) .NE. 2.) fail = .TRUE.
      IF (a(2) .NE. 4.) fail = .TRUE.

      IF (fail) THEN
         WRITE (*, *) 'TEST_REMOVE_ARRAY_ELEMENT:', a
         STOP 'TEST_REMOVE_ARRAY_ELEMENT: Failed'
      ELSE
         WRITE (*, *) 'TEST_REMOVE_ARRAY_ELEMENT: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_remove_array_element

   SUBROUTINE test_repeated_entries(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      INTEGER, ALLOCATABLE :: i(:)
      INTEGER :: j, n

      fail = .FALSE.

      DO j = 1, 2

         IF (j == 1) THEN
            n = 3
            ALLOCATE (i(n))
            i(1) = 3
            i(2) = 2
            i(3) = 1
         ELSE IF (j == 2) THEN
            n = 3
            ALLOCATE (i(n))
            i(1) = 3
            i(2) = 2
            i(3) = 2
         END IF

         IF (verbose) THEN
            CALL write_array_list(i)
         END IF

         IF (j == 1) THEN
            IF (repeated_entries(i)) fail = .TRUE.
         ELSE IF (j == 2) THEN
            IF (.NOT. repeated_entries(i)) fail = .TRUE.
         ELSE
            STOP 'TEST_REPEATED_ENTIRES: Error, something went wrong'
         END IF

         DEALLOCATE (i)

      END DO

      IF (fail) THEN
         STOP 'TEST_REPEATED_ENTRIES: Failed'
      ELSE
         WRITE (*, *) 'TEST_REPEATED_ENTRIES: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_repeated_entries

   SUBROUTINE test_number_of_appearances(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      INTEGER, ALLOCATABLE :: i(:)
      INTEGER :: j, n

      fail = .FALSE.

      DO j = 1, 2

         IF (j == 1) THEN
            n = 3
            ALLOCATE (i(n))
            i(1) = 3
            i(2) = 2
            i(3) = 1
         ELSE IF (j == 2) THEN
            n = 3
            ALLOCATE (i(n))
            i(1) = 3
            i(2) = 2
            i(3) = 2
         END IF

         IF (verbose) THEN
            CALL write_array_list(i)
         END IF

         IF (verbose) THEN
            WRITE (*, *)
         END IF

         IF (j == 1) THEN
            IF (number_of_appearances(1, i) .NE. 1) fail = .TRUE.
            IF (number_of_appearances(2, i) .NE. 1) fail = .TRUE.
            IF (number_of_appearances(3, i) .NE. 1) fail = .TRUE.
         ELSE IF (j == 2) THEN
            IF (number_of_appearances(1, i) .NE. 0) fail = .TRUE.
            IF (number_of_appearances(2, i) .NE. 2) fail = .TRUE.
            IF (number_of_appearances(3, i) .NE. 1) fail = .TRUE.
         ELSE
            STOP 'TEST_NUMBER_OF_APPEARANCES: Error, something went wrong'
         END IF

         DEALLOCATE (i)

      END DO

      IF (fail) THEN
         STOP 'TEST_NUMBER_OF_APPEARANCES: Failed'
      ELSE
         WRITE (*, *) 'TEST_NUMBER_OF_APPEARANCES: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_number_of_appearances

   SUBROUTINE test_array_positions(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      INTEGER, ALLOCATABLE :: i(:), loc(:)
      INTEGER :: j, k, n, m

      fail = .FALSE.

      DO j = 1, 2

         IF (j == 1) THEN
            n = 3
            ALLOCATE (i(n))
            i(1) = 3
            i(2) = 2
            i(3) = 1
         ELSE IF (j == 2) THEN
            n = 3
            ALLOCATE (i(n))
            i(1) = 3
            i(2) = 2
            i(3) = 2
         ELSE
            STOP 'TEST_ARRAY_POSITIONS: Error, something went wrong'
         END IF

         CALL array_positions(2, i, loc, m)

         IF (verbose) THEN
            CALL write_array_list(i)
            WRITE (*, *) 'TEST_ARRAY_POSITIONS: Subtest:', j
            WRITE (*, *) 'TEST_ARRAY_POSITIONS: Locations of ''2'':', (loc(k), k=1, m)
            WRITE (*, *)
         END IF

         IF (j == 1) THEN
            IF (m .NE. 1) fail = .TRUE.
            IF (size(loc) .NE. 1) fail = .TRUE.
            IF (loc(1) .NE. 2) fail = .TRUE.
         ELSE IF (j == 2) THEN
            IF (m .NE. 2) fail = .TRUE.
            IF (size(loc) .NE. 2) fail = .TRUE.
            IF (loc(1) .NE. 2) fail = .TRUE.
            IF (loc(2) .NE. 3) fail = .TRUE.
         ELSE
            STOP 'TEST_ARRAY_POSITIONS: Error, something went wrong'
         END IF

         DEALLOCATE (i, loc)

      END DO

      IF (fail) THEN
         STOP 'TEST_ARRAY_POSITIONS: Failed'
      ELSE
         WRITE (*, *) 'TEST_ARRAY_POSITIONS: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_array_positions

   SUBROUTINE test_unique_entries(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      INTEGER, ALLOCATABLE :: i(:)
      INTEGER :: j, n

      fail = .FALSE.

      DO j = 1, 3

         IF (j == 1) THEN
            n = 5
            ALLOCATE (i(n))
            i(1) = 1
            i(2) = 2
            i(3) = 4
            i(4) = 3
            i(5) = 5
         ELSE IF (j == 2) THEN
            n = 5
            ALLOCATE (i(n))
            i(1) = 1
            i(2) = 1
            i(3) = 2
            i(4) = 3
            i(5) = 2
         ELSE IF (j == 3) THEN
            n = 5
            ALLOCATE (i(n))
            i(1) = 2
            i(2) = 2
            i(3) = 2
            i(4) = 2
            i(5) = 2
         ELSE
            STOP 'TEST_UNIQUE_ENTRIES: Error, something went wrong'
         END IF

         IF (verbose) THEN
            CALL write_array_list(i)
            WRITE (*, *) 'TEST_UNIQUE_ENTRIES: Unique entries:', unique_entries(i)
            WRITE (*, *)
         END IF

         IF (j == 1) THEN
            IF (unique_entries(i) .NE. 5) fail = .TRUE.
         ELSE IF (j == 2) THEN
            IF (unique_entries(i) .NE. 3) fail = .TRUE.
         ELSE IF (j == 3) THEN
            IF (unique_entries(i) .NE. 1) fail = .TRUE.
         ELSE
            STOP 'TEST_UNIQUE_ENTRIES: Error, something went wrong'
         END IF

         DEALLOCATE (i)

      END DO

      IF (fail) THEN
         STOP 'TEST_UNIQUE_ENTRIES: Failed'
      ELSE
         WRITE (*, *) 'TEST_UNIQUE_ENTRIES: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_unique_entries

   SUBROUTINE test_amputate_array(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:)
      INTEGER :: n

      n = 5
      ALLOCATE (a(n))
      a(1) = 1.
      a(2) = 2.
      a(3) = 3.
      a(4) = 4.
      a(5) = 5.

      IF (verbose) CALL write_array_list(a)

      CALL amputate_array(a, 2, 4)

      IF (verbose) CALL write_array_list(a)

      fail = .FALSE.
      IF (size(a) .NE. 3) fail = .TRUE.
      IF (a(1) .NE. 2.) fail = .TRUE.
      IF (a(2) .NE. 3.) fail = .TRUE.
      IF (a(3) .NE. 4.) fail = .TRUE.

      IF (fail) THEN
         STOP 'TEST_AMPUTATE_ARRAY: Failed'
      ELSE
         WRITE (*, *) 'TEST_AMPUTATE_ARRAY: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_amputate_array

   SUBROUTINE test_unique_index(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      INTEGER, ALLOCATABLE :: i(:), match(:), unique(:)
      INTEGER :: n, m, j

      n = 5
      ALLOCATE (i(n), match(n))
      i(1) = 1 ! First unique entry
      i(2) = 1
      i(3) = 3 ! Second unique entry
      i(4) = 2 ! Third unique entry
      i(5) = 3

      CALL unique_index(i, n, unique, m, match)

      IF (verbose) THEN
         WRITE (*, *) 'Original array'
         CALL write_array_list(i)
         WRITE (*, *) 'Unique array'
         CALL write_array_list(unique)
         WRITE (*, *) 'Matching indices'
         CALL write_array_list(match)
      END IF

      fail = .FALSE.
      IF (m .NE. 3) fail = .TRUE.
      IF (size(unique) .NE. 3) fail = .TRUE.
      IF (unique(1) .NE. 1) fail = .TRUE. ! First unique entry
      IF (unique(2) .NE. 3) fail = .TRUE. ! Second unique entry
      IF (unique(3) .NE. 2) fail = .TRUE. ! Third unique entry
      IF (match(1) .NE. 1) fail = .TRUE.  ! Should always be one
      IF (match(2) .NE. 1) fail = .TRUE.  ! Second entry is the same as first
      IF (match(3) .NE. 2) fail = .TRUE.  ! Third entry the second unique entry
      IF (match(4) .NE. 3) fail = .TRUE.  ! Fourth entry is the third unique entry
      IF (match(5) .NE. 2) fail = .TRUE.  ! Fifth entry is the same as second
      DO j = 1, n
         IF (unique(match(j)) .NE. i(j)) fail = .TRUE. ! This is what the routine is supposed to accomplish
      END DO

      IF (fail) THEN
         STOP 'TEST_UNIQUE_INDEX: Failed'
      ELSE
         WRITE (*, *) 'TEST_UNIQUE_INDEX: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_unique_index

   SUBROUTINE test_fill_array(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:)
      REAL :: min, max
      INTEGER :: n
      REAL, PARAMETER :: eps = 1e-6 ! Tolerance for real-number equal

      min = 1.
      max = 10.
      n = 10
      CALL fill_array(min, max, a, n)

      IF (verbose) CALL write_array_list(a)

      fail = .FALSE.
      IF (.NOT. requal(a(1), min, eps)) fail = .TRUE.
      IF (.NOT. requal(a(n), max, eps)) fail = .TRUE.
      IF (.NOT. requal(a(1), 1., eps)) fail = .TRUE.
      IF (.NOT. requal(a(2), 2., eps)) fail = .TRUE.
      IF (.NOT. requal(a(3), 3., eps)) fail = .TRUE.
      IF (.NOT. requal(a(4), 4., eps)) fail = .TRUE.
      IF (.NOT. requal(a(5), 5., eps)) fail = .TRUE.
      IF (.NOT. requal(a(6), 6., eps)) fail = .TRUE.
      IF (.NOT. requal(a(7), 7., eps)) fail = .TRUE.
      IF (.NOT. requal(a(8), 8., eps)) fail = .TRUE.
      IF (.NOT. requal(a(9), 9., eps)) fail = .TRUE.
      IF (.NOT. requal(a(10), 10., eps)) fail = .TRUE.

      IF (fail) THEN
         STOP 'TEST_FILL_ARRAY: Failed'
      ELSE
         WRITE (*, *) 'TEST_FILL_ARRAY: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_fill_array

   SUBROUTINE test_splay(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:), c(:, :, :)
      INTEGER :: i, j, k, n, ii

      n = 2
      ALLOCATE (c(n, n, n))

      IF (verbose) WRITE (*, *) 'Original cube'
      ii = 0
      DO i = 1, n
         DO j = 1, n
            DO k = 1, n
               ii = ii+1
               c(i, j, k) = real(ii)
               IF (verbose) WRITE (*, *) i, j, k, c(i, j, k)
            END DO
         END DO
      END DO

      a = splay(c, n, n, n)

      IF (verbose) THEN
         WRITE (*, *) 'Splayed result'
         CALL write_array_list(a)
      END IF

      fail = .FALSE.
      DO i = 1, n
         IF (a(i) .NE. real(i)) fail = .TRUE.
      END DO

      IF (fail) THEN
         STOP 'TEST_SPLAY: Failed'
      ELSE
         WRITE (*, *) 'TEST_SPLAY: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_splay

   SUBROUTINE test_reduce_array(fail, verbose)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:), b(:)
      INTEGER :: i, n, m

      n = 5
      m = 3
      ALLOCATE (a(n), b(m))
      DO i = 1, n
         a(i) = real(i)
      END DO

      IF (verbose) THEN
         WRITE (*, *) 'Original'
         CALL write_array_list(a)
      END IF

      CALL reduce_array(a, b)

      IF (verbose) THEN
         WRITE (*, *) 'Reduced'
         CALL write_array_list(b)
      END IF

      fail = .FALSE.
      IF (b(1) .NE. a(1)) fail = .TRUE.
      IF (b(m) .NE. a(n)) fail = .TRUE.
      IF (b(1) .NE. 1.) fail = .TRUE.
      IF (b(2) .NE. 3.) fail = .TRUE.
      IF (b(3) .NE. 5.) fail = .TRUE.

      IF (fail) THEN
         STOP 'TEST_REDUCE_ARRAY: Failed'
      ELSE
         WRITE (*, *) 'TEST_REDUCE_ARRAY: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_reduce_array

   SUBROUTINE test_regular_spacing(fail)

      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: fail
      REAL, ALLOCATABLE :: a(:)
      INTEGER :: itest
      LOGICAL :: result
      INTEGER, PARAMETER :: n = 5
      INTEGER, PARAMETER :: ntest = 2

      DO itest = 1, ntest

         ALLOCATE (a(n))
         IF (itest == 1) THEN
            a(1) = 1.
            a(2) = 2.
            a(3) = 3.
            a(4) = 4.
            a(5) = 5.
            result = .TRUE.
         ELSE IF (itest == 2) THEN
            a(1) = 1.
            a(2) = 2.
            a(3) = 3.1
            a(4) = 4.
            a(5) = 5.
            result = .FALSE.
         ELSE
            STOP 'TEST_REGULAR_SPACING: Error, itest specified incorrecly'
         END IF

         IF (result .NEQV. regular_spacing(a)) THEN
            WRITE (*, *) 'TEST_REGULAR_SPACING: Array:', a
            WRITE (*, *) 'TEST_REGULAR_SPACING: Expected result:', result
            fail = .TRUE.
            STOP 'TEST_REGULAR_SPACING: Test failed'
         END IF

         DEALLOCATE (a)

      END DO

      IF (.NOT. fail) THEN
         WRITE (*, *) 'TEST_REGULAR_SPACING: Pass'
         WRITE (*, *)
      END IF

   END SUBROUTINE test_regular_spacing

END PROGRAM array_operations_test
