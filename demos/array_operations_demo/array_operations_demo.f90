PROGRAM array_operations_test

   USE basic_operations
   USE array_operations

   IMPLICIT NONE
   REAL, ALLOCATABLE :: a(:), b(:), c(:, :, :), d(:)
   REAL :: min, max
   INTEGER :: na, i, nb, j, k, itest, nc, n, ilog, m
   CHARACTER(len=256) :: test
   INTEGER, ALLOCATABLE :: idx(:), match(:), unique(:), loc(:)
   INTEGER, PARAMETER :: ntest = 13

   CALL get_command_argument(1, test)
   IF (test == '') THEN
      itest = -1
   ELSE
      READ (test, *) itest
   END IF

   WRITE (*, *)

   IF (itest == -1) THEN
      WRITE (*, *) 'Array tests'
      WRITE (*, *) '==========='
      WRITE (*, *) ' 1 - Test reverse'
      WRITE (*, *) ' 2 - Test reduce'
      WRITE (*, *) ' 3 - Test splay'
      WRITE (*, *) ' 4 - Test binning'
      WRITE (*, *) ' 5 - Test fill'
      WRITE (*, *) ' 6 - Test amputate'
      WRITE (*, *) ' 7 - Test unique index'
      WRITE (*, *) ' 8 - Test unique entries'
      WRITE (*, *) ' 9 - '
      WRITE (*, *) '10 - Test number of appearances, repeated entries, and locations'
      WRITE (*, *) '11 - Test remove array element'
      WRITE (*, *) '12 - Test remove repeated array elements'
      WRITE (*, *) '13 - Test swap arrays'
      READ (*, *) itest
      WRITE (*, *)
   END IF

   IF (itest == 1) THEN

      WRITE (*, *) 'Dimensions of array to reverse (filled up by i)'
      READ (*, *) na

      ALLOCATE (a(na))

      WRITE (*, *) 'Original'
      DO i = 1, na
         a(i) = i
         WRITE (*, *) a(i)
      END DO

      CALL reverse_array(a)

      WRITE (*, *) 'Reversed'
      DO i = 1, na
         WRITE (*, *) a(i)
      END DO

   ELSE IF (itest == 2) THEN

      WRITE (*, *) 'Original array size'
      READ (*, *) na
      WRITE (*, *) 'Size to be reduced to'
      READ (*, *) nb

      ALLOCATE (a(na), b(nb))

      WRITE (*, *) 'Original'
      DO i = 1, na
         a(i) = i
         WRITE (*, *) a(i)
      END DO

      CALL reduce_array(a, b)

      WRITE (*, *) 'Reduced'
      DO i = 1, nb
         WRITE (*, *) b(i)
      END DO
      WRITE (*, *)

   ELSE IF (itest == 3) THEN

      WRITE (*, *) 'Splay array dimensions (it is a cube)'
      READ (*, *) nc

      ALLOCATE (c(nc, nc, nc))

      WRITE (*, *) 'Original cube'
      DO i = 1, nc
         DO j = 1, nc
            DO k = 1, nc
               c(i, j, k) = i*j*k
               WRITE (*, *) i, j, k, i*j*k
            END DO
         END DO
      END DO

      d = splay(c, nc, nc, nc)

      WRITE (*, *) 'Splayed result'
      DO i = 1, nc**3
         WRITE (*, *) i, d(i)
      END DO
      WRITE (*, *)

   ELSE IF (itest == 4) THEN

      WRITE (*, *) 'Dimension of array to bin'
      READ (*, *) na

      ALLOCATE (a(na))

      DO i = 1, na
         a(i) = exp(-float(i))
      END DO

      WRITE (*, *) 'Number of bins:'
      READ (*, *) n

      ALLOCATE (b(n), d(n))
      CALL binning(a, minval(a), maxval(a), b, d, ilog=0)

      DO i = 1, n
         WRITE (*, *) b(i), d(i)
      END DO

   ELSE IF (itest == 5) THEN

      WRITE (*, *) 'min'
      READ (*, *) min
      WRITE (*, *) 'max'
      READ (*, *) max
      WRITE (*, *) 'number'
      READ (*, *) n
      WRITE (*, *) '0 - Linear fill'
      WRITE (*, *) '1 - Log fill'
      READ (*, *) ilog
      WRITE (*, *)

      IF (ilog == 0) THEN
         CALL fill_array(min, max, a, n)
      ELSE IF (ilog == 1) THEN
         CALL fill_array(log(min), log(max), a, n)
         a = exp(a)
      ELSE
         STOP
      END IF

      DO i = 1, n
         WRITE (*, *) i, a(i)
      END DO
      WRITE (*, *)

   ELSE IF (itest == 6) THEN

      n = 10
      CALL fill_array(1., 10., a, n)

      DO i = 1, n
         WRITE (*, *) i, a(i)
      END DO

      m = 5
      CALL amputate_array(a, 1, m)

      DO i = 1, m
         WRITE (*, *) i, a(i)
      END DO

   ELSE IF (itest == 7) THEN

      n = 5
      ALLOCATE (idx(n), match(n))
      idx(1) = 1
      idx(2) = 1
      idx(3) = 3
      idx(4) = 2
      idx(5) = 3

      CALL unique_index(idx, n, unique, m, match)

      WRITE (*, *) 'Original array'
      DO i = 1, n
         WRITE (*, *) idx(i)
      END DO
      WRITE (*, *)

      WRITE (*, *) 'Unique array'
      DO i = 1, m
         WRITE (*, *) unique(i)
      END DO
      WRITE (*, *)

      WRITE (*, *) 'Matching indices'
      DO i = 1, n
         WRITE (*, *) match(i)
      END DO
      WRITE (*, *)

   ELSE IF (itest == 8) THEN

      DO j = 1, 3

         IF (j == 1) THEN
            n = 5
            ALLOCATE (idx(n))
            idx(1) = 1
            idx(2) = 2
            idx(3) = 4
            idx(4) = 3
            idx(5) = 5
         ELSE IF (j == 2) THEN
            n = 5
            ALLOCATE (idx(n))
            idx(1) = 1
            idx(2) = 1
            idx(3) = 2
            idx(4) = 3
            idx(5) = 2
         ELSE IF (j == 3) THEN
            n = 5
            ALLOCATE (idx(n))
            idx(1) = 2
            idx(2) = 2
            idx(3) = 2
            idx(4) = 2
            idx(5) = 2
         END IF

         WRITE (*, *) 'Array:'
         DO i = 1, n
            WRITE (*, *) idx(i)
         END DO
         WRITE (*, *) 'Unique entries:', unique_entries(idx)
         WRITE (*, *)

         DEALLOCATE (idx)

      END DO

   ELSE IF (itest == 10) THEN

      DO j = 1, 2

         IF (j == 1) THEN
            n = 3
            ALLOCATE (idx(n))
            idx(1) = 3
            idx(2) = 2
            idx(3) = 1
         ELSE IF (j == 2) THEN
            n = 3
            ALLOCATE (idx(n))
            idx(1) = 3
            idx(2) = 2
            idx(3) = 2
         END IF

         WRITE (*, *) 'Array:'
         DO i = 1, n
            WRITE (*, *) idx(i)
         END DO
         WRITE (*, *) 'Repeated entries:', repeated_entries(idx)
         WRITE (*, *) 'Number of apperances of ''2'':', number_of_appearances(2, idx)
         CALL array_positions(2, idx, loc, m)
         WRITE (*, *) 'Locations of ''2'':', (loc(k), k=1, m)
         WRITE (*, *)

         DEALLOCATE (idx, loc)

      END DO

   ELSE IF (itest == 11) THEN

      ! 11 - Remove array element

      n = 3
      ALLOCATE (a(n))
      a(1) = 2.
      a(2) = 3.
      a(3) = 4.

      CALL write_array_list(a)

      CALL remove_array_element(2, a)
      WRITE (*, *) 'Removed second element'
      WRITE (*, *)
      n = n-1

      CALL write_array_list(a)

   ELSE IF (itest == 12) THEN

      ! 12 - Remove repeated entries
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

      CALL write_array_list(a)

      CALL remove_repeated_array_elements(a, m)
      WRITE (*, *) 'Removed repeated entries'
      WRITE (*, *)

      CALL write_array_list(a)

   ELSE IF (itest == 13) THEN

      ! 13 - Test swap arrays

      n = 3
      ALLOCATE (a(n), b(n))
      a(1) = 1.
      a(2) = 2.
      a(3) = 3.
      WRITE (*, *) 'Array a'
      CALL write_array_list(a)

      b(1) = 4.
      b(2) = 5.
      b(3) = 6.
      WRITE (*, *) 'Array b'
      CALL write_array_list(b)

      WRITE (*, *) 'Swapping arrays'
      CALL swap_arrays(a, b)
      WRITE (*, *) 'Arrays swapped'
      WRITE (*, *)

      WRITE (*, *) 'Array a'
      CALL write_array_list(a)

      WRITE (*, *) 'Array b'
      CALL write_array_list(b)

   ELSE

      STOP 'Test not specified correctly'

   END IF

END PROGRAM array_operations_test
