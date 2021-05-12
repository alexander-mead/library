MODULE array_operations

   USE basic_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: insert_in_array
   PUBLIC :: within_array
   PUBLIC :: swap_arrays
   PUBLIC :: append
   PUBLIC :: add_to_array
   PUBLIC :: splay
   PUBLIC :: write_array_list
   PUBLIC :: repeated_entries
   PUBLIC :: array_position
   PUBLIC :: reverse_array
   PUBLIC :: remove_array_element
   PUBLIC :: array_positions
   PUBLIC :: amputate_array
   PUBLIC :: unique_index
   PUBLIC :: reduce_array
   PUBLIC :: fill_pixels
   PUBLIC :: binning
   PUBLIC :: merge_arrays
   PUBLIC :: create_mask
   PUBLIC :: apply_mask
   PUBLIC :: smooth_array_tophat
   PUBLIC :: smooth_array_Gaussian
   PUBLIC :: if_allocated_deallocate
   PUBLIC :: regular_spacing
   PUBLIC :: safe_allocate
   PUBLIC :: greater_than_all
   PUBLIC :: greater_than_any
   PUBLIC :: is_in_array
   PUBLIC :: integer_sequence
   PUBLIC :: unique_entries
   PUBLIC :: number_of_appearances
   PUBLIC :: remove_repeated_array_elements
   PUBLIC :: remove_repeated_two_array_elements
   PUBLIC :: fill_array
   PUBLIC :: fill_array_log ! TODO: Remove
   PUBLIC :: outside_array_range
   PUBLIC :: sum_logical
   PUBLIC :: monotonic_array

   INTERFACE is_in_array
      MODULE PROCEDURE is_in_array_integer
      MODULE PROCEDURE is_in_array_character
   END INTERFACE is_in_array

   INTERFACE array_position
      MODULE PROCEDURE array_position_int
      MODULE PROCEDURE array_position_real
   END INTERFACE array_position

   INTERFACE add_to_array
      MODULE PROCEDURE add_to_array_2D
      MODULE PROCEDURE add_to_array_3D
   END INTERFACE add_to_array

   INTERFACE reverse_array
      MODULE PROCEDURE reverse_array_1D
      MODULE PROCEDURE reverse_array_2D
      MODULE PROCEDURE reverse_array_3D
   END INTERFACE reverse_array

   INTERFACE splay
      PROCEDURE splay_2D
      PROCEDURE splay_3D
   END INTERFACE splay

   INTERFACE write_array_list
      PROCEDURE write_array_list_real
      PROCEDURE write_2array_list_real
      PROCEDURE write_array2D_list_real
      PROCEDURE write_array_list_int
   END INTERFACE write_array_list

   INTERFACE apply_mask
      PROCEDURE apply_mask_1D
      PROCEDURE apply_mask_2D
   END INTERFACE apply_mask

   INTERFACE if_allocated_deallocate
      PROCEDURE if_allocated_deallocate_real_1D
      PROCEDURE if_allocated_deallocate_real_2D
      PROCEDURE if_allocated_deallocate_real_3D
      PROCEDURE if_allocated_deallocate_integer_1D
      PROCEDURE if_allocated_deallocate_logical_1D
      PROCEDURE if_allocated_deallocate_character_1D
   END INTERFACE if_allocated_deallocate

   INTERFACE safe_allocate
      PROCEDURE safe_allocate_real
      PROCEDURE safe_allocate2_real
      PROCEDURE safe_allocate_integer
      PROCEDURE safe_allocate_logical
      PROCEDURE safe_allocate_character
   END INTERFACE safe_allocate

   INTERFACE remove_array_element
      PROCEDURE remove_array_element_real
      PROCEDURE remove_array_element_int
   END INTERFACE remove_array_element

   INTERFACE amputate_array
      PROCEDURE amputate_array_real
      PROCEDURE amputate_array_int
   END INTERFACE amputate_array

   INTERFACE concatenate_arrays
      PROCEDURE concatenate_arrays_real
      PROCEDURE concatenate_arrays_int
   END INTERFACE concatenate_arrays

CONTAINS

   LOGICAL FUNCTION monotonic_array(x)

      ! Returns true if sequential array entries are monotonically increasing (or equal) from low to high
      REAL, INTENT(IN) :: x(:)
      INTEGER :: i

      monotonic_array = .TRUE.
      DO i = 1, size(x)-1
         IF (x(i+1) < x(i)) THEN
            monotonic_array = .FALSE.
            EXIT
         END IF
      END DO

   END FUNCTION monotonic_array

   LOGICAL FUNCTION outside_array_range(x, xtab)

      ! Assumes xtab is an array of x values going from low to high
      ! Returns TRUE if x is outside the range of xtab (either low or high) and FALSE if within
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xtab(:)

      IF (x < xtab(1) .OR. x > xtab(size(xtab))) THEN
         outside_array_range = .TRUE.
      ELSE
         outside_array_range = .FALSE.
      END IF

   END FUNCTION outside_array_range

   SUBROUTINE insert_in_array(x, i, a)

      REAL, INTENT(IN) :: x                    ! Value to insert
      INTEGER, INTENT(IN) :: i                 ! Position to insert
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Array in which to insert
      REAL, ALLOCATABLE :: b(:)
      INTEGER :: j, n
      
      n = size(a)
      IF(i < 1 .OR. i > n) THEN
         STOP 'INSERT_IN_ARRAY: Error, integer is outside array'
      END IF
      ALLOCATE(b(n))
      b = a
      DEALLOCATE(a)
      ALLOCATE(a(n+1))

      DO j = 1, n+1
         IF (j < i) THEN
            a(j) = b(j)
         ELSE IF (j == i) THEN
            a(j) = x
         ELSE IF (j > i) THEN
            a(j) = b(j-1)
         END IF
      END DO

   END SUBROUTINE insert_in_array

   LOGICAL FUNCTION greater_than_any(x, a)

      ! Returns TRUE if x is greater than any one of the values in array
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: a(:)
      INTEGER :: i, n

      n = size(a)

      STOP 'GREATER_THAN_ANY: Test this'

      greater_than_any = .FALSE.
      DO i = 1, n
         IF(x > a(i)) THEN
            greater_than_any = .TRUE.
            EXIT
         END IF
      END DO

   END FUNCTION greater_than_any

   LOGICAL FUNCTION greater_than_all(x, a)

      ! Returns TRUE if x is greater than all of the values in array
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: a(:)
      INTEGER :: i, n

      STOP 'GREATER_THAN_ALL: Test this'

      greater_than_all = .TRUE.
      DO i = 1, n
         IF(x <= a(i)) THEN
            greater_than_all = .FALSE.
            EXIT
         END IF
      END DO

   END FUNCTION greater_than_all

   SUBROUTINE safe_allocate_real(x, n)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(x)
      ALLOCATE(x(n))

   END SUBROUTINE safe_allocate_real

   SUBROUTINE safe_allocate2_real(x, nx, ny)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:, :)
      INTEGER, INTENT(IN) :: nx
      INTEGER, INTENT(IN) :: ny

      CALL if_allocated_deallocate(x)
      ALLOCATE(x(nx, ny))

   END SUBROUTINE safe_allocate2_real

   SUBROUTINE safe_allocate_integer(i, n)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: i(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(i)
      ALLOCATE(i(n))

   END SUBROUTINE safe_allocate_integer

   SUBROUTINE safe_allocate_logical(l, n)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      LOGICAL, ALLOCATABLE, INTENT(INOUT) :: l(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(l)
      ALLOCATE(l(n))

   END SUBROUTINE safe_allocate_logical

   SUBROUTINE safe_allocate_character(c, n)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      CHARACTER(len=*), ALLOCATABLE, INTENT(INOUT) :: c(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(c)
      ALLOCATE(c(n))

   END SUBROUTINE safe_allocate_character

   SUBROUTINE if_allocated_deallocate_real_1D(x)

      ! Deallocates an array if it is already allocated
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)

      IF(allocated(x)) DEALLOCATE(x)

   END SUBROUTINE if_allocated_deallocate_real_1D

   SUBROUTINE if_allocated_deallocate_real_2D(x)

      ! Deallocates an array if it is already allocated
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:, :)

      IF(allocated(x)) DEALLOCATE(x)

   END SUBROUTINE if_allocated_deallocate_real_2D

   SUBROUTINE if_allocated_deallocate_real_3D(x)

      ! Deallocates an array if it is already allocated
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:, :, :)

      IF(allocated(x)) DEALLOCATE(x)

   END SUBROUTINE if_allocated_deallocate_real_3D

   SUBROUTINE if_allocated_deallocate_integer_1D(i)

      ! Deallocates an array if it is already allocated
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: i(:)

      IF(allocated(i)) DEALLOCATE(i)

   END SUBROUTINE if_allocated_deallocate_integer_1D

   SUBROUTINE if_allocated_deallocate_logical_1D(l)

      ! Deallocates an array if it is already allocated
      LOGICAL, ALLOCATABLE, INTENT(INOUT) :: l(:)

      IF(allocated(l)) DEALLOCATE(l)

   END SUBROUTINE if_allocated_deallocate_logical_1D

   SUBROUTINE if_allocated_deallocate_character_1D(c)

      ! Deallocates an array if it is already allocated
      CHARACTER(len=*), ALLOCATABLE, INTENT(INOUT) :: c(:)

      IF(allocated(c)) DEALLOCATE(c)

   END SUBROUTINE if_allocated_deallocate_character_1D

   SUBROUTINE smooth_array_tophat(x, f, dx)

      ! Take an array f(x) and smooth it over a fixed width in x
      REAL, INTENT(IN) :: x(:)    ! x coordinates
      REAL, INTENT(INOUT) :: f(:) ! Array to smooth
      REAL, INTENT(IN) :: dx      ! Half-width of top hat
      INTEGER :: i, j, n
      REAL, ALLOCATABLE :: ff(:)
      REAL :: total, weight

      IF (dx .NE. 0.) THEN

         n = size(x)
         IF (n /= size(f)) STOP 'GAUSSIAN_SMOOTH_ARRAY: Error, x and y should be the same size'

         ! Save the original input array
         ff = f

         ! Delete the original array
         f = 0.

         ! Apply top-hat smoothing
         DO i = 1, n
            total = 0.
            IF (abs(x(i)-x(1)) < dx .OR. abs(x(i)-x(n)) < dx) THEN
               f(i) = ff(i)
            ELSE
               DO j = 1, n
                  IF (abs(x(i)-x(j)) < dx) THEN
                     weight = 1.
                  ELSE
                     weight = 0.
                  END IF
                  f(i) = f(i)+ff(j)*weight
                  total = total+weight
               END DO
               f(i) = f(i)/total
            END IF
         END DO

         DEALLOCATE(ff)

      END IF

   END SUBROUTINE smooth_array_tophat

   SUBROUTINE smooth_array_Gaussian(x, f, sigma)

      ! Smooth an array f(x) using a Gaussian kernel
      REAL, INTENT(IN) :: x(:)    ! x coordinates
      REAL, INTENT(INOUT) :: f(:) ! Array to smooth
      REAL, INTENT(IN) :: sigma   ! Width of smoothing Gaussian
      INTEGER :: i, j, n
      REAL, ALLOCATABLE :: ff(:)
      REAL :: weight, total
      REAL, PARAMETER :: nsig = 3. ! Do not smooth if point lies within this number of sigma from edge

      IF (sigma  .NE. 0.) THEN

         n = size(x)
         IF (n /= size(f)) STOP 'GAUSSIAN_SMOOTH_ARRAY: Error, x and y should be the same size'

         ! Save the original input array
         ff = f

         ! Delete the original array
         f = 0.

         ! Apply Gaussian smoothing
         DO i = 1, n
            total = 0.
            IF (abs(x(i)-x(1)) < nsig*sigma .OR. abs(x(i)-x(n)) < nsig*sigma) THEN
               f(i) = ff(i)
            ELSE
               DO j = 1, n
                  weight = exp(-(x(i)-x(j))**2/(2.*sigma**2))
                  f(i) = f(i)+ff(j)*weight
                  total = total+weight
               END DO
               f(i) = f(i)/total
            END IF
         END DO

         DEALLOCATE(ff)

      END IF

   END SUBROUTINE smooth_array_Gaussian

   LOGICAL FUNCTION within_array(x, a)

      ! Checks to see if x falls within the range of values in array a
      REAL, INTENT(IN) :: x    ! Value to check
      REAL, INTENT(IN) :: a(:) ! Array (of x values, presumably)
      INTEGER :: n

      n = size(a)

      IF (x >= minval(a) .AND. x <= maxval(a)) THEN
         within_array = .TRUE.
      ELSE
         within_array = .FALSE.
      END IF

   END FUNCTION within_array

   SUBROUTINE swap_arrays(x, y)

      ! Swap arrays x<->y in a memory-efficient way
      ! Only one excess real number is ever stored
      REAL, INTENT(INOUT) :: x(:) ! Array 1
      REAL, INTENT(INOUT) :: y(:) ! Array 2
      INTEGER :: i, n
      REAL :: h

      n = size(x)
      IF (n /= size(y)) STOP 'SWAP_ARRAYS: Error, input arrays are different sizes'

      ! Loop over array and swap element by element
      DO i = 1, n
         h = x(i)
         x(i) = y(i)
         y(i) = h
      END DO

   END SUBROUTINE swap_arrays

   SUBROUTINE append(x, a)

      ! Append value b to the end of array a(n) to make a new array a(n+1)
      INTEGER, INTENT(IN) :: x
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:)    
      INTEGER, ALLOCATABLE :: hold(:)
      INTEGER :: i, n

      n = size(a)
      ALLOCATE(hold(n))

      hold = a
      DEALLOCATE (a)
      ALLOCATE (a(n+1))

      DO i = 1, n
         a(i) = hold(i)
      END DO
      a(n+1) = x

   END SUBROUTINE append

   SUBROUTINE add_to_array_2D(a, v, i)

      ! Add value 'v' to array 'a' at location 'i' in array
      ! If 'i' is outside the array range then this routine does nothing
      REAL, INTENT(INOUT) :: a(:, :) 
      REAL, INTENT(IN) :: v
      INTEGER, INTENT(IN) :: i(2)
      INTEGER :: j
      LOGICAL :: bin
      INTEGER, PARAMETER :: dim = 2

      bin = .TRUE.
      DO j = 1, dim
         IF (i(j) < 1 .OR. i(j) > size(a, j)) THEN
            bin = .FALSE.
            EXIT
         END IF
      END DO

      IF (bin) a(i(1), i(2)) = a(i(1), i(2))+v

   END SUBROUTINE add_to_array_2D

   SUBROUTINE add_to_array_3D(a, v, i)

      ! Add value 'v' to array 'a' at location 'i' in array
      ! If 'i' is outside the array range then this routine does nothing
      REAL, INTENT(INOUT) :: a(:, :, :)
      REAL, INTENT(IN) :: v
      INTEGER, INTENT(IN) :: i(3)
      INTEGER :: j
      LOGICAL :: bin
      INTEGER, PARAMETER :: dim = 3

      bin = .TRUE.
      DO j = 1, dim
         IF (i(j) < 1 .OR. i(j) > size(a, j)) THEN
            bin = .FALSE.
            EXIT
         END IF
      END DO

      IF (bin) a(i(1), i(2), i(3)) = a(i(1), i(2), i(3))+v

   END SUBROUTINE add_to_array_3D

   INTEGER FUNCTION array_position_int(x, a)

      ! Returns the location in a(n) of value x
      ! If x is not in array then returns zero
      INTEGER, INTENT(IN) :: x    ! Value to check if it is in array
      INTEGER, INTENT(IN) :: a(:) ! Array to check
      INTEGER :: i, n

      n = size(a)

      array_position_int = 0

      DO i = 1, n
         IF (a(i) == x) THEN
            array_position_int = i
            EXIT
         END IF
      END DO

   END FUNCTION array_position_int

   INTEGER FUNCTION array_position_real(x, a, eps)

      ! Returns the location in a(n) of value x
      ! If x is not in array then returns zero
      REAL, INTENT(IN) :: x    ! Value to check if it is in array
      REAL, INTENT(IN) :: a(:) ! Array to check
      REAL, INTENT(IN) :: eps  ! Difference to tolerate
      INTEGER :: i, n

      n = size(a)

      array_position_real = 0

      DO i = 1, n
         IF(requal(x, a(i), eps)) THEN
            array_position_real = i
            EXIT
         END IF
      END DO

   END FUNCTION array_position_real

   INTEGER FUNCTION number_of_appearances(x, a)

      ! Returns the number of appearances in a(n) of value x
      ! If x is not in a(n) then returns zero
      INTEGER, INTENT(IN) :: x    ! Value to check if it is in array
      INTEGER, INTENT(IN) :: a(:) ! Array to check
      INTEGER :: i, n

      n = size(a)

      number_of_appearances = 0

      DO i = 1, n
         IF (a(i) == x) THEN
            number_of_appearances = number_of_appearances+1
         END IF
      END DO

   END FUNCTION number_of_appearances

   SUBROUTINE array_positions(x, a, b, m)

      ! Returns the locations in the array of integer value x
      ! If x is not in array then returns zero
      INTEGER, INTENT(IN) :: x    ! Value to check if it is in array
      INTEGER, INTENT(IN) :: a(:) ! Array to check
      INTEGER, ALLOCATABLE, INTENT(OUT) :: b(:)
      INTEGER, INTENT(OUT) :: m
      INTEGER :: i, p, n

      n = size(a)

      m = number_of_appearances(x, a)
      ALLOCATE (b(m))

      p = 0
      DO i = 1, n
         IF (a(i) == x) THEN
            p = p+1
            b(p) = i
         END IF
      END DO

   END SUBROUTINE array_positions

   SUBROUTINE amputate_array_real(a, i1, i2)

      ! Chop an array of size a(n) down to a smaller size demarked by indices i1, i2
      ! If i1=1 and i2=n then this does nothing
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: i1
      INTEGER, INTENT(IN) :: i2
      REAL, ALLOCATABLE :: b(:)
      INTEGER :: i, m, n

      IF (.NOT. allocated(a)) STOP 'AMPUTATE_ARRAY: Error, array a must be initially allocated'
      n = size(a)

      IF (i2 < i1) THEN
         STOP 'AMPUTATE: Error, i2 should be greater than i1'
      END IF

      m = i2-i1+1
      IF (n < m) THEN
         STOP 'AMPUTATE: Error, new array should be smaller than the old one'
      END IF

      ! Store input array and then deallocate
      ALLOCATE (b(n))
      b = a
      DEALLOCATE (a)

      ! Allocate new output array
      ALLOCATE (a(m))

      ! Fill new output array
      DO i = 1, m
         a(i) = b(i+i1-1)
      END DO

      ! Deallocate holding array
      DEALLOCATE (b)

   END SUBROUTINE amputate_array_real

   SUBROUTINE amputate_array_int(a, i1, i2)

      ! Chop an array of size a(n) down to a smaller size demarked by indices i1, i2
      ! If i1=1 and i2=n then this does nothing
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: i1
      INTEGER, INTENT(IN) :: i2
      INTEGER, ALLOCATABLE :: b(:)
      INTEGER :: i, m, n

      IF (.NOT. allocated(a)) STOP 'AMPUTATE_ARRAY: Error, array a must be initially allocated'
      n = size(a)

      IF (i2 < i1) THEN
         STOP 'AMPUTATE: Error, i2 should be greater than i1'
      END IF

      m = i2-i1+1
      IF (n < m) THEN
         STOP 'AMPUTATE: Error, new array should be smaller than the old one'
      END IF

      ! Store input array and then deallocate
      ALLOCATE (b(n))
      b = a
      DEALLOCATE (a)

      ! Allocate new output array
      ALLOCATE (a(m))

      ! Fill new output array
      DO i = 1, m
         a(i) = b(i+i1-1)
      END DO

      ! Deallocate holding array
      DEALLOCATE (b)

   END SUBROUTINE amputate_array_int

   SUBROUTINE reduce_array(a1, a2)

      ! Reduces the size of array1 to the size of array2
      ! This will not preserve the spacing of entries in array1 and might be a terrible idea in many cases
      ! TODO: Remove this  
      REAL, INTENT(IN) :: a1(:)
      REAL, INTENT(OUT) :: a2(:)
      INTEGER :: i, j, n1, n2

      n1 = size(a1)
      n2 = size(a2)

      IF (n2 > n1) STOP 'REDUCE_ARRAY: Error, n2 > n1'

      DO i = 1, n2
         j = 1+ceiling(real((n1-1)*(i-1))/real(n2-1))
         a2(i) = a1(j)
      END DO

   END SUBROUTINE reduce_array

   SUBROUTINE reverse_array_1D(a)

      ! Reverses the contents of array
      REAL, INTENT(INOUT) :: a(:)
      INTEGER :: i, n
      REAL, ALLOCATABLE :: hold(:)

      n = size(a)

      hold = a

      DO i = 1, n
         a(i) = hold(n-i+1)
      END DO

   END SUBROUTINE reverse_array_1D

   SUBROUTINE reverse_array_2D(a, dim)

      ! Reverses the contents of array
      REAL, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: dim
      INTEGER :: ix, iy
      INTEGER :: nx, ny
      REAL, ALLOCATABLE :: b(:, :)

      IF ((dim /= 1) .AND. (dim /= 2)) THEN
         STOP 'REVERSE_ARRAY_2D: Error, dim must be either 1 or 2'
      END IF

      b = a

      nx = size(a, 1)
      ny = size(a, 2)

      DO iy = 1, ny
         DO ix = 1, nx
            IF (dim == 1) THEN
               a(ix, iy) = b(nx-ix+1, iy)
            ELSE IF (dim == 2) THEN
               a(ix, iy) = b(ix, ny-iy+1)
            ELSE
               STOP 'REVERSE_ARRAY_2D: Error, something went wrong'
            END IF
         END DO
      END DO

   END SUBROUTINE reverse_array_2D

   SUBROUTINE reverse_array_3D(a, dim)

      ! Reverses the contents of array
      REAL, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: dim
      INTEGER :: ix, iy, iz
      INTEGER :: nx, ny, nz
      REAL, ALLOCATABLE :: b(:, :, :)

      IF ((dim /= 1) .AND. (dim /= 2) .AND. (dim /= 3)) THEN
         STOP 'REVERSE_ARRAY_3D: Error, dim must be either 1, 2, or 3'
      END IF

      b = a

      nx = size(a, 1)
      ny = size(a, 2)
      nz = size(a, 3)

      DO iz = 1, nz
         DO iy = 1, ny
            DO ix = 1, nx
               IF (dim == 1) THEN
                  a(ix, iy, iz) = b(nx-ix+1, iy, iz)
               ELSE IF (dim == 2) THEN
                  a(ix, iy, iz) = b(ix, ny-iy+1, iz)
               ELSE IF (dim == 3) THEN
                  a(ix, iy, iz) = b(ix, iy, nz-iz+1)
               ELSE
                  STOP 'REVERSE_ARRAY_3D: Error, something went wrong'
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE reverse_array_3D

   SUBROUTINE remove_array_element_real(i, a)

      ! Remove element i from array a(n) returning an array of size n-1
      INTEGER, INTENT(IN) :: i                 ! Element to remove
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Input array     
      REAL :: b(size(a)-1)
      INTEGER :: j, jj, n

      n = size(a)

      IF (i < 1 .OR. i > n) THEN
         WRITE (*, *) 'Array size:', n
         WRITE (*, *) 'Element to remove:', i
         STOP 'REMOVE_ARRAY_ELEMENT: Error, element to remove is outside array bounds'
      END IF

      jj = 0
      DO j = 1, n
         IF (j == i) THEN
            CYCLE
         ELSE
            jj = jj+1
            b(jj) = a(j)
         END IF
      END DO

      DEALLOCATE (a)
      ALLOCATE (a(n-1))
      a = b

   END SUBROUTINE remove_array_element_real

   SUBROUTINE remove_array_element_int(i, a)

      ! Remove element i from array a(n) returning an array of size n-1
      INTEGER, INTENT(IN) :: i                    ! Element to remove
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Input array     
      INTEGER :: b(size(a)-1)
      INTEGER :: j, jj, n

      n = size(a)

      IF (i < 1 .OR. i > n) THEN
         WRITE (*, *) 'Array size:', n
         WRITE (*, *) 'Element to remove:', i
         STOP 'REMOVE_ARRAY_ELEMENT: Error, element to remove is outside array bounds'
      END IF

      jj = 0
      DO j = 1, n
         IF (j == i) THEN
            CYCLE
         ELSE
            jj = jj+1
            b(jj) = a(j)
         END IF
      END DO

      DEALLOCATE (a)
      ALLOCATE (a(n-1))
      a = b

   END SUBROUTINE remove_array_element_int

   SUBROUTINE remove_repeated_array_elements(a, m)

      ! Remove any repeated entries from the array
      ! Assumes the array is sorted
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Array to consider
      INTEGER, INTENT(OUT) :: m                ! Final size of array
      INTEGER :: i, n
      LOGICAL :: remove(size(a))

      n = size(a)

      remove = .FALSE.

      m = n
      DO i = 1, n-1
         IF (a(i) == a(i+1)) THEN
            remove(i+1) = .TRUE.
            m = m-1
         END IF
      END DO

      a = PACK(a, .NOT. remove)

   END SUBROUTINE remove_repeated_array_elements

   SUBROUTINE remove_repeated_two_array_elements(a, b, m)

      ! Remove any repeated entries in the array a from both a and b
      ! Assumes the array a is sorted
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: b(:)
      INTEGER, INTENT(OUT) :: m
      INTEGER :: i, n
      LOGICAL :: remove(size(a))

      n = size(a)
      IF (n /= size(b)) STOP 'REMOVE_REPEATED_TWO_ARRAY_ELEMENTS: Error, arrays a and b are differnt sizes'

      remove = .FALSE.

      m = n
      DO i = 1, n-1
         IF (a(i) == a(i+1)) THEN
            remove(i+1) = .TRUE.
            m = m-1
         END IF
      END DO

      a = PACK(a,.NOT. remove)
      b = PACK(b,.NOT. remove)

   END SUBROUTINE remove_repeated_two_array_elements

   SUBROUTINE write_array_list_real(a)

      ! Write out a list of array elements in a neat way
      REAL, INTENT(IN) :: a(:)
      INTEGER :: i, n

      n = size(a)

      WRITE (*, *) 'WRITE_ARRAY_LIST: Writing array'
      DO i = 1, n
         WRITE (*, *) 'WRITE_ARRAY_LIST:', i, a(i)
      END DO
      WRITE (*, *) 'WRITE_ARRAY_LIST: Done'
      WRITE (*, *)

   END SUBROUTINE write_array_list_real

   SUBROUTINE write_array2D_list_real(a)

      ! Write out a list of array elements in a neat way
      REAL, INTENT(IN) :: a(:, :)
      INTEGER :: ix, iy, nx, ny

      nx = size(a, 1)
      ny = size(a, 2)

      WRITE (*, *) 'WRITE_ARRAY_LIST: Writing array'
      DO iy = 1, ny
         DO ix = 1, nx
            WRITE (*, *) 'WRITE_ARRAY_LIST:', ix, iy, a(ix, iy)
         END DO
      END DO
      WRITE (*, *) 'WRITE_ARRAY_LIST: Done'
      WRITE (*, *)

   END SUBROUTINE write_array2D_list_real

   SUBROUTINE write_2array_list_real(a, b)

      ! Write out a list of array elements in a neat way
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: b(:)
      INTEGER :: i, n

      n = size(a)
      IF (n /= size(b)) STOP 'Arrays must be the same size'

      WRITE (*, *) 'WRITE_ARRAY_LIST: Writing array'
      DO i = 1, n
         WRITE (*, *) 'WRITE_ARRAY_LIST:', i, a(i), b(i)
      END DO
      WRITE (*, *) 'WRITE_ARRAY_LIST: Done'
      WRITE (*, *)

   END SUBROUTINE write_2array_list_real

   SUBROUTINE write_array_list_int(a)

      ! Write out a list of array elements in a neat way
      INTEGER, INTENT(IN) :: a(:)
      INTEGER :: i, n

      n = size(a)

      WRITE (*, *) 'WRITE_ARRAY_LIST: Writing array'
      DO i = 1, n
         WRITE (*, *) 'WRITE_ARRAY_LIST:', i, a(i)
      END DO
      WRITE (*, *) 'WRITE_ARRAY_LIST: Done'
      WRITE (*, *)

   END SUBROUTINE write_array_list_int

   FUNCTION splay_2D(a, n1, n2)

      ! This splays out a 2D array 'a' into a 1d array 'b' of the same size (n1*n2)
      ! TODO: Ugly order of arguments (convert to subroutine?)
      INTEGER, INTENT(IN) :: n1, n2
      REAL :: splay_2D(n1*n2)
      REAL, INTENT(IN) :: a(n1, n2)
      INTEGER :: i, j, ii

      ! Set sum integer to zero
      ii = 0

      DO j = 1, n2
         DO i = 1, n1
            ii = ii+1
            splay_2D(ii) = a(i, j)
         END DO
      END DO

   END FUNCTION splay_2D

   FUNCTION splay_3D(a, n1, n2, n3)

      ! This splays out a 3D array 'a' into a 1d array 'b' of the same size (n1*n2*n3)
      ! TODO: Should i, j, k order of loops be reversed?
      ! TODO: Ugly order of arguments (convert to subroutine?)
      INTEGER, INTENT(IN) :: n1, n2, n3
      REAL :: splay_3D(n1*n2*n3)
      REAL, INTENT(IN) :: a(n1, n2, n3)
      INTEGER :: i, j, k, ii

      ! Set sum integer to zero
      ii = 0

      DO i = 1, n1
         DO j = 1, n2
            DO k = 1, n3
               ii = ii+1
               splay_3D(ii) = a(i, j, k)
            END DO
         END DO
      END DO

   END FUNCTION splay_3D

   SUBROUTINE binning(a, a1, a2, b, c, ilog)

      ! Histogram? bin numbers in array a in bins b with height c between a1 and a2
      ! TODO: Is this a repeat of a different subroutine called histogram in statistics.f90?
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a2
      REAL, INTENT(OUT) :: b(:)
      REAL, INTENT(OUT) :: c(:)
      INTEGER, INTENT(IN) :: ilog
      REAL :: min, max
      REAL, ALLOCATABLE :: binlim(:)
      INTEGER :: i, j, n

      n = size(b)

      min = a1
      max = a2

      WRITE (*, *) 'Binning'
      WRITE (*, *) 'Min:', min
      WRITE (*, *) 'Max:', max

      IF (ilog == 1) THEN
         min = log10(min)
         max = log10(max)
      END IF

      ! This sets the limits for the bins!
      CALL fill_array(min, max, binlim, n+1)

      ! This sets the centre value for each bin!
      DO i = 1, n
         b(i) = (binlim(i)+binlim(i+1))/2.
      END DO

      IF (ilog == 1) THEN
         binlim = 10.**binlim
         b = 10.**b
      END IF

      c = 0.

      DO i = 1, size(a)
         DO j = 1, n
            IF (a(i) >= binlim(j) .AND. a(i) <= binlim(j+1)) THEN
               c(j) = c(j)+1.
            END IF
         END DO
      END DO

      WRITE (*, *) 'Binning complete'
      WRITE (*, *)

   END SUBROUTINE binning

   SUBROUTINE merge_arrays(a, b, c, nc)

      ! Takes arrays a and b and merges them together to make c with length SIZE(a)+SIZE(b)
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: b(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: c(:)
      INTEGER, INTENT(OUT) :: nc
      INTEGER :: i, na, nb

      na = size(a)
      nb = size(b)
      nc = na+nb

      IF (allocated(c)) DEALLOCATE (c)
      ALLOCATE (c(nc))

      DO i = 1, na
         c(i) = a(i)
      END DO

      DO i = 1, nb
         c(i+na) = b(i)
      END DO

   END SUBROUTINE merge_arrays

   FUNCTION concatenate_arrays_real(a, b)

      ! Concatenate arrays a and b to form new array with length SIZE(a)+SIZE(b)
      ! TODO: Ugly definition for type of function (convert to subroutine?)
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: b(:)
      REAL :: concatenate_arrays_real(size(a)+size(b))
      INTEGER :: i, na, nb

      na = size(a)
      nb = size(b)

      DO i = 1, na
         concatenate_arrays_real(i) = a(i)
      END DO

      DO i = 1, nb
         concatenate_arrays_real(i+na) = b(i)
      END DO

   END FUNCTION concatenate_arrays_real

   FUNCTION concatenate_arrays_int(a, b)

      ! Concatenate arrays a and b to form new array with length SIZE(a)+SIZE(b)
      ! TODO: Ugly definition for type of function (convert to subroutine?)
      INTEGER, INTENT(IN) :: a(:)
      INTEGER, INTENT(IN) :: b(:)
      INTEGER :: concatenate_arrays_int(size(a)+size(b))
      INTEGER :: i, na, nb

      na = size(a)
      nb = size(b)

      DO i = 1, na
         concatenate_arrays_int(i) = a(i)
      END DO

      DO i = 1, nb
         concatenate_arrays_int(i+na) = b(i)
      END DO

   END FUNCTION concatenate_arrays_int

   SUBROUTINE fill_array(min, max, arr, n, ilog)

      ! Fills array 'arr' in linearly spaced intervals
      ! e.g., 4 values between 0 and 1 would be 0, 1/3, 2/3, and 1
      ! This means that min and max are included in the array
      REAL, INTENT(IN) :: min ! Minimum value for array
      REAL, INTENT(IN) :: max ! Maximum value for array
      REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:) ! Output array
      INTEGER, INTENT(IN) :: n ! Number of entries in array
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      INTEGER :: i

      ! Allocate the array, and deallocate it if it is full
      CALL safe_allocate(arr, n)
      arr = 0.

      IF (n == 1) THEN
         IF (min /= max) STOP 'FILL_ARRAY: If array is to be filled with a single value then min and max must be equal'
         arr(1) = min
      ELSE IF (n > 1) THEN
         DO i = 1, n
            arr(i) = progression(min, max, i, n, ilog)
         END DO
      ELSE
         STOP 'FILL_ARRAY: Error, n cannot be negative'
      END IF

   END SUBROUTINE fill_array

   SUBROUTINE fill_array_log(xmin, xmax, x, nx)

      ! Fill an array x in n log-space intervals between xmin and xmax (inclusive)
      ! TODO: Remove
      REAL, INTENT(IN) :: xmin
      REAL, INTENT(IN) :: xmax
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      INTEGER, INTENT(IN) :: nx

      CALL fill_array(log(xmin), log(xmax), x, nx)
      x = exp(x)

   END SUBROUTINE fill_array_log

   SUBROUTINE integer_sequence(is, i1, i2, n)

      ! Reutrns an array of all the integers between i1 and i2 (inclusive)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: is(:)
      INTEGER, INTENT(IN) :: i1
      INTEGER, INTENT(IN) :: i2
      INTEGER, INTENT(OUT) :: n
      INTEGER :: i, j

      IF(i2 < i1) STOP 'INTEGER_SEQUENCE: Error, i2 is less than i1'

      n = i2-i1+1

      CALL safe_allocate(is, n)
      DO j = 1, n
         i = i1+j-1
         is(j) = i
      END DO

   END SUBROUTINE integer_sequence

   SUBROUTINE fill_pixels(min, max, arr, n)

      ! Fill an array between min and max as if the values should correspond to central pixel values
      ! e.g., 3 pixel values between 0 and 1 would be 1/6, 1/2 and 5/6
      ! This means that min and max are not included in the array
      REAL, INTENT(IN) :: min                  ! Minimum value bordering array
      REAL, INTENT(IN) :: max                  ! Maximum value bordering array
      REAL, ALLOCATABLE, INTENT(OUT) :: arr(:) ! Output array
      INTEGER, INTENT(IN) :: n                 ! Number of pixels in array
      REAL :: mmin, mmax, pixel_size

      pixel_size = (max-min)/real(n) ! Calculate the pixel size
      mmin = min+pixel_size/2.       ! Offset minimum value by half a pixel into the region
      mmax = max-pixel_size/2.       ! Offset maximum value by half a pixel into the region

      CALL fill_array(mmin, mmax, arr, n) ! Linearly fill array

   END SUBROUTINE fill_pixels

   SUBROUTINE create_mask(min, max, m, mask)

      ! Flags objects that make the cut as true in mask
      REAL, INTENT(IN) :: min
      REAL, INTENT(IN) :: max
      REAL, INTENT(IN) :: m(:)
      LOGICAL, INTENT(OUT) :: mask(:)
      INTEGER :: i, o, n

      n = size(m)
      IF (n /= size(mask)) STOP 'CREATE_MASK: Error, real and logical array should be the same size'

      WRITE (*, *) 'MASK: Imposing property cut'
      WRITE (*, *) 'MASK: Minimum value:', min
      WRITE (*, *) 'MASK: Maximum value:', max
      WRITE (*, *) 'MASK: Original number of objects', n

      mask = .FALSE.
      DO i = 1, n
         IF (m(i) >= min .AND. m(i) <= max) mask(i) = .TRUE.
      END DO

      o = count(mask)

      WRITE (*, *) 'MASK: Final number of objects:', o
      WRITE (*, *) 'MASK: Fraction remaining:', real(o)/real(n)
      WRITE (*, *) 'MASK: Fraction culled:', 1.-real(o)/real(n)
      WRITE (*, *) 'MASK: Done'
      WRITE (*, *)

   END SUBROUTINE create_mask

   SUBROUTINE apply_mask_1D(x, a, nx, xmask)

      ! Use the xmask to resize arrays x and a(x) according to only true elements of the mask
      ! TODO: Ugly argument list
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:) ! Input x array
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Input a(x) array
      INTEGER, INTENT(INOUT) :: nx
      LOGICAL, INTENT(IN) :: xmask(nx)
      INTEGER :: i, ii, nx_new
      REAL :: a_old(nx), x_old(nx) ! Not necessary, could use x, but makes easier to read
      INTEGER :: nx_old            ! Not necessary, could use n, but makes easier to read
      REAL, ALLOCATABLE :: a_new(:), x_new(:)

      IF(size(a) .NE. nx) STOP 'APPLY_MASK_1D: Error, you have specified array dimension incorrectly'

      ! Save original arrays (unnecessary)
      x_old = x
      a_old = a
      nx_old = nx

      ! Make new arrays
      nx_new = count(xmask)
      ALLOCATE(x_new(nx_new), a_new(nx_new))
      ii = 0
      DO i = 1, nx
         IF(xmask(i)) THEN
            ii = ii+1
            x_new(ii) = x_old(i)
            a_new(ii) = a_old(i)
         END IF
      END DO

      DEALLOCATE(x, a)
      ALLOCATE(x(nx_new), a(nx_new))
      x = x_new
      a = a_new
      nx = nx_new

      DEALLOCATE(x_new, a_new)

   END SUBROUTINE apply_mask_1D

   SUBROUTINE apply_mask_2D(x, y, a, nx, ny, xmask, ymask)

      ! Use xmask and ymask to resize arrays x, y and a(x, y) according to only true elements of the masks
      ! TODO: Ugly argument list
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: y(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(INOUT) :: nx
      INTEGER, INTENT(INOUT) :: ny
      LOGICAL, INTENT(IN) :: xmask(nx)
      LOGICAL, INTENT(IN) :: ymask(ny)
      INTEGER :: i, ii, j, jj, nx_new, ny_new
      REAL :: x_old(nx), y_old(ny), a_old(nx, ny) ! Not necessary, could use x, but makes easier to read
      INTEGER :: nx_old, ny_old ! Not necessary, could use n, but makes easier to read
      REAL, ALLOCATABLE :: x_new(:), y_new(:), a_new(:, :)

      IF(size(a(:,1)) .NE. nx) THEN
         WRITE(*, *) 'APPLY_MASK_2D: nx:', nx
         WRITE(*, *) 'APPLY_MASK_2D: size of array first dimension:', size(a(:,1))
         STOP 'APPLY_MASK_2D: Error, you have specified x array dimension incorrectly'
      END IF
      IF(size(a(1,:)) .NE. ny) THEN
         WRITE(*, *) 'APPLY_MASK_2D: ny:', ny
         WRITE(*, *) 'APPLY_MASK_2D: size of array second dimension:', size(a(1,:))
         STOP 'APPLY_MASK_2D: Error, you have specified y array dimension incorrectly'
      END IF

      ! Save original array
      x_old = x
      y_old = y
      a_old = a
      nx_old = nx
      ny_old = ny

      ! Make new array
      nx_new = count(xmask)
      ny_new = count(ymask)
      ALLOCATE(x_new(nx_new), y_new(ny_new), a_new(nx_new, ny_new))
      jj=0
      DO j = 1, ny
         IF(ymask(j)) THEN
            jj=jj+1
            y_new(jj) = y_old(j)
            ii = 0
            DO i = 1, nx
               IF(xmask(i)) THEN
                  ii = ii+1
                  x_new(ii) = x_old(i)
                  a_new(ii, jj) = a_old(i, j)
               END IF
            END DO
         END IF
      END DO

      DEALLOCATE(x, y, a)
      ALLOCATE(x(nx_new), y(ny_new), a(nx_new, ny_new))
      x = x_new
      y = y_new
      a = a_new
      nx = nx_new
      ny = ny_new

      DEALLOCATE(x_new, y_new, a_new)

   END SUBROUTINE apply_mask_2D

   INTEGER FUNCTION unique_entries(a)

      ! Counts the total number of unique entries in array 'a'
      INTEGER, INTENT(IN) :: a(:) ! Input array
      INTEGER :: i, j, n

      n = size(a)

      ! Initially assume all entries are unique
      unique_entries = n

      ! Double loop to each each pair against each other
      DO i = 1, n
         DO j = i+1, n
            IF (a(j) == a(i)) THEN
               unique_entries = unique_entries-1 ! A non-unique entry has been discovered, so subtract one
               EXIT
            END IF
         END DO
      END DO

   END FUNCTION unique_entries

   SUBROUTINE unique_index(array, n, unique, m, match)

      ! Create an indexing scheme to only the unique elements
      INTEGER, INTENT(IN) :: n                       ! Number of entries in input array
      INTEGER, INTENT(IN) :: array(n)                ! Array to find the unique indices of
      INTEGER, ALLOCATABLE, INTENT(OUT) :: unique(:) ! Output array of unique indices
      INTEGER, INTENT(OUT) :: m                      ! Number of unique indices
      INTEGER, INTENT(OUT) :: match(n)               ! Array for matching input and unique arrays
      INTEGER :: i, j, p
      LOGICAL :: inc

      ! First count the number of unique entries
      m = unique_entries(array)
      ALLOCATE (unique(m))

      ! Set the first unique entry to be the first entry
      unique(1) = array(1)

      p = 1
      DO i = 1, n
         inc = .FALSE.
         DO j = 1, p
            IF (array(i) /= unique(j)) THEN
               unique(p+1) = array(i)
               inc = .TRUE.
               EXIT
            END IF
         END DO
         IF (inc) p = p+1
         IF (p == m) EXIT ! Added this in haste
      END DO

      ! Now fill the matching array
      DO j = 1, m
         DO i = 1, n
            IF (unique(j) == array(i)) THEN
               match(i) = j
            END IF
         END DO
      END DO

   END SUBROUTINE unique_index

   LOGICAL FUNCTION is_in_array_integer(x, a)

      ! Test to see if integer value x is in a(n)
      ! This is very useful in IF statements where otherwise there would be a long chain of .OR.s
      INTEGER, INTENT(IN) :: x
      INTEGER, INTENT(IN) :: a(:)
      INTEGER :: i

      is_in_array_integer = .FALSE.
      DO i = 1, size(a)
         IF (x == a(i)) THEN
            is_in_array_integer = .TRUE.
            EXIT
         END IF     
      END DO

   END FUNCTION is_in_array_integer

   LOGICAL FUNCTION is_in_array_character(x, a)

      ! Test to see if integer value x is in a(n)
      ! This is very useful in IF statements where otherwise there would be a long chain of .OR.'s
      CHARACTER(len=*), INTENT(IN) :: x
      CHARACTER(len=*), INTENT(IN) :: a(:)
      INTEGER :: i

      is_in_array_character = .FALSE.
      DO i = 1, size(a)
         IF (trim(x) == trim(a(i))) THEN
            is_in_array_character = .TRUE.
            EXIT
         END IF     
      END DO

   END FUNCTION is_in_array_character

   LOGICAL FUNCTION repeated_entries(a)

      ! Checks for repeated entries in a(n) either true or false
      INTEGER, INTENT(IN) :: a(:)
      INTEGER :: i, j, n

      n = size(a)

      repeated_entries = .FALSE.

      DO i = 1, n
         DO j = i+1, n
            IF (a(i) == a(j)) THEN
               repeated_entries = .TRUE.
               EXIT
            END IF
         END DO
         IF (repeated_entries) EXIT
      END DO

   END FUNCTION repeated_entries

   LOGICAL FUNCTION regular_spacing(a)

      ! Returns true if array a is regularly spaced
      REAL, INTENT(IN) :: a(:)
      REAL :: amin, amax, b
      INTEGER :: i, n
      REAL, PARAMETER :: eps = 1e-5

      n = size(a)

      amin = a(1)
      amax = a(n)

      regular_spacing = .TRUE.
      DO i = 1, n
         b = progression(amin, amax, i, n)
         IF (.NOT. requal(b, a(i), eps)) THEN
            regular_spacing = .FALSE.
            EXIT
         END IF
      END DO

   END FUNCTION regular_spacing

   INTEGER FUNCTION sum_logical(a)

      ! Count the number of true entries in a logical array
      LOGICAL, INTENT(IN) :: a(:)
      INTEGER :: i

      sum_logical = 0
      DO i = 1, size(a)
         IF (a(i)) sum_logical = sum_logical+1
      END DO

   END FUNCTION sum_logical

END MODULE array_operations
