MODULE array_operations

   ! These interfaces did not work if -fdefault-real-8 and -fdefault-real-8 were set
!!$  INTERFACE fill_array
!!$     MODULE PROCEDURE fill_array_single
!!$     MODULE PROCEDURE fill_array_double
!!$  END INTERFACE fill_array
!!$
!!$  INTERFACE progression
!!$     MODULE PROCEDURE progression_single
!!$     MODULE PROCEDURE progression_double
!!$  END INTERFACE progression
!!$
!!$  INTERFACE progression_log
!!$     MODULE PROCEDURE progression_log_single
!!$     MODULE PROCEDURE progression_log_double
!!$  END INTERFACE progression_log

   PRIVATE

   PUBLIC :: within_array
   PUBLIC :: swap_arrays
   PUBLIC :: append
   PUBLIC :: add_to_array
   PUBLIC :: splay
   PUBLIC :: write_array_list
   PUBLIC :: maximum
   PUBLIC :: sum_double
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
   PUBLIC :: mask
   PUBLIC :: apply_mask
   PUBLIC :: smooth_array
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
   !PUBLIC :: fill_array_double ! TODO: Delete
   PUBLIC :: fill_array_log

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

   INTERFACE splay
      PROCEDURE splay_2D
      PROCEDURE splay_3D
   END INTERFACE splay

   INTERFACE write_array_list
      PROCEDURE write_array_list_real
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
      PROCEDURE safe_allocate_integer
      PROCEDURE safe_allocate_logical
      PROCEDURE safe_allocate_character
   END INTERFACE safe_allocate

CONTAINS

   LOGICAL FUNCTION greater_than_any(x, a)

      ! Returns TRUE if x is greater than any one of the values in array
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: a(:)
      INTEGER :: i, n

      n = size(array)

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
      IMPLICIT NONE
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
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(x)
      ALLOCATE(x(n))

   END SUBROUTINE safe_allocate_real

   ! SUBROUTINE safe_allocate_double(x, n)

   !    IMPLICIT NONE
   !    DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: x(:)
   !    INTEGER, INTENT(IN) :: n

   !    CALL if_allocated_deallocate(x)
   !    ALLOCATE(x(n))

   ! END SUBROUTINE safe_allocate_double

   SUBROUTINE safe_allocate_integer(i, n)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      IMPLICIT NONE
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: i(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(i)
      ALLOCATE(i(n))

   END SUBROUTINE safe_allocate_integer

   SUBROUTINE safe_allocate_logical(l, n)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      IMPLICIT NONE
      LOGICAL, ALLOCATABLE, INTENT(INOUT) :: l(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(l)
      ALLOCATE(l(n))

   END SUBROUTINE safe_allocate_logical

   SUBROUTINE safe_allocate_character(c, n)

      ! Checks array for allocation status, deallocates if necessary, then allocates
      IMPLICIT NONE
      CHARACTER(len=*), ALLOCATABLE, INTENT(INOUT) :: c(:)
      INTEGER, INTENT(IN) :: n

      CALL if_allocated_deallocate(c)
      ALLOCATE(c(n))

   END SUBROUTINE safe_allocate_character

   SUBROUTINE if_allocated_deallocate_real_1D(x)

      ! Deallocates an array if it is already allocated
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)

      IF(ALLOCATED(x)) DEALLOCATE(x)

   END SUBROUTINE if_allocated_deallocate_real_1D

   SUBROUTINE if_allocated_deallocate_real_2D(x)

      ! Deallocates an array if it is already allocated
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:,:)

      IF(ALLOCATED(x)) DEALLOCATE(x)

   END SUBROUTINE if_allocated_deallocate_real_2D

   SUBROUTINE if_allocated_deallocate_real_3D(x)

      ! Deallocates an array if it is already allocated
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:,:,:)

      IF(ALLOCATED(x)) DEALLOCATE(x)

   END SUBROUTINE if_allocated_deallocate_real_3D

   SUBROUTINE if_allocated_deallocate_integer_1D(i)

      ! Deallocates an array if it is already allocated
      IMPLICIT NONE
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: i(:)

      IF(ALLOCATED(i)) DEALLOCATE(i)

   END SUBROUTINE if_allocated_deallocate_integer_1D

   SUBROUTINE if_allocated_deallocate_logical_1D(l)

      ! Deallocates an array if it is already allocated
      IMPLICIT NONE
      LOGICAL, ALLOCATABLE, INTENT(INOUT) :: l(:)

      IF(ALLOCATED(l)) DEALLOCATE(l)

   END SUBROUTINE if_allocated_deallocate_logical_1D

   SUBROUTINE if_allocated_deallocate_character_1D(c)

      ! Deallocates an array if it is already allocated
      IMPLICIT NONE
      CHARACTER(len=*), ALLOCATABLE, INTENT(INOUT) :: c(:)

      IF(ALLOCATED(c)) DEALLOCATE(c)

   END SUBROUTINE if_allocated_deallocate_character_1D

   SUBROUTINE smooth_array(x, m)

      ! Take an array and smooth it by taking the over a fixed number of values in each direction
      ! For example, if m = 2 then poisition i in the smoothed array is the average over i-2, i-1, i, i+1, i+2 in the old array
      ! It probably only makes sense to do this if x(t) is a function and the t are uniformally distributed
      IMPLICIT NONE
      REAL, INTENT(INOUT) :: x(:) ! Array to smooth
      INTEGER, INTENT(IN) :: m    ! Number of entries to smooth over
      INTEGER :: i, iup, idn, j, mm, n
      REAL, ALLOCATABLE :: x_save(:)

      n = size(x)
      ALLOCATE(x_save(n))

      ! Save the original input array
      x_save = x

      ! Delete the original array
      x = 0.

      ! Make a running average over the function
      DO i = 1, n
         mm = 1
         x(i) = x_save(i)
         DO j = 1, m
            iup = i+j
            IF(iup <= n) THEN
               x(i) = x(i) + x_save(iup)
               mm = mm + 1
            END IF
            idn = i-j
            IF(idn >= 1) THEN
               x(i) = x(i) + x_save(idn)
               mm = mm + 1
            END IF
         END DO
         x(i) = x(i)/real(mm)
      END DO

   END SUBROUTINE smooth_array

   LOGICAL FUNCTION within_array(x, a)

      ! Checks to see if x falls within the range of values in array a
      IMPLICIT NONE
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
      IMPLICIT NONE
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

   SUBROUTINE append(a, b)

      ! Append value b to the end of array a(n) to make a new array a(n+1)
      IMPLICIT NONE
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: b
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
      a(n+1) = b

   END SUBROUTINE append

   SUBROUTINE add_to_array_2D(a, m, v, i)

      ! Add value 'v' to array 'a' at location 'i' in array
      ! If 'i' is outside the array range then this routine does nothing
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(INOUT) :: a(m, m)
      REAL, INTENT(IN) :: v
      INTEGER, INTENT(IN) :: i(2)
      INTEGER :: j
      LOGICAL :: bin
      INTEGER, PARAMETER :: dim = 2

      bin = .TRUE.
      DO j = 1, dim
         IF (i(j) < 1 .OR. i(j) > m) THEN
            bin = .FALSE.
            EXIT
         END IF
      END DO

      IF (bin) a(i(1), i(2)) = a(i(1), i(2))+v

   END SUBROUTINE add_to_array_2D

   SUBROUTINE add_to_array_3D(a, m, v, i)

      ! Add value 'v' to array 'a' at location 'i' in array
      ! If 'i' is outside the array range then this routine does nothing
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(INOUT) :: a(m, m, m)
      REAL, INTENT(IN) :: v
      INTEGER, INTENT(IN) :: i(3)
      INTEGER :: j
      LOGICAL :: bin
      INTEGER, PARAMETER :: dim = 3

      bin = .TRUE.
      DO j = 1, dim
         IF (i(j) < 1 .OR. i(j) > m) THEN
            bin = .FALSE.
            EXIT
         END IF
      END DO

      IF (bin) a(i(1), i(2), i(3)) = a(i(1), i(2), i(3))+v

   END SUBROUTINE add_to_array_3D

   INTEGER FUNCTION array_position_int(x, a)

      ! Returns the location in a(n) of value x
      ! If x is not in array then returns zero
      IMPLICIT NONE
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
      USE basic_operations
      IMPLICIT NONE
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
      IMPLICIT NONE
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
      IMPLICIT NONE
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

   REAL FUNCTION sum_double(a)

      ! Sum using double precision, which is necessary for many array elements
      IMPLICIT NONE
      REAL, INTENT(IN) :: a(:)
      DOUBLE PRECISION :: sum
      INTEGER :: i, n

      n = size(a)

      sum = 0.d0

      DO i = 1, n
         sum = sum+a(i)
      END DO

      sum_double = real(sum)

   END FUNCTION sum_double

   SUBROUTINE amputate_array(a, i1, i2)

      ! Chop an array of size a(n) down to a smaller size demarked by indices i1, i2
      ! If i1=1 and i2=n then this does nothing
      IMPLICIT NONE
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

   END SUBROUTINE amputate_array

   SUBROUTINE reduce_array(a1, a2, n2)

      ! Reduces the size of array1 to the size of array2
      ! This will not preserve the spacing of entries in array1 and might be a terrible idea in many cases
      ! TODO: Remove this
      IMPLICIT NONE   
      REAL, INTENT(IN) :: a1(:)
      REAL, INTENT(OUT) :: a2(n2)
      INTEGER, INTENT(IN) :: n2
      INTEGER :: i, j, n1

      n1 = size(a1)

      DO i = 1, n2
         j = 1+ceiling(real((n1-1)*(i-1))/real(n2-1))
         a2(i) = a1(j)
      END DO

   END SUBROUTINE reduce_array

!!$  SUBROUTINE reduceto(arr1,n)
!!$
!!$    ! Reduces the array from whatever size to size 'n'
!!$    ! TODO: Remove
!!$    IMPLICIT NONE
!!$    REAL, ALLOCATABLE, INTENT(INOUT) :: arr1(:)
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, ALLOCATABLE :: hold(:)
!!$    INTEGER :: i, j
!!$
!!$    ALLOCATE(hold(n))
!!$
!!$    DO i=1,n
!!$       j=1+ceiling(real((n-1)*(i-1))/real(n-1))
!!$       hold(i)=arr1(j)
!!$    END DO
!!$
!!$    DEALLOCATE(arr1)
!!$    ALLOCATE(arr1(n))
!!$
!!$    arr1=hold
!!$
!!$    DEALLOCATE(hold)
!!$
!!$  END SUBROUTINE reduceto

   SUBROUTINE reverse_array(a)

      ! Reverses the contents of array
      IMPLICIT NONE
      REAL, INTENT(INOUT) :: a(:)
      INTEGER :: i, n
      REAL :: hold(size(a))

      n = size(a)

      hold = a

      DO i = 1, n
         a(i) = hold(n-i+1)
      END DO

   END SUBROUTINE reverse_array

   SUBROUTINE remove_array_element(a, i)

      ! Remove element i from array a(n) returning an array of size n-1
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Input array
      INTEGER, INTENT(IN) :: i                 ! Element to remove
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

   END SUBROUTINE remove_array_element

   SUBROUTINE remove_repeated_array_elements(a, m)

      ! Remove any repeated entries from the array
      ! Assumes the array is sorted
      IMPLICIT NONE
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
      IMPLICIT NONE
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
      IMPLICIT NONE
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

   SUBROUTINE write_array_list_int(a)

      ! Write out a list of array elements in a neat way
      IMPLICIT NONE
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
      IMPLICIT NONE
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
      IMPLICIT NONE
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

      IMPLICIT NONE
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(OUT) :: b(:)
      REAL, INTENT(OUT) :: c(:)
      INTEGER, INTENT(IN) :: ilog
      REAL :: a1, a2, min, max
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
      IMPLICIT NONE
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: b(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: c(:)
      INTEGER, INTENT(OUT) :: nc
      INTEGER :: i, na, nb

      na = size(a)
      nb = size(b)
      nc = na+nb

      IF (ALLOCATED(c)) DEALLOCATE (c)
      ALLOCATE (c(nc))

      DO i = 1, na
         c(i) = a(i)
      END DO

      DO i = 1, nb
         c(i+na) = b(i)
      END DO

   END SUBROUTINE merge_arrays

   FUNCTION concatenate_arrays(a, b)

      ! Concatenate arrays a and b to form new array with length SIZE(a)+SIZE(b)
      ! TODO: Ugly definition for type of function
      IMPLICIT NONE   
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: b(:)
      REAL :: concatenate_arrays(size(a)+size(b))
      INTEGER :: i, na, nb

      na = size(a)
      nb = size(b)

      DO i = 1, na
         concatenate_arrays(i) = a(i)
      END DO

      DO i = 1, nb
         concatenate_arrays(i+na) = b(i)
      END DO

   END FUNCTION concatenate_arrays

   SUBROUTINE fill_array(min, max, arr, n)

      ! Fills array 'arr' in linearly spaced intervals
      ! e.g., 4 values between 0 and 1 would be 0, 1/3, 2/3, and 1
      ! This means that min and max are included in the array
      USE basic_operations
      IMPLICIT NONE
      REAL, INTENT(IN) :: min ! Minimum value for array
      REAL, INTENT(IN) :: max ! Maximum value for array
      REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:) ! Output array
      INTEGER, INTENT(IN) :: n ! Number of entries in array
      INTEGER :: i

      ! Allocate the array, and deallocate it if it is full
      !IF (ALLOCATED(arr)) DEALLOCATE (arr)
      !ALLOCATE (arr(n))'
      CALL safe_allocate(arr, n)
      arr = 0.

      IF (n == 1) THEN
         arr(1) = min
      ELSE IF (n > 1) THEN
         DO i = 1, n
            arr(i) = progression(min, max, i, n)
         END DO
      END IF

   END SUBROUTINE fill_array

   ! SUBROUTINE fill_array_double(min, max, arr, n)

   !    ! Fills array 'arr' in equally spaced intervals
   !    ! TODO: Not sure if inputting an array like this is okay
   !    ! TODO: Delete
   !    USE basic_operations
   !    IMPLICIT NONE
   !    DOUBLE PRECISION, INTENT(IN) :: min, max
   !    DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: arr(:)
   !    INTEGER, INTENT(IN) :: n
   !    INTEGER :: i

   !    ! Allocate the array, and deallocate it if it is full
   !    !IF (ALLOCATED(arr)) DEALLOCATE (arr)
   !    !ALLOCATE (arr(n))
   !    CALL safe_allocate(arr, n)
   !    arr = 0.d0

   !    IF (n == 1) THEN
   !       arr(1) = min
   !    ELSE IF (n > 1) THEN
   !       DO i = 1, n
   !          arr(i) = progression_double(min, max, i, n)
   !       END DO
   !    END IF

   ! END SUBROUTINE fill_array_double

   SUBROUTINE fill_array_log(xmin, xmax, x, nx)

      ! Fill an array x in n log-space intervals between xmin and xmax (inclusive)
      IMPLICIT NONE
      REAL, INTENT(IN) :: xmin
      REAL, INTENT(IN) :: xmax
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      INTEGER, INTENT(IN) :: nx

      CALL fill_array(log(xmin), log(xmax), x, nx)
      x = exp(x)

   END SUBROUTINE fill_array_log

   SUBROUTINE integer_sequence(is, i1, i2, n)

      ! Reutrns an array of all the integers between i1 and i2 (inclusive)
      IMPLICIT NONE
      INTEGER, ALLOCATABLE, INTENT(OUT) :: is(:)
      INTEGER, INTENT(IN) :: i1
      INTEGER, INTENT(IN) :: i2
      INTEGER, INTENT(OUT) :: n
      INTEGER :: i, j

      IF(i2 < i1) STOP 'INTEGER_SEQUENCE: Error, i2 is less than i1'

      n = i2-i1+1
      !CALL if_allocated_deallocate(is)
      !ALLOCATE(is(n))
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
      IMPLICIT NONE
      REAL, INTENT(IN) :: min ! Minimum value bordering array
      REAL, INTENT(IN) :: max ! Maximum value bordering array
      REAL, ALLOCATABLE, INTENT(OUT) :: arr(:) ! Output array
      INTEGER, INTENT(IN) :: n ! Number of pixels in array
      REAL :: mmin, mmax, pixel_size

      pixel_size = (max-min)/real(n) ! Calculate the pixel size
      mmin = min+pixel_size/2. ! Offset minimum value by half a pixel into the region
      mmax = max-pixel_size/2. ! Offset maximum value by half a pixel into the region

      CALL fill_array(mmin, mmax, arr, n) ! Linearly fill array

   END SUBROUTINE fill_pixels

   REAL FUNCTION maximum(x, y)

      ! From an array y(x) finds the x location of the maximum treating y(x) as a continuous function
      USE special_functions
      IMPLICIT NONE
      REAL, INTENT(IN) :: x(:)
      REAL, INTENT(IN) :: y(:)
      REAL :: x1, x2, x3, y1, y2, y3, a, b, c
      INTEGER :: imax(1), i, n

      n = size(x)
      IF(n /= size(y)) STOP 'MAXIMUM: Error, x and y must be the same size'

      ! Need this to stop a compile-time warning
      maximum = 0.

      ! Integer maximum location
      imax = maxloc(y)
      i = imax(1)

      IF (i == 1 .OR. i == n) THEN

         STOP 'MAXIMUM: Error, maximum array value is at one end of the array'

      ELSE

         ! Get the x positions
         x1 = x(i-1)
         x2 = x(i)
         x3 = x(i+1)

         ! Get the y values
         y1 = y(i-1)
         y2 = y(i)
         y3 = y(i+1)

         ! Fix a parabola around the maximum
         CALL fix_quadratic(a, b, c, x1, y1, x2, y2, x3, y3)

         ! Read off the maximum from the parabola
         maximum = -b/(2.*a)

      END IF

   END FUNCTION maximum

   SUBROUTINE mask(okay, m, n, min, max)

      ! Flags objects that make the cut as 'okay'
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(OUT) :: okay(n)
      REAL, INTENT(IN) :: m(n)
      REAL, INTENT(IN) :: min
      REAL, INTENT(IN) :: max
      INTEGER :: i, o

      WRITE (*, *) 'MASK: Imposing property cut'
      WRITE (*, *) 'MASK: Minimum value:', min
      WRITE (*, *) 'MASK: Maximum value:', max
      WRITE (*, *) 'MASK: Original number of objects', n

      okay = .FALSE.

      DO i = 1, n
         IF (m(i) >= min .AND. m(i) <= max) okay(i) = .TRUE.
      END DO

      o = count(okay)

      WRITE (*, *) 'MASK: Final number of objects:', o
      WRITE (*, *) 'MASK: Fraction remaining:', real(o)/real(n)
      WRITE (*, *) 'MASK: Fraction culled:', 1.-real(o)/real(n)
      WRITE (*, *) 'MASK: Done'
      WRITE (*, *)

   END SUBROUTINE mask

   SUBROUTINE apply_mask_1D(x, a, nx, xmask)

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(INOUT) :: nx
      LOGICAL, INTENT(IN) :: xmask(nx)
      INTEGER :: i, ii, nx_new
      REAL :: a_old(nx), x_old(nx) ! Not necessary, could use x, but makes easier to read
      INTEGER :: nx_old ! Not necessary, could use n, but makes easier to read
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

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: y(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:,:)
      INTEGER, INTENT(INOUT) :: nx
      INTEGER, INTENT(INOUT) :: ny
      LOGICAL, INTENT(IN) :: xmask(nx)
      LOGICAL, INTENT(IN) :: ymask(ny)
      INTEGER :: i, ii, j, jj, nx_new, ny_new
      REAL :: x_old(nx), y_old(ny), a_old(nx, ny) ! Not necessary, could use x, but makes easier to read
      INTEGER :: nx_old, ny_old ! Not necessary, could use n, but makes easier to read
      REAL, ALLOCATABLE :: x_new(:), y_new(:), a_new(:,:)

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
      IMPLICIT NONE
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

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n                       ! Number of entries in input array
      INTEGER, INTENT(IN) :: array(n)                ! Array to find the unique indices of
      INTEGER, ALLOCATABLE, INTENT(OUT) :: unique(:) ! Output array of unique indices
      INTEGER, INTENT(OUT) :: m                      ! Number of unique indices
      INTEGER, INTENT(OUT) :: match(n)               ! Array for matching input and unique arrays
      INTEGER :: i, j, p
      LOGICAL :: increment

      ! First count the number of unique entries
      m = unique_entries(array)
      ALLOCATE (unique(m))

      ! Set the first unique entry to be the first entry
      unique(1) = array(1)

      p = 1
      DO i = 1, n
         increment = .FALSE.
         DO j = 1, p
            IF (array(i) /= unique(j)) THEN
               unique(p+1) = array(i)
               increment = .TRUE.
               EXIT
            END IF
         END DO
         IF (increment) p = p+1
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
      IMPLICIT NONE
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
      ! This is very useful in IF statements where otherwise there would be a long chain of .OR.s
      IMPLICIT NONE
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

      ! Checks for repeated entries in a(n)
      IMPLICIT NONE
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
      USE basic_operations
      IMPLICIT NONE
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

END MODULE array_operations
