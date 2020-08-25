MODULE sorting

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: index
   PUBLIC :: sort
   PUBLIC :: sorted
   PUBLIC :: sorted_index
   PUBLIC :: reindex
   PUBLIC :: isort_bubble
   PUBLIC :: isort_stupid
   PUBLIC :: isort_QsortC

   INTERFACE index
      MODULE PROCEDURE index_real
      MODULE PROCEDURE index_int
   END INTERFACE index

   INTERFACE bubble_index
      MODULE PROCEDURE bubble_index_real
      MODULE PROCEDURE bubble_index_int
   END INTERFACE bubble_index

   INTERFACE stupid_index
      MODULE PROCEDURE stupid_index_real
      MODULE PROCEDURE stupid_index_int
   END INTERFACE stupid_index

   INTEGER, PARAMETER :: isort_stupid = 1
   INTEGER, PARAMETER :: isort_bubble = 2
   INTEGER, PARAMETER :: isort_QsortC = 3

CONTAINS

   SUBROUTINE sort(a, isort)

      ! Sort arrays in order from lowest to highest values
      !INTEGER, INTENT(IN) :: n
      REAL, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: isort

      IF (isort == isort_stupid) THEN
         CALL stupid_sort(a)
      ELSE IF (isort == isort_bubble) THEN
         CALL bubble_sort(a)
      ELSE IF (isort == isort_QsortC) THEN
         CALL QsortC(a)
      ELSE
         STOP 'SORT: Error, isort not specified correctly'
      END IF

   END SUBROUTINE sort

   SUBROUTINE bubble_sort(a)

      ! Bubble sort array 'a' into lowest to highest value
      REAL, INTENT(INOUT) :: a(:)
      REAL :: hold
      INTEGER :: i, n
      LOGICAL :: is_sorted

      n = size(a)

      DO
         is_sorted = .TRUE.
         DO i = 1, n-1
            IF (a(i) > a(i+1)) THEN
               hold = a(i+1)
               a(i+1) = a(i)
               a(i) = hold
               is_sorted = .FALSE.
            END IF
         END DO
         IF (is_sorted) EXIT
      END DO

   END SUBROUTINE bubble_sort

   SUBROUTINE stupid_sort(a)

      ! I have no idea what this is
      REAL, INTENT(INOUT) :: a(:)
      REAL :: hold, min
      INTEGER :: i, j, minl, n

      n = size(a)

      DO i = 1, n-1
         min = a(i)
         minl = i
         DO j = i+1, n
            IF (a(j) < min) THEN
               min = a(j)
               minl = j
            END IF
         END DO
         hold = a(i)
         a(i) = min
         a(minl) = hold
      END DO

   END SUBROUTINE stupid_sort

   RECURSIVE SUBROUTINE QsortC(A)

      ! Stolen from http://www.fortran.com/qsort_c.f95
      REAL, INTENT(INOUT), DIMENSION(:) :: A
      INTEGER :: iq

      IF (size(A) > 1) THEN
         CALL Partition(A, iq)
         CALL QsortC(A(:iq-1))
         CALL QsortC(A(iq:))
      END IF

   END SUBROUTINE QsortC

   SUBROUTINE Partition(A, marker)

      ! Stolen from http://www.fortran.com/qsort_c.f95
      REAL, INTENT(INOUT), DIMENSION(:) :: A
      INTEGER, INTENT(OUT) :: marker
      INTEGER :: i, j
      REAL :: temp
      REAL :: x      ! pivot point

      x = A(1)
      i = 0
      j = size(A)+1

      do
         j = j-1
         do
            if (A(j) <= x) exit
            j = j-1
         end do
         i = i+1
         do
            if (A(i) >= x) exit
            i = i+1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
         elseif (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         endif
      end do

   END SUBROUTINE Partition

   SUBROUTINE index_real(a, ind, isort)

      ! Index the array 'a' from lowest to highest value
      REAL, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: ind(:)
      INTEGER, INTENT(IN) :: isort

      IF (isort == isort_stupid) THEN
         CALL stupid_index_real(a, ind)
      ELSE IF (isort == isort_bubble) THEN
         CALL bubble_index_real(a, ind)
      ELSE
         STOP 'INDEX_REAL: Error, isort specified incorrectly'
      END IF

   END SUBROUTINE index_real

   SUBROUTINE index_int(a, ind, isort)

      ! Index the array 'a' from lowest to highest value  
      INTEGER, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: ind(:)
      INTEGER, INTENT(IN) :: isort
      
      IF (isort == isort_stupid) THEN
         CALL stupid_index_int(a, ind)
      ELSE IF (isort == isort_bubble) THEN
         CALL bubble_index_int(a, ind)
      ELSE
         STOP 'INDEX_REAL: Error, isort specified incorrectly'
      END IF

   END SUBROUTINE index_int

   SUBROUTINE bubble_index_real(a, ind)

      ! Create an index array for a(:) that indexes from smallest to largest value
      REAL, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: ind(:)
      INTEGER :: i, isort, hold, n

      n = size(a)
      IF (n /= size(ind)) STOP 'BUBBLE_INDEX_REAL: Error, a and ind must have the same size'

      DO i = 1, n
         ind(i) = i
      END DO

      DO
         isort = 0
         DO i = 1, n-1
            IF (a(ind(i)) > a(ind(i+1))) THEN
               hold = ind(i+1)
               ind(i+1) = ind(i)
               ind(i) = hold
               isort = 1
            END IF
         END DO
         IF (isort == 0) EXIT
      END DO

   END SUBROUTINE bubble_index_real

   SUBROUTINE bubble_index_int(a, ind)!,verbose)

      ! Create an index array for integer a(:) that indexes from smallest to largest value
      INTEGER, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: ind(:)
      INTEGER :: i, isort, hold, n

      n = size(a)
      IF (n /= size(ind)) STOP 'BUBBLE_INDEX_INT: Error, a and ind must have the same size'

      DO i = 1, n
         ind(i) = i
      END DO

      DO
         isort = 0
         DO i = 1, n-1
            IF (a(ind(i)) > a(ind(i+1))) THEN
               hold = ind(i+1)
               ind(i+1) = ind(i)
               ind(i) = hold
               isort = 1
            END IF
         END DO
         IF (isort == 0) EXIT
      END DO

   END SUBROUTINE bubble_index_int

   SUBROUTINE stupid_index_real(a, ind)
  
      REAL, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: ind(:)
      INTEGER :: i, j, n
      REAL, ALLOCATABLE :: b(:)

      n = size(a)
      IF (n /= size(ind)) STOP 'STUPID_INDEX_REAL: Error, a and ind must have the same size'
      ALLOCATE(b(n))

      b = a

      ! This is probably stupid
      DO i = 1, n
         j = MINLOC(b, 1)
         ind(i) = j
         b(j) = HUGE(b)
      END DO

   END SUBROUTINE stupid_index_real

   SUBROUTINE stupid_index_int(a, ind)

      INTEGER, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: ind(:)
      INTEGER :: i, j, n
      INTEGER, ALLOCATABLE :: b(:)

      n = size(a)
      IF (n /= size(ind)) STOP 'STUPID_INDEX_INT: Error, a and ind must have the same size'
      ALLOCATE(b(n))

      b = a

      ! This is probably stupid
      DO i = 1, n
         j = MINLOC(b, 1)
         ind(i) = j
         b(j) = HUGE(b)
      END DO

   END SUBROUTINE stupid_index_int

   LOGICAL FUNCTION sorted(a)

      ! Checks if array 'a' is sorted from highest to lowest
      REAL, INTENT(IN) :: a(:) ! Input array to check
      INTEGER :: i, n

      n = size(a)

      sorted = .TRUE.

      DO i = 1, n-1
         IF (a(i) > a(i+1)) THEN
            sorted = .FALSE.
            EXIT
         END IF
      END DO

   END FUNCTION sorted

   LOGICAL FUNCTION sorted_index(a, j)

      ! Checks if array indices for 'a' are sorted from highest to lowest
      REAL, INTENT(IN) :: a(:)    ! Input array to check
      INTEGER, INTENT(IN) :: j(:) ! Input array indices to check
      INTEGER :: i, n

      n = size(a)
      IF (n /= size(j)) STOP 'CHECK_SORTED_INDEX: Error a and j should be same size'

      sorted_index = .TRUE.

      DO i = 1, n-1
         IF (a(j(i)) > a(j(i+1))) THEN
            sorted_index = .FALSE.
            EXIT
         END IF
      END DO

   END FUNCTION sorted_index

   SUBROUTINE reindex(a, j)

      ! Reindex the array 'a' with the new indices 'j'
      REAL, INTENT(INOUT) :: a(:) ! Input array to check
      INTEGER, INTENT(IN) :: j(:) ! Input array indices to check
      REAL, ALLOCATABLE :: b(:)
      INTEGER :: i, n

      n = size(a)
      IF (n /= size(j)) STOP 'REINDEX: Error, a and j must have the same size'
      ALLOCATE(b(n))

      b = a ! Store the input array

      a = 0. ! Delete the input array

      ! Loop over values and reindex
      DO i = 1, n
         a(i) = b(j(i))
      END DO

   END SUBROUTINE reindex

END MODULE sorting
