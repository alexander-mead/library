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
   PUBLIC :: isort_quick

   INTERFACE index
      MODULE PROCEDURE index_real
      MODULE PROCEDURE index_int
   END INTERFACE index

   INTERFACE reindex
      MODULE PROCEDURE reindex_real
      MODULE PROCEDURE reindex_int
   END INTERFACE reindex

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
   INTEGER, PARAMETER :: isort_quick = 4

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

   SUBROUTINE index_real(a, idx, isort)

      ! Index the array 'a' from lowest to highest value
      REAL, INTENT(IN) :: a(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: idx(:)
      INTEGER, INTENT(IN) :: isort

      ALLOCATE(idx(size(a)))
      IF (isort == isort_stupid) THEN
         CALL stupid_index_real(a, idx)
      ELSE IF (isort == isort_bubble) THEN
         CALL bubble_index_real(a, idx)
      ELSE
         STOP 'INDEX_REAL: Error, isort specified incorrectly'
      END IF

   END SUBROUTINE index_real

   SUBROUTINE index_int(a, idx, isort)

      ! Index the array 'a' from lowest to highest value  
      INTEGER, INTENT(IN) :: a(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: idx(:)
      INTEGER, INTENT(IN) :: isort
      
      ALLOCATE(idx(size(a)))
      IF (isort == isort_stupid) THEN
         CALL stupid_index_int(a, idx)
      ELSE IF (isort == isort_bubble) THEN
         CALL bubble_index_int(a, idx)
      ELSE IF (isort == isort_quick) THEN
         CALL quick_index(size(a), a, idx)
      ELSE
         STOP 'INDEX_REAL: Error, isort specified incorrectly'
      END IF

   END SUBROUTINE index_int

   SUBROUTINE bubble_index_real(a, idx)

      ! Create an index array for a(:) that indexes from smallest to largest value
      REAL, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: idx(:)
      INTEGER :: i, isort, hold, n

      n = size(a)
      IF (n /= size(idx)) STOP 'BUBBLE_INDEX_REAL: Error, a and ind must have the same size'

      DO i = 1, n
         idx(i) = i
      END DO

      DO
         isort = 0
         DO i = 1, n-1
            IF (a(idx(i)) > a(idx(i+1))) THEN
               hold = idx(i+1)
               idx(i+1) = idx(i)
               idx(i) = hold
               isort = 1
            END IF
         END DO
         IF (isort == 0) EXIT
      END DO

   END SUBROUTINE bubble_index_real

   SUBROUTINE bubble_index_int(a, idx)!,verbose)

      ! Create an index array for integer a(:) that indexes from smallest to largest value
      INTEGER, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: idx(:)
      INTEGER :: i, isort, hold, n

      n = size(a)
      IF (n /= size(idx)) STOP 'BUBBLE_INDEX_INT: Error, a and idx must have the same size'

      DO i = 1, n
         idx(i) = i
      END DO

      DO
         isort = 0
         DO i = 1, n-1
            IF (a(idx(i)) > a(idx(i+1))) THEN
               hold = idx(i+1)
               idx(i+1) = idx(i)
               idx(i) = hold
               isort = 1
            END IF
         END DO
         IF (isort == 0) EXIT
      END DO

   END SUBROUTINE bubble_index_int

   SUBROUTINE stupid_index_real(a, idx)
  
      REAL, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: idx(:)
      INTEGER :: i, j, n
      REAL, ALLOCATABLE :: b(:)

      n = size(a)
      IF (n /= size(idx)) STOP 'STUPID_INDEX_REAL: Error, a and idx must have the same size'
      ALLOCATE(b(n))

      b = a

      ! This is probably stupid
      DO i = 1, n
         j = MINLOC(b, 1)
         idx(i) = j
         b(j) = HUGE(b)
      END DO

   END SUBROUTINE stupid_index_real

   SUBROUTINE stupid_index_int(a, idx)

      INTEGER, INTENT(IN) :: a(:)
      INTEGER, INTENT(OUT) :: idx(:)
      INTEGER :: i, j, n
      INTEGER, ALLOCATABLE :: b(:)

      n = size(a)
      IF (n /= size(idx)) STOP 'STUPID_INDEX_INT: Error, a and idx must have the same size'
      ALLOCATE(b(n))

      b = a

      ! This is probably stupid
      DO i = 1, n
         j = MINLOC(b, 1)
         idx(i) = j
         b(j) = HUGE(b)
      END DO

   END SUBROUTINE stupid_index_int

   SUBROUTINE quick_index(n, arr, indx)
      
      ! I think this is originally from numerial recipes, not sure
      INTEGER :: n, indx(n), M, NSTACK, arr(n), a
      PARAMETER (M=7, NSTACK=50)
      INTEGER :: i, indxt, ir, itemp, j, jstack, k, l, istack(NSTACK)
      do 11 j=1,n
         indx(j)=j
11    end do
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
         do 13 j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do 12 i=j-1,l,-1
               if(arr(indx(i)).le.a)goto 2
               indx(i+1)=indx(i)
12          end do
            i=l-1
2           indx(i+1)=indxt
13       end do
         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l)).gt.arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1)).gt.arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l)).gt.arr(indx(l+1)))then
            itemp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l+1)
         a=arr(indxt)
3        continue
         i=i+1
         if(arr(indx(i)).lt.a)goto 3
4        continue
         j=j-1
         if(arr(indx(j)).gt.a)goto 4
         if(j.lt.i)goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
5        indx(l+1)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack.gt.NSTACK) STOP 'NSTACK too small in index'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1

   END SUBROUTINE quick_index

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

   SUBROUTINE reindex_real(a, j)

      ! Reindex the array 'a' with the new indices 'j'
      REAL, INTENT(INOUT) :: a(:) ! Input array to check
      INTEGER, INTENT(IN) :: j(:) ! Input array indices to check
      REAL, ALLOCATABLE :: b(:)
      INTEGER :: i, n

      n = size(a)
      IF (n /= size(j)) STOP 'REINDEX: Error, a and j must have the same size'
      ALLOCATE(b(n))

      b = a  ! Store the input array
      a = 0. ! Delete the input values

      ! Loop over values and reindex
      DO i = 1, n
         a(i) = b(j(i))
      END DO

   END SUBROUTINE reindex_real

   SUBROUTINE reindex_int(a, j)

      ! Reindex the array 'a' with the new indices 'j'
      INTEGER, INTENT(INOUT) :: a(:) ! Input array to check
      INTEGER, INTENT(IN) :: j(:) ! Input array indices to check
      INTEGER, ALLOCATABLE :: b(:)
      INTEGER :: i, n

      n = size(a)
      IF (n /= size(j)) STOP 'REINDEX: Error, a and j must have the same size'
      ALLOCATE(b(n))

      b = a  ! Store the input array
      a = 0. ! Delete the input values

      ! Loop over values and reindex
      DO i = 1, n
         a(i) = b(j(i))
      END DO

   END SUBROUTINE reindex_int

END MODULE sorting
