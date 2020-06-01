PROGRAM sorting_test

   USE sorting

   IMPLICIT NONE
   LOGICAL :: ifail

   WRITE(*, *)

   CALL sort_test(ifail)

   CONTAINS

   SUBROUTINE sort_test(fail)

      LOGICAL, INTENT(OUT) :: fail
      REAL :: unsorted(10), input(10), sorted(10)
      INTEGER :: isort
      INTEGER, PARAMETER :: nsort = 3

      unsorted = [1.5, 100., -4., 1., 1., 1.2, 25., 81., 1e-7, 1e-1]
      sorted = [-4., 1e-7, 1e-1, 1., 1., 1.2, 1.5, 25., 81., 100.]

      fail = .FALSE.

      DO isort = 1, nsort

         input = unsorted

         CALL sort(unsorted, isort)

         IF (.NOT. all(unsorted == sorted)) fail = .TRUE.

         IF (fail) THEN
            WRITE(*, *) 'Sort test failed', isort
            EXIT
         END IF

      END DO

      IF (.NOT. fail) THEN
         WRITE(*, *) 'SORT_TEST: Passed' 
         WRITE(*, *)
      END IF

   END SUBROUTINE sort_test

END PROGRAM sorting_test