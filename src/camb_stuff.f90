MODULE camb_stuff

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_CAMB_Pk

CONTAINS

   SUBROUTINE read_CAMB_Pk(k, p, n, infile, verbose)

      ! Read in CAMB format P(k) file
      ! Note that this assumes there is one comment line
      USE file_info
      USE constants
      USE logical_operations
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: p(:)
      INTEGER, INTENT(OUT) :: n
      CHARACTER(len=*), INTENT(IN) :: infile
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i

      n = file_length(infile, verbose=.FALSE.)
      n = n-1
      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_PK: CAMB file: ', trim(infile)
         WRITE (*, *) 'READ_CAMB_PK: Number of points:', n
         WRITE (*, *)
      END IF

      ALLOCATE (k(n), p(n))

      OPEN (7, file=infile)
      DO i = 0, n
         IF (i == 0) THEN
            READ (7, *)
         ELSE
            READ (7, *) k(i), p(i)
         END IF
      END DO
      CLOSE (7)

      ! Convert from P(k) to Delta^2(k)
      p = 4.*pi*p*(k**3)/twopi**3

      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_PK: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE read_CAMB_Pk

END MODULE camb_stuff
