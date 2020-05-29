MODULE camb_stuff

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_CAMB_Pk
   PUBLIC :: read_CAMB_Tk
   PUBLIC :: CAMB_column_Tk_CDM
   PUBLIC :: CAMB_column_Tk_baryon
   PUBLIC :: CAMB_column_Tk_photon
   PUBLIC :: CAMB_column_Tk_nu
   PUBLIC :: CAMB_column_Tk_massive_nu
   PUBLIC :: CAMB_column_Tk_total
   PUBLIC :: CAMB_Pk_comment_lines
   PUBLIC :: CAMB_Tk_comment_lines

   INTEGER, PARAMETER :: CAMB_Pk_comment_lines = 1
   INTEGER, PARAMETER :: CAMB_Tk_comment_lines = 1
   INTEGER, PARAMETER :: CAMB_number_Tk = 6
   INTEGER, PARAMETER :: CAMB_column_Tk_CDM = 1
   INTEGER, PARAMETER :: CAMB_column_Tk_baryon = 2
   INTEGER, PARAMETER :: CAMB_column_Tk_photon = 3
   INTEGER, PARAMETER :: CAMB_column_Tk_nu = 4
   INTEGER, PARAMETER :: CAMB_column_Tk_massive_nu = 5
   INTEGER, PARAMETER :: CAMB_column_Tk_total = 6

CONTAINS

   SUBROUTINE read_CAMB_Pk(k, Pk, nk, infile, verbose)

      ! Read in CAMB format P(k) file
      ! Note that this assumes there is one comment line
      USE file_info
      USE constants
      USE basic_operations
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
      INTEGER, INTENT(OUT) :: nk
      CHARACTER(len=*), INTENT(IN) :: infile
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i

      nk = file_length(infile, verbose=.FALSE.)
      nk = nk-CAMB_Pk_comment_lines
      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_PK: CAMB file: ', trim(infile)
         WRITE (*, *) 'READ_CAMB_PK: Number of points:', nk
         WRITE (*, *)
      END IF

      ALLOCATE (k(nk), Pk(nk))

      OPEN (7, file=infile)
      DO i = 1-CAMB_Pk_comment_lines, nk
         IF (i <= 0) THEN
            READ (7, *)
         ELSE
            READ (7, *) k(i), Pk(i)
         END IF
      END DO
      CLOSE (7)

      ! Convert from P(k) to Delta^2(k)
      Pk = 4.*pi*Pk*(k**3)/twopi**3

      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_PK: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE read_CAMB_Pk

   SUBROUTINE read_CAMB_Tk(k, Tk, nk, nTk, infile, verbose)

      ! Read in CAMB format T(k) file
      ! Note that this assumes there is one comment line
      USE file_info
      USE basic_operations
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Tk(:,:)
      INTEGER, INTENT(OUT) :: nk
      INTEGER, INTENT(OUT) :: nTk
      CHARACTER(len=*), INTENT(IN) :: infile
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i, j

      nk = file_length(infile, verbose=.FALSE.)
      nk = nk-CAMB_Tk_comment_lines
      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_Tcold: CAMB file: ', trim(infile)
         WRITE (*, *) 'READ_CAMB_Tcold: Number of points:', nk
         WRITE (*, *)
      END IF

      nTk = CAMB_number_Tk
      ALLOCATE (k(nk), Tk(nTk, nk))

      OPEN (7, file=infile)
      DO i = 1-CAMB_Tk_comment_lines, nk
         IF (i <= 0 ) THEN
            READ (7, *)
         ELSE
            READ (7, *) k(i), (Tk(j,i), j=1,nTk)
         END IF       
      END DO
      CLOSE (7)

      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_Tcold: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE read_CAMB_Tk

END MODULE camb_stuff
