MODULE camb_stuff

   USE constants
   USE basic_operations
   USE file_info

   IMPLICIT NONE

   PRIVATE

   ! Routines
   PUBLIC :: read_CAMB_Pk
   PUBLIC :: write_CAMB_Pk
   PUBLIC :: read_CAMB_Tk

   ! File columns and additional info
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
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
      INTEGER, INTENT(OUT) :: nk
      CHARACTER(len=*), INTENT(IN) :: infile
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i
      INTEGER :: u

      nk = file_length(infile, verbose=.FALSE.)
      nk = nk-CAMB_Pk_comment_lines
      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_PK: CAMB file: ', trim(infile)
         WRITE (*, *) 'READ_CAMB_PK: Number of points:', nk
      END IF

      ALLOCATE (k(nk), Pk(nk))

      OPEN (newunit=u, file=infile)
      DO i = 1-CAMB_Pk_comment_lines, nk
         IF (i <= 0) THEN
            READ (u, *)
         ELSE
            READ (u, *) k(i), Pk(i)
         END IF
      END DO
      CLOSE (u)

      ! Convert from P(k) [(Mpc/h)^3] to Delta^2(k)
      Pk = 4.*pi*Pk*(k**3)/twopi**3

      IF(present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_CAMB_PK: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE read_CAMB_Pk

   SUBROUTINE write_CAMB_Pk(k, Pk, outfile, header)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: Pk(:)
      CHARACTER(len=256), INTENT(IN) :: outfile
      LOGICAL, INTENT(IN) :: header
      INTEGER :: u
      INTEGER :: ik

      OPEN(newunit=u, file=outfile)
      IF (header) WRITE(u, *) '#           k/h    P'
      DO ik = 1, size(k)
         WRITE(u, *) k(ik), Pk(ik)/(4.*pi*(k(ik)/twopi)**3)
      END DO
      CLOSE(u)

   END SUBROUTINE write_CAMB_Pk

   SUBROUTINE read_CAMB_Tk(k, Tk, nk, nTk, infile, verbose)

      ! Read in CAMB format T(k) file
      ! Note that this assumes there is one comment line
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
