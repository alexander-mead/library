MODULE io

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_VD20_power

   CONTAINS

   SUBROUTINE read_VD20_power(k, z, Pk, nk, nz, name)

      USE file_info
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: z(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:,:)
      INTEGER, INTENT(OUT) :: nk
      INTEGER, INTENT(OUT) :: nz
      CHARACTER(len=*), INTENT(IN) :: name
      CHARACTER(len=128), PARAMETER :: dir = '/Users/Mead/Physics/data/VD20/data'
      CHARACTER(len=256) :: infile
      INTEGER :: i, n, ik, iz
      REAL :: crap
      INTEGER, PARAMETER :: header_size = 1
      INTEGER, PARAMETER :: expected_number_of_lines = 5280
      INTEGER, PARAMETER :: nk_expected = 352
      INTEGER, PARAMETER :: nz_expected = 15

      ! Create the input file name
      infile = trim(dir)//'/'//trim(name)//'.dat'

      ! Count the length of the file
      n = file_length(infile, verbose=.TRUE.)
      n = n-header_size

      ! Check that the file length is consistent with expectations
      IF (n < expected_number_of_lines) THEN
         WRITE(*,*) 'READ_VD20_POWER: Name: ', trim(name)
         WRITE(*,*) 'READ_VD20_POWER: File length: ', n
         STOP 'READ_VD20_POWER: Error, file size is unexpected'
      END IF

      ! This should be true unless the files are fucked
      nk = nk_expected
      nz = nz_expected

      ! Allocate arrays for k, z and power
      ALLOCATE(k(nk), z(nz), Pk(nk, nz))

      ! Fill arrays
      ik = 0
      iz = 0
      OPEN(7, file=infile)
      READ(7, *)
      DO i = 1, n       
         IF(MOD(i-1, nk) == 0) THEN
            ik = 0
            iz = iz+1
         END IF
         ik = ik+1
         READ(7, *) z(iz), k(ik), crap, Pk(ik, iz)
      END DO
      CLOSE(7)

   END SUBROUTINE read_VD20_power

END MODULE io