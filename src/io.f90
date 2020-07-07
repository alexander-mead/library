MODULE io

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: VD20_read_power

   CONTAINS

   SUBROUTINE VD20_read_power(k, z, Pk, nk, nz, name)

      USE file_info
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: z(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      INTEGER, INTENT(OUT) :: nk
      INTEGER, INTENT(OUT) :: nz
      CHARACTER(len=*), INTENT(IN) :: name
      CHARACTER(len=128), PARAMETER :: dir = '/Users/Mead/Physics/data/VD20'
      CHARACTER(len=256) :: infile
      INTEGER :: i, n, ik, iz
      INTEGER :: u
      REAL :: zin, zrem
      REAL :: crap
      INTEGER, PARAMETER :: header_size = 1
      LOGICAL, PARAMETER :: verbose = .FALSE.

      ! Create the input file name
      infile = trim(dir)//'/'//trim(name)//'.dat'
      CALL check_file_exists(infile)

      ! Count the length of the file
      n = file_length(infile, verbose)
      n = n-header_size

      ! Should be initialised
      zrem = -1.

      ! Calculate how many k and z values there are from the file
      OPEN(newunit=u, file=infile, status='old')
      READ(u, *)
      nk = 0
      DO
         nk = nk+1
         READ(u, *) zin
         IF (nk == 1) THEN
            zrem = zin
         ELSE IF (zin .NE. zrem) THEN
            nk = nk-1
            EXIT
         END IF
      END DO
      CLOSE(u)
      IF (verbose) WRITE(*, *) 'VD20_READ_POWER: nk:', nk     
      IF(mod(n, nk) .NE. 0) STOP 'VD20_READ_POWER: Error, could not figure out how many k, z values there were'
      nz = n/nk
      IF (verbose) WRITE(*, *) 'VD20_READ_POWER: nk:', nz

      ! Allocate arrays for k, z and power
      ALLOCATE(k(nk), z(nz), Pk(nk, nz))

      ! Fill arrays
      ik = 0
      iz = 0
      OPEN(newunit=u, file=infile, status='old')
      READ(u, *) ! Read header
      DO i = 1, n       
         IF(MOD(i-1, nk) == 0) THEN
            ik = 0
            iz = iz+1
         END IF
         ik = ik+1
         READ(u, *) z(iz), k(ik), crap, Pk(ik, iz)
      END DO
      CLOSE(u)

   END SUBROUTINE VD20_read_power

END MODULE io