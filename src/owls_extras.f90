MODULE owls_extras

   USE owls

CONTAINS

   CHARACTER(len=256) FUNCTION BAHAMAS_power_file_name(model, m, z, ip)

      USE HMx
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: model
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: z
      INTEGER, INTENT(IN) :: ip(2)
      CHARACTER(len=64) :: dir
      CHARACTER(len=32) :: snap, field(2), f1, f2, mesh
      LOGICAL :: lexist
      INTEGER :: j
      INTEGER, PARAMETER :: computer = 1

      ! Directory containing everything
      IF (computer == 1) dir = '/Users/Mead/Physics/BAHAMAS/power/'!M1536'
      IF (computer == 2) dir = '/home/amead/BAHAMAS/power/'!M1536'

      ! Convert the mesh size to a string and append to directory
      WRITE (mesh, *) m
      dir = TRIM(dir)//'M'//trim(adjustl(mesh))

      snap = BAHAMAS_snapshot(z)

      ! Set the fields
      DO j = 1, 2
         IF (ip(j) == field_matter) THEN
            field(j) = 'all'
         ELSE IF (ip(j) == field_cdm) THEN
            field(j) = 'dm'
         ELSE IF (ip(j) == field_gas) THEN
            field(j) = 'gas'
         ELSE IF (ip(j) == field_star) THEN
            field(j) = 'stars'
         ELSE IF (ip(j) == field_electron_pressure) THEN
            field(j) = 'epressure'
         ELSE
            WRITE (*, *) 'BAHAMAS_POWER_FILE_NAME: Field number', j
            WRITE (*, *) 'BAHAMAS_POWER_FILE_NAME: Halo type', ip(j)
            STOP 'BAHAMAS_POWER_FILE_NAME: Error, ip specified incorrectly'
         END IF
      END DO

      DO j = 1, 2

         IF (j == 1) THEN
            f1 = field(1)
            f2 = field(2)
         ELSE IF (j == 2) THEN
            f1 = field(2)
            f2 = field(1)
         ELSE
            STOP 'BAHAMAS_POWER_FILE_NAME: Error, something went wrong'
         END IF

         ! File name
         BAHAMAS_power_file_name = &
            trim(dir)//'/'//trim(model)//'_L400N1024_WMAP9_'//trim(snap)//'_'//trim(f1)//'_'//trim(f2)//'_power.dat'

         ! Check it exists
         INQUIRE (file=BAHAMAS_power_file_name, exist=lexist)

         IF (lexist) THEN
            EXIT
         ELSE IF (j == 2) THEN
            WRITE (*, *) 'BAHAMAS_POWER_FILE_NAME: ', trim(BAHAMAS_power_file_name)
            STOP 'BAHAMAS_POWER_FILE_NAME: Error, file does not exist'
         END IF

      END DO

   END FUNCTION BAHAMAS_power_file_name

   SUBROUTINE read_BAHAMAS_power(k, Pk, nk, z, name, mesh, field, cosm, kmin, kmax, cut_nyquist, subtract_shot, response, verbose)

      USE basic_operations
      USE cosmology_functions
      USE HMx

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
      INTEGER, INTENT(OUT) :: nk
      REAL, INTENT(IN) :: z
      CHARACTER(len=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: mesh
      INTEGER, INTENT(IN) :: field(2)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, OPTIONAL, INTENT(IN) :: kmin, kmax
      LOGICAL, OPTIONAL, INTENT(IN) :: cut_nyquist
      LOGICAL, OPTIONAL, INTENT(IN) :: subtract_shot
      LOGICAL, OPTIONAL, INTENT(IN) :: response
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: Pk_DM(:), Pk_HMcode(:,:)
      REAL :: a(1)
      INTEGER, PARAMETER :: na=1
      CHARACTER(len=256) :: infile, dmonly
      INTEGER, PARAMETER :: field_all_matter(2) = field_matter
      CHARACTER(len=256) :: response_infile = 'DMONLY_2fluid_nu0'
      
      infile = BAHAMAS_power_file_name(name, mesh, z, field)
      CALL read_simulation_power_spectrum(k, Pk, nk, infile, kmin, kmax, cut_nyquist, subtract_shot, verbose)

      IF (present_and_correct(response)) THEN
         dmonly = BAHAMAS_power_file_name(response_infile, mesh, z, field_all_matter)
         CALL read_simulation_power_spectrum(k, Pk_DM, nk, dmonly, kmin, kmax, cut_nyquist, subtract_shot, verbose)
         Pk = Pk/Pk_DM
         a = scale_factor_z(z)
         ALLOCATE (Pk_HMcode(nk,na))
         CALL calculate_HMcode(k, a, Pk_HMcode, nk, na, cosm)
         Pk = Pk*Pk_HMcode(:,1)
      END IF

   END SUBROUTINE read_BAHAMAS_power

   SUBROUTINE read_simulation_power_spectrum(k, Pk, n, infile, kmin, kmax, cut_nyquist, subtract_shot, verbose)

      USE file_info
      USE array_operations
      USE basic_operations

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:), Pk(:) ! Output simulation k and power
      INTEGER, INTENT(OUT) :: n ! Number of output k values
      CHARACTER(len=*), INTENT(IN) :: infile ! Input file location
      REAL, OPTIONAL, INTENT(IN) :: kmin, kmax ! Minimum and maximum k values to cut at
      LOGICAL, OPTIONAL, INTENT(IN) :: cut_nyquist ! Logical to cut Nyquist or not
      LOGICAL, OPTIONAL, INTENT(IN) :: subtract_shot ! Logical to subtract shot noise or not
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose ! Logical verbose
      INTEGER :: i, j
      REAL :: shot, kbig
      LOGICAL :: lexist

      ! Deallocate arrays if they are already allocated
      IF (ALLOCATED(k)) DEALLOCATE (k)
      IF (ALLOCATED(Pk)) DEALLOCATE (Pk)

      ! Check file exists
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) THEN
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: File: ', trim(infile)
         STOP 'READ_SIMULATION_POWER_SPECTRUM: File does not exist'
      END IF

      ! Get file length and allocate arrays for output
      n = file_length(infile, verbose=.FALSE.)
      ALLOCATE (k(n), Pk(n))

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Reading in data'
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: File: ', trim(infile)
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Initial nk:', n
      END IF

      ! Read in data from file
      OPEN (9, file=infile, status='old')
      DO i = 1, n
         READ (9, *) k(i), Pk(i), shot
         IF (present_and_correct(subtract_shot)) Pk(i) = Pk(i)-shot
      END DO
      CLOSE (9)

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: kmin [h/Mpc]:', k(1)
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: kmax [h/Mpc]:', k(n)
      END IF

      IF (present_and_correct(cut_nyquist)) THEN
         ! Find position in array of half-Nyquist
         kbig = k(n)
         DO i = 1, n
            IF (k(i) > kbig/2.) EXIT
         END DO
         ! Cut arrays down to half-Nyquist
         CALL amputate_array(k, n, 1, i)
         CALL amputate_array(Pk, n, 1, i)
         n = i
         ! Write to screen
         IF (present_and_correct(verbose)) THEN
            WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Trimmed to Nyquist frequency'
            WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: New kmax [h/Mpc]:', k(n)
         END IF
      END IF

      IF (PRESENT(kmin)) THEN
         j = 0
         DO i = 1, n-1
            IF (k(i) < kmin .AND. k(i+1) > kmin) THEN
               j = i
            END IF
         END DO
         IF (j .NE. 0) THEN
            CALL amputate_array(k, n, j, n)
            CALL amputate_array(Pk, n, j, n)
            n = n-j+1
         END IF
         IF (present_and_correct(verbose)) THEN
            WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Trimmed to new kmin'
            WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: New kmin [h/Mpc]:', k(1)
         END IF
      END IF

      IF (PRESENT(kmax)) THEN
         j = n
         DO i = 1, n-1
            IF (k(i) < kmax .AND. k(i+1) > kmax) THEN
               j = i
            END IF
         END DO
         IF (j .NE. n) THEN
            CALL amputate_array(k, n, 1, j)
            CALL amputate_array(Pk, n, 1, j)
            !n=m
            n = j
         END IF
         IF (present_and_correct(verbose)) THEN
            WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Trimmed to new kmax'
            WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: New kmax [h/Mpc]:', k(n)
         END IF
      END IF

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         !DO i=1,n
         !   WRITE(*,*) k(i), Pk(i)
         !END DO
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Final nk:', n
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE read_simulation_power_spectrum

END MODULE owls_extras
