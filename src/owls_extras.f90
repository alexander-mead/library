MODULE owls_extras

   USE owls

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BAHAMAS_read_power
   PUBLIC :: BAHAMAS_get_ks
   PUBLIC :: VD20_get_power
   PUBLIC :: VD20_get_response

   INTEGER, PARAMETER :: computer_mac = 1
   INTEGER, PARAMETER :: computer_linux = 2

   ! Set computer
   INTEGER, PARAMETER :: computer = computer_mac

CONTAINS

   SUBROUTINE BAHAMAS_get_ks(mesh, k, nk)

      USE file_info
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mesh
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      INTEGER, INTENT(OUT) :: nk
      INTEGER :: i
      CHARACTER(len=256) :: dir, infile
      CHARACTER(len=256), PARAMETER :: kfile = 'DMONLY_nu0_L400N1024_WMAP9_snap32_all_all_power.dat'

      dir = BAHAMAS_dir(mesh)
      infile = trim(dir)//'/'//trim(kfile)
      nk = file_length(infile, verbose=.FALSE.)
      
      ALLOCATE(k(nk))
      OPEN(7, file=infile)
      DO i = 1, nk
         READ(7, *) k(i)
      END DO
      CLOSE(7)

   END SUBROUTINE BAHAMAS_get_ks

   CHARACTER(len=64) FUNCTION BAHAMAS_dir(m)

      USE string_operations
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      CHARACTER(len=8) :: mesh
      CHARACTER(len=64) :: dir

      ! Directory containing everything
      IF (computer == computer_mac)   dir = '/Users/Mead/Physics/BAHAMAS/power/'
      IF (computer == computer_linux) dir = '/home/amead/BAHAMAS/power/'

      ! Convert the mesh size to a string and append to directory
      mesh = integer_to_string(m)
      BAHAMAS_dir = TRIM(dir)//'M'//trim(mesh)

   END FUNCTION BAHAMAS_dir

   CHARACTER(len=32) FUNCTION BAHAMAS_field_name(i)

      USE HMx
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i

      IF (i == field_matter) THEN
         BAHAMAS_field_name = 'all'
      ELSE IF (i == field_cdm) THEN
         BAHAMAS_field_name = 'dm'
      ELSE IF (i == field_gas) THEN
         BAHAMAS_field_name = 'gas'
      ELSE IF (i == field_stars) THEN
         BAHAMAS_field_name = 'stars'
      ELSE IF (i == field_electron_pressure) THEN
         BAHAMAS_field_name = 'epressure'
      ELSE
         WRITE (*, *) 'BAHAMAS_FIELD_NAME: Field integer:', i
         STOP 'BAHAMAS_FIELD_NAME: Error, field integer specified incorrectly'
      END IF

   END FUNCTION BAHAMAS_field_name

   CHARACTER(len=256) FUNCTION BAHAMAS_power_file_name(model, m, z, ip)

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: model
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: z
      INTEGER, INTENT(IN) :: ip(2)
      CHARACTER(len=64) :: dir
      CHARACTER(len=32) :: snap, field(2), f1, f2
      LOGICAL :: lexist
      INTEGER :: j

      dir = BAHAMAS_dir(m)
      snap = BAHAMAS_snapshot(z)
      DO j = 1, 2
         field(j) = BAHAMAS_field_name(ip(j))
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

   CHARACTER(len=256) FUNCTION BAHAMAS_error_file_name(m, z, ip)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: z
      INTEGER, INTENT(IN) :: ip(2)
      CHARACTER(len=64) :: dir
      CHARACTER(len=32) :: snap, field(2), f1, f2
      LOGICAL :: lexist
      INTEGER :: j

      dir = BAHAMAS_dir(m)
      snap = BAHAMAS_snapshot(z)
      DO j = 1, 2
         field(j) = BAHAMAS_field_name(ip(j))
      END DO

      DO j = 1, 2

         IF (j == 1) THEN
            f1 = field(1)
            f2 = field(2)
         ELSE IF (j == 2) THEN
            f1 = field(2)
            f2 = field(1)
         ELSE
            STOP 'BAHAMAS_ERROR_FILE_NAME: Error, something went wrong'
         END IF

         ! File name
         BAHAMAS_error_file_name = &
            trim(dir)//'/L400N1024_WMAP9_'//trim(snap)//'_'//trim(f1)//'_'//trim(f2)//'_error.dat'

         ! Check it exists
         INQUIRE (file=BAHAMAS_error_file_name, exist=lexist)

         IF (lexist) THEN
            EXIT
         ELSE IF (j == 2) THEN
            WRITE (*, *) 'BAHAMAS_ERROR_FILE_NAME: ', trim(BAHAMAS_error_file_name)
            STOP 'BAHAMAS_ERROR_FILE_NAME: Error, file does not exist'
         END IF

      END DO

   END FUNCTION BAHAMAS_error_file_name

   SUBROUTINE BAHAMAS_read_power(k, Pk, Er, nk, z, name, mesh, field, cosm, kmin, kmax, &
      cut_nyquist, subtract_shot, realisation_errors, response, verbose)

      USE basic_operations
      USE cosmology_functions
      USE HMx

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Er(:)
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
      LOGICAL, OPTIONAL, INTENT(IN) :: realisation_errors
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: Pk_DM(:), crap(:), Pk_HMcode(:,:)
      REAL :: a(1)
      INTEGER, PARAMETER :: na=1
      CHARACTER(len=256) :: infile, dmonly
      INTEGER, PARAMETER :: field_all_matter(2) = field_matter
      CHARACTER(len=256) :: response_infile = 'DMONLY_2fluid_nu0'
      
      infile = BAHAMAS_power_file_name(name, mesh, z, field)
      CALL read_simulation_power_spectrum(k, Pk, Er, nk, infile, kmin, kmax, cut_nyquist, subtract_shot, verbose)
      IF(present_and_correct(realisation_errors)) THEN
         infile = BAHAMAS_error_file_name(mesh, z, field)
         CALL read_simulation_power_spectrum(k, crap, Er, nk, infile, kmin, kmax, cut_nyquist, subtract_shot, verbose)
      END IF

      IF (present_and_correct(response)) THEN
         dmonly = BAHAMAS_power_file_name(response_infile, mesh, z, field_all_matter)
         CALL read_simulation_power_spectrum(k, Pk_DM, crap, nk, dmonly, kmin, kmax, cut_nyquist, subtract_shot, verbose)
         Pk = Pk/Pk_DM
         a = scale_factor_z(z)
         ALLOCATE (Pk_HMcode(nk, na))
         CALL calculate_HMcode(k, a, Pk_HMcode, nk, na, cosm)
         Pk = Pk*Pk_HMcode(:,1)
      END IF

   END SUBROUTINE BAHAMAS_read_power

   SUBROUTINE read_simulation_power_spectrum(k, Pk, Er, nk, infile, kmin, kmax, cut_nyquist, subtract_shot, verbose)

      USE file_info
      USE basic_operations

      ! TODO: Move this to simulations.f90
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)  ! Output simulation k and power
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:) ! Output simulation k and power
      REAL, ALLOCATABLE, INTENT(OUT) :: Er(:) ! Output simulation k and power
      INTEGER, INTENT(OUT) :: nk              ! Number of output k values
      CHARACTER(len=*), INTENT(IN) :: infile ! Input file location
      REAL, OPTIONAL, INTENT(IN) :: kmin, kmax ! Minimum and maximum k values to cut at
      LOGICAL, OPTIONAL, INTENT(IN) :: cut_nyquist ! Logical to cut Nyquist or not
      LOGICAL, OPTIONAL, INTENT(IN) :: subtract_shot ! Logical to subtract shot noise or not
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose ! Logical verbose
      INTEGER :: i
      REAL :: shot, crap

      ! Deallocate arrays if they are already allocated
      IF (ALLOCATED(k))  DEALLOCATE (k)
      IF (ALLOCATED(Pk)) DEALLOCATE (Pk)
      IF (ALLOCATED(Er)) DEALLOCATE (Er)

      CALL check_file_exists(infile)

      ! Get file length and allocate arrays for output
      nk = file_length(infile, verbose=.FALSE.)
      ALLOCATE (k(nk), Pk(nk), Er(nk))

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Reading in data'
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: File: ', trim(infile)
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Initial nk:', nk
      END IF

      ! Read in data from file
      OPEN (9, file=infile, status='old')
      DO i = 1, nk
         READ (9, *) k(i), Pk(i), shot, crap, Er(i)
         IF (present_and_correct(subtract_shot)) Pk(i) = Pk(i)-shot
      END DO
      CLOSE (9)

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: kmin [h/Mpc]:', k(1)
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: kmax [h/Mpc]:', k(nk)
      END IF

      IF (present_and_correct(cut_nyquist)) CALL cut_Nyquist_frequency(k, Pk, Er, nk, verbose)
      IF (PRESENT(kmin)) CALL cut_kmin(kmin, k, Pk, Er, nk, verbose)
      IF (PRESENT(kmax)) CALL cut_kmax(kmax, k, Pk, Er, nk, verbose)

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Final nk:', nk
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE read_simulation_power_spectrum

   SUBROUTINE cut_Nyquist_frequency(k, Pk, Er, nk, verbose)

      USE basic_operations
      USE array_operations
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Pk(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Er(:)
      INTEGER, INTENT(INOUT) :: nk
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i
      REAL :: kbig

      ! Find position in array of half-Nyquist
      kbig = k(nk)
      DO i = 1, nk
         IF (k(i) > kbig/2.) EXIT
      END DO

      ! Cut arrays down to half-Nyquist
      CALL amputate_array(k, nk, 1, i)
      CALL amputate_array(Pk, nk, 1, i)
      CALL amputate_array(Er, nk, 1, i)
      nk = i

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: Trimmed to Nyquist frequency'
         WRITE (*, *) 'READ_SIMULATION_POWER_SPECTRUM: New kmax [h/Mpc]:', k(nk)
      END IF

   END SUBROUTINE cut_Nyquist_frequency

   SUBROUTINE cut_kmin(kmin, k, Pk, Er, nk, verbose)

      USE basic_operations
      USE array_operations
      IMPLICIT NONE
      REAL, INTENT(IN) :: kmin
      REAL, ALLOCATABLE, INTENT(INOUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Pk(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Er(:)
      INTEGER, INTENT(INOUT) :: nk
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i, j

      j = 0
      DO i = 1, nk-1
         IF (k(i) < kmin .AND. k(i+1) > kmin) THEN
            j = i
         END IF
      END DO
      IF (j .NE. 0) THEN
         CALL amputate_array(k, nk, j, nk)
         CALL amputate_array(Pk, nk, j, nk)
         CALL amputate_array(Er, nk, j, nk)
         nk = nk-j+1
      END IF
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'CUT_KMIN: Trimmed to new kmin'
         WRITE (*, *) 'CUT_KMIN: New kmin [h/Mpc]:', k(1)
         WRITE (*, *)
      END IF

   END SUBROUTINE cut_kmin

   SUBROUTINE cut_kmax(kmax, k, Pk, Er, nk, verbose)

      USE basic_operations
      USE array_operations
      IMPLICIT NONE
      REAL, INTENT(IN) :: kmax
      REAL, ALLOCATABLE, INTENT(INOUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Pk(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Er(:)
      INTEGER, INTENT(INOUT) :: nk
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i, j
      
      j = nk
      DO i = 1, nk-1
         IF (k(i) < kmax .AND. k(i+1) > kmax) THEN
            j = i
         END IF
      END DO
      IF (j .NE. nk) THEN
         CALL amputate_array(k, nk, 1, j)
         CALL amputate_array(Pk, nk, 1, j)
         CALL amputate_array(Er, nk, 1, j)
         nk = j
      END IF
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'CUT_KMAX: Trimmed to new kmax'
         WRITE (*, *) 'CUT_KMAX: New kmax [h/Mpc]:', k(nk)
      END IF  

   END SUBROUTINE cut_kmax

   SUBROUTINE VD20_get_power(k, Pk, nk, z, name)

      USE io
      USE array_operations
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
      INTEGER, INTENT(OUT) :: nk
      REAL, INTENT(IN) :: z
      CHARACTER(len=*), INTENT(IN) :: name
      REAL, ALLOCATABLE :: zs(:), P(:, :)
      INTEGER :: iz, nz
      REAL, PARAMETER :: eps = 1e-4

      CALL read_VD20_power(k, zs, P, nk, nz, name)

      iz = array_position(z, zs, nz, eps)

      ALLOCATE(Pk(nk))
      Pk = P(:, iz)

   END SUBROUTINE VD20_get_power

   SUBROUTINE VD20_get_response(k, Rk, nk, z, name)

      USE io
      USE array_operations
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Rk(:)
      INTEGER, INTENT(OUT) :: nk
      REAL, INTENT(IN) :: z
      CHARACTER(len=*), INTENT(IN) :: name
      REAL, ALLOCATABLE :: Pk(:), Pk_dmonly(:)
      CHARACTER(len=256) :: name_dmonly

      name_dmonly = VD20_dmonly_counterpart(name)

      CALL get_VD20_power(k, Pk, nk, z, name)
      CALL get_VD20_power(k, Pk_dmonly, nk, z, name_dmonly)

      ALLOCATE(Rk(nk))

      Rk = Pk/Pk_dmonly

   END SUBROUTINE VD20_get_response

   CHARACTER(len=256) FUNCTION VD20_dmonly_counterpart(name)

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: name

      IF (trim(name) == 'BAHAMAS_Theat7.6_nu0_WMAP9') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_Theat8.0_nu0_WMAP9') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.06_Planck2015') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.06_Planck2015_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.06_WMAP9') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.06_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.12_Planck2015') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.12_Planck2015_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.12_WMAP9') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.12_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.24_Planck2015') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.24_Planck2015_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.24_WMAP9') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.24_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.48_Planck2015') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.48_Planck2015_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0.48_WMAP9') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0.48_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0_BAO_L200N512') THEN
         VD20_dmonly_counterpart = 'DMONLY_nu0_BAO_L200N512'
      ELSE IF (trim(name) == 'BAHAMAS_nu0_Planck2013') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0_Planck2013_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0_WMAP9') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0_WMAP9_L100N512') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0_WMAP9_L100N512'
      ELSE IF (trim(name) == 'BAHAMAS_nu0_WMAP9_v2') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0_v2_WMAP9_L400N1024'
      ELSE IF (trim(name) == 'BAHAMAS_nu0_WMAP9_v3') THEN
         VD20_dmonly_counterpart = 'DMONLY_2fluid_nu0_v3_WMAP9_L400N1024'
      ELSE
         WRITE (*, *) 'VD20_DMONLY_COUNTERPART: Simulation: ', trim(name)
         STOP 'VD20_DMONLY_COUNTERPART: Error, no counterpart found'
      END IF

   END FUNCTION VD20_dmonly_counterpart

END MODULE owls_extras
