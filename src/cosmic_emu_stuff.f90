MODULE cosmic_emu_stuff

   USE file_info
   USE array_operations
   USE cosmology_functions
   USE interpolate
   USE basic_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: get_CosmicEmu_power_z
   PUBLIC :: get_FrankenEmu_power_z
   PUBLIC :: get_MiraTitan_power_z
   PUBLIC :: get_emulator_power
   PUBLIC :: emulator_CosmicEmu
   PUBLIC :: emulator_FrankenEmu
   PUBLIC :: emulator_MiraTitan

   INTEGER, PARAMETER :: emulator_CosmicEmu = 1
   INTEGER, PARAMETER :: emulator_FrankenEmu = 2
   INTEGER, PARAMETER :: emulator_MiraTitan = 3

CONTAINS

   SUBROUTINE get_emulator_power(k, a, Pk, nk, cosm, emulator_version, rebin, verbose)

      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      INTEGER, INTENT(OUT) :: nk
      TYPE(cosmology), INTENT(IN) :: cosm   
      INTEGER, INTENT(IN) :: emulator_version
      LOGICAL, OPTIONAL, INTENT(IN) :: rebin
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: ia, na
      REAL :: z, kmax_rebin
      REAL, ALLOCATABLE :: k_emu(:), Pk_emu(:)
      REAL, PARAMETER :: kmin_rebin = 1e-2
      REAL, PARAMETER :: kmax_CosmicEmu_rebin = 1.   ! Maximum k if rebinnning
      REAL, PARAMETER :: kmax_FrankenEmu_rebin = 10. ! Maximum k if rebinnning
      REAL, PARAMETER :: kmax_MiraTitan_rebin = 7.   ! Maximum k if rebinnning
      INTEGER, PARAMETER :: nk_rebin = 128           ! Number of k points if rebinning
      INTEGER, PARAMETER :: iorder_rebin = 3
      INTEGER, PARAMETER :: ifind_rebin = 3
      INTEGER, PARAMETER :: iinterp_rebin = 2

      IF (present_and_correct(rebin)) THEN
         IF (emulator_version == emulator_CosmicEmu) THEN
            kmax_rebin = kmax_CosmicEmu_rebin
         ELSE IF (emulator_version == emulator_FrankenEmu) THEN
            kmax_rebin = kmax_FrankenEmu_rebin
         ELSE IF (emulator_version == emulator_MiraTitan) THEN
            kmax_rebin = kmax_MiraTitan_rebin
         ELSE
            STOP 'GET_EMULATOR_POWER: Error, emulator_version set incorrectly'
         END IF
      END IF

      na = size(a)

      DO ia = 1, na

         z = redshift_a(a(ia))

         IF (emulator_version == emulator_CosmicEmu) THEN
            CALL get_CosmicEmu_power_z(k_emu, Pk_emu, nk, z, cosm, verbose)
         ELSE IF (emulator_version == emulator_FrankenEmu) THEN
            CALL get_FrankenEmu_power_z(k_emu, Pk_emu, nk, z, cosm, verbose)
         ELSE IF (emulator_version == emulator_MiraTitan) THEN
            CALL get_MiraTitan_power_z(k_emu, Pk_emu, nk, z, cosm, verbose)
         ELSE
            STOP 'GET_EMULATOR_POWER: Error, emulator_version not specified correctly'
         END IF

         IF (present_and_correct(rebin)) THEN
            CALL rebin_array(kmin_rebin, kmax_rebin, nk_rebin, k_emu, Pk_emu, &
               iorder_rebin, &
               ifind_rebin, &
               iinterp_rebin, &
               logx=.TRUE., &
               logf=.TRUE. &
               )
            nk = nk_rebin
         END IF

         IF (.NOT. allocated(k))  ALLOCATE(k(nk))
         IF (.NOT. allocated(Pk)) ALLOCATE(Pk(nk, na))

         k = k_emu
         Pk(:, ia) = Pk_emu

      END DO

   END SUBROUTINE get_emulator_power

   SUBROUTINE get_CosmicEmu_power_z(k, P, n, z, cosm, verbose)
    
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i
      CHARACTER(len=256) :: crap
      REAL :: h
      INTEGER, PARAMETER :: nh = 10  ! Length of header
      INTEGER, PARAMETER :: h_li = 8 ! Line that h in on
      CHARACTER(len=256), PARAMETER :: params = 'emu_params.txt'
      CHARACTER(len=256), PARAMETER :: output = 'emu_power.dat'
      CHARACTER(len=256), PARAMETER :: exe = '/Users/Mead/Physics/CosmicEmu/emu.exe'
      REAL, PARAMETER :: eps_h = 3e-3

      ! Remove previous parameter and power file
      CALL EXECUTE_COMMAND_LINE('rm '//trim(params))
      CALL EXECUTE_COMMAND_LINE('rm '//trim(output))

      ! Write a new parameter file
      OPEN (7, file=params)
      WRITE (7, fmt='(A20,7F10.5)') trim(output), cosm%om_m*cosm%h**2, cosm%om_b*cosm%h**2, cosm%ns, cosm%sig8, cosm%w, z
      CLOSE (7)

      ! Run emu
      CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params))! > /dev/null')

      ! Get length of emu file
      n = file_length(output)
      n = n-nh
      IF (ALLOCATED(k)) DEALLOCATE (k)
      IF (ALLOCATED(P)) DEALLOCATE (P)
      ALLOCATE (k(n), P(n))

      ! Write useful things to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_COSMICEMU_POWER_Z: z:', z
         WRITE (*, *) 'GET_COSMICEMU_POWER_Z: P(k) file length:', n
         WRITE (*, *) 'GET_COSMICEMU_POWER_Z: Reading in P(k): ', trim(output)
      END IF

      ! Read in data file
      OPEN (7, file=output)
      DO i = 1-nh, n
         IF (i == h_li-nh) THEN
            READ (7, *) crap, crap, crap, crap, crap, crap, crap, crap, h
            IF (present_and_correct(verbose)) THEN
               WRITE (*, *) 'GET_COSMICEMU_POWER_Z: CMB derived h:', h
               WRITE (*, *) 'GET_COSMICEMU_POWER_Z: Cosmology h:', cosm%h
               WRITE (*, *) 'GET_COSMICEMU_POWER_Z: h ratio:', cosm%h/h
            END IF
            IF (.NOT. requal(h, cosm%h, eps_h)) STOP 'GET_COSMIC_EMU_POWER_Z: Error, h values differ'
         ELSE IF (i < 1) THEN
            READ (7, *)
         ELSE
            READ (7, *) k(i), P(i)
         END IF
      END DO
      CLOSE (7)

      ! Convert k to k/h
      k = k/cosm%h

      ! Done
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_COSMICEMU_POWER_Z: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_CosmicEmu_power_z

   SUBROUTINE get_FrankenEmu_power_z(k, P, n, z, cosm, verbose)

      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i
      INTEGER, PARAMETER :: nh = 5   ! Length of header
      CHARACTER(len=256), PARAMETER :: params = 'emu_params.txt'
      CHARACTER(len=256), PARAMETER :: output = 'emu_power.dat'
      CHARACTER(len=256), PARAMETER :: exe = '/Users/Mead/Physics/FrankenEmu/emu.exe'

      ! Remove previous parameter and power file
      CALL EXECUTE_COMMAND_LINE('rm '//trim(params))
      CALL EXECUTE_COMMAND_LINE('rm '//trim(output))

      ! Write a new parameter file
      OPEN (7, file=params)
      WRITE (7, fmt='(A20,7F10.5)') trim(output), cosm%om_m*cosm%h**2, cosm%om_b*cosm%h**2, cosm%ns, -cosm%w, cosm%sig8, cosm%h, z
      CLOSE (7)

      ! Run emu
      CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params))! > /dev/null')

      ! Get length of emu file
      n = file_length(output, verbose=.FALSE.)
      n = n-nh
      IF (ALLOCATED(k)) DEALLOCATE (k)
      IF (ALLOCATED(P)) DEALLOCATE (P)
      ALLOCATE (k(n), P(n))

      ! Write useful things to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_FRANKENEMU_POWER_Z: z:', z
         WRITE (*, *) 'GET_FRANKENEMU_POWER_Z: P(k) file length:', n
         WRITE (*, *) 'GET_FRANKENEMU_POWER_Z: Reading in P(k): ', trim(output)
      END IF

      ! Read in data file
      OPEN (7, file=output)
      DO i = 1-nh, n
         IF (i < 1) THEN
            READ (7, *)
         ELSE
            READ (7, *) k(i), P(i)
         END IF
      END DO
      CLOSE (7)

      ! Convert k to k/h
      k = k/cosm%h

      ! Done
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_FRANKENEMU_POWER_Z: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_FrankenEmu_power_z

   SUBROUTINE get_MiraTitan_power_z(k, P, n, z, cosm, verbose)

      USE constants
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      CHARACTER(len=256) :: output
      INTEGER :: i

      CALL EXECUTE_COMMAND_LINE('rm xstar.dat')
      CALL EXECUTE_COMMAND_LINE('rm EMU0.txt')

      OPEN (7, file='xstar.dat')
      WRITE (7, *) (cosm%Om_m*cosm%h**2), (cosm%Om_b*cosm%h**2), cosm%sig8, &
         cosm%h, cosm%ns, cosm%w, cosm%wa, (cosm%om_nu*cosm%h**2), z
      CLOSE (7)

      CALL EXECUTE_COMMAND_LINE('/Users/Mead/Physics/MiraTitan/P_tot/emu.exe')
      output = 'EMU0.txt'

      n = file_length(output, verbose=.FALSE.)
      IF (ALLOCATED(k)) DEALLOCATE (k)
      IF (ALLOCATED(P)) DEALLOCATE (P)
      ALLOCATE (k(n), P(n))

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_MIRATITAN_POWER_Z: z:', z
         WRITE (*, *) 'GET_MIRATITAN_POWER_Z: P(k) file length:', n
         WRITE (*, *) 'GET_MIRATITAN_POWER_Z: Reading in P(k)'
      END IF

      OPEN (7, file=output)
      DO i = 1, n
         READ (7, *) k(i), P(i)
      END DO
      CLOSE (7)

      ! Convert P(k) to Delta^2(k)
      P = P*(k**3)*4.*pi/(2.*pi)**3

      ! Convert k to k/h
      k = k/cosm%h

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_MIRATITAN_POWER_Z: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_MiraTitan_power_z

END MODULE cosmic_emu_stuff
