MODULE cosmic_emu_stuff

   USE file_info
   USE array_operations
   USE cosmology_functions
   USE interpolate
   USE basic_operations
   USE table_integer

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

   REAL, PARAMETER :: eps_h_CosmicEmu = 3e-3

   REAL, PARAMETER :: kmin_rebin_CosmicEmu = 1e-2
   REAL, PARAMETER :: kmax_rebin_CosmicEmu = 1.
   !REAL, PARAMETER :: kmax_rebin_CosmicEmu = 2.7 ! This is the maximum
   !REAL, PARAMETER :: kmax_rebin_CosmicEmu = 0.2 ! Perturbation-theory studies
   INTEGER, PARAMETER :: nk_rebin_CosmicEmu = 128
   CHARACTER(len=256), PARAMETER :: params_CosmicEmu = 'emu_params.txt'
   CHARACTER(len=256), PARAMETER :: output_CosmicEmu = 'emu_power.dat'
   CHARACTER(len=256), PARAMETER :: exe_CosmicEmu = '/Users/Mead/Physics/CosmicEmu/emu.exe'
   INTEGER, PARAMETER :: nh_CosmicEmu = 10  ! Length of header
   INTEGER, PARAMETER :: hl_CosmicEmu = 8   ! Line that h in on

   REAL, PARAMETER :: kmin_rebin_FrankenEmu = 1e-2
   REAL, PARAMETER :: kmax_rebin_FrankenEmu = 10.
   INTEGER, PARAMETER :: nk_rebin_FrankenEmu = 128
   CHARACTER(len=256), PARAMETER :: params_FrankenEmu = 'emu_params.txt'
   CHARACTER(len=256), PARAMETER :: output_FrankenEmu = 'emu_power.dat'
   CHARACTER(len=256), PARAMETER :: exe_FrankenEmu = '/Users/Mead/Physics/FrankenEmu/emu.exe'

   REAL, PARAMETER :: kmin_rebin_MiraTitan = 1e-2
   REAL, PARAMETER :: kmax_rebin_MiraTitan = 7.
   INTEGER, PARAMETER :: nk_rebin_MiraTitan = 128
   CHARACTER(len=256), PARAMETER :: params_MiraTitan = 'xstar.dat'
   CHARACTER(len=256), PARAMETER :: output_MiraTitan = 'EMU0.txt'
   CHARACTER(len=256), PARAMETER :: exe_MiraTitan = '/Users/Mead/Physics/MiraTitan/P_tot/emu.exe'

   INTEGER, PARAMETER :: iorder_rebin = 3
   INTEGER, PARAMETER :: ifind_rebin = ifind_split
   INTEGER, PARAMETER :: iinterp_rebin = iinterp_Lagrange

CONTAINS

   SUBROUTINE get_emulator_power(k, a, Pk, nk, cosm, version, rebin, verbose)

      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      INTEGER, INTENT(OUT) :: nk
      TYPE(cosmology), INTENT(IN) :: cosm   
      INTEGER, INTENT(IN) :: version
      LOGICAL, OPTIONAL, INTENT(IN) :: rebin
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: ia, na
      REAL :: z
      REAL, ALLOCATABLE :: k_emu(:), Pk_emu(:)

      na = size(a)

      DO ia = 1, na

         z = redshift_a(a(ia))

         IF (version == emulator_CosmicEmu) THEN
            CALL get_CosmicEmu_power_z(k_emu, Pk_emu, nk, z, cosm, rebin, verbose)
         ELSE IF (version == emulator_FrankenEmu) THEN
            CALL get_FrankenEmu_power_z(k_emu, Pk_emu, nk, z, cosm, rebin, verbose)
         ELSE IF (version == emulator_MiraTitan) THEN
            CALL get_MiraTitan_power_z(k_emu, Pk_emu, nk, z, cosm, rebin, verbose)
         ELSE
            STOP 'GET_EMULATOR_POWER: Error, version not specified correctly'
         END IF

         IF (.NOT. allocated(k))  ALLOCATE(k(nk))
         IF (.NOT. allocated(Pk)) ALLOCATE(Pk(nk, na))

         k = k_emu
         Pk(:, ia) = Pk_emu

      END DO

   END SUBROUTINE get_emulator_power

   SUBROUTINE rebin_emulator(kmin, kmax, nk, k, Pk)

      REAL, INTENT(IN) :: kmin
      REAL, INTENT(IN) :: kmax
      INTEGER, INTENT(IN) :: nk
      REAL, ALLOCATABLE, INTENT(INOUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Pk(:)
      INTEGER, PARAMETER :: iorder = iorder_rebin
      INTEGER, PARAMETER :: ifind = ifind_rebin
      INTEGER, PARAMETER :: iinterp = iinterp_rebin

      CALL rebin_array(kmin, kmax, nk, k, Pk, &
         iorder_rebin, &
         ifind_rebin, &
         iinterp_rebin, &
         logx=.TRUE., &
         logf=.TRUE. &
         )

   END SUBROUTINE rebin_emulator

   SUBROUTINE get_CosmicEmu_power_z(k, P, n, z, cosm, rebin, verbose)
    
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, OPTIONAL, INTENT(IN) :: rebin
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i
      CHARACTER(len=256) :: crap
      REAL :: h
      INTEGER, PARAMETER :: nh = nh_CosmicEmu ! Length of header
      INTEGER, PARAMETER :: hl = hl_CosmicEmu ! Line that h in on
      CHARACTER(len=256), PARAMETER :: params = params_CosmicEmu
      CHARACTER(len=256), PARAMETER :: output = output_CosmicEmu
      CHARACTER(len=256), PARAMETER :: exe = exe_CosmicEmu
      REAL, PARAMETER :: eps_h = eps_h_CosmicEmu
      REAL, PARAMETER :: kmin_rebin = kmin_rebin_CosmicEmu
      REAL, PARAMETER :: kmax_rebin = kmax_rebin_CosmicEmu   
      INTEGER, PARAMETER :: nk_rebin = nk_rebin_CosmicEmu

      ! Remove previous parameter and power file
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(params))
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(output))

      ! Write a new parameter file
      OPEN (7, file=params)
      WRITE (7, fmt='(A20,7F10.5)') trim(output), cosm%om_m*cosm%h**2, cosm%om_b*cosm%h**2, cosm%ns, cosm%sig8, cosm%w, z
      CLOSE (7)

      IF (present_and_correct(verbose)) THEN
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params))
      ELSE
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params)//' > /dev/null')
      END IF

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
         IF (i == hl-nh) THEN
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

      IF (present_and_correct(rebin)) THEN
         IF (present_and_correct(verbose)) THEN
            WRITE (*, *) 'GET_COSMICEMU_POWER_Z: Rebinning power'
         END IF
         CALL rebin_emulator(kmin_rebin, kmax_rebin, nk_rebin, k, P)
         n = nk_rebin
      END IF

      ! Done
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_COSMICEMU_POWER_Z: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_CosmicEmu_power_z

   SUBROUTINE get_FrankenEmu_power_z(k, P, n, z, cosm, rebin, verbose)

      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, OPTIONAL, INTENT(IN) :: rebin
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i
      INTEGER, PARAMETER :: nh = 5   ! Length of header
      CHARACTER(len=256), PARAMETER :: params = params_FrankenEmu
      CHARACTER(len=256), PARAMETER :: output = output_FrankenEmu
      CHARACTER(len=256), PARAMETER :: exe = exe_FrankenEmu
      REAL, PARAMETER :: kmin_rebin = kmin_rebin_FrankenEmu
      REAL, PARAMETER :: kmax_rebin = kmax_rebin_FrankenEmu
      INTEGER, PARAMETER :: nk_rebin = nk_rebin_FrankenEmu

      ! Remove previous parameter and power file
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(params))
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(output))

      ! Write a new parameter file
      OPEN (7, file=params)
      WRITE (7, fmt='(A20,7F10.5)') trim(output), cosm%om_m*cosm%h**2, cosm%om_b*cosm%h**2, cosm%ns, -cosm%w, cosm%sig8, cosm%h, z
      CLOSE (7)

      IF (present_and_correct(verbose)) THEN
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params))
      ELSE
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params)//' > /dev/null')
      END IF

      ! Get length of emu file
      n = file_length(output)
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

      IF (present_and_correct(rebin)) THEN
         IF (present_and_correct(verbose)) THEN
            WRITE (*, *) 'GET_COSMICEMU_POWER_Z: Rebinning power'
         END IF
         CALL rebin_emulator(kmin_rebin, kmax_rebin, nk_rebin, k, P)
         n = nk_rebin
      END IF

      ! Done
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_FRANKENEMU_POWER_Z: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_FrankenEmu_power_z

   SUBROUTINE get_MiraTitan_power_z(k, P, n, z, cosm, rebin, verbose)

      USE constants
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, OPTIONAL, INTENT(IN) :: rebin
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i
      REAL, PARAMETER :: kmin_rebin = kmin_rebin_MiraTitan
      REAL, PARAMETER :: kmax_rebin = kmax_rebin_MiraTitan
      INTEGER, PARAMETER :: nk_rebin = nk_rebin_MiraTitan
      CHARACTER(len=256), PARAMETER :: params = params_MiraTitan
      CHARACTER(len=256), PARAMETER :: output = output_MiraTitan
      CHARACTER(len=256), PARAMETER :: exe = exe_MiraTitan

      ! Remove previous parameter and power file
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(params))
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(output))

      OPEN (7, file=params)
      WRITE (7, *) (cosm%Om_m*cosm%h**2), (cosm%Om_b*cosm%h**2), cosm%sig8, &
         cosm%h, cosm%ns, cosm%w, cosm%wa, (cosm%om_nu*cosm%h**2), z
      CLOSE (7)

      IF (present_and_correct(verbose)) THEN
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params))
      ELSE
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params)//' > /dev/null')
      END IF

      n = file_length(output)
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
      P = Delta_Pk(P, k)

      ! Convert k to k/h
      k = k/cosm%h

      IF (present_and_correct(rebin)) THEN
         IF (present_and_correct(verbose)) THEN
            WRITE (*, *) 'GET_COSMICEMU_POWER_Z: Rebinning power'
         END IF
         CALL rebin_emulator(kmin_rebin, kmax_rebin, nk_rebin, k, P)
         n = nk_rebin
      END IF

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'GET_MIRATITAN_POWER_Z: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_MiraTitan_power_z

END MODULE cosmic_emu_stuff
