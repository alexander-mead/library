MODULE cosmic_emu_stuff

   USE file_info
   USE array_operations
   USE string_operations
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
   PUBLIC :: calculate_Emulator_power

   PUBLIC :: emulator_CosmicEmu
   PUBLIC :: emulator_FrankenEmu
   PUBLIC :: emulator_MiraTitan
   PUBLIC :: emulator_Euclid

   INTEGER, PARAMETER :: emulator_CosmicEmu = 1
   INTEGER, PARAMETER :: emulator_FrankenEmu = 2
   INTEGER, PARAMETER :: emulator_MiraTitan = 3
   INTEGER, PARAMETER :: emulator_Euclid = 4

   ! Cosmic Emu
   REAL, PARAMETER :: kmin_rebin_CosmicEmu = 1e-2
   REAL, PARAMETER :: kmax_rebin_CosmicEmu = 1.
   !REAL, PARAMETER :: kmax_rebin_CosmicEmu = 2.7 ! This is the maximum
   !REAL, PARAMETER :: kmax_rebin_CosmicEmu = 0.2 ! Perturbation-theory studies
   INTEGER, PARAMETER :: nk_rebin_CosmicEmu = 128
   CHARACTER(len=256), PARAMETER :: params_CosmicEmu = 'emu_params.txt'
   CHARACTER(len=256), PARAMETER :: output_CosmicEmu = 'emu_power.dat'
   CHARACTER(len=256), PARAMETER :: exe_CosmicEmu = '/Users/Mead/Physics/CosmicEmu/emu.exe'
   INTEGER, PARAMETER :: nh_CosmicEmu = 10   ! Length of header
   INTEGER, PARAMETER :: hl_CosmicEmu = 8    ! Line that h in on'
   REAL, PARAMETER :: eps_h_CosmicEmu = 3e-3 ! Tolerance for how close 'h' needs to be to that derived from r_sound

   ! Franken Emu
   REAL, PARAMETER :: kmin_rebin_FrankenEmu = 1e-2
   REAL, PARAMETER :: kmax_rebin_FrankenEmu = 10.
   INTEGER, PARAMETER :: nk_rebin_FrankenEmu = 128
   CHARACTER(len=256), PARAMETER :: params_FrankenEmu = 'emu_params.txt'
   CHARACTER(len=256), PARAMETER :: output_FrankenEmu = 'emu_power.dat'
   CHARACTER(len=256), PARAMETER :: exe_FrankenEmu = '/Users/Mead/Physics/FrankenEmu/emu.exe'

   ! Mira Titan
   REAL, PARAMETER :: kmin_rebin_MiraTitan = 1e-2
   REAL, PARAMETER :: kmax_rebin_MiraTitan = 7.
   INTEGER, PARAMETER :: nk_rebin_MiraTitan = 128
   CHARACTER(len=256), PARAMETER :: params_MiraTitan = 'xstar.dat'
   CHARACTER(len=256), PARAMETER :: output_MiraTitan = 'EMU0.txt'
   CHARACTER(len=256), PARAMETER :: exe_MiraTitan = '/Users/Mead/Physics/MiraTitan/P_tot/emu.exe'

   ! Euclid Emulator
   CHARACTER(len=256), PARAMETER :: exe_Euclid = '/Users/Mead/Physics/EuclidEmulator2/ee2.exe'
   CHARACTER(len=256), PARAMETER :: dir_Euclid = '/Users/Mead/Physics/EuclidEmulator2/'
   CHARACTER(len=256), PARAMETER :: outdir_Euclid = trim(dir_Euclid)//'Mead/'
   CHARACTER(len=256), PARAMETER :: outbase_Euclid = 'temp'
   CHARACTER(len=256), PARAMETER :: outfile_Euclid = trim(outbase_Euclid)//'0.dat'

   ! Rebinnig
   LOGICAL, PARAMETER :: logk_interp_emulator = .TRUE.
   LOGICAL, PARAMETER :: logPk_interp_emulator = .TRUE.
   LOGICAL, PARAMETER :: logBk_interp_emulator = .FALSE.
   INTEGER, PARAMETER :: iorder_rebin = 3
   INTEGER, PARAMETER :: ifind_rebin = ifind_split
   INTEGER, PARAMETER :: iinterp_rebin = iinterp_Lagrange

CONTAINS

   SUBROUTINE calculate_Emulator_power(k, a, Pk, cosm, version, verbose)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, INTENT(IN) :: version
      LOGICAL, INTENT(IN) :: verbose

      IF (is_in_array(version, [emulator_CosmicEmu, emulator_FrankenEmu, emulator_MiraTitan])) THEN
         CALL calculate_CosmicEmulator_power(k, a, Pk, cosm, version, verbose)
      ELSE IF (version == emulator_Euclid) THEN
         CALL calculate_EuclidEmulator_power(k, a, Pk, cosm)
      ELSE
         STOP 'CALCULATE_EMULATOR_POWER: Error, emulator version not recognised'
      END IF

   END SUBROUTINE calculate_Emulator_power

   SUBROUTINE calculate_EuclidEmulator_power(k, a, Pk, cosm)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)    
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: k_EE(:), Bk_EE(:), Pk_lin(:, :)
      INTEGER :: ia, na, nk
      CHARACTER(len=256) :: corrfile = trim(outdir_Euclid)//trim(outfile_Euclid)
      INTEGER, PARAMETER :: iorder = iorder_rebin
      INTEGER, PARAMETER :: ifind = ifind_rebin
      INTEGER, PARAMETER :: iinterp = iinterp_rebin
      LOGICAL, PARAMETER :: logk_interp = logk_interp_emulator
      LOGICAL, PARAMETER :: logBk_interp = logBk_interp_emulator

      ! Allocate array for the output power spectrumm
      nk = size(k)
      na = size(a)
      ALLOCATE(Pk(nk, na))

      ! Calculate the linear power
      ! NOTE: Do this before calling the emulator to make sure that As is set
      CALL calculate_plin(k, a, Pk_lin, cosm)

      ! Run the emulator and interpolate the correction from the EE k to the input k
      DO ia = 1, na
         CALL run_EuclidEumulator(a(ia), cosm)
         CALL read_EuclidEmulator_correction(k_EE, Bk_EE, corrfile)
         CALL interpolate_array(k_EE, Bk_EE, k, Pk(:, ia), iorder, ifind, iinterp, logk_interp, logBk_interp)
      END DO

      ! Multiply the correction by the linear power to get the non-linear spectrum
      Pk = Pk_lin*Pk

   END SUBROUTINE calculate_EuclidEmulator_power

   SUBROUTINE run_EuclidEumulator(a, cosm)

      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      CHARACTER(len=256) :: command, b, m, s, n, h, w, wa, As, z, dir, out
      CHARACTER(len=256), PARAMETER :: exe = exe_Euclid
      CHARACTER(len=256), PARAMETER :: basedir = dir_Euclid
      CHARACTER(len=256), PARAMETER :: outdir = outdir_Euclid
      CHARACTER(len=256), PARAMETER :: outbase = outbase_Euclid
      CHARACTER(len=256), PARAMETER :: outfile = outfile_Euclid

      ! Checks
      IF (cosm%k /= 0.) STOP 'RUN_EUCLIDEMULATOR: Error, does not support curved models'

      ! Remove previous data files
      CALL EXECUTE_COMMAND_LINE('rm -f '//trim(outdir)//trim(outfile))

      ! Commands
      b = ' -b '//trim(real_to_string(cosm%Om_b, 1, 5))
      m = ' -m '//trim(real_to_string(cosm%Om_m, 1, 5))
      s = ' -s '//trim(real_to_string(cosm%m_nu, 1, 5))
      n = ' -n '//trim(real_to_string(cosm%ns, 1, 5))
      h = ' -H '//trim(real_to_string(cosm%h, 1, 5))
      w = ' -W '//trim(real_to_string(cosm%w, 1, 5))
      wa = ' -w '//trim(real_to_string(cosm%wa, 1, 5))
      As = ' -A '//trim(exp_to_string(cosm%As, 1, 5, -9))
      z = ' -z '//trim(real_to_string(redshift_a(a), 2, 5))
      dir = ' -d '//trim(outdir)
      out = ' -o '//trim(outbase)
      command = trim(exe)//trim(b)//trim(m)//trim(s)//trim(n)//trim(h)//trim(w)//trim(wa)//trim(As)//trim(z)//trim(dir)//trim(out)

      ! Run the emulator
      CALL EXECUTE_COMMAND_LINE('cd '//trim(basedir)//' && '//trim(command)//' > /dev/null')

   END SUBROUTINE run_EuclidEumulator

   SUBROUTINE read_EuclidEmulator_correction(k, Bk, infile)

      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Bk(:)
      CHARACTER(len=*), INTENT(IN) :: infile
      INTEGER :: ik, nk, u

      ! Allocate arrays (remember one-line header)
      nk = file_length(infile)-1
      ALLOCATE(k(nk), Bk(nk))

      ! Read in non-linear correction data
      OPEN(newunit=u, file=infile)
      READ(u, *) ! Bypass one-line header
      DO ik = 1, nk
         READ(u, *) k(ik), Bk(ik)
      END DO
      CLOSE(u)

   END SUBROUTINE read_EuclidEmulator_correction

   SUBROUTINE calculate_CosmicEmulator_power(k, a, Pk, cosm, version, verbose)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)    
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, INTENT(IN) :: version
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: k_emu(:), Pk_emu(:, :)
      INTEGER :: ia, nk, na
      INTEGER, PARAMETER :: iorder = iorder_rebin
      INTEGER, PARAMETER :: ifind = ifind_rebin
      INTEGER, PARAMETER :: iinterp = iinterp_rebin
      LOGICAL, PARAMETER :: logk_interp = logk_interp_emulator
      LOGICAL, PARAMETER :: logPk_interp = logPk_interp_emulator

      nk = size(k)
      na = size(a)
      ALLOCATE(Pk(nk, na))

      CALL get_emulator_power(k_emu, a, Pk_emu, nk, cosm, version, rebin=.FALSE., verbose=verbose)

      DO ia = 1, na
         CALL interpolate_array(k_emu, Pk_emu(:, ia), k, Pk(:, ia), iorder, ifind, iinterp, logk_interp, logPk_interp)
      END DO

   END SUBROUTINE calculate_CosmicEmulator_power

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
      INTEGER, PARAMETER :: nh = 5 ! Length of header
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
