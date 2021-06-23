MODULE cosmic_emu_stuff

   ! TODO: Could remove depdence on coosmology_functions and have internal cosmology type for emulator

   USE file_info
   USE basic_operations
   USE array_operations
   USE string_operations
   USE cosmology_functions
   USE interpolate
   USE table_integer

   IMPLICIT NONE

   PRIVATE

   ! Public functions
   PUBLIC :: get_CosmicEmu_power_z
   PUBLIC :: get_FrankenEmu_power_z
   PUBLIC :: get_MiraTitan_power_z
   PUBLIC :: get_emulator_power
   PUBLIC :: calculate_Emulator_power

   ! Emulator versions
   PUBLIC :: emulator_CosmicEmu
   PUBLIC :: emulator_FrankenEmu
   PUBLIC :: emulator_MiraTitan
   PUBLIC :: emulator_Euclid
   PUBLIC :: emulator_BACCO
   PUBLIC :: emulator_NGenHALOFIT

   ! Emulator versions
   INTEGER, PARAMETER :: emulator_CosmicEmu = 1
   INTEGER, PARAMETER :: emulator_FrankenEmu = 2
   INTEGER, PARAMETER :: emulator_MiraTitan = 3
   INTEGER, PARAMETER :: emulator_Euclid = 4
   INTEGER, PARAMETER :: emulator_BACCO = 5
   INTEGER, PARAMETER :: emulator_NGenHALOFIT = 6

   ! Cosmic Emu
   REAL, PARAMETER :: kmin_rebin_CosmicEmu = 1e-2
   REAL, PARAMETER :: kmax_rebin_CosmicEmu = 1.
   !REAL, PARAMETER :: kmax_rebin_CosmicEmu = 2.7 ! This is the maximum
   !REAL, PARAMETER :: kmax_rebin_CosmicEmu = 0.2 ! Perturbation-theory studies
   INTEGER, PARAMETER :: nk_rebin_CosmicEmu = 129
   CHARACTER(len=256), PARAMETER :: params_CosmicEmu = 'emu_params.txt'
   CHARACTER(len=256), PARAMETER :: output_CosmicEmu = 'emu_power.dat'
   CHARACTER(len=256), PARAMETER :: exe_CosmicEmu = '/Users/Mead/Physics/CosmicEmu/emu.exe'
   INTEGER, PARAMETER :: nh_CosmicEmu = 10   ! Length of header
   INTEGER, PARAMETER :: hl_CosmicEmu = 8    ! Line that h in on'
   REAL, PARAMETER :: eps_h_CosmicEmu = 3e-3 ! Tolerance for how close 'h' needs to be to that derived from r_sound

   ! Franken Emu
   REAL, PARAMETER :: kmin_rebin_FrankenEmu = 1e-2
   REAL, PARAMETER :: kmax_rebin_FrankenEmu = 10.
   INTEGER, PARAMETER :: nk_rebin_FrankenEmu = 129
   CHARACTER(len=256), PARAMETER :: params_FrankenEmu = 'emu_params.txt'
   CHARACTER(len=256), PARAMETER :: output_FrankenEmu = 'emu_power.dat'
   CHARACTER(len=256), PARAMETER :: exe_FrankenEmu = '/Users/Mead/Physics/FrankenEmu/emu.exe'

   ! Mira Titan
   REAL, PARAMETER :: kmin_rebin_MiraTitan = 1e-2
   REAL, PARAMETER :: kmax_rebin_MiraTitan = 7.
   INTEGER, PARAMETER :: nk_rebin_MiraTitan = 129
   CHARACTER(len=256), PARAMETER :: params_MiraTitan = 'xstar.dat'
   CHARACTER(len=256), PARAMETER :: output_MiraTitan = 'EMU0.txt'
   CHARACTER(len=256), PARAMETER :: exe_MiraTitan = '/Users/Mead/Physics/MiraTitan/P_tot/emu.exe'

   ! Euclid Emulator
   CHARACTER(len=256), PARAMETER :: exe_Euclid = '/Users/Mead/Physics/EuclidEmulator2/ee2.exe'
   CHARACTER(len=256), PARAMETER :: dir_Euclid = '/Users/Mead/Physics/EuclidEmulator2/'
   CHARACTER(len=256), PARAMETER :: outdir_Euclid = trim(dir_Euclid)//'Mead/'
   CHARACTER(len=256), PARAMETER :: outbase_Euclid = 'temp'
   CHARACTER(len=256), PARAMETER :: outfile_Euclid = trim(outbase_Euclid)//'0.dat'

   ! BACCO
   CHARACTER(len=256), PARAMETER :: exe_BACCO = '/Users/Mead/Physics/BACCO/run_BACCO.py'
   CHARACTER(len=256), PARAMETER :: outfile_BACCO = '/Users/Mead/Physics/BACCO/results.dat'

    ! NGenHALOFIT
   CHARACTER(len=256), PARAMETER :: dir_NGenHALOFIT = '/Users/Mead/Physics/NGenHalofit/'
   CHARACTER(len=256), PARAMETER :: exe_NGenHALOFIT = trim(dir_NGenHALOFIT)//'NGenHalofit.exe' 
   CHARACTER(len=256), PARAMETER :: dir_temp_NGenHALOFIT = trim(dir_NGenHALOFIT)//'Mead/'
   CHARACTER(len=256), PARAMETER :: linfile_NGenHALOFIT = 'linear_power.dat'
   CHARACTER(len=256), PARAMETER :: linfilefull_NGenHALOFIT = trim(dir_temp_NGenHALOFIT)//trim(linfile_NGenHALOFIT)
   CHARACTER(len=256), PARAMETER :: inifile_NGenHALOFIT = trim(dir_temp_NGenHALOFIT)//'input_file.dat'
   CHARACTER(len=256), PARAMETER :: expfile_NGenHALOFIT = trim(dir_temp_NGenHALOFIT)//'expansion_file.dat'
   CHARACTER(len=256), PARAMETER :: powbase_NGenHALOFIT = 'HALOFIT_power'
   CHARACTER(len=256), PARAMETER :: MPTbase_NGenHALOFIT = 'MPT_power'
   CHARACTER(len=256), PARAMETER :: vardir_NGenHALOFIT = trim(dir_NGenHALOFIT)//'DATAPlinCAMB/'
   CHARACTER(len=256), PARAMETER :: varbase_NGenHALOFIT = 'Planck2013.Step_ByHand.HighAcc_matterpower'
   REAL, PARAMETER :: kmin_lin_NGenHALOFIT = 1e-3
   REAL, PARAMETER :: kmax_lin_NGenHALOFIT = 1e2
   INTEGER, PARAMETER :: nk_plin_NGenHALOFIT = 129
   REAL, PARAMETER :: As_NGenHALOFIT = 2.14485e-9
   REAL, PARAMETER :: kpiv_noh_default_NGenHALOFIT = 0.05
   REAL, PARAMETER :: eps_kpiv_NGenHALOFIT = 1e-4
   REAL, PARAMETER :: alin_NGenHALOFIT = 1.
   INTEGER, PARAMETER :: flag_NGenHALOFIT = flag_matter
   INTEGER, PARAMETER :: iorder_interp_NGenHALOFIT = 3
   INTEGER, PARAMETER :: ifind_interp_NGenHALOFIT = ifind_split
   INTEGER, PARAMETER :: iinterp_interp_NGenHALOFIT = iinterp_Lagrange
   REAL, PARAMETER :: w_min_NGenHALOFIT = -1.05
   REAL, PARAMETER :: w_max_NGenHALOFIT = -0.95
   REAL, PARAMETER :: wa_min_NGenHALOFIT = -0.4
   REAL, PARAMETER :: wa_max_NGenHALOFIT = 0.4
   REAL, PARAMETER :: Om_m_min_NGenHALOFIT = 0.21
   REAL, PARAMETER :: Om_m_max_NGenHALOFIT = 0.4
   REAL, PARAMETER :: wc_min_NGenHALOFIT = 0.1
   REAL, PARAMETER :: wc_max_NGenHALOFIT = 0.13
   REAL, PARAMETER :: wb_min_NGenHALOFIT = 0.02
   REAL, PARAMETER :: wb_max_NGenHALOFIT = 0.04
   REAL, PARAMETER :: ns_min_NGenHALOFIT = 0.85
   REAL, PARAMETER :: ns_max_NGenHALOFIT = 1.05
   REAL, PARAMETER :: As_min_NGenHALOFIT = 1.72e-9
   REAL, PARAMETER :: As_max_NGenHALOFIT = 2.58e-9
   REAL, PARAMETER :: nrun_min_NGenHALOFIT = -0.2
   REAL, PARAMETER :: nrun_max_NGenHALOFIT = 0.2

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
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose

      IF (is_in_array(version, [emulator_CosmicEmu, emulator_FrankenEmu, emulator_MiraTitan])) THEN
         CALL calculate_CosmicEmulator_power(k, a, Pk, cosm, version, verbose)
      ELSE IF (version == emulator_Euclid) THEN
         CALL calculate_EuclidEmulator_power(k, a, Pk, cosm)
      ELSE IF (version == emulator_BACCO) THEN
         CALL calculate_BACCO_power(k, a, Pk, cosm)
      ELSE IF (version == emulator_NGenHALOFIT) THEN
         CALL calculate_NGenHALOFIT(k, a, Pk, cosm)
      ELSE
         STOP 'CALCULATE_EMULATOR_POWER: Error, emulator version not recognised'
      END IF

   END SUBROUTINE calculate_Emulator_power

   SUBROUTINE calculate_NGenHALOFIT(k, a, Pk, cosm)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: nk, na, ia
      REAL :: kmin, kmax
      REAL, ALLOCATABLE :: k_HF(:), Pk_HF(:, :)
      INTEGER, PARAMETER :: iorder = iorder_interp_NGenHALOFIT
      INTEGER, PARAMETER :: ifind = ifind_interp_NGenHALOFIT
      INTEGER, PARAMETER :: iinterp = iinterp_interp_NGenHALOFIT

      ! Array sizes
      nk = size(k)
      na = size(a)

      ! Run N-Gen HALOFIT
      ! TODO: Should nk here be the same as in the desired k array?
      kmin = k(1)
      kmax = k(nk)
      CALL run_NgenHALOFIT(kmin, kmax, nk, k_HF, a, Pk_HF, cosm)

      ! Interpolate results on to my k array
      ALLOCATE(Pk(nk, na))
      DO ia = 1, na
         CALL interpolate_array(k_HF, Pk_HF(:, ia), k, Pk(:, ia), iorder, ifind, iinterp, logx=.TRUE., logy=.TRUE.)
      END DO

   END SUBROUTINE calculate_NGenHALOFIT

   SUBROUTINE run_NGenHALOFIT(kmin, kmax, nk, k, a, Pk, cosm)

      USE camb_stuff
      REAL, INTENT(IN) :: kmin
      REAL, INTENT(IN) :: kmax
      INTEGER, INTENT(IN) :: nk
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: ik
      CHARACTER(len=256), PARAMETER :: exe = exe_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: dir = dir_temp_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: linfile = linfilefull_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: inifile = inifile_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: expfile = expfile_NGenHALOFIT
      REAL, PARAMETER :: kmin_lin = kmin_lin_NGenHALOFIT
      REAL, PARAMETER :: kmax_lin = kmax_lin_NGenHALOFIT
      INTEGER, PARAMETER :: nk_lin = nk_plin_NGenHALOFIT
      REAL, PARAMETER :: alin = alin_NGenHALOFIT
      INTEGER, PARAMETER :: flag = flag_NGenHALOFIT
      REAL, PARAMETER :: kpiv_noh_default = kpiv_noh_default_NGenHALOFIT
      REAL, PARAMETER :: eps_kpiv = eps_kpiv_NGenHALOFIT
      REAL, PARAMETER :: w_min = w_min_NGenHALOFIT
      REAL, PARAMETER :: w_max = w_max_NGenHALOFIT
      REAL, PARAMETER :: wa_min = wa_min_NGenHALOFIT
      REAL, PARAMETER :: wa_max = wa_max_NGenHALOFIT
      REAL, PARAMETER :: Om_m_min = Om_m_min_NGenHALOFIT
      REAL, PARAMETER :: Om_m_max = Om_m_max_NGenHALOFIT
      REAL, PARAMETER :: wc_min = wc_min_NGenHALOFIT
      REAL, PARAMETER :: wc_max = wc_max_NGenHALOFIT
      REAL, PARAMETER :: wb_min = wb_min_NGenHALOFIT
      REAL, PARAMETER :: wb_max = wb_max_NGenHALOFIT
      REAL, PARAMETER :: ns_min = ns_min_NGenHALOFIT
      REAL, PARAMETER :: ns_max = ns_max_NGenHALOFIT
      REAL, PARAMETER :: As_min = As_min_NGenHALOFIT
      REAL, PARAMETER :: As_max = As_max_NGenHALOFIT
      REAL, PARAMETER :: nrun_min = nrun_min_NGenHALOFIT
      REAL, PARAMETER :: nrun_max = nrun_max_NGenHALOFIT

      ! Checks
      IF (.NOT. cosm%is_init) STOP 'RUN_NGENHALOFIT: Error, cosmology is not initialised'
      IF (cosm%k /= 0.) WRITE(*, *) 'RUN_NGENHALOFIT: Warning, NGenHalofit only supports flat cosmologies'
      IF (cosm%m_nu /= 0.) WRITE(*, *) 'RUN_NGENHALOFIT: Warning, NGenHalofit does not support massive neutrino cosmologies'
      IF (.NOT. requal(cosm%kpiv*cosm%h, kpiv_noh_default, eps_kpiv)) THEN
        WRITE(*, *) 'RUN_NGENHALOFIT: kpiv [h/Mpc]:', cosm%kpiv
        WRITE(*, *) 'RUN_NGENHALOFIT: kpiv [1/Mpc]:', cosm%kpiv*cosm%h
        STOP 'RUN_NGENHALOFIT: Error, NGenHalofit assumes pivot scale of 0.05 Mpc^-1'
      END IF

      ! Needs to be called for As to be correct
      IF (.NOT. cosm%is_normalised) CALL normalise_power(cosm)

      ! Check parameter range
      IF (.NOT. between(cosm%w, w_min, w_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, w is outside boundary'
      IF (.NOT. between(cosm%wa, wa_min, wa_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, wa is outside boundary'
      IF (.NOT. between(cosm%Om_m, Om_m_min, Om_m_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, Omega_m is outside boundary'
      IF (.NOT. between(cosm%Om_c*cosm%h**2, wc_min, wc_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, omega_c is outside boundary'
      IF (.NOT. between(cosm%Om_b*cosm%h**2, wb_min, wb_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, omega_b is outside boundary'
      IF (.NOT. between(cosm%ns, ns_min, ns_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, ns is outside boundary'
      IF (.NOT. between(cosm%As, As_min, As_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, As is outside boundary'
      IF (.NOT. between(cosm%nrun, nrun_min, nrun_max)) WRITE (*, *) 'RUN_NGENHALOFIT: Warning, alpha is outside boundary'

      ! Remove any previous files
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(dir)//'*')

      ! Write the ini file and list of expansion factors
      CALL write_NGenHALOFIT_ini_file(kmin, kmax, nk, cosm, inifile, dir)
      CALL write_NGenHALOFIT_exp_file(a, expfile)

      ! Write a CAMB format linear power spectrum
      ! TODO: It is a bit lazy to re-use the k, Pk arrays here
      CALL fill_array(kmin_lin, kmax_lin, k, nk_lin, ilog=.TRUE.)
      ALLOCATE(Pk(nk_lin, 1))
      DO ik = 1, nk_lin
         Pk(ik, 1) = plin(k(ik), alin, flag, cosm)
      END DO
      CALL write_CAMB_Pk(k, Pk(:, 1), linfile, header=.FALSE.)
      DEALLOCATE(k, Pk)

      ! Run NGenHALOFIT
      CALL EXECUTE_COMMAND_LINE('cd '//trim(dir_NGenHALOFIT)//' && '//trim(exe)//' '//trim(inifile)//' '//trim(expfile)//' > /dev/null')

      ! Read in the power data
      CALL read_NGenHALOFIT_power(k, a, Pk)

   END SUBROUTINE run_NGenHALOFIT

   SUBROUTINE write_NGenHALOFIT_ini_file(kmin, kmax, nk, cosm, outfile, dir)

      REAL, INTENT(IN) :: kmin
      REAL, INTENT(IN) :: kmax
      INTEGER, INTENT(IN) :: nk
      TYPE(cosmology), INTENT(IN) :: cosm
      CHARACTER(len=256), INTENT(IN) :: outfile
      CHARACTER(len=256), INTENT(IN) :: dir
      INTEGER :: u
      LOGICAL, PARAMETER :: logk = .TRUE.
      REAL, PARAMETER :: As = As_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: powbase = powbase_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: MPTbase = MPTbase_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: linfile = linfile_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: vardir = vardir_NGenHALOFIT
      CHARACTER(len=256), PARAMETER :: varbase = varbase_NGenHALOFIT

      OPEN(newunit=u, file=outfile)

      ! Cosmology
      WRITE (u, *) 'aexp 1.0'
      WRITE (u, *) 'w0', cosm%w
      WRITE (u, *) 'w1', cosm%wa
      WRITE (u, *) 'om_ch20', cosm%Om_c*cosm%h**2
      WRITE (u, *) 'om_bh20', cosm%Om_b*cosm%h**2
      WRITE (u, *) 'om_DE0', 1.-cosm%Om_m
      WRITE (u, *) 'As', cosm%As/As
      WRITE (u, *) 'pindex', cosm%ns
      WRITE (u, *) 'running', cosm%nrun

      ! Input parameters
      WRITE (u, *) 'nPkOut', nk
      WRITE (u, *) 'rkOutMIN', kmin
      WRITE (u, *) 'rkOutMAX', kmax
      IF (logk) THEN
         WRITE (u, *) 'iLogOrLin 1'
      ELSE
         WRITE (u, *) 'iLogOrLin 0'
      END IF
      WRITE (u, *) 'iGenEffSpecTarget	1'
      WRITE (u, *) 'iGenEffSpecVar	0'

      ! File handlers
      WRITE (u, *) 'OutputDir ', trim(dir)
      WRITE (u, *) 'OutputFileBase ', trim(powbase)
      WRITE (u, *) 'OutputMPTFileBase ', trim(MPTbase)
      WRITE (u, *) 'PowDirTarget ', trim(dir)
      WRITE (u, *) 'PowFileTarget ', trim(linfile)
      WRITE (u, *) 'PowDirVar ', trim(vardir)
      WRITE (u, *) 'PowFileBaseVar ', trim(varbase)

      CLOSE(u)

   END SUBROUTINE write_NGenHALOFIT_ini_file

   SUBROUTINE write_NGenHALOFIT_exp_file(a, outfile)

      REAL, INTENT(IN) :: a(:)
      CHARACTER(len=256), INTENT(IN) :: outfile
      INTEGER :: u
      INTEGER :: ia

      ! Write a file of scale factors with one scale factor per line
      OPEN(newunit=u, file=outfile)
      DO ia = 1, size(a)
         WRITE (u, *) a(ia)
      END DO
      CLOSE(u)

   END SUBROUTINE write_NGenHALOFIT_exp_file

   SUBROUTINE read_NGenHALOFIT_power(k, a, Pk)

      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      REAL :: crap
      INTEGER :: ik, ia, nk, na
      INTEGER :: u
      CHARACTER(len=256) :: infile
      CHARACTER(len=256), PARAMETER :: inbase = trim(dir_temp_NGenHALOFIT)//trim(powbase_NGenHALOFIT)

      ! Loop over scale factors
      na = size(a)
      DO ia = 1, na

         ! Get the infile name and allocate arrays if necessary
         infile = trim(inbase)//'.'//trim(integer_to_string(ia-1))//'.dat'
         IF (ia == 1) THEN
            nk = file_length(infile)
            ALLOCATE(k(nk), Pk(nk, na))
         END IF

         ! Read in NGenHalofit power spectrum
         OPEN(newunit=u, file=infile)
         DO ik = 1, size(k)
            !READ(u, *) k(ik), crap, Pk(ik, ia) ! Column 3 is Takahashi HALOFIT
            READ(u, *) k(ik), crap, crap, Pk(ik, ia) ! Column 4 is NGenHALOFIT
         END DO
         CLOSE(u)

         ! Convert P(k) to Delta^2(k)
         Pk(:, ia) = Delta_Pk(Pk(:, ia), k)

      END DO

   END SUBROUTINE read_NGenHALOFIT_power

   SUBROUTINE calculate_BACCO_power(k, a, Pk, cosm)

      ! NOTE: BACCO calculates things for the cold matter power spectrum
      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)    
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: k_emu(:), Bk_emu(:), Pk_mm_lin(:, :), Pk_cc_lin(:, :)
      INTEGER :: ik, ia, nk, na
      CHARACTER(len=256) :: corrfile = trim(outfile_BACCO)
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
      ! NOTE: This calculates the linear matter-matter power, but BACCO only considers cold matter
      CALL calculate_plin(k, a, Pk_mm_lin, flag_matter, cosm)
      CALL calculate_plin(k, a, Pk_cc_lin, flag_cold, cosm)

      ! Run the emulator and interpolate the correction from the EE k to the input k
      DO ia = 1, na
         CALL run_BACCO(a(ia), cosm)
         CALL read_BACCO_correction(k_emu, Bk_emu, corrfile)
         CALL interpolate_array(k_emu, Bk_emu, k, Pk(:, ia), iorder, ifind, iinterp, logk_interp, logBk_interp)
         DO ik = 1, nk
            IF (k(ik) < k_emu(1)) Pk(ik, ia) = 1.
         END DO
      END DO

      ! Multiply the correction by the linear power to get the non-linear spectrum
      ! NOTE: Apply correction to cold spectrum only, then add residual from neutrinos
      ! NOTE: P_mm = P_cc + 2P_nc + P_nn; for both non-linear and linear
      ! NOTE: P_mm - P_cc = 2P_nc + P_nn; so making assumption that P_nc and P_nn are both linear
      Pk = Pk_cc_lin*Pk+Pk_mm_lin-Pk_cc_lin

   END SUBROUTINE calculate_BACCO_power

   SUBROUTINE run_BACCO(a, cosm)

      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      CHARACTER(len=256) :: command, c, s8, b, ns, h, nu, w, wa, aexp
      CHARACTER(len=256), PARAMETER :: exe = exe_BACCO
      CHARACTER(len=256), PARAMETER :: outfile = outfile_BACCO

      ! Checks
      IF (cosm%k /= 0.) STOP 'RUN_BACCO: Error, does not support curved models'

      ! Remove previous data files
      CALL EXECUTE_COMMAND_LINE('rm -f '//trim(outfile))

      ! Commands
      c = ' '//trim(real_to_string(cosm%Om_c+cosm%Om_b, 1, 5))
      s8 = ' '//trim(real_to_string(cosm%sig8*(1.-cosm%f_nu), 1, 5))
      b = ' '//trim(real_to_string(cosm%Om_b, 1, 5))
      ns = ' '//trim(real_to_string(cosm%ns, 1, 5))
      h = ' '//trim(real_to_string(cosm%h, 1, 5))
      nu = ' '//trim(real_to_string(cosm%m_nu, 1, 5))
      w = ' '//trim(real_to_string(cosm%w, 1, 5))
      wa = ' '//trim(real_to_string(cosm%wa, 1, 5))
      aexp = ' '//trim(real_to_string(a, 1, 5))
      command = 'python3 '//trim(exe)//trim(c)//trim(s8)//trim(b)//trim(ns)//trim(h)//trim(nu)//trim(w)//trim(wa)//trim(aexp)

      ! Run the emulator
      CALL EXECUTE_COMMAND_LINE(trim(command)//' > /dev/null')

   END SUBROUTINE run_BACCO

   SUBROUTINE read_BACCO_correction(k, Bk, infile)

      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Bk(:)
      CHARACTER(len=*), INTENT(IN) :: infile
      INTEGER :: ik, nk, u

      ! Allocate arrays (remember one-line header)
      nk = file_length(infile)
      ALLOCATE(k(nk), Bk(nk))

      ! Read in non-linear correction data
      OPEN(newunit=u, file=infile)
      DO ik = 1, nk
         READ(u, *) k(ik), Bk(ik)
      END DO
      CLOSE(u)

   END SUBROUTINE read_BACCO_correction

   SUBROUTINE calculate_EuclidEmulator_power(k, a, Pk, cosm)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)    
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: k_emu(:), Bk_emu(:), Pk_lin(:, :)
      INTEGER :: ik, ia, nk, na
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
      CALL calculate_plin(k, a, Pk_lin, flag_matter, cosm)

      ! Run the emulator and interpolate the correction from the EE k to the input k
      DO ia = 1, na
         CALL run_EuclidEumulator(a(ia), cosm)
         CALL read_EuclidEmulator_correction(k_emu, Bk_emu, corrfile)
         CALL interpolate_array(k_emu, Bk_emu, k, Pk(:, ia), iorder, ifind, iinterp, logk_interp, logBk_interp)
         DO ik = 1, nk
            IF (k(ik) < k_emu(1)) Pk(ik, ia) = 1.
         END DO
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
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
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
