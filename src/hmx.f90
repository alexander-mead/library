MODULE HMx

   ! Module usage statements
   USE constants
   USE array_operations
   USE solve_equations
   USE special_functions
   USE string_operations
   USE calculus_table
   USE cosmology_functions
   USE basic_operations
   USE interpolate
   USE file_info
   USE table_integer

   IMPLICIT NONE

   PRIVATE

   ! Type
   PUBLIC :: halomod

   ! Main routines
   PUBLIC :: assign_halomod
   PUBLIC :: init_halomod
   PUBLIC :: print_halomod
   PUBLIC :: assign_init_halomod

   ! Write routines
   PUBLIC :: write_power
   PUBLIC :: write_power_a
   PUBLIC :: write_power_a_multiple
   PUBLIC :: write_power_fields
   PUBLIC :: read_power_a

   ! Halo-model calculations  
   PUBLIC :: calculate_HMx
   PUBLIC :: calculate_HMx_full
   PUBLIC :: calculate_HMx_old ! TODO: Replace sometime, but Tilman uses this
   PUBLIC :: calculate_HMx_a   ! TODO: Remove from public (hard task)
   PUBLIC :: calculate_HMcode
   PUBLIC :: calculate_HMcode_full
   PUBLIC :: calculate_halomod
   PUBLIC :: calculate_halomod_full
   PUBLIC :: calculate_halomod_response

   ! Halo-model things
   PUBLIC :: set_halo_type
   PUBLIC :: halo_type
   PUBLIC :: M_nu
   PUBLIC :: b_nu
   PUBLIC :: g_nu
   PUBLIC :: nu_M
   PUBLIC :: virial_radius
   PUBLIC :: convert_mass_definitions
   PUBLIC :: win
   PUBLIC :: UPP              ! TODO: Retire
   PUBLIC :: p_1void          ! TODO: Retire
   PUBLIC :: halo_HI_fraction ! TODO: Retire
   PUBLIC :: T_1h
   PUBLIC :: BNL
   PUBLIC :: simple_twohalo
   PUBLIC :: simple_onehalo

   ! Mean things
   PUBLIC :: mean_bias_number_weighted
   PUBLIC :: mean_nu_number_weighted
   PUBLIC :: mean_halo_mass_number_weighted
   PUBLIC :: mean_bias_mass_weighted
   PUBLIC :: mean_nu_mass_weighted
   PUBLIC :: mean_halo_mass_mass_weighted
   PUBLIC :: mean_halo_number_density

   ! Effective power spectrum things
   PUBLIC :: effective_index
   PUBLIC :: effective_curvature

   ! Diagnostics
   PUBLIC :: halo_definitions
   PUBLIC :: halo_diagnostics
   PUBLIC :: halo_properties
   PUBLIC :: write_halo_profiles
   PUBLIC :: write_halo_fractions

   ! Winint functions
   PUBLIC :: winint_diagnostics
   PUBLIC :: winint_speed_tests

   ! Mass functions and bias
   PUBLIC :: mass_function
   PUBLIC :: multiplicity_function
   PUBLIC :: halo_bias

   ! HMcode versions
   PUBLIC :: HMcode2015
   PUBLIC :: HMcode2015_CAMB
   PUBLIC :: HMcode2016
   PUBLIC :: HMcode2016_CAMB
   PUBLIC :: HMcode2016_neutrinofix
   PUBLIC :: HMcode2016_neutrinofix_CAMB
   PUBLIC :: HMcode2018
   PUBLIC :: HMcode2019
   PUBLIC :: HMcode2020
   PUBLIC :: HMcode2020_CAMB
   PUBLIC :: HMcode2020_feedback
   PUBLIC :: HMcode2020_feedback_CAMB

   ! HMx versions
   PUBLIC :: HMx2020_matter_w_temp_scaling
   PUBLIC :: HMx2020_matter_pressure_w_temp_scaling

   ! HMx functions
   PUBLIC :: HMx_alpha
   PUBLIC :: HMx_beta
   PUBLIC :: HMx_eps
   PUBLIC :: HMx_Gamma
   PUBLIC :: HMx_M0
   PUBLIC :: HMx_Astar
   PUBLIC :: HMx_Twhim

   ! Fields
   PUBLIC :: field_dmonly
   PUBLIC :: field_matter
   PUBLIC :: field_cdm
   PUBLIC :: field_gas
   PUBLIC :: field_stars
   PUBLIC :: field_bound_gas
   PUBLIC :: field_ejected_gas
   PUBLIC :: field_electron_pressure
   PUBLIC :: field_void
   PUBLIC :: field_compensated_void
   PUBLIC :: field_central_galaxies
   PUBLIC :: field_satellite_galaxies
   PUBLIC :: field_galaxies
   PUBLIC :: field_HI
   PUBLIC :: field_cold_bound_gas
   PUBLIC :: field_hot_bound_gas
   PUBLIC :: field_normal_bound_gas
   PUBLIC :: field_central_stars
   PUBLIC :: field_satellite_stars
   PUBLIC :: field_CIB_353
   PUBLIC :: field_CIB_545
   PUBLIC :: field_CIB_857
   PUBLIC :: field_halo_11p0_11p5
   PUBLIC :: field_halo_11p5_12p0
   PUBLIC :: field_halo_12p0_12p5
   PUBLIC :: field_halo_12p5_13p0
   PUBLIC :: field_halo_13p0_13p5
   PUBLIC :: field_halo_13p5_14p0
   PUBLIC :: field_halo_14p0_14p5
   PUBLIC :: field_halo_14p5_15p0
   PUBLIC :: field_neutrino
   PUBLIC :: field_n ! Total number of fields

   ! Total number of fitting parameters
   PUBLIC :: param_n

   ! Fitting parameters - hydro
   PUBLIC :: param_alpha
   PUBLIC :: param_beta
   PUBLIC :: param_eps
   PUBLIC :: param_eps2
   PUBLIC :: param_Gamma
   PUBLIC :: param_Zamma
   PUBLIC :: param_M0
   PUBLIC :: param_Astar
   PUBLIC :: param_Twhim
   PUBLIC :: param_cstar
   PUBLIC :: param_sstar
   PUBLIC :: param_mstar
   PUBLIC :: param_fcold
   PUBLIC :: param_fhot
   PUBLIC :: param_eta
   PUBLIC :: param_ibeta
   PUBLIC :: param_gbeta

   ! Fitting parameters - hydro mass indices
   PUBLIC :: param_alphap
   PUBLIC :: param_betap
   PUBLIC :: param_Gammap
   PUBLIC :: param_Zammap
   PUBLIC :: param_cstarp
   PUBLIC :: param_ibetap

   ! Fitting parameters - hydro z indicies
   PUBLIC :: param_alphaz
   PUBLIC :: param_betaz
   PUBLIC :: param_epsz
   PUBLIC :: param_eps2z
   PUBLIC :: param_Gammaz
   PUBLIC :: param_Zammaz
   PUBLIC :: param_M0z
   PUBLIC :: param_Astarz
   PUBLIC :: param_Twhimz
   PUBLIC :: param_cstarz
   PUBLIC :: param_Mstarz
   PUBLIC :: param_etaz
   PUBLIC :: param_ibetaz
   PUBLIC :: param_gbetaz

   ! Fitting parameters - HMcode
   PUBLIC :: param_HMcode_Dv0
   PUBLIC :: param_HMcode_Dv1
   PUBLIC :: param_HMcode_dc0
   PUBLIC :: param_HMcode_dc1
   PUBLIC :: param_HMcode_eta0
   PUBLIC :: param_HMcode_eta1
   PUBLIC :: param_HMcode_f0
   PUBLIC :: param_HMcode_f1
   PUBLIC :: param_HMcode_kstar
   PUBLIC :: param_HMcode_As
   PUBLIC :: param_HMcode_alpha0
   PUBLIC :: param_HMcode_alpha1
   PUBLIC :: param_HMcode_alpha2
   PUBLIC :: param_HMcode_Dvnu
   PUBLIC :: param_HMcode_dcnu
   PUBLIC :: param_HMcode_mbar
   PUBLIC :: param_HMcode_nbar
   PUBLIC :: param_HMcode_Amf
   PUBLIC :: param_HMcode_Amfz
   PUBLIC :: param_HMcode_sbar
   PUBLIC :: param_HMcode_STp
   PUBLIC :: param_HMcode_STq
   PUBLIC :: param_HMcode_kd
   PUBLIC :: param_HMcode_kdp
   PUBLIC :: param_HMcode_nd
   PUBLIC :: param_HMcode_Ap
   PUBLIC :: param_HMcode_Ac
   PUBLIC :: param_HMcode_kp
   PUBLIC :: param_HMcode_mbarz
   PUBLIC :: param_HMcode_sbarz
   PUBLIC :: param_HMcode_fM
   PUBLIC :: param_HMcode_fD
   PUBLIC :: param_HMcode_zD
   PUBLIC :: param_HMcode_Az

   ! Fitting parameters - PT
   PUBLIC :: param_PT_A
   PUBLIC :: param_PT_alpha
   PUBLIC :: param_PT_beta

   ! Halo-model stuff that needs to be recalculated for each new z
   TYPE halomod

      ! Redshift and scale factor
      REAL :: z, a

      ! Switches
      INTEGER :: ip2h, ip1h, ibias, imf, iconc, iDolag, iAs, i2hcor
      INTEGER :: idc, iDv, ieta, i2hdamp, i1hdamp, itrans, ikdamp
      
      ! Flags for sigma 
      INTEGER :: flag_sigma

      ! Void stuff
      LOGICAL :: add_voids

      ! Spherical collapse parameters
      REAL :: dc, Dv
      LOGICAL :: mass_dependent_Dv

      ! HMx baryon parameters
      LOGICAL :: fix_star_concentration, different_Gammas
      REAL :: Theat_array(3), pivot_mass
      REAL :: alpha, alphap, alphaz, alpha_array(3), alphap_array(3), alphaz_array(3)
      REAL :: beta, betap, betaz
      REAL :: eps, epsz, eps_array(3), epsz_array(3)
      REAL :: eps2, eps2z, eps2_array(3), eps2z_array(3)
      REAL :: Gamma, Gammap, Gammaz, Gamma_array(3), Gammap_array(3), Gammaz_array(3)
      REAL :: Zamma, Zammap, Zammaz, Zamma_array(3), Zammap_array(3), Zammaz_array(3)
      REAL :: M0, M0z, M0_array(3), M0z_array(3) 
      REAL :: Astar, Astarz, Astar_array(3), Astarz_array(3)
      REAL :: Twhim, Twhimz, Twhim_array(3), Twhimz_array(3)
      REAL :: cstar, cstarp, cstarz, cstar_array(3), cstarp_array(3), cstarz_array(3)
      REAL :: sstar
      REAL :: Mstar, Mstarz, Mstar_array(3), Mstarz_array(3)
      REAL :: fcold
      REAL :: fhot
      REAL :: eta, etaz, eta_array(3), etaz_array(3)
      REAL :: ibeta, ibetap, ibetaz
      REAL :: gbeta, gbetaz

      ! Tilman 2018 work
      ! TODO: Probably redundant
      !REAL :: Theat
      REAL :: A_alpha, B_alpha, C_alpha, D_alpha, E_alpha      ! Tilman alpha parameters
      REAL :: A_eps, B_eps, C_eps, D_eps                       ! Tilman eps parameters
      REAL :: A_Gamma, B_Gamma, C_Gamma, D_Gamma, E_gamma      ! Tilman Gamma parameters
      REAL :: A_M0, B_M0, C_M0, D_M0, E_M0                     ! Tilman M0 parameters
      REAL :: A_Astar, B_Astar, C_Astar, D_Astar               ! Tilman Astar parameters
      REAL :: A_Twhim, B_Twhim, C_Twhim, D_Twhim               ! Tilman Twhim parameters

      ! Look-up tables
      REAL :: mmin, mmax
      REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:), log_m(:)!, mr(:)
      REAL, ALLOCATABLE :: r500(:), m500(:), c500(:), r200(:), m200(:), c200(:)
      REAL, ALLOCATABLE :: r500c(:), m500c(:), c500c(:), r200c(:), m200c(:), c200c(:)
      INTEGER :: n

      ! Window-function (not used?)
      REAL, ALLOCATABLE :: k(:), wk(:, :, :)
      INTEGER :: nk

      ! HMcode parameters and experimental parameters
      REAL :: knl, rnl, mnl, Rh, Mh, Mp, sigv, Rhh, neff!, ncur

      ! Saturation parameters (e.g., WDM)
      REAL :: nu_saturation
      LOGICAL :: saturation

      ! HOD etc.
      REAL :: mhalo_min, mhalo_max, HImin, HImax, rcore, hmass
      REAL :: n_c, n_s, n_g, rho_HI, sigma_conc

      ! Mass function integrals
      REAL :: gmin, gmax, gbmin, gbmax, gnorm

      ! HMcode parameters
      REAL :: Dv0, Dv1, dc0, dc1, eta0, eta1, f0, f1, ks, As, alp0, alp1, Dvnu, dcnu
      REAL :: HMcode_kstar, HMcode_fdamp, HMcode_kdamp, HMcode_alpha, HMcode_eta, HMcode_A
      LOGICAL :: DMONLY_baryon_recipe, DMONLY_neutrino_halo_mass_correction
      LOGICAL :: DMONLY_neutrinos_affect_virial_radius

      ! HMcode (2015) one-parameter baryon model (A, eta0)
      LOGICAL :: one_parameter_baryons

      ! HMcode (2020) parameters
      REAL :: kd, kdp, Ap, Ac, kp, nd, alp2, fM, fD, zD, Az
      REAL :: mbar, nbar, sbar, mbarz, sbarz
      REAL :: mbar_T, As_T, sbar_T, mbarz_T, sbarz_T, Az_T

      ! Halo types
      INTEGER :: halo_DMONLY, halo_CDM, halo_normal_bound_gas, halo_cold_bound_gas, halo_hot_bound_gas, halo_ejected_gas
      INTEGER :: halo_central_stars, halo_satellite_stars, halo_HI, halo_neutrino, halo_satellites
      INTEGER :: halo_void, halo_compensated_void, electron_pressure

      ! Halo components
      INTEGER :: normalise_baryons
      INTEGER :: frac_central_stars, frac_stars, frac_HI
      INTEGER :: frac_bound_gas, frac_cold_bound_gas, frac_hot_bound_gas

      ! Different 'has' logicals
      LOGICAL :: has_HI, has_galaxies, has_mass_conversions, safe_negative

      LOGICAL :: simple_pivot
      INTEGER :: response_baseline, response_denominator
      LOGICAL :: response_matter_only
      REAL :: acc, small_nu, large_nu
      CHARACTER(len=256) :: name
      INTEGER :: HMx_mode

      ! Halo mass function and linear halo bias
      REAL :: Tinker_bigA, Tinker_bigB, Tinker_bigC, Tinker_a, Tinker_b, Tinker_c
      REAL :: Tinker_alpha, Tinker_beta, Tinker_gamma, Tinker_phi, Tinker_eta
      REAL :: alpha_numu
      REAL :: ST_p, ST_q, ST_A, Amf, Amfz
      LOGICAL :: has_mass_function

      ! Non-linear halo bias
      TYPE(interpolator3D) :: bnl, bnl_lownu, bnl_lowk, bnl_lownu_lowk
      INTEGER :: iextrap_bnl
      LOGICAL :: has_bnl, stitch_bnl_nu
      CHARACTER(len=256) :: bnl_cat, bnl_path

      ! Perturbation theory
      REAL :: PT_alpha, PT_beta, PT_A

      ! HALOFIT
      REAL :: HALOFIT_knl, HALOFIT_neff, HALOFIT_ncur

   END TYPE halomod

   ! General
   INTEGER, PARAMETER :: nhalomod_large = 1000 ! Needs to be larger than the total number of defined 'halo models' TODO: Remove

   ! Halo window function integration
   ! NOTE: acc_win governs the speed of calculations with non-analytic halo profile Fourier transforms
   ! NOTE: It seems that acc_win can be downgraded significantly without impacting accuracy too much
   REAL, PARAMETER :: acc_win = 1e-3           ! Window-function integration accuracy parameter
   INTEGER, PARAMETER :: imeth_win = 12        ! Window-function integration method
   INTEGER, PARAMETER :: winint_order = 3      ! Window-function integration order
   REAL, PARAMETER :: winint_test_seconds = 1. ! Approximately how many seconds should each timing test take
   INTEGER, PARAMETER :: nmeth_win = 13        ! Total number of different winint methods
   INTEGER, PARAMETER :: nlim_bumps = 2        ! Do the bumps approximation after this number of bumps
   LOGICAL, PARAMETER :: winint_exit = .FALSE. ! Exit when the contributions to the integral become small

   ! Delta_v
   REAL, PARAMETER :: M0_Dv_default = 1e14 ! If Delta_v is halo-mass dependent (e.g., MG) then write it at this mass

   ! Mass function
   INTEGER, PARAMETER :: iorder_derivative_mass_function = 3          ! Polynomial order to calculate halo-mass function via numerical derivative
   INTEGER, PARAMETER :: ifind_derivative_mass_function = ifind_split ! Finding scheme for derivative

   ! Diagnostics
   REAL, PARAMETER :: mmin_diag = 1e10 ! Minimum halo mass for diagnostic tests [Msun/h]
   REAL, PARAMETER :: mmax_diag = 1e16 ! Maximum halo mass for diagnostic tests [Msun/h]
   INTEGER, PARAMETER :: n_diag = 101  ! Number of masses to output during diagnostics

   ! Halomodel
   LOGICAL, PARAMETER :: verbose_convert_mass = .FALSE.  ! Verbosity when running the mass-definition-conversion routines
   REAL, PARAMETER :: eps_p2h = 1e-3                     ! Tolerance for requal for lower limit of two-halo power
   REAL, PARAMETER :: g_integral_limit = -1e-4           ! Mininum allowed value for g(nu) integral
   REAL, PARAMETER :: gb_integral_limit = -1e-4          ! Mininum allowed value for g(nu)b(nu) integral
   LOGICAL, PARAMETER :: check_mass_function = .FALSE.   ! Check properties of the mass function
   INTEGER, PARAMETER :: iorder_integration = 3          ! Order for standard integration
   INTEGER, PARAMETER :: iorder_2halo_integration = 3    ! Order for 2-halo integration
   INTEGER, PARAMETER :: iorder_1halo_integration = 1    ! Order for 1-halo integration (basic because of wiggles)
   LOGICAL, PARAMETER :: calculate_Omega_stars = .FALSE. ! Calculate Omega_* in halomod_init

   ! Non-linear radius: solves nu(R_nl)=1
   INTEGER, PARAMETER :: iorder_inversion_rnl = 3                 ! Order for interpolation inversion
   INTEGER, PARAMETER :: ifind_inversion_rnl = ifind_split        ! Finding scheme for table
   INTEGER, PARAMETER :: iinterp_inversion_rnl = iinterp_Lagrange ! Interpolation method
   LOGICAL, PARAMETER :: rnl_extrap = .FALSE.                     ! Extrapolate to find r_nl or use nu(1)? REVERT

   ! Voids
   REAL, PARAMETER :: void_underdensity = -1.      ! Void underdensity
   REAL, PARAMETER :: void_compensation = 1.1      ! How much larger than the void is the compensation region?
   REAL, PARAMETER :: simple_void_radius = 10.     ! Void radius [Mpc/h]
   LOGICAL, PARAMETER :: compensate_voids = .TRUE. ! Do we compenate voids with ridges?
   LOGICAL, PARAMETER :: simple_voids = .FALSE.    ! Use simple voids or not

   ! HMcode
   REAL, PARAMETER :: HMcode_fdamp_min = 0.00 ! Minimum value for f_damp parameter
   REAL, PARAMETER :: HMcode_fdamp_max = 0.99 ! Maximum value for f_damp parameter
   REAL, PARAMETER :: HMcode_alpha_min = 0.50 ! Minimum value for alpha transition parameter
   REAL, PARAMETER :: HMcode_alpha_max = 2.00 ! Maximum value for alpha transition parameter
   REAL, PARAMETER :: HMcode_ks_limit = 7.    ! Limit for (k/k*)^2 in one-halo term

   ! HMcode baryon model
   LOGICAL, PARAMETER :: renormalise_onehalo = .FALSE. ! Enforce that one-halo amplitude is unaffected by mass loss

   ! HMx
   REAL, PARAMETER :: HMx_alpha_min = 1e-2 ! Minimum alpha parameter; needs to be set at not zero
   REAL, PARAMETER :: HMx_beta_min = 1e-2  ! Minimum alpha parameter; needs to be set at not zero
   REAL, PARAMETER :: HMx_Gamma_min = 1.10 ! Minimum polytropic index
   REAL, PARAMETER :: HMx_Gamma_max = 2.00 ! Maximum polytropic index
   REAL, PARAMETER :: HMx_eps_min = -0.99  ! Minimum concentration modification
   REAL, PARAMETER :: HMx_Astar_min = 0.   ! Minimum halo star fraction (I think this can be zero)3
   REAL, PARAMETER :: HMx_eta_min = 0.     ! Minimum value for star index thing
   REAL, PARAMETER :: frac_min = -1e-5     ! Fractions cannot be less than frac min (slightly below zero for roundoff)
   REAL, PARAMETER :: frac_max = 1.        ! Fractions cannot be greater than frac max (exactly unity)

   ! Minimum fraction before haloes are treated as delta functions (dangerous)
   ! TODO: Hydro tests passed if 1e-4, but I worry a bit, also it was a fairly minor speed up (25%)
   ! TODO: Also, not obvious that a constant frac_min is correct because some species have low abundance
   ! but we ususally care about perturbations to this abundance
   REAL, PARAMETER :: frac_min_delta = 0. ! Be very careful with this
   INTEGER, PARAMETER :: iorder_delta = 3 !
   INTEGER, PARAMETER :: ifind_delta = 3  !
   INTEGER, PARAMETER :: imeth_delta = 2  !

   ! Halo types
   LOGICAL, PARAMETER :: verbose_galaxies = .FALSE. ! Verbosity when doing the galaxies initialisation
   LOGICAL, PARAMETER :: verbose_HI = .TRUE.        ! Verbosity when doing the HI initialisation

   ! Non-linear halo bias
   INTEGER, PARAMETER :: bnl_method_a = 1           ! Assign B_NL based on scale factor
   INTEGER, PARAMETER :: bnl_method_sig8 = 2        ! Assign B_NL based on sigma_8
   INTEGER, PARAMETER :: bnl_method_rescale = 3     ! Assign B_NL based on AW10 rescaling
   INTEGER, PARAMETER :: bnl_method = bnl_method_rescale  ! Choose method to assign the Multidark measured B_NL
   LOGICAL, PARAMETER :: fixed_bnl = .FALSE.        ! Is the non-linear bias correction taken at fixed z = 0?
   LOGICAL, PARAMETER :: add_I_11 = .TRUE.          ! Add integral below numin, numin in halo model calculation
   LOGICAL, PARAMETER :: add_I_12_and_I_21 = .TRUE. ! Add integral below numin in halo model calculation
   REAL, PARAMETER :: kmin_bnl = 8e-2               ! Below this wavenumber force B_NL to zero
   REAL, PARAMETER :: numin_bnl = 0.                ! Below this halo mass force  B_NL to zero
   REAL, PARAMETER :: numax_bnl = 10.               ! Above this halo mass force  B_NL to zero
   LOGICAL, PARAMETER :: exclusion_bnl = .FALSE.    ! Attempt to manually include damping from halo exclusion
   LOGICAL, PARAMETER :: fix_minimum_bnl = .FALSE.  ! Force a minimum value for B_NL
   REAL, PARAMETER :: min_bnl = -1.                 ! Minimum value for B_NL (could be below -1; cross spectra can be negative)
   INTEGER, PARAMETER :: iorder_bnl = 1             ! 1 - Linear interpolation
   INTEGER, PARAMETER :: iextrap_bnl = iextrap_lin  ! Linear extrapolation
   LOGICAL, PARAMETER :: store_bnl = .FALSE.        ! Storage interpolator mode?
   REAL, PARAMETER :: bnl_rescale_R1 = 1.           ! Minimum R for rescaling [Mpc/h]
   REAL, PARAMETER :: bnl_rescale_R2 = 10.          ! Maximum R for rescaling [Mpc/h]
   INTEGER, PARAMETER :: icos_bnl_rescale = 37      ! 37 - Multidark
   CHARACTER(len=256), PARAMETER :: dir_bnl = 'BNL' ! Directory containing BNL measurements
   !CHARACTER(len=256), PARAMETER :: cat = 'rockstar' ! Halo catalogue
   !CHARACTER(len=256), PARAMETER :: base_bnl_lownu = '/Users/Mead/Physics/Multidark/data/'//trim(dir_bnl)//'/M512/Bolshoi_'//trim(cat_bnl)
   !CHARACTER(len=256), PARAMETER :: base_bnl = '/Users/Mead/Physics/Multidark/data/'//trim(dir_bnl)//'/M512/MDR1_'//trim(cat_bnl)

   ! HALOFIT (can be used as two-halo term)
   INTEGER, PARAMETER :: HALOFIT_twohalo = HALOFIT_Takahashi

   ! Field types
   INTEGER, PARAMETER :: field_dmonly = 1
   INTEGER, PARAMETER :: field_matter = 2
   INTEGER, PARAMETER :: field_cdm = 3
   INTEGER, PARAMETER :: field_gas = 4
   INTEGER, PARAMETER :: field_stars = 5
   INTEGER, PARAMETER :: field_bound_gas = 6
   INTEGER, PARAMETER :: field_ejected_gas = 7
   INTEGER, PARAMETER :: field_electron_pressure = 8
   INTEGER, PARAMETER :: field_void = 9
   INTEGER, PARAMETER :: field_compensated_void = 10
   INTEGER, PARAMETER :: field_central_galaxies = 11
   INTEGER, PARAMETER :: field_satellite_galaxies = 12
   INTEGER, PARAMETER :: field_galaxies = 13
   INTEGER, PARAMETER :: field_HI = 14
   INTEGER, PARAMETER :: field_cold_bound_gas = 15
   INTEGER, PARAMETER :: field_hot_bound_gas = 16
   INTEGER, PARAMETER :: field_normal_bound_gas = 17
   INTEGER, PARAMETER :: field_central_stars = 18
   INTEGER, PARAMETER :: field_satellite_stars = 19
   INTEGER, PARAMETER :: field_CIB_353 = 20
   INTEGER, PARAMETER :: field_CIB_545 = 21
   INTEGER, PARAMETER :: field_CIB_857 = 22
   INTEGER, PARAMETER :: field_halo_11p0_11p5 = 23
   INTEGER, PARAMETER :: field_halo_11p5_12p0 = 24
   INTEGER, PARAMETER :: field_halo_12p0_12p5 = 25
   INTEGER, PARAMETER :: field_halo_12p5_13p0 = 26
   INTEGER, PARAMETER :: field_halo_13p0_13p5 = 27
   INTEGER, PARAMETER :: field_halo_13p5_14p0 = 28
   INTEGER, PARAMETER :: field_halo_14p0_14p5 = 29
   INTEGER, PARAMETER :: field_halo_14p5_15p0 = 30
   INTEGER, PARAMETER :: field_neutrino = 31
   INTEGER, PARAMETER :: field_n = 31

   ! Parameters to pass to minimization routines
   INTEGER, PARAMETER :: param_alpha = 1
   INTEGER, PARAMETER :: param_eps = 2
   INTEGER, PARAMETER :: param_gamma = 3
   INTEGER, PARAMETER :: param_M0 = 4
   INTEGER, PARAMETER :: param_Astar = 5
   INTEGER, PARAMETER :: param_Twhim = 6
   INTEGER, PARAMETER :: param_cstar = 7
   INTEGER, PARAMETER :: param_fcold = 8
   INTEGER, PARAMETER :: param_mstar = 9
   INTEGER, PARAMETER :: param_sstar = 10
   INTEGER, PARAMETER :: param_alphap = 11
   INTEGER, PARAMETER :: param_Gammap = 12
   INTEGER, PARAMETER :: param_cstarp = 13
   INTEGER, PARAMETER :: param_fhot = 14
   INTEGER, PARAMETER :: param_alphaz = 15
   INTEGER, PARAMETER :: param_Gammaz = 16
   INTEGER, PARAMETER :: param_M0z = 17
   INTEGER, PARAMETER :: param_Astarz = 18
   INTEGER, PARAMETER :: param_Twhimz = 19
   INTEGER, PARAMETER :: param_eta = 20
   INTEGER, PARAMETER :: param_HMcode_Dv0 = 21
   INTEGER, PARAMETER :: param_HMcode_Dv1 = 22
   INTEGER, PARAMETER :: param_HMcode_dc0 = 23
   INTEGER, PARAMETER :: param_HMcode_dc1 = 24
   INTEGER, PARAMETER :: param_HMcode_eta0 = 25
   INTEGER, PARAMETER :: param_HMcode_eta1 = 26
   INTEGER, PARAMETER :: param_HMcode_f0 = 27
   INTEGER, PARAMETER :: param_HMcode_f1 = 28
   INTEGER, PARAMETER :: param_HMcode_kstar = 29
   INTEGER, PARAMETER :: param_HMcode_As = 30
   INTEGER, PARAMETER :: param_HMcode_alpha0 = 31
   INTEGER, PARAMETER :: param_HMcode_alpha1 = 32
   INTEGER, PARAMETER :: param_epsz = 33
   INTEGER, PARAMETER :: param_beta = 34
   INTEGER, PARAMETER :: param_betap = 35
   INTEGER, PARAMETER :: param_betaz = 36
   INTEGER, PARAMETER :: param_etaz = 37
   INTEGER, PARAMETER :: param_cstarz = 38
   INTEGER, PARAMETER :: param_mstarz = 39
   INTEGER, PARAMETER :: param_ibeta = 40
   INTEGER, PARAMETER :: param_ibetap = 41
   INTEGER, PARAMETER :: param_ibetaz = 42
   INTEGER, PARAMETER :: param_gbeta = 43
   INTEGER, PARAMETER :: param_gbetaz = 44
   INTEGER, PARAMETER :: param_HMcode_Dvnu = 45
   INTEGER, PARAMETER :: param_HMcode_dcnu = 46
   INTEGER, PARAMETER :: param_eps2 = 47
   INTEGER, PARAMETER :: param_eps2z = 48
   INTEGER, PARAMETER :: param_Zamma = 49
   INTEGER, PARAMETER :: param_Zammap = 50
   INTEGER, PARAMETER :: param_Zammaz = 51
   INTEGER, PARAMETER :: param_HMcode_mbar = 52
   INTEGER, PARAMETER :: param_HMcode_nbar = 53
   INTEGER, PARAMETER :: param_HMcode_Amf = 54
   INTEGER, PARAMETER :: param_HMcode_sbar = 55
   INTEGER, PARAMETER :: param_HMcode_STp = 56
   INTEGER, PARAMETER :: param_HMcode_STq = 57
   INTEGER, PARAMETER :: param_HMcode_kd = 58
   INTEGER, PARAMETER :: param_HMcode_Ap = 59
   INTEGER, PARAMETER :: param_HMcode_Ac = 60
   INTEGER, PARAMETER :: param_HMcode_kp = 61
   INTEGER, PARAMETER :: param_HMcode_kdp = 62
   INTEGER, PARAMETER :: param_HMcode_mbarz = 63
   INTEGER, PARAMETER :: param_HMcode_sbarz = 64
   INTEGER, PARAMETER :: param_HMcode_Amfz = 65
   INTEGER, PARAMETER :: param_HMcode_nd = 66
   INTEGER, PARAMETER :: param_HMcode_alpha2 = 67
   INTEGER, PARAMETER :: param_HMcode_fM = 68
   INTEGER, PARAMETER :: param_HMcode_fD = 69
   INTEGER, PARAMETER :: param_HMcode_zD = 70
   INTEGER, PARAMETER :: param_HMcode_Az = 71
   INTEGER, PARAMETER :: param_PT_A = 72
   INTEGER, PARAMETER :: param_PT_alpha = 73
   INTEGER, PARAMETER :: param_PT_beta = 74
   INTEGER, PARAMETER :: param_n = 74

   ! HMcode versions
   INTEGER, PARAMETER :: HMcode2015 = 7
   INTEGER, PARAMETER :: HMcode2015_CAMB = 66 
   INTEGER, PARAMETER :: HMcode2016 = 1
   INTEGER, PARAMETER :: HMcode2016_CAMB = 51
   INTEGER, PARAMETER :: HMcode2016_neutrinofix = 92
   INTEGER, PARAMETER :: HMcode2016_neutrinofix_CAMB = 95
   INTEGER, PARAMETER :: HMcode2018 = 53
   INTEGER, PARAMETER :: HMcode2019 = 15
   INTEGER, PARAMETER :: HMcode2020 = 79
   INTEGER, PARAMETER :: HMcode2020_CAMB = 123
   INTEGER, PARAMETER :: HMcode2020_feedback = 103
   INTEGER, PARAMETER :: HMcode2020_feedback_CAMB = 124

   ! HMx versions
   INTEGER, PARAMETER :: HMx2020_matter_w_temp_scaling = 60
   INTEGER, PARAMETER :: HMx2020_matter_pressure_w_temp_scaling = 61

CONTAINS

   SUBROUTINE assign_halomod(ihm, hmod, verbose)

      ! TODO: nhalomod should not be necessary
      INTEGER, INTENT(INOUT) :: ihm
      TYPE(halomod), INTENT(OUT) :: hmod
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i
      INTEGER, PARAMETER :: nhalomod = nhalomod_large ! Some large-enough integer
      CHARACTER(len=256), ALLOCATABLE :: names(:)

      ! Allocate array for names (warning from gfortran 10)
      ! TODO: This is lazy, there must be a better way
      ALLOCATE(names(nhalomod))
      names = ''

      names(1) =  'HMcode (2016)'
      names(2) =  'Basic halomodel (Two-halo term is linear; Dv=200; dc=1.686; c(M) Bullock)'
      names(3) =  'Standard halomodel (Seljak 2000)'
      names(4) =  'Standard halomodel but with HMcode (2015) smoothed transition'
      names(5) =  'Standard halomodel but with Dv=200 and dc=1.686 and Bullock c(M)'
      names(6) =  'Half-accurate HMcode (HMcode 2015, 2016)'
      names(7) =  'HMcode (2015)'
      names(8) =  'Including scatter in halo properties at fixed mass'
      names(9) =  'Parameters for CCL tests (high accuracy)'
      names(10) = 'Comparison of mass conversions with Wayne Hu code'
      names(11) = 'UPP for electron pressure'
      names(12) = 'Spherical collapse calculation for Mead (2017) results'
      names(13) = 'Experimental sigmoid transition'
      names(14) = 'Experimental scale-dependent halo bias'
      names(15) = 'HMcode (2019)'
      names(16) = 'Halo-void model'
      names(17) = 'Tilman HMx - AGN 7.6'
      names(18) = 'Tilman HMx - AGN tuned'
      names(19) = 'Tilman HMx - AGN 8.0'
      names(20) = 'Standard halo-model (Seljak 2000) in response'
      names(21) = 'Cored profile model'
      names(22) = 'Delta function-NFW star profile model response'
      names(23) = 'Tinker (20xx) mass function; virial mass'
      names(24) = 'Non-linear halo bias for M200c haloes'
      names(25) = 'Villaescusa-Navarro HI halomodel'
      names(26) = 'Delta-function mass function'
      names(27) = 'Press & Schecter mass function'
      names(28) = 'HMcode (2016) with one-parameter OWLS AGN model'
      names(29) = 'Adding in cold gas'
      names(30) = 'Adding in hot gas'
      names(31) = 'HMcode (2016) with damped BAO'
      names(32) = 'HMx: AGN 7.6; fixed z; simple pivot'
      names(33) = 'HMx: AGN tuned; fixed z; simple pivot'
      names(34) = 'HMx: AGN 8.0; fixed z; simple pivot'
      names(35) = 'HMx: AGN 7.6; fixed z; clever pivot; f_hot'
      names(36) = 'HMx: AGN tuned; fixed z; clever pivot; f_hot'
      names(37) = 'HMx: AGN 8.0; fixed z; clever pivot; f_hot'
      names(38) = 'HMx: AGN 7.6'
      names(39) = 'HMx: AGN tuned'
      names(40) = 'HMx: AGN 8.0'
      names(41) = 'Put some galaxy mass in the halo/satellites'
      names(42) = 'Tinker (20xx) mass function and bias; M200c'
      names(43) = 'Standard halo-model (Seljak 2000) in matter response'
      names(44) = 'Tinker (20xx) mass function and bias; M200'
      names(45) = 'No stars'
      names(46) = 'Isothermal beta model for gas'
      names(47) = 'Isothermal beta model for gas in response'
      names(48) = 'Non-linear halo bias for standard model'
      names(49) = 'Non-linear halo bias with Tinker (2010) with bias from peak-background split'
      names(50) = 'HMcode (2016) with Dolag pow=1 bug'
      names(51) = 'HMcode (2016) with CAMB parameters'
      names(52) = 'Standard but with Mead (2017) spherical collapse'
      names(53) = 'HMcode (2018)'
      names(54) = 'Matter masquerading as DMONLY'
      names(55) = 'HMx2020: Baseline'
      names(56) = 'HMx2020: Model that fits stars-stars AGN 7.6'
      names(57) = 'HMx2020: Model that fits stars-stars AGN 7.8'
      names(58) = 'HMx2020: Model that fits stars-stars AGN 8.0'   
      names(59) = 'HMx2020: Temperature-dependent model for stars'
      names(60) = 'HMx2020: Temperature-dependent model for matter'
      names(61) = 'HMx2020: Temperature-dependent model for matter, pressure'
      names(62) = 'HMx2020: Temperature-dependent model for matter, CDM, gas, stars'
      names(63) = 'HMx2020: Temperature-dependent model for everything'
      names(64) = 'HMcode (2016) but with HMcode (2020) baryon model'
      names(65) = 'HMx2020: Baseline, no response'
      names(66) = 'HMcode (2015) with CAMB parameters'
      names(67) = 'Dv=200 and dc=1.686'
      names(68) = 'Bullock c(M)'
      names(69) = 'Simple Bullock c(M)'
      names(70) = 'No Dolag correction'
      names(71) = 'Dolag with 1.5 exponent'
      names(72) = 'Linear two-halo term'
      names(73) = 'Linear two-halo term with damped wiggles'
      names(74) = 'dc=1.686 fixed'
      names(75) = 'dc from Mead (2017) fit'
      names(76) = 'Dv from Mead (2017) fit'
      names(77) = 'HMcode (2020) unfitted'
      names(78) = 'HMcode (2020) fitted to Cosmic Emu k<1'
      names(79) = 'HMcode (2020) fitted'
      names(80) = 'Jenkins mass function' 
      names(81) = 'HMx2020: matter-only response; Model that fits stars-stars AGN 7.6'
      names(82) = 'HMx2020: matter-only response; Model that fits stars-stars AGN 7.8'
      names(83) = 'HMx2020: matter-only response; Model that fits stars-stars AGN 8.0'
      names(84) = 'HMx2020: different Gammas; Model that fits stars-stars AGN 7.6'
      names(85) = 'HMx2020: different Gammas; Model that fits stars-stars AGN 7.8'
      names(86) = 'HMx2020: different Gammas; Model that fits stars-stars AGN 8.0'
      names(87) = 'Despali et al. (2016) mass function'
      names(88) = 'Child et al. (2018) concentration-mass relation'
      names(89) = 'All halo gas bound'
      names(90) = 'z-dependent Dolag correction'
      names(91) = 'Tinker (2008) mass function for virial haloes'
      names(92) = 'HMcode (2016) with neutrino halo-mass correction'
      names(93) = 'Warren (2006) mass function'
      names(94) = 'Reed (2007) mass function'
      names(95) = 'HMcode (2016) with neutrino halo-mass correction and CAMB mass range'
      names(96) = 'Philcox (2020) mass function'
      names(97) = 'Sigma calculated for all matter'
      names(98) = 'No DMONLY correction for neutrinos'
      names(99) = 'sigma calculated for normalised cold matter'
      names(100) = 'Massive neutrinos affect virial radius'
      names(101) = 'Spherical collapse calculation'
      names(102) = 'HMcode (2020) baryon model'
      names(103) = 'HMcode (2020) baryon model response'
      names(104) = 'HMcode (2020) baryon model using HMx language'
      names(105) = 'Tweaked baryon parameters'
      names(106) = 'Non-linear halo bias with tweaked baryon parameters'
      names(107) = 'Non-linear halo bias'
      names(108) = 'Bhattacharya et al. (2011) mass function'
      names(109) = 'Perturbation theory inspired two-halo term'
      names(110) = 'One-loop SPT for two-halo term'
      names(111) = 'Perturbation theory inspired and fitted two-halo term'
      names(112) = 'Tinker (2008) mass function with no z evolution'
      names(113) = 'Tinker (2010) mass function and peak-background split bias with no z evolution'
      names(114) = 'Non-linear halo bias'
      names(115) = 'Tinker (2010) mass function'
      names(116) = 'Non-linear halo bias with no extrapolation'
      names(117) = 'Non-linear halo bias with Bolshoi'
      names(118) = 'Non-linear halo bias with Bolshoi and no extrapolation'
      names(119) = 'Non-linear power in two-halo term'
      names(120) = 'Quasi-linear power in two-halo term'
      names(121) = 'Non-linear power in two-halo term with Tinker (2010) mass function'
      names(122) = 'Quasi-linear power in two-halo term with Tinker (2010) mass function'
      names(123) = 'HMcode (2020) with extended mass range'
      names(124) = 'HMcode (2020) feedback with extended mass range'
      names(125) = 'HMcode (2020) unfitted with extended mass range'

      IF (verbose) WRITE (*, *) 'ASSIGN_HALOMOD: Assigning halomodel'

      ! Default options
      hmod%mmin = 1e7  ! Lower mass limit for integration [Msun/h]
      hmod%mmax = 1e17 ! Upper mass limit for integration [Msun/h]
      hmod%n = 128     ! Number of points in integration (128 is okay, 1024 is better)
      hmod%acc = 1e-4  ! Accuracy for continuous integrals (1e-3 is okay, 1e-4 is better)

      ! Small and large values for nu (6 is okay, corrections are suppressed by exp(-large_nu^2)
      ! TODO: It would be nice to be able to fix small_nu = 0.
      hmod%small_nu = 1e-6
      hmod%large_nu = 6.

      ! Two-halo term
      ! 0 - No two-halo term
      ! 1 - Linear theory
      ! 2 - Standard from Seljak (2000)
      ! 3 - Linear theory with damped wiggles
      ! 4 - One-loop SPT, dewiggle and damp
      ! 5 - One-loop SPT
      ! 6 - Non-linear theory from HALOFIT
      ! 7 - Quasi-linear term from HALOFIT
      hmod%ip2h = 2

      ! One-halo term
      ! 0 - No one-halo term
      ! 1 - Standard one-halo term
      hmod%ip1h = 1

      ! Method to correct the two-halo integral
      ! 0 - Do nothing
      ! 1 - Add value of missing integral assuming that W(k)=1 (only works for matter fields)
      ! 2 - Put the missing part of the integrand as a delta function at lower mass limit
      hmod%i2hcor = 2

      ! Order of halo bias to go to
      ! 1 - Linear order (standard)
      ! 2 - Second order
      ! 3 - Full non-linear halo bias
      ! 4 - Fudge from Fedeli (2014b)
      ! 5 - Fudge from me
      hmod%ibias = 1

      ! One-halo term large-scale damping
      ! 0 - No damping
      ! 1 - HMcode (2015, 2016) exponential damping
      ! 2 - k^4 at large scales
      ! 3 - HMcode (2020) k^4 at large scales with sigma8 dependence on wavenumber
      hmod%i1hdamp = 0

      ! Halo mass function and bias
      !  1 - Press & Schecter (1974)
      !  2 - Sheth & Tormen (1999)
      !  3 - Tinker et al. (2010) with z dependence
      !  4 - Delta function in mass
      !  5 - Jenkins et al. (2001)
      !  6 - Despali et al. (2016)
      !  7 - Tinker et al. (2008) with z dependence
      !  8 - Warren et al. (2006; astro-ph/0506395)
      !  9 - Reed et al. (2007; astro-ph/0607150)
      ! 10 - Bhattacharya et al. (2011; 1005.2239)
      ! 11 - Tinker et al. (2010) but with no halo bias
      ! 12 - Sheth & Tormen (2002)
      ! 13 - Sheth & Tormen form, but with free parameters
      ! 14 - Philcox et al. (2020)
      ! 15 - Tinker et al. (2008) no z dependence
      ! 16 - Tinker et al. (2010) no z dependence
      ! 17 - Tinker et al. (2010) with calibrated bias
      ! 18 - Tinker et al. (2010) with calibrated bias and no z dependence
      ! 19 - Sheth & Tormen (1999) with calibrated bias
      hmod%imf = 2

      ! Concentration-mass relation
      ! 1 - Full Bullock et al. (2001; astro-ph/9909159)
      ! 2 - Simple Bullock et al. (2001; astro-ph/9909159)
      ! 3 - Duffy et al. (2008; astro-ph/0804.2486): full 200m
      ! 4 - Duffy et al. (2008; astro-ph/0804.2486): full virial
      ! 5 - Duffy et al. (2008; astro-ph/0804.2486): full 200c
      ! 6 - Duffy et al. (2008; astro-ph/0804.2486): relaxed 200m
      ! 7 - Duffy et al. (2008; astro-ph/0804.2486): relaxed virial
      ! 8 - Duffy et al. (2008; astro-ph/0804.2486): relaxed 200c
      hmod%iconc = 4

      ! How to calculate sigma(R); either cold matter (two definitions) or all matter
      hmod%flag_sigma = flag_ucold

      ! Linear collapse threshold delta_c
      ! 1 - Fixed 1.686
      ! 2 - Nakamura & Suto (1997) fitting function
      ! 3 - HMcode (2016)
      ! 4 - Mead (2017) fitting function
      ! 5 - Spherical-collapse calculation
      ! 6 - HMcode (2015)
      hmod%idc = 2

      ! Virial density Delta_v
      ! 1 - Fixed 200
      ! 2 - Bryan & Norman (1998) fitting function
      ! 3 - HMcode (2016)
      ! 4 - Mead (2017) fitting function
      ! 5 - Spherical-collapse calculation
      ! 6 - Fixed to unity to give Lagrangian radius
      ! 7 - For M200c
      ! 8 - HMcode (2015)
      ! 9 - Fixed to 18pi^2 ~ 178
      hmod%iDv = 2

      ! eta for halo window function
      ! 0 - No
      ! 1 - HMcode (2015, 2016)
      ! 2 - Same as 1 but with sigma_cold dependence
      ! 3 - HMcode (2020)
      hmod%ieta = 1

      ! Concentration-mass rescaling
      ! 0 - No
      ! 1 - HMcode (2015, 2016)
      ! 2 - HMcode (2020) with sigma8 dependence
      hmod%iAs = 0

      ! Two-halo term damping
      ! 0 - No
      ! 1 - HMcode (2015)
      ! 2 - HMcode (2016)
      ! 3 - HMcode (2020) perturbation theory and exponential damping
      hmod%i2hdamp = 0

      ! alpha for two- to one-halo transition region
      ! 0 - No
      ! 1 - Smoothed transition with alpha
      ! 2 - ?
      ! 3 - New HMx transition
      ! 4 - Tanh transition
      ! 5 - HMcode (2020) alpha
      hmod%itrans = 0

      ! Use the Dolag et al. (2004; astro-ph/0309771) c(M) correction for dark energy?
      ! 0 - No
      ! 1 - Exactly as in Dolag et al. (2004)
      ! 2 - As in Dolag et al. (2004) but with a 1.5 exponent
      ! 3 - Dolag (2004) with redshift dependence
      ! 4 - Dolag (2004) with redshift dependence and neutrinos removed
      hmod%iDolag = 1

      ! How to ensure that we do not have more baryons than Om_b
      ! 1 - Remove excess mass from stars
      ! 2 - Limit baryon reserves for gas (Arico et al. 2020; 1911.08471)
      ! 3 - Remove excess mass from bound gas
      hmod%normalise_baryons = 1

      ! Halo bound gas fraction
      ! 1 - Fedeli (2014a) bound gas model
      ! 2 - Schneider & Teyssier (2015) bound gas
      ! 3 - Universal baryon fraction
      hmod%frac_bound_gas = 2

      ! Halo cold-bound gas fraction
      ! 1 - Constant fraction of bound gas is cold gas
      hmod%frac_cold_bound_gas = 1

      ! Halo hot-bound gas fraction
      ! 1 - Constant fraction of bound gas is hot gas
      hmod%frac_hot_bound_gas = 1

      ! Halo star fraction
      ! 1 - Fedeli (2014)
      ! 2 - Constant stellar fraction
      ! 3 - Fedeli (2014) but saturates at high halo mass
      hmod%frac_stars = 3

      ! Halo central-galaxy star fraction
      ! NOTE: If eta=0 then all stars are automatically in the central galaxies
      ! 1 - All stars in central galaxy
      ! 2 - Schneider et al. (2018) split between central galaxy and satellite galaxies
      hmod%frac_central_stars = 2

      ! Halo HI fraction
      ! 1 - Simple model
      ! 2 - Villaescusa-Navarro et al. (2018; 1804.09180)
      hmod%frac_HI = 1

      ! DMONLY halo profile
      ! 1 - Analyical NFW
      ! 2 - Non-analytical NFW (for testing W(k) numerical integration)
      ! 3 - Tophat
      ! 4 - Delta function
      ! 5 - Cored NFW
      ! 6 - Isothermal
      ! 7 - Shell
      hmod%halo_DMONLY = 1

      ! CDM halo profile
      ! 1 - NFW
      hmod%halo_CDM = 1

      ! Bound static gas halo profile
      ! 1 - Simplified Komatsu & Seljak (2001) gas model
      ! 2 - Isothermal beta model
      ! 3 - Full Komatsu & Seljak (2001) gas model
      ! 4 - NFW
      hmod%halo_normal_bound_gas = 1

      ! Bound cold gas halo profile
      ! 1 - Delta function
      hmod%halo_cold_bound_gas = 1

      ! Bound hot gas halo profile
      ! 1 - Isothermal
      hmod%halo_hot_bound_gas = 1

      ! Free gas halo profile
      ! 1 - Isothermal model (out to 2rv)
      ! 2 - Ejected gas model from Schneider & Teyssier (2015)
      ! 3 - Isothermal shell that connects electron pressure and density to bound gas at rv
      ! 4 - Komatsu-Seljak continuation
      ! 5 - Power-law continuation
      ! 6 - Cubic profile
      ! 7 - Smoothly distributed
      ! 8 - Delta function
      hmod%halo_ejected_gas = 7

      ! Different Gamma parameter for pressure and density
      hmod%different_Gammas = .FALSE.

      ! Central stars halo profile
      ! 1 - Fedeli (2014) stellar distribution
      ! 2 - Schneider & Teyssier (2015) stellar distribution
      ! 3 - Delta function
      ! 4 - Delta function at low mass and NFW at high mass
      hmod%halo_central_stars = 1

      ! Satellite stars halo profile
      ! 1 - NFW
      ! 2 - Fedeli (2014) stellar distribution
      hmod%halo_satellite_stars = 1

      ! Do stars see the original halo concentration or that modified by the HMx parameters?
      hmod%fix_star_concentration = .FALSE.

      ! Neutrino halo profile
      ! 1 - Smooth neutrino distribution
      ! 2 - NFW
      hmod%halo_neutrino = 1

      ! HI halo profile
      ! 1 - NFW
      ! 2 - Delta function
      ! 3 - Polynomial with internal exponential hole
      ! 4 - NFW with internal exponential hole
      ! 5 - Modified NFW
      hmod%halo_HI = 1

      ! Satellite galaxy halo profile
      ! 1 - NFW
      ! 2 - Isothermal
      hmod%halo_satellites = 1

      ! Electron pressure
      ! 1 - UPP (Arnaud et al. 2010)
      ! 2 - From gas profiles
      hmod%electron_pressure = 2

      ! Do voids?
      hmod%add_voids = .FALSE.

      ! Set the void model
      ! 1 - Top-hat void
      hmod%halo_void = 1

      ! Set the compensated void model
      ! 1 - Top-hat void
      hmod%halo_compensated_void = 1

      ! Safeguard against negative terms in cross correlations
      hmod%safe_negative = .FALSE.

      ! HMcode parameter defaults
      hmod%Dv0 = 200. ! This should really be from the Bryan & Norman formula
      hmod%Dv1 = 0.
      hmod%dc0 = 1.686
      hmod%dc1 = 0.
      hmod%eta0 = 0.
      hmod%eta1 = 0.
      hmod%f0 = 0.
      hmod%f1 = 0.
      hmod%ks = 0.
      hmod%As = 4.
      hmod%alp0 = 1.
      hmod%alp1 = 0.
      hmod%Dvnu = 0.
      hmod%dcnu = 0.

      ! HMcode (2020) additional parameters
      hmod%DMONLY_neutrino_halo_mass_correction = .TRUE.
      hmod%DMONLY_neutrinos_affect_virial_radius = .FALSE.
      hmod%kd = 0.
      hmod%kdp = 0.
      hmod%kp = 0.
      hmod%nd = 0.
      hmod%Ap = 0.
      hmod%Ac = 0.
      hmod%Az = 0.
      hmod%alp2 = 0.
      hmod%fM = 0.01
      hmod%fD = 1.
      hmod%zD = 100.   

      ! HMcode (2020) baryon recipe
      hmod%DMONLY_baryon_recipe = .FALSE.     
      hmod%mbar = 1e14
      hmod%nbar = 2.
      hmod%sbar = 1e-3 ! Problems with fitting if set to zero due to log
      hmod%mbarz = 0.
      hmod%sbarz = 0.
      hmod%Amfz = 0.
      hmod%mbar_T = 0.
      hmod%As_T = 0.
      hmod%sbar_T = 0.
      hmod%mbarz_T = 0.
      hmod%sbarz_T = 0.
      hmod%Az_T = 0.

      ! ~infinite redshift for Dolag correction
      !hmod%zinf_Dolag = 10.

      ! Baryon model (relates eta0 and A)
      hmod%one_parameter_baryons = .FALSE.

      ! HMx parameters dependence on redshift and AGN temperature
      ! 1 - HMcode 2020
      ! 2 - NOT USED
      ! 3 - Standard
      ! 4 - NOT SUPPORTED
      ! 5 - HMx2020 redshift scalings
      ! 6 - HMx2020 redshift and temperature scalings
      hmod%HMx_mode = 3

      ! Pivot is held fixed or is the non-linear mass (function of cosmology/z)?
      hmod%simple_pivot = .FALSE.
      hmod%pivot_mass = 1e14

      !!! HMx hydro parameters !!!

      ! Variations of HMx parameters with AGN temperature
      hmod%Theat_array(1) = 10**7.6
      hmod%Theat_array(2) = 10**7.8
      hmod%Theat_array(3) = 10**8.0

      ! Non-virial temperature correction for static gas
      hmod%alpha = 1.0      
      hmod%alphap = 0.0
      hmod%alphaz = 0.0
      hmod%alpha_array = hmod%alpha
      hmod%alphap_array = hmod%alphap
      hmod%alphaz_array = hmod%alphaz
      hmod%A_alpha = -0.005
      hmod%B_alpha = 0.022
      hmod%C_alpha = 0.865
      hmod%D_alpha = -13.565
      hmod%E_alpha = 52.516

      ! Non-virial temperature correction for hot gas
      hmod%beta = 1.0
      hmod%betap = 0.0
      hmod%betaz = 0.0

      ! Low-gas-content concentration modification
      hmod%eps = 0.0
      hmod%epsz = 0.0
      hmod%eps_array = hmod%eps
      hmod%epsz_array = hmod%epsz
      hmod%A_eps = -0.289
      hmod%B_eps = 2.147
      hmod%C_eps = 0.129
      hmod%D_eps = -0.867

      ! High-gas-content concentration modification
      hmod%eps2 = 0.        
      hmod%eps2z = 0.
      hmod%eps2_array = hmod%eps2
      hmod%eps2z_array = hmod%eps2z 

      ! Polytropic gas index
      hmod%Gamma = 1.17
      hmod%Gammap = 0.0
      hmod%Gammaz = 0.
      hmod%Gamma_array = hmod%Gamma
      hmod%Gammap_array = hmod%Gammap
      hmod%Gammaz_array = hmod%Gammaz
      hmod%A_Gamma = 0.026
      hmod%B_Gamma = -0.064
      hmod%C_Gamma = 1.150
      hmod%D_Gamma = -17.011
      hmod%E_Gamma = 66.289

      ! Pressure polytropic gas index
      hmod%Zamma = hmod%Gamma
      hmod%Zammap = hmod%Gammap
      hmod%Zammaz = hmod%Gammaz
      hmod%Zamma_array = hmod%Zamma
      hmod%Zammap_array = hmod%Zammap
      hmod%Zammaz_array = hmod%Zammaz

      ! Halo mass that has lost half gas
      hmod%M0 = 1e14     
      hmod%M0z = 0.0   
      hmod%M0_array = hmod%M0
      hmod%M0z_array = hmod%M0z
      hmod%A_M0 = -0.007
      hmod%B_M0 = 0.018
      hmod%C_M0 = 6.788
      hmod%D_M0 = -103.453
      hmod%E_M0 = 406.705

      ! Maximum star-formation efficiency
      hmod%Astar = 0.03
      hmod%Astarz = 0.0
      hmod%Astar_array = hmod%Astar
      hmod%Astarz_array = hmod%Astarz
      hmod%A_Astar = 0.004
      hmod%B_Astar = -0.034
      hmod%C_Astar = -0.017
      hmod%D_Astar = 0.165

      ! WHIM temperature [K]
      hmod%Twhim = 10**6.5
      hmod%Twhimz = 0.0
      hmod%Twhim_array = hmod%Twhim
      hmod%Twhimz_array = hmod%Twhimz
      hmod%A_Twhim = -0.024
      hmod%B_Twhim = 0.077
      hmod%C_Twhim = 0.454
      hmod%D_Twhim = 1.717

      ! Stellar concentration c_* = rv/r_*
      hmod%cstar = 10.      
      hmod%cstarp = 0.0
      hmod%cstarz = 0.0
      hmod%cstar_array = hmod%cstar
      hmod%cstarz_array = hmod%cstarz

      ! sigma_* for f_* distribution
      hmod%sstar = 1.2      

      ! M* for most efficient halo mass for star formation
      hmod%Mstar = 10**12.5
      hmod%Mstarz = 0.0
      hmod%Mstar_array = hmod%Mstar
      hmod%Mstarz_array = hmod%Mstarz

      ! Fraction of bound gas that is cold
      hmod%fcold = 0.0     

      ! Fraction of bound gas that is hot
      hmod%fhot = 0.0       

      ! Power-law for central galaxy mass fraction
      hmod%eta = 0.0
      hmod%etaz = 0.0
      hmod%eta_array = hmod%eta
      hmod%etaz_array = hmod%etaz

      ! Isothermal beta power index
      hmod%ibeta = 2./3.
      hmod%ibetap = 0.0
      hmod%ibetaz = 0.0

      ! Gas fraction power-law
      hmod%gbeta = 0.6
      hmod%gbetaz = 0.0

      !!! !!!

      ! Do we treat the halomodel as a response model (multiply by HMcode) or not?
      ! 0 - No, otherwise choose an ihm for the baseline and denominator
      ! If the denominator remains zero then it is taken to be the DMONLY model from the current ihm
      hmod%response_baseline = 0
      hmod%response_denominator = 0

      ! Do any response calculations apply to the matter spectra only, or to all spectra?
      hmod%response_matter_only = .FALSE.

      ! Halo mass if the mass function is a delta function
      hmod%hmass = 1e13

      ! Scatter in c/c
      hmod%sigma_conc = 0.

      ! HOD parameters
      hmod%mhalo_min = 1e12
      hmod%mhalo_max = 1e16

      ! HI parameters
      hmod%HImin = 1e9
      hmod%HImax = 1e12

      ! NFW core radius
      hmod%rcore = 0.1

      ! Sheth & Tormen (1999) mass function parameters
      hmod%ST_p = 0.3
      hmod%ST_q = 0.707 ! 0.75 is proposed in Sheth & Tormen (2002)

      ! Mass function amplitude (not very physical)
      hmod%Amf = 1.0    

      ! Index to make the mass function integration easier
      ! Should be related to how the mass function diverges at low nu
      ! e.g., ST diverges like nu^-0.6, alpha = 1/(1-0.6) = 2.5
      hmod%alpha_numu = 2.5

      ! Set flags to false
      hmod%has_galaxies = .FALSE.
      hmod%has_HI = .FALSE.
      hmod%has_mass_conversions = .FALSE.
      hmod%has_mass_function = .FALSE.
      hmod%has_bnl = .FALSE.

      ! Saturation (e.g., WDM)
      hmod%saturation = .FALSE.
      hmod%nu_saturation = 0.

      ! Non-linear halo bias
      hmod%iextrap_bnl = iextrap_bnl
      hmod%stitch_bnl_nu = .FALSE.
      hmod%bnl_cat = 'rockstar'
      hmod%bnl_path = '/Users/Mead/Physics/Multidark/data/'

      ! Perturbation(ish) theory
      hmod%PT_A = 1.
      hmod%PT_alpha = 1.
      hmod%PT_beta = 1.

      IF (ihm == -1) THEN
         WRITE (*, *) 'ASSIGN_HALOMOD: Choose your halo model'
         DO i = 1, nhalomod
            IF(trim(names(i)) == '') EXIT
            WRITE (*, *) i, trim(names(i))
         END DO
         READ (*, *) ihm
         WRITE (*, *)
      END IF

      IF (is_in_array(ihm, [1, 7, 15, 28, 31, 50, 51, 53, 64, 66, 92, 95])) THEN
         !  1 - HMcode (2016)
         !  7 - HMcode (2015)
         ! 15 - HMcode (2019)
         ! 28 - HMcode (2016) w/ one parameter baryon model
         ! 31 - HMcode (2016) w/ additional BAO damping
         ! 50 - HMcode (2016) w/ pow=1 bug in Dolag that was present in CAMB
         ! 51 - HMcode (2016) w/ CAMB mass range points
         ! 53 - HMcode (2018)
         ! 64 - HMcode (2016) w/ 2020 baryon recipe
         ! 66 - HMcode (2015) w/ CAMB mass range and points
         ! 92 - HMcode (2016) w/ neutrino halo-mass fix
         ! 95 - HMcode (2016) w/ neutrino halo-mass fix and CAMB mass range and points
         hmod%ip2h = 1
         hmod%i1hdamp = 1
         hmod%iconc = 1
         hmod%idc = 3
         hmod%iDv = 3
         hmod%ieta = 1
         hmod%iAs = 1
         hmod%i2hdamp = 2
         hmod%itrans = 1
         hmod%iDolag = 2
         hmod%zD = 10.
         hmod%flag_sigma = flag_ucold
         hmod%DMONLY_neutrino_halo_mass_correction = .FALSE.
         hmod%Dv0 = 418.
         hmod%Dv1 = -0.352
         hmod%dc0 = 1.59
         hmod%dc1 = 0.0314
         hmod%eta0 = 0.603
         hmod%eta1 = 0.300
         hmod%f0 = 0.0095
         hmod%f1 = 1.37
         hmod%ks = 0.584
         hmod%As = 3.13
         hmod%alp0 = 3.24
         hmod%alp1 = 1.85
         hmod%Dvnu = 0.916
         hmod%dcnu = 0.262
         IF (ihm == 7 .OR. ihm == 66) THEN
            !  7 - HMcode (2015)
            ! 66 - HMcode (2015) but with CAMB halo mass range and number of mass points
            hmod%idc = 6
            hmod%i2hdamp = 1
            hmod%itrans = 1
            hmod%iDolag = 1
            hmod%f0 = 0.188
            hmod%f1 = 4.29
            hmod%alp0 = 2.93
            hmod%alp1 = 1.77       
            hmod%Dvnu = 0.
            hmod%dcnu = 0.
            hmod%flag_sigma = flag_matter
         ELSE IF (ihm == 15) THEN
            ! HMcode (2019)
            hmod%i1hdamp = 2 ! 2 - k^4 at large scales for one-halo term
            hmod%ip2h = 3    ! 3 - Linear theory with damped wiggles
            hmod%i2hdamp = 1 ! 1 - Change back to Mead (2015) model for two-halo damping
            hmod%DMONLY_neutrino_halo_mass_correction = .TRUE. ! Neutrino halo-mass fix
            ! Model 1: Fits at 0.0196 to FrankenEmu nodes
            hmod%Dv0 = 411.
            hmod%Dv1 = -0.333
            hmod%dc0 = 1.631
            hmod%dc1 = 0.0187
            hmod%eta0 = 0.523
            hmod%eta1 = 0.259
            hmod%f0 = 0.0615
            hmod%f1 = 1.607
            hmod%ks = 0.136
            hmod%As = 3.23
            hmod%alp0 = 3.17
            hmod%alp1 = 1.88
            hmod%Dvnu = 0.
            hmod%dcnu = 0.
         ELSE IF (ihm == 28) THEN
            ! HMcode (2016) with one-parameter baryon model
            hmod%one_parameter_baryons = .TRUE.
            !hmod%As = 3.13 ! Default value
            hmod%As = 2.00 ! AGN-ish
         ELSE IF (ihm == 31) THEN
            ! HMcode (2016) with damped BAO
            hmod%ip2h = 3
         ELSE IF (ihm == 50) THEN
            ! HMcode (2016) with bug present in CAMB when Dolag power is set to 1 rather than 1.5
            hmod%iDolag = 1
         ELSE IF (ihm == 53) THEN
            ! HMcode (2018) with improved fitted parameters and neutrino halo-mass fix compared to 2016
            hmod%Dv0 = 416.
            hmod%Dv1 = -0.346
            hmod%dc0 = 1.661
            hmod%dc1 = 0.0131
            hmod%eta0 = 0.501
            hmod%eta1 = 0.273
            hmod%f0 = 0.00612
            hmod%f1 = 1.607
            hmod%ks = 0.927
            hmod%As = 3.09
            hmod%alp0 = 3.11
            hmod%alp1 = 1.91
            hmod%Dvnu = 0. ! Set to zero because not necessary after halo-mass fix
            hmod%dcnu = 0. ! Set to zero because not necessary after halo-mass fix
            hmod%DMONLY_neutrino_halo_mass_correction = .TRUE.
         ELSE IF (ihm == 64) THEN
            ! 64 - HMcode 2016 but with 2020 baryon recipe
            hmod%DMONLY_baryon_recipe = .TRUE.
         ELSE IF (ihm == 92 .OR. ihm == 95) THEN
            ! 92, 95 - HMcode (2016) with neutrino halo-mass fix
            hmod%DMONLY_neutrino_halo_mass_correction = .TRUE. ! Neutrino halo-mass fix
         END IF
         IF (ihm == 51 .OR. ihm == 66 .OR. ihm == 95) THEN
            ! CAMB mass range
            hmod%mmin = 1e0
            hmod%mmax = 1e18
            hmod%n = 256
         END IF
      ELSE IF (ihm == 2) THEN
         ! Basic halo model with linear two halo term (Delta_v = 200, delta_c = 1.686)
         ! This is what the original HMcode model was based off
         ! Note that Delta_v = 200 is not consistent with the ST mass function
         hmod%ip2h = 1
         hmod%idc = 1
         hmod%iDv = 1
         hmod%iconc = 1
      ELSE IF (ihm == 3) THEN
         ! Standard halo-model calculation (Seljak 2000)
         ! This is the default, so do nothing here
      ELSE IF (ihm == 4) THEN
         ! Standard halo-model calculation but with HMcode (2015) smoothed two- to one-halo transition and one-halo damping
         hmod%itrans = 2
         hmod%i1hdamp = 2
         hmod%safe_negative = .TRUE.
      ELSE IF (ihm == 5) THEN
         ! Standard halo-model calculation but with Delta_v = 200 and delta_c = 1.686 fixed and Bullock c(M)
         hmod%idc = 1
         hmod%iDv = 1
         hmod%iconc = 1
      ELSE IF (ihm == 6) THEN
         ! Half-accurate halo-model calculation, inspired by (HMcode 2015, 2016)
         hmod%i1hdamp = 2
         hmod%i2hdamp = 2
         hmod%safe_negative = .TRUE.
         hmod%idc = 3
         hmod%iDv = 3
         hmod%ieta = 1
         hmod%iAs = 1
         hmod%itrans = 1
         hmod%iconc = 1
         hmod%iDolag = 2
      ELSE IF (ihm == 8) THEN
         ! Include scatter in halo properties
         hmod%sigma_conc = 0.25
      ELSE IF (ihm == 9) THEN
         ! For CCL comparison and benchmark generation
         hmod%n = 2048   ! Increase accuracy for the CCL benchmarks
         hmod%acc = 1e-5 ! Increase accuracy for the CCL benchmarks
         hmod%ip2h = 2
         hmod%i2hcor = 2
         hmod%ibias = 1
         hmod%i1hdamp = 0
         hmod%imf = 2
         hmod%iconc = 4 ! 4 - Virial Duffy relation for full sample
         hmod%idc = 2   ! 2 - Virial dc
         hmod%iDv = 2   ! 2 - Virial Dv
         hmod%ieta = 0  ! 0 - No
         hmod%iAs = 0
         hmod%i2hdamp = 0
         hmod%itrans = 0
         hmod%iDolag = 0
      ELSE IF (ihm == 10) THEN
         ! For mass conversions comparison with Wayne Hu's code
         hmod%iconc = 2
         hmod%idc = 1
         hmod%iDv = 1
      ELSE IF (ihm == 11) THEN
         ! UPP
         hmod%electron_pressure = 1
      ELSE IF (ihm == 12) THEN
         ! Spherical-collapse model to produce Mead (2017) results
         hmod%idc = 5    ! Spherical-collapse calculation for delta_c
         hmod%iDv = 5    ! Spherical-collapse calculation for Delta_v
         hmod%imf = 2    ! Sheth & Tormen mass function
         hmod%iconc = 1  ! Bullock et al. c(M) relation
         hmod%iDolag = 1 ! This is important for the accuracy of the z=0 results presented in Mead (2017)
      ELSE IF (ihm == 13) THEN
         ! Experimental sigmoid transition
         hmod%itrans = 4
      ELSE IF (ihm == 14) THEN
         ! Experimental scale-dependent halo bias
         hmod%ibias = 5
         hmod%i1hdamp = 2
      ELSE IF (ihm == 16) THEN
         ! Halo-void model
         hmod%add_voids = .TRUE.
      ELSE IF (ihm == 17 .OR. ihm == 18 .OR. ihm == 19) THEN
         ! 17 - HMx AGN 7.6
         ! 18 - HMx AGN tuned
         ! 19 - HMx AGN 8.0
         hmod%HMx_mode = 4
         hmod%itrans = 3
         !hmod%itrans = 1
         hmod%i1hdamp = 2
         hmod%safe_negative = .TRUE.
         hmod%response_baseline = HMcode2016
         hmod%alp0 = 3.24
         hmod%alp1 = 1.85
         STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
      ELSE IF (ihm == 20) THEN
         ! Standard halo model but as response with HMcode
         hmod%response_baseline = HMcode2016
      ELSE IF (ihm == 21) THEN
         ! Cored NFW halo profile model
         hmod%halo_DMONLY = 5 ! Cored profile
      ELSE IF (ihm == 22) THEN
         ! Different stellar profile
         hmod%halo_central_stars = 2 ! Schneider & Teyssier (2015)
         hmod%response_baseline = HMcode2016
      ELSE IF (ihm == 23) THEN
         ! Tinker (2010) mass function and bias; virial halo mass
         hmod%imf = 3 ! Tinker mass function and bias
      ELSE IF (ihm == 25) THEN
         ! Villaescusa-Navarro HI halo model
         hmod%imf = 3     ! Tinker mass function
         hmod%frac_HI = 3 ! HI mass fraction with z evolution (Villaescusa-Navarro et al. 1804.09180)
         !hmod%halo_HI = 3 ! Exponentially cored polynomial (Villaescusa-Navarro et al. 1804.09180)
         !hmod%halo_HI = 2 ! Delta function
         hmod%halo_HI = 5 ! Modified NFW (Padmanabhan & Refregier 2017; 1607.01021)
         hmod%ibias = 3   ! Non-linear halo bias
      ELSE IF (ihm == 26) THEN
         ! Delta function mass function
         hmod%ip2h = 2        ! Standard two-halo term
         hmod%iDv = 6         ! Lagrangian radius haloes
         hmod%imf = 4         ! Delta function mass function
         hmod%halo_DMONLY = 3 ! Top-hat halo profile
      ELSE IF (ihm == 27) THEN
         hmod%imf = 1 ! Press & Schecter (1974) mass function
      ELSE IF (ihm == 29) THEN
         ! Adding some cold gas
         hmod%fcold = 0.1
      ELSE IF (ihm == 30) THEN
         ! Adding some hot gas
         hmod%fhot = 0.2
      ELSE IF (ihm == 32 .OR. ihm == 33 .OR. ihm == 34 .OR. ihm == 35 .OR. ihm == 36 .OR. ihm == 37) THEN
         ! 32 - Response model for AGN 7.6;   z = 0.0; simple pivot
         ! 33 - Response model for AGN tuned; z = 0.0; simple pivot
         ! 34 - Response model for AGN 8.0;   z = 0.0; simple pivot
         ! 35 - Response model for AGN 7.6;   z = 0.0; clever pivot; f_hot
         ! 36 - Response model for AGN tuned; z = 0.0; clever pivot; f_hot
         ! 37 - Response model for AGN 8.0;   z = 0.0; clever pivot; f_hot
         hmod%response_baseline = HMcode2016
         IF (ihm == 32) THEN
            ! AGN 7.6
            ! Mpiv = 1e14; z = 0.0
            hmod%simple_pivot = .TRUE.
            hmod%alpha = 0.709
            hmod%eps = 0.12
            hmod%Gamma = 1.236
            hmod%M0 = 10.**13.6
            hmod%Astar = 0.045
            hmod%Twhim = 10.**6.07
            hmod%cstar = 11.2
            hmod%fcold = 0.021
            hmod%Mstar = 10.**11.5
            hmod%sstar = 0.85
            hmod%alphap = -0.479
            hmod%Gammap = -0.020
            hmod%cstarp = -0.459
         ELSE IF (ihm == 33) THEN
            ! AGN tuned
            ! Mpiv = 1e14; z = 0.0
            hmod%simple_pivot = .TRUE.
            hmod%alpha = 0.802
            hmod%eps = 0.06
            hmod%Gamma = 1.268
            hmod%M0 = 10.**13.9
            hmod%Astar = 0.043
            hmod%Twhim = 10.**6.09
            hmod%cstar = 12.3
            hmod%fcold = 0.019
            hmod%Mstar = 10.**10.6
            hmod%sstar = 1.00
            hmod%alphap = -0.518
            hmod%Gammap = -0.033
            hmod%cstarp = -0.523
         ELSE IF (ihm == 34) THEN
            ! AGN 8.0
            ! Mpiv = 1e14; z = 0.0
            hmod%simple_pivot = .TRUE.
            hmod%alpha = 0.769
            hmod%eps = -0.09
            hmod%Gamma = 1.274
            hmod%M0 = 10.**14.3
            hmod%Astar = 0.039
            hmod%Twhim = 10.**6.11
            hmod%cstar = 12.7
            hmod%fcold = 0.012
            hmod%Mstar = 10.**10.
            hmod%sstar = 1.16
            hmod%alphap = -0.409
            hmod%Gammap = -0.029
            hmod%cstarp = -0.510
         ELSE IF (ihm == 35) THEN
            ! AGN 7.6
            ! Mpiv = Mh; z = 0.0; f_hot
            hmod%alpha = 1.247
            hmod%eps = 0.087
            hmod%Gamma = 1.253
            hmod%M0 = 10.**13.63
            hmod%Astar = 0.042
            hmod%Twhim = 10.**6.08
            hmod%cstar = 10.80
            hmod%fcold = 0.0153
            hmod%Mstar = 10.**12.12
            hmod%sstar = 0.916
            hmod%alphap = -0.488
            hmod%Gammap = -0.0164
            hmod%cstarp = -0.222
            hmod%fhot = 0.0323
         ELSE IF (ihm == 36) THEN
            ! AGN tuned
            ! Mpiv = Mh; z = 0.0; f_hot
            hmod%alpha = 1.251
            hmod%eps = -0.291
            hmod%Gamma = 1.257
            hmod%M0 = 10.**13.856
            hmod%Astar = 0.0413
            hmod%Twhim = 10.**6.09
            hmod%cstar = 11.77
            hmod%fcold = 0.0066
            hmod%Mstar = 10.**12.12
            hmod%sstar = 0.856
            hmod%alphap = -0.481
            hmod%Gammap = -0.019
            hmod%cstarp = -0.303
            hmod%fhot = 0.0282
         ELSE IF (ihm == 37) THEN
            ! AGN 8.0
            ! Mpiv = Mh; z = 0.0; f_hot
            hmod%alpha = 1.024
            hmod%eps = -0.179
            hmod%Gamma = 1.264
            hmod%M0 = 10.**14.33
            hmod%Astar = 0.0383
            hmod%Twhim = 10.**6.10
            hmod%cstar = 17.09
            hmod%fcold = 0.00267
            hmod%Mstar = 10.**11.52
            hmod%sstar = 0.99
            hmod%alphap = -0.342
            hmod%Gammap = -0.018
            hmod%cstarp = -0.413
            hmod%fhot = 0.0031
         ELSE
            STOP 'ASSIGN_HALOMOD: Error, ihm specified incorrectly'
         END IF
         hmod%beta = hmod%alpha
         hmod%betap = hmod%alphap
      ELSE IF (ihm == 38 .OR. ihm == 39 .OR. ihm == 40) THEN
         ! 38 - AGN 7p6
         ! 39 - AGN tuned
         ! 49 - AGN 8p0
         hmod%response_baseline = HMcode2016
         IF (ihm == 38) THEN
            ! AGN 7p6
            hmod%alpha = 1.52016437
            hmod%eps = 0.06684244
            hmod%Gamma = 1.24494147
            hmod%M0 = 2.84777660E+13
            hmod%Astar = 3.79242003E-02
            hmod%Twhim = 987513.562
            hmod%cstar = 9.78754425
            hmod%fcold = 1.83899999E-02
            hmod%mstar = 2.68670049E+12
            hmod%sstar = 1.13123488
            hmod%alphap = -0.531507611
            hmod%Gammap = -6.00000005E-03
            hmod%cstarp = -0.113370098
            hmod%fhot = 1.08200002E-04
            hmod%alphaz = 0.346108794
            hmod%Gammaz = 0.231210202
            hmod%M0z = -9.32227001E-02
            hmod%Astarz = -0.536369383
            hmod%Twhimz = -0.139042795
            hmod%eta = -0.222550094
         ELSE IF (ihm == 39) THEN
            ! AGN tuned
            hmod%alpha = 1.54074240
            hmod%eps = 0.01597583
            hmod%Gamma = 1.24264216
            hmod%M0 = 6.51173707E+13
            hmod%Astar = 3.45238000E-02
            hmod%Twhim = 1389362.50
            hmod%cstar = 10.4484358
            hmod%fcold = 6.32560020E-03
            hmod%mstar = 2.34850982E+12
            hmod%sstar = 1.25916803
            hmod%alphap = -0.513900995
            hmod%Gammap = -5.66239981E-03
            hmod%cstarp = -6.11630008E-02
            hmod%fhot = 1.26100000E-04
            hmod%alphaz = 0.340969592
            hmod%Gammaz = 0.295596898
            hmod%M0z = -8.13719034E-02
            hmod%Astarz = -0.545276821
            hmod%Twhimz = -0.122411400
            hmod%eta = -0.266318709
         ELSE IF (ihm == 40) THEN
            ! AGN 8p0
            hmod%alpha = 1.45703220
            hmod%eps = -0.128
            hmod%Gamma = 1.24960959
            hmod%M0 = 1.15050950E+14
            hmod%Astar = 3.86818014E-02
            hmod%Twhim = 1619561.00
            hmod%cstar = 19.4119701
            hmod%fcold = 1.10999999E-05
            hmod%mstar = 2.18769510E+12
            hmod%sstar = 0.803893507
            hmod%alphap = -0.528370678
            hmod%Gammap = -3.31420009E-03
            hmod%cstarp = -0.355121315
            hmod%fhot = 3.90500005E-04
            hmod%alphaz = 0.740169585
            hmod%Gammaz = 0.354409009
            hmod%M0z = -2.40819994E-02
            hmod%Astarz = -0.425019890
            hmod%Twhimz = -8.60318989E-02
            hmod%eta = -0.243649304
         ELSE
            STOP 'ASSIGN_HALOMOD: Error, ihm specified incorrectly'
         END IF
         hmod%beta = hmod%alpha
         hmod%betap = hmod%alphap
         hmod%betaz = hmod%alphaz
      ELSE IF (ihm == 41) THEN
         ! Some stellar mass in satellite galaxies
         hmod%eta = -0.3
      ELSE IF (ihm == 42) THEN
         ! Tinker and stuff apprpriate for M200c
         hmod%imf = 3   ! Tinker mass function and bias
         hmod%iDv = 7   ! M200c
         hmod%iconc = 5 ! Duffy for M200c for full sample
         hmod%idc = 1   ! Fixed to 1.686
      ELSE IF (ihm == 43) THEN
         ! Standard halo model but as response with HMcode (2016) but only for matter spectra
         hmod%response_baseline = HMcode2016
         hmod%response_matter_only = .TRUE.
      ELSE IF (ihm == 44) THEN
         ! Tinker and stuff apprpriate for M200
         hmod%imf = 3   ! Tinker mass function and bias
         hmod%iDv = 1   ! M200
         hmod%iconc = 3 ! Duffy for M200 for full sample
         hmod%idc = 1   ! Fixed to 1.686
      ELSE IF (ihm == 45) THEN
         ! No stars
         hmod%Astar = 0.
      ELSE IF (ihm == 46) THEN
         ! Isothermal beta model
         hmod%halo_normal_bound_gas = 2
      ELSE IF (ihm == 47) THEN
         ! Isothermal beta model, response
         hmod%halo_normal_bound_gas = 2
         hmod%response_baseline = HMcode2016
      ELSE IF (is_in_array(ihm, [24, 48, 49, 107, 114, 116, 117, 118])) THEN
         ! Non-linear halo bias
         !  24 - Tinker 2010 with M200c haloes
         !  48 - Sheth & Tormen (1999)
         !  49 - Tinker 2010 mass function with peak-background split halo bias
         ! 107 - Tinker 2010 mass function with peak-background split halo bias and no extrapolation
         ! 114 - DEFAULT: Tinker 2010
         ! 116 - Tinker 2010 and no extrapolation
         ! 117 - Tinker 2010 and adding Bolshoi for low-mass haloes
         ! 118 - Tinker 2010 and adding Bolshoi for low-mass haloes and no extrapolation
         hmod%ibias = 3 ! Non-linear halo bias
         !hmod%i1hdamp = 3 ! One-halo damping like k^4         
         IF (ihm == 48) THEN
            hmod%imf = 2 ! Sheth & Tormen (1999)
         ELSE IF (ihm == 49 .OR. ihm == 107) THEN
            hmod%imf = 3  ! Tinker 2010 with peak-background split bias
         ELSE
            hmod%imf = 17 ! Tinker 2010 mass function and calibrated bias
         END IF  
         IF (ihm == 24) THEN
            hmod%iDv = 7   ! M200c
            hmod%iconc = 5 ! Duffy M200c concentrations for full sample
         END IF
         IF (ihm == 117 .OR. ihm == 118) THEN
            hmod%stitch_bnl_nu = .TRUE.
            hmod%bnl_cat = 'BDMV'
         END IF
         IF (ihm == 107 .OR. ihm == 116 .OR. ihm == 118) hmod%iextrap_bnl = iextrap_zero
      ELSE IF(ihm == 52) THEN
         ! Standard halo model but with Mead (2017) spherical-collapse fitting function
         hmod%idc = 4 ! Mead (2017) fitting function for delta_c
         hmod%iDv = 4 ! Mead (2017) fitting function for Delta_v
      ELSE IF (ihm == 54) THEN
         ! Such that hydro model masquerades as DMONLY
         hmod%frac_bound_gas = 3     ! Universal baryon fraction as bound gas
         hmod%halo_normal_bound_gas = 4    ! NFW profile for bound gas
         hmod%Astar = 0.             ! No stars
         hmod%frac_central_stars = 1 ! All stars are central stars (not necessary, but maybe speeds up)
         hmod%frac_stars = 2         ! Constant star fraction (not necessary, but maybe speeds up)
      ELSE IF (is_in_array(ihm, [55, 56, 57, 58, 59, 60, 61, 62, 63, 65, 81, 82, 83, 84, 85, 86])) THEN
         ! HMx2020: Baseline
         hmod%response_baseline = HMcode2016 ! Model should be calculated as a response
         hmod%halo_central_stars = 3   ! 3 - Delta function for central stars
         hmod%halo_satellite_stars = 1 ! 1 - NFW
         hmod%eta = -0.3               ! eta to split central and satellitee galaxies
         hmod%HMx_mode = 5             ! HMx2020 possible M and z dependence of parameters
         IF (ihm == 65) hmod%response_baseline = 0 ! No response
         IF (is_in_array(ihm, [81, 82, 83])) hmod%response_matter_only = .TRUE.
         IF (is_in_array(ihm, [84, 85, 86])) hmod%different_Gammas = .TRUE.
         IF (ihm == 56 .OR. ihm == 81 .OR. ihm == 84) THEN
            ! AGN 7.6
            hmod%fix_star_concentration = .TRUE.
            hmod%Astar = 0.03477
            hmod%mstar = 10**12.46201
            hmod%Astarz = -0.00926
            hmod%eta = -0.34285
            hmod%mstarz = -0.36641
         ELSE IF (ihm == 57 .OR. ihm == 82 .OR. ihm == 85) THEN
            ! AGN 7.8
            hmod%fix_star_concentration = .TRUE.
            hmod%Astar = 0.03302
            hmod%mstar = 10**12.44789
            hmod%Astarz = -0.00881
            hmod%eta = -0.35556
            hmod%mstarz = -0.35206
         ELSE IF (ihm == 58 .OR. ihm == 83 .OR. ihm == 86) THEN
            ! AGN 8.0
            hmod%fix_star_concentration = .TRUE.
            hmod%Astar = 0.03093
            hmod%mstar = 10**12.39230
            hmod%Astarz = -0.00818
            hmod%eta = -0.35052
            hmod%mstarz = -0.30727
         ELSE IF (ihm == 59 &
                  .OR. ihm == HMx2020_matter_w_temp_scaling &
                  .OR. ihm == HMx2020_matter_pressure_w_temp_scaling) THEN      
            ! 59 - HMx2020 with temperature scaling that fits stars
            ! 60 - HMx2020 with temperature scaling that fits matter (stars fixed)
            ! 61 - HMx2020 with temperature scaling that fits matter, pressure (stars fixed)
            hmod%fix_star_concentration = .TRUE.
            hmod%HMx_mode = 6 ! Scaling with temperature
            hmod%Astar_array = [0.0348, 0.0330, 0.0309]          
            hmod%Mstar_array = 10**[12.4620, 12.4479, 12.3923]
            hmod%Astarz_array = [-0.0093, -0.0088, -0.0082]
            hmod%eta_array = [-0.3428, -0.3556, -0.3505]
            hmod%Mstarz_array = [-0.3664, -0.3521, -0.3073]           
            IF(ihm == HMx2020_matter_w_temp_scaling) THEN
               ! 60 - Additionally fits matter-matter
               hmod%eps_array = [0.2841, 0.2038, 0.0526]
               hmod%Gamma_array = 1.+[0.2363, 0.3376, 0.6237]
               hmod%M0_array = 10**[13.0020, 13.3658, 14.0226]
               hmod%epsz_array = [-0.0046, -0.0047, 0.0365]
            END IF
            IF(ihm == HMx2020_matter_pressure_w_temp_scaling) THEN
               ! 61 - Additionally fits matter, pressure
               hmod%alpha_array = [0.76425, 0.84710, 1.03136]
               hmod%eps_array = [-0.10017, -0.10650, -0.12533]
               hmod%Gamma_array = 1.+[0.16468, 0.17702, 0.19657]
               hmod%M0_array = 10**[13.19486, 13.59369, 14.24798]
               hmod%Twhim_array = 10**[6.67618, 6.65445, 6.66146] ! Weirdly similar (z=0 Twhim all the same... ?)
               hmod%Twhimz_array = [-0.55659, -0.36515, -0.06167]
               hmod%epsz_array = [-0.04559, -0.10730, -0.01107] ! Non monotonic
            END IF
         ELSE IF (ihm == 62) THEN
            ! 62 - Model for matter, CDM, gas, stars
            hmod%HMx_mode = 6 ! Scaling with temperature
            hmod%eps_array = [0.40207, 0.12360, -0.11578]
            hmod%Gamma_array = 1.+[0.27626, 0.29555, 0.28608] ! Weirdly similar, non monotonic
            hmod%M0_array = 10**[13.09777, 13.48537, 14.12539]
            hmod%Astar_array = [0.03460, 0.03417, 0.03205]
            hmod%Mstar_array = 10**[12.55064, 12.37147, 12.30323]
            hmod%Gammaz_array = [-0.05539, -0.09365, -0.13816]
            hmod%Astarz_array = [-0.00917, -0.01046, -0.00935] ! Non monotonic
            hmod%eta_array = [-0.49697, -0.40520, -0.34431]
            hmod%epsz_array = [0.04350, -0.01872, 0.14079] ! Non monotonic
            hmod%Mstarz_array = [-0.46153, 0.01489, -0.08174] ! Non monotonic
         ELSE IF (ihm == 63) THEN  
            ! 62 - Model for matter, CDM, gas, stars, electron pressure
            ! TODO: Not applied becuase it did not work
            hmod%HMx_mode = 6 ! Scaling with temperature
         END IF
      ELSE IF (ihm == 67) THEN
         ! Standard halo-model calculation but with Delta_v = 200 and delta_c = 1.686
         hmod%idc = 1
         hmod%iDv = 1
      ELSE IF (ihm == 68) THEN
         ! Standard halo-model calculation but with Bullock c(M)
         hmod%iconc = 1
      ELSE IF (ihm == 69) THEN
         ! Standard halo-model calculation but with simple Bullock c(M)
         hmod%iconc = 2
      ELSE IF (ihm == 70) THEN
         ! Standard but with no Dolag correction
         hmod%iDolag = 0
      ELSE IF (ihm == 71) THEN
         ! Standard but with Dolag correction with 1.5 exponent
         hmod%iDolag = 2
      ELSE IF (ihm == 72) THEN
         ! Standard but with linear two-halo term
         hmod%ip2h = 1
      ELSE IF (ihm == 73) THEN
         ! Standard but with two-halo term with damped wiggles
         hmod%ip2h = 3
      ELSE IF (ihm == 74) THEN
         ! Standard but with dc = 1.686 fixed
         hmod%idc = 1
      ELSE IF (ihm == 75) THEN
         ! Standard but with dc from Mead (2017) fit
         hmod%idc = 4
      ELSE IF (ihm == 76) THEN
         ! Standard but with Dv from Mead (2017) fit
         hmod%iDv = 4
      ELSE IF (is_in_array(ihm, [77, 78, 79, 102, 103, 104, 123, 124, 125])) THEN
         ! HMcode (2020)
         !  77 - Unfitted
         !  78 - Fitted to Cosmic Emu for k<1
         !  79 - Fitted and published model
         ! 102 - Fitted with baryon model
         ! 103 - Baryon response model
         ! 104 - Baryon model using HMx language
         ! 123 - Extended mass range
         ! 124 - Baryon response model and extended mass range
         ! 125 - Unfitted with extended mass range
         hmod%ip2h = 3    ! 3 - Linear two-halo term with damped wiggles
         hmod%i1hdamp = 3 ! 3 - k^4 at large scales for one-halo term
         hmod%itrans = 1  ! 1 - HMcode alpha-neff smoothing
         hmod%i2hdamp = 3 ! 3 - fdamp for perturbation theory      
         hmod%idc = 4     ! 4 - delta_c from Mead (2017) fit
         hmod%iDv = 4     ! 4 - Delta_v from Mead (2017) fit
         hmod%iconc = 1   ! 1 - Bullock c(M) relation   
         hmod%iDolag = 3  ! 3 - Dolag c(M) correction with sensible z evolution
         hmod%iAs = 2     ! 2 - Vary c(M) relation prefactor with sigma8 dependence
         hmod%ieta = 3    ! 3 - HMcode 2020 eta bloating
         hmod%zD = 10.    ! 10 vs 100 makes a difference for EDE-type cosmologies
         hmod%flag_sigma = flag_ucold ! Cold un-normalised produces better massive-neutrino results
         hmod%DMONLY_neutrino_halo_mass_correction = .TRUE. ! Correct haloes for missing neutrino mass
         IF (ihm == 123 .OR. ihm == 124 .OR. ihm == 125) THEN
            hmod%mmin = 1e0  ! Reduced lower-mass limit
            hmod%mmax = 1e18 ! Increased upper-mass limit
            hmod%n = 256     ! Increase number of points in mass
         END IF
         IF (ihm == 78) THEN
            ! Model 3: 0.00926 for Cosmic Emu
            hmod%f0 = 0.1995332
            hmod%f1 = 0.4271555
            hmod%ks = 0.0300269
            hmod%kp = -1.3693749
            hmod%alp0 = 1.1786775
            hmod%alp1 = 0.2668017
            hmod%alp2 = 0.4473270
            hmod%kd = 0.0414780
            hmod%kdp = -1.0887362
            hmod%nd = 2.7259666
         ELSE IF (ihm == 79 .OR. ihm == 102 .OR. ihm == 123) THEN
            ! Published model: 0.0168 to Franken Emu
            hmod%f0 = 0.2695822
            hmod%f1 = 0.9403087
            hmod%ks = 0.0561778
            hmod%kp = -1.0131066
            hmod%kd = 0.0569871
            hmod%kdp = -1.0890162
            hmod%nd = 2.8534197
            hmod%eta0 = 0.1281210
            hmod%eta1 = -0.3643963
            hmod%alp0 = 1.8751465
            hmod%alp1 = 1.6029913
            hmod%As = 5.1958429
         END IF
         IF (ihm == 102) THEN
            ! 102 - HMcode baryon recipe
            hmod%DMONLY_baryon_recipe = .TRUE.
            hmod%As = 3.08320
            hmod%mbar = 10**13.61510
            hmod%mbarz = -0.178117
            hmod%sbar = 0.02996309
            hmod%sbarz = 0.9240190
            hmod%As_T =  -0.47455
            hmod%mbar_T = 1.830053
            hmod%mbarz_T = -0.202781
            hmod%sbar_T = -0.0050202
            hmod%sbarz_T = -0.12670
         END IF
         IF (ihm == 103 .OR. ihm == 124) THEN
            ! 103 - HMcode response baryon recipe
            ! 124 - As above but with extended mass range
            hmod%DMONLY_baryon_recipe = .TRUE.
            IF (ihm == 103) THEN
               hmod%response_baseline = HMcode2020
               hmod%response_denominator = 77
            ELSE IF (ihm == 124) THEN
               hmod%response_baseline = HMcode2020_CAMB
               hmod%response_denominator = 125
            ELSE
               STOP
            END IF         
            hmod%As = 3.44
            hmod%As_T = -0.496
            hmod%Az = -0.0671
            hmod%Az_T = -0.0371
            hmod%sbar = 0.0201
            hmod%sbar_T = -0.0030
            hmod%sbarz = 0.409
            hmod%sbarz_T = 0.0224
            hmod%mbar = 10**13.87
            hmod%mbar_T = 1.81
            hmod%mbarz = -0.108
            hmod%mbarz_T = 0.195
         END IF
         IF (ihm == 104) THEN
            ! 104 - HMx response baryon model
            hmod%response_baseline = HMcode2020   
            hmod%HMx_mode = 1              ! 1 - HMcode 2020 scalings with z and T_AGN         
            hmod%frac_bound_gas = 2        ! 2 - Schneider & Teyssier 2015 baryon fraction
            hmod%frac_stars = 2            ! 2 - Constant fraction of halo in stars
            hmod%frac_central_stars = 1    ! 1 - All stars are central stars
            hmod%halo_normal_bound_gas = 4 ! 4 - Bound gas is NFW
            hmod%halo_ejected_gas = 7      ! 7 - Ejected gas is smooth
            hmod%halo_central_stars = 3    ! 3 - Central stars are delta function
            hmod%normalise_baryons = 2     ! 2 - Limit baryon reserves for gas
            hmod%eps = -1.+3.42/4.
            hmod%gbeta = 2.
            hmod%M0 = 10**13.83
            hmod%Astar = 0.0231
            hmod%M0z = -0.254
            hmod%Astarz = 1.076
         END IF
      ELSE IF (ihm == 80) THEN
         ! Jenkins mass function (defined for FoF 0.2 haloes)
         hmod%iDv = 9 ! 9 - ~ 178
         hmod%imf = 5
      ELSE IF (ihm == 87) THEN
         ! Despali et al. (2016) mass function
         hmod%imf = 6
      ELSE IF (ihm == 88) THEN
         ! Child et al. (2018) concentration-mass relation
         hmod%iDv = 7 ! M200c definitions
         hmod%imf = 3 ! Tinker (for M200c)
         hmod%iconc = 9 ! Child et al.
      ELSE IF (ihm == 89) THEN
         ! All gas is bound gas
         hmod%frac_bound_gas = 3
      ELSE IF (ihm == 90) THEN
         ! Redshift-dependent Dolag
         hmod%iDolag = 3
      ELSE IF (ihm == 91) THEN
         ! Tinker (2008) mass function
         hmod%imf = 7
      ELSE IF (ihm == 93) THEN
         ! Warren (2006) mass function
         ! Approrpriate for FoF = 0.2 haloes so should change c(M) too
         hmod%imf = 8
      ELSE IF (ihm == 94) THEN
         ! Reed (2007) mass function
         hmod%imf = 9
      ELSE IF (ihm == 96) THEN
         ! Philcox et al. (2020) mass function
         ! Approrpriate for FoF = 0.2 haloes so should change c(M) too
         hmod%imf = 14
      ELSE IF (ihm == 97) THEN
         ! Sigma calculated from total matter power
         hmod%flag_sigma = flag_matter
      ELSE IF (ihm == 98) THEN
         ! No correction for neutrino mass applied to DMONLY haloes
         hmod%DMONLY_neutrino_halo_mass_correction = .FALSE.        
      ELSE IF (ihm == 99) THEN
         ! Cold sigma calculated using 1+delta_c = rho_c/mean_rho_m
         hmod%flag_sigma = flag_cold
      ELSE IF (ihm == 100) THEN
         ! Massive neutrinos affect the calculation of the virial radius
         hmod%DMONLY_neutrinos_affect_virial_radius = .TRUE.
      ELSE IF(ihm == 101) THEN
         ! Standard halo model but with spherical-collapse calculation
         hmod%idc = 5 ! Collapse calculation for delta_c
         hmod%iDv = 5 ! Collapse calculation for Delta_v
      ELSE IF (ihm == 105 .OR. ihm == 106) THEN
         ! Slightly adjusted standard and non-linear bias model
         ! Adjusted to get hydro power large-scale amplitude more correct
         hmod%imf = 17       ! Tinker 2010 mass function and bias         
         hmod%Astar = 0.0320 ! Stellar mass fraction
         IF (ihm == 106) hmod%ibias = 3 ! Non-linear halo bias
      ELSE IF (ihm == 108) THEN
         ! Bhattacharya et al. (2011) mass function
         ! Approrpriate for FoF = 0.2 haloes so should change c(M) too
         hmod%imf = 10
      ELSE IF (ihm == 109 .OR. ihm == 111) THEN 
         ! 109 - PT inspired thing
         ! 111 - Fitted PT inspired thing
         hmod%ip2h = 5   ! 5 - One-loop SPT, dewiggled and damped
         hmod%ks = 0.03  ! One-halo damping wavenumber
         hmod%kd = 0.03  ! Two-halo damping wavenumber
         hmod%itrans = 0 ! For some reason two-halo term is negative at high z and high k
         hmod%ip1h = 0   ! No one-halo term
         IF (ihm == 111) THEN
            hmod%PT_A = 0.781
            hmod%PT_alpha = 0.414
            hmod%PT_beta = 0.479
         END IF
      ELSE IF (ihm == 110) THEN
         ! One-loop SPT
         hmod%ip2h = 4 ! 4 - One-loop SPT in two-halo term
      ELSE IF (ihm == 112) THEN  
         hmod%imf = 15 ! 15 - Tinker (2008) mass function with no z dependence
      ELSE IF (ihm == 113) THEN  
         hmod%imf = 16 ! 16 - Tinker (2010) mass function and peak-background split bias with no z dependence
      ELSE IF (ihm == 115) THEN
         hmod%imf = 17 ! 17 - Tinker (2010) mass function and bias
      ELSE IF (ihm == 119) THEN
         hmod%ip2h = 6 ! 6 - Non-linear power in two-halo term
      ELSE IF (ihm == 120) THEN
         hmod%ip2h = 7 ! 7 - Quasi-linear power in two-halo term
      ELSE IF (ihm == 121) THEN
         hmod%ip2h = 6 ! 6 - Non-linear power in two-halo term
         hmod%imf = 17 ! 17 - Tinker (2010) mass function and bias
      ELSE IF (ihm == 122) THEN
         hmod%ip2h = 7 ! 7 - Quasi-linear power in two-halo term
         hmod%imf = 17 ! 17 - Tinker (2010) mass function and bias
      ELSE
         STOP 'ASSIGN_HALOMOD: Error, ihm specified incorrectly'
      END IF
      hmod%name = names(ihm)

      IF (verbose) THEN
         WRITE (*, *) 'ASSIGN_HALOMOD: ', TRIM(names(ihm))
         WRITE (*, *) 'ASSIGN_HALOMOD: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE assign_halomod

   SUBROUTINE init_halomod(a, hmod, cosm, verbose)

      ! Halo-model initialisation routine
      ! Computes other tables necessary for the one-halo integral
      REAL, INTENT(IN) :: a
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i
      REAL :: z, Om_stars
      REAL, PARAMETER :: glim = g_integral_limit
      REAL, PARAMETER :: gblim = gb_integral_limit
      LOGICAL, PARAMETER :: check_mf = check_mass_function
      LOGICAL, PARAMETER :: calculate_stars = calculate_Omega_stars
      INTEGER, PARAMETER :: flag_sigmaV = flag_matter

      IF (verbose) WRITE (*, *) 'INIT_HALOMOD: Initialising calculation'

      ! Set the redshift (this routine needs to be called anew for each z)
      hmod%a = a
      z = redshift_a(a)
      hmod%z = z
      IF (verbose) THEN  
         WRITE (*, *) 'INIT_HALOMOD: Tables being filled at redshift:', REAL(z)
         WRITE (*, *) 'INIT_HALOMOD: Tables being filled at scale-factor:', REAL(a)
      END IF

      ! Reset flags to false
      hmod%has_galaxies = .FALSE.
      hmod%has_HI = .FALSE.
      hmod%has_mass_conversions = .FALSE.
      hmod%has_mass_function = .FALSE.
      hmod%has_bnl = .FALSE.

      ! Calculate sigma_v
      ! TODO: This is really a 'cosmology' thing
      ! TODO: This is only necessary for some halo models
      hmod%sigv = sigmaV(0., a, flag_sigmaV, cosm)    
      IF (verbose) WRITE (*, *) 'INIT_HALOMOD: sigma_V [Mpc/h]:', REAL(hmod%sigv)

      IF (cosm%img == img_fR) THEN
         hmod%mass_dependent_Dv = .TRUE.
      ELSE
         hmod%mass_dependent_Dv = .FALSE.
      END IF

      ! Calcuate delta_c and Delta_v
      ! Sometimes these calculations rely on the sigmas that should be calculated before this
      hmod%dc = delta_c(hmod, cosm)
      hmod%Dv = Delta_v(M0_Dv_default, hmod, cosm)
      IF (verbose) THEN         
         WRITE (*, *) 'INIT_HALOMOD: Delta_v:', REAL(hmod%Dv)
         WRITE (*, *) 'INIT_HALOMOD: delta_c:', REAL(hmod%dc)
      END IF

      ! Deallocate and reallocate arrays if necessary
      IF (ALLOCATED(hmod%log_m)) CALL deallocate_halomod(hmod)
      CALL allocate_halomod(hmod)

      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: Filling look-up tables'
         WRITE (*, *) 'INIT_HALOMOD: Number of entries:', hmod%n
      END IF

      ! Halo masses
      CALL fill_array_log(hmod%mmin, hmod%mmax, hmod%m, hmod%n)
      hmod%log_m = log(hmod%m)

      ! Lagrangian and virial radius
      hmod%rr = Lagrangian_radius(hmod%m, cosm)
      IF (hmod%DMONLY_neutrinos_affect_virial_radius) THEN
         hmod%rr = hmod%rr*(1.-cosm%f_nu)**(1./3.)
      END IF
      IF(hmod%mass_dependent_Dv) THEN
         DO i = 1, hmod%n
            hmod%rv(i) = virial_radius(hmod%m(i), hmod, cosm)
         END DO
      ELSE
         hmod%rv = hmod%rr/hmod%Dv**(1./3.)
      END IF

      ! Sigma and nu
      DO i = 1, hmod%n
         hmod%sig(i) = sigma(hmod%rr(i), a, hmod%flag_sigma, cosm)
      END DO
      hmod%nu = hmod%dc/hmod%sig

      IF (verbose) WRITE (*, *) 'INIT_HALOMOD: M, R, rv, sigma, nu tables filled'
      IF (hmod%nu(1) > 1.) WRITE(*, *) 'INIT_HALOMOD: WARNING: lowest nu value stored is greater than unity, z:', real(hmod%z)

      ! For spectra with finite variance we have cut-off masses etc.
      IF(cosm%warm) THEN
         hmod%saturation = .TRUE.
         hmod%nu_saturation = hmod%nu(1)
         IF(verbose) THEN
            WRITE(*, *) 'INIT_HALOMOD: Saturation nu:', hmod%nu_saturation
         END IF
      END IF

      ! Write some useful information to the screen
      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: Minimum nu:', REAL(hmod%nu(1))
         WRITE (*, *) 'INIT_HALOMOD: Maximum nu:', REAL(hmod%nu(hmod%n))
         WRITE (*, *) 'INIT_HALOMOD: Minimum sigma:', REAL(hmod%sig(hmod%n))
         WRITE (*, *) 'INIT_HALOMOD: Maximum sigma:', REAL(hmod%sig(1))
         WRITE (*, *) 'INIT_HALOMOD: Minimum R [Mpc/h]:', REAL(hmod%rr(1))
         WRITE (*, *) 'INIT_HALOMOD: Maximum R [Mpc/h]:', REAL(hmod%rr(hmod%n))
         WRITE (*, *) 'INIT_HALOMOD: Minimum R_v [Mpc/h]:', REAL(hmod%rv(1))
         WRITE (*, *) 'INIT_HALOMOD: Maximum R_v [Mpc/h]:', REAL(hmod%rv(hmod%n))
         WRITE (*, *) 'INIT_HALOMOD: Minimum log10(M) [Msun/h]:', REAL(log10(hmod%m(1)))
         WRITE (*, *) 'INIT_HALOMOD: Maximum log10(M) [Msun/h]:', REAL(log10(hmod%m(hmod%n)))
      END IF

      ! Calculate missing mass things if necessary
      IF ((hmod%imf /= 4) .AND. (hmod%i2hcor /= 0)) THEN

         IF (check_mf) hmod%gmin = 1.-integrate_g_nu(hmod%nu(1), hmod%large_nu, hmod)
         IF (check_mf) hmod%gmax = integrate_g_nu(hmod%nu(hmod%n), hmod%large_nu, hmod)
         hmod%gbmin = 1.-integrate_gb_nu(hmod%nu(1), hmod%large_nu, hmod)
         IF (check_mf) hmod%gbmax = integrate_gb_nu(hmod%nu(hmod%n), hmod%large_nu, hmod)
         IF (check_mf) hmod%gnorm = integrate_g_mu(hmod%small_nu, hmod%large_nu, hmod)

         IF (verbose) THEN
            IF (check_mf) WRITE (*, *) 'INIT_HALOMOD: Missing g(nu) at low end:', REAL(hmod%gmin)
            IF (check_mf) WRITE (*, *) 'INIT_HALOMOD: Missing g(nu) at high end:', REAL(hmod%gmax)
            WRITE (*, *) 'INIT_HALOMOD: Missing g(nu)b(nu) at low end:', REAL(hmod%gbmin)
            IF (check_mf) WRITE (*, *) 'INIT_HALOMOD: Missing g(nu)b(nu) at high end:', REAL(hmod%gbmax)
            IF (check_mf) WRITE (*, *) 'INIT_HALOMOD: Total g(nu) integration:', REAL(hmod%gnorm)
         END IF

         IF (check_mf) THEN
            IF (hmod%gmin < 0.) STOP 'INIT_HALOMOD: Error, missing g(nu) at low end is less than zero'
            IF (hmod%gmax < glim) STOP 'INIT_HALOMOD: Error, missing g(nu) at high end is less than zero'
            IF (hmod%gbmin < 0.) STOP 'INIT_HALOMOD: Error, missing g(nu)b(nu) at low end is less than zero'
            IF (hmod%gbmax < gblim) STOP 'INIT_HALOMOD: Error, missing g(nu)b(nu) at high end is less than zero'
         END IF

      END IF

      ! Find non-linear Lagrangian radius and associated scale
      ! This is defined as nu(M_star)=1 *not* sigma(M_star)=1, so depends on delta_c
      ! TODO: Is this necessar for all halo models?
      hmod%rnl = r_nl(hmod)
      hmod%mnl = Lagrangian_mass(hmod%rnl, cosm)
      hmod%knl = 1./hmod%rnl ! Note that there are no factors of 2pi here
      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: Non-linear mass [log10(M*) [Msun/h]]:', REAL(log10(hmod%mnl))
         WRITE (*, *) 'INIT_HALOMOD: Non-linear halo virial radius [Mpc/h]:', REAL(virial_radius(hmod%mnl, hmod, cosm))
         WRITE (*, *) 'INIT_HALOMOD: Non-linear Lagrangian radius [Mpc/h]:', REAL(hmod%rnl)
         WRITE (*, *) 'INIT_HALOMOD: Non-linear wavenumber [h/Mpc]:', REAL(hmod%knl)
      END IF

      ! Calculate the effective spectral index at the collapse scale
      ! TODO: Not necessary for some halo models; does this affect speed?
      hmod%neff = effective_index(hmod%flag_sigma, hmod, cosm)
      !hmod%ncur = effective_curvature(hmod%flag_sigma, hmod, cosm)
      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: Collapse n_eff:', REAL(hmod%neff)
         !WRITE (*, *) 'INIT_HALOMOD: Collapse curvature:', REAL(hmod%ncur)
      END IF

      ! Calculate the amplitude of the the one-halo term
      ! TODO: Not necessary for some halo models
      hmod%Rh = one_halo_amplitude(hmod, cosm)
      hmod%Mh = hmod%Rh*comoving_matter_density(cosm)
      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: One-halo amplitude [Mpc/h]^3:', hmod%Rh
         WRITE (*, *) 'INIT_HALOMOD: One-halo amplitude [log10(M) [Msun/h]]:', log10(hmod%Mh)
      END IF

      ! Applies if we enforce the amplitude of the one-halo term not to be affected by mass loss
      IF (hmod%DMONLY_baryon_recipe .AND. renormalise_onehalo) THEN
         hmod%Rhh = one_halo_modified_amplitude(hmod, cosm)
         IF (verbose) THEN
            WRITE (*, *) 'INIT_HALOMOD: One-halo modified amplitude [Mpc/h]^3:', hmod%Rhh
         END IF
      END IF

      ! Calculate the total stellar mass fraction
      ! TODO: Not necessary for some halo models
      IF (calculate_stars .AND. verbose) THEN
         Om_stars = Omega_stars(hmod, cosm)
         WRITE (*, *) 'INIT_HALOMOD: Omega_*:', Om_stars
         WRITE (*, *) 'INIT_HALOMOD: Omega_* / Omega_m:', Om_stars/cosm%Om_m
      END IF

      ! Calculate HMcode things
      ! TODO: Not necessary for some halo models
      hmod%HMcode_kstar = HMcode_kstar(hmod, cosm)
      hmod%HMcode_fdamp = HMcode_fdamp(hmod, cosm)
      hmod%HMcode_kdamp = HMcode_kdamp(hmod, cosm)
      hmod%HMcode_alpha = HMcode_alpha(hmod, cosm)
      hmod%HMcode_eta = HMcode_eta(hmod, cosm)
      hmod%HMcode_A = HMcode_A(hmod, cosm)
      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: HMcode: k*:', hmod%HMcode_kstar
         WRITE (*, *) 'INIT_HALOMOD: HMcode: fdamp:', hmod%HMcode_fdamp
         WRITE (*, *) 'INIT_HALOMOD: HMcode: kdamp [h/Mpc]:', hmod%HMcode_kdamp
         WRITE (*, *) 'INIT_HALOMOD: HMcode: alpha:', hmod%HMcode_alpha
         WRITE (*, *) 'INIT_HALOMOD: HMcode: eta:', hmod%HMcode_eta
         WRITE (*, *) 'INIT_HALOMOD: HMcode: A:', hmod%HMcode_A
      END IF

      ! Calculate the concentration-mass relation
      ! Note that this requires mnl etc. for some c(M) relations. 
      ! This should be the final thing in init_halomod. Do not move!
      CALL fill_halo_concentration(hmod, cosm)
      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: Halo concentration tables filled'
         WRITE (*, *) 'INIT_HALOMOD: Minimum concentration:', REAL(hmod%c(hmod%n))
         WRITE (*, *) 'INIT_HALOMOD: Maximum concentration:', REAL(hmod%c(1))
      END IF

      ! Initialise HALOFIT if required
      IF (hmod%ip2h == 6 .OR. hmod%ip2h == 7) THEN
         CALL HALOFIT_init(hmod%HALOFIT_knl, hmod%HALOFIT_neff, hmod%HALOFIT_ncur, hmod%a, cosm, verbose)
      END IF

      ! Check for second-order halo bias
      IF (hmod%ibias == 2) STOP 'INIT_HALOMOD: Error, check the second-order halo bias calculation very carefully'

      ! Check that combination of two-halo term and halo bias is valid
      IF ((hmod%ip2h .NE. 2) .AND. is_in_array(hmod%ibias, [3, 4, 5])) THEN
         STOP 'INIT_HALOMOD: Error, your combination of two-halo term and halo bias is not supported'
      END IF

      ! Finish
      IF (verbose) THEN
         WRITE (*, *) 'INIT_HALOMOD: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_halomod

   SUBROUTINE print_halomod(hmod, cosm, verbose)

      ! This subroutine writes out the physical halo-model parameters at some redshift
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      CHARACTER(len=12), PARAMETER :: fmt = '(A30,F10.5)'
      CHARACTER(len=43), PARAMETER :: dashes = '==========================================='

      IF (verbose) THEN

         WRITE (*, *) dashes
         WRITE (*, *) 'HALOMODEL: ', TRIM(hmod%name)
         WRITE (*, *) dashes
         WRITE (*, *) 'HALOMODEL: Accuracy parameters'
         WRITE (*, *) dashes
         WRITE (*, *) 'HALOMODEL: Lower mass for integration [log10(Msun/h)]:', log10(hmod%mmin)
         WRITE (*, *) 'HALOMODEL: Upper mass for integration [log10(Msun/h)]:', log10(hmod%mmax)
         WRITE (*, *) 'HALOMODEL: Number of points in look-up tables:', hmod%n
         WRITE (*, *) 'HALOMODEL: Halo-model accuracy parameter:', hmod%acc
         WRITE (*, *) 'HALOMODEL: Large value of nu:', hmod%large_nu
         WRITE (*, *) dashes

         ! delta_c
         IF (hmod%idc == 1) WRITE (*, *) 'HALOMODEL: Fixed delta_c = 1.686'
         IF (hmod%idc == 2) WRITE (*, *) 'HALOMODEL: delta_c from Nakamura & Suto (1997) fitting function'
         IF (hmod%idc == 3) WRITE (*, *) 'HALOMODEL: delta_c from HMcode (2016) power spectrum fit'
         IF (hmod%idc == 4) WRITE (*, *) 'HALOMODEL: delta_c from Mead (2017) fitting function'
         IF (hmod%idc == 5) WRITE (*, *) 'HALOMODEL: delta_c from spherical-collapse calculation'
         IF (hmod%idc == 6) WRITE (*, *) 'HALOMODEL: delta_c from HMcode (2015) power spectrum fit'

         ! Delta_v
         IF (hmod%iDv == 1)  WRITE (*, *) 'HALOMODEL: Delta_v = 200 fixed'
         IF (hmod%iDv == 2)  WRITE (*, *) 'HALOMODEL: Delta_v from Bryan & Norman (1998) fitting function'
         IF (hmod%iDv == 3)  WRITE (*, *) 'HALOMODEL: Delta_v from HMcode (2016) power spectrum fit'
         IF (hmod%iDv == 4)  WRITE (*, *) 'HALOMODEL: Delta_v from Mead (2017) fitting function'
         IF (hmod%iDv == 5)  WRITE (*, *) 'HALOMODEL: Delta_v from spherical-collapse calculation'
         IF (hmod%iDv == 6)  WRITE (*, *) 'HALOMODEL: Delta_v to give haloes Lagrangian radius'
         IF (hmod%iDv == 7)  WRITE (*, *) 'HALOMODEL: Delta_v for M200c'
         IF (hmod%iDv == 8)  WRITE (*, *) 'HALOMODEL: Delta_v from HMcode (2015) power spectrum fit'
         IF (hmod%iDv == 9)  WRITE (*, *) 'HALOMODEL: Delta_v = 178 fixed'
         IF (hmod%iDv == 10) WRITE (*, *) 'HALOMODEL: Delta_v user defined'

         ! Form of the two-halo term
         IF (hmod%ip2h == 0) WRITE (*, *) 'HALOMODEL: No two-halo term'
         IF (hmod%ip2h == 1) WRITE (*, *) 'HALOMODEL: Linear two-halo term'
         IF (hmod%ip2h == 2) WRITE (*, *) 'HALOMODEL: Standard two-halo term (Seljak 2000)'
         IF (hmod%ip2h == 3) WRITE (*, *) 'HALOMODEL: Linear two-halo term with damped wiggles' 
         IF (hmod%ip2h == 4) WRITE (*, *) 'HALOMODEL: One-loop SPT'  
         IF (hmod%ip2h == 5) WRITE (*, *) 'HALOMODEL: One-loop SPT, dewiggle and damp'
         IF (hmod%ip2h == 6) WRITE (*, *) 'HALOMODEL: Non-linear two-halo term from HALOFIT'
         IF (hmod%ip2h == 7) WRITE (*, *) 'HALOMODEL: Quasi-linear two-halo term from HALOFIT'

         ! Form of the one-halo term
         IF (hmod%ip1h == 0) WRITE (*, *) 'HALOMODEL: No one-halo term'
         IF (hmod%ip1h == 1) WRITE (*, *) 'HALOMODEL: Standard one-halo term'
         
         ! Order to go to in halo bias
         IF (hmod%ibias == 1) WRITE (*, *) 'HALOMODEL: Linear halo bias'
         IF (hmod%ibias == 2) WRITE (*, *) 'HALOMODEL: Second-order halo bias'
         IF (hmod%ibias == 3) WRITE (*, *) 'HALOMODEL: Full non-linear halo bias'
         IF (hmod%ibias == 4) WRITE (*, *) 'HALOMODEL: Scale-dependent halo bias fudge from Fedeli (2014b)'
         IF (hmod%ibias == 5) WRITE (*, *) 'HALOMODEL: Scale-dependent halo bias fudge from me'

         ! Correction for missing low-mass haloes
         IF (hmod%i2hcor == 0) WRITE (*, *) 'HALOMODEL: No two-halo correction applied for missing low-mass haloes'
         IF (hmod%i2hcor == 1) WRITE (*, *) 'HALOMODEL: Two-halo term corrected by adding missing g(nu)b(nu)'
         IF (hmod%i2hcor == 2) WRITE (*, *) 'HALOMODEL: Two-halo term corrected via delta function at low mass end'

         ! Halo mass function
         IF (hmod%imf == 1)  WRITE (*, *) 'HALOMODEL: Press & Schecter (1974) mass function'
         IF (hmod%imf == 2)  WRITE (*, *) 'HALOMODEL: Sheth & Tormen (1999) mass function'
         IF (hmod%imf == 3)  WRITE (*, *) 'HALOMODEL: Tinker et al. (2010) mass function with peak-background split bias'
         IF (hmod%imf == 4)  WRITE (*, *) 'HALOMODEL: Delta function mass function'
         IF (hmod%imf == 5)  WRITE (*, *) 'HALOMODEL: Jenkins et al. (2001) mass function'
         IF (hmod%imf == 6)  WRITE (*, *) 'HALOMODEL: Despali et al. (2016) mass function'
         IF (hmod%imf == 7)  WRITE (*, *) 'HALOMODEL: Tinker et al. (2008) mass function'
         IF (hmod%imf == 8)  WRITE (*, *) 'HALOMODEL: Warren et al. (2006) mass function'
         IF (hmod%imf == 9)  WRITE (*, *) 'HALOMODEL: Reed et al. (2007) mass function'
         IF (hmod%imf == 10) WRITE (*, *) 'HALOMODEL: Bhattacharya et al. (2011) mass function'
         IF (hmod%imf == 11) WRITE (*, *) 'HALOMODEL: Tinker et al. (2010) mass function with no halo bias and z dependence'
         IF (hmod%imf == 12) WRITE (*, *) 'HALOMODEL: Sheth & Tormen (2002) mass function'
         IF (hmod%imf == 13) WRITE (*, *) 'HALOMODEL: Sheth & Tormen form of mass function'
         IF (hmod%imf == 14) WRITE (*, *) 'HALOMODEL: Philcox et al. (2020) mass function'
         IF (hmod%imf == 15) WRITE (*, *) 'HALOMODEL: Tinker et al. (2008) mass function with no z dependence'
         IF (hmod%imf == 16) WRITE (*, *) 'HALOMODEL: Tinker et al. (2010) mass function with no z dependence'
         IF (hmod%imf == 17) WRITE (*, *) 'HALOMODEL: Tinker et al. (2010) mass function with calibrated halo bias'
         IF (hmod%imf == 18) WRITE (*, *) 'HALOMODEL: Tinker et al. (2010) mass function with calibrated halo bias and no z dependence'
         IF (hmod%imf == 19) WRITE (*, *) 'HALOMODEL: Sheth & Tormen (1999) mass function with calibrated halo bias'

         ! Concentration-mass relation
         IF (hmod%iconc == 1)  WRITE (*, *) 'HALOMODEL: Full Bullock et al. (2001) concentration-mass relation'
         IF (hmod%iconc == 2)  WRITE (*, *) 'HALOMODEL: Simple Bullock et al. (2001) concentration-mass relation'
         IF (hmod%iconc == 3)  WRITE (*, *) 'HALOMODEL: Full sample for M200 Duffy et al. (2008) concentration-mass relation'
         IF (hmod%iconc == 4)  WRITE (*, *) 'HALOMODEL: Full sample for Mv Duffy et al. (2008) concentration-mass relation'
         IF (hmod%iconc == 5)  WRITE (*, *) 'HALOMODEL: Full sample for M200c Duffy et al. (2008) concentration-mass relation'
         IF (hmod%iconc == 6)  WRITE (*, *) 'HALOMODEL: Relaxed sample for M200 Duffy et al. (2008) concentration-mass relation'
         IF (hmod%iconc == 7)  WRITE (*, *) 'HALOMODEL: Relaxed sample for Mv Duffy et al. (2008) concentration-mass relation'
         IF (hmod%iconc == 8)  WRITE (*, *) 'HALOMODEL: Relaxed sample for M200c Duffy et al. (2008) concentration-mass relation'
         IF (hmod%iconc == 9)  WRITE (*, *) 'HALOMODEL: Child et al. (2018) M200c concentration-mass relation'
         IF (hmod%iconc == 10) WRITE (*, *) 'HALOMODEL: Diemer & Kravstov (2019) M200c concentration-mass relation'

         ! Concentration-mass relation correction
         IF (hmod%iDolag == 0) WRITE (*, *) 'HALOMODEL: No concentration-mass correction for dark energy'
         IF (hmod%iDolag == 1) WRITE (*, *) 'HALOMODEL: Dolag (2004) concentration correction'
         IF (hmod%iDolag == 2) WRITE (*, *) 'HALOMODEL: Dolag (2004) concentration correction with 1.5 exponent'
         IF (hmod%iDolag == 3) WRITE (*, *) 'HALOMODEL: Dolag (2004) concentration correction with z dependence'
         IF (hmod%iDolag == 4) WRITE (*, *) 'HALOMODEL: Dolag (2004) concentration correction with neutrinos removed and z dependence'
         IF (hmod%iDolag == 5) WRITE (*, *) 'HALOMODEL: Dolag (2004) concentration correction with collapse redshift dependence'

         ! Bound gas fraction
         IF (hmod%frac_bound_gas == 1) WRITE (*, *) 'HALOMODEL: Halo bound gas fraction: Fedeli (2014)'
         IF (hmod%frac_bound_gas == 2) WRITE (*, *) 'HALOMODEL: Halo bound gas fraction: Schneider & Teyssier (2015)'
         IF (hmod%frac_bound_gas == 3) WRITE (*, *) 'HALOMODEL: Halo bound gas fraction: Universal baryon fraction'

         ! Cold gas fraction
         IF (hmod%frac_cold_bound_gas == 1) WRITE (*, *) 'HALOMODEL: Cold gas fraction: Constant fraction of bound halo gas'

         ! Hot gas fraction
         IF (hmod%frac_hot_bound_gas == 1) WRITE (*, *) 'HALOMODEL: Hot gas fraction: Constant fraction of bound halo gas'

         ! Central star fraction
         IF (hmod%frac_stars == 1) WRITE (*, *) 'HALOMODEL: Halo star fraction: Fedeli (2014)'
         IF (hmod%frac_stars == 2) WRITE (*, *) 'HALOMODEL: Halo star fraction: Constant'
         IF (hmod%frac_stars == 3) WRITE (*, *) 'HALOMODEL: Halo star fraction: Fedeli (2014) but saturated at high halo mass'

         ! Satellite star fraction
         IF (hmod%frac_central_stars == 1) WRITE (*, *) 'HALOMODEL: Central star fraction: All stars are central stars'
         IF (hmod%frac_central_stars == 2) WRITE (*, *) 'HALOMODEL: Central star fraction: Schneider et al. (2018)'

         ! HI fraction
         IF (hmod%frac_HI == 1) WRITE (*, *) 'HALOMODEL: Halo HI fraction: Hard truncation at high and low masses'
         IF (hmod%frac_HI == 2) WRITE (*, *) 'HALOMODEL: Halo HI fraction: Villaescusa-Navarro et al. (2018; 1804.09180)'

         ! DMONLY halo model
         IF (hmod%halo_DMONLY == 1) WRITE (*, *) 'HALOMODEL: DMONLY halo profile: NFW'
         IF (hmod%halo_DMONLY == 2) WRITE (*, *) 'HALOMODEL: DMONLY halo profile: Non-analytical NFW (for testing W(k) functions)'
         IF (hmod%halo_DMONLY == 3) WRITE (*, *) 'HALOMODEL: DMONLY halo profile: Tophat'
         IF (hmod%halo_DMONLY == 4) WRITE (*, *) 'HALOMODEL: DMONLY halo profile: Delta function'
         IF (hmod%halo_DMONLY == 5) WRITE (*, *) 'HALOMODEL: DMONLY halo profile: Cored NFW'
         IF (hmod%halo_DMONLY == 6) WRITE (*, *) 'HALOMODEL: DMONLY halo profile: Isothermal'
         IF (hmod%halo_DMONLY == 7) WRITE (*, *) 'HALOMODEL: DMONLY halo profile: Shell'

         ! CDM halo profile
         IF (hmod%halo_CDM == 1) WRITE (*, *) 'HALOMODEL: CDM halo profile: NFW'

         ! Bound gas halo profile
         IF (hmod%halo_normal_bound_gas == 1) WRITE (*, *) 'HALOMODEL: Static gas profile: Simplified Komatsu & Seljak (2001)'
         IF (hmod%halo_normal_bound_gas == 2) WRITE (*, *) 'HALOMODEL: Static gas profile: Isothermal beta profile'
         IF (hmod%halo_normal_bound_gas == 3) WRITE (*, *) 'HALOMODEL: Static gas profile: Full Komatsu & Seljak (2001)'
         IF (hmod%halo_normal_bound_gas == 4) WRITE (*, *) 'HALOMODEL: Static gas profile: NFW'

         ! Cold gas profile
         IF (hmod%halo_cold_bound_gas == 1) WRITE (*, *) 'HALOMODEL: Cold gas profile: Delta function'

         ! Hot gas profile
         IF (hmod%halo_hot_bound_gas == 1) WRITE (*, *) 'HALOMODEL: Hot gas profile: Isothermal'

         ! Free gas halo profile
         IF (hmod%halo_ejected_gas == 1) WRITE (*, *) 'HALOMODEL: Free gas profile: Isothermal model (out to 2rv)'
         IF (hmod%halo_ejected_gas == 2) WRITE (*, *) 'HALOMODEL: Free gas profile: Ejected gas model (Schneider & Teyssier 2015)'
         IF (hmod%halo_ejected_gas == 3) WRITE (*, *) 'HALOMODEL: Free gas profile: Isothermal shell that connects at rv'
         IF (hmod%halo_ejected_gas == 4) WRITE (*, *) 'HALOMODEL: Free gas profile: Komatsu-Seljak continuation'
         IF (hmod%halo_ejected_gas == 5) WRITE (*, *) 'HALOMODEL: Free gas profile: Power-law continuation'
         IF (hmod%halo_ejected_gas == 6) WRITE (*, *) 'HALOMODEL: Free gas profile: Cubic profile'
         IF (hmod%halo_ejected_gas == 7) WRITE (*, *) 'HALOMODEL: Free gas profile: Smoothly distributed'
         IF (hmod%halo_ejected_gas == 8) WRITE (*, *) 'HALOMODEL: Free gas profile: Delta function'

         ! Central star halo profile
         IF (hmod%halo_central_stars == 1) WRITE (*, *) 'HALOMODEL: Central star profile: Fedeli (2014)'
         IF (hmod%halo_central_stars == 2) WRITE (*, *) 'HALOMODEL: Central star profile: Schneider & Teyssier (2015)'
         IF (hmod%halo_central_stars == 3) WRITE (*, *) 'HALOMODEL: Central star profile: Delta function'
         IF (hmod%halo_central_stars == 4) WRITE (*, *) 'HALOMODEL: Central star profile: Delta at low mass and NFW at high mass'

         ! Satellite star halo profile
         IF (hmod%halo_satellite_stars == 1) WRITE (*, *) 'HALOMODEL: Satellite star profile: NFW'
         IF (hmod%halo_satellite_stars == 2) WRITE (*, *) 'HALOMODEL: Satellite star profile: Schneider & Teyssier (2015)'

         ! Neurtino halo profile
         IF (hmod%halo_neutrino == 1) WRITE (*, *) 'HALOMODEL: Neutrino profile: Smoothly distributed'
         IF (hmod%halo_neutrino == 2) WRITE (*, *) 'HALOMODEL: Neutrino profile: NFW'

         ! HI halo profile
         IF (hmod%halo_HI == 1) WRITE (*, *) 'HALOMODEL: HI profile: NFW'
         IF (hmod%halo_HI == 2) WRITE (*, *) 'HALOMODEL: HI profile: Delta function'
         IF (hmod%halo_HI == 3) WRITE (*, *) 'HALOMODEL: HI profile: Polynomial with hole (Villaescusa-Navarro et al. 2018)'
         IF (hmod%halo_HI == 4) WRITE (*, *) 'HALOMODEL: HI profile: Modified NFW with hole (Villaescusa-Navarro et al. 2018)'
         IF (hmod%halo_HI == 5) WRITE (*, *) 'HALOMODEL: HI profile: Modified NFW (Padmanabhan & Refregier 2017)'

         ! Satellite galaxy halo profile
         IF (hmod%halo_satellites == 1) WRITE (*, *) 'HALOMODEL: Satellite galaxy profile: NFW'
         IF (hmod%halo_satellites == 2) WRITE (*, *) 'HALOMODEL: Satellite galaxy profile: Isothermal'

         ! Electron pressure profile
         IF (hmod%electron_pressure == 1) WRITE (*, *) 'HALOMODEL: Electron pressure: Using UPP'
         IF (hmod%electron_pressure == 2) WRITE (*, *) 'HALOMODEL: Electron pressure: Using gas profiles'

         ! Are voids being done?
         IF (hmod%add_voids) THEN

            WRITE (*, *) 'HALOMODEL: Considering voids'

            ! Void model
            IF (hmod%halo_void == 1) WRITE (*, *) 'HALOMODEL: Void model: Top-hat'

            ! Compensated void model
            IF (hmod%halo_compensated_void == 1) WRITE (*, *) 'HALOMODEL: Compensated void model: Top-hat'

         END IF

         ! sigma(R) type
         IF (hmod%flag_sigma == flag_cold) THEN
            WRITE (*, *) 'HALOMODEL: Sigma being calculated using cold matter'
         ELSE IF (hmod%flag_sigma == flag_ucold) THEN
            WRITE (*, *) 'HALOMODEL: Sigma being calculated using un-normalised cold matter'
         ELSE IF (hmod%flag_sigma == flag_matter) THEN
            WRITE (*, *) 'HALOMODEL: Sigma being calculated using total matter'
         ELSE
            STOP 'HALOMODEL: Error, flag for cold/matter specified incorreclty'
         END IF

         ! eta for halo window function
         IF (hmod%ieta == 0) WRITE (*, *) 'HALOMODEL: eta = 0 fixed'
         IF (hmod%ieta == 1) WRITE (*, *) 'HALOMODEL: eta from HMcode (2015, 2016) power spectrum fit'
         IF (hmod%ieta == 2) WRITE (*, *) 'HALOMODEL: eta from HMcode (2015, 2016) but with cold dependence'
         IF (hmod%ieta == 3) WRITE (*, *) 'HALOMODEL: eta from HMcode (2020)'

         ! Small-scale two-halo term damping coefficient
         IF (hmod%i2hdamp == 0) WRITE (*, *) 'HALOMODEL: No two-halo term damping at small scales'
         IF (hmod%i2hdamp == 1) WRITE (*, *) 'HALOMODEL: Two-halo term damping from HMcode (2015)'
         IF (hmod%i2hdamp == 2) WRITE (*, *) 'HALOMODEL: Two-halo term damping from HMcode (2016)'
         IF (hmod%i2hdamp == 3) WRITE (*, *) 'HALOMODEL: Two-halo term damping from HMcode (2020)'

         ! Large-scale one-halo term damping function
         IF (hmod%i1hdamp == 0) WRITE (*, *) 'HALOMODEL: No damping in one-halo term at large scales'
         IF (hmod%i1hdamp == 1) WRITE (*, *) 'HALOMODEL: One-halo large-scale damping via an exponential'
         IF (hmod%i1hdamp == 2) WRITE (*, *) 'HALOMODEL: One-halo large-scale P(k) ~ (k/k*)^4 with simple k*'
         IF (hmod%i1hdamp == 3) WRITE (*, *) 'HALOMODEL: One-halo large-scale P(k) ~ (k/k*)^4 with parameterised k*'

         ! Concentration-mass scaling
         IF (hmod%iAs == 0) WRITE (*, *) 'HALOMODEL: No rescaling of concentration-mass relation'
         IF (hmod%iAs == 1) WRITE (*, *) 'HALOMODEL: Concentration-mass rescaled mass independently (HMcode 2015, 2016)'
         IF (hmod%iAs == 2) WRITE (*, *) 'HALOMODEL: Concentration-mass rescaled as function of sigma (HMcode 2020)'

         ! Two- to one-halo transition region
         IF (hmod%itrans == 0) WRITE (*, *) 'HALOMODEL: Standard sum of two- and one-halo terms'
         IF (hmod%itrans == 1) WRITE (*, *) 'HALOMODEL: Smoothed transition using alpha'
         IF (hmod%itrans == 2) WRITE (*, *) 'HALOMODEL: Experimental smoothed transition'
         IF (hmod%itrans == 3) WRITE (*, *) 'HALOMODEL: Experimental smoothed transition for HMx (2.5 exponent)'
         IF (hmod%itrans == 4) WRITE (*, *) 'HALOMODEL: Tanh transition with k_nl'
         IF (hmod%itrans == 5) WRITE (*, *) 'HALOMODEL: Experimental smoothed transition 2'

         ! Response
         IF (hmod%response_baseline .NE. 0) THEN
            WRITE (*, *) 'HALOMODEL: Halo model is response with baseline:', hmod%response_baseline
            IF (hmod%response_denominator .NE. 0) THEN
               WRITE (*, *) 'HALOMODEL: Response denominator:', hmod%response_denominator
            ELSE
               WRITE (*, *) 'HALOMODEL: Response denominator is DMONLY model'
            END IF
            IF (hmod%response_matter_only) THEN
               WRITE (*, *) 'HALOMODEL: Response applies to matter spectra'
            ELSE
               WRITE (*, *) 'HALOMODEL: Response applies to all spectra'
            END IF
         END IF
         WRITE (*, *) dashes

         ! Numerical parameters
         WRITE (*, *) 'HALOMODEL: Standard parameters'
         WRITE (*, *) dashes
         WRITE (*, fmt=fmt) 'redshift:', hmod%z
         WRITE (*, fmt=fmt) 'scale factor:', hmod%a
         WRITE (*, *) dashes

         ! Perturbation theory
         WRITE (*, *) 'HALOMODEL: Perturbation theory'
         WRITE (*, *) dashes
         WRITE (*, fmt=fmt) 'A:', hmod%PT_A
         WRITE (*, fmt=fmt) 'alpha:', hmod%PT_alpha
         WRITE (*, fmt=fmt) 'beta:', hmod%PT_beta
         WRITE (*, *) dashes

         ! Halo-mass function
         WRITE (*, *) 'HALOMODEL: Mass function parameters'
         WRITE (*, *) dashes
         WRITE (*, fmt=fmt) 'Delta_v:', hmod%Dv
         WRITE (*, fmt=fmt) 'delta_c:', hmod%dc
         WRITE (*, fmt=fmt) 'Amplitude:', hmod%Amf
         WRITE (*, fmt=fmt) 'Amplitude z:', hmod%Amfz
         IF (is_in_array(hmod%imf, [2, 6, 12, 13, 19])) THEN
            WRITE (*, fmt=fmt) 'Sheth & Tormen p:', hmod%ST_p
            WRITE (*, fmt=fmt) 'Sheth & Tormen q:', hmod%ST_q
            WRITE (*, fmt=fmt) 'Sheth & Tormen A:', hmod%ST_A
         ELSE IF (is_in_array(hmod%imf, [3, 11, 16, 17, 18])) THEN
            ! Tinker (2010) mass function
            WRITE (*, fmt=fmt) 'Tinker alpha:', hmod%Tinker_alpha
            WRITE (*, fmt=fmt) 'Tinker beta:', hmod%Tinker_beta
            WRITE (*, fmt=fmt) 'Tinker gamma:', hmod%Tinker_gamma
            WRITE (*, fmt=fmt) 'Tinker eta:', hmod%Tinker_eta
            WRITE (*, fmt=fmt) 'Tinker phi:', hmod%Tinker_phi
         ELSE IF (hmod%imf == 4) THEN
            WRITE (*, *) 'Halo mass [log10(Msun/h)]:', log10(hmod%hmass)
         ELSE IF (hmod%imf == 7 .OR. hmod%imf == 15) THEN
            ! Tinker (2008)
            WRITE (*, fmt=fmt) 'Tinker A:', hmod%Tinker_bigA
            WRITE (*, fmt=fmt) 'Tinker a:', hmod%Tinker_a
            WRITE (*, fmt=fmt) 'Tinker b:', hmod%Tinker_b
            WRITE (*, fmt=fmt) 'Tinker c:', hmod%Tinker_c
         END IF
         IF (hmod%imf == 17 .OR. hmod%imf == 18) THEN
            ! Tinker (2010) calibrated halo bias
            WRITE (*, fmt=fmt) 'Tinker A:', hmod%Tinker_bigA
            WRITE (*, fmt=fmt) 'Tinker B:', hmod%Tinker_bigB
            WRITE (*, fmt=fmt) 'Tinker C:', hmod%Tinker_bigC
            WRITE (*, fmt=fmt) 'Tinker a:', hmod%Tinker_a
            WRITE (*, fmt=fmt) 'Tinker b:', hmod%Tinker_b
            WRITE (*, fmt=fmt) 'Tinker c:', hmod%Tinker_c
         END IF
         WRITE (*, *) dashes

         ! Halo parameters
         WRITE (*, *) 'HALOMODEL: Halo parameters'
         WRITE (*, *) dashes       
         WRITE (*, fmt=fmt) 'Dolag zc:', hmod%zD
         WRITE (*, *) dashes

         ! HMcode parameters
         WRITE (*, *) 'HALOMODEL: HMcode parameters'
         WRITE (*, *) dashes
         WRITE (*, fmt=fmt) 'Dv0:', hmod%Dv0
         WRITE (*, fmt=fmt) 'Dv1:', hmod%Dv1
         WRITE (*, fmt=fmt) 'dc0:', hmod%dc0
         WRITE (*, fmt=fmt) 'dc1:', hmod%dc1
         WRITE (*, fmt=fmt) 'eta0:', hmod%eta0
         WRITE (*, fmt=fmt) 'eta1:', hmod%eta1
         WRITE (*, fmt=fmt) 'f0:', hmod%f0
         WRITE (*, fmt=fmt) 'f1:', hmod%f1
         WRITE (*, fmt=fmt) 'ks:', hmod%ks
         WRITE (*, fmt=fmt) 'kp:', hmod%kp
         WRITE (*, fmt=fmt) 'kd:', hmod%kd
         WRITE (*, fmt=fmt) 'kdp:', hmod%kdp
         WRITE (*, fmt=fmt) 'nd:', hmod%nd
         WRITE (*, fmt=fmt) 'As:', hmod%As
         WRITE (*, fmt=fmt) 'Ap:', hmod%Ap
         WRITE (*, fmt=fmt) 'Ac:', hmod%Ac
         WRITE (*, fmt=fmt) 'Az:', hmod%Az
         WRITE (*, fmt=fmt) 'alpha0:', hmod%alp0
         WRITE (*, fmt=fmt) 'alpha1:', hmod%alp1
         WRITE (*, fmt=fmt) 'alpha2:', hmod%alp2
         WRITE (*, fmt=fmt) 'fM:', hmod%fM
         WRITE (*, fmt=fmt) 'fD:', hmod%fD
         WRITE (*, fmt=fmt) 'Dvnu:', hmod%Dvnu
         WRITE (*, fmt=fmt) 'dcnu:', hmod%dcnu
         IF (hmod%DMONLY_baryon_recipe) THEN
            IF (hmod%mbar == 0.) THEN
               WRITE (*, fmt=fmt) 'M_bar [Msun/h]:', hmod%mbar
            ELSE
               WRITE (*, fmt=fmt) 'log10(M_bar) [Msun/h]:', log10(hmod%mbar)
            END IF
            WRITE (*, fmt=fmt) 'm_bar z:', hmod%mbarz
            WRITE (*, fmt=fmt) 'n_bar:', hmod%nbar
            WRITE (*, fmt=fmt) 's_bar:', hmod%sbar
            WRITE (*, fmt=fmt) 's_bar z:', hmod%sbarz
         END IF
         WRITE (*, *) dashes

         ! HMcode variables
         WRITE (*, *) 'HALOMODEL: HMcode variables'
         WRITE (*, *) dashes
         WRITE (*, fmt=fmt) 'k*:', hmod%HMcode_kstar
         WRITE (*, fmt=fmt) 'fdamp:', hmod%HMcode_fdamp
         WRITE (*, fmt=fmt) 'kdamp:', hmod%HMcode_kdamp
         WRITE (*, fmt=fmt) 'alpha:', hmod%HMcode_alpha
         WRITE (*, fmt=fmt) 'eta:', hmod%HMcode_eta
         WRITE (*, fmt=fmt) 'A:', hmod%HMcode_A
         WRITE (*, *) dashes

         ! HMx
         WRITE (*, *) 'HALOMODEL: HMx parameters'
         WRITE (*, *) dashes
         IF (hmod%HMx_mode == 4) STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
         IF (hmod%HMx_mode == 1) THEN
            WRITE (*, fmt=fmt) 'eps:', HMx_eps(hmod, cosm)
            WRITE (*, fmt=fmt) 'log10(M0) [Msun/h]:', log10(HMx_M0(hmod, cosm))
            WRITE (*, fmt=fmt) 'beta:', HMx_beta(hmod%Mh, hmod, cosm)
            WRITE (*, fmt=fmt) 'A*:', HMx_Astar(hmod, cosm)
            WRITE (*, *) dashes
            WRITE (*, fmt=fmt) 'eps:', hmod%eps  
            IF (hmod%M0 == 0.) THEN 
               WRITE (*, fmt=fmt) 'M0 [Msun/h]:', hmod%M0
            ELSE   
               WRITE (*, fmt=fmt) 'log10(M0) [Msun/h]:', log10(hmod%M0)
            END IF
            WRITE (*, fmt=fmt) 'M0 z index:', hmod%M0z         
            WRITE (*, fmt=fmt) 'beta:', hmod%beta         
            WRITE (*, fmt=fmt) 'A*:', hmod%Astar
            WRITE (*, fmt=fmt) 'A* z index:', hmod%Astarz       
         ELSE IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
            WRITE (*, fmt=fmt) 'alpha:', hmod%alpha
            WRITE (*, fmt=fmt) 'alpha mass index:', hmod%alphap
            WRITE (*, fmt=fmt) 'alpha z index:', hmod%alphaz
            WRITE (*, fmt=fmt) 'beta:', hmod%beta
            WRITE (*, fmt=fmt) 'beta mass index:', hmod%betap
            WRITE (*, fmt=fmt) 'beta z index:', hmod%betaz
            WRITE (*, fmt=fmt) 'epsilon:', hmod%eps
            WRITE (*, fmt=fmt) 'epsilon z index:', hmod%epsz
            WRITE (*, fmt=fmt) 'epsilon2:', hmod%eps2
            WRITE (*, fmt=fmt) 'epsilon2 z index:', hmod%eps2z
            WRITE (*, fmt=fmt) 'Gamma:', hmod%Gamma
            WRITE (*, fmt=fmt) 'Gamma mass index:', hmod%Gammap
            WRITE (*, fmt=fmt) 'Gamma z index:', hmod%Gammaz
            IF (hmod%M0 == 0.) THEN
               WRITE (*, fmt=fmt) 'M0 [Msun/h]:', hmod%M0
            ELSE
               WRITE (*, fmt=fmt) 'log10(M0) [Msun/h]:', log10(hmod%M0)
            END IF
            WRITE (*, fmt=fmt) 'M0 z index:', hmod%M0z
            WRITE (*, fmt=fmt) 'A*:', hmod%Astar
            WRITE (*, fmt=fmt) 'A* z index:', hmod%Astarz
            IF (hmod%Twhim == 0.) THEN
               WRITE (*, fmt=fmt) 'T_WHIM [K]:', hmod%Twhim
            ELSE
               WRITE (*, fmt=fmt) 'log10(T_WHIM) [K]:', log10(hmod%Twhim)
            END IF
            WRITE (*, fmt=fmt) 'T_WHIM z index:', hmod%Twhimz
            WRITE (*, fmt=fmt) 'c*:', hmod%cstar
            WRITE (*, fmt=fmt) 'c* mass index:', hmod%cstarp
            WRITE (*, fmt=fmt) 'c* z index:', hmod%cstarz
            WRITE (*, fmt=fmt) 'sigma*:', hmod%sstar
            WRITE (*, fmt=fmt) 'log10(M*) [Msun/h]:', log10(hmod%Mstar)
            WRITE (*, fmt=fmt) 'M* z index:', hmod%Mstarz
            WRITE (*, fmt=fmt) 'f_cold:', hmod%fcold
            WRITE (*, fmt=fmt) 'f_hot:', hmod%fhot
            WRITE (*, fmt=fmt) 'eta:', hmod%eta
            WRITE (*, fmt=fmt) 'eta z index:', hmod%etaz      
            WRITE (*, fmt=fmt) 'iso beta:', hmod%ibeta
            WRITE (*, fmt=fmt) 'iso beta mass index:', hmod%ibetap
            WRITE (*, fmt=fmt) 'iso beta z index:', hmod%ibetaz
            WRITE (*, fmt=fmt) 'beta gas:', hmod%gbeta   
            WRITE (*, fmt=fmt) 'beta gas z index:', hmod%gbetaz
         ELSE IF (hmod%HMx_mode == 6) THEN         
            WRITE (*, fmt=fmt) 'log10(T_heat) [K]:', log10(cosm%Theat)
            WRITE (*, fmt=fmt) 'alpha:', HMx_alpha(hmod%Mh, hmod, cosm)
            WRITE (*, fmt=fmt) 'beta:', HMx_beta(hmod%Mh, hmod, cosm)
            WRITE (*, fmt=fmt) 'epsilon:', HMx_eps(hmod, cosm)
            WRITE (*, fmt=fmt) 'epsilon2:', HMx_eps2(hmod, cosm)
            WRITE (*, fmt=fmt) 'Gamma:', HMx_Gamma(hmod%Mh, hmod, cosm)
            WRITE (*, fmt=fmt) 'log10(M0) [Msun/h]:', log10(HMx_M0(hmod, cosm))
            WRITE (*, fmt=fmt) 'A*:', HMx_Astar(hmod, cosm)
            WRITE (*, fmt=fmt) 'log10(T_WHIM) [K]:', log10(HMx_Twhim(hmod, cosm))          
            WRITE (*, fmt=fmt) 'c*:', HMx_cstar(hmod%Mh, hmod, cosm)
            WRITE (*, fmt=fmt) 'sigma*:', HMx_sstar(hmod, cosm)
            WRITE (*, fmt=fmt) 'log10(M*):', log10(HMx_Mstar(hmod, cosm))
            WRITE (*, fmt=fmt) 'frac cold:', HMx_fcold(hmod, cosm)
            WRITE (*, fmt=fmt) 'frac hot:', HMx_fhot(hmod, cosm)
            WRITE (*, fmt=fmt) 'eta:', HMx_eta(hmod, cosm)
            WRITE (*, fmt=fmt) 'isothermal beta:', HMx_ibeta(hmod%Mh, hmod, cosm)
            WRITE (*, fmt=fmt) 'beta gas frac:', HMx_gbeta(hmod, cosm)
         ELSE
            STOP 'PRINT_HALOMOD: Something went wrong'
         END IF
         WRITE (*, *) dashes

         ! HOD
         WRITE (*, *) 'HALOMODEL: HOD parameters'
         WRITE (*, *) dashes
         WRITE (*, fmt=fmt) 'log10(M_halo_min) [Msun/h]:', log10(hmod%mhalo_min)
         WRITE (*, fmt=fmt) 'log10(M_halo_max) [Msun/h]:', log10(hmod%mhalo_max)
         WRITE (*, fmt=fmt) 'log10(M_HI_min) [Msun/h]:', log10(hmod%HImin)
         WRITE (*, fmt=fmt) 'log10(M_HI_max) [Msun/h]:', log10(hmod%HImax)
         WRITE (*, *) dashes
         IF (hmod%halo_DMONLY == 5) WRITE (*, fmt=fmt) 'r_core [Mpc/h]:', hmod%rcore
         IF (hmod%sigma_conc /= 0.) WRITE (*, fmt=fmt) 'sigma(c)/c:', hmod%sigma_conc
         WRITE (*, *)

      END IF

   END SUBROUTINE print_halomod

   SUBROUTINE assign_init_halomod(ihm, a, hmod, cosm, verbose)

      ! Both assigns and initialises the halo model
      INTEGER, INTENT(INOUT) :: ihm
      REAL, INTENT(IN) :: a
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose

      CALL assign_halomod(ihm, hmod, verbose)
      CALL init_halomod(a, hmod, cosm, verbose)
      CALL print_halomod(hmod, cosm, verbose)

   END SUBROUTINE assign_init_halomod

   REAL FUNCTION integrate_g_mu(nu1, nu2, hmod)

      ! Integrate g(nu) between nu1 and nu2; this is the fraction of mass in the universe in haloes between nu1 and nu2
      ! Integrating this over all nu should give unity
      ! TODO: Explicity remove divergence by making functions where alpha*mu^(alpha-1) is multiplied by g(nu)
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model    
      REAL :: mu1, mu2
      REAL :: alpha
      INTEGER, PARAMETER :: iorder = iorder_integration

      ! Relation between nu and mu: nu=mu^alpha
      alpha = hmod%alpha_numu
      mu1 = nu1**(1./alpha)
      mu2 = nu2**(1./alpha)

      integrate_g_mu = integrate_hmod(mu1, mu2, g_mu, hmod, hmod%acc, iorder)

   END FUNCTION integrate_g_mu

   REAL FUNCTION integrate_g_nu(nu1, nu2, hmod)

      ! Integrate g(nu) between nu1 and nu2; this is the fraction of mass in the Universe in haloes between nu1 and nu2
      ! Previously called 'mass_interval'
      ! Integrating this over all nu should give unity
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      integrate_g_nu = integrate_hmod(nu1, nu2, g_nu, hmod, hmod%acc, iorder)

   END FUNCTION integrate_g_nu

   REAL FUNCTION integrate_g_nu_on_M(nu1, nu2, hmod)

      ! Integrate g(nu)/M between nu1 and nu2; this is the number density of haloes in the Universe between nu1 and nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      integrate_g_nu_on_M = integrate_hmod(nu1, nu2, g_nu_on_M, hmod, hmod%acc, iorder)

   END FUNCTION integrate_g_nu_on_M

   REAL FUNCTION integrate_gb_nu(nu1, nu2, hmod)

      ! Integrate b(nu)*g(nu) between nu1 and nu2
      ! Previously called 'bias_interval'
      ! This is the unnormalised mass-weighted mean bias in the range nu1 to nu2
      ! Integrating this over all nu should give unity
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      integrate_gb_nu = integrate_hmod(nu1, nu2, gb_nu, hmod, hmod%acc, iorder)

   END FUNCTION integrate_gb_nu

   REAL FUNCTION integrate_gb_nu_on_M(nu1, nu2, hmod)

      ! Integrate b(nu)*g(nu) between nu1 and nu2
      ! Previously called 'bias_interval'
      ! This is the unnormalised mass-weighted mean bias in the range nu1 to nu2
      ! Integrating this over all nu should give unity
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      integrate_gb_nu_on_M = integrate_hmod(nu1, nu2, gb_nu_on_M, hmod, hmod%acc, iorder)

   END FUNCTION integrate_gb_nu_on_M

   REAL FUNCTION integrate_nug_nu(nu1, nu2, hmod)

      ! Integrate nu*g(nu) between nu1 and nu2
      ! Previously called 'nu interval'
      ! This is the unnormalised mass-weighted nu in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      integrate_nug_nu = integrate_hmod(nu1, nu2, nug_nu, hmod, hmod%acc, iorder)

   END FUNCTION integrate_nug_nu

   REAL FUNCTION integrate_Mg_nu(nu1, nu2, hmod)

      ! Integrate nu*g(nu) between nu1 and nu2
      ! Previously called 'nu interval'
      ! This is the unnormalised mass-weighted nu in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      integrate_Mg_nu = integrate_hmod(nu1, nu2, Mg_nu, hmod, hmod%acc, iorder)

   END FUNCTION integrate_Mg_nu

   REAL FUNCTION integrate_nug_nu_on_M(nu1, nu2, hmod)

      ! Integrate nu*g(nu) between nu1 and nu2
      ! Previously called 'nu interval'
      ! This is the unnormalised mass-weighted nu in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      integrate_nug_nu_on_M = integrate_hmod(nu1, nu2, nug_nu_on_M, hmod, hmod%acc, iorder)

   END FUNCTION integrate_nug_nu_on_M

   REAL FUNCTION mean_bias_number_weighted(nu1, nu2, hmod)

      ! Calculate the mean number-weighted bias in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      mean_bias_number_weighted = integrate_gb_nu_on_M(nu1, nu2, hmod)/integrate_g_nu_on_M(nu1, nu2, hmod)

   END FUNCTION mean_bias_number_weighted

   REAL FUNCTION mean_nu_number_weighted(nu1, nu2, hmod)

      ! Calculate the mean number-weighted nu in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      mean_nu_number_weighted = integrate_nug_nu_on_M(nu1, nu2, hmod)/integrate_g_nu_on_M(nu1, nu2, hmod)

   END FUNCTION mean_nu_number_weighted

   REAL FUNCTION mean_halo_mass_number_weighted(nu1, nu2, hmod)

      ! Calculate the mean number-weighted halo mass in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      mean_halo_mass_number_weighted = integrate_g_nu(nu1, nu2, hmod)/integrate_g_nu_on_M(nu1, nu2, hmod)

   END FUNCTION mean_halo_mass_number_weighted

   REAL FUNCTION mean_bias_mass_weighted(nu1, nu2, hmod)

      ! Calculate the mean mass-weighted bias in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      mean_bias_mass_weighted = integrate_gb_nu(nu1, nu2, hmod)/integrate_g_nu(nu1, nu2, hmod)

   END FUNCTION mean_bias_mass_weighted

   REAL FUNCTION mean_nu_mass_weighted(nu1, nu2, hmod)

      ! Calculate the mean mass-weighted nu in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      mean_nu_mass_weighted = integrate_nug_nu(nu1, nu2, hmod)/integrate_g_nu(nu1, nu2, hmod)

   END FUNCTION mean_nu_mass_weighted

   REAL FUNCTION mean_halo_mass_mass_weighted(nu1, nu2, hmod)

      ! Calculate the mean mass-weighted halo mass in the range nu1 to nu2
      REAL, INTENT(IN) :: nu1, nu2         ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      mean_halo_mass_mass_weighted = integrate_Mg_nu(nu1, nu2, hmod)/integrate_g_nu(nu1, nu2, hmod)

   END FUNCTION mean_halo_mass_mass_weighted

   REAL FUNCTION mean_halo_number_density(nu1, nu2, hmod, cosm)

      ! Calculate N(m) where N is the number density of haloes above mass m
      ! Obtained by integrating the mass function
      REAL, INTENT(IN) :: nu1, nu2           ! Range in nu
      TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      INTEGER, PARAMETER :: iorder = iorder_integration ! Order for integration

      mean_halo_number_density = integrate_hmod(nu1, nu2, g_nu_on_M, hmod, hmod%acc, iorder)
      mean_halo_number_density = mean_halo_number_density*comoving_matter_density(cosm)

   END FUNCTION mean_halo_number_density

   FUNCTION halo_type(i)

      ! Name function for halo types
      ! TODO: This must be able to be combined with set_halo_type
      ! TODO: Can this be in the header? Is this even necessary?
      CHARACTER(len=256) :: halo_type
      INTEGER :: i

      halo_type = ''
      IF (i == field_dmonly) halo_type = 'DMONLY'
      IF (i == field_matter) halo_type = 'Matter'
      IF (i == field_cdm) halo_type = 'CDM'
      IF (i == field_gas) halo_type = 'Gas'
      IF (i == field_stars) halo_type = 'Stars'
      IF (i == field_bound_gas) halo_type = 'Bound gas'
      IF (i == field_ejected_gas) halo_type = 'Free gas'
      IF (i == field_electron_pressure) halo_type = 'Electron pressure'
      IF (i == field_void) halo_type = 'Void'
      IF (i == field_compensated_void) halo_type = 'Compensated void'
      IF (i == field_central_galaxies) halo_type = 'Central galaxies or haloes'
      IF (i == field_satellite_galaxies) halo_type = 'Satellite galaxies'
      IF (i == field_galaxies) halo_type = 'Galaxies'
      IF (i == field_HI) halo_type = 'HI'
      IF (i == field_cold_bound_gas) halo_type = 'Cold gas'
      IF (i == field_hot_bound_gas) halo_type = 'Hot gas'
      IF (i == field_normal_bound_gas) halo_type = 'Static gas'
      IF (i == field_central_stars) halo_type = 'Central stars'
      IF (i == field_satellite_stars) halo_type = 'Satellite stars'
      IF (i == field_CIB_353) halo_type = 'CIB 353 GHz'
      IF (i == field_CIB_545) halo_type = 'CIB 545 GHz'
      IF (i == field_CIB_857) halo_type = 'CIB 857 GHz'
      IF (i == field_halo_11p0_11p5) halo_type = 'Haloes 10^11.0 -> 10^11.5'
      IF (i == field_halo_11p5_12p0) halo_type = 'Haloes 10^11.5 -> 10^12.0'
      IF (i == field_halo_12p0_12p5) halo_type = 'Haloes 10^12.0 -> 10^12.5'
      IF (i == field_halo_12p5_13p0) halo_type = 'Haloes 10^12.5 -> 10^13.0'
      IF (i == field_halo_13p0_13p5) halo_type = 'Haloes 10^13.0 -> 10^13.5'
      IF (i == field_halo_13p5_14p0) halo_type = 'Haloes 10^13.5 -> 10^14.0'
      IF (i == field_halo_14p0_14p5) halo_type = 'Haloes 10^14.0 -> 10^14.5'
      IF (i == field_halo_14p5_15p0) halo_type = 'Haloes 10^14.5 -> 10^15.0'
      IF (i == field_neutrino) halo_type = 'Neutrino'
      IF (halo_type == '') STOP 'HALO_TYPE: Error, i not specified correctly'

   END FUNCTION halo_type

   SUBROUTINE set_halo_type(ip)

      ! Set the halo types
      ! TODO: This must be able to be combined with halo_type
      INTEGER, INTENT(OUT) :: ip
      INTEGER :: i

      WRITE (*, *) 'SET_HALO_TYPE: Choose halo type'
      WRITE (*, *) '==============================='
      DO i = 1, field_n
         WRITE (*, fmt='(I3,A3,A26)') i, '- ', TRIM(halo_type(i))
      END DO
      READ (*, *) ip
      WRITE (*, *) '==============================='
      WRITE (*, *)

      IF (ip < 1 .OR. ip > field_n) STOP 'SET_HALO_TYPE: Error, you have chosen a halo type that does not exist'

   END SUBROUTINE set_halo_type

   SUBROUTINE calculate_HMcode(k, a, Pk, nk, na, cosm, version)

      ! Get the HMcode prediction for a cosmology for a range of k and a
      INTEGER, INTENT(IN) :: nk                  ! Number of wavenumbers
      INTEGER, INTENT(IN) :: na                  ! Number of scale factors
      REAL, INTENT(IN) :: k(nk)                  ! Array of wavenumbers [h/Mpc]
      REAL, INTENT(IN) :: a(na)                  ! Array of scale factors
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      TYPE(cosmology), INTENT(INOUT) :: cosm     ! Cosmology
      INTEGER, OPTIONAL, INTENT(IN) :: version   ! The ihm corresponding to the HMcode version
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :), pow_1h(:, :)

      CALL calculate_HMcode_full(k, a, pow_li, pow_2h, pow_1h, Pk, nk, na, cosm, version)

   END SUBROUTINE calculate_HMcode

   SUBROUTINE calculate_HMcode_full(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, cosm, version)

      ! Get the HMcode prediction for a cosmology for a range of k and a
      INTEGER, INTENT(IN) :: nk                  ! Number of wavenumbers
      INTEGER, INTENT(IN) :: na                  ! Number of scale factors
      REAL, INTENT(IN) :: k(nk)                  ! Array of wavenumbers [h/Mpc]
      REAL, INTENT(IN) :: a(na)                  ! Array of scale factors
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_li(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_2h(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_1h(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_hm(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      TYPE(cosmology), INTENT(INOUT) :: cosm     ! Cosmology
      INTEGER, OPTIONAL, INTENT(IN) :: version   ! The ihm corresponding to the HMcode version
      INTEGER :: ihm

      IF(.NOT. present(version)) THEN
         ihm = HMcode2016
      ELSE
         ihm = version
      END IF

      IF(.NOT. is_HMcode(ihm)) STOP 'CALCULATE_HMCODE_FULL: Error, your HMcode version is not recognized'   

      CALL calculate_halomod_full(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, cosm, ihm)

   END SUBROUTINE calculate_HMcode_full

   LOGICAL FUNCTION is_HMcode(i)

      ! Returns TRUE if integer corresponds to a recognised HMcode version
      INTEGER, INTENT(IN) :: i

      IF(is_in_array(i, [&
         HMcode2015, &
         HMcode2015_CAMB, &
         HMcode2016, &
         HMcode2016_CAMB, &
         HMcode2016_neutrinofix, &
         HMcode2016_neutrinofix_CAMB, &
         HMcode2018, &
         HMcode2019, &
         HMcode2020, &
         HMcode2020_CAMB, &
         HMcode2020_feedback, &
         HMcode2020_feedback_CAMB &
         ])) THEN
            is_HMcode = .TRUE.
         ELSE
            is_HMcode = .FALSE.
      END IF

   END FUNCTION is_HMcode

   SUBROUTINE calculate_halomod(k, a, Pk, nk, na, cosm, ihm)

      ! Get the halo model prediction for matter--matter for cosmology for a range of k and a
      ! Assumes DMONLY halo profiles etc.
      INTEGER, INTENT(IN) :: nk                  ! Number of wavenumbers
      INTEGER, INTENT(IN) :: na                  ! Number of scale factors
      REAL, INTENT(IN) :: k(nk)                  ! Array of wavenumbers [h/Mpc]
      REAL, INTENT(IN) :: a(na)                  ! Array of scale factors
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      TYPE(cosmology), INTENT(INOUT) :: cosm     ! Cosmology
      INTEGER, INTENT(INOUT) :: ihm              ! The ihm corresponding to the HMcode version
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :), pow_1h(:, :)

      CALL calculate_halomod_full(k, a, pow_li, pow_2h, pow_1h, Pk, nk, na, cosm, ihm)

   END SUBROUTINE calculate_halomod

   SUBROUTINE calculate_halomod_response(k, a, Pk, nk, na, cosm, ihm_num, ihm_den, ihm_base)

      ! Get the halo model prediction for matter--matter for cosmology for a range of k and a
      ! Assumes DMONLY halo profiles etc.
      INTEGER, INTENT(IN) :: nk                  ! Number of wavenumbers
      INTEGER, INTENT(IN) :: na                  ! Number of scale factors
      REAL, INTENT(IN) :: k(nk)                  ! Array of wavenumbers [h/Mpc]
      REAL, INTENT(IN) :: a(na)                  ! Array of scale factors
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      TYPE(cosmology), INTENT(INOUT) :: cosm     ! Cosmology
      INTEGER, INTENT(IN) :: ihm_num             ! The ihm corresponding to the numerator
      INTEGER, INTENT(IN) :: ihm_den             ! The ihm corresponding to the denominator
      INTEGER, INTENT(IN) :: ihm_base            ! The ihm corresponding to the HMcode version (probably)
      REAL, ALLOCATABLE :: Pk_num(:, :), Pk_den(:, :)
      INTEGER :: ihm

      ihm = ihm_num
      CALL calculate_halomod(k, a, Pk_num, nk, na, cosm, ihm)

      ihm = ihm_den
      CALL calculate_halomod(k, a, Pk_den, nk, na, cosm, ihm)

      ihm = ihm_base
      CALL calculate_halomod(k, a, Pk, nk, na, cosm, ihm)

      Pk = Pk*Pk_num/Pk_den

   END SUBROUTINE calculate_halomod_response

   SUBROUTINE calculate_halomod_full(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, cosm, ihm)

      ! Get the halo model prediction for matter--matter for cosmology for a range of k and a
      ! Assumes DMONLY halo profiles etc.
      ! TODO: Can I directly put in 2D arrays into routine that expects 4D if some dimensions are length 1?
      INTEGER, INTENT(IN) :: nk                  ! Number of wavenumbers
      INTEGER, INTENT(IN) :: na                  ! Number of scale factors
      REAL, INTENT(IN) :: k(nk)                  ! Array of wavenumbers [h/Mpc]
      REAL, INTENT(IN) :: a(na)                  ! Array of scale factors
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_li(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_2h(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_1h(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_hm(:, :) ! Output power array, note that this is Delta^2(k), not P(k)
      TYPE(cosmology), INTENT(INOUT) :: cosm     ! Cosmology
      INTEGER, INTENT(INOUT) :: ihm              ! The ihm corresponding to the HMcode version
      REAL, ALLOCATABLE :: pows_2h(:, :, :, :), pows_1h(:, :, :, :), pows_hm(:, :, :, :)
      INTEGER, PARAMETER :: field(1) = [field_DMONLY]
      INTEGER, PARAMETER :: nf = 1

      CALL calculate_HMx_full(field, k, a, pow_li, pows_2h, pows_1h, pows_hm, nf, nk, na, cosm, ihm)

      ALLOCATE(pow_2h(nk, na), pow_1h(nk, na), pow_hm(nk, na))
      pow_2h = pows_2h(1, 1, :, :)
      pow_1h = pows_1h(1, 1, :, :)
      pow_hm = pows_hm(1, 1, :, :)

   END SUBROUTINE calculate_halomod_full

   SUBROUTINE calculate_HMx(ifield, k, a, pow, nf, nk, na, cosm, ihm)

      ! Public facing function, calculates the halo model power for k and a range
      ! TODO: Change (:,:,k,a) to (k,a,:,:) for speed?
      INTEGER, INTENT(IN) :: nf         ! Number of different fields
      INTEGER, INTENT(IN) :: ifield(nf) ! Indices for different fields
      INTEGER, INTENT(IN) :: nk         ! Number of k points
      REAL, INTENT(IN) :: k(nk)         ! k array [h/Mpc]
      INTEGER, INTENT(IN) :: na         ! Number of a points
      REAL, INTENT(IN) :: a(na)         ! a array
      REAL, ALLOCATABLE, INTENT(OUT) :: pow(:, :, :, :) ! Pow(f1,f2,k,a)
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      INTEGER, INTENT(INOUT) :: ihm
      REAL, ALLOCATABLE :: pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :)

      CALL calculate_HMx_full(ifield, k, a, pow_li, pow_2h, pow_1h, pow, nf, nk, na, cosm, ihm)

   END SUBROUTINE calculate_HMx

   SUBROUTINE calculate_HMx_full(ifield, k, a, pow_li, pow_2h, pow_1h, pow_hm, nf, nk, na, cosm, ihm)

      ! Public facing function, calculates the halo model power for k and a range
      ! TODO: Change (:,:,k,a) to (k,a,:,:) for speed?
      INTEGER, INTENT(IN) :: nf         ! Number of different fields
      INTEGER, INTENT(IN) :: ifield(nf) ! Indices for different fields
      INTEGER, INTENT(IN) :: nk         ! Number of k points
      REAL, INTENT(IN) :: k(nk)         ! k array [h/Mpc]
      INTEGER, INTENT(IN) :: na         ! Number of a points
      REAL, INTENT(IN) :: a(na)         ! a array
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_li(:, :)       ! Pow(k,a)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_2h(:, :, :, :) ! Pow(f1,f2,k,a)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_1h(:, :, :, :) ! Pow(f1,f2,k,a)
      REAL, ALLOCATABLE, INTENT(OUT) :: pow_hm(:, :, :, :) ! Pow(f1,f2,k,a)
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      INTEGER, INTENT(INOUT) :: ihm
      INTEGER :: ia
      TYPE(halomod) :: hmod
      LOGICAL :: verbose
      LOGICAL, PARAMETER :: verbose_here = .FALSE.

      CALL assign_halomod(ihm, hmod, verbose_here)

      ALLOCATE(pow_li(nk, na), pow_2h(nf, nf, nk, na), pow_1h(nf, nf, nk, na), pow_hm(nf, nf, nk, na))

      DO ia = 1, na

         ! So as not to make it overly verbose
         IF (ia == na) THEN
            verbose = verbose_here
         ELSE
            verbose = .FALSE.
         END IF 

         CALL init_halomod(a(ia), hmod, cosm, verbose)

         CALL calculate_HMx_a(ifield, nf, k, nk, &
            pow_li(:, ia), &
            pow_2h(:, :, :, ia), &
            pow_1h(:, :, :, ia), &
            pow_hm(:, :, :, ia), &
            hmod, cosm, verbose)

         CALL print_halomod(hmod, cosm, verbose)

      END DO

   END SUBROUTINE calculate_HMx_full

   SUBROUTINE calculate_HMx_old(ifield, nf, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

      ! Public facing function, calculates the halo model power for k and a range
      ! TODO: Change (:, :, k, a) to (k, a, :, :) for speed?
      ! TODO: Change name to calculate_HMx_manual? Tilman needs this routine still so do not remove it
      INTEGER, INTENT(IN) :: nf               ! Number of different fields
      INTEGER, INTENT(IN) :: ifield(nf)       ! Indices for different fields
      INTEGER, INTENT(IN) :: nk               ! Number of k points
      REAL, INTENT(IN) :: k(nk)               ! k array [h/Mpc]
      INTEGER, INTENT(IN) :: na               ! Number of a points
      REAL, INTENT(IN) :: a(na)               ! a array
      REAL, INTENT(OUT) :: pow_li(:, :)       ! Pow(k, a)
      REAL, INTENT(OUT) :: pow_2h(:, :, :, :) ! Pow(f1, f2, k, a)
      REAL, INTENT(OUT) :: pow_1h(:, :, :, :) ! Pow(f1, f2, k, a)
      REAL, INTENT(OUT) :: pow_hm(:, :, :, :) ! Pow(f1, f2, k, a)
      TYPE(halomod), INTENT(INOUT) :: hmod    ! Halo model
      TYPE(cosmology), INTENT(INOUT) :: cosm  ! Cosmology
      LOGICAL, INTENT(IN) :: verbose          ! Verbosity
      REAL :: z
      INTEGER :: i, j
      LOGICAL :: verbose2

      ! To avoid splurge of stuff printed to screen
      verbose2 = verbose

      ! Do the halo-model calculation by looping over scale factor index
      ! TODO: Why does this loop backwards. Is it to do with the write-to-screen command?
      DO i = na, 1, -1

         z = redshift_a(a(i))
         CALL init_halomod(a(i), hmod, cosm, verbose2)
         CALL print_halomod(hmod, cosm, verbose2)
         CALL calculate_HMx_a(ifield, nf, k, nk, pow_li(:,i), pow_2h(:,:,:,i), pow_1h(:,:,:,i), pow_hm(:,:,:,i),&
            hmod, cosm, verbose2)

         IF (i == na .and. verbose) THEN
            WRITE (*, *) 'CALCULATE_HMx: Doing calculation'
            DO j = 1, nf
               WRITE (*, *) 'CALCULATE_HMx: Haloes:', ifield(j), TRIM(halo_type(ifield(j)))
            END DO
            WRITE (*, *) '======================================='
            WRITE (*, *) '                            a         z'
            WRITE (*, *) '======================================='
         END IF

         IF (verbose) WRITE (*, fmt='(A15,I5,2F10.3)') 'CALCULATE_HMx:', i, a(i), z
         verbose2 = .FALSE.

      END DO

      IF (verbose) THEN
         WRITE (*, *) '======================================='
         WRITE (*, *) 'CALCULATE_HMx: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE calculate_HMx_old

   SUBROUTINE calculate_HMx_a(ifield, nf, k, nk, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)

      ! Calculate halo model Pk(a) for a k range

      INTEGER, INTENT(IN) :: nf
      INTEGER, INTENT(IN) :: ifield(nf)
      INTEGER, INTENT(IN) :: nk
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(OUT) :: pow_li(nk)
      REAL, INTENT(OUT) :: pow_2h(nf, nf, nk)
      REAL, INTENT(OUT) :: pow_1h(nf, nf, nk)
      REAL, INTENT(OUT) :: pow_hm(nf, nf, nk)
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      REAL :: powg_2h(nk), powg_1h(nk), powg_hm(nk)
      REAL, ALLOCATABLE :: upow_2h(:, :, :), upow_1h(:, :, :), upow_hm(:, :, :)
      REAL, ALLOCATABLE :: wk0(:, :)
      REAL :: wk0_den(hmod%n, 1), wk0_base(hmod%n, 1)
      REAL :: base_2h(nk), base_1h(nk), base_hm(nk)
      INTEGER :: i, j, ii, jj, match(nf), nnf
      INTEGER :: ihm
      INTEGER, ALLOCATABLE :: iifield(:)
      TYPE(halomod) :: hmod_base, hmod_den
      INTEGER, PARAMETER :: dmonly(1) = field_dmonly ! Needed because it needs to be an array(1)
      !INTEGER, PARAMETER :: ihm_den = 77
      !INTEGER :: ihmcode = HMcode2016 ! Would like to be a parameter

      ! Make a new indexing scheme for only the unique arrays
      CALL unique_index(ifield, nf, iifield, nnf, match)
      ALLOCATE (upow_2h(nnf, nnf, nk), upow_1h(nnf, nnf, nk), upow_hm(nnf, nnf, nk))

      ! Write to screen
      IF (verbose) THEN
         DO i = 1, nnf
            WRITE (*, *) 'CALCULATE_HMX_A: Halo type:', iifield(i), TRIM(halo_type(iifield(i)))
         END DO
         WRITE (*, *) 'CALCULATE_HMX_A: k min [h/Mpc]:', REAL(k(1))
         WRITE (*, *) 'CALCULATE_HMX_A: k max [h/Mpc]:', REAL(k(nk))
         WRITE (*, *) 'CALCULATE_HMX_A: number of k:', nk
         WRITE (*, *) 'CALCULATE_HMX_A: a:', REAL(hmod%a)
         WRITE (*, *) 'CALCULATE_HMX_A: z:', REAL(hmod%z)
         WRITE (*, *) 'CALCULATE_HMX_A: Calculating halo-model power spectrum'
         WRITE (*, *)
      END IF

      ! Do an HMcode calculation for multiplying the response
      ! Also calculate the k->0 window functions
      IF (hmod%response_baseline .NE. 0) THEN
         ihm = hmod%response_baseline
         CALL assign_init_halomod(ihm, hmod%a, hmod_base, cosm, verbose=.FALSE.)
         CALL init_windows(zero, dmonly, 1, wk0_base, hmod_base%n, hmod_base, cosm)
         IF (hmod%response_denominator .NE. 0) THEN
            ihm = hmod%response_denominator
            CALL assign_init_halomod(ihm, hmod%a, hmod_den, cosm, verbose=.FALSE.)
            CALL init_windows(zero, dmonly, 1, wk0_den, hmod_den%n, hmod_den, cosm)
         ELSE
            CALL init_windows(zero, dmonly, 1, wk0_den, hmod%n, hmod, cosm)
         END IF
      END IF

      ! Get the k->0 window functions
      ALLOCATE(wk0(hmod%n, nnf))
      CALL init_windows(zero, ifield, nnf, wk0, hmod%n, hmod, cosm)

      ! Loop over k values
      ! TODO: add OMP support properly. What is private and what is shared? CHECK THIS!
!!$OMP PARALLEL DO DEFAULT(SHARED)!, private(k,plin,pow_2h,pow_1h,pow,pow_lin)
!!$OMP PARALLEL DO DEFAULT(PRIVATE)
!!$OMP PARALLEL DO FIRSTPRIVATE(nk,cosm,compute_p_lin,k,a,pow_lin,plin,itype1,itype2,z,pow_2h,pow_1h,pow,hmod)
!!$OMP PARALLEL DO
      DO i = 1, nk

         ! Get the linear power
         pow_li(i) = plin(k(i), hmod%a, flag_matter, cosm)

         ! Do the halo model calculation
         CALL calculate_HMx_ka(iifield, wk0, nnf, k(i), pow_li(i), upow_2h(:, :, i), upow_1h(:, :, i), upow_hm(:, :, i), hmod, cosm)

         IF (hmod%response_baseline .NE. 0) THEN

            ! If doing a response then calculate a DMONLY prediction for both the current halomodel and HMcode
            CALL calculate_HMx_ka(dmonly, wk0_base, 1, k(i), pow_li(i), base_2h(i), base_1h(i), base_hm(i), hmod_base, cosm)
            IF (hmod%response_denominator .NE. 0) THEN
               CALL calculate_HMx_ka(dmonly, wk0_den, 1, k(i), pow_li(i), powg_2h(i), powg_1h(i), powg_hm(i), hmod_den, cosm)
            ELSE
               CALL calculate_HMx_ka(dmonly, wk0_den, 1, k(i), pow_li(i), powg_2h(i), powg_1h(i), powg_hm(i), hmod, cosm)
            END IF  

            IF (hmod%response_matter_only) THEN

               ! Only apply response to matter fields
               ! TODO: Slow array accessing here? Change ii and jj loop order?
               DO ii = 1, nf
                  DO jj = 1, nf
                     IF(is_matter_field(iifield(ii)) .AND. is_matter_field(iifield(jj))) THEN
                        upow_2h(ii, jj, i) = upow_2h(ii, jj, i)*base_2h(i)/powg_2h(i) ! Two-halo response times HMcode
                        upow_1h(ii, jj, i) = upow_1h(ii, jj, i)*base_1h(i)/powg_1h(i) ! One-halo response times HMcode
                        upow_hm(ii, jj, i) = upow_hm(ii, jj, i)*base_hm(i)/powg_hm(i) ! Full model response times HMcode
                     END IF
                  END DO
               END DO

            ELSE

               ! Calculate the response for everything
               upow_2h(:, :, i) = upow_2h(:, :, i)*base_2h(i)/powg_2h(i) ! Two-halo response times HMcode
               upow_1h(:, :, i) = upow_1h(:, :, i)*base_1h(i)/powg_1h(i) ! One-halo response times HMcode
               upow_hm(:, :, i) = upow_hm(:, :, i)*base_hm(i)/powg_hm(i) ! Full model response times HMcode

            END IF

         END IF

      END DO
!!$OMP END PARALLEL DO

      ! Now fill the full arrays from the unique arrays
      DO i = 1, nf
         DO j = 1, nf
            ii = match(i)
            jj = match(j)
            pow_1h(i, j, :) = upow_1h(ii, jj, :)
            pow_2h(i, j, :) = upow_2h(ii, jj, :)
            pow_hm(i, j, :) = upow_hm(ii, jj, :)
         END DO
      END DO

   END SUBROUTINE calculate_HMx_a

   LOGICAL FUNCTION is_matter_field(i)

      ! TRUE if the input field is matter and corresponds to dmonly, matter, CDM, gas or star (or consituents)
      INTEGER, INTENT(IN) :: i

      IF (is_in_array(i, [&
            field_dmonly, &
            field_matter, &
            field_cdm, &
            field_gas, &
            field_bound_gas, &
            field_normal_bound_gas, &
            field_hot_bound_gas, &
            field_cold_bound_gas, &
            field_ejected_gas, &      
            field_stars, &    
            field_central_stars, &
            field_satellite_stars, &
            field_neutrino & 
            ])) THEN
         is_matter_field = .TRUE.
      ELSE
         is_matter_field = .FALSE.
      END IF

   END FUNCTION is_matter_field

   SUBROUTINE calculate_HMx_ka(ifield, wk0, nf, k, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm)

      ! Gets the two- and one-halo terms and combines them
      ! TODO: Include profile scatter in two-halo term
      INTEGER, INTENT(IN) :: nf
      INTEGER, INTENT(IN) :: ifield(nf)
      REAL, INTENT(IN) :: wk0(:, :)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: pow_li
      REAL, INTENT(OUT) :: pow_2h(nf, nf)
      REAL, INTENT(OUT) :: pow_1h(nf, nf)
      REAL, INTENT(OUT) :: pow_hm(nf, nf)
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: wk(hmod%n, nf), wk2(hmod%n, 2), wk_product(hmod%n)
      INTEGER :: f1, f2, ih(2)

      ! Calls expressions for two- and one-halo terms and then combines to form the full power spectrum
      IF (k == 0.) THEN

         ! This should really never be called for k=0
         pow_2h = 0.
         pow_1h = 0.

      ELSE

         ! Get the window functions
         CALL init_windows(k, ifield, nf, wk, hmod%n, hmod, cosm)

         ! Loop over fields and get the one-halo term for each pair
         DO f1 = 1, nf
            DO f2 = f1, nf
               IF (hmod%sigma_conc == 0.) THEN
                  wk_product = wk(:, f1)*wk(:, f2)
               ELSE
                  ih(1) = ifield(f1)
                  ih(2) = ifield(f2)
                  wk_product = wk_product_scatter(hmod%n, ih, k, hmod, cosm)
               END IF
               pow_1h(f1, f2) = p_1h(wk_product, hmod%n, k, hmod, cosm)
            END DO
         END DO

         ! If linear theory is used for two-halo term we must recalculate the window functions
         ! For the two-halo term we need the k->0 limit
         IF (is_in_array(hmod%ip2h, [1, 3, 5, 6])) THEN
            wk = wk0
         ELSE
            ! If we are worrying about halo profiles that somehow change between one- and two-halo evaulations then we need this
            ! This applies to ejected halo gas, but also to massive neutrinos and simple baryon feedback models
            CALL add_smooth_component_to_windows(ifield, nf, wk, hmod%n, hmod, cosm)
         END IF

         ! Get the two-halo term
         DO f1 = 1, nf
            DO f2 = f1, nf
               ih(1) = ifield(f1)
               ih(2) = ifield(f2)
               wk2(:, 1) = wk(:, f1)
               wk2(:, 2) = wk(:, f2)
               pow_2h(f1, f2) = p_2h(ih, wk2, hmod%n, k, pow_li, hmod, cosm)
            END DO
         END DO

      END IF

      ! Loop over fields and get the total halo-model power
      DO f1 = 1, nf
         DO f2 = f1, nf
            pow_hm(f1, f2) = p_hm(k, pow_2h(f1, f2), pow_1h(f1, f2), hmod)
         END DO
      END DO

      ! Construct symmetric parts using ij=ji symmetry of spectra
      DO f1 = 1, nf
         DO f2 = f1, nf
            pow_1h(f2, f1) = pow_1h(f1, f2)
            pow_2h(f2, f1) = pow_2h(f1, f2)
            pow_hm(f2, f1) = pow_hm(f1, f2)
         END DO
      END DO

   END SUBROUTINE calculate_HMx_ka

   SUBROUTINE init_windows(k, fields, nf, wk, nm, hmod, cosm)

      ! Fill the window functions for all the different fields
      REAL, INTENT(IN) :: k
      INTEGER, INTENT(IN) :: nf
      INTEGER, INTENT(IN) :: fields(nf)
      INTEGER, INTENT(IN) :: nm
      REAL, INTENT(OUT) :: wk(nm, nf)
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i, j
      REAL :: m, rv, c, rs, nu, eta
      INTEGER :: i_all, i_cdm, i_gas, i_sta, i_neu
      LOGICAL :: quick_matter
      LOGICAL, PARAMETER :: real_space = .FALSE. ! Fourier profiles
      LOGICAL, PARAMETER :: use_quick_matter = .TRUE.

      ! This should be set to false initially
      quick_matter = .FALSE.

      IF (use_quick_matter) THEN

         IF (repeated_entries(fields)) STOP 'INIT_WINDOWS: Error, repeated fields'

         ! Get the array positions corresponding to all, cdm, gas, stars, (neutrinos) if they exist
         i_all = array_position(field_matter, fields)
         i_cdm = array_position(field_cdm, fields)
         i_gas = array_position(field_gas, fields)
         i_sta = array_position(field_stars, fields)
         i_neu = array_position(field_neutrino, fields)

         ! If all, cdm, gas and stars exist then activate the quick-matter mode
         IF (cosm%f_nu == 0.) THEN
            IF ((i_all /= 0) .AND. (i_cdm /= 0) .AND. (i_gas /= 0) .AND. (i_sta /= 0)) THEN
               quick_matter = .TRUE.
            END IF
         ELSE 
            IF ((i_all /= 0) .AND. (i_cdm /= 0) .AND. (i_gas /= 0) .AND. (i_sta /= 0) .AND. (i_neu /= 0)) THEN
               quick_matter = .TRUE.
            END IF
         END IF

      END IF

      ! Get eta
      eta = hmod%HMcode_eta

      ! Calculate the halo window functions for each field
      DO j = 1, nf

         IF (quick_matter .AND. j == i_all) CYCLE

         ! Loop over masses to fill window-function array
         DO i = 1, nm

            m = hmod%m(i)
            rv = hmod%rv(i)
            c = hmod%c(i)
            rs = rv/c
            nu = hmod%nu(i)

            wk(i, j) = win(real_space, fields(j), k*nu**eta, m, rv, rs, hmod, cosm)

         END DO

      END DO

      ! If quick-matter mode is active then create the total matter window by summing contributions
      IF (use_quick_matter .AND. quick_matter) THEN       
            wk(:, i_all) = wk(:, i_cdm)+wk(:, i_gas)+wk(:, i_sta)
         IF(cosm%f_nu /= 0.) THEN
            wk(:, i_all) = wk(:, i_cdm)+wk(:,i_neu)
         END IF
      END IF

   END SUBROUTINE init_windows

   SUBROUTINE add_smooth_component_to_windows(fields, nf, wk, nm, hmod, cosm)

      ! Refills the window functions for the two-halo term if this is necessary
      ! This is for contributions due to unbound gas, and the effect of this on electron pressure
      ! TODO: Have I inculded the halo bias corresponding to the free component correctly?
      INTEGER, INTENT(IN) :: nf
      INTEGER, INTENT(IN) :: fields(nf)
      INTEGER, INTENT(IN) :: nm
      REAL, INTENT(INOUT) :: wk(nm, nf)
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: m, rv, c, rs, nu, fc, pc
      REAL :: rho0, T0
      INTEGER :: i, i_matter, i_gas, i_freegas, i_pressure, i_dmonly, i_neutrino

      ! Check for repeated fields
      IF (repeated_entries(fields)) STOP 'ADD_SMOOTH_COMPONENT_TO_WINDOWS: Error, repeated fields'

      ! If working with DMONLY haloes and an effective neutrino model
      IF (hmod%DMONLY_baryon_recipe) THEN
         i_dmonly = array_position(field_dmonly, fields)
         IF (i_dmonly /= 0) THEN
            DO i = 1, hmod%n
               wk(i, i_dmonly) = unbaryonify_wk(wk(i, i_dmonly), hmod%m(i), hmod, cosm)
            END DO
         END IF
      ELSE IF (hmod%DMONLY_neutrino_halo_mass_correction) THEN
         i_dmonly = array_position(field_dmonly, fields)
         IF (i_dmonly /= 0) wk(:, i_dmonly) = wk(:, i_dmonly)/(1.-cosm%f_nu)
      END IF

      ! Add in neutrinos for matter/neutrino profiles
      IF (hmod%halo_neutrino == 1 .AND. cosm%f_nu /= 0.) THEN
         i_matter = array_position(field_matter, fields)
         i_neutrino = array_position(field_neutrino, fields)
         IF (i_matter == 0 .AND. i_neutrino == 0) THEN
            ! Do nothing
         ELSE
            ! Loop over mass and apply corrections
            DO i = 1, hmod%n
               m = hmod%m(i)
               fc = halo_neutrino_fraction(m, hmod, cosm)*m/comoving_matter_density(cosm)
               IF(i_matter /= 0) wk(i, i_matter) = wk(i, i_matter)+fc
               IF(i_neutrino /= 0) wk(i, i_neutrino) = wk(i, i_neutrino)+fc
            END DO
         END IF
      END IF

      IF (hmod%halo_ejected_gas == 7) THEN

         ! Get the array positions corresponding to all, gas, pressure if they exist
         i_matter = array_position(field_matter, fields)
         i_gas = array_position(field_gas, fields)
         i_freegas = array_position(field_ejected_gas, fields)
         i_pressure = array_position(field_electron_pressure, fields)

         IF (i_matter == 0 .AND. i_gas == 0 .AND. i_pressure == 0 .AND. i_freegas == 0) THEN

            ! Do nothing

         ELSE

            ! Loop over mass and apply corrections
            DO i = 1, hmod%n

               ! Halo variables
               m = hmod%m(i)
               rv = hmod%rv(i)
               c = hmod%c(i)
               rs = rv/c
               nu = hmod%nu(i)

               ! Correction factor for the gas density profiles
               fc = halo_ejected_gas_fraction(m, hmod, cosm)*m/comoving_matter_density(cosm)
               IF (i_matter /= 0)  wk(i, i_matter) =  wk(i, i_matter)+fc
               IF (i_gas /= 0)     wk(i, i_gas) = wk(i, i_gas)+fc
               IF (i_freegas /= 0) wk(i, i_freegas) = wk(i, i_freegas)+fc
               IF (i_pressure /= 0) THEN

                  ! Calculate the value of the density profile prefactor [(Msun/h)/(Mpc/h)^3] and change units from cosmological to SI
                  rho0 = m*halo_ejected_gas_fraction(m, hmod, cosm) ! rho0 in [(Msun/h)/(Mpc/h)^3]
                  rho0 = rho0*msun/Mpc/Mpc/Mpc                   ! Overflow with REAL(4) if you use Mpc**3, this converts to SI units [h^2 kg/m^3]
                  rho0 = rho0*cosm%h**2                          ! Absorb factors of h, so now [kg/m^3]

                  ! This is the total thermal pressure of the WHIM
                  T0 = HMx_Twhim(hmod, cosm) ! [K]

                  ! Factors to convert from Temp x density -> electron pressure (Temp x n; n is all particle number density)
                  pc = (rho0/(mp*cosm%mup))*(kb*T0) ! Multiply window by *number density* (all particles) times temperature time k_B [J/m^3]
                  pc = pc/(eV*(0.01)**(-3))         ! Change units to pressure in [eV/cm^3]
                  pc = pc*cosm%mup/cosm%mue         ! Convert from total thermal pressure to electron pressure

                  ! Add correction to 'electron pressure' haloes
                  wk(i, i_pressure) = wk(i, i_pressure)+pc

               END IF

            END DO

         END IF

      END IF

   END SUBROUTINE add_smooth_component_to_windows

   FUNCTION wk_product_scatter(n, ifield, k, hmod, cosm)

      ! Calculate W(k)^2 while assuming some scatter in c(M)
      INTEGER, INTENT(IN) :: n
      REAL :: wk_product_scatter(n)
      INTEGER, INTENT(IN) :: ifield(2)
      REAL, INTENT(IN) :: k
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i
      REAL :: m, rv, c

      DO i = 1, n
         m = hmod%m(i)
         rv = hmod%rv(i)
         c = hmod%c(i)
         wk_product_scatter(i) = integrate_scatter(c, c*hmod%sigma_conc, ifield, k, m, rv, hmod, cosm, hmod%acc, 3)
      END DO

   END FUNCTION wk_product_scatter

   REAL FUNCTION scatter_integrand(c, mean_c, sigma_c, ih, k, m, rv, hmod, cosm)

      ! Integrand for computing halo profiles with scatter
      REAL, INTENT(IN) :: c
      REAL, INTENT(IN) :: mean_c
      REAL, INTENT(IN) :: sigma_c
      INTEGER, INTENT(IN) :: ih(2)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: wk(2), pc, rs
      INTEGER :: j
      LOGICAL, PARAMETER :: real_space = .FALSE. ! Fourier profiles

      ! Halo profiles
      DO j = 1, 2
         rs = rv/c
         wk(j) = win(real_space, ih(j), k, m, rv, rs, hmod, cosm)
      END DO

      ! Probability distribution
      pc = lognormal_distribution(c, mean_c, sigma_c)

      ! The full integrand
      scatter_integrand = wk(1)*wk(2)*pc

   END FUNCTION scatter_integrand

   REAL FUNCTION p_2h(ih, wk, n, k, pli, hmod, cosm)

      ! Calculates the 'two-halo' power
      INTEGER, INTENT(IN) :: ih(2)
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: wk(n, 2)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: pli
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sigv, fdamp, rhom, kdamp, ndamp
      REAL :: m0, b0, wk0(2), nu0
      REAL :: m1, m2, rv1, rv2, b1, b2, g1, g2, u1, u2, nu1, nu2
      REAL :: Inl, Inl_11, Inl_21, Inl_12, Inl_22, B_NL
      REAL :: I_11, I_12(n), I_21(n), I_22(n, n)
      REAL :: I2h, I2hs(2)
      REAL :: pq, ph, pnl
      INTEGER :: i, j
      INTEGER, PARAMETER :: iorder = iorder_delta
      INTEGER, PARAMETER :: ifind = ifind_delta
      INTEGER, PARAMETER :: imeth = imeth_delta
      INTEGER, PARAMETER :: HALOFIT_version = HALOFIT_twohalo

      ! Necessary to prevent warning for some reason
      I_11 = 0.

      ! Background matter density
      rhom = comoving_matter_density(cosm)

      IF (hmod%ip2h == 0) THEN

         ! No two-halo term
         p_2h = 0.

      ELSE IF (hmod%ip2h == 1) THEN

         ! Linear theory
         p_2h = pli 

      ELSE IF (hmod%ip2h == 2) THEN

         ! Standard two-halo term

         IF (hmod%imf == 4) THEN

            ! In this case the mass function is a delta function...
            m0 = hmod%hmass
            nu0 = nu_M(m0, hmod, cosm)
            b0 = b_nu(nu0, hmod)
            DO j = 1, 2
               wk0(j) = find(log(m0), hmod%log_m, wk(:, j), n, iorder, ifind, imeth)
               I2hs(j) = b0*rhom*wk0(j)/m0
            END DO

         ELSE

            ! ... otherwise we need to do an integral
            DO j = 1, 2
               CALL I_2h(ih(j), I2h, wk(:, j), n, hmod, cosm, ibias=1)
               I2hs(j) = I2h
            END DO

         END IF

         p_2h = pli*I2hs(1)*I2hs(2)

         IF (hmod%ibias == 2) THEN

            ! Second-order bias correction
            ! This needs to have the property that \int f(nu)b2(nu) du = 0
            ! This means it is hard to check that the normalisation is correct
            ! e.g., how much doees the second-order bias of low mass haloes matter?
            DO j = 1, 2
               CALL I_2h(ih(j), I2h, wk(:, j), n, hmod, cosm, ibias=2)
               I2hs(j) = I2h
            END DO
            p_2h = p_2h+(pli**2)*I2hs(1)*I2hs(2)*rhom**2 ! TODO: Should rhom^2 really be here?

         ELSE IF (hmod%ibias == 3) THEN

            ! Full non-linear bias correction
            DO j = 1, n

               ! Halo 2 parameters
               m2 = hmod%m(j)
               nu2 = hmod%nu(j)
               rv2 = hmod%rv(j)
               b2 = b_nu(nu2, hmod)
               g2 = g_nu(nu2, hmod)
               u2 = rhom*wk(j, 2)/m2

               DO i = 1, n

                  ! Halo 1 parameters
                  m1 = hmod%m(i)
                  nu1 = hmod%nu(i)
                  rv1 = hmod%rv(i)
                  b1 = b_nu(nu1, hmod)
                  g1 = g_nu(nu1, hmod)
                  u1 = rhom*wk(i, 1)/m1

                  ! Get the non-linear bias
                  B_NL = BNL(k, nu1, nu2, rv1, rv2, hmod, cosm)

                  ! Bottom-right quadrant (later multiplied by missing part of integrand)
                  IF (i == 1 .AND. j == 1) THEN
                     I_11 = B_NL*u1*u2
                  END IF

                  ! Integrand for upper-left quadrant. Only single integral over nu2. Only nu2 properties varying.
                  IF (i == 1) THEN
                     I_12(j) = B_NL*b2*g2*u2
                  END IF

                  ! Integrand for bottom-right quadrant. Only single integral over nu1. Only nu1 properties varying.
                  IF (j == 1) THEN
                     I_21(i) = B_NL*b1*g1*u1
                  END IF

                  ! Integrand for upper-right quadrant. Main double integral over nu1 and nu2
                  I_22(i, j) = B_NL*b1*b2*g1*g2*u1*u2

               END DO

            END DO

            ! Calculate boundary terms
            Inl_11 = I_11*hmod%gbmin**2
            Inl_12 = integrate_table(hmod%nu, I_12, 1, n, iorder=1)*hmod%gbmin*(wk(1, 1)*rhom/hmod%m(1))
            Inl_21 = integrate_table(hmod%nu, I_21, 1, n, iorder=1)*hmod%gbmin*(wk(1, 2)*rhom/hmod%m(1))
            Inl_22 = integrate_table(hmod%nu, hmod%nu, I_22)

            ! Sum all terms to get total Inl
            Inl = Inl_22
            IF (add_I_11)          Inl = Inl+Inl_11
            IF (add_I_12_and_I_21) Inl = Inl+Inl_12+Inl_21

            ! This is then the complete two-halo term for non-linear bias
            p_2h = p_2h+pli*Inl

         ELSE IF (hmod%ibias == 4) THEN

            ! Fedeli (2014b) 'scale-dependent halo bias'
            p_2h = p_2h*(1.+k)**0.54

         ELSE IF (hmod%ibias == 5) THEN

            ! My own experimental non-linear halo bias model
            p_2h = p_2h*(1.+(k/hmod%knl))**0.5

         END IF

      ELSE IF (hmod%ip2h == 3) THEN 
         
         ! Damped BAO linear theory
         p_2h = P_dwg(k, hmod%a, hmod%PT_beta*hmod%sigv, cosm) 

      ELSE IF (hmod%ip2h == 4) THEN    
         
         ! One-loop SPT
         p_2h = pli+P_SPT(k, hmod%a, cosm) 

      ELSE IF (hmod%ip2h == 5) THEN

         ! IR resummation ish thing
         p_2h = P_PTish(k, hmod%a, hmod%PT_beta*hmod%sigv, hmod%PT_A, hmod%PT_alpha*hmod%sigv, cosm, approx=.FALSE.)

      ELSE IF (hmod%ip2h == 6 .OR. hmod%ip2h == 7) THEN

         ! HALOFIT either non-linear or quasi-linear term
         CALL calculate_HALOFIT_ka(k, hmod%HALOFIT_neff, hmod%HALOFIT_ncur, hmod%HALOFIT_knl, pli, pnl, pq, ph, hmod%a, cosm, ihf=HALOFIT_version)
         IF (hmod%ip2h == 6) THEN
            p_2h = pnl  
         ELSE IF (hmod%ip2h == 7) THEN
            p_2h = pq 
         ELSE
            STOP 'P_2H: Error, something went wrong with HALOFIT'
         END IF

      ELSE

         STOP 'P_2H: Error, ip2h not specified correclty'

      END IF

      ! Get the normalisation correct for non-matter fields if you are using linear theory or damped BAO
      ! TODO: Does not seem to work with smoothed components (what are smoothed components?)
      IF (hmod%ip2h .NE. 2) THEN
         DO j = 1, 2
            IF (ih(j) == field_dmonly .OR. ih(j) == field_matter) THEN
               I2hs(j) = 1.
            ELSE
               CALL I_2h(ih(j), I2h, wk(:, j), n, hmod, cosm, ibias=1)
               I2hs(j) = I2h
            END IF
         END DO
         p_2h = p_2h*I2hs(1)*I2hs(2)
      END IF

      ! Apply damping to the two-halo term
      IF (hmod%i2hdamp == 1 .OR. hmod%i2hdamp == 2) THEN

         ! HMcode 2015, 2016 damping
         sigv = hmod%sigv
         fdamp = hmod%HMcode_fdamp
         IF (fdamp .NE. 0.) THEN
            p_2h = p_2h*(1.-fdamp*(tanh(k*sigv/sqrt(abs(fdamp))))**2)
         END IF

      ELSE IF (hmod%i2hdamp == 3) THEN

         ! HMcode 2020 damping
         fdamp = hmod%HMcode_fdamp
         kdamp = hmod%HMcode_kdamp
         ndamp = hmod%nd
         IF (kdamp .NE. 0.) THEN
            p_2h = p_2h*(1.-fdamp*((k/kdamp)**ndamp)/(1.+(k/kdamp)**ndamp))
         END IF

      END IF

   END FUNCTION p_2h

   SUBROUTINE I_2h(ih, int, wk, n, hmod, cosm, ibias)

      ! Integral for the two-halo power
      INTEGER, INTENT(IN) :: ih
      REAL, INTENT(OUT) :: int
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: wk(n)
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, INTENT(IN) :: ibias
      REAL :: rhom, b, m, nu, integrand(n)
      INTEGER :: i
      INTEGER, PARAMETER :: iorder = iorder_2halo_integration

      rhom = comoving_matter_density(cosm)

      ! ...otherwise you need to do an integral

      DO i = 1, n

         ! Some variables to make equations cleaner below
         m = hmod%m(i)
         nu = hmod%nu(i)

         ! Linear bias term, standard two-halo term integral
         IF (ibias == 1) THEN
            b = b_nu(nu, hmod)
         ELSE IF (ibias == 2) THEN
            b = b2_nu(nu, hmod)
         ELSE
            STOP 'I_2H: Error, ibias specified incorrectly'
         END IF

         integrand(i) = g_nu(nu, hmod)*b*rhom*wk(i)/m

      END DO

      ! Evaluate these integrals from the tabulated values
      int = integrate_table(hmod%nu, integrand, 1, n, iorder)

      IF (ibias == 1) THEN

         IF (hmod%i2hcor == 0) THEN
            ! Do nothing in this case
            ! There will be large errors if there should be any signal from low-mass haloes
            ! e.g., for the matter power spectrum
         ELSE IF (hmod%i2hcor == 1) THEN
            ! Add on the value of integral b(nu)*g(nu) assuming W(k)=1
            ! Advised by Yoo et al. (2006) and Cacciato et al. (2012)
            ! TODO: This will not work for fields that do not have mass fractions defined
            ! TODO: Could make this correct for fields with no mass fraction, see Appendix B of HMx paper
            int = int+hmod%gbmin*halo_fraction(ih, m, hmod, cosm)
         ELSE IF (hmod%i2hcor == 2) THEN
            ! Put the missing part of the integrand as a delta function at the low-mass limit of the integral
            ! Advised by Schmidt (2016)
            int = int+hmod%gbmin*rhom*wk(1)/hmod%m(1)
         ELSE
            STOP 'P_2h: Error, i2hcor not specified correctly'
         END IF

      END IF

   END SUBROUTINE I_2h

   REAL FUNCTION p_1h(wk_product, n, k, hmod, cosm)

      ! Calculates the one-halo term
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: wk_product(n)
      REAL, INTENT(IN) :: k
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: m, g, fac, ks, wk0_product, m0, rhom, integrand(n)
      INTEGER :: i
      INTEGER, PARAMETER :: iorder_hm = iorder_1halo_integration
      INTEGER, PARAMETER :: iorder_df = iorder_delta
      INTEGER, PARAMETER :: ifind_df = ifind_delta
      INTEGER, PARAMETER :: imeth_df = imeth_delta

      IF (hmod%ip1h == 0) THEN

         ! Set the one-halo term to zero
         p_1h = 0.

      ELSE IF (hmod%ip1h == 1) THEN

         ! Matter density
         rhom = comoving_matter_density(cosm)

         IF (hmod%imf == 4) THEN

            ! In this case the mass function is a delta function...
            m0 = hmod%hmass
            wk0_product = find(log(m0), hmod%log_m, wk_product, n, iorder_delta, ifind_delta, imeth_delta)
            p_1h = rhom*wk0_product/m0

         ELSE

            ! ...otherwise you need to do an integral

            ! Calculates the value of the integrand at all nu values!
            DO i = 1, n
               g = g_nu(hmod%nu(i), hmod)
               m = hmod%m(i)
               integrand(i) = g*wk_product(i)/m
            END DO

            ! Carries out the integration
            p_1h = rhom*integrate_table(hmod%nu, integrand, 1, n, iorder_hm)

         END IF

         ! Convert from P(k) -> Delta^2(k)
         p_1h = p_1h*(4.*pi)*(k/twopi)**3

         IF (hmod%i1hdamp == 0) THEN
            ! Do nothing
         ELSE IF (hmod%i1hdamp == 1) THEN
            ! Damping of the 1-halo term at very large scales
            ks = hmod%HMcode_kstar
            IF (ks == 0. .OR. ((k/ks)**2 > HMcode_ks_limit)) THEN
               ! Prevents problems if k/ks is very large
               fac = 0.
            ELSE
               fac = exp(-((k/ks)**2))
            END IF
            p_1h = p_1h*(1.-fac)
         ELSE IF (hmod%i1hdamp == 2 .OR. hmod%i1hdamp == 3) THEN
            ! Note that the power here should be 4 because it multiplies Delta^2(k) ~ k^3 at low k (NOT 7)
            ! Want f(k<<ks) ~ k^4; f(k>>ks) = 1
            ks = hmod%HMcode_kstar
            IF (ks == 0.) THEN
               fac = 1.
            ELSE
               fac = (k/ks)**4/(1.+(k/ks)**4)
            END IF
            p_1h = p_1h*fac
         ELSE
            STOP 'P_1H: Error, i1hdamp not specified correctly'
         END IF

         ! Renormalise the one-halo term if that is desirous
         IF(hmod%DMONLY_baryon_recipe .AND. renormalise_onehalo) THEN
            p_1h = p_1h*hmod%Rh/hmod%Rhh
         END IF

      ELSE
         STOP 'P_1H: Error, ip1h in halomod not set correctly'
      END IF

   END FUNCTION p_1h

   REAL FUNCTION p_hm(k, pow_2h, pow_1h, hmod)

      ! Combines the two- and one-halo terms into the full halo-model power spectrum
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: pow_2h
      REAL, INTENT(IN) :: pow_1h
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: alpha, pow, con

      ! alpha is set to one sometimes, which is just the standard halo-model sum of terms
      ! TODO: No need to have an IF statement around this?
      IF (hmod%itrans == 1 .OR. hmod%itrans == 2 .OR. hmod%itrans == 3 .OR. hmod%itrans == 5) THEN

         ! If either term is less than zero then we need to be careful
         IF (pow_2h < 0. .OR. pow_1h < 0.) THEN

            ! Either the safe option and just sum the components
            IF (hmod%safe_negative) THEN
               p_hm = pow_2h+pow_1h
            ELSE
               WRITE (*, *) 'P_HM: k [h/Mpc]:', k
               WRITE (*, *) 'P_HM: a:', hmod%a
               WRITE (*, *) 'P_HM: Two-halo term:', pow_2h
               WRITE (*, *) 'P_HM: One-halo term:', pow_1h
               STOP 'P_HM: Error, either two- or one-halo term negative, which is a problem for smoothed transition'
            END IF

         ELSE

            ! Do the standard smoothed transition
            alpha = hmod%HMcode_alpha
            p_hm = (pow_2h**alpha+pow_1h**alpha)**(1./alpha)

         END IF

      ELSE IF (hmod%itrans == 4) THEN

         ! Sigmoid transition
         con = 1. ! Constant for the transition
         pow = 1. ! Power for the transition
         p_hm = pow_2h+sigmoid_log(con*(1.*k/hmod%knl), pow)*(pow_1h-pow_2h)

      ELSE

         ! Standard sum
         p_hm = pow_2h+pow_1h

      END IF

      ! If we are adding in power from voids
      IF (hmod%add_voids) THEN
         p_hm = p_hm+p_1void(k, hmod)
      END IF

   END FUNCTION p_hm

   REAL FUNCTION p_1void(k, hmod)

      ! Calculate the one-void power spectrum
      REAL, INTENT(IN) :: k
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: dc, wk, V, rvoid, rcomp, nu
      REAL :: integrand(hmod%n)
      INTEGER :: i, n

      REAL, PARAMETER :: dv = void_underdensity
      REAL, PARAMETER :: fvoid = void_compensation
      REAL, PARAMETER :: rvoid_simple = simple_void_radius
      LOGICAL, PARAMETER :: compensate = compensate_voids
      LOGICAL, PARAMETER :: simple = simple_voids

      IF (simple) THEN
         n = 1
      ELSE
         n = hmod%n
      END IF

      DO i = 1, n

         !Get the void radius and compensation radius
         IF (simple) THEN
            rvoid = rvoid_simple
         ELSE
            rvoid = hmod%rr(i)
            nu = hmod%nu(i)
         END IF
         rcomp = fvoid*rvoid

         !Calculate the compensation over-density
         dc = -dv*rvoid**3/(rcomp**3-rvoid**3)

         !Calculate the void Fourier transform
         IF (compensate) THEN
            wk = (4.*pi/3.)*((dv-dc)*wk_tophat(k*rvoid)*rvoid**3+dc*wk_tophat(k*rcomp)*rcomp**3)
         ELSE
            wk = (4.*pi/3.)*dv*wk_tophat(k*rvoid)*rvoid**3
         END IF

         !Calculate the void volume
         IF (compensate) THEN
            V = rcomp**3
         ELSE
            V = rvoid**3
         END IF

         IF (simple .EQV. .FALSE.) THEN
            integrand(i) = g_nu(nu, hmod)*wk**2/V
         END IF

      END DO

      !Calculate the void one-halo term
      IF (simple) THEN
         p_1void = wk**2/V
      ELSE
         p_1void = integrate_table(hmod%nu, integrand, 1, n, 1)
      END IF

      p_1void = p_1void*(4.*pi)*(k/twopi)**3

   END FUNCTION p_1void

   REAL FUNCTION BNL(k, nu1, nu2, rv1, rv2, hmod, cosm)

      ! Non-linear halo bias function
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: nu1
      REAL, INTENT(IN) :: nu2
      REAL, INTENT(IN) :: rv1
      REAL, INTENT(IN) :: rv2
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: kk
      REAL, PARAMETER :: kmin = kmin_bnl        ! Below this wavenumber set BNL to zero
      REAL, PARAMETER :: numin = numin_bnl      ! Below this halo mass set BNL to zero
      REAL, PARAMETER :: numax = numax_bnl      ! Above this halo mass set BNL to zero
      REAL, PARAMETER :: min_value = min_bnl    ! Minimum value that BNL is allowed to be (could be below -1)
      LOGICAL, PARAMETER :: halo_exclusion = exclusion_bnl
      LOGICAL, PARAMETER :: fix_min = fix_minimum_bnl

      IF (.NOT. hmod%has_bnl) CALL init_BNL(hmod, cosm)

      ! Ensure that k is not outside array boundary at high end
      kk = k
      CALL fix_maximum(kk, hmod%bnl%xmax)

      ! Decide on which interpolator to evalute or if to set BNL to zero
      IF (kk < kmin) THEN
         BNL = 0.
      ELSE IF (nu1 < numin .OR. nu2 < numin) THEN
         BNL = 0.
      ELSE IF (nu1 > numax .OR. nu2 > numax) THEN
         BNL = 0.
      ELSE
         IF (hmod%stitch_bnl_nu .AND. ((nu1 < hmod%Bnl%ymin) .OR. (nu2 < hmod%Bnl%ymin))) THEN
            BNL = evaluate_interpolator(kk, nu1, nu2, hmod%Bnl_lownu)
         ELSE
            BNL = evaluate_interpolator(kk, nu1, nu2, hmod%Bnl)
         END IF
      END IF

      ! Halo exclusion
      !IF(halo_exclusion) BNL=BNL*exp(-((kk*(rv1+rv2))**2)/4.)
      IF(halo_exclusion) BNL=BNL-sigmoid_tanh(kk**2-1./(rv1+rv2)**2)*(BNL-min_value)

      ! Fix a minimum value for BNL
      IF(fix_min) CALL fix_minimum(BNL, min_value)

   END FUNCTION BNL

   SUBROUTINE init_BNL(hmod, cosm)

      ! Initialisation for the non-linear halo bias term
      ! TODO: Interpolate between redshifts
      USE multidark_stuff
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      TYPE(interpolator3D) :: bnl_interp
      INTEGER :: i, j, ibin, jbin, ik
      INTEGER :: nbin, nk
      INTEGER :: snap
      INTEGER :: u
      INTEGER :: icos
      REAL :: crap, sig8, s, a_MD
      REAL, ALLOCATABLE :: k(:), nu(:), B(:, :, :)
      TYPE(cosmology) :: cosm_MD
      CHARACTER(len=256) :: infile, inbase, base
      CHARACTER(len=256) :: base_bnl, base_bnl_lownu
      LOGICAL, PARAMETER :: verbose = .TRUE.
      INTEGER, PARAMETER :: flag = flag_ucold
      INTEGER, PARAMETER :: icos_MD = icos_bnl_rescale
      REAL, PARAMETER :: R1 = bnl_rescale_R1
      REAL, PARAMETER :: R2 = bnl_rescale_R2
      INTEGER, PARAMETER :: irescale = irescale_sigma

      ! Multidark scale factors
      REAL, PARAMETER :: as_MD(35) = &
         [0.257, 0.287, 0.318, 0.348, 0.378, 0.409, 0.439, &
         0.470, 0.500, 0.530, 0.561, 0.591, 0.621, 0.652, &
         0.682, 0.713, 0.728, 0.743, 0.758, 0.773, 0.788, &
         0.804, 0.819, 0.834, 0.849, 0.864, 0.880, 0.895, &
         0.910, 0.925, 0.940, 0.956, 0.971, 0.986, 1.001]

      ! Snapshots corresponding to scale factors
      INTEGER, PARAMETER :: snaps(35) = &
         [36, 38, 40, 42, 44, 46, 48, &
         50, 52, 54, 56, 58, 60, 62, &
         64, 66, 67, 68, 69, 70, 71, &
         72, 73, 74, 75, 76, 77, 78, &
         79, 80, 81, 82, 83, 84, 85]

      ! Input files
      base_bnl = trim(hmod%bnl_path)//trim(dir_bnl)//'/M512/MDR1_'//trim(hmod%bnl_cat)
      base_bnl_lownu = trim(hmod%bnl_path)//trim(dir_bnl)//'/M512/Bolshoi_'//trim(hmod%bnl_cat)

      IF (verbose) WRITE (*, *) 'INIT_BNL: Running'

      ! Loop over possible different BNL inputs
      DO j = 1, 2

         ! Choose input file type
         IF (j == 1) THEN
            base = base_bnl
         ELSE IF (j == 2 .AND. hmod%stitch_bnl_nu) THEN
            base = base_bnl_lownu
         ELSE
            CYCLE
         END IF

         ! Get the required snapshot
         IF (fixed_bnl) THEN
            IF(snippet_in_string('Bolshoi', base)) THEN
               snap = 416
            ELSE         
               snap = 85
            END IF
         ELSE
            IF(snippet_in_string('Bolshoi', base)) THEN
               STOP 'INIT_BNL: Searching BNL does not work with Bolshoi yet'
            ELSE
               IF (bnl_method == bnl_method_a) THEN
                  snap = nearest_Multidark_snapshot(hmod%a)
               ELSE IF (bnl_method == bnl_method_sig8) THEN
                  sig8 = sigma(8., hmod%a, flag_ucold, cosm)
                  snap = nearest_Multidark_snapshot_sig8(sig8)
               ELSE IF (bnl_method == bnl_method_rescale) THEN
                  icos = icos_MD
                  CALL assign_init_cosmology(icos, cosm_MD, verbose=.FALSE.)
                  CALL calculate_rescaling_parameters(R1, R2, as_MD, s, a_MD, hmod%a, cosm_MD, cosm, irescale, verbose=.FALSE.)
                  snap = nearest_Multidark_snapshot(a_MD)
                  !snap = snaps(array_position(a_MD, as_MD))
               ELSE
                  STOP 'INIT_BNL: Error, bnl_method is specified incorrectly'     
               END IF
            END IF
         END IF
         inbase = trim(base)//'_'//integer_to_string(snap)

         ! Read in the nu values of the bins
         infile = trim(inbase)//'_binstats.dat'
         IF (verbose) WRITE (*, *) 'INIT_BNL: Input file: ', trim(infile)
         nbin = file_length(infile)
         IF (verbose) WRITE (*, *) 'INIT_BNL: Number of nu bins:', nbin
         IF(ALLOCATED(nu)) DEALLOCATE(nu)
         ALLOCATE(nu(nbin))
         OPEN (newunit=u, file=infile)
         DO i = 1, nbin
            READ (u, *) crap, crap, crap, crap, crap, nu(i)
            IF (verbose) WRITE (*, *) 'INIT_BNL: nu bin', i, nu(i)
         END DO
         CLOSE (u)
         IF (verbose) WRITE (*, *) 'INIT_BNL: Done with nu'

         ! Read in k and Bnl(k, nu1, nu2)
         infile = trim(inbase)//'_bnl.dat'
         IF (verbose) WRITE (*, *) 'INIT_BNL: Input file: ', trim(infile)
         nk = file_length(infile)/nbin**2
         IF (verbose) WRITE (*, *) 'INIT_BNL: Number of k values: ', nk          
         IF(ALLOCATED(k)) DEALLOCATE(k)
         IF(ALLOCATED(B)) DEALLOCATE(B)
         ALLOCATE(k(nk), B(nk, nbin, nbin))
         OPEN(newunit=u, file=infile)
         DO ibin = 1, nbin
            DO jbin = 1, nbin
               DO ik = 1, nk
                  READ (u, *) k(ik), B(ik, ibin, jbin)
               END DO
            END DO
         END DO
         CLOSE(u)        
         B = B-1. ! Convert from 1+B_NL to B_NL

         ! Rescale k values if necessary
         IF (bnl_method == bnl_method_rescale) k = s*k

         ! Initialise interpolator
         CALL init_interpolator(k, nu, nu, B, Bnl_interp, &
            iorder = iorder_bnl, &
            iextrap = hmod%iextrap_bnl, &
            store = store_bnl, &
            logx = .TRUE., &
            logy = .FALSE., &
            logz = .FALSE., &
            logf = .FALSE.)

         ! Fill the correct interpolator
         IF (j == 1) THEN
            hmod%Bnl = Bnl_interp
         ELSE IF (j == 2) THEN
            hmod%Bnl_lownu = Bnl_interp
         ELSE IF (j == 3) THEN
            hmod%Bnl_lowk = Bnl_interp
         ELSE IF (j == 4) THEN
            hmod%Bnl_lownu_lowk = Bnl_interp
         ELSE
            STOP 'INIT_BNL: Error, something went wrong with init loop'
         END IF

      END DO

      ! Change flag
      hmod%has_bnl = .TRUE.

      IF (verbose) THEN
         WRITE (*, *) 'INIT_BNL: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_BNL

   REAL FUNCTION T_1h(k1, k2, ih, hmod, cosm)

      ! Halo model one-halo trispectrum
      REAL, INTENT(IN) :: k1, k2
      INTEGER, INTENT(IN) :: ih(2)
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: wk1(hmod%n, 2), wk2(hmod%n, 2)
      REAL :: uk1(hmod%n, 2), uk2(hmod%n, 2), uk_quad(hmod%n), uk0_quad
      REAL :: rhom, g, m, m0, integrand(hmod%n)
      INTEGER :: i
      INTEGER, PARAMETER :: iorder_hm = iorder_1halo_integration
      INTEGER, PARAMETER :: iorder_tr = iorder_delta
      INTEGER, PARAMETER :: ifind_tr = ifind_delta
      INTEGER, PARAMETER :: imeth_tr = imeth_delta

      ! Matter density
      rhom = comoving_matter_density(cosm)

      ! Get the window functions for the two different k values
      CALL init_windows(k1, ih, 2, wk1, hmod%n, hmod, cosm)
      CALL init_windows(k2, ih, 2, wk2, hmod%n, hmod, cosm)

      ! Create normalised wk
      DO i = 1, hmod%n
         m = hmod%m(i)
         uk1(i, :) = wk1(i, :)*rhom/m
         uk2(i, :) = wk2(i, :)*rhom/m
      END DO

      ! TODO: Not sure if this is the correct k and field combinations
      uk_quad = uk1(:, 1)*uk1(:, 2)*uk2(:, 1)*uk2(:, 2)
      !uk_quad=uk1(:,1)*uk1(:,1)*uk2(:,2)*uk2(:,2) ! Could be this
      !uk_quad=uk1(:,1)*uk2(:,2)

      IF (hmod%imf == 4) THEN

         ! In this case the mass function is a delta function...
         m0 = hmod%hmass
         uk0_quad = find(log(m0), hmod%log_m, uk_quad, hmod%n, iorder_tr, ifind_tr, imeth_tr)
         T_1h = rhom*uk0_quad/m0

      ELSE

         ! ...otherwise you need to do an integral

         ! Calculates the value of the integrand at all nu values!
         DO i = 1, hmod%n
            g = g_nu(hmod%nu(i), hmod)
            m = hmod%m(i)
            integrand(i) = g*uk_quad(i)*rhom/m
         END DO

         ! Carries out the integration
         T_1h = integrate_table(hmod%nu, integrand, 1, hmod%n, iorder_hm)

      END IF

   END FUNCTION T_1h

   SUBROUTINE halo_diagnostics(hmod, cosm, dir)

      ! Writes out to file a whole set of halo diagnostics
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=*), INTENT(IN) :: dir
      REAL :: mass
      CHARACTER(len=64) :: ext
      CHARACTER(len=512) :: base
      CHARACTER(len=1024) :: outfile
      INTEGER :: m, mi1, mi2

      ! Integer 10^m to produce haloes between
      REAL, PARAMETER :: m1 = mmin_diag
      REAL, PARAMETER :: m2 = mmax_diag

      WRITE (*, *) 'HALO_DIAGNOSTICS: Outputting diagnostics'

      outfile = TRIM(dir)//'/mass_fractions.dat'
      WRITE (*, *) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
      CALL write_halo_fractions(hmod, cosm, outfile)

      IF (hmod%z == 0.0) THEN
         ext = '_z0.0.dat'
      ELSE IF (hmod%z == 0.5) THEN
         ext = '_z0.5.dat'
      ELSE IF (hmod%z == 1.0) THEN
         ext = '_z1.0.dat'
      ELSE IF (hmod%z == 2.0) THEN
         ext = '_z2.0.dat'
      ELSE
         STOP 'HALO_DIAGNOSTICS: Error, need to make this better with z'
      END IF

      mi1 = NINT(log10(m1))
      mi2 = NINT(log10(m2))

      DO m = mi1, mi2

         mass = 10.**m

         base = TRIM(dir)//'/halo_profile_m'
         outfile = number_file(base, m, ext)
         WRITE (*, *) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
         CALL write_halo_profiles(mass, hmod, cosm, outfile)

         base = TRIM(dir)//'/halo_window_m'
         outfile = number_file(base, m, ext)
         WRITE (*, *) 'HALO_DIAGNOSTICS: ', TRIM(outfile)
         CALL write_halo_transforms(mass, hmod, cosm, outfile)

      END DO

      WRITE (*, *) 'HALO_DIAGNOSTICS: Done'
      WRITE (*, *)

   END SUBROUTINE halo_diagnostics

   SUBROUTINE halo_definitions(hmod, cosm, dir)

      ! Writes out to files the different halo definitions
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=*), INTENT(IN) :: dir
      CHARACTER(len=256) :: fradius, fmass, fconc
      CHARACTER(len=64) :: ext
      INTEGER :: i

      WRITE (*, *) 'HALO_DEFINITIONS: Outputting definitions'

      IF (hmod%z == 0.0) THEN
         ext = '_z0.0.dat'
      ELSE IF (hmod%z == 0.5) THEN
         ext = '_z0.5.dat'
      ELSE IF (hmod%z == 1.0) THEN
         ext = '_z1.0.dat'
      ELSE IF (hmod%z == 2.0) THEN
         ext = '_z2.0.dat'
      ELSE
         STOP 'HALO_DEFINITIONS: Error, need to make this better with z'
      END IF

      fradius = TRIM(dir)//'/radius'//TRIM(ext)
      fmass = TRIM(dir)//'/mass'//TRIM(ext)
      fconc = TRIM(dir)//'/concentration'//TRIM(ext)

      WRITE (*, *) 'HALO_DEFINITIONS: ', TRIM(fradius)
      WRITE (*, *) 'HALO_DEFINITIONS: ', TRIM(fmass)
      WRITE (*, *) 'HALO_DEFINITIONS: ', TRIM(fconc)

      IF (hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod, cosm)

      OPEN (7, file=fradius)
      OPEN (8, file=fmass)
      OPEN (9, file=fconc)
      DO i = 1, hmod%n
         WRITE (7, *) hmod%rv(i), hmod%r200(i), hmod%r500(i), hmod%r200c(i), hmod%r500c(i)
         WRITE (8, *) hmod%m(i), hmod%m200(i), hmod%m500(i), hmod%m200c(i), hmod%m500c(i)
         WRITE (9, *) hmod%c(i), hmod%c200(i), hmod%c500(i), hmod%c200c(i), hmod%c500c(i)
      END DO
      CLOSE (7)
      CLOSE (8)
      CLOSE (9)

      WRITE (*, *) 'HALO_DEFINITIONS: Done'
      WRITE (*, *)

   END SUBROUTINE halo_definitions

   SUBROUTINE halo_properties(hmod, dir)

      ! Writes out the halo properties
      TYPE(halomod), INTENT(INOUT) :: hmod
      CHARACTER(len=*), INTENT(IN) :: dir
      CHARACTER(len=256) :: output
      CHARACTER(len=64) :: ext
      INTEGER :: i

      WRITE (*, *) 'HALO_PROPERTIES: Outputting definitions'

      IF (hmod%z == 0.0) THEN
         ext = '_z0.0.dat'
      ELSE IF (hmod%z == 0.5) THEN
         ext = '_z0.5.dat'
      ELSE IF (hmod%z == 1.0) THEN
         ext = '_z1.0.dat'
      ELSE IF (hmod%z == 2.0) THEN
         ext = '_z2.0.dat'
      ELSE
         STOP 'HALO_PROPERTIES: Error, need to make this better with z'
      END IF

      output = TRIM(dir)//'/properies'//TRIM(ext)
      WRITE (*, *) 'HALO_PROPERTIES: ', TRIM(output)

      OPEN (7, file=output)
      DO i = 1, hmod%n
         WRITE (7, *) hmod%m(i), hmod%rr(i), hmod%rv(i), hmod%nu(i), hmod%c(i), hmod%sig(i), hmod%sigf(i), hmod%zc(i)
      END DO
      CLOSE (7)

      WRITE (*, *) 'HALO_PROPERTIES: Done'
      WRITE (*, *)

   END SUBROUTINE halo_properties

   SUBROUTINE write_halo_fractions(hmod, cosm, outfile)

      ! Writes out the halo mass fractions
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL :: m
      INTEGER :: i
      REAL, PARAMETER :: mmin = mmin_diag
      REAL, PARAMETER :: mmax = mmax_diag
      INTEGER, PARAMETER :: n = n_diag

      OPEN (7, file=outfile)
      DO i = 1, n
         m = exp(progression(log(mmin), log(mmax), i, n))
         WRITE (7, *) m, &
            halo_fraction(field_cdm, m, hmod, cosm), &
            halo_fraction(field_gas, m, hmod, cosm), &
            halo_fraction(field_bound_gas, m, hmod, cosm), &
            halo_fraction(field_ejected_gas, m, hmod, cosm), &
            halo_fraction(field_stars, m, hmod, cosm), &
            halo_fraction(field_central_stars, m, hmod, cosm), &
            halo_fraction(field_satellite_stars, m, hmod, cosm), &
            halo_fraction(field_neutrino, m, hmod, cosm)
      END DO
      CLOSE (7)

   END SUBROUTINE write_halo_fractions

   SUBROUTINE write_halo_profiles(m, hmod, cosm, outfile)

      ! Writes out the halo density profiles
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL :: r, rv, rs, c
      INTEGER :: i, j, nf
      INTEGER, ALLOCATABLE :: fields(:)
      REAL, PARAMETER :: rmin = 1e-3  ! Mininum r/rv
      REAL, PARAMETER :: rmax = 1.1e0 ! Maximum r/rv
      INTEGER, PARAMETER :: n = 512   ! Number of points
      LOGICAL, PARAMETER :: real_space = .TRUE. ! Real profiles
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      ! Calculate halo attributes
      rv = exp(find(log(m), hmod%log_m, log(hmod%rv), hmod%n, iorder, ifind, imeth))
      c = find(log(m), hmod%log_m, hmod%c, hmod%n, iorder, ifind, imeth)
      rs = rv/c

      ! Field types
      nf = 5
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_stars
      fields(5) = field_electron_pressure

      ! Write file
      OPEN (7, file=outfile)
      DO i = 1, n
         r = exp(progression(log(rmin), log(rmax), i, n))
         r = r*rv
         WRITE (7, *) r/rv, (win(real_space, fields(j), r, m, rv, rs, hmod, cosm)*rv**3, j=1, nf) ! rv**3 here is from r^2 dr
      END DO
      CLOSE (7)

   END SUBROUTINE write_halo_profiles

   SUBROUTINE write_halo_transforms(m, hmod, cosm, outfile)

      ! Writes out the Fourier transform of the halo density profiles
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL :: x, rv, c, rs, k, rhobar
      INTEGER :: i, j, nf
      INTEGER, ALLOCATABLE :: fields(:)
      REAL, PARAMETER :: xmin = 1e-1      ! Minimum r/rv
      REAL, PARAMETER :: xmax = 1e2       ! Maximum r/rv
      INTEGER, PARAMETER :: n = 512       ! Number of points
      LOGICAL, PARAMETER :: rsp = .FALSE. ! Fourier profiles
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      ! Calculate halo attributes
      rv = exp(find(log(m), hmod%log_m, log(hmod%rv), hmod%n, iorder, ifind, imeth))
      c = find(log(m), hmod%log_m, hmod%c, hmod%n, iorder, ifind, imeth)
      rs = rv/c

      ! Field types
      nf = 5
      ALLOCATE (fields(nf))
      fields(1) = field_matter
      fields(2) = field_cdm
      fields(3) = field_gas
      fields(4) = field_stars
      fields(5) = field_electron_pressure

      ! Need mean density
      rhobar = comoving_matter_density(cosm)

      ! Write file
      OPEN (7, file=outfile)
      DO i = 1, n
         x = exp(progression(log(xmin), log(xmax), i, n))
         k = x/rv
         WRITE (7, *) x, (win(rsp, fields(j), k, m, rv, rs, hmod, cosm)*rhobar/m, j=1, nf)
      END DO
      CLOSE (7)

   END SUBROUTINE write_halo_transforms

   REAL FUNCTION delta_c(hmod, cosm)

      ! Linear collapse density
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: a, sig

      ! Scale factor
      a = hmod%a

      IF (hmod%idc == 1) THEN
         ! Fixed value
         delta_c = 1.686
      ELSE IF (hmod%idc == 2) THEN
         ! From Nakamura & Suto (1997) LCDM fitting function
         delta_c = dc_NakamuraSuto(a, cosm)
      ELSE IF (hmod%idc == 3 .OR. hmod%idc == 6) THEN
         ! From HMcode (2015, 2016)
         sig = sigma(8., a, flag_matter, cosm)
         delta_c = hmod%dc0+hmod%dc1*log(sig)
         IF (hmod%idc == 3) THEN
            ! HMcode(2016) addition of small cosmology and explicit neutrino dependence          
            delta_c = delta_c*(dc_NakamuraSuto(a, cosm)/dc0)
            IF (cosm%img == img_nDGP) THEN
               delta_c = delta_c*(1.+0.008*sigmoid_tanh(log10(2.33*a*cosm%H0rc**(-0.65))/0.4))
            END IF
         END IF
      ELSE IF (hmod%idc == 4) THEN
         ! From Mead (2017) fitting function
         delta_c = dc_Mead(a, cosm)
      ELSE IF (hmod%idc == 5) THEN
         ! From spherical-collapse calculation
         delta_c = dc_spherical(a, cosm)
      ELSE
         WRITE (*, *) 'DELTA_C: idc:', hmod%idc
         STOP 'DELTA_C: Error, idc defined incorrectly'
      END IF

      ! Massive neutrino dependence
      delta_c = delta_c*(1.+hmod%dcnu*cosm%f_nu)

   END FUNCTION delta_c

   REAL FUNCTION Delta_v(M, hmod, cosm)

      ! Virialised overdensity
      REAL, INTENT(IN) :: M
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: a

      ! Scale factor
      a = hmod%a

      IF (hmod%iDv == 1) THEN
         ! Fixed value; M200
         Delta_v = 200.
      ELSE IF (hmod%iDv == 2) THEN
         ! From Bryan & Norman (1998; arXiv:astro-ph/9710107) fitting functions
         Delta_v = Dv_BryanNorman(a, cosm)
      ELSE IF (hmod%iDv == 3 .OR. hmod%iDv == 8) THEN
         ! From HMcode(2015, 2016)  
         ! This has Omega_m(a) dependence in HMcode (2016), maybe cold would have been better?
         Delta_v = hmod%Dv0*Omega_m(a, cosm)**hmod%Dv1 
         IF (hmod%iDv == 3) THEN
            ! HMcode(2016) additions           
            IF (cosm%img == img_nDGP) THEN
               Delta_v = Delta_v*(1.-a*0.0388*cosm%H0rc**(-0.7))
            ELSE IF (cosm%img == img_fR) THEN
               Delta_v = Delta_v*Geff_fR(M, a, cosm)**(-0.5)
            END IF
         END IF
      ELSE IF (hmod%iDv == 4) THEN
         ! From Mead (2017) fitting function
         Delta_v = Dv_Mead(a, cosm)
      ELSE IF (hmod%iDv == 5) THEN
         ! From spheircal-collapse calculation
         Delta_v = Dv_spherical(a, cosm)
      ELSE IF (hmod%iDv == 6) THEN
         ! Lagrangian radius
         Delta_v = 1.
      ELSE IF (hmod%iDv == 7) THEN
         ! M200c
         Delta_v = 200./Omega_m(a, cosm)
      ELSE IF (hmod%iDv == 9) THEN
         ! 18pi^2 ~178
         Delta_v = Dv0
      ELSE IF (hmod%iDv == 10) THEN
         Delta_v = hmod%Dv0
      ELSE
         STOP 'DELTA_V: Error, iDv defined incorrectly'
      END IF

      ! Massive neutrino dependence
      Delta_v = Delta_v*(1.+hmod%Dvnu*cosm%f_nu)

   END FUNCTION Delta_v

   REAL FUNCTION Geff_fR(M, a, cosm)

      ! This is an approximate toy screening model
      REAL, INTENT(IN) :: M
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: M0
      REAL, PARAMETER :: Gmin = 1.
      REAL, PARAMETER :: Gmax = 4./3.

      ! Screening mass in f(R) models
      M0 = (10**26.4)*(abs(fr_a(a, cosm))**1.5)*(cosm%Om_m**(-0.5))
      Geff_fR = Gmin+(Gmax-Gmin)*sigmoid_tanh(log10(M/M0)/0.38)
  
   END FUNCTION Geff_fR

   REAL FUNCTION HMcode_kstar(hmod, cosm)

      ! Calculates the one-halo damping wave number [h/Mpc]
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sigv, sig8

      IF (hmod%i1hdamp == 0) THEN
         HMcode_kstar = 0.
      ELSE IF (hmod%i1hdamp == 1) THEN
         !sigv = sigmaV(0., hmod%a, flag_matter, cosm)
         sigv = hmod%sigv
         HMcode_kstar = hmod%ks/sigv
      ELSE IF(hmod%i1hdamp == 2) THEN
         HMcode_kstar = hmod%ks
      ELSE IF (hmod%i1hdamp == 3) THEN
         sig8 = sigma(8., hmod%a, flag_ucold, cosm)
         HMcode_kstar = hmod%ks*sig8**hmod%kp
      ELSE
         STOP 'HMCODE_KSTAR: Error, i1hdamp  specified incorrectly'
      END IF

   END FUNCTION HMcode_kstar

   REAL FUNCTION HMcode_fdamp(hmod, cosm)

      ! Calculates the linear-theory damping factor
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sig, sigv

      IF (hmod%i2hdamp == 0) THEN
         ! Set to 0 for the standard linear theory two halo term
         HMcode_fdamp = 0.
      ELSE IF (hmod%i2hdamp == 1) THEN
         ! HMcode (2015)
         sig = sigma(8., hmod%a, flag_matter, cosm)
         HMcode_fdamp = hmod%f0*sig**hmod%f1
      ELSE IF (hmod%i2hdamp == 2) THEN
         ! HMcode (2016)
         sigv = sigmaV(100., hmod%a, flag_matter, cosm)
         HMcode_fdamp = hmod%f0*sigv**hmod%f1
      ELSE IF (hmod%i2hdamp == 3) THEN
         ! HMcode (2020)
         sig = sigma(8., hmod%a, flag_ucold, cosm)
         !sig = sigma(8., hmod%a, flag_matter, cosm)
         HMcode_fdamp = hmod%f0*sig**hmod%f1
      ELSE
         STOP 'HMcode_FDAMP: Error, i2hdamp defined incorrectly'
      END IF

      ! Catches extreme values of fdamp that occur for ridiculous cosmologies
      IF (HMcode_fdamp < HMcode_fdamp_min) HMcode_fdamp = HMcode_fdamp_min
      IF (HMcode_fdamp > HMcode_fdamp_max) HMcode_fdamp = HMcode_fdamp_max

   END FUNCTION HMcode_fdamp

   REAL FUNCTION HMcode_kdamp(hmod, cosm)

      ! 2-halo term damping wavenumber
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sig8

      sig8 = sigma(8., hmod%a, flag_ucold, cosm)
      HMcode_kdamp = hmod%kd*sig8**hmod%kdp

   END FUNCTION HMcode_kdamp

   REAL FUNCTION HMcode_alpha(hmod, cosm)

      ! Calculates the alpha to smooth the two- to one-halo transition
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: n_eff!, x, sig8, y
      REAL :: crap

      ! To prevent compile-time warnings
      crap = cosm%A

      IF (hmod%itrans == 1) THEN
         ! From HMcode (2015, 2016)
         n_eff = hmod%neff
         !n_eff = effective_index(hmod%flag_sigma, hmod, cosm)
         IF(hmod%alp1 == 0.) THEN
            HMcode_alpha = hmod%alp0
         ELSE
            HMcode_alpha = hmod%alp0*(hmod%alp1**n_eff)
         END IF
      ELSE IF (hmod%itrans == 2) THEN
         ! Specially for HMx, exponentiated HMcode (2016) result
         n_eff = hmod%neff
         !n_eff = effective_index(hmod%flag_sigma, hmod, cosm)
         IF(hmod%alp1 == 0.) THEN
            HMcode_alpha = hmod%alp0**1.5
         ELSE
            HMcode_alpha = (hmod%alp0*hmod%alp1**n_eff)**1.5
         END IF
      ELSE IF (hmod%itrans == 3) THEN
         ! Exponentiated HMcode (2016) result
         n_eff = hmod%neff
         !n_eff = effective_index(hmod%flag_sigma, hmod, cosm)
         IF(hmod%alp1 == 0.) THEN
            HMcode_alpha = hmod%alp0**2.5
         ELSE
            HMcode_alpha = (hmod%alp0*hmod%alp1**n_eff)**2.5
         END IF
      ELSE IF (hmod%itrans == 5) THEN
         ! Experimental smoothed transition           
         !x = sigma(8., hmod%a, flag_ucold, cosm)
         !n_eff = effective_index(hmod%flag_sigma, hmod, cosm)
         n_eff = hmod%neff
         !y = effective_curvature(hmod%flag_sigma, hmod, cosm)
         IF (hmod%alp1 == 0.) THEN
            HMcode_alpha = hmod%alp0
         ELSE
            HMcode_alpha = hmod%alp0*(hmod%alp1**n_eff)    ! Relation 1
         END IF
         !HMcode_alpha = hmod%alp0+x*hmod%alp1       ! Relation 2
         !HMcode_alpha = hmod%alp0*abs(x)**hmod%alp1 ! Relation 3
         !HMcode_alpha = hmod%alp0+hmod%alp1*x+hmod%alp2*x**2
         !HMcode_alpha = hmod%alp0+hmod%alp1*x+hmod%alp2*y
         !HMcode_alpha = hmod%alp0*(hmod%alp1**x)+hmod%alp2*y
      ELSE
         HMcode_alpha = 1.
      END IF

      ! Catches values of alpha that are crazy
      IF (HMcode_alpha < HMcode_alpha_min) HMcode_alpha = HMcode_alpha_min
      IF (HMcode_alpha > HMcode_alpha_max) HMcode_alpha = HMcode_alpha_max

   END FUNCTION HMcode_alpha

   REAL FUNCTION HMcode_eta(hmod, cosm)

      ! Calculates the eta that comes into the bastardised one-halo term
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: eta0, sig

      IF (hmod%ieta == 0) THEN
         HMcode_eta = 0.
      ELSE IF (hmod%ieta == 1) THEN
         ! From HMcode(2015; arXiv 1505.07833, 2016)
         IF (hmod%one_parameter_baryons) THEN
            eta0 = 0.98-hmod%As*0.12
         ELSE
            eta0 = hmod%eta0
         END IF
         sig = sigma(8., hmod%a, flag_matter, cosm)
         HMcode_eta = eta0-hmod%eta1*sig
      ELSE IF (hmod%ieta == 2) THEN
         sig = sigma(8., hmod%a, flag_ucold, cosm)
         HMcode_eta = hmod%eta0-hmod%eta1*sig
      ELSE IF (hmod%ieta == 3) THEN
         sig = sigma(8., hmod%a, flag_ucold, cosm)
         HMcode_eta = hmod%eta0*sig**hmod%eta1
      ELSE
         STOP 'HMcode_ETA: Error, ieta defined incorrectly'
      END IF

   END FUNCTION HMcode_eta

   REAL FUNCTION HMcode_A(hmod, cosm)

      ! Halo concentration pre-factor
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sig, As, Az

      IF (hmod%iAs == 0) THEN
         ! Set to 4 for the standard Bullock value
         HMcode_A = 4.
      ELSE IF (hmod%iAs == 1) THEN
         ! This is the 'A' halo-concentration parameter in HMcode(2015; 2016)
         HMcode_A = hmod%As
      ELSE IF (hmod%iAs == 2) THEN
         ! HMcode (2020)
         sig = sigma(8., hmod%a, flag_ucold, cosm)
         IF (hmod%DMONLY_baryon_recipe) THEN
            As = HMcode_DMONLY_baryon_model(cosm%Theat, hmod%As_T, hmod%As)
            Az = HMcode_DMONLY_baryon_model(cosm%Theat, hmod%Az_T, hmod%Az)
         ELSE
            As = hmod%As
            Az = hmod%Az
         END IF
         HMcode_A = (As*(sig**hmod%Ap)+hmod%Ac)*10**(hmod%z*Az)
      ELSE
         STOP 'HMcode_A: Error, iAs defined incorrectly'
      END IF

      ! Now this is divided by 4 so as to be relative to the Bullock base result
      HMcode_A = HMcode_A/4.

   END FUNCTION HMcode_A

   REAL FUNCTION pivot_mass(hmod)

      ! Returns the 'pivot mass' [Msun/h]
      TYPE(halomod), INTENT(IN) :: hmod
      !REAL, PARAMETER :: simple_pivot_mass = 1e14

      IF (hmod%simple_pivot) THEN
         !pivot_mass = simple_pivot_mass
         pivot_mass = hmod%pivot_mass
      ELSE
         pivot_mass = hmod%Mh
      END IF

   END FUNCTION pivot_mass

   REAL FUNCTION HMx_alpha(m, hmod, cosm)

      ! Static gas temperature
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, T, A, B, C, D, E, Mpiv, alpha, alphap, alphaz

      IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN   
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
            alpha = hmod%alpha
            alphap = hmod%alphap
            alphaz = hmod%alphaz
         ELSE IF (hmod%HMx_mode == 6) THEN
            alpha = HMx2020_Temperature_scaling(hmod%alpha_array, hmod, cosm)
            alphap = HMx2020_Temperature_scaling(hmod%alphap_array, hmod, cosm)
            alphaz = HMx2020_Temperature_scaling(hmod%alphaz_array, hmod, cosm)
         ELSE
            STOP 'HMx_ALPHA: Error, HMx_mode not specified correctly'
         END IF
         Mpiv = pivot_mass(hmod)
         z = hmod%z
         IF (hmod%HMx_mode == 3) THEN
            HMx_alpha = alpha*((m/Mpiv)**alphap)*(1.+z)**alphaz
         ELSE IF (is_in_array(hmod%HMx_mode, [5, 6])) THEN
            HMx_alpha = (alpha+z*alphaz)*((m/Mpiv)**alphap)
         ELSE
            STOP 'HMx_ALPHA: Error, HMx_mode not specified correctly'
         END IF
      ELSE IF (hmod%HMx_mode == 4) THEN
         STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
         A = hmod%A_alpha
         B = hmod%B_alpha
         C = hmod%C_alpha
         D = hmod%D_alpha
         E = hmod%E_alpha
         z = hmod%z
         !T = log10(hmod%Theat)
         !HMx_alpha=A*(1.+z)*T+B*(1.+z)+C*T+D
         HMx_alpha = (A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E
         HMx_alpha = HMx_alpha/(1.+z) ! NEW: (1+z) accounts for previous wrong definition of temperature
      ELSE
         HMx_alpha = hmod%alpha
      END IF

      CALL fix_minimum(HMx_alpha, HMx_alpha_min)

   END FUNCTION HMx_alpha

   REAL FUNCTION HMx_beta(m, hmod, cosm)

      ! Hot gas temperature
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, Mpiv, beta, betaz, betap
    
      IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
         beta = hmod%beta
         betap = hmod%betap
         betaz = hmod%betaz
         Mpiv = pivot_mass(hmod)
         z = hmod%z
         HMx_beta = beta*((m/Mpiv)**betap)*(1.+z)**betaz
      ELSE IF (hmod%HMx_mode == 4) THEN
         HMx_beta = HMx_alpha(m, hmod, cosm)
      ELSE
         HMx_beta = hmod%beta
      END IF

      CALL fix_minimum(HMx_beta, HMx_beta_min)

   END FUNCTION HMx_beta

   REAL FUNCTION HMx_eps(hmod, cosm)

      ! Halo concentration epsilon
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, T, A, B, C, D, eps, epsz

      IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
            eps = hmod%eps
            epsz = hmod%epsz
         ELSE IF (hmod%HMx_mode == 6) THEN
            eps = HMx2020_Temperature_scaling(hmod%eps_array, hmod, cosm)
            epsz = HMx2020_Temperature_scaling(hmod%epsz_array, hmod, cosm)
         ELSE
            STOP 'HMx_EPS: Error, HMx_mode not specified correctly'
         END IF
         z = hmod%z
         IF(hmod%HMx_mode == 3) THEN
            HMx_eps = eps*(1.+z)**epsz
         ELSE IF(hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            HMx_eps = eps+z*epsz
         ELSE
            STOP 'HMx_EPS: Error, HMx_mode not specified correctly'
         END IF 
      ELSE IF (hmod%HMx_mode == 4) THEN
         STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
         A = hmod%A_eps
         B = hmod%B_eps
         C = hmod%C_eps
         D = hmod%D_eps
         z = hmod%z
         !T = log10(hmod%Theat)
         HMx_eps = A*(1.+z)*T+B*(1.+z)+C*T+D
         HMx_eps = 10**HMx_eps-1.
      ELSE
         HMx_eps = hmod%eps
      END IF

      CALL fix_minimum(HMx_eps, HMx_eps_min)

   END FUNCTION HMx_eps

   REAL FUNCTION HMx_eps2(hmod, cosm)

      ! Halo concentration epsilon
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, eps2, eps2z

      IF (hmod%HMx_mode == 1) THEN
         HMx_eps2 = HMx_eps(hmod, cosm)
      ELSE IF (is_in_array(hmod%HMx_mode, [3, 4, 5, 6])) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 4 .OR. hmod%HMx_mode == 5) THEN
            eps2 = hmod%eps2
            eps2z = hmod%eps2z
         ELSE IF (hmod%HMx_mode == 6) THEN
            eps2 = HMx2020_Temperature_scaling(hmod%eps2_array, hmod, cosm)
            eps2z = HMx2020_Temperature_scaling(hmod%eps2z_array, hmod, cosm)
         ELSE
            STOP 'HMX_EPS2: Error, something went wrong'
         END IF        
         z = hmod%z
         HMx_eps2 = eps2+z*eps2z
      ELSE
         HMx_eps2 = hmod%eps2
      END IF

   END FUNCTION HMx_eps2

   REAL FUNCTION HMx_Gamma(m, hmod, cosm)

      ! Komatsu-Seljak gas profile index
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, T, A, B, C, D, E, Mpiv, Gamma, Gammap, Gammaz

      IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
            Gamma = hmod%Gamma
            Gammap = hmod%Gammap
            Gammaz = hmod%Gammaz
         ELSE IF (hmod%HMx_mode == 6) THEN
            Gamma = HMx2020_Temperature_scaling(hmod%Gamma_array, hmod, cosm)
            Gammap = HMx2020_Temperature_scaling(hmod%Gammap_array, hmod, cosm)
            Gammaz = HMx2020_Temperature_scaling(hmod%Gammaz_array, hmod, cosm)
         ELSE
            STOP 'HMx_GAMMA: Error, HMx_mode not specified correctly'
         END IF 
         Mpiv = pivot_mass(hmod)
         z = hmod%z
         IF (hmod%HMx_mode == 3) THEN
            HMx_Gamma = Gamma*((m/Mpiv)**Gammap)*(1.+z)**Gammaz
         ELSE IF (hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            HMx_Gamma = (Gamma+z*Gammaz)*((m/Mpiv)**Gammap)
         ELSE
            STOP 'HMx_GAMMA: Error, HMx_mode not specified correctly'
         END IF
      ELSE IF (hmod%HMx_mode == 4) THEN
         STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
         A = hmod%A_Gamma
         B = hmod%B_Gamma
         C = hmod%C_Gamma
         D = hmod%D_Gamma
         E = hmod%E_Gamma
         z = hmod%z
         !T = log10(hmod%Theat)
         !HMx_Gamma=A*(1.+z)*T+B*(1.+z)+C*T+D
         HMx_Gamma = (A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E
      ELSE
         HMx_Gamma = hmod%Gamma
      END IF

      CALL fix_minimum(HMx_Gamma, HMx_Gamma_min)
      CALL fix_maximum(HMx_Gamma, HMx_Gamma_max)

   END FUNCTION HMx_Gamma

   REAL FUNCTION HMx_Zamma(m, hmod, cosm)

      ! Komatsu-Seljak index for pressure profile if different from gas density
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, Mpiv, Zamma, Zammap, Zammaz

      IF (.NOT. hmod%different_Gammas) THEN
         HMx_Zamma = HMx_Gamma(m, hmod, cosm)
      ELSE
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
               Zamma = hmod%Zamma
               Zammap = hmod%Zammap
               Zammaz = hmod%Zammaz
            ELSE IF (hmod%HMx_mode == 6) THEN
               Zamma = HMx2020_Temperature_scaling(hmod%Zamma_array, hmod, cosm)
               Zammap = HMx2020_Temperature_scaling(hmod%Zammap_array, hmod, cosm)
               Zammaz = HMx2020_Temperature_scaling(hmod%Zammaz_array, hmod, cosm)
            ELSE
               STOP 'HMx_GAMMA_PRESSURE: Error, HMx_mode not specified correctly'
            END IF 
            Mpiv = pivot_mass(hmod)
            z = hmod%z
            IF (hmod%HMx_mode == 3) THEN
               HMx_Zamma = Zamma*((m/Mpiv)**Zammap)*(1.+z)**Zammaz
            ELSE IF (hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
               HMx_Zamma = (Zamma+z*Zammaz)*((m/Mpiv)**Zammap)
            ELSE
               STOP 'HMx_GAMMA_PRESSURE: Error, HMx_mode not specified correctly'
            END IF
         ELSE
            HMx_Zamma = hmod%Zamma
         END IF

         CALL fix_minimum(HMx_Zamma, HMx_Gamma_min)
         CALL fix_maximum(HMx_Zamma, HMx_Gamma_max)

      END IF

   END FUNCTION HMx_Zamma

   REAL FUNCTION HMx_M0(hmod, cosm)

      ! Gas fraction turn-over mass
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, T, A, B, C, D, E, M0, M0z

      IF (hmod%HMx_Mode == 1) THEN
         HMx_M0 = hmod%M0*10**(hmod%M0z*hmod%z)
      ELSE IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
            M0 = hmod%M0
            M0z = hmod%M0z
         ELSE IF (hmod%HMx_mode == 6) THEN 
            M0 = exp(HMx2020_Temperature_scaling(log(hmod%M0_array), hmod, cosm))
            M0z = HMx2020_Temperature_scaling(hmod%M0z_array, hmod, cosm)
         ELSE
            STOP 'HMx_M0: Error, HMx_mode not specified correctly'
         END IF
         z = hmod%z
         IF(hmod%HMx_mode == 3) THEN
            HMx_M0 = M0**((1.+z)**M0z)
         ELSE IF (hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            HMx_M0 = M0*exp(M0z*z)
         ELSE
            STOP 'HMx_M0: Error, HMx_mode not specified correctly'
         END IF
      ELSE IF (hmod%HMx_mode == 4) THEN
         STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
         A = hmod%A_M0
         B = hmod%B_M0
         C = hmod%C_M0
         D = hmod%D_M0
         E = hmod%E_M0
         z = hmod%z
         HMx_M0 = (A*(1.+z)**2+B*(1.+z)+C)*T**2+D*T+E
         HMx_M0 = 10**HMx_M0
      ELSE
         HMx_M0 = hmod%M0
      END IF

   END FUNCTION HMx_M0

   REAL FUNCTION HMx_Astar(hmod, cosm)

      ! Star fraction ampltiude
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, T, A, B, C, D, Astar, Astarz

      IF (hmod%HMx_mode == 1) THEN
         HMx_Astar = hmod%Astar*(1.+hmod%Astarz*hmod%z**2)
      ELSE IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
            Astar = hmod%Astar
            Astarz = hmod%Astarz     
         ELSE IF (hmod%HMx_mode == 6) THEN
            Astar = HMx2020_Temperature_scaling(hmod%Astar_array, hmod, cosm)
            Astarz = HMx2020_Temperature_scaling(hmod%Astarz_array, hmod, cosm)      
         ELSE
            STOP 'HMx_ASTAR: Error, HMx_mode not specified correctly'
         END IF
         z = hmod%z
         IF (hmod%HMx_mode == 3) THEN
            HMx_Astar = Astar*(1.+z)**Astarz
         ELSE IF (hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            HMx_Astar = Astar+z*Astarz
         ELSE
            STOP 'HMx_ASTAR: Error, HMx_mode not specified correctly'
         END IF
      ELSE IF (hmod%HMx_mode == 4) THEN
         STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
         A = hmod%A_Astar
         B = hmod%B_Astar
         C = hmod%C_Astar
         D = hmod%D_Astar
         z = hmod%z
         !T = log10(hmod%Theat)
         HMx_Astar = A*(1.+z)*T+B*(1.+z)+C*T+D
      ELSE
         HMx_Astar = hmod%Astar
      END IF

      CALL fix_minimum(HMx_Astar, HMx_Astar_min)

   END FUNCTION HMx_Astar

   REAL FUNCTION HMx_Twhim(hmod, cosm)

      ! WHIM temperature
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, T, A, B, C, D, Twhim, Twhimz

      IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 5) THEN
            Twhim = hmod%Twhim
            Twhimz = hmod%Twhimz
         ELSE IF (hmod%HMx_mode == 6) THEN
            Twhim = exp(HMx2020_Temperature_scaling(log(hmod%Twhim_array), hmod, cosm))
            Twhimz = HMx2020_Temperature_scaling(hmod%Twhimz_array, hmod, cosm)
         ELSE
            STOP 'HMx_TWHIM: Error, HMx_mode not specified correctly'
         END IF
         z = hmod%z
         IF (hmod%HMx_mode == 3) THEN
            HMx_Twhim = Twhim**((1.+z)**Twhimz)
         ELSE IF (is_in_array(hmod%HMx_mode, [5, 6])) THEN
            !HMx_Twhim = Twhim*(exp(Twhimz)**z)
            HMx_Twhim = Twhim*exp(Twhimz*z)
         ELSE
            STOP 'HMx_TWHIM: Error, HMx_mode not specified correctly'
         END IF
      ELSE IF (hmod%HMx_mode == 4) THEN
         STOP 'ASSIGN_HALOMOD: Halomodel Theat no longer supported'
         A = hmod%A_Twhim
         B = hmod%B_Twhim
         C = hmod%C_Twhim
         D = hmod%D_Twhim
         z = hmod%z
         !T = log10(hmod%Theat)
         !HMx_Twhim=A*(1.+z)*T+B*(1.+z)+C*T+D
         HMx_Twhim = (A*(1.+z)**2+B*(1.+z)+C)*T+D
         HMx_Twhim = 10**HMx_Twhim
      ELSE
         HMx_Twhim = hmod%Twhim
      END IF

   END FUNCTION HMx_Twhim

   REAL FUNCTION HMx_cstar(m, hmod, cosm)

      ! Stellar-density concentration
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: Mpiv, z, cstar, cstarp, cstarz

      IF (is_in_array(hmod%HMx_mode, [3, 4, 5, 6])) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 4 .OR. hmod%HMx_mode == 5) THEN
            cstar = hmod%cstar
            cstarp = hmod%cstarp
            cstarz = hmod%cstarz
         ELSE IF (hmod%HMx_mode == 6) THEN
            cstar = HMx2020_Temperature_scaling(hmod%cstar_array, hmod, cosm)
            cstarp = HMx2020_Temperature_scaling(hmod%cstarp_array, hmod, cosm)
            cstarz = HMx2020_Temperature_scaling(hmod%cstarz_array, hmod, cosm)
         ELSE
            STOP 'HMx_CSTAR: Error, HMx_mode not specified correctly'
         END IF
         Mpiv = pivot_mass(hmod)
         z = hmod%z
         HMx_cstar = cstar*((m/Mpiv)**cstarp)*(1.+z)**cstarz
      ELSE
         HMx_cstar = hmod%cstar
      END IF

   END FUNCTION HMx_cstar

   REAL FUNCTION HMx_sstar(hmod, cosm)

      ! Star fraction width
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: crap

      ! Prevent warning
      crap = cosm%A

      HMx_sstar = hmod%sstar

   END FUNCTION HMx_sstar

   REAL FUNCTION HMx_Mstar(hmod, cosm)

      ! Peak halo mass for star-formation efficiency
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: Mstar, Mstarz, z

      IF (is_in_array(hmod%HMx_mode, [3, 4, 5, 6])) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 4 .OR. hmod%HMx_mode == 5) THEN
            Mstar = hmod%mstar
            Mstarz = hmod%mstarz
         ELSE IF (hmod%HMx_mode == 6) THEN
            Mstar = exp(HMx2020_Temperature_scaling(log(hmod%Mstar_array), hmod, cosm))
            Mstarz = HMx2020_Temperature_scaling(hmod%Mstarz_array, hmod, cosm)
         ELSE
            STOP 'HMx_MSTAR: Error, HMx_mode not specified correctly'
         END IF
         z = hmod%z
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 4) THEN
            HMx_Mstar = Mstar**((1.+z)**Mstarz)
         ELSE IF(hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            HMx_Mstar = Mstar*exp(Mstarz*z)
         ELSE  
            STOP 'HMx_MSTAR: Error, HMx_mode not specified correctly'
         END IF
      ELSE
         HMx_Mstar = hmod%Mstar
      END IF

   END FUNCTION HMx_Mstar

   REAL FUNCTION HMx_fcold(hmod, cosm)

      ! Cold-gas fraction
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: crap

      ! Prevent warning
      crap = cosm%A

      HMx_fcold = hmod%fcold

   END FUNCTION HMx_fcold

   REAL FUNCTION HMx_fhot(hmod, cosm)

      ! Hot-gas fraction
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: crap

      ! Prevent warning
      crap = cosm%A

      HMx_fhot = hmod%fhot

   END FUNCTION HMx_fhot

   REAL FUNCTION HMx_eta(hmod, cosm)

      ! Satellite galaxy fraction
      TYPE(cosmology), INTENT(IN) :: cosm
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: eta, etaz, z

      IF (is_in_array(hmod%HMx_mode, [3, 4, 5, 6])) THEN
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 4 .OR. hmod%HMx_mode == 5) THEN
            eta = hmod%eta
            etaz = hmod%etaz
         ELSE IF (hmod%HMx_mode == 6) THEN
            eta = HMx2020_Temperature_scaling(hmod%eta_array, hmod, cosm)
            etaz = HMx2020_Temperature_scaling(hmod%etaz_array, hmod, cosm)        
         ELSE
            STOP 'HMx_ETA: Error, HMx_mode not specified correctly'
         END IF
         z = hmod%z
         IF (hmod%HMx_mode == 3 .OR. hmod%HMx_mode == 4) THEN
            HMx_eta = eta*(1.+z)**etaz
         ELSE IF (hmod%HMx_mode == 5 .OR. hmod%HMx_mode == 6) THEN
            HMx_eta = eta+z*etaz
         ELSE
            STOP 'HMx_ETA: Error, HMx_mode not specified correctly'
         END IF
      ELSE
         HMx_eta = hmod%eta
      END IF

      IF (HMx_eta > HMx_eta_min) HMx_eta = HMx_eta_min

   END FUNCTION HMx_eta

   REAL FUNCTION HMx_ibeta(m, hmod, cosm)

      ! Isothermal-beta profile index
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: z, Mpiv, ibeta, ibetap, ibetaz, crap

      ! Prevent warning
      crap = cosm%A

      ibeta = hmod%ibeta
      ibetap = hmod%ibetap
      ibetaz = hmod%ibetaz

      Mpiv = pivot_mass(hmod)
      z = hmod%z      
      HMx_ibeta = ibeta*((m/Mpiv)**ibetap)*((1.+z)**ibetaz)

   END FUNCTION HMx_ibeta

   REAL FUNCTION HMx_gbeta(hmod, cosm)

      ! Power-law decline for halo gas fractions
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: gbeta, gbetaz, z, crap

      ! Prevent warning
      crap = cosm%A
 
      gbeta = hmod%gbeta
      gbetaz = hmod%gbetaz

      z = hmod%z
      HMx_gbeta = gbeta*(1.+z)**gbetaz

   END FUNCTION HMx_gbeta

   REAL FUNCTION HMx2020_Temperature_scaling(array, hmod, cosm)

      REAL, INTENT(IN) :: array(3) ! Array containing three values of the parameter to be found   
      TYPE(halomod), INTENT(IN) :: hmod ! Halomodel
      TYPE(cosmology), INTENT(IN) :: cosm
      INTEGER, PARAMETER :: iorder = 1 ! Linear interpolation here (only 3 points)
      INTEGER, PARAMETER :: ifind = ifind_split
      INTEGER, PARAMETER :: iinterp = iinterp_Lagrange
      INTEGER, PARAMETER :: n = 3
      REAL :: logT

      logT = log(cosm%Theat)
      HMx2020_Temperature_scaling = find(logT, log(hmod%Theat_array), array, n, iorder, ifind, iinterp)

   END FUNCTION HMx2020_Temperature_scaling

   REAL FUNCTION r_nl(hmod)

      ! Calculates Lagrangian R_nl where nu(R_nl)=1.
      TYPE(halomod), INTENT(INOUT) :: hmod
      INTEGER, PARAMETER :: iorder = iorder_inversion_rnl
      INTEGER, PARAMETER :: ifind = ifind_inversion_rnl
      INTEGER, PARAMETER :: iinterp = iinterp_inversion_rnl

      IF ((.NOT. rnl_extrap) .AND. hmod%nu(1) > 1.) THEN
         ! This catches some very strange values for cosmologies where the non-linear radius can be very large
         ! This can be a very bad idea since r_nl is wrong and depends on the lower-mass limit in hmod
         r_nl = hmod%rr(1)
      ELSE
         r_nl = exp(find(log(1.), log(hmod%nu), log(hmod%rr), hmod%n, iorder, ifind, iinterp))
      END IF

   END FUNCTION r_nl

   SUBROUTINE allocate_halomod(hmod)

      ! Allocates memory for the look-up tables at sets array to zero
      TYPE(halomod), INTENT(INOUT) :: hmod
      INTEGER :: n

      ! Number of entries in look-up table
      n = hmod%n

      ! Allocate all arrays
      ALLOCATE (hmod%log_m(n))
      ALLOCATE (hmod%zc(n), hmod%m(n), hmod%c(n), hmod%rv(n))!, hmod%mr(n))
      ALLOCATE (hmod%nu(n), hmod%rr(n), hmod%sigf(n), hmod%sig(n))
      ALLOCATE (hmod%m500(n), hmod%r500(n), hmod%c500(n))
      ALLOCATE (hmod%m500c(n), hmod%r500c(n), hmod%c500c(n))
      ALLOCATE (hmod%m200(n), hmod%r200(n), hmod%c200(n))
      ALLOCATE (hmod%m200c(n), hmod%r200c(n), hmod%c200c(n))

      hmod%log_m = 0.
      hmod%zc = 0.
      hmod%m = 0.
      hmod%c = 0.
      hmod%rv = 0.
      hmod%nu = 0.
      hmod%rr = 0.
      hmod%sigf = 0.
      hmod%sig = 0.

      hmod%m500 = 0.
      hmod%r500 = 0.
      hmod%c500 = 0.

      hmod%m500c = 0.
      hmod%r500c = 0.
      hmod%c500c = 0.

      hmod%m200 = 0.
      hmod%r200 = 0.
      hmod%c200 = 0.

      hmod%m200c = 0.
      hmod%r200c = 0.
      hmod%c200c = 0.

   END SUBROUTINE allocate_halomod

   SUBROUTINE deallocate_halomod(hmod)

      ! Deallocates the look-up tables
      TYPE(halomod) :: hmod

      ! Deallocates look-up tables
      DEALLOCATE (hmod%log_m)
      DEALLOCATE (hmod%zc, hmod%m, hmod%c, hmod%rv)!, hmod%mr)
      DEALLOCATE (hmod%nu, hmod%rr, hmod%sigf, hmod%sig)
      DEALLOCATE (hmod%m500, hmod%r500, hmod%c500, hmod%m500c, hmod%r500c, hmod%c500c)
      DEALLOCATE (hmod%m200, hmod%r200, hmod%c200, hmod%m200c, hmod%r200c, hmod%c200c)

   END SUBROUTINE deallocate_halomod

   REAL FUNCTION Omega_stars(hmod, cosm)

      ! Calculate the cosmological density in star mass
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (hmod%imf == 4) THEN
         Omega_stars = halo_star_fraction(hmod%hmass, hmod, cosm)
         Omega_stars = Omega_stars*comoving_matter_density(cosm)/comoving_critical_density(hmod%a, cosm)
      ELSE
         Omega_stars = rhobar_tracer(hmod%nu(1), hmod%large_nu, rhobar_star_integrand, hmod, cosm)
         Omega_stars = Omega_stars/comoving_critical_density(hmod%a, cosm)
      END IF

   END FUNCTION Omega_stars

   SUBROUTINE init_galaxies(hmod, cosm)

      ! Calculate the number densities of galaxies
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: nu_min, nu_max

      nu_min = nu_M(hmod%mhalo_min, hmod, cosm)
      nu_max = nu_M(hmod%mhalo_max, hmod, cosm)
      hmod%n_c = rhobar_tracer(nu_min, nu_max, rhobar_central_integrand, hmod, cosm)
      hmod%n_s = rhobar_tracer(nu_min, nu_max, rhobar_satellite_integrand, hmod, cosm)
      hmod%n_g = hmod%n_c+hmod%n_s
      IF (verbose_galaxies) THEN
         WRITE (*, *) 'INIT_GALAXIES: Comoving number density of central galaxies [(Mpc/h)^-3]:', REAL(hmod%n_c)
         WRITE (*, *) 'INIT_GALAXIES: Comoving number density of satellite galaxies [(Mpc/h)^-3]:', REAL(hmod%n_s)
         WRITE (*, *) 'INIT_GALAXIES: Comoving number density of all galaxies [(Mpc/h)^-3]:', REAL(hmod%n_g)
         WRITE (*, *)
      END IF

      hmod%has_galaxies = .TRUE.

   END SUBROUTINE init_galaxies

   SUBROUTINE init_HI(hmod, cosm)

      ! Calculates the background HI density by integrating the HI mass function
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: nu_min, nu_max

      nu_min = hmod%nu(1)
      nu_max = hmod%nu(hmod%n)
      hmod%rho_HI = rhobar_tracer(nu_min, hmod%large_nu, rhobar_HI_integrand, hmod, cosm)
      IF (verbose_HI) THEN
         WRITE (*, *) 'INIT_HI: z:', hmod%z
         WRITE (*, *) 'INIT_HI: HI density [log10(rho/(Msun/h)/(Mpc/h)^3)]:', real(log10(hmod%rho_HI))
         WRITE (*, *) 'INIT_HI: Omega_HI(z):', hmod%rho_HI/comoving_critical_density(hmod%a, cosm)
         WRITE (*, *) 'INIT_HI: Omega_HI(z) relative to z=0 critical density:', hmod%rho_HI/comoving_critical_density(1., cosm)
         WRITE (*, *) 'INIT_HI: rho_HI / rho_matter:', hmod%rho_HI/comoving_matter_density(cosm)
         WRITE (*, *)
      END IF

      hmod%has_HI = .TRUE.

   END SUBROUTINE init_HI

   REAL FUNCTION nu_R(R, hmod, cosm)

      ! Calculates nu(R) where R is the comoving Lagrangian halo radius
      REAL, INTENT(IN) :: R
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      nu_R = hmod%dc/sigma(R, hmod%a, hmod%flag_sigma, cosm)

   END FUNCTION nu_R

   REAL FUNCTION nu_M(M, hmod, cosm)

      ! Calculates nu(M) where M is the halo mass
      REAL, INTENT(IN) :: M
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: R

      R = Lagrangian_radius(M, cosm)
      nu_M = nu_R(R, hmod, cosm)

   END FUNCTION nu_M

   REAL FUNCTION M_nu(nu, hmod)

      ! Calculates M(nu) where M is the halo mass and nu is the peak height
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      IF(hmod%saturation .AND. nu<hmod%nu_saturation) THEN
         M_nu = 0.
      ELSE
         M_nu = exp(find(nu, hmod%nu, hmod%log_m, hmod%n, iorder, ifind, imeth))
      END IF

   END FUNCTION M_nu

   REAL FUNCTION rhobar_central_integrand(nu, hmod, cosm)

      ! Integrand for the number density of central galaxies
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Could remove
      REAL :: M, crap

      crap = cosm%A

      M = M_nu(nu, hmod)
      rhobar_central_integrand = N_centrals(M, hmod)*g_nu(nu, hmod)/M

   END FUNCTION rhobar_central_integrand

   REAL FUNCTION rhobar_satellite_integrand(nu, hmod, cosm)

      ! Integrand for the number density of satellite galaxies
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Could remove
      REAL :: M, crap

      crap = cosm%A

      M = M_nu(nu, hmod)
      rhobar_satellite_integrand = N_satellites(M, hmod)*g_nu(nu, hmod)/M

   END FUNCTION rhobar_satellite_integrand

   REAL FUNCTION rhobar_star_integrand(nu, hmod, cosm)

      ! Integrand for the matter density of stars
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: M

      M = M_nu(nu, hmod)
      rhobar_star_integrand = halo_star_fraction(M, hmod, cosm)*g_nu(nu, hmod)

   END FUNCTION rhobar_star_integrand

   REAL FUNCTION rhobar_HI_integrand(nu, hmod, cosm)

      ! Integrand for the HI mass density
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: M

      M = M_nu(nu, hmod)
      rhobar_HI_integrand = halo_HI_fraction(M, hmod, cosm)*g_nu(nu, hmod)

   END FUNCTION rhobar_HI_integrand

   REAL FUNCTION rhobar_tracer(nu_min, nu_max, integrand, hmod, cosm)

      ! Calculate the mean density of a tracer
      ! Integrand here is a function of mass, i.e. I(M); R = rho * Int I(M)dM
      ! TODO: This function is comparatively slow and should be accelerated somehow
      ! TODO: This uses integrate_hmod_cosm_exp, which is weird, surely can use some transformed integrand instead?
      REAL, INTENT(IN) :: nu_min, nu_max
      REAL, EXTERNAL :: integrand
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      INTERFACE
         FUNCTION integrand(M, hmod_interface, cosm_interface)
            IMPORT :: halomod, cosmology   
            REAL, INTENT(IN) :: M
            TYPE(halomod), INTENT(INOUT) :: hmod_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION integrand
      END INTERFACE

      rhobar_tracer=integrate_hmod_cosm_exp(log(nu_min), log(nu_max), integrand, hmod, cosm, hmod%acc, 3)
      rhobar_tracer=rhobar_tracer*comoving_matter_density(cosm)

   END FUNCTION rhobar_tracer

   REAL FUNCTION one_halo_amplitude(hmod, cosm)

      ! Calculates the amplitude of the shot-noise plateau of the one-halo term [Mpc/h]^3
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: integrand(hmod%n), g, m
      INTEGER :: i

      IF (hmod%imf == 4) THEN
         ! Special case for delta-function mass function
         one_halo_amplitude = hmod%hmass
      ELSE
         ! Calculates the value of the integrand at all nu values!
         DO i = 1, hmod%n
            g = g_nu(hmod%nu(i), hmod)
            m = hmod%m(i)
            integrand(i) = g*m
         END DO
         one_halo_amplitude = integrate_table(hmod%nu, integrand, 1, hmod%n, 1)
      END IF

      one_halo_amplitude = one_halo_amplitude/comoving_matter_density(cosm)

   END FUNCTION one_halo_amplitude

   REAL FUNCTION one_halo_modified_amplitude(hmod, cosm)

      ! Calculates the amplitude of the shot-noise plateau of the one-halo term [Mpc/h]^3
      ! TODO: Is this necesarry? How is this different from one_halo_amplitude?
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: integrand(hmod%n), g, m
      INTEGER :: i

      IF (hmod%imf == 4) THEN
         ! Special case for delta-function mass function
         one_halo_modified_amplitude = hmod%hmass*DMONLY_halo_mass_fraction(hmod%hmass, hmod, cosm)**2
      ELSE
         ! Calculates the value of the integrand at all nu values!
         DO i = 1, hmod%n
            g = g_nu(hmod%nu(i), hmod)
            m = hmod%m(i)
            integrand(i) = g*m*DMONLY_halo_mass_fraction(m, hmod, cosm)**2
         END DO
         one_halo_modified_amplitude = integrate_table(hmod%nu, integrand, 1, hmod%n, 1)
      END IF

      one_halo_modified_amplitude = one_halo_modified_amplitude/comoving_matter_density(cosm)

   END FUNCTION one_halo_modified_amplitude

   REAL FUNCTION mass_function(m, hmod, cosm)

      ! Calculates the halo mass function, what I call n(M); some people call dn/dM 
      ! In either case the units are [(Msun/h)^-1(Mpc/h)^-3]
      ! TODO: Relies on a potentially unstable numerical derivative
      REAL, INTENT(IN) :: m                  ! Halo mass [Msun/h]
      TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

      IF(m == 0.) THEN
         mass_function = 0.
      ELSE
         mass_function = multiplicity_function(m, hmod, cosm)*comoving_matter_density(cosm)/m**2
      END IF

   END FUNCTION mass_function

   REAL FUNCTION multiplicity_function(m, hmod, cosm)

      ! Returns the dimensionless multiplicity function: M^2n(M)/rho; n(M) = dn/dM sometimes
      ! TODO: Is there a way to avoid an unstable numerical derivative here? Need dnu/dln(M)
      REAL, INTENT(IN) :: m                  ! Halo mass [Msun/h]
      TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL :: nu, dnu_dlnm
      INTEGER, PARAMETER :: iorder = iorder_derivative_mass_function
      INTEGER, PARAMETER :: ifind = ifind_derivative_mass_function

      IF(m == 0.) THEN
         multiplicity_function = 0.
      ELSE
         nu = nu_M(m, hmod, cosm)
         dnu_dlnm = derivative_table(log(m), log(hmod%m), hmod%nu, iorder, ifind)
         multiplicity_function = g_nu(nu, hmod)*dnu_dlnm
      END IF

   END FUNCTION multiplicity_function

   REAL FUNCTION halo_bias(m, hmod, cosm)

      ! Calculate b(M)
      REAL, INTENT(IN) :: m                  ! Halo mass [Msun/h]
      TYPE(halomod), INTENT(INOUT) :: hmod   ! Halo model
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL :: nu

      nu = nu_M(m, hmod, cosm)
      halo_bias = b_nu(nu, hmod)

   END FUNCTION halo_bias

   SUBROUTINE convert_mass_definitions(hmod, cosm)

      ! Create look-up tables for conversions between mass definitions
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: rhom, rhoc, Dv, a

      ! Scale factor
      a = hmod%a

      ! Get the densities
      rhom = comoving_matter_density(cosm)
      rhoc = comoving_critical_density(a, cosm)
      IF (hmod%mass_dependent_Dv) THEN
         STOP 'CONVERT_MASS_DEFINITIONS: Error, this does not work for a mass-dependent Delta_v'
      ELSE
         Dv = Delta_v(M0_Dv_default, hmod, cosm)
      END IF

      ! Calculate Delta = 200, 500 and Delta_c = 200, 500 quantities
      CALL convert_mass_definition(hmod%rv, hmod%c, hmod%m, Dv, 1., hmod%r500, hmod%c500, hmod%m500, 500., 1., hmod%n)
      CALL convert_mass_definition(hmod%rv, hmod%c, hmod%m, Dv, 1., hmod%r200, hmod%c200, hmod%m200, 200., 1., hmod%n)
      CALL convert_mass_definition(hmod%rv, hmod%c, hmod%m, Dv, rhom, hmod%r500c, hmod%c500c, hmod%m500c, 500., rhoc, hmod%n)
      CALL convert_mass_definition(hmod%rv, hmod%c, hmod%m, Dv, rhom, hmod%r200c, hmod%c200c, hmod%m200c, 200., rhoc, hmod%n)

      hmod%has_mass_conversions = .TRUE.

   END SUBROUTINE convert_mass_definitions

   SUBROUTINE convert_mass_definition(r1, c1, m1, D1, rho1, r2, c2, m2, D2, rho2, n)

      !Converts mass definition from Delta_1 rho_1 overdense to Delta_2 rho_2 overdense
      INTEGER, INTENT(IN) :: n   ! Number of entries in tables
      REAL, INTENT(IN) :: r1(n)  ! Array of initial virial radii [Mpc/h]
      REAL, INTENT(IN) :: c1(n)  ! Array of initial halo concentration
      REAL, INTENT(IN) :: m1(n)  ! Array of initial halo mass [Msun/h]
      REAL, INTENT(IN) :: D1     ! Initial halo overdensity definition (e.g., 200, 500)
      REAL, INTENT(IN) :: rho1   ! Initial halo overdensity defintion (critical or mass)
      REAL, INTENT(OUT) :: r2(n) ! Output array of virial radii with new definition [Mpc/h]
      REAL, INTENT(OUT) :: c2(n) ! Output array of halo concentration with new definition
      REAL, INTENT(OUT) :: m2(n) ! Output array of halo mass with new definition [Msun/h]
      REAL, INTENT(IN) :: D2     ! Final halo overdensity definition (e.g., 200, 500)
      REAL, INTENT(IN) :: rho2   ! Final halo overdensity defintion (critical or mass)
      REAL :: f(n)
      REAL :: rmin, rmax, rs
      REAL, ALLOCATABLE :: r(:)
      INTEGER :: i, j

      IF (verbose_convert_mass) THEN
         WRITE (*, *) 'CONVERT_MASS_DEFINITION: Converting mass definitions:'
         WRITE (*, *) 'CONVERT_MASS_DEFINITION: Initial overdensity:', D1
         WRITE (*, *) 'CONVERT_MASS_DEFINITION: Final overdensity:', D2
      END IF

      ! Make an array of general 'r' values for the solution later
      rmin = r1(1)/10. ! Should be sufficient for reasonable cosmologies
      rmax = r1(n)*10. ! Should be sufficient for reasonable cosmologies
      CALL fill_array(log(rmin), log(rmax), r, n) ! Necessary to log space
      r = exp(r)

      IF (verbose_convert_mass) THEN
         WRITE (*, *) 'CONVERT_MASS_DEFINITION: rmin:', rmin
         WRITE (*, *) 'CONVERT_MASS_DEFINITION: rmax:', rmax
         WRITE (*, *) 'CONVERT_MASS_DEFINITION: nr:', n
      END IF

      ! Now use the find algorithm to invert L(r_i)=R(r_j) so that r_j=R^{-1}[L(r_i)]
      IF (verbose_convert_mass) THEN
         WRITE (*, *) '========================================================================================================'
         WRITE (*, *) '         M_old         rv_old          c_old    M_new/M_old    r_new/r_old    c_new/c_old  M`_new/M`_old'
         WRITE (*, *) '========================================================================================================'
      END IF
      DO i = 1, n

         ! Calculate the halo scale radius
         rs = r1(i)/c1(i)

         ! Fill up the f(r) table which needs to be solved for R | f(R)=0
         DO j = 1, n
            f(j) = r1(i)**3*rho1*D1/NFW_factor(r1(i)/rs)-r(j)**3*rho2*D2/NFW_factor(r(j)/rs) ! NFW
            !f(j)=r1(i)**2*rho1*D1-r(j)**2*rho2*D2 !Isothemral sphere
         END DO

         ! First find the radius R | f(R)=0; I am fairly certain that I can use log on 'r' here
         r2(i) = exp(solve_find(log(r), f))

         ! Now do the concentration and mass conversions
         c2(i) = r2(i)/rs
         m2(i) = m1(i)*(rho2*D2*r2(i)**3)/(rho1*D1*r1(i)**3)
         !m2(i)=m1(i)*NFW_factor(c2(i))/NFW_factor(c1(i))

         IF (verbose_convert_mass) THEN
            WRITE (*, fmt='(ES15.7,6F15.7)') m1(i), r1(i), c1(i), m2(i)/m1(i), r2(i)/r1(i), c2(i)/c1(i), &
               NFW_factor(c2(i))/NFW_factor(c1(i))
         END IF

      END DO

      IF(verbose_convert_mass) THEN
         WRITE(*,*) '========================================================================================================'
      END IF

      IF (verbose_convert_mass) THEN
         WRITE (*, *) 'CONVERT_MASS_DEFINITION: Done'
         WRITE (*, *)
         STOP
      END IF

   END SUBROUTINE convert_mass_definition

   REAL FUNCTION NFW_factor(x)

      ! The NFW 'mass' factor that crops up all the time
      ! This is X(c) in M(r) = M X(r/rs) / X(c)
      REAL, INTENT(IN) :: x

      NFW_factor = log(1.+x)-x/(1.+x)

   END FUNCTION NFW_factor

   REAL FUNCTION virial_radius(m, hmod, cosm)

      ! Comoving halo virial radius
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      virial_radius = (3.*m/(4.*pi*comoving_matter_density(cosm)*Delta_v(m, hmod, cosm)))**(1./3.)

   END FUNCTION virial_radius

   REAL FUNCTION effective_index(flag_sigma, hmod, cosm)

      ! Power spectrum effective slope at the non-linear scale
      ! Defined as -3. + d ln sigma^2 / d ln r, so pertains to P(k) not Delta^2(k)
      INTEGER, INTENT(IN) :: flag_sigma
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      
      ! Find effective index at collapse
      effective_index = neff(hmod%rnl, hmod%a, flag_sigma, cosm)

      ! For some bizarre cosmologies r_nl is very small, so almost no collapse has occured
      ! In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
      IF (effective_index < cosm%ns-4.) effective_index = cosm%ns-4.
      IF (effective_index > cosm%ns)    effective_index = cosm%ns

   END FUNCTION effective_index

   REAL FUNCTION effective_curvature(flag_sigma, hmod, cosm)

      ! Power spectrum effective curvature at the non-linear scale
      ! Defined as d^2 ln sigma^2 / d ln r^2, so pertains to P(k) not Delta^2(k)
      INTEGER, INTENT(IN) :: flag_sigma
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      
      ! Effective curvature at the collapse scale
      effective_curvature = ncur(hmod%rnl, hmod%a, flag_sigma, cosm)

   END FUNCTION effective_curvature

   SUBROUTINE fill_halo_concentration(hmod, cosm)

      ! Fill look-up table with halo concentrations as a function of halo mass
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: mnl, m, zc, z, a, fg, rr, nu
      INTEGER :: i

      ! Get the redshift
      z = hmod%z
      a = hmod%a

      ! Any initialisation for the c(M) relation goes here
      IF (hmod%iconc == 1) THEN
         ! Fill the collapse z look-up table
         CALL zcoll_Bullock(z, hmod, cosm)
      ELSE IF (hmod%iconc == 2 .OR. hmod%iconc == 9) THEN
         mnl = hmod%mnl
      ELSE IF (hmod%iconc == 10) THEN
         fg = growth_rate(hmod%a, cosm)
      END IF

      ! Fill concentration-mass for all halo masses
      DO i = 1, hmod%n

         ! Choose concentration-mass relation
         IF (hmod%iconc == 1) THEN
            zc = hmod%zc(i)
            hmod%c(i) = conc_Bullock(z, zc)
         ELSE IF (hmod%iconc == 2) THEN
            m = hmod%m(i)
            hmod%c(i) = conc_Bullock_simple(m, mnl)
         ELSE IF (hmod%iconc == 3) THEN
            m = hmod%m(i)
            hmod%c(i) = conc_Duffy_full_M200(m, z)
         ELSE IF (hmod%iconc == 4) THEN
            m = hmod%m(i)
            hmod%c(i) = conc_Duffy_full_virial(m, z)
         ELSE IF (hmod%iconc == 5) THEN
            hmod%c(i) = conc_Duffy_full_M200c(m, z)
         ELSE IF (hmod%iconc == 6) THEN
            m = hmod%m(i)
            hmod%c(i) = conc_Duffy_relaxed_M200(m, z)
         ELSE IF (hmod%iconc == 7) THEN
            m = hmod%m(i)
            hmod%c(i) = conc_Duffy_relaxed_virial(m, z)
         ELSE IF (hmod%iconc == 8) THEN
            m = hmod%m(i)
            hmod%c(i) = conc_Duffy_relaxed_M200c(m, z)
         ELSE IF (hmod%iconc == 9) THEN
            m = hmod%m(i)
            hmod%c(i) = conc_Child(m, mnl)
         ELSE IF (hmod%iconc == 10) THEN
            nu = hmod%nu(i)
            rr = hmod%rr(i)
            hmod%c(i) = conc_Diemer(nu, rr, fg, a, cosm)
         ELSE
            STOP 'FILL_HALO_CONCENTRATION: Error, iconc specified incorrectly'
         END IF

         ! Rescale halo concentrations via the 'A' HMcode parameter
         hmod%c(i) = hmod%c(i)*hmod%HMcode_A

      END DO

      ! Dolag (2004) prescription for adding DE dependence
      IF (hmod%iDolag /= 0) CALL Dolag_correction(hmod, cosm)

      ! Rescale the concentration-mass relation for gas the epsilon parameter
      ! This only rescales the concentrations of haloes that *contain* substantial amounts of gas
      ! TODO: I feel this should be removed. It has the possibility to go very wrong, also unnecessary for non-hydro
      DO i = 1, hmod%n
         m = hmod%m(i)
         hmod%c(i) = hmod%c(i)*hydro_concentration_modification(m, hmod, cosm)
      END DO

   END SUBROUTINE fill_halo_concentration

   REAL FUNCTION hydro_concentration_modification(m, hmod, cosm)

      ! Fixed so that concentration modification applies to low-mass haloes
      ! Note very carefully that all modifications to the halo concentration must go in here
      ! This is very important because they need to be undone for DMONLY haloes
      ! Please be extremely careful
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: gas_fraction, eps1, eps2
      REAL, PARAMETER :: halo_conc_mod_min = 1e-2

      ! Fraction of original gas content remaining in halo
      gas_fraction = halo_bound_gas_fraction(m, hmod, cosm)/(cosm%Om_b/cosm%Om_m)
      eps1 = HMx_eps(hmod, cosm)  ! Gas-poor limit: c -> c*(1+e1)
      eps2 = HMx_eps2(hmod, cosm) ! Gas-rich limit: c -> c*(1+e2)
      hydro_concentration_modification = 1.+eps1+(eps2-eps1)*gas_fraction
      CALL fix_minimum(hydro_concentration_modification, halo_conc_mod_min)

   END FUNCTION hydro_concentration_modification

   SUBROUTINE Dolag_correction(hmod, cosm)

      ! Applies the Dolag et al. (2004; astro-ph/0309771) correction to concentration-mass relation
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: g_LCDM, g_wCDM, gc_LCDM, gc_wCDM, ginf_LCDM, ginf_wCDM
      REAL :: a, ac, zc, ainf, rc, rinf, r
      INTEGER :: i
      TYPE(cosmology) :: cosm_LCDM
      LOGICAL :: remove_neutrinos
      LOGICAL, PARAMETER :: make_lambda = .TRUE.
      LOGICAL, PARAMETER :: make_flat = .TRUE.

      IF (hmod%iDolag == 5) THEN

         IF (.NOT. allocated(hmod%zc)) STOP 'DOLAG_CORRECTION: Error, this will only work if collapse redshifts have been defined'

         ! Set the LCDM cosmology
         cosm_LCDM = convert_cosmology(cosm, make_lambda, make_flat, remove_neutrinos=.FALSE.)

         ! Calculate the current growth factors
         a = hmod%a
         g_wCDM = grow(a, cosm)
         g_LCDM = grow(a, cosm_LCDM)
         r = g_wCDM/g_LCDM

         ! Apply the correction using the collapse redshift of each halo individually
         DO i = 1, hmod%n
            zc = hmod%zc(i)*hmod%fD
            IF (zc > hmod%z) THEN
               ac = scale_factor_z(zc)
               gc_wCDM = grow(ac, cosm)
               gc_LCDM = grow(ac, cosm_LCDM)
               rc = gc_wCDM/gc_LCDM
               hmod%c(i) = hmod%c(i)*rc/r
            END IF
         END DO

      ELSE

         IF (hmod%z < hmod%zD) THEN

            ! The 'infinite' scale factor
            ainf = scale_factor_z(hmod%zD)

            ! Save the growth function in the current cosmology
            ginf_wCDM = grow(ainf, cosm)

            ! Make a flat LCDM cosmology and calculate growth
            IF (hmod%iDolag == 1 .OR. hmod%iDolag == 2 .OR. hmod%iDolag == 3) THEN
               remove_neutrinos = .FALSE.
            ELSE IF (hmod%iDolag == 4) THEN
               remove_neutrinos = .TRUE.
            ELSE
               STOP 'DOLAG_CORRECTION: Error, iDolag not specified correctly'
            END IF
            cosm_LCDM = convert_cosmology(cosm, make_lambda, make_flat, remove_neutrinos)

            ! Growth factor in LCDM at 'infinity' calculated using Linder approximation
            ginf_LCDM = grow_Linder(ainf, cosm_LCDM)

            ! Fractional difference compared to LCDM
            rinf = ginf_wCDM/ginf_LCDM

            IF (hmod%iDolag == 1) THEN
               ! Standard correction (HMcode 2015)
               hmod%c = hmod%c*rinf
            ELSE IF (hmod%iDolag == 2) THEN
               ! Changed this to a power of 1.5 in HMcode 2016, produces more accurate results for extreme DE
               hmod%c = hmod%c*rinf**1.5
            ELSE IF (hmod%iDolag == 3 .OR. hmod%iDolag == 4) THEN
               ! Correction with a sensible redshift dependence, which is otherwise missing from Dolag
               a = hmod%a
               g_wCDM = grow(a, cosm)
               g_LCDM = grow_Linder(a, cosm_LCDM)
               r = g_wCDM/g_LCDM
               hmod%c = hmod%c*rinf/r
            ELSE
               STOP 'DOLAG_CORRECTION: Error, iDolag specified incorrectly'
            END IF

         END IF

   END IF

   END SUBROUTINE Dolag_correction

   REAL FUNCTION conc_Bullock(z, zc)

      ! Calculate the Bullock (2001; astro-ph/9908159) concentration
      REAL, INTENT(IN) :: z
      REAL, INTENT(IN) :: zc
      REAL, PARAMETER :: A = 4. ! Pre-factor for Bullock relation

      conc_Bullock = A*(1.+zc)/(1.+z)

   END FUNCTION conc_Bullock

   SUBROUTINE zcoll_Bullock(z, hmod, cosm)

      ! This fills up the halo collapse redshift table as per Bullock relations
      ! TODO: Convert to using solve root-finding routines
      REAL, INTENT(IN) :: z
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: dc, af, zf, RHS, a, rf, sig, growz, f
      INTEGER :: i

      a = scale_factor_z(z)

      ! Fills up a table for sigma(fM) for Bullock c(m) relation
      IF (hmod%iconc == 1) THEN
         f = hmod%fM**(1./3.)
         DO i = 1, hmod%n
            rf = hmod%rr(i)*f
            sig = sigma(rf, a, hmod%flag_sigma, cosm)
            hmod%sigf(i) = sig
         END DO   
      ELSE
         STOP 'ZCOLL_BULLOCK: Error, halo concentration model specified incorrectly'
      END IF

      ! I don't think this is really consistent with dc varying as a function of z
      ! but the change will *probably* be *very* small
      dc = hmod%dc

      ! Calculate the growth function at the current scale factor
      growz = grow(a, cosm)

      ! Do numerical inversion
      DO i = 1, hmod%n

         ! TODO: Probably more consistent to use sigma(R,a) here
         RHS = dc*growz/hmod%sigf(i)

         IF (RHS > growz) THEN
            ! If the halo forms 'in the future' then set the formation z to the current z
            zf = z
         ELSE
            af = inverse_interpolator(RHS, cosm%grow)
            zf = redshift_a(af)
         END IF

         hmod%zc(i) = zf

      END DO

   END SUBROUTINE zcoll_Bullock

   FUNCTION conc_Bullock_simple(m, mstar)

      ! The simple concentration-mass relation from Bullock et al. (2001; astro-ph/9908159 equation 18)
      REAL :: conc_Bullock_simple
      REAL, INTENT(IN) :: m, mstar

      conc_Bullock_simple = 9.*(m/mstar)**(-0.13)

   END FUNCTION conc_Bullock_simple

   REAL FUNCTION conc_Duffy_full_M200c(m, z)

      ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: z

      REAL, PARAMETER :: m_piv = 2e12 ! Pivot mass [Msun/h]
      REAL, PARAMETER :: A = 5.71
      REAL, PARAMETER :: B = -0.084
      REAL, PARAMETER :: C = -0.47

      ! Equation (4) in 0804.2486, parameters from 4th row of Table 1
      conc_Duffy_full_M200c = A*(m/m_piv)**B*(1.+z)**C

   END FUNCTION conc_Duffy_full_M200c

   REAL FUNCTION conc_Duffy_relaxed_M200c(m, z)

      ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: z

      REAL, PARAMETER :: m_piv = 2e12 ! Pivot mass [Msun/h]
      REAL, PARAMETER :: A = 6.71
      REAL, PARAMETER :: B = -0.091
      REAL, PARAMETER :: C = -0.44

      ! Equation (4) in 0804.2486, parameters from 4th row of Table 1
      conc_Duffy_relaxed_M200c = A*(m/m_piv)**B*(1.+z)**C

   END FUNCTION conc_Duffy_relaxed_M200c

   REAL FUNCTION conc_Duffy_full_virial(m, z)

      ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: z

      REAL, PARAMETER :: m_piv = 2e12 ! Pivot mass [Msun/h]
      REAL, PARAMETER :: A = 7.85
      REAL, PARAMETER :: B = -0.081
      REAL, PARAMETER :: C = -0.71

      ! Equation (4) in 0804.2486, parameters from 6th row of Table 1
      conc_Duffy_full_virial = A*(m/m_piv)**B*(1.+z)**C

   END FUNCTION conc_Duffy_full_virial

   REAL FUNCTION conc_Duffy_relaxed_virial(m, z)

      ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: z

      REAL, PARAMETER :: m_piv = 2e12 ! Pivot mass [Msun/h]
      REAL, PARAMETER :: A = 9.23
      REAL, PARAMETER :: B = -0.090
      REAL, PARAMETER :: C = -0.69

      ! Equation (4) in 0804.2486, parameters from 6th row of Table 1
      conc_Duffy_relaxed_virial = A*(m/m_piv)**B*(1.+z)**C

   END FUNCTION conc_Duffy_relaxed_virial

   REAL FUNCTION conc_Duffy_full_M200(m, z)

      ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: z

      REAL, PARAMETER :: m_piv = 2e12 ! Pivot mass [Msun/h]
      REAL, PARAMETER :: A = 10.14
      REAL, PARAMETER :: B = -0.081
      REAL, PARAMETER :: C = -1.01

      ! Equation (4) in 0804.2486, parameters from 10th row of Table 1
      conc_Duffy_full_M200 = A*(m/m_piv)**B*(1.+z)**C

   END FUNCTION conc_Duffy_full_M200

   REAL FUNCTION conc_Duffy_relaxed_M200(m, z)

      ! Duffy et al (2008; 0804.2486) c(M) relation for WMAP5, See Table 1
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: z

      REAL, PARAMETER :: m_piv = 2e12 ! Pivot mass [Msun/h]
      REAL, PARAMETER :: A = 11.93
      REAL, PARAMETER :: B = -0.090
      REAL, PARAMETER :: C = -0.99

      ! Equation (4) in 0804.2486, parameters from 10th row of Table 1
      conc_Duffy_relaxed_M200 = A*(m/m_piv)**B*(1.+z)**C

   END FUNCTION conc_Duffy_relaxed_M200

   REAL FUNCTION conc_Child(m, mnl)

      ! Concentration-mass relation from Child et al. (2018; 1804.10199); equation (18)
      ! Applicable for M200c
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: mnl
      REAL :: f
      REAL, PARAMETER :: A = 4.61
      REAL, PARAMETER :: b = 638.65
      REAL, PARAMETER :: p = -0.07
      REAL, PARAMETER :: c0 = 3.59

      f = (m/mnl)/b

      conc_Child = A*((f/(1.+f))**p-1.)+c0

   END FUNCTION conc_Child

   REAL FUNCTION conc_Diemer(nu, r, fg, a, cosm)

      ! Diemer & Joyce (2019) concentration-mass relation
      ! Parameters taken from 'mean' column of Table 2
      ! Relation is for M200c halo-boundary definition
      REAL, INTENT(IN) :: nu
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: fg
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: bigA, bigB, bigC, alpha, n_eff
      REAL, PARAMETER :: kappa = 0.42
      REAL, PARAMETER :: a0 = 2.37
      REAL, PARAMETER :: a1 = 1.74
      REAL, PARAMETER :: b0 = 3.39
      REAL, PARAMETER :: b1 = 1.82
      REAL, PARAMETER :: ca = 0.20
      INTEGER, PARAMETER :: flag_power = flag_matter

      ! Suppress warnings
      alpha = nu

      STOP 'CONC_DIEMER: Need to implement this effectively, need neff(M)'
      n_eff = neff(kappa*r, a, flag_power, cosm)

      bigA = a0*(1.+a1*(n_eff+3.))
      bigB = b0*(1.+b1*(n_eff+3.))
      bigC = 1.-ca*(1.-fg)
      
      conc_Diemer = 1.

   END FUNCTION conc_Diemer

   REAL FUNCTION baryonify_wk(wk, m, hmod, cosm)

      ! Change DMONLY halo window to account for feedback mass loss and star formation
      REAL, INTENT(IN) :: wk
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: wkn, sb

      wkn = wk
      sb = HMcode_sbar(cosm, hmod)
      wkn = wkn*DMONLY_halo_mass_fraction(m, hmod, cosm) ! Account for gas expulsion
      wkn = wkn+sb*m/comoving_matter_density(cosm)       ! Account for star formation

      baryonify_wk = wkn

   END FUNCTION baryonify_wk

   REAL FUNCTION unbaryonify_wk(wk, m, hmod, cosm)

      ! Inverse of baryonify_wk
      ! This is only used if windows need to be rebuilt for the two-halo term
      ! It will not be used if the two-halo term is linear theory, or dewiggled linear
      REAL, INTENT(IN) :: wk
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: wko, sb

      wko = wk
      sb = HMcode_sbar(cosm, hmod)
      wko = (wko-sb*m/comoving_matter_density(cosm))     ! Subtract stars
      wko = wko/DMONLY_halo_mass_fraction(m, hmod, cosm) ! Remove gas expulsion

      unbaryonify_wk = wko

   END FUNCTION unbaryonify_wk

   REAL FUNCTION DMONLY_halo_mass_fraction(m, hmod, cosm)

      ! Simple baryon model where high-mass haloes have a mass fraction of 1 and low-mass haloes have Omega_c/Omega_m
      ! This also accounts for massive neutrinos since it is written in terms of Om_c and Om_b (rather than Om_m)
      ! TODO: This is independent of k, so probably should be computed outside and k dependent function for speed
      ! TODO: Could just precompute this once for each M in the halomod init function
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(IN) :: hmod
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: r, fc, fb, fs, mb, n

      mb = HMcode_mbar(cosm, hmod)
      n = hmod%nbar ! Power-law index
      r = (m/mb)**n ! If m>>m0 then r becomes large, if m<<m0 then r=0  
      fb = cosm%Om_b/cosm%Om_m!-HMcode_sbar(cosm, hmod) ! Halo baryon fraction                              
      fc = cosm%Om_c/cosm%Om_m ! Halo CDM fraction
      fs = HMcode_sbar(cosm, hmod) ! Halo star fraction
      DMONLY_halo_mass_fraction = fc+(fb-fs)*r/(1.+r) ! Remaining halo mass fraction

   END FUNCTION DMONLY_halo_mass_fraction

   REAL FUNCTION HMcode_DMONLY_baryon_model(T, a1, a0)

      ! Scaling for the HMcode-2020 temperature dependence of the baryon parameters
      REAL, INTENT(IN) :: T
      REAL, INTENT(IN) :: a1
      REAL, INTENT(IN) :: a0
      REAL, PARAMETER :: Tmid = 7.8

      HMcode_DMONLY_baryon_model = a1*(log10(T)-Tmid)+a0

   END FUNCTION HMcode_DMONLY_baryon_model

   REAL FUNCTION HMcode_mbar(cosm, hmod)

      ! HMcode-2020 baryonic halo mass
      TYPE(cosmology), INTENT(IN) :: cosm
      TYPE(halomod) :: hmod
      REAL :: mbar, mbarz, z

      mbar = 10.**HMcode_DMONLY_baryon_model(cosm%Theat, hmod%mbar_T, log10(hmod%mbar))
      mbarz = HMcode_DMONLY_baryon_model(cosm%Theat, hmod%mbarz_T, hmod%mbarz)
      z = hmod%z
      HMcode_mbar = mbar*10**(mbarz*z)

   END FUNCTION HMcode_mbar

   REAL FUNCTION HMcode_sbar(cosm, hmod)

      ! HMcode-2020 star formation fraction
      TYPE(cosmology), INTENT(IN) :: cosm
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: sbar, sbarz, z

      sbar = HMcode_DMONLY_baryon_model(cosm%Theat, hmod%sbar_T, hmod%sbar)
      sbarz = HMcode_DMONLY_baryon_model(cosm%Theat, hmod%sbarz_T, hmod%sbarz)
      z = hmod%z
      HMcode_sbar = sbar*10**(z*sbarz)

      CALL fix_maximum(HMcode_sbar, cosm%Om_b/cosm%Om_m)

   END FUNCTION HMcode_sbar

   REAL FUNCTION win(real_space, ifield, k, m, rv, rs, hmod, cosm)

      ! Selects the halo profile type
      LOGICAL, INTENT(IN) :: real_space
      INTEGER, INTENT(IN) :: ifield
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: nu, mmin, mmax

      IF (ifield == field_dmonly) THEN
         win = win_DMONLY(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_neutrino) THEN
         win = win_neutrino(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_matter) THEN
         win = win_matter(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_cdm) THEN
         win = win_CDM(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_gas) THEN
         win = win_gas(real_space, ifield, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_stars) THEN
         win = win_stars(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_bound_gas) THEN
         win = win_bound_gas(real_space, ifield, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_ejected_gas) THEN
         win = win_ejected_gas(real_space, ifield, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_electron_pressure) THEN
         win = win_electron_pressure(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_void) THEN
         win = win_void(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_compensated_void) THEN
         win = win_compensated_void(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_central_galaxies) THEN
         win = win_centrals(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_satellite_galaxies) THEN
         win = win_satellites(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_galaxies) THEN
         win = win_galaxies(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_HI) THEN
         win = win_HI(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_cold_bound_gas) THEN
         win = win_cold_bound_gas(real_space, ifield, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_hot_bound_gas) THEN
         win = win_hot_bound_gas(real_space, ifield, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_normal_bound_gas) THEN
         win = win_normal_bound_gas(real_space, ifield, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_central_stars) THEN
         win = win_central_stars(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_satellite_stars) THEN
         win = win_satellite_stars(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_CIB_353 .OR. ifield == field_CIB_545 .OR. ifield == field_CIB_857) THEN
         IF (ifield == field_CIB_353) THEN
            nu = 353.e9 ! Frequency [Hz]
         ELSE IF (ifield == field_CIB_545) THEN
            nu = 545.e9 ! Frequency [Hz]
         ELSE IF (ifield == field_CIB_857) THEN
            nu = 857.e9 ! Frequency [Hz]
         ELSE
            STOP 'WIN_TYPE: Error, ifield specified incorrectly'
         END IF
         win = win_CIB(real_space, nu, k, m, rv, rs, hmod, cosm)
      ELSE IF (ifield == field_halo_11p0_11p5 .OR. ifield == field_halo_11p5_12p0 .OR. &
               ifield == field_halo_12p0_12p5 .OR. ifield == field_halo_12p5_13p0 .OR. &
               ifield == field_halo_13p0_13p5 .OR. ifield == field_halo_13p5_14p0 .OR. &
               ifield == field_halo_14p0_14p5 .OR. ifield == field_halo_14p5_15p0) THEN
         IF (ifield == field_halo_11p0_11p5) THEN
            mmin = 10**11.0 ! Minimum halo mass [Msun/h]
            mmax = 10**11.5 ! Maximum halo mass [Msun/h]
         ELSE IF (ifield == field_halo_11p5_12p0) THEN
            mmin = 10**11.5 ! Minimum halo mass [Msun/h]
            mmax = 10**12.0 ! Maximum halo mass [Msun/h]
         ELSE IF (ifield == field_halo_12p0_12p5) THEN
            mmin = 10**12.0 ! Minimum halo mass [Msun/h]
            mmax = 10**12.5 ! Maximum halo mass [Msun/h]
         ELSE IF (ifield == field_halo_12p5_13p0) THEN
            mmin = 10**12.5 ! Minimum halo mass [Msun/h]
            mmax = 10**13.0 ! Maximum halo mass [Msun/h]
         ELSE IF (ifield == field_halo_13p0_13p5) THEN
            mmin = 10**13.0 ! Minimum halo mass [Msun/h]
            mmax = 10**13.5 ! Maximum halo mass [Msun/h]
         ELSE IF (ifield == field_halo_13p5_14p0) THEN
            mmin = 10**13.5 ! Minimum halo mass [Msun/h]
            mmax = 10**14.0 ! Maximum halo mass [Msun/h]
         ELSE IF (ifield == field_halo_14p0_14p5) THEN
            mmin = 10**14.0 ! Minimum halo mass [Msun/h]
            mmax = 10**14.5 ! Maximum halo mass [Msun/h]
         ELSE IF (ifield == field_halo_14p5_15p0) THEN
            mmin = 10**14.5 ! Minimum halo mass [Msun/h]
            mmax = 10**15.0 ! Maximum halo mass [Msun/h]
         ELSE
            STOP 'WIN_TYPE: Error, ifield specified incorrectly'
         END IF
         win = win_haloes(real_space, mmin, mmax, k, m, rv, rs, hmod, cosm)
      ELSE
         WRITE (*, *) 'WIN_TYPE: ifield:', ifield
         STOP 'WIN_TYPE: Error, ifield specified incorreclty'
      END IF

   END FUNCTION win

   REAL FUNCTION win_DMONLY(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for all matter under the assumption that it is all CDM
      ! TODO: Possibly remove this since DMONLY should really be set by the halo_model ihm
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, rss, p1, p2

      ! Set additional halo parameters to zero
      p1 = 0.
      p2 = 0.

      rmin = 0.
      rmax = rv

      IF (hmod%halo_DMONLY == 1) THEN
         ! Analytical NFW
         irho = 5
      ELSE IF (hmod%halo_DMONLY == 2) THEN
         ! Non-analyical NFW
         irho = 4
      ELSE IF (hmod%halo_DMONLY == 3) THEN
         ! Tophat
         irho = 2
      ELSE IF (hmod%halo_DMONLY == 4) THEN
         ! Delta function
         irho = 0
      ELSE IF (hmod%halo_DMONLY == 5) THEN
         ! Cored NFW
         irho = 24
         p1 = hmod%rcore
      ELSE IF (hmod%halo_DMONLY == 6) THEN
         ! Isothermal
         irho = 1
      ELSE IF (hmod%halo_DMONLY == 7) THEN
         ! Shell
         irho = 28
      ELSE
         STOP 'WIN_DMONLY: Error, halo_DMONLY specified incorrectly'
      END IF

      ! Force it to use the gravity-only concentration relation (unapply gas correction)
      ! TODO: This is super ugly and should be improved somehow; also is unncessary calculations so slow
      ! TODO: Somehow should store modified and unmodified halo concentrations
      ! TODO: Or maybe the DMONLY option should be removed? This would probably be the more sensible thing to do.
      rss = rs*hydro_concentration_modification(m, hmod, cosm) ! This is *very* important

      IF (real_space) THEN
         r = k
         win_DMONLY = rho(r, rmin, rmax, rv, rss, p1, p2, irho)
         win_DMONLY = win_DMONLY/normalisation(rmin, rmax, rv, rss, p1, p2, irho)
      ELSE
         !Properly normalise and convert to overdensity
         win_DMONLY = m*win_norm(k, rmin, rmax, rv, rss, p1, p2, irho)/comoving_matter_density(cosm)
      END IF

      ! Fudges for DMONLY halo models to account for baryons and/or massive neutrinos
      ! These are not k dependent so should probably not be computed here
      IF (hmod%DMONLY_baryon_recipe) THEN
         ! This automatically accounts for massive neutrinos too
         win_DMONLY = baryonify_wk(win_DMONLY, m, hmod, cosm)
      ELSE IF (hmod%DMONLY_neutrino_halo_mass_correction) THEN
         win_DMONLY = win_DMONLY*(1.-cosm%f_nu)
      END IF

   END FUNCTION win_DMONLY

   REAL FUNCTION win_matter(real_space, k, m, rv, rs, hmod, cosm)

      ! The halo profile of all the matter density
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: CDM, gas, stars, neutrino

      CDM = win_CDM(real_space, k, m, rv, rs, hmod, cosm)
      gas = win_gas(real_space, field_gas, k, m, rv, rs, hmod, cosm)
      stars = win_stars(real_space, k, m, rv, rs, hmod, cosm)
      neutrino = win_neutrino(real_space, k, m, rv, rs, hmod, cosm)

      win_matter = CDM+gas+stars+neutrino

   END FUNCTION win_matter

   REAL FUNCTION win_CDM(real_space, k, m, rv, rs, hmod, cosm)

      ! The halo profile for CDM density
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2, frac

      frac = halo_CDM_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_CDM = 0.

      ELSE

         ! Default maximum and minimum radii
         rmin = 0.
         rmax = rv

         ! Default additional halo parameters
         p1 = 0.
         p2 = 0.

         IF (hmod%halo_CDM == 1) THEN
            irho = 5 ! Analytical NFW
         ELSE
            STOP 'WIN_CDM: Error, halo_CDM specified incorrectly'
         END IF

         IF (real_space) THEN
            r = k
            win_CDM = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
            win_CDM = win_CDM/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
         ELSE
            ! Properly normalise and convert to overdensity
            win_CDM = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)/comoving_matter_density(cosm)
         END IF

         win_CDM = frac*win_CDM

      END IF

   END FUNCTION win_CDM

   REAL FUNCTION win_gas(real_space, itype, k, m, rv, rs, hmod, cosm)

      ! Halo profile for gas density
      LOGICAL, INTENT(IN) :: real_space
      INTEGER, INTENT(IN) :: itype
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: win_bound, win_free

      win_bound = win_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)
      win_free = win_ejected_gas(real_space, itype, k, m, rv, rs, hmod, cosm)
      win_gas = win_bound+win_free

   END FUNCTION win_gas

   REAL FUNCTION win_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)

      ! Halo profile for totality of bound-gas density
      LOGICAL, INTENT(IN) :: real_space
      INTEGER, INTENT(IN) :: itype
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: win_static, win_cold, win_hot

      win_static = win_normal_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)
      win_cold = win_cold_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)
      win_hot = win_hot_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)

      win_bound_gas = win_static+win_cold+win_hot

   END FUNCTION win_bound_gas

   REAL FUNCTION win_normal_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)

      ! Halo profile for the standard bound-gas density
      LOGICAL, INTENT(IN) :: real_space
      INTEGER, INTENT(IN) :: itype
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: rho0, T0, r, a
      REAL :: rmin, rmax, p1, p2, frac
      INTEGER :: irho_density, irho_pressure

      frac = halo_normal_bound_gas_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_normal_bound_gas = 0.

      ELSE

         ! Initially set the halo parameters to zero
         p1 = 0.
         p2 = 0.

         ! Set maximum and minimum integration radius
         rmin = 0.
         rmax = rv

         IF (frac < frac_min_delta) THEN
            ! Treat as delta functions if there is not much abundance
            irho_density = 0
            irho_pressure = 0
         ELSE IF (hmod%halo_normal_bound_gas == 1 .OR. hmod%halo_normal_bound_gas == 3) THEN
            ! Komatsu & Seljak (2001) profile
            IF (hmod%halo_normal_bound_gas == 1) THEN
               ! Simplified KS model
               irho_density = 11
               irho_pressure = 13
            ELSE IF (hmod%halo_normal_bound_gas == 3) THEN
               ! Full KS model
               irho_density = 21
               irho_pressure = 23
            END IF
            p1 = HMx_Gamma(m, hmod, cosm)
            p2 = HMx_Zamma(m, hmod, cosm)
         ELSE IF (hmod%halo_normal_bound_gas == 2) THEN
            ! Set cored isothermal profile
            !irho_density=6   ! Isothermal beta model with beta=2/3
            irho_density = 15 ! Isothermal beta model with general beta
            p1 = HMx_ibeta(m, hmod, cosm)
            p2 = p1
            irho_pressure = irho_density ! Okay to use density for pressure because temperature is constant (isothermal)
         ELSE IF (hmod%halo_normal_bound_gas == 4) THEN
            ! NFW
            irho_density = 5 ! Analytical NFW
            irho_pressure = 0
         ELSE
            STOP 'WIN_STATIC_GAS: Error, halo_normal_bound_gas not specified correctly'
         END IF

         IF (itype == field_matter .OR. itype == field_gas .OR. itype == field_bound_gas .OR. itype == field_normal_bound_gas) THEN

            ! Density profile of bound gas
            IF (real_space) THEN
               r = k
               win_normal_bound_gas = rho(r, rmin, rmax, rv, rs, p1, p2, irho_density)
               win_normal_bound_gas = win_normal_bound_gas/normalisation(rmin, rmax, rv, rs, p1, p2, irho_density)
            ELSE
               ! Properly normalise and convert to overdensity
               win_normal_bound_gas = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho_density)/comoving_matter_density(cosm)
            END IF

            win_normal_bound_gas = frac*win_normal_bound_gas

         ELSE IF (itype == field_electron_pressure) THEN

            ! Electron pressure profile of bound gas
            ! NOTE: Swapped p1, p2 to p2, p1 here to get the gas index correct in the case that the Gammas differ
            ! NOTE: This is a bit of a hack and there is probably a much cuter solution 
            IF (real_space) THEN
               r = k
               win_normal_bound_gas = rho(r, rmin, rmax, rv, rs, p2, p1, irho_pressure)
            ELSE
               ! The electron pressure window is T(r) x rho_e(r), we want unnormalised, so multiply through by normalisation
               ! TODO: Can I make the code more efficient here by having an unnorm window function?
               win_normal_bound_gas = win_norm(k, rmin, rmax, rv, rs, p2, p1, irho_pressure)
               win_normal_bound_gas = win_normal_bound_gas*normalisation(rmin, rmax, rv, rs, p2, p1, irho_pressure)
            END IF

            ! Calculate the value of the density profile prefactor and change units from cosmological to SI
            rho0 = m*frac/normalisation(rmin, rmax, rv, rs, p1, p2, irho_density)
            rho0 = rho0*msun/mpc/mpc/mpc ! Overflow with 4-byte real numbers if you use mpc**3
            rho0 = rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

            ! Calculate the value of the temperature prefactor [K]
            a = hmod%a
            T0 = HMx_alpha(m, hmod, cosm)*virial_temperature(m, rv, hmod%a, cosm)

            ! Convert from Temp x density -> electron prsessure (Temp x n; n is all particle number density)
            win_normal_bound_gas = win_normal_bound_gas*(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply by *number density* times temperature k_B [J/m^3]
            win_normal_bound_gas = win_normal_bound_gas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
            win_normal_bound_gas = win_normal_bound_gas*cosm%mup/cosm%mue ! Convert from total thermal pressure to electron pressure

         ELSE

            STOP 'WIN_STATIC_GAS: Error, itype not specified correctly'

         END IF

      END IF

   END FUNCTION win_normal_bound_gas

   REAL FUNCTION win_cold_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)

      ! Halo profile for the cold-gas density
      LOGICAL, INTENT(IN) :: real_space
      INTEGER, INTENT(IN) :: itype
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r, rmin, rmax, p1, p2, frac
      INTEGER :: irho

      frac = halo_cold_bound_gas_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_cold_bound_gas = 0.

      ELSE

         ! Initially set the halo parameters to zero
         p1 = 0.
         p2 = 0.

         ! Set maximum and minimum integration radius
         rmin = 0.
         rmax = rv

         IF (frac < frac_min_delta) THEN
            ! Treat as delta functions if there is not much abundance
            irho = 0
         ELSE IF (hmod%halo_cold_bound_gas == 1) THEN
            ! Delta function
            irho = 0
         ELSE
            STOP 'WIN_COLD_GAS: Error, halo_cold_bound_gas not specified correctly'
         END IF

         IF (itype == field_matter .OR. itype == field_gas .OR. itype == field_bound_gas .OR. itype == field_cold_bound_gas) THEN

            ! Density profile of cold gas
            IF (real_space) THEN
               r = k
               win_cold_bound_gas = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
               win_cold_bound_gas = win_cold_bound_gas/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
            ELSE
               ! Properly normalise and convert to overdensity
               win_cold_bound_gas = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)/comoving_matter_density(cosm)
            END IF

            win_cold_bound_gas = frac*win_cold_bound_gas

         ELSE IF (itype == field_electron_pressure) THEN

            ! No electron-pressure contribution from the cold gas
            win_cold_bound_gas = 0.

         ELSE

            STOP 'WIN_COLD_GAS: Error, itype not specified correctly'

         END IF

      END IF

   END FUNCTION win_cold_bound_gas

   REAL FUNCTION win_hot_bound_gas(real_space, itype, k, m, rv, rs, hmod, cosm)

      ! Halo profile for hot-bound gas density
      LOGICAL, INTENT(IN) :: real_space
      INTEGER, INTENT(IN) :: itype
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r, rmin, rmax, p1, p2, frac
      INTEGER :: irho_density, irho_pressure
      REAL :: rho0, a, T0

      frac = halo_hot_bound_gas_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_hot_bound_gas = 0.

      ELSE

         ! Initially set the halo parameters to zero
         p1 = 0.
         p2 = 0.

         ! Set maximum and minimum integration radius
         rmin = 0.
         rmax = rv

         IF (frac < frac_min_delta) THEN
            ! Treat as delta functions if there is not much abundance
            irho_density = 0
            irho_pressure = 0
         ELSE IF (hmod%halo_hot_bound_gas == 1) THEN
            ! Isothermal
            irho_density = 1
            irho_pressure = 1
         ELSE
            STOP 'WIN_HOT_GAS: Error, halo_hot_bound_gas not specified correctly'
         END IF

         IF (itype == field_matter .OR. itype == field_gas .OR. itype == field_bound_gas .OR. itype == field_hot_bound_gas) THEN

            ! Density profile of hot gas
            IF (real_space) THEN
               r = k
               win_hot_bound_gas = rho(r, rmin, rmax, rv, rs, p1, p2, irho_density)
               win_hot_bound_gas = win_hot_bound_gas/normalisation(rmin, rmax, rv, rs, p1, p2, irho_density)
            ELSE
               ! Properly normalise and convert to overdensity
               win_hot_bound_gas = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho_density)/comoving_matter_density(cosm)
            END IF

            win_hot_bound_gas = frac*win_hot_bound_gas

         ELSE IF (itype == field_electron_pressure) THEN

            ! Electron-pressure profile of bound gas
            IF (real_space) THEN
               r = k
               win_hot_bound_gas = rho(r, rmin, rmax, rv, rs, p1, p2, irho_pressure)
            ELSE
               ! The electron pressure window is T(r) x rho_e(r), we want unnormalised, so multiply through by normalisation
               ! TODO: Can I make the code more efficient here by having an un-normalised window function?
               win_hot_bound_gas = win_norm(k, rmin, rmax, rv, rs, p1, p2, irho_pressure)
               win_hot_bound_gas = win_hot_bound_gas*normalisation(rmin, rmax, rv, rs, p1, p2, irho_pressure)
            END IF

            ! Calculate the value of the density profile prefactor and change units from cosmological to SI
            rho0 = m*frac/normalisation(rmin, rmax, rv, rs, p1, p2, irho_density)
            rho0 = rho0*msun/mpc/mpc/mpc ! Overflow with REAL*4 if you use mpc**3
            rho0 = rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

            ! Calculate the value of the temperature prefactor [K]
            a = hmod%a
            T0 = HMx_beta(m, hmod, cosm)*virial_temperature(m, rv, hmod%a, cosm)

            ! Convert from Temp x density -> electron pressure (Temp x n; n is all particle number density)
            win_hot_bound_gas = win_hot_bound_gas*(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply by number density times temperature k_B [J/m^3]
            win_hot_bound_gas = win_hot_bound_gas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
            win_hot_bound_gas = win_hot_bound_gas*cosm%mup/cosm%mue ! Convert from total thermal pressure to electron pressure

         ELSE

            STOP 'WIN_HOT_GAS: Error, itype not specified correctly'

         END IF

      END IF

   END FUNCTION win_hot_bound_gas

   REAL FUNCTION win_ejected_gas(real_space, itype, k, m, rv, rs, hmod, cosm)

      ! Halo profile for the ejected gas component
      LOGICAL, INTENT(IN) :: real_space
      INTEGER, INTENT(IN) :: itype
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: re, rmin, rmax, r, A, rho0, rhov, T0, p1, p2, beta, c, thing, m0, frac
      INTEGER :: irho_density, irho_pressure

      ! Enable to force the electron pressure to be matched at the virial radius
      ! This is enabled by default for some halo gas/pressure models
      LOGICAL :: match_electron_pressure = .FALSE.

      frac = halo_ejected_gas_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         ! Sometimes the free-gas fraction will be zero, in which case this avoids problems
         win_ejected_gas = 0.

      ELSE

         ! Set the halo 'parameter' variables to zero initially
         p1 = 0.
         p2 = 0.

         ! Set default min and max radii
         rmin = 0.
         rmax = rv

         IF (frac < frac_min_delta) THEN

            ! Treat as delta functions if there is not much abundance
            irho_density = 0
            irho_pressure = 0

         ELSE IF (hmod%halo_ejected_gas == 1) THEN

            ! Simple isothermal model, motivated by constant velocity and rate expulsion
            irho_density = 1
            irho_pressure = irho_density ! Okay because T is constant
            rmin = 0.
            rmax = 2.*rv

         ELSE IF (hmod%halo_ejected_gas == 2) THEN

            ! Ejected gas model from Schneider & Teyssier (2015)
            irho_density = 10
            irho_pressure = irho_density ! Okay because T is constant
            rmin = rv
            re = rv
            p1 = re
            rmax = 15.*re ! Needs to be such that integral converges (15rf seems okay)

         ELSE IF (hmod%halo_ejected_gas == 3) THEN

            ! Now do isothermal shell connected to the KS profile continuously
            irho_density = 16
            irho_pressure = irho_density ! Okay because T is constant

            ! Isothermal model with continuous link to KS
            rhov = win_normal_bound_gas(.TRUE., 2, rv, m, rv, rs, hmod, cosm) ! value of the density at the halo boundary for bound gas
            A = rhov/rho(rv, 0., rv, rv, rs, p1, p2, irho_density) ! This is A, as in A/r^2

            rmin = rv
            rmax = rv+frac/(4.*pi*A) ! This ensures density continuity and mass conservation

            c = 10. ! How many times larger than the virial radius can the gas cloud go?
            IF (rmax > c*rv) rmax = c*rv ! This needs to be set otherwise get huge decrement in gas power at large scales
            match_electron_pressure = .TRUE. ! Match the electron pressure at the boundary

         ELSE IF (hmod%halo_ejected_gas == 4) THEN

            ! Ejected gas is a continuation of the KS profile
            irho_density = 11 ! KS
            irho_pressure = 13 ! KS
            rmin = rv
            rmax = 2.*rv
            p1 = HMx_Gamma(m, hmod, cosm)

         ELSE IF (hmod%halo_ejected_gas == 5) THEN

            m0 = 1e14

            IF (m < m0) THEN

               irho_density = 0
               irho_pressure = irho_density
               rmin = 0.
               rmax = rv

            ELSE

               ! Set the density profile to be the power-law profile
               irho_density = 17
               irho_pressure = irho_density ! Not okay

               ! Calculate the KS index at the virial radius
               c = rv/rs
               beta = (c-(1.+c)*log(1.+c))/((1.+c)*log(1.+c))
               beta = beta/(HMx_Gamma(m, hmod, cosm)-1.) ! This is the power-law index at the virial radius for the KS gas profile
               p1 = beta
               IF (beta <= -3.) beta = -2.9 ! If beta<-3 then there is only a finite amount of gas allowed in the free component

               ! Calculate the density at the boundary of the KS profile
               rhov = win_normal_bound_gas(.TRUE., 2, rv, m, rv, rs, hmod, cosm)

               ! Calculate A as in rho(r)=A*r**beta
               A = rhov/rho(rv, 0., rv, rv, rs, p1, p2, irho_density)

               ! Set the minimum radius for the power-law to be the virial radius
               rmin = rv

               ! Set the maximum radius so that it joins to KS profile seamlessly
               thing = (beta+3.)*frac/(4.*pi*A)+(rhov*rv**3)/A
               IF (thing > 0.) THEN
                  ! This then fixes the condition of contiunity in amplitude and gradient
                  rmax = thing**(1./(beta+3.))
               ELSE
                  ! If there are no sohmodions then fix to 10rv and accept discontinuity
                  ! There may be no sohmodion if there is a lot of free gas and if beta<-3
                  rmax = 10.*rv
               END IF

            END IF

         ELSE IF (hmod%halo_ejected_gas == 6) THEN

            ! Cubic profile
            rmin = rv
            rmax = 3.*rv
            irho_density = 18
            irho_pressure = irho_density

         ELSE IF (hmod%halo_ejected_gas == 7) THEN

            ! Smooth profile (rho=0)
            rmin = 0.
            rmax = rv
            irho_density = 19
            irho_pressure = irho_density

         ELSE IF (hmod%halo_ejected_gas == 8) THEN

            ! Delta function
            rmin = 0.
            rmax = rv
            irho_density = 0
            irho_pressure = irho_density

         ELSE
            STOP 'WIN_FREE_GAS: Error, halo_ejected_gas specified incorrectly'
         END IF

         ! Density profile
         IF (itype == field_matter .OR. itype == field_gas .OR. itype == field_ejected_gas) THEN

            ! Density profile of free gas
            IF (real_space) THEN
               r = k
               win_ejected_gas = rho(r, rmin, rmax, rv, rs, p1, p2, irho_density)
               win_ejected_gas = win_ejected_gas/normalisation(rmin, rmax, rv, rs, p1, p2, irho_density)
            ELSE
               ! Properly normalise and convert to overdensity
               win_ejected_gas = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho_density)/comoving_matter_density(cosm)
            END IF

            win_ejected_gas = frac*win_ejected_gas

            ! Electron pressure profile
         ELSE IF (itype == field_electron_pressure) THEN

            ! If we are applying a pressure-matching condition
            IF (match_electron_pressure) THEN

               STOP 'WIN_FREE_GAS: Check the units and stuff here *very* carefully'

               r = k
               IF (r > rmin .AND. r < rmax) THEN
                  ! Only works for isothermal profile
                  win_ejected_gas = win_normal_bound_gas(.TRUE., 6, rv, m, rv, rs, hmod, cosm)*(r/rv)**(-2)
               ELSE
                  win_ejected_gas = 0.
               END IF

            ELSE

               ! Electron pressure profile of free gas
               IF (real_space) THEN
                  r = k
                  win_ejected_gas = rho(r, rmin, rmax, rv, rs, p1, p2, irho_pressure)
               ELSE
                  win_ejected_gas = win_norm(k, rmin, rmax, rv, rs, p1, p2, irho_pressure)
                  win_ejected_gas = win_ejected_gas*normalisation(rmin, rmax, rv, rs, p1, p2, irho_pressure)
               END IF

               ! Calculate the value of the density profile prefactor [(Msun/h)/(Mpc/h)^3] and change units from cosmological to SI
               rho0 = m*frac/normalisation(rmin, rmax, rv, rs, p1, p2, irho_density) ! rho0 in [(Msun/h)/(Mpc/h)^3]
               rho0 = rho0*msun/Mpc/Mpc/Mpc ! Overflow with REAL(4) if you use Mpc**3, this converts to SI units [h^2 kg/m^3]
               rho0 = rho0*cosm%h**2 ! Absorb factors of h, so now [kg/m^3]

               ! This is the total thermal pressure of the WHIM
               T0 = HMx_Twhim(hmod, cosm) ! [K]

               ! Factors to convert from Temp x density -> electron pressure (Temp x n; n is all particle number density)
               win_ejected_gas = win_ejected_gas*(rho0/(mp*cosm%mup))*(kb*T0) ! Multiply by number density times temperature k_B [J/m^3]
               win_ejected_gas = win_ejected_gas/(eV*(0.01)**(-3)) ! Change units to pressure in [eV/cm^3]
               win_ejected_gas = win_ejected_gas*cosm%mup/cosm%mue ! Convert from total thermal pressure to electron pressure

            END IF

         ELSE

            STOP 'WIN_FREE_GAS: Error, itype not specified correctly'

         END IF

      END IF

   END FUNCTION win_ejected_gas

   REAL FUNCTION win_stars(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for the stellar density
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: win_central, win_satellite

      win_central = win_central_stars(real_space, k, m, rv, rs, hmod, cosm)
      win_satellite = win_satellite_stars(real_space, k, m, rv, rs, hmod, cosm)

      win_stars = win_central+win_satellite

   END FUNCTION win_stars

   REAL FUNCTION win_central_stars(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for central stellar density
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: rstar, r, rmin, rmax, p1, p2, frac, rss

      frac = halo_central_star_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_central_stars = 0.

      ELSE

         ! Initially set p1, p2
         p1 = 0.
         p2 = 0.

         ! Set default maximuim and minimum radii
         rmin = 0.
         rmax = rv

         IF (frac < frac_min_delta) THEN
            ! Treat as delta functions if there is not much abundance
            irho = 0
         ELSE IF (hmod%halo_central_stars == 1) THEN
            ! Fedeli (2014)
            irho = 7
            rstar = rv/HMx_cstar(m, hmod, cosm)
            p1 = rstar
            rmax = rv ! Set so that not too much bigger than rstar, otherwise bumps integration goes mad
         ELSE IF (hmod%halo_central_stars == 2) THEN
            ! Schneider & Teyssier (2015), following Mohammed (2014)
            irho = 9
            rstar = rv/HMx_cstar(m, hmod, cosm)
            p1 = rstar
            rmax = min(10.*rstar, rv) ! Set so that not too much bigger than rstar, otherwise bumps integration goes crazy
         ELSE IF (hmod%halo_central_stars == 3) THEN
            ! Delta function
            irho = 0
         ELSE IF (hmod%halo_central_stars == 4) THEN
            ! Transition mass between NFW and delta function
            ! TODO: mstar here is the same as in the stellar halo-mass fraction. It should probably not be this
            IF (m < HMx_Mstar(hmod, cosm)) THEN
               irho = 0 ! Delta function
            ELSE
               irho = 5 ! NFW
            END IF
         ELSE
            STOP 'WIN_CENTRAL_STARS: Error, halo_central_stars specified incorrectly'
         END IF

         ! So that the stars see the original halo
         IF (hmod%fix_star_concentration) THEN
            rss = rs*hydro_concentration_modification(m, hmod, cosm)
         ELSE
            rss = rs
         END IF

         IF (real_space) THEN
            r = k
            win_central_stars = rho(r, rmin, rmax, rv, rss, p1, p2, irho)
            win_central_stars = win_central_stars/normalisation(rmin, rmax, rv, rss, p1, p2, irho)
         ELSE
            ! Properly normalise and convert to overdensity
            win_central_stars = m*win_norm(k, rmin, rmax, rv, rss, p1, p2, irho)/comoving_matter_density(cosm)
         END IF

         win_central_stars = frac*win_central_stars

      END IF

   END FUNCTION win_central_stars

   REAL FUNCTION win_satellite_stars(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for the satellite stellar density
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: rstar, r, rmin, rmax, p1, p2, frac, rss

      frac = halo_satellite_star_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_satellite_stars = 0.

      ELSE

         ! Initially set p1, p2
         p1 = 0.
         p2 = 0.

         ! Set default maximuim and minimum radii
         rmin = 0.
         rmax = rv

         IF (frac < frac_min_delta) THEN
            ! Treat as delta functions if there is not much abundance
            irho = 0
         ELSE IF (hmod%halo_satellite_stars == 1) THEN
            ! NFW
            irho = 5
         ELSE IF (hmod%halo_satellite_stars == 2) THEN
            ! Fedeli (2014)
            irho = 7
            rstar = rv/HMx_cstar(m, hmod, cosm)
            p1 = rstar
            rmax = rv ! Set so that not too much bigger than rstar, otherwise bumps integration misbehaves
         ELSE
            STOP 'WIN_SATELLITE_STARS: Error, halo_satellite_stars specified incorrectly'
         END IF

         ! So that the stars see the original halo
         IF (hmod%fix_star_concentration) THEN
            rss = rs*hydro_concentration_modification(m, hmod, cosm)
         ELSE
            rss = rs
         END IF

         IF (real_space) THEN
            r = k
            win_satellite_stars = rho(r, rmin, rmax, rv, rss, p1, p2, irho)
            win_satellite_stars = win_satellite_stars/normalisation(rmin, rmax, rv, rss, p1, p2, irho)
         ELSE
            ! Properly normalise and convert to overdensity
            win_satellite_stars = m*win_norm(k, rmin, rmax, rv, rss, p1, p2, irho)/comoving_matter_density(cosm)
         END IF

         win_satellite_stars = frac*win_satellite_stars

      END IF

   END FUNCTION win_satellite_stars

   REAL FUNCTION win_neutrino(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for the collapsed halo neutrino density
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r, frac, rmin, rmax, p1, p2
      INTEGER :: irho

      frac = halo_neutrino_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_neutrino = 0.

      ELSE

         ! Set some default values
         p1 = 0.
         p2 = 0.
         rmin = 0.
         rmax = rv

         ! Smooth profile (rho=0)
         IF (hmod%halo_neutrino == 1) THEN
            irho = 19 ! Smooth
         ELSE IF (hmod%halo_neutrino == 2) THEN
            irho = 5 ! NFW 
         ELSE  
            STOP 'WIN_NEUTRINO: Error, neutrino profile not specified correctly'
         END IF

         IF (real_space) THEN
            r = k
            win_neutrino = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
            win_neutrino = win_neutrino/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
         ELSE
            ! Properly normalise and convert to overdensity
            win_neutrino = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)/comoving_matter_density(cosm)
         END IF

         win_neutrino = frac*win_neutrino

      END IF

   END FUNCTION win_neutrino

   REAL FUNCTION win_electron_pressure(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo electron pressure profile function for the sum of bound + unbound electron gas
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (hmod%electron_pressure == 1) THEN
         ! This overrides everything and just uses the UPP
         win_electron_pressure = UPP(real_space, k, m, rv, rs, hmod, cosm)
      ELSE IF (hmod%electron_pressure == 2) THEN
         ! Otherwise use...
         win_electron_pressure = win_gas(real_space, field_electron_pressure, k, m, rv, rs, hmod, cosm)
      ELSE
         STOP 'WIN_ELECTRON_PRESSURE: Error, electron_pressure specified incorrectly'
      END IF

   END FUNCTION win_electron_pressure

   REAL FUNCTION win_void(real_space, k, m, rv, rs, hmod, cosm)

      ! Void profile
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2

      ! Initially set p1, p2
      p1 = 0.
      p2 = 0.

      ! Set default maximuim and minimum radii
      rmin = 0.
      rmax = rv

      IF (hmod%halo_void == 1) THEN
         ! Top-hat void
         irho = 2
         rmin = 0.
         rmax = 10.*rv
      ELSE
         STOP 'WIN_VOID: Error, halo_void specified incorrectly'
      END IF

      IF (real_space) THEN
         r = k
         win_void = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
         win_void = win_void/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
      ELSE
         win_void = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)/comoving_matter_density(cosm)
      END IF

   END FUNCTION win_void

   REAL FUNCTION win_compensated_void(real_space, k, m, rv, rs, hmod, cosm)

      ! Profile for compensated voids
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2

      ! Initially set p1, p2
      p1 = 0.
      p2 = 0.

      ! Set default maximuim and minimum radii
      rmin = 0.
      rmax = rv

      IF (hmod%halo_compensated_void == 1) THEN
         ! Top-hat
         irho = 2
         rmin = 0.
         rmax = 10.*rv
      ELSE
         STOP 'WIN_COMPENSATED_VOID: Error, halo_compensated_void specified incorrectly'
      END IF

      IF (real_space) THEN
         r = k
         win_compensated_void = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
         win_compensated_void = win_compensated_void/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
      ELSE
         win_compensated_void = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)/comoving_matter_density(cosm)
      END IF

   END FUNCTION win_compensated_void

   REAL FUNCTION win_haloes(real_space, mmin, mmax, k, m, rv, rs, hmod, cosm)

      ! Halo profile function for haloes
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: mmin
      REAL, INTENT(IN) :: mmax
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2, N, nhalo
      REAL :: nu1, nu2

      IF (m < mmin .OR. m > mmax) THEN
         N = 0.
      ELSE
         N = 1.
      END IF

      IF (N == 0.) THEN

         win_haloes = 0.

      ELSE

         ! Default minimum and maximum radii
         rmin = 0.
         rmax = rv

         ! This calculation really only needs to be done once
         ! This could slow down the calculation by a large amount
         ! TODO: Ensure this is only done once; see init_galaxies
         nu1 = nu_M(mmin, hmod, cosm)
         nu2 = nu_M(mmax, hmod, cosm)
         nhalo = mean_halo_number_density(nu1, nu2, hmod, cosm)

         ! Default additional halo parameters
         p1 = 0.
         p2 = 0.

         ! Delta function
         irho = 0

         IF (real_space) THEN
            r = k
            win_haloes = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
            win_haloes = win_haloes/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
         ELSE
            win_haloes = win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)/nhalo
         END IF

      END IF

   END FUNCTION win_haloes

   REAL FUNCTION win_galaxies(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for all galaxies
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: wc, ws

      IF (.NOT. hmod%has_galaxies) CALL init_galaxies(hmod, cosm)

      wc = win_centrals(real_space, k, m, rv, rs, hmod, cosm)*hmod%n_c
      ws = win_satellites(real_space, k, m, rv, rs, hmod, cosm)*hmod%n_s
      win_galaxies = (wc+ws)/hmod%n_g

   END FUNCTION win_galaxies

   REAL FUNCTION win_centrals(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for central galaxies
      ! TODO: Where does the normalisation division by hmod%n_c go?
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2, N

      IF (.NOT. hmod%has_galaxies) CALL init_galaxies(hmod, cosm)

      N = N_centrals(m, hmod)

      IF (N == 0.) THEN

         win_centrals = 0.

      ELSE

         ! Default minimum and maximum radii
         rmin = 0.
         rmax = rv

         ! Default additional halo parameters
         p1 = 0.
         p2 = 0.

         ! Delta functions
         irho = 0

         IF (real_space) THEN
            r = k
            win_centrals = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
            win_centrals = win_centrals/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
         ELSE
            win_centrals = win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)!/hmod%n_c
         END IF
         win_centrals = N*win_centrals/hmod%n_c

      END IF

   END FUNCTION win_centrals

   REAL FUNCTION win_satellites(real_space, k, m, rv, rs, hmod, cosm)

      ! Halo profile for satellite galaxies
      ! TODO: Where does the normalisation division by hmod%n_s go?
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2, N

      IF (.NOT. hmod%has_galaxies) CALL init_galaxies(hmod, cosm)

      N = N_satellites(m, hmod)

      IF (N == 0.) THEN

         win_satellites = 0.

      ELSE

         ! Initially set p1, p2
         p1 = 0.
         p2 = 0.

         ! Default maximum and minimum radii
         rmin = 0.
         rmax = rv
         
         IF (hmod%halo_satellites == 1) THEN
            irho = 5 ! 5 - NFW profile
         ELSE IF (hmod%halo_satellites == 2) THEN
            irho = 1 ! 1 - Isothermal profile
         ELSE
            STOP 'WIN_SATELLITES: Error, halo profile not set correctly'
         END IF

         IF (real_space) THEN
            r = k
            win_satellites = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
            win_satellites = win_satellites/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
         ELSE
            win_satellites = win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)!/hmod%n_s
         END IF
         win_satellites = N*win_satellites/hmod%n_s

      END IF

   END FUNCTION win_satellites

   REAL FUNCTION win_CIB(real_space, nu, k, m, rv, rs, hmod, cosm)

      ! Halo profile for all matter under the assumption that it is all CDM
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: nu
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2, z
      REAL :: crap

      REAL, PARAMETER :: a = 1. ! Dust blob size relative to halo virial radius
      REAL, PARAMETER :: T = 15. ! Dust temperature [K]
      REAL, PARAMETER :: beta = 1.6 ! Grey-body power-law index

      ! Set additional halo parameters to zero
      p1 = 0.
      p2 = 0.

      rmin = 0.
      rmax = rv

      ! Prevent compile warnings
      crap = m
      crap = cosm%Om_m

      ! Halo type
      ! 0 - Delta function
      ! 1 - Isothermal
      ! 5 - NFW
      irho = 1

      IF (real_space) THEN
         r = k
         win_CIB = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
         win_CIB = win_CIB/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
      ELSE
         !Properly normalise and convert to overdensity
         win_CIB = win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)
      END IF

      z = hmod%z
      win_CIB = grey_body_nu((1.+z)*nu, T, beta) ! Get the black-body radiance [W m^-2 Sr^-1 Hz^-1]
      win_CIB = win_CIB/Jansky     ! Convert units to Jansky [Jy Sr^-1]
      win_CIB = win_CIB*(a*rv)**2  ! [Jy Sr^-1 (Mpc/h)^2]
      !win_CIB=win_CIB/(1e-3+luminosity_distance(a,cosm))**2 ! Bad idea because divide by zero when z=0

   END FUNCTION win_CIB

   REAL FUNCTION grey_body_nu(nu, T, beta)

      ! Grey body irradiance [W m^-2 Hz^-1 Sr^-1]
      USE physics
      REAL, INTENT(IN) :: nu          ! Observing frequency [Hz]
      REAL, INTENT(IN) :: T           ! Grey body temperature [K]
      REAL, INTENT(IN) :: beta        ! Power-law index
      REAL, PARAMETER :: tau = 1.     ! Emissivity (degenerate with R in CIB work)
      REAL, PARAMETER :: nu0 = 545.e9 ! Reference frequency [Hz]

      grey_body_nu = tau*((nu/nu0)**beta)*black_body_nu(nu, T)

   END FUNCTION grey_body_nu

   REAL FUNCTION N_centrals(m, hmod)

      ! The number of central galaxies as a function of halo mass
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod

      IF (m > hmod%mhalo_min .AND. m < hmod%mhalo_max) THEN
         N_centrals = 1
      ELSE
         N_centrals = 0
      END IF

   END FUNCTION N_centrals

   REAL FUNCTION N_satellites(m, hmod)

      ! The number of satellite galxies as a function of halo mass
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod

      IF (m < hmod%mhalo_min .OR. m > hmod%mhalo_max) THEN
         N_satellites = 0
      ELSE
         N_satellites = floor(m/hmod%mhalo_min)-1
      END IF

   END FUNCTION N_satellites

   REAL FUNCTION N_galaxies(m, hmod)

      ! The number of central galaxies as a function of halo mass
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod

      N_galaxies = N_centrals(m, hmod)+N_satellites(m, hmod)

   END FUNCTION N_galaxies

   REAL FUNCTION win_HI(real_space, k, m, rv, rs, hmod, cosm)

      ! Returns the real or Fourier space HI halo profile
      LOGICAL, INTENT(IN) :: real_space ! Real space or Fourier space
      REAL, INTENT(IN) :: k  ! Comoving wave vector (or radius)
      REAL, INTENT(IN) :: m  ! Halo mass
      REAL, INTENT(IN) :: rv ! Halo virial radius
      REAL, INTENT(IN) :: rs ! Halo scale radius
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halomodel
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      INTEGER :: irho
      REAL :: r, rmin, rmax, p1, p2, frac
      REAL :: r0, alpha, c_HI, r_HI, z

      IF (hmod%has_HI .EQV. .FALSE.) CALL init_HI(hmod, cosm)

      frac = halo_HI_fraction(m, hmod, cosm)

      IF (frac == 0.) THEN

         win_HI = 0.

      ELSE

         ! Default minimum and maximum radii
         rmin = 0.
         rmax = rv

         ! Default additional halo parameters
         p1 = 0.
         p2 = 0.

         IF (hmod%halo_HI == 1) THEN
            ! NFW profile
            irho = 5
         ELSE IF (hmod%halo_HI == 2) THEN
            ! Delta function
            irho = 0
         ELSE IF (hmod%halo_HI == 3) THEN
            ! Polynomial with exponential cut off Villaescusa-Navarro et al. (1804.09180)
            irho = 25
            r0 = 10**(-2.5) ! 0.003 Mpc/h (really small)
            alpha = 3.00 ! Tending to homogeneity
            p1 = r0
            p2 = alpha
         ELSE IF (hmod%halo_HI == 4) THEN
            ! Modified NFW with exponential cut off Villaescusa-Navarro et al. (1804.09180)
            irho = 26
            r0 = 10**(-2.5) ! 0.003 Mpc/h (really small)
            r_HI = 10**(-3.0)
            p1 = r0
            p2 = r_HI
         ELSE IF (hmod%halo_HI == 5) THEN
            ! Modified NFW from Padmanabhan & Refreiger (1607.01021)
            irho = 27
            z = hmod%z
            c_HI = 130.
            c_HI = 4.*c_HI*((M/1e11)**(-0.109))/(1.+z)
            r_HI = rv/c_HI
            p1 = r_HI
         ELSE
            STOP 'win_HI: Error, halo_HI not specified correctly'
         END IF

         IF (real_space) THEN
            r = k
            win_HI = rho(r, rmin, rmax, rv, rs, p1, p2, irho)
            win_HI = win_HI/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
         ELSE
            !win_HI=win_norm(k,rmin,rmax,rv,rs,p1,p2,irho) ! Wrong, but is what I first sent Richard and Kiyo
            win_HI = m*win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)/hmod%rho_HI
         END IF

         win_HI = frac*win_HI

      END IF

   END FUNCTION win_HI

   REAL FUNCTION virial_temperature(M, rv, a, cosm)

      ! Halo physical virial temperature [K]
      ! Calculates the temperature as if pristine gas falls into the halo
      ! Energy is equally distributed between the particles
      REAL, INTENT(IN) :: M  ! virial mass
      REAL, INTENT(IN) :: rv ! virial radius
      REAL, INTENT(IN) :: a  ! scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm ! cosmology

      REAL, PARAMETER :: modes = 3. ! 1/2 k_BT per mode, 3 modes for 3 dimensions

      virial_temperature = bigG*((M*msun)*mp*cosm%mup)/(a*rv*mpc) ! NEW: a in denominator converts comoving->physical halo radius
      virial_temperature = virial_temperature/(kb*modes/2.) ! Convert to temperature from energy

   END FUNCTION virial_temperature

   REAL FUNCTION UPP(real_space, k, m, rv, rs, hmod, cosm)

      ! Universal electron pressure profile (Arnaud et al. 2010; arxiv:0910.1234)
      ! Note *very* well that this is for *electron* pressure (see Arnaud 2010 and my notes on this)
      ! Note that is is also the physical pressure, and relates to the comoving pressure via (1+z)^3
      ! The units of the pressure profile are [eV/cm^3]
      LOGICAL, INTENT(IN) :: real_space
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r500c, rmin, rmax, a, z, r, m500c, E
      REAL, PARAMETER :: alphap = 0.12 ! Exponent correction
      REAL, PARAMETER :: b = 0.        ! Hydrostatic mass bias
      INTEGER, PARAMETER :: irho = 14  ! Set UPP profile
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      IF (hmod%has_mass_conversions .EQV. .FALSE.) CALL convert_mass_definitions(hmod, cosm)

      ! Get r500 for UPP
      r500c = exp(find(log(m), hmod%log_m, log(hmod%r500c), hmod%n, iorder, ifind, imeth)) ! [Mpc/h]

      ! Set the radius range for the profile
      rmin = 0.
      rmax = rv

      IF (real_space) THEN
         r = k
         UPP = rho(r, rmin, rmax, rv, rs, r500c, zero, irho)
      ELSE
         UPP = winint(k, rmin, rmax, rv, rs, r500c, zero, irho, imeth_win)
      END IF

      ! UPP, P(x), equation 4.1 in Ma et al. (2015)
      m500c = exp(find(log(m), hmod%log_m, log(hmod%m500c), hmod%n, iorder, ifind, imeth)) ![Msun/h]
      m500c = m500c*(1.-b) ![Msun/h]

      ! Dimensionless Hubble parameter
      a = hmod%a ! Scale-factor appears explicitly here
      E = sqrt(Hubble2(a, cosm))

      ! Pre-factors from equation 4.1 in Ma et al. (2015) [eV cm^-3 - no h factors]
      UPP = UPP*((m500c/2.1e14)**(alphap+2./3.))*(E**(8./3.))*1.65*(cosm%h/0.7)**2

      ! The standard UPP is written is in physical units
      ! scale to comoving using (1+z)^3 because pressure ~energy density
      z = hmod%z ! Redshift appears explicitly here
      UPP = UPP/(1.+z)**3

   END FUNCTION UPP

   REAL FUNCTION rho(r, rmin, rmax, rv, rs, p1, p2, irho)

      ! This is an UNNORMALISED halo profile of any sort

      ! Types of profile
      ! ================
      !  0 - Delta function at r=0
      !  1 - Isothermal: r^-2
      !  2 - Top hat: rho ~ constant
      !  3 - Moore (1999)
      !  4 - NFW (1997)
      !  5 - Analytic NFW
      !  6 - Beta model with beta=2/3
      !  7 - Star profile
      !  8 - Komatsu & Seljak (2001) according to Schneider & Teyssier (2015)
      !  9 - Stellar profile from Schneider & Teyssier (2015)
      ! 10 - Ejected gas profile (Schneider & Teyssier 2015)
      ! 11 - Simplified Komatsu & Seljak (2001) density
      ! 12 - Simplified Komatsu & Seljak (2001) temperature
      ! 13 - Simplified Komatsu & Seljak (2001) pressure
      ! 14 - Universal pressure profile
      ! 15 - Isothermal beta model, beta=0.86 (Ma et al. 2015)
      ! 16 - Isothermal exterior
      ! 17 - Power-law profile
      ! 18 - Cubic profile: r^-3
      ! 19 - Smooth profile (rho = 0, not really physical)
      ! 20 - Exponential profile
      ! 21 - Full Komatsu & Seljak (2001) density
      ! 22 - Full Komatsu & Seljak (2001) temperature
      ! 23 - Full Komatsu & Seljak (2001) pressure
      ! 24 - Cored NFW profile (Copeland, Taylor & Hall 2018)
      ! 25 - Polynomial with central hole (Villaescusa-Navarro et al. 2018)
      ! 26 - Modified NFW with central hole (Villaescusa-Navarro et al. 2018)
      ! 27 - Modified NFW (Padmanabhan & Refregier 2018)
      ! 28 - Shell

      REAL, INTENT(IN) :: r, rmin, rmax, rv, rs, p1, p2 ! Standard profile parameters
      INTEGER, INTENT(IN) :: irho
      REAL :: y, ct, t, c, beta, Gamma, r500c, rt, A, re, rstar, B, rb, r0, alpha, rh, eta0
      REAL :: f1, f2
      REAL :: crap

      ! UPP parameters
      REAL, PARAMETER :: P0 = 6.41
      REAL, PARAMETER :: c500 = 1.81
      REAL, PARAMETER :: alpha_UPP = 1.33
      REAL, PARAMETER :: beta_UPP = 4.13
      REAL, PARAMETER :: gamma_UPP = 0.31

      ! To stop compile-time warnings
      crap = p2

      IF (r < rmin .OR. r > rmax) THEN
         ! The profile is considered to be zero outside this region
         rho = 0.
      ELSE
         IF (irho == 0) THEN
            ! Delta function
            ! Do not assign any value to rho as this gets handled properly elsewhere
            rho = 0.
         ELSE IF (irho == 1 .OR. irho == 16) THEN
            ! Isothermal
            rho = 1./r**2
         ELSE IF (irho == 2) THEN
            ! Top hat
            rho = 1.
         ELSE IF (irho == 3) THEN
            ! Moore (1999)
            y = r/rs
            rho = 1./((y**1.5)*(1.+y**1.5))
         ELSE IF (irho == 4 .OR. irho == 5) THEN
            ! NFW (1997)
            y = r/rs
            rho = 1./(y*(1.+y)**2)
         ELSE IF (irho == 6) THEN
            ! Isothermal beta model (X-ray gas; SZ profiles; beta=2/3 fixed)
            ! Also known as 'cored isothermal profile'
            ! In general this is (1+(r/rs)**2)**(-3*beta/2); beta=2/3 cancels the last power
            rho = 1./(1.+(r/rs)**2)
         ELSE IF (irho == 7) THEN
            ! Stellar profile from Fedeli (2014a)
            rstar = p1
            y = r/rstar
            rho = (1./y)*exp(-y)
         ELSE IF (irho == 8) THEN
            ! Komatsu & Seljak (2001) profile with NFW transition radius
            ! VERY slow to calculate the W(k) for some reason
            ! Also creates a weird upturn in P(k) that I do not think can be correct
            STOP 'RHO: This profile is very difficult for some reason'
            t = sqrt(5.)
            rt = rv/t
            y = r/rs
            c = rs/rv
            ct = c/t
            Gamma = (1.+3.*ct)*log(1.+ct)/((1.+ct)*log(1.+ct)-ct)
            IF (r <= rt) THEN
               ! Komatsu Seljak in the interior
               rho = (log(1.+y)/y)**Gamma
            ELSE
               ! NFW in the outskirts
               A = ((rt/rs)*(1.+rt/rs)**2)*(log(1.+rt/rs)/(rt/rs))**Gamma
               rho = A/(y*(1.+y)**2)
            END IF
         ELSE IF (irho == 9) THEN
            ! Stellar profile from Schneider & Teyssier (2015) via Mohammed (2014)
            rstar = p1
            rho = exp(-(r/(2.*rstar))**2)/r**2
            ! Converting to y caused the integration to crash for some reason !?!
            !y=r/rs
            !rho=exp(-(y/2.)**2.)/y**2.
         ELSE IF (irho == 10) THEN
            ! Ejected gas profile from Schneider & Teyssier (2015)
            re = p1
            rho = exp(-0.5*(r/re)**2)
         ELSE IF (irho == 11 .OR. irho == 12 .OR. irho == 13 .OR. irho == 21 .OR. irho == 22 .OR. irho == 23) THEN
            ! Komatsu & Seljak (2001) profiles for density, temperature and pressure
            IF (irho == 11 .OR. irho == 12 .OR. irho == 13) THEN
               Gamma = p1
               y = r/rs
               rho = log(1.+y)/y
            ELSE IF (irho == 21 .OR. irho == 22 .OR. irho == 23) THEN
               c = rv/rs
               Gamma = p1+0.01*(c-6.5)
               eta0 = 0.00676*(c-6.5)**2+0.206*(c-6.5)+2.48
               f1 = (3./eta0)*(Gamma-1.)/Gamma
               f2 = c/NFW_factor(c)
               B = f1*f2
               IF (B > 1.) B = 1.
               y = r/rs
               rho = 1.-B*(1.-log(1.+y)/y)
            ELSE
               STOP 'RHO: Error, irho specified incorrectly'
            END IF
            IF (irho == 11 .OR. irho == 21) THEN
               ! KS density profile
               rho = rho**(1./(Gamma-1.))
            ELSE IF (irho == 12 .OR. irho == 22) THEN
               ! KS temperature profile (no Gamma dependence)
               rho = rho
            ELSE IF (irho == 13 .OR. irho == 23) THEN
               ! KS pressure profile
               rho = rho**(Gamma/(Gamma-1.))
            END IF
         ELSE IF (irho == 14) THEN
            ! UPP is in terms of r500c, not rv
            r500c = p1
            ! UPP funny-P(x), equation 4.2 in Ma et al. (2015)
            f1 = (c500*r/r500c)**gamma_UPP
            f2 = (1.+(c500*r/r500c)**alpha_UPP)**((beta_UPP-gamma_UPP)/alpha_UPP)
            rho = P0/(f1*f2)
         ELSE IF (irho == 15) THEN
            ! Isothermal beta model with general beta
            !beta=0.86 in Ma et al. (2015)
            beta = p1
            !WRITE(*,*) 'Beta:', beta
            rho = (1.+(r/rs)**2)**(-3.*beta/2.)
         ELSE IF (irho == 16) THEN
            ! Isothermal exterior
            rho = 1./r**2
         ELSE IF (irho == 17) THEN
            ! Power-law profile
            beta = p1
            rho = r**beta
         ELSE IF (irho == 18) THEN
            ! Cubic profile
            rho = r**(-3)
         ELSE IF (irho == 19) THEN
            ! Smooth profile
            rho = 0.
         ELSE IF (irho == 20) THEN
            ! Exponential profile (HI from Padmanabhan et al. 2017)
            re = p1
            rho = exp(-r/re)
         ELSE IF (irho == 24) THEN
            ! Cored NFW (Copeland, Taylor & Hall 2018)
            rb = p1
            f1 = (r+rb)/rs
            f2 = (1.+r/rs)**2
            rho = 1./(f1*f2)
         ELSE IF (irho == 25) THEN
            ! polynomial with central exponential hole
            r0 = p1
            alpha = p2
            rho = (r**(-alpha))*exp(-r0/r)
         ELSE IF (irho == 26) THEN
            ! modified NFW with central exponential hole
            r0 = p1
            rh = p2
            rho = (1./((0.75+r/rh)*(1.+r/rh)**2))*exp(-r0/r)
         ELSE IF (irho == 27) THEN
            ! modified NFW from Padmanabhan & Refregier (2017; 1607.01021)
            rh = p1
            rho = (1./((0.75+r/rh)*(1.+r/rh)**2))
         ELSE IF (irho == 28) THEN
            ! Shell
            rho = 0.
         ELSE
            STOP 'RHO: Error, irho not specified correctly'
         END IF

      END IF

   END FUNCTION rho

   REAL FUNCTION rhor2at0(irho)

      ! This is the value of rho(r)*r^2 at r=0
      ! For most profiles this is zero, BUT not if rho(r->0) -> r^-2
      ! Note if rho(r->0) -> r^n with n<-2 then the profile mass would diverge!
      INTEGER, INTENT(IN) :: irho

      IF (irho == 0) THEN
         STOP 'RHOR2AT0: You should not be here for a delta-function profile'
      ELSE IF (irho == 1 .OR. irho == 9) THEN
         !1 - Isothermal
         !9 - Stellar profile from Schneider & Teyssier (2015)
         rhor2at0 = 1.
      ELSE IF (irho == 18) THEN
         STOP 'RHOR2AT0: Error, profile diverges at the origin'
      ELSE
         rhor2at0 = 0.
      END IF

   END FUNCTION rhor2at0

   REAL FUNCTION normalisation(rmin, rmax, rv, rs, p1, p2, irho)

      ! This calculates the normalisation of a halo
      ! This is the integral of 4pir^2*rho(r)*dr between rmin and rmax
      ! This is the total profile 'mass' if rmin=0 and rmax=rv

      ! Profile results
      !  0 - Delta function (M = 1)
      !  1 - Isothermal (M = 4pi*rv)
      !  2 - Top hat (M = (4pi/3)*rv^3)
      !  3 - Moore (M = (8pi/3)*rv^3*ln(1+c^1.5)/c^3)
      !  4 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
      !  5 - NFW (M = 4pi*rs^3*[ln(1+c)-c/(1+c)])
      !  6 - Beta model with beta=2/3 (M = 4pi*rs^3*(rv/rs-atan(rv/rs)))
      !  7 - Fedeli stellar model (M = 4pi*rstar^2 * [1-exp(-rmax/rstar)*(1.+rmax/rstar)]
      !  8 - No
      !  9 - Stellar profile (Schneider & Teyssier 2015; M = 4(pi^1.5)*rstar; assumed rmax ~ infinity)
      ! 10 - Ejected gas profile (Schneider & Teyssier 2015; M = 4pi*sqrt(pi/2)*re^3; assumed rmax ~ infinity)
      ! 11 - No
      ! 12 - No
      ! 13 - No
      ! 14 - No
      ! 15 - No
      ! 16 - Isothermal shell (M = 4pi*(rmax-rmin))
      ! 17 - No
      ! 18 - Cubic profile (M = 4pi*log(rmax/rmin))
      ! 19 - Smooth profile (M = 1; prevents problems)
      ! 20 - No
      ! 21 - No
      ! 22 - No
      ! 23 - No
      ! 24 - No
      ! 25 - No
      ! 26 - No
      ! 27 - No
      ! 28 - Shell (M = 4pi*rv^3)

      REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
      INTEGER, INTENT(IN) :: irho
      REAL :: cmax, re, rstar, beta, rb, c, b

      IF (irho == 0) THEN
         ! Delta function
         normalisation = 1.
      ELSE IF (irho == 1 .OR. irho == 16) THEN
         ! Isothermal
         normalisation = 4.*pi*(rmax-rmin)
      ELSE IF (irho == 2) THEN
         ! Top hat
         normalisation = 4.*pi*(rmax**3-rmin**3)/3.
      ELSE IF (irho == 3 .AND. (rmin /= 0)) THEN
         ! Moore et al. (1999)
         cmax = rmax/rs
         normalisation = (2./3.)*4.*pi*(rs**3)*log(1.+cmax**1.5)
      ELSE IF (irho == 4 .OR. irho == 5 .AND. (rmin /= 0)) THEN
         ! NFW (1997)
         cmax = rmax/rs
         normalisation = norm_NFW(rs, cmax)
      ELSE IF (irho == 6 .AND. (rmin /= 0)) THEN
         ! Isothermal beta model with beta=2/3
         cmax = rmax/rs
         normalisation = 4.*pi*(rs**3)*(cmax-atan(cmax))
      ELSE IF (irho == 7 .AND. (rmin /= 0)) THEN
         ! Fedeli (2014) stellar mode
         ! This would be even easier if rmax -> infinity (just 4*pi*rstar^2)
         rstar = p1
         normalisation = 4.*pi*(rstar**3)*(1.-exp(-rmax/rstar)*(1.+rmax/rstar))
         !normalisation=4.*pi*rstar**3 ! rmax/rstar -> infinity limit (rmax >> rstar)
      ELSE IF (irho == 9 .AND. (rmin /= 0)) THEN
         ! Stellar profile from Schneider & Teyssier (2015)
         ! CARE: Assumed to go on to r -> infinity
         rstar = p1
         normalisation = 4.*(pi**(3./2.))*rstar
      ELSE IF (irho == 10 .AND. (rmin /= 0)) THEN
         ! Ejected gas profile from Schneider & Teyssier (2015)
         ! CARE: Assumed to go on to r -> infinity
         re = p1
         normalisation = 4.*pi*sqrt(pi/2.)*re**3
      ELSE IF (irho == 17) THEN
         ! Power-law profile
         beta = p1
         normalisation = (4.*pi/(beta+3.))*(rmax**(beta+3.)-rmin**(beta+3.))
      ELSE IF (irho == 18) THEN
         ! Cubic profile
         normalisation = 4.*pi*log(rmax/rmin)
      ELSE IF (irho == 19) THEN
         ! Smooth profile, needs a normalisation because divided by
         ! CARE: This cannot be set to zero
         normalisation = 1.
      ELSE IF (irho == 24 .AND. (rv == rmax)) THEN
         ! Cored NFW profile
         rb = p1
         IF (rb == 0.) THEN
            ! This is then the standard NFW case
            c = rv/rs
            normalisation = norm_NFW(rs, c)
         ELSE
            ! Otherwise there is actually a core
            b = rv/rb
            c = rv/rs
            normalisation = (4.*pi*rs**3)/(b-c)**2
            normalisation = normalisation*(b*(b-2.*c)*NFW_factor(c)+(log(1.+b)-b/(1.+c))*c**2)
         END IF
      ELSE IF (irho == 28) THEN
         ! Shell
         normalisation = 4.*pi*rv**3
      ELSE
         ! Otherwise need to do the integral numerically
         ! k=0 gives normalisation
         normalisation = winint(zero, rmin, rmax, rv, rs, p1, p2, irho, imeth_win)
      END IF

   END FUNCTION normalisation

   REAL FUNCTION norm_NFW(rs, c)

      REAL, INTENT(IN) :: rs
      REAL, INTENT(IN) :: c

      norm_NFW = 4.*pi*(rs**3)*NFW_factor(c)

   END FUNCTION norm_NFW

   REAL FUNCTION win_norm(k, rmin, rmax, rv, rs, p1, p2, irho)

      ! Normalised halo window functions

      ! Types of profile
      ! ================
      !  0 - Delta function at r=0
      !  1 - Isothermal: r^-2
      !  2 - Top hat: constant
      !  3 - No
      !  4 - No
      !  5 - NFW
      !  6 - No
      !  7 - Star profile
      !  8 - No
      !  9 - Stellar profile from Schneider & Teyssier (2015)
      ! 10 - Ejected gas profile (Schneider & Teyssier 2015)
      ! 11 - No
      ! 12 - No
      ! 13 - No
      ! 14 - No
      ! 15 - No
      ! 16 - Isothermal shell
      ! 17 - No
      ! 18 - No
      ! 19 - Smooth profile
      ! 20 - Exponential profile
      ! 21 - No
      ! 22 - No
      ! 23 - No
      ! 24 - Cored NFW profile (Copeland, Taylor & Hall 2018)
      ! 25 - No
      ! 26 - No
      ! 27 - No
      ! 28 - Shell

      ! Calculates the normalised spherical Fourier Transform of the density profile
      ! Note that this means win_norm(k->0)=1
      ! and that win must be between 0 and 1
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: rmin
      REAL, INTENT(IN) :: rmax
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      REAL, INTENT(IN) :: p1
      REAL, INTENT(IN) :: p2
      INTEGER, INTENT(IN) :: irho
      REAL :: re, f1, f2, rstar, kstar, rb

      IF (k == 0.) THEN

         ! If called for the zero mode (e.g. for the normalisation)
         win_norm = 1.

      ELSE

         IF (irho == 0) THEN
            ! Delta function profile is not localised in Fourier Space
            win_norm = 1.
         ELSE IF (irho == 1) THEN
            win_norm = wk_isothermal(k*rmax)
         ELSE IF (irho == 2) THEN
            ! Analytic for top hat
            win_norm = wk_tophat(k*rmax)
         ELSE IF (irho == 5) THEN
            ! Analytic for NFW
            win_norm = win_NFW(k, rmax, rs)
         ELSE IF (irho == 7) THEN
            ! Analytic for Fedeli (2014) stellar profile
            rstar = p1
            kstar = k*rstar
            f1 = kstar-exp(-rmax/rstar)*(sin(k*rmax)+kstar*cos(k*rmax))
            f2 = kstar*(1.+kstar**2)
            win_norm = f1/f2
            !win_norm=1./(1.+kstar**2) !bigRstar -> infinity limit (rmax >> rstar)
         ELSE IF (irho == 9) THEN
            ! Normalisation is only valid if rmin=0 and rmax=inf
            rstar = p1
            win_norm = (sqrt(pi)/2.)*erf(k*rstar)/(k*rstar)
         ELSE IF (irho == 10) THEN
            ! Ejected gas profile
            re = p1
            win_norm = exp(-1.5*(k*re)**2.)
         ELSE IF (irho == 16) THEN
            ! Isothermal shells
            win_norm = wk_isothermal_2(k*rmax, k*rmin)
         ELSE IF (irho == 19) THEN
            ! Smooth profile
            win_norm = 0.
         ELSE IF (irho == 20) THEN
            ! Exponential profile
            re = p1
            win_norm = 1./(1.+(k*re)**2)**2
         ELSE IF (irho == 24) THEN
            ! Cored NFW profile
            rb = p1
            IF (rb == 0.) THEN
               ! In this case there is no core
               win_norm = win_NFW(k, rmax, rs)
            ELSE
               ! Otherwise there is a core
               win_norm = win_cored_NFW(k, rmax, rs, rb)
            END IF
         ELSE IF (irho == 28) THEN
            ! Shell
            win_norm = sinc(k*rv)
         ELSE
            ! Numerical integral over the density profile (slower)
            win_norm = winint(k, rmin, rmax, rv, rs, p1, p2, irho, imeth_win)/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
         END IF

      END IF

   END FUNCTION win_norm

   SUBROUTINE winint_speed_tests(k, nk, rmin, rmax, rv, rs, p1, p2, irho, base, ext)

      INTEGER, INTENT(IN) :: nk, irho
      REAL, INTENT(IN) :: k(nk), rmin, rmax, rv, rs, p1, p2
      CHARACTER(len=*), INTENT(IN) :: base
      CHARACTER(len=*), INTENT(IN) :: ext
      CHARACTER(len=256) :: outfile
      INTEGER :: j, i, ii, imeth, n, ntime
      LOGICAL :: timing
      REAL :: t1, t2, w

      ! j = 1: Time calculation
      ! j = 2: Write calculation out
      DO j = 1, 2

         timing = .FALSE.
         IF (j == 2) THEN
            timing = .TRUE.
            WRITE (*, *) 'WININT_SPEED_TESTS: Doing this many evaluations for timing test:', ntime
            WRITE (*, *)
         END IF

         DO imeth = 1, nmeth_win

            WRITE (*, *) 'WININT_SPEED_TESTS: Method:', imeth
            IF (imeth == 1) WRITE (*, *) 'WININT_SPEED_TESTS: Simple integration'
            IF (imeth == 2) WRITE (*, *) 'WININT_SPEED_TESTS: Simple integration over bumps'
            IF (imeth == 3) WRITE (*, *) 'WININT_SPEED_TESTS: Standard intergation'
            IF (imeth == 4) WRITE (*, *) 'WININT_SPEED_TESTS: Standard integration over bumps'
            IF (imeth == 5) WRITE (*, *) 'WININT_SPEED_TESTS: Constant approximation for bumps'
            IF (imeth == 6) WRITE (*, *) 'WININT_SPEED_TESTS: Linear approximation for bumps'
            IF (imeth == 7) WRITE (*, *) 'WININT_SPEED_TESTS: Quadratic approximation for bumps'
            IF (imeth == 8) WRITE (*, *) 'WININT_SPEED_TESTS: Cubic approximation for bumps'
            IF (imeth == 9) WRITE (*, *) 'WININT_SPEED_TESTS: Hybrid with constant approximation for bumps'
            IF (imeth == 10) WRITE (*, *) 'WININT_SPEED_TESTS: Hybrid with linear approximation for bumps'
            IF (imeth == 11) WRITE (*, *) 'WININT_SPEED_TESTS: Hybrid with quadratic approximation for bumps'
            IF (imeth == 12) WRITE (*, *) 'WININT_SPEED_TESTS: Hybrid with cubic approximation for bumps'
            IF (imeth == 13) WRITE (*, *) 'WININT_SPEED_TESTS: Full hybrid'

            IF (timing .EQV. .FALSE.) THEN
               outfile = number_file(base, imeth, ext)
               WRITE (*, *) 'WININT_SPEED_TESTS: Writing data: ', trim(outfile)
               OPEN (7, file=outfile)
               n = 1
               IF (imeth == 1) CALL cpu_time(t1)
            ELSE
               CALL cpu_time(t1)
               n = ntime
            END IF

            ! Loop over number of iterations
            DO ii = 1, n

               ! Loop over wave number and do integration
               DO i = 1, nk
                  w = winint(k(i), rmin, rmax, rv, rs, p1, p2, irho, imeth)/normalisation(rmin, rmax, rv, rs, p1, p2, irho)
                  IF (.NOT. timing) THEN
                     WRITE (7, *) k(i), w
                  END IF
               END DO

            END DO

            IF (.NOT. timing) THEN
               CLOSE (7)
               IF (imeth == 1) THEN
                  CALL cpu_time(t2)
                  ntime = ceiling(winint_test_seconds/(t2-t1))
               END IF
            ELSE
               CALL cpu_time(t2)
               WRITE (*, fmt='(A30,F10.5)') 'WININT_SPEED_TESTS: Time [s]:', t2-t1
            END IF

            WRITE (*, *)

         END DO

      END DO

   END SUBROUTINE winint_speed_tests

   SUBROUTINE winint_diagnostics(rmin, rmax, rv, rs, p1, p2, irho, outfile)

      ! Write out the winint integrand as a function of k
      REAL, INTENT(IN) :: rmin, rmax, rv, rs, p1, p2
      INTEGER, INTENT(IN) :: irho
      CHARACTER(len=256), INTENT(IN) :: outfile
      INTEGER :: i, j
      REAL :: r, k
      REAL, ALLOCATABLE :: integrand(:)

      REAL, PARAMETER :: kmin = 1e-1
      REAL, PARAMETER :: kmax = 1e2
      INTEGER, PARAMETER :: nr = 256 ! Number of points in r
      INTEGER, PARAMETER :: nk = 16  ! Number of points in k

      WRITE (*, *) 'WININT_DIAGNOSTICS: Doing these'
      WRITE (*, *) 'WININT_DIAGNOSTICS: minimum r [Mpc/h]:', REAL(rmin)
      WRITE (*, *) 'WININT_DIAGNOSTICS: maximum r [Mpc/h]:', REAL(rmax)
      WRITE (*, *) 'WININT_DIAGNOSTICS: virial radius [Mpc/h]:', REAL(rv)
      WRITE (*, *) 'WININT_DIAGNOSTICS: scale radius [Mpc/h]:', REAL(rs)
      WRITE (*, *) 'WININT_DIAGNOSTICS: concentration:', REAL(rv/rs)
      WRITE (*, *) 'WININT_DIAGNOSTICS: halo parameter 1:', p1
      WRITE (*, *) 'WININT_DIAGNOSTICS: halo parameter 2:', p2
      WRITE (*, *) 'WININT_DIAGNOSTICS: profile number:', irho
      WRITE (*, *) 'WININT_DIAGNOSTICS: outfile: ', TRIM(outfile)

      ALLOCATE (integrand(nk))

      OPEN (7, file=outfile)
      DO i = 1, nr
         r = progression(0., rmax, i, nr)
         DO j = 1, nk
            k = exp(progression(log(kmin), log(kmax), j, nk))
            integrand(j) = winint_integrand(r, rmin, rmax, rv, rs, p1, p2, irho)*sinc(r*k)
         END DO
         WRITE (7, *) r/rv, (integrand(j), j=1, nk)
      END DO
      CLOSE (7)

      WRITE (*, *) 'WININT_DIAGNOSTICS: Done'
      WRITE (*, *)

   END SUBROUTINE winint_diagnostics

   REAL FUNCTION winint(k, rmin, rmax, rv, rs, p1, p2, irho, imeth)

      ! Calculates W(k,M)
      REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
      INTEGER, INTENT(IN) :: irho, imeth

      ! Integration method
      ! imeth =  1: Simple integration
      ! imeth =  2: Bumps with simple integration
      ! imeth =  3: Standard integration
      ! imeth =  4: Bumps with standard integration
      ! imeth =  5: Constant approximation for bumps
      ! imeth =  6: Linear approximation for bumps
      ! imeth =  7: Quadratic approximation for bumps
      ! imeth =  8: Cubic approximation for bumps
      ! imeth =  9: Hybrid with constant approximation for bumps
      ! imeth = 10: Hybrid with linear approximation for bumps
      ! imeth = 11: Hybrid with quadratic approximation for bumps
      ! imeth = 12: Hybrid with cubic approximation for bumps

      ! Bump methods go crazy with some star profiles (those that drop too fast)
      ! You need to make sure that the rmax for the integration does not extend too far out

      IF (imeth == 1) THEN
         winint = integrate_window_normal(rmin, rmax, k, rmin, rmax, rv, rs, p1, p2, irho, winint_order, acc_win)
      ELSE IF (imeth == 3) THEN
         winint = integrate_window_store(rmin, rmax, k, rmin, rmax, rv, rs, p1, p2, irho, winint_order, acc_win)
      ELSE
         winint = winint_bumps(k, rmin, rmax, rv, rs, p1, p2, irho, winint_order, acc_win, imeth)
      END IF

   END FUNCTION winint

   REAL FUNCTION integrate_window_normal(a, b, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)

      ! Integration routine using 'normal' method to calculate the normalised halo FT
      REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
      INTEGER, INTENT(IN) :: irho
      INTEGER, INTENT(IN) :: iorder
      REAL, INTENT(IN) :: acc
      REAL, INTENT(IN) :: a, b
      REAL :: sum
      REAL :: r, dr, winold, weight
      INTEGER :: n, i, j

      INTEGER, PARAMETER :: jmin = 5
      INTEGER, PARAMETER :: jmax = 30
      INTEGER, PARAMETER :: ninit = 2

      winold = 0.

      IF (a == b) THEN

         integrate_window_normal = 0.

      ELSE

         ! Integrates to required accuracy!
         DO j = 1, jmax

            ! Increase the number of integration points each go until convergence
            n = ninit*(2**(j-1))

            ! Set the integration sum variable to zero
            sum = 0.

            DO i = 1, n

               ! Get the weights
               IF (iorder == 1) THEN
                  ! Composite trapezium weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.5
                  ELSE
                     weight = 1.
                  END IF
               ELSE IF (iorder == 2) THEN
                  ! Composite extended formula weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.416666666666
                  ELSE IF (i == 2 .OR. i == n-1) THEN
                     weight = 1.083333333333
                  ELSE
                     weight = 1.
                  END IF
               ELSE IF (iorder == 3) THEN
                  ! Composite Simpson weights
                  IF (i == 1 .OR. i == n) THEN
                     weight = 0.375
                  ELSE IF (i == 2 .OR. i == n-1) THEN
                     weight = 1.166666666666
                  ELSE IF (i == 3 .OR. i == n-2) THEN
                     weight = 0.958333333333
                  ELSE
                     weight = 1.
                  END IF
               ELSE
                  STOP 'INTEGRATE_WINDOW_NORMAL: Error, order specified incorrectly'
               END IF

               ! Now get r and do the function evaluations
               r = progression(a, b, i, n)
               sum = sum+weight*winint_integrand(r, rmin, rmax, rv, rs, p1, p2, irho)*sinc(r*k)

            END DO

            ! The dr are all equally spaced
            dr = (b-a)/REAL(n-1)

            integrate_window_normal = REAL(sum)*dr

            IF ((j > jmin) .AND. integrate_window_normal == 0.) THEN
               EXIT
            ELSE IF ((j > jmin) .AND. (ABS(-1.+integrate_window_normal/winold) < acc)) THEN
               EXIT
            ELSE
               winold = integrate_window_normal
            END IF

         END DO

      END IF

   END FUNCTION integrate_window_normal

   REAL FUNCTION integrate_window_store(a, b, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder, irho
      REAL, INTENT(IN) :: a, b
      INTEGER :: i, j, n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = 5
      INTEGER, PARAMETER :: jmax = 30

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_window_store = 0.

      ELSE

         ! Reset the sum variables for the integration
         sum_2n = 0.d0
         sum_n = 0.d0
         sum_old = 0.d0
         sum_new = 0.d0

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/REAL(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = winint_integrand(a, rmin, rmax, rv, rs, p1, p2, irho)*sinc(a*k)
               f2 = winint_integrand(b, rmin, rmax, rv, rs, p1, p2, irho)*sinc(b*k)
               sum_2n = 0.5*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = progression(a, b, i, n)
                  fx = winint_integrand(x, rmin, rmax, rv, rs, p1, p2, irho)*sinc(x*k)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_WINDOW_STORE: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_WINDOW_STORE: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO

         integrate_window_store = sum_new

      END IF

   END FUNCTION integrate_window_store

   REAL FUNCTION winint_bumps(k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc, imeth)

      ! Integration routine to calculate the normalised halo FT
      REAL, INTENT(IN) :: k, rmin, rmax, rv, rs, p1, p2
      INTEGER, INTENT(IN) :: irho
      INTEGER, INTENT(IN) :: iorder, imeth
      REAL, INTENT(IN) :: acc
      REAL :: sum, w
      REAL :: r1, r2
      INTEGER :: i, n

      ! This MUST be set to zero for this routine
      IF (rmin /= 0.) STOP 'WININT_BUMPS: Error, rmin must be zero'

      ! Calculate the number of nodes of sinc(k*rmax) for 0<=r<=rmax
      n = floor(k*rmax/pi)

      ! Set the sum variable to zero
      sum = 0.

      ! Integrate over each chunk between nodes separately
      DO i = 0, n

         ! Set the lower integration limit
         IF (k == 0.) THEN
            ! Special case when k=0 to avoid division by zero
            r1 = 0.
         ELSE
            r1 = i*pi/k
         END IF

         ! Set the upper integration limit
         IF (k == 0. .OR. i == n) THEN
            ! Special case when on last section because end is rmax, not a node!
            r2 = rmax
         ELSE
            r2 = (i+1)*pi/k
         END IF

         ! Now do the integration along a section
         IF (imeth == 2) THEN
            w = integrate_window_normal(r1, r2, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)
         ELSE IF (i == 0 .OR. i == n .OR. k == 0. .OR. imeth == 4) THEN
            w = integrate_window_store(r1, r2, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)
         ELSE IF (imeth == 13) THEN
            w = winint_hybrid(r1, r2, i, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)
         ELSE
            IF (imeth == 5 .OR. (imeth == 9 .AND. n > nlim_bumps)) THEN
               w = winint_approx(r1, r2, i, k, rmin, rmax, rv, rs, p1, p2, irho, iorder=0)
            ELSE IF (imeth == 6 .OR. (imeth == 10 .AND. n > nlim_bumps)) THEN
               w = winint_approx(r1, r2, i, k, rmin, rmax, rv, rs, p1, p2, irho, iorder=1)
            ELSE IF (imeth == 7 .OR. (imeth == 11 .AND. n > nlim_bumps)) THEN
               w = winint_approx(r1, r2, i, k, rmin, rmax, rv, rs, p1, p2, irho, iorder=2)
            ELSE IF (imeth == 8 .OR. (imeth == 12 .AND. n > nlim_bumps)) THEN
               w = winint_approx(r1, r2, i, k, rmin, rmax, rv, rs, p1, p2, irho, iorder=3)
            ELSE
               w = integrate_window_store(r1, r2, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)
            END IF
         END IF

         sum = sum+w

         ! Exit if the contribution to the sum is very tiny, this seems to be necessary to prevent crashes
         IF (winint_exit .AND. ABS(w) < acc*ABS(sum)) EXIT

      END DO

      winint_bumps = sum

   END FUNCTION winint_bumps

   REAL FUNCTION winint_approx(rn, rm, i, k, rmin, rmax, rv, rs, p1, p2, irho, iorder)

      ! Approximate forms for the integral over the sine bump times a polynomial
      USE special_functions
      REAL, INTENT(IN) :: rn, rm
      INTEGER, INTENT(IN) :: i
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: rmin, rmax
      REAL, INTENT(IN) :: rv, rs, p1, p2
      INTEGER, INTENT(IN) :: irho
      INTEGER, INTENT(IN) :: iorder
      REAL :: x0, x1, x2, x3
      REAL :: y0, y1, y2, y3
      REAL :: a0, a1, a2, a3
      REAL :: w

      IF (iorder == 0) THEN
         x0 = (rn+rm)/2. ! Middle
      ELSE IF (iorder == 1) THEN
         x0 = rn ! Beginning
         x1 = rm ! End
      ELSE IF (iorder == 2) THEN
         x0 = rn ! Beginning
         x1 = (rn+rm)/2. ! Middle
         x2 = rm ! End
      ELSE IF (iorder == 3) THEN
         x0 = rn ! Beginning
         x1 = rn+1.*(rm-rn)/3. ! Middle
         x2 = rn+2.*(rm-rn)/3. ! Middle
         x3 = rm ! End
      ELSE
         STOP 'WININT_APPROX: Error, iorder specified incorrectly'
      END IF

      y0 = winint_integrand(x0, rmin, rmax, rv, rs, p1, p2, irho)/x0

      IF (iorder == 0) THEN
         a0 = y0
      ELSE IF (iorder == 1) THEN
         y1 = winint_integrand(x1, rmin, rmax, rv, rs, p1, p2, irho)/x1
         CALL fix_polynomial(a1, a0, [x0, x1], [y0, y1])
      ELSE IF (iorder == 2) THEN
         y1 = winint_integrand(x1, rmin, rmax, rv, rs, p1, p2, irho)/x1
         y2 = winint_integrand(x2, rmin, rmax, rv, rs, p1, p2, irho)/x2
         CALL fix_polynomial(a2, a1, a0, [x0, x1, x2], [y0, y1, y2])
      ELSE IF (iorder == 3) THEN
         y1 = winint_integrand(x1, rmin, rmax, rv, rs, p1, p2, irho)/x1
         y2 = winint_integrand(x2, rmin, rmax, rv, rs, p1, p2, irho)/x2
         y3 = winint_integrand(x3, rmin, rmax, rv, rs, p1, p2, irho)/x3
         CALL fix_polynomial(a3, a2, a1, a0, [x0, x1, x2, x3], [y0, y1, y2, y3])
      ELSE
         STOP 'WININT_APPROX: Error, iorder specified incorrectly'
      END IF

      w = 2.*a0
      IF (iorder == 1 .OR. iorder == 2 .OR. iorder == 3) w = w+(rm+rn)*a1
      IF (iorder == 2 .OR. iorder == 3) w = w+(rm**2+rn**2-4./k**2)*a2
      IF (iorder == 3) w = w+(rm**3+rn**3-6.*(rm+rn)/k**2)*a3

      w = ((-1)**i)*w/k**2

      winint_approx = w

   END FUNCTION winint_approx

   REAL FUNCTION winint_hybrid(rn, rm, i, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)

      ! An attempt to do an automatic combination of approximation and proper integration
      ! It turned out to be very slow
      USE special_functions
      REAL, INTENT(IN) :: rn, rm
      INTEGER, INTENT(IN) :: i
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: rmin, rmax
      REAL, INTENT(IN) :: rv, rs, p1, p2
      INTEGER, INTENT(IN) :: irho
      INTEGER, INTENT(IN) :: iorder
      REAL, INTENT(IN) :: acc
      REAL :: x0, x1, x2, x3
      REAL :: y0, y1, y2, y3
      REAL :: a0, a1, a2, a3
      REAL :: w, rmid, epsa0

      REAL, PARAMETER :: eps = 1e-2

      rmid = (rn+rm)/2.

      x0 = rn
      x1 = rn+1.*(rm-rn)/3.
      x2 = rn+2.*(rm-rn)/3.
      x3 = rm

      y0 = winint_integrand(x0, rmin, rmax, rv, rs, p1, p2, irho)/x0
      y1 = winint_integrand(x1, rmin, rmax, rv, rs, p1, p2, irho)/x1
      y2 = winint_integrand(x2, rmin, rmax, rv, rs, p1, p2, irho)/x2
      y3 = winint_integrand(x3, rmin, rmax, rv, rs, p1, p2, irho)/x3

      CALL fix_polynomial(a3, a2, a1, a0, [x0, x1, x2, x3], [y0, y1, y2, y3])

      epsa0 = eps*a0
      IF (ABS(a3*rmid**3) < epsa0 .AND. ABS(a2*rmid**2) < epsa0 .AND. ABS(a1*rmid) < epsa0) THEN
         w = (rm**3+rn**3-6.*(rm+rn)/k**2)*a3+(rm**2+rn**2-4./k**2)*a2+(rm+rn)*a1+2.*a0
         w = w*(-1)**i
      ELSE
         w = integrate_window_store(rn, rm, k, rmin, rmax, rv, rs, p1, p2, irho, iorder, acc)
      END IF

      winint_hybrid = w

   END FUNCTION winint_hybrid

   REAL FUNCTION winint_integrand(r, rmin, rmax, rv, rs, p1, p2, irho)

      ! The integrand for the W(k) integral
      ! Note that the sinc function is *not* included
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: rmin
      REAL, INTENT(IN) :: rmax
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      REAL, INTENT(IN) :: p1
      REAL, INTENT(IN) :: p2
      INTEGER, INTENT(IN) :: irho

      IF (r == 0.) THEN
         winint_integrand = 4.*pi*rhor2at0(irho)
      ELSE
         winint_integrand = 4.*pi*(r**2)*rho(r, rmin, rmax, rv, rs, p1, p2, irho)
      END IF

   END FUNCTION winint_integrand

   REAL FUNCTION win_NFW(k, rv, rs)

      ! The analytic normalised (W(k=0)=1) Fourier Transform of the NFW profile
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      REAL :: c, ks
      REAL :: p1, p2, p3
      REAL :: rmin, rmax
      INTEGER, PARAMETER :: irho = 4

      c = rv/rs
      ks = k*rv/c

      p1 = cos(ks)*F_NFW(k, rv, c)
      p2 = sin(ks)*G_NFW(k, rv, c)
      p3 = sin(ks*c)/(ks*(1.+c))

      rmin = 0.
      rmax = rv
      !win_NFW = 4.*pi*(rs**3)*(p1+p2-p3)/normalisation(rmin, rmax, rv, rs, zero, zero, irho)
      win_NFW = 4.*pi*(rs**3)*(p1+p2-p3)/norm_NFW(rs, c)

   END FUNCTION win_NFW

   REAL FUNCTION win_cored_NFW(k, rv, rs, rb)

      ! The analytic normalised (W(k=0)=1) Fourier Transform of the cored NFW profile
      ! Appendix A of Copeland, Taylor & Hall (1712.07112)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      REAL, INTENT(IN) :: rb
      REAL :: b, c
      REAL :: rmin, rmax
      INTEGER, PARAMETER :: irho = 24

      b = rv/rb
      c = rv/rs

      rmin = 0.
      rmax = rv

      win_cored_NFW = NFW_factor(c)*win_NFW(k, rv, rs)
      win_cored_NFW = win_cored_NFW+(c/(b-c))*(1./(k*rs))*(G_NFW(k, rv, c)*cos(k*rs)-F_NFW(k, rv, c)*sin(k*rs))
      win_cored_NFW = win_cored_NFW-(c/(b-c))*(1./(k*rs))*(G_NFW(k, rv, b)*cos(k*rb)-F_NFW(k, rv, b)*sin(k*rb))
      win_cored_NFW = 4.*pi*(rs**3)*(b/(b-c))*win_cored_NFW/normalisation(rmin, rmax, rv, rs, rb, zero, irho)

   END FUNCTION win_cored_NFW

   REAL FUNCTION F_NFW(k, rv, c)

      ! Equation (A3; top) from Copeland, Taylor & Hall (1712.07112)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: c
      REAL :: ks

      ks = k*rv/c

      F_NFW = Ci(ks*(1.+c))-Ci(ks)

   END FUNCTION F_NFW

   REAL FUNCTION G_NFW(k, rv, c)

      ! Equation (A3; bottom) from Copeland, Taylor & Hall (1712.07112)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: c
      REAL :: ks

      ks = k*rv/c

      G_NFW = Si(ks*(1.+c))-Si(ks)

   END FUNCTION G_NFW

   REAL FUNCTION b_nu(nu, hmod)

      ! Bias function selection
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod

      IF (.NOT. hmod%has_mass_function) CALL init_mass_function(hmod)

      IF (hmod%imf == 1) THEN
         b_nu = b_PS(nu, hmod)
      ELSE IF (is_in_array(hmod%imf, [2, 6, 12, 13])) THEN
         b_nu = b_ST(nu, hmod)
      ELSE IF (hmod%imf == 3 .OR. hmod%imf == 16) THEN
         b_nu = b_Tinker2010_peakbackground(nu, hmod)
      ELSE IF (hmod%imf == 4 .OR. hmod%imf == 11) THEN
         b_nu = 1.
      ELSE IF (hmod%imf == 5) THEN
         b_nu = b_Jenkins(nu, hmod)
      ELSE IF (hmod%imf == 7 .OR. hmod%imf == 15) THEN
         b_nu = b_Tinker2008(nu, hmod)
      ELSE IF (hmod%imf == 8) THEN
         b_nu = b_Warren(nu, hmod)
      ELSE IF (hmod%imf == 9) THEN
         b_nu = b_Reed(nu, hmod)
      ELSE IF (hmod%imf == 10 .OR. hmod%imf == 14) THEN
         b_nu = b_Bhattacharya(nu, hmod)
      ELSE IF (hmod%imf == 17 .OR. hmod%imf == 18) THEN
         b_nu = b_Tinker2010_calibrated(nu, hmod)
      ELSE IF (hmod%imf == 19) THEN
         b_nu = b_SMT(nu, hmod)
      ELSE
         STOP 'B_NU: Error, imf not specified correctly'
      END IF

   END FUNCTION b_nu

   REAL FUNCTION b_PS(nu, hmod)

      ! Press & Scheter (1974) halo bias
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod

      b_PS = 1.+(nu**2-1.)/hmod%dc

   END FUNCTION b_PS

   REAL FUNCTION b_ST(nu, hmod)

      ! Sheth & Tormen (1999) halo bias (equation 12 in 9901122)
      ! Comes from peak-background split
      ! Haloes defined with SO relative to mean matter density with SC Delta_v relation
      ! A redshift dependent delta_c is used for barrier height, again from SC
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: dc, p, q

      p = hmod%ST_p
      q = hmod%ST_q
      dc = hmod%dc

      b_ST = 1.+(q*(nu**2)-1.+2.*p/(1.+(q*nu**2)**p))/dc

   END FUNCTION b_ST

   REAL FUNCTION b_SMT(nu, hmod)

      ! Sheth, Mo & Tormen (2001) halo bias
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: dc, anu2, f1, f2, f3, f4
      REAL, PARAMETER :: a = 0.707
      REAL, PARAMETER :: b = 0.5
      REAL, PARAMETER :: c = 0.6

      dc = hmod%dc

      anu2 = a*nu**2

      f1 = sqrt(a)*anu2
      f2 = sqrt(a)*b*anu2**(1.-c)
      f3 = anu2**c
      f4 = anu2**c+b*(1.-c)*(1-c/2.)

      b_SMT = 1.+(f1+f2-f3/f4)/(dc*sqrt(a))

   END FUNCTION b_SMT

   REAL FUNCTION b_Tinker2008(nu, hmod)

      ! Peak-background split applied to Tinker et al. (2008) mass function
      ! Note that the 2008 mass function is not normalised correctly, so this will not have the integral constraint
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: dc, sig
      REAL :: a, b, c

      a = hmod%Tinker_a
      b = hmod%Tinker_b
      c = hmod%Tinker_c

      dc = hmod%dc
      sig = dc/nu

      b_Tinker2008 = 1.-(a*b**a/(sig**a+b**a)-2.*c/sig**2)/dc

   END FUNCTION b_Tinker2008

   REAL FUNCTION b_Tinker2010_calibrated(nu, hmod)

      ! Equation (6) from Tinker et al. (2010; 1001.3162)
      ! Parameter values from Table 2
      ! TODO: Could do an init function?
      ! NOTE: Tinker takes dc = 1.686 and this does not vary with cosmology
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: bigA, bigB, bigC, a, b, c, dc

      bigA = hmod%Tinker_bigA
      bigB = hmod%Tinker_bigB
      bigC = hmod%Tinker_bigC
      a = hmod%Tinker_a
      b = hmod%Tinker_b
      c = hmod%Tinker_c

      dc = hmod%dc

      b_Tinker2010_calibrated = 1.-bigA*nu**a/(nu**a+dc**a)+bigB*nu**b+bigC*nu**c

   END FUNCTION b_Tinker2010_calibrated

   REAL FUNCTION b_Tinker2010_peakbackground(nu, hmod)

      ! Tinker et al. (2010; 1001.3162) halo bias derived using the peak-background split
      ! Equation (15)
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: dc
      REAL :: beta, gamma, phi, eta
      REAL :: f1, f2

      beta = hmod%Tinker_beta
      gamma = hmod%Tinker_gamma
      phi = hmod%Tinker_phi
      eta = hmod%Tinker_eta

      dc = hmod%dc

      f1 = (gamma*nu**2-(1.+2.*eta))/dc
      f2 = (2.*phi/dc)/(1.+(beta*nu)**(2.*phi))
      b_Tinker2010_peakbackground = 1.+f1+f2

   END FUNCTION b_Tinker2010_peakbackground

   REAL FUNCTION b_Jenkins(nu, hmod)

      ! Bias from applying peak-background split to Jenkins (2001) mass function
      ! Haloes defined with FoF with b = 0.2
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: x

      x = log(nu)-0.0876

      IF (x > 0.) THEN
         b_Jenkins = 1.+3.8*abs(x)**2.8/hmod%dc
      ELSE
         b_Jenkins = 1.-3.8*abs(x)**2.8/hmod%dc
         IF (b_Jenkins < 0.) b_Jenkins = 0.
      END IF

   END FUNCTION b_Jenkins

   REAL FUNCTION b_Warren(nu, hmod)

      ! Bias from peak-backgound split applied to Warren et al. (2006) mass function
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: sig, dc
      REAL, PARAMETER :: b = 0.2538
      REAL, PARAMETER :: c = 1.1982

      dc = hmod%dc
      sig = dc/nu

      b_Warren = 1.-(b/(1.+b*sig**b)-2.*c/sig**2)/dc

   END FUNCTION b_Warren

   REAL FUNCTION b_Reed(nu, hmod)

      ! Halo bias from peak-backgound split applied to Reed et al. (2007) mass function
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: crap

      ! Suppress warnings
      crap = nu
      crap = hmod%a

      b_Reed = 1.

   END FUNCTION b_Reed

   REAL FUNCTION b_Bhattacharya(nu, hmod)

      ! Bhattacharya et al. (2011; 1005.2239) mass function for FoF b = 0.2 haloes
      ! Bias derived using the peak-background split
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: f1, f2, f3, dc
      REAL :: a, p, q

      IF (hmod%imf == 10) THEN
         ! Bhattacharya et al. (2011)
         STOP 'G_BHATTACHARYA: Check this carefully'
         a = 0.788/(1.+hmod%z)**0.01
         p = 0.807
         q = 1.795
      ELSE IF (hmod%imf == 14) THEN
         ! Philcox et al. (2020)
         a = 0.774
         p = 0.637
         q = 1.663
      ELSE
         STOP 'G_BHATTACHARYA: Mass function specified incorrectly'
      END IF

      f1 = a*nu**2
      f2 = 2*p/(1.+(a*nu**2)**p)
      f3 = q
      dc = hmod%dc

      b_Bhattacharya = (f1-f2+f3)/dc

   END FUNCTION b_Bhattacharya

   REAL FUNCTION b2_nu(nu, hmod)

      ! Second-order bias function selection
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod

      IF (is_in_array(hmod%imf, [1, 2, 6, 12, 13])) THEN
         b2_nu = b2_ST(nu, hmod)
      ELSE
         STOP 'B2_NU: Error, second-order bias not specified this mass function'
      END IF

   END FUNCTION b2_nu

   REAL FUNCTION b2_ST(nu, hmod)

      ! Sheth, Mo & Tormen (2001) second-order bias
      ! Notation follows from Cooray & Sheth (2002) pp 25-26
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: eps1, eps2, E1, E2, dc, p, q
      REAL, PARAMETER :: a2 = -17./21.

      IF (hmod%imf == 1) THEN
         ! Press & Schecter (1974)
         p = 0.
         q = 1.
         STOP 'B2_ST: Check this very carefully for Press & Schecter form'
      ELSE IF (is_in_array(hmod%imf, [2, 6, 12, 13])) THEN
         p = hmod%ST_p
         q = hmod%ST_q      
      ELSE
         STOP 'B2_ST: Error, incorrect mass function'
      END IF

      dc = hmod%dc

      eps1 = (q*nu**2-1.)/dc
      eps2 = (q*nu**2)*(q*nu**2-3.)/dc**2
      E1 = (2.*p)/(dc*(1.+(q*nu**2)**p))
      E2 = ((1.+2.*p)/dc+2.*eps1)*E1

      b2_ST = 2.*(1.+a2)*(eps1+E1)+eps2+E2

   END FUNCTION b2_ST

   REAL RECURSIVE FUNCTION g_nu(nu, hmod)

      ! Mass function
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod

      IF (.NOT. hmod%has_mass_function) CALL init_mass_function(hmod)

      IF (hmod%imf == 1) THEN
         g_nu = g_PS(nu, hmod)
      ELSE IF (is_in_array(hmod%imf, [2, 6, 12, 13, 19])) THEN
         g_nu = g_ST(nu, hmod)
      ELSE IF (hmod%imf == 3 .OR. hmod%imf == 11 .OR. hmod%imf == 16 .OR. hmod%imf == 17 .OR. hmod%imf == 18) THEN
         g_nu = g_Tinker2010(nu, hmod)
      ELSE IF (hmod%imf == 4) THEN
         STOP 'G_NU: Error, this function should not be used for delta-mass-function'
      ELSE IF (hmod%imf == 5) THEN
         g_nu = g_Jenkins(nu, hmod)
      ELSE IF (hmod%imf == 7 .OR. hmod%imf == 15) THEN
         g_nu = g_Tinker2008(nu, hmod)
      ELSE IF (hmod%imf == 8) THEN
         g_nu = g_Warren(nu, hmod)
      ELSE IF (hmod%imf == 9) THEN
         g_nu = g_Reed(nu, hmod)
      ELSE IF (hmod%imf == 10 .OR. hmod%imf == 14) THEN
         g_nu = g_Bhattacharya(nu, hmod)
      ELSE
         STOP 'G_NU: Error, imf specified incorrectly'
      END IF

      g_nu=g_nu*HMcode_AMF(hmod)

   END FUNCTION g_nu

   REAL FUNCTION HMcode_AMF(hmod)

      TYPE(halomod), INTENT(IN) :: hmod
  
      HMcode_AMF = hmod%Amf+hmod%Amfz*hmod%z

   END FUNCTION HMcode_AMF

   REAL FUNCTION g_mu(mu, hmod)

      ! Mass function g(nu) times dnu/dmu = alpha*mu^(alpha-1); nu=mu^alpha
      ! This transformation makes the integral easier at low nu
      ! TODO: Implement a Taylor expansion for low nu/mu
      REAL, INTENT(IN) :: mu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: nu, alpha

      ! Relation between nu and mu: nu=mu**alpha
      alpha = hmod%alpha_numu
      nu = mu**alpha

      g_mu = g_nu(nu, hmod)*alpha*mu**(alpha-1.)

   END FUNCTION g_mu

   REAL FUNCTION g_PS(nu, hmod)

      ! Press & Scheter (1974) mass function!
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: crap
      REAL, PARAMETER :: A = sqrt(2./pi)

      ! Stop compile-time warnings
      crap = hmod%a

      g_PS = A*exp(-(nu**2)/2.)

   END FUNCTION g_PS

   REAL FUNCTION g_ST(nu, hmod)

      ! Sheth & Tormen (1999) mass function, equation (10) in arXiv:9901122
      ! Note I use nu=dc/sigma(M) while Sheth & Tormen (1999) use nu=(dc/sigma)^2
      ! Haloes defined with SO relative to mean matter density with SC Delta_v relation
      ! A redshift dependent delta_c is used for barrier height, again from SC
      ! TODO: Implement a Taylor expansion for low nu
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: p, q, A

      ! Mass-function parameters
      A = hmod%ST_A
      p = hmod%ST_p
      q = hmod%ST_q

      g_st = A*(1.+((q*nu**2)**(-p)))*exp(-q*nu**2/2.)

   END FUNCTION g_ST

   REAL FUNCTION g_Tinker2008(nu, hmod)

      ! Tinker et al. (2008; 0803.2706) mass function
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: bigA, a, b, c, sig

      bigA = hmod%Tinker_bigA
      a = hmod%Tinker_a
      b = hmod%Tinker_b
      c = hmod%Tinker_c

      sig = hmod%dc/nu

      ! The actual mass function
      g_Tinker2008 = (bigA/nu)*((sig/b)**(-a)+1.)*exp(-c/sig**2)

   END FUNCTION g_Tinker2008

   REAL FUNCTION g_Tinker2010(nu, hmod)

      ! Tinker et al. (2010; 1001.3162) mass function
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: alpha, beta, gamma, phi, eta

      alpha = hmod%Tinker_alpha
      beta = hmod%Tinker_beta
      gamma = hmod%Tinker_gamma
      phi = hmod%Tinker_phi
      eta = hmod%Tinker_eta

      ! The actual mass function
      g_Tinker2010 = alpha*(1.+(beta*nu)**(-2.*phi))*nu**(2.*eta)*exp(-0.5*gamma*nu**2)

   END FUNCTION g_Tinker2010

   REAL FUNCTION g_Jenkins(nu, hmod)

      ! Mass function from astro-ph/0005260
      ! Haloes defined with FoF with b = 0.2
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: x, crap

      ! Prevent compile-time warning
      crap = hmod%a

      x = log(nu)-0.0876

      g_Jenkins = 0.315*exp(-abs(x)**3.8)/nu

   END FUNCTION g_Jenkins

   REAL FUNCTION g_Warren(nu, hmod)

      ! Warren et al. (2006; astro-ph/0506395) mass function for FoF = 0.2 haloes
      ! The paper claims that this relates to ~ 280x mean matter density haloes
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: sig
      REAL, PARAMETER :: bigA = 0.7234 ! Parameters from Equation (8)
      REAL, PARAMETER :: a = 1.625
      REAL, PARAMETER :: b = 0.2538
      REAL, PARAMETER :: c = 1.1982

      sig = hmod%dc/nu

      ! Equation (5), conversion from f(sigma) to g(nu) via 1/nu factor
      g_Warren = (bigA/nu)*(sig**(-a)+b)*exp(-c/sig**2) 

   END FUNCTION g_Warren

   REAL FUNCTION g_Reed(nu, hmod)

      ! Reed et al. (2007; astro-ph/0607150) mass function
      ! NOTE: Paper contains typos in formulae in equations (10) and (12) be very careful
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: omeg, sig, G1, G2, f1, f2, f3, ne
      REAL, PARAMETER :: bigA = 0.310 ! Below equation (12)
      REAL, PARAMETER :: c = 1.08     ! Below equation (11)
      REAL, PARAMETER :: ca = 0.764   ! Below equation (12)
      REAL, PARAMETER :: a = ca/c
      REAL, PARAMETER :: p = 0.3      ! Below equation (5) from Sheth & Tormen

      sig = hmod%dc/nu
      omeg = sqrt(ca)*nu

      ! Equation (12)
      G1 = exp(-(log(omeg)-0.788)**2/(2.*0.6**2))
      G2 = exp(-(log(omeg)-1.138)**2/(2.*0.2**2))

      ! TODO: Implement this
      STOP 'G_REED: Need to implement neff(M) here'
      ne = 0. 

      ! Equation (12)
      ! Note that there are typos in the paper, see my emails
      f1 = 1.+1.02*omeg**(-2.*p)+0.6*G1+0.4*G2 ! TYPO: 2p -> -2p
      f2 = bigA*omeg*sqrt(2./pi)
      f3 = exp(-omeg**2/2.-0.0325*omeg**p/(3.+ne)**2) ! TYPO: omega -> omega^2
      
      g_Reed = f1*f2*f3

   END FUNCTION g_Reed

   REAL FUNCTION g_Bhattacharya(nu, hmod)

      ! Bhattacharya et al. (2011; 1005.2239) mass function for FoF b = 0.2 haloes
      ! Normalised such that all mass is in haloes
      REAL, INTENT(IN) :: nu
      TYPE(halomod), INTENT(IN) :: hmod
      REAL :: f1, f2, f3
      REAL :: bigA, a, p, q, z

      IF (hmod%imf == 10) THEN
         ! Bhattacharya et al. (2011)
         STOP 'G_BHATTACHARYA: Check this carefully'
         z = hmod%z
         bigA = 0.333/(1.+z)**0.11
         a = 0.788/(1.+z)**0.01
         p = 0.807
         q = 1.795
      ELSE IF (hmod%imf == 14) THEN
         ! Philcox et al. (2020)
         bigA = 0.33873 ! Where exactly does this come from?
         a = 0.774
         p = 0.637
         q = 1.663
      ELSE
         STOP 'G_BHATTACHARYA: Mass function specified incorrectly'
      END IF

      f1 = bigA*sqrt(2./pi)*exp(-a*nu**2/2.)
      f2 = 1.+(a*nu**2)**(-p)
      f3 = (a*nu**2)**(q/2.)
      
      g_Bhattacharya = f1*f2*f3/nu

   END FUNCTION g_Bhattacharya

   SUBROUTINE init_mass_function(hmod)

      ! Initialise anything to do with the halo mass function
      TYPE(halomod), INTENT(INOUT) :: hmod

      IF (is_in_array(hmod%imf, [2, 6, 12, 13, 19])) THEN
         CALL init_ST(hmod)
      ELSE IF (is_in_array(hmod%imf, [3, 11, 16, 17, 18])) THEN
         CALL init_Tinker2010(hmod)
      ELSE IF (hmod%imf == 7 .OR. hmod%imf == 15) THEN
         CALL init_Tinker2008(hmod)
      END IF

      hmod%has_mass_function = .TRUE.

   END SUBROUTINE init_mass_function

   SUBROUTINE init_ST(hmod)

      ! Normalises ST mass function
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: p, q, A

      ! ST parameters
      IF (hmod%imf == 2 .OR. hmod%imf == 19) THEN
         ! Sheth & Tormen (1999)
         p = 0.3
         q = 0.707
      ELSE IF (hmod%imf == 6) THEN
         ! Despali et al. (2016)      
         p = 0.2579
         q = 0.7663
         A = 0.3298
      ELSE IF (hmod%imf == 12) THEN
         ! Sheth & Tormen (2002)
         p = 0.3
         q = 0.75
      ELSE IF (hmod%imf == 13) THEN
         ! Sheth & Tormen form
         p = hmod%ST_p
         q = hmod%ST_q
      ELSE
         STOP 'INIT_ST: Error, imf set incorrectly'
      END IF

      ! Normalisation of ST mass function (involves Gamma function)
      IF (is_in_array(hmod%imf, [2, 6, 13])) THEN
         !hmod%ST_A = 1./(sqrt(pi/(2.*q))+(1./sqrt(q))*(2.**(-p-0.5))*Gamma(0.5-p))
         A = sqrt(2.*q)/(sqrt(pi)+Gamma(0.5-p)/2**p)
      END IF

      ! Set the parameters
      hmod%ST_A = A
      hmod%ST_p = p
      hmod%ST_q = q

   END SUBROUTINE init_ST

   SUBROUTINE init_Tinker2008(hmod)

      ! Initialise the parameters of the Tinker et al. (2010) mass function and bias
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: bigA, a, b, c
      REAL :: z, log_Dv, alpha
      LOGICAL :: z_dependence

      ! Parameter arrays from Tinker (2008): Table 2
      INTEGER, PARAMETER :: n = 9 ! Number of entries in parameter lists
      REAL, PARAMETER :: log_Deltav(n) = log([200., 300., 400., 600., 800., 1200., 1600., 2400., 3200.])
      REAL, PARAMETER :: bigAs(n)= [0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260]
      REAL, PARAMETER :: as(n) = [1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66]
      REAL, PARAMETER :: bs(n) = [2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41]
      REAL, PARAMETER :: cs(n) = [1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44]
      REAL, PARAMETER :: bigA_z_exp = -0.14
      REAL, PARAMETER :: a_z_exp = -0.06

      INTEGER, PARAMETER :: iorder = 1 ! Order for interpolation
      INTEGER, PARAMETER :: ifind = 3  ! Scheme for finding (3 - Mid-point splitting)
      INTEGER, PARAMETER :: imeth = 2  ! Method for interpolation (2 - Lagrange polynomial)

      ! Choose redshift dependence
      IF (hmod%imf == 15) THEN
         z_dependence = .FALSE.
      ELSE
         z_dependence = .TRUE.
      END IF

      ! Get these from the halo-model structure
      z = hmod%z
      IF (hmod%mass_dependent_Dv) THEN
         STOP 'INIT_TINKER2008: Error, this does not work with mass-dependent Delta_v'
      ELSE
         log_Dv = log(hmod%Dv)
      END IF

      ! Delta_v dependence (changed to log Dv finding)
      bigA = find(log_Dv, log_Deltav, bigAs, n, iorder, ifind, imeth)
      a = find(log_Dv, log_Deltav, as, n, iorder, ifind, imeth)
      b = find(log_Dv, log_Deltav, bs, n, iorder, ifind, imeth)
      c = find(log_Dv, log_Deltav, cs, n, iorder, ifind, imeth)

      ! Redshift dependence
      IF (z_dependence) THEN
         bigA = bigA*(1.+z)**bigA_z_exp                    ! Equation (5)
         a = a*(1.+z)**a_z_exp                             ! Equation (6)
         alpha = 10**(-(0.75/log10(exp(log_Dv)/75.))**1.2) ! Equation (8)
         b = b*(1.+z)**(-alpha)                            ! Equation (7)
      END IF

      hmod%Tinker_bigA = bigA
      hmod%Tinker_a = a
      hmod%Tinker_b = b
      hmod%Tinker_c = c

   END SUBROUTINE init_Tinker2008

   SUBROUTINE init_Tinker2010(hmod)

      ! Initialise the parameters of the Tinker et al. (2010) mass function and bias
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL :: alpha, beta, Gamma, phi, eta
      REAL :: z, log_Dv, y
      LOGICAL :: z_dependence

      ! Parameter arrays from Tinker (2010): Table 4
      ! NOTE: alpha parameters below are for normalisation assuming no redshift evolution
      ! NOTE: If parameters have explicit redshift evolution then alpha should be calculated
      INTEGER, PARAMETER :: n = 9 ! Number of entries in parameter lists
      REAL, PARAMETER :: log_Deltav(n) = log([200., 300., 400., 600., 800., 1200., 1600., 2400., 3200.])
      REAL, PARAMETER :: alpha0(n) = [0.368, 0.363, 0.385, 0.389, 0.393, 0.365, 0.379, 0.355, 0.327]
      REAL, PARAMETER :: beta0(n) = [0.589, 0.585, 0.544, 0.543, 0.564, 0.623, 0.637, 0.673, 0.702]
      REAL, PARAMETER :: gamma0(n) = [0.864, 0.922, 0.987, 1.09, 1.20, 1.34, 1.50, 1.68, 1.81]
      REAL, PARAMETER :: phi0(n) = [-0.729, -0.789, -0.910, -1.05, -1.20, -1.26, -1.45, -1.50, -1.49]
      REAL, PARAMETER :: eta0(n) = [-0.243, -0.261, -0.261, -0.273, -0.278, -0.301, -0.301, -0.319, -0.336]
      REAL, PARAMETER :: beta_z_exp = 0.20
      REAL, PARAMETER :: gamma_z_exp = -0.01
      REAL, PARAMETER :: phi_z_exp = -0.08
      REAL, PARAMETER :: eta_z_exp = 0.27

      INTEGER, PARAMETER :: iorder = 1 ! Order for interpolation
      INTEGER, PARAMETER :: ifind = 3  ! Scheme for finding (3 - Mid-point splitting)
      INTEGER, PARAMETER :: imeth = 2  ! Method for interpolation (2 - Lagrange polynomial)

      ! Choose redshift dependence
      ! NOTE: Redshift-dependent parameters are only given for M200 haloes in Tinker (2010) paper
      ! TODO: Should the Dv = 200 condition be enforced?
      IF (hmod%Dv == 200. .AND. (hmod%imf == 3 .OR. hmod%imf == 11 .OR. hmod%imf == 17)) THEN
         z_dependence = .TRUE.
      ELSE
         z_dependence = .FALSE.
      END IF

      ! Get these from the halo-model structure
      z = hmod%z
      IF (z > 3.) z = 3. ! Recommendation from Tinker et al. (2010)
      IF (hmod%mass_dependent_Dv) THEN
         STOP 'INIT_TINKER2010: Error, this does not work with mass-dependent Delta_v'
      ELSE
         log_Dv = log(hmod%Dv)
      END IF

      ! Delta_v dependence (changed to log Dv finding)
      alpha = find(log_Dv, log_Deltav, alpha0, n, iorder, ifind, imeth)
      beta = find(log_Dv, log_Deltav, beta0, n, iorder, ifind, imeth)
      gamma = find(log_Dv, log_Deltav, gamma0, n, iorder, ifind, imeth)
      phi = find(log_Dv, log_Deltav, phi0, n, iorder, ifind, imeth)
      eta = find(log_Dv, log_Deltav, eta0, n, iorder, ifind, imeth)

      ! Redshift dependence
      IF (z_dependence) THEN
         beta = beta*(1.+z)**beta_z_exp    ! Equation (9)
         gamma = gamma*(1.+z)**gamma_z_exp ! Equation (12)
         phi = phi*(1.+z)**phi_z_exp       ! Equation (10)
         eta = eta*(1.+z)**eta_z_exp       ! Equation (11)
      END IF

      ! Set the Tinker parameters
      IF (z_dependence) THEN
         hmod%Tinker_alpha = 1. ! Fix to unity before normalisation
      ELSE
         hmod%Tinker_alpha = alpha
      END IF
      hmod%Tinker_beta = beta
      hmod%Tinker_gamma = gamma
      hmod%Tinker_phi = phi
      hmod%Tinker_eta = eta

      ! Set the flag
      hmod%has_mass_function = .TRUE.

      ! Explicitly normalise
      IF (z_dependence) THEN
         hmod%Tinker_alpha = hmod%Tinker_alpha/integrate_g_mu(hmod%small_nu, hmod%large_nu, hmod)
      END IF

      ! Calibrated bias initialisation
      IF (hmod%imf == 17 .OR. hmod%imf == 18) THEN

         y = log10(hmod%Dv)
         
         hmod%Tinker_bigA = 1.+0.24*y*exp(-(4./y)**4)
         hmod%Tinker_a = 0.44*y-0.88
         hmod%Tinker_bigB = 0.183
         hmod%Tinker_b = 1.5
         hmod%Tinker_bigC = 0.019+0.107*y+0.19*exp(-(4./y)**4)
         hmod%Tinker_c = 2.4 

      END IF

   END SUBROUTINE init_Tinker2010

   REAL FUNCTION gb_nu(nu, hmod)

      ! g(nu)*b(nu); useful as an integrand
      REAL, INTENT(IN) :: nu               ! Proxy mass variable
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

      gb_nu = g_nu(nu, hmod)*b_nu(nu, hmod)

   END FUNCTION gb_nu

   REAL FUNCTION gb_nu_on_M(nu, hmod)

      ! g(nu)*b(nu)/M(nu); useful as an integrand
      REAL, INTENT(IN) :: nu               ! Proxy mass variable
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

      gb_nu_on_M = g_nu(nu, hmod)*b_nu(nu, hmod)/M_nu(nu, hmod)

   END FUNCTION gb_nu_on_M

   REAL FUNCTION nug_nu(nu, hmod)

      ! nu*g(nu); useful as an integrand
      REAL, INTENT(IN) :: nu               ! Proxy mass variable
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

      nug_nu = nu*g_nu(nu, hmod)

   END FUNCTION nug_nu

   REAL FUNCTION Mg_nu(nu, hmod)

      ! M(nu)*g(nu); useful as an integrand
      REAL, INTENT(IN) :: nu               ! Proxy mass variable
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

      Mg_nu = M_nu(nu, hmod)*g_nu(nu, hmod)

   END FUNCTION Mg_nu

   REAL FUNCTION nug_nu_on_M(nu, hmod)

      ! nu*g(nu)/M(nu); useful as an integrand
      REAL, INTENT(IN) :: nu               ! Proxy mass variable
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

      nug_nu_on_M = nu*g_nu(nu, hmod)/M_nu(nu, hmod)

   END FUNCTION nug_nu_on_M

   REAL FUNCTION g_nu_on_M(nu, hmod)

      ! g(nu)/M(nu); useful as an integrand
      REAL, INTENT(IN) :: nu               ! Proxy mass variable
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model

      g_nu_on_M = g_nu(nu, hmod)/M_nu(nu, hmod)

   END FUNCTION g_nu_on_M

   FUNCTION wk_isothermal(x)

      ! The normlaised Fourier Transform of an isothermal profile
      REAL :: wk_isothermal
      REAL, INTENT(IN) :: x

      REAL, PARAMETER :: dx = 1e-3

      ! Taylor expansion used for low |x| to avoid cancellation problems

      IF (ABS(x) < ABS(dx)) THEN
         ! Taylor series at low x
         wk_isothermal = 1.-(x**2)/18.
      ELSE
         wk_isothermal = Si(x)/x
      END IF

   END FUNCTION wk_isothermal

   FUNCTION wk_isothermal_2(x, y)

      ! The normlaised Fourier Transform of an isothemral profile from x -> y
      REAL :: wk_isothermal_2
      REAL, INTENT(IN) :: x, y

      wk_isothermal_2 = (Si(x)-Si(y))/(x-y)

   END FUNCTION wk_isothermal_2

   REAL FUNCTION halo_fraction(itype, m, hmod, cosm)

      ! Mass fraction of a type within a halo
      INTEGER, INTENT(IN) :: itype
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (itype == field_dmonly .OR. itype == field_matter) THEN
         ! TODO: Is this necessary? Is this ever called?
         halo_fraction = 1.
      ELSE IF (itype == field_cdm) THEN
         halo_fraction = halo_CDM_fraction(m, hmod, cosm)
      ELSE IF (itype == field_gas) THEN
         halo_fraction = halo_gas_fraction(m, hmod, cosm)
      ELSE IF (itype == field_stars) THEN
         halo_fraction = halo_star_fraction(m, hmod, cosm)
      ELSE IF (itype == field_bound_gas) THEN
         halo_fraction = halo_bound_gas_fraction(m, hmod, cosm)
      ELSE IF (itype == field_ejected_gas) THEN
         halo_fraction = halo_ejected_gas_fraction(m, hmod, cosm)
      ELSE IF (itype == field_cold_bound_gas) THEN
         halo_fraction = halo_cold_bound_gas_fraction(m, hmod, cosm)
      ELSE IF (itype == field_hot_bound_gas) THEN
         halo_fraction = halo_hot_bound_gas_fraction(m, hmod, cosm)
      ELSE IF (itype == field_normal_bound_gas) THEN
         halo_fraction = halo_normal_bound_gas_fraction(m, hmod, cosm)
      ELSE IF (itype == field_central_stars) THEN
         halo_fraction = halo_central_star_fraction(m, hmod, cosm)
      ELSE IF (itype == field_satellite_stars) THEN
         halo_fraction = halo_satellite_star_fraction(m, hmod, cosm)
      ELSE IF (itype == field_neutrino) THEN
         halo_fraction = halo_neutrino_fraction(m, hmod, cosm)
      ELSE
         STOP 'HALO_FRACTION: Error, itype not specified correcntly'
      END IF

   END FUNCTION halo_fraction

   REAL FUNCTION halo_CDM_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in CDM
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: crap

      ! To prevent compile-time warning
      crap = m
      crap = hmod%a

      ! Always the universal value
      halo_CDM_fraction = cosm%om_c/cosm%om_m

   END FUNCTION halo_CDM_fraction

   REAL FUNCTION halo_neutrino_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in neutrinos
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: crap

      ! To prevent compile-time warning
      crap = m
      crap = hmod%a

      ! Always the universal value
      halo_neutrino_fraction = cosm%f_nu

   END FUNCTION halo_neutrino_fraction

   REAL FUNCTION halo_gas_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in gas
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      halo_gas_fraction = halo_bound_gas_fraction(m, hmod, cosm)+halo_ejected_gas_fraction(m, hmod, cosm)

      IF(halo_gas_fraction < frac_min .OR. halo_gas_fraction > frac_max) THEN
         WRITE(*, *) 'HALO_GAS_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_GAS_FRACTION: Halo gas fraction:', halo_gas_fraction
         STOP 'HALO_GAS_FRACTION: Error, halo gas fraction is outside sensible range'
      END IF

   END FUNCTION halo_gas_fraction

   REAL FUNCTION halo_bound_gas_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in bound gas
      ! This function should NOT reference halo_star_fraction() otherwise infinite recursion generated
      ! TODO: may be wise to have the prefactor in all of these equations be: (Om_b/Om_m - f_*)
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: m0, sig, beta, x, maxgas, f_star, f_total

      IF (hmod%normalise_baryons == 1 .OR. hmod%normalise_baryons == 3) THEN
         maxgas = cosm%om_b/cosm%om_m
      ELSE IF (hmod%normalise_baryons == 2) THEN
         maxgas = cosm%om_b/cosm%om_m - halo_star_fraction(m, hmod, cosm)
      ELSE
         STOP 'HALO_BOUND_GAS_FRACTION: Error, normalise_baryons not specified correctly'
      END IF

      IF (hmod%frac_bound_gas == 1) THEN
         ! From Fedeli (2014a)
         m0 = 1e12
         sig = 3.
         IF (m < m0) THEN
            halo_bound_gas_fraction = 0.
         ELSE
            halo_bound_gas_fraction = erf(log10(m/m0)/sig)*maxgas
         END IF
      ELSE IF (hmod%frac_bound_gas == 2) THEN
         ! From Schneider & Teyssier (2015)
         M0 = HMx_M0(hmod, cosm)
         beta = HMx_gbeta(hmod, cosm)
         x = (m/M0)**beta
         halo_bound_gas_fraction = maxgas*x/(1.+x)
      ELSE IF (hmod%frac_bound_gas == 3) THEN
         ! Universal baryon fraction model
         ! Be *very* careful with generating an infinite recursion here
         ! DO NOT ACCOUNT FOR STAR FRACTION HERE DUE TO INFINITE RECURSION
         ! Stars are subtracted from elsewhere if this is a problem
         !halo_bound_gas_fraction = cosm%om_b/cosm%om_m-halo_star_fraction(m, hmod, cosm) ! ABSOLUTELY DO NOT DO THIS
         halo_bound_gas_fraction = maxgas
      ELSE
         STOP 'HALO_BOUND_GAS_FRACTION: Error, frac_bound_gas not specified correctly'
      END IF

      ! Remove bound gas if the total gas+star fraction would be greater than universal baryon fraction
      IF (hmod%normalise_baryons == 3) THEN
         f_star = halo_star_fraction(m, hmod, cosm) ! Star fraction in halo
         f_total = f_star+halo_bound_gas_fraction   ! Total baryon content inside halo
         IF(f_total > maxgas) THEN
            halo_bound_gas_fraction = halo_bound_gas_fraction-(f_total-maxgas)
         END IF
      END IF

      IF (halo_bound_gas_fraction < frac_min .OR. halo_bound_gas_fraction > frac_max) THEN
         WRITE(*, *) 'HALO_BOUND_GAS_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_BOUND_GAS_FRACTION:', halo_bound_gas_fraction
         STOP 'HALO_BOUND_GAS_FRACTION: Error, halo bound gas fraction is outside sensible range'
      END IF

   END FUNCTION halo_bound_gas_fraction

   REAL FUNCTION halo_normal_bound_gas_fraction(m, hmod, cosm)

      ! Mass fraction of total gas that is static and bound in haloes
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: bound, cold, hot

      bound = halo_bound_gas_fraction(m, hmod, cosm)
      cold = halo_cold_bound_gas_fraction(m, hmod, cosm)
      hot = halo_hot_bound_gas_fraction(m, hmod, cosm)

      halo_normal_bound_gas_fraction = bound-cold-hot

      IF(halo_normal_bound_gas_fraction < frac_min .OR. halo_normal_bound_gas_fraction > frac_max) THEN
         WRITE(*, *) 'HALO_STATIC_GAS_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_STATIC_GAS_FRACTION:', halo_normal_bound_gas_fraction
         STOP 'HALO_STATIC_GAS_FRACTION: Error, static gas fraction is outside sensible range'
      END IF

   END FUNCTION halo_normal_bound_gas_fraction

   REAL FUNCTION halo_cold_bound_gas_fraction(m, hmod, cosm)

      ! Mass fraction of total gas that is cold and bound in haloes
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r

      IF (hmod%frac_cold_bound_gas == 1) THEN
         ! Constant fraction of cold halo gas
         r = HMx_fcold(hmod, cosm)
      ELSE
         STOP 'HALO_COLD_GAS_FRACTION: Error, frac_cold_bound_gas not specified correctly'
      END IF

      IF(r < frac_min .OR. r > frac_max) THEN
         WRITE(*, *) 'HALO_COLD_GAS_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_COLD_GAS_FRACTION: r:', r
         STOP 'HALO_COLD_GAS_FRACTION: Error, fraction of bound gas that is cold is outside sensible range'
      END IF
      halo_cold_bound_gas_fraction = r*halo_bound_gas_fraction(m, hmod, cosm)

   END FUNCTION halo_cold_bound_gas_fraction

   REAL FUNCTION halo_hot_bound_gas_fraction(m, hmod, cosm)

      ! Mass fraction of total gas that is hot and bound in haloes
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r

      IF (hmod%frac_hot_bound_gas == 1) THEN
         ! Constant fraction of hot halo gas
         r = HMx_fhot(hmod, cosm)
      ELSE
         STOP 'HALO_HOT_GAS_FRACTION: Error, frac_hot_bound_gas not specified correctly'
      END IF

      IF(r < frac_min .OR. r > frac_max) THEN
         WRITE(*, *) 'HALO_HOT_GAS_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_HOT_GAS_FRACTION: r:', r
         STOP 'HALO_HOT_GAS_FRACTION: Error, fraction of bound gas that is hot must be between zero and one'
      END IF
      halo_hot_bound_gas_fraction = r*halo_bound_gas_fraction(m, hmod, cosm)

   END FUNCTION halo_hot_bound_gas_fraction

   REAL FUNCTION halo_ejected_gas_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in free gas
      ! This fucntion can reference halo_bound_gas_fraction() and halo_star_fraction()
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      ! This is necessarily all the gas that is not bound or in stars
      halo_ejected_gas_fraction = cosm%om_b/cosm%om_m-halo_star_fraction(m, hmod, cosm)-halo_bound_gas_fraction(m, hmod, cosm)

      IF (halo_ejected_gas_fraction < frac_min .OR. halo_ejected_gas_fraction > frac_max) THEN
         WRITE(*, *) 'HALO_FREE_GAS_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_FREE_GAS_FRACTION: Baryon fraction:', cosm%om_b/cosm%om_m
         WRITE(*, *) 'HALO_FREE_GAS_FRACTION: Halo star fraction:', halo_star_fraction(m, hmod, cosm)
         WRITE(*, *) 'HALO_FREE_GAS_FRACTION: Halo bound-gas fraction:', halo_bound_gas_fraction(m, hmod, cosm)        
         WRITE(*, *) 'HALO_FREE_GAS_FRACTION: Halo free-gas fraction:', halo_ejected_gas_fraction
         STOP 'HALO_FREE_GAS_FRACTION: Error, free-gas fraction is outside sensible range'
      END IF

   END FUNCTION halo_ejected_gas_fraction

   REAL FUNCTION halo_star_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in stars
      ! This fucntion can reference halo_bound_gas_fraction()
      ! TODO: Remove calls to halo_bound_gas fraction to remove infinite recursions
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: m0, sig, A
      REAL :: f_bound, f_total

      IF (hmod%frac_stars == 1 .OR. hmod%frac_stars == 3) THEN
         ! Fedeli (2014)
         A = HMx_Astar(hmod, cosm)
         m0 = HMx_Mstar(hmod, cosm)
         sig = HMx_sstar(hmod, cosm)
         halo_star_fraction = A*exp(-((log10(m/m0))**2)/(2.*sig**2))
         IF (hmod%frac_stars == 3) THEN
            ! Suggested by Ian, the relation I have is for the central stellar mass
            ! in reality this saturates for high-mass haloes (due to satellite contribution)
            IF (halo_star_fraction < A/3. .AND. m > m0) halo_star_fraction = A/3.
         END IF
      ELSE IF (hmod%frac_stars == 2) THEN
         ! Constant star fraction
         halo_star_fraction = HMx_Astar(hmod, cosm)
      ELSE
         STOP 'HALO_STAR_FRACTION: Error, frac_stars specified incorrectly'
      END IF

      ! Take away from the stars if the sum of bound gas and stars in the halo would be greater than universal baryon fraction
      ! Be careful with generating an infinite recursion here
      ! TODO: Could remove this an subtract excess stars from prefactor in bound gas
      IF (hmod%normalise_baryons == 1) THEN
         f_bound = halo_bound_gas_fraction(m, hmod, cosm) ! Bound gas in halo
         f_total = f_bound+halo_star_fraction             ! Total baryon content inside halo
         IF(f_total > cosm%om_b/cosm%om_m) THEN
            halo_star_fraction = halo_star_fraction-(f_total-cosm%om_b/cosm%om_m)
         END IF
      END IF

      IF(halo_star_fraction < frac_min .OR. halo_star_fraction > frac_max) THEN
         WRITE(*, *) 'HALO_STAR_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_STAR_FRACTION: Halo star fraction', halo_star_fraction
         STOP 'HALO_STAR_FRACTION: Error, star fraction must be between zero and one'
      END IF

   END FUNCTION halo_star_fraction

   REAL FUNCTION halo_central_star_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in central stars
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r

      IF (hmod%frac_central_stars == 1) THEN
         ! All stellar mass in centrals
         r = 1.
      ELSE IF (hmod%frac_central_stars == 2) THEN
         ! Schnedier et al. (2018)
         IF (M <= HMx_Mstar(hmod, cosm)) THEN
            ! For M<M* all galaxy mass is in centrals
            r = 1.
         ELSE
            ! Otherwise there is a power law, eta ~ -0.3, higher mass haloes have more mass in satellites
            r = (m/HMx_Mstar(hmod, cosm))**HMx_eta(hmod, cosm)
         END IF       
      ELSE
         STOP 'HALO_CENTRAL_STAR_FRACTION: Error, frac_central_stars specified incorrectly'
      END IF

      IF(r < frac_min .OR. r > frac_max) THEN
         WRITE(*, *) 'HALO_CENTRAL_STAR_FRACTION: Halo mass [log10(Msun/h)]:', log10(m)
         WRITE(*, *) 'HALO_CENTRAL_STAR_FRACTION: r:', r
         STOP 'HALO_CENTRAL_STAR_FRACTION: Error, fraction of stars that are central must be between zero and one'
      END IF
      halo_central_star_fraction = r*halo_star_fraction(m, hmod, cosm)

   END FUNCTION halo_central_star_fraction

   REAL FUNCTION halo_satellite_star_fraction(m, hmod, cosm)

      ! Mass fraction of a halo in satellite stars
      REAL, INTENT(IN) :: m
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm

      halo_satellite_star_fraction = halo_star_fraction(m, hmod, cosm)-halo_central_star_fraction(m, hmod, cosm)

   END FUNCTION halo_satellite_star_fraction

   REAL FUNCTION halo_HI_fraction(M, hmod, cosm)

      ! Dimensionless M_HI/M_halo
      ! TODO: This should probably be related to the cold gas fraction
      REAL, INTENT(IN) :: M ! Host-halo mass
      TYPE(halomod), INTENT(INOUT) :: hmod ! Halo model
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL :: alpha, M0, Mmin
      REAL :: crap

      ! Prevent compile-time warnings
      crap = cosm%A

      IF (hmod%frac_HI == 1) THEN
         ! Simple model with hard truncation
         IF (m >= hmod%HImin .AND. m <= hmod%HImax) THEN
            halo_HI_fraction = 1.
         ELSE
            halo_HI_fraction = 0.
         END IF
      ELSE IF (hmod%frac_HI == 2 .OR. hmod%frac_HI == 3) THEN
         ! From Villaescusa-Navarro et al. (2018; 1804.09180)
         ! 2 - Just z=0 results
         ! 3 - Using z evolution
         IF (hmod%frac_HI == 2 .OR. hmod%z == 0.) THEN
            ! FoF values from Table 1
            !alpha=0.24
            !M0=4.3e10
            !Mmin=2e12
            ! FoF-SO values from Table 1
            alpha = 0.16
            M0 = 4.1e10
            Mmin = 2.4e12
         ELSE IF (hmod%z == 1.) THEN
            ! FoF values from Table 1
            !alpha=0.53
            !M0=1.5e10
            !Mmin=6e11
            ! FoF-SO values from Table 1
            alpha = 0.43
            M0 = 1.8e10
            Mmin = 8.6e11
         ELSE IF (hmod%z == 2.) THEN
            ! FoF values from Table 1
            !alpha=0.60
            !M0=1.3e10
            !Mmin=3.6e11
            ! FoF-SO values from Table 1
            alpha = 0.51
            M0 = 1.5e10
            Mmin = 4.6e11
         ELSE IF (hmod%z == 3.) THEN
            ! FoF values from Table 1
            !alpha=0.76
            !M0=2.9e9
            !Mmin=6.7e10
            ! FoF-SO values from Table 1
            alpha = 0.69
            M0 = 3.7e9
            Mmin = 9.6e10
         ELSE IF (hmod%z == 4.) THEN
            ! FoF values from Table 1
            !alpha=0.79
            !M0=1.4e9
            !Mmin=2.1e10
            ! FoF-SO values from Table 1
            alpha = 0.61
            M0 = 4.5e9
            Mmin = 7.6e10
         ELSE IF (hmod%z == 5.) THEN
            ! FoF values from Table 1
            !alpha=0.74
            !M0=1.9e9
            !Mmin=2e10
            ! FoF-SO values from Table 1
            alpha = 0.59
            M0 = 4.1e9
            Mmin = 5.4e10
         ELSE
            STOP 'HALO_HI_FRACTION: Error, redshift not supported here'
         END IF
         halo_HI_fraction = (M0/M)*((M/Mmin)**alpha)*exp(-(M/Mmin)**(-0.35))
      ELSE
         STOP 'HALO_HI_FRACTION: Error, frac_HI specified incorrectly'
      END IF

   END FUNCTION halo_HI_fraction

   REAL RECURSIVE FUNCTION integrate_hmod(a, b, f, hmod, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, EXTERNAL :: f
      TYPE(halomod), INTENT(INOUT) :: hmod
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = 5
      INTEGER, PARAMETER :: jmax = 30

      INTERFACE
         FUNCTION f(x_interface, hmod_interface)
            IMPORT :: halomod
            REAL, INTENT(IN) :: x_interface
            TYPE(halomod), INTENT(INOUT) :: hmod_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_hmod = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.d0
         sum_n = 0.d0
         sum_old = 0.d0
         sum_new = 0.d0

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/REAL(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, hmod)
               f2 = f(b, hmod)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*REAL(i-1)/REAL(n-1)
                  fx = f(x, hmod)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.d0+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_HMOD: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_HMOD: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO
         
         integrate_hmod = real(sum_new)

      END IF

   END FUNCTION integrate_hmod

   REAL RECURSIVE FUNCTION integrate_hmod_cosm(a, b, f, hmod, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, EXTERNAL :: f
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = 5
      INTEGER, PARAMETER :: jmax = 30

      INTERFACE
         FUNCTION f(x_interface, hmod_interface, cosm_interface)
            IMPORT :: halomod
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            TYPE(halomod), INTENT(INOUT) :: hmod_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_hmod_cosm = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.d0
         sum_n = 0.d0
         sum_old = 0.d0
         sum_new = 0.d0

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/REAL(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, hmod, cosm)
               f2 = f(b, hmod, cosm)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*REAL(i-1)/REAL(n-1)
                  fx = f(x, hmod, cosm)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.d0+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_HMOD_COSM: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_HMOD_COSM: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO
         
         integrate_hmod_cosm = REAL(sum_new)

      END IF

   END FUNCTION integrate_hmod_cosm

   REAL RECURSIVE FUNCTION integrate_hmod_cosm_exp(a, b, f, hmod, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, EXTERNAL :: f
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = 5
      INTEGER, PARAMETER :: jmax = 30

      INTERFACE
         FUNCTION f(nu, hmod_interface, cosm_interface)
            IMPORT :: halomod
            IMPORT :: cosmology
            REAL, INTENT(IN) :: nu
            TYPE(halomod), INTENT(INOUT) :: hmod_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_hmod_cosm_exp = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.d0
         sum_n = 0.d0
         sum_old = 1.d0 ! Should not be zero
         sum_new = 0.d0

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/REAL(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(exp(a), hmod, cosm)*exp(a)
               f2 = f(exp(b), hmod, cosm)*exp(b)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*REAL(i-1)/REAL(n-1)
                  fx = f(exp(x), hmod, cosm)*exp(x)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.d0+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_HMOD_COSM_EXP: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_HMOD_COSM_EXP: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO

         integrate_hmod_cosm_exp = REAL(sum_new)

      END IF

   END FUNCTION integrate_hmod_cosm_exp

   REAL FUNCTION integrate_scatter(c, dc, ih, k, m, rv, hmod, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: c
      REAL, INTENT(IN) :: dc
      INTEGER, INTENT(IN) :: ih(2)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: rv
      TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      REAL :: a, b
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = 5
      INTEGER, PARAMETER :: jmax = 30
      REAL, PARAMETER :: nsig = 5

      a = c/(1.+nsig*dc)
      b = c*(1.+nsig*dc)

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_scatter = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.d0
         sum_n = 0.d0
         sum_old = 1.d0 ! Should not be zero
         sum_new = 0.d0

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/REAL(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = scatter_integrand(a, c, dc, ih, k, m, rv, hmod, cosm)
               f2 = scatter_integrand(b, c, dc, ih, k, m, rv, hmod, cosm)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*REAL(i-1)/REAL(n-1)
                  fx = scatter_integrand(x, c, dc, ih, k, m, rv, hmod, cosm)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.d0+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_SCATTER: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_SCATTER: Integration timed out'
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO

         integrate_scatter = REAL(sum_new)

      END IF

   END FUNCTION integrate_scatter

   SUBROUTINE write_power(k, pow_li, pow_2h, pow_1h, pow_hm, nk, output, verbose)

      ! Write P(k) of all halo-model terms
      INTEGER, INTENT(IN) :: nk
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(IN) :: pow_li(nk)
      REAL, INTENT(IN) :: pow_2h(nk)
      REAL, INTENT(IN) :: pow_1h(nk)
      REAL, INTENT(IN) :: pow_hm(nk)
      CHARACTER(len=*), INTENT(IN) :: output   
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i

      IF (verbose) THEN
         WRITE (*, *) 'WRITE_POWER: Writing power to ', TRIM(output)
      END IF

      ! Loop over k values
      OPEN (7, file=output)
      DO i = 1, nk
         WRITE (7, fmt='(5ES20.10)') k(i), pow_li(i), pow_2h(i), pow_1h(i), pow_hm(i)
      END DO
      CLOSE (7)

      IF (verbose) THEN
         WRITE (*, *) 'WRITE_POWER: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE write_power

   SUBROUTINE write_power_fields(k, pow_li, pow_2h, pow_1h, pow_hm, nk, fields, nf, base, verbose)

      ! Write P(k) for all halo-model components with a different file for each field pair
      INTEGER, INTENT(IN) :: nk
      INTEGER, INTENT(IN) :: nf
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(IN) :: pow_li(nk)
      REAL, INTENT(IN) :: pow_2h(nf, nf, nk)
      REAL, INTENT(IN) :: pow_1h(nf, nf, nk)
      REAL, INTENT(IN) :: pow_hm(nf, nf, nk)
      INTEGER, INTENT(IN) :: fields(nf)
      CHARACTER(len=*), INTENT(IN) :: base   
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: j1, j2
      CHARACTER(len=256) :: outfile
      CHARACTER(len=256) :: ext='.dat'

      IF(verbose) WRITE(*, *) 'WRITE_POWER_FIELDS: Writing all field power'

      DO j1 = 1, nf
         DO j2 = 1, nf
            outfile = number_file2(base, fields(j1), '', fields(j2), ext)
            IF(verbose) WRITE(*,*) 'WRITE_POWER_FIELDS: Output file: ', trim(outfile)
            CALL write_power(k, pow_li, pow_2h(j1, j2, :), pow_1h(j1, j2, :), pow_hm(j1, j2, :), nk, outfile, verbose=.FALSE.)
         END DO
      END DO

      IF(verbose) THEN
         WRITE(*, *) 'WRITE_POWER_FIELDS: Done'
         WRITE(*, *)
      END IF
   
   END SUBROUTINE write_power_fields

   SUBROUTINE write_power_a_multiple(k, a, pow_li, pow_2h, pow_1h, pow_hm, nk, na, base, verbose)

      ! Write P(k,a) for all halo-model components with a separte file for each component
      INTEGER, INTENT(IN) :: nk
      INTEGER, INTENT(IN) :: na
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(IN) :: a(na)
      REAL, INTENT(IN) :: pow_li(nk, na)
      REAL, INTENT(IN) :: pow_2h(nk, na)
      REAL, INTENT(IN) :: pow_1h(nk, na)
      REAL, INTENT(IN) :: pow_hm(nk, na)
      CHARACTER(len=*), INTENT(IN) :: base
      LOGICAL, INTENT(IN) :: verbose
      REAL :: pow(nk, na)
      INTEGER :: i
      CHARACTER(len=512) :: output
      LOGICAL :: verbose2

      DO i = 1, 4
         IF (i == 1) THEN
            output = TRIM(base)//'_linear.dat'
            pow = pow_li
         ELSE IF (i == 2) THEN
            output = TRIM(base)//'_2h.dat'
            pow = pow_2h
         ELSE IF (i == 3) THEN
            output = TRIM(base)//'_1h.dat'
            pow = pow_1h
         ELSE IF (i == 4) THEN
            output = TRIM(base)//'_hm.dat'
            pow = pow_hm
         ELSE
            STOP 'WRITE_POWER_A_MULTIPLE: Error, something went wrong'
         END IF
         IF (i == 1) THEN
            verbose2 = verbose
         ELSE
            verbose2 = .FALSE.
         END IF
         CALL write_power_a(k, a, pow, nk, na, output, verbose2)
      END DO

   END SUBROUTINE write_power_a_multiple

   SUBROUTINE write_power_a(k, a, pow, nk, na, output, verbose)

      ! Write P(k,a) 
      INTEGER, INTENT(IN) :: nk
      INTEGER, INTENT(IN) :: na
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(IN) :: a(na)     
      REAL, INTENT(IN) :: pow(nk, na)
      CHARACTER(len=*), INTENT(IN) :: output     
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i, j

      ! Print to screen
      IF (verbose) THEN
         WRITE (*, *) 'WRITE_POWER_A: The first entry of the file is hashes - #####'
         WRITE (*, *) 'WRITE_POWER_A: The remainder of the first row are the scale factors - a'
         WRITE (*, *) 'WRITE_POWER_A: The remainder of the first column are the wave numbers - k'
         WRITE (*, *) 'WRITE_POWER_A: Each row then gives the power at that k and a'
         WRITE (*, *) 'WRITE_POWER_A: Output: ', TRIM(output)
      END IF

      ! Write out data to files
      OPEN (7, file=output)
      DO i = 0, nk
         IF (i == 0) THEN
            WRITE (7, fmt='(A20,40F20.10)') '#####', (a(j), j=1, na)
         ELSE
            WRITE (7, fmt='(F20.10,40E20.10)') k(i), (pow(i, j), j=1, na)
         END IF
      END DO
      CLOSE (7)

      ! Print to screen
      IF (verbose) THEN
         WRITE (*, *) 'WRITE_POWER_A: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE write_power_a

   SUBROUTINE read_power_a(k, a, pow, infile, na, verbose)

      ! Write P(k,a) 
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: a(:)     
      REAL, ALLOCATABLE, INTENT(OUT) :: pow(:, :)
      CHARACTER(len=*), INTENT(IN) :: infile
      INTEGER, INTENT(IN) :: na
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: ik, ia, u, nk
      CHARACTER(len=256) :: crap

      ! Print to screen
      IF (verbose) WRITE (*, *) 'READ_POWER_A: Input file: ', TRIM(infile)

      nk = file_length(infile)-1
      ALLOCATE(k(nk), a(na), pow(nk, na))

      ! Write data to files
      OPEN (newunit=u, file=infile, status='old')
      DO ik = 0, nk
         IF (ik == 0) THEN
            READ (u, *) crap, (a(ia), ia=1,na)
         ELSE
            READ (u, *) k(ik), (pow(ik, ia), ia=1,na)
         END IF
      END DO
      CLOSE (u)

      ! Print to screen
      IF (verbose) THEN
         WRITE (*, *) 'READ_POWER_A: Minimum k [h/Mpc]:', k(1)
         WRITE (*, *) 'READ_POWER_A: Maximum k [h/Mpc]:', k(nk)
         WRITE (*, *) 'READ_POWER_A: Number of points in k:', nk
         WRITE (*, *) 'READ_POWER_A: Minimum a:', a(1)
         WRITE (*, *) 'READ_POWER_A: Maximum a:', a(na)
         WRITE (*, *) 'READ_POWER_A: Number of points in a:', na
         WRITE (*, *) 'READ_POWER_A: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE read_power_a

   REAL FUNCTION P_PTish(k, a, sigv, Amp, R, cosm, approx)

      USE basic_operations
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: sigv
      REAL, INTENT(IN) :: R
      REAL, INTENT(IN) :: Amp
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, OPTIONAL, INTENT(IN) :: approx
      REAL :: P_dw, P_lin, P_loop
      INTEGER, PARAMETER :: flag = flag_matter

      P_dw = p_dwg(k, a, sigv, cosm)
      P_lin = plin(k, a, flag, cosm)
      IF (present_and_correct(approx)) THEN
         P_loop = P_SPT_approx(k, a, cosm)
      ELSE  
         P_loop = P_SPT(k, a, cosm)
      END IF
      P_loop = Amp*P_loop*exp(-(k*R)**2)
      P_PTish = P_dw*(P_lin+P_loop)/P_lin

   END FUNCTION P_PTish

   REAL FUNCTION simple_twohalo(k, Pk, rv)

      ! A simple two-halo term with haloes all of the same mass and are unbiased
      USE special_functions
      REAL, INTENT(IN) :: k  ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: Pk ! Power spectrum (usually linear)
      REAL, INTENT(IN) :: rv ! Virial radius [Mpc/h]

      simple_twohalo = Pk*wk_tophat(k*rv)**2

   END FUNCTION simple_twohalo

   REAL FUNCTION simple_onehalo(k, rv, M, rhobar)

      ! Does a simple halo-model calculation with haloes all of the same mass
      USE special_functions
      REAL, INTENT(IN) :: k      ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: rv     ! Virial radius [Mpc/h]
      REAL, INTENT(IN) :: M      ! Halo mass [Msun/h]
      REAL, INTENT(IN) :: rhobar ! Mean matter density [Msun/h / (Mpc/h)^3]

      simple_onehalo = 4.*pi*(k/twopi)**3*(wk_tophat(k*rv)**2)*M/rhobar

   END FUNCTION simple_onehalo

END MODULE HMx
