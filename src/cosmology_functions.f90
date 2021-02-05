MODULE cosmology_functions

   USE basic_operations
   USE array_operations
   USE string_operations
   USE interpolate
   USE constants
   USE file_info 
   USE table_integer

   IMPLICIT NONE

   PRIVATE

   ! Type
   PUBLIC :: cosmology

   ! Basic routines
   PUBLIC :: assign_cosmology
   PUBLIC :: init_cosmology
   PUBLIC :: print_cosmology
   PUBLIC :: assign_init_cosmology
   PUBLIC :: convert_cosmology

   ! Scale factor and z
   PUBLIC :: scale_factor_z
   PUBLIC :: redshift_r
   PUBLIC :: redshift_a
   PUBLIC :: scale_factor_r

   ! Background
   PUBLIC :: Hubble2
   PUBLIC :: Omega_m ! TODO: Retire
   PUBLIC :: Omega_c
   PUBLIC :: Omega_b
   PUBLIC :: Omega_r ! TODO: Retire
   PUBLIC :: Omega_g
   PUBLIC :: Omega_nu
   PUBLIC :: Omega_v
   PUBLIC :: Omega_w
   PUBLIC :: Omega
   PUBLIC :: w_de
   PUBLIC :: w_de_total
   PUBLIC :: w_eff

   ! Background without radiation
   PUBLIC :: Omega_m_norad

   ! Dark energy models
   PUBLIC :: iw_LCDM
   PUBLIC :: iw_wCDM
   PUBLIC :: iw_waCDM
   PUBLIC :: iw_QUICC
   PUBLIC :: iw_IDE1
   PUBLIC :: iw_IDE2
   PUBLIC :: iw_IDE3
   PUBLIC :: iw_BDE

   ! Modified gravity models
   PUBLIC :: img_none
   PUBLIC :: img_nDGP
   PUBLIC :: img_fR
   PUBLIC :: img_nDGP_lin
   PUBLIC :: img_fR_lin

   ! Modified gravity things
   PUBLIC :: fR_a

   ! Distances and times
   PUBLIC :: f_k
   PUBLIC :: fdash_k
   PUBLIC :: comoving_distance
   PUBLIC :: physical_distance
   PUBLIC :: comoving_particle_horizon
   PUBLIC :: physical_particle_horizon
   PUBLIC :: comoving_angular_distance
   PUBLIC :: physical_angular_distance
   PUBLIC :: luminosity_distance
   PUBLIC :: cosmic_time
   PUBLIC :: look_back_time
   PUBLIC :: conformal_time

   ! Densities
   PUBLIC :: comoving_critical_density
   PUBLIC :: comoving_matter_density
   PUBLIC :: physical_critical_density
   PUBLIC :: physical_matter_density

   ! Linear growth
   PUBLIC :: grow
   PUBLIC :: grow_Linder
   PUBLIC :: grow_CPT
   PUBLIC :: ungrow
   PUBLIC :: ungrow_approx
   PUBLIC :: ungrow_integral
   PUBLIC :: growth_rate
   PUBlIC :: growth_rate_Linder
   PUBLIC :: growth_index
   PUBLIC :: growth_index_Linder
   PUBLIC :: growth_index_default
   PUBLIC :: acc_growth

   ! Spherical collapse
   PUBLIC :: Dv_BryanNorman
   PUBLIC :: Dv_Mead
   PUBLIC :: Dv_Spherical
   PUBLIC :: dc_NakamuraSuto
   PUBLIC :: dc_Mead
   PUBLIC :: dc_Spherical

   ! Power spectrum
   PUBLIC :: Pk_Delta
   PUBLIC :: Delta_Pk
   PUBLIC :: Plin ! TODO: P_lin ?
   PUBLIC :: P_dwg
   PUBLIC :: P_smt
   PUBLIC :: calculate_Plin
   PUBLIC :: calculate_P_dwg
   PUBLIC :: calculate_P_smt
   
   ! sigma et al.
   PUBLIC :: sigma
   PUBLIC :: sigmaV
   PUBLIC :: neff
   PUBLIC :: ncur
   PUBLIC :: Lagrangian_mass
   PUBLIC :: Lagrangian_radius

   ! Normalisation
   PUBLIC :: As_norm
   PUBLIC :: sigma8_norm
   PUBLIC :: pval_norm

   ! Correlation function
   PUBLIC :: xi_lin

   ! Power spectrum normalisation
   PUBLIC :: normalise_power
   PUBLIC :: norm_sigma8
   PUBLIC :: norm_none
   PUBLIC :: norm_As
   PUBLIC :: norm_pval

   ! Flags for power type
   PUBLIC :: flag_matter
   PUBLIC :: flag_cold
   PUBLIC :: flag_ucold

   ! Transfer functions
   PUBLIC :: iTk_none
   PUBLIC :: iTk_EH
   PUBLIC :: iTk_nw
   PUBLIC :: iTk_DEFW
   PUBLIC :: iTk_CAMB
   PUBLIC :: iTk_external

   PUBLIC :: init_external_linear_power_tables

   ! CAMB interface
   PUBLIC :: get_CAMB_power
   PUBLIC :: CAMB_HALOFIT_Smith
   PUBLIC :: CAMB_HALOFIT_Bird
   PUBLIC :: CAMB_HALOFIT_Takahashi
   PUBLIC :: CAMB_HMcode2015
   PUBLIC :: CAMB_HMcode2016
   PUBLIC :: CAMB_HMcode2020
   PUBLIC :: CAMB_HMcode2020_feedback
   PUBLIC :: CAMB_HALOFIT_Smith_paper
   PUBLIC :: CAMB_HALOFIT_Bird_paper
   PUBLIC :: CAMB_HALOFIT_Takahashi_paper

   ! HALOFIT
   PUBLIC :: calculate_HALOFIT
   PUBLIC :: calculate_HALOFIT_a
   PUBLIC :: calculate_HALOFIT_ka
   PUBLIC :: HALOFIT_init
   PUBLIC :: HALOFIT_Smith_code
   PUBLIC :: HALOFIT_Smith_paper
   PUBLIC :: HALOFIT_Bird_code
   PUBLIC :: HALOFIT_Bird_paper
   PUBLIC :: HALOFIT_Bird_CAMB
   PUBLIC :: HALOFIT_Takahashi
   PUBLIC :: HALOFIT_CAMB
   PUBLIC :: HALOFIT_CLASS

   ! Perturbation theory
   PUBLIC :: P_SPT
   PUBLIC :: P_SPT_approx
   PUBLIC :: P_SPT_sum
   PUBLIC :: SPT_integrand

   ! Rescaling
   PUBLIC :: calculate_rescaling_parameters
   PUBLIC :: irescale_sigma
   PUBLIC :: irescale_power
   
   INTERFACE integrate_cosm
      MODULE PROCEDURE integrate_1_cosm
      MODULE PROCEDURE integrate_2_cosm
      MODULE PROCEDURE integrate_3_cosm
      MODULE PROCEDURE integrate_3_flag_cosm
      MODULE PROCEDURE integrate_4_flag_cosm
      MODULE PROCEDURE integrate_4_flag_2cosm
   END INTERFACE integrate_cosm

   ! Contains cosmological parameters that only need to be calculated once
   TYPE cosmology

      ! Primary parameters
      CHARACTER(len=256) :: name            ! Name for cosmological model
      REAL :: Om_m, Om_b, Om_v, Om_w        ! Densities
      REAL :: ns, nrun, nrunrun             ! Primordial spectrum
      REAL :: h, w, wa, m_wdm, YH           ! Cosmological parameters
      REAL :: a1, a2, nstar, ws, am, dm, wm ! Dark-energy parameters
      REAL :: T_CMB                         ! CMB temperature
      REAL :: neff, m_nu                    ! Neutrinos
      REAL :: H0rc, fr0, nfR                ! Modified gravity parameters
      REAL :: Om_m_pow, Om_b_pow, h_pow     ! Cosmological parameters used for P(k) if different from background
      REAL :: b0, b1, b2, b3, b4            ! BDE parameters
      REAL :: A_bump, k_bump, sigma_bump    ! Power-spectrum bump   
      REAL :: Theat                         ! AGN temperature
      REAL :: Lbox                          ! Box size [Mpc/h]
      INTEGER :: n_nu                       ! Number of massive neutrinos
      INTEGER :: bump                       ! Type of bump to add to P(k)
      INTEGER :: norm_method                ! Power normalisation scheme
      INTEGER :: iw                         ! Switch for dark-energy type
      INTEGER :: img                        ! Switch for modified gravity
      INTEGER :: iTk                        ! Switch for transfer function type
      INTEGER :: iTc                        ! Switch for cold transfer function type
      LOGICAL :: box                        ! Constrain the calculation to take place in a box?    
      LOGICAL :: warm                       ! Is DM warm?
      LOGICAL :: power_Omegas               ! Are the Omegas for the background different from those for the perturbations?
      LOGICAL :: derive_gas_numbers         ! Should mu_e and mu_p be derived or not?

      ! Variables that might be primary or dervied
      REAL :: kpiv, As, kval, pval, sig8 ! Power spectrum normalisation   
      REAL :: mue, mup                   ! Gas parameters

      ! Derived parameters
      REAL :: A                          ! Power spectrum amplitude
      REAL :: Om, Om_k, Om_c, Om_g, Om_r ! Derived Omegas
      REAL :: Om_nu, f_nu, a_nu          ! Neutrinos
      REAL :: Om_nu_rad, omega_nu, T_nu  ! Neutrinos
      REAL :: omega_m, omega_b, omega_c  ! Physical densities
      REAL :: k                          ! Curvature [(Mpc/h)^-2]
      REAL :: z_CMB                      ! Distance to the CMB
      REAL :: Om_c_pow                   ! Cosmological parameters used for P(k) if different from background
      REAL :: age, horizon               ! Derived distance/time
      REAL :: YHe                        ! Derived thermal parameters
      REAL :: Om_ws, astar, a1n, a2n     ! Derived DE parameters
      REAL :: gnorm                      ! Growth-factor normalisation
      REAL :: kbox                       ! Wavenumber of box mode
      LOGICAL :: scale_dependent_growth  ! Is the linear growth scale dependent in this cosmology?
      
      ! Look-up tables that are filled during a calculation
      REAL, ALLOCATABLE :: k_Plin(:), a_Plin(:), Plin_array(:, :) ! Arrays for input linear P(k) TODO: Remove
      INTEGER :: nk_plin, na_plin ! Number of array entries TODO: Remove
      TYPE(interpolator1D) :: grow, grate, agrow, dc, Dv, dist, time, Xde, rspt ! 1D interpolators
      TYPE(interpolator2D) :: sigma, plin, wiggle, Tk_matter, Tk_cold ! 2D interpolators
      LOGICAL :: analytical_Tk                                                          
      LOGICAL :: has_distance, has_growth, has_sigma, has_spherical, has_power ! Interpolator status
      LOGICAL :: has_wiggle, has_SPT, has_time, has_Xde ! Interpolator status
      LOGICAL :: is_init, is_normalised ! Flags to check if things have been done 

      ! CAMB
      CHARACTER(len=256) :: CAMB_exe, CAMB_temp_dir

      ! Verbose
      LOGICAL :: verbose

      ! Error handling
      INTEGER :: status

   END TYPE cosmology

   ! Global parameters
   REAL, PARAMETER :: acc_cosm = 1e-4        ! Global accuacy for the cosmological integrations
   INTEGER, PARAMETER :: ncosmo_large = 1000 ! Needs to be larger than the total number of defined cosmologies TODO: Remove

   ! Writing to screen parameters
   REAL, PARAMETER :: small_curve = 1e-5 ! Used to decide if writing curvature to screen to avoid silly numbers

   ! Equation of state parameters for cosmological fluids
   REAL, PARAMETER :: w_c = 0.    ! CDM
   REAL, PARAMETER :: w_b = 0.    ! Baryons
   REAL, PARAMETER :: w_g = 1./3. ! Photons
   REAL, PARAMETER :: w_v = -1.   ! Vacuum energy

   ! Neutrino energy-density evolution methods
   INTEGER, PARAMETER :: neutrino_basic = 1                 ! Basic approximation that ensures relativistic and non-relativisit limits
   INTEGER, PARAMETER :: neutrino_Komatsu = 2               ! Accurate approximation from WMAP7 paper
   INTEGER, PARAMETER :: neutrino_method = neutrino_Komatsu ! Select neutrino method

   ! Parameters from WMAP7 Komatsu et al. 2011 paper
   REAL, PARAMETER :: A_Komatsu = 0.3173 ! Acutally 180*zeta(3)/7*pi^4
   REAL, PARAMETER :: p_Komatsu = 1.83   ! Fitted

   ! Neutrinos
   REAL, PARAMETER :: f_nu_limit = 0.5 ! If f_nu is larger than this then stop the calculation
   REAL, PARAMETER :: a_nu_limit = 0.2 ! If neutrinos are too hot then stop the calculation

   ! Dark energy density
   INTEGER, PARAMETER :: iw_LCDM = 1  ! Vacuum dark energy
   INTEGER, PARAMETER :: iw_QUICC = 2 ! QUICC dark energy
   INTEGER, PARAMETER :: iw_waCDM = 3 ! w(a) dark energy
   INTEGER, PARAMETER :: iw_wCDM = 4  ! Constant w dark energy
   INTEGER, PARAMETER :: iw_IDE1 = 5  ! Intermediate dark energy model 1
   INTEGER, PARAMETER :: iw_IDE2 = 6  ! Intermediate dark energy model 2
   INTEGER, PARAMETER :: iw_IDE3 = 7  ! Intermediate dark energy model 3
   INTEGER, PARAMETER :: iw_BDE = 8   ! Bound dark energy (1812.01133)  

   ! Dark energy integration and interpolation
   REAL, PARAMETER :: amin_Xde = 1e-4                ! Minimum scale factor for direction integration to get Xde
   REAL, PARAMETER :: amax_Xde = 1.                  ! Maximum scale factor for direction integration to get Xde
   INTEGER, PARAMETER :: n_Xde = 128                 ! Number of points for Xde interpolation
   LOGICAL, PARAMETER :: tabulate_Xde = .TRUE.       ! Tabulate Xde for interpolation
   REAL, PARAMETER :: acc_integration_Xde = acc_cosm ! Accuracy for direct integration of dark energy density
   INTEGER, PARAMETER :: iorder_integration_Xde = 3  ! Polynomial order for time integration
   INTEGER, PARAMETER :: iorder_interp_Xde = 3       ! Polynomial order for time interpolation
   INTEGER, PARAMETER :: iextrap_Xde = iextrap_lin   ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_Xde = .TRUE.          ! Storage in interpolator

   ! Modified gravity
   INTEGER, PARAMETER :: img_none = 0     ! Standard gravity
   INTEGER, PARAMETER :: img_nDGP = 1     ! Normal-branch DGP gravity with a LCDM background
   INTEGER, PARAMETER :: img_fR = 2       ! f(R) gravity with a LCDM background
   INTEGER, PARAMETER :: img_nDGP_lin = 3 ! Linearised nDGP (only affects spherical model and HMcode)
   INTEGER, PARAMETER :: img_fR_lin = 4   ! Linearised f(R) (only affects spherical model and HMcode)

   ! Distance
   REAL, PARAMETER :: amin_distance = 1e-4               ! Minimum scale factor in look-up table
   REAL, PARAMETER :: amax_distance = 1.                 ! Maximum scale factor in look-up table
   INTEGER, PARAMETER :: n_distance = 128                ! Number of scale factor entries in look-up table
   REAL, PARAMETER :: atay_distance = 1e-5               ! Below this do a Taylor expansion to avoid divergence
   INTEGER, PARAMETER :: iorder_integration_distance = 3 ! Polynomial order for distance integration
   INTEGER, PARAMETER :: iorder_interp_distance = 3      ! Polynomial order for distance interpolation
   INTEGER, PARAMETER :: iextrap_distance = iextrap_std  ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_distance = .TRUE.         ! Pre-calculate interpolation coefficients?

   ! Time
   REAL, PARAMETER :: amin_time = 1e-4               ! Minimum scale factor in look-up table
   REAL, PARAMETER :: amax_time = 1.                 ! Maximum scale factor in look-up table
   INTEGER, PARAMETER :: n_time = 128                ! Number of scale factor entries in look-up table
   REAL, PARAMETER :: atay_time = 1e-5               ! Below this do a Taylor expansion to avoid divergence
   INTEGER, PARAMETER :: iorder_integration_time = 3 ! Polynomial order for time integration
   INTEGER, PARAMETER :: iorder_interp_time = 3      ! Polynomial order for time interpolation
   INTEGER, PARAMETER :: iextrap_time = iextrap_std  ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_time = .TRUE.         ! Pre-calculate interpolation coefficients?

   ! Linear power
   REAL, PARAMETER :: kmin_abs_plin = 0.       ! Power below this wavenumber is set to zero [h/Mpc]
   REAL, PARAMETER :: kmax_abs_plin = 1e8      ! Power above this wavenumber is set to zero [h/Mpc]
   LOGICAL, PARAMETER :: plin_extrap = .FALSE. ! Extrapolate high-k power assuming P(k) ~ ln(k)^2 k^(n-3)?
   INTEGER, PARAMETER :: iTk_none = 0          ! Pure power-law spectrum
   INTEGER, PARAMETER :: iTk_EH = 1            ! Eisenstein & Hu linear spectrum
   INTEGER, PARAMETER :: iTk_CAMB = 2          ! CAMB linear spectrum
   INTEGER, PARAMETER :: iTk_DEFW = 3          ! DEFW linear spectrum
   INTEGER, PARAMETER :: iTk_external = 4      ! External linear spectrum
   INTEGER, PARAMETER :: iTk_nw = 5            ! No-wiggle Eisenstein & Hu linear spectrum
   INTEGER, PARAMETER :: norm_none = 0         ! Power spectrum does not need to be normalised
   INTEGER, PARAMETER :: norm_sigma8 = 1       ! Normalise power spectrum via sigma8 value
   INTEGER, PARAMETER :: norm_pval = 2         ! Normalise power spectrum via specifying a value at a k 
   INTEGER, PARAMETER :: norm_As = 3           ! Normalise power spectrum vis As value as in CAMB   
   INTEGER, PARAMETER :: flag_matter = 1       ! Flag to get the total matter power spectrum
   INTEGER, PARAMETER :: flag_cold = 2         ! Flag to get the cold (CDM+baryons) power spectrum with 1+delta = rho_cold/mean_rho_matter
   INTEGER, PARAMETER :: flag_ucold = 3        ! Flag to get the cold (CDM+baryons) power spectrum with 1+delta = rho_cold/mean_rho_cold

   ! Linear power interpolation
   ! NOTE: In CAMB HMcode different accuracy parameters are used, for example linear interpolation
   LOGICAL, PARAMETER :: interp_all_power = .FALSE.      ! Create interpolators for all linear power spectra, even analytical ones
   REAL, PARAMETER :: kmin_plin = 1e-3                   ! Minimum wavenumber used [h/Mpc]
   REAL, PARAMETER :: kmax_plin = 1e2                    ! Maximum wavenumber used [h/Mpc]
   INTEGER, PARAMETER :: nk_plin = 128                   ! Number of k points to use
   REAL, PARAMETER :: amin_plin = 0.1                    ! Minimum a value for Pk growth if scale dependent
   REAL, PARAMETER :: amax_plin = 1.0                    ! Maximum a value for Pk growth if scale dependent
   INTEGER, PARAMETER :: na_plin = 16                    ! Number of a values if growth is scale dependent
   INTEGER, PARAMETER :: iorder_interp_plin = 3          ! Polynomial order
   INTEGER, PARAMETER :: ifind_interp_plin = ifind_split ! Finding scheme in table (only linear if rebinning)
   INTEGER, PARAMETER :: iinterp_plin = iinterp_Lagrange ! Method for interpolation polynomials
   INTEGER, PARAMETER :: iextrap_plin = iextrap_lin      ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_plin = .TRUE.             ! Pre-calculate interpolation coefficients

   ! CAMB interface
   ! NOTE: rebin_CAMB makes a big difference in comparisons with CAMB HMcode; in CAMB P(k) *is* generally rebinned
   REAL, PARAMETER :: pk_min_CAMB = 1e-10                      ! Minimum value of power at low k (remove k with less than this) 
   REAL, PARAMETER :: nmax_CAMB = 2.                           ! How many times more to go than kmax due to inaccuracy near k limit
   CHARACTER(len=256), PARAMETER :: de_CAMB = 'ppf'            ! CAMB dark-energy prescription, either 'ppf' or 'fluid'
   LOGICAL, PARAMETER :: rebin_CAMB = .FALSE.                  ! Should we rebin CAMB or just use default k spacing?
   INTEGER, PARAMETER :: iorder_rebin_CAMB = 3                 ! Polynomial order for interpolation if rebinning P(k)
   INTEGER, PARAMETER :: ifind_rebin_CAMB = ifind_split        ! Finding scheme for interpolation if rebinning P(k) (*definitely* not linear)
   INTEGER, PARAMETER :: iinterp_rebin_CAMB = iinterp_Lagrange ! Interpolation scheme if rebinning P(k)

   ! Cold transfer function methods
   INTEGER, PARAMETER :: iTc_none = 0  ! Assume cold power is indentical to matter
   INTEGER, PARAMETER :: iTc_total = 1 ! Assume neutrinos are completely hot
   INTEGER, PARAMETER :: iTc_EH = 2    ! Eisenstein & Hu approximation
   INTEGER, PARAMETER :: iTc_CAMB = 3  ! Taken from CAMB

   ! EH cold transfer function
   LOGICAL, PARAMETER :: Tk_cold_EdS_growth = .FALSE. ! Use the (incorrect) EdS growth function in the fitting function

   ! CAMB cold transfer function
   INTEGER, PARAMETER :: iextrap_Tk = iextrap_lin ! Extrapolation scheme for cold interpolation
   INTEGER, PARAMETER :: iorder_interp_Tk = 3     ! Order for cold interpolatin
   LOGICAL, PARAMETER :: store_Tk = .TRUE.        ! Storage for cold interpolation

   ! Wiggle smoothing extraction methods
   INTEGER, PARAMETER :: dewiggle_tophat = 1   ! Top-hat smoothing
   INTEGER, PARAMETER :: dewiggle_Gaussian = 2 ! Gaussian smoothing

   ! Linear power spectrum smoothing methods  
   ! TODO: Setting scale_grow_wiggle to .FALSE. may save time for massive-nu models
   REAL, PARAMETER :: wiggle_dx = 0.20                     ! Smoothing half-width if using top-hat smoothing
   REAL, PARAMETER :: wiggle_sigma = 0.25                  ! Smoothing width if using Gaussian smoothing
   INTEGER, PARAMETER :: wiggle_smooth = dewiggle_Gaussian ! Type of smoothing to use
   LOGICAL, PARAMETER :: divide_by_nowiggle = .TRUE.       ! Should we reduce dynamic range with EH no-wiggle?
   REAL, PARAMETER :: knorm_nowiggle = 0.03                ! Wavenumber at which to force linear and nowiggle to be identical [Mpc/h]
   LOGICAL, PARAMETER :: scale_grow_wiggle = .TRUE.        ! Treat the wiggle as being different at different 'a'

   ! Wiggle extraction and interpolation
   REAL, PARAMETER :: kmin_wiggle = 5e-3               ! Minimum wavenumber to calulate wiggle [Mpc/h]
   REAL, PARAMETER :: kmax_wiggle = 5.                 ! Maximum wavenumber to calulate wiggle [Mpc/h]
   INTEGER, PARAMETER :: nk_wiggle = 512               ! Number of k points to store wiggle
   LOGICAL, PARAMETER :: store_wiggle = .TRUE.         ! Pre-calculate interpolation coefficients 
   INTEGER, PARAMETER :: iorder_interp_wiggle = 3      ! Order for wiggle interpolator
   INTEGER, PARAMETER :: iextrap_wiggle = iextrap_zero ! Should be zeros because interpolator stores only wiggle

   ! Correlation function
   ! TODO: This works very poorly
   INTEGER, PARAMETER :: method_xi = 1         ! Method for xi integration
   INTEGER, PARAMETER :: iorder_xi = 3         ! Polynomial order for xi(r) integration
   INTEGER, PARAMETER :: min_humps_xi = 5      ! Minimum number of humps fox xi(r) integration if using humps integration
   INTEGER, PARAMETER :: max_humps_xi = 100000 ! Maximum number of humps fox xi(r) integration if using humps integration
   REAL, PARAMETER :: rsplit_xi = 10.          ! Value of r to split alpha values for conventional integration
   REAL, PARAMETER :: alpha_lo_xi = 2.         ! Low r value of alpha for conventional integration
   REAL, PARAMETER :: alpha_hi_xi = 1.5        ! High r value of alpha for conventional integration

   ! Growth ODE
   REAL, PARAMETER :: aini_growth = 1e-4                     ! Starting value for growth integratiton (should start | Omega_m(a)=1)
   REAL, PARAMETER :: afin_growth = 2.                       ! Finishing value for growth integratiton (CARE: changed from a=1 to a=2)
   REAL, PARAMETER :: acc_ODE_growth = acc_cosm              ! Accuracy parameter for growth ODE solving
   INTEGER, PARAMETER :: imeth_ODE_growth = 3                ! Method for solving growth ODE
   INTEGER, PARAMETER :: iorder_ODE_interpolation_growth = 3 ! Polynomial order for growth interpolation for ODE solution
   INTEGER, PARAMETER :: ifind_ODE_interpolation_growth = 3  ! Finding scheme for growth interpolation for ODE solution
   INTEGER, PARAMETER :: imeth_ODE_interpolation_growth = 2  ! Method for growth interpolation for ODE solution
   LOGICAL, PARAMETER :: cold_growth = .FALSE.               ! Should smooth neutrinos be accounted for in growth calculations?
   LOGICAL, PARAMETER :: EDE_growth_ics = .TRUE.             ! Should we try to account for EDE in growth initial conditions?

   ! Growth integrals (approximate)
   REAL, PARAMETER :: acc_integral_grow = acc_cosm ! Accuracy parameter for growth integral
   INTEGER, PARAMETER :: iorder_integral_grow = 3  ! Polynomial order for growth integral
   REAL, PARAMETER :: aeps_integral_grow = 1e-3    ! Minimum scale factor below which to approximate integrand

   ! Growth interpolation
   REAL, PARAMETER :: amin_growth = 1e-3            ! Minimum value to store
   REAL, PARAMETER :: amax_growth = afin_growth     ! Maximum value to store
   INTEGER, PARAMETER :: na_growth = 128            ! Number of entries for interpolation tables
   INTEGER, PARAMETER :: iorder_interp_grow = 3     ! Polynomial order for growth interpolation
   INTEGER, PARAMETER :: iextrap_grow = iextrap_lin ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_grow = .TRUE.        ! Pre-calculate interpolation coefficients?

   ! Growth rate interpolation
   INTEGER, PARAMETER :: iorder_interp_rate = 3     ! Polynomial order for growth rate interpolation for ODE solution
   INTEGER, PARAMETER :: iextrap_rate = iextrap_lin ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_rate = .TRUE.        ! Pre-calculate interpolation coefficients?

   ! Growth rate index
   REAL, PARAMETER :: growth_index_default = 6./11. ! Default index value (perturbation theory for LCDM)
   REAL, PARAMETER :: growth_index_limit = 0.01     ! 1-Omega_m value below which to switch to the default

   ! Accumualted growth integration and interpolation
   INTEGER, PARAMETER :: iorder_integration_agrow = 3 ! Polynomial order for accumulated growth integration
   INTEGER, PARAMETER :: iorder_interp_agrow = 3      ! Polynomial order for interpolation of accumulated growth
   INTEGER, PARAMETER :: iextrap_agrow = iextrap_lin  ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_agrow = .TRUE.         ! Pre-calculate interpolation coefficients?

   ! sigma(R) integration
   REAL, PARAMETER :: alpha_sigma = 3.     ! Exponent to increase speed (1 is terrible, 2, 3, 4 all okay)
   REAL, PARAMETER :: acc_sigma = acc_cosm ! Accuracy parameter for sigma(R) integration
   INTEGER, PARAMETER :: iorder_sigma = 3  ! Polynomial order for sigma(R) integration
   
   ! sigma(R) tabulation and interpolation
   REAL, PARAMETER :: rmin_sigma = 1e-4                   ! Minimum r value (NB. sigma(R) needs to be power-law below) [Mpc/h]
   REAL, PARAMETER :: rmax_sigma = 1e3                    ! Maximum r value (NB. sigma(R) needs to be power-law above) [Mpc/h]
   INTEGER, PARAMETER :: nr_sigma = 128                   ! Number of r entries for sigma(R) tables
   REAL, PARAMETER :: amin_sigma = amin_plin              ! Minimum a value for sigma(R,a) tables when growth is scale dependent
   REAL, PARAMETER :: amax_sigma = amax_plin              ! Maximum a value for sigma(R,a) tables when growth is scale dependent
   INTEGER, PARAMETER :: na_sigma = 16                    ! Number of a values for sigma(R,a) tables
   INTEGER, PARAMETER :: iorder_interp_sigma = 3          ! Polynomial order for sigma(R) interpolation 
   INTEGER, PARAMETER :: ifind_interp_sigma = ifind_split ! Finding scheme for sigma(R) interpolation (changing to linear not speedy)
   INTEGER, PARAMETER :: sigma_store = flag_ucold         ! Which version of sigma should be tabulated (0 for none)
   INTEGER, PARAMETER :: iextrap_sigma = iextrap_lin      ! Extrapolation for sigma(R) interpolator
   LOGICAL, PARAMETER :: store_sigma = .TRUE.             ! Pre-calculate interpolation coefficients?

   ! sigma_v(R) integration
   REAL, PARAMETER :: alpha_sigmaV = 3.     ! Exponent to increase integration speed
   REAL, PARAMETER :: acc_sigmaV = acc_cosm ! Accuracy parameter for sigma(R) integration
   INTEGER, PARAMETER :: iorder_sigmaV = 3  ! Polynomial order for sigmaV(R) integration

   ! neff integration
   REAL, PARAMETER :: alpha_dsigma = 2.     ! Exponent to increase integration speed (3 seemed to give inaccuracies; no idea why)
   REAL, PARAMETER :: acc_dsigma = acc_cosm ! Accuracy parameter for neff(R) integration
   INTEGER, PARAMETER :: iorder_dsigma = 3  ! Polynomial order for neff(R) integration

   ! ncur integration
   REAL, PARAMETER :: alpha_ddsigma = 2.     ! Exponent to increase integration speed
   REAL, PARAMETER :: acc_ddsigma = acc_cosm ! Accuracy parameter for ncur(R) integration
   INTEGER, PARAMETER :: iorder_ddsigma = 3  ! Polynomial order for ncur(R) integration

   ! Spherical collapse
   INTEGER, PARAMETER :: imeth_ODE_spherical = 3 ! Method for spherical collapse ODE solving
   REAL, PARAMETER :: amax_spherical = 2.        ! Maximum scale factor to consider
   REAL, PARAMETER :: dmin_spherical = 1e-7      ! Minimum starting value for perturbation
   REAL, PARAMETER :: dmax_spherical = 1e-3      ! Maximum starting value for perturbation
   INTEGER, PARAMETER :: m_spherical = 128       ! Number of collapse scale-factors to try to calculate
   INTEGER, PARAMETER :: n_spherical = 100000    ! Number of points for ODE calculations
   REAL, PARAMETER :: dinf_spherical = 1e8       ! Value considered to be 'infinite' for the perturbation

   ! delta_c
   INTEGER, PARAMETER :: iorder_interp_dc = 3     ! Polynomial order for delta_c interpolation
   INTEGER, PARAMETER :: ifind_interp_dc = 3      ! Finding scheme for delta_c interpolation
   INTEGER, PARAMETER :: imeth_interp_dc = 2      ! Method for delta_c interpolation
   INTEGER, PARAMETER :: iextrap_dc = iextrap_std ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_dc = .TRUE.        ! Pre-calculate interpolation coefficients?
   LOGICAL, PARAMETER :: cold_Nakamura = .FALSE.  ! Use cold Omega_m in Nakamura & Suto formula?

   ! Delta_v
   INTEGER, PARAMETER :: iorder_interp_Dv = 3     ! Polynomial order for Delta_v interpolation
   INTEGER, PARAMETER :: ifind_interp_Dv = 3      ! Finding scheme for Delta_v interpolation
   INTEGER, PARAMETER :: imeth_interp_Dv = 2      ! Method for Delta_v interpolation
   INTEGER, PARAMETER :: iextrap_Dv = iextrap_std ! Extrapolation scheme
   LOGICAL, PARAMETER :: store_Dv = .TRUE.        ! Pre-calculate interpolation coefficients?
   LOGICAL, PARAMETER :: cold_Bryan = .FALSE.     ! Use cold Omega_m in Bryan & Norman formula?

   ! HALOFIT versions
   INTEGER, PARAMETER :: HALOFIT_Smith_code = 1  ! Smith et al. (2003; https://www.roe.ac.uk/~jap/haloes/) with numbers as in online code
   INTEGER, PARAMETER :: HALOFIT_Bird_code = 2   ! Bird et al. (2012; https://arxiv.org/abs/1109.4416) with numbers as in online code
   INTEGER, PARAMETER :: HALOFIT_Takahashi = 3   ! Takahashi et al. (2012; https://arxiv.org/abs/1208.2701)
   INTEGER, PARAMETER :: HALOFIT_CAMB = 4        ! Version as used in CAMB (2020) an unpublished hybrid of Takahashi and Bird
   INTEGER, PARAMETER :: HALOFIT_CLASS = 5       ! Version as used in CLASS (2020) ?
   INTEGER, PARAMETER :: HALOFIT_Smith_paper = 6 ! Smith et al. (2003; https://arxiv.org/abs/astro-ph/0207664) with numbers as in paper
   INTEGER, PARAMETER :: HALOFIT_Bird_paper = 7  ! Bird et al. (2012; https://arxiv.org/abs/1109.4416) with numbers as in paper
   INTEGER, PARAMETER :: HALOFIT_Bird_CAMB = 8   ! Bird et al. (2012; https://arxiv.org/abs/1109.4416) with numbers as in CAMB

   ! HALOFIT parameters
   REAL, PARAMETER :: HALOFIT_acc = 1e-3  ! Accuracy for HALOFIT calculation of knl, neff, ncur (1e-3 is in public code and CAMB)
   REAL, PARAMETER :: HALOFIT_logr1 = -3. ! log10(R / Mpc/h) minimum filter size to start search (note this is -2 in CAMB)
   REAL, PARAMETER :: HALOFIT_logr2 = 3.5 ! log10(R / Mpc/h) maximum filter size to start search

   ! CAMB non-linear numbering schemes
   INTEGER, PARAMETER :: CAMB_HALOFIT_Smith = 1            ! Smith et al. (2003; https://www.roe.ac.uk/~jap/haloes/)
   INTEGER, PARAMETER :: CAMB_HALOFIT_Bird = 2             ! Bird et al. (2012; https://arxiv.org/abs/1109.4416)
   INTEGER, PARAMETER :: CAMB_HALOFIT_Takahashi = 4        ! Takahashi et al. (2012; https://arxiv.org/abs/1208.2701)
   INTEGER, PARAMETER :: CAMB_HALOFIT_Smith_paper = 11     ! Smith et al. (2003; https://www.roe.ac.uk/~jap/haloes/)
   INTEGER, PARAMETER :: CAMB_HALOFIT_Bird_paper = 12      ! Bird et al. (2012; https://arxiv.org/abs/1109.4416)
   INTEGER, PARAMETER :: CAMB_HALOFIT_Takahashi_paper = 13 ! Takahashi et al. (2012; https://arxiv.org/abs/1208.2701)
   INTEGER, PARAMETER :: CAMB_HMcode2015 = 8               ! Mead et al. (2015; https://arxiv.org/abs/1505.07833)
   INTEGER, PARAMETER :: CAMB_HMcode2016 = 5               ! Mead et al. (2016; https://arxiv.org/abs/1602.02154)
   INTEGER, PARAMETER :: CAMB_HMcode2020 = 9               ! Mead et al. (2020; https://arxiv.org/abs/2009.01858)
   INTEGER, PARAMETER :: CAMB_HMcode2020_feedback = 10     ! Mead et al. (2020; https://arxiv.org/abs/2009.01858)

   ! Standard perturbation theory (SPT)
   REAL, PARAMETER :: kmin_integrate_SPT = 1e-4    ! Minimum wavenumber for integration [h/Mpc]: -9.2 in Komatsu code
   REAL, PARAMETER :: kmax_integrate_SPT = 20.     ! Maximum wavenumber for integration [h/Mpc]: e^3 ~ 20. in Komatsu code
   REAL, PARAMETER :: eps_multiple_SPT = 0.17      ! Multiple of wavenumber to take care in P_22 around pole: 0.17 in Komatsu code
   REAL, PARAMETER :: q_on_k_limit_F3 = 100.       ! Limit of q/k for Taylor expansion in F3 kernel: 1e2 in Komatsu code
   REAL, PARAMETER :: acc_integrate_SPT = acc_cosm ! Accuracy for integration
   INTEGER, PARAMETER :: iorder_integrate_SPT = 3  ! Order for integration
   LOGICAL, PARAMETER :: interp_SPT = .TRUE.       ! Use interpolator for SPT?
   REAL, PARAMETER :: kmin_interpolate_SPT = 5e-3  ! Minimum wavenumber for interpolator [h/Mpc]
   REAL, PARAMETER :: kmax_interpolate_SPT = 0.5   ! Maximum wavenumber for interpolator [h/Mpc]
   INTEGER, PARAMETER :: nk_interpolate_SPT = 64   ! Number of wavenumber points in interpolator
   INTEGER, PARAMETER :: iorder_interp_SPT = 3     ! Order for interpolation
   INTEGER, PARAMETER :: iextrap_SPT = iextrap_std ! Extrapolation scheme for interpolation
   LOGICAL, PARAMETER :: store_interp_SPT = .TRUE. ! Store coefficients for interpolation?

   ! Rescaling
   INTEGER, PARAMETER :: irescale_sigma = 1          ! Rescale by matching sigma(R)
   INTEGER, PARAMETER :: irescale_power = 2          ! Rescale by matching P(k)
   INTEGER, PARAMETER :: flag_rescale = flag_ucold   ! Type of power to match by rescaling
   REAL, PARAMETER :: acc_rescale = acc_cosm         ! Rescaling accuracy
   INTEGER, PARAMETER :: iorder_rescale_integral = 3 ! Order for integration
   REAL, PARAMETER :: smin_rescale = 0.33            ! Minimum s (R -> sR)
   REAL, PARAMETER :: smax_rescale = 3.00            ! Maximum s (R -> sR)
   INTEGER, PARAMETER :: ns_rescale = 268            ! Number of points in s
   !INTEGER, PARAMETER :: ns_rescale = nint(1+100.*(smax_rescale-smin_rescale)) ! Number of points in s

   ! Cosmic Emu (for finding h)
   ! TODO: Remove
   CHARACTER(len=256), PARAMETER :: params_CosmicEmu = 'emu_params.txt'
   CHARACTER(len=256), PARAMETER :: output_CosmicEmu = 'emu_power.dat'
   CHARACTER(len=256), PARAMETER :: exe_CosmicEmu = '/Users/Mead/Physics/CosmicEmu/emu.exe'
   INTEGER, PARAMETER :: nh_CosmicEmu = 10 ! Length of header
   INTEGER, PARAMETER :: hl_CosmicEmu = 8  ! Line that h in on

   ! General cosmological integrations
   INTEGER, PARAMETER :: jmin_integration = 5  ! Minimum number of points: 2^(j-1)
   INTEGER, PARAMETER :: jmax_integration = 30 ! Maximum number of points: 2^(j-1) TODO: Could lower to make time-out faster
   INTEGER, PARAMETER :: jmax_ode = 20 ! Maximum number of ODE attempts
   INTEGER, PARAMETER :: nmin_ode = 10 ! Minimum number of points for ODE to solve

CONTAINS

   SUBROUTINE assign_cosmology(icosmo, cosm, verbose)

      ! Assigns the 'primary' cosmological parameters (primary according to my definition)
      ! This routine *only* assigns parameters, it does and should not do *any* calculations
      INTEGER, INTENT(INOUT) :: icosmo
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i
      REAL :: Xe, Xi
      REAL :: wm, wb, wc

      ! Names of pre-defined cosmologies
      ! TODO: This is lazy, there must be a better way
      INTEGER, PARAMETER :: ncosmo = ncosmo_large
      CHARACTER(len=256), ALLOCATABLE :: names(:)

      ! Allocate array for names (warning from gfortran 10)
      ALLOCATE(names(ncosmo))
      names = ''

      names(1)  = 'Boring'
      names(2)  = 'WMAP7 (cosmo-OWLS)'
      names(3)  = 'Planck 2013 (cosmo-OWLS/BAHAMAS)'
      names(4)  = 'WMAP9 (BAHAMAS; EH)'
      names(5)  = 'Open'
      names(6)  = 'Einstein-de Sitter'
      names(7)  = 'IDE I'
      names(8)  = 'IDE II'
      names(9)  = 'IDE III'
      names(10) = 'Random neutrinoless Mira Titan cosmology'
      names(11) = 'nDGP - Strong'
      names(12) = 'nDGP - Medium'
      names(13) = 'nDGP - Weak'
      names(14) = 'WDM'
      names(15) = 'TCDM'
      names(16) = 'Boring: w = -0.7'
      names(17) = 'Boring: w = -1.3'
      names(18) = 'Boring: w = -1.0; wa =  0.5'
      names(19) = 'Boring: w = -1.0; wa = -0.5'
      names(20) = 'Boring: w = -0.7; wa = -1.5'
      names(21) = 'Boring: w = -1.3; wa =  0.5'
      names(22) = 'IDE3'
      names(23) = 'IDE10'
      names(24) = 'Random Mira Titan cosmology'
      names(25) = 'Random Fanken Emu cosmology'
      names(26) = 'Boring: CAMB linear spectrum'
      names(27) = 'Illustris TNG 75'
      names(28) = 'Boring: Finite box'
      names(29) = 'Boring LCDM: Open; z=1 normalisation'
      names(30) = 'Boring LCDM: w = -0.7; z=1 normalisation'
      names(31) = 'Boring LCDM: w = -1.3; z=1 normalisation'
      names(32) = 'Boring LCDM: w = -1.0; wa =  0.5; z=1 normalisation'
      names(33) = 'Boring LCDM: w = -1.0; wa = -0.5; z=1 normalisation'
      names(34) = 'Boring LCDM: w = -0.7; wa = -1.5; z=1 normalisation'
      names(35) = 'Boring LCDM: w = -1.3; wa =  0.5; z=1 normalisation'
      names(36) = 'Einsten-de Sitter LCDM; z=1 normalisation'
      names(37) = 'Multidark: WMAP 5'
      names(38) = 'Random Cosmic Emu cosmology'
      names(39) = 'Random LCDM cosmology (EH linear)'
      names(40) = 'Random nu-w(a)CDM cosmology'
      names(41) = 'SCDM with some CDM replaced by 4eV neutrinos'
      names(42) = 'Planck 2018 (no neutrinos)'
      names(43) = 'Multidark: WMAP 5 with lower sigma_8'
      names(44) = 'Multidark: WMAP 5 with Eisenstein & Hu'
      names(45) = 'Random w(a)CDM cosmology'
      names(46) = 'Random wCDM cosmology'
      names(47) = 'Random LCDM cosmology'
      names(48) = 'Random nu-LCDM cosmology'
      names(49) = 'CFHTLenS (Heymans et al. 2013)'
      names(50) = 'CFHTLenS (Kilbinger et al. 2013)'
      names(51) = 'CAMB difference cosmology 1'
      names(52) = 'CAMB difference cosmology 2'
      names(53) = 'CAMB difference cosmology 3'
      names(54) = 'CAMB difference cosmology 4'
      names(55) = 'CAMB difference cosmology 5'
      names(56) = 'Planck 2018'
      names(57) = 'Bound dark energy'
      names(58) = 'Extreme bound dark energy'
      names(59) = 'Power bump'
      names(60) = 'Boring but with some CDM replaced with 0.3eV neutrinos'
      names(61) = 'WMAP9 (BAHAMAS AGN 7.6)'
      names(62) = 'WMAP9 (BAHAMAS AGN 8.0)'
      names(63) = 'Planck 2015 0.06eV (BAHAMAS)'
      names(64) = 'Planck 2015 0.06eV (BAHAMAS; AGN 7.6)'
      names(65) = 'Planck 2015 0.06eV (BAHAMAS; AGN 8.0)'
      names(66) = 'Power bump: A = 0.15; k = 0.05h/Mpc; sigma = 1.0'
      names(67) = 'Power bump: A = 0.15; k = 0.10h/Mpc; sigma = 1.0'
      names(68) = 'Power bump: A = 0.15; k = 1.00h/Mpc; sigma = 1.0'
      names(69) = 'Axel power no bump'
      names(70) = 'WMAP9 (Extreme low AGN temperature)'
      names(71) = 'WMAP9 (Extreme high AGN temperature)'
      names(72) = 'Power bump: A = 0.15; k = 0.05h/Mpc; sigma = 0.1'
      names(73) = 'Power bump: A = 0.15; k = 0.10h/Mpc; sigma = 0.1'
      names(74) = 'Power bump: A = 0.15; k = 1.00h/Mpc; sigma = 0.1'
      names(75) = 'WMAP9 with 0.06eV neutrinos (BAHAMAS)'
      names(76) = 'WMAP9 with 0.12eV neutrinos (BAHAMAS)'
      names(77) = 'WMAP9 with 0.24eV neutrinos (BAHAMAS)'
      names(78) = 'WMAP9 with 0.48eV neutrinos (BAHAMAS)'
      names(79) = 'Power bump: A = 0.15; k = 0.05h/Mpc; sigma = 0.3'
      names(80) = 'Power bump: A = 0.15; k = 0.10h/Mpc; sigma = 0.3'
      names(81) = 'Power bump: A = 0.15; k = 1.00h/Mpc; sigma = 0.3'
      names(82) = 'Harrison-Zeldovich'
      names(83) = 'Boring but with some CDM replaced with 1.00eV neutrinos'
      names(84) = 'f(R) with F4, n=1'
      names(85) = 'f(R) with F5, n=1'
      names(86) = 'f(R) with F6, n=1'
      names(87) = 'Linearised f(R) with F4, n=1'
      names(88) = 'Linearised f(R) with F5, n=1'
      names(89) = 'Linearised f(R) with F6, n=1'
      names(90) = 'Linearised nDGP - Strong'
      names(91) = 'Linearised nDGP - Medium'
      names(92) = 'Linearised nDGP - Weak'
      names(93) = 'Planck 2015 0.12eV (BAHAMAS)'
      names(94) = 'Planck 2015 0.24eV (BAHAMAS)'
      names(95) = 'Planck 2015 0.48eV (BAHAMAS)'
      names(96) = 'WMAP9 (CAMB)'
      names(97) = 'WMAP9 (BAHAMAS AGN 7.6; CAMB)'
      names(98) = 'WMAP9 (BAHAMAS AGN 8.0; CAMB)'
      names(99) = 'Boring but with no-wiggle linear power'
      names(238) = 'Random Mira Titan cosmology with constant dark energy'
      names(239) = 'Bolshoi: Planck'
      names(240) = 'HMcode test: Open'
      names(241) = 'HMcode test: Low w'
      names(242) = 'HMcode test: High w'
      names(243) = 'HMcode test: Medium-mass neutrinos'
      names(244) = 'HMcode test: High-mass neutrinos'
      names(245) = 'HMcode test: Low spectral index'
      names(246) = 'HMcode test: High spectral index'
      names(247) = 'HMcode test: Low baryon fraction'
      names(248) = 'HMcode test: High baryon fraction'
      names(249) = 'HMcode test: Early dark energy'
      names(250) = 'HMcode test: Low AGN temperature'
      names(251) = 'HMcode test: High AGN temperature'
      names(252) = 'HMcode test:'
      names(253) = 'Random nu-w(a)CDM cosmology with random AGN temperature'
      names(254) = 'Random LCDM cosmology with random AGN temperature'
      names(255) = 'Random nu-LCDM cosmology with random AGN temperature'
      names(256) = 'Random w(a)CDM cosmology with random AGN temperature'
      names(257) = 'Random wCDM cosmology with random AGN temperature'
      names(258) = 'CAMB test cosmology'
      names(259) = 'Boring but with low sigma_8 = 0.5'
      names(260) = 'Boring bit with high sigma_8 = 1.2'
      names(261) = 'Millenium WMAP 1'
      names(262) = 'TCDM - from rescaling'
      names(263) = 'WMAP3'
      names(264) = 'Boring normalised via As'
      names(265) = 'Fiducial Smith & Angulo (2019)'
      names(266) = 'Random Euclid cosmology'
      names(267) = 'Random w(a)CDM Euclid cosmology'
      names(268) = 'Random nu-wCDM Euclid cosmology'
      names(269) = 'Random wCDM Euclid cosmology'
      names(270) = 'Random BACCO cosmology'
      names(271) = 'Random NGen-HALOFIT cosmology'
      names(272) = 'Random NGen-HALOFIT cosmology without running index'
      names(273) = 'Random Dark Quest cosmology'
      names(274) = 'Rounded Planck (2018)'

      names(100) = 'Mira Titan M000'
      names(101) = 'Mira Titan M001'
      names(102) = 'Mira Titan M002'
      names(103) = 'Mira Titan M003'
      names(104) = 'Mira Titan M004'
      names(105) = 'Mira Titan M005'
      names(106) = 'Mira Titan M006'
      names(107) = 'Mira Titan M007'
      names(108) = 'Mira Titan M008'
      names(109) = 'Mira Titan M009'
      names(110) = 'Mira Titan M010'
      names(111) = 'Mira Titan M011'
      names(112) = 'Mira Titan M012'
      names(113) = 'Mira Titan M013'
      names(114) = 'Mira Titan M014'
      names(115) = 'Mira Titan M015'
      names(116) = 'Mira Titan M016'
      names(117) = 'Mira Titan M017'
      names(118) = 'Mira Titan M018'
      names(119) = 'Mira Titan M019'
      names(120) = 'Mira Titan M020'
      names(121) = 'Mira Titan M021'
      names(122) = 'Mira Titan M022'
      names(123) = 'Mira Titan M023'
      names(124) = 'Mira Titan M024'
      names(125) = 'Mira Titan M025'
      names(126) = 'Mira Titan M026'
      names(127) = 'Mira Titan M027'
      names(128) = 'Mira Titan M028'
      names(129) = 'Mira Titan M029'
      names(130) = 'Mira Titan M030'
      names(131) = 'Mira Titan M031'
      names(132) = 'Mira Titan M032'
      names(133) = 'Mira Titan M033'
      names(134) = 'Mira Titan M034'
      names(135) = 'Mira Titan M035'
      names(136) = 'Mira Titan M036'

      names(200) = 'Franken Emu M000'
      names(201) = 'Franken Emu M001'
      names(202) = 'Franken Emu M002'
      names(203) = 'Franken Emu M003'
      names(204) = 'Franken Emu M004'
      names(205) = 'Franken Emu M005'
      names(206) = 'Franken Emu M006'
      names(207) = 'Franken Emu M007'
      names(208) = 'Franken Emu M008'
      names(209) = 'Franken Emu M009'
      names(210) = 'Franken Emu M010'
      names(211) = 'Franken Emu M011'
      names(212) = 'Franken Emu M012'
      names(213) = 'Franken Emu M013'
      names(214) = 'Franken Emu M014'
      names(215) = 'Franken Emu M015'
      names(216) = 'Franken Emu M016'
      names(217) = 'Franken Emu M017'
      names(218) = 'Franken Emu M018'
      names(219) = 'Franken Emu M019'
      names(220) = 'Franken Emu M020'
      names(221) = 'Franken Emu M021'
      names(222) = 'Franken Emu M022'
      names(223) = 'Franken Emu M023'
      names(224) = 'Franken Emu M024'
      names(225) = 'Franken Emu M025'
      names(226) = 'Franken Emu M026'
      names(227) = 'Franken Emu M027'
      names(228) = 'Franken Emu M028'
      names(229) = 'Franken Emu M029'
      names(230) = 'Franken Emu M030'
      names(231) = 'Franken Emu M031'
      names(232) = 'Franken Emu M032'
      names(233) = 'Franken Emu M033'
      names(234) = 'Franken Emu M034'
      names(235) = 'Franken Emu M035'
      names(236) = 'Franken Emu M036'
      names(237) = 'Franken Emu M037'

      names(300) = 'Cosmic Emu M000'
      names(301) = 'Cosmic Emu M001'
      names(302) = 'Cosmic Emu M002'
      names(303) = 'Cosmic Emu M003'
      names(304) = 'Cosmic Emu M004'
      names(305) = 'Cosmic Emu M005'
      names(306) = 'Cosmic Emu M006'
      names(307) = 'Cosmic Emu M007'
      names(308) = 'Cosmic Emu M008'
      names(309) = 'Cosmic Emu M009'
      names(310) = 'Cosmic Emu M010'
      names(311) = 'Cosmic Emu M011'
      names(312) = 'Cosmic Emu M012'
      names(313) = 'Cosmic Emu M013'
      names(314) = 'Cosmic Emu M014'
      names(315) = 'Cosmic Emu M015'
      names(316) = 'Cosmic Emu M016'
      names(317) = 'Cosmic Emu M017'
      names(318) = 'Cosmic Emu M018'
      names(319) = 'Cosmic Emu M019'
      names(320) = 'Cosmic Emu M020'
      names(321) = 'Cosmic Emu M021'
      names(322) = 'Cosmic Emu M022'
      names(323) = 'Cosmic Emu M023'
      names(324) = 'Cosmic Emu M024'
      names(325) = 'Cosmic Emu M025'
      names(326) = 'Cosmic Emu M026'
      names(327) = 'Cosmic Emu M027'
      names(328) = 'Cosmic Emu M028'
      names(329) = 'Cosmic Emu M029'
      names(330) = 'Cosmic Emu M030'
      names(331) = 'Cosmic Emu M031'
      names(332) = 'Cosmic Emu M032'
      names(333) = 'Cosmic Emu M033'
      names(334) = 'Cosmic Emu M034'
      names(335) = 'Cosmic Emu M035'
      names(336) = 'Cosmic Emu M036'
      names(337) = 'Cosmic Emu M037'

      IF (verbose) WRITE (*, *) 'ASSIGN_COSMOLOGY: Assigning cosmological model parameters'

      IF (icosmo == -1) THEN
         WRITE (*, *) 'ASSIGN_COSMOLOGY: Choose cosmological model'
         WRITE (*, *) '==========================================='
         DO i = 1, size(names)-1
            IF (i >= 100 .AND. i <= 136) THEN
               ! Do nothing
            ELSE IF (i >= 200 .AND. i <= 237) THEN
               ! Do nothing
            ELSE IF (i >= 300 .AND. i <= 337) THEN
               ! Do nothing
            ELSE IF (names(i) .NE. '') THEN
               WRITE (*, *) i, '- ', trim(names(i))
            END IF
         END DO
         WRITE (*, *) ' 100 -> 136 - Mira Titan M000 -> M036'
         WRITE (*, *) ' 200 -> 237 - Franken Emu M000 -> M037'
         WRITE (*, *) ' 300 -> 337 - Cosmic Emu M000 -> M037'
         READ (*, *) icosmo
         WRITE (*, *) '==========================================='
      END IF

      ! Set verbosity
      cosm%verbose = verbose

      ! Set the name of the cosmological model
      cosm%name = names(icosmo)

      ! Linear power spectrum
      cosm%iTk = iTk_EH ! Default to Eisenstein & Hu (1998)
      cosm%iTc = iTc_EH ! Default to Eisenstein & Hu (1999)

      ! Boring default cosmology
      cosm%Om_m = 0.3    ! Total matter (CDM + baryons + massive neutrino) density
      cosm%Om_b = 0.05   ! Baryon density
      cosm%Om_v = 0.7    ! Vacuum density
      cosm%Om_w = 0.     ! Dark-energy density (in addition to vacuum density)
      cosm%m_nu = 0.     ! Neutrino mass
      cosm%h = 0.7       ! Dimensionless Hubble parameter
      cosm%ns = 0.96     ! Spectral index
      cosm%nrun = 0.     ! Spectral tilt
      cosm%w = -1.       ! Dark energy equation of state
      cosm%wa = 0.       ! Dark energy time-varying equation of state
      cosm%T_CMB = 2.725 ! CMB temperature [K]
      cosm%neff = 3.046  ! Effective number of relativistic neutrinos
      cosm%YH = 0.76     ! Hydrogen mass fraction
      cosm%N_nu = 3      ! Integer number of (currently degenerate) massive neutrinos

      ! Dark energy
      cosm%iw = iw_LCDM ! Dark energy type

      ! Modified gravity
      cosm%img = img_none ! Modified gravity type
      cosm%H0rc = 0.      ! Note that zero is very strong
      cosm%fR0 = 0.       ! f(R) amplitude parameter
      cosm%nfR = 0        ! f(R) index parameter

      ! Warm dark matter
      cosm%warm = .FALSE. ! Is CDM actually WDM?
      cosm%m_wdm = 0.     ! WDM particle mass [keV]

      ! AGN
      cosm%Theat = 10**7.8 ! AGN temperature for feedback model

      !! Normalisation !!

      ! Should initially be set to 1 and then will be changed later by normalisation methods
      cosm%A = 1.

      ! Power spectrum normalisation
      !cosm%norm_method = norm_none   ! No normalisation, just take value of A from below
      cosm%norm_method = norm_sigma8 ! Normalise using sigma_8
      !cosm%norm_method = norm_pval   ! Large-scale structure normalisation at a specific wavenumber at z=0
      !cosm%norm_method = norm_As     ! As normalisation   
      cosm%sig8 = 0.8                ! norm_sigma8: sigma(R=8, a=1) normalisation
      cosm%kval = 0.001              ! norm_pval: Wavenumber for normalisation [h/Mpc]
      cosm%pval = 0.1973236854e-06   ! norm_pval: Delta^2 value to get sig8 = 0.8 for a boring cosmology
      cosm%kpiv = 0.05/0.7           ! norm_As: Wavenumber at which to define As normalisation NOTE: [h/Mpc]
      !cosm%As = 2.1e-9               ! norm_As: 2.1e-9 is a sensible value
      cosm%As = 1.97269e-9           ! norm_As: This gives sigma_8 = 0.8 in the boring cosmology with kpiv = 0.05/Mpc
      
      !! !!

      ! Power bump
      cosm%bump = 0
      cosm%A_bump = 0.
      cosm%k_bump = 0.
      cosm%sigma_bump = 0.

      ! Alternative dark energy models
      cosm%a1 = 0.
      cosm%a2 = 0.
      cosm%nstar = 0.
      cosm%ws = 0.
      cosm%am = 0.
      cosm%dm = 0.
      cosm%wm = 0.
      cosm%b0 = 0.
      cosm%b1 = 0.
      cosm%b2 = 0.
      cosm%b3 = 0.
      cosm%b4 = 0.

      ! Consider box size
      cosm%box = .FALSE.
      cosm%Lbox = 100. ! Box size [Mpc/h]

      ! Set is flags to negative
      cosm%is_init = .FALSE.
      cosm%is_normalised = .FALSE.
      cosm%analytical_Tk = .FALSE.

      ! Interpolators
      CALL reset_interpolator_status(cosm)

      ! Omegas for power spectrum if different from background cosmological parameters; false by default
      cosm%Om_m_pow = 0.
      cosm%Om_b_pow = 0.
      cosm%h_pow = 0.
      cosm%power_Omegas = .FALSE. 

      ! Gas options
      cosm%derive_gas_numbers = .TRUE.

      ! CAMB
      cosm%CAMB_exe = 'camb'
      cosm%CAMB_temp_dir = '/Users/Mead/Physics/CAMB_files/tmp/'

      IF (icosmo == 1) THEN
         ! Boring - do nothing
      ELSE IF (icosmo == 2) THEN
         ! cosmo-OWLS - WMAP7 (1312.5462)
         cosm%iTk = iTk_CAMB
         cosm%Om_m = 0.272
         cosm%Om_b = 0.0455
         cosm%Om_v = 1.-cosm%Om_m
         cosm%h = 0.704
         cosm%sig8 = 0.81
         cosm%ns = 0.967
      ELSE IF (icosmo == 3) THEN
         ! Planck 2013 (cosmo-OWLS/BAHAMAS; 1312.5462/1603.02702; no neutrinos)
         cosm%iTk = iTk_CAMB
         cosm%Om_m = 0.3175
         cosm%Om_b = 0.0490
         cosm%Om_v = 1.-cosm%Om_m
         cosm%h = 0.6711
         cosm%ns = 0.9624
         cosm%sig8 = 0.8341
      ELSE IF (is_in_array(icosmo, [4, 61, 62, 70, 71, 75, 76, 77, 78, 96, 97, 98])) THEN
         ! BAHAMAS - WMAP 9
         !  4 - WMAP9 (1712.02411; no neutrinos; EH Tk)
         ! 61 - WMAP9 with T_AGN = 7.6 (low; EH Tk)
         ! 62 - WMAP9 with T_AGN = 8.0 (high; EH Tk)
         ! 70 - WMAP9 with T_AGN = 7.0 (extremely low; EH Tk)
         ! 71 - WMAP9 with T_AGN = 8.6 (extrmely high; EH Tk)
         ! 75 - WMAP9 with 0.06eV neutrinos (CAMB Tk)
         ! 76 - WMAP9 with 0.12eV neutrinos (CAMB Tk)
         ! 77 - WMAP9 with 0.24eV neutrinos (CAMB Tk)
         ! 78 - WMAP9 with 0.48eV neutrinos (CAMB Tk)
         ! 96 - WMAP9 (as 4: T_AGN 7.8 but with CAMB Tk)
         ! 97 - WMAP9 (as 4: T_AGN 7.6 but with CAMB Tk)
         ! 98 - WMAP9 (as 4: T_AGN 8.0 but with CAMB Tk)
         IF (is_in_array(icosmo, [4, 61, 62, 70, 71])) THEN
            cosm%iTk = iTk_EH
         ELSE
            cosm%iTk = iTk_CAMB
         END IF
         cosm%h = 0.7000
         cosm%Om_b = 0.0463
         cosm%Om_m = 0.2330+cosm%Om_b
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns = 0.9720
         cosm%sig8 = 0.8211
         cosm%derive_gas_numbers = .FALSE.
         cosm%mup = 0.61
         Xi = 1.08
         Xe = 1.17
         cosm%mue = cosm%mup*(Xe+Xi)/Xe
         IF (icosmo == 61 .OR. icosmo == 97) THEN
            cosm%Theat = 10**7.6 ! Low AGN temperature
         ELSE IF (icosmo == 62 .OR. icosmo == 98) THEN
            cosm%Theat = 10**8.0 ! High AGN temperature
         ELSE IF (icosmo == 70) THEN
            cosm%Theat = 10**7.0 ! Extremely low AGN temperature
         ELSE IF (icosmo == 71) THEN
            cosm%Theat = 10**8.6 ! Extremely high AGN temperature
         END IF
         IF (icosmo == 75) THEN
            cosm%m_nu = 0.06 ! 0.06eV (minimal) neutrino mass
            cosm%sig8 = 0.8069
         ELSE IF (icosmo == 76) THEN
            cosm%m_nu = 0.12 ! 0.12eV neutrinos
            cosm%sig8 = 0.7924
         ELSE IF (icosmo == 77) THEN
            cosm%m_nu = 0.24 ! 0.24eV neutrinos
            cosm%sig8 = 0.7600
         ELSE IF (icosmo == 78) THEN
            cosm%m_nu = 0.48 ! 0.48eV neutrinos
            cosm%sig8 = 0.7001
         END IF
      ELSE IF (icosmo == 5) THEN
         !  5 - Open model
         cosm%Om_v = 0.
      ELSE IF (icosmo == 6 .OR. icosmo == 15) THEN
         !  6 - Einstein-de Sitter (SCDM)
         ! 15 - Einstein-de Sitter (TCDM)
         IF (icosmo == 15) THEN
            cosm%power_Omegas = .TRUE.
            cosm%Om_m_pow = cosm%Om_m
            cosm%Om_b_pow = cosm%Om_b
            cosm%h_pow = cosm%h
         END IF
         cosm%Om_m = 1.
         cosm%Om_v = 0.
      ELSE IF (icosmo == 7) THEN
         ! IDE I
         cosm%iw = iw_IDE1
         cosm%astar = 0.1
         cosm%nstar = 3.
         cosm%Om_ws = 0.3
         cosm%Om_m = 0.3
         cosm%Om_v = 0.7
      ELSE IF (icosmo == 8) THEN
         ! IDE II model
         cosm%iw = iw_IDE2
         cosm%astar = 0.1
         cosm%nstar = 3.    
         cosm%Om_ws = 0.3
         cosm%Om_m = 0.3
         cosm%Om_w = 0.7
         cosm%Om_v = 0. ! No vacuum necessary here
      ELSE IF (icosmo == 9) THEN
         ! IDE III model
         cosm%iw = iw_IDE3
         cosm%astar = 0.1
         cosm%Om_ws = 0.3
         cosm%ws = 0.5
         cosm%Om_m = 0.3
         cosm%Om_w = 0.7
         cosm%Om_v = 0.
      ELSE IF (is_in_array(icosmo, [11, 12, 13, 90, 91, 92])) THEN
         ! nDGP
         IF (is_in_array(icosmo, [11, 12, 13])) THEN
            cosm%img = img_nDGP
         ELSE
            cosm%img = img_nDGP_lin
         END IF
         IF (icosmo == 11 .OR. icosmo == 90) THEN
            ! Strong         
            cosm%H0rc = 0.1
            cosm%sig8 = 0.8*(1.0176/0.7790) ! Normalise to boring at a<<1
         ELSE IF (icosmo == 12 .OR. icosmo == 91) THEN
            ! Medium
            cosm%H0rc = 0.5
            cosm%sig8 = 0.8*(0.8717/0.7790) ! Normalise to boring at a<<1
         ELSE IF (icosmo == 13 .OR. icosmo == 92) THEN
            ! Weak
            cosm%H0rc = 2.0
            cosm%sig8 = 0.8*(0.8093/0.7790) ! Normalise to boring at a<<1
         ELSE
            STOP 'ASSIGN_COSMOLOGY: Error, something went wrong with nDGP'
         END IF
      ELSE IF (icosmo == 14) THEN
         ! WDM
         cosm%warm = .TRUE.
         cosm%m_wdm = 1.
      ELSE IF (icosmo == 16) THEN
         ! w = -0.7
         cosm%iw = iw_wCDM
         cosm%w = -0.7
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
      ELSE IF (icosmo == 17) THEN
         ! w = -1.3
         cosm%iw = iw_wCDM
         cosm%w = -1.3
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
      ELSE IF (icosmo == 18) THEN
         ! wa = 0.5
         cosm%iw = iw_waCDM
         cosm%wa = 0.5
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
      ELSE IF (icosmo == 19) THEN
         ! wa = -0.5
         cosm%iw = iw_waCDM
         cosm%wa = -0.5
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
      ELSE IF (icosmo == 20) THEN
         ! w = -0.7; wa = -1.5
         cosm%iw = iw_waCDM
         cosm%w = -0.7
         cosm%wa = -1.5
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
      ELSE IF (icosmo == 21) THEN
         ! w = -1.3; wa = 0.5
         cosm%iw = iw_waCDM
         cosm%w = -1.3
         cosm%wa = 0.5
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
      ELSE IF (icosmo == 22 .OR. icosmo == 23) THEN
         ! IDE II models
         cosm%iw = iw_IDE2
         cosm%Om_m = 0.3
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0. ! No vacuum necessary here
         IF (icosmo == 22) THEN
            ! IDE2 with n* = 3
            cosm%nstar = 3.
            cosm%astar = 0.01
            cosm%Om_ws = 0.1
         ELSE IF (icosmo == 23) THEN
            ! IDE2 with n* = 10
            cosm%nstar = 10.
            cosm%astar = 0.1
            cosm%Om_ws = 0.02
         END IF
      ELSE IF (icosmo == 10 .OR. icosmo == 24 .OR. icosmo == 238) THEN
         ! Random Mira Titan cosmology
         CALL random_Mira_Titan_cosmology(cosm)
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
         IF (icosmo == 10)  cosm%m_nu = 0. ! No massive neutrinos
         IF (icosmo == 238) cosm%wa = 0.   ! No time-varying dark energy
      ELSE IF (icosmo == 25) THEN
         ! Random Franken Emu cosmology
         CALL random_Franken_Emu_cosmology(cosm)
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo == 26) THEN
         ! Boring with CAMB linear spectrum
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo == 27) THEN
         ! Illustris; L = 75 Mpc/h
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
         cosm%Om_m = 0.3089
         cosm%Om_b = 0.0486
         cosm%Om_v = 1.-cosm%Om_m
         cosm%h = 0.6774
         cosm%ns = 0.9667
         cosm%sig8 = 0.8159
         cosm%box = .TRUE.
         cosm%Lbox = 75. ! 75 Mpc/h box
      ELSE IF (icosmo == 28) THEN
         ! Finite box
         cosm%box = .TRUE.
      ELSE IF (icosmo == 29) THEN
         ! Boring: Open; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.88397
      ELSE IF (icosmo == 30) THEN
         ! Boring: w = -0.7; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.83253
      ELSE IF (icosmo == 31) THEN
         ! Boring: w = -1.3; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.77462
      ELSE IF (icosmo == 32) THEN
         ! Boring: w = -1.0; wa =  0.5; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.81022
      ELSE IF (icosmo == 33) THEN
         ! Boring: w = -1.0; wa = -0.5; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.79116
      ELSE IF (icosmo == 34) THEN
         ! Boring: w = -0.7; wa = -1.5; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.80090
      ELSE IF (icosmo == 35) THEN
         ! Boring: w = -1.3; wa =  0.5; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.78197
      ELSE IF (icosmo == 36) THEN
         ! Boring; EdS; z=1 normalisation for Mead 2017; LCDM
         cosm%sig8 = 0.65380
      ELSE IF (icosmo == 37 .OR. icosmo == 43 .OR. icosmo == 44) THEN
         ! 37 - Multidark: WMAP5 (also Bolshoi)
         ! 43 - Multidark: WMAP5 with lowered sigma_8
         ! 44 - Multidark: WMAP5 with Eisenstein & Hu transfer function
         cosm%h = 0.70
         cosm%Om_b = 0.0469
         cosm%Om_m = 0.27
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns = 0.95
         cosm%sig8 = 0.82 ! Seems wrong at z=0, data more like sigma_8 = 0.80
         cosm%iTk = iTk_CAMB ! CAMB
         IF(icosmo == 43) cosm%sig8 = 0.80 ! Check to see if better matches with lower sigma_8
         IF(icosmo == 44) cosm%iTk = iTk_EH ! Eisenstein & Hu T(k)
      ELSE IF (icosmo == 38) THEN
         ! Random cosmic emu model
         CALL random_Cosmic_Emu_cosmology(cosm)
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
      ELSE IF (is_in_array(icosmo, [39, 40, 45, 46, 47, 48])) THEN
         ! Random cosmologies
         IF(icosmo == 39) THEN
            CALL random_LCDM_cosmology(cosm)
         ELSE IF(icosmo == 40) THEN
            CALL random_cosmology(cosm)
         ELSE IF(icosmo == 45) THEN
            CALL random_waCDM_cosmology(cosm)
         ELSE IF(icosmo == 46) THEN
            CALL random_wCDM_cosmology(cosm)
         ELSE IF(icosmo == 47) THEN
            CALL random_LCDM_cosmology(cosm)
         ELSE IF(icosmo == 48) THEN
            CALL random_nuLCDM_cosmology(cosm)
         ELSE
            STOP 'ASSIGN_COSMOLOGY: Error, something went wrong with random cosmology'
         END IF
         IF(icosmo == 39) THEN
            cosm%iTk = iTk_EH
         ELSE
            cosm%iTk = iTk_CAMB
         END IF
      ELSE IF (icosmo == 41) THEN
         ! SCDM with high neutrino mass
         cosm%Om_m = 1.
         cosm%Om_v = 0.
         cosm%m_nu = 4.
      ELSE IF (icosmo == 56 .OR. icosmo == 42) THEN
         ! 56 - Planck 2018 (Plik, from Table 1 of https://arxiv.org/abs/1807.06209)
         ! 42 - Same, but with neutrino mass fixed to zero and nothing else changed
         cosm%iTk = iTk_CAMB
         cosm%h = 0.6732
         cosm%ns = 0.96605
         cosm%m_nu = 0.06
         IF(icosmo == 42) cosm%m_nu = 0.
         cosm%Om_m = 0.3158
         cosm%Om_b = 0.022383/cosm%h**2
         cosm%Om_v = 1.-cosm%Om_m
         cosm%sig8 = 0.8120
         IF (icosmo == 56) THEN
            cosm%As = 2.08924e-9
            cosm%kpiv = 0.05/cosm%h
         END IF
      ELSE IF (icosmo == 274) THEN
         ! 274 - Rounded Planck 2018
         cosm%iTk = iTk_CAMB
         cosm%h = 0.67
         cosm%ns = 0.97
         cosm%m_nu = 0.06
         cosm%Om_m = 0.32
         cosm%Om_b = 0.049
         cosm%Om_v = 1.-cosm%Om_m
         cosm%sig8 = 0.81
      ELSE IF(icosmo == 49) THEN
         ! CFHTLenS best-fitting cosmology (Heymans 2013; combined with WMAP 7)
         ! From first line of Table 3 of https://arxiv.org/pdf/1303.1808.pdf
         ! CFHTLenS combined with WMAP7 and R11
         cosm%Om_m = 0.255
         cosm%Om_b = 0.0437
         cosm%Om_v = 1.-cosm%om_m
         cosm%sig8 = 0.794
         cosm%ns = 0.967
         cosm%h = 0.717
      ELSE IF(icosmo == 50) THEN
         ! CFHTLenS best-fitting cosmology (Kilbinger 2013; combined with WMAP 7)
         ! From first line in Table 3 of https://arxiv.org/pdf/1212.3338.pdf
         ! CFHTLenS combined with WMAP7
         cosm%Om_m = 0.274
         cosm%Om_b = 0.0456
         cosm%Om_v = 1.-cosm%Om_m
         cosm%sig8 = 0.815
         cosm%ns = 0.966
         cosm%h = 0.702
      ELSE IF (is_in_array(icosmo, [51, 52, 53, 54, 55])) THEN
         ! CAMB fail cosmologies for HMcode comparisons
         cosm%iTk = iTk_CAMB
         cosm%iw = iw_waCDM
         cosm%Om_v = 0.
         IF (icosmo == 51) THEN
            ! Weird wiggle
            ! HMcode (2016, 2020) both fail
            ! Solved by increasing accuracy of sigma(R) integration; interaction between BAO and T(k)?
            cosm%Om_m = 0.16793
            cosm%Om_w = 1.-cosm%Om_m
            cosm%h = 0.69341
            cosm%Om_b = 0.05574
            cosm%ns = 0.81903
            cosm%sig8 = 0.73592
         ELSE IF (icosmo == 52) THEN
            ! Quite weird wiggle
            ! HMcode (2016, 2020) both fail
            ! Solved by increasing accuracy of sigma(R) integration; interaction between BAO and T(k)?
            cosm%Om_m = 0.16644
            cosm%Om_w = 1.-cosm%Om_m
            cosm%h = 0.82246
            cosm%Om_b = 0.054119
            cosm%ns = 0.84478
            cosm%sig8 = 0.62890
            cosm%w = -0.91698
            cosm%wa = -0.36060
         ELSE IF (icosmo == 53) THEN
            ! BAO wiggle damping slightly wrong; high neutrino and baryon fractions
            ! HMcode (2020) fails but (2016) passes
            ! Solved by making the BAO wiggle extraction scale-dependent in HMx
            cosm%Om_m = 0.19165
            cosm%Om_w = 1.-cosm%Om_m
            cosm%h = 0.50576
            cosm%Om_b = 0.05983
            cosm%ns = 0.79467
            cosm%sig8 = 0.8839
            cosm%w = -1.28965
            cosm%wa = -0.71800
            cosm%f_nu = 0.11069
            cosm%m_nu = neutrino_constant(cosm)*cosm%f_nu*cosm%Om_m*cosm%h**2
         ELSE IF (icosmo == 54) THEN
            ! Weird wiggle; low h; high baryon fraction
            ! HMcode (2016, 2020) both fail
            ! Solved by increasing accuracy of sigma(R) integration; interaction between BAO and T(k)?
            cosm%Om_m = 0.18504
            cosm%Om_b = 0.06696
            cosm%Om_w = 1.-cosm%Om_m
            cosm%sig8 = 0.77520
            cosm%ns = 0.98360
            cosm%h = 0.42278
            cosm%m_nu = 0.07639
            cosm%w = -1.26075
            cosm%wa = -0.07341
         ELSE IF (icosmo == 55) THEN
            ! Small-scale suppression, EDE-ish
            ! HMcode (2020) fails but (2016) passes
            ! Solved by synchronising initial time for growth ODE 
            cosm%Om_m = 0.28330
            cosm%Om_b = 0.04034
            cosm%Om_w = 1.-cosm%Om_m
            cosm%sig8 = 0.7213
            cosm%ns = 1.18159
            cosm%h = 0.99249
            cosm%f_nu = 0.010366
            cosm%m_nu = neutrino_constant(cosm)*cosm%f_nu*cosm%Om_m*cosm%h**2
            cosm%w = -0.73926
            cosm%wa = 0.61511
         ELSE
            STOP 'ASSIGN_COSMOLOGY: Error, something went wrong with CAMB difference cosmologies'
         END IF
      ELSE IF(icosmo == 57 .OR. icosmo == 58) THEN
         ! Bound dark energy (1812.01133)
         cosm%iw = iw_BDE
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
         IF(icosmo == 57) THEN
            ! Model from (1812.01133)
            cosm%b0 = -0.9296
            cosm%b1 = -3.752
            cosm%b2 = -5.926
            cosm%b3 = -4.022
            cosm%b4 = -0.999
         ELSE IF(icosmo == 58) THEN
            ! Extreme model
            cosm%b0 = 0.
            cosm%b1 = 0.
            cosm%b2 = 0.
            cosm%b3 = 0.
            cosm%b4 = -1.
         ELSE
            STOP 'ASSIGN_COSMOLOGY: Error, something went wrong with BDE'
         END IF 
      ELSE IF (icosmo == 59) THEN
         ! Bump in power
         cosm%norm_method = norm_pval ! Normalise like this to prevent bump annoying sigma8        
         cosm%bump = 1
         cosm%A_bump = 0.08
         cosm%k_bump = 5.
         cosm%sigma_bump = 0.5      
      ELSE IF (is_in_array(icosmo, [66, 67, 68, 69, 72, 73, 74, 79, 80, 81])) THEN
         ! Axel bump cosmologies
         cosm%h = 0.7
         cosm%Om_b = 0.05
         cosm%Om_m = 0.30
         cosm%ns = 0.96
         cosm%Om_v = 1.-cosm%Om_m
         cosm%norm_method = norm_pval
         cosm%pval = 1.995809e-7 ! Gives sigma8 = 0.8 for no bump   
         cosm%iTk = iTk_CAMB
         IF (is_in_array(icosmo, [66, 67, 68, 72, 73, 74, 79, 80, 81])) THEN
            cosm%bump = 2
            cosm%A_bump = 0.15
         END IF
         IF (icosmo == 66 .OR. icosmo == 67 .OR. icosmo == 68) cosm%sigma_bump = 1.
         IF (icosmo == 72 .OR. icosmo == 73 .OR. icosmo == 74) cosm%sigma_bump = 0.1
         IF (icosmo == 79 .OR. icosmo == 80 .OR. icosmo == 81) cosm%sigma_bump = 0.3
         IF (icosmo == 66 .OR. icosmo == 72 .OR. icosmo == 79) cosm%k_bump = 0.05
         IF (icosmo == 67 .OR. icosmo == 73 .OR. icosmo == 80) cosm%k_bump = 0.1
         IF (icosmo == 68 .OR. icosmo == 74 .OR. icosmo == 81) cosm%k_bump = 1.0
      ELSE IF (icosmo == 60) THEN
         ! Boring cosmology but with exciting neutrino mass
         cosm%m_nu = 0.3
         cosm%iTk = iTk_CAMB
      ELSE IF (is_in_array(icosmo, [63, 64, 65])) THEN
         ! BAHAMAS Planck 2015 cosmologies
         ! Note well that these cosmologies all have a  neutrino mass
         ! 63 - Planck 2015 with 0.06eV neutrinos (BAHAMAS; Table 1 of 1712.02411)
         ! 64 - Planck 2015 with 0.06eV neutrinos but with 10^7.6 AGN temperature
         ! 65 - Planck 2015 with 0.06eV neutrinos but with 10^8.0 AGN temperature
         cosm%iTk = iTk_CAMB
         cosm%m_nu = 0.06
         cosm%h = 0.6787
         cosm%Om_b = 0.0482
         cosm%Om_m = cosm%Om_b+0.2571+0.0014
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns = 0.9701
         cosm%sig8 = 0.8085
         IF (icosmo == 64) cosm%Theat = 10**7.6
         IF (icosmo == 65) cosm%Theat = 10**8.0
      ELSE IF (icosmo == 93) THEN
         ! 93 - BAHAMAS Planck 2015 but with 0.12eV neutrinos (other parameters changed too)
         cosm%iTk = iTk_CAMB
         cosm%m_nu = 0.12
         cosm%h = 0.6768
         cosm%Om_b = 0.0488
         cosm%Om_m = cosm%Om_b+0.2574+0.0029
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns = 0.9693
         cosm%sig8 = 0.7943
      ELSE IF (icosmo == 94) THEN
         ! 94 - BAHAMAS Planck 2015 but with 0.24eV neutrinos (other parameters changed too)
         cosm%iTk = iTk_CAMB
         cosm%m_nu = 0.24
         cosm%h = 0.6723
         cosm%Om_b = 0.0496
         cosm%Om_m = cosm%Om_b+0.2576+0.0057
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns = 0.9733
         cosm%sig8 = 0.7664
      ELSE IF (icosmo == 95) THEN
         ! 95 - BAHAMAS Planck 2015 but with 0.48eV neutrinos (other parameters changed too)
         cosm%iTk = iTk_CAMB
         cosm%m_nu = 0.48
         cosm%h = 0.6643
         cosm%Om_b = 0.0513
         cosm%Om_m = cosm%Om_b+0.2567+0.0117
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns = 0.9811
         cosm%sig8 = 0.7030
      ELSE IF (icosmo == 82) THEN
         ! Harrison - Zel'dovich
         cosm%iTk = iTk_none
         cosm%ns = 1.
      ELSE IF (icosmo == 83) THEN
         ! Boring cosmology but with very exciting neutrino mass
         cosm%m_nu = 1.
         cosm%iTk = iTk_CAMB
      ELSE IF (is_in_array(icosmo, [84, 85, 86, 87, 88, 89])) THEN
         ! f(R) models
         IF (is_in_array(icosmo, [84, 85, 86])) THEN
            cosm%img = img_fR
         ELSE
            cosm%img = img_fR_lin
         END IF
         cosm%nfR = 1
         IF (icosmo == 84 .OR. icosmo == 87) THEN
            cosm%fR0 = -1e-4
            cosm%sig8 = 0.8*sqrt(1.9732/1.6810) ! Normalise to boring at k<<1
         ELSE IF (icosmo == 85 .OR. icosmo == 88) THEN
            cosm%fR0 = -1e-5
            cosm%sig8 = 0.8*sqrt(1.9732/1.82991) ! Normalise to boring at k<<1
         ELSE IF (icosmo == 86 .OR. icosmo == 89) THEN
            cosm%fR0 = -1e-6
            cosm%sig8 = 0.8*sqrt(1.9732/1.9364) ! Normalise to boring at k<<1
         ELSE
            STOP 'ASSIGN_COSMOLOGY: Something went wrong with f(R) models'
         END IF
      ELSE IF (icosmo == 99) THEN
         ! No wiggle linear power
         cosm%iTk = iTk_nw
      ELSE IF (icosmo == 239) THEN
         ! Bolshoi: Planck (https://www.cosmosim.org/cms/simulations/bolshoip/)
         cosm%Om_m = 0.30711
         cosm%Om_v = 1.-cosm%Om_m
         cosm%h = 0.70
         cosm%Om_b = 0.048
         cosm%ns = 0.96
         cosm%sig8 = 0.82
      ELSE IF (is_in_array(icosmo, [240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252])) THEN
         ! HMcode test cosmologies
         ! 240 - Open
         ! 241 - Low w
         ! 242 - High w
         ! 243 - Medium-mass neutrinos
         ! 244 - High-mass neutrinos
         ! 245 - Low spectral index
         ! 246 - High spectral index
         ! 247 - Low baryon fraction
         ! 248 - High baryon fraction
         ! 249 - Early dark energy
         ! 250 - Low feedback temperature
         ! 251 - High feedback temperature
         ! 252 - 
         cosm%iTk = iTk_CAMB
         IF (icosmo == 240) THEN
            ! 240 -  Open
            cosm%Om_v = 0.
         ELSE IF (icosmo == 241 .OR. icosmo == 242) THEN
            ! 241 - Low w
            ! 242 - High w
            cosm%iw = iw_wCDM
            cosm%Om_w = cosm%Om_v
            cosm%Om_v = 0.
            IF (icosmo == 241) THEN
               cosm%w = -0.7
               cosm%As = 2.46981e-9
            END IF 
            IF (icosmo == 242) THEN
               cosm%w = -1.3
               cosm%As = 1.75070e-9
            END IF
         ELSE IF (icosmo == 243) THEN
            ! 243 - Medium-mass neutrinos
            cosm%m_nu = 0.3
            cosm%As = 2.35868e-9
         ELSE IF (icosmo == 244) THEN
            ! 244 - High-mass neutrinos
            cosm%m_nu = 0.9
            cosm%As = 3.43507e-9
         ELSE IF (icosmo == 245) THEN
            ! 245 - Low spectral index
            cosm%ns = 0.7
            cosm%As = 2.38515e-9
         ELSE IF (icosmo == 246) THEN
            ! 246 - High spectral index
            cosm%ns = 1.3
            cosm%As = 1.47601e-9
         ELSE IF (icosmo == 247) THEN
            ! 247 - Low baryon fraction
            cosm%Om_b = 0.01
            cosm%As = 1.20028e-9
         ELSE IF (icosmo == 248) THEN
            ! 248 - High baryon fraction
            cosm%Om_b = 0.1
            cosm%As = 3.88822e-9
         ELSE IF (icosmo == 249) THEN
            ! 249 - Early dark energy
            cosm%iw = iw_waCDM
            cosm%wa = 0.9
            cosm%Om_w = cosm%Om_v
            cosm%Om_v = 0.
            cosm%As = 3.35538e-9
         ELSE IF (icosmo == 250) THEN
            ! 250 - Low AGN temperature
            cosm%Theat = 10**7.6
         ELSE IF (icosmo == 251) THEN
            ! 251 - High AGN temperature
            cosm%Theat = 10**8.0
         END IF
      ELSE IF (is_in_array(icosmo, [253, 254, 255, 256, 257])) THEN
         ! Random nu-w(a)CDM with random feedback temperature
         IF(icosmo == 253) CALL random_cosmology(cosm)
         IF(icosmo == 254) CALL random_LCDM_cosmology(cosm)
         IF(icosmo == 255) CALL random_nuLCDM_cosmology(cosm)
         IF(icosmo == 256) CALL random_waCDM_cosmology(cosm)
         IF(icosmo == 257) CALL random_wCDM_cosmology(cosm)
         CALL random_AGN_temperature(cosm)
         cosm%iTk = iTk_CAMB
      ELSE IF (icosmo == 258) THEN
         ! Test cosmology from camb_test.py
         cosm%iTk = iTk_CAMB
         cosm%norm_method = norm_As
         cosm%h = 0.675
         cosm%Om_b = 0.022/cosm%h**2
         cosm%Om_c = 0.122/cosm%h**2
         cosm%m_nu = 0.07
         cosm%N_nu = 3
         cosm%Om_nu = cosm%m_nu/(neutrino_constant(cosm)*cosm%h**2)
         cosm%Om_m = cosm%Om_b+cosm%Om_c+cosm%Om_nu
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns = 0.965
         cosm%As = 2e-9
      ELSE IF (icosmo ==  259) THEN
         ! Boring but with low sigma_8
         cosm%sig8 = 0.5
      ELSE IF (icosmo ==  260) THEN
         ! Boring but with high sigma_8
         cosm%sig8 = 1.2
      ELSE IF (icosmo == 261) THEN
         ! Millenium - WMAP 1
         cosm%Om_m = 0.25
         cosm%Om_v = 1.-cosm%Om_m
         cosm%Om_b = 0.045
         cosm%h = 0.73
         cosm%sig8 = 0.9
         cosm%ns = 1.
         cosm%iTk = iTK_CAMB
      ELSE IF (icosmo == 262) THEN
         ! TCDM (from Mead & Peacock 2014; 1308.5183; spectral shape similar to LCDM)
         cosm%Om_m = 1.
         cosm%Om_v = 0.
         cosm%Om_b = 0.
         cosm%h = 0.5
         cosm%sig8 = 0.8
         cosm%ns = 1.
         cosm%power_Omegas = .TRUE.
         cosm%Om_m_pow = 0.3 ! Only need Gamma = Om_m*h to be 0.21
         cosm%Om_b_pow = 0.1 ! Only need Gamma = Om_m*h to be 0.21
         cosm%h_pow = 0.7    ! Only need Gamma = Om_m*h to be 0.21
         cosm%iTk = iTk_DEFW
      ELSE IF (icosmo == 263) THEN
         ! WMAP3 (from Angulo & White 2010; 0912.4277)
         cosm%Om_m = 0.238
         cosm%Om_v = 1.-cosm%Om_m
         cosm%Om_b = 0.0416
         cosm%h = 0.732
         cosm%sig8 = 0.761
         cosm%ns = 1.
         cosm%iTk = iTk_CAMB  
      ELSE IF (icosmo == 264) THEN
         ! Boring, but normalised via As
         cosm%norm_method = norm_As
      ELSE IF (icosmo == 265) THEN
         ! Fiducial from Smith & Angulo (2019)
         cosm%iTk = iTk_CAMB
         cosm%norm_method = norm_As
         wc = 0.11889
         wb = 0.022161
         wm = wb+wc
         !cosm%Om_v = 0.6914   ! This is the value written in the paper, which is incorrect
         cosm%Om_v = 0.6928849 ! This is the correct value according to Robert Smith
         cosm%Om_m = 1.-cosm%Om_v
         cosm%h = sqrt(wm/cosm%Om_m)
         cosm%Om_b = wb/cosm%h**2
         cosm%ns = 0.9611
         !cosm%sig8 = 0.8279
         cosm%As = 2.14818e-9
         cosm%kpiv = 0.05/cosm%h
         cosm%T_CMB = 2.726
         cosm%neff = 3.040
      ELSE IF (icosmo == 266 .OR. icosmo == 267 .OR. icosmo == 268 .OR. icosmo == 269) THEN
         ! Random Euclid cosmology
         CALL random_Euclid_cosmology(cosm)
         cosm%iTk = iTk_CAMB
         IF (icosmo == 267 .OR. icosmo == 269) cosm%m_nu = 0.
         IF (icosmo == 268 .OR. icosmo == 269) cosm%wa = 0.
      ELSE IF (icosmo == 270) THEN
         ! Random BACCO cosmology
         CALL random_BACCO_cosmology(cosm)
         cosm%iTk = iTk_CAMB
      ELSE IF (icosmo == 271 .OR. icosmo == 272) THEN
         ! Random NGenHALOFIT cosmology
         CALL random_NGenHALOFIT_cosmology(cosm)
         cosm%iTk = iTk_CAMB
         IF (icosmo == 272) cosm%nrun = 0.
      ELSE IF (icosmo == 273) THEN
         ! Random Dark Quest cosmology
         CALL random_Dark_Quest_cosmology(cosm)
         cosm%iTk = iTk_CAMB
      ELSE IF (icosmo >= 100 .AND. icosmo <= 137) THEN
         ! Mira Titan nodes
         CALL Mira_Titan_node_cosmology(icosmo-100, cosm)
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo >= 200 .AND. icosmo <= 237) THEN
         ! Franken Emu nodes (which are the same as Franken Emu nodes)
         CALL Franken_Emu_node_cosmology(icosmo-200, cosm)
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo >= 300 .AND. icosmo <= 337) THEN
         ! Cosmic Emu nodes (which are the same as Franken Emu nodes)
         CALL Cosmic_Emu_node_cosmology(icosmo-300, cosm)
         cosm%iTk = iTk_CAMB ! Set to CAMB linear power
      ELSE
         STOP 'ASSIGN_COSMOLOGY: Error, icosmo not specified correctly'
      END IF

      IF (cosm%verbose) THEN
         WRITE (*, *) 'ASSIGN_COSMOLOGY: Cosmology: ', trim(cosm%name)
         WRITE (*, *) 'ASSIGN_COSMOLOGY: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE assign_cosmology

   SUBROUTINE init_cosmology(cosm)

      ! Calcualtes derived parameters
      ! TODO: Force the cold transfer function to come from CAMB if possible?
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Xs, f1, f2, neff_constant
      REAL :: rho_g, Om_g_h2, f_nu_rad
      REAL, PARAMETER :: small = small_curve ! Some small number for writing curvature things

      ! Is statements
      cosm%is_init = .FALSE.
      cosm%is_normalised = .FALSE.

      ! TODO: Force the cold transfer function to come from CAMB if possible?
      !IF (cosm%iTk == iTk_CAMB) cosm%iTc = iTc_CAMB
      
      ! Overall power normalisaiton, should initially be unity
      cosm%A = 1.

      ! Things to do with finite box
      IF (cosm%box) cosm%kbox = twopi/cosm%Lbox

      IF (cosm%verbose) WRITE (*, *) 'INIT_COSMOLOGY: Calculating derived parameters'

      ! Calculate radiation density (includes photons and neutrinos at recombination)
      cosm%T_nu = cosm%T_CMB*(4./11.)**(1./3.)           ! Neutrino temperature [K]
      rho_g = 4.*SBconst*cosm%T_CMB**4/c_light**3        ! Photon physical density at z=0 from CMB temperature [kg/m^3]
      Om_g_h2 = rho_g/critical_density                   ! Photon cosmological density [h^2]
      cosm%Om_g = Om_g_h2/cosm%h**2                      ! Photon density parameter
      neff_constant = (7./8.)*(4./11.)**(4./3.)          ! Number that converted photon density to neutrino density (~ 0.2271)
      cosm%Om_nu_rad = cosm%Om_g*cosm%neff*neff_constant ! Relativisitic neutrino density (assuming they never behave like radiation always)
      cosm%Om_r = cosm%Om_g+cosm%Om_nu_rad               ! Radiation is sum of photon and neutrino densities
      f_nu_rad = cosm%Om_nu_rad/cosm%Om_r                ! Fraction of radiation that is neutrinos (~0.40)

      ! Redshift of the last-scatting surface (depends on Omega_m and Omega_b)
      cosm%z_CMB = z_CMB(cosm)
      IF (cosm%verbose) WRITE(*, *) 'INIT_COSMOLOGY: z_LSS:', cosm%z_CMB

      ! Information about how radiation density is calculated
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: Omega_g:', cosm%Om_g
         WRITE (*, *) 'INIT_COSMOLOGY: Omega_nu (radiation):', cosm%Om_nu_rad
         WRITE (*, *) 'INIT_COSMOLOGY: Omega_r:', cosm%Om_r
         WRITE (*, *) 'INIT_COSMOLOGY: f_nu (radiation):', f_nu_rad
      END IF

      ! Check that radiation density is not absurd
      IF (cosm%Om_r > 1e-3) THEN
         STOP 'INIT_COSMOLOGY: Error, radiation density is too high'
      END IF

!!$    ! Correction to vacuum density in order for radiation to maintain flatness
!!$    cosm%Om_v_mod=cosm%Om_v-cosm%Om_r
!!$    If(cosm%verbose) THEN
!!$       WRITE(*,*) 'INIT_COSMOLOGY: Altering vacuum density to account for radiation and maintain flatness'
!!$       WRITE(*,*) 'INIT_COSMOLOGY: Omega_v prior to change:', cosm%Om_v
!!$       WRITE(*,*) 'INIT_COSMOLOGY: Omega_v post change:', cosm%Om_v_mod
!!$    END IF

      ! Massive neutrinos
      ! TODO: Add support for non-degenerate neutrino species and masses
      IF (cosm%m_nu == 0.) THEN
         cosm%Om_nu = cosm%Om_nu_rad
         cosm%a_nu = 1. ! TODO: Should this be larger?
         cosm%f_nu = 0.
      ELSE
         IF (neutrino_method == neutrino_basic) THEN
            cosm%Om_nu = cosm%m_nu/(neutrino_constant(cosm)*cosm%h**2)
            IF (cosm%Om_nu >= cosm%Om_nu_rad) THEN
               cosm%a_nu = cosm%Om_nu_rad/cosm%Om_nu
            ELSE
               cosm%Om_nu = cosm%Om_nu_rad
               cosm%a_nu = 1. ! TODO: Should this be larger?
            END IF
         ELSE IF (neutrino_method == neutrino_Komatsu) THEN
            ! This does not use the 94.1eV approximation
            ! TODO: Should this be a division by neff or by 3?
            cosm%a_nu = cosm%T_nu/(cosm%m_nu/cosm%neff)*(kB/eV) 
            cosm%Om_nu = cosm%Om_nu_rad*Komatsu_nu(1./cosm%a_nu)
         ELSE
            STOP 'INIT_COSMOLOGY: Error, neutrino method not recognised'
         END IF
         cosm%f_nu = cosm%Om_nu/cosm%Om_m
      END IF

      ! Write neutrino information to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: Omega_nu:', cosm%Om_nu
         WRITE (*, *) 'INIT_COSMOLOGY: a_nu:', cosm%a_nu
         WRITE (*, *) 'INIT_COSMOLOGY: z_nu:', redshift_a(cosm%a_nu)
         WRITE (*, *) 'INIT_COSMOLOGY: f_nu:', cosm%f_nu
      END IF

      ! Check neutrino mass fraction is not too high
      ! TODO: Improve modelling for very light neutrinos
      IF (cosm%f_nu > f_nu_limit) STOP 'INIT_COSMOLOGY: Error, neutrino mass fraction is too high'
      IF ((cosm%m_nu .NE. 0.) .AND. (cosm%a_nu > a_nu_limit)) THEN
         WRITE(*, *) 'INIT_COSMOLOGY: Neutrino mass [eV]:', cosm%m_nu
         STOP 'INIT_COSMOLOGY: Error, neutrinos are too light'
      END IF

      ! Decide on scale-dependent growth
      ! NOTE: This is the only place the scale_dependent_growth variable should be set (!!!)
      IF ((cosm%m_nu .NE. 0.) .OR. (cosm%img == img_fR) .OR. (cosm%img == img_fR_lin)) THEN
         cosm%scale_dependent_growth = .TRUE.
      ELSE
         cosm%scale_dependent_growth = .FALSE.
      END IF

      ! Decide on triviality of the cold spectrum
      IF (cosm%m_nu == 0.) cosm%iTc = iTc_none

      ! Write to screen
      IF (cosm%verbose) THEN
         IF (cosm%scale_dependent_growth) WRITE (*, *) 'INIT_COSMOLOGY: Scale-dependent growth'
         IF (cosm%iTc == iTc_none)        WRITE (*, *) 'INIT_COSMOLOGY: Trivial cold spectrum'
      END IF

      ! Would need to include MGCAMB to make this work
      IF ((cosm%img .NE. img_none) .AND. (cosm%m_nu .NE. 0.)) THEN
         STOP 'INIT_COSMOLOGY: Error, modified gravity not compatible with massive neutrinos'
      END IF

      ! Would need to modify formulas to make this compatable
      IF (((cosm%img == img_fR) .OR. (cosm%img == img_fR_lin)) .AND. (cosm%Om_w /= 0.)) THEN
         STOP 'INIT_COSMOLOGY: f(R) gravity not currently compatible with dark energy'
      END IF

      ! Derived cosmological parameters
      cosm%Om_c = cosm%Om_m-cosm%Om_b ! Omega_m defined to include CDM, baryons and massive neutrinos
      IF (cosm%m_nu .NE. 0.) cosm%Om_c = cosm%Om_c-cosm%Om_nu
      cosm%Om = cosm%Om_m+cosm%Om_v+cosm%Om_w    ! Ignore radiation here
      cosm%Om_k = 1.-cosm%Om
      cosm%k = (cosm%Om-1.)/(Hdist**2)
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: Omega_c:', cosm%Om_c
         WRITE (*, *) 'INIT_COSMOLOGY: Omega:', cosm%Om
         WRITE (*, *) 'INIT_COSMOLOGY: Omega_k:', cosm%Om_k
         WRITE (*, *) 'INIT_COSMOLOGY: k [Mpc/h]^-2:', cosm%k
         IF (abs(cosm%k) > small) THEN
            WRITE (*, *) 'INIT_COSMOLOGY: k_rad [Mpc/h]:', 1./sqrt(abs(cosm%k)) ! Curvature radius
         END IF
      END IF

      IF (cosm%Om_c < 0.) STOP 'INIT_COSMOLOGY: Error, CDM density is negative'

      ! Physical density parameters
      cosm%omega_m = cosm%Om_m*cosm%h**2   ! Physical matter density
      cosm%omega_b = cosm%Om_b*cosm%h**2   ! Physical baryon density
      cosm%omega_c = cosm%Om_c*cosm%h**2   ! Physical CDM density
      cosm%omega_nu = cosm%Om_nu*cosm%h**2 ! Physical neutrino density

      ! Write physical densities to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: omega_m:', cosm%omega_m
         WRITE (*, *) 'INIT_COSMOLOGY: omega_b:', cosm%omega_b
         WRITE (*, *) 'INIT_COSMOLOGY: omega_c:', cosm%omega_c
         WRITE (*, *) 'INIT_COSMOLOGY: omega_nu:', cosm%omega_nu
      END IF

      ! Using different background Omegas compared to power Omegas
      IF (cosm%power_Omegas) THEN
         cosm%Om_c_pow = cosm%Om_m_pow-cosm%Om_b_pow
      ELSE
         cosm%Om_c_pow = 0.
      END IF

      ! Gas parameters
      IF (cosm%derive_gas_numbers) THEN
         ! Mean mass per gas particle divided by proton mass
         ! ~0.588 if fH=0.76, gas is ionised and H and He only; 0.61 in BAHAMAS
         cosm%mup = 4./(5.*cosm%YH+3.)
         ! Mean mass per gas electron divided by proton mass
         ! ~1.136 if fH=0.76, gas is ionised and H and He only; 1.17 in BAHAMAS
         cosm%mue = 2./(1.+cosm%YH)       
      END IF
      cosm%YHe = 1.-cosm%YH ! Helium mass fraction

      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: mu_p:', cosm%mup
         WRITE (*, *) 'INIT_COSMOLOGY: mu_e:', cosm%mue
      END IF

      ! Cosmology is now initialised
      cosm%is_init = .TRUE.

      ! Checks for dark energy models
      IF ((cosm%iw .NE. iw_waCDM) .AND. (cosm%wa .NE. 0.)) STOP 'INIT_COSMOLOGY: wa is set but waCDM is not being used'
      IF ((cosm%iw == iw_LCDM) .AND. (cosm%w .NE. -1.)) STOP 'INIT_COSMOLOGY: Dark energy is set to LCDM but w is not -1'

      ! Extra initialisation for dark-energy models
      IF (cosm%iw == iw_IDE1) THEN
         ! IDE I
         !Om_w=Om_w*(Om_m*astar**(-3)+Om_v)/(X(astar)*(1.-Om_w))
         f1 = cosm%Om_ws*X_de(cosm%astar, cosm)+cosm%Om_ws*cosm%astar**(-2)
         f2 = X_de(cosm%astar, cosm)*(1.-cosm%Om_ws)+cosm%Om_ws*cosm%astar**(-2)
         cosm%Om_w = cosm%Om_ws*(Hubble2(cosm%a, cosm)-f1/f2)
      ELSE IF (cosm%iw == iw_IDE2) THEN
         ! IDE II
         ! Define a1^n
         cosm%a1n = cosm%astar**cosm%nstar
         ! Necessary for first step below
         cosm%a2n = cosm%a1n
         ! All neccessary to convert parameters to a1,a2
         f1 = cosm%Om_ws*(Hubble2(cosm%astar, cosm)-cosm%Om_w*X_de(cosm%astar, cosm))
         f2 = cosm%Om_w*(1.-cosm%Om_ws)
         Xs = f1/f2
         Xs = Xs**(cosm%nstar/6.)
         ! Top and bottom of fraction
         f1 = cosm%a1n*(2.*Xs-(1.+cosm%a1n))
         f2 = (1.+cosm%a1n)-2.*Xs*cosm%a1n
         cosm%a2n = f1/f2 ! Finally! a2
         !IF(a2<a1) a2=a1
      ELSE IF (cosm%iw == iw_IDE3) THEN
         ! IDE III
         ! Scale-factor at which Om_w(a*) is most important
         cosm%a1 = cosm%astar
         ! Needs to be set for X(a*) and H2(a*) below (which cancel each other)
         cosm%a2 = cosm%astar
         f1 = cosm%Om_ws*(Hubble2(cosm%astar, cosm)-cosm%Om_w*X_de(cosm%astar, cosm))
         f2 = cosm%Om_w*(1.-cosm%Om_ws)
         cosm%a2 = cosm%astar*(f1/f2)**(1./(3.*(1.+cosm%ws)))
      END IF

      ! Useful variables
      cosm%age = 0.      
      cosm%gnorm = 0.
      cosm%horizon = 0.

      ! Interpolators
      CALL reset_interpolator_status(cosm)

      ! Switch analytical transfer function
      IF (is_in_array(cosm%iTk, [iTk_EH, iTk_DEFW, iTk_none, iTk_nw])) THEN
         cosm%analytical_Tk = .TRUE.
      ELSE
         cosm%analytical_Tk = .FALSE.
      END IF

      ! Write finishing message to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_cosmology

   SUBROUTINE print_cosmology(cosm)

      ! Prints the cosmological parameters to the screen
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: small = small_curve
      CHARACTER(len=16), PARAMETER :: format = '(A11,A16,F11.5)'
      CHARACTER(len=37), PARAMETER :: dashes = '====================================='

      IF (cosm%verbose) THEN
         WRITE (*, *) dashes
         WRITE (*, *) 'COSMOLOGY: ', trim(cosm%name)
         WRITE (*, *) dashes
         WRITE (*, *) 'COSMOLOGY: Background parameters'
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_m:', cosm%Om_m
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_b:', cosm%Om_b
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_v:', cosm%Om_v
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_w:', cosm%Om_w
         WRITE (*, fmt=format) 'COSMOLOGY:', 'h:', cosm%h
         WRITE (*, fmt=format) 'COSMOLOGY:', 'T_CMB [K]:', cosm%T_CMB
         WRITE (*, fmt=format) 'COSMOLOGY:', 'n_eff:', cosm%neff
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Y_H:', cosm%YH
         !WRITE(*,fmt=format) 'COSMOLOGY:', 'm_nu 1 [eV]:', cosm%m_nu(1)
         !WRITE(*,fmt=format) 'COSMOLOGY:', 'm_nu 2 [eV]:', cosm%m_nu(2)
         !WRITE(*,fmt=format) 'COSMOLOGY:', 'm_nu 3 [eV]:', cosm%m_nu(3)
         WRITE (*, fmt=format) 'COSMOLOGY:', 'M_nu [eV]:', cosm%m_nu
         IF (cosm%m_nu /= 0.) WRITE (*, fmt='(A11,A16,I11)') 'COSMOLOGY:', 'N_nu:', cosm%N_nu
         WRITE (*, *) dashes
         IF (cosm%iw == iw_LCDM) THEN
            WRITE (*, *) 'COSMOLOGY: Dark energy: Vacuum'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'w:', -1.
         ELSE IF (cosm%iw == iw_QUICC) THEN
            WRITE (*, *) 'COSMOLOGY: Dark energy: QUICC'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'w0:', cosm%w
            WRITE (*, fmt=format) 'COSMOLOGY:', 'wm:', cosm%wm
            WRITE (*, fmt=format) 'COSMOLOGY:', 'am:', cosm%am
            WRITE (*, fmt=format) 'COSMOLOGY:', 'dm:', cosm%dm
         ELSE IF (cosm%iw == iw_waCDM) THEN
            WRITE (*, *) 'COSMOLOGY: Dark energy: w(a)CDM'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'w0:', cosm%w
            WRITE (*, fmt=format) 'COSMOLOGY:', 'wa:', cosm%wa
         ELSE IF (cosm%iw == iw_wCDM) THEN
            WRITE (*, *) 'COSMOLOGY: Dark energy: Constant w'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'w:', cosm%w
         ELSE IF (cosm%iw == iw_IDE1) THEN
            WRITE (*, *) 'COSMOLOGY: Dark energy: IDE I'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'a*:', cosm%astar
            WRITE (*, fmt=format) 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
            WRITE (*, fmt=format) 'COSMOLOGY:', 'n*:', cosm%nstar
         ELSE IF (cosm%iw == iw_IDE2) THEN
            WRITE (*, *) 'COSMOLOGY: Dark energy: IDE II'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'a*:', cosm%astar
            WRITE (*, fmt=format) 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
            WRITE (*, fmt=format) 'COSMOLOGY:', 'n*:', cosm%nstar
         ELSE IF (cosm%iw == iw_IDE3) THEN
            WRITE (*, *) 'COSMOLOGY: Dark energy: IDE III'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'a*:', cosm%a1
            WRITE (*, fmt=format) 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
            WRITE (*, fmt=format) 'COSMOLOGY:', 'w*:', cosm%ws
         ELSE IF (cosm%iw == iw_BDE) THEN
            WRITE(*, *) 'COSMOLOGY: Dark energy: Bound'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'b0:', cosm%b0
            WRITE (*, fmt=format) 'COSMOLOGY:', 'b1:', cosm%b1
            WRITE (*, fmt=format) 'COSMOLOGY:', 'b2:', cosm%b2
            WRITE (*, fmt=format) 'COSMOLOGY:', 'b3:', cosm%b3
            WRITE (*, fmt=format) 'COSMOLOGY:', 'b4:', cosm%b4
         END IF
         WRITE (*, *) dashes
         IF (cosm%img .NE. img_none) THEN
            IF(cosm%img == img_nDGP .OR. cosm%img == img_nDGP_lin) THEN
               IF (cosm%img == img_nDGP) WRITE(*, *) 'COSMOLOGY: nDGP modified gravity'
               IF (cosm%img == img_nDGP_lin) WRITE(*, *) 'COSMOLOGY: Linearised nDGP modified gravity'
               WRITE(*, fmt=format) 'COSMOLOGY:', 'H0rc:', cosm%H0rc
            ELSE IF (cosm%img == img_fR .OR. cosm%img == img_fR_lin) THEN
               IF (cosm%img == img_fR) WRITE(*, *) 'COSMOLOGY: f(R) with LCDM background'
               IF (cosm%img == img_fR_lin) WRITE(*, *) 'COSMOLOGY: Linearised f(R) with LCDM background'
               WRITE(*, fmt=format) 'COSMOLOGY:', 'log10(-fR0):', log10(-cosm%fR0)
               WRITE(*, fmt=format) 'COSMOLOGY:', 'nfR:', cosm%nfR
            END IF
            WRITE (*, *) dashes
         END IF
         IF(cosm%iTk == iTk_none) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: Pure power law'
         ELSE IF(cosm%iTk == iTk_EH) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: Eisenstein & Hu'
         ELSE IF(cosm%iTk == iTk_nw) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: No-wiggle'
         ELSE IF(cosm%iTk == iTk_CAMB) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: CAMB'
         ELSE IF(cosm%iTk == iTk_DEFW) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: DEFW'
         ELSE IF(cosm%iTk == iTk_external) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: External'
         ELSE
            STOP 'COSMOLOGY: Error, iTk not set properly'
         END IF   
         WRITE (*, fmt=format) 'COSMOLOGY:', 'n_s:', cosm%ns
         WRITE (*, fmt=format) 'COSMOLOGY:', 'n_run:', cosm%nrun
         WRITE (*, fmt=format) 'COSMOLOGY:', 'n_runrun:', cosm%nrunrun
         WRITE (*, fmt=format) 'COSMOLOGY:', 'kpiv [h/Mpc]:', cosm%kpiv
         WRITE (*, fmt=format) 'COSMOLOGY:', 'kpiv [1/Mpc]:', cosm%kpiv*cosm%h
         IF(cosm%norm_method == norm_sigma8) THEN
            WRITE (*, *) 'COSMOLOGY: Normalisation: sigma_8'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'sigma_8:', cosm%sig8
         ELSE IF(cosm%norm_method == norm_pval) THEN
            WRITE (*, *) 'COSMOLOGY: Normalisation: Power'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'k [h/Mpc]:', cosm%kval
            WRITE (*, fmt=format) 'COSMOLOGY:', 'Delta^2:', cosm%pval
         ELSE IF(cosm%norm_method == norm_As) THEN
            WRITE (*, *) 'COSMOLOGY: Normalisation: A_s'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'As [10^-9]:', cosm%As*1e9
         END IF      
         IF (cosm%box) WRITE (*, fmt=format) 'COSMOLOGY:', 'L_box [Mpc/h]:', cosm%Lbox
         IF (cosm%warm) THEN
            WRITE (*, fmt=format) 'COSMOLOGY:', 'm_wdm [keV]:', cosm%m_wdm
         END IF
         IF (cosm%bump .NE. 0) THEN
            WRITE (*, fmt=format) 'COSMOLOGY:', 'A bu:', cosm%A_bump
            WRITE (*, fmt=format) 'COSMOLOGY:', 'k_bu [Mpc/h]:', cosm%k_bump
            WRITE (*, fmt=format) 'COSMOLOGY:', 'sigma bu:', cosm%sigma_bump
         END IF
         WRITE (*, *) dashes
         WRITE (*, *) 'COSMOLOGY: Baryon feedback'
         WRITE (*, fmt=format) 'COSMOLOGY:', 'log10(T_AGN/K):', log10(cosm%Theat)
         WRITE (*, *) dashes
         WRITE (*, *) 'COSMOLOGY: Derived parameters'
         WRITE (*, fmt=format) 'COSMOLOGY:', 'omega_m:', cosm%omega_m
         WRITE (*, fmt=format) 'COSMOLOGY:', 'omega_c:', cosm%omega_c
         WRITE (*, fmt=format) 'COSMOLOGY:', 'omega_b:', cosm%omega_b
         WRITE (*, fmt=format) 'COSMOLOGY:', 'omega_nu:', cosm%omega_nu
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_g:', cosm%Om_g
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_r:', cosm%Om_r
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega:', cosm%Om
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_c:', cosm%Om_c
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_nu:', cosm%Om_nu
         WRITE (*, fmt=format) 'COSMOLOGY:', 'f_nu:', cosm%f_nu
!!$       WRITE(*,fmt=format) 'COSMOLOGY:', 'Omega_v'':', cosm%Om_v_mod
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Omega_k:', cosm%Om_k
         WRITE (*, fmt=format) 'COSMOLOGY:', 'k [Mpc/h]^-2:', cosm%k
         IF (abs(cosm%k) > small) THEN
            WRITE (*, fmt=format) 'COSMOLOGY:', 'k_rad [Mpc/h]:', 1./sqrt(abs(cosm%k))
         END IF
         WRITE (*, fmt=format) 'COSMOLOGY:', 'mu_p:', cosm%mup
         WRITE (*, fmt=format) 'COSMOLOGY:', 'mu_e:', cosm%mue
         WRITE (*, fmt=format) 'COSMOLOGY:', 'z_CMB:', cosm%z_CMB
         IF (cosm%iw == iw_IDE2) THEN
            WRITE (*, *) dashes
            WRITE (*, *) 'COSMOLOGY: IDE II'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'a1^n:', cosm%a1n
            WRITE (*, fmt=format) 'COSMOLOGY:', 'a2^n:', cosm%a2n
         END IF
         WRITE (*, *) dashes
         WRITE (*, *)
      END IF

   END SUBROUTINE print_cosmology

   SUBROUTINE assign_init_cosmology(icosmo, cosm, verbose)

      ! Both assigns and initialises the cosmological model
      INTEGER, INTENT(INOUT) :: icosmo
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose

      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

   END SUBROUTINE assign_init_cosmology

   SUBROUTINE reset_interpolator_status(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      cosm%has_distance = .FALSE.
      cosm%has_growth = .FALSE.
      cosm%has_sigma = .FALSE.
      cosm%has_spherical = .FALSE.
      cosm%has_power = .FALSE.
      cosm%has_time = .FALSE.
      cosm%has_Xde = .FALSE.
      cosm%has_wiggle = .FALSE.
      cosm%has_SPT = .FALSE.

   END SUBROUTINE reset_interpolator_status

   REAL FUNCTION neutrino_constant(cosm)

      ! Critical mass for neutrino density to close Universe [eV] 
      ! Roughly 94.1 eV, or is it 93.03 eV, or 93.14 eV?; https://arxiv.org/pdf/1812.02102.pdf
      ! Not really a constant because it depends on T_CMB, and also maybe Neff?
      ! TODO: Should there be a factor of Neff/N (~3.046/3)^(3/4) here (converts 94.14 -> 93.14 eV)?
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: C

      C = (14./11.)*180.*Riemann_3/(7.*pi**4)
      C = C*SBconst*cosm%T_CMB**3/(kB*critical_density*c_light)
      C = C*eV/c_light**2  
      C = 1./C
      !C = C*(cosm%neff/3.)**0.75 ! This converts 94.1 -> 93.14 eV
      neutrino_constant = C

   END FUNCTION neutrino_constant

   REAL FUNCTION Komatsu_nu(y)

      ! Equation (26) in Komatsu et al. (2011; https://arxiv.org/pdf/1001.4538.pdf)
      ! When f(y->0) = 1, f(y->infinity) = Ay
      REAL, INTENT(IN) :: y
      REAL, PARAMETER :: A = A_Komatsu
      REAL, PARAMETER :: p = p_Komatsu

      Komatsu_nu = (1.+(A*y)**p)**(1./p)

   END FUNCTION Komatsu_nu
   
   REAL FUNCTION xi_lin(r, a, flag, cosm)

      ! Computes the 3D linear matter correlation function by integrating over P(k)
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i
      REAL :: k1, k2, xi_bit, xi
      INTEGER, PARAMETER :: min_humps = min_humps_xi
      INTEGER, PARAMETER :: max_humps = max_humps_xi
      INTEGER, PARAMETER :: method = method_xi
      INTEGER, PARAMETER :: iorder = iorder_xi

      STOP 'XI_LIN: This is ridiculuously slow for large R'

      IF (method == 1) THEN

         xi_lin = integrate_cosm(0., 1., xi_integrand_transformed, r, a, flag, cosm, acc_cosm, iorder)

      ELSE IF (method == 2) THEN   

         ! Loop over humps
         xi = 0. ! Set summation variable to zero
         DO i = 0, max_humps

            k1 = i*pi/r
            k2 = (i+1)*pi/r

            xi_bit = integrate_cosm(k1, k2, xi_integrand, r, a, flag, cosm, acc_cosm, iorder)

            xi = xi+xi_bit

            IF (i > min_humps) THEN
               IF (abs(xi_bit/xi) < acc_cosm) THEN
                  EXIT
               END IF
            END IF

            IF (i == max_humps) THEN
               WRITE (*, *) 'XI_LIN: r [Mpc/h]:', r
               WRITE (*, *) 'XI_LIN: Minimum number of humps:', min_humps
               WRITE (*, *) 'XI_LIN: Maximum number of humps:', max_humps
               STOP 'XI_LIN: Error, maximum number of humps exceeded'
            END IF

         END DO

         xi_lin = xi

      ELSE

         STOP 'XI_LIN: Error, method specified incorrectly'

      END IF

   END FUNCTION xi_lin

   REAL FUNCTION xi_integrand(k, r, a, flag, cosm)

      ! Integrand for the 3D linear matter correlation function
      USE special_functions
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc] (integrated over)
      REAL, INTENT(IN) :: r ! Separation [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

      IF (k == 0.) THEN
         xi_integrand = 0.
      ELSE
         xi_integrand = sinc(k*r)*plin(k, a, flag, cosm)/k
      END IF

   END FUNCTION xi_integrand

   REAL FUNCTION xi_integrand_transformed(t, r, a, flag, cosm)

      ! Integrand for the 3D linear matter correlation function
      ! TODO: Optimise alpha(r)
      USE special_functions
      REAL, INTENT(IN) :: t ! kr=(-1+1/t)^alpha (integrated over 0:1)
      REAL, INTENT(IN) :: r ! Separation [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL :: kr, k, alpha
      REAL, PARAMETER :: alpha_lo = alpha_lo_xi
      REAL, PARAMETER :: alpha_hi = alpha_hi_xi
      REAL, PARAMETER :: rsplit = rsplit_xi

      IF (r < rsplit) THEN
         alpha = alpha_lo
      ELSE
         alpha = alpha_hi
      END IF

      IF (t == 0. .OR. t == 1.) THEN
         xi_integrand_transformed = 0.
      ELSE
         kr = (-1.+1./t)**alpha
         k = kr/r
         xi_integrand_transformed = sinc(kr)*plin(k, a, flag, cosm)*alpha/(t*(1.-t))
      END IF

   END FUNCTION xi_integrand_transformed

   SUBROUTINE normalise_power(cosm)

      ! Normalise the power spectrum using whatever scheme you choose to do this with
      TYPE(cosmology), INTENT(INOUT) :: cosm

      ! Check that transfer function is okay for massive neutrinos
      IF ((cosm%m_nu /= 0.) .AND. is_in_array(cosm%iTk, [iTk_none, iTk_DEFW, iTk_EH, iTk_nw])) THEN
         STOP 'INIT_COSMOLOGY: You cannot use a linear power fitting function for massive neutrino cosmologies'
      END IF

      ! Check that transfer function is okay for modified gravity
      IF ((cosm%img .NE. img_none) .AND. cosm%iTk == iTk_CAMB) THEN
         STOP 'INIT_COSMOLOGY: Modified gravity not compatible with using a CAMB transfer function'
      END IF

      ! Get the CAMB power if necessary
      ! TODO: Should all init_power stuff go here or in init_cosmology? It should all be kept together.
      IF (cosm%iTk == iTk_CAMB) THEN
         CALL init_CAMB_linear(cosm)
      ELSE IF (cosm%iTk == iTk_external) THEN
         CALL init_external_linear(cosm)
      ELSE IF (cosm%img == img_fR .OR. cosm%img == img_fR_lin) THEN
         CALL init_fR_linear(cosm)
      ELSE IF (cosm%analytical_Tk .AND. interp_all_power) THEN
         CALL init_analytical_linear(cosm)
      END IF

      ! Change the flag *before* doing the normalisation calculation because they call power
      cosm%is_normalised = .TRUE.

      ! Normalise the linear spectrum
      IF (cosm%norm_method == norm_none) THEN
         ! No need to do anything
      ELSE IF ((cosm%iTk == iTk_external) .AND. (cosm%norm_method /= norm_none)) THEN
         ! TODO: External linear could be renormalised here if necessary
         STOP 'NORMALISE_POWER: Error, external power should not be re-normalised'
      ELSE IF (cosm%norm_method == norm_sigma8) THEN
         CALL normalise_power_sigma8(cosm)
      ELSE IF (cosm%norm_method == norm_pval) THEN
         CALL normalise_power_pval(cosm)
      ELSE IF (cosm%norm_method == norm_As) THEN
         CALL normalise_power_As(cosm) 
      ELSE
         STOP 'NORMALISE_POWER: Error, normalisation method not specified correctly'
      END IF

      ! If normalisation is not done via sigma8 then calculate the correct sigma8
      IF (cosm%norm_method .NE. norm_sigma8) CALL reset_sigma8(cosm)
      IF (cosm%norm_method .NE. norm_pval)   CALL reset_pval(cosm)
      IF (cosm%norm_method .NE. norm_As)     CALL reset_As(cosm) ! TODO: Make this work
      !WRITE(*, fmt='(A5, F10.5)') 'As:', cosm%As/1e-9 ! TODO: Remove
      !STOP ! TODO: Remove

   END SUBROUTINE normalise_power

   SUBROUTINE normalise_power_sigma8(cosm)

      ! Normalising the power spectrum using sigma8
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sigma8_initial, sigma8_final, kbox_save

      IF (cosm%verbose) WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Normalising power to get correct sigma_8'

      ! Remove the k-cut for the normalisation     
      IF (cosm%box) THEN
         kbox_save = cosm%kbox
         cosm%kbox = 0.
      ELSE
         kbox_save = 0. ! Need to give this a value otherwise get a warning in debug mode
      END IF

      ! Calculate the initial sigma_8 value (will not be correct)
      sigma8_initial = sigma8_norm(cosm)
      IF (cosm%verbose) WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Initial sigma_8:', real(sigma8_initial)

      ! Normalisation factor 
      cosm%A = cosm%A*cosm%sig8/sigma8_initial

      ! Replace the k-cut if necessary
      IF (cosm%box) THEN
         cosm%kbox = kbox_save
      END IF

      ! Check that the normalisation has been done correctly
      sigma8_final = sigma8_norm(cosm)

      ! Write to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Target sigma_8:', real(cosm%sig8)
         WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Final sigma_8 (calculated):', real(sigma8_final)
         WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE normalise_power_sigma8

   SUBROUTINE normalise_power_pval(cosm)

      ! Normalise the power spectrum by fixing the power at some wavenumber at a=1
      TYPE(cosmology), INTENT(INOUT) :: cosm
   
      cosm%A = cosm%A*sqrt(cosm%pval/pval_norm(cosm))

   END SUBROUTINE normalise_power_pval

   SUBROUTINE normalise_power_As(cosm)

      ! Normalise the power spectrum by fixing As
      ! Only do this for analytical transfer functions
      ! Otherwise do nothing because the normalisation will already be correct
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%analytical_Tk) THEN
         cosm%A = cosm%A*sqrt(cosm%As/As_norm(cosm))
      ELSE
         cosm%A = 1.
      END IF

   END SUBROUTINE normalise_power_As

   SUBROUTINE reset_sigma8(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      cosm%sig8 = sigma8_norm(cosm)

   END SUBROUTINE reset_sigma8

   SUBROUTINE reset_pval(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      cosm%pval = pval_norm(cosm)

   END SUBROUTINE reset_pval

   SUBROUTINE reset_As(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%analytical_Tk) THEN
         cosm%As = As_norm(cosm)       
      ELSE
         cosm%As = cosm%As*cosm%A**2
      END IF

   END SUBROUTINE reset_As

   RECURSIVE REAL FUNCTION sigma8_norm(cosm)

      ! Calculate the value of sigma8 from the linear power spectrum
      ! TODO: Needs to call sigma_intergral, rather than sigma, not sure why, maybe regression?
      TYPE(cosmology), INTENT(INOUT) :: cosm   ! Cosmology
      REAL, PARAMETER :: R = 8.                ! Because we are doing sigma(R = 8 Mpc/h) normalisation
      REAL, PARAMETER :: a = 1.                ! Because we are doing sigma(R = 8 Mpc/h, a = 1) normalisation
      INTEGER, PARAMETER :: flag = flag_matter ! sigma8 is defined for linear matter power, not cold matter

      sigma8_norm = sigma_integral(R, a, flag, cosm)
      !sigma8 = sigma(R, a, flag_matter, cosm) ! Why not? 

   END FUNCTION sigma8_norm

   RECURSIVE REAL FUNCTION pval_norm(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm   ! Cosmology
      REAL, PARAMETER :: a = 1.                ! Because we are doing sigma(R = 8 Mpc/h, a = 1) normalisation
      INTEGER, PARAMETER :: flag = flag_matter ! sigma8 is defined for linear matter power, not cold matter

      pval_norm = plin(cosm%kval, a, flag, cosm)

   END FUNCTION pval_norm

   RECURSIVE REAL FUNCTION As_norm(cosm)

      ! Calculate the value of A_s from the linear power spectrum; defined using kpiv
      ! See equation (10) of https://arxiv.org/pdf/1807.00040.pdf
      TYPE(cosmology), INTENT(INOUT) :: cosm   ! Cosmology
      REAL :: kpiv, Tk, g
      REAL, PARAMETER :: a = 1.                ! Normalisation is at a=1
      INTEGER, PARAMETER :: flag = flag_matter ! A_s is defined from linear total matter

      kpiv = cosm%kpiv              ! Pivot wavenumber: As = As(kp)
      Tk = Tk_matter(kpiv, a, cosm) ! Matter transfer function at a = 1
      g = ungrow(a, cosm)           ! Growth factor at a=1 normalised such that g(a<<1) = a
      As_norm = (25./4.)*((cosm%Om_m/g)**2)*((kpiv*Hdist)**(-4))*plin(kpiv, a, flag, cosm)/Tk**2

   END FUNCTION As_norm

   SUBROUTINE init_fR_linear(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: dinit, vinit
      REAL, ALLOCATABLE :: k(:), gk(:, :), Pk(:, :)
      REAL, ALLOCATABLE :: d(:), v(:)
      REAL, ALLOCATABLE :: a_ode(:), a_lin(:)
      INTEGER :: ik, ia
      REAL :: norm, f, g
      REAL, PARAMETER :: kmin = kmin_plin
      REAL, PARAMETER :: kmax = kmax_plin
      INTEGER, PARAMETER :: nk = nk_plin
      REAL, PARAMETER :: amin_lin = amin_plin
      REAL, PARAMETER :: amax_lin = amax_plin
      INTEGER, PARAMETER :: na_lin = na_plin
      REAL, PARAMETER :: aini_ode = aini_growth
      REAL, PARAMETER :: afin_ode = amax_plin
      INTEGER, PARAMETER :: na_ode = na_growth
      INTEGER, PARAMETER :: imeth_ode = imeth_ODE_growth
      REAL, PARAMETER :: acc_ode = acc_ODE_growth
      LOGICAL, PARAMETER :: ilog_k = .TRUE.
      LOGICAL, PARAMETER :: ilog_a = .TRUE.

      IF(cosm%verbose) WRITE (*, *) 'INIT_FR_LINEAR: Starting'

      CALL fill_array(kmin, kmax, k, nk, ilog=ilog_k)
      ALLOCATE (gk(nk, na_ode))

      ! Set the initial conditions to be in the cold matter growing mode
      ! NOTE: For massive neutrinos or EDE there is no asymptotic g(a) ~ a limit
      IF (EDE_growth_ics) THEN
         IF (cold_growth) THEN
            f = 1.-Omega_cold_norad(aini_ode, cosm)+cosm%f_nu
         ELSE
            f = 1.-Omega_m_norad(aini_ode, cosm)
         END IF
      ELSE
         IF (cold_growth) THEN
            f = cosm%f_nu
         ELSE
            f = 0.
         END IF
      END IF
      dinit = aini_ode**(1.-3.*f/5.)
      vinit = (1.-3.*f/5.)*aini_ode**(-3.*f/5.)

      ! Write to screen
      IF(cosm%verbose) THEN
         WRITE (*, *) 'INIT_FR_LINEAR: Solving scale-dependent growth ODE'
         WRITE (*, *) 'INIT_FR_LINEAR: kmin [h/Mpc]:', kmin
         WRITE (*, *) 'INIT_FR_LINEAR: kmax [h/Mpc]:', kmax
         WRITE (*, *) 'INIT_FR_LINEAR: nk:', nk
         WRITE (*, *) 'INIT_FR_LINEAR: Linear power amin:', amin_lin
         WRITE (*, *) 'INIT_FR_LINEAR: Linear power amax:', amax_lin
         WRITE (*, *) 'INIT_FR_LINEAR: Linear power na:', na_lin
         WRITE (*, *) 'INIT_FR_LINEAR: Scale-growth amin:', aini_ode
         WRITE (*, *) 'INIT_FR_LINEAR: Scale-growth amax:', afin_ode
         WRITE (*, *) 'INIT_FR_LINEAR: Scale-growth na:', na_ode
      END IF

      ! Loop over all wavenumbers and solve g(k) equation
      DO ik = 1, nk
         CALL ODE2_adaptive_cosm(d, v, k(ik), a_ode, cosm, &
            aini_ode, afin_ode, &
            dinit, vinit, &
            dvda, &
            acc_ode, na_ode, imeth_ode, ilog=ilog_a)
         gk(ik, :) = d
      END DO

      ! Get normalisation for g(k) 
      ! NOTE: Only will work if init_growth has run, so call ungrow rather than cosm%gnorm, norm = gk(1, na_ode) is lazy too
      norm = ungrow(1., cosm)
      gk = gk/norm
      DO ia = 1, na_ode     
         gk(:, ia) = gk(:, ia)/grow(a_ode(ia), cosm) ! Isolate the pure scale-dependent part of the growth, as g(a) already in P_lin(k)
      END DO

      ! Write to screen
      IF(cosm%verbose) THEN
         WRITE (*, *) 'INIT_FR_LINEAR: ODE solved'
         WRITE (*, *) 'INIT_FR_LINEAR: Calculating power'
      END IF

      ! Need to set is_normalised now so that linear power (from EH) can be got below
      cosm%is_normalised = .TRUE.

      ! Get the linear power, will not be correct at this stage
      CALL fill_array(amin_lin, amax_lin, a_lin, na_lin, ilog=ilog_a)
      CALL calculate_Plin(k, a_lin, Pk, flag_matter, cosm)

      ! Mutliply through by scale-dependent growth squared to get f(R) linear shape
      DO ia = 1, na_lin
         DO ik = 1, nk
            g = find(a_lin(ia), a_ode, gk(ik, :), na_ode, iorder=3, ifind=ifind_split, iinterp=iinterp_Lagrange)
            Pk(ik, ia) = Pk(ik, ia)*g**2
         END DO
      END DO
      CALL init_linear(k, a_lin, Pk, cosm)

      ! Write to screen
      IF(cosm%verbose) THEN
         WRITE (*, *) 'INIT_FR_LINEAR: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_fR_linear

   REAL FUNCTION comoving_critical_density(a, cosm)

      ! Comoving critical density [(Msun/h) / (Mpc/h)^3]
      ! For LCDM, ignoring radiation, this is constant in the past, increases like a^3 in the future
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      comoving_critical_density = physical_critical_density(a, cosm)*a**3

   END FUNCTION comoving_critical_density

   REAL FUNCTION physical_critical_density(a, cosm)

      ! Physical critical density [(Msun/h) / (Mpc/h)^3]
      ! For LCDM, ignoring radiation, tends to a constant in the future, behaves like a^-3 in the past
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_critical_density = critical_density_cos*Hubble2(a, cosm)

   END FUNCTION physical_critical_density

   PURE REAL FUNCTION comoving_matter_density(cosm)

      ! Comoving matter density [(Msun/h) / (Mpc/h)^3]
      ! Not a function of redshift, constant value throughout time
      TYPE(cosmology), INTENT(IN) :: cosm

      comoving_matter_density = critical_density_cos*cosm%Om_m

   END FUNCTION comoving_matter_density

   PURE REAL FUNCTION physical_matter_density(a, cosm)

      ! Physical matter density [(Msun/h) / (Mpc/h)^3]
      ! Proportional to a^-3 always
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm

      physical_matter_density = comoving_matter_density(cosm)*a**(-3)

   END FUNCTION physical_matter_density

   REAL FUNCTION Hubble2(a, cosm)

      ! Calculates Hubble^2 in units such that H^2(a=1)=1
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'HUBBLE2: Error, cosmology is not initialised'
      Hubble2 = &
         cosm%Om_c*X_c(a)+ &
         cosm%Om_b*X_b(a)+ &
         cosm%Om_g*X_r(a)+ &
         cosm%Om_nu*X_nu(a, cosm)+ &
         cosm%Om_v*X_v(a)+ &
         cosm%Om_w*X_de(a, cosm)+ &
         (1.-cosm%Om)*a**(-2)

   END FUNCTION Hubble2

   REAL FUNCTION Hubble2_norad(a, cosm)

      ! Squared Hubble parameter but ignoring the radiation component
      ! Units such that Hubble^2(a=1)=1 (as long as radiation is not important at a=1)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (a > cosm%a_nu) THEN
         Hubble2_norad = Hubble2(a, cosm)-cosm%Om_g*a**(-4)
      ELSE
         Hubble2_norad = Hubble2(a, cosm)-cosm%Om_r*a**(-4)
      END IF

   END FUNCTION Hubble2_norad

   REAL FUNCTION AH(a, cosm)

      ! Acceleration function: \ddot{a}/a
      ! Some people call this the Raychaudhuri equation
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'AH: Error, cosmology is not initialised'
      AH = &
         cosm%Om_c*(1.+3.*w_c)*X_c(a)+ &
         cosm%Om_b*(1.+3.*w_b)*X_b(a)+ &
         cosm%Om_g*(1.+3.*w_g)*X_g(a)+ &
         cosm%Om_nu*(1.+3.*w_nu(a, cosm))*X_nu(a, cosm)+ &
         cosm%Om_v*(1.+3.*w_v)*X_v(a)+ &
         cosm%Om_w*(1.+3.*w_de(a, cosm))*X_de(a, cosm)
      AH = -AH/2.

   END FUNCTION AH

   REAL FUNCTION AH_norad(a, cosm)

      ! Acceleration function without the radiation contribution
      ! NOTE: It is correct that we should add Om_r*a^-4 here because of the -(1/2) factor in AH
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (a > cosm%a_nu) THEN
         AH_norad = AH(a, cosm)+cosm%Om_g*a**(-4)
      ELSE
         AH_norad = AH(a, cosm)+cosm%Om_r*a**(-4)
      END IF

   END FUNCTION AH_norad

   REAL FUNCTION Omega_m(a, cosm)

      ! This calculates Omega_m variations with scale factor (note this is not proportional to a^-3 always)
      ! TODO: Maybe retire eventually, since Omega_m is not a real thing due to neutrinos
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_M: Error, cosmology is not initialised'
      Omega_m = cosm%Om_m*X_m(a)/Hubble2(a, cosm)

   END FUNCTION Omega_m

   REAL FUNCTION Omega_m_norad(a, cosm)

      ! This calculates Omega_m variations with scale factor, but ignoring photon contribution
      ! This ensures that Omega_m_norad(a->0) -> 1
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_M_NORAD: Error, cosmology is not initialised'
      IF (a > cosm%a_nu) THEN
         Omega_m_norad = cosm%Om_m*X_m(a)/Hubble2_norad(a, cosm)
      ELSE
         Omega_m_norad = (cosm%Om_c+cosm%Om_b)*X_m(a)/Hubble2_norad(a, cosm)
      END IF

   END FUNCTION Omega_m_norad

   REAL FUNCTION Omega_c(a, cosm)

      ! This calculates Omega_c variations with scale factor (note this is not proportional to a^-3 always)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_C: Error, cosmology is not initialised'
      Omega_c = cosm%Om_c*X_c(a)/Hubble2(a, cosm)

   END FUNCTION Omega_c

   REAL FUNCTION Omega_b(a, cosm)

      ! This calculates Omega_b variations with scale factor (note this is not proportional to a^-3 always)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_B: Error, cosmology is not initialised'
      Omega_b = cosm%Om_b*X_b(a)/Hubble2(a, cosm)

   END FUNCTION Omega_b

   REAL FUNCTION Omega_cold_norad(a, cosm)

      ! This calculates Omega_c variations with scale factor, but ignoring photon and neutrino components
      ! This ensures that Omega_m_norad(a->0) -> 1
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_COLD_NORAD: Error, cosmology is not initialised'
      Omega_cold_norad = (cosm%Om_c*X_c(a)+cosm%Om_b*X_b(a))/Hubble2_norad(a, cosm)

   END FUNCTION Omega_cold_norad

   REAL FUNCTION Omega_r(a, cosm)

      ! This calculates Omega_r variations with scale factor (note this is *not* proportional to a^-4 always)
      ! TODO: Maybe retire eventually, since Omega_r is not a real thing due to neutrinos
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_R: Error, cosmology is not initialised'
      Omega_r = cosm%Om_r*X_r(a)/Hubble2(a, cosm)

   END FUNCTION Omega_r

   REAL FUNCTION Omega_g(a, cosm)

      ! This calculates Omega_g variations with scale factor (this is proportional to a^-4 always)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_G: Error, cosmology is not initialised'
      Omega_g = cosm%Om_g*X_g(a)/Hubble2(a, cosm)

   END FUNCTION Omega_g

   REAL FUNCTION Omega_nu(a, cosm)

      ! This calculates Omega_nu variations with scale factor (changes between a^-4 and a^-3)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_NU: Error, cosmology is not initialised'
      Omega_nu = cosm%Om_nu*X_nu(a, cosm)/Hubble2(a, cosm)

   END FUNCTION Omega_nu

   REAL FUNCTION Omega_v(a, cosm)

      ! This calculates Omega_v variations with scale factor (note this is not constant)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_V: Error, cosmology is not initialised'
      Omega_v = cosm%Om_v*X_v(a)/Hubble2(a, cosm)

   END FUNCTION Omega_v

   REAL FUNCTION Omega_w(a, cosm)

      ! This calculates Omega_w variations with scale factor
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_W: Error, cosmology is not initialised'
      Omega_w = cosm%Om_w*X_de(a, cosm)/Hubble2(a, cosm)

   END FUNCTION Omega_w

   REAL FUNCTION Omega(a, cosm)

      ! This calculates total Omega variations with a
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA: Error, cosmology is not initialised'
      Omega = &
         Omega_c(a, cosm)+ &
         Omega_b(a, cosm)+ &
         Omega_g(a, cosm)+ &
         Omega_nu(a, cosm)+ &
         Omega_v(a, cosm)+ &
         Omega_w(a, cosm)

   END FUNCTION Omega

   REAL FUNCTION w_nu(a, cosm)

      ! Scaling of neutrino density
      ! TODO: Account for radiation -> matter transition properly
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: x, y
      REAL, PARAMETER :: bigA = A_Komatsu
      REAL, PARAMETER :: p = p_Komatsu

      IF (neutrino_method == neutrino_basic) THEN
         IF (a > cosm%a_nu) THEN
            w_nu = 0.
         ELSE
            w_nu = 1./3.
         END IF
      ELSE IF (neutrino_method == neutrino_Komatsu) THEN
         y = a/cosm%a_nu
         x = (bigA*y)**p
         w_nu = (1.-x/(1.+x))/3.
      ELSE
         STOP 'W_NU: Error, neutrino method not recognised'
      END IF

   END FUNCTION w_nu

   REAL FUNCTION w_de(a, cosm)

      ! Variations of the dark energy equation-of-state parameter w(a)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: p1, p2, p3, p4
      REAL :: f1, f2, f3, f4
      REAL :: z

      IF (cosm%iw == iw_LCDM) THEN
         ! LCDM
         w_de = -1.
      ELSE IF (cosm%iw == iw_QUICC) THEN
         ! QUICC parameterisation
         p1 = 1.+exp(cosm%am/cosm%dm)
         p2 = 1.-exp(-(a-1.)/cosm%dm)
         p3 = 1.+exp(-(a-cosm%am)/cosm%dm)
         p4 = 1.-exp(1./cosm%dm)
         w_de = cosm%w+(cosm%wm-cosm%w)*p1*p2/(p3*p4)
      ELSE IF (cosm%iw == iw_waCDM) THEN
         ! w(a)CDM
         w_de = cosm%w+(1.-a)*cosm%wa
      ELSE IF (cosm%iw == iw_wCDM) THEN
         ! wCDM
         w_de = cosm%w
      ELSE IF (cosm%iw == iw_IDE1) THEN
         ! IDE I
         w_de = ((a/cosm%astar)**cosm%nstar-1.)/((a/cosm%astar)**cosm%nstar+1.)
      ELSE IF (cosm%iw == iw_IDE2) THEN
         ! IDE II
         f1 = a**cosm%nstar-cosm%a1n
         f2 = a**cosm%nstar+cosm%a1n
         f3 = a**cosm%nstar-cosm%a2n
         f4 = a**cosm%nstar+cosm%a2n
         w_de = -1.+real(f1/f2-f3/f4)
      ELSE IF (cosm%iw == iw_IDE3) THEN
         ! IDE III
         IF (a < cosm%a1) THEN
            w_de = -1.
         ELSE IF (cosm%a1 <= a .AND. a < cosm%a2) THEN
            w_de = cosm%ws
         ELSE IF (a >= cosm%a2) THEN
            w_de = -1.
         ELSE
            STOP 'W_DE: Error, something went wrong'
         END IF
      ELSE IF (cosm%iw == iw_BDE) THEN
         ! Bound dark energy
         z = redshift_a(a)
         w_de = 0.
         w_de = w_de + &
            cosm%b0*z**0 + &
            cosm%b1*z**1 + &
            cosm%b2*z**2 + &
            cosm%b3*z**3 + &
            cosm%b4*z**4
         w_de = w_de*a**4
      ELSE
         STOP 'W_DE: Error, value of iw set incorrectly'
      END IF

   END FUNCTION w_de

   REAL FUNCTION w_de_total(a, cosm)

      ! Average equation-of-state over the dark energy components (vacuum and w currently)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%Om_v == 0. .AND. cosm%Om_w == 0.) THEN
         w_de_total = -1.
      ELSE
         w_de_total = w_de(a, cosm)*Omega_w(a, cosm)-Omega_v(a, cosm)
         w_de_total = w_de_total/(Omega_w(a, cosm)+Omega_v(a, cosm))
      END IF

   END FUNCTION w_de_total

   REAL FUNCTION w_eff(a, cosm)

      ! Average equation-of-state over all components
      ! Note that matter does not enter in this equation because w_matter = 0
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      w_eff = &
         w_de(a, cosm)*Omega_w(a, cosm)+ &
         w_v*Omega_v(a, cosm)+ &
         w_g*Omega_g(a, cosm)+ &
         w_nu(a, cosm)*Omega_nu(a, cosm)
      w_eff = w_eff/Omega(a, cosm)

   END FUNCTION w_eff

   REAL FUNCTION X_m(a)

      ! Scaling of matter density
      ! TODO: Retire because matter is not really a thing due to neutrinos
      REAL, INTENT(IN) :: a

      X_m = a**(-3)

   END FUNCTION X_m

   REAL FUNCTION X_c(a)

      ! Scaling of CDM density
      REAL, INTENT(IN) :: a

      X_c = a**(-3)

   END FUNCTION X_c

   REAL FUNCTION X_b(a)

      ! Scaling of baryon density
      REAL, INTENT(IN) :: a

      X_b = a**(-3)

   END FUNCTION X_b

   REAL FUNCTION X_r(a)

      ! Scaling of radiation density
      ! TODO: Retire because radiation is not really a thing due to neutrinos
      REAL, INTENT(IN) :: a

      X_r = a**(-4)

   END FUNCTION X_r

   REAL FUNCTION X_g(a)

      ! Scaling of photon density
      REAL, INTENT(IN) :: a

      X_g = a**(-4)

   END FUNCTION X_g

   REAL FUNCTION X_nu(a, cosm)

      ! Scaling of neutrino density
      ! TODO: Account for radiation -> matter transition properly
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (neutrino_method == neutrino_basic) THEN
         IF (a > cosm%a_nu) THEN
            X_nu = a**(-3)
         ELSE
            X_nu = (cosm%Om_nu_rad/cosm%Om_nu)/a**4
         END IF
      ELSE IF (neutrino_method == neutrino_Komatsu) THEN
         X_nu = (cosm%Om_nu_rad/cosm%Om_nu)*Komatsu_nu(a/cosm%a_nu)/a**4
      ELSE
         STOP 'X_NU: Error, neutrino method not recognised'
      END IF

   END FUNCTION X_nu

   REAL FUNCTION X_v(a)

      ! Scaling of vacuum density
      REAL, INTENT(IN) :: a
      REAL :: crap

      crap = a

      X_v = 1.

   END FUNCTION X_v

   REAL FUNCTION X_de(a, cosm)

      ! Scaling for dark energy density (i.e., if w=0 X(a)=a^-3, if w=-1 X(a)=const etc.)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: f1, f2, f3, f4

      IF (cosm%iw == iw_LCDM) THEN
         ! LCDM
         X_de = 1.
      ELSE IF (cosm%iw == iw_waCDM) THEN
         ! w(a)CDM
         X_de = (a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))
      ELSE IF (cosm%iw == iw_wCDM) THEN
         ! wCDM
         X_de = a**(-3.*(1.+cosm%w))
      ELSE IF (cosm%iw == iw_IDE1) THEN
         ! IDE I
         X_de = ((1.+(a/cosm%astar)**cosm%nstar)/(1.+(1./cosm%astar)**cosm%nstar))**(-6./cosm%nstar)
      ELSE IF (cosm%iw == iw_IDE2) THEN
         ! IDE II
         f1 = a**cosm%nstar+cosm%a1n
         f2 = 1.+cosm%a1n
         f3 = 1.+cosm%a2n
         f4 = a**cosm%nstar+cosm%a2n
         X_de = real(f1*f3/(f2*f4))**(-6./cosm%nstar)
      ELSE IF (cosm%iw == iw_IDE3) THEN
         ! IDE III
         IF (a < cosm%a1) THEN
            X_de = (cosm%a1/cosm%a2)**(-3.*(1.+cosm%ws))
         ELSE IF (cosm%a1 <= a .AND. a < cosm%a2) THEN
            X_de = (a/cosm%a2)**(-3.*(1.+cosm%ws))
         ELSE IF (a >= cosm%a2) THEN
            X_de = 1.
         ELSE
            STOP 'X_DE: Error, something went wrong'
         END IF
      ELSE
         ! Generally true, but this integration can make calculations very slow
         IF (tabulate_Xde) THEN
            IF(.NOT. cosm%has_Xde) CALL init_Xde(cosm)
            X_de = evaluate_interpolator(a, cosm%Xde)
         ELSE
            X_de = Xde_integral(a, cosm)
         END IF
      END IF

   END FUNCTION X_de

   SUBROUTINE init_Xde(cosm)

      ! Initialise an interpolator for X_de(a)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: amin = amin_Xde
      REAL, PARAMETER :: amax = amax_Xde
      INTEGER, PARAMETER :: n = n_Xde
      INTEGER :: i
      REAL, ALLOCATABLE :: a(:), Xde(:)

      CALL fill_array(amin, amax, a, n)
      ALLOCATE(Xde(n))

      IF (cosm%verbose) THEN
         WRITE(*, *) 'INIT_XDE: minimum a:', amin
         WRITE(*, *) 'INIT_XDE: maximum a:', amax
         WRITE(*, *) 'INIT_XDE: number of a:', n
      END IF

      DO i = 1, n
         Xde(i) = Xde_integral(a(i), cosm)
      END DO

      IF (cosm%verbose) THEN
         WRITE(*, *) 'INIT_XDE: minimum X_de:', Xde(n)
         WRITE(*, *) 'INIT_XDE: maximum X_de:', Xde(1)
      END IF

      CALL init_interpolator(a, Xde, cosm%Xde, &
         iorder = iorder_interp_Xde, &
         iextrap = iextrap_Xde, &
         store = store_Xde, &
         logx = .TRUE., &
         logf = .TRUE.)
      cosm%has_Xde = .TRUE.

      IF (cosm%verbose) THEN
         WRITE(*, *) 'INIT_XDE: Done'
         WRITE(*, *)
      END IF

   END SUBROUTINE init_Xde

   REAL FUNCTION Xde_integral(a, cosm)

      ! Integral to calculate X_de
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: acc = acc_integration_Xde
      INTEGER, PARAMETER :: iorder = iorder_integration_Xde

      Xde_integral = (a**(-3))*exp(3.*integrate_cosm(a, 1., integrand_de, cosm, acc, iorder))

   END FUNCTION Xde_integral

   REAL FUNCTION integrand_de(a, cosm)

      ! The integrand for the X_de(a) integral
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      integrand_de = w_de(a, cosm)/a

   END FUNCTION integrand_de

   TYPE(cosmology) FUNCTION convert_cosmology(cosm, make_lambda, make_flat, remove_neutrinos)

      ! Make a vanilla LCDM version of an input cosmology
      ! This will be a flat cosmology with standard Lambda dark energy
      ! It will also have zero neutrino mass
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, INTENT(IN) :: make_lambda
      LOGICAL, INTENT(IN) :: make_flat
      LOGICAL, INTENT(IN) :: remove_neutrinos

      ! Initially set the new cosmology to the input cosmology
      convert_cosmology = cosm

      ! Remove dark energy
      IF (make_lambda) THEN
         convert_cosmology%iw = iw_LCDM
         convert_cosmology%w = -1.
         convert_cosmology%wa = 0.
         convert_cosmology%Om_w = 0.
         convert_cosmology%Om_v = cosm%Om_v+cosm%Om_w
      END IF

      ! Flatten by forcing the dark-energy or vacuum density to sum with matter to unity
      IF (make_flat) THEN
         IF (make_lambda) THEN
            convert_cosmology%Om_v = 1.-cosm%Om_m
         ELSE
            IF ((cosm%Om_v .NE. 0.) .AND. (cosm%Om_w .NE. 0.)) THEN
               STOP 'CONVERT_COSMOLOGY: Error, no unique way to flatten a cosmology with Omeega_v and Omega_w'
            ELSE IF (cosm%Om_w .NE. 0.) THEN
               convert_cosmology%Om_w = 1.-cosm%Om_m
            ELSE
               convert_cosmology%Om_v = 1.-cosm%Om_m
            END IF
         END IF
      END IF

      ! Remove neutrinos will convert nu to CDM since Omega_c is a derived parameter
      IF (remove_neutrinos) convert_cosmology%m_nu = 0.
    
      ! Ensure the new cosmology is not verbose
      convert_cosmology%verbose = .FALSE.

      ! Initialise
      CALL init_cosmology(convert_cosmology)

   END FUNCTION convert_cosmology

   ELEMENTAL REAL FUNCTION scale_factor_z(z)

      ! Scale factor corresponding to redshift 'z'
      REAL, INTENT(IN) :: z

      scale_factor_z = 1./(1.+z)

   END FUNCTION scale_factor_z

   ELEMENTAL REAL FUNCTION redshift_a(a)

      ! Redshift corresponding to scale-factor 'a'
      REAL, INTENT(IN) :: a

      redshift_a = -1.+1./a

   END FUNCTION redshift_a

   REAL FUNCTION scale_factor_r(r, cosm)

      ! Scale factor corresponding to comoving distance 'r'
      REAL, INTENT(IN) :: r
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: p

      IF (.NOT. cosm%has_distance) CALL init_distance(cosm)
      IF (r == 0.) THEN
         scale_factor_r = 1.
      ELSE
         p = cosm%horizon-r
         scale_factor_r = inverse_interpolator(p, cosm%dist)
      END IF

   END FUNCTION scale_factor_r

   REAL FUNCTION redshift_r(r, cosm)

      ! The redshift corresponding to comoving distance 'r'
      REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
      TYPE(cosmology), INTENT(INOUT) :: cosm

      redshift_r = redshift_a(scale_factor_r(r, cosm))

   END FUNCTION redshift_r

   ELEMENTAL REAL FUNCTION f_k(r, cosm)

      ! Curvature function, also comoving angular-diameter distance [Mpc/h]
      REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (cosm%k > 0.) THEN
         f_k = sin(sqrt(cosm%k)*r)/sqrt(cosm%k)
      ELSE IF (cosm%k < 0.) THEN
         f_k = sinh(sqrt(-cosm%k)*r)/sqrt(-cosm%k)
      ELSE
         f_k = r
      END IF

   END FUNCTION f_k

   ELEMENTAL REAL FUNCTION fdash_k(r, cosm)

      ! Derivative of curvature function df_k(r)/dr
      REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (cosm%k > 0.) THEN       
         fdash_k = cos(sqrt(cosm%k)*r)
      ELSE IF (cosm%k < 0.) THEN
         fdash_k = cosh(sqrt(-cosm%k)*r)
      ELSE
         fdash_k = 1.
      END IF

   END FUNCTION fdash_k

   REAL FUNCTION comoving_particle_horizon(a, cosm)

      ! The comoving particle horizon [Mpc/h]
      ! Related to the conformal time via a change in dimension
      ! This is the furthest distance a particle can have travelled since a=0
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (.NOT. cosm%has_distance) CALL init_distance(cosm)

      IF (a == 0.) THEN
         comoving_particle_horizon = 0.
      ELSE IF (a > 1.) THEN
         WRITE (*, *) 'COMOVING_PARTICLE_HORIZON: a:', a
         STOP 'COMOVING_PARTICLE_HORIZON: Error, tried to calculate particle horizon in the future'
      ELSE
         comoving_particle_horizon = evaluate_interpolator(a, cosm%dist)
      END IF

   END FUNCTION comoving_particle_horizon

   REAL FUNCTION physical_particle_horizon(a, cosm)

      ! The physical particle horizon [Mpc/h]
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_particle_horizon = comoving_particle_horizon(a, cosm)*a

   END FUNCTION physical_particle_horizon

   REAL FUNCTION comoving_distance(a, cosm)

      ! The comoving distance [Mpc/h]
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: p

      ! Ensures that init_distance is run and therefore that horizon is calculated
      p = comoving_particle_horizon(a, cosm) 

      comoving_distance = cosm%horizon-p

   END FUNCTION comoving_distance

   REAL FUNCTION physical_distance(a, cosm)

      ! The physical distance [Mpc/h]
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_distance = comoving_distance(a, cosm)*a

   END FUNCTION physical_distance

   REAL FUNCTION physical_angular_distance(a, cosm)

      ! The physical angular-diameter distance [Mpc/h]
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_angular_distance = a*comoving_angular_distance(a, cosm)

   END FUNCTION physical_angular_distance

   REAL FUNCTION comoving_angular_distance(a, cosm)

      ! The comoving angular-diameter distance [Mpc/h]
      ! Some people call this the 'effective distance'
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      comoving_angular_distance = f_k(comoving_distance(a, cosm), cosm)

   END FUNCTION comoving_angular_distance

   REAL FUNCTION luminosity_distance(a, cosm)

      ! The luminosity distance [Mpc/h]
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      luminosity_distance = f_k(comoving_distance(a, cosm), cosm)/a

   END FUNCTION luminosity_distance

   REAL FUNCTION conformal_time(a, cosm)

      ! Conformal time at scale factor a [Gyrs/h]
      ! Sometimes denoted eta, with eta(t) = int dt/a between t=0 and t
      ! Same as the particle horizon to within a factor of c
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      conformal_time = (Htime/Hdist)*comoving_particle_horizon(a, cosm)

   END FUNCTION conformal_time

   SUBROUTINE init_distance(cosm)

      ! Fill up tables of a vs. p(a) (comoving particle horizon)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: zmin, zmax, b
      REAL, ALLOCATABLE :: a(:), r(:)
      INTEGER :: i
      REAL, PARAMETER :: amin = amin_distance
      REAL, PARAMETER :: amax = amax_distance
      INTEGER, PARAMETER :: iorder = iorder_integration_distance ! Order for integration
      INTEGER, PARAMETER :: n = n_distance

      ! Calculate redshifts
      zmin = redshift_a(amax)
      zmax = redshift_a(amin)
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_DISTANCE: Redshift range for distance tables'
         WRITE (*, *) 'INIT_DISTANCE: minimum z:', real(zmin)
         WRITE (*, *) 'INIT_DISTANCE: maximum z:', real(zmax)
         WRITE (*, *) 'INIT_DISTANCE: minimum a:', real(amin)
         WRITE (*, *) 'INIT_DISTANCE: maximum a:', real(amax)
      END IF

      ! Fill array of 'a' in log space
      CALL fill_array_log(amin, amax, a, n)
      ALLOCATE(r(n))

      ! Now do the r(a) calculation
      DO i = 1, n
         b = sqrt(a(i)) ! Parameter to make the integrand not diverge for small values (a=b^2)
         r(i) = integrate_cosm(0., b, distance_integrand, cosm, acc_cosm, iorder)
      END DO
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_DISTANCE: minimum r [Mpc/h]:', real(r(1))
         WRITE (*, *) 'INIT_DISTANCE: maximum r [Mpc/h]:', real(r(n))
      END IF

      CALL init_interpolator(a, r, cosm%dist, &
         iextrap = iextrap_distance, &
         iorder = iorder_interp_distance, &
         store = store_distance, &
         logx = .TRUE., &
         logf = .TRUE. &
         )

      ! Find the horizon distance in your cosmology
      ! exp(log) ensures the value is exactly the same as what comes out of the (log) interpolator
      cosm%horizon = exp(log(integrate_cosm(0., 1., distance_integrand, cosm, acc_cosm, iorder)))
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_DISTANCE: Horizon distance [Mpc/h]:', real(cosm%horizon)
         WRITE (*, *) 'INIT_DISTANCE: Done'
         WRITE (*, *)
      END IF

      cosm%has_distance = .TRUE.

   END SUBROUTINE init_distance

   REAL FUNCTION distance_integrand(b, cosm)

      ! Integrand for the cosmic-distance calculation [Mpc/h]
      ! Cast in terms of b=sqrt(a) to remove integrand divergence at a=0
      ! This means that the integrand is 2b/H(a)a^2, rather than 1/H(a)a^2
      REAL, INTENT(IN) :: b
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: a

      ! Relation beween a and b
      a = b**2

      IF (a < atay_distance) THEN
         IF (cosm%Om_r == 0.) THEN
            distance_integrand = 2.*Hdist/sqrt(cosm%Om_m)
         ELSE
            distance_integrand = 2.*Hdist*b*(1.-0.5*(cosm%Om_m/cosm%Om_r)*b**2)/sqrt(cosm%Om_r)
         END IF
      ELSE
         distance_integrand = 2.*Hdist*b/(sqrt(Hubble2(a, cosm))*a**2)
      END IF

   END FUNCTION distance_integrand

   REAL FUNCTION cosmic_time(a, cosm)

      ! Age of the universe [Gyr/h]
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (.NOT. cosm%has_time) CALL init_time(cosm)

      IF (a == 0.) THEN
         cosmic_time = 0.
      ELSE
         cosmic_time = evaluate_interpolator(a, cosm%time)
      END IF

   END FUNCTION cosmic_time

   REAL FUNCTION look_back_time(a, cosm)

      ! The time in the past corresponding to scale factor 'a' [Gyr/h]
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: t

      ! Ensures that init_time is run and therefore that age is calculated
      t = cosmic_time(a, cosm) 
      
      look_back_time = cosm%age-t

   END FUNCTION look_back_time

   SUBROUTINE init_time(cosm)

      ! Initialises interpolator of a vs. r(a) (comoving particle horizon)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: zmin, zmax
      REAL, ALLOCATABLE :: a(:), t(:)
      INTEGER :: i
      REAL, PARAMETER :: amin = amin_time
      REAL, PARAMETER :: amax = amax_time
      INTEGER, PARAMETER :: iorder = iorder_integration_time
      INTEGER, PARAMETER :: n = n_time

      ! Calculate redshifts
      zmin = redshift_a(amax)
      zmax = redshift_a(amin)
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_TIME: Redshift range for time tables'
         WRITE (*, *) 'INIT_TIME: minimum z:', real(zmin)
         WRITE (*, *) 'INIT_TIME: maximum z:', real(zmax)
         WRITE (*, *) 'INIT_TIME: minimum a:', real(amin)
         WRITE (*, *) 'INIT_TIME: maximum a:', real(amax)
      END IF

      ! Fill array of 'a' in log space
      CALL fill_array_log(amin, amax, a, n)
      ALLOCATE(t(n))

      ! Now do the t(a) calculation
      DO i = 1, n
         t(i) = integrate_cosm(0., a(i), time_integrand, cosm, acc_cosm, iorder)
      END DO
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_TIME: minimum t [Gyr/h]:', real(t(1))
         WRITE (*, *) 'INIT_TIME: maximum t [Gyr/h]:', real(t(n))
      END IF

      CALL init_interpolator(a, t, cosm%time, &
         iorder = iorder_interp_time, &
         iextrap = iextrap_time, &
         store = store_time, &
         logx = .TRUE., &
         logf = .TRUE.)

      ! Find the horizon distance in your cosmology
      ! exp(log) ensures the value is exactly the same as what comes out of the (log) look-up tables
      cosm%age = exp(log(integrate_cosm(0., 1., time_integrand, cosm, acc_cosm, iorder)))
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_TIME: Age [Gyr/h]:', real(cosm%age)
         WRITE (*, *) 'INIT_TIME: Done'
         WRITE (*, *)
      END IF

      cosm%has_time = .TRUE.

   END SUBROUTINE init_time

   REAL FUNCTION time_integrand(a, cosm)

      ! The integrand for the cosmic-distance calculation [Gyr/h]
      ! This is 1/H(a)*a
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (a < atay_time) THEN
         ! Taylor expansions
         IF (cosm%Om_r == 0.) THEN
            time_integrand = Htime*sqrt(a/cosm%Om_m)
         ELSE
            time_integrand = Htime*a*(1.-0.5*cosm%Om_m*a/cosm%Om_r)/sqrt(cosm%Om_r)
         END IF
      ELSE
         time_integrand = Htime/(a*sqrt(Hubble2(a, cosm)))
      END IF

   END FUNCTION time_integrand

   REAL FUNCTION Tk_matter(k, a, cosm)

      ! Transfer function selection
      ! Note that my definition is that T(k<<1, a) = unity for all a
      REAL, INTENT(IN) :: k               ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: a               ! Scale factor
      TYPE(cosmology), INTENT(IN) :: cosm ! Cosmology

      IF (cosm%analytical_Tk) THEN
         IF (cosm%iTk == iTk_none) THEN
            Tk_matter = 1.
         ELSE IF (cosm%iTk == iTk_EH) THEN
            Tk_matter = Tk_EH(k, cosm)
         ELSE IF (cosm%iTk == iTk_DEFW) THEN
            Tk_matter = Tk_DEFW(k, cosm)
         ELSE IF (cosm%iTk == iTk_nw) THEN
            Tk_matter = Tk_nw(k, cosm)
         ELSE
            WRITE (*, *) 'TK_MATTER: iTk:', cosm%iTk
            STOP 'TK_MATTER: Error, iTk specified incorrectly'
         END IF
      ELSE
         Tk_matter = evaluate_interpolator(k, a, cosm%Tk_matter)
      END IF

      ! Additional weirdness
      Tk_matter = Tk_matter*Tk_factor(k, cosm)

   END FUNCTION Tk_matter

   REAL FUNCTION Tk_cold(k, a, cosm)

      ! Ratio of transfer function for cold matter relative to all matter
      ! Note that my definition is that T(k<<1, a) = (Om_c+Om_b)/Om_m (not unity) for all a
      ! TODO: Force the cold transfer function to come from CAMB if possible?
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: a ! Scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

      IF ((cosm%iTc == iTc_none) .OR. (cosm%m_nu == 0.)) THEN
         ! Assuming the cold spectrum is exactly the matter spectrum
         Tk_cold = 1. 
      ELSE IF (cosm%iTc == iTc_total) THEN
         ! This approximation assumes that the neutrinos are as clustered as the rest of the mass
         ! This is only true on scales greater than the neutrino free-streaming scale
         Tk_cold = (cosm%Om_c+cosm%Om_b)/cosm%Om_m 
      ELSE IF (cosm%iTc == iTc_EH) THEN
         ! Use the Eisenstein and Hu approximation
         Tk_cold = Tk_cold_EH(k, a, cosm)
      ELSE IF (cosm%iTc == iTc_CAMB) THEN
         ! Use look-up tables from CAMB transfer functions
         Tk_cold = evaluate_interpolator(k, a, cosm%Tk_cold)
      ELSE
         STOP 'TCOLD: Error, cold transfer function method not specified correctly'
      END IF

      ! Additional weirdness
      Tk_cold = Tk_cold*Tk_factor(k, cosm)

   END FUNCTION Tk_cold

   REAL FUNCTION Tk_factor(k, cosm)

      ! Additional weird things to put in the transfer function
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm

      Tk_factor = 1.
      IF (cosm%warm) Tk_factor = Tk_factor*Tk_WDM(k, cosm) ! Damp transfer function if considering WDM
      IF (cosm%bump == 1) THEN
         Tk_factor = Tk_factor*Tk_bump(k, cosm)
      ELSE IF (cosm%bump == 2) THEN
         Tk_factor = Tk_factor*Tk_bump_Mexico(k, cosm)
      END IF

   END FUNCTION Tk_factor

   REAL FUNCTION Tk_DEFW(k, cosm)

      ! The DEFW transfer function approximation (astro-ph/xxx.xxxx)
      ! Relies on the power-spectrum scale parameter Gamma=Omega_m*h
      ! This function was written by John Peacock
      ! NOTE: I removed double precision for q8 and Tk8 from this
      ! TODO: Why the 1e-20 in the q4 and q8 definitions?
      REAL, INTENT(IN) :: k               ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm ! Cosmology
      REAL :: keff, q4, tk4, q8, tk8, Gamma

      ! Calculate shape parameter
      IF (cosm%power_Omegas) THEN
         Gamma = cosm%Om_m_pow*cosm%h_pow
      ELSE
         Gamma = cosm%Om_m*cosm%h
      END IF

      keff = 0.172+0.011*log(Gamma/0.36)*log(Gamma/0.36)
      q4 = 1.e-20+k/Gamma    
      q8 = 1.e-20+keff/Gamma
      tk4 = 1./(1.+(6.4*q4+(3.0*q4)**1.5+(1.7*q4)**2)**1.13)**(1./1.13)
      tk8 = 1./(1.+(6.4*q8+(3.0*q8)**1.5+(1.7*q8)**2)**1.13)**(1./1.13)

      tk_defw = tk4/tk8

   END FUNCTION Tk_DEFW

   REAL FUNCTION Tk_EH(k, cosm)

      USE special_functions

      ! Eisenstein & Hu fitting function (arXiv: 9709112)
      ! JP: the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
      ! JP: remember I use k/h, whereas they use pure k, Om_m is cdm + baryons
      ! TODO: Could have an init for this as many things only need to be calculated once
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: rk, e, thet, b1, b2, zd, ze, rd, re, rke, s, rks
      REAL :: q
      REAL :: y, g, ab
      REAL :: a1, a2, ac
      REAL :: bc
      REAL :: f, fac
      REAL :: c1, c2, tc
      REAL :: bb, bn, ss, tb
      REAL :: Om_m, Om_b, h, Omh2, Obh2

      ! Define some useful variables
      IF (cosm%power_Omegas) THEN
         Om_m = cosm%Om_m_pow
         Om_b = cosm%Om_b_pow
         h = cosm%h_pow
      ELSE
         Om_m = cosm%Om_m
         Om_b = cosm%Om_b
         h = cosm%h
      END IF

      ! Physical densities
      Omh2 = Om_m*h**2
      Obh2 = Om_b*h**2

      ! Convert wave-number from h/Mpc to 1/Mpc
      rk = k*h

      ! 2.718...
      e = exp(1.)

      ! CMB temperature (Section 2)
      thet = cosm%T_CMB/2.7

      b1 = 0.313*(Omh2)**(-0.419)*(1.+0.607*(Omh2)**0.674) ! Equation (4)
      b2 = 0.238*(Omh2)**0.223 ! Equation (4)
      zd = 1291.*(1.+b1*(Obh2)**b2)*(Omh2)**0.251/(1.+0.659*(Omh2)**0.828) ! Drag redshift; equation (4)
      ze = 2.50e4*Omh2/thet**4 ! z_eq; equation (2)

      ! Equation (5); changed from /zd -> /(1+zd)
      ! Now lines up with http://background.uchicago.edu/~whu/transfer/tf_fit.c (thanks Steven Murray)
      ! Difference between paper and public code
      rd = 31500.*Obh2/thet**4/(1.+zd)

      re = 31500.*Obh2/thet**4/ze ! Equation (5)
      rke = 7.46e-2*Omh2/thet**2 ! k_eq; equation (3)
      s = (2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1.+sqrt(re))) ! Sound horizon at drag; equation (6)
      rks = 1.6*((Obh2)**0.52)*((Omh2)**0.73)*(1.+(10.4*Omh2)**(-0.95)) ! Silk k; equation(7)

      q = rk/13.41/rke ! Equation (10)

      y = (1.+ze)/(1.+zd) ! y that enters in equation G(y) in equations (14) and (15)
      g = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.))) ! Equation (15)
      ab = g*2.07*rke*s/(1.+rd)**(0.75) ! Equation (14)

      a1 = (46.9*Omh2)**0.670*(1+(32.1*Omh2)**(-0.532)) ! Equation (11)
      a2 = (12.0*Omh2)**0.424*(1+(45.0*Omh2)**(-0.582)) ! Equation (11)
      ac = (a1**(-Om_b/Om_m))*(a2**(-(Om_b/Om_m)**3)) ! Equation (11)

      b1 = 0.944/(1.+(458.*Omh2)**(-0.708)) ! Equation (12)
      b2 = (0.395*Omh2)**(-0.0266) ! Equation (12)
      bc = 1./(1.+b1*((1.-Om_b/Om_m)**b2-1.)) ! Equation (12)

      f = 1./(1.+(rk*s/5.4)**4) ! Equation (18)

      c1 = 14.2+386./(1.+69.9*q**1.08) ! Equation (20) without alpha_c in the 14.2/alpha_c first bit
      c2 = 14.2/ac+386./(1.+69.9*q**1.08) ! Equation (20) (C function should have explicity alpha_c dependence in paper)
      tc = f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q**2)+(1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q**2) ! Equation (17)

      bb = 0.5+(Om_b/Om_m)+(3.-2.*Om_b/Om_m)*sqrt((17.2*Omh2)**2+1.) ! Equation (24)
      bn = 8.41*(Omh2)**0.435 ! Equation (23)
      ss = s/(1.+(bn/rk/s)**3)**(1./3.) ! Equation (22)
      tb = log(e+1.8*q)/(log(e+1.8*q)+c1*q**2)/(1.+(rk*s/5.2)**2) ! First term in equation (21)

      ! Removed this IF statement as it produced a discontinuity in P_lin(k) as cosmology
      ! was varied - thanks David Copeland for pointing this out
      !IF((rk/rks**1.4)>7.) THEN
      !   fac=0.
      !ELSE
      fac = exp(-(rk/rks)**1.4) ! Silk-damping factor from equation (21)
      !END IF

      tb = (tb+ab*fac/(1.+(bb/rk/s)**3))*sinc(rk*ss) ! Equation (21)

      tk_eh = real((Om_b/Om_m)*tb+(1.-Om_b/Om_m)*tc) ! The weighted mean of baryon and CDM transfer functions

   END FUNCTION Tk_EH

   REAL FUNCTION Tk_nw(k, cosm)

      ! No-wiggle transfer function from astro-ph:9709112
      ! TODO: Could have an init for this as many things only need to be calculated once
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: q, L, C, Gamma, wm, wb, s, h, alpha, rb
      REAL, PARAMETER :: e = exp(1.)

      ! Useful parameters to make equations shorter
      IF (cosm%power_Omegas) THEN
         h = cosm%h_pow               ! Hubble factor
         wm = cosm%Om_m_pow*cosm%h**2 ! Real matter density
         wb = cosm%Om_b_pow*cosm%h**2 ! Real baryon density      
      ELSE
         h = cosm%h               ! Hubble factor
         wm = cosm%Om_m*cosm%h**2 ! Real matter density
         wb = cosm%Om_b*cosm%h**2 ! Real baryon density       
      END IF
      rb = wb/wm ! Baryon ratio

      ! These only needs to be calculated once
      s = 44.5*log(9.83/wm)/sqrt(1.+10.*wb**0.75)              ! Equation (26)
      alpha = 1.-0.328*log(431.*wm)*rb+0.38*log(22.3*wm)*rb**2 ! Equation (31)

      ! Functions of k
      Gamma = (wm/h)*(alpha+(1.-alpha)/(1.+(0.43*k*s*h)**4)) ! Equation (30)
      q = k*(cosm%T_CMB/2.7)**2/Gamma ! Equation (28)
      L = log(2.*e+1.8*q)             ! Equation (29)
      C = 14.2+731./(1.+62.5*q)       ! Equation (29)
      Tk_nw = L/(L+C*q**2)            ! Equation (29)

   END FUNCTION Tk_nw

   REAL FUNCTION Tk_WDM(k, cosm)

      ! Warm dark matter 'correction' to the standard transfer function
      ! This version and equation references were taken from arxiv:1605.05973
      ! Originally from Bode et al. (2001; arixv:0010389)
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: alpha, mu

      alpha = 0.074*0.7*cosm%m_wdm**(-1.15) ! alpha from equation (5), units Mpc/h
      mu = 1.12                             ! mu from equation (4), dimensionless

      Tk_wdm = (1.+(alpha*k)**(2.*mu))**(-5./mu) ! Equation (2)

   END FUNCTION Tk_WDM

   REAL FUNCTION Tk_bump(k, cosm)

      ! Put a Gaussian bump in a linear power spectrum
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm

      Tk_bump = 1.+cosm%A_bump*exp(-(log(k/cosm%k_bump)**2/(2.*cosm%sigma_bump**2)))

   END FUNCTION Tk_bump

   REAL FUNCTION Tk_bump_Mexico(k, cosm)

      ! Put a Gaussian bump in a linear power spectrum, no factor of 2 in the exponential
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm

      Tk_bump_Mexico = sqrt(1.+cosm%A_bump*exp(-(log(k/cosm%k_bump)**2/cosm%sigma_bump**2)))

   END FUNCTION Tk_bump_Mexico

   ! REAL FUNCTION Tnu_approx(cosm)
   !
   !    ! How the matter power spectrum would be changed if some fraction of the mass is converted to massive neutrinos
   !    ! Approximation for how power is suppressed by massive nu at small scales
   !    ! Calculated assuming perturbation grow from z~1000 and that neutrinos are hot and therefore completely smooth
   !    ! Related to the growth-function approximation: g(a) = a^(1-3f_nu/5)
   !    ! TODO: This is NEVER used and is potentially VERY confusing
   !    ! TODO: This is NOT the transfer function relating the matter spectrum to the cold spectrum
   !    ! TODO: This IS the transfer function relation matter power spectra in different models
   !    ! TODO: Take EXTREME caution here
   !    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
   !
   !    Tnu_approx = sqrt(1.-8.*cosm%f_nu)
   !
   ! END FUNCTION Tnu_approx

   REAL FUNCTION Tk_cold_EH(k, a, cosm)

      ! Calculates the ratio of T(k) for cold vs. all matter
      ! Cold perturbation defined such that 1+delta = rho_cold/rho_matter
      ! Uses approximations from Eisenstein & Hu (1999; astro-ph/9710252)
      ! Note that this assumes that the neutrino mass is split evenly between the number of massive species
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc 
      REAL, INTENT(IN) :: a ! Scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL :: D, Dcb, Dcbnu, pcb, zeq, q, yfs, z
      REAL :: BigT
      LOGICAL, PARAMETER :: EdS_growth = Tk_cold_EdS_growth

      IF (cosm%m_nu == 0.) THEN

         ! Fix to unity if there are no neutrinos
         Tk_cold_EH = 1.

      ELSE

         ! Get the redshift
         z = redshift_a(a)

         ! Growth exponent under the assumption that neutrinos are completely unclustered (equation 11)
         pcb = (5.-sqrt(1.+24.*(1.-cosm%f_nu)))/4.

         ! Theta for temperature (BigT=T/2.7 K)
         BigT = cosm%T_CMB/2.7

         ! The matter-radiation equality redshift
         zeq = (2.5e4)*cosm%Om_m*(cosm%h**2)*BigT**(-4)

         ! The growth function normalised such that D=(1.+z_eq)/(1+z) at early times (when Omega_m ~ 1)
         ! Note that I previously used the 'wrong' growth rate below in that I took the EdS result
         ! NOTE: The two different options below gives per-cent level differences for CAMB vs HMx comparison
         IF (EdS_growth) THEN
            D = a ! Einstein-de Sitter growth function
         ELSE
            D = ungrow(a, cosm) ! Standard growth function normalised g(a<<1) = a
         END IF
         D = (1.+zeq)*D

         ! Wave number relative to the horizon scale at equality (equation 5)
         ! Extra factor of h becauase all my k are in units of h/Mpc
         q = k*cosm%h*BigT**2/(cosm%om_m*cosm%h**2)

         ! Free streaming scale (equation 14)
         ! Note that Eisenstein & Hu (1999) only consider the case of 3 neutrinos
         ! with Nnu of these being massive with the mass split evenly between Nnu species.
         yfs = 17.2*cosm%f_nu*(1.+0.488*cosm%f_nu**(-7./6.))*(cosm%N_nu*q/cosm%f_nu)**2

         ! These are (almost) the scale-dependent growth functions for each component in Eisenstein & Hu (1999)
         ! Some part is missing, but this cancels when they are divided by each other, which is all I need them for.
         ! Equations (12) and (13)
         Dcb = (1.+(D/(1.+yfs))**0.7)**(pcb/0.7)
         Dcbnu = ((1.-cosm%f_nu)**(0.7/pcb)+(D/(1.+yfs))**0.7)**(pcb/0.7)

         ! Finally, the ratio
         Tk_cold_EH = Dcb*(1.-cosm%f_nu)/Dcbnu

      END IF

   END FUNCTION Tk_cold_EH

   REAL FUNCTION primordial_spectrum(k, cosm)

      ! The unnormalised, but dimensionless, primordial spectrum
      ! TODO: Add running of running, ...
      REAL, INTENT(IN) :: k
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: pow, lnkpiv

      lnkpiv = log(k/cosm%kpiv)
      pow = cosm%ns+3.+cosm%nrun*lnkpiv/2.!+cosm%nrunrun*lnkpiv**2/6.
      primordial_spectrum = (k/cosm%kpiv)**pow

   END FUNCTION primordial_spectrum

   REAL RECURSIVE FUNCTION plin(k, a, flag, cosm)

      ! Linear matter power spectrum
      ! Must be recursive function because normalise_power calls this function again
      ! TODO: Should not access interpolator internals
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: kmax, pmax
      REAL, PARAMETER :: kmin_zero = kmin_abs_plin ! Below this wavenumber the power is fixed to zero
      REAL, PARAMETER :: kmax_zero = kmax_abs_plin ! Above this wavenumber the power is fixed to zero

      ! This line generates recursion
      IF (.NOT. cosm%is_normalised) CALL normalise_power(cosm)

      IF (k <= kmin_zero) THEN
         ! If plin happens to be foolishly called for very low k
         ! This call should never happen, but may in integrals
         plin = 0.
      ELSE IF (k > kmax_zero) THEN
         ! Avoids some issues if plin is called for absurdly high k values
         ! For some reason crashes can occur if this is the case
         plin = 0.
      ELSE IF (cosm%box .AND. k < cosm%kbox) THEN
         ! If investigating effects caused by a finite box size
         plin = 0.
      ELSE
         IF (cosm%has_power) THEN
            kmax = cosm%plin%xmax ! TODO: Should not be accessing array internals like this
            IF (cosm%scale_dependent_growth) THEN
               IF (plin_extrap .AND. k > kmax) THEN
                  pmax = evaluate_interpolator(kmax, a, cosm%plin)
                  plin = plin_extrapolation(k, kmax, pmax, cosm%ns)
               ELSE
                  plin = evaluate_interpolator(k, a, cosm%plin)
               END IF
            ELSE
               IF (plin_extrap .AND. k > kmax) THEN
                  pmax = evaluate_interpolator(kmax, 1., cosm%plin)
                  plin = plin_extrapolation(k, kmax, pmax, cosm%ns)
               ELSE
                  plin = evaluate_interpolator(k, 1., cosm%plin)
               END IF
               plin = plin*grow(a, cosm)**2
            END IF
         ELSE
            ! In this case get the power from the transfer function
            plin = primordial_spectrum(k, cosm)*(Tk_matter(k, a, cosm)*grow(a, cosm))**2
         END IF
      END IF
      plin = plin*cosm%A**2

      IF (flag == flag_cold .OR. flag == flag_ucold) THEN
         plin = plin*Tk_cold(k, a, cosm)**2
         IF (flag == flag_ucold) plin = plin/(1.-cosm%f_nu)**2
      END IF

   END FUNCTION plin

   REAL FUNCTION plin_extrapolation(k,kmax,pmax,ns)

      ! Extrapolate linear power at small scales assuming Delta^2(k) goes like ln(k)^2 k^(n-1)
      ! This works really badly if kmax is not high enough; maybe best just to use power-law
      ! TODO: It is really weird that log(k) appears, rather than log(k/kmax), this must be wrong!
      ! TODO: Check ln(k)^2 k^(3+n) at small scales (massive neutrinos; also k has dimensions?)
      REAL, INTENT(IN) :: k    ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: kmax ! Maximum wavenumber [h/Mpc]
      REAL, INTENT(IN) :: pmax ! Delta^2(k) at maximum wavenumber
      REAL, INTENT(IN) :: ns   ! Specral index

      plin_extrapolation = pmax*((log(k)/log(kmax))**2)*(k/kmax)**(ns-1.)

   END FUNCTION plin_extrapolation

   SUBROUTINE calculate_Plin(k, a, Pk, flag, cosm)

      ! Fill array P(k, a) from input arrays of k and a
      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology) :: cosm
      INTEGER :: ik, ia
      INTEGER :: nk, na

      nk = size(k)
      na = size(a)
      ALLOCATE(Pk(nk, na))
      DO ia = 1, na
         DO ik = 1, nk
            Pk(ik, ia) = plin(k(ik), a(ia), flag, cosm)
         END DO
      END DO

   END SUBROUTINE calculate_Plin

   SUBROUTINE init_sigma(cosm)

      ! This fills up interpolator of r vs. sigma(r) across a range in r
      ! It is used only in look-up for further calculations of sigma(r) and not otherwise
      ! This prevents a large number of calls to the sigma integration functions in HMx
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: R(:), a(:), sig(:, :)
      INTEGER :: ir, ia, na
      REAL :: amin, amax
      REAL, PARAMETER :: rmin = rmin_sigma
      REAL, PARAMETER :: rmax = rmax_sigma
      INTEGER, PARAMETER :: nr = nr_sigma

      ! Write to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_SIGMA: Filling sigma(R) interpolation tables'
         WRITE (*, *) 'INIT_SIGMA: R minimum [Mpc/h]:', real(rmin)
         WRITE (*, *) 'INIT_SIGMA: R maximum [Mpc/h]:', real(rmax)
         WRITE (*, *) 'INIT_SIGMA: Number of R points:', nr
      END IF

      ! Allocate and fill array of R values
      CALL fill_array_log(rmin, rmax, R, nr)

      ! This does not need to be evaulated at multiple a unless growth is scale dependent
      IF (cosm%scale_dependent_growth) THEN

         IF (cosm%has_power) THEN
            ! TODO: Should probably not be accessing interpolator internals like this
            amin = cosm%plin%ymin
            amax = cosm%plin%ymax
         ELSE
            amin = amin_sigma
            amax = amax_sigma
         END IF

         na = na_sigma
         CALL fill_array_log(amin, amax, a, na)
         IF (cosm%verbose) THEN
            WRITE (*, *) 'INIT_SIGMA: amin:', real(amin)
            WRITE (*, *) 'INIT_SIGMA: amax:', real(amax)
            WRITE (*, *) 'INIT_SIGMA: Number of a points:', na
         END IF

      ELSE
         na = 1
         ALLOCATE(a(na))
         a = 1.
      END IF

      ALLOCATE(sig(nr, na))     

      ! Loop over R and a and calculate sigma(R,a)
      IF (cosm%verbose) WRITE (*, *) 'INIT_SIGMA: Calculating sigma(R)'
      DO ia = 1, na
         DO ir = 1, nr
            sig(ir, ia) = sigma_integral(R(ir), a(ia), sigma_store, cosm)
         END DO
      END DO

      ! Initialise interpolator
      CALL init_interpolator(R, a, sig, cosm%sigma, &
         iorder_interp_sigma, &
         iextrap_sigma, &
         store = store_sigma, &
         logx = .TRUE., &
         logy = .TRUE., &
         logf = .TRUE. &
         )

      ! Change flag so that it is known that the look-up tables are filled
      cosm%has_sigma = .TRUE.

      ! Write useful information to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_SIGMA: Minimum sigma (a=1):', sig(nr, na)
         WRITE (*, *) 'INIT_SIGMA: Maximum sigma (a=1):', sig(1, na)
         WRITE (*, *) 'INIT_SIGMA: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_sigma

   REAL FUNCTION find_sigma(R, a, cosm)

      ! Evaluate interpolator for sigma(R)
      REAL, INTENT(IN) :: R ! Smoothing scale to calculate sigma [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
   
      IF (cosm%scale_dependent_growth) THEN
         find_sigma = evaluate_interpolator(R, a, cosm%sigma)      
      ELSE
         find_sigma = grow(a, cosm)*evaluate_interpolator(R, 1., cosm%sigma)
      END IF

   END FUNCTION find_sigma

   RECURSIVE REAL FUNCTION sigma(R, a, flag, cosm)

      ! Either calculates sigma(R) or evaluates interpolator
      REAL, INTENT(IN) :: R ! Smoothing scale to calculate sigma [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

      IF (.NOT. cosm%has_sigma) CALL init_sigma(cosm)
      IF(flag == sigma_store) THEN
         sigma = find_sigma(R, a, cosm)
      ELSE
         sigma = sigma_integral(R, a, flag, cosm)
      END IF

   END FUNCTION sigma

   RECURSIVE REAL FUNCTION sigma_integral(r, a, flag, cosm)

      ! Calculates sigma(R) by intergration
      REAL, INTENT(IN) :: r ! Smoothing scale to calculate sigma [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: acc = acc_sigma ! TODO: Because this calculates sigma^2(R), but we only care about sigma(R)?
      INTEGER, PARAMETER :: iorder = iorder_sigma
      REAL, PARAMETER :: tmin = 0. ! Minimum value for integration (t=0 => k->infinity)
      REAL, PARAMETER :: tmax = 1. ! Minimum value for integration (t=1 => k->0)

      sigma_integral = integrate_cosm(tmin, tmax, sigma2_integrand, r, a, flag, cosm, acc, iorder)
      sigma_integral = sqrt(sigma_integral)

   END FUNCTION sigma_integral

   RECURSIVE REAL FUNCTION sigma2_integrand(t, R, a, flag, cosm)

      ! The integrand for the sigma(R) integrals
      USE special_functions
      REAL, INTENT(IN) :: t
      REAL, INTENT(IN) :: R
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: k, kR, w_hat
      REAL, PARAMETER :: alpha = alpha_sigma

      ! Integrand to the sigma integral in terms of t. Defined by kR=(1/t-1)**alpha
      ! alpha can be any positive number, can even be a function of R
      IF (t <= 0. .OR. t >= 1.) THEN
         ! t=0 corresponds to k=infintiy when W(kR)=0
         ! t=1 corresponds to k=0 when P(k)=0
         sigma2_integrand = 0.
      ELSE
         kR = (-1.+1./t)**alpha
         k = kR/R
         w_hat = wk_tophat(kR)
         sigma2_integrand = plin(k, a, flag, cosm)*(w_hat**2)*alpha/(t*(1.-t))
      END IF

   END FUNCTION sigma2_integrand

   REAL FUNCTION sigmaV(R, a, flag, cosm)

      ! Standard deviation in the displacement field [Mpc/h]
      ! Integrand changed from k [0,inf] to t[0,1] via kR = (1/t-1)**alpha
      REAL, INTENT(IN) :: R
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: tmin = 0.
      REAL, PARAMETER :: tmax = 1.
      REAL, PARAMETER :: acc = acc_sigmaV
      INTEGER, PARAMETER :: iorder = iorder_sigmaV

      sigmaV = integrate_cosm(tmin, tmax, sigmaV2_integrand, R, a, flag, cosm, acc, iorder)

      ! Convert 3D sigmaV^2 to 1D sigmaV
      sigmaV = sqrt(sigmaV/3.)

   END FUNCTION sigmaV

   REAL FUNCTION sigmaV2_integrand(t, R, a, flag, cosm)

      ! This is the integrand for the velocity dispersion integral
      USE special_functions
      REAL, INTENT(IN) :: t
      REAL, INTENT(IN) :: R
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: k, kR, w_hat
      REAL, PARAMETER :: alpha = alpha_sigmaV

      IF (t == 0. .OR. t == 1.) THEN
         ! t = 0 corresponds to k = infinity when W(kR) = 0
         ! t = 1 corresponds to k = 0 when P(k) = 0
         sigmaV2_integrand = 0.
      ELSE
         IF (R == 0.) THEN
            kR = 0.
            k = (-1.+1./t)**alpha
         ELSE
            kR = (-1.+1./t)**alpha
            k = kR/R
         END IF
         w_hat = wk_tophat(kR)
         sigmaV2_integrand = (plin(k, a, flag, cosm)/k**2)*(w_hat**2)*alpha/(t*(1.-t))
      END IF

   END FUNCTION sigmaV2_integrand

   REAL FUNCTION neff(r, a, flag, cosm)

      ! Effective spectral index on scale R
      ! 3 + neff = -dln(sigma^2)/dlnR and ranges from ~1 on large scales to ~-3 on small scales
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm

      neff = -3.-dsigma(r, a, flag, cosm)

   END FUNCTION neff

   REAL FUNCTION dsigma(r, a, flag, cosm)

      ! dln(sigma^2)/dlnR
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sig
      REAL, PARAMETER :: tmin = 0.
      REAL, PARAMETER :: tmax = 1.
      REAL, PARAMETER :: acc = acc_dsigma
      INTEGER, PARAMETER :: iorder = iorder_dsigma

      sig = sigma(r, a, flag, cosm)
      dsigma = 2.*integrate_cosm(tmin, tmax, dsigma_integrand, r, a, flag, cosm, acc, iorder)/sig**2

   END FUNCTION dsigma

   REAL FUNCTION dsigma_integrand(t, R, a, flag, cosm)

      ! Integral for calculating dln(sigma^2)/dlnR
      ! Transformation is kR = (1/t-1)**alpha
      USE special_functions
      REAL, INTENT(IN) :: t
      REAL, INTENT(IN) :: R
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: w_hat, w_hat_deriv, k, kR
      REAL, PARAMETER :: alpha = alpha_dsigma

      IF (t <= 0. .OR. t >= 1.) THEN
         dsigma_integrand = 0.
      ELSE
         kR = (-1.+1./t)**alpha
         k = kR/R
         w_hat = wk_tophat(kR)
         w_hat_deriv = wk_tophat_deriv(kR)
         dsigma_integrand = plin(k, a, flag, cosm)*w_hat*w_hat_deriv*kR*alpha/(t*(1.-t))     
      END IF

   END FUNCTION dsigma_integrand

   REAL FUNCTION ncur(r, a, flag, cosm)

      ! -d^2ln(sigma^2)/dlnR^2
      ! Spectral curvature in linear spectrum on scale 'r'
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm

      ncur = -ddsigma(r, a, flag, cosm)

   END FUNCTION ncur

   REAL FUNCTION ddsigma(r, a, flag, cosm)

      ! Integral for calculating dln(sigma^2)/dlnR
      ! Transformation is kR = (1/t-1)**alpha
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sig, dsig
      REAL, PARAMETER :: tmin = 0.
      REAL, PARAMETER :: tmax = 1.
      REAL, PARAMETER :: acc = acc_ddsigma
      INTEGER, PARAMETER :: iorder = iorder_ddsigma

      sig = sigma(r, a, flag, cosm)
      dsig = dsigma(r, a, flag, cosm)
      ddsigma = 2.*integrate_cosm(tmin, tmax, ddsigma_integrand, r, a, flag, cosm, acc, iorder)/sig**2
      ddsigma = ddsigma+dsig-dsig**2

   END FUNCTION ddsigma

   REAL FUNCTION ddsigma_integrand(t, R, a, flag, cosm)

      ! Transformation is kR = (1/t-1)**alpha
      USE special_functions
      REAL, INTENT(IN) :: t
      REAL, INTENT(IN) :: R
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: w_hat, dw_hat, ddw_hat, k, kR
      REAL, PARAMETER :: alpha = alpha_ddsigma

      IF (t <= 0. .OR. t >= 1.) THEN
         ddsigma_integrand = 0.
      ELSE
         kR = (-1.+1./t)**alpha
         k = kR/R
         w_hat = wk_tophat(kR)
         dw_hat = wk_tophat_deriv(kR)
         ddw_hat = wk_tophat_dderiv(kR)
         ddsigma_integrand = plin(k, a, flag, cosm)*(w_hat*ddw_hat+dw_hat**2)*(kR**2)*alpha/(t*(1.-t))     
      END IF

   END FUNCTION ddsigma_integrand

   ELEMENTAL REAL FUNCTION Lagrangian_mass(R, cosm)

      ! Lagrangian mass associated with comoving scale R
      REAL, INTENT(IN) :: R
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: rho

      rho = comoving_matter_density(cosm)
      Lagrangian_mass = 4.*pi*R**3*rho/3.

   END FUNCTION Lagrangian_mass

   ELEMENTAL REAL FUNCTION Lagrangian_radius(M, cosm)

      ! Comoving Lagranigan radius associated with mass M
      USE special_functions
      REAL, INTENT(IN) :: M
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: rho

      rho = comoving_matter_density(cosm)
      Lagrangian_radius = cbrt(M/(4.*pi*rho/3.))

   END FUNCTION Lagrangian_radius

   REAL RECURSIVE FUNCTION grow(a, cosm)

      ! Scale-independent growth function, normalised | g(a=1)=1
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
      IF (a == 1.) THEN
         grow = 1.
      ELSE
         grow = evaluate_interpolator(a, cosm%grow)
      END IF

   END FUNCTION grow

   REAL FUNCTION grow_Linder(a, cosm)

      ! Calculate the growth function from the Linder growth rate via integration
      ! Defined such that g(a=1) = 1
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: acc = acc_integral_grow
      INTEGER, PARAMETER :: iorder = iorder_integral_grow

      grow_Linder = exp(-integrate_cosm(a, 1., grow_Linder_integrand, cosm, acc, iorder))

   END FUNCTION grow_Linder

   REAL FUNCTION grow_Linder_integrand(a, cosm)

      ! Integrand for the approximate growth integral using Linder approximate growth rate
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      grow_Linder_integrand = growth_rate_Linder(a, cosm)/a

   END FUNCTION grow_Linder_integrand

   REAL FUNCTION grow_CPT(a, cosm)

      ! Carroll, Press & Turner (1992) approximation to growth function (good to 5%)
      ! https://ui.adsabs.harvard.edu/abs/1992ARA%26A..30..499C/abstract
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Om_mz, Om_vz, Om_m, Om_v

      ! Get all necessary Omega values
      Om_mz = Omega_m_norad(a, cosm)
      Om_m = cosm%Om_m
      Om_vz = Omega_v(a, cosm)+Omega_w(a, cosm)    
      Om_v = cosm%Om_v+cosm%Om_w

      ! Now call CPT twice, second time to normalise it
      grow_CPT = CPT(a, Om_mz, Om_vz)/CPT(1., Om_m, Om_v)

   END FUNCTION grow_CPT

   REAL FUNCTION CPT(a, Om_m, Om_v)

      ! The CPT growth function approximation from 1992
      REAL, INTENT(IN) :: a    ! Scale factor
      REAL, INTENT(IN) :: Om_m ! Matter-density parameter
      REAL, INTENT(IN) :: Om_v ! Vacuum-density parameter

      CPT = a*Om_m/((Om_m**(4./7.))-Om_v+(1.+Om_m/2.)*(1.+Om_v/70.))

   END FUNCTION CPT

   REAL RECURSIVE FUNCTION ungrow(a, cosm)

      ! Scale-independent growth function normalised such that g(a) = a at early (matter-dominated) times
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
      ungrow = cosm%gnorm*grow(a, cosm)

   END FUNCTION ungrow

   REAL FUNCTION ungrow_approx(a, cosm)

      ! Approximate scale-independent growth function
      ! Obtained by integrating growth rate of Omega_m^gamma(a) normalised | g(a->0) = a
      ! TODO: Will not be correct for mixed Omega_v, Omega_w models
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: f1, f2, f3, w

      IF (cosm%iw == iw_LCDM) THEN
         ungrow_approx = a*(1.-(2./11.)*(cosm%Om_v/cosm%Om_m)*a**3)
      ELSE
         w = w_de(a, cosm)
         f1 = (w-1.)/(w*(5.-6.*w))
         f2 = cosm%Om_w/cosm%Om_m
         f3 = a**(-3.*w)
         ungrow_approx = a*(1.+f1*f2*f3)
      END IF

   END FUNCTION ungrow_approx

   REAL FUNCTION ungrow_integral(a, cosm)

      ! Integral solution to growth equation from Heath (1977) for pressureless cosmologies
      ! Only exactly correct if dark energy is w = -1 or w = -1/3
      ! Not sure how well it does for more general dark-energy models
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: acc = acc_integral_grow
      INTEGER, PARAMETER :: iorder = iorder_integral_grow
      REAL :: f

      f = 2.5*cosm%Om_m*sqrt(Hubble2_norad(a, cosm))
      ungrow_integral = f*integrate_cosm(0., a, ungrow_integrand, cosm, acc, iorder)

   END FUNCTION ungrow_integral

   REAL FUNCTION ungrow_integrand(a, cosm)

      ! Integrand for integral solution to growth equation
      ! Split at some scale factor to avoid infinities
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: aeps = aeps_integral_grow

      IF (a < aeps) THEN
         ungrow_integrand = a**1.5*cosm%Om_m**(-1.5)
      ELSE
         ungrow_integrand = a**(-3)*Hubble2_norad(a, cosm)**(-1.5)
      END IF

   END FUNCTION ungrow_integrand

   REAL RECURSIVE FUNCTION growth_rate(a, cosm)

      ! Growth rate: dln(g) / dln(a) ~ Omega_m(a)^0.55 for LCDM
      ! Transitions from 1 at high z to zero at high z when DE stops growth
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
      growth_rate = evaluate_interpolator(a, cosm%grate)

   END FUNCTION growth_rate

   REAL FUNCTION growth_rate_Linder(a, cosm)

      ! Approximation for the growth rate from Linder astro-ph/0507263
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: gamma, Om_m

      IF ((cosm%Om_v .NE. 0.) .AND. (cosm%Om_w .NE. 0.)) STOP 'GROWTH_RATE_LINDER: Error, does not work if Omega_v and Omega_w both non zero'

      gamma = growth_index_Linder(cosm)

      IF (cold_growth) THEN
         Om_m = Omega_cold_norad(a, cosm)
      ELSE
         Om_m = Omega_m_norad(a, cosm)
      END IF

      growth_rate_Linder = Om_m**gamma

   END FUNCTION growth_rate_Linder

   REAL RECURSIVE FUNCTION growth_index(a, cosm)

      ! Calculates gamma in f = Omega_m(a)^gamma
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Om_m
      REAL, PARAMETER :: gamma_default = growth_index_default
      REAL, PARAMETER :: gamma_limit = growth_index_limit

      IF (cold_growth) THEN
         Om_m = Omega_cold_norad(a, cosm)
      ELSE
         Om_m = Omega_m_norad(a, cosm)
      END IF

      IF (abs(1.-Om_m) < gamma_limit) THEN
         growth_index = gamma_default
      ELSE
         growth_index = log(growth_rate(a, cosm))/log(Om_m)
      END IF

   END FUNCTION growth_index

   REAL FUNCTION growth_index_Linder(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: weff

      IF (cosm%iw == iw_LCDM) THEN
         growth_index_Linder = 0.55
      ELSE
         ! Evaluate the equation of state at z=1
         ! Bizarre discontinuous slope
         weff = w_de(0.5, cosm) 
         IF (weff < -1.) THEN
            growth_index_Linder = 0.55+0.02*(1.+weff)
         ELSE
            growth_index_Linder = 0.55+0.05*(1.+weff)
         END IF
      END IF

   END FUNCTION growth_index_Linder

   REAL RECURSIVE FUNCTION acc_growth(a, cosm)

      ! Accumulated growth function: int_0^a g(a)/a da
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
      acc_growth = evaluate_interpolator(a, cosm%agrow)

   END FUNCTION acc_growth

   SUBROUTINE init_growth(cosm)

      ! Fills look-up tables for scale-dependent growth: a vs. g(a), f(a) and G(a)
      ! TODO: Figure out why if I set amax=10, rather than amax=1, I start getting weird f(a) around a=0.001
      USE calculus_table
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i
      REAL, ALLOCATABLE :: a(:), growth(:), rate(:), agrow(:)
      REAL, ALLOCATABLE :: d(:), v(:)
      REAL :: dinit, vinit, f
      REAL :: g0, f0, bigG0
      REAL, PARAMETER :: k = 0. ! Scale-indepdent so large-scale limit fixed
      REAL, PARAMETER :: aini = aini_growth
      REAL, PARAMETER :: afin = afin_growth
      REAL, PARAMETER :: amin = amin_growth
      REAL, PARAMETER :: amax = amax_growth
      INTEGER, PARAMETER :: ng = na_growth
      REAL, PARAMETER :: acc_ODE = acc_ODE_growth
      INTEGER, PARAMETER :: imeth_ODE = imeth_ODE_growth
      INTEGER, PARAMETER :: iorder_interp = iorder_ODE_interpolation_growth
      INTEGER, PARAMETER :: ifind_interp = ifind_ODE_interpolation_growth
      INTEGER, PARAMETER :: imeth_interp = imeth_ODE_interpolation_growth
      INTEGER, PARAMETER :: iorder_agrow = iorder_integration_agrow

      ! Set the initial conditions to be in the cold matter growing mode
      ! Note that for massive neutrinos or EDE there is no asymptotic g(a) ~ a limit
      IF (EDE_growth_ics) THEN
         IF (cold_growth) THEN
            f = 1.-Omega_cold_norad(aini, cosm)+cosm%f_nu
         ELSE
            f = 1.-Omega_m_norad(aini, cosm)
         END IF
      ELSE
         IF (cold_growth) THEN
            f = cosm%f_nu
         ELSE
            f = 0.
         END IF
      END IF
      dinit = aini**(1.-3.*f/5.)
      vinit = (1.-3.*f/5.)*aini**(-3.*f/5.)

      ! Write some useful information to the screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_GROWTH: Solving growth equation' 
         WRITE (*, *) 'INIT_GROWTH: Minimum scale factor:', aini
         WRITE (*, *) 'INIT_GROWTH: Maximum scale factor:', afin
         WRITE (*, *) 'INIT_GROWTH: Fraction of missing mass at initial time:', f
         WRITE (*, *) 'INIT_GROWTH: Initial delta:', dinit
         WRITE (*, *) 'INIT_GROWTH: Initial delta derivative:', vinit
         WRITE (*, *) 'INIT_GROWTH: Number of points for look-up tables:', ng
      END IF

      ! Solve the growth ODE
      CALL ODE2_adaptive_cosm(d, v, k, a, cosm, aini, afin, &
         dinit, vinit, dvda, acc_ODE, ng, imeth_ODE, ilog=.TRUE.)
      IF (cosm%verbose) WRITE (*, *) 'INIT_GROWTH: ODE done'    
      rate = v*a/d ! Calculate growth rate

      ! Normalise so that g(z=0)=1 and store the normalising factor
      cosm%gnorm = find(1., a, d, ng, iorder_interp, ifind_interp, imeth_interp)
      IF (cosm%verbose) WRITE (*, *) 'INIT_GROWTH: Unnormalised growth at z=0:', real(cosm%gnorm)
      growth = d/cosm%gnorm

      CALL init_interpolator(a, growth, cosm%grow, &
         iorder = iorder_interp_grow, &
         iextrap = iextrap_grow, &
         store = store_grow, &
         logx = .TRUE., &
         logf = .TRUE. &
         )

      CALL init_interpolator(a, rate, cosm%grate, &
         iorder = iorder_interp_rate, &
         iextrap = iextrap_rate, &
         store = store_rate, &
         logx = .TRUE., &
         logf = .FALSE. &
         )

      !! Table integration to calculate G(a)=int_0^a g(a')/a' da' !!

      ! Set to zero, because I have an x=x+y thing later on
      ALLOCATE(agrow(ng))
      agrow = 0.

      ! Do the integral up to table position i, which fills the accumulated growth table
      DO i = 1, ng
         ! Do the integral using the arrays
         IF (i > 1) THEN
            agrow(i) = integrate_table(a, cosm%gnorm*growth/a, 1, i, iorder_agrow)
         END IF
         ! Add missing section; g(a=0)/0 = 1, so you just add on a rectangle of height g*a/a=g
         agrow(i) = agrow(i)+cosm%gnorm*growth(1)
      END DO

      CALL init_interpolator(a, agrow, cosm%agrow, &
         iorder = iorder_interp_agrow, &
         iextrap = iextrap_agrow, &
         store = store_agrow, &
         logx = .TRUE., &
         logf = .TRUE.  &
         )

      ! Set the flag to true so that this subroutine is only called once
      cosm%has_growth = .TRUE.

      ! Write stuff about growth parameter at a=1 to the screen
      ! Note that has_growth = .TRUE. must have been set before to avoid a recursion
      IF (cosm%verbose) THEN
         g0 = grow(1., cosm)
         f0 = growth_rate(1., cosm)
         bigG0 = acc_growth(1., cosm)
         WRITE (*, *) 'INIT_GROWTH: Normalised growth at z=0:', g0
         WRITE (*, *) 'INIT_GROWTH: Growth rate at z=0:', f0
         WRITE (*, *) 'INIT_GROWTH: Integrated growth at z=0:', bigG0
         WRITE (*, *) 'INIT_GROWTH: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_growth

   REAL FUNCTION velocity(d, v, k, a, cosm)

      ! This is the dd/da in \dot{\delta}/da = dd/da; needed for growth function solution
      ! TODO: dd/da = v could just be built into the ODE solver somehow?
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: crap

      ! Prevent compile-time warnings
      crap = d
      crap = k
      crap = cosm%A
      crap = a

      velocity = v

   END FUNCTION velocity

   REAL FUNCTION dvda(d, v, k, a, cosm)

      ! This is the dv in \ddot{\delta} = dv; needed for growth function ODEsolution
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Om_m, f1, f2

      IF (cold_growth) THEN
         Om_m = Omega_cold_norad(a, cosm)
      ELSE
         Om_m = Omega_m_norad(a, cosm)
      END IF

      f1 = 1.5*Om_m*G_lin(d, v, k, a, cosm)*d/a**2
      f2 = -(2.+AH_norad(a, cosm)/Hubble2_norad(a, cosm))*v/a
      dvda = f1+f2

   END FUNCTION dvda

   REAL FUNCTION dvdanl(d, v, k, a, cosm)

      ! Function used for ODE solver in non-linear growth calculation
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Om_m, f1, f2, f3

      IF (cold_growth) THEN
         Om_m = Omega_cold_norad(a, cosm)
      ELSE
         Om_m = Omega_m_norad(a, cosm)
      END IF

      f1 = 1.5*Om_m*G_nl(d, v, k, a, cosm)*d*(1.+d)/a**2
      f2 = -(2.+AH_norad(a, cosm)/Hubble2_norad(a, cosm))*v/a
      f3 = (4./3.)*(v**2)/(1.+d)

      dvdanl = f1+f2+f3

   END FUNCTION dvdanl

   REAL FUNCTION G_lin(d, v, k, a, cosm)

      ! Linear effective gravitational constant
      ! FEB 2021: Fixed 1+mu -> mu bug in f(R) gravity, gave double gravitational constant
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: crap

      ! Prevent compile-time warnings
      crap = d
      crap = v

      IF (cosm%img == img_none) THEN
         G_lin = 1.
      ELSE IF (cosm%img == img_nDGP .OR. cosm%img == img_nDGP_lin) THEN
         G_lin = 1.+1./(3.*beta_dgp(a, cosm))
      ELSE IF (cosm%img == img_fR .OR. cosm%img == img_fR_lin) THEN
         G_lin = mu_fR(k, a, cosm)
      ELSE
         STOP 'G_LIN: Error, img not specified correctly'
      END IF

   END FUNCTION G_lin

   REAL FUNCTION beta_DGP(a, cosm)

      ! DGP beta function
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
  
      beta_DGP = 1.+(2./3.)*sqrt(Hubble2(a, cosm))*cosm%H0rc*(2.+AH(a,cosm)/Hubble2(a, cosm))
  
   END FUNCTION beta_DGP

   REAL FUNCTION mu_fR(k, a, cosm)

      ! Linear gravity modification for f(R)
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
  
      mu_fR = 4.-1./(1.+Compton2(a, cosm)*(k/a)**2)
      mu_fR = mu_fR/3.
  
   END FUNCTION mu_fR

   REAL FUNCTION Compton2(a, cosm)

      ! Squared Compton wavelength for f(R)
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: Rbar0

      Rbar0 = Rbar(1., cosm)
      Compton2 = -3.*(cosm%nfR+1.)*(cosm%fR0/Rbar0)*(Rbar0/Rbar(a, cosm))**(cosm%nfR+2.)

    END FUNCTION Compton2
  
   REAL FUNCTION Rbar(a, cosm)

      ! Background R value for f(R)
      ! TODO: Include dark energy properly via Om_w
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: Om_v

      Om_v = cosm%Om_v+cosm%Om_w
      Rbar=3.*(cosm%Om_m*(a**(-3))+4.*Om_v)/Hdist**2

   END FUNCTION Rbar

   REAL FUNCTION fR_a(a, cosm)

      ! TODO: Include dark energy via Om_w
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: c1, c2, Om_v

      Om_v = cosm%Om_v+cosm%Om_w
      c1 = 1.+4.*Om_v/cosm%Om_m
      c2 = (a**(-3))+4.*Om_v/cosm%Om_m

      fR_a = cosm%fR0*((c1/c2)**(cosm%nfR+1.))
  
   END FUNCTION fR_a

   REAL FUNCTION G_nl(d, v, k, a, cosm)

      ! Non-linear effective gravitational constant modification
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%img == img_none) THEN
         G_nl = 1.
      ELSE IF (cosm%img == img_nDGP) THEN  
         G_nl = G_DGP(d, a, cosm)
      ELSE IF (cosm%img == img_nDGP_lin .OR. cosm%img == img_fR_lin) THEN
         G_nl = G_lin(d, v, k, a, cosm)
      ELSE IF (cosm%img == img_fR) THEN
         STOP 'G_NL: Error, non-linear calculation not supported for f(R)'
      ELSE
         STOP 'G_NL: Error, img not specified correctly'
      END IF

   END FUNCTION G_nl

   REAL FUNCTION G_DGP(d, a, cosm)

      ! Non-linear effective gravitational constant for DGP models in spherical case
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: M, r, r3
      REAL, PARAMETER :: eps = 1e-2

      M = 1. ! Mass perturbation can be anything, it cancels out in the end  
      !r = (3.*M/(4.*pi*comoving_matter_density(cosm)*d))**(1./3.) 
      r = Lagrangian_radius(d, cosm) ! Convert the mass perturbation to a comoving radius (M=4*pi*r^3*delta/3)
      r = r*a ! Convert comoving -> physical radius      
      r3 = (r/r_Vainshtein_DGP(M, a, cosm))**3 ! G_nl depends on r3 only
      IF((1./r3) < eps) THEN
         ! High r3 expansion to avoid cancellation problems
         G_DGP = 1.+(1./(3.*beta_DGP(a, cosm)))*(1.-1./(4.*r3))
      ELSE
         ! Standard formula
         G_DGP = 1.+(2./(3.*beta_DGP(a, cosm)))*r3*(sqrt(1.+1./r3)-1.)
      END IF

   END FUNCTION G_DGP

   REAL FUNCTION r_Vainshtein_DGP(M, a, cosm)

      ! nDGP Vainshtein radius in physical coordinates
      ! This Vainshtein radius is in physical coordinates
      ! M(r) is the mass *perturbation* interior to r
      REAL, INTENT(IN) :: M
      REAL, INTENT(IN) :: a 
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: GN = bigG_cos/H0_cos**2 ! G/H0^2 in units (Mpc/h)^3 (M_sun/h)^-1
      
      r_Vainshtein_DGP = (16.*GN*M*cosm%H0rc**2)/(9.*beta_DGP(a, cosm)**2)
      r_Vainshtein_DGP = r_Vainshtein_DGP**(1./3.)
  
    END FUNCTION r_Vainshtein_DGP

   REAL FUNCTION dc_NakamuraSuto(a, cosm)

      ! Nakamura & Suto (1997; arXiv:astro-ph/9612074) fitting formula for spherical-collapse in LCDM
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Om_mz
      LOGICAL :: cold = cold_Nakamura

      IF (cold) THEN
         Om_mz = Omega_cold_norad(a, cosm)
      ELSE
         Om_mz = Omega_m_norad(a, cosm)
      END IF
      dc_NakamuraSuto = dc0*(1.+0.012299*log10(Om_mz))

   END FUNCTION dc_NakamuraSuto

   REAL FUNCTION Dv_BryanNorman(a, cosm)

      ! Bryan & Norman (1998; arXiv:astro-ph/9710107) spherical over-density fitting function
      ! Here overdensity is defined relative to the background matter density, rather than the critical density
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: x, Om_mz
      LOGICAL :: cold = cold_Bryan

      IF (cold) THEN
         Om_mz = Omega_cold_norad(a, cosm)
      ELSE
         Om_mz = Omega_m_norad(a, cosm)
      END IF
      x = Om_mz-1.

      IF (cosm%Om_v == 0. .AND. cosm%Om_w == 0.) THEN
         ! Open model results
         Dv_BryanNorman = Dv0+60.*x-32.*x**2
         Dv_BryanNorman = Dv_BryanNorman/Om_mz
      ELSE
         ! LCDM results
         Dv_BryanNorman = Dv0+82.*x-39.*x**2
         Dv_BryanNorman = Dv_BryanNorman/Om_mz
      END IF

   END FUNCTION Dv_BryanNorman

   REAL FUNCTION dc_Mead(a, cosm)

      ! delta_c fitting function from Mead (2017; 1606.05345)
      REAL, INTENT(IN) :: a ! scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: lg, bG, Om_m, ai
      TYPE(cosmology) :: cosm_LCDM

      ! See Appendix A of Mead (2017) for naming convention
      REAL, PARAMETER :: p10 = -0.0069
      REAL, PARAMETER :: p11 = -0.0208
      REAL, PARAMETER :: p12 = 0.0312
      REAL, PARAMETER :: p13 = 0.0021
      INTEGER, PARAMETER :: a1 = 1
      REAL, PARAMETER :: p20 = 0.0001
      REAL, PARAMETER :: p21 = -0.0647
      REAL, PARAMETER :: p22 = -0.0417
      REAL, PARAMETER :: p23 = 0.0646
      INTEGER, PARAMETER :: a2 = 0

      IF (cosm%m_nu /= 0.) THEN
         cosm_LCDM = convert_cosmology(cosm, make_lambda=.FALSE., make_flat=.FALSE., remove_neutrinos=.TRUE.)
      ELSE
         cosm_LCDM = cosm
      END IF

      lg = ungrow(a, cosm_LCDM)
      bG = acc_growth(a, cosm_LCDM)
      Om_m = Omega_m_norad(a, cosm_LCDM)
      ai = a
 
      dc_Mead = 1.  
      dc_Mead = dc_Mead+f_Mead(lg/ai, bG/ai, p10, p11, p12, p13)*log10(Om_m)**a1
      dc_Mead = dc_Mead+f_Mead(lg/ai, bG/ai, p20, p21, p22, p23)
      dc_Mead = dc_Mead*dc0*(1.-0.041*cosm%f_nu)

   END FUNCTION dc_Mead

   REAL FUNCTION Dv_Mead(a, cosm)

      ! Delta_v fitting function from Mead (2017; 1606.05345)
      REAL, INTENT(IN) :: a !scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: lg, bG, Om_m, ai
      TYPE(cosmology) :: cosm_LCDM

      ! See Appendix A of Mead (2017) for naming convention
      REAL, PARAMETER :: p30 = -0.79
      REAL, PARAMETER :: p31 = -10.17
      REAL, PARAMETER :: p32 = 2.51
      REAL, PARAMETER :: p33 = 6.51
      INTEGER, PARAMETER :: a3 = 1
      REAL, PARAMETER :: p40 = -1.89
      REAL, PARAMETER :: p41 = 0.38
      REAL, PARAMETER :: p42 = 18.8
      REAL, PARAMETER :: p43 = -15.87
      INTEGER, PARAMETER :: a4 = 2

      IF (cosm%m_nu /= 0.) THEN
         cosm_LCDM = convert_cosmology(cosm, make_lambda=.FALSE., make_flat=.FALSE., remove_neutrinos=.TRUE.)
      ELSE
         cosm_LCDM = cosm
      END IF

      lg = ungrow(a, cosm_LCDM)
      bG = acc_growth(a, cosm_LCDM)
      Om_m = Omega_m_norad(a, cosm_LCDM)
      ai = a

      Dv_Mead = 1.    
      Dv_Mead = Dv_Mead+f_Mead(lg/ai, bG/ai, p30, p31, p32, p33)*log10(Om_m)**a3
      Dv_Mead = Dv_Mead+f_Mead(lg/ai, bG/ai, p40, p41, p42, p43)*log10(Om_m)**a4
      Dv_Mead = Dv_Mead*Dv0*(1.+0.763*cosm%f_nu)

   END FUNCTION Dv_Mead

   REAL FUNCTION f_Mead(x, y, p0, p1, p2, p3)

      ! Equation (A3) in Mead (2017)
      REAL, INTENT(IN) :: x, y
      REAL, INTENT(IN) :: p0, p1, p2, p3

      f_Mead = p0+p1*(1.-x)+p2*(1.-x)**2+p3*(1.-y)

   END FUNCTION f_Mead

   REAL FUNCTION dc_spherical(a, cosm)

      ! Get delta_c from a spherical-collapse calculation or look-up table
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interp_dc
      INTEGER, PARAMETER :: ifind = ifind_interp_dc
      INTEGER, PARAMETER :: imeth = imeth_interp_dc

      IF (cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)

      IF (a < cosm%dc%xmin) THEN
         dc_spherical = dc0
      ELSE
         dc_spherical = evaluate_interpolator(a, cosm%dc)
      END IF

   END FUNCTION dc_spherical

   REAL FUNCTION Dv_spherical(a, cosm)

      ! Get Delta_v from a spherical-collapse calculation or look-up table
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interp_Dv
      INTEGER, PARAMETER :: ifind = ifind_interp_Dv
      INTEGER, PARAMETER :: imeth = imeth_interp_Dv

      IF (cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)

      IF (a < cosm%Dv%xmin) THEN
         Dv_spherical = Dv0
      ELSE
         Dv_spherical = evaluate_interpolator(a, cosm%Dv)
      END IF

   END FUNCTION Dv_spherical

   SUBROUTINE init_spherical_collapse(cosm)

      ! Initialise the spherical-collapse calculation
      ! TODO: Care with initial conditions and neutrinos/EDE here
      USE table_integer
      USE minimization
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: dinit, ainit, vinit, ac, f
      REAL :: av, a_rmax, d_rmax, rmax, rv
      REAL, ALLOCATABLE :: d(:), a(:), v(:), dc(:), Dv(:), aa(:)
      REAL, ALLOCATABLE :: dnl(:), vnl(:), rnl(:)
      REAL, ALLOCATABLE :: a_coll(:), r_coll(:)
      INTEGER :: i, j, k, k2, nn
      REAL, PARAMETER :: amax = amax_spherical
      INTEGER, PARAMETER :: imeth_ODE = imeth_ODE_spherical
      REAL, PARAMETER :: dmin = dmin_spherical
      REAL, PARAMETER :: dmax = dmax_spherical
      INTEGER, PARAMETER :: n = n_spherical
      INTEGER, PARAMETER :: m = m_spherical

      IF (cosm%verbose) WRITE (*, *) 'SPHERICAL_COLLAPSE: Doing integration'

      IF (cosm%verbose) THEN
         WRITE (*, *) 'SPHERICAL_COLLAPSE: delta min', dmin
         WRITE (*, *) 'SPHERICAL_COLLAPSE: delta max', dmax
         WRITE (*, *) 'SPHERICAL_COLLAPSE: number of collapse points attempted', m
      END IF
      ALLOCATE(a(m), dc(m), Dv(m))
      a = 0.
      dc = 0.
      Dv = 0.

      ! BCs for integration. Note ainit=dinit means that collapse should occur around a=1 for dmin
      ! amax should be slightly greater than 1 to ensure at least a few points for a>0.9 (i.e not to miss out a=1)
      f = cosm%f_nu
      ainit = dmin**(1./(1.-3.*f/5.))
      vinit = (1.-3.*f/5.)*ainit**(-3.*f/5.) ! vinit=1 is EdS growing mode solution

      ! Now loop over all initial density fluctuations
      DO j = 1, m

         ! log range of initial delta
         dinit = exp(progression(log(dmin), log(dmax), j, m))

         ! Do both with the same a1 and a2 and using the same number of time steps
         ! This means that arrays a, and anl will be identical, which simplifies calculation
         CALL ODE2_spherical(dnl, vnl, 0., aa, cosm, ainit, amax, dinit, vinit, dvdanl, n, imeth_ODE, .TRUE.)
         DEALLOCATE (aa)
         CALL ODE2_spherical(d, v, 0., aa, cosm, ainit, amax, dinit, vinit, dvda, n, imeth_ODE, .TRUE.)

         ! If this condition is met then collapse occured some time a<amax
         IF (dnl(n) == 0.) THEN

            !! delta_c calcualtion !!

            ALLOCATE (rnl(n))

            rnl = aa*(1.+dnl)**(-1./3.)

            ! Find the collapse point (very crude)
            ! More accurate calculations seem to be worse
            ! I think this is due to the fact that delta spikes very quickly
            DO i = 1, n
               IF (dnl(i) == 0.) THEN
                  ! k is the new maxium size of the arrays
                  k = i-1
                  EXIT
               END IF
            END DO

            ! Cut away parts of the arrays for a>ac
            CALL amputate_array(aa, 1, k)
            CALL amputate_array(d, 1, k)
            CALL amputate_array(dnl, 1, k)
            CALL amputate_array(rnl, 1, k)

            ! Collapse has occured so use previous a as ac and d as dc
            ac = aa(k)
            dc(j) = d(k)

            !! !!

            !! Now to Delta_v calculation !!

            ! Find the a values when the perturbation is maximum size
            a_rmax = find_array_maximum(aa, rnl)

            ! Find the over-density at this point
            d_rmax = exp(find(log(a_rmax), log(aa), log(dnl), SIZE(aa), &
               iorder=1, ifind=ifind_split, iinterp=iinterp_Lagrange))

            ! Find the maximum radius
            rmax = find(log(a_rmax), log(aa), rnl, SIZE(aa), &
               iorder=1, ifind=ifind_split, iinterp=iinterp_Lagrange)

            ! The radius of the perturbation when it is virialised is half maximum
            ! This might not be appropriate for LCDM models (or anything with DE)
            rv = rmax/2.

            ! Need to assign new arrays for the collapse branch of r such that it is monotonic
            !k2=int_split(d_rmax,dnl,k)
            k2 = find_table_integer(d_rmax, dnl, ifind_split)

            ! Allocate collapse branch arrays
            ALLOCATE (a_coll(k-k2+1), r_coll(k-k2+1))

            ! Fill collapse branch arrays
            DO i = k2, k
               a_coll(i-k2+1) = aa(i)
               r_coll(i-k2+1) = rnl(i)
            END DO

            ! Find the scale factor when the perturbation has reached virial radius
            av = exp(find(rv, r_coll, log(a_coll), SIZE(r_coll), &
               iorder=3, ifind=ifind_split, iinterp=iinterp_Lagrange))

            ! Deallocate collapse branch arrays
            DEALLOCATE (a_coll, r_coll)

            ! Spherical model approximation is that perturbation is at virial radius when
            ! 'collapse' is considered to have occured, which has already been calculated
            Dv(j) = exp(find(log(av), log(aa), log(dnl), SIZE(aa), &
               iorder=1, ifind=ifind_split, iinterp=iinterp_Lagrange))
            Dv(j) = Dv(j)*(ac/av)**3
            Dv(j) = Dv(j)+1.

            !!

            a(j) = ac

            DEALLOCATE (rnl)

         END IF

         ! Deallocate arrays ready for next calculation
         DEALLOCATE (d, v, aa)
         DEALLOCATE (dnl, vnl)

      END DO

      IF (cosm%verbose) WRITE (*, *) 'SPHERICAL COLLAPSE: calculation complete'

      ! Reverse the arrays so that they run lowest a to highest a
      CALL reverse_array(a)
      CALL reverse_array(dc)
      CALL reverse_array(Dv)

      IF (cosm%verbose) THEN
         WRITE (*, *) '===================================='
         WRITE (*, *) 'Point  scalefactor  delta_c  Delta_v'
         WRITE (*, *) '===================================='
         DO i = 1, m
            IF (a(i) == 0.) EXIT
            WRITE (*, fmt='(I5,F13.4,F9.4,F9.1)') i, a(i), dc(i), Dv(i)
         END DO
         WRITE (*, *) '===================================='
      END IF

      ! Calculate the maximum sizes for these new arrays
      DO i = 1, m
         IF (a(i) == 0.) EXIT
      END DO
      nn = i-1

      IF (cosm%verbose) THEN
         WRITE (*, *) 'SPHERICAL_COLLAPSE: number of collapse points:', nn
         WRITE (*, *)
      END IF

      ! Remove bits of the array that are unnecessary
      CALL amputate_array(a, 1, nn)
      CALL amputate_array(dc, 1, nn)
      CALL amputate_array(Dv, 1, nn)

      CALL init_interpolator(a, dc, cosm%dc, &
         iorder = iorder_interp_dc, &
         iextrap = iextrap_dc, &
         store = store_dc, &
         logx = .TRUE., &
         logf = .FALSE. &
         )

      CALL init_interpolator(a, Dv, cosm%Dv, &
         iorder = iorder_interp_Dv, &
         iextrap = iextrap_Dv, &
         store = store_Dv, &
         logx = .TRUE., &
         logf = .FALSE. &
         )

      ! Set the flag
      cosm%has_spherical = .TRUE.

   END SUBROUTINE init_spherical_collapse

   SUBROUTINE get_CAMB_power(a, na, k_Pk, Pk, nk_Pk, k_Tk, Tk_m, Tk_c, nk_Tk, non_linear, halofit_version, cosm)

      ! Runs CAMB to get a power spectrum
      ! TODO: Could this be moved to CAMB stuff? Not easily, because it requires cosmology type
      ! TODO: New CAMB cosmology class?
      ! TODO: Could split this up into run_CAMB etc. etc.
      USE CAMB_stuff
      INTEGER, INTENT(IN) :: na
      REAL, INTENT(IN) :: a(na)
      REAL, ALLOCATABLE, INTENT(OUT) :: k_Pk(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      INTEGER, INTENT(OUT) :: nk_Pk
      REAL, ALLOCATABLE, INTENT(OUT) :: k_Tk(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Tk_m(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: Tk_c(:, :)
      INTEGER, INTENT(OUT) :: nk_Tk
      LOGICAL, INTENT(IN) :: non_linear
      INTEGER, INTENT(IN) :: halofit_version
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i, j
      INTEGER :: nTk
      REAL, ALLOCATABLE :: Pk_CAMB(:), Tk_CAMB(:,:)
      REAL :: Om_c, Om_b, Om_nu, h, ombh2, omch2, omnuh2
      CHARACTER(len=256) :: infile
      CHARACTER(len=256) :: camb, dir, root, matterpower, transfer, params
      REAL, PARAMETER :: kmax = kmax_plin ! Maximum wavenumber to get the power to [h/Mpc]
      REAL, PARAMETER :: nmax = nmax_CAMB ! Multiplicative factor to go beyond kmax
      CHARACTER(len=256) :: de = de_CAMB

      ! Sort executable, directories and files
      camb = cosm%CAMB_exe
      dir = cosm%CAMB_temp_dir
      root = trim(dir)//'temp'
      matterpower = trim(root)//'_matterpower_'
      transfer = trim(root)//'_transfer_'
      params = trim(root)//'_create_params.ini'

      IF(cosm%iTk .NE. iTk_CAMB) STOP 'GET_CAMB_POWER: Normalisation will not work unless iTk_CAMB is set'

      IF (cosm%verbose) THEN
         WRITE(*,*) 'GET_CAMB_POWER: Running CAMB'
         WRITE(*,*) 'GET_CAMB_POWER: kpiv [h/Mpc]:', real(cosm%kpiv)
         WRITE(*,*) 'GET_CAMB_POWER: kpiv [1/Mpc]:', real(cosm%kpiv*cosm%h)
         WRITE(*,*) 'GET_CAMB_POWER: As:', real(cosm%As)
         WRITE(*,*) 'GET_CAMB_POWER: Minimum a:', a(1)
         WRITE(*,*) 'GET_CAMB_POWER: Maximum a:', a(na)
         WRITE(*,*) 'GET_CAMB_POWER: Number of a:', na
      END IF

      ! Set the cosmological parameters for CAMB
      IF (cosm%power_Omegas) THEN
         Om_c = cosm%Om_c_pow
         Om_b = cosm%Om_b_pow
         Om_nu = 0.
         h = cosm%h_pow
      ELSE
         Om_c = cosm%Om_c
         Om_b = cosm%Om_b
         IF (cosm%m_nu == 0.) THEN
            Om_nu = 0.
         ELSE
            Om_nu = cosm%Om_nu
         END IF
         h = cosm%h
      END IF

      IF (.NOT. is_in_array(cosm%iw, [iw_LCDM, iw_wCDM, iw_waCDM])) STOP 'GET_CAMB_POWER: Can only cope with LCDM, wCDM or w(a)CDM'

      ! Physical density parameters that CAMB requires
      ombh2 = Om_b*h**2
      omch2 = Om_c*h**2
      omnuh2 = Om_nu*h**2

      ! Remove previous parameters file
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(params))

      ! Open new parameters file for writing
      OPEN (7, file=params, status='replace')

      ! Output root
      WRITE (7, *) 'output_root = ', trim(root)

      ! Things to get
      WRITE (7, *) 'get_scalar_cls = F'
      WRITE (7, *) 'get_vector_cls = F'
      WRITE (7, *) 'get_tensor_cls = F'
      WRITE (7, *) 'get_transfer = T'

      ! Lensing
      WRITE (7, *) 'do_lensing = F'

      IF (non_linear) THEN
         WRITE (7, *) 'do_nonlinear = 1'
      ELSE
         WRITE (7, *) 'do_nonlinear = 0'
      END IF
      WRITE (7, *) 'halofit_version =', halofit_version

      ! Standard cosmological parameters
      WRITE (7, *) 'ombh2 =', ombh2
      WRITE (7, *) 'omch2 =', omch2
      WRITE (7, *) 'omnuh2 =', omnuh2
      WRITE (7, *) 'omk = ', cosm%Om_k
      WRITE (7, *) 'hubble =', 100.*h

      ! Dark energy
      WRITE (7, *) 'dark_energy_model = ', trim(de)
      WRITE (7, *) 'w =', cosm%w
      WRITE (7, *) 'wa =', cosm%wa
      WRITE (7, *) 'cs2_lam = 1'

      ! CMB temperature and Helium fraction
      WRITE (7, *) 'temp_cmb = ', cosm%T_CMB
      WRITE (7, *) 'helium_fraction = ', cosm%YHe

      ! Neutrinos
      WRITE (7, *) 'share_delta_neff = T'
      IF (omnuh2 == 0.) THEN
         ! Massless neutrinos
         WRITE (7, *) 'massless_neutrinos = ', cosm%neff
         WRITE (7, *) 'nu_mass_eigenstates = 0'
         WRITE (7, *) 'massive_neutrinos = 0'
      ELSE
         WRITE (7, *) 'massless_neutrinos = ', cosm%neff-cosm%N_nu
         WRITE (7, *) 'massive_neutrinos = ', cosm%N_nu
         WRITE (7, *) 'nu_mass_eigenstates = 1'         
         WRITE (7, *) 'nu_mass_fractions = 1'
      END IF

      ! Primordial power spectrum properties
      WRITE (7, *) 'initial_power_num = 1'
      WRITE (7, *) 'scalar_spectral_index(1) =', cosm%ns
      WRITE (7, *) 'scalar_nrun(1) =', cosm%nrun
      WRITE (7, *) 'scalar_amp(1) =', cosm%As
      WRITE (7, *) 'pivot_scalar =', cosm%kpiv*cosm%h ! Note that CAMB uses 1/Mpc whereas I use h/Mpc
      WRITE (7, *) 'pivot_tensor =', cosm%kpiv*cosm%h ! Note that CAMB uses 1/Mpc whereas I use h/Mpc

      ! Reionisation
      WRITE (7, *) 'reionization = F'
      WRITE (7, *) 're_use_optical_depth = F'
      WRITE (7, *) 're_optical_depth = 0.09'
      WRITE (7, *) 're_delta_redshift = 1.5'
      WRITE (7, *) 're_ionization_frac = -1'

      ! RECFAST
      WRITE (7, *) 'recombination_model = Recfast'
      WRITE (7, *) 'RECFAST_fudge = 1.14'
      WRITE (7, *) 'RECFAST_fudge_He = 0.86'
      WRITE (7, *) 'RECFAST_Heswitch = 6'
      WRITE (7, *) 'RECFAST_Hswitch = T'

      ! Adiabatic/isocurvature etc.
      WRITE (7, *) 'initial_condition = 1' ! 1 - Adiabatic initial condtions

      ! ?
      WRITE (7, *) 'COBE_normalize = F'
      WRITE (7, *) 'CMB_outputscale = 7.4311e12'

      ! Power spectrum requirements
      WRITE (7, *) 'transfer_high_precision = T'
      WRITE (7, *) 'transfer_kmax = ', nmax*kmax
      WRITE (7, *) 'transfer_k_per_logint = 0'
      WRITE (7, *) 'transfer_interp_matterpower = T'
      WRITE (7, *) 'transfer_power_var = 7'
      WRITE (7, *) 'transfer_num_redshifts =', na
      DO i = 1, na
         WRITE (7, *) trim(number_file(trim('transfer_redshift('), i, trim(') ='))), redshift_a(a(i))
      END DO

      ! Bispectra
      WRITE (7, *) 'do_lensing_bispectrum = F'
      WRITE (7, *) 'do_primordial_bispectrum = F'

      ! Verbosity (0 writes nothing; 1 writes some useful stuff; 2 writes HMcode stuff)
      IF (cosm%verbose) THEN
         WRITE (7, *) 'feedback_level = 1'
      ELSE
         WRITE (7, *) 'feedback_level = 0'
      END IF

      ! Print out a comment describing file headers
      WRITE (7, *) 'output_file_headers = T'

      ! Write out derived parameters
      WRITE (7, *) 'derived_parameters = T'

      ! ?
      WRITE (7, *) 'massive_nu_approx = 1'

      ! Accuracy parameters for polarization and reionization
      WRITE (7, *) 'accurate_polarization = T'
      WRITE (7, *) 'accurate_reionization = T'

      WRITE (7, *) 'lensing_method = 1'
      WRITE (7, *) 'accurate_BB = F'

      ! ?
      WRITE (7, *) 'do_tensor_neutrinos = T'

      ! ?
      WRITE (7, *) 'do_late_rad_truncation = T'

      ! OMP threads
      WRITE (7, *) 'number_of_threads = 4'

      ! Accuracy parameters
      WRITE (7, *) 'accuracy_boost = 1'
      WRITE (7, *) 'l_accuracy_boost = 1'
      WRITE (7, *) 'l_sample_boost = 1'

      ! HMcode baryonic feedback
      WRITE(7, *) 'HMcode_logT_AGN  =', log10(cosm%Theat)

      ! Close file
      CLOSE (7)

      ! Remove old files and run CAMB
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(transfer)//'*')
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(matterpower)//'*')
      IF (cosm%verbose) THEN
         CALL EXECUTE_COMMAND_LINE(trim(camb)//' '//trim(params))
      ELSE
         CALL EXECUTE_COMMAND_LINE(trim(camb)//' '//trim(params)//' > /dev/null')
      END IF
      IF (cosm%verbose) WRITE (*, *) 'GET_CAMB_POWER: CAMB run complete'

      ! Loop over redshifts and read CAMB power into arrays
      DO j = 1, na
         infile = number_file(matterpower, j, trim('.dat'))
         CALL read_CAMB_Pk(k_Pk, Pk_CAMB, nk_Pk, infile)
         IF (j == 1) THEN
            ALLOCATE(Pk(nk_Pk, na))
         END IF
         Pk(:, j) = Pk_CAMB
      END DO
      DEALLOCATE(Pk_CAMB)

      ! Do pruning
      IF (cosm%verbose) WRITE (*, *) 'GET_CAMB_POWER: nk before pruning:', nk_Pk
      CALL prune_CAMB(k_Pk, a, Pk, nk_Pk, na)
      IF (cosm%verbose) THEN
         WRITE (*, *) 'GET_CAMB_POWER: nk after pruning:', nk_Pk
         WRITE (*, *) 'GET_CAMB_POWER: Getting transfer functions'
      END IF

      ! Loop over redshifts and read CAMB transfer functions into arrays
      DO j = 1, na

         ! Read in the raw CAMB data 
         infile = number_file(transfer, j, trim('.dat'))
         CALL read_CAMB_Tk(k_Tk, Tk_CAMB, nk_Tk, nTk, infile)
         IF (j == 1) THEN
            ALLOCATE(Tk_m(nk_Tk, na), Tk_c(nk_Tk, na))
         END IF

         ! Normalise according to my definitions
         Tk_CAMB(CAMB_column_Tk_CDM, :) = cosm%Om_c*Tk_CAMB(CAMB_column_Tk_CDM, :)
         Tk_CAMB(CAMB_column_Tk_baryon, :) = cosm%Om_b*Tk_CAMB(CAMB_column_Tk_baryon, :)
         Tk_CAMB(CAMB_column_Tk_total, :) = cosm%Om_m*Tk_CAMB(CAMB_column_Tk_total, :)

         ! Normalise both matter and cold matter
         Tk_m(:, j) = Tk_CAMB(CAMB_column_Tk_total, :)/Tk_CAMB(CAMB_column_Tk_total, 1)
         Tk_c(:, j) = (Tk_CAMB(CAMB_column_Tk_CDM, :)+Tk_CAMB(CAMB_column_Tk_baryon, :))/Tk_CAMB(CAMB_column_Tk_total, :)

      END DO

      ! Write to screen
      IF(cosm%verbose) THEN
         WRITE (*, *) 'GET_CAMB_POWER: Transfer function kmin [h/Mpc]:', k_Tk(1)
         WRITE (*, *) 'GET_CAMB_POWER: Transfer function kmax [h/Mpc]:', k_Tk(nk_Tk)
         WRITE (*, *) 'GET_CAMB_POWER: Matter transfer function (kmin):', Tk_m(1, 1)
         WRITE (*, *) 'GET_CAMB_POWER: Cold transfer function (kmin):', Tk_c(1, 1)
         WRITE (*, *) 'GET_CAMB_POWER: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_CAMB_power

   SUBROUTINE prune_CAMB(k, a, Pk, nk, na)

      ! Remove some k values from the CAMB calculation of P_lin(k)
      REAL, ALLOCATABLE, INTENT(INOUT) :: k(:)
      INTEGER, INTENT(IN) :: na
      REAL, INTENT(IN) :: a(na)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Pk(:,:)
      INTEGER, INTENT(INOUT) :: nk
      INTEGER :: i, naa
      LOGICAL :: kkeep(nk), akeep(na)
      REAL, ALLOCATABLE :: aa(:)
      REAL, PARAMETER :: kmax = kmax_plin     ! Maximum wavenumber to get the power to
      REAL, PARAMETER :: pk_min = pk_min_CAMB ! Minimum power to consider

      ! Initially assume everything is being kept
      kkeep = .TRUE.
      akeep = .TRUE. 

      ! Mask array elements that are either too high k or have too low power
      DO i = 1, nk
         IF(k(i) > kmax .OR. Pk(i,na) < pk_min) THEN
            kkeep(i) = .FALSE.
         END IF
      END DO

      ALLOCATE(aa(na)) ! Define because a is INTENT(IN)
      aa=a             ! Define because a is INTENT(IN)
      naa=na           ! Define because na is INTENT(IN)
      CALL apply_mask(k, aa, Pk, nk, naa, kkeep, akeep) ! This will change nk
      IF((naa .NE. na) .OR. (size(aa) .NE. size(a))) STOP 'PRUNE_CAMB: Error, something went wrong'

   END SUBROUTINE prune_CAMB

   SUBROUTINE init_CAMB_linear(cosm)

      ! Initialise the CAMB linear power spectrum calculation
      USE CAMB_stuff
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: a(:), kPk(:), Pk(:, :), kTk(:), Tm(:, :), Tc(:, :), kkPk(:), Pkk(:, :)
      REAL :: k, fac
      INTEGER :: i, j, na, nPk, nTk, ik
      CHARACTER(len=256) :: camb, dir, root, matterpower, transfer, params
      LOGICAL, PARAMETER :: non_linear = .FALSE. ! Should not use non-linear when trying to get linear theory
      INTEGER, PARAMETER :: halofit_version = 0  ! Irrelevant here because we are getting linear spectrum
      REAL, PARAMETER :: kmin_rebin = kmin_plin  ! Minimum k if rebinning [h/Mpc]
      REAL, PARAMETER :: kmax_rebin = kmax_plin  ! Maximum k if rebinning [h/Mpc] 
      INTEGER, PARAMETER :: nk = nk_plin         ! Number of k if rebinning
      REAL, PARAMETER :: amin = amin_plin        ! Minimum scale factor to get from CAMB
      REAL, PARAMETER :: amax = amax_plin        ! Maximum scale factor to get from CAMB
      LOGICAL, PARAMETER :: rebin = rebin_CAMB   ! Should we rebin CAMB input P(k)?
      LOGICAL, PARAMETER :: fill_Tk_interpolators = .TRUE.

      ! Sort executable, directories and files
      camb = cosm%CAMB_exe
      dir = cosm%CAMB_temp_dir
      root = trim(dir)//'temp'
      matterpower = trim(root)//'_matterpower_'
      transfer = trim(root)//'_transfer_'
      params = trim(root)//'_create_params.ini'

      IF (cosm%verbose) THEN
         WRITE(*,*) 'INIT_CAMB_LINEAR: Getting linear power from CAMB'
      END IF

      ! For scale-dependent growth fill the scale-factor array first
      IF (cosm%scale_dependent_growth) THEN
         na = na_plin
         CALL fill_array_log(amin, amax, a, na)          
      ELSE
         ! For non-scale-dependent growth need a size-one array of a(1)=1.
         na = 1
         ALLOCATE(a(na))
         a = 1.
      END IF

      ! Get the CAMB P(k) at a series of 'a' values
      CALL get_CAMB_power(a, na, kPk, Pk, nPk, kTk, Tm, Tc, nTk, non_linear, halofit_version, cosm)

      ! Apply non-CAMB transfer functions to power (spikes, WDM etc.)
      DO ik = 1, nPk
         k = kPk(ik)
         fac = Tk_factor(k, cosm)
         Pk(ik, :) = Pk(ik, :)*fac**2
      END DO

      ! Initialise transfer function interpolation
      ! Transfer function needed for As
      ! TODO: Could save time here if total matter and cold matter identical
      ! TODO: This could be time consuming if interpolators filled but never used?
      IF (fill_Tk_interpolators .OR. cosm%iTc == iTc_CAMB) THEN
      
         ! Initialise total matter transfer functions for interpolation
         CALL init_interpolator(kTk, a, Tm, cosm%Tk_matter, &
            iorder = iorder_interp_Tk, &
            iextrap = iextrap_Tk, &
            store = store_Tk, &
            logx = .TRUE., &
            logy = .TRUE., &
            logf = .TRUE.)

         ! Initialise cold matter transfer function interpolation
         CALL init_interpolator(kTk, a, Tc, cosm%Tk_cold, &
            iorder = iorder_interp_Tk, &
            iextrap = iextrap_Tk, &
            store = store_Tk, &
            logx = .TRUE., &
            logy = .TRUE., &
            logf = .TRUE.)

      END IF

      DEALLOCATE(kTk, Tm, Tc)    

      IF (cosm%verbose) THEN
         WRITE(*,*) 'INIT_CAMB_LINEAR: kmin [h/Mpc]:', kPk(1)
         WRITE(*,*) 'INIT_CAMB_LINEAR: kmax [h/Mpc]:', kPk(nPk)
         WRITE(*,*) 'INIT_CAMB_LINEAR: Number of points in k:', nPk
         IF(na .NE. 1) THEN
            WRITE(*,*) 'INIT_CAMB_LINEAR: amin:', a(1)
            WRITE(*,*) 'INIT_CAMB_LINEAR: amax:', a(na)
         END IF
         WRITE(*,*) 'INIT_CAMB_LINEAR: Number of points in a:', na
      END IF

      IF (rebin) THEN

         ! TODO: Dedicated rebinnig 

         CALL fill_array_log(kmin_rebin, kmax_rebin, kkPk, nk)
         ALLOCATE(Pkk(nk, na))

         kPk = log(kPk)
         Pk = log(Pk)
        
         DO j = 1, na
            DO i = 1, nk
               Pkk(i, j) = find(log(kkPk(i)), kPk, Pk(:, j), nPk, &
                  iorder_rebin_CAMB, &
                  ifind_rebin_CAMB, &
                  iinterp_rebin_CAMB)
            END DO
         END DO

         DEALLOCATE(kPk)
         kPk = kkPk
         DEALLOCATE(Pk)
         Pk = exp(Pkk)

         IF (cosm%verbose) THEN
            WRITE(*,*) 'INIT_CAMB_LINEAR: Rebinning kmin [h/Mpc]:', kPk(1)
            WRITE(*,*) 'INIT_CAMB_LINEAR: Rebinning kmax [h/Mpc]:', kPk(nk)
            WRITE(*,*) 'INIT_CAMB_LINEAR: Number of rebinned points in k:', nk
         END IF

      END IF

      CALL init_linear(kPk, a, Pk, cosm)

      DEALLOCATE(kPk, Pk)     

      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_CAMB_LINEAR: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_CAMB_linear

   SUBROUTINE init_analytical_linear(cosm)

      ! Initialise an interpolator even if using an analyitcal power spectrum
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: k(:), a(:), Pk(:, :)
      INTEGER :: i
      REAL, PARAMETER :: kmin = kmin_plin
      REAL, PARAMETER :: kmax = kmax_plin
      INTEGER, PARAMETER :: nk = nk_plin
      INTEGER, PARAMETER :: na = 1

      CALL fill_array_log(kmin, kmax, k, nk)
      ALLOCATE(a(na), Pk(nk, na))
      a = 1.

      DO i = 1, nk
         Pk(i, 1) = plin(k(i), a(1), flag_matter, cosm)
      END DO
      
      CALL init_linear(k, a, Pk, cosm)

   END SUBROUTINE init_analytical_linear

   SUBROUTINE init_linear(k, a, Pk, cosm)

      ! Initialise an interpolator with values for the linear power spectrum
      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: Pk(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: nk, na

      nk = size(k)
      IF (nk /= size(Pk, 1)) THEN
         WRITE(*, *) 'INIT_LINEAR: Sizes', nk, size(Pk, 1)
         STOP 'INIT_LINEAR: Error, k and Pk arrays are not the same size'
      END IF    

      na = size(a)
      IF (na /= size(Pk, 2)) THEN
         WRITE(*, *) 'INIT_LINEAR: Sizes', na, size(Pk, 2)
         STOP 'INIT_LINEAR: Error, a and Pk arrays are not the same size'
      END IF
      
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_LINEAR: kmin [h/Mpc]:', k(1)
         WRITE (*, *) 'INIT_LINEAR: kmax [h/Mpc]:', k(nk)
         WRITE (*, *) 'INIT_LINEAR: nk:', nk
         WRITE (*, *) 'INIT_LINEAR: amin:', a(1)
         WRITE (*, *) 'INIT_LINEAR: amax:', a(na)
         WRITE (*, *) 'INIT_LINEAR: na:', na
         WRITE (*, *) 'INIT_LINEAR: Initialising interpolator'
      END IF

      CALL init_interpolator(k, a, Pk, cosm%plin, &
         iorder = iorder_interp_plin, &
         iextrap = iextrap_plin, &
         store = store_plin, &
         logx = .TRUE., &
         logy = .TRUE., &
         logf = .TRUE.)

      IF (cosm%is_normalised) cosm%A = 1.

      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_LINEAR: Done'
         WRITE (*, *)
      END IF

      ! Now the linear power interpolator is filled
      cosm%has_power = .TRUE.

   END SUBROUTINE init_linear

   SUBROUTINE init_external_linear_power_tables(cosm, k, a, Pk)

      ! TILMAN: Wrote this
      ! TODO: Rename to 'read_external_linear_power'
      ! TODO: This really is P(k) here, not Delta^2(k). Should always be Delta^2(k) despite variable being called Pk   
      ! TODO: Change cosm to be final argument to be consistent with other routines
      ! TODO: Really should remove cosm%k_plin, cosm%a_plin, cosm%Plin_array
      ! TODO: MUST talk to Tilman before changing this at all
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmolgy
      REAL, INTENT(IN) :: k(:)     ! Array of wavenumbers
      REAL, INTENT(IN) :: a(:)     ! Array of scale factors
      REAL, INTENT(IN) :: Pk(:, :) ! Array of P(k, a) (NOTE: Really P(k) here, not Delta^2(k))
      INTEGER :: nk, na, i

      ! Get sizes of arrays
      nk = size(k)
      na = size(a)
      IF (nk /= size(Pk, 1) .OR. na /= size(Pk, 2)) THEN
         WRITE(*, *) 'INIT_EXTERNAL_LINEAR_POWER_TABLES: Sizes of k, a, and Pk are inconsistent:', nk, na, shape(Pk)
         cosm%status = 1
         RETURN
      END IF
      IF (na == 1 .AND. a(na) /= 1.) THEN
         STOP 'INIT_EXTERNAL_LINEAR_POWER_TABLES: Error, if linear power is provided at a single scale factor then input must be at a=1.'
      ELSE IF (na > 1 .AND. na < 4) THEN
         STOP 'INIT_EXTERNAL_LINEAR_POWER_TABLES: Error, input linear power needs at least 4 scale factors for interpolation to work'
      END IF
      cosm%nk_plin = nk
      cosm%na_plin = na

      ! Fill internal cosm arrays from external power spectrum
      IF (.NOT. allocated(cosm%k_Plin)) ALLOCATE(cosm%k_Plin(nk))
      IF (.NOT. allocated(cosm%a_Plin)) ALLOCATE(cosm%a_Plin(na))
      IF (.NOT. allocated(cosm%Plin_array)) ALLOCATE(cosm%Plin_array(nk, na))
      cosm%k_Plin = k
      cosm%a_Plin = a
      FORALL (i = 1:nk) cosm%Plin_array(i, :) = Pk(i, :)*k(i)**3/(2.*pi**2)

      ! Set flags, assume this is called after assign_cosmology, but before init_cosmology
      cosm%itk = iTk_external
      cosm%norm_method = norm_none

   END SUBROUTINE init_external_linear_power_tables

   SUBROUTINE init_external_linear(cosm)

      ! TILMAN: Wrote this
      ! The purpose of this is *only* to init interpolators and set has_power
      ! TODO: Really should remove cosm%k_plin, cosm%a_plin, cosm%plin_array
      ! TODO: Could this replace or use init_linear?
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: nk, na, nk_pk, na_pk!, nk_pka, na_pka
      INTEGER :: plin_shape(2)

      ! Check k array is allocated properly
      IF (allocated(cosm%k_Plin)) THEN
         nk = size(cosm%k_Plin)
      ELSE
         WRITE (*, *) 'INIT_EXTERNAL_LINEAR: cosmology%k_Plin has not been allocated!'
         cosm%status = 1
         RETURN
      END IF

      ! Check a array is allocated properly
      IF (allocated(cosm%a_Plin)) THEN
         na = size(cosm%a_Plin)
      ELSE
         WRITE (*, *) 'INIT_EXTERNAL_LINEAR: cosmology%a_Plin has not been allocated!'
         cosm%status = 1
         RETURN
      END IF

      ! Check P(k,a) array is allocated properly
      IF (allocated(cosm%Plin_array)) THEN
         plin_shape = shape(cosm%Plin_array)
         nk_pk = plin_shape(1)
         na_pk = plin_shape(2)
      ELSE
         WRITE (*, *) 'INIT_EXTERNAL_LINEAR: cosmology%Plin_array has not been allocated!'
         cosm%status = 1
         RETURN
      END IF

      ! Check k sizes of arrays agree
      IF (nk /= nk_pk .OR. nk /= cosm%nk_plin) THEN
         WRITE (*,*) 'INIT_EXTERNAL_LINEAR: Sizes of cosmology%Plin_array, cosmology%k_plin, or cosmology%nk_plin are inconsistent:', nk_pk, nk, cosm%nk_plin
         cosm%status = 1
         RETURN
      END IF

      ! Check a sizes of arrays agree
      IF(na /= na_pk .OR. na /= cosm%na_plin) THEN
         WRITE (*, *) 'INIT_EXTERNAL_LINEAR: Sizes of cosmology%Plin_array, cosmology%a_plin, or cosmology%na_plin are inconsistent:', na_pk, na, cosm%na_plin
         cosm%status = 1
         RETURN
      END IF

      IF (cosm%scale_dependent_growth .AND. na == 1) THEN
         STOP 'INIT_EXTERNAL_LINEAR: Error, growth is scale dependent but linear spectrum only provided at a single scale factor'
      END IF

      CALL init_linear(cosm%k_plin, cosm%a_plin, cosm%Plin_array, cosm)
      DEALLOCATE(cosm%k_plin, cosm%a_plin, cosm%Plin_array)

   END SUBROUTINE init_external_linear

   SUBROUTINE random_AGN_temperature(cosm)

      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: logT_AGN
      REAL, PARAMETER :: logT_AGN_min = 7.5
      REAL, PARAMETER :: logT_AGN_max = 8.1

      logT_AGN = random_uniform(logT_AGN_min, logT_AGN_max)
      cosm%Theat = 10**logT_AGN

   END SUBROUTINE random_AGN_temperature

   SUBROUTINE random_cosmology(cosm)

      ! Generate some random cosmological parameters
      ! TODO: Use wb rather than wa? See https://arxiv.org/abs/2003.12116
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: Om_c_min = 0.1
      REAL, PARAMETER :: Om_c_max = 0.5
      REAL, PARAMETER :: Om_b_min = 0.02
      REAL, PARAMETER :: Om_b_max = 0.07
      REAL, PARAMETER :: h_min = 0.4
      REAL, PARAMETER :: h_max = 1.0
      REAL, PARAMETER :: n_min = 0.7
      REAL, PARAMETER :: n_max = 1.3
      REAL, PARAMETER :: w_min = -1.3
      REAL, PARAMETER :: w_max = -0.7
      REAL, PARAMETER :: wa_min=-1.
      REAL, PARAMETER :: wa_max=1.
      REAL, PARAMETER :: sig8_min = 0.6
      REAL, PARAMETER :: sig8_max = 0.9
      REAL, PARAMETER :: mnu_min = 0.01
      REAL, PARAMETER :: mnu_max = 0.60
      REAL, PARAMETER :: w_lim = 0.

      !CALL RNG_set(seed=0)

      cosm%h = random_uniform(h_min, h_max)

      cosm%Om_c = random_uniform(Om_c_min, Om_c_max)

      cosm%Om_b = random_uniform(Om_b_min, Om_b_max)

      cosm%Om_m = cosm%Om_c+cosm%Om_b

      ! Enforce flatness
      ! Note - need to have Om_w for dark enegry
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m

      cosm%ns = random_uniform(n_min, n_max)

      cosm%w = random_uniform(w_min, w_max)

      ! Ensure that there is no early dark energy (require w+wa ~< 0.)
      DO
         cosm%wa = random_uniform(wa_min, wa_max)
         IF((cosm%w+cosm%wa)<w_lim) EXIT
      END DO

      cosm%iw = iw_waCDM ! Fix iw = 3 for w(a)CDM

      cosm%sig8 = random_uniform(sig8_min, sig8_max)

      cosm%m_nu = random_uniform(mnu_min, mnu_max)
 
   END SUBROUTINE random_cosmology

   SUBROUTINE random_waCDM_cosmology(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_cosmology(cosm)
      cosm%m_nu = 0. ! Fix massive neutrinos to zero

   END SUBROUTINE random_waCDM_cosmology

   SUBROUTINE random_wCDM_cosmology(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_waCDM_cosmology(cosm)
      cosm%wa = 0.      ! Time-independent w
      cosm%iw = iw_wCDM ! Fix for constant-w wCDM

   END SUBROUTINE random_wCDM_cosmology

   SUBROUTINE random_LCDM_cosmology(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_wCDM_cosmology(cosm)
      cosm%w = -1.      ! Fix w = -1
      cosm%iw = iw_LCDM ! Fix to vacuum dark energy

   END SUBROUTINE random_LCDM_cosmology

   SUBROUTINE random_nuLCDM_cosmology(cosm)

      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_cosmology(cosm)
      cosm%w = -1.
      cosm%wa = 0.
      cosm%iw = iw_LCDM

   END SUBROUTINE random_nuLCDM_cosmology

   SUBROUTINE get_CosmicEmu_h(cosm)

      ! TODO: This should be removed eventually
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: u
      INTEGER :: i
      REAL :: h
      CHARACTER(len=256) :: crap
      CHARACTER(len=256), PARAMETER :: params = params_CosmicEmu
      CHARACTER(len=256), PARAMETER :: output = output_CosmicEmu
      CHARACTER(len=256), PARAMETER :: exe = exe_CosmicEmu
      INTEGER, PARAMETER :: hl = hl_CosmicEmu
      INTEGER, PARAMETER :: nh = nh_CosmicEmu
      REAL, PARAMETER :: z = 0.

      ! Remove previous parameter and power file
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(params))
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(output))

      ! Write a new parameter file
      OPEN (newunit=u, file=params)
      WRITE (u, fmt='(A20,7F10.5)') trim(output), cosm%Om_m*cosm%h**2, cosm%Om_b*cosm%h**2, cosm%ns, cosm%sig8, cosm%w, z
      CLOSE (u)

      ! Run emu
      IF(present_and_correct(cosm%verbose)) THEN
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params))
      ELSE
         CALL EXECUTE_COMMAND_LINE(trim(exe)//' '//trim(params)//' > /dev/null')
      END IF

      ! Read in data file
      OPEN (newunit=u, file=output)
      DO i = 1, nh
         IF (i == hl) THEN
            READ (u, *) crap, crap, crap, crap, crap, crap, crap, crap, h
            cosm%h = h
            EXIT
         ELSE
            READ(u, *)
         END IF
      END DO
      CLOSE (u)

   END SUBROUTINE get_CosmicEmu_h

   SUBROUTINE random_Cosmic_Emu_cosmology(cosm)

      ! Generate some random cosmological parameters
      ! h is fixed by the distance to the CMB last-scattering surface
      ! Described below equation (8) of https://arxiv.org/pdf/0902.0429.pdf
      ! TODO: Should replace get_CosmicEmu_h with my own calculation
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: dlss, zlss, alss, rs 
      REAL :: wm, wb
      REAL, PARAMETER :: om_m_min = 0.120
      REAL, PARAMETER :: om_m_max = 0.155
      REAL, PARAMETER :: om_b_min = 0.0215
      REAL, PARAMETER :: om_b_max = 0.0235
      REAL, PARAMETER :: n_min = 0.85
      REAL, PARAMETER :: n_max = 1.05
      REAL, PARAMETER :: w_min = -1.3
      REAL, PARAMETER :: w_max = -0.7
      REAL, PARAMETER :: sig8_min = 0.616
      REAL, PARAMETER :: sig8_max = 0.9
      REAL, PARAMETER :: wm_mid = 0.135
      REAL, PARAMETER :: la = 302.4
      LOGICAL, PARAMETER :: calculate_h = .FALSE.

      ! Uniform-random sampling of cosmic-emu parameter space
      wm = random_uniform(om_m_min, om_m_max)
      wb = random_uniform(om_b_min, om_b_max)
      cosm%ns = random_uniform(n_min, n_max)
      cosm%w = random_uniform(w_min, w_max)
      cosm%sig8 = random_uniform(sig8_min, sig8_max)

      ! Deal with Hubble parameter separately
      IF (calculate_h) THEN
         zlss = z_CMB(cosm)
         alss = scale_factor_z(zlss)
         dlss = comoving_distance(alss, cosm)
         rs = r_sound(cosm)
         WRITE(*, *) 'l_A (target):', la
         WRITE(*, *) 'l_A (this cosmology)', pi*dlss/rs
         STOP 'RANDOM_COSMIC_EMU_COSMOLOGY: Error, this needs to be solved numerically'
      ELSE
         cosm%h = 1. ! Need to set to unity here
         cosm%Om_m = wm
         cosm%Om_b = wb
         CALL get_CosmicEmu_h(cosm) ! Reset h
      END IF

      ! Convert to HMcode primary parameters
      cosm%Om_m = wm/cosm%h**2
      cosm%Om_b = wb/cosm%h**2

      ! Dark energy stuff
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m
      cosm%iw = iw_wCDM ! Constant w models only in Cosmic Emu

      ! Consistent normalisation for As
      cosm%kpiv = 0.05/cosm%h

   END SUBROUTINE random_Cosmic_Emu_cosmology

   SUBROUTINE random_Franken_Emu_cosmology(cosm)

      ! Generate some random cosmological parameters within the FrankenEmu cube
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: wm, wb
      REAL, PARAMETER :: om_m_min = 0.120
      REAL, PARAMETER :: om_m_max = 0.155
      REAL, PARAMETER :: om_b_min = 0.0215
      REAL, PARAMETER :: om_b_max = 0.0235
      REAL, PARAMETER :: n_min = 0.85
      REAL, PARAMETER :: n_max = 1.05
      REAL, PARAMETER :: h_min = 0.55
      REAL, PARAMETER :: h_max = 0.85
      REAL, PARAMETER :: w_min = -1.3
      REAL, PARAMETER :: w_max = -0.7
      REAL, PARAMETER :: sig8_min = 0.616
      REAL, PARAMETER :: sig8_max = 0.9

      ! Uniform-random sampling of FrankenEmu parameter space
      cosm%h = random_uniform(h_min, h_max)
      wm = random_uniform(om_m_min, om_m_max)
      wb = random_uniform(om_b_min, om_b_max)
      cosm%ns = random_uniform(n_min, n_max)
      cosm%w = random_uniform(w_min, w_max)
      cosm%sig8 = random_uniform(sig8_min, sig8_max)

      ! Convert to HMcode primary parameters
      cosm%Om_m = wm/cosm%h**2
      cosm%Om_b = wb/cosm%h**2

      ! Enforce flatness
      ! Note - need to have Om_w for dark enegry
      cosm%iw = iw_wCDM ! Constant w models only in Franken Emu
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m

      ! For consistent As normalisation
      cosm%kpiv = 0.05/cosm%h

   END SUBROUTINE random_Franken_Emu_cosmology

   SUBROUTINE random_Mira_Titan_cosmology(cosm)

      ! Generate some random cosmological parameters for the Mira Titan hypercube
      ! TODO: Use wb rather than wa? See https://arxiv.org/abs/2003.12116
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: om_m, om_b, om_nu
      REAL, PARAMETER :: om_m_min = 0.120
      REAL, PARAMETER :: om_m_max = 0.155
      REAL, PARAMETER :: om_b_min = 0.0215
      REAL, PARAMETER :: om_b_max = 0.0235
      REAL, PARAMETER :: om_nu_min = 1e-4 ! Not 0 to avoid too light neutrinos
      REAL, PARAMETER :: om_nu_max = 1e-2
      REAL, PARAMETER :: n_min = 0.85
      REAL, PARAMETER :: n_max = 1.05
      REAL, PARAMETER :: h_min = 0.55
      REAL, PARAMETER :: h_max = 0.85
      REAL, PARAMETER :: w_min = -1.3
      REAL, PARAMETER :: w_max = -0.7
      REAL, PARAMETER :: wa_min = -1.73
      REAL, PARAMETER :: wa_max = 1.28
      REAL, PARAMETER :: sig8_min = 0.7
      REAL, PARAMETER :: sig8_max = 0.9

      ! Uniform random sampling
      cosm%h = random_uniform(h_min, h_max)
      om_m = random_uniform(om_m_min, om_m_max)
      om_b = random_uniform(om_b_min, om_b_max)
      om_nu = random_uniform(om_nu_min, om_nu_max)
      cosm%ns = random_uniform(n_min, n_max)
      cosm%w = random_uniform(w_min, w_max)
      cosm%sig8 = random_uniform(sig8_min, sig8_max) 
   
      ! Enforce 0.3 <= (-w0-wa)^(1/4) as in Mira Titan paper
      ! TODO: Use wb rather than wa?
      DO
         cosm%wa = random_uniform(wa_min, wa_max)
         IF (0.0081 <= -cosm%w-cosm%wa .AND. 2.769 >= -cosm%w-cosm%wa) EXIT
      END DO

      ! Convert to my primary parameters
      cosm%Om_m = om_m/cosm%h**2   
      cosm%Om_b = om_b/cosm%h**2 
      cosm%m_nu = neutrino_constant(cosm)*om_nu

      ! Enforce flatness, ensure Omega_w is used for dark energy, Omega_v = 0
      cosm%iw = iw_waCDM ! w(a)CDM models in Mira Titan
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m  

      ! For consistent As normalisation
      cosm%kpiv = 0.05/cosm%h

   END SUBROUTINE random_Mira_Titan_cosmology

   SUBROUTINE random_Euclid_cosmology(cosm)

      ! Generate random cosmological parameters for the Euclid Emulator hypercube
      ! Table (2) from https://arxiv.org/pdf/2010.11288.pdf
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: Om_b_min = 0.04
      REAL, PARAMETER :: Om_b_max = 0.06
      REAL, PARAMETER :: Om_m_min = 0.24
      REAL, PARAMETER :: Om_m_max = 0.40
      REAL, PARAMETER :: m_nu_min = 0.01 ! Not 0 to avoid too-light neutrinos
      REAL, PARAMETER :: m_nu_max = 0.15
      REAL, PARAMETER :: ns_min = 0.92
      REAL, PARAMETER :: ns_max = 1.00
      REAL, PARAMETER :: h_min = 0.61
      REAL, PARAMETER :: h_max = 0.73
      REAL, PARAMETER :: w_min = -1.3
      REAL, PARAMETER :: w_max = -0.7
      REAL, PARAMETER :: wa_min = -0.7
      REAL, PARAMETER :: wa_max = 0.5
      REAL, PARAMETER :: As_min = 1.7e-9
      REAL, PARAMETER :: As_max = 2.5e-9

      ! Randomly generate primary parameters
      cosm%Om_b = random_uniform(Om_b_min, Om_b_max)
      cosm%Om_m = random_uniform(Om_m_min, Om_m_max)
      cosm%m_nu = random_uniform(m_nu_min, m_nu_max)
      cosm%h = random_uniform(h_min, h_max)  
      cosm%ns = random_uniform(ns_min, ns_max)
      cosm%w = random_uniform(w_min, w_max)
      cosm%wa = random_uniform(wa_min, wa_max)
      cosm%As = random_uniform(As_min, As_max)
      
      ! Enforce flatness, ensure Omega_w is used for w(a)CDM dark energy, Omega_v = 0
      cosm%iw = iw_waCDM 
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m     

      ! CMB temperature [K]
      cosm%T_CMB = 2.7255 

      ! Normalisation; Ensure kpiv = 0.05/Mpc; NOTE: My units are h/Mpc
      cosm%norm_method = norm_As
      cosm%kpiv = 0.05/cosm%h

   END SUBROUTINE random_Euclid_cosmology

   SUBROUTINE random_BACCO_cosmology(cosm)

      ! Generate random cosmological parameters for the BACCO Emulator hypercube
      ! http://www.dipc.org/bacco/emulator.html
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: Om_cold_min = 0.23
      REAL, PARAMETER :: Om_cold_max = 0.40
      REAL, PARAMETER :: sig8_cold_min = 0.73
      REAL, PARAMETER :: sig8_cold_max = 0.90
      REAL, PARAMETER :: Om_b_min = 0.04
      REAL, PARAMETER :: Om_b_max = 0.06
      REAL, PARAMETER :: ns_min = 0.92
      REAL, PARAMETER :: ns_max = 1.01
      REAL, PARAMETER :: h_min = 0.6
      REAL, PARAMETER :: h_max = 0.8
      REAL, PARAMETER :: m_nu_min = 0.01 ! Not 0 to avoid too light neutrinos
      REAL, PARAMETER :: m_nu_max = 0.40     
      REAL, PARAMETER :: w_min = -1.15
      REAL, PARAMETER :: w_max = -0.85
      REAL, PARAMETER :: wa_min = -0.3
      REAL, PARAMETER :: wa_max = 0.3
      REAL :: sig8_cold, Om_cold, f_nu, Om_nu

      ! Randomly generate BACCO primary parameters
      Om_cold = random_uniform(Om_cold_min, Om_cold_max)
      sig8_cold = random_uniform(sig8_cold_min, sig8_cold_max)
      cosm%Om_b = random_uniform(Om_b_min, Om_b_max)
      cosm%ns = random_uniform(ns_min, ns_max)
      cosm%h = random_uniform(h_min, h_max)
      cosm%m_nu = random_uniform(m_nu_min, m_nu_max)
      cosm%w = random_uniform(w_min, w_max)
      
      ! Avoid dark energy equation of state crossing -1   
      !IF (cosm%w == -1.) THEN
      !   cosm%wa = 0.
      !ELSE IF (cosm%w < -1.) THEN
      !   cosm%wa = random_uniform(wa_min, -cosm%w-1.)
      !ELSE
      !   cosm%wa = random_uniform(-1.-cosm%w, wa_max)
      !END IF
      IF (cosm%w <= -1.) THEN
         cosm%wa = 0.
      ELSE
         cosm%wa = random_uniform(-1.-cosm%w, wa_max)
      END IF

      ! CMB temperature [K]
      cosm%T_CMB = 2.7255

      ! Convert to my primary parameters
      Om_nu = cosm%m_nu/(neutrino_constant(cosm)*cosm%h**2) ! Change CMB temperature before this
      f_nu = Om_nu/(Om_nu+Om_cold)
      cosm%Om_m = Om_cold/(1.-f_nu)
      cosm%sig8 = sig8_cold/(1.-f_nu)
        
      ! Enforce flatness, ensure Omega_w is used for w(a)CDM dark energy, Omega_v = 0
      cosm%iw = iw_waCDM
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m

   END SUBROUTINE random_BACCO_cosmology

   SUBROUTINE random_NGenHALOFIT_cosmology(cosm)

      ! Generate random cosmological parameters for the NGenHALOFIT hypercube
      ! From abstract of https://arxiv.org/pdf/1807.00040.pdf
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: wm, wc, wb
      REAL, PARAMETER :: w_min = -1.05
      REAL, PARAMETER :: w_max = -0.95
      REAL, PARAMETER :: wa_min = -0.4
      REAL, PARAMETER :: wa_max = 0.4
      REAL, PARAMETER :: Om_m_min = 0.21
      REAL, PARAMETER :: Om_m_max = 0.40
      REAL, PARAMETER :: wc_min = 0.1
      REAL, PARAMETER :: wc_max = 0.13
      REAL, PARAMETER :: wb_min = 0.02
      REAL, PARAMETER :: wb_max = 0.024
      REAL, PARAMETER :: ns_min = 0.85
      REAL, PARAMETER :: ns_max = 1.05
      REAL, PARAMETER :: As_min = 1.72e-9
      REAL, PARAMETER :: As_max = 2.58e-9
      REAL, PARAMETER :: nrun_min = -0.2
      REAL, PARAMETER :: nrun_max = 0.2

      ! Randomly generate primary parameters
      cosm%w = random_uniform(w_min, w_max)
      cosm%wa = random_uniform(wa_min, wa_max)
      cosm%Om_m = random_uniform(Om_m_min, Om_m_max)
      wc = random_uniform(wc_min, wc_max)
      wb = random_uniform(wb_min, wb_max)
      cosm%ns = random_uniform(ns_min, ns_max)
      cosm%As = random_uniform(As_min, As_max)
      cosm%nrun = random_uniform(nrun_min, nrun_max)
   
      ! Convert to my primary parameters
      wm = wc+wb
      cosm%h = sqrt(wm/cosm%Om_m)
      cosm%Om_b = wb/cosm%h**2
      
      ! Enforce flatness, ensure Omega_w is used for w(a)CDM dark energy, Omega_v = 0
      cosm%iw = iw_waCDM 
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m

      ! Normalisation; Ensure kpiv = 0.05/Mpc; NOTE: My units are h/Mpc
      cosm%norm_method = norm_As
      cosm%kpiv = 0.05/cosm%h

      ! Other cosmological parameters
      cosm%T_CMB = 2.726
      cosm%neff = 3.040

   END SUBROUTINE random_NGenHALOFIT_cosmology

   SUBROUTINE random_Dark_Quest_cosmology(cosm)

      ! Generate random cosmological parameters for the Dark Quest hypercube
      ! Equations (25) from https://arxiv.org/pdf/1811.09504.pdf
      USE random_numbers
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: lnAs, wb, wc, wnu, wm
      REAL, PARAMETER :: wb_min = 0.0211375
      REAL, PARAMETER :: wb_max = 0.0233625
      REAL, PARAMETER :: wc_min = 0.10782
      REAL, PARAMETER :: wc_max = 0.13178
      REAL, PARAMETER :: Om_w_min = 0.54752
      REAL, PARAMETER :: Om_w_max = 0.82128
      REAL, PARAMETER :: lnAs_min = 2.4752
      REAL, PARAMETER :: lnAs_max = 3.7128
      REAL, PARAMETER :: ns_min = 0.916275
      REAL, PARAMETER :: ns_max = 1.012725
      REAL, PARAMETER :: w_min = -1.2
      REAL, PARAMETER :: w_max = -0.8

      ! Randomly generate primary parameters
      wb = random_uniform(wb_min, wb_max)
      wc = random_uniform(wc_min, wc_max)
      cosm%Om_w = random_uniform(Om_w_min, Om_w_max) 
      lnAs = random_uniform(lnAs_min, lnAs_max) 
      cosm%ns = random_uniform(ns_min, ns_max)
      cosm%w = random_uniform(w_min, w_max)

      ! Fixed neutrino density
      wnu = 0.00064
      
      ! Enforce flatness, ensure Omega_w is used for wCDM dark energy, Omega_v = 0
      cosm%iw = iw_wCDM 
      cosm%Om_v = 0.
      cosm%Om_m = 1.-cosm%Om_w
      wm = wc+wb+wnu
      cosm%h = sqrt(wm/cosm%Om_m)
      cosm%As = exp(lnAs)/1e10
      cosm%Om_b = wb/cosm%h**2  

      ! CMB temperature [K]
      !cosm%T_CMB = 2.7255 

      ! Neutrino mass (only after T_CMB has been set)
      cosm%m_nu = wnu*neutrino_constant(cosm)

      ! Normalisation; Ensure kpiv = 0.05/Mpc; NOTE: My units are h/Mpc
      cosm%norm_method = norm_As
      cosm%kpiv = 0.05/cosm%h

   END SUBROUTINE random_Dark_Quest_cosmology

   SUBROUTINE Cosmic_Emu_node_cosmology(node, cosm)

      INTEGER, INTENT(IN) :: node
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: om_m, om_b

      ! M000 -> M037 of Cosmic_Emu (Table 1 in 1304.7849)
      IF (node == 0) THEN
         ! M000 (not included in emulator construction)
         om_m = 0.1296
         om_b = 0.0224
         cosm%ns = 0.9700
         cosm%w = 1.000
         cosm%sig8 = 0.8000
         cosm%h = 0.7200
      ELSE IF (node == 1) THEN
         ! M001
         om_m = 0.1539
         om_b = 0.0231
         cosm%ns = 0.9468
         cosm%w = 0.816
         cosm%sig8 = 0.8161
         cosm%h = 0.5977
      ELSE IF (node == 2) THEN
         ! M002
         om_m = 0.1460
         om_b = 0.0227
         cosm%ns = 0.8952
         cosm%w = 0.758
         cosm%sig8 = 0.8548
         cosm%h = 0.5970
      ELSE IF (node == 3) THEN
         ! M003
         om_m = 0.1324
         om_b = 0.0235
         cosm%ns = 0.9984
         cosm%w = 0.874
         cosm%sig8 = 0.8484
         cosm%h = 0.6763
      ELSE IF (node == 4) THEN
         ! M004
         om_m = 0.1381
         om_b = 0.0227
         cosm%ns = 0.9339
         cosm%w = 1.087
         cosm%sig8 = 0.7000
         cosm%h = 0.7204
      ELSE IF (node == 5) THEN
         ! M005
         om_m = 0.1358
         om_b = 0.0216
         cosm%ns = 0.9726
         cosm%w = 1.242
         cosm%sig8 = 0.8226
         cosm%h = 0.7669
      ELSE IF (node == 6) THEN
         ! M006
         om_m = 0.1516
         om_b = 0.0229
         cosm%ns = 0.9145
         cosm%w = 1.223
         cosm%sig8 = 0.6705
         cosm%h = 0.7040
      ELSE IF (node == 7) THEN
         ! M007
         om_m = 0.1268
         om_b = 0.0223
         cosm%ns = 0.9210
         cosm%w = 0.70001 ! Changed to avoid problems
         cosm%sig8 = 0.7474
         cosm%h = 0.6189
      ELSE IF (node == 8) THEN
         ! M008
         om_m = 0.1448
         om_b = 0.0223
         cosm%ns = 0.9855
         cosm%w = 1.203
         cosm%sig8 = 0.8090
         cosm%h = 0.7218
      ELSE IF (node == 9) THEN
         ! M009
         om_m = 0.1392
         om_b = 0.0234
         cosm%ns = 0.9790
         cosm%w = 0.739
         cosm%sig8 = 0.6692
         cosm%h = 0.6127
      ELSE IF (node == 10) THEN
         ! M010
         om_m = 0.1403
         om_b = 0.0218
         cosm%ns = 0.8565
         cosm%w = 0.990
         cosm%sig8 = 0.7556
         cosm%h = 0.6695
      ELSE IF (node == 11) THEN
         ! M011
         om_m = 0.1437
         om_b = 0.0234
         cosm%ns = 0.8823
         cosm%w = 1.126
         cosm%sig8 = 0.7276
         cosm%h = 0.7177
      ELSE IF (node == 12) THEN
         ! M012
         om_m = 0.1223
         om_b = 0.0225
         cosm%ns = 1.0048
         cosm%w = 0.971
         cosm%sig8 = 0.6271
         cosm%h = 0.7396
      ELSE IF (node == 13) THEN
         ! M013
         om_m = 0.1482
         om_b = 0.0221
         cosm%ns = 0.9597
         cosm%w = 0.855
         cosm%sig8 = 0.6508
         cosm%h = 0.6107
      ELSE IF (node == 14) THEN
         ! M014
         om_m = 0.1471
         om_b = 0.0233
         cosm%ns = 1.0306
         cosm%w = 1.010
         cosm%sig8 = 0.7075
         cosm%h = 0.6688
      ELSE IF (node == 15) THEN
         ! M015
         om_m = 0.1415
         om_b = 0.0230
         cosm%ns = 1.0177
         cosm%w = 1.281
         cosm%sig8 = 0.7692
         cosm%h = 0.7737
      ELSE IF (node == 16) THEN
         ! M016
         om_m = 0.1245
         om_b = 0.0218
         cosm%ns = 0.9403
         cosm%w = 1.145
         cosm%sig8 = 0.7437
         cosm%h = 0.7929
      ELSE IF (node == 17) THEN
         ! M017
         om_m = 0.1426
         om_b = 0.0215
         cosm%ns = 0.9274
         cosm%w = 0.893
         cosm%sig8 = 0.6865
         cosm%h = 0.6305
      ELSE IF (node == 18) THEN
         ! M018
         om_m = 0.1313
         om_b = 0.0216
         cosm%ns = 0.8887
         cosm%w = 1.029
         cosm%sig8 = 0.6440
         cosm%h = 0.7136
      ELSE IF (node == 19) THEN
         ! M019
         om_m = 0.1279
         om_b = 0.0232
         cosm%ns = 0.8629
         cosm%w = 1.184
         cosm%sig8 = 0.6159
         cosm%h = 0.8120
      ELSE IF (node == 20) THEN
         ! M020
         om_m = 0.1290
         om_b = 0.0220
         cosm%ns = 1.0242
         cosm%w = 0.797
         cosm%sig8 = 0.7972
         cosm%h = 0.6442
      ELSE IF (node == 21) THEN
         ! M021
         om_m = 0.1335
         om_b = 0.0221
         cosm%ns = 1.0371
         cosm%w = 1.165
         cosm%sig8 = 0.6563
         cosm%h = 0.7601
      ELSE IF (node == 22) THEN
         ! M022
         om_m = 0.1505
         om_b = 0.0225
         cosm%ns = 1.0500
         cosm%w = 1.107
         cosm%sig8 = 0.7678
         cosm%h = 0.6736
      ELSE IF (node == 23) THEN
         ! M023
         om_m = 0.1211
         om_b = 0.0220
         cosm%ns = 0.9016
         cosm%w = 1.261
         cosm%sig8 = 0.6664
         cosm%h = 0.8694
      ELSE IF (node == 24) THEN
         ! M024
         om_m = 0.1302
         om_b = 0.0226
         cosm%ns = 0.9532
         cosm%w = 1.300
         cosm%sig8 = 0.6644
         cosm%h = 0.8380
      ELSE IF (node == 25) THEN
         ! M025
         om_m = 0.1494
         om_b = 0.0217
         cosm%ns = 1.0113
         cosm%w = 0.719
         cosm%sig8 = 0.7398
         cosm%h = 0.5724
      ELSE IF (node == 26) THEN
         ! M026
         om_m = 0.1347
         om_b = 0.0232
         cosm%ns = 0.9081
         cosm%w = 0.952
         cosm%sig8 = 0.7995
         cosm%h = 0.6931
      ELSE IF (node == 27) THEN
         ! M027
         om_m = 0.1369
         om_b = 0.0224
         cosm%ns = 0.8500
         cosm%w = 0.836
         cosm%sig8 = 0.7111
         cosm%h = 0.6387
      ELSE IF (node == 28) THEN
         ! M028
         om_m = 0.1527
         om_b = 0.0222
         cosm%ns = 0.8694
         cosm%w = 0.932
         cosm%sig8 = 0.8068
         cosm%h = 0.6189
      ELSE IF (node == 29) THEN
         ! M029
         om_m = 0.1256
         om_b = 0.0228
         cosm%ns = 1.0435
         cosm%w = 0.913
         cosm%sig8 = 0.7087
         cosm%h = 0.7067
      ELSE IF (node == 30) THEN
         ! M030
         om_m = 0.1234
         om_b = 0.0230
         cosm%ns = 0.8758
         cosm%w = 0.777
         cosm%sig8 = 0.6739
         cosm%h = 0.6626
      ELSE IF (node == 31) THEN
         ! M031
         om_m = 0.1550
         om_b = 0.0219
         cosm%ns = 0.9919
         cosm%w = 1.068
         cosm%sig8 = 0.7041
         cosm%h = 0.6394
      ELSE IF (node == 32) THEN
         ! M032
         om_m = 0.1200
         om_b = 0.0229
         cosm%ns = 0.9661
         cosm%w = 1.048
         cosm%sig8 = 0.7556
         cosm%h = 0.7901
      ELSE IF (node == 33) THEN
         ! M033
         om_m = 0.1399
         om_b = 0.0225
         cosm%ns = 1.0407
         cosm%w = 1.147
         cosm%sig8 = 0.8645
         cosm%h = 0.7286
      ELSE IF (node == 34) THEN
         ! M034
         om_m = 0.1497
         om_b = 0.0227
         cosm%ns = 0.9239
         cosm%w = 1.000
         cosm%sig8 = 0.8734
         cosm%h = 0.6510
      ELSE IF (node == 35) THEN
         ! M035
         om_m = 0.1485
         om_b = 0.0221
         cosm%ns = 0.9604
         cosm%w = 0.853
         cosm%sig8 = 0.8822
         cosm%h = 0.6100
      ELSE IF (node == 36) THEN
         ! M036
         om_m = 0.1216
         om_b = 0.0233
         cosm%ns = 0.9387
         cosm%w = 0.706
         cosm%sig8 = 0.8911
         cosm%h = 0.6421
      ELSE IF (node == 37) THEN
         ! M037
         om_m = 0.1495
         om_b = 0.0228
         cosm%ns = 1.0233
         cosm%w = 1.294
         cosm%sig8 = 0.8999 ! Moved off boundary
         cosm%h = 0.7313
      ELSE
         STOP 'COSMIC_EMU_NODE_COSMOLOGY: Error, node specified incorrectly'
      END IF

      cosm%w = -cosm%w
      cosm%iw = iw_wCDM
      cosm%Om_m = om_m/cosm%h**2
      cosm%Om_b = om_b/cosm%h**2
      cosm%Om_w = 1.-cosm%Om_m
      cosm%Om_v = 0.
      cosm%kpiv = 0.05/cosm%h

   END SUBROUTINE Cosmic_Emu_node_cosmology

   SUBROUTINE Franken_Emu_node_cosmology(node, cosm)

      ! Node cosmolgies from Franken Emu, almost the same as Cosmic Emu apart from M023
      INTEGER, INTENT(IN) :: node
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: om_m, om_b

      IF (node == 23) THEN
         ! M023
         om_m = 0.1211
         om_b = 0.0220
         cosm%ns = 0.9016
         cosm%w = -1.261
         cosm%sig8 = 0.6664
         cosm%h = 0.8500 ! Have to round down from 0.8694 to 0.85 for FrankenEmu shrunken space
         cosm%Om_m = om_m/cosm%h**2
         cosm%Om_b = om_b/cosm%h**2
         cosm%Om_w = 1.-cosm%Om_m
         cosm%Om_v = 0.
         cosm%iw = iw_wCDM
         cosm%kpiv = 0.05/cosm%h
      ELSE
         CALL Cosmic_Emu_node_cosmology(node, cosm)
      END IF

   END SUBROUTINE Franken_Emu_node_cosmology

   SUBROUTINE Mira_Titan_node_cosmology(node, cosm)

      ! Node cosmologies from Mira Titan
      INTEGER, INTENT(IN) :: node
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: om_m, om_b, om_nu

      ! M000 -> M036 of Mira Titan (m_nu = 0 eV for M000 -> M010; Table 3 in 1705.03388)
      IF (node == 0) THEN
         ! M000 (WMAP7 - not included in emulator construction)
         om_m = 0.1335
         om_b = 0.02258
         cosm%sig8 = 0.8
         cosm%h = 0.71
         cosm%ns = 0.963
         cosm%w = -1.0
         cosm%wa = 0.0
         om_nu = 0.
      ELSE IF (node == 1) THEN
         ! M001
         ! Something seems funnny with this cosmology, the halo-model prediction is awful
         ! This model has w0 ~ -wa ~ 0.7, so w(a=0) = w0+wa ~ 0 so DE scales almost like matter at early times
         om_m = 0.1472
         om_b = 0.02261
         cosm%sig8 = 0.8778
         cosm%h = 0.6167
         cosm%ns = 0.9611
         cosm%w = -0.7000
         cosm%wa = 0.67220
         om_nu = 0.
      ELSE IF (node == 2) THEN
         ! M002
         om_m = 0.1356
         om_b = 0.02328
         cosm%sig8 = 0.8556
         cosm%h = 0.7500
         cosm%ns = 1.0500
         cosm%w = -1.0330
         cosm%wa = 0.91110
         om_nu = 0.
      ELSE IF (node == 3) THEN
         ! M003
         om_m = 0.1550
         om_b = 0.02194
         cosm%sig8 = 0.9000
         cosm%h = 0.7167
         cosm%ns = 0.8944
         cosm%w = -1.1000
         cosm%wa = -0.28330
         om_nu = 0.
      ELSE IF (node == 4) THEN
         ! M004
         om_m = 0.1239
         om_b = 0.02283
         cosm%sig8 = 0.7889
         cosm%h = 0.5833
         cosm%ns = 0.8722
         cosm%w = -1.1670
         cosm%wa = 1.15000
         om_nu = 0.
      ELSE IF (node == 5) THEN
         ! M005
         om_m = 0.1433
         om_b = 0.02350
         cosm%sig8 = 0.7667
         cosm%h = 0.8500
         cosm%ns = 0.9833
         cosm%w = -1.2330
         cosm%wa = -0.04445
         om_nu = 0.
      ELSE IF (node == 6) THEN
         ! M006
         om_m = 0.1317
         om_b = 0.021501 ! Moved off boundary
         cosm%sig8 = 0.8333
         cosm%h = 0.5500
         cosm%ns = 0.9167
         cosm%w = -0.7667
         cosm%wa = 0.19440
         om_nu = 0.
      ELSE IF (node == 7) THEN
         ! M007
         om_m = 0.1511
         om_b = 0.02217
         cosm%sig8 = 0.8111
         cosm%h = 0.8167
         cosm%ns = 1.0280
         cosm%w = -0.8333
         cosm%wa = -1.0000
         om_nu = 0.
      ELSE IF (node == 8) THEN
         ! M008
         om_m = 0.1200
         om_b = 0.02306
         cosm%sig8 = 0.7000
         cosm%h = 0.6833
         cosm%ns = 1.0060
         cosm%w = -0.9000
         cosm%wa = 0.43330
         om_nu = 0.
      ELSE IF (node == 9) THEN
         ! M009
         om_m = 0.1394
         om_b = 0.02172
         cosm%sig8 = 0.7444
         cosm%h = 0.6500
         cosm%ns = 0.8500
         cosm%w = -0.9667
         cosm%wa = -0.76110
         om_nu = 0.
      ELSE IF (node == 10) THEN
         ! M010 (last cosmology with massless neutrinos)
         om_m = 0.1278
         om_b = 0.02239
         cosm%sig8 = 0.7222
         cosm%h = 0.7833
         cosm%ns = 0.9389
         cosm%w = -1.3000
         cosm%wa = -0.52220
         om_nu = 0.
      ELSE IF (node == 11) THEN
         ! M011 (first cosmology with massive neutrinos)
         om_m = 0.1227
         om_b = 0.0220
         cosm%sig8 = 0.7151
         cosm%h = 0.5827
         cosm%ns = 0.9357
         cosm%w = -1.0821
         cosm%wa = 1.0646
         om_nu = 0.000345
      ELSE IF (node == 12) THEN
         ! M012
         om_m = 0.1241
         om_b = 0.0224
         cosm%sig8 = 0.7472
         cosm%h = 0.8315
         cosm%ns = 0.8865
         cosm%w = -1.2325
         cosm%wa = -0.7646
         om_nu = 0.001204
      ELSE IF (node == 13) THEN
         ! M013
         om_m = 0.1534
         om_b = 0.0232
         cosm%sig8 = 0.8098
         cosm%h = 0.7398
         cosm%ns = 0.8706
         cosm%w = -1.2993
         cosm%wa = 1.2236
         om_nu = 0.003770
      ELSE IF (node == 14) THEN
         ! M014
         om_m = 0.1215
         om_b = 0.0215
         cosm%sig8 = 0.8742
         cosm%h = 0.5894
         cosm%ns = 1.0151
         cosm%w = -0.7281
         cosm%wa = -0.2088
         om_nu = 0.001752
      ELSE IF (node == 15) THEN
         ! M015
         om_m = 0.1250
         om_b = 0.0224
         cosm%sig8 = 0.8881
         cosm%h = 0.6840
         cosm%ns = 0.8638
         cosm%w = -1.0134
         cosm%wa = 0.0415
         om_nu = 0.002789
      ELSE IF (node == 16) THEN
         ! M016
         om_m = 0.1499
         om_b = 0.0223
         cosm%sig8 = 0.7959
         cosm%h = 0.6452
         cosm%ns = 1.0219
         cosm%w = -1.0139
         cosm%wa = 0.9434
         om_nu = 0.002734
      ELSE IF (node == 17) THEN
         ! M017
         om_m = 0.1206
         om_b = 0.0215
         cosm%sig8 = 0.7332
         cosm%h = 0.7370
         cosm%ns = 1.0377
         cosm%w = -0.9472
         cosm%wa = -0.9897
         om_nu = 0.000168
      ELSE IF (node == 18) THEN
         ! M018
         om_m = 0.1544
         om_b = 0.0217
         cosm%sig8 = 0.7982
         cosm%h = 0.6489
         cosm%ns = 0.9026
         cosm%w = -0.7091
         cosm%wa = 0.6409
         om_nu = 0.006419
      ELSE IF (node == 19) THEN
         ! M019
         om_m = 0.1256
         om_b = 0.0222
         cosm%sig8 = 0.8547
         cosm%h = 0.8251
         cosm%ns = 1.0265
         cosm%w = -0.9813
         cosm%wa = -0.3393
         om_nu = 0.004673
      ELSE IF (node == 20) THEN
         ! M020
         om_m = 0.1514
         om_b = 0.0225
         cosm%sig8 = 0.7561
         cosm%h = 0.6827
         cosm%ns = 0.9913
         cosm%w = -1.0101
         cosm%wa = -0.7778
         om_nu = 0.009777
      ELSE IF (node == 21) THEN
         ! M021
         om_m = 0.1472
         om_b = 0.0221
         cosm%sig8 = 0.8475
         cosm%h = 0.6583
         cosm%ns = 0.9613
         cosm%w = -0.9111
         cosm%wa = -1.5470
         om_nu = 0.000672
      ELSE IF (node == 22) THEN
         ! M022
         om_m = 0.1384
         om_b = 0.0231
         cosm%sig8 = 0.8328
         cosm%h = 0.8234
         cosm%ns = 0.9739
         cosm%w = -0.9312
         cosm%wa = 0.5939
         om_nu = 0.008239
      ELSE IF (node == 23) THEN
         ! M023
         om_m = 0.1334
         om_b = 0.0225
         cosm%sig8 = 0.7113
         cosm%h = 0.7352
         cosm%ns = 0.9851
         cosm%w = -0.8971
         cosm%wa = 0.3247
         om_nu = 0.003733
      ELSE IF (node == 24) THEN
         ! M024
         om_m = 0.1508
         om_b = 0.0229
         cosm%sig8 = 0.7002
         cosm%h = 0.7935
         cosm%ns = 0.8685
         cosm%w = -1.0322
         cosm%wa = 1.0220
         om_nu = 0.003063
      ELSE IF (node == 25) THEN
         ! M025
         om_m = 0.1203
         om_b = 0.0230
         cosm%sig8 = 0.8773
         cosm%h = 0.6240
         cosm%ns = 0.9279
         cosm%w = -0.8282
         cosm%wa = -1.5005
         om_nu = 0.007024
      ELSE IF (node == 26) THEN
         ! M026
         om_m = 0.1224
         om_b = 0.0222
         cosm%sig8 = 0.7785
         cosm%h = 0.7377
         cosm%ns = 0.8618
         cosm%w = -0.7463
         cosm%wa = 0.3647
         om_nu = 0.002082
      ELSE IF (node == 27) THEN
         ! M027
         om_m = 0.1229
         om_b = 0.0234
         cosm%sig8 = 0.8976
         cosm%h = 0.8222
         cosm%ns = 0.9698
         cosm%w = -1.0853
         cosm%wa = 0.8683
         om_nu = 0.002902
      ELSE IF (node == 28) THEN
         ! M028
         om_m = 0.1229
         om_b = 0.0231
         cosm%sig8 = 0.8257
         cosm%h = 0.6109
         cosm%ns = 0.9885
         cosm%w = -0.9311
         cosm%wa = 0.8693
         om_nu = 0.009086
      ELSE IF (node == 29) THEN
         ! M029
         om_m = 0.1274
         om_b = 0.0228
         cosm%sig8 = 0.8999
         cosm%h = 0.8259
         cosm%ns = 0.8505
         cosm%w = -0.7805
         cosm%wa = 0.5688
         om_nu = 0.006588
      ELSE IF (node == 30) THEN
         ! M030
         om_m = 0.1404
         om_b = 0.0222
         cosm%sig8 = 0.8232
         cosm%h = 0.6852
         cosm%ns = 0.8679
         cosm%w = -0.8594
         cosm%wa = -0.4637
         om_nu = 0.008126
      ELSE IF (node == 31) THEN
         ! M031
         om_m = 0.1386
         om_b = 0.0229
         cosm%sig8 = 0.7693
         cosm%h = 0.6684
         cosm%ns = 1.0478
         cosm%w = -1.2670
         cosm%wa = 1.2536
         om_nu = 0.006502
      ELSE IF (node == 32) THEN
         ! M032
         om_m = 0.1369
         om_b = 0.021501 ! Moved off boundary
         cosm%sig8 = 0.8812
         cosm%h = 0.8019
         cosm%ns = 1.0005
         cosm%w = -0.7282
         cosm%wa = -1.6927
         om_nu = 0.000905
      ELSE IF (node == 33) THEN
         ! M033
         om_m = 0.1286
         om_b = 0.0230
         cosm%sig8 = 0.7005
         cosm%h = 0.6752
         cosm%ns = 1.0492
         cosm%w = -0.7119
         cosm%wa = -0.8184
         om_nu = 0.007968
      ELSE IF (node == 34) THEN
         ! M034
         om_m = 0.1354
         om_b = 0.0216
         cosm%sig8 = 0.7018
         cosm%h = 0.5970
         cosm%ns = 0.8791
         cosm%w = -0.8252
         cosm%wa = -1.1148
         om_nu = 0.003602
      ELSE IF (node == 35) THEN
         ! M035
         om_m = 0.1359
         om_b = 0.0228
         cosm%sig8 = 0.8210
         cosm%h = 0.6815
         cosm%ns = 0.9872
         cosm%w = -1.1642
         cosm%wa = -0.1801
         om_nu = 0.004440
      ELSE IF (node == 36) THEN
         ! M036
         om_m = 0.1390
         om_b = 0.0220
         cosm%sig8 = 0.8631
         cosm%h = 0.6477
         cosm%ns = 0.8985
         cosm%w = -0.8632
         cosm%wa = 0.8285
         om_nu = 0.001082
      ELSE
         STOP 'MIRA_TITAN_NODE_COSMOLOGY: Error, i specified incorrectly'
      END IF

      cosm%iw = iw_waCDM
      cosm%Om_m = om_m/cosm%h**2
      cosm%Om_b = om_b/cosm%h**2
      cosm%m_nu = neutrino_constant(cosm)*om_nu
      cosm%N_nu = 3
      cosm%Om_w = 1.-cosm%Om_m
      cosm%Om_v = 0.
      cosm%kpiv = 0.05/cosm%h

   END SUBROUTINE Mira_Titan_node_cosmology

   SUBROUTINE HALOFIT_init(rknl, rneff, rncur, a, cosm, verbose)

      ! Halofit initialisation routine taken from https://www.roe.ac.uk/~jap/haloes/
      ! Calculate the non-linear wavenumber (rknl), effective spectral index (rneff) and curvature (rncur) 
      ! of the power spectrum at the desired redshift, using the method described in Smith et al (2003).
      REAL, INTENT(OUT) :: rknl
      REAL, INTENT(OUT) :: rneff
      REAL, INTENT(OUT) :: rncur
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      REAL :: xlogr1, xlogr2, rmid, sig, d1, d2, diff
      REAL, PARAMETER :: diff_limit = HALOFIT_acc
      REAL, PARAMETER :: logr1_init = HALOFIT_logr1
      REAL, PARAMETER :: logr2_init = HALOFIT_logr2

      IF (verbose) WRITE (*, *) 'HALOFIT_INIT: computing effective spectral quantities:'

      xlogr1 = logr1_init
      xlogr2 = logr2_init
      DO
         rmid = 10**((xlogr2+xlogr1)/2.)
         CALL wint_HALOFIT(rmid, sig, d1, d2, a, cosm)
         diff = sig-1.0
         IF (abs(diff) <= diff_limit) THEN
            rknl = 1./rmid
            rneff = -3.-d1
            rncur = -d2
            EXIT
         ELSE IF (diff > diff_limit) THEN
            xlogr1 = log10(rmid)
         ELSE IF (diff < -diff_limit) THEN
            xlogr2 = log10(rmid)
         END IF
      END DO

      IF (verbose) THEN
         WRITE (*, *) 'HALOFIT_INIT: z =', redshift_a(a)
         WRITE (*, *) 'HALOFIT_INIT: rknl [h/Mpc] =', rknl
         WRITE (*, *) 'HALOFIT_INIT: rneff =', rneff
         WRITE (*, *) 'HALOFIT_INIT: rncur =', rncur
         WRITE (*, *) 'HALOFIT_INIT: initialised'
         WRITE (*, *)
      END IF

   END SUBROUTINE HALOFIT_init

   SUBROUTINE calculate_HALOFIT(k, a, Pk, cosm, version)

      ! Fill array P(k,a) with HALOFIT power spectrum for input arrays of k and a
      REAL, INTENT(IN) :: k(:) ! Array of wavenumbers
      REAL, INTENT(IN) :: a(:) ! Scale factor
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :) ! Output power array
      TYPE(cosmology), INTENT(INOUT) :: cosm   ! Cosmology
      INTEGER, OPTIONAL, INTENT(IN) :: version ! HALOFIT version
      REAL, ALLOCATABLE :: Pli(:), Pq(:), Ph(:), Pnl(:)
      INTEGER :: ia, na, nk
      LOGICAL, PARAMETER :: verbose = .FALSE.

      nk = size(k)
      na = size(a)
      ALLOCATE(Pli(nk), Pq(nk), Ph(nk), Pnl(nk))
      ALLOCATE(Pk(nk, na))

      DO ia = 1, na
         CALL calculate_HALOFIT_a(k, a(ia), Pli, Pq, Ph, Pnl, nk, cosm, verbose, version)
         Pk(:, ia) = Pnl
      END DO

   END SUBROUTINE calculate_HALOFIT

   SUBROUTINE calculate_HALOFIT_a(k, a, Pli, Pq, Ph, Pnl, n, cosm, verbose, ihf)

      ! Get a HALOFIT P(k,a) prediction for an array of k values
      INTEGER, INTENT(IN) :: n       ! Size of input k array
      REAL, INTENT(IN) :: k(n)       ! Array of wavenumbers
      REAL, INTENT(IN) :: a          ! Scale factor
      REAL, INTENT(OUT) :: Pli(n)   ! Output array of linear power
      REAL, INTENT(OUT) :: Pq(n)     ! Output array of quasi-linear power
      REAL, INTENT(OUT) :: Ph(n)     ! Output array of halo power
      REAL, INTENT(OUT) :: Pnl(n)    ! Output array of HALOFIT power
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmological model
      LOGICAL, INTENT(IN) :: verbose ! Verbosity
      INTEGER, INTENT(IN) :: ihf     ! Interger for HALOFIT version
      REAL :: rknl, rneff, rncur
      INTEGER :: i

      CALL HALOFIT_init(rknl, rneff, rncur, a, cosm, verbose)

      DO i = 1, n
         Pli(i) = plin(k(i), a, flag_matter, cosm)
         CALL calculate_HALOFIT_ka(k(i), rneff, rncur, rknl, Pli(i), Pnl(i), Pq(i), Ph(i), a, cosm, ihf)
      END DO

   END SUBROUTINE calculate_HALOFIT_a

   SUBROUTINE wint_HALOFIT(r, sig, d1, d2, a, cosm)

      ! HALOFIT window integration routine
      REAL, INTENT(IN) :: r
      REAL, INTENT(OUT) :: sig
      REAL, INTENT(OUT) :: d1
      REAL, INTENT(OUT) :: d2
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sum1, sum2, sum3, t, y, x, w1, w2, w3, rk
      INTEGER :: i
      INTEGER, PARAMETER :: nint = 3000

      sum1 = 0.d0
      sum2 = 0.d0
      sum3 = 0.d0
      DO i = 1, nint
         t = (float(i)-0.5)/float(nint)
         y = -1.+1./t
         rk = y
         d2 = plin(rk, a, flag_matter, cosm)
         x = y*r
         w1 = exp(-x*x)
         w2 = 2*x*x*w1
         w3 = 4*x*x*(1-x*x)*w1
         sum1 = sum1+w1*d2/y/t/t
         sum2 = sum2+w2*d2/y/t/t
         sum3 = sum3+w3*d2/y/t/t
      END DO
      sum1 = sum1/float(nint)
      sum2 = sum2/float(nint)
      sum3 = sum3/float(nint)
      sig = sqrt(sum1)
      d1 = -sum2/sum1
      d2 = -sum2*sum2/sum1/sum1-sum3/sum1

   END SUBROUTINE wint_HALOFIT

   SUBROUTINE calculate_HALOFIT_ka(rk, rn, rncur, rknl, pli, pnl, pq, ph, a, cosm, ihf)

      ! Calculates the HALOFIT power spectrum after rn, rncur and rknl have been pre-calculated
      ! Smith et al. (2003 ;astro-ph/0207664)
      ! Bird et al. (2012; arXiv:)
      ! Takahashi et al. (2012; arXiv:)
      ! Takahashi et al. (2012) but taken from CAMB with some neutrino stuff
      ! Takahashi et al. (2012) but taken from CLASS with some neutrino stuff (https://github.com/cmbant/CAMB/issues/44)
      REAL, INTENT(IN) :: rk
      REAL, INTENT(IN) :: rn
      REAL, INTENT(IN) :: rncur
      REAL, INTENT(IN) :: rknl
      REAL, INTENT(IN) :: pli
      REAL, INTENT(OUT) :: pnl
      REAL, INTENT(OUT) :: pq
      REAL, INTENT(OUT) :: ph
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, INTENT(IN) :: ihf
      REAL :: gam, aa, bb, cc, mu, nu, alpha, beta, f1, f2, f3, y, P, Q, fy
      REAL :: om_mz, om_vz, fnu, om_m, wz
      real :: f1a, f2a, f3a, f1b, f2b, f3b, frac

      ! Necessary cosmological parameters
      Om_m = cosm%Om_m        
      Om_mz = Omega_m(a, cosm)
      Om_vz = Omega_v(a, cosm)+Omega_w(a, cosm) ! Note this well
      wz = w_de(a, cosm) ! Choice here; do you use w or w(z)? w(z) is better I think; CAMB makes w(z) choice too
      fnu = cosm%Om_nu/cosm%Om_m

      IF (ihf == HALOFIT_Smith_paper) THEN
         ! Smith et al. (2003); the numbers here are EXACTLY those quoted in the paper!
         aa = 10**(1.4861+1.8369*rn+1.6762*rn**2+0.7940*rn**3+0.1670*rn**4-0.6206*rncur) ! Smith equation (C9)
         bb = 10**(0.9463+0.9466*rn+0.3084*rn**2-0.9400*rncur) ! Smith equation (C10)
         cc = 10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur) ! Smith equation (C11)
         gam = 0.8649+0.2989*rn+0.1631*rncur ! Smith equation (C12)
         alpha = 1.3884+0.3700*rn-0.1452*rn**2 ! Smith equation (C13)
         beta = 0.8291+0.9854*rn+0.3401*rn**2 ! Smith equation (C14)
         mu = 10**(-3.5442+0.1908*rn) ! Smith equation (C15)
         nu = 10**(0.9589+1.2857*rn) ! Smith equation (C16)
      ELSE IF (ihf == HALOFIT_Smith_code) THEN
         ! Smith et al. (2003); the numbers here are EXACTLY those from the online code and the same numbers as in CAMB
         aa = 10**(1.4861+1.83693*rn+1.67618*rn**2+0.7940*rn**3+0.1670756*rn**4-0.620695*rncur)
         bb = 10**(0.9463+0.9466*rn+0.3084*rn**2-0.940*rncur)
         cc = 10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur)
         gam = 0.86485+0.2989*rn+0.1631*rncur
         alpha = 1.38848+0.3701*rn-0.1452*rn**2
         beta = 0.8291+0.9854*rn+0.3400*rn**2
         mu = 10**(-3.54419+0.19086*rn)
         nu = 10**(0.95897+1.2857*rn)
      ELSE IF (ihf == HALOFIT_Bird_paper) THEN
         ! Bird et al. (2012); based off Smith et al. (2003) with numbers taken directly from the papers
         aa = 10**(1.4861+1.8369*rn+1.6762*rn**2+0.7940*rn**3+0.1670*rn**4-0.6206*rncur)
         bb = 10**(0.9463+0.9466*rn+0.3084*rn**2-0.9400*rncur)
         cc = 10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur)
         gam = 0.8649+0.2989*rn+0.1631*rncur+0.316-0.0765*rn-0.835*rncur ! Bird equation (A5), last 3 terms are new
         alpha = 1.3884+0.3700*rn-0.1452*rn**2
         beta = 0.8291+0.9854*rn+0.3401*rn**2+fnu*(-6.49+1.44*rn**2) ! Bird equation (A10), fnu term is new
         mu = 10**(-3.5442+0.1908*rn)
         nu = 10**(0.9589+1.2857*rn)
      ELSE IF (ihf == HALOFIT_Bird_code) THEN
         ! Bird et al. (2012); based off Smith et al. (2003) with numbers taken from the online code
         aa = 10**(1.4861+1.83693*rn+1.67618*rn**2+0.7940*rn**3+0.1670756*rn**4-0.620695*rncur)
         bb = 10**(0.9463+0.9466*rn+0.3084*rn**2-0.940*rncur)
         cc = 10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur)
         gam = 0.86485+0.2989*rn+0.1631*rncur+0.316-0.0765*rn-0.835*rncur ! Bird equation (A5), last 3 terms are new       
         alpha = 1.38848+0.3701*rn-0.1452*rn**2
         beta = 0.8291+0.9854*rn+0.3400*rn**2+fnu*(-6.49+1.44*rn**2) ! Bird equation (A10), fnu term is new  
         mu = 10**(-3.54419+0.19086*rn)
         nu = 10**(0.95897+1.2857*rn)
      ELSE IF (ihf == HALOFIT_Bird_CAMB) THEN
         ! Bird et al. (2012); based off Smith et al. (2003) with numbers taken from CAMB
         aa = 10**(1.4861+1.83693*rn+1.67618*rn**2+0.7940*rn**3+0.1670756*rn**4-0.620695*rncur)
         bb = 10**(0.9463+0.9466*rn+0.3084*rn**2-0.940*rncur)
         cc = 10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur)
         gam = 0.86485+0.2989*rn+0.1631*rncur+0.3159-0.0765*rn-0.8350*rncur ! Bird equation (A5), last 3 terms are new
         alpha = 1.38848+0.3701*rn-0.1452*rn**2
         beta = 0.8291+0.9854*rn+0.3400*rn**2+fnu*(-6.4868+1.4373*rn**2) ! Bird equation (A10), fnu term is new
         mu = 10**(-3.54419+0.19086*rn)
         nu = 10**(0.95897+1.2857*rn)
      ELSE IF (ihf == HALOFIT_Takahashi) THEN
         ! Takahashi et al. (2012); complete refit, all parameters different from Smith et al. (2003)
         aa = 10**(1.5222+2.8553*rn+2.3706*rn**2+0.9903*rn**3+0.2250*rn**4-0.6038*rncur+0.1749*Om_vz*(1.+wz)) ! Takahashi equation (A6)
         bb = 10**(-0.5642+0.5864*rn+0.5716*rn**2-1.5474*rncur+0.2279*Om_vz*(1.+wz)) ! Takahashi equation (A7)
         cc = 10**(0.3698+2.0404*rn+0.8161*rn**2+0.5869*rncur) ! Takahashi equation (A8)
         gam = 0.1971-0.0843*rn+0.8460*rncur ! Takahashi equation (A9)
         alpha = abs(6.0835+1.3373*rn-0.1959*rn**2-5.5274*rncur) ! Takahashi equation (A10; note the ABS)
         beta = 2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+0.3980*rn**4-0.1682*rncur ! Takahashi equation (A11)
         mu = 0. ! Takahashi equation (A12)
         nu = 10**(5.2105+3.6902*rn) ! Takahashi equation (A13)
      ELSE IF (ihf == HALOFIT_CAMB) THEN
         ! Unpublished CAMB from halofit_ppf.f90; based on Takahashi et al. (2012) plus Bird et al. (2012)
         aa = 10**(1.5222+2.8553*rn+2.3706*rn**2+0.9903*rn**3+0.2250*rn**4-0.6038*rncur+0.1749*Om_vz*(1.+wz)) ! Same as Takahashi
         bb = 10**(-0.5642+0.5864*rn+0.5716*rn**2-1.5474*rncur+0.2279*Om_vz*(1.+wz)) ! Same as Takahashi
         cc = 10**(0.3698+2.0404*rn+0.8161*rn**2+0.5869*rncur) ! Same as Takahashi
         gam = 0.1971-0.0843*rn+0.8460*rncur ! Same as Takahashi; no Bird modification
         alpha = abs(6.0835+1.3373*rn-0.1959*rn**2-5.5274*rncur) ! Same as Takahashi
         beta = 2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+0.3980*rn**4-0.1682*rncur+fnu*(1.081+0.395*rn**2) ! CAMB; halofit_ppf.f90
         mu = 0. ! Same as Takahashi
         nu = 10**(5.2105+3.6902*rn) ! Same as Takahashi
      ELSE
         STOP 'HALOFIT: Error, ihf specified incorrectly'
      END IF

      ! Here we need Omega at redshift of interest and calculate the Omega evolution
      ! As far as I know, no HALOFIT extension has modified these numbers
      IF (abs(1.-Om_mz) > 0.01) THEN

         ! Open model
         f1a = Om_mz**(-0.0732)
         f2a = Om_mz**(-0.1423)
         f3a = Om_mz**(0.0725)

         ! Flat LCDM
         f1b = Om_mz**(-0.0307)
         f2b = Om_mz**(-0.0585)
         f3b = Om_mz**(0.0743)

         ! Linearly interpolate between LCDM and open case for mixed model
         frac = Om_vz/(1.-Om_mz)
         f1 = frac*f1b+(1.-frac)*f1a
         f2 = frac*f2b+(1.-frac)*f2a
         f3 = frac*f3b+(1.-frac)*f3a

      ELSE

         ! Set to one for EdS
         f1 = 1.
         f2 = 1.
         f3 = 1.

      END IF

      ! Ratio of current wave number to the non-linear wave number
      y = rk/rknl

      ! Quasi-linear term    
      IF (ihf == HALOFIT_Bird_code .OR. ihf == HALOFIT_Bird_paper .OR. ihf == HALOFIT_Bird_CAMB) THEN
         ! Bird et al. (2012)
         P = 26.3*fnu*rk**2/(1.+1.5*rk**2)
      ELSE IF (ihf == HALOFIT_CAMB) THEN
         ! Unpublished CAMB stuff from halofit.f90; note that this is used by ALL HALOFIT versions in CAMB, including Smith and Bird versions
         P = 47.48*fnu*rk**2/(1.+1.5*rk**2)
      ELSE
         P = 0.
      END IF
      pq = pli*(1.+P)
      fy = y/4.+y**2/8. ! Smith (below C2)
      pq = pli*(1.+pq)**beta/(1.+pq*alpha)*exp(-fy)

      ! Halo term
      ph = aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam)) ! Smith (C4); Bird (A2); Takahashi (A3ii)
      ph = ph/(1.+mu*y**(-1)+nu*y**(-2))                   ! Smith (C3); Bird (A1); Takahashi (A3i)
      IF (ihf == HALOFIT_Bird_code .OR. ihf == HALOFIT_Bird_paper .OR. ihf == HALOFIT_Bird_CAMB) THEN
         Q = fnu*(2.080-12.4*(Om_m-0.3))/(1.+1.2e-3*y**3) ! Bird equation (A6); note Omega_m term
      ELSE IF (ihf == HALOFIT_CAMB) THEN
         ! Unpublished CAMB stuff from halofit.f90; note that this is used by ALL HALOFIT versions in CAMB, including Smith and Bird versions   
         Q = fnu*0.977 ! CAMB; halofit_ppf.f90; halofit; note no Omega_m term
      ELSE
         Q = 0.
      END IF
      ph = ph*(1.+Q)

      ! Sum the quasi-linear and halo terms to get the total
      pnl = pq+ph

   END SUBROUTINE calculate_HALOFIT_ka

   REAL FUNCTION z_CMB(cosm)

      ! Fit to the redshift of last scattering, used in the original Cosmic Emu for h determination
      ! Equation (23) from Hu & White (1997) https://arxiv.org/pdf/astro-ph/9609079.pdf
      ! TODO: What does this assume about dark energy?
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: b1, b2
      REAL :: om_b, om_m

      om_m = cosm%Om_m*cosm%h**2
      om_b = cosm%Om_b*cosm%h**2

      b1 = (0.0783*om_b**(-0.238))/(1.+39.5*om_b**0.763)
      b2 = 0.560/(1.+21.1*om_b**1.81)

      z_CMB = 1048.*(1.+0.00124*om_b**(-0.738))*(1.+b1*om_m**b2)

   END FUNCTION z_CMB

   REAL FUNCTION r_sound(cosm)

      ! Sound horizon [Mpc/h] calculated in WKB approximation, used in the original Cosmic Emu for h determination
      ! Equation (B6) from Hu & Sugiyama (1995) https://arxiv.org/pdf/astro-ph/9407093.pdf
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: a_eq, R, R_eq

      ! Scale factor at matter-radiation equality
      a_eq = scale_factor_eq(cosm)

      ! Scale factors normalised to 0.75 at baryon-photon equality
      R = 0.75*(1./(cosm%Om_g/cosm%Om_b))
      R_eq = 0.75*(a_eq/(cosm%Om_g/cosm%Om_b))

      ! Sound horizon approximation
      r_sound = (2./3.)*sqrt(6./R_eq)*log((sqrt(1.+R)+sqrt(R+R_eq))/(1.+sqrt(R_eq)))
      r_sound = r_sound/k_eq(cosm)

      STOP 'R_SOUND: Check this, I think it is wrong'

   END FUNCTION r_sound

   REAL FUNCTION scale_factor_eq(cosm)

      ! Scale factor at matter-radiation equality
      TYPE(cosmology), INTENT(IN) :: cosm

      scale_factor_eq = cosm%Om_r/cosm%Om_m

   END FUNCTION scale_factor_eq

   REAL FUNCTION k_eq(cosm)

      ! Wavenumber at matter-radiation equality [h/Mpc]
      ! TODO: How exactly is this defined?
      ! TODO: Why is it 4-2*sqrt(2), rather than 1, pi or 2pi?
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: a_eq

      a_eq = scale_factor_eq(cosm)
      k_eq = (4.-2.*sqrt(2.))/comoving_particle_horizon(a_eq, cosm)

      STOP 'K_EQ: Check this'

   END FUNCTION k_eq

   ELEMENTAL REAL FUNCTION Pk_Delta(Delta, k)

      ! Converts dimensionless Delta^2(k) to P(k) [Mpc/h]^3
      REAL, INTENT(IN) :: Delta ! Power spectrum in Delta^2(k) dimensionless form
      REAL, INTENT(IN) :: k     ! Wavenumber [h/Mpc]

      Pk_Delta = Delta/(k/twopi)**3
      Pk_Delta = Pk_Delta/(4.*pi)

   END FUNCTION Pk_Delta

   ELEMENTAL REAL FUNCTION Delta_Pk(Pk, k)

      ! Converts P(k) [Mpc/h]^3 to dimensionless Delta^2(k) 
      REAL, INTENT(IN) :: Pk ! Power spectrum in P(k) [Mpc/h]^3
      REAL, INTENT(IN) :: k  ! Wavenumber [h/Mpc]

      Delta_Pk = (4.*pi)*((k/twopi)**3)*Pk

   END FUNCTION Delta_Pk

   SUBROUTINE init_wiggle(cosm)

      ! Isolate the power spectrum wiggle
      USE special_functions
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: k(:), a(:), Pk(:, :)
      REAL, ALLOCATABLE :: Pk_smooth(:, :), Pk_wiggle(:, :)
      INTEGER  :: na
      REAL, PARAMETER :: amin = amin_plin
      REAL, PARAMETER :: amax = amax_plin
      REAL, PARAMETER :: kmin = kmin_wiggle
      REAL, PARAMETER :: kmax = kmax_wiggle
      INTEGER, PARAMETER  :: nk = nk_wiggle

      ! Words
      IF (cosm%verbose) WRITE(*, *) 'INIT_WIGGLE: Starting'

      ! Checks
      IF (cosm%box) STOP 'INIT_WIGGLE: Error, cannot extract wiggle from truncated linear power'
      IF (cosm%scale_dependent_growth .AND. (.NOT. cosm%has_power)) THEN
         STOP 'INIT_WIGGLE: This will not work for scale-dependent growth unless you have tabulated power'
      END IF

      ! Allocate array for k
      CALL fill_array_log(kmin, kmax, k, nk)
      IF (scale_grow_wiggle .AND. cosm%scale_dependent_growth) THEN
         na = na_plin
         CALL fill_array_log(amin, amax, a, na)
      ELSE
         na = 1
         a = [1.]      
      END IF

      ! Get the linear power spectrum
      IF(.NOT. cosm%is_normalised) STOP 'INIT_WIGGLE: Error, linear power must be normalised'
      CALL calculate_Plin(k, a, Pk, flag_matter, cosm)

      ! Write details to screen
      IF (cosm%verbose) THEN
         WRITE(*, *) 'INIT_WIGGLE: kmin [h/Mpc]:', k(1)
         WRITE(*, *) 'INIT_WIGGLE: kmax [h/Mpc]:', k(nk)
         WRITE(*, *) 'INIT_WIGGLE: nk:', nk
         WRITE(*, *) 'INIT_WIGGLE: Splitting into wiggle and broad-band'
      END IF

      ! Calculate the smooth power spectrum
      IF (cosm%verbose) WRITE(*, *) 'INIT_WIGGLE: Calculating smooth power spectrum'
      CALL smooth_power(k, a, Pk, Pk_smooth, cosm)

      ! Isolate the wiggle
      IF (cosm%verbose) WRITE(*, *) 'INIT_WIGGLE: Isolating wiggle'   
      Pk_wiggle = Pk-Pk_smooth

      ! Init interpolator
      IF (cosm%verbose) WRITE(*, *) 'INIT_WIGGLE: Initialising interpolator'
      CALL init_interpolator(k, a, Pk_wiggle, cosm%wiggle, &
         iorder = iorder_interp_wiggle, &
         iextrap = iextrap_wiggle, &
         store = store_wiggle, &
         logx = .TRUE., &
         logy = .TRUE., &
         logf = .FALSE.)

      ! Set the flag
      cosm%has_wiggle = .TRUE.

      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_WIGGLE: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_wiggle

   REAL FUNCTION P_wig(k, a, cosm)

      ! Calculate the wiggle part of the power spectrum
      ! Defined like P_wiggle = P_lin - P_broadband
      ! Note that this scales with the growth factor as the right-hand side does
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: g

      IF (.NOT. cosm%has_wiggle) CALL init_wiggle(cosm)
      P_wig = evaluate_interpolator(k, a, cosm%wiggle)
      IF (scale_grow_wiggle .AND. cosm%scale_dependent_growth) THEN       
         g = 1. ! TODO: This should be if cosm%wiggle is really 1D and a = 1 only
      ELSE
         g = grow(a, cosm) 
      END IF
      P_wig = P_wig*g**2

   END FUNCTION P_wig

   REAL FUNCTION P_smt(k, a, cosm)

      ! The smooth 'broad band' power spectrum, which is the linear minus the entire wiggle
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: flag = flag_matter

      P_smt = Plin(k, a, flag, cosm)-P_wig(k, a, cosm)

   END FUNCTION P_smt

   REAL FUNCTION P_dwg(k, a, sigv, cosm)

      ! Calculate the dewiggled power spectrum, which is linear but with damped wiggles
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: sigv
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: f, P_linear
      INTEGER, PARAMETER :: flag = flag_matter

      IF (.NOT. cosm%is_normalised) CALL normalise_power(cosm) 
      P_linear = plin(k, a, flag, cosm) ! Needed here to make sure it is init before init_wiggle   
      f = exp(-(k*sigv)**2)
      P_dwg = P_linear+(f-1.)*P_wig(k, a, cosm)

   END FUNCTION P_dwg

   ! SUBROUTINE calculate_P_wig(k, a, Pk, cosm)

   !    REAL, INTENT(IN) :: k(:)
   !    REAL, INTENT(IN) :: a(:)
   !    REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
   !    TYPE(cosmology), INTENT(INOUT) :: cosm
   !    INTEGER :: ik, ia
   !    INTEGER :: nk, na

   !    nk = size(k)
   !    na = size(a)
   !    ALLOCATE(Pk(nk, na))
   !    DO ia = 1, na
   !       DO ik = 1, nk
   !          Pk(ik, ia) = P_wig(k(ik), a(ia), cosm)
   !       END DO
   !    END DO

   ! END SUBROUTINE calculate_P_wig

   SUBROUTINE calculate_P_smt(k, a, Pk, cosm)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: ik, ia
      INTEGER :: nk, na

      nk = size(k)
      na = size(a)
      ALLOCATE(Pk(nk, na))
      DO ia = 1, na
         DO ik = 1, nk
            Pk(ik, ia) = P_smt(k(ik), a(ia), cosm)
         END DO
      END DO

   END SUBROUTINE calculate_P_smt

   SUBROUTINE calculate_P_dwg(k, a, Pk, cosm)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: ik, ia
      INTEGER :: nk, na
      REAL :: sigv
      INTEGER, PARAMETER :: flag = flag_matter

      nk = size(k)
      na = size(a)
      ALLOCATE(Pk(nk, na))
      DO ia = 1, na
         sigv = sigmaV(0., a(ia), flag, cosm)
         DO ik = 1, nk
            Pk(ik, ia) = P_dwg(k(ik), a(ia), sigv, cosm)
         END DO
      END DO

   END SUBROUTINE calculate_P_dwg

   REAL FUNCTION P_nw(k, a, cosm)

      ! Calculates the un-normalised no-wiggle power spectrum 
      ! Comes from the Eisenstein & Hu approximation
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      P_nw = ((cosm%A*grow(a, cosm))**2)*primordial_spectrum(k, cosm)*Tk_nw(k, cosm)**2

   END FUNCTION P_nw

   SUBROUTINE calculate_nowiggle(k, a, Pk_nw, cosm)

      ! Calculate the normalised no wiggle power spectrum at a range of k and a
      ! Comes from the Eisenstein & Hu approximation
      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk_nw(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: ik, ia, nk, na
      REAL :: Pk_lin, Pk_nw_norm
      REAL, PARAMETER :: knorm = knorm_nowiggle
      INTEGER, PARAMETER :: iorder = 3 ! TODO: Move to header
      INTEGER, PARAMETER :: ifind = ifind_split ! TODO: Move to header
      INTEGER, PARAMETER :: iinterp = iinterp_Lagrange ! TODO: Move to header

      ! Allocate arrays
      nk = size(k)
      na = size(a)
      ALLOCATE(Pk_nw(nk, na))

      ! Get the Eisenstein & Hu no-wiggle power spectrum
      DO ia = 1, na
         DO ik = 1, nk
            Pk_nw(ik, :) = P_nw(k(ik), a(ia), cosm)
         END DO
      END DO

      ! Force spectra to agree at the minimum wavenumber
      DO ia = 1, na
         Pk_lin = Plin(knorm, a(ia), flag_matter, cosm)
         Pk_nw_norm = find(knorm, k, Pk_nw(:, ia), nk, iorder, ifind, iinterp)
         Pk_nw(:, ia) = Pk_nw(:, ia)*Pk_lin/Pk_nw_norm
      END DO

   END SUBROUTINE calculate_nowiggle

   SUBROUTINE smooth_power(k, a, Pk, Pk_smt, cosm)

      ! Calculate the normalised smoothed power spectrum at a range of k
      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: a(:)
      REAL, INTENT(IN) :: Pk(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk_smt(:, :)
      REAL, ALLOCATABLE :: Pk_nwg(:, :)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: ia, na
      REAL, PARAMETER :: dx = wiggle_dx
      REAL, PARAMETER :: sig = wiggle_sigma
      LOGICAL, PARAMETER :: divide = divide_by_nowiggle

      ! Reduce dynamic range
      IF (divide) THEN
         CALL calculate_nowiggle(k, a, Pk_nwg, cosm)
         Pk_smt = Pk/Pk_nwg
      ELSE   
         Pk_smt = log(Pk)
      END IF

      ! Smooth linear power
      na = size(a)
      DO ia = 1, na
         IF (wiggle_smooth == dewiggle_Gaussian) THEN
            CALL smooth_array_Gaussian(log(k), Pk_smt(:, ia), sig)
         ELSE IF (wiggle_smooth == dewiggle_tophat) THEN
            CALL smooth_array_tophat(log(k), Pk_smt(:, ia), dx)
         ELSE
            STOP 'SMOOTH_POWER: Error, smoothing method not recognised'
         END IF
      END DO

      ! Return dynamic range
      IF (divide) THEN
         Pk_smt = Pk_smt*Pk_nwg
      ELSE
         Pk_smt = exp(Pk_smt)
      END IF

   END SUBROUTINE smooth_power

   SUBROUTINE init_SPT(cosm)

      ! Initialise interpolator for SPT
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: k(:), Rk(:)
      INTEGER :: ik
      REAL, PARAMETER :: kmin = kmin_interpolate_SPT
      REAL, PARAMETER :: kmax = kmax_interpolate_SPT
      INTEGER, PARAMETER :: nk = nk_interpolate_SPT
      REAL, PARAMETER :: a = 1.
      INTEGER, PARAMETER :: flag = flag_matter
      INTEGER, PARAMETER :: iorder = iorder_interp_SPT
      INTEGER, PARAMETER :: iextrap = iextrap_SPT
      LOGICAL, PARAMETER :: store = store_interp_SPT

      IF(cosm%verbose) THEN
         WRITE(*, *) 'INIT_SPT: Initialising interpolator for SPT'
         WRITE(*, *) 'INIT_SPT: Calculating ratio with linear power at a =', a
         WRITE(*, *) 'INIT_SPT: Minimum wavenumber [h/Mpc]:', kmin
         WRITE(*, *) 'INIT_SPT: Minimum wavenumber [h/Mpc]:', kmax
         WRITE(*, *) 'INIT_SPT: Number of k points:', nk
      END IF

      ! Fill wavenumber array
      CALL fill_array_log(kmin, kmax, k, nk)

      ! Calculate SPT power and create ratio with respect to linear power
      ALLOCATE(Rk(nk))
      DO ik = 1, nk
         Rk(ik) = P_SPT_sum(k(ik), a, flag, cosm)/plin(k(ik), a, flag, cosm)
      END DO

      ! Initialise the interpolator
      CALL init_interpolator(k, Rk, cosm%rspt, iorder, iextrap, store, logx=.TRUE., logf=.FALSE.)
      cosm%has_SPT = .TRUE.

      IF(cosm%verbose) THEN
         WRITE(*, *) 'INIT_SPT: Done'
         WRITE(*, *)
      END IF

   END SUBROUTINE init_SPT

   RECURSIVE REAL FUNCTION P_SPT(k, a, cosm)

      ! Returns the one-loop SPT power spectrum in Delta^2(k) form
      ! Scales exactly as g(a)^4
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: kmin = kmin_integrate_SPT
      REAL, PARAMETER :: kmax = kmax_integrate_SPT
      INTEGER, PARAMETER :: flag = flag_matter
      LOGICAL, PARAMETER :: interp = interp_SPT

      ! Fill interpolator if necessary
      IF (interp_SPT .AND. (.NOT. cosm%has_SPT)) CALL init_SPT(cosm)

      IF (cosm%has_SPT) THEN
         IF (k < cosm%rspt%xmin) THEN
            P_SPT = 0.
         ELSE
            P_SPT = (grow(a, cosm)**2)*plin(k, a, flag, cosm)*evaluate_interpolator(k, cosm%rspt)
         END IF
      ELSE
         IF (k < kmin) THEN
            P_SPT = 0.
         ELSE IF (k > kmax) THEN
            STOP 'P_SPT: You are trying to evalulate SPT for too large a wavenumber'
         ELSE
            P_SPT = P_SPT_sum(k, a, flag, cosm)
         END IF
      END IF

   END FUNCTION P_SPT

   REAL FUNCTION P_SPT_sum(k, a, flag, cosm)

      ! One-loop power spectrum from summing P_13 and P_22; Delta^2(k)
      ! Scales exactly as g(a)^4
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/compute_pkd.f90 
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm

      ! Note factor of 2 in front of P_13
      P_SPT_sum = P_22_SPT(k, a, flag, cosm)+2.*P_13_SPT(k, a, flag, cosm)

   END FUNCTION P_SPT_sum

   REAL FUNCTION P_SPT_approx(k, a, cosm)

      ! An approximation to the SPT result, valid for k < 0.15 h/Mpc ish
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: ratio
      REAL, PARAMETER :: A1 = -0.0374
      REAL, PARAMETER :: A2 = -0.0054
      REAL, PARAMETER :: A3 = 0.0643
      REAL, PARAMETER :: kn = 0.1
      REAL, PARAMETER :: sig8 = 0.8
      INTEGER, PARAMETER :: flag = flag_matter

      !ratio = A1*(k/k1)+A3*(k/k3)**3
      ratio = A1*(k/kn)+A2*(k/kn)**2+A3*(k/kn)**3
      ratio = ratio*(grow(a, cosm)*cosm%sig8/sig8)**2
      P_SPT_approx = ratio*plin(k, a, flag, cosm)

   END FUNCTION P_SPT_approx

   REAL FUNCTION SPT_integrand(q, k, a, flag, cosm)

      REAL, INTENT(IN) :: q
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: I_22, I_13
      REAL :: lnq, epsi, epsk, m
      REAL, PARAMETER :: eps_multiple = eps_multiple_SPT

      epsk = k*eps_multiple

      ! This fixes the upper limit on the mu integral within I_22
      IF (q > k-epsk .AND. q < k+epsk) THEN
         epsi = epsk
      ELSE
         epsi = 0.
      END IF

      ! This fixes the multiple for the different parts of I_22
      ! See equation (22) of https://arxiv.org/pdf/astro-ph/9311070.pdf
      IF (q < epsk) THEN
         m = 2.
      ELSE
         m = 1.
      END IF

      ! Evaluate the integrands I_22 and I_13
      lnq = log(q)
      I_22 = 2.*m*Delta_Pk(P_22_integrand(lnq, k, a, epsi, flag, cosm)/twopi**2, k)
      I_13 = 0.5*Plin(k, a, flag, cosm)*P_13_integrand(lnq, k, a, flag, cosm)/twopi**2

      ! Sum to get the total
      SPT_integrand = I_22+I_13

   END FUNCTION SPT_integrand

   REAL FUNCTION P_22_SPT(k, a, flag, cosm)

      ! The P_22 contribution to the one-loop power from SPT
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/pkd.f90
      ! Routine called PkD in Komatsu code
      ! See https://arxiv.org/pdf/astro-ph/9311070.pdf for details
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: epsk, lnkmin, lnkmax, epsi, f(4)
      INTEGER :: i
      REAL, PARAMETER :: kmin = kmin_integrate_SPT
      REAL, PARAMETER :: kmax = kmax_integrate_SPT
      REAL, PARAMETER :: eps_multiple = eps_multiple_SPT
      REAL, PARAMETER :: acc = acc_integrate_SPT
      INTEGER, PARAMETER :: iorder = iorder_integrate_SPT

      ! A wavenumber close to k
      epsk = k*eps_multiple
      
      ! Loop over the 4 segments of the P_22 integral
      DO i = 1, 4

         ! Set boundary conditions for each segment
         IF (i == 1) THEN
            lnkmin = log(kmin)
            lnkmax = log(epsk)
            epsi = 0.
         ELSE IF (i == 2) THEN
            lnkmin = lnkmax
            lnkmax = log(k-epsk)
            epsi = 0.
         ELSE IF (i == 3) THEN
            lnkmin = lnkmax
            lnkmax = log(k+epsk)
            epsi = epsk
         ELSE IF (i == 4) THEN
            lnkmin = lnkmax
            lnkmax = log(kmax)
            epsi = 0.
         ELSE
            STOP 'P_22_SPT: Something went very wrong'
         END IF

         ! Do the integration (note that this is a nested integral)
         f(i) = integrate_cosm(lnkmin, lnkmax, P_22_integrand, k, a, epsi, flag, cosm, acc, iorder)

      END DO

      ! Sum contributions to create final result
      ! Extra factor of 2 before f(1) comes from excising pole (Jain & Bertschinger 1993; https://arxiv.org/pdf/astro-ph/9311070.pdf)
      P_22_SPT = 2.*(2.*f(1)+f(2)+f(3)+f(4))/twopi**2 ! Line 71 in Komatsu code
      P_22_SPT = Delta_Pk(P_22_SPT, k)                ! Convert from P(k) to Delta^2(k)

   END FUNCTION P_22_SPT

   REAL FUNCTION P_22_integrand(lnq, k, a, eps, flag, cosm)

      ! Integrand for the P_22 term in standard perturbation theory
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/pkd.f90
      ! Function called dp22dq_d in Komatsu code
      ! TODO: Clean up fudge for mu integration limits
      REAL, INTENT(IN) :: lnq
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: eps
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: q, mumax
      REAL, PARAMETER :: mumin = -1.
      REAL, PARAMETER :: acc = acc_integrate_SPT
      INTEGER, PARAMETER :: iorder = iorder_integrate_SPT

      ! Convert from log(q) integration variable
      q = exp(lnq)

      ! Set the maximum integration limit on mu (bit of a fudge)
      IF (eps == 0.) THEN
         mumax = 1.
      ELSE
         mumax = (k**2+q**2-eps**2)/(2.*k*q)
      END IF

      ! Do the integration and correct by factors
      ! In Komatsu code: dp22dq_d = q**3.*linear_pk(q)*ss but he uses P(k) whereas I use Delta^2(k) hence the 2pi^2 below
      P_22_integrand = integrate_cosm(mumin, mumax, P_22_inner_integrand, q, k, a, flag, cosm, acc, iorder)
      P_22_integrand = P_22_integrand*plin(q, a, flag, cosm)*2.*pi**2 

   END FUNCTION P_22_integrand

   REAL FUNCTION P_22_inner_integrand(mu, q, k, a, flag, cosm)

      ! Inner integral over mu that appears in P_22 from SPT
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/pkd.f90
      ! Function called dFFdmu in Komatsu code
      REAL, INTENT(IN) :: mu
      REAL, INTENT(IN) :: q
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: p

      p = sqrt(k**2+q**2-2.*k*q*mu)
      P_22_inner_integrand = Pk_Delta(plin(p, a, flag, cosm), p)*F2_SPT(k, q, mu)**2

   END FUNCTION P_22_inner_integrand

   REAL FUNCTION P_13_SPT(k, a, flag, cosm)

      ! The P_13 contribution to the non-linear power from SPT
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/pkd.f90     
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm 
      REAL, PARAMETER :: kmin = kmin_integrate_SPT
      REAL, PARAMETER :: kmax = kmax_integrate_SPT
      REAL, PARAMETER :: acc = acc_integrate_SPT
      INTEGER, PARAMETER :: iorder = iorder_integrate_SPT

      P_13_SPT = integrate_cosm(log(kmin), log(kmax), P_13_integrand, k, a, flag, cosm, acc, iorder)
      P_13_SPT = 0.5*P_13_SPT*Plin(k, a, flag, cosm)/twopi**2

   END FUNCTION P_13_SPT

   REAL FUNCTION P_13_integrand(lnq, k, a, flag, cosm)

      ! Integrand for the P_13 term in standard perturbation theory
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/pkd.f90
      ! Called dp13dq_d in Komatsu code
      REAL, INTENT(IN) :: lnq
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: q

      q = exp(lnq)
      P_13_integrand = q*F3_SPT(k, q)*Pk_Delta(plin(q, a, flag, cosm), q)

   END FUNCTION P_13_integrand

   REAL FUNCTION F2_SPT(k, q, mu)

      ! F2 kernel from standard perturbation theory
      ! vec{p} = vec{k} - vec{q}
      ! F2(vec{q},vec{k}-vec{q}) in Eq.(52) of Smith, Scoccimarro & Sheth (2006)
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/pkd.f90
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: q
      REAL, INTENT(IN) :: mu
      REAL :: p, kp, kq, pq

      kq = k*q*mu
      kp = k**2.-kq
      pq = kq-q**2 
      p  = sqrt(k**2+q**2-2*kq)
      F2_SPT = 5./7. &
         + (1./2.)*pq/(p*q)*((p/q)+(q/p)) &
         + (2./7.)*pq**2/(p*q)**2

   END FUNCTION F2_SPT

   REAL FUNCTION F3_SPT(k, q)

      ! F3 kernel from standard perturbation theory
      ! Taken from https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/powerspectrum/density3pt/pkd/pkd.f90
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: q
      REAL, PARAMETER :: q_on_k_limit = q_on_k_limit_F3

      ! Ref: Eq.[2.25] of Makino, Sasaki and Suto, PRD, 46, 585 (1992)
      IF (q/k < q_on_k_limit) THEN
         F3_SPT = k**2./252.*( &
            12.*(k/q)**2-158.+100.*(q/k)**2-42*(q/k)**4 &
            +(3./k**5/q**3)*(q**2-k**2)**3*(7.*q**2+2.*k**2) &
            *log((k+q)/abs(k-q)))
      ELSE
         ! Asymptotic value for q/k->infinity, accurate to 2e-5 at q/k=100
         F3_SPT = (-122./315.)*k**2
      END IF

   END FUNCTION F3_SPT

   SUBROUTINE calculate_rescaling_parameters(x1_tgt, x2_tgt, as_ini, s, a, a_tgt, cosm_ini, cosm_tgt, irescale, verbose)

      ! Calculates the AW10 (s, a) rescaling coefficients to go from cosm_ini to cosm_tgt
      ! cosm_tgt will be approximated by cosm_ini with R -> R/s (or k -> s*k) at scale factor a
      ! Note *NOT* a_tgt, which is the desired scale factor in the target cosmology
      REAL, INTENT(IN) :: x1_tgt                 ! Minimum value for minimization of (s, a) pair (either R or k)
      REAL, INTENT(IN) :: x2_tgt                 ! Maximum value for minimization of (s, a) pair (either R or k)
      REAL, INTENT(IN) :: as_ini(:)              ! Array of possible a values from the initial cosmology
      REAL, INTENT(OUT) :: s                     ! AW10 physical rescaling factor (R -> R/s)
      REAL, INTENT(OUT) :: a                     ! AW10 rescaling scale factor
      REAL, INTENT(IN) :: a_tgt                  ! Desired scale factor in the target cosmology
      TYPE(cosmology), INTENT(INOUT) :: cosm_ini ! Initial cosmmology
      TYPE(cosmology), INTENT(INOUT) :: cosm_tgt ! Target cosmology
      INTEGER, INTENT(IN) :: irescale            ! Rescaling type (either sigma(R) or P(k))
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose   ! Verbosity
      REAL, ALLOCATABLE :: s_tab(:)
      INTEGER :: ia, is, ia_best, is_best, na
      REAL :: a_best, s_best, cost, cost_best
      REAL :: k1_tgt, k2_tgt, R1_tgt, R2_tgt
      REAL, PARAMETER :: smin = smin_rescale
      REAL, PARAMETER :: smax = smax_rescale
      INTEGER, PARAMETER :: ns = ns_rescale
      !REAL, PARAMETER :: amin = amin_rescale
      !REAL, PARAMETER :: amax = amax_rescale
      !INTEGER, PARAMETER :: na = na_rescale

      ! Total number of a values
      na = size(as_ini)

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Starting'
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Minimum s:', smin
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Maximum s:', smax
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Number of s values:', ns
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Minimum a:', as_ini(1)
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Maximum a:', as_ini(na)
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Number of s values:', na
         IF (irescale == irescale_sigma) THEN
            WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Target R1 [Mpc/h]:', x1_tgt
            WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Target R2 [Mpc/h]:', x2_tgt
         ELSE IF (irescale == irescale_power) THEN
            WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Target k1 [h/Mpc]:', x1_tgt
            WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Target k2 [h/Mpc]:', x2_tgt
         ELSE
            STOP 'CALCULATE_RESCALING_PARAMETERS: Error, irescale not specified correctly'
         END IF
      END IF

      ! Set the range of either R or k over which to minimize
      IF (irescale == irescale_sigma) THEN
         R1_tgt = x1_tgt
         R2_tgt = x2_tgt
      ELSE IF (irescale == irescale_power) THEN
         k1_tgt = x1_tgt
         k2_tgt = x2_tgt
      ELSE
         STOP 'CALCULATE_RESCALING_PARAMETERS: Error, irescale not specified correctly'
      END IF

      ! Fill arrays of s, a values
      CALL fill_array(smin, smax, s_tab, ns)

      ! Loop over grid of values to find minimum cost
      cost_best = huge(cost_best)
      DO is = 1, ns
         DO ia = 1, na
      
            ! Calculate the cost associated with s, a values
            ! Note well that the integration variable is in the target cosmology
            s = s_tab(is)
            a = as_ini(ia)
            IF (irescale == irescale_sigma) THEN     
               cost = rescaling_cost_sigma(R1_tgt, R2_tgt, s, a, a_tgt, cosm_ini, cosm_tgt)
            ELSE IF (irescale == irescale_power) THEN
               cost = rescaling_cost_power(k1_tgt, k2_tgt, s, a, a_tgt, cosm_ini, cosm_tgt)
            ELSE
               STOP 'CALCULATE_RESCALING_PARAMETERS: Error, irescale not specified correctly'
            END IF

            ! If the cost is minimum then remember it
            IF (cost < cost_best) THEN
               s_best = s
               a_best = a
               is_best = is
               ia_best = ia
               cost_best = cost
            END IF

         END DO
      END DO

      ! Fix to minimum-cost values
      s = s_best
      a = a_best

      ! Write to screen and report error if best-fit is on the boundary
      IF (present_and_correct(verbose)) THEN
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Best-fitting s:', s
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Best-fitting a:', a
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Best-fitting z:', redshift_a(a)
      END IF
      IF (is_best == 1 .OR. is_best == ns) STOP 'CALCULATE_RESCALING_PARAMETERS: Error, best-fitting s is on boundary'
      IF (present_and_correct(verbose)) THEN
         WRITE(*, *) 'CALCULATE_RESCALING_PARAMETERS: Done'
         WRITE(*, *)
      END IF
   
   END SUBROUTINE calculate_rescaling_parameters

   REAL FUNCTION rescaling_cost_sigma(R1_tgt, R2_tgt, s, a_ini, a_tgt, cosm_ini, cosm_tgt)

      ! Rescaling cost function for minimizing sigma(R, a) differences
      ! Note that this integrates over a ln(R) range
      REAL, INTENT(IN) :: R1_tgt
      REAL, INTENT(IN) :: R2_tgt
      REAL, INTENT(IN) :: s
      REAL, INTENT(IN) :: a_ini
      REAL, INTENT(IN) :: a_tgt
      TYPE(cosmology), INTENT(INOUT) :: cosm_ini
      TYPE(cosmology), INTENT(INOUT) :: cosm_tgt
      INTEGER, PARAMETER :: flag = flag_rescale
      REAL, PARAMETER :: acc = acc_rescale
      INTEGER, PARAMETER :: iorder = iorder_rescale_integral

      ! Integrate and normalise 
      rescaling_cost_sigma = integrate_cosm(log(R1_tgt), log(R2_tgt), rescaling_cost_sigma_integrand, &
         s, a_ini, a_tgt, flag, cosm_ini, cosm_tgt, acc, iorder)
      rescaling_cost_sigma = rescaling_cost_sigma/log(R2_tgt/R1_tgt)

   END FUNCTION rescaling_cost_sigma

   REAL FUNCTION rescaling_cost_sigma_integrand(lnR, s, a_ini, a_tgt, flag, cosm_ini, cosm_tgt)

      ! Integrand for rescaling cost function in terms of sigma
      ! Note that this integrates over ln(R) in the target cosmology
      REAL, INTENT(IN) :: lnR
      REAL, INTENT(IN) :: s
      REAL, INTENT(IN) :: a_ini
      REAL, INTENT(IN) :: a_tgt
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm_ini
      TYPE(cosmology), INTENT(INOUT) :: cosm_tgt
      REAL :: R

      ! Change input from ln(R) to R and evaluate integrand
      R = exp(lnR)
      rescaling_cost_sigma_integrand = (1.-sigma(R/s, a_ini, flag, cosm_ini)/sigma(R, a_tgt, flag, cosm_tgt))**2

   END FUNCTION rescaling_cost_sigma_integrand

   REAL FUNCTION rescaling_cost_power(k1_tgt, k2_tgt, s, a_ini, a_tgt, cosm_ini, cosm_tgt)

      ! Rescaling cost function for minimizing Delta^2(k, a) differences
      ! Note that this integrates over a ln(k) range
      REAL, INTENT(IN) :: k1_tgt
      REAL, INTENT(IN) :: k2_tgt
      REAL, INTENT(IN) :: s
      REAL, INTENT(IN) :: a_ini
      REAL, INTENT(IN) :: a_tgt
      TYPE(cosmology), INTENT(INOUT) :: cosm_ini
      TYPE(cosmology), INTENT(INOUT) :: cosm_tgt
      INTEGER, PARAMETER :: flag = flag_rescale
      REAL, PARAMETER :: acc = acc_rescale
      INTEGER, PARAMETER :: iorder = iorder_rescale_integral

      ! Integrate and then normalise
      rescaling_cost_power = integrate_cosm(log(k1_tgt), log(k2_tgt), rescaling_cost_power_integrand, &
         s, a_ini, a_tgt, flag, cosm_ini, cosm_tgt, acc, iorder)
      rescaling_cost_power = rescaling_cost_power/log(k2_tgt/k1_tgt)

   END FUNCTION rescaling_cost_power

   REAL FUNCTION rescaling_cost_power_integrand(lnk, s, a_ini, a_tgt, flag, cosm_ini, cosm_tgt)

      ! Integrand for rescaling cost function in terms of power
      ! Note that this integrates over ln(k) in the target cosmology
      REAL, INTENT(IN) :: lnk
      REAL, INTENT(IN) :: s
      REAL, INTENT(IN) :: a_ini
      REAL, INTENT(IN) :: a_tgt
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm_ini
      TYPE(cosmology), INTENT(INOUT) :: cosm_tgt
      REAL :: k

      ! Change input from ln(k) to k and evaluate integrand
      k = exp(lnk)
      rescaling_cost_power_integrand = (1.-plin(k*s, a_ini, flag, cosm_ini)/plin(k, a_tgt, flag, cosm_tgt))**2

   END FUNCTION rescaling_cost_power_integrand

   SUBROUTINE ODE2_spherical(x, v, k, t, cosm, ti, tf, xi, vi, fv, n, imeth, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and creates arrays of x, v, t values
      ! ODE solver has a fixed number of time steps
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, INTENT(IN) :: k
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: imeth
      LOGICAL, INTENT(IN) :: ilog
      INTEGER :: i

      ! v' = fv
      INTERFACE      
         FUNCTION fv(x_interface, v_interface, k_interface, t_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: v_interface
            REAL, INTENT(IN) :: k_interface
            REAL, INTENT(IN) :: t_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION fv
      END INTERFACE

      ! Allocate arrays
      ALLOCATE (x(n), v(n), t(n))
      x = 0.
      v = 0.
      t = 0.

      ! xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
      x(1) = xi
      v(1) = vi

      ! Fill time array
      CALL fill_array(ti, tf, t, n, ilog)

      ! Loop over all time steps
      DO i = 1, n-1
         CALL ODE2_advance_cosm(x(i), x(i+1), v(i), v(i+1), k, t(i), t(i+1), cosm, velocity, fv, imeth)     
         IF (x(i+1) > dinf_spherical) EXIT ! Escape from the ODE solver when the perturbation is ~collapsed
      END DO

   END SUBROUTINE ODE2_spherical

   SUBROUTINE ODE2_adaptive_cosm(x, v, k, t, cosm, ti, tf, xi, vi, fv, acc, n, iode, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
      ! Adaptive, such that time steps are increased until convergence is achieved
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, INTENT(IN) :: k
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      REAL, INTENT(IN) :: acc ! Desired accuracy
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      REAL, ALLOCATABLE :: xh(:), th(:), vh(:)
      INTEGER :: j, in, ip, nn, np
      LOGICAL :: fail
      INTEGER, PARAMETER :: jmax = jmax_ode
      INTEGER, PARAMETER :: nmin = nmin_ode

      ! v' = fv
      INTERFACE       
         FUNCTION fv(x_interface, v_interface, k_interface, t_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: v_interface
            REAL, INTENT(IN) :: k_interface
            REAL, INTENT(IN) :: t_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION fv
      END INTERFACE

      ! Loop over attempts
      DO j = 1, jmax

         ! Number of points for this attempt
         nn = 1+(n-1)*2**(j-1)

         ! Solve the ODE with 'n' points
         CALL ODE2_cosm(x, v, k, t, cosm, ti, tf, xi, vi, fv, nn, iode, ilog)

         ! Check to see if integration has succeeded
         IF (j == 1 .OR. nn < nmin) THEN           
            fail = .TRUE.
         ELSE
            fail = .FALSE.
            DO ip = 1, np
               in = 2*ip-1
               IF ((.NOT. requal(xh(ip), x(in), acc)) .OR. (.NOT. requal(vh(ip), v(in), acc))) THEN
                  fail = .TRUE.
                  DEALLOCATE (xh, th, vh)
                  EXIT
               END IF
            END DO
         END IF

         ! If integration failed then store and move on to next attempt, otherwise exit
         IF (fail) THEN           
            xh = x
            vh = v
            th = t
            np = nn
         ELSE
            EXIT
         END IF

      END DO

      IF (fail) THEN
         STOP 'ODE_ADAPTIVE_COSMOLOGY: Error, failed to converge'
      ELSE
         xh = x
         vh = v
         th = t
         DEALLOCATE(x, v, t)
         ALLOCATE(x(n), v(n), t(n))
         DO ip = 1, n
            in = 1+(ip-1)*(nn-1)/(n-1)
            x(ip) = xh(in)
            v(ip) = vh(in)
            t(ip) = th(in)
         END DO
      END IF

   END SUBROUTINE ODE2_adaptive_cosm

   SUBROUTINE ODE2_cosm(x, v, k, t, cosm, ti, tf, xi, vi, fv, n, iode, ilog)

      ! Integrates second-order ODE: x'' + ax' + bx = c from ti to tf
      ! Outputs array of x(t) and v(t)
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, INTENT(IN) :: k
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      INTEGER :: i

      ! x'' = v' = fv(x, x', t)
      INTERFACE      
         FUNCTION fv(x_in, v_in, k_in, t_in, cosm_in)
            IMPORT :: cosmology     
            REAL, INTENT(IN) :: x_in, v_in, k_in, t_in
            TYPE(cosmology), INTENT(INOUT) :: cosm_in
         END FUNCTION fv
      END INTERFACE
   
      ! Allocate and set xi and vi to the initial values of x and v
      ALLOCATE(x(n), v(n))
      x(1) = xi
      v(1) = vi

      ! Fill the time array
      CALL fill_array(ti, tf, t, n, ilog)

      ! Advance the system through all n-1 time steps
      DO i = 1, n-1
         CALL ODE2_advance_cosm(x(i), x(i+1), v(i), v(i+1), k, t(i), t(i+1), cosm, velocity, fv, iode)
      END DO

   END SUBROUTINE ODE2_cosm

   SUBROUTINE ODE2_advance_cosm(x1, x2, v1, v2, k, t1, t2, cosm, fx, fv, imeth)

      ! Advance the ODE system from t1 to t2
      ! ODE integration methods
      ! imeth = 1: Crude
      ! imeth = 2: Mid-point
      ! imeth = 3: Fourth order Runge-Kutta
      REAL, INTENT(IN) :: x1
      REAL, INTENT(OUT) :: x2
      REAL, INTENT(IN) :: v1
      REAL, INTENT(OUT) :: v2
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: t1
      REAL, INTENT(IN) :: t2
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: imeth  
      REAL :: dt
      REAL :: kx1, kx2, kx3, kx4
      REAL :: kv1, kv2, kv3, kv4

      ! x' = fx; v' = fv
      INTERFACE     
         FUNCTION fx(x_interface, v_interface, k_interface, t_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: v_interface
            REAL, INTENT(IN) :: k_interface
            REAL, INTENT(IN) :: t_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION fx
         FUNCTION fv(x_interface, v_interface, k_interface, t_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: v_interface
            REAL, INTENT(IN) :: k_interface
            REAL, INTENT(IN) :: t_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION fv
      END INTERFACE

      ! Time step
      dt = real(t2-t1)

      IF (imeth == 1) THEN

         ! Crude method!
         kx1 = dt*fx(x1, v1, k, t1, cosm)
         kv1 = dt*fv(x1, v1, k, t1, cosm)

         x2 = x1+kx1
         v2 = v1+kv1

      ELSE IF (imeth == 2) THEN

         ! Mid-point method!
         kx1 = dt*fx(x1, v1, k, t1, cosm)
         kv1 = dt*fv(x1, v1, k, t1, cosm)
         kx2 = dt*fx(x1+kx1/2., v1+kv1/2., k, t1+dt/2., cosm)
         kv2 = dt*fv(x1+kx1/2., v1+kv1/2., k, t1+dt/2., cosm)

         x2 = x1+kx2
         v2 = v1+kv2

      ELSE IF (imeth == 3) THEN

         ! RK4 (Holy Christ, this is so fast compared to above methods)!
         kx1 = dt*fx(x1, v1, k, t1, cosm)
         kv1 = dt*fv(x1, v1, k, t1, cosm)
         kx2 = dt*fx(x1+kx1/2., v1+kv1/2., k, t1+dt/2., cosm)
         kv2 = dt*fv(x1+kx1/2., v1+kv1/2., k, t1+dt/2., cosm)
         kx3 = dt*fx(x1+kx2/2., v1+kv2/2., k, t1+dt/2., cosm)
         kv3 = dt*fv(x1+kx2/2., v1+kv2/2., k, t1+dt/2., cosm)
         kx4 = dt*fx(x1+kx3, v1+kv3, k, t1+dt, cosm)
         kv4 = dt*fv(x1+kx3, v1+kv3, k, t1+dt, cosm)

         x2 = x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.
         v2 = v1+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.

      ELSE

         STOP 'ODE_ADVANCE: Error, imeth specified incorrectly'

      END IF

   END SUBROUTINE ODE2_advance_cosm

   REAL RECURSIVE FUNCTION integrate_1_cosm(a, b, f, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a ! Integration lower limit
      REAL, INTENT(IN) :: b ! Integration upper limit
      REAL, EXTERNAL :: f
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL, INTENT(IN) :: acc ! Accuracy
      INTEGER, INTENT(IN) :: iorder ! Order for integration
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_1_cosm = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.d0
         sum_n = 0.d0
         sum_old = 0.d0
         sum_new = 0.d0

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, cosm)
               f2 = f(b, cosm)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*real(i-1)/real(n-1)
                  fx = f(x, cosm)
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
                  STOP 'INTEGRATE_COSM_1: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_COSM_1: Integration timed out'
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

         integrate_1_cosm = real(sum_new)

      END IF

   END FUNCTION integrate_1_cosm

   REAL RECURSIVE FUNCTION integrate_2_cosm(a, b, f, y, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a ! Integration lower limit for first arguement in 'f'
      REAL, INTENT(IN) :: b ! Integration upper limit for first arguement in 'f'
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: y ! Second argument in 'f'
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL, INTENT(IN) :: acc ! Accuracy
      INTEGER, INTENT(IN) :: iorder ! Order for integration
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x_interface, y_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: y_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_2_cosm = 0.

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
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, y, cosm)
               f2 = f(b, y, cosm)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*real(i-1)/real(n-1)
                  fx = f(x, y, cosm)
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
                  STOP 'INTEGRATE_COSM_2: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_COSM_2: Integration timed out'
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
         
         integrate_2_cosm = real(sum_new)

      END IF

   END FUNCTION integrate_2_cosm

   REAL RECURSIVE FUNCTION integrate_3_cosm(a, b, f, y, z, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a ! Integration lower limit for first argument in 'f'
      REAL, INTENT(IN) :: b ! Integration upper limit for first argument in 'f'
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: y ! Second argument in 'f'
      REAL, INTENT(IN) :: z ! Third argument in 'f'
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL, INTENT(IN) :: acc ! Accuracy
      INTEGER, INTENT(IN) :: iorder ! Order for integration
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x_interface, y_interface, z_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: y_interface
            REAL, INTENT(IN) :: z_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_3_cosm = 0.

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
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, y, z, cosm)
               f2 = f(b, y, z, cosm)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*real(i-1)/real(n-1)
                  fx = f(x, y, z, cosm)
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
                  STOP 'INTEGRATE_COSM_3: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.     
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_COSM_3: Integration timed out'
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

         integrate_3_cosm = real(sum_new)

      END IF

   END FUNCTION integrate_3_cosm

   REAL RECURSIVE FUNCTION integrate_3_flag_cosm(a, b, f, y, z, flag, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a ! Integration lower limit for first argument in 'f'
      REAL, INTENT(IN) :: b ! Integration upper limit for first argument in 'f'
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: y ! Second argument in 'f'
      REAL, INTENT(IN) :: z ! Third argument in 'f'
      INTEGER, INTENT(IN) :: flag ! Flag argument in 'f'
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL, INTENT(IN) :: acc ! Accuracy
      INTEGER, INTENT(IN) :: iorder ! Order for integration
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x_interface, y_interface, z_interface, flag_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: y_interface
            REAL, INTENT(IN) :: z_interface
            INTEGER, INTENT(IN) :: flag_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_3_flag_cosm = 0.

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
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, y, z, flag, cosm)
               f2 = f(b, y, z, flag, cosm)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*real(i-1)/real(n-1)
                  fx = f(x, y, z, flag, cosm)
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
                  STOP 'INTEGRATE_COSM_4: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_COSM_4: Integration timed out'
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

         integrate_3_flag_cosm = real(sum_new)

      END IF

   END FUNCTION integrate_3_flag_cosm

   REAL RECURSIVE FUNCTION integrate_4_flag_cosm(a, b, f, y, z, z1, flag, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a ! Integration lower limit for first argument in 'f'
      REAL, INTENT(IN) :: b ! Integration upper limit for first argument in 'f'
      REAL, EXTERNAL :: f
      REAL, INTENT(IN) :: y ! Second argument in 'f'
      REAL, INTENT(IN) :: z ! Third argument in 'f'
      REAL, INTENT(IN) :: z1 ! Fourth argument in 'f'
      INTEGER, INTENT(IN) :: flag ! Flag argument in 'f'
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL, INTENT(IN) :: acc ! Accuracy
      INTEGER, INTENT(IN) :: iorder ! Order for integration
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x_interface, y_interface, z_interface, z1_interface, flag_interface, cosm_interface)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_interface
            REAL, INTENT(IN) :: y_interface
            REAL, INTENT(IN) :: z_interface
            REAL, INTENT(IN) :: z1_interface
            INTEGER, INTENT(IN) :: flag_interface
            TYPE(cosmology), INTENT(INOUT) :: cosm_interface
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_4_flag_cosm = 0.

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
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, y, z, z1, flag, cosm)
               f2 = f(b, y, z, z1, flag, cosm)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*real(i-1)/real(n-1)
                  fx = f(x, y, z, z1, flag, cosm)
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
                  STOP 'INTEGRATE_COSM_4: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_COSM_4: Integration timed out'
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

         integrate_4_flag_cosm = real(sum_new)

      END IF

   END FUNCTION integrate_4_flag_cosm

   REAL RECURSIVE FUNCTION integrate_4_flag_2cosm(a, b, f, y, z, z1, flag, cosm1, cosm2, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      REAL, INTENT(IN) :: a                   ! Integration lower limit for first argument in 'f'
      REAL, INTENT(IN) :: b                   ! Integration upper limit for first argument in 'f'
      REAL, EXTERNAL :: f                     ! Function to integrate over
      REAL, INTENT(IN) :: y                   ! Second argument in 'f'
      REAL, INTENT(IN) :: z                   ! Third argument in 'f'
      REAL, INTENT(IN) :: z1                  ! Fourth argument in 'f'
      INTEGER, INTENT(IN) :: flag             ! Flag argument in 'f'
      TYPE(cosmology), INTENT(INOUT) :: cosm1 ! Cosmology
      TYPE(cosmology), INTENT(INOUT) :: cosm2 ! Cosmology
      REAL, INTENT(IN) :: acc                 ! Accuracy
      INTEGER, INTENT(IN) :: iorder           ! Order for integration
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x_int, y_int, z_int, z1_int, flag_int, cosm1_int, cosm2_int)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x_int
            REAL, INTENT(IN) :: y_int
            REAL, INTENT(IN) :: z_int
            REAL, INTENT(IN) :: z1_int
            INTEGER, INTENT(IN) :: flag_int
            TYPE(cosmology), INTENT(INOUT) :: cosm1_int
            TYPE(cosmology), INTENT(INOUT) :: cosm2_int
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_4_flag_2cosm = 0.

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
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = f(a, y, z, z1, flag, cosm1, cosm2)
               f2 = f(b, y, z, z1, flag, cosm1, cosm2)
               sum_2n = 0.5d0*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = a+(b-a)*real(i-1)/real(n-1)
                  fx = f(x, y, z, z1, flag, cosm1, cosm2)
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
                  STOP 'INTEGRATE_COSM_4: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE IF (j == jmax) THEN
               pass = .FALSE.
               STOP 'INTEGRATE_COSM_4: Integration timed out'
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

         integrate_4_flag_2cosm = real(sum_new)

      END IF

   END FUNCTION integrate_4_flag_2cosm

END MODULE cosmology_functions