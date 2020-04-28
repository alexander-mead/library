MODULE cosmology_functions

   USE interpolate
   USE constants
   USE file_info
   USE array_operations

   IMPLICIT NONE

   PRIVATE

   ! Type
   PUBLIC :: cosmology

   ! Basic routines
   PUBLIC :: assign_cosmology
   PUBLIC :: init_cosmology
   PUBLIC :: print_cosmology
   PUBLIC :: assign_init_cosmology

   ! Scale factor and z
   PUBLIC :: scale_factor_z
   PUBLIC :: redshift_r
   PUBLIC :: redshift_a
   PUBLIC :: scale_factor_r

   ! Friedmann
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

   ! Densities
   PUBLIC :: comoving_critical_density
   PUBLIC :: comoving_matter_density
   PUBLIC :: physical_critical_density
   PUBLIC :: physical_matter_density

   ! Linear growth
   PUBLIC :: ungrow
   PUBLIC :: grow
   PUBLIC :: grow_CPT
   PUBLIC :: grow_Linder
   PUBLIC :: growth_rate
   PUBlIC :: growth_rate_Linder
   PUBLIC :: acc_growth

   ! Spherical collapse
   PUBLIC :: Dv_BryanNorman
   PUBLIC :: Dv_Mead
   PUBLIC :: Dv_Spherical
   PUBLIC :: dc_NakamuraSuto
   PUBLIC :: dc_Mead
   PUBLIC :: dc_Spherical

   ! Power and correlation
   PUBLIC :: Pk_Delta
   PUBLIC :: Delta_Pk
   PUBLIC :: p_lin
   PUBLIC :: calculate_plin
   PUBLIC :: p_dewiggle
   PUBLIC :: sigma8
   PUBLIC :: sigma
   PUBLIC :: sigmaV
   PUBLIC :: neff
   PUBLIC :: ncur
   PUBLIC :: xi_lin
   PUBLIC :: flag_power_total
   PUBLIC :: flag_power_cold
   PUBLIC :: flag_power_cold_unorm
   PUBLIC :: norm_sigma8
   PUBLIC :: norm_none
   PUBLIC :: itk_none
   PUBLIC :: itk_EH
   PUBLIC :: itk_DEFW
   PUBLIC :: itk_CAMB
   PUBLIC :: itk_external

   ! CAMB interface
   PUBLIC :: get_CAMB_power
   PUBLIC :: CAMB_nonlinear_HALOFIT_Smith
   PUBLIC :: CAMB_nonlinear_HALOFIT_Bird
   PUBLIC :: CAMB_nonlinear_HALOFIT_Takahashi 
   PUBLIC :: CAMB_nonlinear_HMcode2015
   PUBLIC :: CAMB_nonlinear_HMcode2016

   ! HALOFIT
   PUBLIC :: calculate_HALOFIT_a
   PUBLIC :: calculate_HALOFIT
   PUBLIC :: HALOFIT_init
   PUBLIC :: HALOFIT_Smith
   PUBLIC :: HALOFIT_Smith_paper
   PUBLIC :: HALOFIT_Bird
   PUBLIC :: HALOFIT_Bird_paper
   PUBLIC :: HALOFIT_Takahashi
   PUBLIC :: HALOFIT_CAMB
   PUBLIC :: HALOFIT_CLASS
   
   INTERFACE integrate_cosm
      MODULE PROCEDURE integrate_cosm_1
      MODULE PROCEDURE integrate_cosm_2
      MODULE PROCEDURE integrate_cosm_3
      MODULE PROCEDURE integrate_cosm_4
   END INTERFACE integrate_cosm

   ! Contains cosmological parameters that only need to be calculated once
   TYPE cosmology

      ! Primary parameters
      CHARACTER(len=256) :: name            ! Name for cosmological model
      REAL :: Om_m, Om_b, Om_v, Om_w, m_nu  ! Densities
      REAL :: h, ns, w, wa, m_wdm, YH       ! Cosmological parameters
      REAL :: a1, a2, nstar, ws, am, dm, wm ! Dark-energy parameters
      REAL :: z_CMB, T_CMB, neff            ! CMB/radiation parameters
      REAL :: Om_m_pow, Om_b_pow, h_pow     ! Cosmological parameters used for P(k) if different from background
      REAL :: b0, b1, b2, b3, b4            ! BDE parameters
      REAL :: A_bump, k_bump, sigma_bump    ! Power-spectrum bump   
      REAL :: Theat                         ! AGN temperature
      REAL :: Lbox                          ! Box size [Mpc/h]
      INTEGER :: bump                       ! Type of bump to add to P(k)
      INTEGER :: norm_method                ! Power normalisation scheme
      INTEGER :: iw                         ! Switch for dark-energy type
      INTEGER :: itk                        ! Switch for transfer function type
      LOGICAL :: box                        ! Constrain the calculation to take place in a box?    
      LOGICAL :: warm                       ! Is DM warm?
      LOGICAL :: power_Omegas               ! Are the Omegas for the background different from those for the perturbations?
      LOGICAL :: derive_gas_numbers         ! Should mu_e and mu_p be derived or not?

      ! Variables that might be primary or secondary depening on the power normalisation
      REAL :: kpiv, As, kval, pval, sig8    ! Power spectrum normalisation   
      REAL :: mue, mup                      ! Gas parameters

      ! Derived parameters
      REAL :: A, Gamma, k                ! Power spectrum amplitude and shape parameter for DEFW
      REAL :: Om, Om_k, Om_c, Om_g, Om_r ! Derived Omegas
      REAL :: Om_nu, f_nu, a_nu, z_nu    ! Neutrinos
      REAL :: Om_nu_rad, omega_nu        ! Neutrinos
      REAL :: omega_m, omega_b, omega_c  ! Physical densities
      REAL :: Om_c_pow                   ! Cosmological parameters used for P(k) if different from background
      REAL :: age, horizon               ! Derived distance/time
      REAL :: YHe                        ! Derived thermal parameters
      REAL :: Om_ws, astar, a1n, a2n     ! Derived DE parameters
      REAL :: gnorm                      ! Growth-factor normalisation
      REAL :: kbox                       ! Wavenumber of box mode
      LOGICAL :: scale_dependent_growth  ! Is the linear growth scale dependent in this cosmology?
      LOGICAL :: trivial_cold            ! Is the cold spectrum trivially related to the matter spectrum?    
      
      ! Look-up tables that are filled during a calculation
      REAL, ALLOCATABLE :: log_plin(:), log_k_plin(:), log_plina(:, :), log_a_plin(:)        ! Arrays for input linear P(k)
      REAL, ALLOCATABLE :: log_k_Tcold(:), Tcold(:,:), p_wiggle(:)                           ! Arrays for cold T(k) and wiggle P(k)
      REAL, ALLOCATABLE :: log_r_sigma(:), log_a_sigma(:)                                    ! Arrays for R and a for sigma(R, a)
      REAL, ALLOCATABLE :: log_sigma(:), log_sigmaa(:, :), log_sigmac(:), log_sigmaca(:,:)   ! Arrays for sigma(R, a)
      REAL, ALLOCATABLE :: log_a_growth(:), log_growth(:), growth_rate(:), log_acc_growth(:) ! Arrays for growth
      REAL, ALLOCATABLE :: log_p(:), log_a_p(:)        ! Arrays for distance (particle horizon)
      REAL, ALLOCATABLE :: log_t(:), log_a_t(:)        ! Arrays for time
      REAL, ALLOCATABLE :: log_a_dcDv(:), dc(:), Dv(:) ! Arrays for spherical-collapse parameters
      REAL, ALLOCATABLE :: log_a_Xde(:), log_Xde(:)    ! Arrays for dark-energy density
      INTEGER :: n_growth, n_p, n_t, n_dcDv, nr_sigma, na_sigma, nk_plin, nk_Tcold, na_plin, n_Xde ! Number of array entries
      REAL :: amin_sigma, amax_sigma                                                               ! Ranges of arrays
      LOGICAL :: has_distance, has_growth, has_sigma, has_spherical, has_power, has_time, has_Xde  ! What has been calculated
      LOGICAL :: has_wiggle
      LOGICAL :: is_init, is_normalised ! Flags to check if things have been done  

      ! For random cosmologies
      !LOGICAL :: seeded
      !INTEGER :: iseed

      ! Verbose
      LOGICAL :: verbose

      ! Error handling
      INTEGER :: status

   END TYPE cosmology

   ! Global parameters
   REAL, PARAMETER :: acc_cosm = 1e-4 ! Global accuacy for the cosmological integrations

   ! Writing to screen parameters
   REAL, PARAMETER :: small_curve = 1e-5 ! Used to decide if writing curvature to screen to avoid silly numbers

   ! Equation of state parameters for cosmological fluids
   REAL, PARAMETER :: w_c = 0.    ! CDM
   REAL, PARAMETER :: w_b = 0.    ! Baryons
   REAL, PARAMETER :: w_g = 1./3. ! Photons
   REAL, PARAMETER :: w_v = -1.   ! Vacuum energy

   ! Dark energy density
   INTEGER, PARAMETER :: iw_LCDM = 1     ! Vacuum dark energy
   INTEGER, PARAMETER :: iw_QUICC = 2    ! QUICC dark energy
   INTEGER, PARAMETER :: iw_waCDM = 3    ! w(a) dark energy
   INTEGER, PARAMETER :: iw_wCDM = 4     ! Constant w dark energy
   INTEGER, PARAMETER :: iw_IDE1 = 5     ! Intermediate dark energy model 1
   INTEGER, PARAMETER :: iw_IDE2 = 6     ! Intermediate dark energy model 2
   INTEGER, PARAMETER :: iw_IDE3 = 7     ! Intermediate dark energy model 3
   INTEGER, PARAMETER :: iw_BDE = 8      ! Bound dark energy (1812.01133)  
   REAL, PARAMETER :: amin_Xde = 1e-4    ! Minimum scale factor for direction integration to get Xde
   REAL, PARAMETER :: amax_Xde = 1.      ! Maximum scale factor for direction integration to get Xde
   INTEGER, PARAMETER :: n_Xde = 128     ! Number of points for Xde interpolation
   LOGICAL, PARAMETER :: tabulate_Xde = .TRUE.        ! Tabulate Xde for interpolation
   REAL, PARAMETER :: acc_integration_Xde = acc_cosm  ! Accuracy for direct integration of dark energy density
   INTEGER, PARAMETER :: iorder_integration_Xde = 3   ! Polynomial order for time integration
   INTEGER, PARAMETER :: iorder_interpolation_Xde = 3 ! Polynomial order for time interpolation
   INTEGER, PARAMETER :: ifind_interpolation_Xde = 3  ! Finding scheme for time interpolation
   INTEGER, PARAMETER :: imeth_interpolation_Xde = 2  ! Method for time interpolation

   ! Distance
   ! Changing to linear integer finding provides very little speed increase
   REAL, PARAMETER :: amin_distance = 1e-4                 ! Minimum scale factor in look-up table
   REAL, PARAMETER :: amax_distance = 1.                   ! Maximum scale factor in look-up table
   INTEGER, PARAMETER :: n_distance = 128                  ! Number of scale factor entries in look-up table
   REAL, PARAMETER :: atay_distance = 1e-5                 ! Below this do a Taylor expansion to avoid divergence
   INTEGER, PARAMETER :: iorder_integration_distance = 3   ! Polynomial order for distance integration
   INTEGER, PARAMETER :: iorder_interpolation_distance = 3 ! Polynomial order for distance interpolation
   INTEGER, PARAMETER :: ifind_interpolation_distance = 3  ! Finding scheme for distance interpolation
   INTEGER, PARAMETER :: imeth_interpolation_distance = 2  ! Method for distance interpolation
   INTEGER, PARAMETER :: iorder_inversion_distance = 3     ! Polynomial order for distance inversion
   INTEGER, PARAMETER :: ifind_inversion_distance = 3      ! Finding scheme for distance inversion
   INTEGER, PARAMETER :: imeth_inversion_distance = 2      ! Method for distance inversion

   ! Time
   ! Changing to linear integer finding provides very little speed increase
   REAL, PARAMETER :: amin_time = 1e-4                 ! Minimum scale factor in look-up table
   REAL, PARAMETER :: amax_time = 1.                   ! Maximum scale factor in look-up table
   INTEGER, PARAMETER :: n_time = 128                  ! Number of scale factor entries in look-up table
   REAL, PARAMETER :: atay_time = 1e-5                 ! Below this do a Taylor expansion to avoid divergence
   INTEGER, PARAMETER :: iorder_integration_time = 3   ! Polynomial order for time integration
   INTEGER, PARAMETER :: iorder_interpolation_time = 3 ! Polynomial order for time interpolation
   INTEGER, PARAMETER :: ifind_interpolation_time = 3  ! Finding scheme for time interpolation
   INTEGER, PARAMETER :: imeth_interpolation_time = 2  ! Method for time interpolation

   ! Power
   REAL, PARAMETER :: kmin_plin = 0.           ! Power below this wavenumber is set to zero [h/Mpc]
   REAL, PARAMETER :: kmax_plin = 1e8          ! Power above this wavenumber is set to zero [h/Mpc]
   LOGICAL, PARAMETER :: plin_extrap = .FALSE. ! Extrapolate high-k power assumning P(k) ~ ln(k)^2 k^(n-3) (otherwise power law)?
   INTEGER, PARAMETER :: itk_none = 0          ! Pure power-law spectrum
   INTEGER, PARAMETER :: itk_EH = 1            ! Eisenstein and Hu linear spectrum
   INTEGER, PARAMETER :: itk_CAMB = 2          ! CAMB linear spectrum
   INTEGER, PARAMETER :: itk_DEFW = 3          ! DEFW linear spectrum
   INTEGER, PARAMETER :: itk_external = 4      ! DEFW linear spectrum
   INTEGER, PARAMETER :: norm_sigma8 = 1       ! Normalise power spectrum via sigma8 value
   INTEGER, PARAMETER :: norm_value = 2        ! Normalise power spectrum via specifying a value at a k 
   INTEGER, PARAMETER :: norm_As = 3           ! Normalise power spectrum vis As value as in CAMB
   INTEGER, PARAMETER :: norm_none = 4         ! Power spectrum does not need to be normalised
   INTEGER, PARAMETER :: flag_power_total = 1      ! Flag to get the total matter power spectrum
   INTEGER, PARAMETER :: flag_power_cold = 2       ! Flag to get the cold (CDM+baryons) power spectrum with 1+delta = rho_cold/mean_rho_matter
   INTEGER, PARAMETER :: flag_power_cold_unorm = 3 ! Flag to get the cold (CDM+baryons) power spectrum with 1+delta = rho_cold/mean_rho_cold

   ! Correlation function
   ! TODO: This works very poorly
   INTEGER, PARAMETER :: method_xi = 1         ! Method for xi integration
   INTEGER, PARAMETER :: iorder_xi = 3         ! Polynomial order for xi(r) integration
   INTEGER, PARAMETER :: min_humps_xi = 5      ! Minimum number of humps fox xi(r) integration if using humps integration
   INTEGER, PARAMETER :: max_humps_xi = 100000 ! Maximum number of humps fox xi(r) integration if using humps integration
   REAL, PARAMETER :: rsplit_xi = 10.          ! Value of r to split alpha values for conventional integration
   REAL, PARAMETER :: alpha_lo_xi = 2.         ! Low r value of alpha for conventional integration
   REAL, PARAMETER :: alpha_hi_xi = 1.5        ! High r value of alpha for conventional integration

   ! CAMB interface
   REAL, PARAMETER :: kmin_CAMB = 1e-3         ! Minimum wavenumber to use if rebinning
   REAL, PARAMETER :: kmax_CAMB = 1e2          ! Maximum wavenumber to get (also to use for rebinning)
   INTEGER, PARAMETER :: nk_CAMB = 128         ! Number of k values to use if rebinning
   REAL, PARAMETER :: amin_CAMB = 0.1          ! Minimum a value for P(k,a) tables when linear growth is scale dependent
   REAL, PARAMETER :: amax_CAMB = 1.0          ! Maximum a value for P(k,a) tables when linear growth is scale dependent
   INTEGER, PARAMETER :: na_CAMB = 16          ! Number of a values for linear P(k,a) tables if growth is scale dependent
   REAL, PARAMETER :: pk_min_CAMB = 1e-10      ! Minimum value of power at low k (remove k with less than this) 
   REAL, PARAMETER :: nmax_CAMB = 2.           ! How many times more to go than kmax due to inaccuracy near k limit
   LOGICAL, PARAMETER :: rebin_CAMB = .FALSE.  ! Should we rebin CAMB or just use default k?
   INTEGER, PARAMETER :: iorder_rebin_CAMB = 3 ! Polynomial order for interpolation on CAMB rebinning
   INTEGER, PARAMETER :: ifind_rebin_CAMB = 3  ! Finding scheme for interpolation on CAMB rebinning (*definitely* not linear)
   INTEGER, PARAMETER :: imeth_rebin_CAMB = 2  ! Method for interpolation on CAMB rebinning

   ! Linear power interpolation
   INTEGER, PARAMETER :: iorder_interpolation_plin = 3 ! Order for interpolation
   INTEGER, PARAMETER :: ifind_interpolation_plin = 3  ! Finding scheme in table (only linear if rebinning)
   INTEGER, PARAMETER :: imeth_interpolation_plin = 2  ! Method for polynomials

   ! Growth ODE
   REAL, PARAMETER :: ainit_growth = 1e-3                    ! Starting value for growth integratiton (should start | Omega_m(a)=1)
   REAL, PARAMETER :: amax_growth = 1.                       ! Finishing value for growth integratiton (should be a=1)
   INTEGER, PARAMETER :: n_growth = 128                      ! Number of entries for growth look-up tables
   REAL, PARAMETER :: acc_ODE_growth = acc_cosm              ! Accuracy parameter for growth ODE solving
   INTEGER, PARAMETER :: imeth_ODE_growth = 3                ! Method for solving growth ODE
   INTEGER, PARAMETER :: iorder_ODE_interpolation_growth = 3 ! Polynomial order for growth interpolation for ODE solution
   INTEGER, PARAMETER :: ifind_ODE_interpolation_growth = 3  ! Finding scheme for growth interpolation for ODE solution
   INTEGER, PARAMETER :: imeth_ODE_interpolation_growth = 2  ! Method for growth interpolation for ODE solution

   ! Growth integral (LCDM only)
   REAL, PARAMETER :: acc_integral_grow = acc_cosm ! Accuracy parameter for growth integral solving (wCDM only)
   INTEGER, PARAMETER :: iorder_integral_grow = 3  ! Polynomial order for growth integral solving (wCDM only)

   ! Growth interpolation
   ! Changing to linear integer finding is wrong because tables not stored lin-log
   INTEGER, PARAMETER :: iorder_interpolation_grow = 3 ! Polynomial order for growth interpolation
   INTEGER, PARAMETER :: ifind_interpolation_grow = 3  ! Finding scheme for growth interpolation
   INTEGER, PARAMETER :: imeth_interpolation_grow = 2  ! Method for growth interpolation

   ! Growth rate interpolation
   ! Changing to linear integer finding is wrong because tables not stored lin-log
   INTEGER, PARAMETER :: iorder_interpolation_rate = 3 ! Polynomial order for growth rate interpolation for ODE solution
   INTEGER, PARAMETER :: ifind_interpolation_rate = 3  ! Finding scheme for growth rate interpolation
   INTEGER, PARAMETER :: imeth_interpolation_rate = 2  ! Method for growth rate interpolation

   ! Accumulated growth integration
   INTEGER, PARAMETER :: iorder_integration_agrow = 3 ! Polynomial order for growth interpolation for ODE solution

   ! Accumualted growth interpolation
   ! Changing to linear integer finding is wrong because tables not stored lin-log
   INTEGER, PARAMETER :: iorder_interpolation_agrow = 3 ! Polynomial order for growth interpolation for ODE solution
   INTEGER, PARAMETER :: ifind_interpolation_agrow = 3  ! Finding scheme for accumulated growth rate interpolation
   INTEGER, PARAMETER :: imeth_interpolation_agrow = 2  ! Method for accumulated growth rate interpolation

   ! sigma(R) integration
   REAL, PARAMETER :: alpha_sigma = 3.     ! Exponent to increase speed (1 is terrible, 2, 3, 4 all okay)
   REAL, PARAMETER :: acc_sigma = acc_cosm ! Accuracy parameter for sigma(R) integration
   INTEGER, PARAMETER :: iorder_sigma = 3  ! Polynomial order for sigma(R) integration
   
   ! sigma(R) tabulation and interpolation
   ! Changing to linear integer finding provides very little speed increase
   REAL, PARAMETER :: rmin_sigma = 1e-4                 ! Minimum r value (NB. sigma(R) needs to be power-law below)
   REAL, PARAMETER :: rmax_sigma = 1e3                  ! Maximum r value (NB. sigma(R) needs to be power-law above)
   INTEGER, PARAMETER :: nr_sigma = 128                 ! Number of r entries for sigma(R) tables
   REAL, PARAMETER :: amin_sigma = 0.1                  ! Minimum a value for sigma(R,a) tables when growth is scale dependent
   REAL, PARAMETER :: amax_sigma = 1.0                  ! Maximum a value for sigma(R,a) tables when growth is scale dependent
   INTEGER, PARAMETER :: na_sigma = 16                  ! Number of a values for sigma(R,a) tables
   INTEGER, PARAMETER :: iorder_interpolation_sigma = 3 ! Polynomial order for sigma(R) interpolation 
   INTEGER, PARAMETER :: ifind_interpolation_sigma = 3  ! Finding scheme for sigma(R) interpolation (changing to linear not speedy)
   INTEGER, PARAMETER :: imeth_interpolation_sigma = 2  ! Method for sigma(R) interpolation
   INTEGER, PARAMETER :: sigma_cold_store = flag_power_cold_unorm ! Which version of the cold spectrum should be tabulated (can be 0 for none)
   !INTEGER, PARAMETER :: sigma_cold_store = 0           ! Which version of the cold spectrum should be tabulated (can be 0 for none)

   ! sigma_v(R) integration
   REAL, PARAMETER :: alpha_sigmaV = 3.     ! Exponent to increase integration speed
   REAL, PARAMETER :: acc_sigmaV = acc_cosm ! Accuracy parameter for sigma(R) integration
   INTEGER, PARAMETER :: iorder_sigmaV = 3  ! Polynomial order for sigmaV(R) integration

   ! neff integration
   REAL, PARAMETER :: alpha_dsigma = 2.     ! Exponent to increase integration speed (3 seemed to give inaccuracies; no idea why)
   REAL, PARAMETER :: acc_dsigma = acc_cosm ! Accuracy parameter for sigma(R) integration
   INTEGER, PARAMETER :: iorder_dsigma = 3  ! Polynomial order for neff(R) integration

   ! ncur integration
   REAL, PARAMETER :: alpha_ddsigma = 2.     ! Exponent to increase integration speed (3 seemed to give inaccuracies; no idea why)
   REAL, PARAMETER :: acc_ddsigma = acc_cosm ! Accuracy parameter for sigma(R) integration
   INTEGER, PARAMETER :: iorder_ddsigma = 3  ! Polynomial order for neff(R) integration

   ! Spherical collapse
   INTEGER, PARAMETER :: imeth_ODE_spherical = 3 ! Method for spherical collapse ODE solving
   REAL, PARAMETER :: amax_spherical = 2.        ! Maximum scale factor to consider
   REAL, PARAMETER :: dmin_spherical = 1e-7      ! Minimum starting value for perturbation
   REAL, PARAMETER :: dmax_spherical = 1e-3      ! Maximum starting value for perturbation
   INTEGER, PARAMETER :: m_spherical = 128       ! Number of collapse scale-factors to try to calculate
   INTEGER, PARAMETER :: n_spherical = 100000    ! Number of points for ODE calculations
   REAL, PARAMETER :: dinf_spherical = 1e8       ! Value considered to be 'infinite' for the perturbation

   ! delta_c
   INTEGER, PARAMETER :: iorder_interpolation_dc = 3 ! Polynomial order for delta_c interpolation
   INTEGER, PARAMETER :: ifind_interpolation_dc = 3  ! Finding scheme for delta_c interpolation
   INTEGER, PARAMETER :: imeth_interpolation_dc = 2  ! Method for delta_c interpolation

   ! Delta_v
   INTEGER, PARAMETER :: iorder_interpolation_Dv = 3 ! Polynomial order for Delta_v interpolation
   INTEGER, PARAMETER :: ifind_interpolation_Dv = 3  ! Finding scheme for Delta_v interpolation
   INTEGER, PARAMETER :: imeth_interpolation_Dv = 2  ! Method for Delta_v interpolation

   ! Halofit
   INTEGER, PARAMETER :: HALOFIT_Smith = 1       ! Smith et al. (https://www.roe.ac.uk/~jap/haloes/)
   INTEGER, PARAMETER :: HALOFIT_Bird = 2        ! Bird et al. (2011; https://arxiv.org/abs/1109.4416)
   INTEGER, PARAMETER :: HALOFIT_Takahashi = 3   ! Takahashi et al. (2012; https://arxiv.org/abs/1208.2701)
   INTEGER, PARAMETER :: HALOFIT_CAMB = 4        ! Version as used in CAMB  (date; ???)
   INTEGER, PARAMETER :: HALOFIT_CLASS = 5       ! Version as used in CLASS (date; ???)
   INTEGER, PARAMETER :: HALOFIT_Smith_paper = 6 ! Smith et al. (2003; https://arxiv.org/abs/astro-ph/0207664)
   INTEGER, PARAMETER :: HALOFIT_Bird_paper = 7  ! Bird et al. ()

   ! CAMB non-linear
   INTEGER, PARAMETER :: CAMB_nonlinear_HALOFIT_Smith = 1
   INTEGER, PARAMETER :: CAMB_nonlinear_HALOFIT_Bird = 2
   INTEGER, PARAMETER :: CAMB_nonlinear_HALOFIT_Takahashi = 4
   INTEGER, PARAMETER :: CAMB_nonlinear_HMcode2015 = 8
   INTEGER, PARAMETER :: CAMB_nonlinear_HMcode2016 = 5

   ! General integrations
   INTEGER, PARAMETER :: jmin_integration = 5
   INTEGER, PARAMETER :: jmax_integration = 30

CONTAINS

   SUBROUTINE assign_cosmology(icosmo, cosm, verbose)

      ! Assigns the 'primary' cosmological parameters (primary according to my definition)
      ! This routine *only* assigns parameters, it does and should not do *any* calculations
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i
      REAL :: Xe, Xi

      ! Names of pre-defined cosmologies
      INTEGER, PARAMETER :: ncosmo = 337
      CHARACTER(len=256) :: names(0:ncosmo)
      names = ''
      names(0)  = 'User defined'
      names(1)  = 'Boring'
      names(2)  = 'WMAP7 (cosmo-OWLS)'
      names(3)  = 'Planck 2013 (cosmo-OWLS/BAHAMAS)'
      names(4)  = 'WMAP9 (BAHAMAS)'
      names(5)  = 'Open'
      names(6)  = 'Einstein-de Sitter'
      names(7)  = 'IDE I (user)'
      names(8)  = 'IDE II (user)'
      names(9)  = 'IDE III (user)'
      names(10) = 'IDE3'
      names(11) = 'IDE10'
      names(12) = 'LCDM (user)'
      names(13) = 'w(a)CDM (user)'
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
      names(63) = 'Planck 2015'
      names(64) = 'Planck 2015 (AGN 7.6)'
      names(65) = 'Planck 2015 (AGN 8.0)'
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
         DO i = 0, size(names)-1
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
      cosm%itk = itk_EH ! Default to Eisenstein & Hu

      ! Boring default cosmology
      cosm%Om_m = 0.3
      cosm%Om_b = 0.05
      cosm%Om_v = 1.-cosm%Om_m
      cosm%Om_w = 0.
      cosm%m_nu = 0.
      cosm%h = 0.7
      cosm%ns =  0.96
      cosm%w = -1.
      cosm%wa = 0.
      cosm%T_CMB = 2.725 ! CMB temperature [K]
      cosm%z_CMB = 1087. ! Redshift of the last-scatting surface
      cosm%neff = 3.046  ! Effective number of relativistic neutrinos
      cosm%YH = 0.76     ! Hydrogen mass fraction

      ! Dark energy
      cosm%iw = iw_LCDM ! Default LCDM

      ! Warm dark matter
      cosm%warm = .FALSE.
      cosm%m_wdm = 0. ! Particle mass [keV]

      ! AGN
      cosm%Theat = 10**7.8 ! AGN temperature

      !! Normalisation !!

      ! Overall power normalisaiton, should always be set to 1 and then will be changed later
      cosm%A = 1.

      ! Large-scale structure normalisation
      cosm%norm_method = norm_sigma8
      cosm%sig8 = 0.8              

      ! Large-scale structure normalisation at a specific wavenumber
      !cosm%norm_method = norm_value
      cosm%kval = 0.001
      cosm%pval = 0.1973236854e-06 ! Power value to get sig8 = 0.8 for a boring cosmology

      ! As CMB normalisation
      !cosm%norm_method = norm_As
      cosm%kpiv = 0.05 ! Wavenumber at which to define the normalisation
      cosm%As = 2.1e-9 ! This is a generally sensible value
      
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

      ! Set is/has flags to negative
      cosm%is_init = .FALSE.
      cosm%is_normalised = .FALSE.
      cosm%has_distance = .FALSE.
      cosm%has_growth = .FALSE.
      cosm%has_sigma = .FALSE.
      cosm%has_spherical = .FALSE.
      cosm%has_power = .FALSE.
      cosm%has_time = .FALSE.
      cosm%has_Xde = .FALSE.
      cosm%has_wiggle = .FALSE.

      ! For random cosmologies
      !cosm%seed = 0
      !cosm%seeded = .FALSE.

      ! Omegas for power spectrum if different from background cosmological parameters; false by default
      cosm%Om_m_pow = 0.
      cosm%Om_b_pow = 0.
      cosm%h_pow = 0.
      cosm%power_Omegas = .FALSE. 

      ! Gas options
      cosm%derive_gas_numbers = .TRUE.

      IF (icosmo == 0) THEN
         STOP 'TODO: implement user decision here'
      ELSE IF (icosmo == 1) THEN
         ! Boring - do nothing
      ELSE IF (icosmo == 2) THEN
         ! cosmo-OWLS - WMAP7 (1312.5462)
         cosm%Om_m = 0.272
         cosm%Om_b = 0.0455
         cosm%Om_v = 1.-cosm%Om_m
         cosm%h = 0.704
         cosm%sig8 = 0.81
         cosm%ns =  0.967
      ELSE IF (icosmo == 3) THEN
         ! Planck 2013 (cosmo-OWLS/BAHAMAS; 1312.5462/1603.02702; no neutrinos)
         cosm%Om_m = 0.3175
         cosm%Om_b = 0.0490
         cosm%Om_v = 1.-cosm%Om_m
         cosm%h = 0.6711
         cosm%ns =  0.9624
         cosm%sig8 = 0.8341
      ELSE IF (is_in_array(icosmo, [4, 61, 62, 70, 71, 75, 76, 77, 78])) THEN
         ! TODO: Should this use CAMB linear spectrum?
         !  4 - WMAP9 (BAHAMAS; 1712.02411; no neutrinos)
         ! 61 - WMAP9 with T_AGN = 10^7.6 K (low)
         ! 62 - WMAP9 with T_AGN = 10^8.0 K (high)
         ! 70 - WMAP9 with T_AGN = 10^7.0 K (extremely low)
         ! 71 - WMAP9 with T_AGN = 10^8.6 K (extrmely high)
         ! 75 - WMAP9 with 0.06eV neutrinos
         ! 76 - WMAP9 with 0.12eV neutrinos
         ! 77 - WMAP9 with 0.24eV neutrinos
         ! 78 - WMAP9 with 0.48eV neutrinos
         cosm%h = 0.7000
         cosm%Om_b = 0.0463
         cosm%Om_m = 0.2330+cosm%Om_b
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns =  0.9720
         cosm%sig8 = 0.8211
         cosm%derive_gas_numbers = .FALSE.
         cosm%mup = 0.61
         Xi = 1.08
         Xe = 1.17
         cosm%mue = cosm%mup*(Xe+Xi)/Xe
         IF (icosmo == 4)  cosm%Theat = 10**7.8 ! Standard AGN temperature
         IF (icosmo == 61) cosm%Theat = 10**7.6 ! Low AGN temperature
         IF (icosmo == 62) cosm%Theat = 10**8.0 ! High AGN temperature
         IF (icosmo == 70) cosm%Theat = 10**7.0 ! Extremely low AGN temperature
         IF (icosmo == 71) cosm%Theat = 10**8.6 ! Extremely high AGN temperature
         IF (icosmo == 75) cosm%m_nu = 0.06 ! 0.06eV (minimal) neutrino mass
         IF (icosmo == 76) cosm%m_nu = 0.12 ! 0.12eV neutrinos
         IF (icosmo == 77) cosm%m_nu = 0.24 ! 0.24eV neutrinos
         IF (icosmo == 78) cosm%m_nu = 0.48 ! 0.48eV neutrinos
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
         WRITE (*, *) 'a*:'
         READ (*, *) cosm%astar
         WRITE (*, *) 'n*:'
         READ (*, *) cosm%nstar
         WRITE (*, *) 'Om_w(a*):'
         READ (*, *) cosm%Om_ws
         cosm%Om_m = 0.3
         cosm%Om_v = 0.7
      ELSE IF (icosmo == 8) THEN
         ! IDE II model
         cosm%iw = iw_IDE2
         WRITE (*, *) 'n*:'
         READ (*, *) cosm%nstar
         WRITE (*, *) 'a*:'
         READ (*, *) cosm%astar
         WRITE (*, *) 'Om_w(a*):'
         READ (*, *) cosm%Om_ws
         cosm%Om_m = 0.3
         cosm%Om_w = 0.7
         cosm%Om_v = 0. !No vacuum necessary here
      ELSE IF (icosmo == 9) THEN
         ! IDE III model
         cosm%iw = iw_IDE3
         WRITE (*, *) 'a*:'
         READ (*, *) cosm%astar
         WRITE (*, *) 'Om_w(a*):'
         READ (*, *) cosm%Om_ws
         WRITE (*, *) 'w*:'
         READ (*, *) cosm%ws
         cosm%Om_m = 0.3
         cosm%Om_w = 0.7
         cosm%Om_v = 0.
      ELSE IF (icosmo == 10 .OR. icosmo == 11) THEN
         ! Fixed IDE II models
         cosm%iw = iw_IDE2
         cosm%Om_w = cosm%Om_v
         cosm%Om_v = 0.
         IF (icosmo == 10) cosm%nstar = 3.
         IF (icosmo == 11) cosm%nstar = 10.
         cosm%astar = 0.01
         cosm%Om_ws = 0.2
      ELSE IF (icosmo == 12) THEN
         ! Curved
         WRITE (*, *) 'Om_m:'
         READ (*, *) cosm%Om_m
         WRITE (*, *) 'Om_v:'
         READ (*, *) cosm%Om_v
      ELSE IF (icosmo == 13) THEN
         ! w(a)CDM
         cosm%iw = iw_waCDM
         WRITE (*, *) 'w0:'
         READ (*, *) cosm%w
         WRITE (*, *) 'wa:'
         READ (*, *) cosm%wa
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
      ELSE IF (icosmo == 24) THEN
         ! Random Mira Titan cosmology
         CALL random_Mira_Titan_cosmology(cosm)
         cosm%itk = itk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo == 25) THEN
         ! Random Franken Emu cosmology
         CALL random_Franken_Emu_cosmology(cosm)
         cosm%itk = itk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo == 26) THEN
         ! Boring with CAMB linear spectrum
         cosm%itk = itk_CAMB   ! Set to CAMB linear power
         cosm%Om_w = cosm%Om_v ! Necessary for CAMB
         cosm%Om_v = 0.        ! Necessary for CAMB
      ELSE IF (icosmo == 27) THEN
         ! Illustris; L = 75 Mpc/h
         cosm%itk = itk_CAMB ! Set to CAMB linear power
         cosm%Om_m = 0.3089
         cosm%Om_b = 0.0486
         cosm%Om_w = 1.-cosm%Om_m ! Necessary to be Omega_w for CAMB
         cosm%h = 0.6774
         cosm%ns =  0.9667
         cosm%sig8 = 0.8159
         cosm%Om_v = 0. ! Necessary for CAMB
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
         ! 37 - Multidark: WMAP5
         ! 43 - Multidark: WMAP5 with altered sigma_8
         ! 44 - Multidark: WMAP5 with Eisenstein & Hu transfer function
         cosm%h = 0.70
         cosm%Om_b = 0.0469
         cosm%Om_m = 0.27
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns =  0.95
         cosm%sig8 = 0.82 ! Seems wrong at z=0, data more like sigma_8 = 0.80
         cosm%itk = itk_CAMB ! CAMB
         IF(icosmo == 43) cosm%sig8 = 0.80 ! Check to see if better matches with lower sigma_8
         IF(icosmo == 44) cosm%itk = itk_EH ! Eisenstein & Hu T(k)
      ELSE IF (icosmo == 38) THEN
         ! Random cosmic emu model
         CALL random_Cosmic_Emu_cosmology(cosm)
         cosm%itk = itk_CAMB ! Set to CAMB linear power
         cosm%iw = iw_wCDM    ! Set to constant w dark energy
         cosm%Om_v = 0.      ! Necessary for CAMB
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
            cosm%itk = itk_EH
         ELSE
            cosm%itk = itk_CAMB
         END IF
      ELSE IF (icosmo == 41) THEN
         ! SCDM with high neutrino mass
         cosm%Om_m = 1.
         cosm%Om_v = 0.
         cosm%m_nu = 4.
      ELSE IF (icosmo == 56 .OR. icosmo == 42) THEN
         ! 56 - Planck 2018 (Plik, from Table 1 of https://arxiv.org/abs/1807.06209)
         ! 42 - Same, but with neutrino mass fixed to zero and nothing else changed
         cosm%itk = itk_CAMB
         cosm%h = 0.6732
         cosm%ns =  0.96605
         cosm%m_nu = 0.06
         IF(icosmo == 42) cosm%m_nu = 0.
         cosm%Om_m = 0.3158
         cosm%Om_b = 0.022383/cosm%h**2
         cosm%Om_v = 1.-cosm%Om_m
         cosm%sig8 = 0.8120
      ELSE IF(icosmo == 49) THEN
         ! CFHTLenS best-fitting cosmology (Heymans 2013; combined with WMAP 7)
         ! From first line of Table 3 of https://arxiv.org/pdf/1303.1808.pdf
         ! CFHTLenS combined with WMAP7 and R11
         cosm%Om_m = 0.255
         cosm%Om_b = 0.0437
         cosm%Om_v = 1.-cosm%om_m
         cosm%sig8 = 0.794
         cosm%ns =  0.967
         cosm%h = 0.717
         cosm%w = -1.
      ELSE IF(icosmo == 50) THEN
         ! CFHTLenS best-fitting cosmology (Kilbinger 2013; combined with WMAP 7)
         ! From first line in Table 3 of https://arxiv.org/pdf/1212.3338.pdf
         ! CFHTLenS combined with WMAP7
         cosm%Om_m = 0.274
         cosm%Om_b = 0.0456
         cosm%Om_v = 1.-cosm%om_m
         cosm%sig8 = 0.815
         cosm%ns =  0.966
         cosm%h = 0.702
         cosm%w = -1.
      ELSE IF (is_in_array(icosmo, [51, 52, 53, 54, 55])) THEN
         ! CAMB test problem cosmologies
         cosm%itk = itk_CAMB
         cosm%iw = iw_waCDM
         cosm%Om_v = 0.
         IF(icosmo == 51) THEN
            ! CAMB difference cosmology 1 (low sig8; low ns; low h)
            cosm%Om_m = 0.16947
            cosm%Om_b = 0.03972
            cosm%Om_w = 1.-cosm%Om_m
            cosm%sig8 = 0.61272
            cosm%ns =  0.70918
            cosm%h = 0.54980
            cosm%M_nu = 0.
            cosm%w = -1.
            cosm%wa = 0.
         ELSE IF(icosmo == 52) THEN
            ! CAMB difference cosmology 2 (low sig8; low ns; low h)
            cosm%Om_m = 0.16021
            cosm%Om_b = 0.05953
            cosm%Om_w = 1.-cosm%Om_m
            cosm%sig8 = 0.69869
            cosm%ns =  0.70630
            cosm%h = 0.52845
            cosm%M_nu = 0.0
            cosm%w = -1.
            cosm%wa = 0.
         ELSE IF(icosmo == 53) THEN
            ! CAMB difference cosmology 3 (low sig8; low ns; low h)
            cosm%Om_m = 0.17575
            cosm%Om_b = 0.03993
            cosm%Om_w = 1.-cosm%Om_m
            cosm%sig8 = 0.61299
            cosm%ns =  0.72517
            cosm%h = 0.44980
            cosm%M_nu = 0.
            cosm%w = -1.
            cosm%wa = 0.
         ELSE IF(icosmo == 54) THEN
            ! CAMB difference cosmology 4 (low sig8; low ns; low h)
            cosm%Om_m = 0.15586
            cosm%Om_b = 0.05109
            cosm%Om_w = 1.-cosm%Om_m
            cosm%sig8 = 0.65709
            cosm%ns =  0.70027
            cosm%h = 0.44819
            cosm%M_nu = 0.20378
            cosm%w = -1.01372
            cosm%wa = -0.90823
         ELSE IF(icosmo == 55) THEN
            ! CAMB difference cosmology 5 (low h; high Om_b/Om_m)
            cosm%Om_m = 0.18504
            cosm%Om_b = 0.06696
            cosm%Om_w = 1.-cosm%Om_m
            cosm%sig8 = 0.77520
            cosm%ns =  0.98360
            cosm%h = 0.42278
            cosm%M_nu = 0.07639
            cosm%w = -1.26075
            cosm%wa = -0.07341
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
         cosm%norm_method = norm_value ! Normalise like this to prevent bump annoying sigma8        
         cosm%bump = 1
         cosm%A_bump = 0.08
         cosm%k_bump = 5.
         cosm%sigma_bump = 0.5      
      ELSE IF (is_in_array(icosmo, [66, 67, 68, 69, 72, 73, 74, 79, 80, 81])) THEN
         ! Axel bump cosmologies
         cosm%h = 0.7
         cosm%Om_b = 0.05
         cosm%Om_m = 0.30
         cosm%ns =  0.96
         cosm%Om_v = 1.-cosm%Om_m
         cosm%norm_method = norm_value
         cosm%pval = 1.995809e-7 ! Gives sigma8 = 0.8 for no bump   
         cosm%itk = itk_CAMB
         IF (cosm%itk == itk_CAMB) THEN
            cosm%Om_w = cosm%Om_v
            cosm%Om_v = 0.
         END IF
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
         cosm%itk = itk_CAMB
      !ELSE IF (icosmo == 63 .OR. icosmo == 64 .OR. icosmo == 65) THEN
      ELSE IF (is_in_array(icosmo, [63, 64, 65, 79, 80, 81])) THEN
         ! 63 - Planck 2015 (BAHAMAS; Table 1 of 1712.02411)
         ! 64 - Planck 2015 but with 10^7.6 AGN temperature
         ! 65 - Planck 2015 but with 10^8.0 AGN temperature
         ! 79 - Planck 2015 but with 0.12eV neutrinos
         ! 80 - Planck 2015 but with 0.12eV neutrinos
         ! 81 - Planck 2015 but with 0.12eV neutrinos
         ! Note well that this has a neutrino mass
         cosm%itk = itk_CAMB
         cosm%m_nu = 0.06
         cosm%h = 0.6787
         cosm%Om_b = 0.0482
         cosm%Om_m = 0.2571+cosm%Om_b
         cosm%Om_v = 1.-cosm%Om_m
         cosm%ns =  0.9701
         cosm%sig8 = 0.8085
         IF (icosmo == 64) cosm%Theat = 10**7.6
         IF (icosmo == 65) cosm%Theat = 10**8.0
         IF (icosmo == 79) cosm%m_nu = 0.12
         IF (icosmo == 80) cosm%m_nu = 0.24
         IF (icosmo == 81) cosm%m_nu = 0.48
      ELSE IF (icosmo >= 100 .AND. icosmo <= 137) THEN
         ! Mira Titan nodes
         CALL Mira_Titan_node_cosmology(icosmo-100, cosm)
         cosm%itk = itk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo >= 200 .AND. icosmo <= 237) THEN
         ! Franken Emu nodes (which are the same as Franken Emu nodes)
         CALL Franken_Emu_node_cosmology(icosmo-200, cosm)
         cosm%itk = itk_CAMB ! Set to CAMB linear power
      ELSE IF (icosmo >= 300 .AND. icosmo <= 337) THEN
         ! Cosmic Emu nodes (which are the same as Franken Emu nodes)
         CALL Cosmic_Emu_node_cosmology(icosmo-300, cosm)
         cosm%itk = itk_CAMB ! Set to CAMB linear power
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
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Xs, f1, f2
      REAL :: rho_g, Om_g_h2, rho_crit, f_nu_rad
      REAL, PARAMETER :: small = small_curve ! Some small number for writing curvature things

      ! Set all 'has/is' to false
      cosm%has_distance = .FALSE.
      cosm%has_growth = .FALSE.
      cosm%has_sigma = .FALSE.
      cosm%has_spherical = .FALSE.
      cosm%has_power = .FALSE.
      cosm%has_Xde = .FALSE.
      cosm%is_init = .FALSE.
      cosm%is_normalised = .FALSE.
      cosm%has_wiggle = .FALSE.
      cosm%A = 1. ! Overall power normalisaiton, should always be unity

      ! Things to do with finite box
      IF (cosm%box) cosm%kbox = twopi/cosm%Lbox

      IF (cosm%verbose) WRITE (*, *) 'INIT_COSMOLOGY: Calculating derived parameters'

      ! Calculate radiation density (includes photons and neutrinos at recombination)
      rho_g = (4.*SBconst*cosm%T_CMB**4/c_light**3) ! Photon physical density at z=0 from CMB temperature [kg/m^3]
      rho_crit = 3.*H0**2/(8.*pi*bigG) ! Critical density [h^2 kg/m^3]
      Om_g_h2 = rho_g/rho_crit ! Photon cosmological density [h^2]
      cosm%Om_g = Om_g_h2/cosm%h**2 ! Photon density parameter
      cosm%Om_nu_rad = cosm%Om_g*neff_constant*cosm%neff ! Relativisitic neutrino density
      cosm%Om_r = cosm%Om_g+cosm%Om_nu_rad ! Radiation is sum of photon and neutrino densities
      f_nu_rad = cosm%Om_nu_rad/cosm%Om_r ! Fraction of radiation that is in neutrinos (~0.40)

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
      cosm%Om_nu = cosm%m_nu/(neutrino_constant*cosm%h**2)
      cosm%f_nu = cosm%Om_nu/cosm%Om_m
      IF (cosm%Om_nu >= cosm%Om_nu_rad) THEN
         !cosm%a_nu=1./(1900.*cosm%m_nu)
         cosm%a_nu = cosm%Om_nu_rad/cosm%Om_nu
      ELSE
         cosm%a_nu = 1.
      END IF
      cosm%z_nu = redshift_a(cosm%a_nu)
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: Omega_nu (matter):', cosm%Om_nu
         WRITE (*, *) 'INIT_COSMOLOGY: a_nu:', cosm%a_nu
         WRITE (*, *) 'INIT_COSMOLOGY: z_nu:', cosm%z_nu
         WRITE (*, *) 'INIT_COSMOLOGY: f_nu:', cosm%f_nu
      END IF
      IF ((cosm%m_nu .NE. 0) .AND. is_in_array(cosm%itk, [itk_none, itk_DEFW, itk_EH])) THEN
         STOP 'INIT_COSMOLOGY: You cannot use a linear power fitting function for massive neutrino cosmologies'
      END IF

      ! Check neutrino mass fraction is not too high
      IF (cosm%f_nu > 0.5) STOP 'INIT_COSMOLOGY: Error, neutrino mass fraction is too high'
      IF ((cosm%m_nu .NE. 0.) .AND. cosm%a_nu > 0.2) THEN
         WRITE(*, *) 'INIT_COSMOLOGY: Neutrino mass [eV]:', cosm%m_nu
         STOP 'INIT_COSMOLOGY: Error, neutrinos are too light'
      END IF

      ! Decide on scale-dependent growth and cold spectrum
      cosm%scale_dependent_growth = .FALSE.
      cosm%trivial_cold = .TRUE.
      IF (cosm%m_nu .NE. 0.) THEN
         cosm%scale_dependent_growth = .TRUE.
         cosm%trivial_cold = .FALSE.
      END IF

      IF (cosm%verbose .AND. cosm%scale_dependent_growth) WRITE (*, *) 'INIT_COSMOLOGY: Scale-dependent growth'
      IF (cosm%verbose .AND. (.NOT. cosm%trivial_cold))   WRITE (*, *) 'INIT_COSMOLOGY: Non-trivial cold spectrum'

      ! Derived cosmological parameters
      cosm%Om_c = cosm%Om_m-cosm%Om_b-cosm%Om_nu ! Omega_m defined to include CDM, baryons and massive neutrinos
      cosm%Om = cosm%Om_m+cosm%Om_v+cosm%Om_w ! Ignore radiation here
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

      ! Gamma for DEFW
      cosm%Gamma = cosm%Om_m*cosm%h
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: Gamma:', cosm%Gamma
      END IF

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

      ! Ensure deallocate distances
      cosm%has_distance = .FALSE.
      CALL if_allocated_deallocate(cosm%log_p)
      CALL if_allocated_deallocate(cosm%log_a_p)
      cosm%horizon = 0.

      ! Ensure deallocate time
      cosm%has_time = .FALSE.
      CALL if_allocated_deallocate(cosm%log_t)
      CALL if_allocated_deallocate(cosm%log_a_t)
      cosm%age = 0.

      ! Ensure deallocate growth
      cosm%has_growth = .FALSE.
      CALL if_allocated_deallocate(cosm%log_a_growth)
      CALL if_allocated_deallocate(cosm%log_growth)
      CALL if_allocated_deallocate(cosm%growth_rate)
      CALL if_allocated_deallocate(cosm%log_acc_growth)
      cosm%gnorm = 0.

      ! Ensure deallocate sigma
      cosm%has_sigma = .FALSE.
      CALL if_allocated_deallocate(cosm%log_r_sigma)
      CALL if_allocated_deallocate(cosm%log_a_sigma)
      CALL if_allocated_deallocate(cosm%log_sigma)
      CALL if_allocated_deallocate(cosm%log_sigmaa)
      CALL if_allocated_deallocate(cosm%log_sigmac)
      CALL if_allocated_deallocate(cosm%log_sigmaca)

      ! Switch for power?
      IF (cosm%itk == itk_EH .OR. cosm%itk == itk_DEFW .OR. cosm%itk == itk_none) THEN
         ! Default to use internal linear P(k) from Eisenstein & Hu
         cosm%has_power = .FALSE.
      ELSE
         cosm%has_power = .TRUE.
      END IF

      IF (cosm%itk /= itk_external) THEN
         ! Ensure deallocate linear-power tables if not using external tables
         CALL if_allocated_deallocate(cosm%log_k_plin)
         CALL if_allocated_deallocate(cosm%log_a_plin)
         CALL if_allocated_deallocate(cosm%log_plin)
         CALL if_allocated_deallocate(cosm%log_plina)
      ENDIF

      ! TILMAN: Added these
      cosm%nr_sigma = 0
      cosm%na_sigma = 0

      ! TILMAN: Added these
      cosm%amin_sigma = 0.0
      cosm%amax_sigma = 0.0

      ! Ensure delloacte spherical-collapse arrays
      cosm%has_spherical = .FALSE.
      CALL if_allocated_deallocate(cosm%log_a_dcDv)
      CALL if_allocated_deallocate(cosm%dc)
      CALL if_allocated_deallocate(cosm%Dv)

      ! Write finishing message to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_COSMOLOGY: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_cosmology

   SUBROUTINE print_cosmology(cosm)

      ! Prints the cosmological parameters to the screen
      IMPLICIT NONE
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
         WRITE (*, fmt=format) 'COSMOLOGY:', 'z_CMB:', cosm%z_CMB
         WRITE (*, fmt=format) 'COSMOLOGY:', 'n_eff:', cosm%neff
         WRITE (*, fmt=format) 'COSMOLOGY:', 'Y_H:', cosm%YH
         !WRITE(*,fmt=format) 'COSMOLOGY:', 'm_nu 1 [eV]:', cosm%m_nu(1)
         !WRITE(*,fmt=format) 'COSMOLOGY:', 'm_nu 2 [eV]:', cosm%m_nu(2)
         !WRITE(*,fmt=format) 'COSMOLOGY:', 'm_nu 3 [eV]:', cosm%m_nu(3)
         WRITE (*, fmt=format) 'COSMOLOGY:', 'M_nu [eV]:', cosm%m_nu
         WRITE (*, fmt=format) 'COSMOLOGY:', 'log10(T_AGN/K):', log10(cosm%Theat)
         WRITE (*, *) dashes
         !WRITE (*, *) 'COSMOLOGY: Dark energy'
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
            WRITE (*, *) 'COSMOLOGY: Dark energy: CPL'
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
         IF(cosm%itk == itk_none) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: Pure power law'
         ELSE IF(cosm%itk == itk_EH) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: Eisenstein & Hu'
         ELSE IF(cosm%itk == itk_CAMB) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: CAMB'
         ELSE IF(cosm%itk == itk_DEFW) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: DEFW'
         ELSE IF(cosm%itk == 4) THEN
            WRITE(*,*) 'COSMOLOGY: Linear: External'
         ELSE
            STOP 'COSMOLOGY: Error, itk not set properly'
         END IF   
         WRITE (*, fmt=format) 'COSMOLOGY:', 'n_s:', cosm%ns  
         IF(cosm%norm_method == norm_sigma8) THEN
            WRITE (*, *) 'COSMOLOGY: Normalisation: sigma_8'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'sigma_8:', cosm%sig8
         ELSE IF(cosm%norm_method == norm_value) THEN
            WRITE (*, *) 'COSMOLOGY: Normalisation: Power'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'k [h/Mpc]:', cosm%kval
            WRITE (*, fmt=format) 'COSMOLOGY:', 'D^2 [h/Mpc]:', cosm%pval
         ELSE IF(cosm%norm_method == norm_As) THEN
            WRITE (*, *) 'COSMOLOGY: Normalisation: A_s'
            WRITE (*, fmt=format) 'COSMOLOGY:', 'ks [1/Mpc]:', cosm%kpiv
            WRITE (*, fmt=format) 'COSMOLOGY:', 'As:', cosm%As
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
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: icosmo
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose

      CALL assign_cosmology(icosmo, cosm, verbose)
      CALL init_cosmology(cosm)
      CALL print_cosmology(cosm)

   END SUBROUTINE assign_init_cosmology
   
   REAL FUNCTION xi_lin(r, a, flag, cosm)

      ! Computes the 3D linear matter correlation function by integrating over P(k)
      IMPLICIT NONE
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i
      REAL :: k1, k2, xi_bit
      DOUBLE PRECISION :: xi8
      INTEGER, PARAMETER :: min_humps = min_humps_xi
      INTEGER, PARAMETER :: max_humps = max_humps_xi
      INTEGER, PARAMETER :: method = method_xi
      INTEGER, PARAMETER :: iorder = iorder_xi

      STOP 'XI_LIN: This is ridiculuously slow for large R'

      IF (method == 1) THEN

         xi_lin = integrate_cosm(0., 1., xi_integrand_transformed, r, a, flag, cosm, acc_cosm, iorder)

      ELSE IF (method == 2) THEN

         ! Set summation variable to zero
         xi8 = 0.

         ! Loop over humps
         DO i = 0, max_humps

            k1 = i*pi/r
            k2 = (i+1)*pi/r

            xi_bit = integrate_cosm(k1, k2, xi_integrand, r, a, flag, cosm, acc_cosm, iorder)

            xi8 = xi8+xi_bit

            IF (i > min_humps) THEN
               IF (abs(xi_bit/real(xi8)) < acc_cosm) THEN
                  EXIT
               END IF
            END IF

            IF (i == max_humps) THEN
               WRITE (*, *) 'XI_LIN: r [Mpc/h]:', r
               WRITE (*, *) 'XI_LIN: Minimum number of humps:', min_humps
               WRITE (*, *) 'XI_LIN: Maximum number of humps:', max_humps
               WRITE (*, *) 'XI_LIN: Warning, maximum number of humps exceeded'
               STOP
            END IF

         END DO

         xi_lin = real(xi8)

      ELSE

         STOP 'XI_LIN: Error, method specified incorrectly'

      END IF

   END FUNCTION xi_lin

   REAL FUNCTION xi_integrand(k, r, a, flag, cosm)

      ! Integrand for the 3D linear matter correlation function
      USE special_functions
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc] (integrated over)
      REAL, INTENT(IN) :: r ! Separation [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

      IF (k == 0.) THEN
         xi_integrand = 0.
      ELSE
         xi_integrand = sinc(k*r)*p_lin(k, a, flag, cosm)/k
      END IF

   END FUNCTION xi_integrand

   REAL FUNCTION xi_integrand_transformed(t, r, a, flag, cosm)

      ! Integrand for the 3D linear matter correlation function
      ! TODO: Optimise alpha(r)
      USE special_functions
      IMPLICIT NONE
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
         xi_integrand_transformed = sinc(kr)*p_lin(k, a, flag, cosm)*alpha/(t*(1.-t))
      END IF

   END FUNCTION xi_integrand_transformed

   SUBROUTINE normalise_power(cosm)

      ! Normalise the power spectrum using whatever scheme you choose to do this with
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm  

      ! Get the CAMB power if necessary
      IF (cosm%itk == itk_CAMB)     CALL init_CAMB_linear(cosm)
      IF (cosm%itk == itk_external) CALL init_external_linear(cosm)

      ! Change the flag *before* doing the normalisation calculation because it calls power
      cosm%is_normalised = .TRUE.

      ! Normalise the linear spectrum
      IF (cosm%norm_method == norm_sigma8) THEN
         CALL normalise_power_sigma8(cosm)
      ELSE IF (cosm%norm_method == norm_value) THEN
         CALL normalise_power_value(cosm)
      ELSE IF (cosm%itk == itk_CAMB .AND. cosm%norm_method == norm_As) THEN
         ! Do nothing because the CAMB will have done the normalisation correctly
      ELSE IF (cosm%itk == itk_external .AND. cosm%norm_method == norm_none) THEN
         ! No need to do anything
      ELSE
         STOP 'NORMALISE_POWER: Error, normalisation method not specified correctly'
      END IF

      ! If normalisation is not done via sigma8 then calculate sigma8
      IF (cosm%norm_method .NE. norm_sigma8) THEN
         CALL reset_sigma8(cosm)
      END IF

   END SUBROUTINE normalise_power

   SUBROUTINE normalise_power_sigma8(cosm)

      ! Normalising the power spectrum using sigma8
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: sigma8_initial, sigma8_final, kbox_save
      LOGICAL, PARAMETER :: run_CAMB_twice = .FALSE. ! This is almost always a stupid thing to do

      IF (cosm%verbose) WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Normalising power to get correct sigma_8'

      ! Remove the k-cut for the normalisation          
      kbox_save = 0. ! Need to give this a value otherwise get a warning in debug mode
      IF (cosm%box) THEN
         kbox_save = cosm%kbox
         cosm%kbox = 0.
      END IF

      ! Calculate the initial sigma_8 value (will not be correct)
      sigma8_initial = sigma8(cosm)

      IF (cosm%verbose) WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Initial sigma_8:', real(sigma8_initial)

      IF (cosm%itk == itk_EH .OR. cosm%itk == itk_DEFW .OR. cosm%itk == itk_none) THEN
         cosm%A = cosm%sig8/sigma8_initial
         !cosm%A=391.0112 ! Appropriate for sigma_8=0.8 in the boring model (for tests)
      ELSE IF(cosm%itk == itk_CAMB .OR. cosm%itk == itk_external) THEN
         ! TODO: Resetting As is not really necesary and might not be logical to do here
         ! TODO: Need to think about normalisation parameters as primary vs. seconary
         cosm%As = cosm%As*(cosm%sig8/sigma8_initial)**2 
         ! Normalisation
         IF (run_CAMB_twice) THEN
            ! Run again to normalise..., almost certainly foolish
            CALL init_CAMB_linear(cosm)
         ELSE
            ! ... or normalise using sigma8 and rescaling linear power
            IF (cosm%scale_dependent_growth) THEN
               cosm%log_plina = cosm%log_plina+2.*log(cosm%sig8/sigma8_initial)
            ELSE
               cosm%log_plin = cosm%log_plin+2.*log(cosm%sig8/sigma8_initial)
            END IF
         END IF
      ELSE
         STOP 'NORMALISE_POWER_SIGMA8: Error, cannot normalise with this itk'
      END IF

      ! Replace the k-cut if necessary
      IF (cosm%box) THEN
         cosm%kbox = kbox_save
      END IF

      ! Check that the normalisation has been done correctly
      sigma8_final = sigma8(cosm)

      ! Write to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Target sigma_8:', real(cosm%sig8)
         WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Final sigma_8 (calculated):', real(sigma8_final)
         WRITE (*, *) 'NORMALISE_POWER_SIGMA8: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE normalise_power_sigma8

   SUBROUTINE normalise_power_value(cosm)

      ! Normalise the power spectrum by fixing the power at some wavenumber at a=1
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: fac
      REAL, PARAMETER :: a = 1.

      IF (cosm%itk == itk_EH .OR. cosm%itk == itk_DEFW .OR. cosm%itk == itk_none) THEN
         cosm%A = cosm%A*sqrt(cosm%pval/p_lin(cosm%kval, a, flag_power_total, cosm))
      ELSE
         fac = cosm%pval/p_lin(cosm%kval, a, flag_power_total, cosm)
         IF (cosm%scale_dependent_growth) THEN
            cosm%log_plina = cosm%log_plina+log(fac)
         ELSE
            cosm%log_plin = cosm%log_plin+log(fac)
         END IF
      END IF

   END SUBROUTINE normalise_power_value

   SUBROUTINE reset_sigma8(cosm)

      ! Set the value of sigma8 to be correct in the cosmology if it is not used in normalisation
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm

      cosm%sig8 = sigma8(cosm)

   END SUBROUTINE reset_sigma8

   RECURSIVE REAL FUNCTION sigma8(cosm)

      ! Calculate the value of sigma8 from the linear power spectrum
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: R = 8. ! Because we are doing sigma(R = 8 Mpc/h) normalisation
      REAL, PARAMETER :: a = 1. ! Because we are doing simga(R = 8 Mpc/h, a = 1) normalisation

      sigma8 = sigma_integral(R, a, flag_power_total, cosm)

   END FUNCTION sigma8

   REAL FUNCTION comoving_critical_density(a, cosm)

      ! Comoving critical density [(Msun/h) / (Mpc/h)^3]
      ! For LCDM this is constant in the past, increases like a^3 in the future
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      comoving_critical_density = physical_critical_density(a, cosm)*a**3

   END FUNCTION comoving_critical_density

   REAL FUNCTION physical_critical_density(a, cosm)

      ! Physical critical density [(Msun/h) / (Mpc/h)^3]
      ! For LCDM tends to a constant in the future, behaves like a^-3 in the past
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_critical_density = critical_density_cos*Hubble2(a, cosm)

   END FUNCTION physical_critical_density

   REAL FUNCTION comoving_matter_density(cosm)

      ! Comoving matter density [(Msun/h) / (Mpc/h)^3]
      ! Not a function of redshift, constant value throughout time
      IMPLICIT NONE
      TYPE(cosmology), INTENT(IN) :: cosm

      comoving_matter_density = critical_density_cos*cosm%Om_m

   END FUNCTION comoving_matter_density

   REAL FUNCTION physical_matter_density(a, cosm)

      ! Physical matter density [(Msun/h) / (Mpc/h)^3]
      ! Proportional to a^-3 always
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm

      physical_matter_density = comoving_matter_density(cosm)*a**(-3)

   END FUNCTION physical_matter_density

   REAL FUNCTION Hubble2(a, cosm)

      ! Calculates Hubble^2 in units such that H^2(a=1)=1
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'HUBBLE2: Error, cosmology is not initialised'
      Hubble2 = &
         cosm%Om_c*X_c(a)+ &
         cosm%Om_b*X_b(a)+ &
         cosm%Om_g*X_r(a)+ &
         Omega_nu_0(a, cosm)*X_nu(a, cosm)+ &
         cosm%Om_v*X_v(a)+ &
         cosm%Om_w*X_de(a, cosm)+ &
         (1.-cosm%Om)*a**(-2)

   END FUNCTION Hubble2

   REAL FUNCTION Hubble2_norad(a, cosm)

      ! Squared Hubble parameter but ignoring the radiation component
      ! Units such that Hubble^2(a=1)=1 (as long as radiation is not important at a=1)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      Hubble2_norad = Hubble2(a, cosm)-cosm%Om_r*a**(-4)

   END FUNCTION Hubble2_norad

   REAL FUNCTION AH(a, cosm)

      ! Acceleration function: \ddot{a}/a
      ! Some people call this the Raychaudhuri equation
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'AH: Error, cosmology is not initialised'
      AH = &
         cosm%Om_c*(1.+3*w_c)*X_c(a)+ &
         cosm%Om_b*(1.+3*w_b)*X_b(a)+ &
         cosm%Om_g*(1.+3*w_g)*X_g(a)+ &
         Omega_nu_0(a, cosm)*(1.+3*w_nu(a, cosm))*X_nu(a, cosm)+ &
         cosm%Om_v*(1.+3*w_v)*X_v(a)+ &
         cosm%Om_w*(1.+3.*w_de(a, cosm))*X_de(a, cosm)
      AH = -AH/2.

   END FUNCTION AH

   REAL FUNCTION AH_norad(a, cosm)

      ! Acceleration function without the radiation contribution
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      AH_norad = AH(a, cosm)+cosm%Om_r*a**(-4)

   END FUNCTION AH_norad

   REAL FUNCTION Omega_m(a, cosm)

      ! This calculates Omega_m variations with scale factor (note this is not proportional to a^-3 always)
      ! TODO: Maybe retire eventually, since Omega_m is not a real thing due to neutrinos
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_M: Error, cosmology is not initialised'
      Omega_m = cosm%Om_m*X_m(a)/Hubble2(a, cosm)

   END FUNCTION Omega_m

   REAL FUNCTION Omega_cold_norad(a, cosm)

      ! This calculates Omega_cold variations with scale factor, but ignoring photon and neutrino components
      ! This ensures that Omega_m_norad(a->0) -> 1
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_COLD_NORAD: Error, cosmology is not initialised'
      Omega_cold_norad = (cosm%Om_c*X_c(a)+cosm%Om_b*X_b(a))/Hubble2_norad(a, cosm)

   END FUNCTION Omega_cold_norad

   REAL FUNCTION Omega_c(a, cosm)

      ! This calculates Omega_c variations with scale factor (note this is not proportional to a^-3 always)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_C: Error, cosmology is not initialised'
      Omega_c = cosm%Om_c*X_c(a)/Hubble2(a, cosm)

   END FUNCTION Omega_c

   REAL FUNCTION Omega_b(a, cosm)

      ! This calculates Omega_b variations with scale factor (note this is not proportional to a^-3 always)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_B: Error, cosmology is not initialised'
      Omega_b = cosm%Om_b*X_b(a)/Hubble2(a, cosm)

   END FUNCTION Omega_b

   REAL FUNCTION Omega_r(a, cosm)

      ! This calculates Omega_r variations with scale factor (note this is *not* proportional to a^-4 always)
      ! TODO: Maybe retire eventually, since Omega_r is not a real thing due to neutrinos
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_R: Error, cosmology is not initialised'
      Omega_r = cosm%Om_r*X_r(a)/Hubble2(a, cosm)

   END FUNCTION Omega_r

   REAL FUNCTION Omega_g(a, cosm)

      ! This calculates Omega_g variations with scale factor (note this is *not* proportional to a^-4 always)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_G: Error, cosmology is not initialised'
      Omega_g = cosm%Om_g*X_g(a)/Hubble2(a, cosm)

   END FUNCTION Omega_g

   REAL FUNCTION Omega_nu(a, cosm)

      ! This calculates Omega_nu variations with scale factor
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_NU: Error, cosmology is not initialised'
      Omega_nu = Omega_nu_0(a, cosm)*X_nu(a, cosm)/Hubble2(a, cosm)

   END FUNCTION Omega_nu

   REAL FUNCTION Omega_nu_0(a, cosm)

      ! Calcuate the Omega_nu0 value to use in calculations
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (a > cosm%a_nu) THEN
         Omega_nu_0 = cosm%Om_nu
      ELSE
         Omega_nu_0 = cosm%Om_nu_rad
      END IF

   END FUNCTION Omega_nu_0

   REAL FUNCTION Omega_v(a, cosm)

      ! This calculates Omega_v variations with scale factor (note this is not constant)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_V: Error, cosmology is not initialised'
      Omega_v = cosm%Om_v*X_v(a)/Hubble2(a, cosm)

   END FUNCTION Omega_v

   REAL FUNCTION Omega_w(a, cosm)

      ! This calculates Omega_w variations with scale factor
      IMPLICIT NONE
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
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (a > cosm%a_nu) THEN
         w_nu = 0.
      ELSE
         w_nu = 1./3.
      END IF

   END FUNCTION w_nu

   REAL FUNCTION w_de(a, cosm)

      ! Variations of the dark energy equation-of-state parameter w(a)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: p1, p2, p3, p4
      DOUBLE PRECISION :: f1, f2, f3, f4
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
      IMPLICIT NONE
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
      IMPLICIT NONE
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
      IMPLICIT NONE
      REAL, INTENT(IN) :: a

      X_m = a**(-3)

   END FUNCTION X_m

   REAL FUNCTION X_c(a)

      ! Scaling of CDM density
      IMPLICIT NONE
      REAL, INTENT(IN) :: a

      X_c = a**(-3)

   END FUNCTION X_c

   REAL FUNCTION X_b(a)

      ! Scaling of baryon density
      IMPLICIT NONE
      REAL, INTENT(IN) :: a

      X_b = a**(-3)

   END FUNCTION X_b

   REAL FUNCTION X_r(a)

      ! Scaling of radiation density
      ! TODO: Retire because radiation is not really a thing due to neutrinos
      IMPLICIT NONE
      REAL, INTENT(IN) :: a

      X_r = a**(-4)

   END FUNCTION X_r

   REAL FUNCTION X_g(a)

      ! Scaling of photon density
      IMPLICIT NONE
      REAL, INTENT(IN) :: a

      X_g = a**(-4)

   END FUNCTION X_g

   REAL FUNCTION X_nu(a, cosm)

      ! Scaling of neutrino density
      ! TODO: Account for radiation -> matter transition properly
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (a > cosm%a_nu) THEN
         X_nu = a**(-3)
      ELSE
         X_nu = a**(-4)
      END IF

   END FUNCTION X_nu

   REAL FUNCTION X_v(a)

      ! Scaling of vacuum density
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      REAL :: crap

      crap = a

      X_v = 1.

   END FUNCTION X_v

   REAL FUNCTION X_de(a, cosm)

      ! Scaling for dark energy density (i.e., if w=0 x(a)=a^-3, if w=-1 x(a)=const etc.)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      DOUBLE PRECISION :: f1, f2, f3, f4
      !REAL, PARAMETER :: acc = acc_Xde          ! Accuracy for direct integration
      !INTEGER, PARAMETER :: iorder = iorder_Xde ! Order for direct integration

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
         IF(tabulate_Xde) THEN
            IF(.NOT. cosm%has_Xde) CALL init_Xde(cosm)
            X_de = exp(find(log(a), cosm%log_a_Xde, cosm%log_Xde, cosm%n_Xde, &
               iorder_interpolation_Xde, &
               ifind_interpolation_Xde, &
               imeth_interpolation_Xde))
         ELSE
            X_de = Xde_integral(a, cosm)
         END IF
      END IF

   END FUNCTION X_de

   SUBROUTINE init_Xde(cosm)

      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: amin = amin_Xde
      REAL, PARAMETER :: amax = amax_Xde
      INTEGER, PARAMETER :: n = n_Xde
      INTEGER :: i
      REAL, ALLOCATABLE :: a(:), Xde(:)

      CALL fill_array(amin, amax, a, n)
      ALLOCATE(Xde(n))

      WRITE(*, *) 'INIT_XDE: minimum a:', amin
      WRITE(*, *) 'INIT_XDE: maximum a:', amax
      WRITE(*, *) 'INIT_XDE: number of a:', n

      DO i = 1, n
         Xde(i) = Xde_integral(a(i), cosm)
      END DO

      WRITE(*, *) 'INIT_XDE: minimum X_de:', Xde(n)
      WRITE(*, *) 'INIT_XDE: maximum X_de:', Xde(1)

      cosm%log_a_Xde = log(a)
      cosm%log_Xde = log(Xde)
      cosm%n_Xde = n
      cosm%has_Xde = .TRUE.

      WRITE(*, *) 'INIT_XDE: Done'
      WRITE(*, *)

   END SUBROUTINE init_Xde

   REAL FUNCTION Xde_integral(a, cosm)

      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: acc = acc_integration_Xde
      INTEGER, PARAMETER :: iorder = iorder_integration_Xde

      Xde_integral = (a**(-3))*exp(3.*integrate_cosm(a, 1., integrand_de, cosm, acc, iorder))

   END FUNCTION

   REAL FUNCTION integrand_de(a, cosm)

      ! The integrand for the X_de(a) integral
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      integrand_de = w_de(a, cosm)/a

   END FUNCTION integrand_de

   ELEMENTAL REAL FUNCTION scale_factor_z(z)

      ! The scale factor corresponding to redshift 'z'
      IMPLICIT NONE
      REAL, INTENT(IN) :: z

      ! IF (z <= -1.) THEN
      !    WRITE (*, *) 'SCALE_FACTOR_Z: z:', z
      !    STOP 'SCALE_FACTOR_Z: Error, routine called for z<=-1'
      ! END IF

      scale_factor_z = 1./(1.+z)

   END FUNCTION scale_factor_z

   ELEMENTAL REAL FUNCTION redshift_a(a)

      ! The redshift corresponding to scale-factor 'a'
      IMPLICIT NONE
      REAL, INTENT(IN) :: a

      ! IF (a == 0.) THEN
      !    WRITE (*, *) 'REDSHIFT_A: a:', a
      !    STOP 'REDSHIFT_A: Error, routine called with a=0'
      ! END IF

      redshift_a = -1.+1./a

   END FUNCTION redshift_a

   REAL FUNCTION scale_factor_r(r, cosm)

      ! The scale factor corresponding to comoving distance 'r'
      IMPLICIT NONE
      REAL, INTENT(IN) :: r
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: p
      INTEGER, PARAMETER :: iorder = iorder_inversion_distance
      INTEGER, PARAMETER :: ifind = ifind_inversion_distance
      INTEGER, PARAMETER :: imeth = imeth_inversion_distance

      IF (cosm%has_distance .EQV. .FALSE.) CALL init_distance(cosm)
      IF (r == 0.) THEN
         scale_factor_r = 1.
      ELSE
         p = cosm%horizon-r
         scale_factor_r = exp(find(log(p), cosm%log_p, cosm%log_a_p, cosm%n_p, iorder, ifind, imeth))
      END IF

   END FUNCTION scale_factor_r

   REAL FUNCTION redshift_r(r, cosm)

      ! The redshift corresponding to comoving distance 'r'
      IMPLICIT NONE
      REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
      TYPE(cosmology), INTENT(INOUT) :: cosm

      redshift_r = redshift_a(scale_factor_r(r, cosm))

   END FUNCTION redshift_r

   ELEMENTAL REAL FUNCTION f_k(r, cosm)

      ! Curvature function, also comoving angular-diameter distance [Mpc/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (cosm%k == 0.) THEN
         f_k = r
      ELSE IF (cosm%k < 0.) THEN
         f_k = sinh(sqrt(-cosm%k)*r)/sqrt(-cosm%k)
      !ELSE IF (cosm%k > 0.) THEN
      ELSE
         f_k = sin(sqrt(cosm%k)*r)/sqrt(cosm%k)
      ! ELSE
      !    STOP 'F_K: Something went wrong'
      END IF

   END FUNCTION f_k

   ELEMENTAL REAL FUNCTION fdash_k(r, cosm)

      ! Derivative of curvature function df_k(r)/dr
      IMPLICIT NONE
      REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
      TYPE(cosmology), INTENT(IN) :: cosm

      IF (cosm%k == 0.) THEN
         fdash_k = 1.
      ELSE IF (cosm%k < 0.) THEN
         fdash_k = cosh(sqrt(-cosm%k)*r)
      !ELSE IF (cosm%k > 0.) THEN
      ELSE
         fdash_k = cos(sqrt(cosm%k)*r)
      ! ELSE
      !    STOP 'FDASH_K: Something went wrong'
      END IF

   END FUNCTION fdash_k

   REAL FUNCTION comoving_particle_horizon(a, cosm)

      ! The comoving particle horizon [Mpc/h]
      ! This is the furthest distance a particle can have travelled since a=0
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interpolation_distance
      INTEGER, PARAMETER :: ifind = ifind_interpolation_distance
      INTEGER, PARAMETER :: imeth = imeth_interpolation_distance

      IF (cosm%has_distance .EQV. .FALSE.) CALL init_distance(cosm)

      IF (a == 0.) THEN
         comoving_particle_horizon = 0.
      ELSE IF (a > 1.) THEN
         WRITE (*, *) 'COMOVING_PARTICLE_HORIZON: a:', a
         STOP 'COMOVING_PARTICLE_HORIZON: Error, tried to calculate particle horizon in the future'
      ELSE
         comoving_particle_horizon = exp(find(log(a), cosm%log_a_p, cosm%log_p, cosm%n_p, iorder, ifind, imeth))
      END IF

   END FUNCTION comoving_particle_horizon

   REAL FUNCTION physical_particle_horizon(a, cosm)

      ! The physical particle horizon [Mpc/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_particle_horizon = comoving_particle_horizon(a, cosm)*a

   END FUNCTION physical_particle_horizon

   REAL FUNCTION comoving_distance(a, cosm)

      ! The comoving distance [Mpc/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: p

      ! Ensures that init_distance is run and therefore that horizon is calculated
      p = comoving_particle_horizon(a, cosm)

      ! Now calculate the comoving distance
      comoving_distance = cosm%horizon-p

   END FUNCTION comoving_distance

   REAL FUNCTION physical_distance(a, cosm)

      ! The physical distance [Mpc/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_distance = comoving_distance(a, cosm)*a

   END FUNCTION physical_distance

   REAL FUNCTION physical_angular_distance(a, cosm)

      ! The physical angular-diameter distance [Mpc/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      physical_angular_distance = a*comoving_angular_distance(a, cosm)

   END FUNCTION physical_angular_distance

   REAL FUNCTION comoving_angular_distance(a, cosm)

      ! The comoving angular-diameter distance [Mpc/h]
      ! Some people call this the 'effective distance'
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      comoving_angular_distance = f_k(comoving_distance(a, cosm), cosm)

   END FUNCTION comoving_angular_distance

   REAL FUNCTION luminosity_distance(a, cosm)

      ! The luminosity distance [Mpc/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      luminosity_distance = f_k(comoving_distance(a, cosm), cosm)/a

   END FUNCTION luminosity_distance

   SUBROUTINE init_distance(cosm)

      ! Fill up tables of a vs. p(a) (comoving particle horizon)
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: zmin, zmax, amin, amax, a, b, r
      INTEGER :: i
      INTEGER, PARAMETER :: iorder = iorder_integration_distance ! Order for integration

      ! Get from PARAMETERS
      amin = amin_distance
      amax = amax_distance
      cosm%n_p = n_distance

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
      CALL fill_array(log(amin), log(amax), cosm%log_a_p, cosm%n_p)
      IF (ALLOCATED(cosm%log_p)) DEALLOCATE (cosm%log_p)
      ALLOCATE (cosm%log_p(cosm%n_p))

      ! Now do the r(a) calculation
      DO i = 1, cosm%n_p
         !WRITE(*,*) i
         a = exp(cosm%log_a_p(i))
         !WRITE(*,*) i, a
         b = sqrt(a) ! Parameter to make the integrand not diverge for small values (a=b^2)
         !WRITE(*,*) i, a, b
         r = integrate_cosm(0., b, distance_integrand, cosm, acc_cosm, iorder)
         !WRITE(*,*) i, a, b, r
         cosm%log_p(i) = log(r)
         !STOP
      END DO
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_DISTANCE: minimum r [Mpc/h]:', real(exp(cosm%log_p(1)))
         WRITE (*, *) 'INIT_DISTANCE: maximum r [Mpc/h]:', real(exp(cosm%log_p(cosm%n_p)))
      END IF

      ! Find the horizon distance in your cosmology
      ! exp(log) ensures the value is the same as what comes out of the (log) look-up tables
      cosm%horizon = exp(log(integrate_cosm(0., 1., distance_integrand, cosm, acc_cosm, iorder)))
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_DISTANCE: Horizon distance [Mpc/h]:', real(cosm%horizon)
         WRITE (*, *) 'INIT_DISTANCE: Done'
         WRITE (*, *)
      END IF

      cosm%has_distance = .TRUE.

   END SUBROUTINE init_distance

   REAL FUNCTION distance_integrand(b, cosm)

      ! The integrand for the cosmic-distance calculation [Mpc/h]
      ! Cast in terms of b=sqrt(a) to remove integrand divergence at a=0
      ! This means that the integrand is 2b/H(a)a^2, rather than 1/H(a)a^2
      IMPLICIT NONE
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

      ! The age of the universe [Gyr/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interpolation_time
      INTEGER, PARAMETER :: ifind = ifind_interpolation_time
      INTEGER, PARAMETER :: imeth = imeth_interpolation_time

      IF (.NOT. cosm%has_time) CALL init_time(cosm)

      IF (a == 0.) THEN
         cosmic_time = 0.
      ELSE
         cosmic_time = exp(find(log(a), cosm%log_a_t, cosm%log_t, cosm%n_p, iorder, ifind, imeth))
      END IF

   END FUNCTION cosmic_time

   REAL FUNCTION look_back_time(a, cosm)

      ! The time in the past [Gyr/h]
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: t

      t = cosmic_time(a, cosm) ! Ensures that init_time is run and therefore that age is calculated
      look_back_time = cosm%age-t

   END FUNCTION look_back_time

   SUBROUTINE init_time(cosm)

      ! Fill up tables of a vs. r(a) (comoving particle horizon)
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: zmin, zmax, amin, amax, a, t
      INTEGER :: i
      INTEGER, PARAMETER :: iorder = iorder_integration_time

      ! Get from PARAMETERS
      amin = amin_time
      amax = amax_time
      cosm%n_t = n_time

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
      CALL fill_array(log(amin), log(amax), cosm%log_a_t, cosm%n_t)
      IF (ALLOCATED(cosm%log_t)) DEALLOCATE (cosm%log_t)
      ALLOCATE (cosm%log_t(cosm%n_t))

      ! Now do the r(a) calculation
      DO i = 1, cosm%n_t
         a = exp(cosm%log_a_t(i))
         t = integrate_cosm(0., a, time_integrand, cosm, acc_cosm, iorder)
         cosm%log_t(i) = log(t)
      END DO
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_TIME: minimum t [Gyr/h]:', real(exp(cosm%log_t(1)))
         WRITE (*, *) 'INIT_TIME: maximum t [Gyr/h]:', real(exp(cosm%log_t(cosm%n_t)))
      END IF

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
      IMPLICIT NONE
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

   REAL FUNCTION Tk_matter(k, cosm)

      ! Transfer function selection
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(INOUT) :: cosm

      IF (cosm%itk == itk_none) THEN
         Tk_matter = 1.
      ELSE IF (cosm%itk == itk_EH) THEN
         Tk_matter = Tk_EH(k, cosm)
      ELSE IF (cosm%itk == itk_DEFW) THEN
         Tk_matter = Tk_DEFW(k, cosm)
      ELSE
         WRITE (*, *) 'TK: itk:', cosm%itk
         STOP 'TK: Error, itk specified incorrectly'
      END IF

      ! Damp transfer function if considering WDM
      !IF (cosm%warm)      Tk_matter = Tk_matter*Tk_WDM(k, cosm)
      !IF (cosm%bump == 1) Tk_matter = Tk_matter*Tk_bump(k, cosm)
      !IF (cosm%bump == 2) Tk_matter = Tk_matter*Tk_bump_Mexico(k, cosm)
      Tk_matter = Tk_matter*Tk_factor(k, cosm)

   END FUNCTION Tk_matter

   REAL FUNCTION Tk_factor(k, cosm)

      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(INOUT) :: cosm

      ! Damp transfer function if considering WDM
      Tk_factor = 1.
      IF (cosm%warm) Tk_factor = Tk_factor*Tk_WDM(k, cosm)
      IF (cosm%bump == 1) THEN
         Tk_factor = Tk_factor*Tk_bump(k, cosm)
      ELSE IF (cosm%bump == 2) THEN
         Tk_factor = Tk_factor*Tk_bump_Mexico(k, cosm)
      END IF

   END FUNCTION Tk_factor

   REAL FUNCTION Tk_DEFW(k, cosm)

      ! The DEFW transfer function approximation
      ! Relies on the power-spectrum scale parameter Gamma=Omega_m*h
      ! This function was written by John Peacock
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: keff, q, tk4
      DOUBLE PRECISION :: q8, tk8

      keff = 0.172+0.011*log(cosm%Gamma/0.36)*log(cosm%Gamma/0.36)
      q = 1.e-20+k/cosm%Gamma
      q8 = 1.e-20+keff/cosm%Gamma
      tk4 = 1./(1.+(6.4*q+(3.0*q)**1.5+(1.7*q)**2)**1.13)**(1./1.13)
      tk8 = 1./(1.+(6.4*q8+(3.0*q8)**1.5+(1.7*q8)**2)**1.13)**(1./1.13)

      tk_defw = tk4/real(tk8)

   END FUNCTION Tk_DEFW

   REAL FUNCTION Tk_EH(k, cosm)

      USE special_functions

      ! Eisenstein & Hu fitting function (arXiv: 9709112)
      ! JP: the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
      ! JP: remember I use k/h, whereas they use pure k, Om_m is cdm + baryons
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(INOUT) :: cosm
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

      ! Wave-number
      rk = k*h ! Convert to [1/Mpc]

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

      !tb=(tb+ab*fac/(1.+(bb/rk/s)**3))*sin(rk*ss)/rk/ss ! Equation (21)
      tb = (tb+ab*fac/(1.+(bb/rk/s)**3))*sinc(rk*ss) ! Equation (21)

      tk_eh = real((Om_b/Om_m)*tb+(1.-Om_b/Om_m)*tc) ! The weighted mean of baryon and CDM transfer functions

   END FUNCTION Tk_EH

   REAL FUNCTION Tk_WDM(k, cosm)

      ! Warm dark matter 'correction' to the standard transfer function
      ! This version and equation references were taken from arxiv:1605.05973
      ! Originally from Bode et al. (2001; arixv:0010389)
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: alpha, mu

      alpha = 0.074*0.7*cosm%m_wdm**(-1.15) ! alpha from equation (5), units Mpc/h
      mu = 1.12                             ! mu from equation (4), dimensionless

      Tk_wdm = (1.+(alpha*k)**(2.*mu))**(-5./mu) ! Equation (2)

   END FUNCTION Tk_WDM

   REAL FUNCTION Tk_bump(k, cosm)

      ! Put a Gaussian bump in a linear power spectrum
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(INOUT) :: cosm

      Tk_bump = 1.+cosm%A_bump*exp(-(log(k/cosm%k_bump)**2/(2.*cosm%sigma_bump**2)))

   END FUNCTION Tk_bump

   REAL FUNCTION Tk_bump_Mexico(k, cosm)

      ! Put a Gaussian bump in a linear power spectrum
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      TYPE(cosmology), INTENT(INOUT) :: cosm

      Tk_bump_Mexico = sqrt(1.+cosm%A_bump*exp(-(log(k/cosm%k_bump)**2/cosm%sigma_bump**2)))

   END FUNCTION Tk_bump_Mexico

   REAL FUNCTION Tcold(k, a, cosm)

      ! Ratio of transfer function for cold matter relative to all matter
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: a ! Scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      INTEGER, PARAMETER :: method_cold = 3

      IF (cosm%trivial_cold .OR. method_cold == 1 .OR. cosm%f_nu == 0.) THEN
         ! Assuming the cold spectrum is exactly the matter spectrum
         Tcold = 1. 
      ELSE IF (method_cold == 2) THEN
         ! This approximation assumes that the neutrinos are as clustered as the rest of the mass
         ! This is only true on scales greater than the neutrino free-streaming scale
         Tcold = (cosm%Om_c+cosm%Om_b)/cosm%Om_m 
      ELSE IF (method_cold == 3) THEN
         ! Use the Eisenstein and Hu approximation
         Tcold = Tcold_EH(k, a, cosm)
      ELSE IF (method_cold == 4) THEN
         ! Use look-up tables from CAMB transfer functions
         Tcold = find(log(k), cosm%log_k_Tcold, log(a), cosm%log_a_plin, cosm%Tcold, cosm%nk_Tcold, cosm%na_plin, &
            iorder = iorder_interpolation_plin, &
            ifind = ifind_interpolation_plin, &
            iinterp = iinterp_polynomial) ! No Lagrange polynomials 
      ELSE
         STOP 'TCOLD: Error, method not specified correctly'
      END IF

   END FUNCTION Tcold

   ! REAL FUNCTION Tcold_approx(cosm)

   !    ! How the matter power spectrum would be changed if some fraction of the mass is converted to massive neutrinos
   !    ! Approximation for how power is suppressed by massive nu at small scales
   !    ! Calculated assuming perturbation grow from z~1000 and that neutrinos are hot and therefore completely smooth
   !    ! Related to the growth-function approximation: g(a) = a^(1-3f_nu/5)
   !    ! TODO: This is NEVER used and is potentially VERY confusing
   !    ! TODO: This is NOT the transfer function relating the matter spectrum to the cold spectrum
   !    ! TODO: This IS the transfer function relation matter power spectra in different models
   !    ! TODO: Take EXTREME caution here
   !    IMPLICIT NONE
   !    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

   !    Tcold_approx = sqrt(1.-8.*cosm%f_nu)

   ! END FUNCTION Tcold_approx

   REAL FUNCTION Tcold_EH(k, a, cosm)

      ! Calculates the ratio of T(k) for cold vs. all matter
      ! Cold perturbation defined such that 1+delta = rho_cold/rho_matter
      ! Uses approximations in Eisenstein & Hu (1999; arXiv 9710252)
      ! Note that this assumes that there are exactly 3 species of neutrinos
      ! Nnu<=3 of these being massive, and with the mass split evenly between the number of massive species.
      IMPLICIT NONE
      REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: a ! Scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL :: D, Dcb, Dcbnu, pcb, zeq, q, yfs, z
      REAL :: BigT
      INTEGER, PARAMETER :: N_massive_nu = 3

      IF (cosm%f_nu == 0.) THEN

         ! Fix to unity if there are no neutrinos
         Tcold_EH = 1.

      ELSE

         ! Get the redshift
         z = redshift_a(a)

         ! Growth exponent under the assumption that neutrinos are completely unclustered (equation 11)
         pcb = (5.-sqrt(1.+24.*(1.-cosm%f_nu)))/4.

         ! Theta for temperature (BigT=T/2.7 K)
         BigT = cosm%T_CMB/2.7

         ! The matter-radiation equality redshift
         zeq = (2.5e4)*cosm%om_m*(cosm%h**2.)*(BigT**(-4.))

         ! The growth function normalised such that D=(1.+z_eq)/(1+z) at early times (when Omega_m \approx 1)
         ! For my purpose (just the ratio) seems to work better using the EdS growth function result, \propto a .
         ! In any case, can't use grow at the moment because that is normalised by default.
         !D = (1.+zeq)/(1.+z) ! TODO: Could update this
         D = (1.+zeq)*ungrow(a, cosm)

         ! Wave number relative to the horizon scale at equality (equation 5)
         ! Extra factor of h becauase all my k are in units of h/Mpc
         q = k*cosm%h*BigT**2./(cosm%om_m*cosm%h**2.)

         ! Free streaming scale (equation 14)
         ! Note that Eisenstein & Hu (1999) only consider the case of 3 neutrinos
         ! with Nnu of these being massive with the mass split evenly between Nnu species.
         yfs = 17.2*cosm%f_nu*(1.+0.488*cosm%f_nu**(-7./6.))*(real(N_massive_nu)*q/cosm%f_nu)**2.

         ! These are (almost) the scale-dependent growth functions for each component in Eisenstein & Hu (1999)
         ! Some part is missing, but this cancels when they are divided by each other, which is all I need them for.
         ! Equations (12) and (13)
         Dcb = (1.+(D/(1.+yfs))**0.7)**(pcb/0.7)
         Dcbnu = ((1.-cosm%f_nu)**(0.7/pcb)+(D/(1.+yfs))**0.7)**(pcb/0.7)

         ! Finally, the ratio
         Tcold_EH = Dcb*(1.-cosm%f_nu)/Dcbnu

      END IF

   END FUNCTION Tcold_EH

   SUBROUTINE calculate_plin(k, a, Pk, nk, na, cosm)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nk
      INTEGER, INTENT(IN) :: na
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(IN) :: a(na)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      TYPE(cosmology) :: cosm
      INTEGER :: ik, ia
      INTEGER, PARAMETER :: flag = flag_power_total

      ALLOCATE(Pk(nk, na))

      DO ia = 1, na
         DO ik = 1, nk
            Pk(ik, ia) = p_lin(k(ik), a(ia), flag, cosm)
         END DO
      END DO

   END SUBROUTINE calculate_plin

   REAL RECURSIVE FUNCTION p_lin(k, a, flag, cosm)

      ! Linear matter power spectrum
      ! Must be a recursive function because normalise_power calls this function again
      IMPLICIT NONE
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: kmax, pmax
      INTEGER :: nk
      INTEGER, PARAMETER :: iorder = iorder_interpolation_plin ! Order for interpolation (3 - cubic)
      INTEGER, PARAMETER :: ifind = ifind_interpolation_plin   ! Finding scheme in table (3 - Mid-point method)
      INTEGER, PARAMETER :: imeth = imeth_interpolation_plin   ! Method for polynomials (2 - Lagrange polynomials)

      ! Using init_power seems to provide no significant speed improvements to HMx
      !IF(cosm%has_power .EQV. .FALSE.) CALL init_power(cosm)

      ! This line generates a recursion
      IF (.NOT. cosm%is_normalised) THEN
         !cosm%is_normalised = .TRUE. ! TODO: Is this line necessary?         
         CALL normalise_power(cosm)
      END IF

      IF (k <= kmin_plin) THEN
         ! If p_lin happens to be foolishly called for 0 mode
         ! This call should never happen, but may in integrals
         p_lin = 0.
      ELSE IF (k > kmax_plin) THEN
         ! Avoids some issues if p_lin is called for absurdly high k values
         ! For some reason crashes can occur if this is the case
         p_lin = 0.
      ELSE IF (cosm%box .AND. k < cosm%kbox) THEN
         ! If investigating effects caused by a finite box size
         p_lin = 0.
      ELSE
         IF (cosm%has_power) THEN
            IF(plin_extrap) THEN
               nk = cosm%nk_plin
               kmax = exp(cosm%log_k_plin(nk))
            END IF
            IF (cosm%scale_dependent_growth) THEN
               IF(plin_extrap .AND. k > kmax) THEN
                  pmax = exp(find(log(a), cosm%log_a_plin, cosm%log_plina(nk,:), cosm%na_plin, &
                     iorder, ifind, imeth))
                  p_lin = p_lin_extrapolation(k, kmax, pmax, cosm%ns)
               ELSE
                  p_lin = exp(find(log(k), cosm%log_k_plin, log(a), cosm%log_a_plin,&
                     cosm%log_plina, cosm%nk_plin, cosm%na_plin, iorder, ifind, iinterp_polynomial))
               END IF
            ELSE
               IF(plin_extrap .AND. k > kmax) THEN
                  pmax = exp(cosm%log_plin(nk))
                  p_lin = p_lin_extrapolation(k, kmax, pmax, cosm%ns)
               ELSE
                  p_lin = exp(find(log(k), cosm%log_k_plin, cosm%log_plin, cosm%nk_plin, iorder, ifind, imeth))
               END IF
               p_lin = (grow(a, cosm)**2)*p_lin
            END IF
         ELSE
            ! In this case get the power from the transfer function
            p_lin = (cosm%A**2)*(grow(a, cosm)**2)*(Tk_matter(k, cosm)**2)*(k**(cosm%ns+3.))
         END IF
      END IF

      IF (flag == flag_power_cold .OR. flag == flag_power_cold_unorm) THEN
         p_lin = p_lin*Tcold(k, a, cosm)**2
         IF (flag == flag_power_cold_unorm) p_lin = p_lin/(1.-cosm%f_nu)**2
      END IF

   END FUNCTION p_lin

   REAL FUNCTION p_lin_extrapolation(k,kmax,pmax,ns)

      ! Extrapolate linear power at small scales assuming Delta^2(k) goes like ln(k)^2 k^(n-1)
      ! This works really badly if kmax is not high enough; maybe best just to use power-law
      ! TODO: It is really weird that log(k) appears, rather than log(k/kmax), this must be wrong!
      ! TODO: Check ln(k)^2 k^(3+n) at small scales (massive neutrinos; also k has dimensions?)
      IMPLICIT NONE
      REAL, INTENT(IN) :: k    ! Wavenumber [h/Mpc]
      REAL, INTENT(IN) :: kmax ! Maximum wavenumber [h/Mpc]
      REAL, INTENT(IN) :: pmax ! Delta^2(k) at maximum wavenumber
      REAL, INTENT(IN) :: ns   ! Specral index

      p_lin_extrapolation = pmax*((log(k)/log(kmax))**2)*(k/kmax)**(ns-1.)

   END FUNCTION p_lin_extrapolation

   SUBROUTINE init_sigma(cosm)

      ! This fills up tables of r vs. sigma(r) across a range in r
      ! It is used only in look-up for further calculations of sigma(r) and not otherwise
      ! This prevents a large number of calls to the sigma integration functions in HMx
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: R, a, sig
      INTEGER :: i, j, iflag, ntype, flag

      ! Deallocate tables if they are already allocated
      CALL if_allocated_deallocate(cosm%log_r_sigma)
      CALL if_allocated_deallocate(cosm%log_a_sigma)
      CALL if_allocated_deallocate(cosm%log_sigma)
      CALL if_allocated_deallocate(cosm%log_sigmac)
      CALL if_allocated_deallocate(cosm%log_sigmaa)
      CALL if_allocated_deallocate(cosm%log_sigmaca)

      ! Set nr_sigma, na_sigma, amin_sigma, and amax_sigma to defaults if not set by the user
      IF(cosm%nr_sigma == 0) cosm%nr_sigma = nr_sigma
      IF(cosm%na_sigma == 0) cosm%na_sigma = na_sigma
      IF(cosm%amin_sigma == 0.0) cosm%amin_sigma = amin_sigma
      IF(cosm%amax_sigma == 0.0) cosm%amax_sigma = amax_sigma

      ! Write to screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_SIGMA: Filling sigma(R) interpolation tables'
         WRITE (*, *) 'INIT_SIGMA: R minimum [Mpc/h]:', real(rmin_sigma)
         WRITE (*, *) 'INIT_SIGMA: R maximum [Mpc/h]:', real(rmax_sigma)
         WRITE (*, *) 'INIT_SIGMA: Number of points:', cosm%nr_sigma
      END IF

      ! Allocate and fill array of R values
      CALL fill_array(log(rmin_sigma), log(rmax_sigma), cosm%log_r_sigma, cosm%nr_sigma)

      ! Allocate carrays
      IF (cosm%scale_dependent_growth) THEN
         CALL fill_array(log(cosm%amin_sigma), log(cosm%amax_sigma), cosm%log_a_sigma, cosm%na_sigma)
         ALLOCATE (cosm%log_sigmaa(cosm%nr_sigma, cosm%na_sigma))
         IF (sigma_cold_store .NE. 0) ALLOCATE (cosm%log_sigmaca(nr_sigma, cosm%na_sigma))   
      ELSE
         ALLOCATE (cosm%log_sigma(cosm%nr_sigma))
         IF (sigma_cold_store .NE. 0) ALLOCATE (cosm%log_sigmac(cosm%nr_sigma))
      END IF

      ! Do we need to loop over both matter and cold or not?
      IF(cosm%trivial_cold) THEN
         ntype = 1
      ELSE
         ntype = 2
      END IF

      DO iflag = 1, ntype

         IF (iflag == 1) flag = flag_power_total
         IF (iflag == 2) flag = sigma_cold_store

         IF (flag == 0) CYCLE

         ! Do the calculations to fill the look-up tables
         ! Deal with scale-dependent/independent growth separately
         IF (cosm%scale_dependent_growth) THEN

            ! Loop over R and a and calculate sigma(R,a)
            DO j = 1, cosm%na_sigma
               a = exp(cosm%log_a_sigma(j))
               DO i = 1, cosm%nr_sigma
                  R = exp(cosm%log_r_sigma(i))
                  sig = sigma_integral(R, a, flag, cosm)
                  IF (flag == flag_power_total) THEN
                     cosm%log_sigmaa(i, j) = log(sig)
                  ELSE IF (flag == sigma_cold_store) THEN
                     cosm%log_sigmaca(i, j) = log(sig)
                  END IF
               END DO
            END DO

         ELSE

            ! Loop over R values and calculate sigma(R)
            a = 1. ! Fill this table at a=1
            DO i = 1, cosm%nr_sigma
               R = exp(cosm%log_r_sigma(i))
               sig = sigma_integral(R, a, flag, cosm)
               IF (flag == flag_power_total) THEN
                  cosm%log_sigma(i) = log(sig)
               ELSE IF (flag == sigma_cold_store) THEN
                  cosm%log_sigmac(i) = log(sig)
               END IF
            END DO

         END IF

      END DO

      ! If the cold spectrum is trivial then the cold sigma is equal to the matter sigma
      IF ((sigma_cold_store .NE. 0) .AND. cosm%trivial_cold) THEN
         IF (cosm%scale_dependent_growth) THEN
            cosm%log_sigmaca = cosm%log_sigmaa
         ELSE
            cosm%log_sigmac = cosm%log_sigma
         END IF
      END IF

      ! Write useful information to screen
      IF (cosm%verbose) THEN
         IF (cosm%scale_dependent_growth) THEN
            WRITE (*, *) 'INIT_SIGMA: Minimum sigma (a=1):', exp(cosm%log_sigmaa(cosm%nr_sigma,cosm%na_sigma))
            WRITE (*, *) 'INIT_SIGMA: Maximum sigma (a=1):', exp(cosm%log_sigmaa(1,cosm%na_sigma))
         ELSE
            WRITE (*, *) 'INIT_SIGMA: Minimum sigma (a=1):', exp(cosm%log_sigma(cosm%nr_sigma))
            WRITE (*, *) 'INIT_SIGMA: Maximum sigma (a=1):', exp(cosm%log_sigma(1))
         END IF
         WRITE (*, *) 'INIT_SIGMA: Done'
         WRITE (*, *)
      END IF

      ! Change flag so that it is known that the look-up tables are filled
      cosm%has_sigma = .TRUE.

   END SUBROUTINE init_sigma

   REAL FUNCTION find_sigma(R, a, flag, cosm)

      ! Search look-up tables for sigma(R)
      IMPLICIT NONE
      REAL, INTENT(IN) :: R ! Smoothing scale to calculate sigma [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      INTEGER, PARAMETER :: iorder = iorder_interpolation_sigma ! Order for interpolation 
      INTEGER, PARAMETER :: ifind = ifind_interpolation_sigma   ! Finding scheme in table
      INTEGER, PARAMETER :: imeth = imeth_interpolation_sigma   ! Method for polynomials

      IF (cosm%scale_dependent_growth) THEN
         IF (flag == flag_power_total) THEN
            find_sigma = exp(find(log(R), cosm%log_r_sigma, log(a), cosm%log_a_sigma, cosm%log_sigmaa,&
               cosm%nr_sigma, cosm%na_sigma, &
               iorder, ifind, iinterp=iinterp_polynomial)) ! No Lagrange polynomials for 2D interpolation
         ELSE IF (flag == sigma_cold_store) THEN
            find_sigma = exp(find(log(R), cosm%log_r_sigma, log(a), cosm%log_a_sigma, cosm%log_sigmaca,&
               cosm%nr_sigma, cosm%na_sigma, &
               iorder, ifind, iinterp=iinterp_polynomial)) ! No Lagrange polynomials for 2D interpolation
         ELSE
            STOP 'FIND_SIGMA: Error, flag not specified correctly'
         END IF
      ELSE
         IF (flag == flag_power_total) THEN
            find_sigma = exp(find(log(R), cosm%log_r_sigma, cosm%log_sigma, cosm%nr_sigma, iorder, ifind, imeth))
         ELSE IF (flag == sigma_cold_store) THEN
            find_sigma = exp(find(log(R), cosm%log_r_sigma, cosm%log_sigmac, cosm%nr_sigma, iorder, ifind, imeth))
         ELSE
            STOP 'FIND_SIGMA: Error, flag not specified correctly'
         END IF
         find_sigma = grow(a, cosm)*find_sigma
      END IF

   END FUNCTION find_sigma

   RECURSIVE REAL FUNCTION sigma(R, a, flag, cosm)

      ! Finds sigma_all from look-up tables
      IMPLICIT NONE
      REAL, INTENT(IN) :: R ! Smoothing scale to calculate sigma [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

      IF (.NOT. cosm%has_sigma) CALL init_sigma(cosm)
      IF(flag == flag_power_total .OR. flag == sigma_cold_store) THEN
         sigma = find_sigma(R, a, flag, cosm)
      ELSE
         sigma = sigma_integral(R, a, flag, cosm)
      END IF

   END FUNCTION sigma

   RECURSIVE REAL FUNCTION sigma_integral(r, a, flag, cosm)

      ! Calculates sigma(R) by intergration
      IMPLICIT NONE
      REAL, INTENT(IN) :: r ! Smoothing scale to calculate sigma [Mpc/h]
      REAL, INTENT(IN) :: a ! Scale factor
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: acc = acc_sigma ! Because this calculates sigma^2(R), but we only care about sigma(R)
      INTEGER, PARAMETER :: iorder = iorder_sigma
      REAL, PARAMETER :: tmin=0.
      REAL, PARAMETER :: tmax=1.

      sigma_integral = integrate_cosm(tmin, tmax, sigma2_integrand, r, a, flag, cosm, acc, iorder)
      sigma_integral = sqrt(sigma_integral)

   END FUNCTION sigma_integral

   RECURSIVE REAL FUNCTION sigma2_integrand(t, R, a, flag, cosm)

      ! The integrand for the sigma(R) integrals
      USE special_functions
      IMPLICIT NONE
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
         sigma2_integrand = p_lin(k, a, flag, cosm)*(w_hat**2)*alpha/(t*(1.-t))
      END IF

   END FUNCTION sigma2_integrand

   REAL FUNCTION sigmaV(R, a, flag, cosm)

      ! Standard deviation in the displacement field [Mpc/h]
      ! Integrand changed from k [0,inf] to t[0,1] via kR = (1/t-1)**alpha
      IMPLICIT NONE
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
      IMPLICIT NONE
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
         sigmaV2_integrand = (p_lin(k, a, flag, cosm)/k**2)*(w_hat**2)*alpha/(t*(1.-t))
      END IF

   END FUNCTION sigmaV2_integrand

   REAL FUNCTION neff(r, a, flag, cosm)

      ! Effective spectral index on scale R
      ! 3 + neff = -dln(sigma^2)/dlnR and ranges from ~1 on large scales to ~-3 on small scales
      IMPLICIT NONE
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm

      neff = -3.-dsigma(r, a, flag, cosm)

   END FUNCTION neff

   REAL FUNCTION dsigma(r, a, flag, cosm)

      ! dln(sigma^2)/dlnR
      IMPLICIT NONE
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
      IMPLICIT NONE
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
         dsigma_integrand = p_lin(k, a, flag, cosm)*w_hat*w_hat_deriv*kR*alpha/(t*(1.-t))     
      END IF

   END FUNCTION dsigma_integrand

   REAL FUNCTION ncur(r, a, flag, cosm)

      ! Effective spectral index on scale R
      ! 3 + neff = -dln(sigma^2)/dlnR and ranges from ~1 on large scales to ~-3 on small scales
      IMPLICIT NONE
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      TYPE(cosmology), INTENT(INOUT) :: cosm

      ncur = -ddsigma(r, a, flag, cosm)

   END FUNCTION ncur

   REAL FUNCTION ddsigma(r, a, flag, cosm)

      ! Integral for calculating the effective linear P(k) power-law index
      ! This is dlnP(k)/dlnk and ranges from ~1 on large scales to ~-3 on small scales
      IMPLICIT NONE
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
      IMPLICIT NONE
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
         ddsigma_integrand = p_lin(k, a, flag, cosm)*(w_hat*ddw_hat+dw_hat**2)*(kR**2)*alpha/(t*(1.-t))     
      END IF

   END FUNCTION ddsigma_integrand

   REAL RECURSIVE FUNCTION grow(a, cosm)

      ! Scale-independent growth function | normalised g(a=1)=1
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interpolation_grow
      INTEGER, PARAMETER :: ifind = ifind_interpolation_grow
      INTEGER, PARAMETER :: imeth = imeth_interpolation_grow

      IF (cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
      IF (a == 1.) THEN
         grow = 1.
      ELSE
         grow = exp(find(log(a), cosm%log_a_growth, cosm%log_growth, cosm%n_growth, iorder, ifind, imeth))
      END IF

   END FUNCTION grow

   REAL RECURSIVE FUNCTION ungrow(a, cosm)

      ! Growth function normalised such that g(a)~a at early (matter-dominated) times
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      ungrow = cosm%gnorm*grow(a, cosm)

   END FUNCTION ungrow

   REAL RECURSIVE FUNCTION growth_rate(a, cosm)

      ! Growth rate: dln(g) / dln(a)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interpolation_rate
      INTEGER, PARAMETER :: ifind = ifind_interpolation_rate
      INTEGER, PARAMETER :: imeth = imeth_interpolation_rate

      IF (cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
      growth_rate = find(log(a), cosm%log_a_growth, cosm%growth_rate, cosm%n_growth, iorder, ifind, imeth)

   END FUNCTION growth_rate

   REAL RECURSIVE FUNCTION acc_growth(a, cosm)

      ! Accumulated growth function: int_0^a g(a)/a da
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interpolation_agrow
      INTEGER, PARAMETER :: ifind = ifind_interpolation_agrow
      INTEGER, PARAMETER :: imeth = imeth_interpolation_agrow

      IF (cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
      acc_growth = exp(find(log(a), cosm%log_a_growth, cosm%log_acc_growth, cosm%n_growth, iorder, ifind, imeth))

   END FUNCTION acc_growth

   REAL FUNCTION growth_rate_Linder(a, cosm)

      ! Approximation for the growth rate from Linder astro-ph/0507263
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: gam

      IF (cosm%w < -1.) THEN
         gam = 0.55+0.02*(1.+cosm%w)
      ELSE IF (cosm%w > -1) THEN
         gam = 0.55+0.05*(1.+cosm%w)
      ELSE
         gam = 0.55
      END IF

      growth_rate_Linder = Omega_cold_norad(a, cosm)**gam

   END FUNCTION growth_rate_Linder

   REAL FUNCTION grow_Linder_integrand(a, cosm)

      ! Integrand for the approximate growth integral using Linder approximate growth rate
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm

      grow_Linder_integrand = growth_rate_Linder(a, cosm)/a

   END FUNCTION grow_Linder_integrand

   REAL FUNCTION grow_Linder(a, cosm)

      ! Calculate the growth function from the Linder growth rate via integration
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, PARAMETER :: acc = acc_integral_grow
      INTEGER, PARAMETER :: iorder = iorder_integral_grow

      grow_Linder = exp(-integrate_cosm(a, 1., grow_Linder_integrand, cosm, acc, iorder))

   END FUNCTION grow_Linder

   REAL FUNCTION grow_CPT(a, cosm)

      ! Carroll, Press & Turner (1992) approximation to growth function (good to 5%)
      ! https://ui.adsabs.harvard.edu/abs/1992ARA%26A..30..499C/abstract
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Om_mz, Om_vz, Om_m, Om_v

      ! Get all necessary Omega values
      Om_mz = Omega_cold_norad(a, cosm)
      Om_vz = Omega_v(a, cosm)+Omega_w(a, cosm)
      Om_m = cosm%Om_b+cosm%Om_c
      Om_v = cosm%Om_v+cosm%Om_w

      ! Now call CPT twice, second time to normalise it
      grow_CPT = CPT(a, Om_mz, Om_vz)/CPT(1., Om_m, Om_v)

   END FUNCTION grow_CPT

   REAL FUNCTION CPT(a, Om_m, Om_v)

      ! The CPT growth function approximation from 1992
      IMPLICIT NONE
      REAL, INTENT(IN) :: a    ! Scale factor
      REAL, INTENT(IN) :: Om_m ! Matter-density parameter
      REAL, INTENT(IN) :: Om_v ! Vacuum-density parameter

      CPT = a*Om_m/((Om_m**(4./7.))-Om_v+(1.+Om_m/2.)*(1.+Om_v/70.))

   END FUNCTION CPT

   SUBROUTINE init_growth(cosm)

      ! Fills a table of the growth function vs. a
      ! TODO: Figure out why if I set amax=10, rather than amax=1, I start getting weird f(a) around a=0.001
      ! TODO: Should look-up tables be stored log or lin? Currently stored log but not even log spacing in a
      USE calculus_table
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i, na
      REAL :: a
      REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
      REAL :: dinit, vinit
      REAL :: g0, f0, bigG0
      REAL, PARAMETER :: ainit = ainit_growth
      REAL, PARAMETER :: amax = amax_growth
      INTEGER, PARAMETER :: ng = n_growth
      REAL, PARAMETER :: acc_ODE = acc_ODE_growth
      INTEGER, PARAMETER :: imeth_ODE = imeth_ODE_growth
      INTEGER, PARAMETER :: iorder_int = iorder_ODE_interpolation_growth
      INTEGER, PARAMETER :: ifind_int = ifind_ODE_interpolation_growth
      INTEGER, PARAMETER :: imeth_int = imeth_ODE_interpolation_growth
      INTEGER, PARAMETER :: iorder_agrow = iorder_integration_agrow

      !! First do growth factor and growth rate !!

      ! Set the initial conditions to be in the cold matter growing mode
      dinit = ainit**(1.-3.*cosm%f_nu/5.)
      vinit = (1.-3.*cosm%f_nu/5.)*ainit**(-3.*cosm%f_nu/5.)

      ! Write some useful information to the screen
      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_GROWTH: Solving growth equation'
         WRITE (*, *) 'INIT_GROWTH: Minimum scale factor:', ainit
         WRITE (*, *) 'INIT_GROWTH: Maximum scale factor:', amax
         WRITE (*, *) 'INIT_GROWTH: Initial delta:', dinit
         WRITE (*, *) 'INIT_GROWTH: Initial delta derivative:', vinit
         WRITE (*, *) 'INIT_GROWTH: Number of points for look-up tables:', ng
      END IF

      ! Solve the ODE
      CALL ODE_adaptive_cosmology(d_tab, v_tab, 0., a_tab, cosm, ainit, amax, &
         dinit, vinit, ddda, dvda, acc_ODE, imeth_ODE, .FALSE.)
      IF (cosm%verbose) WRITE (*, *) 'INIT_GROWTH: ODE done'
      na = SIZE(a_tab)

      ! Convert dv/da to f = dlng/dlna for later, so v_tab should really be f_tab from now on
      v_tab = v_tab*a_tab/d_tab

      ! Normalise so that g(z=0)=1
      cosm%gnorm = find(1., a_tab, d_tab, na, iorder_int, ifind_int, imeth_int)
      IF (cosm%verbose) WRITE (*, *) 'INIT_GROWTH: unnormalised growth at z=0:', real(cosm%gnorm)
      d_tab = d_tab/cosm%gnorm

      ! Allocate arrays
      IF (ALLOCATED(cosm%log_a_growth)) DEALLOCATE (cosm%log_a_growth)
      IF (ALLOCATED(cosm%log_growth)) DEALLOCATE (cosm%log_growth)
      IF (ALLOCATED(cosm%growth_rate)) DEALLOCATE (cosm%growth_rate)
      IF (ALLOCATED(cosm%log_acc_growth)) DEALLOCATE (cosm%log_acc_growth)
      cosm%n_growth = ng

      ! This downsamples the tables that come out of the ODE solver (which can be a bit long)
      ! Could use some table-interpolation routine here to save time
      ALLOCATE (cosm%log_a_growth(ng))
      ALLOCATE (cosm%log_growth(ng))
      ALLOCATE (cosm%growth_rate(ng))
      DO i = 1, ng
         a = progression(ainit, amax, i, ng) ! TODO: Should this be log?
         cosm%log_a_growth(i) = a
         cosm%log_growth(i) = exp(find(log(a), log(a_tab), log(d_tab), na, iorder_int, ifind_int, imeth_int))
         cosm%growth_rate(i) = find(log(a), log(a_tab), v_tab, na, iorder_int, ifind_int, imeth_int)
      END DO

      !! !!

      !! Table integration to calculate G(a)=int_0^a g(a')/a' da' !!

      ! Allocate array
      ALLOCATE (cosm%log_acc_growth(ng))

      ! Set to zero, because I have an x=x+y thing later on
      cosm%log_acc_growth = 0.

      ! Do the integral up to table position i, which fills the accumulated growth table
      DO i = 1, ng

         ! Do the integral using the arrays
         IF (i > 1) THEN
            cosm%log_acc_growth(i) = integrate_table(cosm%log_a_growth, cosm%gnorm*cosm%log_growth/cosm%log_a_growth, &
               ng, 1, i, iorder_agrow)
         END IF

         ! Then add on the section that is missing from the beginning
         ! NB. g(a=0)/0 = 1, so you just add on a rectangle of height g*a/a=g
         cosm%log_acc_growth(i) = cosm%log_acc_growth(i)+cosm%gnorm*cosm%log_growth(1)

      END DO

      !! !!

      ! Make the some of the tables log for easier interpolation
      cosm%log_a_growth = log(cosm%log_a_growth)
      cosm%log_growth = log(cosm%log_growth)
      cosm%log_acc_growth = log(cosm%log_acc_growth)

      ! Set the flag to true so that this subroutine is only called once
      cosm%has_growth = .TRUE.

      ! Write stuff about growth parameter at a=1 to the screen
      IF (cosm%verbose) THEN
         g0 = grow(1., cosm)
         f0 = growth_rate(1., cosm)
         bigG0 = acc_growth(1., cosm)
         WRITE (*, *) 'INIT_GROWTH: normalised growth at z=0:', g0
         WRITE (*, *) 'INIT_GROWTH: growth rate at z=0:', f0
         WRITE (*, *) 'INIT_GROWTH: integrated growth at z=0:', bigG0
      END IF

      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_GROWTH: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_growth

   REAL FUNCTION ddda(d, v, k, a, cosm)

      ! Needed for growth function solution
      ! This is the dd in \dot{\delta}=dd
      IMPLICIT NONE
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: crap

      ! TODO: prevent compile-time warnings
      crap = d
      crap = k
      crap = cosm%A
      crap = a

      ddda = v

   END FUNCTION ddda

   REAL FUNCTION dvda(d, v, k, a, cosm)

      ! Needed for growth function solution
      ! This is the dv in \ddot{\delta}=dv
      IMPLICIT NONE
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: f1, f2
      REAL :: crap

      ! TODO: prevent compile-time warning
      crap = k

      f1 = 1.5*Omega_cold_norad(a, cosm)*d/(a**2)
      f2 = -(2.+AH_norad(a, cosm)/Hubble2_norad(a, cosm))*(v/a)
      dvda = f1+f2

   END FUNCTION dvda

   REAL FUNCTION dvdanl(d, v, k, a, cosm)

      ! Function used for ODE solver in non-linear growth calculation
      IMPLICIT NONE
      REAL, INTENT(IN) :: d
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: f1, f2, f3
      REAL :: crap

      ! TODO: prevent compile-time warning
      crap = k

      f1 = 1.5*Omega_cold_norad(a, cosm)*d*(1.+d)/(a**2)
      f2 = -(2.+AH_norad(a, cosm)/Hubble2_norad(a, cosm))*(v/a)
      f3 = 4.*(v**2)/(3.*(1.+d))

      dvdanl = f1+f2+f3

   END FUNCTION dvdanl

   REAL FUNCTION dc_NakamuraSuto(a, cosm)

      ! Nakamura & Suto (1997; arXiv:astro-ph/9612074) fitting formula for spherical-collapse in LCDM
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: Om_mz

      Om_mz = Omega_cold_norad(a, cosm)
      dc_NakamuraSuto = dc0*(1.+0.012299*log10(Om_mz))

   END FUNCTION dc_NakamuraSuto

   REAL FUNCTION Dv_BryanNorman(a, cosm)

      ! Bryan & Norman (1998; arXiv:astro-ph/9710107) spherical over-density fitting function
      ! Here overdensity is defined relative to the background matter density, rather than the critical density
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: x, Om_mz

      Om_mz = Omega_cold_norad(a, cosm)
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
      IMPLICIT NONE
      REAL, INTENT(IN) :: a ! scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: lg, bG, Om_m

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

      lg = ungrow(a, cosm)
      bG = acc_growth(a, cosm)
      Om_m = Omega_cold_norad(a, cosm)

      dc_Mead = 1.
      dc_Mead = dc_Mead+f_Mead(lg/a, bG/a, p10, p11, p12, p13)*log10(Om_m)**a1
      dc_Mead = dc_Mead+f_Mead(lg/a, bG/a, p20, p21, p22, p23)
      dc_Mead = dc_Mead*dc0

   END FUNCTION dc_Mead

   REAL FUNCTION Dv_Mead(a, cosm)

      ! Delta_v fitting function from Mead (2017; 1606.05345)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a !scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: lg, bG, Om_m

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

      lg = ungrow(a, cosm)
      bG = acc_growth(a, cosm)
      Om_m = Omega_cold_norad(a, cosm)

      Dv_Mead = 1.
      Dv_Mead = Dv_Mead+f_Mead(lg/a, bG/a, p30, p31, p32, p33)*log10(Om_m)**a3
      Dv_Mead = Dv_Mead+f_Mead(lg/a, bG/a, p40, p41, p42, p43)*log10(Om_m)**a4
      Dv_Mead = Dv_Mead*Dv0

   END FUNCTION Dv_Mead

   REAL FUNCTION f_Mead(x, y, p0, p1, p2, p3)

      ! Equation A3 in Mead (2017)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x, y
      REAL, INTENT(IN) :: p0, p1, p2, p3

      f_Mead = p0+p1*(1.-x)+p2*(1.-x)**2+p3*(1.-y)

   END FUNCTION f_Mead

   REAL FUNCTION dc_spherical(a, cosm)

      ! Get delta_c from a spherical-collapse calculation or look-up table
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interpolation_dc
      INTEGER, PARAMETER :: ifind = ifind_interpolation_dc
      INTEGER, PARAMETER :: imeth = imeth_interpolation_dc

      IF (cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)

      IF (log(a) < cosm%log_a_dcDv(1)) THEN
         dc_spherical = dc0
      ELSE
         dc_spherical = find(log(a), cosm%log_a_dcDv, cosm%dc, cosm%n_dcDv, iorder, ifind, imeth)
      END IF

   END FUNCTION dc_spherical

   REAL FUNCTION Dv_spherical(a, cosm)

      ! Get Delta_v from a spherical-collapse calculation or look-up table
      IMPLICIT NONE
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, PARAMETER :: iorder = iorder_interpolation_Dv
      INTEGER, PARAMETER :: ifind = ifind_interpolation_Dv
      INTEGER, PARAMETER :: imeth = imeth_interpolation_Dv

      IF (cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)

      IF (log(a) < cosm%log_a_dcDv(1)) THEN
         Dv_spherical = Dv0
      ELSE
         Dv_spherical = find(log(a), cosm%log_a_dcDv, cosm%Dv, cosm%n_dcDv, iorder, ifind, imeth)
      END IF

   END FUNCTION Dv_spherical

   SUBROUTINE init_spherical_collapse(cosm)

      ! Initialise the spherical-collapse calculation
      USE table_integer
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: dinit, ainit, vinit, ac, dc
      REAL :: av, a_rmax, d_rmax, Dv, rmax, rv
      REAL, ALLOCATABLE :: d(:), a(:), v(:)
      REAL, ALLOCATABLE :: dnl(:), vnl(:), rnl(:)
      REAL, ALLOCATABLE :: a_coll(:), r_coll(:)
      INTEGER :: i, j, k, k2
      REAL, PARAMETER :: amax = amax_spherical
      INTEGER, PARAMETER :: imeth_ODE = imeth_ODE_spherical
      REAL, PARAMETER :: dmin = dmin_spherical
      REAL, PARAMETER :: dmax = dmax_spherical
      INTEGER, PARAMETER :: n = n_spherical
      INTEGER, PARAMETER :: m = m_spherical

      IF (cosm%verbose) WRITE (*, *) 'SPHERICAL_COLLAPSE: Doing integration'

      ! Allocate arrays
      IF (ALLOCATED(cosm%log_a_dcDv)) DEALLOCATE (cosm%log_a_dcDv)
      IF (ALLOCATED(cosm%dc)) DEALLOCATE (cosm%dc)
      IF (ALLOCATED(cosm%Dv)) DEALLOCATE (cosm%Dv)
      ALLOCATE (cosm%log_a_dcDv(m))
      ALLOCATE (cosm%dc(m))
      ALLOCATE (cosm%Dv(m))
      cosm%log_a_dcDv = 0.
      cosm%dc = 0.
      cosm%Dv = 0.

      IF (cosm%verbose) THEN
         WRITE (*, *) 'SPHERICAL_COLLAPSE: delta min', dmin
         WRITE (*, *) 'SPHERICAL_COLLAPSE: delta max', dmax
         WRITE (*, *) 'SPHERICAL_COLLAPSE: number of collapse points attempted', m
      END IF

      ! BCs for integration. Note ainit=dinit means that collapse should occur around a=1 for dmin
      ! amax should be slightly greater than 1 to ensure at least a few points for a>0.9 (i.e not to miss out a=1)
      ainit = dmin**(1./(1.-3.*cosm%f_nu/5.))
      vinit = (1.-3.*cosm%f_nu/5.)*ainit**(-3.*cosm%f_nu/5.) ! vinit=1 is EdS growing mode solution

      ! Now loop over all initial density fluctuations
      DO j = 1, m

         ! log range of initial delta
         dinit = exp(progression(log(dmin), log(dmax), j, m))

         ! Do both with the same a1 and a2 and using the same number of time steps
         ! This means that arrays a, and anl will be identical, which simplifies calculation
         CALL ODE_spherical(dnl, vnl, 0., a, cosm, ainit, amax, dinit, vinit, ddda, dvdanl, n, imeth_ODE, .TRUE.)
         DEALLOCATE (a)
         CALL ODE_spherical(d, v, 0., a, cosm, ainit, amax, dinit, vinit, ddda, dvda, n, imeth_ODE, .TRUE.)

         ! If this condition is met then collapse occured some time a<amax
         IF (dnl(n) == 0.) THEN

            !! delta_c calcualtion !!

            ALLOCATE (rnl(n))

            rnl = a*(1.+dnl)**(-1./3.)

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
            CALL amputate_array(a, n, 1, k)
            CALL amputate_array(d, n, 1, k)
            CALL amputate_array(dnl, n, 1, k)
            CALL amputate_array(rnl, n, 1, k)

            ! Collapse has occured so use previous a as ac and d as dc
            ac = a(k)
            dc = d(k)

            !! !!

            !! Now to Delta_v calculation !!

            ! Find the a values when the perturbation is maximum size
            a_rmax = maximum(a, rnl, k)

            ! Find the over-density at this point
            d_rmax = exp(find(log(a_rmax), log(a), log(dnl), SIZE(a), &
               iorder=1, ifind=ifind_split, iinterp=iinterp_Lagrange))

            ! Find the maximum radius
            rmax = find(log(a_rmax), log(a), rnl, SIZE(a), &
               iorder=1, ifind=ifind_split, iinterp=iinterp_Lagrange)

            ! The radius of the perturbation when it is virialised is half maximum
            ! This might not be appropriate for LCDM models (or anything with DE)
            rv = rmax/2.

            ! Need to assign new arrays for the collapse branch of r such that it is monotonic
            !k2=int_split(d_rmax,dnl,k)
            k2 = find_table_integer(d_rmax, dnl, k, ifind_split)

            ! Allocate collapse branch arrays
            ALLOCATE (a_coll(k-k2+1), r_coll(k-k2+1))

            ! Fill collapse branch arrays
            DO i = k2, k
               a_coll(i-k2+1) = a(i)
               r_coll(i-k2+1) = rnl(i)
            END DO

            ! Find the scale factor when the perturbation has reached virial radius
            av = exp(find(rv, r_coll, log(a_coll), SIZE(r_coll), &
               iorder=3, ifind=ifind_split, iinterp=iinterp_Lagrange))

            ! Deallocate collapse branch arrays
            DEALLOCATE (a_coll, r_coll)

            ! Spherical model approximation is that perturbation is at virial radius when
            ! 'collapse' is considered to have occured, which has already been calculated
            Dv = exp(find(log(av), log(a), log(dnl), SIZE(a), &
               iorder=1, ifind=ifind_split, iinterp=iinterp_Lagrange))
            Dv = Dv*(ac/av)**3
            Dv = Dv+1.

            !!

            cosm%log_a_dcDv(j) = ac
            cosm%dc(j) = dc
            cosm%Dv(j) = Dv

            DEALLOCATE (rnl)

         END IF

         ! Deallocate arrays ready for next calculation
         DEALLOCATE (d, v, a)
         DEALLOCATE (dnl, vnl)

      END DO

      IF (cosm%verbose) WRITE (*, *) 'SPHERICAL COLLAPSE: calculation complete'

      ! Reverse the arrays so that they run lowest a to highest a
      CALL reverse_array(cosm%log_a_dcDv, m)
      CALL reverse_array(cosm%dc, m)
      CALL reverse_array(cosm%Dv, m)

      IF (cosm%verbose) THEN
         WRITE (*, *) '===================================='
         WRITE (*, *) 'Point  scalefactor  delta_c  Delta_v'
         WRITE (*, *) '===================================='
         DO i = 1, m
            IF (cosm%log_a_dcDv(i) == 0.) EXIT
            WRITE (*, fmt='(I5,F13.4,F9.4,F9.1)') i, cosm%log_a_dcDv(i), cosm%dc(i), cosm%Dv(i)
         END DO
         WRITE (*, *) '===================================='
      END IF

      ! Calculate the maximum sizes for these new arrays
      DO i = 1, m
         IF (cosm%log_a_dcDv(i) == 0.) EXIT
      END DO
      cosm%n_dcDv = i-1

      IF (cosm%verbose) THEN
         WRITE (*, *) 'SPHERICAL_COLLAPSE: number of collapse points:', cosm%n_dcDv
         WRITE (*, *)
      END IF

      ! Remove bits of the array that are unnecessary
      CALL amputate_array(cosm%log_a_dcDv, m, 1, cosm%n_dcDv)
      CALL amputate_array(cosm%dc, m, 1, cosm%n_dcDv)
      CALL amputate_array(cosm%Dv, m, 1, cosm%n_dcDv)

      ! Take a logarithm
      cosm%log_a_dcDv = log(cosm%log_a_dcDv)

      ! Set the flag
      cosm%has_spherical = .TRUE.

   END SUBROUTINE init_spherical_collapse

   SUBROUTINE ODE_spherical(x, v, kk, t, cosm, ti, tf, xi, vi, fx, fv, n, imeth, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and creates arrays of x, v, t values
      ! ODE solver has a fixed number of time steps
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, INTENT(IN) :: kk
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
      DOUBLE PRECISION, ALLOCATABLE :: x8(:), v8(:), t8(:)
      INTEGER :: i

      ! imeth sets ODE solving method
      ! imeth = 1: Crude method
      ! imeth = 2: Mid-point method
      ! imeth = 3: Runge-Kutta

      INTERFACE

         !fx is what x' is equal to
         FUNCTION fx(x, v, k, t, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x, v, k, t
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION fx

         ! fv is what v' is equal to
         FUNCTION fv(x, v, k, t, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x, v, k, t
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION fv

      END INTERFACE

      ! Allocate arrays
      ALLOCATE (x8(n), v8(n), t8(n))

      ! Need to be set to zero for this to work in the spherical-collapse case
      x8 = 0.d0
      v8 = 0.d0
      t8 = 0.d0

      ! xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
      x8(1) = xi
      v8(1) = vi

      ! Fill time array
      IF (ilog) THEN
         CALL fill_array_double(dble(log(ti)), dble(log(tf)), t8, n)
         t8 = exp(t8)
      ELSE
         CALL fill_array_double(dble(ti), dble(tf), t8, n)
      END IF

      DO i = 1, n-1

         CALL ODE_advance_cosmology(x8(i), x8(i+1), v8(i), v8(i+1), t8(i), t8(i+1), fx, fv, imeth, kk, cosm)

         ! Needed to escape from the ODE solver when the perturbation is ~collapsed
         IF (x8(i+1) > dinf_spherical) EXIT

      END DO

      IF (ALLOCATED(x)) DEALLOCATE (x)
      IF (ALLOCATED(v)) DEALLOCATE (v)
      IF (ALLOCATED(t)) DEALLOCATE (t)
      ALLOCATE (x(n), v(n), t(n))
      x = real(x8)
      v = real(v8)
      t = real(t8)

   END SUBROUTINE ODE_spherical

   SUBROUTINE ODE_adaptive_cosmology(x, v, kk, t, cosm, ti, tf, xi, vi, fx, fv, acc, imeth, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
      ! Adaptive, such that time steps are increased until convergence is achieved
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, INTENT(IN) :: kk
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      REAL, INTENT(IN) :: acc ! Desired accuracy
      INTEGER, INTENT(IN) :: imeth
      LOGICAL, INTENT(IN) :: ilog
      DOUBLE PRECISION, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
      INTEGER :: i, j, n, k, np, ifail, kn
      INTEGER, PARAMETER :: jmax = 30
      INTEGER, PARAMETER :: ninit = 100

      ! imeth sets ODE solving method
      ! imeth = 1: Crude method
      ! imeth = 2: Mid-point method
      ! imeth = 3: Runge-Kutta

      INTERFACE

         ! fx is what x' is equal to
         FUNCTION fx(x, v, k, t, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x, v, k, t
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION fx

         ! fv is what v' is equal to
         FUNCTION fv(x, v, k, t, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x, v, k, t
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION fv

      END INTERFACE

      DO j = 1, jmax

         n = 1+ninit*(2**(j-1))

         ALLOCATE (x8(n), v8(n), t8(n))

         x8 = 0.d0
         v8 = 0.d0
         t8 = 0.d0

         ! xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
         x8(1) = xi
         v8(1) = vi

         ! Fill time array
         IF (ilog) THEN
            CALL fill_array_double(dble(log(ti)), dble(log(tf)), t8, n)
            t8 = exp(t8)
         ELSE
            CALL fill_array_double(dble(ti), dble(tf), t8, n)
         END IF

         ifail = 0

         DO i = 1, n-1
            CALL ODE_advance_cosmology(x8(i), x8(i+1), v8(i), v8(i+1), t8(i), t8(i+1), fx, fv, imeth, kk, cosm)
         END DO

         IF (j == 1) ifail = 1

         IF (j .NE. 1) THEN

            np = 1+(n-1)/2

            DO k = 1, 1+(n-1)/2

               kn = 2*k-1

               IF (ifail == 0) THEN

                  IF (xh(k) > acc .AND. x8(kn) > acc .AND. (abs(xh(k)/x8(kn))-1.) > acc) ifail = 1
                  IF (vh(k) > acc .AND. v8(kn) > acc .AND. (abs(vh(k)/v8(kn))-1.) > acc) ifail = 1

                  IF (ifail == 1) THEN
                     DEALLOCATE (xh, th, vh)
                     EXIT
                  END IF

               END IF
            END DO

         END IF

         IF (ifail == 0) THEN
            IF (ALLOCATED(x)) DEALLOCATE (x)
            IF (ALLOCATED(v)) DEALLOCATE (v)
            IF (ALLOCATED(t)) DEALLOCATE (t)
            ALLOCATE (x(n), v(n), t(n))
            x = real(x8)
            v = real(v8)
            t = real(t8)
            EXIT
         END IF

         ALLOCATE (xh(n), th(n), vh(n))
         xh = x8
         vh = v8
         th = t8
         DEALLOCATE (x8, t8, v8)

      END DO

   END SUBROUTINE ODE_adaptive_cosmology

   SUBROUTINE ODE_advance_cosmology(x1, x2, v1, v2, t1, t2, fx, fv, imeth, k, cosm)

      ! Advance the ODE system from t1 to t2
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: x1
      DOUBLE PRECISION, INTENT(OUT) :: x2
      DOUBLE PRECISION, INTENT(IN) :: v1
      DOUBLE PRECISION, INTENT(OUT) :: v2
      DOUBLE PRECISION, INTENT(IN) :: t1
      DOUBLE PRECISION, INTENT(IN) :: t2
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: imeth
      REAL, INTENT(IN) :: k
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: x, v, t, dt
      REAL :: kx1, kx2, kx3, kx4
      REAL :: kv1, kv2, kv3, kv4

      INTERFACE

         ! fx is what x' is equal to
         FUNCTION fx(x, v, k, t, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x, v, k, t
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION fx

         ! fv is what v' is equal to
         FUNCTION fv(x, v, k, t, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x, v, k, t
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION fv

      END INTERFACE

      ! TODO: Is this necessary?
      x = real(x1)
      v = real(v1)
      t = real(t1)

      ! Time step
      dt = real(t2-t1)

      IF (imeth == 1) THEN

         ! Crude method!
         kx1 = dt*fx(x, v, k, t, cosm)
         kv1 = dt*fv(x, v, k, t, cosm)

         x2 = x1+kx1
         v2 = v1+kv1

      ELSE IF (imeth == 2) THEN

         ! Mid-point method!
         kx1 = dt*fx(x, v, k, t, cosm)
         kv1 = dt*fv(x, v, k, t, cosm)
         kx2 = dt*fx(x+kx1/2., v+kv1/2., k, t+dt/2., cosm)
         kv2 = dt*fv(x+kx1/2., v+kv1/2., k, t+dt/2., cosm)

         x2 = x1+kx2
         v2 = v1+kv2

      ELSE IF (imeth == 3) THEN

         ! RK4 (Holy Christ, this is so fast compared to above methods)!
         kx1 = dt*fx(x, v, k, t, cosm)
         kv1 = dt*fv(x, v, k, t, cosm)
         kx2 = dt*fx(x+kx1/2., v+kv1/2., k, t+dt/2., cosm)
         kv2 = dt*fv(x+kx1/2., v+kv1/2., k, t+dt/2., cosm)
         kx3 = dt*fx(x+kx2/2., v+kv2/2., k, t+dt/2., cosm)
         kv3 = dt*fv(x+kx2/2., v+kv2/2., k, t+dt/2., cosm)
         kx4 = dt*fx(x+kx3, v+kv3, k, t+dt, cosm)
         kv4 = dt*fv(x+kx3, v+kv3, k, t+dt, cosm)

         x2 = x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.d0
         v2 = v1+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.d0

      ELSE

         STOP 'ODE_ADVANCE: Error, imeth specified incorrectly'

      END IF

   END SUBROUTINE ODE_advance_cosmology

   REAL RECURSIVE FUNCTION integrate_cosm_1(a, b, f, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      IMPLICIT NONE
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
      DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_cosm_1 = 0.

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

         integrate_cosm_1 = real(sum_new)

      END IF

   END FUNCTION integrate_cosm_1

   REAL RECURSIVE FUNCTION integrate_cosm_2(a, b, f, y, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      IMPLICIT NONE
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
      DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x, y, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x
            REAL, INTENT(IN) :: y
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_cosm_2 = 0.

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
         
         integrate_cosm_2 = real(sum_new)

      END IF

   END FUNCTION integrate_cosm_2

   REAL RECURSIVE FUNCTION integrate_cosm_3(a, b, f, y, z, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      IMPLICIT NONE
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
      DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x, y, z, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x
            REAL, INTENT(IN) :: y
            REAL, INTENT(IN) :: z
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_cosm_3 = 0.

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

         integrate_cosm_3 = real(sum_new)

      END IF

   END FUNCTION integrate_cosm_3

   REAL RECURSIVE FUNCTION integrate_cosm_4(a, b, f, y, z, flag, cosm, acc, iorder)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      IMPLICIT NONE
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
      DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      INTEGER, PARAMETER :: jmin = jmin_integration
      INTEGER, PARAMETER :: jmax = jmax_integration

      INTERFACE
         FUNCTION f(x, y, z, flag, cosm)
            IMPORT :: cosmology
            REAL, INTENT(IN) :: x
            REAL, INTENT(IN) :: y
            REAL, INTENT(IN) :: z
            INTEGER, INTENT(IN) :: flag
            TYPE(cosmology), INTENT(INOUT) :: cosm
         END FUNCTION f
      END INTERFACE

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_cosm_4 = 0.

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

         integrate_cosm_4 = real(sum_new)

      END IF

   END FUNCTION integrate_cosm_4

   SUBROUTINE get_CAMB_power(a, na, k_Pk, Pk, nkPk, k_Tc, Tc, nkTc, non_linear, halofit_version, cosm)

      ! Runs CAMB to get a power spectrum
      ! TODO: Could this be moved to CAMB stuff? Not easily, because it requires cosmology type
      ! TODO: New CAMB cosmology class?
      ! TODO: Could split this up into run_CAMB etc. etc.
      USE CAMB_stuff
      USE string_operations
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: na
      REAL, INTENT(IN) :: a(na)
      REAL, ALLOCATABLE, INTENT(OUT) :: k_Pk(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      INTEGER, INTENT(OUT) :: nkPk
      REAL, ALLOCATABLE, INTENT(OUT) :: k_Tc(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Tc(:, :)
      INTEGER, INTENT(OUT) :: nkTc
      LOGICAL, INTENT(IN) :: non_linear
      INTEGER, INTENT(IN) :: halofit_version
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i, j
      INTEGER :: nTk
      REAL, ALLOCATABLE :: Pk_CAMB(:), Tk_CAMB(:,:)
      REAL :: Om_c, Om_b, Om_nu, h, ombh2, omch2, omnuh2
      CHARACTER(len=256) :: infile
      REAL, PARAMETER :: kmax = kmax_CAMB ! Maximum wavenumber to get the power to
      REAL, PARAMETER :: nmax = nmax_CAMB ! Multiplicative factor to go beyond kmax
      CHARACTER(len=256), PARAMETER :: camb = 'camb'
      CHARACTER(len=256), PARAMETER :: dir = '/Users/Mead/Physics/CAMB_files/tmp/'
      CHARACTER(len=256), PARAMETER :: root = trim(dir)//'temp'
      CHARACTER(len=256), PARAMETER :: matterpower = trim(root)//'_matterpower_'
      CHARACTER(len=256), PARAMETER :: transfer = trim(root)//'_transfer_'
      CHARACTER(len=256), PARAMETER :: params = trim(root)//'_params.ini' 

      IF (cosm%verbose) THEN
         WRITE(*,*) 'GET_CAMB_POWER: Running CAMB'
         WRITE(*,*) 'GET_CAMB_POWER: kpiv [1/Mpc]:', real(cosm%kpiv)
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
         Om_nu = cosm%Om_nu
         h = cosm%h
      END IF

      IF (cosm%Om_v .NE. 0) STOP 'GET_CAMB_POWER: Need to provide Omega_w to CAMB, not Omega_v'

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
      WRITE (7, *) 'dark_energy_model = ppf'
      WRITE (7, *) 'w =', cosm%w
      WRITE (7, *) 'wa =', cosm%wa
      WRITE (7, *) 'cs2_lam = 1'

      ! CMB temperature and Helium fraction
      WRITE (7, *) 'temp_cmb = ', cosm%T_CMB
      WRITE (7, *) 'helium_fraction = ', cosm%YHe

      ! Neutrinos
      !WRITE(7,*) 'massless_neutrinos = ', cosm%neff-real(cosm%N_massive_nu)
      WRITE (7, *) 'share_delta_neff = T'
      IF (omnuh2 == 0.) THEN
         ! Massless neutrinos
         WRITE (7, *) 'massless_neutrinos = ', cosm%neff
         WRITE (7, *) 'nu_mass_eigenstates = 0'
         WRITE (7, *) 'massive_neutrinos = 0'
      ELSE
         ! Three equal-mass neutrinos
         WRITE (7, *) 'massless_neutrinos = ', cosm%neff-3.
         WRITE (7, *) 'nu_mass_eigenstates = 1'
         WRITE (7, *) 'massive_neutrinos = 3'
         WRITE (7, *) 'nu_mass_fractions = 1'
      END IF

      ! Primoridial power spectrum properties
      WRITE (7, *) 'initial_power_num = 1'
      WRITE (7, *) 'scalar_spectral_index(1) =', cosm%ns
      WRITE (7, *) 'scalar_nrun(1) = 0'
      WRITE (7, *) 'scalar_amp(1) =', cosm%As
      WRITE (7, *) 'pivot_scalar =', cosm%kpiv
      WRITE (7, *) 'pivot_tensor =', cosm%kpiv

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

      ! Verbosity
      WRITE (7, *) 'feedback_level = 2'

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

      ! Close file
      CLOSE (7)

      ! Remove old files and run CAMB
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(transfer)//'*')
      CALL EXECUTE_COMMAND_LINE('rm -rf '//trim(matterpower)//'*')
      IF (cosm%verbose) THEN
         CALL EXECUTE_COMMAND_LINE(trim(camb)//' '//trim(params)//' > /dev/null')
      ELSE
         CALL EXECUTE_COMMAND_LINE(trim(camb)//' '//trim(params))
      END IF
      IF (cosm%verbose) WRITE (*, *) 'GET_CAMB_POWER: CAMB run complete'

      ! Loop over redshifts and read CAMB power into arrays
      DO j = 1, na
         infile = number_file(matterpower, j, trim('.dat'))
         CALL read_CAMB_Pk(k_Pk, Pk_CAMB, nkPk, infile)
         IF (j == 1) THEN
            ALLOCATE(Pk(nkPk, na))
         END IF
         Pk(:,j) = Pk_CAMB
      END DO
      DEALLOCATE(Pk_CAMB)

      ! Do pruning
      IF (cosm%verbose) WRITE (*, *) 'GET_CAMB_POWER: nk before pruning:', nkPk
      CALL prune_CAMB(k_Pk, a, Pk, nkPk, na)
      IF (cosm%verbose) THEN
         WRITE (*, *) 'GET_CAMB_POWER: nk after pruning:', nkPk
         WRITE (*, *) 'GET_CAMB_POWER: Getting transfer functions'
      END IF

      ! Loop over redshifts and read CAMB transfer functions into arrays
      DO j = 1, na
         infile = number_file(transfer, j, trim('.dat'))
         CALL read_CAMB_Tk(k_Tc, Tk_CAMB, nkTc, nTk, infile)
         IF (j == 1) THEN
            ALLOCATE(Tc(nkTc, na))
         END IF
         Tk_CAMB(CAMB_column_Tk_CDM,:) = cosm%Om_c*Tk_CAMB(CAMB_column_Tk_CDM,:)
         Tk_CAMB(CAMB_column_Tk_baryon,:) = cosm%Om_b*Tk_CAMB(CAMB_column_Tk_baryon,:)
         Tk_CAMB(CAMB_column_Tk_total,:) = cosm%Om_m*Tk_CAMB(CAMB_column_Tk_total,:)  
         Tc(:,j) = (Tk_CAMB(CAMB_column_Tk_CDM,:)+Tk_CAMB(CAMB_column_Tk_baryon,:))/Tk_CAMB(CAMB_column_Tk_total,:)
      END DO

      IF(cosm%verbose) THEN
         WRITE (*, *) 'GET_CAMB_POWER: Transfer function kmin [h/Mpc]:', k_Tc(1)
         WRITE (*, *) 'GET_CAMB_POWER: Transfer function kmax [h/Mpc]:', k_Tc(nkTc)
         WRITE (*, *) 'GET_CAMB_POWER: Transfer function T(kmin):', Tc(1,1)
         WRITE (*, *) 'GET_CAMB_POWER: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE get_CAMB_power

   SUBROUTINE prune_CAMB(k, a, Pk, nk, na)

      ! Remove some k values from the CAMB calculation of P_lin(k)
      USE array_operations
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(INOUT) :: k(:)
      INTEGER, INTENT(IN) :: na
      REAL, INTENT(IN) :: a(na)
      REAL, ALLOCATABLE, INTENT(INOUT) :: Pk(:,:)
      INTEGER, INTENT(INOUT) :: nk
      INTEGER :: i, naa
      LOGICAL :: kkeep(nk), akeep(na)
      REAL, ALLOCATABLE :: aa(:)
      REAL, PARAMETER :: kmax = kmax_CAMB     ! Maximum wavenumber to get the power to
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

      ! Prune arrays
      !nkk=nk ! Define because this will be overwritten and we do not want to change nk yet
      !CALL apply_mask(k, nkk, kkeep)

      ALLOCATE(aa(na)) ! Define because a is INTENT(IN)
      aa=a   ! Define because a is INTENT(IN)
      naa=na ! Define because na is INTENT(IN)
      CALL apply_mask(k, aa, Pk, nk, naa, kkeep, akeep) ! This will change nk
      IF((naa .NE. na) .OR. (size(aa) .NE. size(a))) STOP 'PRUNE_CAMB: Error, something went wrong'

   END SUBROUTINE prune_CAMB

   SUBROUTINE init_CAMB_linear(cosm)

      ! Initialise the CAMB linear power spectrum calculation
      USE CAMB_stuff
      USE string_operations
      USE array_operations
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL, ALLOCATABLE :: a(:), kPk(:), Pk(:, :), kTc(:), Tc(:, :)
      REAL :: k, fac
      INTEGER :: i, j, na, nPk, nTc, ik
      CHARACTER(len=256), PARAMETER :: camb = 'camb'
      CHARACTER(len=256), PARAMETER :: dir = '/Users/Mead/Physics/CAMB_files/tmp/'
      CHARACTER(len=256), PARAMETER :: root = trim(dir)//'temp'
      CHARACTER(len=256), PARAMETER :: matterpower = trim(root)//'_matterpower_'
      CHARACTER(len=256), PARAMETER :: transfer = trim(root)//'_transfer_'
      CHARACTER(len=256), PARAMETER :: params = trim(root)//'_params.ini'
      LOGICAL, PARAMETER :: non_linear = .FALSE.       ! Should not use non-linear when trying to get linear theory
      INTEGER, PARAMETER :: halofit_version = 0        ! Irrelevant here
      REAL, PARAMETER :: amin = amin_CAMB              ! Minimum scale factor to get from CAMB
      REAL, PARAMETER :: amax = amax_CAMB              ! Maximum scale factor to get from CAMB
      LOGICAL, PARAMETER :: rebin = rebin_CAMB         ! Should we rebin CAMB input P(k)?
      INTEGER, PARAMETER :: iorder = iorder_rebin_CAMB ! Order for interpolation if rebinning
      INTEGER, PARAMETER :: ifind = ifind_rebin_CAMB   ! Finding scheme for interpolation if rebinning
      INTEGER, PARAMETER :: imeth = imeth_rebin_CAMB   ! Method for interpolation if rebinning

      IF (cosm%verbose) THEN
         WRITE(*,*) 'INIT_CAMB_LINEAR: Getting linear power from CAMB'
      END IF

      ! For scale-dependent growth fill the scale-factor array first
      IF (cosm%scale_dependent_growth) THEN
         na = na_CAMB
         CALL fill_array(log(amin), log(amax), cosm%log_a_plin, na)          
      ELSE
         ! For non-scale-dependent growth need a size-one array of a(1)=1.
         na = 1
      END IF

      ! Fill arrays
      cosm%na_plin = na
      ALLOCATE(a(na))
      IF (cosm%scale_dependent_growth) THEN
         a = exp(cosm%log_a_plin)
      ELSE
         a = 1.
      END IF

      ! Get the CAMB P(k) at a series of 'a' values
      CALL get_CAMB_power(a, na, kPk, Pk, nPk, kTc, Tc, nTc, non_linear, halofit_version, cosm)

      ! Apply non-CAMB transfer functions
      DO ik = 1, nPk
         k = kPk(ik)
         fac = Tk_factor(k, cosm)
         Pk(ik,:) = Pk(ik, :)*fac**2
      END DO
      DO ik = 1, nTc
         k = kTc(ik)
         fac = Tk_factor(k, cosm)
         Tc(ik,:) = Tc(ik, :)*fac
      END DO

      CALL if_allocated_deallocate(cosm%log_k_Tcold)
      CALL if_allocated_deallocate(cosm%Tcold)
      ALLOCATE(cosm%log_k_Tcold(nTc), cosm%Tcold(nTc, na))
      cosm%log_k_Tcold = log(kTc)
      cosm%Tcold = Tc
      cosm%nk_Tcold = nTc
      DEALLOCATE(kTc, Tc)

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

      IF(rebin) THEN

         CALL fill_array(log(kmin_CAMB), log(kmax_CAMB), cosm%log_k_plin, nk_CAMB)

         IF(cosm%scale_dependent_growth) THEN         
            IF(ALLOCATED(cosm%log_plina)) DEALLOCATE(cosm%log_plina)
            ALLOCATE(cosm%log_plina(nk_CAMB, na))
            DO j = 1, na
               DO i= 1, nk_CAMB
                  cosm%log_plina(i, j) = find(cosm%log_k_plin(i), log(kPk), log(Pk(:,j)), nPk, iorder, ifind, imeth)
               END DO
            END DO
         ELSE
            IF(ALLOCATED(cosm%log_plin)) DEALLOCATE(cosm%log_plin)
            ALLOCATE(cosm%log_plin(nk_CAMB))
            DO i = 1, nk_CAMB
               cosm%log_plin(i) = find(cosm%log_k_plin(i), log(kPk), log(Pk(:,1)), nPk, iorder, ifind, imeth)
            END DO
         END IF
         cosm%nk_plin = nk_CAMB

         IF (cosm%verbose) THEN
            WRITE(*,*) 'INIT_CAMB_LINEAR: Rebinning kmin [h/Mpc]:', exp(cosm%log_k_plin(1))
            WRITE(*,*) 'INIT_CAMB_LINEAR: Rebinning kmax [h/Mpc]:', exp(cosm%log_k_plin(cosm%nk_plin))
            WRITE(*,*) 'INIT_CAMB_LINEAR: Number of rebinned points in k:', cosm%nk_plin
         END IF

      ELSE  

         IF (ALLOCATED(cosm%log_k_plin)) DEALLOCATE (cosm%log_k_plin)
         ALLOCATE (cosm%log_k_plin(nPk))
         cosm%nk_plin = nPk
         cosm%log_k_plin = log(kPk)

         IF (cosm%scale_dependent_growth) THEN      
            IF (ALLOCATED(cosm%log_plina))  DEALLOCATE (cosm%log_plina)
            ALLOCATE (cosm%log_plina(nPk, na))
            cosm%log_plina = log(Pk)
         ELSE
            IF (ALLOCATED(cosm%log_plin))   DEALLOCATE (cosm%log_plin)
            ALLOCATE (cosm%log_plin(nPk))
            cosm%log_plin = log(Pk(:,1))
         END IF

      END IF

      DEALLOCATE(kPk, Pk)

      ! Now the linear power arrays are filled
      cosm%has_power = .TRUE.

      IF (cosm%verbose) THEN
         WRITE (*, *) 'INIT_CAMB_LINEAR: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE init_CAMB_linear

   SUBROUTINE init_external_linear(cosm)
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: nk, nk_plin, na, na_plin
      INTEGER :: plina_shape(2)

      IF(ALLOCATED(cosm%log_k_plin)) THEN
         nk = SIZE(cosm%log_k_plin)
      ELSE
         write(*,*) "cosmology%log_k_plin has not been allocated!"
         cosm%status = 1
         RETURN
      ENDIF

      IF(ALLOCATED(cosm%log_a_plin)) THEN
         na = SIZE(cosm%log_a_plin)
      ELSE
         write(*,*) "cosmology%log_a_plin has not been allocated!"
         cosm%status = 1
         RETURN
      ENDIF

      IF(ALLOCATED(cosm%log_plin)) THEN
         nk_plin = SIZE(cosm%log_plin)
      ELSE
         write(*,*) "cosmology%log_plin has not been allocated!"
         cosm%status = 1
         RETURN
      ENDIF

      IF(ALLOCATED(cosm%log_plina)) THEN
         plina_shape = SHAPE(cosm%log_plina)
         na_plin = plina_shape(2)
      ELSE
         write(*,*) "cosmology%log_plina has not been allocated!"
         cosm%status = 1
         RETURN
      ENDIF

      IF(nk /= nk_plin .OR. nk /= cosm%nk_plin) THEN
         write(*,*) "Sizes of cosmology%log_plin, cosmology%log_k_plin, or cosmology%nk_plin are inconsistent:", nk_plin, nk, cosm%nk_plin
         cosm%status = 1
         RETURN
      ENDIF

      IF(na /= cosm%na_plin .OR. na /= na_plin) THEN
         write(*,*) "Sizes of cosmology%log_plina, cosmology%log_a_plin, or cosmology%na_plin are inconsistent:", na_plin, na, cosm%na_plin
         cosm%status = 1
         RETURN
      ENDIF

      cosm%amin_sigma = MINVAL(EXP(cosm%log_a_plin))
      cosm%amax_sigma = MAXVAL(EXP(cosm%log_a_plin))

      cosm%has_power = .TRUE.

   END SUBROUTINE init_external_linear

   SUBROUTINE random_cosmology(cosm)

      ! Generate some random cosmological parameters
      USE random_numbers
      IMPLICIT NONE
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
      !cosm%Om_b = cosm%Om_c/5.

      cosm%Om_m = cosm%Om_c+cosm%Om_b

      ! Enforce flatness
      ! Note - need to have Om_w for dark enegry
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m

      cosm%ns =  random_uniform(n_min, n_max)

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

      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_cosmology(cosm)
      cosm%m_nu = 0. ! Fix massive neutrinos to zero

   END SUBROUTINE random_waCDM_cosmology

   SUBROUTINE random_wCDM_cosmology(cosm)

      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_waCDM_cosmology(cosm)
      cosm%wa = 0.      ! Time-independent w
      cosm%iw = iw_wCDM ! iFix w = 4 for wCDM

   END SUBROUTINE random_wCDM_cosmology

   SUBROUTINE random_LCDM_cosmology(cosm)

      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_wCDM_cosmology(cosm)
      cosm%m_nu = 0.
      cosm%w = -1.
      cosm%iw = iw_LCDM ! Fix iw = 1 for LCDM

   END SUBROUTINE random_LCDM_cosmology

   SUBROUTINE random_nuLCDM_cosmology(cosm)

      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm

      CALL random_cosmology(cosm)
      cosm%w = -1.
      cosm%wa = 0.
      cosm%iw = iw_LCDM ! Fix iw = 1 for LCDM

   END SUBROUTINE random_nuLCDM_cosmology

   SUBROUTINE random_Cosmic_Emu_cosmology(cosm)

      ! Generate some random cosmological parameters
      USE random_numbers
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm

      STOP 'RANDOM_COSMIC_EMU_COSMOLOGY: Need to implement the CMB distance condition on h'
      CALL random_Franken_Emu_cosmology(cosm)

   END SUBROUTINE random_Cosmic_Emu_cosmology

   SUBROUTINE random_Franken_Emu_cosmology(cosm)

      ! Generate some random cosmological parameters within the FrankenEmu cube
      USE random_numbers
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm

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

      cosm%h = random_uniform(h_min, h_max)

      cosm%Om_m = random_uniform(om_m_min, om_m_max)/cosm%h**2

      cosm%Om_b = random_uniform(om_b_min, om_b_max)/cosm%h**2

      ! Enforce flatness
      ! Note - need to have Om_w for dark enegry
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m

      cosm%ns =  random_uniform(n_min, n_max)

      cosm%w = random_uniform(w_min, w_max)

      cosm%sig8 = random_uniform(sig8_min, sig8_max)

      cosm%iw = iw_wCDM ! Constant w models only in Franken Emu

   END SUBROUTINE random_Franken_Emu_cosmology

   SUBROUTINE random_Mira_Titan_cosmology(cosm)

      ! Generate some random cosmological parameters for the Mira Titan hypercube
      USE random_numbers
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: om_m, om_b, om_nu

      REAL, PARAMETER :: om_m_min = 0.120
      REAL, PARAMETER :: om_m_max = 0.155

      REAL, PARAMETER :: om_b_min = 0.0215
      REAL, PARAMETER :: om_b_max = 0.0235

      REAL, PARAMETER :: om_nu_min = 0.00
      REAL, PARAMETER :: om_nu_max = 0.01

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

      cosm%h = random_uniform(h_min, h_max)

      om_m = random_uniform(om_m_min, om_m_max)
      cosm%Om_m = om_m/cosm%h**2

      om_b = random_uniform(om_b_min, om_b_max)
      cosm%Om_b = om_b/cosm%h**2

      om_nu = random_uniform(om_nu_min, om_nu_max)
      cosm%m_nu = neutrino_constant*om_nu

      ! Enforce flatness, ensure Omega_w is used for dark energy, Omega_v = 0
      cosm%Om_v = 0.
      cosm%Om_w = 1.-cosm%Om_m

      cosm%ns =  random_uniform(n_min, n_max)

      cosm%w = random_uniform(w_min, w_max)

      ! Enforce 0.3 <= (-w0-wa)^(1/4) as in Mira Titan paper
      DO
         cosm%wa = random_uniform(wa_min, wa_max)
         IF (0.0081 <= -cosm%w-cosm%wa .AND. 2.769 >= -cosm%w-cosm%wa) EXIT
      END DO

      cosm%sig8 = random_uniform(sig8_min, sig8_max)

      cosm%iw = iw_waCDM ! w(a)CDM models in Mira Titan

   END SUBROUTINE random_Mira_Titan_cosmology

   SUBROUTINE Cosmic_Emu_node_cosmology(node, cosm)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: node
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: om_m, om_b

      ! M000 -> M037 of Cosmic_Emu (Table 1 in 1304.7849)
      IF (node == 0) THEN
         ! M000 (not included in emulator construction)
         om_m = 0.1296
         om_b = 0.0224
         cosm%ns =  0.9700
         cosm%w = 1.000
         cosm%sig8 = 0.8000
         cosm%h = 0.7200
      ELSE IF (node == 1) THEN
         ! M001
         om_m = 0.1539
         om_b = 0.0231
         cosm%ns =  0.9468
         cosm%w = 0.816
         cosm%sig8 = 0.8161
         cosm%h = 0.5977
      ELSE IF (node == 2) THEN
         ! M002
         om_m = 0.1460
         om_b = 0.0227
         cosm%ns =  0.8952
         cosm%w = 0.758
         cosm%sig8 = 0.8548
         cosm%h = 0.5970
      ELSE IF (node == 3) THEN
         ! M003
         om_m = 0.1324
         om_b = 0.0235
         cosm%ns =  0.9984
         cosm%w = 0.874
         cosm%sig8 = 0.8484
         cosm%h = 0.6763
      ELSE IF (node == 4) THEN
         ! M004
         om_m = 0.1381
         om_b = 0.0227
         cosm%ns =  0.9339
         cosm%w = 1.087
         cosm%sig8 = 0.7000
         cosm%h = 0.7204
      ELSE IF (node == 5) THEN
         ! M005
         om_m = 0.1358
         om_b = 0.0216
         cosm%ns =  0.9726
         cosm%w = 1.242
         cosm%sig8 = 0.8226
         cosm%h = 0.7669
      ELSE IF (node == 6) THEN
         ! M006
         om_m = 0.1516
         om_b = 0.0229
         cosm%ns =  0.9145
         cosm%w = 1.223
         cosm%sig8 = 0.6705
         cosm%h = 0.7040
      ELSE IF (node == 7) THEN
         ! M007
         om_m = 0.1268
         om_b = 0.0223
         cosm%ns =  0.9210
         cosm%w = 0.70001 ! Changed to avoid problems
         cosm%sig8 = 0.7474
         cosm%h = 0.6189
      ELSE IF (node == 8) THEN
         ! M008
         om_m = 0.1448
         om_b = 0.0223
         cosm%ns =  0.9855
         cosm%w = 1.203
         cosm%sig8 = 0.8090
         cosm%h = 0.7218
      ELSE IF (node == 9) THEN
         ! M009
         om_m = 0.1392
         om_b = 0.0234
         cosm%ns =  0.9790
         cosm%w = 0.739
         cosm%sig8 = 0.6692
         cosm%h = 0.6127
      ELSE IF (node == 10) THEN
         ! M010
         om_m = 0.1403
         om_b = 0.0218
         cosm%ns =  0.8565
         cosm%w = 0.990
         cosm%sig8 = 0.7556
         cosm%h = 0.6695
      ELSE IF (node == 11) THEN
         ! M011
         om_m = 0.1437
         om_b = 0.0234
         cosm%ns =  0.8823
         cosm%w = 1.126
         cosm%sig8 = 0.7276
         cosm%h = 0.7177
      ELSE IF (node == 12) THEN
         ! M012
         om_m = 0.1223
         om_b = 0.0225
         cosm%ns =  1.0048
         cosm%w = 0.971
         cosm%sig8 = 0.6271
         cosm%h = 0.7396
      ELSE IF (node == 13) THEN
         ! M013
         om_m = 0.1482
         om_b = 0.0221
         cosm%ns =  0.9597
         cosm%w = 0.855
         cosm%sig8 = 0.6508
         cosm%h = 0.6107
      ELSE IF (node == 14) THEN
         ! M014
         om_m = 0.1471
         om_b = 0.0233
         cosm%ns =  1.0306
         cosm%w = 1.010
         cosm%sig8 = 0.7075
         cosm%h = 0.6688
      ELSE IF (node == 15) THEN
         ! M015
         om_m = 0.1415
         om_b = 0.0230
         cosm%ns =  1.0177
         cosm%w = 1.281
         cosm%sig8 = 0.7692
         cosm%h = 0.7737
      ELSE IF (node == 16) THEN
         ! M016
         om_m = 0.1245
         om_b = 0.0218
         cosm%ns =  0.9403
         cosm%w = 1.145
         cosm%sig8 = 0.7437
         cosm%h = 0.7929
      ELSE IF (node == 17) THEN
         ! M017
         om_m = 0.1426
         om_b = 0.0215
         cosm%ns =  0.9274
         cosm%w = 0.893
         cosm%sig8 = 0.6865
         cosm%h = 0.6305
      ELSE IF (node == 18) THEN
         ! M018
         om_m = 0.1313
         om_b = 0.0216
         cosm%ns =  0.8887
         cosm%w = 1.029
         cosm%sig8 = 0.6440
         cosm%h = 0.7136
      ELSE IF (node == 19) THEN
         ! M019
         om_m = 0.1279
         om_b = 0.0232
         cosm%ns =  0.8629
         cosm%w = 1.184
         cosm%sig8 = 0.6159
         cosm%h = 0.8120
      ELSE IF (node == 20) THEN
         ! M020
         om_m = 0.1290
         om_b = 0.0220
         cosm%ns =  1.0242
         cosm%w = 0.797
         cosm%sig8 = 0.7972
         cosm%h = 0.6442
      ELSE IF (node == 21) THEN
         ! M021
         om_m = 0.1335
         om_b = 0.0221
         cosm%ns =  1.0371
         cosm%w = 1.165
         cosm%sig8 = 0.6563
         cosm%h = 0.7601
      ELSE IF (node == 22) THEN
         ! M022
         om_m = 0.1505
         om_b = 0.0225
         cosm%ns =  1.0500
         cosm%w = 1.107
         cosm%sig8 = 0.7678
         cosm%h = 0.6736
      ELSE IF (node == 23) THEN
         ! M023
         om_m = 0.1211
         om_b = 0.0220
         cosm%ns =  0.9016
         cosm%w = 1.261
         cosm%sig8 = 0.6664
         cosm%h = 0.8694
      ELSE IF (node == 24) THEN
         ! M024
         om_m = 0.1302
         om_b = 0.0226
         cosm%ns =  0.9532
         cosm%w = 1.300
         cosm%sig8 = 0.6644
         cosm%h = 0.8380
      ELSE IF (node == 25) THEN
         ! M025
         om_m = 0.1494
         om_b = 0.0217
         cosm%ns =  1.0113
         cosm%w = 0.719
         cosm%sig8 = 0.7398
         cosm%h = 0.5724
      ELSE IF (node == 26) THEN
         ! M026
         om_m = 0.1347
         om_b = 0.0232
         cosm%ns =  0.9081
         cosm%w = 0.952
         cosm%sig8 = 0.7995
         cosm%h = 0.6931
      ELSE IF (node == 27) THEN
         ! M027
         om_m = 0.1369
         om_b = 0.0224
         cosm%ns =  0.8500
         cosm%w = 0.836
         cosm%sig8 = 0.7111
         cosm%h = 0.6387
      ELSE IF (node == 28) THEN
         ! M028
         om_m = 0.1527
         om_b = 0.0222
         cosm%ns =  0.8694
         cosm%w = 0.932
         cosm%sig8 = 0.8068
         cosm%h = 0.6189
      ELSE IF (node == 29) THEN
         ! M029
         om_m = 0.1256
         om_b = 0.0228
         cosm%ns =  1.0435
         cosm%w = 0.913
         cosm%sig8 = 0.7087
         cosm%h = 0.7067
      ELSE IF (node == 30) THEN
         ! M030
         om_m = 0.1234
         om_b = 0.0230
         cosm%ns =  0.8758
         cosm%w = 0.777
         cosm%sig8 = 0.6739
         cosm%h = 0.6626
      ELSE IF (node == 31) THEN
         ! M031
         om_m = 0.1550
         om_b = 0.0219
         cosm%ns =  0.9919
         cosm%w = 1.068
         cosm%sig8 = 0.7041
         cosm%h = 0.6394
      ELSE IF (node == 32) THEN
         ! M032
         om_m = 0.1200
         om_b = 0.0229
         cosm%ns =  0.9661
         cosm%w = 1.048
         cosm%sig8 = 0.7556
         cosm%h = 0.7901
      ELSE IF (node == 33) THEN
         ! M033
         om_m = 0.1399
         om_b = 0.0225
         cosm%ns =  1.0407
         cosm%w = 1.147
         cosm%sig8 = 0.8645
         cosm%h = 0.7286
      ELSE IF (node == 34) THEN
         ! M034
         om_m = 0.1497
         om_b = 0.0227
         cosm%ns =  0.9239
         cosm%w = 1.000
         cosm%sig8 = 0.8734
         cosm%h = 0.6510
      ELSE IF (node == 35) THEN
         ! M035
         om_m = 0.1485
         om_b = 0.0221
         cosm%ns =  0.9604
         cosm%w = 0.853
         cosm%sig8 = 0.8822
         cosm%h = 0.6100
      ELSE IF (node == 36) THEN
         ! M036
         om_m = 0.1216
         om_b = 0.0233
         cosm%ns =  0.9387
         cosm%w = 0.706
         cosm%sig8 = 0.8911
         cosm%h = 0.6421
      ELSE IF (node == 37) THEN
         ! M037
         om_m = 0.1495
         om_b = 0.0228
         cosm%ns =  1.0233
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

   END SUBROUTINE Cosmic_Emu_node_cosmology

   SUBROUTINE Franken_Emu_node_cosmology(node, cosm)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: node
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: om_m, om_b

      IF (node == 23) THEN
         ! M023
         om_m = 0.1211
         om_b = 0.0220
         cosm%ns =  0.9016
         cosm%w = -1.261
         cosm%sig8 = 0.6664
         cosm%h = 0.8500 ! Have to round down from 0.8694 to 0.85 for FrankenEmu shrunken space
         cosm%Om_m = om_m/cosm%h**2
         cosm%Om_b = om_b/cosm%h**2
         cosm%Om_w = 1.-cosm%Om_m
         cosm%Om_v = 0.
         cosm%iw = iw_wCDM
      ELSE
         CALL Cosmic_Emu_node_cosmology(node, cosm)
      END IF

   END SUBROUTINE Franken_Emu_node_cosmology

   SUBROUTINE Mira_Titan_node_cosmology(node, cosm)

      IMPLICIT NONE
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
         cosm%ns =  0.963
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
         cosm%ns =  0.9611
         cosm%w = -0.7000
         cosm%wa = 0.67220
         om_nu = 0.
      ELSE IF (node == 2) THEN
         ! M002
         om_m = 0.1356
         om_b = 0.02328
         cosm%sig8 = 0.8556
         cosm%h = 0.7500
         cosm%ns =  1.0500
         cosm%w = -1.0330
         cosm%wa = 0.91110
         om_nu = 0.
      ELSE IF (node == 3) THEN
         ! M003
         om_m = 0.1550
         om_b = 0.02194
         cosm%sig8 = 0.9000
         cosm%h = 0.7167
         cosm%ns =  0.8944
         cosm%w = -1.1000
         cosm%wa = -0.28330
         om_nu = 0.
      ELSE IF (node == 4) THEN
         ! M004
         om_m = 0.1239
         om_b = 0.02283
         cosm%sig8 = 0.7889
         cosm%h = 0.5833
         cosm%ns =  0.8722
         cosm%w = -1.1670
         cosm%wa = 1.15000
         om_nu = 0.
      ELSE IF (node == 5) THEN
         ! M005
         om_m = 0.1433
         om_b = 0.02350
         cosm%sig8 = 0.7667
         cosm%h = 0.8500
         cosm%ns =  0.9833
         cosm%w = -1.2330
         cosm%wa = -0.04445
         om_nu = 0.
      ELSE IF (node == 6) THEN
         ! M006
         om_m = 0.1317
         om_b = 0.021501 ! Moved off boundary
         cosm%sig8 = 0.8333
         cosm%h = 0.5500
         cosm%ns =  0.9167
         cosm%w = -0.7667
         cosm%wa = 0.19440
         om_nu = 0.
      ELSE IF (node == 7) THEN
         ! M007
         om_m = 0.1511
         om_b = 0.02217
         cosm%sig8 = 0.8111
         cosm%h = 0.8167
         cosm%ns =  1.0280
         cosm%w = -0.8333
         cosm%wa = -1.0000
         om_nu = 0.
      ELSE IF (node == 8) THEN
         ! M008
         om_m = 0.1200
         om_b = 0.02306
         cosm%sig8 = 0.7000
         cosm%h = 0.6833
         cosm%ns =  1.0060
         cosm%w = -0.9000
         cosm%wa = 0.43330
      ELSE IF (node == 9) THEN
         ! M009
         om_m = 0.1394
         om_b = 0.02172
         cosm%sig8 = 0.7444
         cosm%h = 0.6500
         cosm%ns =  0.8500
         cosm%w = -0.9667
         cosm%wa = -0.76110
         om_nu = 0.
      ELSE IF (node == 10) THEN
         ! M010 (last cosmology with massless neutrinos)
         om_m = 0.1278
         om_b = 0.02239
         cosm%sig8 = 0.7222
         cosm%h = 0.7833
         cosm%ns =  0.9389
         cosm%w = -1.3000
         cosm%wa = -0.52220
         om_nu = 0.
      ELSE IF (node == 11) THEN
         ! M011 (first cosmology with massive neutrinos)
         om_m = 0.1227
         om_b = 0.0220
         cosm%sig8 = 0.7151
         cosm%h = 0.5827
         cosm%ns =  0.9357
         cosm%w = -1.0821
         cosm%wa = 1.0646
         om_nu = 0.000345
      ELSE IF (node == 12) THEN
         ! M012
         om_m = 0.1241
         om_b = 0.0224
         cosm%sig8 = 0.7472
         cosm%h = 0.8315
         cosm%ns =  0.8865
         cosm%w = -1.2325
         cosm%wa = -0.7646
         om_nu = 0.001204
      ELSE IF (node == 13) THEN
         ! M013
         om_m = 0.1534
         om_b = 0.0232
         cosm%sig8 = 0.8098
         cosm%h = 0.7398
         cosm%ns =  0.8706
         cosm%w = -1.2993
         cosm%wa = 1.2236
         om_nu = 0.003770
      ELSE IF (node == 14) THEN
         ! M014
         om_m = 0.1215
         om_b = 0.0215
         cosm%sig8 = 0.8742
         cosm%h = 0.5894
         cosm%ns =  1.0151
         cosm%w = -0.7281
         cosm%wa = -0.2088
         om_nu = 0.001752
      ELSE IF (node == 15) THEN
         ! M015
         om_m = 0.1250
         om_b = 0.0224
         cosm%sig8 = 0.8881
         cosm%h = 0.6840
         cosm%ns =  0.8638
         cosm%w = -1.0134
         cosm%wa = 0.0415
         om_nu = 0.002789
      ELSE IF (node == 16) THEN
         ! M016
         om_m = 0.1499
         om_b = 0.0223
         cosm%sig8 = 0.7959
         cosm%h = 0.6452
         cosm%ns =  1.0219
         cosm%w = -1.0139
         cosm%wa = 0.9434
         om_nu = 0.002734
      ELSE IF (node == 17) THEN
         ! M017
         om_m = 0.1206
         om_b = 0.0215
         cosm%sig8 = 0.7332
         cosm%h = 0.7370
         cosm%ns =  1.0377
         cosm%w = -0.9472
         cosm%wa = -0.9897
         om_nu = 0.000168
      ELSE IF (node == 18) THEN
         ! M018
         om_m = 0.1544
         om_b = 0.0217
         cosm%sig8 = 0.7982
         cosm%h = 0.6489
         cosm%ns =  0.9026
         cosm%w = -0.7091
         cosm%wa = 0.6409
         om_nu = 0.006419
      ELSE IF (node == 19) THEN
         ! M019
         om_m = 0.1256
         om_b = 0.0222
         cosm%sig8 = 0.8547
         cosm%h = 0.8251
         cosm%ns =  1.0265
         cosm%w = -0.9813
         cosm%wa = -0.3393
         om_nu = 0.004673
      ELSE IF (node == 20) THEN
         ! M020
         om_m = 0.1514
         om_b = 0.0225
         cosm%sig8 = 0.7561
         cosm%h = 0.6827
         cosm%ns =  0.9913
         cosm%w = -1.0101
         cosm%wa = -0.7778
         om_nu = 0.009777
      ELSE IF (node == 21) THEN
         ! M021
         om_m = 0.1472
         om_b = 0.0221
         cosm%sig8 = 0.8475
         cosm%h = 0.6583
         cosm%ns =  0.9613
         cosm%w = -0.9111
         cosm%wa = -1.5470
         om_nu = 0.000672
      ELSE IF (node == 22) THEN
         ! M022
         om_m = 0.1384
         om_b = 0.0231
         cosm%sig8 = 0.8328
         cosm%h = 0.8234
         cosm%ns =  0.9739
         cosm%w = -0.9312
         cosm%wa = 0.5939
         om_nu = 0.008239
      ELSE IF (node == 23) THEN
         ! M023
         om_m = 0.1334
         om_b = 0.0225
         cosm%sig8 = 0.7113
         cosm%h = 0.7352
         cosm%ns =  0.9851
         cosm%w = -0.8971
         cosm%wa = 0.3247
         om_nu = 0.003733
      ELSE IF (node == 24) THEN
         ! M024
         om_m = 0.1508
         om_b = 0.0229
         cosm%sig8 = 0.7002
         cosm%h = 0.7935
         cosm%ns =  0.8685
         cosm%w = -1.0322
         cosm%wa = 1.0220
         om_nu = 0.003063
      ELSE IF (node == 25) THEN
         ! M025
         om_m = 0.1203
         om_b = 0.0230
         cosm%sig8 = 0.8773
         cosm%h = 0.6240
         cosm%ns =  0.9279
         cosm%w = -0.8282
         cosm%wa = -1.5005
         om_nu = 0.007024
      ELSE IF (node == 26) THEN
         ! M026
         om_m = 0.1224
         om_b = 0.0222
         cosm%sig8 = 0.7785
         cosm%h = 0.7377
         cosm%ns =  0.8618
         cosm%w = -0.7463
         cosm%wa = 0.3647
         om_nu = 0.002082
      ELSE IF (node == 27) THEN
         ! M027
         om_m = 0.1229
         om_b = 0.0234
         cosm%sig8 = 0.8976
         cosm%h = 0.8222
         cosm%ns =  0.9698
         cosm%w = -1.0853
         cosm%wa = 0.8683
         om_nu = 0.002902
      ELSE IF (node == 28) THEN
         ! M028
         om_m = 0.1229
         om_b = 0.0231
         cosm%sig8 = 0.8257
         cosm%h = 0.6109
         cosm%ns =  0.9885
         cosm%w = -0.9311
         cosm%wa = 0.8693
         om_nu = 0.009086
      ELSE IF (node == 29) THEN
         ! M029
         om_m = 0.1274
         om_b = 0.0228
         cosm%sig8 = 0.8999
         cosm%h = 0.8259
         cosm%ns =  0.8505
         cosm%w = -0.7805
         cosm%wa = 0.5688
         om_nu = 0.006588
      ELSE IF (node == 30) THEN
         ! M030
         om_m = 0.1404
         om_b = 0.0222
         cosm%sig8 = 0.8232
         cosm%h = 0.6852
         cosm%ns =  0.8679
         cosm%w = -0.8594
         cosm%wa = -0.4637
         om_nu = 0.008126
      ELSE IF (node == 31) THEN
         ! M031
         om_m = 0.1386
         om_b = 0.0229
         cosm%sig8 = 0.7693
         cosm%h = 0.6684
         cosm%ns =  1.0478
         cosm%w = -1.2670
         cosm%wa = 1.2536
         om_nu = 0.006502
      ELSE IF (node == 32) THEN
         ! M032
         om_m = 0.1369
         om_b = 0.021501 ! Moved off boundary
         cosm%sig8 = 0.8812
         cosm%h = 0.8019
         cosm%ns =  1.0005
         cosm%w = -0.7282
         cosm%wa = -1.6927
         om_nu = 0.000905
      ELSE IF (node == 33) THEN
         ! M033
         om_m = 0.1286
         om_b = 0.0230
         cosm%sig8 = 0.7005
         cosm%h = 0.6752
         cosm%ns =  1.0492
         cosm%w = -0.7119
         cosm%wa = -0.8184
         om_nu = 0.007968
      ELSE IF (node == 34) THEN
         ! M034
         om_m = 0.1354
         om_b = 0.0216
         cosm%sig8 = 0.7018
         cosm%h = 0.5970
         cosm%ns =  0.8791
         cosm%w = -0.8252
         cosm%wa = -1.1148
         om_nu = 0.003602
      ELSE IF (node == 35) THEN
         ! M035
         om_m = 0.1359
         om_b = 0.0228
         cosm%sig8 = 0.8210
         cosm%h = 0.6815
         cosm%ns =  0.9872
         cosm%w = -1.1642
         cosm%wa = -0.1801
         om_nu = 0.004440
      ELSE IF (node == 36) THEN
         ! M036
         om_m = 0.1390
         om_b = 0.0220
         cosm%sig8 = 0.8631
         cosm%h = 0.6477
         cosm%ns =  0.8985
         cosm%w = -0.8632
         cosm%wa = 0.8285
         om_nu = 0.001082
      ELSE
         STOP 'MIRA_TITAN_NODE_COSMOLOGY: Error, i specified incorrectly'
      END IF

      cosm%iw = iw_waCDM
      cosm%Om_m = om_m/cosm%h**2
      cosm%Om_b = om_b/cosm%h**2
      !cosm%m_nu=neutrino_constant*om_nu/3. ! Split equally over three neutrinos
      cosm%m_nu = neutrino_constant*om_nu
      cosm%Om_w = 1.-cosm%Om_m
      cosm%Om_v = 0.

   END SUBROUTINE Mira_Titan_node_cosmology

   SUBROUTINE HALOFIT_init(rknl, rneff, rncur, a, cosm, verbose)

      ! Halofit initialisation routine taken from https://www.roe.ac.uk/~jap/haloes/
      ! Calculate the non-linear wavenumber (rknl), effective spectral index (rneff) and curvature (rncur) 
      ! of the power spectrum at the desired redshift, using the method described in Smith et al (2003).
      IMPLICIT NONE
      REAL, INTENT(OUT) :: rknl
      REAL, INTENT(OUT) :: rneff
      REAL, INTENT(OUT) :: rncur
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      LOGICAL, INTENT(IN) :: verbose
      REAL :: xlogr1, xlogr2, rmid, sig, d1, d2, diff
      REAL, PARAMETER :: diff_limit = 0.001

      IF (verbose) WRITE (*, *) 'HALOFIT_INIT: computing effective spectral quantities:'

      xlogr1 = -3.0
      xlogr2 = 3.5
      DO
         rmid = 10**((xlogr2+xlogr1)/2.0)
         CALL wint_HALOFIT(rmid, sig, d1, d2, a, cosm)
         diff = sig-1.0
         IF (abs(diff) <= diff_limit) THEN
            rknl = 1./rmid
            rneff = -3-d1
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

   SUBROUTINE calculate_HALOFIT(k, a, Pk, nk, na, cosm, version)

      IMPLICIT NONE
      INTEGER :: nk ! Number of points in k
      INTEGER :: na ! Number of points in a
      REAL, INTENT(IN) :: k(nk) ! Array of wavenumbers
      REAL, INTENT(IN) :: a(na) ! Scale factor
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :) ! Output power array
      TYPE(cosmology), INTENT(INOUT) :: cosm   ! Cosmology
      INTEGER, OPTIONAL, INTENT(IN) :: version ! HALOFIT version
      REAL :: Plin(nk), Pq(nk), Ph(nk), Pnl(nk)
      INTEGER :: j
      LOGICAL, PARAMETER :: verbose = .FALSE.

      ALLOCATE(Pk(nk, na))

      DO j = 1, na
         CALL calculate_HALOFIT_a(k, a(j), Plin, Pq, Ph, Pnl, nk, cosm, verbose, version)
         Pk(:, j) = Pnl
      END DO

   END SUBROUTINE calculate_HALOFIT

   SUBROUTINE calculate_HALOFIT_a(k, a, Plin, Pq, Ph, Pnl, n, cosm, verbose, ihf)

      ! Get a HALOFIT P(k,a) prediction
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n       ! Size of input k array
      REAL, INTENT(IN) :: k(n)       ! Array of wavenumbers
      REAL, INTENT(IN) :: a          ! Scale factor
      REAL, INTENT(OUT) :: Plin(n)   ! Output array of linear power
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
         Plin(i) = p_lin(k(i), a, flag_power_total, cosm)
         CALL HALOFIT(k(i), rneff, rncur, rknl, Plin(i), Pnl(i), Pq(i), Ph(i), a, cosm, ihf)
      END DO

   END SUBROUTINE calculate_HALOFIT_a

   SUBROUTINE wint_HALOFIT(r, sig, d1, d2, a, cosm)

      ! Halofit window integration routine
      IMPLICIT NONE
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
         d2 = p_lin(rk, a, flag_power_total, cosm)
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

   SUBROUTINE HALOFIT(rk, rn, rncur, rknl, plin, pnl, pq, ph, a, cosm, ihf)

      ! Calculates the HALOFIT power spectrum after rn, rncur and rknl have been pre-calculated
      ! ihf = 1 - Smith et al. 2003 (astro-ph/0207664)
      ! ihf = 2 - Bird et al. 2012
      ! ihf = 3 - Takahashi et al. 2012
      ! ihf = 4 - Takahashi et al. but taken from CAMB with some neutrino stuff
      ! ihf = 5 - Takahashi et al. but taken from CLASS with some neutrino stuff (https://github.com/cmbant/CAMB/issues/44)
      IMPLICIT NONE
      REAL, INTENT(IN) :: rk
      REAL, INTENT(IN) :: rn
      REAL, INTENT(IN) :: rncur
      REAL, INTENT(IN) :: rknl
      REAL, INTENT(IN) :: plin
      REAL, INTENT(OUT) :: pnl
      REAL, INTENT(OUT) :: pq
      REAL, INTENT(OUT) :: ph
      REAL, INTENT(IN) :: a
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, INTENT(IN) :: ihf
      REAL :: gam, aa, bb, cc, mu, nu, alpha, beta, f1, f2, f3, y, Q, fy
      REAL :: om_mz, om_vz, fnu, om_m, wz
      real :: f1a, f2a, f3a, f1b, f2b, f3b, frac

      ! Necessary cosmological parameters
      Om_m = cosm%Om_m         ! TODO: Should neutrinos be included?
      Om_mz = Omega_m(a, cosm) ! TODO: Should neutrinos be included?
      Om_vz = Omega_v(a, cosm)+Omega_w(a, cosm) ! Note this well
      wz = w_de(a, cosm) ! Choice here; do you use w or w(z)? w(z) is better I think; CAMB makes w(z) choice too
      fnu = cosm%Om_nu/cosm%Om_m ! TODO: Does this mean that neutrinos should not be included in Om_m?

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
      ELSE IF (ihf == HALOFIT_Smith) THEN
         ! Smith et al. (2003); the numbers here are EXACTLY those from the online code
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
         gam = 0.8649+0.2989*rn+0.1631*rncur+0.316-0.0765*rn-0.835*rncur ! Bird equation (A5)
         alpha = 1.3884+0.3700*rn-0.1452*rn**2
         beta = 0.8291+0.9854*rn+0.3401*rn**2+fnu*(-6.49+1.44*rn**2) ! Bird equation (A10)
         mu = 10**(-3.5442+0.1908*rn)
         nu = 10**(0.9589+1.2857*rn)
      ELSE IF (ihf == HALOFIT_Bird) THEN
         ! Bird et al. (2012); based off Smith et al. (2003) with numbers taken from the online code
         aa = 10**(1.4861+1.83693*rn+1.67618*rn**2+0.7940*rn**3+0.1670756*rn**4-0.620695*rncur)
         bb = 10**(0.9463+0.9466*rn+0.3084*rn**2-0.940*rncur)
         cc = 10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur)
         gam = 0.86485+0.2989*rn+0.1631*rncur+0.316-0.0765*rn-0.835*rncur ! Bird equation (A5)
         alpha = 1.38848+0.3701*rn-0.1452*rn**2
         beta = 0.8291+0.9854*rn+0.3400*rn**2+fnu*(-6.49+1.44*rn**2) ! Bird equation (A10)
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
         gam = 0.1971-0.0843*rn+0.8460*rncur ! Same as Takahashi
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

      IF (ihf == HALOFIT_Smith .OR. ihf == HALOFIT_Smith_paper) THEN
         ! Smith et al. (2003)
         fy = y/4.+y**2/8.                                    ! Smith (below C2)
         ph = aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam)) ! Smith (C4)
         ph = ph/(1.+mu*y**(-1)+nu*y**(-2))                   ! Smith (C3)
         pq = plin*(1.+plin)**beta/(1.+plin*alpha)*exp(-fy)   ! Smith (C2)
      ELSE IF (ihf == HALOFIT_Bird .OR. ihf == HALOFIT_Bird_paper) THEN
         ! Bird et al. (2012)
         fy = y/4.+y**2/8.
         ph = aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam)) ! Bird equation (A2)
         ph = ph/(1.+mu*y**(-1)+nu*y**(-2))                   ! Bird equation (A1)
         Q = fnu*(2.080-12.4*(Om_m-0.3))/(1.+(1.2e-3)*y**3)   ! Bird equation (A6)
         ph = ph*(1.+Q)                                       ! Bird equation (A7)
         pq = plin*(1.+(26.3*fnu*rk**2.)/(1.+1.5*rk**2))      ! Bird equation (A9)
         pq = plin*(1.+pq)**beta/(1.+pq*alpha)*exp(-fy)       ! Bird equation (A8)
      ELSE IF (ihf == HALOFIT_Takahashi) THEN
         ! Takahashi et al. (2012)
         fy = y/4.+y**2/8.                                    ! Takahashi equation (below A2)
         ph = aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam)) ! Takahashi equation (A3ii)
         ph = ph/(1.+mu*y**(-1)+nu*y**(-2))                   ! Takahashi equation (A3i)
         pq = plin*(1.+plin)**beta/(1.+plin*alpha)*exp(-fy)   ! Takahashi equation (A2)
      ELSE IF (ihf == HALOFIT_CAMB) THEN
         ! Unpublished CAMB stuff from halofit_ppf.f90
         fy = y/4.+y**2/8.
         ph = aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam))
         ph = ph/(1.+mu*y**(-1)+nu*y**(-2))*(1.+fnu*0.977)    ! CAMB; halofit_ppf.f90; halofit
         pq = plin*(1.+fnu*47.48*rk**2/(1.+1.5*rk**2))        ! CAMB; halofit_ppf.f90; halofit
         pq = plin*(1.+pq)**beta/(1.+pq*alpha)*exp(-fy)
      ELSE
         STOP 'HALOFIT: Error, ihf specified incorrectly'
      END IF

      pnl = pq+ph

   END SUBROUTINE HALOFIT

   REAL FUNCTION z_lss(cosm)

      ! Fit to the redshift of last scattering
      ! Equation (23) from Hu & White (1997) https://arxiv.org/pdf/astro-ph/9609079.pdf
      IMPLICIT NONE
      TYPE(cosmology), INTENT(IN) :: cosm
      REAL :: b1, b2
      REAL :: om_b, om_m

      om_m = cosm%Om_m*cosm%h**2
      om_b = cosm%Om_b*cosm%h**2

      b1 = (0.0783*om_b**(-0.238))/(1.+39.5*om_b**0.763)
      b2 = 0.560/(1.+21.1*om_b**1.81)

      z_lss = 1048.*(1.+0.00124*(om_b)**(-0.738))*(1.+b1*om_m**b2)

   END FUNCTION

   REAL FUNCTION r_sound(cosm)

      ! Sound horizon calculated in WKB approximation
      ! Equation (B6) from Hu & Sugiyama (1995) https://arxiv.org/pdf/astro-ph/9407093.pdf
      ! TODO: Code this up!
      IMPLICIT NONE
      TYPE(cosmology), INTENT(IN) :: cosm

      r_sound = cosm%A ! Warning suppression
      r_sound = 0.
      STOP 'R_SOUND: Code this up!!!'

   END FUNCTION r_sound

   SUBROUTINE init_wiggle(cosm)

      ! Isolate the power spectrum wiggle
      ! TODO: Smoothing should be over one wiggle period, not just fixed ns
      ! TODO: Convert this to some-sort of wiggle T(k), rather than an addition thing?
      ! TODO: This does not work very well
      USE special_functions
      IMPLICIT NONE
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: logkv(4), logpv(4)
      REAL, ALLOCATABLE :: logk(:), Pk(:), logPk(:)
      REAL, ALLOCATABLE :: Pk_smooth(:), logPk_smooth(:), Pk_wiggles(:)
      INTEGER :: i, nk

      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2
      INTEGER, PARAMETER :: wiggle_Lagrange = 1
      INTEGER, PARAMETER :: ns = 10
      INTEGER, PARAMETER :: wiggle_smooth = 2
      INTEGER :: wiggle_method = wiggle_Lagrange
      !INTEGER :: wiggle_method = wiggle_smooth

      IF (.NOT. ALLOCATED(cosm%log_k_plin)) STOP 'WIGGLE_INIT: Error, P(k) needs to be tabulated for this to work'

      nk = cosm%nk_plin

      ! Allocate arrays
      ALLOCATE (logk(nk), Pk(nk), logPk(nk))
      ALLOCATE (Pk_wiggles(nk), Pk_smooth(nk))

      ! Allocate the internal arrays from the cosmology arrays
      !k = exp(cosm%log_k_plin)
      IF(cosm%scale_dependent_growth) THEN
         logPk = cosm%log_plina(:,cosm%na_plin)
      ELSE
         logPk = cosm%log_plin
      END IF
      logk = cosm%log_k_plin
      Pk = exp(logPk)

      IF(wiggle_method == wiggle_Lagrange) THEN

         ! Fixed k values - CAMB
         logkv(1) = log(0.008)
         logkv(2) = log(0.01)
         logkv(3) = log(0.8)
         logkv(4) = log(1.0)

         ! Fix p values
         DO i = 1, 4
            logpv(i) = find(logkv(i), logk, logPk, nk, iorder, ifind, imeth)
         END DO

         ! Create new smooth spectrum that has no wiggles
         DO i = 1, nk
            IF (logk(i) <= logkv(1) .OR. logk(i) >= logkv(4)) THEN
               Pk_smooth(i) = Pk(i)
            ELSE
               Pk_smooth(i) = exp(Lagrange_polynomial(logk(i), 3, logkv, logpv))
            END IF
         END DO

      ELSE IF (wiggle_method == wiggle_smooth) THEN

         ALLOCATE (logPk_smooth(nk))
         logPk_smooth = logPk
         CALL smooth_array(logPk_smooth, nk, ns)
         Pk_smooth = exp(logPk_smooth)

      ELSE
         STOP 'INIT_WIGGLE: Error, wiggle method is not specified correctly'
      END IF

      ! Isolate just the wiggles
      ! It is difficult to make this operation (subtraction) fast given log arrays
      Pk_wiggles = Pk-Pk_smooth

      CALL if_allocated_deallocate(cosm%p_wiggle)
      ALLOCATE(cosm%p_wiggle(nk))
      cosm%p_wiggle = Pk_wiggles

      ! Set the flag
      cosm%has_wiggle = .TRUE.

   END SUBROUTINE init_wiggle

   REAL FUNCTION p_dewiggle(k, a, flag, sigv, cosm)

      ! Call the dewiggled power spectrum
      ! TODO: Convert this to a wiggle T(k)?
      IMPLICIT NONE
      REAL, INTENT(IN) :: k
      REAL, INTENT(IN) :: a
      INTEGER, INTENT(IN) :: flag
      REAL, INTENT(IN) :: sigv
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: p_wiggle, f, p_linear
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      p_linear = p_lin(k, a, flag, cosm) ! Needed here to make sure it is init before init_wiggle
      IF (.NOT. cosm%has_wiggle) CALL init_wiggle(cosm)
      p_wiggle = find(log(k), cosm%log_k_plin, cosm%p_wiggle, cosm%nk_plin, iorder, ifind, imeth)
      f = exp(-(k*sigv)**2)
      p_dewiggle = p_linear+(f-1.)*p_wiggle*grow(a, cosm)**2

   END FUNCTION p_dewiggle

   REAL FUNCTION Pk_Delta(Delta, k)

      ! Converts dimensionless Delta^2(k) to P(k) [Mpc/h]^3
      IMPLICIT NONE
      REAL, INTENT(IN) :: Delta ! Power spectrum in Delta^2(k) dimensionless form
      REAL, INTENT(IN) :: k     ! Wavenumber [h/Mpc]

      Pk_Delta = Delta/(k/twopi)**3
      Pk_Delta = Pk_Delta/(4.*pi)

   END FUNCTION Pk_Delta

   REAL FUNCTION Delta_Pk(Pk, k)

      ! Converts P(k) [Mpc/h]^3 to dimensionless Delta^2(k) 
      IMPLICIT NONE
      REAL, INTENT(IN) :: Pk ! Power spectrum in P(k) [Mpc/h]^3
      REAL, INTENT(IN) :: k  ! Wavenumber [h/Mpc]

      Delta_Pk = (4.*pi)*((k/twopi)**3)*Pk

   END FUNCTION Delta_Pk

END MODULE cosmology_functions