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
  PUBLIC :: p_lin
  PUBLIC :: sigma_all
  PUBLIC :: sigma_cold
  PUBLIC :: sigmaV
  PUBLIC :: xi_lin

  ! Halofit
  PUBLIC :: calculate_halofit_a
  
  INTERFACE integrate_cosm
     MODULE PROCEDURE integrate1_cosm
     MODULE PROCEDURE integrate2_cosm
     MODULE PROCEDURE integrate3_cosm
  END INTERFACE integrate_cosm

  ! Contains cosmological parameters that only need to be calculated once
  TYPE cosmology

     ! Primary parameters
     CHARACTER(len=256) :: name ! Name for cosmological model
     REAL :: Om_m, Om_b, Om_v, Om_w, m_nu      ! Primary parameters
     REAL :: h, n, sig8, w, wa, inv_m_wdm, YH  ! Primary parameters
     REAL :: a1, a2, ns, ws, am, dm, wm        ! DE parameters      
     REAL :: z_CMB, T_CMB, neff                ! Less primary parameters
     REAL :: Om_m_pow, Om_b_pow, h_pow         ! Cosmological parameters used for P(k) if different from background

     ! Derived parameters
     REAL :: A, Gamma, k                ! Power spectrum amplitude and shape parameter for DEFW
     REAL :: Om, Om_k, Om_c, Om_g, Om_r ! Derived Omegas
     REAL :: Om_nu, f_nu, a_nu, z_nu    ! Neutrinos
     REAL :: Om_nu_rad, omega_nu        ! Neutrinos
     REAL :: omega_m, omega_b, omega_c  ! Physical densities
     REAL :: Om_c_pow                   ! Cosmological parameters used for P(k) if different from background
     REAL :: age, horizon               ! Derived distance/time
     REAL :: mue, mup, YHe              ! Derived thermal parameters        
     REAL :: Om_ws, as, a1n, a2n        ! Derived DE parameters
     REAL :: gnorm                      ! Growth-factor normalisation

     ! Box size
     REAL :: Lbox, kbox
     LOGICAL :: box

     ! Switches
     INTEGER :: iw, itk 
     LOGICAL :: power_Omegas, derive_gas_numbers, growk
     
     ! Look-up tables
     REAL, ALLOCATABLE :: log_plin(:), log_k_plin(:), log_plina(:,:), log_a_plin(:)         ! Arrays for input linear P(k)
     REAL, ALLOCATABLE :: log_sigma(:), log_r_sigma(:), log_a_sigma(:), log_sigmaa(:,:)     ! Arrays for sigma(R)
     REAL, ALLOCATABLE :: log_a_growth(:), log_growth(:), growth_rate(:), log_acc_growth(:) ! Arrays for growth
     REAL, ALLOCATABLE :: log_p(:), log_a_p(:)         ! Arrays for distance (particle horizon)
     REAL, ALLOCATABLE :: log_t(:), log_a_t(:)         ! Arrays for time   
     REAL, ALLOCATABLE :: log_a_dcDv(:), dc(:), Dv(:)  ! Arrays for spherical-collapse parameters
     INTEGER :: n_growth, n_p, n_t, n_dcDv, nr_sigma, na_sigma, nk_plin, na_plin ! Array entries  
     LOGICAL :: has_distance, has_growth, has_sigma, has_spherical, has_power, has_time ! What has been calculated   

     ! Have normalisations and initialisations been done?
     LOGICAL :: is_init, is_normalised

     ! Verbose
     LOGICAL :: verbose
     
  END TYPE cosmology

! Global parameters
REAL, PARAMETER :: acc_cosm=1e-4 ! Global accuacy for the cosmological integrations

! Equation of state parameters
REAL, PARAMETER :: w_c=0.    ! CDM
REAL, PARAMETER :: w_b=0.    ! Baryons
REAL, PARAMETER :: w_g=1./3. ! Photons
REAL, PARAMETER :: w_v=-1.   ! Vacuum energy

! Distance
REAL, PARAMETER :: amin_distance=1e-4 ! Minimum scale factor in look-up table
REAL, PARAMETER :: amax_distance=1.   ! Maximum scale factor in look-up table
INTEGER, PARAMETER :: n_distance=128  ! Number of scale factor entries in look-up table
REAL, PARAMETER :: atay_distance=1e-5 ! Below this do a Taylor expansion to avoid divergence

! Time
REAL, PARAMETER :: amin_time=1e-4 ! Minimum scale factor in look-up table
REAL, PARAMETER :: amax_time=1.   ! Maximum scale factor in look-up table
INTEGER, PARAMETER :: n_time=128  ! Number of scale factor entries in look-up table
REAL, PARAMETER :: atay_time=1e-5 ! Below this do a Taylor expansion to avoid divergence

! Growth
REAL, PARAMETER :: ainit_growth=1e-3 ! Starting value for integratiton (should start | Omega_m(a)=1)
REAL, PARAMETER :: amax_growth=1.    ! Finishing value for integratiton (should be a=1)
INTEGER, PARAMETER :: n_growth=128   ! Number of entries for growth tables  

! Spherical collapse
REAL, PARAMETER :: amax_spherical=2.     ! Maximum scale factor to consider
REAL, PARAMETER :: dmin_spherical=1e-7   ! Minimum starting value for perturbation
REAL, PARAMETER :: dmax_spherical=1e-3   ! Maximum starting value for perturbation
INTEGER, PARAMETER :: m_spherical=128    ! Number of collapse scale-factors to try to calculate
INTEGER, PARAMETER :: n_spherical=100000 ! Number of points for ODE calculations
REAL, PARAMETER :: dinf_spherical=1e8    ! Value considered to be 'infinite' for the perturbation

! Power
REAL, PARAMETER :: kmin_plin=0.  ! Power below this wavenumber is set to zero [h/Mpc]
REAL, PARAMETER :: kmax_plin=1e8 ! Power above this wavenumber is set to zero [h/Mpc]
REAL, PARAMETER :: amin_plin=0.1 ! Minimum a value for P(k,a) tables when linear growth is scale dependent
REAL, PARAMETER :: amax_plin=1.0 ! Maximum a value for P(k,a) tables when linear growth is scale dependent
INTEGER, PARAMETER :: na_plin=16 ! Number of a values for sigma(R,a) tables

! sigma(R)
REAL, PARAMETER :: alpha_sigma=3.    ! I have made no attempt to optimise this number, nor tried alpha(R)
REAL, PARAMETER :: sigma_out=10.     ! How far out to go in 1/R units for sigma^2 split integral
REAL, PARAMETER :: Rsplit_sigma=1e-2 ! R value at which to split between the integration methods
REAL, PARAMETER :: rmin_sigma=1e-4   ! Minimum r value (NB. sigma(R) needs to be power-law below)
REAL, PARAMETER :: rmax_sigma=1e3    ! Maximum r value (NB. sigma(R) needs to be power-law above)
INTEGER, PARAMETER :: nr_sigma=128   ! Number of r entries for sigma(R) tables
REAL, PARAMETER :: amin_sigma=0.1    ! Minimum a value for sigma(R,a) tables when linear growth is scale dependent
REAL, PARAMETER :: amax_sigma=1.0    ! Maximum a value for sigma(R,a) tables when linear growth is scale dependent
INTEGER, PARAMETER :: na_sigma=16    ! Number of a values for sigma(R,a) tables

! sigma_v(R)
REAL, PARAMETER :: alpha_sigmaV=3.

CONTAINS

   SUBROUTINE assign_cosmology(icosmo,cosm,verbose)

    ! Assigns the 'primary' cosmological parameters (primary according to my definition)
    ! This routine *only* assigns parameters, it does and should not do *any* calculations
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(INOUT) :: icosmo
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i
    REAL :: Xe, Xi

    ! Names of pre-defined cosmologies    
    INTEGER, PARAMETER :: ncosmo=337
    CHARACTER(len=256) :: names(0:ncosmo)
    names=''
    names(0)='User defined'
    names(1)='Boring'
    names(2)='WMAP7 (cosmo-OWLS version; 1312.5462)'
    names(3)='Planck 2013 (cosmo-OWLS version; 1312.5462)'
    names(4)='WMAP9 (BAHAMAS version: 1712.02411)'
    names(5)='Open'
    names(6)='Einstein-de Sitter'
    names(7)='IDE I (user)'
    names(8)='IDE II (user)'
    names(9)='IDE III (user)'
    names(10)='IDE3'
    names(11)='IDE10'
    names(12)='LCDM (user)'
    names(13)='w(a)CDM (user)'
    names(14)='WDM'
    names(15)='TCDM'
    names(16)='Boring: w = -0.7'
    names(17)='Boring: w = -1.3'
    names(18)='Boring: w = -1.0; wa =  0.5'
    names(19)='Boring: w = -1.0; wa = -0.5'
    names(20)='Boring: w = -0.7; wa = -1.5'
    names(21)='Boring: w = -1.3; wa =  0.5'
    names(22)='IDE3'
    names(23)='IDE10'
    names(24)='Random Mira Titan cosmology'
    names(25)='Random Fanken Emu cosmology'
    names(26)='Boring: CAMB linear spectrum'
    names(27)='Illustris TNG 75'
    names(28)='Boring: Finite box'
    names(29)='Boring LCDM: Open; z=1 normalisation'
    names(30)='Boring LCDM: w = -0.7; z=1 normalisation'
    names(31)='Boring LCDM: w = -1.3; z=1 normalisation'
    names(32)='Boring LCDM: w = -1.0; wa =  0.5; z=1 normalisation'
    names(33)='Boring LCDM: w = -1.0; wa = -0.5; z=1 normalisation'
    names(34)='Boring LCDM: w = -0.7; wa = -1.5; z=1 normalisation'
    names(35)='Boring LCDM: w = -1.3; wa =  0.5; z=1 normalisation'
    names(36)='Einsten-de Sitter LCDM; z=1 normalisation'
    names(37)='Multidark: WMAP 5'
    names(38)='Random Cosmic Emu cosmology'
    names(39)='Random cosmology'
    names(40)='Random CAMB cosmology'
    names(41)='SCDM with high neutrino mass'
    
    names(100)='Mira Titan M000'
    names(101)='Mira Titan M001'
    names(102)='Mira Titan M002'
    names(103)='Mira Titan M003'
    names(104)='Mira Titan M004'
    names(105)='Mira Titan M005'
    names(106)='Mira Titan M006'
    names(107)='Mira Titan M007'
    names(108)='Mira Titan M008'
    names(109)='Mira Titan M009'
    names(110)='Mira Titan M010'
    names(111)='Mira Titan M011'
    names(112)='Mira Titan M012'
    names(113)='Mira Titan M013'
    names(114)='Mira Titan M014'
    names(115)='Mira Titan M015'
    names(116)='Mira Titan M016'
    names(117)='Mira Titan M017'
    names(118)='Mira Titan M018'
    names(119)='Mira Titan M019'
    names(120)='Mira Titan M020'
    names(121)='Mira Titan M021'
    names(122)='Mira Titan M022'
    names(123)='Mira Titan M023'
    names(124)='Mira Titan M024'
    names(125)='Mira Titan M025'
    names(126)='Mira Titan M026'
    names(127)='Mira Titan M027'
    names(128)='Mira Titan M028'
    names(129)='Mira Titan M029'
    names(130)='Mira Titan M030'
    names(131)='Mira Titan M031'
    names(132)='Mira Titan M032'
    names(133)='Mira Titan M033'
    names(134)='Mira Titan M034'
    names(135)='Mira Titan M035'
    names(136)='Mira Titan M036'
    
    names(200)='Franken Emu M000'
    names(201)='Franken Emu M001'
    names(202)='Franken Emu M002'
    names(203)='Franken Emu M003'
    names(204)='Franken Emu M004'
    names(205)='Franken Emu M005'
    names(206)='Franken Emu M006'
    names(207)='Franken Emu M007'
    names(208)='Franken Emu M008'
    names(209)='Franken Emu M009'
    names(210)='Franken Emu M010'
    names(211)='Franken Emu M011'
    names(212)='Franken Emu M012'
    names(213)='Franken Emu M013'
    names(214)='Franken Emu M014'
    names(215)='Franken Emu M015'
    names(216)='Franken Emu M016'
    names(217)='Franken Emu M017'
    names(218)='Franken Emu M018'
    names(219)='Franken Emu M019'
    names(220)='Franken Emu M020'
    names(221)='Franken Emu M021'
    names(222)='Franken Emu M022'
    names(223)='Franken Emu M023'
    names(224)='Franken Emu M024'
    names(225)='Franken Emu M025'
    names(226)='Franken Emu M026'
    names(227)='Franken Emu M027'
    names(228)='Franken Emu M028'
    names(229)='Franken Emu M029'
    names(230)='Franken Emu M030'
    names(231)='Franken Emu M031'
    names(232)='Franken Emu M032'
    names(233)='Franken Emu M033'
    names(234)='Franken Emu M034'
    names(235)='Franken Emu M035'
    names(236)='Franken Emu M036'
    names(237)='Franken Emu M037'

    names(300)='Cosmic Emu M000'
    names(301)='Cosmic Emu M001'
    names(302)='Cosmic Emu M002'
    names(303)='Cosmic Emu M003'
    names(304)='Cosmic Emu M004'
    names(305)='Cosmic Emu M005'
    names(306)='Cosmic Emu M006'
    names(307)='Cosmic Emu M007'
    names(308)='Cosmic Emu M008'
    names(309)='Cosmic Emu M009'
    names(310)='Cosmic Emu M010'
    names(311)='Cosmic Emu M011'
    names(312)='Cosmic Emu M012'
    names(313)='Cosmic Emu M013'
    names(314)='Cosmic Emu M014'
    names(315)='Cosmic Emu M015'
    names(316)='Cosmic Emu M016'
    names(317)='Cosmic Emu M017'
    names(318)='Cosmic Emu M018'
    names(319)='Cosmic Emu M019'
    names(320)='Cosmic Emu M020'
    names(321)='Cosmic Emu M021'
    names(322)='Cosmic Emu M022'
    names(323)='Cosmic Emu M023'
    names(324)='Cosmic Emu M024'
    names(325)='Cosmic Emu M025'
    names(326)='Cosmic Emu M026'
    names(327)='Cosmic Emu M027'
    names(328)='Cosmic Emu M028'
    names(329)='Cosmic Emu M029'
    names(330)='Cosmic Emu M030'
    names(331)='Cosmic Emu M031'
    names(332)='Cosmic Emu M032'
    names(333)='Cosmic Emu M033'
    names(334)='Cosmic Emu M034'
    names(335)='Cosmic Emu M035'
    names(336)='Cosmic Emu M036'
    names(337)='Cosmic Emu M037'

    IF(verbose) WRITE(*,*) 'ASSIGN_COSMOLOGY: Assigning cosmological model parameters'

    IF(icosmo==-1) THEN
       WRITE(*,*) 'ASSIGN_COSMOLOGY: Choose cosmological model'
       WRITE(*,*) '==========================================='
       DO i=0,SIZE(names)-1
          IF(i>=100 .AND. i<=136) THEN
             ! Do nothing
          ELSE IF(i>=200 .AND. i<=237) THEN
             ! Do nothing
          ELSE IF(i>=300 .AND. i<=337) THEN
             ! Do nothing
          ELSE IF(names(i) .NE. '') THEN
             WRITE(*,*) i, '- ', trim(names(i))
          END IF
       END DO
       WRITE(*,*) ' 100 -> 136 - Mira Titan M000 -> M036'
       WRITE(*,*) ' 200 -> 237 - Franken Emu M000 -> M037'
       WRITE(*,*) ' 300 -> 337 - Cosmic Emu M000 -> M037'
       READ(*,*) icosmo
       WRITE(*,*) '==========================================='
    END IF

    ! Set verbosity
    cosm%verbose=verbose

    ! Set the name of the cosmological model
    cosm%name=names(icosmo)

    ! Linear power spectrum
    ! 1 - Eisenstein & Hu
    ! 2 - CAMB
    ! 3 - DEFW
    ! 4 - External
    cosm%itk=1

    ! Boring default cosmology
    cosm%Om_m=0.3
    cosm%Om_b=0.05
    cosm%Om_v=1.-cosm%Om_m
    cosm%Om_w=0.
    cosm%m_nu=0.
    cosm%h=0.7
    cosm%sig8=0.8
    cosm%n=0.96
    cosm%w=-1.
    cosm%wa=0.
    cosm%T_CMB=2.725 ! CMB temperature [K]
    cosm%z_CMB=1087. ! Redshift of the last-scatting surface
    cosm%neff=3.046  ! Effective number of relativistic neutrinos
    cosm%YH=0.76     ! Hydrogen mass fraction

    ! Default dark energy is Lambda
    cosm%iw=1

    ! Default to have no WDM
    cosm%inv_m_wdm=0. ! Inverse WDM mass [1/keV]

    ! Consider box size
    cosm%box=.FALSE.
    cosm%Lbox=100.

    ! Set flags
    cosm%is_init=.FALSE.
    cosm%is_normalised=.FALSE.

    ! Omegas for power spectrum if different from background cosmological parameters
    cosm%Om_m_pow=0.
    cosm%Om_b_pow=0.
    cosm%h_pow=0.
    cosm%power_Omegas=.FALSE.

    ! Gas options
    cosm%derive_gas_numbers=.TRUE.

    IF(icosmo==0) THEN
       STOP 'TODO: implement user decision here'
    ELSE IF(icosmo==1) THEN
       ! Boring - do nothing
    ELSE IF(icosmo==2) THEN
       ! cosmo-OWLS - WMAP7 (1312.5462)
       cosm%Om_m=0.272
       cosm%Om_b=0.0455
       cosm%Om_v=1.-cosm%Om_m
       cosm%h=0.704
       cosm%sig8=0.81
       cosm%n=0.967
    ELSE IF(icosmo==3) THEN
       ! cosmo-OWLS - Planck 2013 (1312.5462)
       cosm%Om_m=0.3175
       cosm%Om_b=0.0490
       cosm%Om_v=1.-cosm%Om_m
       cosm%h=0.6711
       cosm%n=0.9624
       cosm%sig8=0.834
    ELSE IF(icosmo==4) THEN
       ! BAHAMAS - WMAP9 (1712.02411)
       cosm%h=0.7000
       cosm%Om_b=0.0463
       cosm%Om_m=0.2330+cosm%Om_b
       cosm%Om_v=1.-cosm%Om_m       
       cosm%n=0.9720
       cosm%sig8=0.8211
       cosm%derive_gas_numbers=.FALSE.
       cosm%mup=0.61
       Xi=1.08
       Xe=1.17
       cosm%mue=cosm%mup*(Xe+Xi)/Xe
    ELSE IF(icosmo==5) THEN
       !  5 - Open model
       ! 29 - Open model; normalised for Mead 2017 z=1 response
       cosm%Om_v=0.
    ELSE IF(icosmo==6 .OR. icosmo==15) THEN
       !  6 - Einstein-de Sitter (SCDM)
       ! 15 - Einstein-de Sitter (TCDM)
       IF(icosmo==15) THEN
          cosm%power_Omegas=.TRUE.
          cosm%Om_m_pow=cosm%Om_m
          cosm%Om_b_pow=cosm%Om_b
          cosm%h_pow=cosm%h
       END IF
       cosm%Om_m=1.
       cosm%Om_v=0.
    ELSE IF(icosmo==7) THEN
       ! IDE I
       cosm%iw=5
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%as
       WRITE(*,*) 'n*:'
       READ(*,*) cosm%ns
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       cosm%Om_m=0.3
       cosm%Om_v=0.7
    ELSE IF(icosmo==8) THEN
       ! IDE II model
       cosm%iw=6      
       WRITE(*,*) 'n*:'
       READ(*,*) cosm%ns
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%as
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       cosm%Om_m=0.3
       cosm%Om_w=0.7
       cosm%Om_v=0. !No vacuum necessary here
    ELSE IF(icosmo==9) THEN
       ! IDE III model
       cosm%iw=7
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%as
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       WRITE(*,*) 'w*:'
       READ(*,*) cosm%ws
       cosm%Om_m=0.3
       cosm%Om_w=0.7
       cosm%Om_v=0.
    ELSE IF(icosmo==10 .OR. icosmo==11) THEN
       cosm%iw=6
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
       IF(icosmo==10) cosm%ns=3
       IF(icosmo==11) cosm%ns=10
       cosm%as=0.01
       cosm%Om_ws=0.2
    ELSE IF(icosmo==12) THEN
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Om_v:'
       READ(*,*) cosm%Om_v
    ELSE IF(icosmo==13) THEN
       cosm%iw=3
       WRITE(*,*) 'w0:'
       READ(*,*) cosm%w
       WRITE(*,*) 'wa:'
       READ(*,*) cosm%wa
    ELSE IF(icosmo==14) THEN
       ! WDM
       cosm%inv_m_wdm=1.
    ELSE IF(icosmo==16) THEN
       ! w = -0.7
       cosm%iw=4
       cosm%w=-0.7
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
    ELSE IF(icosmo==17) THEN
       ! w = -1.3
       cosm%iw=4
       cosm%w=-1.3
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
    ELSE IF(icosmo==18) THEN
       ! wa = 0.5
       cosm%iw=3
       cosm%wa=0.5
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
    ELSE IF(icosmo==19) THEN
       ! wa = -0.5
       cosm%iw=3
       cosm%wa=-0.5
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
    ELSE IF(icosmo==20) THEN
       ! w = -0.7; wa = -1.5
       cosm%iw=3
       cosm%w=-0.7
       cosm%wa=-1.5
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
    ELSE IF(icosmo==21) THEN
       ! w = -1.3; wa = 0.5
       cosm%iw=3
       cosm%w=-1.3
       cosm%wa=0.5
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
    ELSE IF(icosmo==22 .OR. icosmo==23) THEN
       ! IDE II models
       cosm%iw=6
       cosm%Om_m=0.3
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0. ! No vacuum necessary here
       IF(icosmo==22) THEN
          ! IDE 3
          cosm%ns=3
          cosm%as=0.01
          cosm%Om_ws=0.1
       ELSE IF(icosmo==23) THEN
          ! IDE 10
          cosm%ns=10
          cosm%as=0.1
          cosm%Om_ws=0.02
       END IF
    ELSE IF(icosmo==24) THEN
       ! Random Mira Titan cosmology
       CALL random_Mira_Titan_cosmology(cosm)
       cosm%itk=2 ! Set to CAMB linear power
       cosm%iw=3 ! Set to w(a) dark energy
       cosm%Om_v=0. ! Necessary for CAMB
    ELSE IF(icosmo==25) THEN
       ! Random Franken Emu cosmology
       CALL random_Franken_Emu_cosmology(cosm)   
       cosm%itk=2 ! Set to CAMB linear power
       cosm%iw=4 ! Set to constant w dark energy
       cosm%Om_v=0. ! Necessary for CAMB
    ELSE IF(icosmo==26) THEN
       ! Boring with CAMB linear spectrum
       cosm%itk=2 ! Set to CAMB linear power
       cosm%Om_w=cosm%Om_v ! Necessary for CAMB
       cosm%Om_v=0. ! Necessary for CAMB
    ELSE IF(icosmo==27) THEN
       ! Illustris; L = 75 Mpc/h
       cosm%itk=1 ! Set to CAMB linear power
       cosm%Om_m=0.3089
       cosm%Om_b=0.0486
       cosm%Om_w=1.-cosm%Om_m ! Necessary to be Omega_w for CAMB
       cosm%h=0.6774
       cosm%n=0.9667
       cosm%sig8=0.8159       
       cosm%Om_v=0. ! Necessary for CAMB
       cosm%box=.TRUE.
       cosm%Lbox=75. ! 75 Mpc/h box
    ELSE IF(icosmo==28) THEN
       ! Finite box
       cosm%box=.TRUE.
    ELSE IF(icosmo==29) THEN
       ! Boring: Open; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.88397
    ELSE IF(icosmo==30) THEN
       ! Boring: w = -0.7; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.83253
    ELSE IF(icosmo==31) THEN
       ! Boring: w = -1.3; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.77462
    ELSE IF(icosmo==32) THEN
       ! Boring: w = -1.0; wa =  0.5; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.81022
    ELSE IF(icosmo==33) THEN
       ! Boring: w = -1.0; wa = -0.5; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.79116
    ELSE IF(icosmo==34) THEN
       ! Boring: w = -0.7; wa = -1.5; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.80090
    ELSE IF(icosmo==35) THEN
       ! Boring: w = -1.3; wa =  0.5; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.78197
    ELSE IF(icosmo==36) THEN
       ! Boring; EdS; z=1 normalisation for Mead 2017; LCDM
       cosm%sig8=0.65380
    ELSE IF(icosmo==37) THEN
       ! Multidark: WMAP5
       cosm%h=0.70
       cosm%Om_b=0.0469
       cosm%Om_m=0.27
       cosm%Om_v=1.-cosm%Om_m
       cosm%n=0.95
       cosm%sig8=0.82
    ELSE IF(icosmo==38) THEN
       ! Random cosmic emu model
       CALL random_Cosmic_Emu_cosmology(cosm)   
       cosm%itk=2 ! Set to CAMB linear power
       cosm%iw=4 ! Set to constant w dark energy
       cosm%Om_v=0. ! Necessary for CAMB
    ELSE IF(icosmo==39) THEN
       ! Random cosmology
       CALL random_cosmology(cosm)
       cosm%itk=1
       cosm%Om_v=0.
       cosm%iw=4
    ELSE IF(icosmo==40) THEN
       ! Random cosmology
       CALL random_cosmology(cosm)
       cosm%itk=2
       cosm%Om_v=0.
       cosm%iw=4
    ELSE IF(icosmo==41) THEN
       ! SCDM with high neutrino mass
       cosm%Om_m=1.
       cosm%Om_v=0.
       cosm%m_nu=4.
    ELSE IF(icosmo>=100 .AND. icosmo<=137) THEN
       ! Mira Titan nodes
       CALL Mira_Titan_node_cosmology(icosmo-100,cosm)
       cosm%itk=2   ! Set to CAMB linear power
       cosm%iw=3    ! Set to w(a) dark energy
       cosm%Om_v=0. ! Necessary for CAMB
    ELSE IF(icosmo>=200 .AND. icosmo<=237) THEN
       ! Franken Emu nodes (which are the same as Franken Emu nodes)
       CALL Franken_Emu_node_cosmology(icosmo-200,cosm)
       cosm%itk=2   ! Set to CAMB linear power
       cosm%iw=4    ! Set to constant w dark energy
       cosm%Om_v=0. ! Necessary for CAMB
    ELSE IF(icosmo>=300 .AND. icosmo<=337) THEN
       ! Cosmic Emu nodes (which are the same as Franken Emu nodes)
       CALL Cosmic_Emu_node_cosmology(icosmo-300,cosm)
       cosm%itk=2   ! Set to CAMB linear power
       cosm%iw=4    ! Set to constant w dark energy
       cosm%Om_v=0. ! Necessary for CAMB
    ELSE
       STOP 'ASSIGN_COSMOLOGY: Error, icosmo not specified correctly'
    END IF

    IF(cosm%verbose) THEN
       WRITE(*,*) 'ASSIGN_COSMOLOGY: Cosmology: ', trim(cosm%name)
       WRITE(*,*) 'ASSIGN_COSMOLOGY: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE assign_cosmology

  SUBROUTINE init_cosmology(cosm)

    ! Calcualtes derived parameters
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: Xs, f1, f2
    REAL :: rho_g, Om_g_h2, rho_crit, f_nu_rad
    REAL, PARAMETER :: small=1e-5 ! Some small number for writing curvature things

    ! Set all 'has/is' to false
    cosm%has_distance=.FALSE.
    cosm%has_growth=.FALSE.
    cosm%has_sigma=.FALSE.
    cosm%has_spherical=.FALSE.
    cosm%has_power=.FALSE.
    cosm%is_init=.FALSE.
    cosm%is_normalised=.FALSE.

    ! Things to do with finite box
    IF(cosm%box) cosm%kbox=twopi/cosm%Lbox

    IF(cosm%verbose) WRITE(*,*) 'INIT_COSMOLOGY: Calculating derived parameters'

    ! Calculate radiation density (includes photons and neutrinos at recombination)
    rho_g=(4.*SBconst*cosm%T_CMB**4/c_light**3) ! Photon physical density at z=0 from CMB temperature [kg/m^3]
    rho_crit=3.*H0**2/(8.*pi*bigG) ! Critical density [h^2 kg/m^3] TODO: Constants?
    Om_g_h2=rho_g/rho_crit ! Photon cosmological density [h^2]
    cosm%Om_g=Om_g_h2/cosm%h**2 ! Photon density parameter
    cosm%Om_nu_rad=cosm%Om_g*neff_constant*cosm%neff ! Relativisitic neutrino density
    cosm%Om_r=cosm%Om_g+cosm%Om_nu_rad ! Radiation is sum of photon and neutrino densities
    f_nu_rad=cosm%Om_nu_rad/cosm%Om_r ! Fraction of radiation that is in neutrinos (~0.40)

    ! Information about how radiation density is calculated
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_g:', cosm%Om_g
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_nu (radiation):', cosm%Om_nu_rad
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_r:', cosm%Om_r
       WRITE(*,*) 'INIT_COSMOLOGY: f_nu (radiation):', f_nu_rad          
    END IF

    ! Check that radiation density is not absurd
    IF(cosm%Om_r>1e-3) THEN
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
    cosm%Om_nu=cosm%m_nu/(neutrino_constant*cosm%h**2)
    cosm%f_nu=cosm%Om_nu/cosm%Om_m
    IF(cosm%m_nu .NE. 0.) THEN
       !cosm%a_nu=1./(1900.*cosm%m_nu)
       cosm%a_nu=cosm%Om_nu_rad/cosm%Om_nu
    ELSE
       cosm%a_nu=1.
    END IF
    cosm%z_nu=redshift_a(cosm%a_nu)
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_nu (matter):', cosm%Om_nu
       WRITE(*,*) 'INIT_COSMOLOGY: a_nu:', cosm%a_nu
       WRITE(*,*) 'INIT_COSMOLOGY: z_nu:', cosm%z_nu
       WRITE(*,*) 'INIT_COSMOLOGY: f_nu:', cosm%f_nu
    END IF

    ! Check neutrino mass fraction is not too high
    IF(cosm%f_nu>0.5) STOP 'INIT_COSMOLOGY: Error, neutrino mass fraction is too high'
    IF((cosm%m_nu .NE. 0.) .AND. cosm%a_nu>0.2) STOP 'INIT_COSMOLOGY: Error, neutrinos are too light'

    ! Decide on scale-dependent growth
    IF(cosm%m_nu .NE. 0.) THEN
       cosm%growk=.TRUE.
       IF(cosm%verbose) WRITE(*,*) 'INIT_COSMOLOGY: Scale-dependent growth'
    ELSE
       cosm%growk=.FALSE.
    END IF

    ! Check that we are able to cope with scale-dependent growth
    !IF(cosm%growk .AND. (cosm%itk .NE. 2)) STOP 'INIT_COSMOLOGY: Error, scale-dependent growth requires CAMB for linear spectra'

    ! Derived cosmological parameters    
    cosm%Om_c=cosm%Om_m-cosm%Om_b-cosm%Om_nu ! Omega_m defined to include CDM, baryons and massive neutrinos
    !cosm%Om=cosm%Om_m+cosm%Om_v_mod+cosm%Om_r+cosm%Om_w
    cosm%Om=cosm%Om_m+cosm%Om_v+cosm%Om_w ! Ignore radiation here
    cosm%Om_k=1.-cosm%Om
    cosm%k=(cosm%Om-1.)/(Hdist**2)
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_c:', cosm%Om_c
       WRITE(*,*) 'INIT_COSMOLOGY: Omega:', cosm%Om
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_k:', cosm%Om_k
       WRITE(*,*) 'INIT_COSMOLOGY: k [Mpc/h]^-2:', cosm%k
       IF(abs(cosm%k)>small) THEN
          WRITE(*,*) 'INIT_COSMOLOGY: k_rad [Mpc/h]:', 1./sqrt(abs(cosm%k)) ! Curvature radius
       END IF
    END IF

    IF(cosm%Om_c<0.) STOP 'INIT_COSMOLOGY: Error, CDM density is negative'

    ! Physical density parameters
    cosm%omega_m=cosm%Om_m*cosm%h**2   ! Physical matter density
    cosm%omega_b=cosm%Om_b*cosm%h**2   ! Physical baryon density
    cosm%omega_c=cosm%Om_c*cosm%h**2   ! Physical CDM density
    cosm%omega_nu=cosm%Om_nu*cosm%h**2 ! Physical neutrino density

    ! Write physical densities to screen
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: omega_m:', cosm%omega_m
       WRITE(*,*) 'INIT_COSMOLOGY: omega_b:', cosm%omega_b
       WRITE(*,*) 'INIT_COSMOLOGY: omega_c:', cosm%omega_c
       WRITE(*,*) 'INIT_COSMOLOGY: omega_nu:', cosm%omega_nu 
    END IF

    ! Using different background Omegas compared to power Omegas
    IF(cosm%power_Omegas) THEN
       cosm%Om_c_pow=cosm%Om_m_pow-cosm%Om_b_pow
    ELSE
       cosm%Om_c_pow=0.
    END IF

    ! Gas parameters
    IF(cosm%derive_gas_numbers) THEN       
       ! Mean mass per gas particle divided by proton mass
       ! ~0.588 if fH=0.76, gas is ionised and H and He only; 0.61 in BAHAMAS
       cosm%mup=4./(5.*cosm%YH+3.)      
       ! Mean mass per gas electron divided by proton mass
       ! ~1.136 if fH=0.76, gas is ionised and H and He only; 1.17 in BAHAMAS
       cosm%mue=2./(1.+cosm%YH)       
       cosm%YHe=1.-cosm%YH ! Helium mass fraction      
    END IF

    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: mu_p:', cosm%mup
       WRITE(*,*) 'INIT_COSMOLOGY: mu_e:', cosm%mue
    END IF

    ! Gamma for DEFW
    cosm%Gamma=cosm%Om_m*cosm%h
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Gamma:', cosm%Gamma
    END IF
    
    cosm%is_init=.TRUE.

    ! Dark energy models
    IF(cosm%iw==5) THEN
       !Om_w=Om_w*(Om_m*astar**(-3)+Om_v)/(X(astar)*(1.-Om_w))
      f1=cosm%Om_ws*X_de(cosm%as,cosm)+cosm%Om_ws*cosm%as**(-2)
      f2=X_de(cosm%as,cosm)*(1.-cosm%Om_ws)+cosm%Om_ws*cosm%as**(-2)
       cosm%Om_w=cosm%Om_ws*(Hubble2(cosm%a,cosm)-f1/f2)
       !cosm%Om_w=cosm%Om_ws*(Hubble2(cosm%a,cosm)-cosm%Om_ws*X_de(cosm%as,cosm)+cosm%Om_ws*cosm%as**(-2))/(X_de(cosm%as,cosm)*(1.-cosm%Om_ws)+cosm%Om_ws*cosm%as**(-2))
    ELSE IF(cosm%iw==6) THEN
       ! Define a1^n
       cosm%a1n=cosm%as**cosm%ns
       ! Necessary for first step below
       cosm%a2n=cosm%a1n 
       ! All neccessary to convert parameters to a1,a2
       f1=cosm%Om_ws*(Hubble2(cosm%as,cosm)-cosm%Om_w*X_de(cosm%as,cosm))
       f2=cosm%Om_w*(1.-cosm%Om_ws)
       Xs=f1/f2
       Xs=Xs**(cosm%ns/6.)
       ! Top and bottom of fraction
       f1=cosm%a1n*(2.*Xs-(1.+cosm%a1n))
       f2=(1.+cosm%a1n)-2.*Xs*cosm%a1n
       cosm%a2n=f1/f2 ! Finally! a2
       !IF(a2<a1) a2=a1
    ELSE IF(cosm%iw==7) THEN
       ! Scale-factor at which Om_w(a*) is most important
       cosm%a1=cosm%as
       ! Needs to be set for X(a*) and H2(a*) below (which cancel each other)
       cosm%a2=cosm%as 
       f1=cosm%Om_ws*(Hubble2(cosm%as,cosm)-cosm%Om_w*X_de(cosm%as,cosm))
       f2=cosm%Om_w*(1.-cosm%Om_ws)
       cosm%a2=cosm%as*(f1/f2)**(1./(3.*(1.+cosm%ws)))
    END IF

    ! Ensure deallocate distances
    cosm%has_distance=.FALSE.
    IF(ALLOCATED(cosm%log_p))   DEALLOCATE(cosm%log_p)
    IF(ALLOCATED(cosm%log_a_p)) DEALLOCATE(cosm%log_a_p)
    cosm%horizon=0.

    ! Ensure deallocate time
    cosm%has_time=.FALSE.
    IF(ALLOCATED(cosm%log_t))   DEALLOCATE(cosm%log_t)
    IF(ALLOCATED(cosm%log_a_t)) DEALLOCATE(cosm%log_a_t)
    cosm%age=0.

    ! Ensure deallocate growth
    cosm%has_growth=.FALSE.
    IF(ALLOCATED(cosm%log_a_growth))    DEALLOCATE(cosm%log_a_growth)
    IF(ALLOCATED(cosm%log_growth))      DEALLOCATE(cosm%log_growth)
    IF(ALLOCATED(cosm%growth_rate))     DEALLOCATE(cosm%growth_rate)
    IF(ALLOCATED(cosm%log_acc_growth))  DEALLOCATE(cosm%log_acc_growth)
    cosm%gnorm=0.

    ! Ensure deallocate sigma
    cosm%has_sigma=.FALSE.
    IF(ALLOCATED(cosm%log_r_sigma)) DEALLOCATE(cosm%log_r_sigma)
    IF(ALLOCATED(cosm%log_a_sigma)) DEALLOCATE(cosm%log_a_sigma)
    IF(ALLOCATED(cosm%log_sigma))   DEALLOCATE(cosm%log_sigma)
    IF(ALLOCATED(cosm%log_sigmaa))  DEALLOCATE(cosm%log_sigmaa)

    ! Switch for power?
    IF(cosm%itk==1 .OR. cosm%itk==3) THEN
       ! Default to use internal linear P(k) from Eisenstein & Hu
       cosm%has_power=.FALSE.
    ELSE
       cosm%has_power=.TRUE.
    END IF
    
    ! Ensure deallocate linear-power tables
    !IF(cosm%has_power .EQV. .FALSE.) THEN
    IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
    IF(ALLOCATED(cosm%log_a_plin)) DEALLOCATE(cosm%log_a_plin)
    IF(ALLOCATED(cosm%log_plin))   DEALLOCATE(cosm%log_plin)
    IF(ALLOCATED(cosm%log_plina))  DEALLOCATE(cosm%log_plina)   
    !END IF

    ! Ensure delloacte spherical-collapse arrays
    cosm%has_spherical=.FALSE.
    IF(ALLOCATED(cosm%log_a_dcDv)) DEALLOCATE(cosm%log_a_dcDv)
    IF(ALLOCATED(cosm%dc))         DEALLOCATE(cosm%dc)
    IF(ALLOCATED(cosm%Dv))         DEALLOCATE(cosm%Dv)

    ! Write finishing message to screen
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Done'
       WRITE(*,*)
    END IF
    
  END SUBROUTINE init_cosmology

  SUBROUTINE print_cosmology(cosm)

    ! Prints the cosmological parameters to the screen
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, PARAMETER :: small=1e-5

    IF(cosm%verbose) THEN
       WRITE(*,*) '===================================='
       WRITE(*,*) 'COSMOLOGY: ', trim(cosm%name)
       WRITE(*,*) '===================================='
       WRITE(*,*) 'COSMOLOGY: Standard parameters'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_m:', cosm%Om_m
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_b:', cosm%Om_b  
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_v:', cosm%Om_v
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_w:', cosm%Om_w
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'h:', cosm%h
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'sigma_8:', cosm%sig8
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n_s:', cosm%n
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'T_CMB [K]:', cosm%T_CMB
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'z_CMB:', cosm%z_CMB
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n_eff:', cosm%neff
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Y_H:', cosm%YH
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'm_nu [eV]:', cosm%m_nu
       !WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'm_nu 1 [eV]:', cosm%m_nu(1)
       !WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'm_nu 2 [eV]:', cosm%m_nu(2)
       !WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'm_nu 3 [eV]:', cosm%m_nu(3)
       IF(cosm%box) WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'L_box [Mpc/h]:', cosm%Lbox
       IF(cosm%inv_m_wdm .NE. 0.) THEN
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'm_wdm [keV]:', 1./cosm%inv_m_wdm
       END IF
       WRITE(*,*) '===================================='
       WRITE(*,*) 'COSMOLOGY: Derived parameters'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'omega_m:', cosm%omega_m
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'omega_c:', cosm%omega_c
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'omega_b:', cosm%omega_b
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'omega_nu:', cosm%omega_nu
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_g:', cosm%Om_g
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_r:', cosm%Om_r
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega:', cosm%Om
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_c:', cosm%Om_c
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_nu:', cosm%Om_nu
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'f_nu:', cosm%f_nu
!!$       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_v'':', cosm%Om_v_mod
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_k:', cosm%Om_k
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'k [Mpc/h]^-2:', cosm%k
       IF(abs(cosm%k)>small) THEN
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'k_rad [Mpc/h]:', 1./sqrt(abs(cosm%k))
       END IF
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'mu_p:', cosm%mup
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'mu_e:', cosm%mue
       WRITE(*,*) '===================================='
       IF(cosm%iw==1) THEN
          WRITE(*,*) 'COSMOLOGY: Vacuum energy'
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w:', -1.
       ELSE IF(cosm%iw==2) THEN
          WRITE(*,*) 'COSMOLOGY: QUICC dark energy prescription'
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w0:', cosm%w
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'wm:', cosm%wm
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'am:', cosm%am
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'dm:', cosm%dm
       ELSE IF(cosm%iw==3) THEN
          WRITE(*,*) 'COSMOLOGY: w(a) = w0+wa(1.-a)'
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w0:', cosm%w
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'wa:', cosm%wa
       ELSE IF(cosm%iw==4) THEN
          WRITE(*,*) 'COSMOLOGY: Constant w'
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w:', cosm%w
       ELSE IF(cosm%iw==5) THEN
          WRITE(*,*) 'COSMOLOGY: IDE I'
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a*:', cosm%as
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n*:', cosm%ns
       ELSE IF(cosm%iw==6) THEN
          WRITE(*,*) 'COSMOLOGY: IDE II'
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a*:', cosm%as
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n*:', cosm%ns
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a1^n (derived):', cosm%a1n
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a2^n (derived):', cosm%a2n
       ELSE IF(cosm%iw==7) THEN
          WRITE(*,*) 'COSMOLOGY: IDE III'
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a*:', cosm%a1
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
          WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w*:', cosm%ws
       END IF
       WRITE(*,*) '===================================='
       WRITE(*,*)
    END IF

  END SUBROUTINE print_cosmology

  REAL FUNCTION xi_lin(r,a,cosm)

    ! Computes the 3D linear matter correlation function by integrating over P(k)
    IMPLICIT NONE
    REAL, INTENT(IN) :: r
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: k1, k2, xi_bit
    DOUBLE PRECISION :: xi8
    INTEGER, PARAMETER :: min_humps=5
    INTEGER, PARAMETER :: max_humps=100000
    INTEGER, PARAMETER :: method=1

    STOP 'XI_LIN: This is ridiculuously slow for large R'

    IF(method==1) THEN

       xi_lin=integrate_cosm(0.,1.,xi_integrand_transformed,r,a,cosm,acc_cosm,3)

    ELSE IF(method==2) THEN

       ! Set summation variable to zero
       xi8=0.

       ! Loop over humps
       DO i=0,max_humps

          k1=i*pi/r
          k2=(i+1)*pi/r

          xi_bit=integrate_cosm(k1,k2,xi_integrand,r,a,cosm,acc_cosm,3)

          xi8=xi8+xi_bit

          IF(i>min_humps) THEN
             IF(abs(xi_bit/real(xi8))<acc_cosm) THEN
                EXIT
             END IF
          END IF

          IF(i==max_humps) THEN
             WRITE(*,*) 'XI_LIN: r [Mpc/h]:', r
             WRITE(*,*) 'XI_LIN: Minimum number of humps:', min_humps
             WRITE(*,*) 'XI_LIN: Maximum number of humps:', max_humps
             WRITE(*,*) 'XI_LIN: Warning, maximum number of humps exceeded'
             STOP
          END IF

       END DO

       xi_lin=real(xi8)

    ELSE

       STOP 'XI_LIN: Error, method specified incorrectly'
       
    END IF
    
  END FUNCTION xi_lin

  REAL FUNCTION xi_integrand(k,r,a,cosm)

    ! Integrand for the 3D linear matter correlation function
    USE special_functions
    IMPLICIT NONE
    REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc] (integrated over)
    REAL, INTENT(IN) :: r ! Separation [Mpc/h]
    REAL, INTENT(IN) :: a ! Scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

    IF(k==0.) THEN
       xi_integrand=0.
    ELSE
       xi_integrand=sinc(k*r)*p_lin(k,a,cosm)/k
    END IF
    
  END FUNCTION xi_integrand

  REAL FUNCTION xi_integrand_transformed(t,r,a,cosm)

    ! Integrand for the 3D linear matter correlation function
    ! TODO: Optimise alpha(r)
    USE special_functions
    IMPLICIT NONE
    REAL, INTENT(IN) :: t ! kr=(-1+1/t)^alpha (integrated over 0:1)
    REAL, INTENT(IN) :: r ! Separation [Mpc/h]
    REAL, INTENT(IN) :: a ! Scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
    REAL :: kr, k, alpha
    REAL, PARAMETER :: rsplit=10.

    IF(r<rsplit) THEN
       alpha=2.
    ELSE
       alpha=1.5
    END IF
    
    IF(t==0. .OR. t==1.) THEN
       xi_integrand_transformed=0.
    ELSE
       kr=(-1.+1./t)**alpha
       k=kr/r
       xi_integrand_transformed=sinc(kr)*p_lin(k,a,cosm)*alpha/(t*(1.-t))
    END IF
    
  END FUNCTION xi_integrand_transformed

  SUBROUTINE normalise_power(cosm)

    ! Get the required sigma_8 by re-normalising the power spectrum
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: sigi, sigf, kbox
    REAL, PARAMETER :: R=8. ! Because we are doing sigma(R=8Mpc/h) normalisation
    REAL, PARAMETER :: a=1. ! Because we are doing simga(R=8Mpc/h,a=1) normalisation
    LOGICAL, PARAMETER :: run_twice=.FALSE. ! This is almost always a stupid thing to do

    ! Need to give this a value otherwise get a warning in debug mode
    kbox=0.

    ! Change the flag *before* doing this calculation because it calls power
    cosm%is_normalised=.TRUE.

    IF(cosm%itk==1) THEN

       ! Remove the k-cut for the normalisation
       IF(cosm%box) THEN
          kbox=cosm%kbox
          cosm%kbox=0.
       END IF

       ! This needs to be set here for the sigma routines below to work
       cosm%A=1.

       IF(cosm%verbose) WRITE(*,*) 'NORMALISE_POWER: Normalising power to get correct sigma_8'

       ! Calculate the initial sigma_8 value (will not be correct)
       sigi=sigma_all_integral(R,a,cosm)

       IF(cosm%verbose) WRITE(*,*) 'NORMALISE_POWER: Initial sigma_8:', real(sigi)

       ! Reset the normalisation to give the correct sigma8
       cosm%A=cosm%sig8/sigi
       !cosm%A=391.0112 ! Appropriate for sigma_8=0.8 in the boring model (for tests)

       ! Recalculate sigma8, should be correct this time
       sigf=sigma_all_integral(R,a,cosm)

       ! Write to screen
       IF(cosm%verbose) THEN
          WRITE(*,*) 'NORMALISE_POWER: Normalisation factor:', real(cosm%A)
          WRITE(*,*) 'NORMALISE_POWER: Target sigma_8:', real(cosm%sig8)
          WRITE(*,*) 'NORMALISE_POWER: Final sigma_8 (calculated):', real(sigf)
          WRITE(*,*) 'NORMALISE_POWER: Done'
          WRITE(*,*)
       END IF

       ! Replace the k-cut
       IF(cosm%box) THEN
          cosm%kbox=kbox
       END IF

    ELSE IF(cosm%itk==2) THEN

       ! Run first time to get power
       cosm%A=2.1e-9
       CALL get_CAMB_power(non_linear=.FALSE.,halofit_version=5,cosm=cosm)
       sigi=sigma_all_integral(R,a,cosm)

       IF(cosm%verbose) THEN
          WRITE(*,*) 'NORMALISE_POWER: Normalising power to get correct sigma_8'
          WRITE(*,*) 'NORMALISE_POWER: Initial As:', real(cosm%A)
          WRITE(*,*) 'NORMALISE_POWER: Initial sigma_8:', real(sigi)          
       END IF

       ! Normalisation
       IF(run_twice) THEN
          ! Run again to normalise
          cosm%A=cosm%A*(cosm%sig8/sigi)**2 
          CALL get_CAMB_power(non_linear=.FALSE.,halofit_version=5,cosm=cosm)          
       ELSE
          ! Normalise using sigma8 and rescaling linear power
         IF(cosm%growk) THEN
            cosm%log_plina=cosm%log_plina+2.*log(cosm%sig8/sigi)
         ELSE
            cosm%log_plin=cosm%log_plin+2.*log(cosm%sig8/sigi)
         END IF
       END IF

       ! Check that the normalisation has been done correctly
       sigf=sigma_all_integral(R,a,cosm)

       ! Write to screen
       IF(cosm%verbose) THEN
          WRITE(*,*) 'NORMALISE_POWER: Target sigma_8:', real(cosm%sig8)
          WRITE(*,*) 'NORMALISE_POWER: Final sigma_8 (calculated):', real(sigf)
          WRITE(*,*) 'NORMALISE_POWER: Done'
          WRITE(*,*)
       END IF

    ELSE

       STOP 'NORMALISE_POWER: Error, cannot normalise with this itk'
       
    END IF

  END SUBROUTINE normalise_power

  REAL FUNCTION comoving_critical_density(a,cosm)

    ! Comoving critical density [(Msun/h) / (Mpc/h)^3]
    ! For LCDM this is constant in the past, increases like a^3 in the future
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    comoving_critical_density=physical_critical_density(a,cosm)*a**3

  END FUNCTION comoving_critical_density

  REAL FUNCTION physical_critical_density(a,cosm)

    ! Physical critical density [(Msun/h) / (Mpc/h)^3]
    ! For LCDM tends to a constant in the future, behaves like a^-3 in the past
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_critical_density=critical_density_cos*Hubble2(a,cosm)

  END FUNCTION physical_critical_density

  REAL FUNCTION comoving_matter_density(cosm)

    ! Comoving matter density [(Msun/h) / (Mpc/h)^3]
    ! Not a function of redshift, constant value throughout time
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm

    comoving_matter_density=critical_density_cos*cosm%Om_m

  END FUNCTION comoving_matter_density

  REAL FUNCTION physical_matter_density(a,cosm)

    ! Physical matter density [(Msun/h) / (Mpc/h)^3]
    ! Proportional to a^-3 always
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_matter_density=comoving_matter_density(cosm)*a**(-3)

  END FUNCTION physical_matter_density

  REAL FUNCTION Hubble2(a,cosm)

    ! Calculates Hubble^2 in units such that H^2(a=1)=1
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'HUBBLE2: Error, cosmology is not initialised'
    !Hubble2=cosm%Om_m*a**(-3)+cosm%Om_r*a**(-4)+cosm%Om_v_mod+cosm%Om_w*X_de(a,cosm)+(1.-cosm%Om)*a**(-2)
    Hubble2=&
         cosm%Om_c*X_c(a)+&
         cosm%Om_b*X_b(a)+&
         cosm%Om_g*X_r(a)+&
         Omega_nu_0(a,cosm)*X_nu(a,cosm)+&
         cosm%Om_v*X_v(a)+&
         cosm%Om_w*X_de(a,cosm)+&
         (1.-cosm%Om)*a**(-2)

  END FUNCTION Hubble2

  REAL FUNCTION Hubble2_norad(a,cosm)

    ! Squared Hubble parameter but ignoring the radiation component
    ! Units such that Hubble^2(a=1)=1 (as long as radiation is not important at a=1)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Hubble2_norad=Hubble2(a,cosm)-cosm%Om_r*a**(-4)
    
  END FUNCTION Hubble2_norad

  REAL FUNCTION AH(a,cosm)

    ! Acceleration function: \ddot{a}/a
    ! Some people call this the Raychaudhuri equation
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'AH: Error, cosmology is not initialised'
    !AH=cosm%Om_m*a**(-3)+2.*cosm%Om_r*a**(-4)-2.*cosm%Om_v_mod+cosm%Om_w*(1.+3.*w_de(a,cosm))*X_de(a,cosm)
    !AH=-AH/2.
    AH=&
         cosm%Om_c*(1.+3*w_c)*X_c(a)+&
         cosm%Om_b*(1.+3*w_b)*X_b(a)+&
         cosm%Om_g*(1.+3*w_g)*X_g(a)+&
         Omega_nu_0(a,cosm)*(1.+3*w_nu(a,cosm))*X_nu(a,cosm)+&
         cosm%Om_v*(1.+3*w_v)*X_v(a)+&
         cosm%Om_w*(1.+3.*w_de(a,cosm))*X_de(a,cosm)
    AH=-AH/2.

  END FUNCTION AH

  REAL FUNCTION AH_norad(a,cosm)

    ! Acceleration function without the radiation contribution
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    AH_norad=AH(a,cosm)+cosm%Om_r*a**(-4)
    
  END FUNCTION AH_norad

  REAL FUNCTION Omega_m(a,cosm)

    ! This calculates Omega_m variations with a (note this is not proportional to a^-3 always)
    ! TODO: Retire
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !Omega_m=cosm%Om_m*a**(-3)/Hubble2(a,cosm)
    Omega_m=cosm%Om_m*X_m(a)/Hubble2(a,cosm)

  END FUNCTION Omega_m

!!$  REAL FUNCTION Omega_m_norad(a,cosm)
!!$
!!$    ! This calculates Omega_m variations with a, but ignoring any radiation component
!!$    ! This ensures that Omega_m_norad(a->0) -> 1
!!$    ! TODO: Retire
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_m_norad: Error, cosmology is not initialised'
!!$    !Omega_m_norad=cosm%Om_m*a**(-3)/Hubble2_norad(a,cosm)
!!$    Omega_m_norad=cosm%Om_m*X_m(a)/Hubble2_norad(a,cosm)
!!$
!!$  END FUNCTION Omega_m_norad

  REAL FUNCTION Omega_cold_norad(a,cosm)

    ! This calculates Omega_cold variations with a, but ignoring photon and neutrino components
    ! This ensures that Omega_m_norad(a->0) -> 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_cold_norad: Error, cosmology is not initialised'
    !Omega_cold_norad=(cosm%Om_c+cosm%Om_b)*a**(-3)/Hubble2_norad(a,cosm)
    Omega_cold_norad=(cosm%Om_c*X_c(a)+cosm%Om_b*X_b(a))/Hubble2_norad(a,cosm)

  END FUNCTION Omega_cold_norad

  REAL FUNCTION Omega_c(a,cosm)

    ! This calculates Omega_c variations with a (note this is not proportional to a^-3 always)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_c: Error, cosmology is not initialised'
    !Omega_c=cosm%Om_c*a**(-3)/Hubble2(a,cosm)
    Omega_c=cosm%Om_c*X_c(a)/Hubble2(a,cosm)

  END FUNCTION Omega_c

  REAL FUNCTION Omega_b(a,cosm)

    ! This calculates Omega_m variations with a (note this is not proportional to a^-3 always)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !Omega_b=cosm%Om_b*a**(-3)/Hubble2(a,cosm)
    Omega_b=cosm%Om_b*X_b(a)/Hubble2(a,cosm)

  END FUNCTION Omega_b

  REAL FUNCTION Omega_r(a,cosm)

    ! This calculates Omega_r variations over time (note this is *not* proportional to a^-4 always)
    ! TODO: Retire
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_r: Error, cosmology is not initialised'
    !Omega_r=cosm%Om_r*a**(-4)/Hubble2(a,cosm)
    Omega_r=cosm%Om_r*X_r(a)/Hubble2(a,cosm)

  END FUNCTION Omega_r

  REAL FUNCTION Omega_g(a,cosm)

    ! This calculates Omega_g variations over time (note this is *not* proportional to a^-4 always)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_g: Error, cosmology is not initialised'
    !Omega_r=cosm%Om_g*a**(-4)/Hubble2(a,cosm)
    Omega_g=cosm%Om_g*X_g(a)/Hubble2(a,cosm)

  END FUNCTION Omega_g

  REAL FUNCTION Omega_nu(a,cosm)

    ! This calculates Omega_nu variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

!!$    ! Scale-factor for transition between radiation and matter
!!$    IF(cosm%m_nu==0.) THEN
!!$       rad=.TRUE.
!!$    ELSE
!!$       anr=1./(1900.*cosm%m_nu)
!!$       IF(a>anr) THEN
!!$          rad=.FALSE.
!!$       ELSE
!!$          rad=.TRUE.
!!$       END IF
!!$    END IF
    
!!$    IF(a<cosm%a_nu) THEN
!!$       Omega_nu=cosm%Om_nu_rad*X_nu(a,cosm)/Hubble2(a,cosm)       
!!$    ELSE
!!$       Omega_nu=cosm%Om_nu*X_nu(a,cosm)/Hubble2(a,cosm)
!!$    END IF

    Omega_nu=Omega_nu_0(a,cosm)*X_nu(a,cosm)/Hubble2(a,cosm)

  END FUNCTION Omega_nu

  REAL FUNCTION Omega_nu_0(a,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(a>cosm%a_nu) THEN
       Omega_nu_0=cosm%Om_nu
    ELSE
       Omega_nu_0=cosm%Om_nu_rad
    END IF
    
  END FUNCTION Omega_nu_0

!!$  REAL FUNCTION Omega_nu(a,cosm)
!!$
!!$    ! This calculates Omega_nu variations over time
!!$    ! TODO: Fix this so that behavior transitions correctly
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_nu: Error, cosmology is not initialised'
!!$    Omega_nu=cosm%Om_nu*a**(-3)/Hubble2(a,cosm)
!!$
!!$  END FUNCTION Omega_nu

  REAL FUNCTION Omega_v(a,cosm)

    ! This calculates Omega_v variations with a (note this is not constant)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !Omega_v=cosm%Om_v_mod/Hubble2(a,cosm)
    !Omega_v=cosm%Om_v_mod*X_v(a)/Hubble2(a,cosm)
    Omega_v=cosm%Om_v*X_v(a)/Hubble2(a,cosm)

  END FUNCTION Omega_v

  REAL FUNCTION Omega_w(a,cosm)

    ! This calculates Omega_w variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA_w: Error, cosmology is not initialised'
    Omega_w=cosm%Om_w*X_de(a,cosm)/Hubble2(a,cosm)

  END FUNCTION Omega_w
 
  REAL FUNCTION Omega(a,cosm)

    ! This calculates total Omega variations with a
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA: Error, cosmology is not initialised'
    !Omega=Omega_m(a,cosm)+Omega_r(a,cosm)+Omega_v(a,cosm)+Omega_w(a,cosm)
    Omega=Omega_c(a,cosm)+Omega_b(a,cosm)+Omega_g(a,cosm)+Omega_nu(a,cosm)+Omega_v(a,cosm)+Omega_w(a,cosm)

  END FUNCTION Omega

!!$  REAL FUNCTION w_m(a)
!!$
!!$    ! Scaling of matter density
!!$    ! TODO: Retire
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    REAL :: crap
!!$
!!$    crap=a
!!$
!!$    w_m=0.
!!$    
!!$  END FUNCTION w_m
!!$
!!$  REAL FUNCTION w_c(a)
!!$
!!$    ! Scaling of CDM density
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    REAL :: crap
!!$
!!$    crap=a
!!$
!!$    w_c=0.
!!$    
!!$  END FUNCTION w_c
!!$
!!$  REAL FUNCTION w_b(a)
!!$
!!$    ! Scaling of baryon density
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    REAL :: crap
!!$
!!$    crap=a
!!$
!!$    w_b=0.
!!$
!!$  END FUNCTION w_b
!!$
!!$  REAL FUNCTION w_r(a)
!!$
!!$    ! Scaling of radiation density
!!$    ! TODO: Retire
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    REAL :: crap
!!$
!!$    crap=a
!!$    
!!$    w_r=1./3.
!!$    
!!$  END FUNCTION w_r
!!$
!!$  REAL FUNCTION w_g(a)
!!$
!!$    ! Scaling of photon density
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    REAL :: crap
!!$
!!$    crap=a
!!$
!!$    w_g=1./3.
!!$    
!!$  END FUNCTION w_g

  REAL FUNCTION w_nu(a,cosm)

    ! Scaling of neutrino density
    ! TODO: Account for radiation -> matter transition properly
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(a>cosm%a_nu) THEN
       w_nu=0.
    ELSE
       w_nu=1./3.
    END IF

  END FUNCTION w_nu

!!$  REAL FUNCTION w_v(a)
!!$
!!$    ! Scaling of vacuum density
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: a
!!$    REAL :: crap
!!$
!!$    crap=a
!!$
!!$    w_v=-1.
!!$    
!!$  END FUNCTION w_v

  REAL FUNCTION w_de(a,cosm)

    ! Variations of the dark energy parameter w(a)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: p1, p2, p3, p4
    DOUBLE PRECISION :: f1, f2, f3, f4

    IF(cosm%iw==1) THEN
       ! LCDM
       w_de=-1.
    ELSE IF(cosm%iw==2) THEN
       ! QUICC parameterisation
       p1=1.+exp(cosm%am/cosm%dm)
       p2=1.-exp(-(a-1.)/cosm%dm)
       p3=1.+exp(-(a-cosm%am)/cosm%dm)
       p4=1.-exp(1./cosm%dm)
       w_de=cosm%w+(cosm%wm-cosm%w)*p1*p2/(p3*p4)
    ELSE IF(cosm%iw==3) THEN
       ! w(a)CDM
       w_de=cosm%w+(1.-a)*cosm%wa
    ELSE IF(cosm%iw==4) THEN
       ! wCDM
       w_de=cosm%w
    ELSE IF(cosm%iw==5) THEN
       ! IDE I
       w_de=((a/cosm%as)**cosm%ns-1.)/((a/cosm%as)**cosm%ns+1.)
    ELSE IF(cosm%iw==6) THEN
       ! IDE II
       f1=a**cosm%ns-cosm%a1n
       f2=a**cosm%ns+cosm%a1n
       f3=a**cosm%ns-cosm%a2n
       f4=a**cosm%ns+cosm%a2n
       w_de=-1.+real(f1/f2-f3/f4)
    ELSE IF(cosm%iw==7) THEN
       ! IDE III
       IF(a<cosm%a1) THEN
          w_de=-1.
       ELSE IF(cosm%a1<=a .AND. a<cosm%a2) THEN
          w_de=cosm%ws
       ELSE IF(a>=cosm%a2) THEN
          w_de=-1.
       ELSE
          STOP 'W_DE: Error, something went wrong'
       END IF
    ELSE
       STOP 'W_DE: Error, value of iw set incorrectly'
    END IF

  END FUNCTION w_de

  REAL FUNCTION w_de_total(a,cosm)

    ! Do an average over the dark energy components
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%Om_v==0. .AND. cosm%Om_w==0.) THEN
       w_de_total=-1.
    ELSE
       w_de_total=w_de(a,cosm)*Omega_w(a,cosm)-Omega_v(a,cosm)
       w_de_total=w_de_total/(Omega_w(a,cosm)+Omega_v(a,cosm))
    END IF

  END FUNCTION w_de_total

  REAL FUNCTION w_eff(a,cosm)

    ! Do an average over all components
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    w_eff=&
         w_de(a,cosm)*Omega_w(a,cosm)+&
         w_v*Omega_v(a,cosm)+&
         w_g*Omega_g(a,cosm)+&
         w_nu(a,cosm)*Omega_nu(a,cosm)
    w_eff=w_eff/Omega(a,cosm)

  END FUNCTION w_eff

  REAL FUNCTION X_m(a)

    ! Scaling of matter density
    ! TODO: Retire
    IMPLICIT NONE
    REAL, INTENT(IN) :: a

    X_m=a**(-3)
    
  END FUNCTION X_m

  REAL FUNCTION X_c(a)

    ! Scaling of CDM density
    IMPLICIT NONE
    REAL, INTENT(IN) :: a

    X_c=a**(-3)
    
  END FUNCTION X_c

  REAL FUNCTION X_b(a)

    ! Scaling of baryon density
    IMPLICIT NONE
    REAL, INTENT(IN) :: a

    X_b=a**(-3)

  END FUNCTION X_b

  REAL FUNCTION X_r(a)

    ! Scaling of radiation density
    ! TODO: Retire
    IMPLICIT NONE
    REAL, INTENT(IN) :: a

    X_r=a**(-4)
    
  END FUNCTION X_r

  REAL FUNCTION X_g(a)

    ! Scaling of photon density
    IMPLICIT NONE
    REAL, INTENT(IN) :: a

    X_g=a**(-4)
    
  END FUNCTION X_g

  REAL FUNCTION X_nu(a,cosm)

    ! Scaling of neutrino density
    ! TODO: Account for radiation -> matter transition properly
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(a>cosm%a_nu) THEN
       X_nu=a**(-3)
    ELSE
       X_nu=a**(-4)
    END IF

  END FUNCTION X_nu

  REAL FUNCTION X_v(a)

    ! Scaling of vacuum density
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    REAL :: crap

    crap=a

    X_v=1.
    
  END FUNCTION X_v

  REAL FUNCTION X_de(a,cosm)

    ! Scaling for dark energy density (i.e., if w=0 x(a)=a^-3, if w=-1 x(a)=const etc.)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    DOUBLE PRECISION :: f1, f2, f3, f4

    IF(cosm%iw==1) THEN
       ! LCDM
       X_de=1.
    ELSE IF(cosm%iw==3) THEN
       ! w(a)CDM
       X_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))
    ELSE IF(cosm%iw==4) THEN
       ! wCDM
       X_de=a**(-3.*(1.+cosm%w))
    ELSE IF(cosm%iw==5) THEN
       ! IDE I
       X_de=((1.+(a/cosm%as)**cosm%ns)/(1.+(1./cosm%as)**cosm%ns))**(-6./cosm%ns)
    ELSE IF(cosm%iw==6) THEN
       ! IDE II
       f1=a**cosm%ns+cosm%a1n
       f2=1.+cosm%a1n
       f3=1.+cosm%a2n
       f4=a**cosm%ns+cosm%a2n
       X_de=real(f1*f3/(f2*f4))**(-6./cosm%ns)
    ELSE IF(cosm%iw==7) THEN
       ! IDE III
       IF(a<cosm%a1) THEN
          X_de=(cosm%a1/cosm%a2)**(-3.*(1.+cosm%ws))
       ELSE IF(cosm%a1<=a .AND. a<cosm%a2) THEN
          X_de=(a/cosm%a2)**(-3.*(1.+cosm%ws))
       ELSE IF(a>=cosm%a2) THEN
          X_de=1.
       ELSE
          STOP 'X_DE: Error, something went wrong'
       END IF
    ELSE
       ! Generally true, doing this integration can make calculations very slow
       STOP 'X_DE: Error, this integration routine has not been tested'
       X_de=(a**(-3))*exp(3.*integrate_cosm(a,1.,integrand_de,cosm,acc_cosm,3))
    END IF

  END FUNCTION X_de

  REAL FUNCTION integrand_de(a,cosm)

    ! The integrand for the X_de(a) integral
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    integrand_de=w_de(a,cosm)/a

  END FUNCTION integrand_de

  REAL FUNCTION scale_factor_z(z)

    ! The scale factor corresponding to redshift 'z'
    IMPLICIT NONE
    REAL, INTENT(IN) :: z

    IF(z<0.) THEN
       WRITE(*,*) 'SCALE_FACTOR_Z: z:', z
       STOP 'SCALE_FACTOR_Z: Error, routine called for z<0'
    END IF

    scale_factor_z=1./(1.+z)

  END FUNCTION scale_factor_z

  REAL FUNCTION redshift_a(a)

    ! The redshift corresponding to scale-factor 'a'
    IMPLICIT NONE
    REAL, INTENT(IN) :: a

    IF(a==0.) THEN
       WRITE(*,*) 'REDSHIFT_A: a:', a
       STOP 'REDSHIFT_A: Error, routine called with a=0'
    ELSE IF(a>1.) THEN
       WRITE(*,*) 'REDSHIFT_A: a:', a
       STOP 'REDSHIFT_A: Error, routine called with a>1'
    END IF

    redshift_a=-1.+1./a

  END FUNCTION redshift_a

  REAL FUNCTION scale_factor_r(r,cosm)

    ! The scale factor corresponding to comoving distance 'r'
    IMPLICIT NONE
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: p  

    IF(cosm%has_distance .EQV. .FALSE.) CALL init_distance(cosm)
    IF(r==0.) THEN
       scale_factor_r=1.
    ELSE
       p=cosm%horizon-r
       scale_factor_r=exp(find(log(p),cosm%log_p,cosm%log_a_p,cosm%n_p,3,3,2))
    END IF

  END FUNCTION scale_factor_r

  REAL FUNCTION redshift_r(r,cosm)

    ! The redshift corresponding to comoving distance 'r'
    IMPLICIT NONE
    REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
    TYPE(cosmology), INTENT(INOUT) :: cosm

    redshift_r=redshift_a(scale_factor_r(r,cosm))

  END FUNCTION redshift_r

  REAL FUNCTION f_k(r,cosm)

    ! Curvature function, also comoving angular-diameter distance [Mpc/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%k==0.) THEN
       f_k=r
    ELSE IF(cosm%k<0.) THEN
       f_k=sinh(sqrt(-cosm%k)*r)/sqrt(-cosm%k)
    ELSE IF(cosm%k>0.) THEN
       f_k=sin(sqrt(cosm%k)*r)/sqrt(cosm%k)
    ELSE
       STOP 'F_K: Something went wrong'
    END IF

  END FUNCTION f_k

  REAL FUNCTION fdash_k(r,cosm)

    ! Derivative of curvature function df_k(r)/dr
    IMPLICIT NONE
    REAL, INTENT(IN) :: r ! Comoving distance [Mpc/h]
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%k==0.) THEN
       fdash_k=1.
    ELSE IF(cosm%k<0.) THEN
       fdash_k=cosh(sqrt(-cosm%k)*r)
    ELSE IF(cosm%k>0.) THEN
       fdash_k=cos(sqrt(cosm%k)*r)
    ELSE
       STOP 'FDASH_K: Something went wrong'
    END IF

  END FUNCTION fdash_k

  REAL FUNCTION comoving_particle_horizon(a,cosm)

    ! The comoving particle horizon [Mpc/h]
    ! This is the furthest distance a particle can have travelled since a=0
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_distance .EQV. .FALSE.) CALL init_distance(cosm)

    IF(a==0.) THEN
       comoving_particle_horizon=0.
    ELSE IF(a>1.) THEN
       WRITE(*,*) 'COMOVING_PARTICLE_HORIZON: a:', a
       STOP 'COMOVING_PARTICLE_HORIZON: Error, tried to calculate particle horizon in the future'
    ELSE       
       comoving_particle_horizon=exp(find(log(a),cosm%log_a_p,cosm%log_p,cosm%n_p,3,3,2))
    END IF

  END FUNCTION comoving_particle_horizon

  REAL FUNCTION physical_particle_horizon(a,cosm)

    ! The physical particle horizon [Mpc/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_particle_horizon=comoving_particle_horizon(a,cosm)*a

  END FUNCTION physical_particle_horizon

  REAL FUNCTION comoving_distance(a,cosm)

    ! The comoving distance [Mpc/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: p

    p=comoving_particle_horizon(a,cosm) ! Ensures that init_distance is run and therefore that horizon is calculated
    comoving_distance=cosm%horizon-p

  END FUNCTION comoving_distance

  REAL FUNCTION physical_distance(a,cosm)

    ! The physical distance [Mpc/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_distance=comoving_distance(a,cosm)*a

  END FUNCTION physical_distance

  REAL FUNCTION physical_angular_distance(a,cosm)

    ! The physical angular-diameter distance [Mpc/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_angular_distance=a*comoving_angular_distance(a,cosm)

  END FUNCTION physical_angular_distance

  REAL FUNCTION comoving_angular_distance(a,cosm)

    ! The comoving angular-diameter distance [Mpc/h]
    ! Some people call this the 'effective distance'
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    comoving_angular_distance=f_k(comoving_distance(a,cosm),cosm)

  END FUNCTION comoving_angular_distance

  REAL FUNCTION luminosity_distance(a,cosm)

    ! The luminosity distance [Mpc/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    luminosity_distance=f_k(comoving_distance(a,cosm),cosm)/a

  END FUNCTION luminosity_distance

  SUBROUTINE init_distance(cosm)

    ! Fill up tables of a vs. p(a) (comoving particle horizon)
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: zmin, zmax, amin, amax, a, b, r
    INTEGER :: i
    INTEGER, PARAMETER :: iorder=3 ! Order for integration

    ! Get from PARAMETERS
    amin=amin_distance
    amax=amax_distance
    cosm%n_p=n_distance

    ! Calculate redshifts
    zmin=redshift_a(amax)
    zmax=redshift_a(amin)
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_DISTANCE: Redshift range for distance tables'
       WRITE(*,*) 'INIT_DISTANCE: minimum z:', real(zmin)
       WRITE(*,*) 'INIT_DISTANCE: maximum z:', real(zmax)
       WRITE(*,*) 'INIT_DISTANCE: minimum a:', real(amin)
       WRITE(*,*) 'INIT_DISTANCE: maximum a:', real(amax)
    END IF

    ! Fill array of 'a' in log space
    CALL fill_array(log(amin),log(amax),cosm%log_a_p,cosm%n_p)
    IF(ALLOCATED(cosm%log_p)) DEALLOCATE(cosm%log_p)
    ALLOCATE(cosm%log_p(cosm%n_p))

    ! Now do the r(a) calculation
    DO i=1,cosm%n_p
       a=exp(cosm%log_a_p(i))
       b=sqrt(a) ! Parameter to make the integrand not diverge for small values (a=b^2)
       r=integrate_cosm(0.,b,distance_integrand,cosm,acc_cosm,iorder)
       cosm%log_p(i)=log(r)
    END DO
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_DISTANCE: minimum r [Mpc/h]:', real(exp(cosm%log_p(1)))
       WRITE(*,*) 'INIT_DISTANCE: maximum r [Mpc/h]:', real(exp(cosm%log_p(cosm%n_p)))
    END IF

    ! Find the horizon distance in your cosmology
    ! exp(log) ensures the value is the same as what comes out of the (log) look-up tables
    cosm%horizon=exp(log(integrate_cosm(0.,1.,distance_integrand,cosm,acc_cosm,iorder))) 
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_DISTANCE: Horizon distance [Mpc/h]:', real(cosm%horizon)
       WRITE(*,*) 'INIT_DISTANCE: Done'
       WRITE(*,*)
    END IF

    cosm%has_distance=.TRUE.

  END SUBROUTINE init_distance

  REAL FUNCTION distance_integrand(b,cosm)

    ! The integrand for the cosmic-distance calculation [Mpc/h]
    ! Cast in terms of b=sqrt(a) to remove integrand divergence at a=0
    ! This means that the integrand is 2b/H(a)a^2, rather than 1/H(a)a^2
    IMPLICIT NONE
    REAL, INTENT(IN) :: b
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: a

    ! Relation beween a and b
    a=b**2

    IF(a<atay_distance) THEN
       IF(cosm%Om_r==0.) THEN
          distance_integrand=2.*Hdist/sqrt(cosm%Om_m)
       ELSE
          distance_integrand=2.*Hdist*b*(1.-0.5*(cosm%Om_m/cosm%Om_r)*b**2)/sqrt(cosm%Om_r)
       END IF
    ELSE
       distance_integrand=2.*Hdist*b/(sqrt(Hubble2(a,cosm))*a**2)
    END IF

  END FUNCTION distance_integrand

!!$  REAL FUNCTION age_of_universe(cosm)
!!$
!!$    ! The total age of the universe [Gyr/h]
!!$    IMPLICIT NONE
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    !age_of_universe=cosmic_time(1.,cosm)
!!$    age_of_universe=cosm%age
!!$
!!$  END FUNCTION age_of_universe

  REAL FUNCTION cosmic_time(a,cosm)

    ! The age of the universe [Gyr/h]
    ! TODO: Look-up table for t(a)?
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(.NOT. cosm%has_time) CALL init_time(cosm)

    !cosmic_time=integrate_cosm(0.,a,time_integrand,cosm,acc_cosm,3)
    IF(a==0.) THEN
       cosmic_time=0.
    ELSE
       cosmic_time=exp(find(log(a),cosm%log_a_t,cosm%log_t,cosm%n_p,3,3,2))
    END IF

  END FUNCTION cosmic_time

  REAL FUNCTION look_back_time(a,cosm)

    ! The time in the past [Gyr/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: t

    !look_back_time=integrate_cosm(a,1.,time_integrand,cosm,acc_cosm,3)
    t=cosmic_time(a,cosm) ! Ensures that init_time is run and therefore that age is calculated
    look_back_time=cosm%age-t

  END FUNCTION look_back_time

  SUBROUTINE init_time(cosm)

    ! Fill up tables of a vs. r(a) (comoving particle horizon)
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: zmin, zmax, amin, amax, a, t
    INTEGER :: i
    INTEGER, PARAMETER :: iorder=3

    ! Get from PARAMETERS
    amin=amin_time
    amax=amax_time
    cosm%n_t=n_time

    ! Calculate redshifts
    zmin=redshift_a(amax)
    zmax=redshift_a(amin)
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_TIME: Redshift range for time tables'
       WRITE(*,*) 'INIT_TIME: minimum z:', real(zmin)
       WRITE(*,*) 'INIT_TIME: maximum z:', real(zmax)
       WRITE(*,*) 'INIT_TIME: minimum a:', real(amin)
       WRITE(*,*) 'INIT_TIME: maximum a:', real(amax)
    END IF

    ! Fill array of 'a' in log space
    CALL fill_array(log(amin),log(amax),cosm%log_a_t,cosm%n_t)
    IF(ALLOCATED(cosm%log_t)) DEALLOCATE(cosm%log_t)
    ALLOCATE(cosm%log_t(cosm%n_t))

    ! Now do the r(a) calculation
    DO i=1,cosm%n_t
       a=exp(cosm%log_a_t(i))
       t=integrate_cosm(0.,a,time_integrand,cosm,acc_cosm,iorder)
       cosm%log_t(i)=log(t)
    END DO
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_TIME: minimum t [Gyr/h]:', real(exp(cosm%log_t(1)))
       WRITE(*,*) 'INIT_TIME: maximum t [Gyr/h]:', real(exp(cosm%log_t(cosm%n_t)))
    END IF

    ! Find the horizon distance in your cosmology
    ! exp(log) ensures the value is the same as what comes out of the (log) look-up tables
    cosm%age=exp(log(integrate_cosm(0.,1.,time_integrand,cosm,acc_cosm,iorder))) 
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_TIME: Age [Gyr/h]:', real(cosm%age)
       WRITE(*,*) 'INIT_TIME: Done'
       WRITE(*,*)
    END IF

    cosm%has_time=.TRUE.

  END SUBROUTINE init_time

  REAL FUNCTION time_integrand(a,cosm)

    ! The integrand for the cosmic-distance calculation [Gyr/h]
    ! This is 1/H(a)*a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(a<atay_time) THEN
       ! Taylor expansions
       IF(cosm%Om_r==0.) THEN
          time_integrand=Htime*sqrt(a/cosm%Om_m)
       ELSE
          time_integrand=Htime*a*(1.-0.5*cosm%Om_m*a/cosm%Om_r)/sqrt(cosm%Om_r)
       END IF
    ELSE
       time_integrand=Htime/(a*sqrt(Hubble2(a,cosm)))
    END IF

  END FUNCTION time_integrand

  REAL FUNCTION Tk(k,cosm)

    ! Transfer function selection
    IMPLICIT NONE
    REAL :: k ! Wavenumber [h/Mpc]
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%itk==1) THEN
       Tk=Tk_eh(k,cosm)
    ELSE IF(cosm%itk==3) THEN
       Tk=Tk_DEFW(k,cosm)
    ELSE
       WRITE(*,*) 'TK: itk:', cosm%itk
       STOP 'TK: Error, itk specified incorrectly'
    END IF

    ! Damp transfer function if considering WDM
    IF(cosm%inv_m_wdm .NE. 0.) Tk=Tk*Tk_wdm(k,cosm)

  END FUNCTION Tk

  REAL FUNCTION Tk_DEFW(k,cosm)

    ! The DEFW transfer function approximation
    ! Relies on the power-spectrum scale parameter Gamma=Omega_m*h
    ! This function was written by John Peacock   
    IMPLICIT NONE
    REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: keff, q, tk4
    DOUBLE PRECISION :: q8, tk8

    keff=0.172+0.011*log(cosm%Gamma/0.36)*log(cosm%Gamma/0.36)
    q=1.e-20 + k/cosm%Gamma
    q8=1.e-20 + keff/cosm%Gamma
    tk4=1./(1.+(6.4*q+(3.0*q)**1.5+(1.7*q)**2)**1.13)**(1./1.13)
    tk8=1./(1.+(6.4*q8+(3.0*q8)**1.5+(1.7*q8)**2)**1.13)**(1./1.13)

    tk_defw=tk4/real(tk8)

  END FUNCTION Tk_DEFW

  REAL FUNCTION Tk_EH(k,cosm)

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
    IF(cosm%power_Omegas) THEN
       Om_m=cosm%Om_m_pow
       Om_b=cosm%Om_b_pow
       h=cosm%h_pow
    ELSE
       Om_m=cosm%Om_m
       Om_b=cosm%Om_b 
       h=cosm%h
    END IF

    ! Physical densities
    Omh2=Om_m*h**2
    Obh2=Om_b*h**2

    ! Wave-number
    rk=k*h ! Convert to [1/Mpc]

    ! 2.718...
    e=exp(1.)

    ! CMB temperature (Section 2)
    thet=cosm%T_CMB/2.7
    
    b1=0.313*(Omh2)**(-0.419)*(1.+0.607*(Omh2)**0.674) ! Equation (4)
    b2=0.238*(Omh2)**0.223 ! Equation (4)
    zd=1291.*(1.+b1*(Obh2)**b2)*(Omh2)**0.251/(1.+0.659*(Omh2)**0.828) ! Drag redshift; equation (4)
    ze=2.50e4*Omh2/thet**4 ! z_eq; equation (2)
    rd=31500.*Obh2/thet**4/(1.+zd) ! Equation (5); changed from /zd -> /(1+zd) because it lines up with http://background.uchicago.edu/~whu/transfer/tf_fit.c (thanks Steven Murray)
    re=31500.*Obh2/thet**4/ze ! Equation (5)
    rke=7.46e-2*Omh2/thet**2 ! k_eq; equation (3)
    s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1.+sqrt(re))) ! Sound horizon at drag; equation (6)
    rks=1.6*( (Obh2)**0.52 ) * ( (Omh2)**0.73 ) * (1.+(10.4*Omh2)**(-0.95)) ! Silk k; equation(7)

    q=rk/13.41/rke ! Equation (10)

    y=(1.+ze)/(1.+zd) ! y that enters in equation G(y) in equations (14) and (15)
    g=y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.))) ! Equation (15)
    ab=g*2.07*rke*s/(1.+rd)**(0.75) ! Equation (14)

    a1=(46.9*Omh2)**0.670*(1+(32.1*Omh2)**(-0.532)) ! Equation (11)
    a2=(12.0*Omh2)**0.424*(1+(45.0*Omh2)**(-0.582)) ! Equation (11)
    ac=(a1**(-Om_b/Om_m)) * (a2**(-(Om_b/Om_m)**3)) ! Equation (11)

    b1=0.944/(1.+(458.*Omh2)**(-0.708)) ! Equation (12)
    b2=(0.395*Omh2)**(-0.0266) ! Equation (12)
    bc=1./(1.+b1*((1.-Om_b/Om_m)**b2-1.)) ! Equation (12)

    f=1./(1.+(rk*s/5.4)**4) ! Equation (18)

    c1=14.2 + 386./(1.+69.9*q**1.08) ! Equation (20) without alpha_c in the 14.2/alpha_c first bit
    c2=14.2/ac + 386./(1.+69.9*q**1.08) ! Equation (20) (C function should have explicity alpha_c dependence in paper)
    tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q**2) +(1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q**2) ! Equation (17)

    bb=0.5+(Om_b/Om_m) + (3.-2.*Om_b/Om_m)*sqrt((17.2*Omh2)**2+1.) ! Equation (24)
    bn=8.41*(Omh2)**0.435 ! Equation (23)
    ss=s/(1.+(bn/rk/s)**3)**(1./3.) ! Equation (22)
    tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q**2)/(1.+(rk*s/5.2)**2) ! First term in equation (21)

    ! Removed this IF statement as it produced a discontinuity in P_lin(k) as cosmology
    ! was varied - thanks David Copeland for pointing this out
    !IF((rk/rks**1.4)>7.) THEN
    !   fac=0.
    !ELSE
    fac=exp(-(rk/rks)**1.4) ! Silk-damping factor from equation (21)
    !END IF
    
    !tb=(tb+ab*fac/(1.+(bb/rk/s)**3))*sin(rk*ss)/rk/ss ! Equation (21)
    tb=(tb+ab*fac/(1.+(bb/rk/s)**3))*sinc(rk*ss) ! Equation (21)

    tk_eh=real((Om_b/Om_m)*tb+(1.-Om_b/Om_m)*tc) ! The weighted mean of baryon and CDM transfer functions

  END FUNCTION Tk_EH

  REAL FUNCTION Tk_WDM(k,cosm)

    ! Warm dark matter 'correction' to the standard transfer function
    ! This version and equation references were taken from arxiv:1605.05973
    ! Originally from Bode et al. (2001; arixv:0010389)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: alpha, mu, m_wdm

    m_wdm=1./cosm%inv_m_wdm

    alpha=0.074*0.7*m_wdm**(-1.15) ! alpha from equation (5), units Mpc/h
    mu=1.12                        ! mu from equation (4), dimensionless

    Tk_wdm=(1.+(alpha*k)**(2.*mu))**(-5./mu) ! Equation (2)
    
  END FUNCTION Tk_WDM

REAL FUNCTION Tcold(k,a,cosm)

! Ratio of transfer function for cold matter relative to all matter
! TODO: Add exact result from CAMB
IMPLICIT NONE
REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
REAL, INTENT(IN) :: a ! Scale factor
TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
INTEGER, PARAMETER :: method=3

IF(method==1) THEN
   Tcold=1.
ELSE IF(method==2) THEN
   Tcold=Tcold_approx(cosm)
ELSE IF(method==3) THEN
   Tcold=Tcold_ratio(k,a,cosm)
ELSE
   STOP 'TCOLD: Error, method not specified correctly'
END IF

END FUNCTION Tcold

REAL FUNCTION Tcold_approx(cosm)

! Approximation for how power is suppressed by massive nu at small scales
! Calculated assuming perturbation grow from z~1000 and that neutrinos are hot and therefore completely smooth
! Related to the growth-function approximation: g(a) = a^(1-3f_nu/5)
IMPLICIT NONE
TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

Tcold_approx=sqrt(1.-8.*cosm%f_nu)

END FUNCTION Tcold_approx

REAL FUNCTION Tcold_ratio(k,a,cosm)

   ! Calculates the ratio of T(k) for cold vs. all matter
   ! Uses approximations in Eisenstein & Hu (1999; arXiv 9710252)
   ! Note that this assumes that there are exactly 3 species of neutrinos with
   ! Nnu<=3 of these being massive, and with the mass split evenly between the number of massive species.
   IMPLICIT NONE
   REAL, INTENT(IN) :: k ! Wavenumber [h/Mpc]
   REAL, INTENT(IN) :: a ! Scale factor
   TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
   REAL :: D, Dcb, Dcbnu, pcb, zeq, q, yfs, z
   REAL :: BigT
   INTEGER, PARAMETER :: N_massive_nu=3

   IF(cosm%f_nu==0.) THEN

      Tcold_ratio=1.

   ELSE

      ! Get the redshift
      z=redshift_a(a)

      ! Growth exponent under the assumption that neutrinos are completely unclustered (equation 11)
      pcb=(5.-sqrt(1.+24.*(1.-cosm%f_nu)))/4.

      ! Theta for temperature (BigT=T/2.7 K)
      BigT=cosm%T_CMB/2.7

      ! The matter-radiation equality redshift
      zeq=(2.5e4)*cosm%om_m*(cosm%h**2.)*(BigT**(-4.))

      ! The growth function normalised such that D=(1.+z_eq)/(1+z) at early times (when Omega_m \approx 1)
      ! For my purpose (just the ratio) seems to work better using the EdS growth function result, \propto a .
      ! In any case, can't use grow at the moment because that is normalised by default.
      D=(1.+zeq)/(1.+z) ! TODO: Could update this

      ! Wave number relative to the horizon scale at equality (equation 5)
      ! Extra factor of h becauase all my k are in units of h/Mpc
      q=k*cosm%h*BigT**2./(cosm%om_m*cosm%h**2.)

      ! Free streaming scale (equation 14)
      ! Note that Eisenstein & Hu (1999) only consider the case of 3 neutrinos
      ! with Nnu of these being massive with the mass split evenly between Nnu species.
      yfs=17.2*cosm%f_nu*(1.+0.488*cosm%f_nu**(-7./6.))*(real(N_massive_nu)*q/cosm%f_nu)**2.

      ! These are (almost) the scale-dependent growth functions for each component in Eisenstein & Hu (1999)
      ! Some part is missing, but this cancels when they are divided by each other, which is all I need them for.
      ! Equations (12) and (13)
      Dcb=(1.+(D/(1.+yfs))**0.7)**(pcb/0.7)
      Dcbnu=((1.-cosm%f_nu)**(0.7/pcb)+(D/(1.+yfs))**0.7)**(pcb/0.7)

      ! Finally, the ratio
      Tcold_ratio=Dcb/Dcbnu

   END IF

   END FUNCTION Tcold_ratio

  REAL RECURSIVE FUNCTION p_lin(k,a,cosm)

    ! Linear matter power spectrum
    ! Must be a recursive function because normalise_power calls this function again
    ! CHANGE: Changed RECURSIVE FUNCTION to REAL RECURSIVE FUNCTION
    IMPLICIT NONE
    REAL, INTENT (IN) :: k, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, PARAMETER :: iorder=3 ! Order for interpolation (3 - cubic)
    INTEGER, PARAMETER :: ifind=3  ! Finding scheme in table (3 - Mid-point method)
    INTEGER, PARAMETER :: imeth=2  ! Method for polynomials (2 - Lagrange polynomials)

    ! Using init_power seems to provide no significant speed improvements to HMx
    !IF(cosm%has_power .EQV. .FALSE.) CALL init_power(cosm)

    ! This line generates a recursion
    IF(.NOT. cosm%is_normalised) THEN
      cosm%is_normalised = .TRUE. ! This line is not necessary
      CALL normalise_power(cosm)
    END IF

    IF(k<=kmin_plin) THEN
       ! If p_lin happens to be foolishly called for 0 mode
       ! This call should never happen, but may in integrals
       p_lin=0.
    ELSE IF(k>kmax_plin) THEN
       ! Avoids some issues if p_lin is called for very (absurdly) high k values
       ! For some reason crashes can occur if this is the case
       p_lin=0.
    ELSE IF(cosm%box .AND. k<cosm%kbox) THEN
       ! If investigating effects caused by a finite box size
       p_lin=0.
    ELSE
       IF(cosm%has_power) THEN
          ! TODO: Do something cleverer here. Could use the ln(k)^2 behaviour at high k, could just truncate...
          IF(cosm%growk) THEN
             p_lin=exp(find(log(k),cosm%log_k_plin,log(a),cosm%log_a_plin,cosm%log_plina,cosm%nk_plin,cosm%na_plin,iorder,ifind,imeth=1))
          ELSE
             p_lin=(grow(a,cosm)**2)*exp(find(log(k),cosm%log_k_plin,cosm%log_plin,cosm%nk_plin,iorder,ifind,imeth))
          END IF
          !p_lin=exp(find(log(k),cosm%log_k_plin,cosm%log_plin,cosm%n_plin,iorder,ifind,imeth))
       ELSE
          ! In this case get the power from the transfer function
          p_lin=(cosm%A**2)*(grow(a,cosm)**2)*(Tk(k,cosm)**2)*(k**(cosm%n+3.))
       END IF
    END IF

  END FUNCTION p_lin

!!$  SUBROUTINE init_power(cosm)
!!$
!!$    ! Fill a look-up table for the linear power spectrum from a fitting function
!!$    IMPLICIT NONE
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    INTEGER :: i
!!$    REAL :: k
!!$
!!$    ! PARAMETERS
!!$    REAL, PARAMETER :: kmin=1e-3
!!$    REAL, PARAMETER :: kmax=1e2
!!$    INTEGER, PARAMETER :: n_plin=256
!!$
!!$    ! Set the numer of points in the look-up tables
!!$    ! Note you need to have enough to resolve the BAO well
!!$    ! Probably some non-log/linear spacing would be best (CAMB does this)
!!$    cosm%n_plin=n_plin
!!$
!!$    ! Allocate arrays
!!$    IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
!!$    IF(ALLOCATED(cosm%log_plin)) DEALLOCATE(cosm%log_plin)
!!$    ALLOCATE(cosm%log_k_plin(cosm%n_plin),cosm%log_plin(cosm%n_plin))
!!$
!!$    ! Get values for the linear power spectrum
!!$    DO i=1,cosm%n_plin
!!$       k=progression_log(kmin,kmax,i,cosm%n_plin)
!!$       cosm%log_k_plin(i)=k
!!$       cosm%log_plin(i)=(Tk(k,cosm)**2)*k**(cosm%n+3.)
!!$    END DO
!!$
!!$    ! Take logarithms
!!$    cosm%log_k_plin=log(cosm%log_k_plin)
!!$    cosm%log_plin=log(cosm%log_plin)
!!$
!!$    ! Change the flag to true
!!$    cosm%has_power=.TRUE.
!!$    
!!$  END SUBROUTINE init_power

  SUBROUTINE init_sigma(cosm)

    ! This fills up tables of r vs. sigma(r) across a range in r
    ! It is used only in look-up for further calculations of sigma(r) and not otherwise
    ! and prevents a large number of calls to the sigma integration functions
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: r, a, sig
    INTEGER :: i, j
    LOGICAL, PARAMETER :: cold=.TRUE.

    !IF(cosm%inv_m_wdm .NE. 0.) STOP 'INIT_SIGMA: This will crash with WDM'

    ! Deallocate tables if they are already allocated
    IF(ALLOCATED(cosm%log_r_sigma)) DEALLOCATE(cosm%log_r_sigma)
    IF(ALLOCATED(cosm%log_sigma))   DEALLOCATE(cosm%log_sigma)
    IF(ALLOCATED(cosm%log_a_sigma)) DEALLOCATE(cosm%log_a_sigma)
    IF(ALLOCATED(cosm%log_sigmaa))  DEALLOCATE(cosM%log_sigmaa)

    ! Write to screen
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_SIGMA: Filling sigma(R) interpolation table'
       WRITE(*,*) 'INIT_SIGMA: R minimum [Mpc/h]:', real(rmin_sigma)
       WRITE(*,*) 'INIT_SIGMA: R maximum [Mpc/h]:', real(rmax_sigma)
       WRITE(*,*) 'INIT_SIGMA: number of points:', nr_sigma
    END IF

    ! Allocate and fill array of R values
    CALL fill_array(log(rmin_sigma),log(rmax_sigma),cosm%log_r_sigma,nr_sigma)
    cosm%nr_sigma=nr_sigma

    IF(cosm%growk) THEN
      cosm%na_sigma=na_sigma
      CALL fill_array(log(amin_sigma),log(amax_sigma),cosm%log_a_sigma,na_sigma)
      ALLOCATE(cosm%log_sigmaa(nr_sigma,na_sigma))
      DO j=1,na_sigma
         a=exp(cosm%log_a_sigma(j))
         DO i=1,nr_sigma
            r=exp(cosm%log_r_sigma(i))
            sig=sigma_cold_integral(r,a,cosm)
            !sig=sigma_all_integral(r,a,cosm)
            cosm%log_sigmaa(i,j)=log(sig)
         END DO
      END DO
    ELSE
      ! Loop over R values and calculate sigma(R)
      ALLOCATE(cosm%log_sigma(nr_sigma))
      a=1. ! Fill this table at a=1
      DO i=1,nr_sigma
         !r=exp(progression(log(rmin_sigma),log(rmax_sigma),i,nsig))
         r=exp(cosm%log_r_sigma(i))
         sig=sigma_cold_integral(r,a,cosm)
         !sig=sigma_all_integral(r,a,cosm)
         cosm%log_sigma(i)=log(sig)
      END DO
   END IF

    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_SIGMA: Done'
       WRITE(*,*)
    END IF

    ! Change flag so that it is known that the look-up tables are filled
    cosm%has_sigma=.TRUE.

  END SUBROUTINE init_sigma

REAL FUNCTION sigma_all_integral(r,a,cosm)

   ! Calculates sigma(R) by intergration
   IMPLICIT NONE
   REAL, INTENT(IN) :: r ! Smoothing scale to calculate sigma [Mpc/h]
   REAL, INTENT(IN) :: a ! Scale factor
   TYPE(cosmology), INTENT(INOUT) :: cosm
   REAL :: tmin, tmax, kmin, kmax, part1, part2

   ! Integration method changes depending on r to make this as fast as possible
   IF(r>=Rsplit_sigma) THEN
      ! Integration upper limit (c = 1 corresponds to k = 0)
      tmin=0.
      tmax=1.
      sigma_all_integral=sqrt(integrate_cosm(tmin,tmax,sigma2_all_integrand_transformed,r,a,cosm,2.*acc_cosm,3))
   ELSE IF(r<Rsplit_sigma) THEN
      ! Integration limits, the split of the integral is done at k = 1/R
      tmin=t_sigma_integrand(ksplit_sigma(R),R)
      tmax=1.  ! Integration limit corresponding to k=0
      kmin=ksplit_sigma(R) ! From the split wavenumber...
      kmax=sigma_out/R     ! ...out to k = inf, but in practice just go out a finite distance in kR
      part1=integrate_cosm(tmin,tmax,sigma2_all_integrand_transformed,r,a,cosm,2.*acc_cosm,3)          
      part2=integrate_cosm(kmin,kmax,sigma2_all_integrand,r,a,cosm,2.*acc_cosm,3)
      sigma_all_integral=sqrt(part1+part2)
   ELSE
      STOP 'INIT_SIGMA: Error, something went wrong'
   END IF

END FUNCTION sigma_all_integral

REAL FUNCTION sigma_cold_integral(r,a,cosm)

   ! Calculates sigma(R) by intergration
   IMPLICIT NONE
   REAL, INTENT(IN) :: r ! Smoothing scale to calculate sigma [Mpc/h]
   REAL, INTENT(IN) :: a ! Scale factor
   TYPE(cosmology), INTENT(INOUT) :: cosm
   REAL :: tmin, tmax, kmin, kmax, part1, part2

   ! Integration method changes depending on r to make this as fast as possible
   IF(r>=Rsplit_sigma) THEN
      ! Integration upper limit (c = 1 corresponds to k = 0)
      tmin=0.
      tmax=1.
      sigma_cold_integral=sqrt(integrate_cosm(tmin,tmax,sigma2_cold_integrand_transformed,r,a,cosm,2.*acc_cosm,3))
   ELSE IF(r<Rsplit_sigma) THEN
      ! Integration limits, the split of the integral is done at k = 1/R
      tmin=t_sigma_integrand(ksplit_sigma(R),R)
      tmax=1.  ! Integration limit corresponding to k=0
      kmin=ksplit_sigma(R) ! From the split wavenumber...
      kmax=sigma_out/R     ! ...out to k = inf, but in practice just go out a finite distance in kR
      part1=integrate_cosm(tmin,tmax,sigma2_cold_integrand_transformed,r,a,cosm,2.*acc_cosm,3)          
      part2=integrate_cosm(kmin,kmax,sigma2_cold_integrand,r,a,cosm,2.*acc_cosm,3)
      sigma_cold_integral=sqrt(part1+part2)
   ELSE
      STOP 'INIT_SIGMA: Error, something went wrong'
   END IF

END FUNCTION sigma_cold_integral

REAL FUNCTION sigma_all(R,a,cosm)

! Gets sigma
! TODO: Could add look-up table

IMPLICIT NONE
REAL, INTENT(IN) :: R ! Smoothing scale to calculate sigma [Mpc/h]
REAL, INTENT(IN) :: a ! Scale factor
TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology

sigma_all=sigma_all_integral(R,a,cosm)

END FUNCTION sigma_all

REAL FUNCTION sigma_cold(R,a,cosm)

   ! Finds sigma_cold from look-up tables
   IMPLICIT NONE
   REAL, INTENT(IN) :: R ! Smoothing scale to calculate sigma [Mpc/h]
   REAL, INTENT(IN) :: a ! Scale factor
   TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
   INTEGER, PARAMETER :: iorder=3 ! Order for interpolation (3 - Cubic)
   INTEGER, PARAMETER :: ifind=3  ! Finding scheme in table (3 - Mid-point method)
   INTEGER, PARAMETER :: imeth=2  ! Method for polynomials (2 - Lagrange polynomails)
   
   IF(cosm%has_sigma .EQV. .FALSE.) CALL init_sigma(cosm)

   IF(cosm%growk) THEN
      sigma_cold=exp(find(log(R),cosm%log_r_sigma,log(a),cosm%log_a_sigma,cosm%log_sigmaa,cosm%nr_sigma,cosm%na_sigma,iorder,ifind,imeth=1))
   ELSE
      sigma_cold=grow(a,cosm)*exp(find(log(R),cosm%log_r_sigma,cosm%log_sigma,cosm%nr_sigma,iorder,ifind,imeth))
   END IF
   !sigma=grow(a,cosm)*exp(find(log(R),cosm%log_r_sigma,cosm%log_sigma,cosm%n_sigma,iorder,ifind,imeth))

END FUNCTION sigma_cold

REAL FUNCTION sigma2_all_integrand(k,R,a,cosm)

   ! The integrand for the sigma(R) integrals
   USE special_functions
   IMPLICIT NONE
   REAL, INTENT(IN) :: k
   REAL, INTENT(IN) :: R
   REAL, INTENT(IN) :: a
   TYPE(cosmology), INTENT(INOUT) :: cosm
   REAL :: w_hat

   IF(k==0.) THEN
      sigma2_all_integrand=0.
   ELSE
      w_hat=wk_tophat(k*R)       
      sigma2_all_integrand=p_lin(k,a,cosm)*(w_hat**2)/k
   END IF

END FUNCTION sigma2_all_integrand

REAL FUNCTION sigma2_cold_integrand(k,R,a,cosm)

IMPLICIT NONE
REAL, INTENT(IN) :: k
REAL, INTENT(IN) :: R
REAL, INTENT(IN) :: a
TYPE(cosmology), INTENT(INOUT) :: cosm

sigma2_cold_integrand=sigma2_all_integrand(k,R,a,cosm)*Tcold(k,a,cosm)**2

END FUNCTION sigma2_cold_integrand

REAL FUNCTION sigma2_all_integrand_transformed(t,R,a,cosm)

   ! The integrand for the sigma(R) integrals
   USE special_functions
   IMPLICIT NONE
   REAL, INTENT(IN) :: t
   REAL, INTENT(IN) :: R
   REAL, INTENT(IN) :: a
   TYPE(cosmology), INTENT(INOUT) :: cosm
   REAL :: k, kR, w_hat

   ! Integrand to the sigma integral in terms of t. Defined by kR=(1/t-1)**alpha
   ! alpha_sigma can be any positive number, can even be a function of R
   IF(t==0.) THEN
      ! t=0 corresponds to k=infintiy when W(kR)=0
      sigma2_all_integrand_transformed=0.
   ELSE IF(t==1.) THEN
      ! t=1 corresponds to k=0 when P(k)=0
      sigma2_all_integrand_transformed=0.
   ELSE
      k=k_sigma_integrand(t,R)
      kR=k*R
      w_hat=wk_tophat(kR)
      sigma2_all_integrand_transformed=p_lin(k,a,cosm)*(w_hat**2)*alpha_sigma/(t*(1.-t))
   END IF

END FUNCTION sigma2_all_integrand_transformed

REAL FUNCTION sigma2_cold_integrand_transformed(t,R,a,cosm)

  IMPLICIT NONE
  REAL, INTENT(IN) :: t
  REAL, INTENT(IN) :: R
  REAL, INTENT(IN) :: a
  TYPE(cosmology), INTENT(INOUT) :: cosm
  REAL :: k

  IF(t==0.) THEN
   sigma2_cold_integrand_transformed=0.
  ELSE IF(t==1.) THEN
   sigma2_cold_integrand_transformed=0.
  ELSE
   k=k_sigma_integrand(t,R)
   sigma2_cold_integrand_transformed=sigma2_all_integrand_transformed(t,R,a,cosm)*Tcold(k,a,cosm)**2
  END IF

END FUNCTION sigma2_cold_integrand_transformed

  REAL FUNCTION ksplit_sigma(R)

    ! Wavenumber at which to split the integral into two [h/Mpc]
    IMPLICIT NONE
    REAL, INTENT(IN) :: R

    ksplit_sigma=1./R
    
  END FUNCTION ksplit_sigma

  REAL FUNCTION k_sigma_integrand(t,R)

    ! How k relates to t and R in the transformed sigma^2(R) integrand
    ! kR=(-1+1/t)^alpha
    ! This function is the inverse of t_sigma_integrand
    IMPLICIT NONE
    REAL, INTENT(IN) :: t
    REAL, INTENT(IN) :: R
    REAL :: kR

    IF(t==0.) THEN
       STOP 'K_SIGMA_INTEGRAND: Error, this should not be called with t = 0'
    ELSE
       kR=(-1.+1./t)**alpha_sigma
       k_sigma_integrand=kR/R
    END IF
    
  END FUNCTION k_sigma_integrand

  REAL FUNCTION t_sigma_integrand(k,R)

    ! How t relates to k and R in the transformed sigma^2(R) integrand
    ! t = 1/[1+(kR)^(1/alpha)]
    ! This function is the inverse of sigma_integrand_k
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: R

    t_sigma_integrand=1./(1.+(k*R)**(1./alpha_sigma))
    
  END FUNCTION t_sigma_integrand

  REAL FUNCTION sigmaV(R,a,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: R
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, PARAMETER :: tmin=0.
    REAL, PARAMETER :: tmax=1.

    sigmaV=integrate_cosm(tmin,tmax,sigmaV2_integrand,R,a,cosm,2.*acc_cosm,3)

    ! Convert 3D sigmaV^2 to 1D sigmaV
    sigmaV=sqrt(sigmaV/3.)

  END FUNCTION sigmaV

  REAL FUNCTION sigmaV2_integrand(t,R,a,cosm)

    ! This is the integrand for the velocity dispersion integral
    ! TODO: Optimize alpha(R); not really a problem when only called for R = 0, 100 Mpc/h
   ! TODO: Should this use cold or all matter power spectrum?
    USE special_functions
    IMPLICIT NONE
    REAL, INTENT(IN) :: t
    REAL, INTENT(IN) :: R
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: k, kR, w_hat, alpha

    IF(t==0. .OR. t==1.) THEN
       ! t = 0 corresponds to k = infinity when W(kR) = 0
       ! t = 1 corresponds to k = 0 when P(k) = 0
       sigmaV2_integrand=0.
    ELSE
       alpha=alpha_sigmaV
       IF(R==0.) THEN
          kR=0.
          k=(-1.+1./t)**alpha
       ELSE
          kR=(-1.+1./t)**alpha
          k=kR/R          
       END IF
       w_hat=wk_tophat(kR)
       !sigmaV2_integrand=(Tcold_ratio(k,a,cosm)**2)*(p_lin(k,a,cosm)/k**2)*(w_hat**2)*alpha/(t*(1.-t))
       !sigmaV2_integrand=(Tcold_approx(cosm)**2)*(p_lin(k,a,cosm)/k**2)*(w_hat**2)*alpha/(t*(1.-t))
       sigmaV2_integrand=(p_lin(k,a,cosm)/k**2)*(w_hat**2)*alpha/(t*(1.-t))
    END IF

  END FUNCTION sigmaV2_integrand

  REAL FUNCTION grow(a,cosm)

    ! Scale-independent growth function | normalised g(z=0)=1
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
    IF(a==1.) THEN
       grow=1.
    ELSE       
       grow=exp(find(log(a),cosm%log_a_growth,cosm%log_growth,cosm%n_growth,3,3,2))
    END IF

  END FUNCTION grow

  REAL FUNCTION ungrow(a,cosm)

    ! Unnormalised growth function
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
   
    ungrow=cosm%gnorm*grow(a,cosm)

  END FUNCTION ungrow

  REAL FUNCTION growth_rate(a,cosm)

    ! Growth rate: dln(g) / dln(a)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)    
    growth_rate=find(log(a),cosm%log_a_growth,cosm%growth_rate,cosm%n_growth,3,3,2)

  END FUNCTION growth_rate

  REAL FUNCTION acc_growth(a,cosm)

    ! Accumulated growth function: int_0^a g(a)/a da
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)    
    acc_growth=exp(find(log(a),cosm%log_a_growth,cosm%log_acc_growth,cosm%n_growth,3,3,2))

  END FUNCTION acc_growth

  REAL FUNCTION growth_rate_Linder(a,cosm)

    ! Approximation for the growth rate from Linder xxxx.xxxx
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: gam

    IF(cosm%w<-1.) THEN
       gam=0.55+0.02*(1.+cosm%w)
    ELSE IF(cosm%w>-1) THEN
       gam=0.55+0.05*(1.+cosm%w)
    ELSE
       gam=0.55
    END IF

    !growth_rate_Linder=Omega_m_norad(a,cosm)**gam
    growth_rate_Linder=Omega_cold_norad(a,cosm)**gam
    
  END FUNCTION growth_rate_Linder

  REAL FUNCTION grow_Linder_integrand(a,cosm)

    ! Integrand for the approximate growth integral using Linder approximate growth rate
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    grow_Linder_integrand=growth_rate_Linder(a,cosm)/a

  END FUNCTION grow_Linder_integrand

  REAL FUNCTION grow_Linder(a,cosm)

    ! Calculate the growth function from the Linder growth rate via integration
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    grow_Linder=exp(-integrate_cosm(a,1.,grow_Linder_integrand,cosm,acc_cosm,3))
    
  END FUNCTION grow_Linder

  REAL FUNCTION grow_CPT(a,cosm)

    ! Carrol, Press & Turner (1992) approximation to growth function (good to 5%)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: Om_mz, Om_vz, Om_m, Om_v

    ! Get all necessary Omega values
    !Om_mz=Omega_m_norad(a,cosm)
    Om_mz=Omega_cold_norad(a,cosm)
    Om_vz=Omega_v(a,cosm)+Omega_w(a,cosm)
    Om_m=cosm%Om_b+cosm%Om_c
    Om_v=cosm%Om_v+cosm%Om_w

    ! Now call CPT twice, second time to normalise it
    grow_CPT=CPT(a,Om_mz,Om_vz)/CPT(1.,Om_m,Om_v)

  END FUNCTION grow_CPT

  REAL FUNCTION CPT(a,Om_m,Om_v)

    ! The CPT growth function approximation from 1992
    IMPLICIT NONE
    REAL, INTENT(IN) :: a    ! Scale factor
    REAL, INTENT(IN) :: Om_m ! Matter-density parameter
    REAL, INTENT(IN) :: Om_v ! Vacuum-density parameter

    CPT=a*Om_m/((Om_m**(4./7.))-Om_v+(1.+Om_m/2.)*(1.+Om_v/70.))
    
  END FUNCTION CPT

  SUBROUTINE init_growth(cosm)

    ! Fills a table of the growth function vs. a
    ! TODO: Figure out why if I set amax=10, rather than amax=1, I start getting weird f(a) around a=0.001
    USE calculus_table
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, na
    REAL :: a
    REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL :: dinit, vinit
    REAL :: g0, f0, bigG0   

    !! First do growth factor and growth rate !!
    
    ! These set the initial conditions to be the Omega_m(a)=1 growing mode
    !dinit=ainit_growth
    !vinit=1.

    ! Set the initial conditions to be in the cold matter growing mode
    dinit=ainit_growth**(1.-3.*cosm%f_nu/5.)
    vinit=(1.-3.*cosm%f_nu/5.)*ainit_growth**(-3.*cosm%f_nu/5.)
    
    ! Write some useful information to the screen
    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_GROWTH: Solving growth equation'
       WRITE(*,*) 'INIT_GROWTH: Minimum scale factor:', ainit_growth
       WRITE(*,*) 'INIT_GROWTH: Maximum scale factor:', amax_growth
       WRITE(*,*) 'INIT_GROWTH: Initial delta:', dinit
       WRITE(*,*) 'INIT_GROWTH: Initial delta derivative:', vinit
       WRITE(*,*) 'INIT_GROWTH: Number of points for look-up tables:', n_growth
    END IF

    ! Solve the ODE
    CALL ODE_adaptive_cosmology(d_tab,v_tab,0.,a_tab,cosm,ainit_growth,amax_growth,dinit,vinit,ddda,dvda,acc_cosm,3,.FALSE.)
    IF(cosm%verbose) WRITE(*,*) 'INIT_GROWTH: ODE done'
    na=SIZE(a_tab)

    ! Convert dv/da to f = dlng/dlna for later, so v_tab should really be f_tab from now on
    v_tab=v_tab*a_tab/d_tab

    ! Normalise so that g(z=0)=1
    cosm%gnorm=find(1.,a_tab,d_tab,na,3,3,2)
    IF(cosm%verbose) WRITE(*,*) 'INIT_GROWTH: unnormalised growth at z=0:', real(cosm%gnorm)
    d_tab=d_tab/cosm%gnorm   

    ! Allocate arrays
    IF(ALLOCATED(cosm%log_a_growth))    DEALLOCATE(cosm%log_a_growth)
    IF(ALLOCATED(cosm%log_growth))      DEALLOCATE(cosm%log_growth)
    IF(ALLOCATED(cosm%growth_rate))     DEALLOCATE(cosm%growth_rate)
    IF(ALLOCATED(cosm%log_acc_growth))  DEALLOCATE(cosm%log_acc_growth)
    cosm%n_growth=n_growth

    ! This downsamples the tables that come out of the ODE solver (which can be a bit long)
    ! Could use some table-interpolation routine here to save time
    ALLOCATE(cosm%log_a_growth(n_growth))
    ALLOCATE(cosm%log_growth(n_growth))
    ALLOCATE(cosm%growth_rate(n_growth))
    DO i=1,n_growth
       a=progression(ainit_growth,amax_growth,i,n_growth)
       cosm%log_a_growth(i)=a
       cosm%log_growth(i)=exp(find(log(a),log(a_tab),log(d_tab),na,3,3,2))
       cosm%growth_rate(i)=find(log(a),log(a_tab),v_tab,na,3,3,2)
    END DO

    !! !!
 
    !! Table integration to calculate G(a)=int_0^a g(a')/a' da' !!

    ! Allocate array
    ALLOCATE(cosm%log_acc_growth(n_growth))

    ! Set to zero, because I have an x=x+y thing later on
    cosm%log_acc_growth=0.
    
    ! Do the integral up to table position i, which fills the accumulated growth table
    DO i=1,n_growth

       ! Do the integral using the arrays
       IF(i>1) THEN
          cosm%log_acc_growth(i)=integrate_table(cosm%log_a_growth,cosm%gnorm*cosm%log_growth/cosm%log_a_growth,n_growth,1,i,3)
       END IF
       
       ! Then add on the section that is missing from the beginning
       ! NB. g(a=0)/0 = 1, so you just add on a rectangle of height g*a/a=g
       cosm%log_acc_growth(i)=cosm%log_acc_growth(i)+cosm%gnorm*cosm%log_growth(1)
       
    END DO

    !! !!

    ! Write stuff about growth parameter at a=1 to the screen    
    IF(cosm%verbose) THEN
       f0=find(1.,cosm%log_a_growth,cosm%growth_rate,n_growth,3,3,2)
       g0=find(1.,cosm%log_a_growth,cosm%log_growth,n_growth,3,3,2)
       bigG0=find(1.,cosm%log_a_growth,cosm%log_acc_growth,n_growth,3,3,2)      
       WRITE(*,*) 'INIT_GROWTH: normalised growth at z=0:', g0
       WRITE(*,*) 'INIT_GROWTH: growth rate at z=0:', f0
       WRITE(*,*) 'INIT_GROWTH: integrated growth at z=0:', bigG0
    END IF

    ! Make the some of the tables log for easier interpolation
    cosm%log_a_growth=log(cosm%log_a_growth)
    cosm%log_growth=log(cosm%log_growth)
    cosm%log_acc_growth=log(cosm%log_acc_growth)

    ! Set the flag to true so that this subroutine is only called once
    cosm%has_growth=.TRUE.

    IF(cosm%verbose) THEN
       WRITE(*,*) 'INIT_GROWTH: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE init_growth

  REAL FUNCTION ddda(d,v,k,a,cosm)

    ! Needed for growth function solution
    ! This is the dd in \dot{\delta}=dd
    IMPLICIT NONE
    REAL, INTENT(IN) :: d
    REAL, INTENT(IN) :: v
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: crap
    
    ! To prevent compile-time warnings
    crap=d
    crap=k
    crap=cosm%A
    crap=a    

    ddda=v

  END FUNCTION ddda

  REAL FUNCTION dvda(d,v,k,a,cosm)

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

    ! To prevent compile-time warning
    crap=k

    !f1=1.5*Omega_m_norad(a,cosm)*d/(a**2)
    f1=1.5*Omega_cold_norad(a,cosm)*d/(a**2)
    f2=-(2.+AH_norad(a,cosm)/Hubble2_norad(a,cosm))*(v/a)
    dvda=f1+f2

  END FUNCTION dvda

  REAL FUNCTION dvdanl(d,v,k,a,cosm)

    ! Function used for ODE solver in non-linear growth calculation
    IMPLICIT NONE
    REAL, INTENT(IN) :: d
    REAL, INTENT(IN) :: v
    REAL, INTENT(IN) :: k
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: f1, f2, f3
    REAL :: crap

    ! To prevent compile-time warning
    crap=k

    !f1=1.5*Omega_m_norad(a,cosm)*d*(1.+d)/(a**2)
    f1=1.5*Omega_cold_norad(a,cosm)*d*(1.+d)/(a**2)
    f2=-(2.+AH_norad(a,cosm)/Hubble2_norad(a,cosm))*(v/a)
    f3=4.*(v**2)/(3.*(1.+d))

    dvdanl=f1+f2+f3

  END FUNCTION dvdanl

  REAL FUNCTION dc_NakamuraSuto(a,cosm)

    ! Nakamura & Suto (1997; arXiv:astro-ph/9612074) fitting formula for LCDM
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: Om_mz

    !Om_mz=Omega_m_norad(a,cosm)
    Om_mz=Omega_cold_norad(a,cosm)
    dc_NakamuraSuto=dc0*(1.+0.012299*log10(Om_mz))

  END FUNCTION dc_NakamuraSuto

  REAL FUNCTION Dv_BryanNorman(a,cosm)

    ! Bryan & Norman (1998; arXiv:astro-ph/9710107) spherical over-density fitting function
    ! Here overdensity is defined relative to the background matter density, rather than the critical density
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: x, Om_mz    

    !Om_mz=Omega_m_norad(a,cosm)
    Om_mz=Omega_cold_norad(a,cosm)
    x=Om_mz-1.

    IF(cosm%Om_v==0. .AND. cosm%Om_w==0.) THEN
       ! Open model results
       Dv_BryanNorman=Dv0+60.*x-32.*x**2
       Dv_BryanNorman=Dv_BryanNorman/Om_mz
    ELSE
       ! LCDM results
       Dv_BryanNorman=Dv0+82.*x-39.*x**2
       Dv_BryanNorman=Dv_BryanNorman/Om_mz
    END IF

  END FUNCTION Dv_BryanNorman

  REAL FUNCTION dc_Mead(a,cosm)

    ! delta_c fitting function from Mead (2017)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a ! scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: lg, bG, Om_m

    ! See Appendix A of Mead (2017) for naming convention
    REAL, PARAMETER :: p10=-0.0069
    REAL, PARAMETER :: p11=-0.0208
    REAL, PARAMETER :: p12=0.0312
    REAL, PARAMETER :: p13=0.0021
    INTEGER, PARAMETER :: a1=1
    REAL, PARAMETER :: p20=0.0001
    REAL, PARAMETER :: p21=-0.0647
    REAL, PARAMETER :: p22=-0.0417
    REAL, PARAMETER :: p23=0.0646
    INTEGER, PARAMETER :: a2=0

    lg=ungrow(a,cosm)
    bG=acc_growth(a,cosm)
    !Om_m=Omega_m_norad(a,cosm)
    Om_m=Omega_cold_norad(a,cosm)

    dc_Mead=1.
    dc_Mead=dc_Mead+f_Mead(lg/a,bG/a,p10,p11,p12,p13)*log10(Om_m)**a1
    dc_Mead=dc_Mead+f_Mead(lg/a,bG/a,p20,p21,p22,p23)
    dc_Mead=dc_Mead*dc0

  END FUNCTION dc_Mead

  REAL FUNCTION Dv_Mead(a,cosm)

    ! Delta_v fitting function from Mead (2017)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a !scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: lg, bG, Om_m

    ! See Appendix A of Mead (2017) for naming convention
    REAL, PARAMETER :: p30=-0.79
    REAL, PARAMETER :: p31=-10.17
    REAL, PARAMETER :: p32=2.51
    REAL, PARAMETER :: p33=6.51
    INTEGER, PARAMETER :: a3=1
    REAL, PARAMETER :: p40=-1.89
    REAL, PARAMETER :: p41=0.38
    REAL, PARAMETER :: p42=18.8
    REAL, PARAMETER :: p43=-15.87
    INTEGER, PARAMETER :: a4=2

    lg=ungrow(a,cosm)
    bG=acc_growth(a,cosm)
    !Om_m=Omega_m_norad(a,cosm)
    Om_m=Omega_cold_norad(a,cosm)

    Dv_Mead=1.
    Dv_Mead=Dv_Mead+f_Mead(lg/a,bG/a,p30,p31,p32,p33)*log10(Om_m)**a3
    Dv_Mead=Dv_Mead+f_Mead(lg/a,bG/a,p40,p41,p42,p43)*log10(Om_m)**a4
    Dv_Mead=Dv_Mead*Dv0

  END FUNCTION Dv_Mead

  REAL FUNCTION f_Mead(x,y,p0,p1,p2,p3)

    ! Equation A3 in Mead (2017)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x, y
    REAL, INTENT(IN) :: p0, p1, p2, p3

    f_Mead=p0+p1*(1.-x)+p2*(1.-x)**2+p3*(1.-y)

  END FUNCTION f_Mead

  REAL FUNCTION dc_spherical(a,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)
    
    IF(log(a)<cosm%log_a_dcDv(1)) THEN
       dc_spherical=dc0
    ELSE
       dc_spherical=find(log(a),cosm%log_a_dcDv,cosm%dc,cosm%n_dcDv,3,3,2)
    END IF

  END FUNCTION dc_spherical

  REAL FUNCTION Dv_spherical(a,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)
    
    IF(log(a)<cosm%log_a_dcDv(1)) THEN
       Dv_spherical=Dv0
    ELSE
       Dv_spherical=find(log(a),cosm%log_a_dcDv,cosm%Dv,cosm%n_dcDv,3,3,2)
    END IF

  END FUNCTION Dv_spherical

  SUBROUTINE init_spherical_collapse(cosm)

    USE table_integer
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: dinit, ainit, vinit, ac, dc
    REAL :: av, a_rmax, d_rmax, Dv, rmax, rv
    REAL, ALLOCATABLE :: d(:), a(:), v(:)
    REAL, ALLOCATABLE :: dnl(:), vnl(:), rnl(:)
    REAL, ALLOCATABLE :: a_coll(:), r_coll(:)
    INTEGER :: i, j, k, k2, n, m

    n=n_spherical ! Number of points in integration
    m=m_spherical ! Number of scale factors to try and calculate (usually you get fewer)

    IF(cosm%verbose) WRITE(*,*) 'SPHERICAL_COLLAPSE: Doing integration'

    ! Allocate arrays
    IF(ALLOCATED(cosm%log_a_dcDv)) DEALLOCATE(cosm%log_a_dcDv)
    IF(ALLOCATED(cosm%dc))         DEALLOCATE(cosm%dc)
    IF(ALLOCATED(cosm%Dv))         DEALLOCATE(cosm%Dv)
    ALLOCATE(cosm%log_a_dcDv(m))
    ALLOCATE(cosm%dc(m))
    ALLOCATE(cosm%Dv(m))
    cosm%log_a_dcDv=0.
    cosm%dc=0.
    cosm%Dv=0.

    IF(cosm%verbose) THEN
       WRITE(*,*) 'SPHERICAL_COLLAPSE: delta min', dmin_spherical
       WRITE(*,*) 'SPHERICAL_COLLAPSE: delta max', dmax_spherical
       WRITE(*,*) 'SPHERICAL_COLLAPSE: number of collapse points attempted', m
    END IF

    ! BCs for integration. Note ainit=dinit means that collapse should occur around a=1 for dmin
    ! amax should be slightly greater than 1 to ensure at least a few points for a>0.9 (i.e not to miss out a=1)
    ! TODO: Change to account for massive neutrinos: g(a)=a^(1-3*f_nu/5)
    !ainit=dmin_spherical
    !vinit=1.*(dmin_spherical/ainit) ! vinit=1 is EdS growing mode solution
    ainit=dmin_spherical**(1./(1.-3.*cosm%f_nu/5.))
    vinit=(1.-3.*cosm%f_nu/5.)*ainit**(-3.*cosm%f_nu/5.) ! vinit=1 is EdS growing mode solution

    ! Now loop over all initial density fluctuations
    DO j=1,m
       
       ! log range of initial delta
       dinit=exp(progression(log(dmin_spherical),log(dmax_spherical),j,m))

       ! Do both with the same a1 and a2 and using the same number of time steps
       ! This means that arrays a, and anl will be identical, which simplifies calculation
       CALL ODE_spherical(dnl,vnl,0.,a,cosm,ainit,amax_spherical,dinit,vinit,ddda,dvdanl,n,3,.TRUE.)
       DEALLOCATE(a)
       CALL ODE_spherical(d,v,0.,a,cosm,ainit,amax_spherical,dinit,vinit,ddda,dvda,n,3,.TRUE.)

       ! If this condition is met then collapse occured some time a<amax
       IF(dnl(n)==0.) THEN

          !! delta_c calcualtion !!

          ALLOCATE(rnl(n))

          rnl=a*(1.+dnl)**(-1./3.)

          ! Find the collapse point (very crude)
          ! More accurate calculations seem to be worse
          ! I think this is due to the fact that delta spikes very quickly
          DO i=1,n
             IF(dnl(i)==0.) THEN
                ! k is the new maxium size of the arrays
                k=i-1
                EXIT
             END IF
          END DO

          ! Cut away parts of the arrays for a>ac
          !CALL amputate(a,n,k)
          !CALL amputate(d,n,k)
          !CALL amputate(dnl,n,k)
          !CALL amputate(rnl,n,k)
          CALL amputate_array(a,n,1,k)
          CALL amputate_array(d,n,1,k)
          CALL amputate_array(dnl,n,1,k)
          CALL amputate_array(rnl,n,1,k)

          ! Collapse has occured so use previous a as ac and d as dc
          ac=a(k)
          dc=d(k)

          !! !!

          !! Now to Delta_v calculation !!

          ! Find the a values when the perturbation is maximum size
          a_rmax=maximum(a,rnl,k)

          ! Find the over-density at this point
          d_rmax=exp(find(log(a_rmax),log(a),log(dnl),SIZE(a),1,3,2))

          ! Find the maximum radius
          rmax=find(log(a_rmax),log(a),rnl,SIZE(a),1,3,2)

          ! The radius of the perturbation when it is virialised is half maximum
          ! This might not be appropriate for LCDM models (or anything with DE)
          rv=rmax/2.

          ! Need to assign new arrays for the collapse branch of r such that it is monotonic
          !k2=int_split(d_rmax,dnl,k)
          k2=select_table_integer(d_rmax,dnl,k,imeth=3)

          ! Allocate collapse branch arrays
          ALLOCATE(a_coll(k-k2+1),r_coll(k-k2+1))

          ! Fill collapse branch arrays
          DO i=k2,k
             a_coll(i-k2+1)=a(i)
             r_coll(i-k2+1)=rnl(i)
          END DO

          ! Find the scale factor when the perturbation has reached virial radius
          av=exp(find(rv,r_coll,log(a_coll),SIZE(r_coll),3,3,2))

          ! Deallocate collapse branch arrays
          DEALLOCATE(a_coll,r_coll)

          ! Spherical model approximation is that perturbation is at virial radius when
          ! 'collapse' is considered to have occured, which has already been calculated
          Dv=exp(find(log(av),log(a),log(dnl),SIZE(a),1,3,2))*(ac/av)**3.
          Dv=Dv+1.

          !!

          cosm%log_a_dcDv(j)=ac
          cosm%dc(j)=dc
          cosm%Dv(j)=Dv

          DEALLOCATE(rnl)

       END IF

       ! Deallocate arrays ready for next calculation
       DEALLOCATE(d,v,a)
       DEALLOCATE(dnl,vnl)

    END DO

    IF(cosm%verbose) WRITE(*,*) 'SPHERICAL COLLAPSE: calculation complete'

    ! Reverse the arrays so that they run lowest a to highest a
    CALL reverse_array(cosm%log_a_dcDv,m)
    CALL reverse_array(cosm%dc,m)
    CALL reverse_array(cosm%Dv,m)

    IF(cosm%verbose) THEN
       WRITE(*,*) '===================================='
       WRITE(*,*) 'Point  scalefactor  delta_c  Delta_v'
       WRITE(*,*) '===================================='
       DO i=1,m
          IF(cosm%log_a_dcDv(i)==0.) EXIT
          WRITE(*,fmt='(I5,F13.4,F9.4,F9.1)') i, cosm%log_a_dcDv(i), cosm%dc(i), cosm%Dv(i)
       END DO
       WRITE(*,*) '===================================='
    END IF

    ! Calculate the maximum sizes for these new arrays
    DO i=1,m
       IF(cosm%log_a_dcDv(i)==0.) EXIT
    END DO
    cosm%n_dcDv=i-1

    IF(cosm%verbose) THEN
       WRITE(*,*) 'SPHERICAL_COLLAPSE: number of collapse points:', cosm%n_dcDv
       WRITE(*,*)
    END IF

    ! Remove bits of the array that are unnecessary
    !CALL amputate(cosm%log_a_dcDv,m,cosm%n_dcDv)
    !CALL amputate(cosm%dc,m,cosm%n_dcDv)
    !CALL amputate(cosm%Dv,m,cosm%n_dcDv)
    CALL amputate_array(cosm%log_a_dcDv,m,1,cosm%n_dcDv)
    CALL amputate_array(cosm%dc,m,1,cosm%n_dcDv)
    CALL amputate_array(cosm%Dv,m,1,cosm%n_dcDv)

    ! Take a logarithm
    cosm%log_a_dcDv=log(cosm%log_a_dcDv)

    ! Set the flag
    cosm%has_spherical=.TRUE.

  END SUBROUTINE init_spherical_collapse

  SUBROUTINE ODE_spherical(x,v,kk,t,cosm,ti,tf,xi,vi,fx,fv,n,imeth,ilog)

    ! Solves 2nd order ODE x''(t) from ti to tf and creates arrays of x, v, t values
    ! I have sometimes called this ODE_crass
    ! It has a fixed number of time steps, n
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
       FUNCTION fx(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fx

       ! fv is what v' is equal to
       FUNCTION fv(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fv

    END INTERFACE

    ! Allocate arrays
    ALLOCATE(x8(n),v8(n),t8(n))

    ! Need to be set to zero for this to work in the spherical-collapse case
    x8=0.d0
    v8=0.d0
    t8=0.d0

    ! xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    x8(1)=xi
    v8(1)=vi

    ! Fill time array
    IF(ilog) THEN
       CALL fill_array_double(dble(log(ti)),dble(log(tf)),t8,n)
       t8=exp(t8)
    ELSE
       CALL fill_array_double(dble(ti),dble(tf),t8,n)
    END IF

    DO i=1,n-1
       
       CALL ODE_advance_cosmology(x8(i),x8(i+1),v8(i),v8(i+1),t8(i),t8(i+1),fx,fv,imeth,kk,cosm)
       
       ! Needed to escape from the ODE solver when the perturbation is ~collapsed
       IF(x8(i+1)>dinf_spherical) EXIT
       
    END DO

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(v)) DEALLOCATE(v)
    IF(ALLOCATED(t)) DEALLOCATE(t)
    ALLOCATE(x(n),v(n),t(n))
    x=real(x8)
    v=real(v8)
    t=real(t8)

  END SUBROUTINE ODE_spherical

  SUBROUTINE ODE_adaptive_cosmology(x,v,kk,t,cosm,ti,tf,xi,vi,fx,fv,acc,imeth,ilog)

    ! Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    ! acc is the desired accuracy across the entire solution
    ! time steps are increased until convergence is achieved
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
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: imeth
    LOGICAL, INTENT(IN) :: ilog
    DOUBLE PRECISION, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)      
    INTEGER :: i, j, n, k, np, ifail, kn
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=100

    ! imeth sets ODE solving method
    ! imeth = 1: Crude method
    ! imeth = 2: Mid-point method
    ! imeth = 3: Runge-Kutta   

    INTERFACE

       ! fx is what x' is equal to
       FUNCTION fx(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fx

       ! fv is what v' is equal to
       FUNCTION fv(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fv

    END INTERFACE

    DO j=1,jmax

       n=1+ninit*(2**(j-1))

       ALLOCATE(x8(n),v8(n),t8(n))

       x8=0.d0
       v8=0.d0
       t8=0.d0

       ! xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
       x8(1)=xi
       v8(1)=vi

       ! Fill time array
       IF(ilog) THEN
          CALL fill_array_double(dble(log(ti)),dble(log(tf)),t8,n)
          t8=exp(t8)
       ELSE
          CALL fill_array_double(dble(ti),dble(tf),t8,n)
       END IF

       ifail=0

       DO i=1,n-1
          CALL ODE_advance_cosmology(x8(i),x8(i+1),v8(i),v8(i+1),t8(i),t8(i+1),fx,fv,imeth,kk,cosm)
       END DO

       IF(j==1) ifail=1

       IF(j .NE. 1) THEN

          np=1+(n-1)/2

          DO k=1,1+(n-1)/2

             kn=2*k-1

             IF(ifail==0) THEN

                IF(xh(k)>acc .AND. x8(kn)>acc .AND. (abs(xh(k)/x8(kn))-1.)>acc) ifail=1
                IF(vh(k)>acc .AND. v8(kn)>acc .AND. (abs(vh(k)/v8(kn))-1.)>acc) ifail=1

                IF(ifail==1) THEN
                   DEALLOCATE(xh,th,vh)
                   EXIT
                END IF

             END IF
          END DO

       END IF

       IF(ifail==0) THEN
          IF(ALLOCATED(x)) DEALLOCATE(x)
          IF(ALLOCATED(v)) DEALLOCATE(v)
          IF(ALLOCATED(t)) DEALLOCATE(t)
          ALLOCATE(x(n),v(n),t(n))
          x=real(x8)
          v=real(v8)
          t=real(t8)
          EXIT
       END IF

       ALLOCATE(xh(n),th(n),vh(n))
       xh=x8
       vh=v8
       th=t8
       DEALLOCATE(x8,t8,v8)

    END DO

  END SUBROUTINE ODE_adaptive_cosmology

  SUBROUTINE ODE_advance_cosmology(x1,x2,v1,v2,t1,t2,fx,fv,imeth,k,cosm)

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
       FUNCTION fx(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fx

       ! fv is what v' is equal to
       FUNCTION fv(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fv

    END INTERFACE

    x=real(x1)
    v=real(v1)
    t=real(t1)

    dt=real(t2-t1)

    IF(imeth==1) THEN

       ! Crude method!
       kx1=dt*fx(x,v,k,t,cosm)
       kv1=dt*fv(x,v,k,t,cosm)

       x2=x1+kx1
       v2=v1+kv1

    ELSE IF(imeth==2) THEN

       ! Mid-point method!
       kx1=dt*fx(x,v,k,t,cosm)
       kv1=dt*fv(x,v,k,t,cosm)
       kx2=dt*fx(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)
       kv2=dt*fv(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)

       x2=x1+kx2
       v2=v1+kv2

    ELSE IF(imeth==3) THEN

       ! RK4 (Holy Christ, this is so fast compared to above methods)!
       kx1=dt*fx(x,v,k,t,cosm)
       kv1=dt*fv(x,v,k,t,cosm)
       kx2=dt*fx(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)
       kv2=dt*fv(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)
       kx3=dt*fx(x+kx2/2.,v+kv2/2.,k,t+dt/2.,cosm)
       kv3=dt*fv(x+kx2/2.,v+kv2/2.,k,t+dt/2.,cosm)
       kx4=dt*fx(x+kx3,v+kv3,k,t+dt,cosm)
       kv4=dt*fv(x+kx3,v+kv3,k,t+dt,cosm)

       x2=x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.d0
       v2=v1+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.d0

    ELSE

       STOP 'ODE_ADVANCE: Error, imeth specified incorrectly'

    END IF

  END SUBROUTINE ODE_advance_cosmology

  REAL FUNCTION integrate1_cosm(a,b,f,cosm,acc,iorder)

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
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       FUNCTION f(x,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate1_cosm=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/real(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=f(a,cosm)
             f2=f(b,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*real(i-1)/real(n-1)
                fx=f(x,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (abs(-1.d0+sum_new/sum_old)<acc)) THEN
             ! jmin avoids spurious early convergence
             !integrate=real(sum_new)
             !WRITE(*,*) 'INTEGRATE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate1_cosm=real(sum_new)

    END IF

  END FUNCTION integrate1_cosm

  REAL FUNCTION integrate2_cosm(a,b,f,y,cosm,acc,iorder)

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
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       FUNCTION f(x,y,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x
         REAL, INTENT(IN) :: y
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate2_cosm=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/real(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=f(a,y,cosm)
             f2=f(b,y,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*real(i-1)/real(n-1)
                fx=f(x,y,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (abs(-1.d0+sum_new/sum_old)<acc)) THEN
             ! jmin avoids spurious early convergence
             !integrate=real(sum_new)
             !WRITE(*,*) 'INTEGRATE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate2_cosm=real(sum_new)

    END IF

  END FUNCTION integrate2_cosm

  REAL FUNCTION integrate3_cosm(a,b,f,y,z,cosm,acc,iorder)

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
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       FUNCTION f(x,y,z,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: x
         REAL, INTENT(IN) :: y
         REAL, INTENT(IN) :: z
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate3_cosm=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/real(n-1)

          IF(j==1) THEN

             ! The first go is just the trapezium of the end points
             f1=f(a,y,z,cosm)
             f2=f(b,y,z,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*real(i-1)/real(n-1)
                fx=f(x,y,z,cosm)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (abs(-1.d0+sum_new/sum_old)<acc)) THEN
             ! jmin avoids spurious early convergence
             !integrate=real(sum_new)
             !WRITE(*,*) 'INTEGRATE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate3_cosm=real(sum_new)

    END IF

  END FUNCTION integrate3_cosm

  SUBROUTINE get_CAMB_power(non_linear,halofit_version,cosm)

    USE CAMB_stuff
    USE string_operations
    IMPLICIT NONE
    LOGICAL,INTENT(IN) :: non_linear
    INTEGER, INTENT(IN) :: halofit_version
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, ALLOCATABLE :: k(:), Pk(:)
    INTEGER :: i, nk
    REAL :: Om_c, Om_b, Om_nu, h, ombh2, omch2, omnuh2
    CHARACTER(len=256), PARAMETER :: camb='camb'
    CHARACTER(len=256), PARAMETER :: dir='/Users/Mead/Physics/CAMB_files/tmp/'
    CHARACTER(len=256), PARAMETER :: root=trim(dir)//'temp'
    CHARACTER(len=256), PARAMETER :: matterpower=trim(root)//'_matterpower_'
    CHARACTER(len=256), PARAMETER :: transfer=trim(root)//'_transfer_'
    CHARACTER(len=256), PARAMETER :: params=trim(root)//'_params.ini'

    ! Set the cosmological parameters for CAMB
    IF(cosm%power_Omegas) THEN
       Om_c=cosm%Om_c_pow
       Om_b=cosm%Om_b_pow
       Om_nu=0.
       h=cosm%h_pow
    ELSE
       Om_c=cosm%Om_c
       Om_b=cosm%Om_b
       Om_nu=cosm%Om_nu       
       h=cosm%h
    END IF

    ! For scale-dependent growth fill the scale-factor array first
    IF(cosm%growk) THEN
      CALL fill_array(log(amin_plin),log(amax_plin),cosm%log_a_plin,na_plin)
      cosm%na_plin=na_plin
    END IF  

    ! Physical density parameters that CAMB requires
    ombh2=Om_b*h**2
    omch2=Om_c*h**2
    omnuh2=Om_nu*h**2

    !IF(cosm%Om_v .NE. 0.) STOP 'GET_CAMB_POWER: Error, Omega_v not zero, should set Omega_w'

    ! Remove previous parameters file
    CALL SYSTEM('rm '//trim(params))

    ! Open new parameters file for writing
    OPEN(7,file=params,status='replace')

    ! Output root
    WRITE(7,*) 'output_root = ', trim(root)

    ! Things to get
    WRITE(7,*) 'get_scalar_cls = F'
    WRITE(7,*) 'get_vector_cls = F'
    WRITE(7,*) 'get_tensor_cls = F'
    WRITE(7,*) 'get_transfer = T'

    ! Lensing
    WRITE(7,*) 'do_lensing = F'
    
    IF(non_linear) THEN
       WRITE(7,*) 'do_nonlinear = 1'
    ELSE
       WRITE(7,*) 'do_nonlinear = 0'
    END IF    
    WRITE(7,*) 'halofit_version =', halofit_version
    
    ! Standard cosmological parameters
    WRITE(7,*) 'ombh2 =', ombh2
    WRITE(7,*) 'omch2 =', omch2
    WRITE(7,*) 'omnuh2 =', omnuh2
    WRITE(7,*) 'omk = ', cosm%Om_k
    WRITE(7,*) 'hubble =', 100.*h

    ! Dark energy
    WRITE(7,*) 'dark_energy_model = ppf'
    WRITE(7,*) 'w =', cosm%w
    WRITE(7,*) 'wa =', cosm%wa
    WRITE(7,*) 'cs2_lam = 1'

    ! CMB temperature and Helium fraction
    WRITE(7,*) 'temp_cmb = ', cosm%T_CMB
    WRITE(7,*) 'helium_fraction = ', cosm%YHe

    ! Neutrinos
    !WRITE(7,*) 'massless_neutrinos = ', cosm%neff-real(cosm%N_massive_nu)   
    WRITE(7,*) 'share_delta_neff = T'   
    IF(omnuh2==0.) THEN
       ! Massless neutrinos
       WRITE(7,*) 'massless_neutrinos = ', cosm%neff
       WRITE(7,*) 'nu_mass_eigenstates = 0'
       WRITE(7,*) 'massive_neutrinos = 0'
    ELSE
       ! Three equal-mass neutrinos
       WRITE(7,*) 'massless_neutrinos = ', cosm%neff-3.
       WRITE(7,*) 'nu_mass_eigenstates = 1'
       WRITE(7,*) 'massive_neutrinos = 3'
       WRITE(7,*) 'nu_mass_fractions = 1'
    END IF

    ! Primoridial power spectrum properties
    WRITE(7,*) 'initial_power_num = 1'
    WRITE(7,*) 'scalar_spectral_index(1) =', cosm%n
    WRITE(7,*) 'scalar_nrun(1) = 0'
    WRITE(7,*) 'scalar_amp(1) =', cosm%A
    WRITE(7,*) 'pivot_scalar = 0.05'
    WRITE(7,*) 'pivot_tensor = 0.05'

    ! Reionisation
    WRITE(7,*) 'reionization = T'
    WRITE(7,*) 're_use_optical_depth = T'
    WRITE(7,*) 're_optical_depth = 0.09'
    WRITE(7,*) 're_delta_redshift = 1.5'
    WRITE(7,*) 're_ionization_frac = -1'
    
    ! RECFAST
    WRITE(7,*) 'recombination_model = Recfast'
    WRITE(7,*) 'RECFAST_fudge = 1.14'
    WRITE(7,*) 'RECFAST_fudge_He = 0.86'
    WRITE(7,*) 'RECFAST_Heswitch = 6'
    WRITE(7,*) 'RECFAST_Hswitch = T'

    ! Adiabatic/isocurvature etc.
    WRITE(7,*) 'initial_condition = 1' ! 1 - Adiabatic initial condtions

    ! ?
    WRITE(7,*) 'COBE_normalize = F'
    WRITE(7,*) 'CMB_outputscale = 7.4311e12'

    ! Power spectrum requirements
    WRITE(7,*) 'transfer_high_precision = F'
    WRITE(7,*) 'transfer_kmax = 10'
    WRITE(7,*) 'transfer_k_per_logint = 0'
    WRITE(7,*) 'transfer_interp_matterpower = T'
    WRITE(7,*) 'transfer_power_var = 7'
    IF(cosm%growk) THEN
      WRITE(7,*) 'transfer_num_redshifts =', na_plin
      DO i=1,na_plin
         !WRITE(7,*) 'transfer_redshift(1) = 0'
         WRITE(7,*) number_file(trim('transfer_redshift('),i,trim(') =')), redshift_a(exp(cosm%log_a_plin(i)))
      END DO 
    ELSE
      WRITE(7,*) 'transfer_num_redshifts = 1'
      WRITE(7,*) 'transfer_redshift(1) = 0'
    END IF

    ! Bispectra
    WRITE(7,*) 'do_lensing_bispectrum = F'
    WRITE(7,*) 'do_primordial_bispectrum = F'

    ! Verbosity
    WRITE(7,*) 'feedback_level = 1'

    ! Print out a comment describing file headers
    WRITE(7,*) 'output_file_headers = T'

    ! Write out derived parameters
    WRITE(7,*) 'derived_parameters = T'

    ! ?
    WRITE(7,*) 'massive_nu_approx = 1'

    ! Accuracy parameters for polarization and reionization
    WRITE(7,*) 'accurate_polarization = T'
    WRITE(7,*) 'accurate_reionization = T'

    WRITE(7,*) 'lensing_method = 1'
    WRITE(7,*) 'accurate_BB = F'

    ! ?
    WRITE(7,*) 'do_tensor_neutrinos = T'

    ! ?
    WRITE(7,*) 'do_late_rad_truncation = T'

    ! OMP threads
    WRITE(7,*) 'number_of_threads = 4'

    ! Accuracy parameters
    WRITE(7,*) 'accuracy_boost = 1'
    WRITE(7,*) 'l_accuracy_boost = 1'
    !WRITE(7,*) 'high_accuracy_default = T' ! Removed in CAMB v1
    !WRITE(7,*) 'use_spline_template =  T'  ! Removed in CAMB v1
    WRITE(7,*) 'l_sample_boost = 1'

    ! Close file
    CLOSE(7)

    ! Remove old files and run CAMB
    IF(cosm%verbose) WRITE(*,*) 'GET_CAMB_POWER: Running CAMB (note weird problems with this function in library)'
    CALL SYSTEM('rm '//trim(transfer)//'*')
    CALL SYSTEM('rm '//trim(matterpower)//'*')
    IF(cosm%verbose) THEN
       CALL SYSTEM(trim(camb)//' '//trim(params)//' > /dev/null')
    ELSE
       CALL SYSTEM(trim(camb)//' '//trim(params))
    END IF
    IF(cosm%verbose) WRITE(*,*) 'GET_CAMB_POWER: CAMB run complete'

   IF(cosm%growk) THEN

      ! Loop over redshifts and read CAMB power into arrays
      DO i=1,na_plin
         CALL read_CAMB_Pk(k,Pk,nk,number_file(matterpower,i,trim('.dat')))
         IF(i==1) THEN
            cosm%nk_plin=nk
            IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
            IF(ALLOCATED(cosm%log_plina))  DEALLOCATE(cosm%log_plina)
            ALLOCATE(cosm%log_k_plin(nk))
            ALLOCATE(cosm%log_plina(nk,na_plin))
            cosm%log_k_plin=log(k)
         END IF
         cosm%log_plina(:,i)=log(Pk)
      END DO

   ELSE 

      ! Read in the power spectrum that has been generated
      CALL read_CAMB_Pk(k,Pk,nk,trim(matterpower)//trim('1.dat'))

      ! Add to cosm arrays and convert to log
      IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
      IF(ALLOCATED(cosm%log_plin))   DEALLOCATE(cosm%log_plin)
      ALLOCATE(cosm%log_plin(nk))
      ALLOCATE(cosm%log_k_plin(nk))
      cosm%log_k_plin=log(k)
      cosm%log_plin=log(Pk)
      cosm%nk_plin=nk
   
   END IF

    cosm%has_power=.TRUE.

    IF(cosm%verbose) THEN
       WRITE(*,*) 'GET_CAMB_POWER: Done'
       WRITE(*,*)
    END IF
    
  END SUBROUTINE get_CAMB_power

  SUBROUTINE random_cosmology(cosm)

    ! Generate some random cosmological parameters
    USE random_numbers
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: Om_c_min=0.001
    REAL, PARAMETER :: Om_c_max=0.5

    REAL, PARAMETER :: Om_b_min=0.01
    REAL, PARAMETER :: Om_b_max=0.07

    REAL, PARAMETER :: h_min=0.4
    REAL, PARAMETER :: h_max=1.0

    REAL, PARAMETER :: n_min=0.7
    REAL, PARAMETER :: n_max=1.3

    REAL, PARAMETER :: w_min=-1.3
    REAL, PARAMETER :: w_max=-0.7

    !REAL, PARAMETER :: wa_min=-1.3
    !REAL, PARAMETER :: wa_max=-0.7

    REAL, PARAMETER :: sig8_min=0.6
    REAL, PARAMETER :: sig8_max=0.9

    !CALL RNG_set(seed=0)

    cosm%h=random_uniform(h_min,h_max)

    cosm%Om_c=random_uniform(Om_c_min,Om_c_max)

    cosm%Om_b=random_uniform(Om_b_min,Om_b_max)

    cosm%Om_m=cosm%Om_c+cosm%Om_b

    ! Enforce flatness
    ! Note - need to have Om_w for dark enegry
    cosm%Om_w=1.-cosm%Om_m

    cosm%n=random_uniform(n_min,n_max)   

    cosm%w=random_uniform(w_min,w_max)

    cosm%sig8=random_uniform(sig8_min,sig8_max)

  END SUBROUTINE random_cosmology

  SUBROUTINE random_Cosmic_Emu_cosmology(cosm)

    ! Generate some random cosmological parameters
    USE random_numbers
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: om_m_min=0.120
    REAL, PARAMETER :: om_m_max=0.155

    REAL, PARAMETER :: om_b_min=0.0215
    REAL, PARAMETER :: om_b_max=0.0235

    REAL, PARAMETER :: n_min=0.85
    REAL, PARAMETER :: n_max=1.05

    REAL, PARAMETER :: w_min=-1.3
    REAL, PARAMETER :: w_max=-0.7

    REAL, PARAMETER :: sig8_min=0.616
    REAL, PARAMETER :: sig8_max=0.9

    STOP 'RANDOM_COSMIC_EMU_COSMOLOGY: Need to implement the CMB distance condition on h'

    cosm%Om_m=random_uniform(om_m_min,om_m_max)/cosm%h**2

    cosm%Om_b=random_uniform(om_b_min,om_b_max)/cosm%h**2

    ! Enforce flatness
    ! Note - need to have Om_w for dark enegry
    cosm%Om_w=1.-cosm%Om_m

    cosm%n=random_uniform(n_min,n_max)   

    cosm%w=random_uniform(w_min,w_max)

    cosm%sig8=random_uniform(sig8_min,sig8_max)

  END SUBROUTINE random_Cosmic_Emu_cosmology

  SUBROUTINE random_Franken_Emu_cosmology(cosm)

    ! Generate some random cosmological parameters
    USE random_numbers
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: om_m_min=0.120
    REAL, PARAMETER :: om_m_max=0.155

    REAL, PARAMETER :: om_b_min=0.0215
    REAL, PARAMETER :: om_b_max=0.0235

    REAL, PARAMETER :: n_min=0.85
    REAL, PARAMETER :: n_max=1.05

    REAL, PARAMETER :: h_min=0.55
    REAL, PARAMETER :: h_max=0.85

    REAL, PARAMETER :: w_min=-1.3
    REAL, PARAMETER :: w_max=-0.7

    REAL, PARAMETER :: sig8_min=0.616
    REAL, PARAMETER :: sig8_max=0.9

    cosm%h=random_uniform(h_min,h_max)

    cosm%Om_m=random_uniform(om_m_min,om_m_max)/cosm%h**2

    cosm%Om_b=random_uniform(om_b_min,om_b_max)/cosm%h**2

    ! Enforce flatness
    ! Note - need to have Om_w for dark enegry
    cosm%Om_w=1.-cosm%Om_m

    cosm%n=random_uniform(n_min,n_max)   

    cosm%w=random_uniform(w_min,w_max)

    cosm%sig8=random_uniform(sig8_min,sig8_max)


  END SUBROUTINE random_Franken_Emu_cosmology

  SUBROUTINE Cosmic_Emu_node_cosmology(node,cosm)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: node
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: om_m, om_b

    ! M000 -> M037 of Cosmic_Emu (Table 1 in 1304.7849)
    IF(node==0) THEN
       ! M000
       om_m=0.1296
       om_b=0.0224
       cosm%n=0.9700
       cosm%w=1.000
       cosm%sig8=0.8000
       cosm%h=0.7200
    ELSE IF(node==1) THEN
       ! M001
       om_m=0.1539
       om_b=0.0231
       cosm%n=0.9468
       cosm%w=0.816
       cosm%sig8=0.8161
       cosm%h=0.5977
    ELSE IF(node==2) THEN
       ! M002
       om_m=0.1460
       om_b=0.0227
       cosm%n=0.8952
       cosm%w=0.758
       cosm%sig8=0.8548
       cosm%h=0.5970
    ELSE IF(node==3) THEN
       ! M003
       om_m=0.1324
       om_b=0.0235
       cosm%n=0.9984
       cosm%w=0.874
       cosm%sig8=0.8484
       cosm%h=0.6763
    ELSE IF(node==4) THEN
       ! M004
       om_m=0.1381
       om_b=0.0227
       cosm%n=0.9339
       cosm%w=1.087
       cosm%sig8=0.7000
       cosm%h=0.7204
    ELSE IF(node==5) THEN
       ! M005
       om_m=0.1358
       om_b=0.0216
       cosm%n=0.9726
       cosm%w=1.242
       cosm%sig8=0.8226
       cosm%h=0.7669
    ELSE IF(node==6) THEN
       ! M006
       om_m=0.1516
       om_b=0.0229
       cosm%n=0.9145
       cosm%w=1.223
       cosm%sig8=0.6705
       cosm%h=0.7040
    ELSE IF(node==7) THEN
       ! M007
       om_m=0.1268
       om_b=0.0223
       cosm%n=0.9210
       cosm%w=0.70001 ! Changed to avoid problems
       cosm%sig8=0.7474
       cosm%h=0.6189
    ELSE IF(node==8) THEN
       ! M008
       om_m=0.1448
       om_b=0.0223
       cosm%n=0.9855
       cosm%w=1.203
       cosm%sig8=0.8090
       cosm%h=0.7218
    ELSE IF(node==9) THEN
       ! M009
       om_m=0.1392
       om_b=0.0234
       cosm%n=0.9790
       cosm%w=0.739
       cosm%sig8=0.6692
       cosm%h=0.6127
    ELSE IF(node==10) THEN
       ! M010
       om_m=0.1403
       om_b=0.0218
       cosm%n=0.8565
       cosm%w=0.990
       cosm%sig8=0.7556
       cosm%h=0.6695
    ELSE IF(node==11) THEN
       ! M011
       om_m=0.1437
       om_b=0.0234
       cosm%n=0.8823
       cosm%w=1.126
       cosm%sig8=0.7276
       cosm%h=0.7177
    ELSE IF(node==12) THEN
       ! M012
       om_m=0.1223
       om_b=0.0225
       cosm%n=1.0048
       cosm%w=0.971
       cosm%sig8=0.6271
       cosm%h=0.7396
    ELSE IF(node==13) THEN
       ! M013
       om_m=0.1482
       om_b=0.0221
       cosm%n=0.9597
       cosm%w=0.855
       cosm%sig8=0.6508
       cosm%h=0.6107
    ELSE IF(node==14) THEN
       ! M014
       om_m=0.1471
       om_b=0.0233
       cosm%n=1.0306
       cosm%w=1.010
       cosm%sig8=0.7075
       cosm%h=0.6688
    ELSE IF(node==15) THEN
       ! M015
       om_m=0.1415
       om_b=0.0230
       cosm%n=1.0177
       cosm%w=1.281
       cosm%sig8=0.7692
       cosm%h=0.7737
    ELSE IF(node==16) THEN
       ! M016
       om_m=0.1245
       om_b=0.0218
       cosm%n=0.9403
       cosm%w=1.145
       cosm%sig8=0.7437
       cosm%h=0.7929
    ELSE IF(node==17) THEN
       ! M017
       om_m=0.1426
       om_b=0.0215
       cosm%n=0.9274
       cosm%w=0.893
       cosm%sig8=0.6865
       cosm%h=0.6305
    ELSE IF(node==18) THEN
       ! M018
       om_m=0.1313
       om_b=0.0216
       cosm%n=0.8887
       cosm%w=1.029
       cosm%sig8=0.6440
       cosm%h=0.7136
    ELSE IF(node==19) THEN
       ! M019
       om_m=0.1279
       om_b=0.0232
       cosm%n=0.8629
       cosm%w=1.184
       cosm%sig8=0.6159
       cosm%h=0.8120
    ELSE IF(node==20) THEN
       ! M020
       om_m=0.1290
       om_b=0.0220
       cosm%n=1.0242
       cosm%w=0.797
       cosm%sig8=0.7972
       cosm%h=0.6442
    ELSE IF(node==21) THEN
       ! M021
       om_m=0.1335
       om_b=0.0221
       cosm%n=1.0371
       cosm%w=1.165
       cosm%sig8=0.6563
       cosm%h=0.7601
    ELSE IF(node==22) THEN
       ! M022
       om_m=0.1505
       om_b=0.0225
       cosm%n=1.0500
       cosm%w=1.107
       cosm%sig8=0.7678
       cosm%h=0.6736
    ELSE IF(node==23) THEN
       ! M023
       om_m=0.1211
       om_b=0.0220
       cosm%n=0.9016
       cosm%w=1.261
       cosm%sig8=0.6664
       cosm%h=0.8694
    ELSE IF(node==24) THEN
       ! M024
       om_m=0.1302
       om_b=0.0226
       cosm%n=0.9532
       cosm%w=1.300
       cosm%sig8=0.6644
       cosm%h=0.8380
    ELSE IF(node==25) THEN
       ! M025
       om_m=0.1494
       om_b=0.0217
       cosm%n=1.0113
       cosm%w=0.719
       cosm%sig8=0.7398
       cosm%h=0.5724
    ELSE IF(node==26) THEN
       ! M026
       om_m=0.1347
       om_b=0.0232
       cosm%n=0.9081
       cosm%w=0.952
       cosm%sig8=0.7995
       cosm%h=0.6931
    ELSE IF(node==27) THEN
       ! M027
       om_m=0.1369
       om_b=0.0224
       cosm%n=0.8500
       cosm%w=0.836
       cosm%sig8=0.7111
       cosm%h=0.6387
    ELSE IF(node==28) THEN
       ! M028
       om_m=0.1527
       om_b=0.0222
       cosm%n=0.8694
       cosm%w=0.932
       cosm%sig8=0.8068
       cosm%h=0.6189
    ELSE IF(node==29) THEN
       ! M029
       om_m=0.1256
       om_b=0.0228
       cosm%n=1.0435
       cosm%w=0.913
       cosm%sig8=0.7087
       cosm%h=0.7067
    ELSE IF(node==30) THEN
       ! M030
       om_m=0.1234
       om_b=0.0230
       cosm%n=0.8758
       cosm%w=0.777
       cosm%sig8=0.6739
       cosm%h=0.6626
    ELSE IF(node==31) THEN
       ! M031
       om_m=0.1550
       om_b=0.0219
       cosm%n=0.9919
       cosm%w=1.068
       cosm%sig8=0.7041
       cosm%h=0.6394
    ELSE IF(node==32) THEN
       ! M032
       om_m=0.1200
       om_b=0.0229
       cosm%n=0.9661
       cosm%w=1.048
       cosm%sig8=0.7556
       cosm%h=0.7901
    ELSE IF(node==33) THEN
       ! M033
       om_m=0.1399
       om_b=0.0225
       cosm%n=1.0407
       cosm%w=1.147
       cosm%sig8=0.8645
       cosm%h=0.7286
    ELSE IF(node==34) THEN
       ! M034
       om_m=0.1497
       om_b=0.0227
       cosm%n=0.9239
       cosm%w=1.000
       cosm%sig8=0.8734
       cosm%h=0.6510
    ELSE IF(node==35) THEN
       ! M035
       om_m=0.1485
       om_b=0.0221
       cosm%n=0.9604
       cosm%w=0.853
       cosm%sig8=0.8822
       cosm%h=0.6100
    ELSE IF(node==36) THEN
       ! M036
       om_m=0.1216
       om_b=0.0233
       cosm%n=0.9387
       cosm%w=0.706
       cosm%sig8=0.8911
       cosm%h=0.6421
    ELSE IF(node==37) THEN
       ! M037
       om_m=0.1495
       om_b=0.0228
       cosm%n=1.0233
       cosm%w=1.294
       cosm%sig8=0.8999 ! Moved off boundary
       cosm%h=0.7313
    ELSE
       STOP 'COSMIC_EMU_NODE_COSMOLOGY: Error, node specified incorrectly'
    END IF

    cosm%w=-cosm%w
    cosm%Om_m=om_m/cosm%h**2
    cosm%Om_b=om_b/cosm%h**2
    cosm%Om_w=1.-cosm%Om_m
    
  END SUBROUTINE Cosmic_Emu_node_cosmology

  SUBROUTINE Franken_Emu_node_cosmology(node,cosm)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: node
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: om_m, om_b

    IF(node==23) THEN
       ! M023
       om_m=0.1211
       om_b=0.0220
       cosm%n=0.9016
       cosm%w=-1.261
       cosm%sig8=0.6664
       cosm%h=0.8500 ! Have to round down from 0.8694 to 0.85 for FrankenEmu shrunken space
       cosm%Om_m=om_m/cosm%h**2
       cosm%Om_b=om_b/cosm%h**2
       cosm%Om_w=1.-cosm%Om_m
    ELSE
       CALL Cosmic_Emu_node_cosmology(node,cosm)
    END IF
    
  END SUBROUTINE Franken_Emu_node_cosmology

  SUBROUTINE random_Mira_Titan_cosmology(cosm)

    !Generate some random cosmological parameters for the Mira Titan hypercube
    USE random_numbers
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: om_m, om_b, om_nu

    REAL, PARAMETER :: om_m_min=0.120
    REAL, PARAMETER :: om_m_max=0.155

    REAL, PARAMETER :: om_b_min=0.0215
    REAL, PARAMETER :: om_b_max=0.0235

    REAL, PARAMETER :: om_nu_min=0.00
    REAL, PARAMETER :: om_nu_max=0.01

    REAL, PARAMETER :: n_min=0.85
    REAL, PARAMETER :: n_max=1.05

    REAL, PARAMETER :: h_min=0.55
    REAL, PARAMETER :: h_max=0.85

    REAL, PARAMETER :: w_min=-1.3
    REAL, PARAMETER :: w_max=-0.7

    REAL, PARAMETER :: wa_min=-1.73
    REAL, PARAMETER :: wa_max=1.28

    REAL, PARAMETER :: sig8_min=0.7
    REAL, PARAMETER :: sig8_max=0.9

    cosm%h=random_uniform(h_min,h_max)

    om_m=random_uniform(om_m_min,om_m_max)
    cosm%Om_m=om_m/cosm%h**2

    om_b=random_uniform(om_b_min,om_b_max)
    cosm%Om_b=om_b/cosm%h**2

    om_nu=random_uniform(om_nu_min,om_nu_max)
    !cosm%m_nu=neutrino_constant*om_nu/3. ! Split equally over 3 neutrinos
    cosm%m_nu=neutrino_constant*om_nu
    cosm%m_nu=0.

    ! Enforce flatness, ensure Omega_w is used for dark energy, Omega_v = 0
    cosm%Om_w=1.-cosm%Om_m

    cosm%n=random_uniform(n_min,n_max)

    cosm%w=random_uniform(w_min,w_max)

    ! Enforce 0.3 <= (-w0-wa)^(1/4)
    DO
       cosm%wa=random_uniform(wa_min,wa_max)
       !IF(0.3<=(-cosm%w-cosm%wa)**(1./4.)) EXIT
       IF(0.0081<=-cosm%w-cosm%wa .AND. 2.769>=-cosm%w-cosm%wa) EXIT
    END DO

    cosm%sig8=random_uniform(sig8_min,sig8_max)

  END SUBROUTINE random_Mira_Titan_cosmology

  SUBROUTINE Mira_Titan_node_cosmology(node,cosm)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: node
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: om_m, om_b, om_nu

    ! M000 -> M036 of Mira Titan (m_nu = 0 eV for M000 -> M010; Table 3 in 1705.03388)
    IF(node==0) THEN
       ! M000
       om_m=0.1335
       om_b=0.02258
       cosm%sig8=0.8
       cosm%h=0.71
       cosm%n=0.963
       cosm%w=-1.0
       cosm%wa=0.0
       om_nu=0.
    ELSE IF(node==1) THEN
      ! M001
      ! Something seems funnny with this cosmology, the halo-model prediction is awful
      ! This model has w0 ~ -wa ~ 0.7, so w(a=0) = w0+wa ~ 0, so DE scales almost like matter at early times
      om_m=0.1472
      om_b=0.02261
      cosm%sig8=0.8778
      cosm%h=0.6167
      cosm%n=0.9611
      cosm%w=-0.7000
      cosm%wa=0.67220
      om_nu=0.
    ELSE IF(node==2) THEN
       ! M002
       om_m=0.1356
       om_b=0.02328
       cosm%sig8=0.8556
       cosm%h=0.7500
       cosm%n=1.0500
       cosm%w=-1.0330
       cosm%wa=0.91110
       om_nu=0.
    ELSE IF(node==3) THEN
       ! M003
       om_m=0.1550
       om_b=0.02194
       cosm%sig8=0.9000
       cosm%h=0.7167
        cosm%n=0.8944
       cosm%w=-1.1000
       cosm%wa=-0.28330
       om_nu=0.
    ELSE IF(node==4) THEN
       ! M004
       om_m=0.1239
       om_b=0.02283
       cosm%sig8=0.7889
       cosm%h=0.5833
       cosm%n=0.8722
       cosm%w=-1.1670
       cosm%wa=1.15000
       om_nu=0.
    ELSE IF(node==5) THEN
       ! M005
       om_m=0.1433
       om_b=0.02350
       cosm%sig8=0.7667
       cosm%h=0.8500
       cosm%n=0.9833
       cosm%w=-1.2330
       cosm%wa=-0.04445
       om_nu=0.
    ELSE IF(node==6) THEN
       ! M006
       om_m=0.1317
       om_b=0.021501 ! Moved off boundary
       cosm%sig8=0.8333
       cosm%h=0.5500
       cosm%n=0.9167
       cosm%w=-0.7667
       cosm%wa=0.19440
       om_nu=0.
    ELSE IF(node==7) THEN
       ! M007
       om_m=0.1511
       om_b=0.02217
       cosm%sig8=0.8111
       cosm%h=0.8167
       cosm%n=1.0280
       cosm%w=-0.8333
       cosm%wa=-1.0000
       om_nu=0.
    ELSE IF(node==8) THEN
       ! M008
       om_m=0.1200
       om_b=0.02306
       cosm%sig8=0.7000
       cosm%h=0.6833
       cosm%n=1.0060
       cosm%w=-0.9000
       cosm%wa=0.43330
    ELSE IF(node==9) THEN
       ! M009
       om_m=0.1394
       om_b=0.02172
       cosm%sig8=0.7444
       cosm%h=0.6500
       cosm%n=0.8500
       cosm%w=-0.9667
       cosm%wa=-0.76110
       om_nu=0.
    ELSE IF(node==10) THEN
       ! M010
       om_m=0.1278
       om_b=0.02239
       cosm%sig8=0.7222
       cosm%h=0.7833
       cosm%n=0.9389
       cosm%w=-1.3000
       cosm%wa=-0.52220
       om_nu=0.
    ELSE IF(node==11) THEN
       ! M011
       om_m=0.1227
       om_b=0.0220
       cosm%sig8=0.7151
       cosm%h=0.5827
       cosm%n=0.9357
       cosm%w=-1.0821
       cosm%wa=1.0646
       om_nu=0.000345
    ELSE IF(node==12) THEN
       ! M012
       om_m=0.1241
       om_b=0.0224
       cosm%sig8=0.7472
       cosm%h=0.8315
       cosm%n=0.8865
       cosm%w=-1.2325
       cosm%wa=-0.7646
       om_nu=0.001204
    ELSE IF(node==13) THEN
       ! M013
       om_m=0.1534
       om_b=0.0232
       cosm%sig8=0.8098
       cosm%h=0.7398
       cosm%n=0.8706
       cosm%w=-1.2993
       cosm%wa=1.2236
       om_nu=0.003770
    ELSE IF(node==14) THEN
       ! M014
       om_m=0.1215
       om_b=0.0215
       cosm%sig8=0.8742
       cosm%h=0.5894
       cosm%n=1.0151
       cosm%w=-0.7281
       cosm%wa=-0.2088
       om_nu=0.001752
    ELSE IF(node==15) THEN
       ! M015
       om_m=0.1250
       om_b=0.0224
       cosm%sig8=0.8881
       cosm%h=0.6840
       cosm%n=0.8638
       cosm%w=-1.0134
       cosm%wa=0.0415
       om_nu=0.002789
    ELSE IF(node==16) THEN
       ! M016
       om_m=0.1499
       om_b=0.0223
       cosm%sig8=0.7959
       cosm%h=0.6452
       cosm%n=1.0219
       cosm%w=-1.0139
       cosm%wa=0.9434
       om_nu=0.002734
    ELSE IF(node==17) THEN
       ! M017
       om_m=0.1206
       om_b=0.0215
       cosm%sig8=0.7332
       cosm%h=0.7370
       cosm%n=1.0377
       cosm%w=-0.9472
       cosm%wa=-0.9897
       om_nu=0.000168
    ELSE IF(node==18) THEN
       ! M018
       om_m=0.1544
       om_b=0.0217
       cosm%sig8=0.7982
       cosm%h=0.6489
       cosm%n=0.9026
       cosm%w=-0.7091
       cosm%wa=0.6409
       om_nu=0.006419
    ELSE IF(node==19) THEN
       ! M019
       om_m=0.1256
       om_b=0.0222
       cosm%sig8=0.8547
       cosm%h=0.8251
       cosm%n=1.0265
       cosm%w=-0.9813
       cosm%wa=-0.3393
       om_nu=0.004673
    ELSE IF(node==20) THEN
       ! M020
       om_m=0.1514
       om_b=0.0225
       cosm%sig8=0.7561
       cosm%h=0.6827
       cosm%n=0.9913
       cosm%w=-1.0101
       cosm%wa=-0.7778
       om_nu=0.009777
    ELSE IF(node==21) THEN
       ! M021
       om_m=0.1472
       om_b=0.0221
       cosm%sig8=0.8475
       cosm%h=0.6583
       cosm%n=0.9613
       cosm%w=-0.9111
       cosm%wa=-1.5470
       om_nu=0.000672
    ELSE IF(node==22) THEN
       ! M022
       om_m=0.1384
       om_b=0.0231
       cosm%sig8=0.8328
       cosm%h=0.8234
       cosm%n=0.9739
       cosm%w=-0.9312
       cosm%wa=0.5939
       om_nu=0.008239
    ELSE IF(node==23) THEN
       ! M023
       om_m=0.1334
       om_b=0.0225
       cosm%sig8=0.7113
       cosm%h=0.7352
       cosm%n=0.9851
       cosm%w=-0.8971
       cosm%wa=0.3247
       om_nu=0.003733
    ELSE IF(node==24) THEN
       ! M024
       om_m=0.1508
       om_b=0.0229
       cosm%sig8=0.7002
       cosm%h=0.7935
       cosm%n=0.8685
       cosm%w=-1.0322
       cosm%wa=1.0220
       om_nu=0.003063
    ELSE IF(node==25) THEN
       ! M025
       om_m=0.1203
       om_b=0.0230
       cosm%sig8=0.8773
       cosm%h=0.6240
       cosm%n=0.9279
       cosm%w=-0.8282
       cosm%wa=-1.5005
       om_nu=0.007024
    ELSE IF(node==26) THEN
       ! M026
       om_m=0.1224
       om_b=0.0222
       cosm%sig8=0.7785
       cosm%h=0.7377
       cosm%n=0.8618
       cosm%w=-0.7463
       cosm%wa=0.3647
       om_nu=0.002082
    ELSE IF(node==27) THEN
       ! M027
       om_m=0.1229
       om_b=0.0234
       cosm%sig8=0.8976
       cosm%h=0.8222
       cosm%n=0.9698
       cosm%w=-1.0853
       cosm%wa=0.8683
       om_nu=0.002902
    ELSE IF(node==28) THEN
       ! M028
       om_m=0.1229
       om_b=0.0231
       cosm%sig8=0.8257
       cosm%h=0.6109
       cosm%n=0.9885
       cosm%w=-0.9311
       cosm%wa=0.8693
       om_nu=0.009086
    ELSE IF(node==29) THEN
       ! M029
       om_m=0.1274
       om_b=0.0228
       cosm%sig8=0.8999
       cosm%h=0.8259
       cosm%n=0.8505
       cosm%w=-0.7805
       cosm%wa=0.5688
       om_nu=0.006588
    ELSE IF(node==30) THEN
       ! M030
       om_m=0.1404
       om_b=0.0222
       cosm%sig8=0.8232
       cosm%h=0.6852
       cosm%n=0.8679
       cosm%w=-0.8594
       cosm%wa=-0.4637
       om_nu=0.008126
    ELSE IF(node==31) THEN
       ! M031
       om_m=0.1386
       om_b=0.0229
       cosm%sig8=0.7693
       cosm%h=0.6684
       cosm%n=1.0478
       cosm%w=-1.2670
       cosm%wa=1.2536
       om_nu=0.006502
    ELSE IF(node==32) THEN
       ! M032
       om_m=0.1369
       om_b=0.021501 ! Moved off boundary
       cosm%sig8=0.8812
       cosm%h=0.8019
       cosm%n=1.0005
       cosm%w=-0.7282
       cosm%wa=-1.6927
       om_nu=0.000905
    ELSE IF(node==33) THEN
       ! M033
       om_m=0.1286
       om_b=0.0230
       cosm%sig8=0.7005
       cosm%h=0.6752
       cosm%n=1.0492
       cosm%w=-0.7119
       cosm%wa=-0.8184
       om_nu=0.007968
    ELSE IF(node==34) THEN
       ! M034
       om_m=0.1354
       om_b=0.0216
       cosm%sig8=0.7018
       cosm%h=0.5970
       cosm%n=0.8791
       cosm%w=-0.8252
       cosm%wa=-1.1148
       om_nu=0.003602
    ELSE IF(node==35) THEN
       ! M035
       om_m=0.1359
       om_b=0.0228
       cosm%sig8=0.8210
       cosm%h=0.6815
       cosm%n=0.9872
       cosm%w=-1.1642
       cosm%wa=-0.1801
       om_nu=0.004440
    ELSE IF(node==36) THEN
       ! M036
       om_m=0.1390
       om_b=0.0220
       cosm%sig8=0.8631
       cosm%h=0.6477
       cosm%n=0.8985
       cosm%w=-0.8632
       cosm%wa=0.8285
       om_nu=0.001082
    ELSE
       STOP 'MIRA_TITAN_NODE_COSMOLOGY: Error, i specified incorrectly'
    END IF

    cosm%Om_m=om_m/cosm%h**2
    cosm%Om_b=om_b/cosm%h**2
    !cosm%m_nu=neutrino_constant*om_nu/3. ! Split equally over three neutrinos
    cosm%m_nu=neutrino_constant*om_nu
    cosm%Om_w=1.-cosm%Om_m
        
  END SUBROUTINE Mira_Titan_node_cosmology

  SUBROUTINE halofit_init(rknl,rneff,rncur,a,cosm,verbose)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: rknl
    REAL, INTENT(OUT) :: rneff
    REAL, INTENT(OUT) :: rncur
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    REAL :: xlogr1, xlogr2, rmid, sig, d1, d2, diff  

    ! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and 
    ! curvature (rncur) of the power spectrum at the desired redshift, using method 
    ! described in Smith et al (2003).

    IF(verbose) WRITE(*,*) 'HALOFIT_INIT: computing effective spectral quantities:'

    xlogr1=-3.0
    xlogr2=3.5
10  rmid=(xlogr2+xlogr1)/2.0
    rmid=10**rmid
    call wint(rmid,sig,d1,d2,a,cosm)
    diff=sig-1.0
    if (abs(diff).le.0.001) then
       rknl=1./rmid
       rneff=-3-d1
       rncur=-d2                  
    elseif (diff.gt.0.001) then
       xlogr1=log10(rmid)
       goto 10
    elseif (diff.lt.-0.001) then
       xlogr2=log10(rmid)
       goto 10
    endif

    !Write out effective quantities!
    !    write(*,20) 'rknl [h/Mpc] =',rknl,'rneff=',rneff, 'rncur=',rncur
    !20  format(a15,f12.6,2x,a6,f12.6,2x,a6,f12.6)

    IF(verbose) THEN
       WRITE(*,*) 'HALOFIT_INIT: rknl [h/Mpc] =', rknl
       WRITE(*,*) 'HALOFIT_INIT: rneff =', rneff
       WRITE(*,*) 'HALOFIT_INIT: rncur =', rncur
       WRITE(*,*) 'HALOFIT_INIT: initialised'
       WRITE(*,*)
    END IF

    !ihf=0

  END SUBROUTINE halofit_init

  SUBROUTINE wint(r,sig,d1,d2,a,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT) :: sig
    REAL, INTENT(OUT) :: d1
    REAL, INTENT(OUT) :: d2
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: sum1, sum2, sum3, t, y, x, w1, w2, w3, rk
    INTEGER :: i,nint    

    nint=3000
    sum1=0.d0
    sum2=0.d0
    sum3=0.d0
    DO i=1,nint
       t=(float(i)-0.5)/float(nint)
       y=-1.+1./t
       rk=y
       d2=p_lin(rk,a,cosm)
       x=y*r
       w1=exp(-x*x)
       w2=2*x*x*w1
       w3=4*x*x*(1-x*x)*w1
       sum1=sum1+w1*d2/y/t/t
       sum2=sum2+w2*d2/y/t/t
       sum3=sum3+w3*d2/y/t/t
    END DO
    sum1=sum1/float(nint)
    sum2=sum2/float(nint)
    sum3=sum3/float(nint)
    sig=sqrt(sum1)
    d1=-sum2/sum1
    d2=-sum2*sum2/sum1/sum1 - sum3/sum1

  END SUBROUTINE wint

  SUBROUTINE calculate_halofit_a(k,a,Plin,Pq,Ph,Pnl,n,cosm,verbose,ihf)

    IMPLICIT NONE
    REAL, INTENT(IN) :: k(n)
    REAL, INTENT(IN) :: a
    REAL, INTENT(OUT) :: Plin(n)
    REAL, INTENT(OUT) :: Pq(n)
    REAL, INTENT(OUT) :: Ph(n)
    REAL, INTENT(OUT) :: Pnl(n)
    INTEGER, INTENT(IN) :: n
    TYPE(cosmology), INTENT(INOUT) :: cosm
    LOGICAL, INTENT(IN) :: verbose
    INTEGER, INTENT(IN) :: ihf
    REAL :: rknl, rneff, rncur
    INTEGER :: i

    CALL halofit_init(rknl,rneff,rncur,a,cosm,verbose)

    DO i=1,n
       Plin(i)=p_lin(k(i),a,cosm)
       CALL halofit(k(i),rneff,rncur,rknl,Plin(i),Pnl(i),Pq(i),Ph(i),a,cosm,ihf)
    END DO
    
  END SUBROUTINE calculate_halofit_a

  SUBROUTINE halofit(rk,rn,rncur,rknl,plin,pnl,pq,ph,a,cosm,ihf)

    ! Calculates the HALOFIT power spectrum after rn, rncur and rknl have been pre-calculated
    ! ihf = 1 - Smith et al. 2003
    ! ihf = 2 - Bird et al. 2011
    ! ihf = 3 - Takahashi et al. 2012
    ! ihf = 4 - Takahashi et al. but taken from CAMB with some neutrino stuff
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
    Om_m=cosm%Om_m        ! TODO: Should neutrinos be included?
    Om_mz=Omega_m(a,cosm) ! TODO: Should neutrinos be included?
    Om_vz=Omega_v(a,cosm)+Omega_w(a,cosm) ! Note this well
    wz=w_de(a,cosm) ! Choice here; do you use w or w(z)? w(z) is better I think
    fnu=cosm%Om_nu/cosm%Om_m ! TODO: Does this mean that neutrinos should not be included in Om_m?

    IF(ihf==1) THEN
       ! Smith et al. (2003)
       aa=10**(1.4861+1.8369*rn+1.6762*rn**2+0.7940*rn**3+0.1670*rn**4-0.6206*rncur) ! Smith equation (C9)
       bb=10**(0.9463+0.9466*rn+0.3084*rn**2-0.9400*rncur) ! Smith equation (C10)
       cc=10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur) ! Smith equation (C11)
       gam=0.8649+0.2989*rn+0.1631*rncur ! Smith equation (C12)
       alpha=1.3884+0.3700*rn-0.1452*rn**2 ! Smith equation (C13)
       beta=0.8291+0.9854*rn+0.3401*rn**2 ! Smith equation (C14)
       mu=10**(-3.5442+0.1908*rn) ! Smith equation (C15)
       nu=10**(0.9589+1.2857*rn) ! Smith equation (C16)
    ELSE IF(ihf==2) THEN
       ! Bird et al. (2012); based off Smith et al. (2003)
       aa=10**(1.4861+1.8369*rn+1.6762*rn**2+0.7940*rn**3+0.1670*rn**4-0.6206*rncur)
       bb=10**(0.9463+0.9466*rn+0.3084*rn**2-0.9400*rncur)
       cc=10**(-0.2807+0.6669*rn+0.3214*rn**2-0.0793*rncur)
       gam=0.8649+0.2989*rn+0.1631*rncur+0.316-0.0765*rn-0.835*rncur ! Bird equation (A5)
       alpha=1.3884+0.3700*rn-0.1452*rn**2
       beta=0.8291+0.9854*rn+0.3401*rn**2+fnu*(-6.49+1.44*rn**2) ! Bird equation (A10)
       mu=10**(-3.5442+0.1908*rn)
       nu=10**(0.9589+1.2857*rn)
    ELSE IF(ihf==3) THEN
       ! Takahashi et al. (2012); complete refit, all parameters different from Smith et al. (2003)
       aa=10**(1.5222+2.8553*rn+2.3706*rn**2+0.9903*rn**3+0.2250*rn**4-0.6038*rncur+0.1749*Om_vz*(1.+wz)) ! Takahashi equation (A6)
       bb=10**(-0.5642+0.5864*rn+0.5716*rn**2-1.5474*rncur+0.2279*Om_vz*(1.+wz)) ! Takahashi equation (A7)
       cc=10**(0.3698+2.0404*rn+0.8161*rn**2+0.5869*rncur) ! Takahashi equation (A8)
       gam=0.1971-0.0843*rn+0.8460*rncur ! Takahashi equation (A9)
       alpha=abs(6.0835+1.3373*rn-0.1959*rn**2-5.5274*rncur) ! Takahashi equation (A10; note the ABS
       beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+0.3980*rn**4-0.1682*rncur ! Takahashi equation (A11)
       mu=0. ! Takahashi equation (A12)
       nu=10**(5.2105+3.6902*rn) ! Takahashi equation (A13)
    ELSE IF(ihf==4) THEN
       ! Unpublished CAMB from halofit_ppf.f90; based on Takahashi et al. (2012) plus Bird et al. (2012)
       aa=10**(1.5222+2.8553*rn+2.3706*rn**2+0.9903*rn**3+0.2250*rn**4-0.6038*rncur+0.1749*Om_vz*(1.+wz))
       bb=10**(-0.5642+0.5864*rn+0.5716*rn**2-1.5474*rncur+0.2279*Om_vz*(1.+wz))
       cc=10**(0.3698+2.0404*rn+0.8161*rn**2+0.5869*rncur)
       gam=0.1971-0.0843*rn+0.8460*rncur
       alpha=abs(6.0835+1.3373*rn-0.1959*rn**2-5.5274*rncur) ! Note ABS
       beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+0.3980*rn**4-0.1682*rncur+fnu*(1.081+0.395*rn**2) ! CAMB; halofit_ppf.f90
       mu=0.
       nu=10**(5.2105+3.6902*rn)
    ELSE
       STOP 'HALOFIT: Error, ihf specified incorrectly'
    END IF

    ! Here we need Omega at redshift of interest and calculate the Omega evolution
    ! No HALOFIT extension has modified these numbers
    IF(abs(1.-Om_mz)>0.01) THEN

       ! Open model
       f1a=Om_mz**(-0.0732) 
       f2a=Om_mz**(-0.1423)
       f3a=Om_mz**(0.0725)

       ! Flat LCDM
       f1b=Om_mz**(-0.0307) 
       f2b=Om_mz**(-0.0585)
       f3b=Om_mz**(0.0743)

       ! Linearly interpolate between LCDM and open case for mixed model
       frac=Om_vz/(1.-Om_mz)
       f1=frac*f1b+(1.-frac)*f1a
       f2=frac*f2b+(1.-frac)*f2a
       f3=frac*f3b+(1.-frac)*f3a
       
    ELSE

       ! Set to one for EdS
       f1=1.
       f2=1.
       f3=1.
       
    END IF

    ! Ratio of current wave number to the non-linear wave number
    y=rk/rknl

    IF(ihf==1) THEN
       ! Smith et al. (2003)
       fy=y/4.+y**2/8.                                    ! Smith (below C2)
       ph=aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam)) ! Smith (C4)
       ph=ph/(1.+mu*y**(-1)+nu*y**(-2))                   ! Smith (C3)
       pq=plin*(1.+plin)**beta/(1.+plin*alpha)*exp(-fy)   ! Smith (C2)
    ELSE IF(ihf==2) THEN
       ! Bird et al. (2012)
       fy=y/4.+y**2/8.
       ph=aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam)) ! Bird equation (A2)
       ph=ph/(1.+mu*y**(-1)+nu*y**(-2))                   ! Bird equation (A1)
       Q=fnu*(2.080-12.4*(Om_m-0.3))/(1.+(1.2e-3)*y**3)   ! Bird equation (A6)
       ph=ph*(1.+Q)                                       ! Bird equation (A7)
       pq=plin*(1.+(26.3*fnu*rk**2.)/(1.+1.5*rk**2))      ! Bird equation (A9)
       pq=plin*(1.+pq)**beta/(1.+pq*alpha)*exp(-fy)       ! Bird equation (A8)
    ELSE IF(ihf==3) THEN
       ! Takahashi et al. (2012)
       fy=y/4.+y**2/8.                                    ! Takahashi equation (below A2)
       ph=aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam)) ! Takahashi equation (A3ii)
       ph=ph/(1.+mu*y**(-1)+nu*y**(-2))                   ! Takahashi equation (A3i)
       pq=plin*(1.+plin)**beta/(1.+plin*alpha)*exp(-fy)   ! Takahashi equation (A2)
    ELSE IF(ihf==4) THEN
       ! Unpublished CAMB stuff from halofit_ppf.f90
       fy=y/4.+y**2/8.
       ph=aa*y**(f1*3.)/(1.+bb*y**f2+(f3*cc*y)**(3.-gam))
       ph=ph/(1.+mu*y**(-1)+nu*y**(-2))*(1.+fnu*0.977)    ! CAMB; halofit_ppf.f90; halofit
       pq=plin*(1.+fnu*47.48*rk**2/(1.+1.5*rk**2))        ! CAMB; halofit_ppf.f90; halofit
       pq=plin*(1.+pq)**beta/(1.+pq*alpha)*exp(-fy)
    ELSE
       STOP 'HALOFIT: Error, ihf specified incorrectly'
    END IF

    pnl=pq+ph

  END SUBROUTINE halofit

END MODULE cosmology_functions