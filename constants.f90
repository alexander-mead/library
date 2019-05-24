MODULE constants

  IMPLICIT NONE

  !!

  ! Mathematical constants
  REAL, PARAMETER :: pi=3.14159265359 ! pi  
  REAL, PARAMETER :: em=0.5772156649  ! Euler–Mascheroni
  REAL, PARAMETER :: zero=0.          ! zero
  REAL, PARAMETER :: one=1.           ! one
  REAL, PARAMETER :: degrees=360.     ! Degrees in a circle

  ! Mathematical derived  
  REAL, PARAMETER :: twopi=2.*pi           ! 2pi or tau ~6.283
  REAL, PARAMETER :: rad2deg=degrees/twopi ! radians-to-degrees conversion
  REAL, PARAMETER :: deg2rad=1./rad2deg    ! degress-to-radians conversion

  !!

  !!

  ! Physical constants
  REAL, PARAMETER :: kB=1.38064852e-23        ! Boltzmann constant [m^2 kg s^-2 K^-1]  
  REAL, PARAMETER :: mp=1.6726219e-27         ! proton mass [kg]
  REAL, PARAMETER :: me=9.10938356e-31        ! electron mass [kg]
  REAL, PARAMETER :: bigG=6.67408e-11         ! Gravitational constant [kg^-1 m^3 s^-2]
  REAL, PARAMETER :: eV=1.60218e-19           ! electronvolt [kg m^2 s^-2]  
  REAL, PARAMETER :: SBconst=5.670367e-8      ! Steffan-Boltzmann constant [kg s^-3 K^-4]
  REAL, PARAMETER :: c_light=2.99792458e8     ! speed of light [m/s]
  REAL, PARAMETER :: sigma_T=6.6524587158e-29 ! Thompson-scatter cross section [m^2]
  REAL, PARAMETER :: h_Planck=6.62607004e-34  ! Planck constant [kg m^2/s]

  ! Derived physical constants
  REAL, PARAMETER :: mp_eV=mp*(c_light)**2/eV ! Proton rest mass in eV [~938 MeV]
  REAL, PARAMETER :: me_eV=me*(c_light)**2/eV ! Electron rest mass in eV [~511 keV]
 
  !!

  !!
  
  ! Astronomy constants

  REAL, PARAMETER :: au=149597870700.               ! Astronomical unit [m] (https://en.wikipedia.org/wiki/Astronomical_unit)  
  REAL, PARAMETER :: H0_cos=100.                    ! Hubble parameter [km s^-1 (Mpc/h)^-1]  
  REAL, PARAMETER :: Msun=1.98847e30                ! Solar mass [kg] (https://en.wikipedia.org/wiki/Solar_mass)
  REAL, PARAMETER :: Jansky=1e-26                   ! [W m^-2 Hz^-1 / Jy]
  !REAL, PARAMETER :: H0=3.243e-18                   ! Hubble parameter, H0 [s^-1]
  !REAL, PARAMETER :: Mpc=3.086e22                   ! Mpc [m]  
  !REAL, PARAMETER :: SI_to_Jansky=1e26              ! [W m^-2 Hz^-1 / Jy]

  ! Derived astronomy constants
  REAL, PARAMETER :: pc=60.*60.*180*au/pi                           ! Parsec [m] ~3.086e16
  REAL, PARAMETER :: kpc=pc*1e3                                     ! kpc [m]
  REAL, PARAMETER :: Mpc=kpc*1e3                                    ! Mpc [m]
  REAL, PARAMETER :: Gpc=Mpc*1e3                                    ! Gpc [m]
  REAL, PARAMETER :: H0=H0_cos*1e3/Mpc                              ! Hubble constant [s^-1] ~3.243e-18
  REAL, PARAMETER :: bigG_cos=bigG*Msun/(1e6*Mpc)                   ! Gravitational constant [(Msun/h)^-1 (km/s)^2 (Mpc/h)]
  REAL, PARAMETER :: critical_density=3.*H0_cos**2/(8.*pi*bigG_cos) ! Universal critical density at (equal to 3*H0^2 / 8piG) [(M_sun/h)/(Mpc/h)^3] ~2.776e11
  REAL, PARAMETER :: Hdist=c_light/(H0_cos*1e3)                     ! Hubble distance (c/H0) [Mpc/h] ~3000
  REAL, PARAMETER :: Htime=1./(H0*60*60*24*365.25*1e9)              ! Hubble time (1/H0) [Gyrs/h]    ~9.78
  REAL, PARAMETER :: yfac=sigma_T/(me*c_light**2)                   ! sigma_T/m_e*c^2 [kg^-1 s^2], prefactor of Compton-y integral over *pressure*
  REAL, PARAMETER :: dc0=(3./20.)*(12.*pi)**(2./3.)                 ! Einstein-de Sitter linear collapse density ~1.686
  REAL, PARAMETER :: Dv0=18.*pi**2                                  ! Einsten-de Sitter virialised collapse threshold ~178

   ! TODO: Relate these to fundamental constants
  REAL, PARAMETER :: neutrino_constant=94.1         ! Critical mass for neutrino density to close Universe [eV] (or is it 93.03 eV, or 93.14 eV; https://arxiv.org/pdf/1812.02102.pdf)
  !REAL, PARAMETER :: neff_contribution=0.227        ! Contribution to Omega_r per n_eff (7/8)*(4/11)**(4/3)
  !REAL, PARAMETER :: critical_density=2.7755e11     ! Universal critical density at (equal to 3*H0^2 / 8piG) [(M_sun/h)/(Mpc/h)^3]
  !REAL, PARAMETER :: Htime=9.7776                   ! Hubble time (1/H0) [Gyrs/h]

  !!
  
END MODULE constants
