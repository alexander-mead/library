MODULE constants

   USE precision
   
   IMPLICIT NONE

   PUBLIC

   !!

   ! Mathematical constants
   REAL, PARAMETER :: pi = 3.14159265359_dp ! pi
   REAL, PARAMETER :: em = 0.5772156649_dp  ! Euler–Mascheroni
   REAL, PARAMETER :: zero = 0._dp          ! zero
   REAL, PARAMETER :: one = 1._dp           ! one
   REAL, PARAMETER :: degrees = 360._dp     ! Degrees in a circle

   ! Mathematical derived constants
   REAL, PARAMETER :: twopi = 2._dp*pi        ! 2pi or tau ~6.283
   REAL, PARAMETER :: rad2deg = degrees/twopi ! radians-to-degrees conversion
   REAL, PARAMETER :: deg2rad = 1._dp/rad2deg ! degress-to-radians conversion
   REAL, PARAMETER :: Riemann_3 = 1.202056903 ! Riemann Zeta function with n=3

   !!

   !!

   ! Physical constants in SI units
   REAL, PARAMETER :: kB = 1.38064852e-23_dp        ! Boltzmann constant [kg m^2 s^-2 K^-1]
   REAL, PARAMETER :: mp = 1.6726219e-27_dp         ! proton mass [kg]
   REAL, PARAMETER :: me = 9.10938356e-31_dp        ! electron mass [kg]
   REAL, PARAMETER :: bigG = 6.67408e-11_dp         ! Gravitational constant [kg^-1 m^3 s^-2]
   REAL, PARAMETER :: eV = 1.60218e-19_dp           ! electronvolt [kg m^2 s^-2]
   REAL, PARAMETER :: SBconst = 5.670367e-8_dp      ! Steffan-Boltzmann constant [kg s^-3 K^-4]
   REAL, PARAMETER :: c_light = 2.99792458e8_dp     ! speed of light [m/s]
   REAL, PARAMETER :: sigma_T = 6.6524587158e-29_dp ! Thompson-scatter cross section [m^2]
   REAL, PARAMETER :: h_Planck = 6.62607004e-34_dp  ! Planck constant [kg m^2/s]

   ! Derived physical constants in eV units
   REAL, PARAMETER :: mp_eV = mp*(c_light)**2/eV ! Proton rest mass in eV [~938 MeV]
   REAL, PARAMETER :: me_eV = me*(c_light)**2/eV ! Electron rest mass in eV [~511 keV]

   !!

   !!

   ! Astronomy constants
   REAL, PARAMETER :: au = 149597870700._dp  ! Astronomical unit [m] (https://en.wikipedia.org/wiki/Astronomical_unit)
   REAL, PARAMETER :: H0_cos = 100._dp       ! Hubble parameter [km s^-1 (Mpc/h)^-1]
   REAL, PARAMETER :: Msun = 1.98847e30_dp   ! Solar mass [kg] (https://en.wikipedia.org/wiki/Solar_mass)
   REAL, PARAMETER :: Jansky = 1e-26_dp      ! [W m^-2 Hz^-1 / Jy]

   ! Derived astronomy constants
   REAL, PARAMETER :: parsec = 60.*60.*180.*au/pi                          ! Parsec [m] ~3.086e16 m
   REAL, PARAMETER :: kpc = parsec*1e3                                     ! kpc [m] ~3.086e19 m
   REAL, PARAMETER :: Mpc = kpc*1e3                                        ! Mpc [m] ~3.086e22 m
   REAL, PARAMETER :: Gpc = Mpc*1e3                                        ! Gpc [m] ~3.086e25 m
   REAL, PARAMETER :: H0 = H0_cos*1e3/Mpc                                  ! Hubble constant [s^-1] ~3.241e-18 s^-1 (1e3 is km -> m)
   REAL, PARAMETER :: bigG_cos = bigG*Msun/(Mpc*1e3**2)                    ! Gravitational constant [(Msun/h)^-1 (km/s)^2 (Mpc/h)] ~4.301e-9 (1e3**2 m -> km)
   REAL, PARAMETER :: critical_density = 3.*H0**2/(8.*pi*bigG)             ! Critical density [h^2 kg/m^3] ~1.878e-26 [h^2 kg/m^3]
   REAL, PARAMETER :: critical_density_cos = 3.*H0_cos**2/(8.*pi*bigG_cos) ! Universal critical density at (equal to 3*H0^2 / 8piG) [(M_sun/h)/(Mpc/h)^3] ~2.776e11
   REAL, PARAMETER :: Hdist = c_light/(H0_cos*1e3)                         ! Hubble distance (c/H0) [Mpc/h] ~3000 Mpc/h (1e3 km -> m)
   REAL, PARAMETER :: Htime = 1./(H0*60.*60.*24.*365.25*1e9)               ! Hubble time (1/H0) [Gyr/h] ~9.78 Gyr/h
   REAL, PARAMETER :: yfac = sigma_T/(me*c_light**2)                       ! sigma_T/m_e*c^2 [kg^-1 s^2], prefactor of Compton-y integral over *pressure*
   REAL, PARAMETER :: dc0 = (3./20.)*(12.*pi)**(2./3.)                     ! Einstein-de Sitter linear collapse density ~1.686
   REAL, PARAMETER :: Dv0 = 18.*pi**2                                      ! Einsten-de Sitter virialised collapse threshold ~178
   
   !!

END MODULE constants
