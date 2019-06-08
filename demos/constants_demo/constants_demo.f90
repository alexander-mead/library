PROGRAM constants_test

  USE constants
  
  IMPLICIT NONE

  WRITE(*,*)

  WRITE(*,*) 'Mathematical constants'
  WRITE(*,*) '======================'
  WRITE(*,*) 'pi:', pi
  WRITE(*,*) 'EM:', em
  WRITE(*,*) 'Zero:', zero
  WRITE(*,*) 'One:', one
  WRITE(*,*) 'Degrees in a circle:', degrees
  WRITE(*,*) 'Two pi:', twopi
  WRITE(*,*) 'Radian [deg]:', rad2deg  
  WRITE(*,*) 'Degree [rad]:', deg2rad
  WRITE(*,*)

  WRITE(*,*) 'Physical constants'
  WRITE(*,*) '=================='
  WRITE(*,*) 'Boltzmann [m^2 kg s^-2 K^-1]:', kB
  WRITE(*,*) 'Proton mass [kg]:', mp
  WRITE(*,*) 'Electron mass [kg]:', me
  WRITE(*,*) 'Graviational constant [kg^-1 m^3 s^-2]:', bigG
  WRITE(*,*) 'Electron volt [kg m^2 s^-2]:', eV
  WRITE(*,*) 'Steffan-Boltzmann [kg s^-3 K^-4]:', SBconst
  WRITE(*,*) 'Speed of light [m/s]:', c_light
  WRITE(*,*) 'Thompson scattering [m^2]:', sigma_T
  WRITE(*,*) 'Planck [kg m^2/s]:', h_Planck
  WRITE(*,*) 'Proton mass [eV]:', mp_eV
  WRITE(*,*) 'Electron mass [eV]:', me_eV
  WRITE(*,*)

  WRITE(*,*) 'Astronomy constants'
  WRITE(*,*) '==================='
  WRITE(*,*) 'Astronomical unit [m]:', au
  WRITE(*,*) 'Hubble constant [km/s (Mpc/h)^-1]:', H0_cos
  WRITE(*,*) 'Solar mass [kg]:', Msun
  WRITE(*,*) 'Jansky [W m^-2 Hz^-1]:', Jansky
  WRITE(*,*) 'Parsec [m]:', pc
  WRITE(*,*) 'Kiloparsec [m]:', kpc
  WRITE(*,*) 'Megaparsec [m]:', Mpc
  WRITE(*,*) 'Gigaparsec [m]:', Gpc
  WRITE(*,*) 'Hubble constant [s^-1]:', H0
  WRITE(*,*) 'Gravitational constant [(Msun/h)^-1 (km/s)^2 (Mpc/h)]:', bigG_cos
  WRITE(*,*) 'Critical density [(Msun/h) (Mpc/h)^-3]:', critical_density
  WRITE(*,*) 'Hubble distance [Mpc/h]:', Hdist
  WRITE(*,*) 'Hubble time [Gyrs/h]:', Htime
  WRITE(*,*) 'y factor [kg^-1 s^2]:', yfac
  WRITE(*,*) 'Linear-collapse threshold:', dc0
  WRITE(*,*) 'Virial overdensity:', Dv0
  WRITE(*,*) 'Neutrino constant [eV]:', neutrino_constant
  WRITE(*,*) 'Neff contribution per nu:', neff_contribution
  WRITE(*,*)
  
END PROGRAM constants_test
