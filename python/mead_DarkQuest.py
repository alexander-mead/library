# Standard imports
import numpy as np
import sys

# My imports
#sys.path.append('/Users/Mead/Physics/library/python')
#import mead_cosmology as cosm

# Constants
rhoc = 2.775373e11 # Critical density [(Msun/h) (Mpc/h)^-3]
dc = 1.686 # Collapse threshold for nu definition
Dv = 200. # Spherical-overdensity halo definition
np_min = 200 # Minimum number of halo particles
npart = 2048 # Cube root of number of simulation particles
Lbox_HR = 1000. # Box size for high-resolution simulations [Mpc/h]
Lbox_LR = 2000. # Box size for low-resolution simulations [Mpc/h]

# Minimum and maximum values of cosmological parameters in the emulator
wb_min = 0.0211375
wb_max = 0.0233625
wc_min = 0.10782
wc_max = 0.13178
Om_w_min = 0.54752
Om_w_max = 0.82128
lnAs_min = 2.4752
lnAs_max = 3.7128
ns_min = 0.916275
ns_max = 1.012725
w_min = -1.2
w_max = -0.8

# Parameters
log_interp_sigma = True

def minimum_halo_mass(emu):

    Mbox_HR = comoving_matter_density(emu)*Lbox_HR**3
    mmin = Mbox_HR*np_min/npart**3
    return mmin

def comoving_matter_density(emu):

   Om_m = emu.cosmo.get_Omega0()
   rhom = rhoc*Om_m
   return rhom

def nu_R(emu, R, z):

   M = Mass_R(emu, R)
   nu = nu_M(emu, M, z)
   return nu

def nu_M(emu, M, z):

   nu = dc/sigma_M(emu, M, z)
   return nu

def radius_M(emu, M):

   rhom = comoving_matter_density(emu)
   radius = (3.*M/(4.*np.pi*rhom))**(1./3.)
   return radius

def Mass_R(emu, R):

   rhom = comoving_matter_density(emu)
   Mass = 4.*np.pi*(R**3)*rhom/3.
   return Mass

def Mass_nu(emu, nu, z):

   # Import
   from scipy.interpolate import InterpolatedUnivariateSpline as ius

   # Options
   log_interp = log_interp_sigma # Should sigma(M) be interpolated logarithmically?
   
   # Get internal M vs sigma arrays
   Ms_internal = emu.massfunc.Mlist
   sig0s_internal = emu.massfunc.sigs0
   sigs_internal = sig0s_internal*emu.Dgrowth_from_z(z)
   nus_internal = dc/sigs_internal 

   # Make an interpolator for sigma(M)  
   if (log_interp):
      mass_interpolator = ius(nus_internal, np.log(Ms_internal))
   else:
      mass_interpolator = ius(nus_internal, Ms_internal)

   # Get sigma(M) from the interpolator at the desired masses
   if (log_interp):
      Mass = np.exp(mass_interpolator(nu))
   else:
      Mass = mass_interpolator(nu)

   return Mass

def sigma_R(emu, R, z):

   M = Mass_R(emu, R)
   sigma = sigma_M(emu, M, z)
   return sigma

def sigma_M(emu, M, z):

   # Import
   from scipy.interpolate import InterpolatedUnivariateSpline as ius

   # Options
   log_interp = log_interp_sigma # Should sigma(M) be interpolated logarithmically?
   
   # Get internal M vs sigma arrays
   Ms_internal = emu.massfunc.Mlist
   sigs_internal = emu.massfunc.sigs0

   # Make an interpolator for sigma(M)  
   if (log_interp):
      sigma_interpolator = ius(np.log(Ms_internal), np.log(sigs_internal))
   else:
      sigma_interpolator = ius(Ms_internal, sigs_internal)
   
   # Get sigma(M) from the interpolator at the desired masses
   if (log_interp):
      sigma0 = np.exp(sigma_interpolator(np.log(M)))
   else:
      sigma0 = sigma_interpolator(M)

   # Growth function (g(z=0)=1)
   g = emu.Dgrowth_from_z(z) 
   sigma = g*sigma0

   # Result assuming scale-independent growth
   return sigma

# Beta_NL function
def beta_NL(emu, vars, ks, z, var='Mass'):
    
   # Parameters
   # ibias = 1: Linear bias from emulator
   # ibias = 2: Linear bias from auto-spectrum at large wavenumber
   # ibias = 3: Linear bias from cross-spectrum at large wavenumber
   ibias = 2              # Choose source of linear bias
   force_BNL_zero = False # Should beta be forced to zero at large scales?
   klin = 0.02            # Large 'linear' scale [h/Mpc]

   # Set array name sensibly
   if (var == 'Mass'):
      Ms = vars
   elif (var == 'Radius'):
      Rs = vars
      Ms = Mass_R(emu, Rs)
   elif (var == 'nu'):
      nus = vars
      Ms = Mass_nu(emu, nus, z)
   else:
      print('Error, variable for beta_NL not recognised')
   
   # klin must be a numpy array for this to work later
   klin = np.array([klin]) 
   
   # Linear power
   Pk_lin = emu.get_pklin_from_z(ks, z)
   
   # Calculate beta_NL by looping over mass arrays
   beta = np.zeros((len(Ms), len(Ms), len(ks)))  
   for im1, M1 in enumerate(Ms):
      for im2, M2 in enumerate(Ms):
         
         # Halo-halo power spectrum
         Pk_hh = emu.get_phh_mass(ks, M1, M2, z)
         
         # Linear halo bias
         if (ibias == 1):
               b1 = emu.get_bias_mass(M1, z)[0]
               b2 = emu.get_bias_mass(M2, z)[0]
         elif (ibias == 2):
               b1 = np.sqrt(emu.get_phh_mass(klin, M1, M1, z)/emu.get_pklin_from_z(klin, z))
               b2 = np.sqrt(emu.get_phh_mass(klin, M2, M2, z)/emu.get_pklin_from_z(klin, z))
         elif (ibias == 3):
               b1 = emu.get_phm_mass(klin, M1, z)/emu.get_pklin_from_z(klin, z)
               b2 = emu.get_phm_mass(klin, M2, z)/emu.get_pklin_from_z(klin, z)
               
         # Create beta_NL
         beta[im1, im2, :] = Pk_hh/(b1*b2*Pk_lin)-1.

         # Force Beta_NL to be zero at large scales if necessary
         if (force_BNL_zero):
            Pk_hh0 = emu.get_phh_mass(klin, M1, M2, z)
            Pk_lin0 = emu.get_pklin_from_z(klin, z)
            db = Pk_hh0/(b1*b2*Pk_lin0)-1.
            #beta[im1, im2, :] = beta[im1, im2, :]-db # Additive correction
            beta[im1, im2, :] = (beta[im1, im2, :]+1.)/(db+1.)-1. # Multiplicative correction
         
   return beta 