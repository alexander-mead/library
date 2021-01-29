# Standard imports
import numpy as np
import sys

# Other imports
sys.path.append('/Users/Mead/Physics/DarkQuest')
import darkemu

# My imports
sys.path.append('/Users/Mead/Physics/library/python')
import mead_constants as const
import mead_cosmology as cosmo

# Constants
dc = 1.686      # Collapse threshold for nu definition
Dv = 200.       # Spherical-overdensity halo definition
np_min = 200    # Minimum number of halo particles
npart = 2048    # Cube root of number of simulation particles
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

class cosmology():

    def __init__(self, wb=0.02225, wc=0.1198, Om_w=0.6844, lnAs=3.094, ns=0.9645, w=-1.):

        # Primary parameters
        self.wb = wb
        self.wc = wc
        self.Om_w = Om_w
        self.lnAs = lnAs
        self.ns = ns
        self.w = w

        # Fixed parameters
        self.wnu = 0.00064
        self.Om_k = 0.

        # Derived parameters
        self.wm = self.wc+self.wb+self.wnu
        self.Om_m = 1.-self.Om_w # Flatness
        self.h = np.sqrt(self.wm/self.Om_m)
        self.Om_b = self.wb/self.h**2
        self.Om_c = self.wc/self.h**2
        self.As = np.exp(self.lnAs)/1e10
        self.m_nu = self.wnu*const.nuconst
        self.Om_nu = self.wnu/self.h**2

    def print(self):

        # Write primary parameters to screen
        print('Dark Quest primary parameters')
        print('omega_b: %1.4f' % (self.wb))
        print('omega_c: %1.4f' % (self.wc))  
        print('Omega_w: %1.4f' % (self.Om_w))
        print('As [1e9]: %1.4f' % (self.As*1e9))
        print('ns: %1.4f' % (self.ns))
        print('w: %1.4f' % (self.w))
        print()

        #print('Dark Quest fixed parameters')
        #print('omega_nu: %1.4f' % (self.wnu))
        #print('Omega_k: %1.4f' % (self.Om_k))
        #print()

        # Write derived parameters to screen
        print('Dark Quest derived parameters')
        print('Omega_m: %1.4f' % (self.Om_m))
        print('Omega_b: %1.4f' % (self.Om_b))      
        print('omega_m: %1.4f' % (self.wm))
        print('h: %1.4f' % (self.h))      
        print('Omegea_c: %1.4f' % (self.Om_c))
        print('Omega_nu: %1.4f' % (self.Om_nu))      
        print('m_nu [eV]: %1.4f' % (self.m_nu))
        print()

# Random cosmological parameters from the Dark Quest hypercube
def random_cosmology():

    wb = np.random.uniform(wb_min, wb_max)
    wc = np.random.uniform(wc_min, wc_max)
    Om_w = np.random.uniform(Om_w_min, Om_w_max)
    lnAs = np.random.uniform(lnAs_min, lnAs_max)
    ns = np.random.uniform(ns_min, ns_max)
    w = np.random.uniform(w_min, w_max)

    cpar = cosmology(wb=wb, wc=wc, Om_w=Om_w, lnAs=lnAs, ns=ns, w=w)

    return cpar

# Create a set of my cosmological parameters from a Dark Quest set
def create_mead_cosmology(cpar, verbose=False):

    # Make a mead cosmology
    cosm = cosmo.cosmology(Om_m=cpar.Om_m, Om_b=cpar.Om_b, Om_w=cpar.Om_w, h=cpar.h, 
                           As=cpar.As, ns=cpar.ns, w=cpar.w, m_nu=cpar.m_nu)

    # Print to screen
    if (verbose):
        cosm.print()

    return cosm

# Convert my cosmology into a Dark Quest cosmology
def convert_mead_cosmology(cosm):

    # Get Dark Quest parameters from my structure
    wb = cosm.w_b
    wc = cosm.w_c
    Om_w = cosm.Om_w
    lnAs = np.log(cosm.As*1e10)
    ns = cosm.ns
    w = cosm.w
    cpar = cosmology(wb=wb, wc=wc, Om_w=Om_w, lnAs=lnAs, ns=ns, w=w)

    return  cpar

# Initialise the emulator for a given set of cosmological parameters
def init_emulator(cpar):

    # Start Dark Quest
    print('Initialize Dark Quest')
    emu = darkemu.base_class()
    print('')

    # Initialise emulator
    cparam = np.array([cpar.wb, cpar.wc, cpar.Om_w, cpar.lnAs, cpar.ns, cpar.w])   
    cpar.print()
    emu.set_cosmology(cparam) # I think that this does a load of emulator init steps

    return emu

# Matter power spectrum
def Pk_mm(emu, ks, zs, nonlinear=False):

    if (isinstance(zs, float)):
        if (nonlinear):
            Pk = emu.get_pmnl(ks, zs)
        else:         
            Pk = emu.get_pklin_from_z(ks, zs)
    else:
        Pk = np.zeros((len(zs), len(ks)))
        for iz, z in enumerate(zs):
            if (nonlinear):
                Pk[iz, :] = emu.get_pmnl(ks, z)
            else:         
                Pk[iz, :] = emu.get_pklin_from_z(ks, z)

    return Pk

# Minimum halo mass for the set of cosmological parameters
def minimum_halo_mass(emu):

    Mbox_HR = comoving_matter_density(emu)*Lbox_HR**3
    mmin = Mbox_HR*np_min/npart**3
    return mmin

    # Comoving matter density
def comoving_matter_density(emu):

    Om_m = emu.cosmo.get_Omega0()
    #rhom = const.rhoc*Om_m
    rhom = cosmo.comoving_matter_density(Om_m)
    return rhom

def nu_R(emu, R, z):

    M = Mass_R(emu, R)
    nu = nu_M(emu, M, z)
    return nu

def nu_M(emu, M, z):

    nu = dc/sigma_M(emu, M, z)
    return nu

def Radius_M(emu, M):

    #rhom = comoving_matter_density(emu)
    #radius = (3.*M/(4.*np.pi*rhom))**(1./3.)
    Om_m = emu.cosmo.get_Omega0()
    radius = cosmo.Radius_M(M, Om_m)
    return radius

def Mass_R(emu, R):

    #rhom = comoving_matter_density(emu)
    #Mass = 4.*np.pi*(R**3)*rhom/3.
    Om_m = emu.cosmo.get_Omega0()
    Mass = cosmo.Mass_R(R, Om_m)
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
    
    # Source of linear halo bias
    # ibias = 1: Linear bias from emulator
    # ibias = 2: Linear bias from auto-spectrum at large wavenumber
    # ibias = 3: Linear bias from cross-spectrum at large wavenumber
    ibias = 2

    # Force Beta_NL to zero?
    # ibias = 0: No
    # ibias = 1: Yes, via addative correction
    # ibias = 3: Yes, via multiplicative correction
    force_BNL_zero = 0

    # Parameters
    klin = 0.02  # Large 'linear' scale [h/Mpc]

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
            if (force_BNL_zero != 0):
                Pk_hh0 = emu.get_phh_mass(klin, M1, M2, z)
                Pk_lin0 = emu.get_pklin_from_z(klin, z)
                db = Pk_hh0/(b1*b2*Pk_lin0)-1.
                if (force_BNL_zero == 1):
                    beta[im1, im2, :] = beta[im1, im2, :]-db # Additive correction
                elif (force_BNL_zero == 2):
                    beta[im1, im2, :] = (beta[im1, im2, :]+1.)/(db+1.)-1. # Multiplicative correction
                else:
                    raise ValueError('force_BNL_zero not set correctly')
         
    return beta 

def calculate_rescaling_params(emu_ori, emu_tgt, z_tgt, M1_tgt, M2_tgt):
    
    R1_tgt = Radius_M(emu_tgt, M1_tgt)
    R2_tgt = Radius_M(emu_tgt, M2_tgt)
    
    s, z = cosmo.calculate_AW10_rescaling_parameters(z_tgt, R1_tgt, R2_tgt, 
                                                     lambda Ri, zi: sigma_R(emu_ori, Ri, zi), 
                                                     lambda Ri, zi: sigma_R(emu_tgt, Ri, zi)
                                                    )
    return s, z