# Import statements
import numpy as np
import mead_constants as constants
from scipy.special import erf as erf

# Mass flux of dark matter through spherical shell of radius R [kg/s]
# mu: dark matter wind speed [km/s] 
# sig: dark matter velocity dispersion [km/s] 
# rho: dark matter mass denisty [GeV/cm^3]
# R: radius of shell [au]
def hit_rate(mu,sig,rho,R):
    
    # Convert units
    mu_SI = mu*1e3
    sig_SI = sig*1e3    
    R_SI = R*constants.au
    rho_SI = rho*1e6*constants.eV_mass*1e9
    
    # Calculation
    x = mu/sig # No need for SI here because ratio
    a = np.pi*(1.+x**2)*erf(x/np.sqrt(2.))
    b = np.sqrt(2.*np.pi)*x*np.exp(-x**2/2.)
    p = (a+b)/x
    
    # Result
    return p*rho_SI*sig_SI*R_SI**2

# Unnormalised velocity distributions for dark matter mass flux through a surface [kg/m^2/s]
# v: speed [no units necessary]
# mu: wind speed [no units necessary]
# sig: one-dimensional velocity dispersion [no units necessary]
def velocity_flux_distribution(v,mu,sig):
        return ((v/sig)**2)*np.sinh(v*mu/sig**2)*np.exp(-v**2/(2.*sig**2))
