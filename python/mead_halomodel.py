import numpy as np
import sys

sys.path.append('/Users/Mead/Physics/library/python')
import mead_general as mead
import mead_cosmology as cosmo

class halomod():

    def __init__(self, ihm=1):

        self.dc = 1.686
        self.Dv = 200.

    def halo_mass_function(self, nu):

        return np.sqrt(2./np.pi)*np.exp(-(nu**2)/2.)

    def linear_halo_bias(self, nu):

        return 1.+(nu**2-1.)/self.dc

def Pk_hm(hmod, Ms, ks, rhou, rhov, Pk_lin, Om_m, sigma=None, low_mass_u=True, low_mass_v=True, Fourier_u=True, Fourier_v=True):

    # hmod - halomodel class
    # Ms - Array of halo masses [Msun/h]
    # ks - Array of wavenumbers [h/Mpc]
    # rhou(Ms, ks/rs) - Either Fourier transform of halo profile 'u' [u(Mpc/h)^3] or real-space profile from 0 to rv [u]
    # rhov(Ms, ks/rs) - Either Fourier transform of halo profile 'v' [v(Mpc/h)^3] or real-space profile from 0 to rv [v]
    # Pk_lin(k) - Function to evaluate the linear power spectrum [(Mpc/h)^3]
    # Om_m - Cosmological matter density
    # sigma(R) - Function to evaluate the linear sigma(R)
    # low_mass_u - Should a correction be made for low-mass haloes for field 'u'?
    # low_mass_v - Should a correction be made for low-mass haloes for field 'v'?
    # Fourier_u - Are haloes for field 'u' provided in Fourier space or real space?
    # Fourier_v - Are haloes for field 'v' provided in Fourier space or real space?

    verbose = True

    # Create arrays of R (Lagrangian) and nu values that correspond to the halo mass
    Rs = cosmo.Radius_M(Ms, Om_m)
    nus = cosmo.nu_R(Rs, Pk_lin, dc=hmod.dc)

    # Calculate the missing halo-bias from the low-mass part of the integral if required
    if low_mass_u or low_mass_v:
        integrand = hmod.halo_mass_function(nus)*hmod.linear_halo_bias(nus)
        A = 1.-np.trapz(integrand, nus)
        if verbose:
            print('Missing halo-bias-mass from two-halo integrand:', A)
            print('')

    #if (sigma == None):
    #    sigma = cosmo.sigma_R(Rs, Pk_lin)

    # Compute the halo window function given a 'density' profile
    def halo_window(Ms, rs, rho):
        W = np.empty_like(rho)
        for iM, M in enumerate(Ms):
            rv = virial_radius(M, hmod.Dv, Om_m)
            rs = np.linspace(0., rv)
            for ik, k in enumerate(ks):
                integrand = rs*np.sin(k*rs)*rho[iM, :]
                W[iM, ik] = (4.*np.pi/k)*np.trapz(integrand, rs)
        return W

    # Calculate the halo profile Fourier transforms if necessary
    if Fourier_u:
        Wu = rhou
    else:
        rs = virial_radius(Ms, hmod.Dv, Om_m)
        Wu = halo_window(Ms, rs, rhou)
        raise ValueError('Non-Fourier profiles not yet supported')

    # Calculate the halo profile Fourier transforms if necessary    
    if Fourier_v:
        Wv = rhov
    else:
        rs = virial_radius(Ms, hmod.Dv, Om_m)
        Wv = halo_window(Ms, rs, rhov)
        raise ValueError('Non-Fourier profiles not yet supported')

    # Integral that appears in the two-halo term
    def I_2h(hmod, Ms, nus, W, Om_m, low_mass=False):
        integrand = W*hmod.linear_halo_bias(nus)*hmod.halo_mass_function(nus)/Ms
        I_2h = np.trapz(integrand, nus)
        if (low_mass):
            I_2h = I_2h+A*W[0]/Ms[0]
        I_2h = I_2h*cosmo.comoving_matter_density(Om_m)     
        return I_2h

    # Evaluate the two-halo term
    def P_2h(hmod, P_lin, k, Ms, nus, Wu, Wv, Om_m, low_mass_u=False, low_mass_v=False):
        return P_lin(k)*I_2h(hmod, Ms, nus, Wu, Om_m, low_mass_u)*I_2h(hmod, Ms, nus, Wv, Om_m, low_mass_v)

    # Evaluate the one-halo term
    def P_1h(hmod, Ms, nus, Wu, Wv, Om_m):
        integrand = Wu*Wv*hmod.halo_mass_function(nus)/Ms
        P_1h = np.trapz(integrand, nus)
        P_1h = P_1h*cosmo.comoving_matter_density(Om_m)
        return P_1h

    Pk_2h_array = np.zeros((len(ks)))
    Pk_1h_array = np.zeros((len(ks)))
    Pk_hm_array = np.zeros((len(ks)))
    for ik, k in enumerate(ks):
        Pk_2h_array[ik] = P_2h(hmod, Pk_lin, k, Ms, nus, Wu[:, ik], Wv[:, ik], Om_m, low_mass_u, low_mass_v)
        Pk_1h_array[ik] = P_1h(hmod, Ms, nus, Wu[:, ik], Wv[:, ik], Om_m)
        Pk_hm_array[ik] = Pk_2h_array[ik]+Pk_1h_array[ik]
    return Pk_2h_array, Pk_1h_array, Pk_hm_array

def virial_radius(M, Dv, Om_m):
    # Calculate the halo virial radius based on the halo mass and overdensity condition
    return cosmo.Radius_M(M, Om_m)/mead.cbrt(Dv)