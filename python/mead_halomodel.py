import numpy as np
import sys

sys.path.append('/Users/Mead/Physics/library/python')
import mead_cosmo as cosmo

class halomod():

    def __init__(self, ihm=1):

        self.dc = 1.686
        self.Dv = 200.

    def halo_mass_function(self, nu):

        return np.sqrt(2./np.pi)*np.exp(-(nu**2)/2.)

    def halo_linear_bias(self, nu):

        return 1.+(nu**2-1.)/self.dc

def Pk_hm(hmod, Ms, ks, Wu, Wv, Pk_lin, Om_m, sigma=None):

    Rs = cosmo.Radius_M(Ms, Om_m)
    nus = cosmo.nu_R(Rs, Pk_lin, dc=hmod.dc)

    #if (sigma == None):
    #    sigma = cosmo.sigma_R(Rs, Pk_lin)

    def P_1h(hmod, Ms, nus, Wu, Wv, Om_m):
        integrand = Wu*Wv*hmod.halo_mass_function(hmod, nus)/Ms
        P_1h = np.trapz(integrand, nus)
        P_1h = P_1h*cosmo.comoving_matter_density(Om_m)
        return P_1h

    def I_2h(hmod, Ms, nus, W, Om_m):
        integrand = W*hmod.linear_halo_bias(nus)*hmod.halo_mass_function(nus)/Ms
        I_2h = np.trapz(integrand, nus)
        I_2h = I_2h*cosmo.comoving_matter_density(Om_m)
        return 1. # Should be I_2h

    def P_2h(hmod, P_lin, k):
        return P_lin(k)*I_2h(hmod)*I_2h(hmod)

    Pk_hm  = np.zeros((len(ks)))
    for ik, k in enumerate(ks):
        Pk_hm[ik] = P_2h(hmod, Pk_lin, k)+P_1h(hmod, Ms, Wu[ik], Wv[ik], Om_m)
    return Pk_hm