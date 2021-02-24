import math
import numpy as np
import sys
import scipy.integrate as integrate

sys.path.append('/Users/Mead/Physics/library/python')
import mead_general as mead
import mead_cosmology as cosmo

# Parameters
Dv0 = 18.*np.pi**2
dc0 = (3./20.)*(12.*np.pi)**(2./3.)
Dv_close_rel_tol = 1e-3
Tinker_z_dep = False

# Halo-model integration scheme
#halo_integration = integrate.trapezoid
halo_integration = integrate.simpson

# W(k) integration scheme
#win_integration = integrate.trapezoid
#win_integration = integrate.simpson
win_integration = integrate.romb

# Halo model types
PS = 'Press & Schecter 1974'
ST = 'Sheth & Tormen 1999'
Tinker2010 = 'Tinker 2010'

# Defaults
hm_def = PS
Dv_def = 200
dc_def = 1.686

class halomod():

    def __init__(self, z, hm=hm_def, Dv=Dv_def, dc=dc_def):

        # Store internal variables
        self.hm = hm
        self.dc = dc
        self.Dv = Dv

        if hm == ST:
            # Sheth-Tormen mass function parameters
            from scipy.special import gamma as Gamma
            p = 0.3
            q = 0.707
            self.p_ST = p
            self.q_ST = q
            self.A_ST = np.sqrt(2.*q)/(np.sqrt(np.pi)+Gamma(0.5-p)/2**p) # A ~ 0.21

        elif hm == Tinker2010:

            # Tinker et al. (2010; https://arxiv.org/abs/1001.3162)

            # Check Delta_v value
            if not math.isclose(Dv, 200., rel_tol=Dv_close_rel_tol):
                print('Warning, Tinker (2010) only coded up for Dv=200')

            # Check delta_c value
            if not math.isclose(dc, 1.686, rel_tol=Dv_close_rel_tol):
                print('Warning, dc = 1.686 assumed in Tinker (2010)')

            # Mass function from Table 4
            alpha = 0.368
            beta = 0.589
            gamma = 0.864
            phi = -0.729
            eta = -0.243
            if Tinker_z_dep:
                beta = beta*(1.+z)**0.20
                gamma = gamma*(1.+z)**-0.01
                phi = phi*(1.+z)**-0.08
                eta = eta*(1.+z)**0.27
            self.alpha_Tinker = alpha
            self.beta_Tinker = beta
            self.gamma_Tinker = gamma
            self.phi_Tinker = phi
            self.eta_Tinker = eta

            # Calibrated halo bias (not from peak-background split) from Table 2
            y = np.log10(self.Dv)
            exp = np.exp(-(4./y)**4)
            self.A_Tinker = 1.+0.24*y*exp
            self.a_Tinker = 0.44*y-0.88
            self.B_Tinker = 0.183
            self.b_Tinker = 1.5
            self.C_Tinker = 0.019+0.107*y+0.19*exp
            self.c_Tinker = 2.4

    def halo_mass_function(self, nu):

        # Halo mass function g(nu) with nu=delta_c/sigma(M)
        # Integral of g(nu) over all nu is unity

        if self.hm == PS:
            return np.sqrt(2./np.pi)*np.exp(-(nu**2)/2.)
        elif self.hm == ST:
            A = self.A_ST
            q = self.q_ST
            p = self.p_ST
            return A*(1.+((q*nu**2)**(-p)))*np.exp(-q*nu**2/2.)
        elif self.hm == Tinker2010:
            alpha = self.alpha_Tinker
            beta = self.beta_Tinker
            gamma = self.gamma_Tinker
            phi = self.phi_Tinker
            eta = self.eta_Tinker
            f1 = (1.+(beta*nu)**(-2.*phi))
            f2 = nu**(2.*eta)
            f3 = np.exp(-gamma*nu**2/2.)
            return alpha*f1*f2*f3
        else:
            raise ValueError('Halo model ihm not recognised in halo_mass_function')

    def linear_halo_bias(self, nu):

        # Halo linear bias b(nu) with nu=delta_c/sigma(M)
        # Integral of b(nu)*g(nu) over all nu is unity

        if self.hm == PS:
            return 1.+(nu**2-1.)/self.dc
        elif self.hm == ST:
            p = self.p_ST
            q = self.q_ST
            return 1.+(q*(nu**2)-1.+2.*p/(1.+(q*nu**2)**p))/self.dc
        elif self.hm == Tinker2010:
            A = self.A_Tinker
            a = self.a_Tinker
            B = self.B_Tinker
            b = self.b_Tinker
            C = self.C_Tinker
            c = self.c_Tinker
            fA = A*nu**a/(nu**a+self.dc**a)
            fB = B*nu**b
            fC = C*nu**c
            return 1.-fA+fB+fC
        else:
            raise ValueError('Halo model ihm not recognised in linear_halo_bias')

def get_nus(Ms, dc, Om_m, sigmas=None, sigma=None, Pk_lin=None):

    # Create arrays of R (Lagrangian) and nu values that correspond to the halo mass
    Rs = cosmo.Radius_M(Ms, Om_m)

    # Convert R values to nu via sigma(R)
    if sigmas is not None:
        nus = dc/sigmas
    elif sigma is not None:
        nus = dc/sigma(Rs) # Else use the provided sigma(R) function
    elif Pk_lin is not None:
        nus = cosmo.nu_R(Rs, Pk_lin, dc) # Uses crappy sigma(R) integration
    else:
        raise ValueError('Error, you need to specify (at least) one of Pk_lin, sigma or sigmas') 
    return nus

def mean_hm(hmod, Ms, Fs, Om_m, sigmas=None, sigma=None, Pk_lin=None):

    # Calculate the mean of some f(M) over halo mass <f>: int f(M)n(M)dM where n(M) = dn/dM in some notations
    # Note that the units of n(M) are [(Msun/h)^{-1} (Mpc/h)^{-3}] so the units of the result are [F (Mpc/h)^{-3}]
    # Common: <M/rho> = 1 over all halo mass (equivalent to int g(nu)dnu = 1)
    # Common: <b(M)M/rho> = 1 over all halo mass (equivalent to int g(nu)b(nu)dnu = 1)
    # Common: <N(M)> with N the number of galaxies in each halo mass; gives mean number density of galaxies
    # Common: <b(M)N(M)>/<N(M)> with N the number of galaxies in each halo mass; gives mean bias of galaxies
    # hmod - halomodel class
    # Ms - Array of halo masses [Msun/h]
    # Fs(Ms) - Array of function to calculate mean density of
    # Om_m - Cosmological matter density at z=0
    # Pk_lin(k) - Function to get linear power at z of interest
    # sigma(R) - Function to get sigma(R) at z of interest
    # nus(M) - Previously calculated array of nu values corresponding to M

    nus = get_nus(Ms, hmod.dc, Om_m, sigmas, sigma, Pk_lin)
    integrand = (Fs/Ms)*hmod.halo_mass_function(nus)
    #return np.trapz(integrand, nus)*cosmo.comoving_matter_density(Om_m)
    return halo_integration(integrand, nus)*cosmo.comoving_matter_density(Om_m)

# Compute the halo window function given a 'density' profile Prho(r) = 4*pi*r^2*rho(r)
# This should almost certainly be done with a dedicated integration routine
def halo_window(ks, rs, Prho):

    # ks: array of wavenumbers [h/Mpc]
    # rs: array of radii (usually going from r=0 to r=rv)
    # Prho[rs]: array of Prho values at different radii

    import scipy.integrate as integrate

    W = np.empty_like(Prho)
    for ik, k in enumerate(ks):
        integrand = np.sinc(k*rs/np.pi)*Prho # Note that numpy sinc function has an odd definition
        if win_integration == integrate.romb:
            dr = (rs[-1]-rs[0])/(len(rs)-1)
            W[ik] = win_integration(integrand, dr)
        elif win_integration in [integrate.trapezoid, integrate.simps]:
            W[ik] = win_integration(integrand, rs)
        else:
            raise ValueError('Halo window function integration method not recognised')
    return W

def Pk_hm(hmod, Ms, ks, rho_uv, Pk_lin, Om_m, beta=None, sigmas=None, sigma=None, low_mass_uv=[False,False], Fourier_uv=[True, True]):

    # TODO: Remove Pk_lin dependence?
    # hmod - halomodel class
    # Ms - Array of halo masses [Msun/h]
    # ks - Array of wavenumbers [h/Mpc]
    # rho_uv(2, Ms, ks/rs) - Array of either Fourier transform of halo profile 'u' and 'v' [u(Mpc/h)^3] or real-space profile from 0 to rv [u]
    # Pk_lin(k) - Function to evaluate the linear power spectrum [(Mpc/h)^3]
    # Om_m - Cosmological matter density
    # beta(M1, M2, k) - Array of beta_NL values at points Ms, Ms, ks
    # sigmas(Ms) - Optional pre-computed array of sigma(M) values corresponding to Ms
    # sigma(R) - Optional function to evaluate the linear sigma(R)
    # low_mass_uv - Should a correction be made for low-mass haloes for field 'u'?
    # Fourier_uv - Are haloes for field 'u' provided in Fourier space or real space?

    # Parameters
    verbose = True

    # Create arrays of R (Lagrangian radius) and nu values that correspond to the halo mass
    nus = get_nus(Ms, hmod.dc, Om_m, sigmas, sigma, Pk_lin)

    # Calculate the missing halo-bias from the low-mass part of the integral if required
    if low_mass_uv[0] or low_mass_uv[1]:
        integrand = hmod.halo_mass_function(nus)*hmod.linear_halo_bias(nus)
        #A = 1.-np.trapz(integrand, nus)
        A = 1.-halo_integration(integrand, nus)
        if verbose:
            print('Missing halo-bias-mass from two-halo integrand:', A)
            print('')

    # Calculate the halo profile Fourier transforms if necessary
    Wuv = np.empty_like(rho_uv)
    for i, _ in enumerate(Fourier_uv):
        if Fourier_uv[i]:
            Wuv[i, :, :] = rho_uv[i, :, :]
        else:
            nr = rho_uv.shape[2] # This is nk, but I suppose it need not be
            for iM, M in enumerate(Ms):
                rv = virial_radius(M, hmod.Dv, Om_m)
                rs = np.linspace(0., rv, nr)
                Wuv[i, iM, :] = halo_window(ks, rs, rho_uv[i, iM, :])

    # Evaluate the integral that appears in the two-halo term
    def I_2h(hmod, Ms, nus, W, Om_m, low_mass):
        integrand = W*hmod.linear_halo_bias(nus)*hmod.halo_mass_function(nus)/Ms
        #I_2h = np.trapz(integrand, nus)
        I_2h = halo_integration(integrand, nus)
        if low_mass:
            I_2h = I_2h+A*W[0]/Ms[0]
        I_2h = I_2h*cosmo.comoving_matter_density(Om_m)
        return I_2h

    # Two-halo term at a specific wavenumber
    def P_2h(hmod, Pk_lin, k, Ms, nus, Wuv, Om_m, low_mass_uv, beta=None):
        if beta is None:
            I_NL = 0.          
        else:
            I_NL = I_beta(hmod, beta, Ms, nus, Wuv, Om_m)
        Iu = I_2h(hmod, Ms, nus, Wuv[0, :], Om_m, low_mass_uv[0])
        Iv = I_2h(hmod, Ms, nus, Wuv[1, :], Om_m, low_mass_uv[1])
        return Pk_lin(k)*(Iu*Iv+I_NL)

    # One-halo term at a specific wavenumber
    def P_1h(hmod, Ms, nus, Wuv, Om_m):
        integrand = Wuv[0, :]*Wuv[1, :]*hmod.halo_mass_function(nus)/Ms
        #P_1h = np.trapz(integrand, nus)
        P_1h = halo_integration(integrand, nus)
        P_1h = P_1h*cosmo.comoving_matter_density(Om_m)
        return P_1h

    # Evaluates the beta_NL integral
    # TODO: Symmetry could be used to speed-up population of integrand by half
    def I_beta(hmod, beta, Ms, nus, Wuv, Om_m):
        integrand = np.zeros((len(nus), len(nus)))
        for iM1, nu1 in enumerate(nus):
            for iM2, nu2 in enumerate(nus):
                M1 = Ms[iM1]
                M2 = Ms[iM2]
                W1 = Wuv[0, iM1]
                W2 = Wuv[1, iM2]
                g1 = hmod.halo_mass_function(nu1)
                g2 = hmod.halo_mass_function(nu2)
                b1 = hmod.linear_halo_bias(nu1)
                b2 = hmod.linear_halo_bias(nu2)
                integrand[iM1, iM2] = beta[iM1, iM2]*W1*W2*g1*g2*b1*b2/(M1*M2)
        return mead.trapz2d(integrand, nus, nus)*cosmo.comoving_matter_density(Om_m)**2

    # Combine everything and return
    Pk_2h_array = np.zeros((len(ks)))
    Pk_1h_array = np.zeros((len(ks)))
    Pk_hm_array = np.zeros((len(ks)))
    for ik, k in enumerate(ks):
        if beta is None:
            Pk_2h_array[ik] = P_2h(hmod, Pk_lin, k, Ms, nus, Wuv[:, :, ik], Om_m, low_mass_uv)
        else:
            Pk_2h_array[ik] = P_2h(hmod, Pk_lin, k, Ms, nus, Wuv[:, :, ik], Om_m, low_mass_uv, beta[:, :, ik])
        Pk_1h_array[ik] = P_1h(hmod, Ms, nus, Wuv[:, :, ik], Om_m)
        Pk_hm_array[ik] = Pk_2h_array[ik]+Pk_1h_array[ik]

    return Pk_2h_array, Pk_1h_array, Pk_hm_array

# Halo virial radius based on the halo mass and overdensity condition
def virial_radius(M, Dv, Om_m):   
    return cosmo.Radius_M(M, Om_m)/mead.cbrt(Dv)

# Isothermal density profile multiplied by 4*pi*r^2
def Prho_isothermal(r, M, rv):
    return M/rv

# NFW density profile multiplied by 4*pi*r^2
def Prho_NFW(r, M, rv, c):
    rs = rv/c
    return M*r/(NFW_factor(c)*(1.+r/rs)**2*rs**2)

def Prho_UPP(r, z, M, r500, cosm):

    alphap = 0.12
    h = cosm.h

    def p(x):

        P0 = 6.41
        c500 = 1.81
        alpha = 1.33
        beta = 4.13
        gamma = 0.31
        y = c500*x

        f1 = y**(2.-gamma)
        f2 = (1.+y**alpha)**(beta-gamma)/alpha
        p = P0*(h/0.7)**(-3./2.)*f1*(r500/c500)**2/f2
        return p

    a = cosmo.scale_factor_z(z)
    H = cosmo.H(cosm, a)

    f1 = 1.65*(h/0.7)**2*H**(8./3.)
    f2 = (M/2.1e14)**(2./3.+alphap)

    return f1*f2*p(r/r500)*4.*np.pi

    


# Converts a Prho profile to a rho profile
# Take care evaluating this at zero (which will give infinity)
def rho_Prho(Prho, r, *args):
    return Prho(r, *args)/(4.*np.pi*r**2)

# Density profile for an isothermal halo
def rho_isothermal(r, M, rv):
    return rho_Prho(Prho_isothermal, r, M, rv)

# Density profile for an NFW halo
def rho_NFW(r, M, rv, c):
    return rho_Prho(Prho_NFW, r, M, rv, c)

# Normalised Fourier tranform for a delta-function profile
def win_delta():
    return 1.

# Normalised Fourier transform for an isothermal profile
def win_isothermal(k, rv):
    from scipy.special import sici as SiCi
    Si, _ = SiCi(k*rv)
    return Si/(k*rv)

# Normalised Fourier transform for an NFW profile
def win_NFW(k, rv, c):
    from scipy.special import sici as SiCi
    rs = rv/c
    kv = k*rv
    ks = k*rs
    Sisv, Cisv = SiCi(ks+kv)
    Sis, Cis = SiCi(ks)
    f1 = np.cos(ks)*(Cisv-Cis)
    f2 = np.sin(ks)*(Sisv-Sis)
    f3 = np.sin(kv)/(ks+kv)
    f4 = NFW_factor(c)
    return (f1+f2-f3)/f4

# Factor from normalisation that always appears in NFW equations
def NFW_factor(c):
    return np.log(1.+c)-c/(1.+c)

def conc_Duffy(M, z):

    M_piv = 2e12 # Pivot mass [Msun/h]
    A = 10.14
    B = -0.081
    C = -1.01

    # Equation (4) in 0804.2486, parameters from 10th row of Table 1
    # Appropriate for the full sample defined via M200
    return A*(M/M_piv)**B*(1.+z)**C

def HOD_Zheng(M, Mmin=1e12, sigma=0.15, M0=1e12, M1=1e13, alpha=1.):
    from scipy.special import erf
    Nc = 0.5*(1.+erf(np.log(M/Mmin)/sigma))
    Ns = Nc*np.heaviside(M-M0, 0.5)*((M-M0)/M1)**alpha
    return Nc, Ns

# LCDM fitting function for the critical linear collapse density from Nakamura & Suto (1997; https://arxiv.org/abs/astro-ph/9612074)
# Cosmology dependence is very weak
def dc_NakamuraSuto(Om_mz):
    return dc0*(1.+0.012299*np.log10(Om_mz))

# LCDM fitting function for virial overdensity from Bryan & Norman (1998; https://arxiv.org/abs/astro-ph/9710107)
# Note that here Dv is defined relative to background matter density, whereas in paper it is relative to critical density
def Dv_BryanNorman(Om_mz):
    x = Om_mz-1.
    Dv = Dv0+82.*x-39.*x**2
    return Dv/Om_mz
