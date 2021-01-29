# Standard
import numpy as np
import sys

# Other imports
import camb

# My imports
sys.path.append('/Users/Mead/Physics/library/python')
import mead_constants as const
import mead_cosmology as cosmo

# Create my cosmological parameters from a CAMB set
# TODO: I guessed at the names of the cosmological parameters here. Probably wrong.
def create_cosmology(pars, verbose=False):

    # CAMB cosmological parameters
    results = camb.get_results(pars)

    # Convert to my cosmology   
    Om_w = results.Params.Om_v
    Om_m = results.Params.Om_m
    h = results.Params.H0/100.
    Om_b = results.Params.Om_b
    As = results.Params.As
    ns = results.Params.ns
    w = results.Params.w
    wa = results.Params.wa
    m_nu = results.Params.m_nu
    cosm = cosmo.cosmology(Om_m=Om_m, Om_b=Om_b, Om_w=Om_w, h=h, As=As, ns=ns, w=w, wa=wa, m_nu=m_nu)

    # Print to screen
    if (verbose):
        cosm.print()

    return cosm

# Convert my cosmology into a CAMB cosmology
# TODO: Neutrino hierarchy? CMB temperature? Helium fraction? Tau?
def convert_mead_cosmology(cosm, verbose=False):

    # Set up a new default set of parameters for CAMB
    pars = camb.CAMBparams()

    # Get CAMB cosmology from my structure
    wb = cosm.w_b
    wc = cosm.w_c
    h = cosm.h
    As = cosm.As
    ns = cosm.ns
    w = cosm.w
    wa = cosm.wa
    Om_k = cosm.Om_k
    m_nu = cosm.m_nu

    # This function sets standard and helium set using BBN consistency
    pars.set_cosmology(ombh2=wb, omch2=wc, H0=100.*h, mnu=m_nu, omk=Om_k)
    pars.set_dark_energy(w=w, wa=wa, dark_energy_model='ppf') 
    pars.InitPower.set_params(As=As, ns=ns, r=0)

    # Print out the full list of cosmological parameters
    if (verbose): 
        print(pars)

    return pars

# Get a 2D P(k,z) array for the linear or non-linear power spectrum
def Pk_mm(ks, zs, pars, nonlinear=None):

    from camb import model

    # Parameters
    kmax = 100.

    if (nonlinear == None):
        nl = False
        nlmod = 'mead2020'
    else:
        nl = True
        nlmod = nonlinear

    # Get non-linear power from CAMB
    pars.set_matter_power(redshifts=zs, kmax=kmax)
    pars.NonLinearModel.set_params(halofit_version=nlmod)
    #results = camb.get_results(pars)
    #ks_CAMB, _, Pk_CAMB = results.get_nonlinear_matter_power_spectrum(params=pars)

    # Create a power spectrum interpolation object
    Pk_CAMB_interp = camb.get_matter_power_interpolator(pars, 
                                                        nonlinear = nl, 
                                                        hubble_units = True, 
                                                        k_hunit = True, 
                                                        kmax = kmax,
                                                        var1 = model.Transfer_tot,
                                                        var2 = model.Transfer_tot, 
                                                        zmax = max(zs),
                                                        )

    # Evaluate interpolator at the desired k and z points
    Pk_CAMB = Pk_CAMB_interp.P(zs, ks)

    return Pk_CAMB

