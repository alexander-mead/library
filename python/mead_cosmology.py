import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import mead_interpolation as interpolation

## Cosmological constants ##

conH0 = 2998. # c/H0 in Mpc/h
invH0 = 9.778 # 1/H0 in Gyr/h
H0 = 100. # H0 in h km/s/Mpc
T0 = 2.73 # CMB temperature in Kelvin
tau = 0. # Optical depth in dimensionless units

## ##

## Assign and initialise ##

def assign_cosmology():

    global Om_r, Om_m, Om_v, Om_w
    global ide, w0, wa, a1, a2, nw

    print('Assign_cosmology: Setting cosmological parameters')

    # Cosmological parameters
    Om_r = 0.0
    Om_m = 0.3
    Om_v = 0.0
    Om_w = 0.7
    ide = 0

    # ide=1,2 (1 - w(a)CDM; 2 - wCDM)
    w0 = -1.
    wa = 0.

    # ide=5 (IDE II)
    a1 = 1e-2
    a2 = 0.0776
    nw = 3

    print('Assign_cosmology: Done')
    print()
    
def initialise_cosmology():

    global Om
    global a_tab, a_lintab, a_logtab
    global amin, amax
    
    # Derived parameters
    Om = Om_r+Om_m+Om_v+Om_w

    # Set 'a' range for integrations
    amin = 1e-5
    amax = 1.
    na = 64

    # Set the range for 'a' (not sure if linear or log is better here)
    a_logtab = np.logspace(np.log10(amin),np.log10(amax),na)
    a_lintab = np.linspace(0.,amax,na)
    #a_lintab_nozero = np.delete(a_lintab,0) #Remove the a=0 entry from the linear 'a' table

    # Set the table of values of 'a' to use for interpolation functions
    a_tab = a_logtab

    print('Initialise_cosmology: min -> max in a_tab:', a_tab[0], '->', a_tab[-1])
    print('Initialise_cosmology: min -> max in a_logtab:', a_logtab[0], '->', a_logtab[-1])
    print('Initialise_cosmology: min -> max in a_lintab:', a_lintab[0], '->', a_lintab[-1])
    #print('min -> max in a_lintab_nozero:', a_lintab_nozero[0],a_lintab_nozero[-1])
    print()

    #Add check for big bang by direct integration of Raychaudhuri equation
    
    #initialise_distances()
    #initialise_growth()

def print_cosmology():

    print('Print_cosmology: Parameters')
    print('Print_cosmology: Om_r:', Om_r)
    print('Print_cosmology: Om_m:', Om_m)
    print('Print_cosmology: Om_v:', Om_v)
    print('Print_cosmology: Om_w:', Om_w)
    print('Print_cosmology: Om:', Om)
    if(ide==0):
        print('Print_cosmology: Dark energy: LCDM')
    elif(ide==1):
        print('Print_cosmology: Dark energy: w(a)CDM')
        print('Print_cosmology: w0:', w0)
        print('Print_cosmology: wa:', wa)
    elif(ide==2):
        print('Print_cosmology: Dark enrgy: wCDM')
        print('Print_cosmology: w:', w0)
    elif(ide==5):
        print('Print_cosmology: Dark energy: IDE II')
        print('Print_cosmology: a1:', a1)
        print('Print_cosmology: a2:', a2)
        print('Print_cosmology: nw:', nw)
    else:
        print('Print_cosmology: Error, ide not specified correctly')
        exit()
    print('Print_cosmology: Done')
    print()

## ##
    
## Definitions of cosmological functions ##

# Hubble function: \dot(a)/a
def H(a):
    H2=(Om_r*a**-4)+(Om_m*a**-3)+Om_w*X_de(a)+Om_v+(1.-Om)*a**-2
    return np.sqrt(H2)

# Acceleration function: \ddot(a)/a
def AH(a):
    AH=-0.5*((2.*Om_r*a**-4)+(Om_m*a**-3)+Om_w*(1.+3*w_de(a)*X_de(a))-2.*Om_v)
    return AH

# Dark energy density as a function of 'a'
def X_de(a):
    if(ide==0):
        return np.full_like(a, 1.) # Make numpy compatible
    if(ide==1):
        return a**(-3.*(1.+w0+wa))*np.exp(-3.*wa*(1.-a))
    elif(ide==2):
        return a**(-3.*(1.+w0))
    elif(ide==5):
        f1=(a/a1)**nw+1.
        f2=(1./a1)**nw+1.
        f3=(1./a2)**nw+1.
        f4=(a/a2)**nw+1.
        return ((f1/f2)*(f3/f4))**(-6./nw)

# Dark energy equation-of-state parameter
def w_de(a):
    if(ide==0):
        return np.full_like(a, -1.) # Make numpy compatible
    elif(ide==1):
        return w0+(1.-a)*wa
    elif(ide==2):
        return np.full_like(a, w0) # Make numpy compatible
    elif(ide==5):
        f1=(a/a1)**nw-1.
        f2=(a/a1)**nw+1.
        f3=(a/a2)**nw-1.
        f4=(a/a2)**nw+1.
        return -1.+f1/f2-f3/f4
    
# Omega_r as a function of 'a'
def Omega_r(a):
    return Om_r*(a**-4)/H(a)**2
    
# Omega_m as a function of 'a'
def Omega_m(a):
    return Om_m*(a**-3)/H(a)**2

# Omega_w as a function of 'a'
def Omega_w(a):
    return Om_w*X_de(a)/H(a)**2

# Omega_v as a function of 'a'
def Omega_v(a):
    return Om_v/H(a)**2

# Total Omega as a function of 'a'
def Omega(a):
    return Omega_r(a)+Omega_m(a)+Omega_w(a)+Omega_v(a)

# Plot Omega_i(a)
def plot_omegas():
    
    plt.figure(1,figsize=(20, 6))

    # Omegas - Linear
    plt.subplot(122)
    plt.plot(a_logtab,Omega_r(a_logtab),label=r'$\Omega_r(a)$')
    plt.plot(a_logtab,Omega_m(a_logtab),label=r'$\Omega_m(a)$')
    plt.plot(a_logtab,Omega_w(a_logtab),label=r'$\Omega_w(a)$')
    plt.plot(a_logtab,Omega_v(a_logtab),label=r'$\Omega_v(a)$')
    plt.xlabel(r'$a$')
    plt.ylabel(r'$\Omega_i(a)$')
    plt.legend()

    # Omegas - Log
    plt.subplot(121)
    plt.semilogx(a_logtab,Omega_r(a_logtab),label=r'$\Omega_r(a)$')
    plt.semilogx(a_logtab,Omega_m(a_logtab),label=r'$\Omega_m(a)$')
    plt.semilogx(a_logtab,Omega_w(a_logtab),label=r'$\Omega_w(a)$')
    plt.semilogx(a_logtab,Omega_v(a_logtab),label=r'$\Omega_v(a)$')
    plt.xlabel(r'$a$')
    plt.ylabel(r'$\Omega_i(a)$')
    plt.legend()

    plt.figure(2,figsize=(20, 6))

    # w(a) - Linear
    plt.subplot(121)
    plt.axhline(0,c='k',ls=':')
    plt.axhline(1,c='k',ls=':')
    plt.axhline(-1,c='k',ls=':')
    plt.plot(a_lintab,w_de(a_lintab))
    plt.xlabel(r'$a$')
    plt.ylabel(r'$w(a)$')
    plt.ylim((-1.05,1.05))

    # w(a) - Log
    plt.subplot(122)
    plt.axhline(0,c='k',ls=':')
    plt.axhline(1,c='k',ls=':')
    plt.axhline(-1,c='k',ls=':')
    plt.semilogx(a_logtab,w_de(a_logtab))
    plt.xlabel(r'$a$')
    plt.ylabel(r'$w(a)$')
    plt.ylim((-1.05,1.05))
    
    plt.show()

# Distance and age integrals
#plot_distances = True

def initialise_distances():

    global r_tab, t_tab, rp_tab
    global r, t, rp
    global r0, t0

    # A small number
    small=0.

    ###

    # Integrand for the r(a) calculations and vectorise
    def r_integrand(a):
        return 1./(H(a)*a**2)
    r_integrand_vec=np.vectorize(r_integrand)

    # Function to integrate to get rp(a) (PARTICLE HORIZON) and vectorise
    def rp_integrate(a):
        rp,_=integrate.quad(r_integrand_vec, 0., a)
        return rp
    rp_integrate_vec=np.vectorize(rp_integrate)

    # Function to integrate to get r(a) and vectorise
    def r_integrate(a):
        r,_=integrate.quad(r_integrand_vec, a, 1.)
        return r
    r_integrate_vec=np.vectorize(r_integrate)

    ###

    # Integrand for the t(a) calculation and vectorise
    def t_integrand(a):
        return 1./(H(a)*a)
    t_integrand_vec=np.vectorize(t_integrand)

    # Function to integrate to get r(a) and vectorise
    def t_integrate(a):
        t,_=integrate.quad(t_integrand_vec, 0., a)
        return t
    t_integrate_vec=np.vectorize(t_integrate)

    ###

    # Call the vectorised integration routine
    r_tab=r_integrate_vec(a_tab)
    rp_tab=rp_integrate_vec(a_tab)
    t_tab=t_integrate_vec(a_tab)

    # Add in values r(a=0) and t(a=0) if the 'nonzero' table has been used
    #r_tab=np.insert(r_tab,0,0.)
    #t_tab=np.insert(t_tab,0,0.)

    # Interpolaton function for r(a)
    r_func=interp1d(a_tab,r_tab,kind='cubic',fill_value='extrapolate') #This needs to be linear, not log
    def r_vectorize(a):
        if(a<=small):
            return rp_integrate(1.)
        elif(a>1.):
            print('Error, r(a>1) called:', a)
            return None
        else:
            return r_func(a)
    r=np.vectorize(r_vectorize)

    # Interpolaton function for rp(a)
    rp_func=interpolation.log_interp1d(a_tab,rp_tab,kind='cubic',fill_value='extrapolate')
    def rp_vectorize(a):
        if(a<=small):
            return 0.
        elif(a>1.):
            print('Error, rp(a>1) called:', a)
            return None
        else:
            return rp_func(a)
    rp=np.vectorize(rp_vectorize)

    # Interpolaton function for t(a)
    t_func=interpolation.log_interp1d(a_tab,t_tab,kind='cubic',fill_value='extrapolate')
    def t_vectorize(a):
        if(a<=small):
            return 0.
        elif(a>1.):
            print('Error, t(a>1) called:', a)
            return None
        else:        
            return t_func(a)
    t=np.vectorize(t_vectorize)

    r0=rp(1.)
    t0=t(1.)
    print('Initialise_distances: Horizon size [dimensionless]:', r0)
    print('Initialise_distances: Horizon size [Mpc/h]:', conH0*r0)
    print('Initialise_distances: Universe age [dimensionless]:', t0)
    print('Initialise_distances: Universe age [Gyr/h]:', invH0*t0)
    print('Initialise_distances: r(0):',  r(0.))
    print('Initialise_distances: rp(0):', rp(0.))
    print('Initialise_distances: t(0):',  t(0.))
    print()

# Plot cosmic distances and times (dimensionless)
def plot_distances():
    
    plt.figure(3,figsize=(20, 6))

    # Plot cosmic distance (dimensionless) vs. a on linear scale
    plt.subplot(121)
    plt.plot(a_tab,r_tab,'go',label=r'$r(a)$')
    plt.plot(a_lintab,r(a_lintab),'g-',label='interpolation')
    plt.plot(a_tab,rp_tab,'bo',label=r'$r_p(a)$')
    plt.plot(a_lintab,rp(a_lintab),'b-',label='interpolation')
    plt.plot(a_tab,t_tab,'ro',label=r'$t(a)$')
    plt.plot(a_lintab,t(a_lintab),'r-',label='interpolation')
    plt.legend()
    plt.xlabel(r'$a$')
    plt.ylabel(r'$r(a)$ or $t(a)$')

    # Plot cosmic distance (dimensionless) vs. a on log scale
    plt.subplot(122)
    plt.loglog(a_tab,r_tab,'go',label=r'$r(a)$')
    plt.loglog(a_logtab,r(a_logtab),'g-',label='interpolation')
    plt.loglog(a_tab,rp_tab,'bo',label=r'$r_p(a)$')
    plt.loglog(a_logtab,rp(a_logtab),'b-',label='interpolation')
    plt.loglog(a_tab,t_tab,'ro',label=r'$t(a)$')
    plt.loglog(a_logtab,t(a_logtab),'r-',label='interpolation')
    plt.legend()
    plt.xlabel(r'$a$')
    plt.ylabel(r'$r(a)$ or $t(a)$')
    plt.show()

def initialise_growth():

    global g_tab, f_tab
    global g, f

    print('Initialise_growth: Solving growth equations')

    # Calculate things associated with the linear growth

    # Set initial conditions for the ODE integration
    #d_init = a_lintab_nozero[0]
    amin=a_tab[0]
    d_init = amin
    v_init = 1.

    # Function to calculate delta'
    def dd(v,d,a):
        dd=1.5*Omega_m(a)*d/a**2-(2.+AH(a)/H(a)**2)*v/a
        return dd

    # Function to get delta' and v' in the correct format for odeint
    # Note that it returns [v',delta'] in the 'wrong' order
    #X=[delta,v]
    def dv(X,t):
        return [X[1],dd(X[1],X[0],t)]   

    # Use odeint to get g(a) and f(a) = d ln(g)/d ln(a)
    #gv=odeint(dv,[d_init,v_init],a_lintab_nozero)
    gv=odeint(dv,[d_init,v_init],a_tab)
    g_tab=gv[:,0]
    f_tab=a_tab*gv[:,1]/gv[:,0]

    print('Initialise_growth: ODE solved')

    # Add in the values g(a=0) and f(a=0) if using the 'nonzero' tab
    #g_tab=np.insert(g_tab,0,0.)
    #f_tab=np.insert(f_tab,0,1.)

    print('Initialise_growth: Creating interpolators')

    # Create interpolation function for g(a)
    g_func=interpolation.log_interp1d(a_tab,g_tab,kind='cubic',fill_value='extrapolate')
    def g_vectorize(a):
        if(a<amin):
            return a
        elif(a>1.):
            print('Error, g(a>1) called:', a)
            return None
        else:
            return g_func(a)
    g=np.vectorize(g_vectorize)

    # Create interpolation function for f(a) = dln(g)/dln(a)
    f_func=interpolation.log_interp1d(a_tab,f_tab,kind='cubic',fill_value='extrapolate')
    def f_vectorize(a):
        if(a<amin):
            return 1.
        elif(a>1.):
            print('Error, f(a>1) called:', a)
            return None
        else:       
            return f_func(a)
    f=np.vectorize(f_vectorize)

    print('Initialise_growth: Interpolators done')

    # Check values
    print('Initialise_growth: g(0):', g(0.))
    print('Initialise_growth: g(amin):', g(amin))
    print('Initialise_growth: g(1):', g(1.))
    print('Initialise_growth: f(0):', f(0.))
    print('Initialise_growth: f(amin):', f(amin))
    print('Initialise_growth: f(1):', f(1.))
    print()

# Plot g(a) and f(a)
def plot_growth():
    
    plt.figure(1,figsize=(20, 6))

    # Linear scale
    plt.subplot(121)
    plt.plot(a_tab,g_tab,'bo',label='g(a)')
    plt.plot(a_tab,f_tab,'ro',label='f(a)')
    plt.plot(a_lintab,g(a_lintab),'b-',label=r'interpolation')
    plt.plot(a_lintab,f(a_lintab),'r-',label=r'interpolation')
    plt.legend()
    plt.xlabel(r'$a$')
    plt.xlim((0,1.0))
    plt.ylabel(r'$g(a)$ or $f(a)$')
    plt.ylim((0,1.05))

    # Log scale
    plt.subplot(122)
    plt.semilogx(a_tab,g_tab,'bo',label=r'$g(a)$')
    plt.semilogx(a_tab,f_tab,'ro',label=r'$f(a)$')
    plt.semilogx(a_logtab,g(a_logtab),'b-',label=r'interpolation')
    plt.semilogx(a_logtab,f(a_logtab),'r-',label=r'interpolation')
    plt.legend()
    plt.xlabel(r'$a$')
    plt.ylabel(r'$g(a)$ or $f(a)$')
    plt.ylim((0,1.05))

    # Show the plot
    plt.show()

# Get P(k) from CAMB file
def read_CAMB(fname):
    kPk=np.loadtxt(fname)
    k=kPk[:,0]
    Pk=kPk[:,1]
    return k, Pk

def create_Pk(k_tab,Pk_tab):

    global Pk

    # Create P(k) interpolation function
    Pk_func=interpolation.log_interp1d(k_tab,Pk_tab,kind='cubic')
    def Pk_vectorize(k):
        if(k<k_tab[0]):
            a=np.log(Pk_tab[1]/Pk_tab[0])/np.log(k_tab[1]/k_tab[0])
            b=np.log(Pk_tab[0])-a*np.log(k_tab[0])
            return np.exp(a*np.log(k)+b)
        elif(k>k_tab[-1]):
            a=np.log(Pk_tab[-2]/Pk_tab[-1])/np.log(k_tab[-2]/k_tab[-1])
            b=np.log(Pk_tab[-1])-a*np.log(k_tab[-1])
            return np.exp(a*np.log(k)+b)
        else:
            return Pk_func(k)
    Pk=np.vectorize(Pk_vectorize)

    return Pk
