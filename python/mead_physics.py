import numpy as np
import mead_constants as const

# Black-body irradiance [W m^-2 Hz^-1 Sr^-1]
# nu - emission frequency [Hz]
# T - black-body temperature [K]
def black_body_nu(nu,T):
    h = const.h
    c = const.c
    kB = const.kB
    a = 2.*(h*nu**3)/c**2
    x = (h*nu)/(kB*T)
    #a = (2.*const.h*nu**3)/const.c**2
    #x = const.h*nu/(const.kB*T)
    return a/(np.exp(x)-1.)

# Returns the frequency corresponding to the black-body maximum (Wein's law) [Hz]
# T - black-body temperature [K]
def Wein_law_nu(T):
    a = 5.879e10
    return T*a
