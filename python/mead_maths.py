import numpy as np

# Tophat Fourier transform function
# Normalised such that T(x=0)=1
def Tophat(x):
    xmin=1e-5
    return np.where(np.abs(x)<xmin, 1.-x**2/10., (3./x**3)*(np.sin(x)-x*np.cos(x)))

# Gaussian function
# Normalised such that G(x=0)=1
# mu - mean
# sig - standard deviation
def Gaussian(x, mu, sig):
    return np.exp(-(x-mu)**2/(2.*sig**2))

# Provides the two real solutions for x for a quadratic a*x^2 + b*x + c = 0 
# TODO: Expand for complex solutions
def solve_quadratic(a, b, c):
    des = b**2-4.*a*c
    if (des > 0.):
        root = np.sqrt(b**2-4.*a*c)
    else:
        print('FUCK')
    f1 = (root-b)/(2.*a)
    f2 = (-root-b)/(2.*a)
    return f1, f2