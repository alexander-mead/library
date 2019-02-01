import numpy as np

# Tophat Fourier transform function
# Normalised such that T(x=0)=1
def Tophat(x):
    xmin=1e-5
    return np.where(np.abs(x)<xmin, 1.-x**2/10., 3.*(np.sin(x)-x*np.cos(x))/x**3)

# Gaussian function
# Normalised such that G(x=0)=1
# mu - mean
# sig - standard deviation
def Gaussian(x, mu, sig):
    return np.exp(-(x-mu)**2/(2.*sig**2))
