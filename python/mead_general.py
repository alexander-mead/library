import numpy as np
from scipy import interpolate
from scipy import integrate

# From https://stackoverflow.com/questions/29346292/logarithmic-interpolation-in-python
#def log_interp1d(xx, yy, kind='cubic'):
#    logx = np.log10(xx)
#    logy = np.log10(yy)
#    lin_interp = interpolate.interp1d(logx, logy, kind=kind)
#    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
#    return log_interp

# My version of: https://stackoverflow.com/questions/29346292/logarithmic-interpolation-in-python
# log10 replaced by log and other options for interp1d added, variables renamed
def log_interp1d(x,y,\
                 kind='linear',\
                 axis=-1,\
                 copy=True,\
                 bounds_error=None,\
                 fill_value=np.nan,\
                 assume_sorted=False):
    logx = np.log(x)
    logy = np.log(y)
    lin_interp = interpolate.interp1d(logx,logy,\
                                      kind=kind,\
                                      axis=axis,\
                                      copy=copy,\
                                      bounds_error=bounds_error,\
                                      fill_value=fill_value,\
                                      assume_sorted=assume_sorted)
    log_interp = lambda z: np.exp(lin_interp(np.log(z)))
    return log_interp

## Logarithmically space a set of points
## What about np.logspace??
#def logspace(xmin,xmax,nx,\
#             endpoint=True,\
#             retstep=False,\
#             dtype=None):
#    return np.exp(np.linspace(np.log(xmin),np.log(xmax),nx,\
#                              endpoint=endpoint,\
#                              retstep=retstep,\
#                              dtype=dtype))


# A routine to integrate in log space. 
# This may actually be pretty useless... not sure. Should do speed tests
def integrate_quad_log(func,a,b,\
                       args=(),\
                       full_output=0,\
                       epsabs=1.49e-08,\
                       epsrel=1.49e-08,\
                       limit=50,\
                       points=None,\
                       weight=None,\
                       wvar=None,\
                       wopts=None,\
                       maxp1=50,\
                       limlst=50):
    loga=np.log(a)
    logb=np.log(b)
    ans=integrate.quad(lambda x, *args: np.exp(x)*func(np.exp(x),*args), loga, logb,\
                       args=args,\
                       full_output=full_output,\
                       epsabs=epsabs,\
                       epsrel=epsrel,\
                       limit=limit,\
                       points=points,\
                       weight=weight,\
                       wvar=wvar,\
                       wopts=wopts,\
                       maxp1=maxp1,\
                       limlst=limlst)
    return ans
