import numpy as np
from scipy import integrate

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

# A very simple rectangular integration that assumes equal sized bins
def integrate_rectangular(fx,x):
    dx = x[1]-x[0]
    return sum(fx)*dx