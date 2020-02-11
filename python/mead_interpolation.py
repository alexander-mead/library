import numpy as np
from scipy import interpolate

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