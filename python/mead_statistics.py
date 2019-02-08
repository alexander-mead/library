# Import statements
import numpy as np
from scipy.integrate import quad

# Calculate the normalisation via integration
def normalisation(f,x1,x2,*args):
    norm,_ = quad(lambda x: f(x,*args), x1, x2)
    return norm

# Cumulative distribution function via integration
def cumulative_simple(x,f,x1,*args):
    C,_ = quad(lambda x: f(x,*args), x1, x)
    return C
cumulative = np.vectorize(cumulative_simple)

# Draw a random number given a cumulative probability
def draw_from_distribution(x,Cx):       
    r = np.random.uniform(Cx[0],Cx[-1])
    xi = np.interp(r, Cx, x)    
    return xi

# Calculate the nth moment of distribution f
def moment(n,f,x1,x2,*args):
    norm = normalisation(f,x1,x2,*args)
    m,_ = quad(lambda x: (x**n)*f(x,*args)/norm, x1, x2)
    return m

# Calculate the variance of distribution f
def variance(f,x1,x2,*args):
    m1 = moment(1,f,x1,x2,*args)
    m2 = moment(2,f,x1,x2,*args)
    return m2-m1**2
