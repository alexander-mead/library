# Import statements
import numpy as np
from scipy.integrate import quad

# Calculate the normalisation via integration
def normalisation(f,x1,x2,*args):
    norm,_ = quad(lambda x: f(x,*args), x1, x2)
    return norm

# Cumulative distribution function via integration
def cumulative(x,f,x1,*args):
    C,_ = quad(lambda x: f(x,*args), x1, x)
    return C
cumulative = np.vectorize(cumulative)

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

# Draw random numbers from a 1D distribution
# n - number of draws to make from f
# f - f(x,y) to draw from
# x1,x2 - limits on x axis
# nx - number of points along x axis
def draw_from_1D(n,f,x1,x2,nx,*args):
    x = np.linspace(x1,x2,nx)    
    C = cumulative(x,f,x1,*args)
    xi = np.zeros(n)
    for i in range(n):
        xi[i] = draw_from_distribution(x,C)
    return xi

# Draw random numbers from a 2D distribution
# n - number of draws to make from f
# f - f(x,y) to draw from
# x1,x2 - limits on x axis
# y1,y2 - limits on y axis
# nx - number of points along x axis
# ny - number of points along y axis
def draw_from_2D(n,f,x1,x2,nx,y1,y2,ny):

    # Pixel sizes in x and y
    dx = (x2-x1)/np.real(nx)
    dy = (y2-y1)/np.real(ny)

    # Linearly spaced arrays of values corresponding to pixel centres
    x = np.linspace(x1+dx/2.,x2-dx/2.,nx)
    y = np.linspace(y1+dy/2.,y2-dy/2.,ny)

    # Make a grid of xy coordinate pairs
    xy = np.array(np.meshgrid(x,y))
    xy = xy.reshape(2,nx*ny) # Reshape the grid (2 here coresponds to 2 coordinates: x,y)
    xy = np.transpose(xy).tolist() # Convert to a long list
    
    # Make array of function values corresponding to the xy coordinates
    X,Y = np.meshgrid(x,y)
    z = f(X,Y)      # Array of function values
    z = z.flatten() # Flatten array to create a long list of function values
    z = z/sum(z)    # Force normalisation

    # Make a list of integers linking xy coordiantes to function values
    i = list(range(z.size)) 

    # Make the random choices with probabilties proportional to the function value
    # The integer chosen by this can then be matched to the xy coordiantes
    j = np.random.choice(i,n,replace=True,p=z) 
    
    # Now match the integers to the xy coordinates
    xs = []
    ys = []
    for i in range(n):
        xi,yi = xy[j[i]]
        xs.append(xi)
        ys.append(yi)

    # Random numbers for inter-pixel displacement
    dxs = np.random.uniform(-dx/2.,dx/2.,n) 
    dys = np.random.uniform(-dy/2.,dy/2.,n)

    # Apply uniform-random displacement within a pixel
    xs = xs+dxs
    ys = ys+dys
        
    return xs,ys
