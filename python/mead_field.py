def make_Gaussian_random_field_2D(mean_value, power_spectrum, map_size, mesh_cells, periodic, *args):

   # mean_value: mean value for the field
   # power_spectrum: P(k, *args) power spectrum for the field [(length units)^2]
   # map_size: side length for the map [length units]
   # mesh_cells: number of mesh cells for the field
   # periodic: should the map be considered to be periodic or not?
   # *args: extra arguments for power spectrum

   import numpy as np
   #from numpy import pi as pi
   #from numpy import fft.fftfreq as fft.fftfreq
   
   # TODO: Enforce Hermitian condition
   # TODO: Use real-to-real FFT
   # TODO: Generalise to nD

   # Parameters
   pad_fraction = 2 # Amount to increase size of array if non-periodic

   if(periodic):
      tran_cells = mesh_cells
      tran_size = map_size
   else:
      tran_cells = pad_fraction*mesh_cells
      tran_size = pad_fraction*map_size
   
   cfield = np.zeros((tran_cells, tran_cells), dtype=complex)
   
   mesh_size = tran_size/tran_cells

   cfield = np.zeros((tran_cells, tran_cells), dtype=complex)
   
   # Look-up arrays for the wavenumbers
   kx = np.fft.fftfreq(tran_cells, d=mesh_size)*2.*np.pi # 2pi needed to convert frequency to angular frequency
   ky = kx # The k in the y direction are identical to those in the x direction

   # Fill the map in Fourier space
   # TODO: Loops in python are bad. Is there a clever way of doing this?
   for i in range(tran_cells):
      for j in range(tran_cells):

            # Absolute wavenumber corresponding to cell in FFT array
            k = np.sqrt(kx[i]**2+ky[j]**2)

            if(i == 0 and j == 0):
               cfield[i, j] = mean_value
               #cfield[i, j] = 0.
            else:
               sigma = np.sqrt(power_spectrum(k, *args))/tran_size
               (x, y) = np.random.normal(0., sigma, size=2)
               cfield[i, j] = complex(x, y)

   # FFT
   cfield = np.fft.ifft2(cfield)
   cfield = cfield*tran_cells**2 # Normalisation
   
   # Convert to real
   field = np.zeros((mesh_cells, mesh_cells))
   field = np.real(cfield[0:mesh_cells, 0:mesh_cells]) # For non-periodic arrays we discard some
   #field = field+mean
   
   return field

# Return a list of the integer coordinates of all local maxima in a 2D field
def find_field_peaks_2D(field, nx, ny, periodic):

   # Initialise an empty list
   peaks = []

   # Loop over all cells in field
   for ix in range(nx):
      for iy in range(ny):

         # Get the central value of tjhe field
         central_value = field[ix, iy]

         # Get the coordinates of the neighbours
         neighbour_cells = neighbour_cells_2D(ix, iy, nx, ny, periodic)
         #print('Neighbour cells:', type(neighbour_cells), neighbour_cells)
         i, j = zip(*neighbour_cells)
         #print(ix, iy, i, j)

         if (all(neighbour_value < central_value for neighbour_value in field[i, j])):
            peaks.append([ix, iy])

   return peaks

#def downward_trajectory(i, j, field, nx, ny, periodic):
#   
#   neighbour_cells = neighbour_cells_2D(i, j, nx, ny, periodic)
#
#   i, j = zip(*neighbour_cells)
#   field_values = field(neighbour_cells)

# Return a list of all cells that are neighbouring some integer cell coordinate
def neighbour_cells_2D(ix, iy, nx, ny, periodic):

   # Check if the cell is actually within the array
   if ((ix < 0) or (ix >= nx) or (iy < 0) or (iy >= ny)):
      raise ValueError('Cell coordinates are wrong')

   # Initialise an empty list
   cells = []

   # Name the cells sensibly
   ix_cell = ix
   ix_mini = ix-1
   ix_maxi = ix+1
   iy_cell = iy
   iy_mini = iy-1
   iy_maxi = iy+1

   # Deal with the edge cases for periodic and non-periodic fields sensibly
   if (ix_mini == -1):
      if (periodic):
         ix_mini = nx-1
      else:
         ix_mini = None
   if (iy_mini == -1):
      if (periodic):
         iy_mini = nx-1
      else:
         iy_mini = None
   if (ix_maxi == nx):
      if (periodic):
         ix_maxi = 0
      else:
         ix_maxi = None
   if (iy_maxi == ny):
      if (periodic):
         iy_maxi = 0
      else:
         iy_maxi = None

   # Add the neighbouring cells to the list
   for ixx in [ix_mini, ix_cell, ix_maxi]:
      for iyy in [iy_mini, iy_cell, iy_maxi]:
         if ((ixx != None) and (iyy != None)):
            cells.append([ixx, iyy])

   # Remove the initial cell from the list
   cells.remove([ix_cell, iy_cell])

   # Return a list of lists
   return cells