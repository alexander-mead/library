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