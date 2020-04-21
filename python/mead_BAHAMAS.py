power_dir = '/Users/Mead/Physics/BAHAMAS/power'

# Measured BAHAMAS power spectra file names
def power_file_name(mesh, model, snap, field_pair):
   import os.path as path
   field1 = field_pair[0]
   field2 = field_pair[1]
   file_name1 = power_dir+'/M'+str(mesh)+'/'+model+'_L400N1024_WMAP9_snap'+str(snap)+'_'+field1+'_'+field2+'_power.dat'
   file_name2 = power_dir+'/M'+str(mesh)+'/'+model+'_L400N1024_WMAP9_snap'+str(snap)+'_'+field2+'_'+field1+'_power.dat'
   if path.isfile(file_name1):
      return file_name1
   else:
      return file_name2

# Measured BAHAMAS errors between different realisations of the AGN_TUNED_nu0 model
def error_file_name(mesh, snap, field_pair):
   field1 = field_pair[0]
   field2 = field_pair[1]
   return power_dir+'/M'+str(mesh)+'/L400N1024_WMAP9_snap'+str(snap)+'_'+field1+'_'+field2+'_error.dat'

# Corresponding HMcode power file
def HMcode_file_name(mesh, snap):
   return power_dir+'/M'+str(mesh)+'/HMcode/snap'+str(snap)+'_HMcode.dat'

# Get the snapshot number corresponding to different BAHAMAS redshifts
def z_to_snap(z):
   if(z == 0.000):
      snap = 32
   elif(z == 0.125):
      snap = 31
   elif(z == 0.250):
      snap = 30
   elif(z == 0.375):
      snap = 29
   elif(z == 0.500):
      snap = 28
   elif(z == 0.750):
      snap = 27
   elif(z == 1.000):
      snap = 26
   elif(z == 1.250):
      snap = 25
   elif(z == 1.500):
      snap = 24
   elif(z == 1.750):
      snap = 23
   elif(z == 2.000):
      snap = 22
   else:
      print('Redshift = ', z)
      raise ValueError('Snapshot not stored corresponding to this z')
   return snap

# Read a BAHAMAS power/error file and output k, power, error
def get_measured_power(mesh, model, z, field_pair, realisation_errors=False, correct_shot_noise=False):

   from numpy import loadtxt

   column_k = 0
   column_power = 1
   column_shot = 2
   column_modes = 3
   column_error = 4

   snap = z_to_snap(z)
   infile = power_file_name(mesh, model, snap, field_pair)
   data = loadtxt(infile)
   k = data[:, column_k]
   power = data[:, column_power]
   modes = data[:, column_modes]
   shot = data[:, column_shot]

   if(realisation_errors):
      infile = error_file_name(mesh, snap, field_pair)
      data = loadtxt(infile)
   error = data[:, column_error]

   # Subtract shot noise
   if (correct_shot_noise):
      power = power - shot

   return k, power, shot, modes, error

def get_measured_response(mesh, model, z, field_pair):

   # Get the power from the hydro model
   k, power, _, _, _ = get_measured_power(mesh, model, z, field_pair,
      realisation_errors=False,
      correct_shot_noise=True
      )

   # Get the power from the DMONLY model
   dmonly_name = 'DMONLY_2fluid_nu0'
   _, dmonly, _, _, _ = get_measured_power(mesh, dmonly_name, z, field_pair=('all', 'all'),
      realisation_errors=False,
      correct_shot_noise=True
      )

   # Make the response
   response = power/dmonly

   return k, response

def get_HMcode_power(mesh, z):
   import numpy as np
   snap = z_to_snap(z)
   infile = HMcode_file_name(mesh, snap)
   data = np.loadtxt(infile)
   k = data[:, 0]
   power = data[:, 1]
   return k, power

def get_HMcode_corrected_power(mesh, model, z, field_pair, realisation_errors=False, correct_shot_noise=False):

   _, _, shot, modes, error = get_measured_power(mesh, model, z, field_pair, realisation_errors, correct_shot_noise)
   _, response = get_measured_response(mesh, model, z, field_pair)
   k, power = get_HMcode_power(mesh, z)

   power = power*response

   return k, power, shot, modes, error

