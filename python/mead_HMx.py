# File names for the halo-model power spectra that are fitted to the BAHAMAS data
def fitted_power_file_name(mode, model, z, field_pair):
   _, f1 = field_name_to_letter_and_integer(field_pair[0])
   _, f2 = field_name_to_letter_and_integer(field_pair[1])
   dir = '/Users/Mead/Physics/HMx/fitting'
   file = dir+'/m'+str(mode)+'/'+model+'/z%1.3f_n3000_c1_power_%d%d.dat' % (z, f1, f2)
   return file

# Read the halo-model power spectra and output the wavenumber and power spectrum
def get_fitted_power(mode, model, z, field_pair):
   from numpy import loadtxt
   infile = fitted_power_file_name(mode, model, z, field_pair)
   data = loadtxt(infile)
   k = data[:,0]
   power = data[:,4]
   return k, power

# Read the halo-model power spectra and output the wavenumber and power spectrum
def get_fitted_response(mode, model, z, field_pair):
   from numpy import loadtxt
   k, power = get_fitted_power(mode, model, z, field_pair)
   _, dmonly = get_fitted_power(mode, model, z, field_pair=('dmonly','dmonly'))
   response = power/dmonly
   return k, response

# Output the correpsonding letter and integer for a given field (e.g., matter -> m, 2)
def field_name_to_letter_and_integer(name):
   if (name == 'dmonly'):
      symbol = 'm'
      integer = 1
   elif (name == 'all'):
      symbol = 'm'
      integer = 2
   elif (name == 'dm'):
      symbol = 'c'
      integer = 3
   elif (name =='gas'):
      symbol = 'g'
      integer = 4
   elif (name == 'stars'):
      symbol = '*'
      integer = 5
   elif (name == 'epressure'):
      symbol = 'p'
      integer = 8
   else:
      print('Field name = ', name)
      raise ValueError('Field name not recognised')
   return symbol, integer

    # Parameter file name
def parameter_file_name(fit, model, z):
   from os.path import isfile
   directory = '/Users/Mead/Physics/HMx/'
   file = directory+'fitting/m%d/%s/z%1.3f_n3000_c1_params.dat' % (fit, model, z)
   if not isfile(file):
      print(file)
      raise NameError('File does not exist')
   return file

# Read best-fitting parameters from file
def read_parameter_file(file):
   # Parameters governing *_params.dat files
   from numpy import genfromtxt, float
   parameter_file_header_length = 4
   parameter_file_footer_length = 1
   parameter_file_names_column = 1
   parameter_file_best_fit_column = 3
   data = genfromtxt(file, dtype = 'str', 
                     skip_header = parameter_file_header_length, 
                     skip_footer = parameter_file_footer_length
                     )
   
   parameter_names = data[:, parameter_file_names_column]
   best_fitting_values = data[:, parameter_file_best_fit_column].astype(float)
   return parameter_names, best_fitting_values

# Read the best-fitting parameters for a range of models and redshifts
def read_parameter_file_models_zs(fit, models, zs):
   from numpy import asarray
   best_fitting_values_models = []
   for model in models:
      best_fitting_values_z = []
      for z in zs:
         infile = parameter_file_name(fit, model, z)
         parameter_names, parameter_values = read_parameter_file(infile)
         best_fitting_values_z.append(parameter_values)
      best_fitting_values_models.append(best_fitting_values_z)
   best_fitting_values_models = asarray(best_fitting_values_models)
   return parameter_names, best_fitting_values_models