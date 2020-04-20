# Read the halo-model power spectra and output the wavenumber and power spectrum
def get_power(infile):
   from numpy import loadtxt
   column_k = 0
   column_power = 4
   data = loadtxt(infile)
   k = data[:, column_k]
   power = data[:, column_power]
   return k, power

# File names for the halo-model power spectra that are fitted to the BAHAMAS data
def fitted_power_file_name(mode, model, z, chain, field_pair, label):
   _, f1 = field_name_to_letter_and_integer(field_pair[0])
   _, f2 = field_name_to_letter_and_integer(field_pair[1])
   dir = '/Users/Mead/Physics/HMx/fitting'
   file = dir+'/m'+str(mode)+'/'+model+'/c%d_z%1.3f_cos1_%s_power_%d%d.dat' % (chain, z, label, f1, f2)
   return file

# Read the halo-model power spectra and output the wavenumber and power spectrum
def get_fitted_power(mode, model, z, chain, field_pair, label):
   #from numpy import loadtxt
   infile = fitted_power_file_name(mode, model, z, chain, field_pair, label)
   k, power = get_power(infile)
   return k, power

# Read the halo-model power spectra and output the wavenumber and power spectrum
def get_fitted_response(mode, model, z, chain, field_pair, label):
   #from numpy import loadtxt
   k, power = get_fitted_power(mode, model, z, chain, field_pair, label)
   _, dmonly = get_fitted_power(mode, model, z, chain, ('dmonly','dmonly'), label)
   response = power/dmonly
   return k, response

def get_HMcode_power(mode, model, z, chain):

   from numpy import loadtxt

   dir = '/Users/Mead/Physics/HMx/fitting'
   infile = dir+'/m'+str(mode)+'/'+model+'/c1_z%1.3f_cos1_HMcode_power.dat' % (z)

   data = loadtxt(infile)

   column_k = 0
   column_power = 1
   
   k = data[:, column_k]
   power = data[:, column_power]

   return k, power

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
def parameter_file_name(fit, model, chain, z=None):
   from os.path import isfile
   directory = '/Users/Mead/Physics/HMx/'
   if (z == None):
      file = directory+'fitting/m%d/%s/c%d_params.dat' % (fit, model, chain)
   else:
      file = directory+'fitting/m%d/%s/c%d_z%1.3f_params.dat' % (fit, model, chain, z)
   if not isfile(file):
      print(file)
      raise NameError('File does not exist')
   return file

# Read best-fitting parameters from file
def read_parameter_file(file):

   # Import statements
   from sys import path 
   from numpy import genfromtxt, float

   path.append('/Users/Mead/Physics/library/python')
   from mead_strings import read_first_number_from_line_in_file

   # Read the file
   parameter_file_figure_of_merit_line = 2
   parameter_file_header_length = 6
   parameter_file_footer_length = 1
   parameter_file_names_column = 1
   parameter_file_best_fit_column = 3

   # Import the table from the text file
   data = genfromtxt(file, dtype = 'str', 
                     skip_header = parameter_file_header_length, 
                     skip_footer = parameter_file_footer_length
                     )
   
   # Make arrays of parameter names and best-fitting values
   parameter_names = data[:, parameter_file_names_column]
   best_fitting_values = data[:, parameter_file_best_fit_column].astype(float)

   # Get the figure of merit
   figure_of_merit = read_first_number_from_line_in_file(file, parameter_file_figure_of_merit_line)

   # Fin 
   return parameter_names, best_fitting_values, figure_of_merit

# Read the best-fitting parameters for a range of models and redshifts
def read_parameter_file_models_zs(fit, models, chain, zs):

   from numpy import asarray 

   # Initialise empty lists for output
   best_fitting_values_models = []
   figures_of_merit = []

   for model in models:

      # Initialise empty lists 
      best_fitting_values_z = []
      figures_of_merit_z = []

      for z in zs:

         # Get the input file name
         infile = parameter_file_name(fit, model, chain, z)

         # Read the data
         parameter_names, parameter_values, figure_of_merit = read_parameter_file(infile)

         # Append information to lists
         best_fitting_values_z.append(parameter_values)
         figures_of_merit_z.append(figure_of_merit)

      # Append lists to lists
      best_fitting_values_models.append(best_fitting_values_z)
      figures_of_merit.append(figures_of_merit_z)

   # Convert lists to arrays
   best_fitting_values_models = asarray(best_fitting_values_models)
   figures_of_merit = asarray(figures_of_merit)

   # Return with the good stuff
   return parameter_names, best_fitting_values_models, figures_of_merit

   # Read the best-fitting parameters for a range of models and redshifts
def read_parameter_file_models(fit, models, chain):

   from numpy import asarray 

   # Initialise empty lists for output
   best_fitting_values_models = []
   figures_of_merit = []

   for model in models:

      # Get the input file name
      infile = parameter_file_name(fit, model, chain)

      # Read the data
      parameter_names, parameter_values, figure_of_merit = read_parameter_file(infile)

      # Append information to lists
      best_fitting_values_models.append(parameter_values)
      figures_of_merit.append(figure_of_merit)

      # Append lists to lists
      #best_fitting_values_models.append(best_fitting_values)
      #figures_of_merit.append(figures_of_merit)

   # Convert lists to arrays
   best_fitting_values_models = asarray(best_fitting_values_models)
   figures_of_merit = asarray(figures_of_merit)

   # Return with the good stuff
   return parameter_names, best_fitting_values_models, figures_of_merit