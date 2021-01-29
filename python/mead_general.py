# https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
def file_length(fname):

    with open(fname) as f:
        for i, _ in enumerate(f):
            pass
    return i + 1

# Returns values of the array at the list of array position integers
def array_values_at_indices(array, list_of_array_positions):

   if (len(array.shape) == 1):
      return array[list_of_array_positions]
   elif (len(array.shape) == 2):
      ix, iy = zip(*list_of_array_positions)
      result = array[ix, iy]
      return result
   else:
      ValueError('Error, this only works in either one or two dimensions at the moment')
      return None

# Return a logarithmically spaced range of numbers
def logspace(xmin, xmax, nx):

   from numpy import logspace, log10
   return logspace(log10(xmin), log10(xmax), nx)

# Mutliply all elements in a list by a constant
def multiply_list_elements(multiple, list):

   return [multiple*x for x in list]

# Sequential colors
def seq_color(i, n, cmap):

    return cmap(i/(n-1))

def colour(i):

    color = 'C%d' % i
    return color



    
