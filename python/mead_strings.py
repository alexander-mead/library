# Extract a number from a 
# Based on https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python/4289415#4289415
def read_first_number_from_line_in_file(file, line_number):

   # Open the file
   with open(file, 'r') as f:

      # Get the first line
      for _ in range(line_number):
         line_text = f.readline()

      # Now find the number within the string
      # This will find the first number I think
      for parts_of_words in line_text.split():
         try:
            number = float(parts_of_words)
         except ValueError:
            pass

   # Close the file
   f.close()

   # Return with the number
   return number