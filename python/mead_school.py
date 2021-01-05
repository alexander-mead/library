# Multiplication practice
def multiplication_practice(rmin, rmax, n):

   import time
   import random

   # Initial white space
   print()

   t1 = time.time()
   nc = 0
   for i in range(n):

      a = random.randint(rmin, rmax)
      b = random.randint(rmin, rmax)

      print('%d times %d' % (a, b))
      c = float(input())

      if (c == a*b):
         nc = nc+1
      else:
         print('Wrong, the correct answer is %d' % (a*b))
      print()

   T = time.time()-t1
   print('You got %d out of %d correct!' % (nc, n))
   print('You took %1.2f seconds' % (T))
   print('That is %1.2f seconds per multiplication' % (T/n))