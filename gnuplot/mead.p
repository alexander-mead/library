## Mead gnuplot library ##

## Commonly used functions ##

# Returns maximum of x, y
max(x, y) = (x > y ? x : y)

# Linear spacing between xmin and xmax
progression(xmin, xmax, i, n) = xmin+(xmax-xmin)*(real(i)-1)/(real(n)-1)

# Colour scheme thing
#color_scheme(c1, c2, c3) = sprintf('%d, %d, %d', c1, c2, c3) # Need to print out integers rather than string
#color_scheme(c) = sprintf('%d, %d, %d', c[1], c[2], c[3]) # Does not seem to work with array argument

## Cosmology ##

khlab = 'k / h Mpc^{-1}'
plab = '{/Symbol D}^2(k)'

## Colour schemes from http://gnuplot.sourceforge.net/demo/pm3dcolors.html ##

# Rainbow
array rainbow[3]
rainbow[1] = 33
rainbow[2] = 13
rainbow[3] = 10

# AFM hot
array AFM_hot[3]
AFM_hot[1] = 34
AFM_hot[2] = 35
AFM_hot[3] = 36

# green-red-violet
array grv[3]
grv[1] = 3 
grv[2] = 11
grv[3] = 6

# ocean (try other permutations too)
array ocean[3]
ocean[1] = 23
ocean[2] = 28
ocean[3] = 3

# hot
array hot[3]
hot[1] = 21
hot[2] = 22
hot[3] = 23

# color printable on black & whitee
array cpg[3]
cpg[1] = 30
cpg[2] = 31
cpg[3] = 32