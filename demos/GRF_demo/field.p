reset

field='data/field_slice.dat'

set size square

set xlabel 'x / h^{-1} Mpc'

set ylabel 'y / h^{-1} Mpc'

set cblabel '{/Symbol d}'

unset key

plot field with image
