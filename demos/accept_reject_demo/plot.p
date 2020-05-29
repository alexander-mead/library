reset

x1 = 0. # where binning starts
x2 = 1. # where binning ends
n = 100 # the number of bins
width = (x2-x1)/n # binwidth; evaluates to 1.0
bin(x) = width*(floor((x-x1)/width)+0.5)+x1

# Number of points generated for accept-reject
np = 100000

unset key

#set xrange [Min:Max+0.01]
set xlabel 'x'

set ylabel 'N'

plot 'results.dat' using (bin($1)):(1.0) smooth freq with boxes,\
   3.*(np/n)*x**2 w l lw 2 ti 'Expectation'
