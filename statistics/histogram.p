reset

set xlabel 'x'
set xrange [-5:5]

set ylabel 'PDF'
unset ytics
set format y ''

x0=0.
sigma=1.

Gaussian(x,x0,sigma)=exp(-(x-x0)**2/(2.*sigma**2))/(sigma*sqrt(2.*pi))

plot Gaussian(x,x0,sigma) w l lc -1 lw 3 ti 'Distribution',\
     'histogram.dat' u 1:2 w l lw 3 ti 'Histogram'
