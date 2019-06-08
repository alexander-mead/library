reset

set term aqua

set xlabel 't'

set ylabel 'x'

plot 'results.dat' u 1:2 w l lw 4 ti 'Analytic Solution',\
     'results.dat' u 1:3 w l lw 2 ti 'ODE solver'
