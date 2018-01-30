reset

set xlabel 'x'

set ylabel 'P(x)'
unset ytics
set format y ''

plot 'rng.dat' u 1:2 w l lw 3 ti 'Generator',\
     'rng.dat' u 1:3 w l lw 3 ti 'Function'
