reset

set term aqua dashed

amin=0.
amax=1.
set xrange [amin:amax]
set xlabel 'a'

Rmin=0.
Rmax=15000.
set yrange [Rmin:Rmax]
set ylabel 'R / h^{-1} Mpc'

set key top right

plot 'distance.dat' u 1:3 w l lw 3 dt 1 ti 'Comoving angular distance',\
     'distance.dat' u 1:4 w l lw 3 dt 1 ti 'Physical angular distance',\
     'distance.dat' u 1:5 w l lw 3 dt 1 ti 'Luminosity distance',\
     'distance.dat' u 1:2 w l lw 3 dt 2 ti 'Comoving distance'