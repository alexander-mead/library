reset

set key top left
set xrange [0:5]
set yrange [0:5]

plot 'points.dat' u 1:2 w p ps 2 pt 7 lc rgb 'black' ti 'Data',\
     'cubic.dat' u 1:2 w l lw 5 lc rgb 'red' ti 'Cubic',\
     'cubic.dat' u 1:3 w l lw 2 lc rgb 'orange' ti 'Lagrange'
