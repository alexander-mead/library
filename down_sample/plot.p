reset

unset key

set xlabel 'x'
set xrange [0:1]

set ylabel 'y'
set yrange [0:1]

set zlabel 'z'
set zrange [0:1]

splot \
'cube.dat' u 1:2:3 w p ps 1 pt 7 lc rgb 'black'

pause -1

splot \
'cube.dat' u 1:2:3 w p ps 1 pt 7 lc rgb 'black',\
'accepted.dat' u 1:2:3 w p ps 1 pt 7 lc rgb 'red'

pause -1

splot \
'accepted.dat' u 1:2:3 w p ps 1 pt 7 lc rgb 'red'