reset

#set term x11

plot \
'results_linear.dat' u 1:2 w l lw 3 lc rgb 'black' ti 'Truth',\
'results_linear.dat' u 1:3 w l lw 3 lc rgb 'red' ti 'Linear',\
'results_quadratic.dat' u 1:3 w l lw 3 lc rgb 'green' ti 'Quadratic',\
'results_cubic.dat' u 1:3 w l lw 3 lc rgb 'blue' ti 'Cubic'

pause -1

set yrange [0.8:1.2]

plot \
1 w l ls -1 noti,\
'results_linear.dat' u 1:4 w l lw 3 lc rgb 'red' ti 'Linear',\
'results_quadratic.dat' u 1:4 w l lw 3 lc rgb 'green' ti 'Quadratic',\
'results_cubic.dat' u 1:4 w l lw 3 lc rgb 'blue' ti 'Cubic'
