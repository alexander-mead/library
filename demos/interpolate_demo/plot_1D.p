unset multiplot
reset

set term qt dashed

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

set ylabel 'y'

plot \
'table.dat' u 1:2 w p ps 3 lc rgb 'black' ti 'Data',\
'results.dat' u 1:2 w l lw 3 lc rgb 'black' ti 'Truth',\
'results.dat' u 1:2 w l lw 3 lc rgb 'red' ti 'Linear',\
'results.dat' u 1:3 w l lw 3 lc rgb 'green' ti 'Quadratic',\
'results.dat' u 1:4 w l lw 3 lc rgb 'blue' ti 'Cubic'

#pause -1

set xlabel 'x'
set format x

set ylabel 'Error'
set yrange [0.97:1.03]

plot 1 w l ls -1 noti,\
     1.01 w l dt 2 lc 'black' noti,\
     0.99 w l dt 2 lc 'black' noti,\
     'ratio.dat' u 1:2 w l lw 3 lc rgb 'red' noti 'Linear',\
     'ratio.dat' u 1:3 w l lw 3 lc rgb 'green' noti 'Quadratic',\
     'ratio.dat' u 1:4 w l lw 3 lc rgb 'blue' noti 'Cubic'

unset multiplot
