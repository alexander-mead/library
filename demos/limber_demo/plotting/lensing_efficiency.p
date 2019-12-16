unset multiplot
reset

if(!exists('print')) {print = 0}
if(print==0) set term qt dashed
if(print==1) set term post enh col; set output 'lensing_efficiency.eps'

file1='data/efficiency1.dat'
file2='data/efficiency2.dat'

set ylabel 'q(r)'
set yrange [0:1]

set multiplot layout 2,1

set xlabel 'r / (Mpc/h)'

plot file1 u 1:2 w l lw 3 ti 'Tracer 1',\
     file2 u 1:2 w l lw 3 ti 'Tracer 2'

set xlabel 'z'

plot file1 u 1:3 w l lw 3 noti,\
     file2 u 1:3 w l lw 3 noti

unset multiplot
