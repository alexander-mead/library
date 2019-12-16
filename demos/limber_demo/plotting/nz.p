reset

if(print==0) set term qt dashed
if(print==1) set term post enh col; set output 'nz.eps'

nz1='data/nz1.dat'
nz2='data/nz2.dat'

set xlabel 'z'
set xrange [0:*]

set ylabel 'n(z)'
set yrange [0:*]

plot nz1 u 1:2 w l lw 3 ti 'Field 1',\
     nz2 u 1:2 w l lw 3 ti 'Field 2'
