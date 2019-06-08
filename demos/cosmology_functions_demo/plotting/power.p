reset

if(!exists('print')){print=0}
if(print==0){set term aqua}
if(print==1){set term post enh col; set output 'power.eps'}

set log x
set xlabel 'k / h Mpc^{-1}'
set format x '10^{%T}'

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

set key top left

plot 'power.dat' u 1:2 w l lw 2 ti 'z = 0',\
     'power.dat' u 1:3 w l lw 2 ti 'z = 1'
