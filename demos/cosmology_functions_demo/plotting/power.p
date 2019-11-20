reset

if(!exists('print')){print=0}
#if(print==0){set term aqua}
if(print==0){set term qt}
if(print==1){set term post enh col; set output 'power.eps'}

power='data/power.dat'

set log x
set xlabel 'k / h Mpc^{-1}'
set format x '10^{%T}'

pmin=1e-10
pmax=1e3
set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'
set yrange [pmin:pmax]

set key top left

plot power u 1:2 w l lw 2 ti 'z = 0',\
     power u 1:3 w l lw 2 ti 'z = 1'
