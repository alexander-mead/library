reset

if(!exists('print')){print=0}
if(print==0){set term aqua}
if(print==1){set term post enh col; set output 'power.eps'}

set log x
set xlabel 'R / h^{-1} Mpc'
set format x '10^{%T}'

set log y
set ylabel '{/Symbol s}(R) or {/Symbol s}_v(R) / h^{-1} Mpc'
set format y '10^{%T}'

set key top right

plot 'sigma.dat' u 1:2 w l lw 2 ti 'sigma: z = 0',\
     'sigma.dat' u 1:3 w l lw 2 ti 'sigma: z = 1',\
     'sigmaV.dat' u 1:2 w l lw 2 ti 'sigmaV: z = 0',\
     'sigmaV.dat' u 1:3 w l lw 2 ti 'sigmaV: z = 1'
