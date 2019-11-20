reset

if(!exists('print')){print=0}
#if(print==0){set term aqua}
if(print==0){set term qt}
if(print==1){set term post enh col; set output 'power.eps'}

sigma='data/sigma.dat'
sigmaV='data/sigmaV.dat'

sigmin=1e-3
sigmax=1e2

set log x
set xlabel 'R / h^{-1} Mpc'
set format x '10^{%T}'

set log y
set ylabel '{/Symbol s}(R) or {/Symbol s}_v(R) / h^{-1} Mpc'
set format y '10^{%T}'
set yrange [sigmin:sigmax]

set key top right

plot sigma u 1:2 w l lw 2 ti 'sigma: z = 0',\
     sigma u 1:3 w l lw 2 ti 'sigma: z = 1',\
     sigmaV u 1:2 w l lw 2 ti 'sigmaV: z = 0',\
     sigmaV u 1:3 w l lw 2 ti 'sigmaV: z = 1'
