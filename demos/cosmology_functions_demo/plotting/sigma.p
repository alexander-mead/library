reset

if(!exists('print')){print=0}
#if(print==0){set term aqua}
if(print==0){set term qt}
if(print==1){set term post enh col; set output 'power.eps'}

sigma='data/sigma.dat'
sigmaV='data/sigmaV.dat'

sigmin=1e-3
sigmax=1e2

rmin=0.97
rmax=1.03

set log x

set multiplot layout 2, 1 margins 0.06,0.98,0.06,0.98

set key bottom left

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol s}(R) or {/Symbol s}_v(R) / h^{-1} Mpc'
set format y '10^{%T}'
set yrange [sigmin:sigmax]

plot sigma  u 1:2 w l lw 3 lc 1 dt 1 ti 'sigma: all matter',\
     sigma  u 1:3 w l lw 2 lc 2 dt 1 ti 'sigma: cold matter',\
     sigma  u 1:4 w l lw 2 lc 3 dt 1 ti 'sigma: cold matter (un-normalised)',\
     sigmaV u 1:2 w l lw 3 lc 1 dt 2 ti 'sigmaV: all matter',\
     sigmaV u 1:3 w l lw 2 lc 2 dt 2 ti 'sigmaV: cold matter',\
     sigmaV u 1:4 w l lw 2 lc 3 dt 2 ti 'sigmaV: cold matter (un-normalised)'

set xlabel 'R / h^{-1} Mpc'
set format x '10^{%T}'

unset log y
set yrange [rmin:rmax]
set format y

plot 1 w l lt -1 noti,\
      sigma  u 1:($3/$2) w l lw 2 lc 2 dt 1 noti,\
      sigma  u 1:($4/$2) w l lw 2 lc 3 dt 1 noti,\
      sigmaV u 1:($3/$2) w l lw 2 lc 2 dt 2 noti,\
      sigmaV u 1:($4/$2) w l lw 2 lc 3 dt 2 noti

unset multiplot
