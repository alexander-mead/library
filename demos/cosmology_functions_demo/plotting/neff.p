reset

#set term aqua dashed
set term qt dashed

file = 'data/neff.dat'

nmin=-3.
nmax=1.

rmin=-0.015
rmax=0.015

set log x

set multiplot layout 2, 1 margins 0.08, 0.98, 0.1, 0.98

set key top left

set xlabel ''
set format x ''

set ylabel 'n_{eff}(R)'
set yrange[nmin:nmax]

plot file u 1:2 w l lw 3 lc 1 ti 'All matter',\
     file u 1:3 w l lw 2 lc 2 ti 'Cold matter'

set xlabel 'R / h^{-1} Mpc'
set format x '10^{%T}'

set ylabel 'n_{eff,cold}(R) - n_{eff,all}(R)'
set yrange [rmin:rmax]

plot 0 w l lt -1 noti,\
   file u 1:($3-$2) w l lc 2 lw 2 noti

unset multiplot