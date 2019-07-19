reset

file = 'data/neff.dat'

set log x
set xlabel 'R / h^{-1} Mpc'
set format x '10^{%T}'

nmin=-3.
nmax=1.
set ylabel 'n_{eff}(R)'
set yrange[nmin:nmax]

plot file u 1:2 w l lw 3 noti,\
   file u 1:3 w l lw 3 noti