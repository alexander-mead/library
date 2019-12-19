reset

file = 'data/dewiggle.dat'

klab = 'k / h Mpc^{-1}'

plab = '{/Symbol D}^2(k)'

rlab = 'P_{no-wiggle}(k) / P_{lin}(k)'
rmin = 0.9
rmax = 1.1

set log x

set multiplot layout 2, 1 margins 0.05,0.98,0.07,0.98

set format x ''

set log y
set format y '10^{%T}'
set ylabel plab

set key top left

plot file u 1:2 w l lc -1 lw 2 ti 'Linear',\
     file u 1:3 w l lc 1  lw 2 ti 'No-wiggle'

set format x
set xlabel klab

unset log y
set format y 
set ylabel rlab
set yrange [rmin:rmax]

plot 1 w l lt -1 noti,\
   file u 1:($2/$3) w l lc 1 lw 2 noti

unset multiplot