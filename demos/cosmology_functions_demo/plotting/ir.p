reset

# Files
file = 'data/ir.dat'

# k axis
klab = 'k / h Mpc^{-1}'

# Power axis
dlab = '{/Symbol D}^2(k)'
dmin = 1e-3
dmax = 1e1

# Ratio axis
rlab = 'P_{IR}(k) / P_{lin}(k)'
rmin = 0.95
rmax = 1.15

# x axis
set log x

set multiplot layout 2, 1 margins 0.1, 0.98, 0.1, 0.98

   set xlabel ''
   set format x ''

   set log y
   set format y '10^{%T}'
   set ylabel dlab
   set yrange [dmin:dmax]

   set key top left

   plot file u 1:2 w l lw 3 ti 'Linear',\
      file u 1:3 w l lw 3 ti 'IR resummation'

   set xlabel klab
   set format x

   unset log y
   set format y
   set ylabel rlab
   set yrange [rmin:rmax]

   unset key

   plot 1 w l lt -1 noti,\
      file u 1:($3/$2) w l lw 3

unset multiplot