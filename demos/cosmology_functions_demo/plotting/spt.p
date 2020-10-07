reset

# Files
file = 'data/spt.dat'

# k axis
klab = 'k / h Mpc^{-1}'
kmin = 1e-2
kmax = 0.5
kmax_fit = 0.15

# Power axis
dlab = '{/Symbol D}^2(k)'
dmin = 1e-3
dmax = 1e1

# Ratio axis
rlab = 'P_{SPT}(k) / P_{lin}(k)'
rmin = 0.95
rmax = 1.15

# x axis
set log x
set xrange [kmin:kmax]

# Fitting function
#R(k) = f*(k-ks)**2/ks**2+1.-f
R(k) = 1.-A*(k/kr)+B*(k/kt)**3

# Fitting function parameters
#f = 0.015
#ks = 0.04
A = 0.035
kr = 0.08
B = 0.5
kt = 0.2

# Fit the fitting function
fit [x=kmin:kmax_fit] R(x) file u 1:(($2+$3)/$2) via A, B

set multiplot layout 2, 1 margins 0.1, 0.98, 0.1, 0.98

   set xlabel ''
   set format x ''

   set log y
   set format y '10^{%T}'
   set ylabel dlab
   set yrange [dmin:dmax]

   set key top left

   plot file u 1:2 w l lc 1 lw 3 ti 'Linear',\
      file u 1:($2+$3) w l lc 2 lw 3 ti 'Linear plus one loop',\
      file u 1:($2+$4) w l lc 3 lw 3 ti 'Approximate linear plus one loop',\
      file u 1:($2*R($1)) w p pt 7 ps .5 lc 3 dt 2 lw 3 ti 'Fitted functional form'

   set xlabel klab
   set format x

   unset log y
   set format y
   set ylabel rlab
   set yrange [rmin:rmax]

   unset key

   plot 1 w l lt -1 noti,\
      file u 1:(($2+$3)/$2) w l lc 2 lw 3,\
      file u 1:(($2+$4)/$2) w l lc 3 lw 3,\
      R(x) w p pt 7 ps .5 lw 3 lc 3 dt 2

unset multiplot