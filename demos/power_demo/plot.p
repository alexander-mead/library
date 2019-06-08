unset multiplot
reset

real_power='power_real.dat'
cplx_power='power_complex.dat'

set log x
set xrange [0.9:*]

set key top left

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

s=1.01

plot cplx_power u (s*$1/(2.*pi)):2:5 w e lw 2 ti 'Complex (shifted in x for clarity)',\
     real_power u ($1/(2.*pi)):2:5   w e lw 2 ti 'Real'

set xlabel 'kL / 2{/Symbol p}'
set format x

unset log y
set ylabel 'P_{complex}(k) / P_{real}(k)'
set format y
dx=1e-3
#set yrange [1-dx:1+dx]

plot '<paste '.cplx_power.' '.real_power.'' u ($1/(2.*pi)):($2/$7) w p noti

unset multiplot
