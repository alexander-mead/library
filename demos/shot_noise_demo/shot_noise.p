reset

if(!exists('print')){print=0}
if(print==0){set term aqua dashed}

set log x
set xrange [0.9:*]

set key top left

#set title 'I do not know why the shot noise is always slightly low'

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

plot 'power.dat' u ($1/(2.*pi)):2:5 w e lw 1 lc 1 dt 1 ti 'Measurements',\
     'power.dat' u ($1/(2.*pi)):3 w l   lw 2 lc 1 dt 2 ti 'Simple shot noise calculation',\
     'power.dat' u ($1/(2.*pi)):6 w l   lw 2 lc 1 dt 3 ti 'Fancy shot noise calculation'

unset log y
set format y
set yrange [0.5:1.5]

set xlabel 'kL / 2{/Symbol p}'
set format x '10^{%T}'

plot 1 w l lt -1 noti,\
     'power.dat' u ($1/(2.*pi)):($2/$6):($5/$6) w e lc 1 noti

unset multiplot
