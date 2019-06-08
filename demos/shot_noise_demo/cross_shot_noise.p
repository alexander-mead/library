unset multiplot
reset

if(!exists('print')){print=0}
if(print==0){set term aqua dashed}
if(print==1){set term post enh col; set output 'cross_shot_noise.eps'}

set log x
set xrange [0.9:*]

set log y
set format y '10^{%T}'
set ylabel '{/Symbol D}^2(k)'

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set key top left

unset format x
set format x ''

plot 'power_11.dat' u ($1/(2.*pi)):2:5 w e lc 1 ti '1-1',\
     'power_12.dat' u ($1/(2.*pi)):2:5 w e lc 2 ti '1-2',\
     'power_22.dat' u ($1/(2.*pi)):2:5 w e lc 3 ti '2-2',\
     'power_11.dat' u ($1/(2.*pi)):3   w l lc 1 noti,\
     'power_12.dat' u ($1/(2.*pi)):3   w l lc 2 noti,\
     'power_22.dat' u ($1/(2.*pi)):3   w l lc 3 noti

set xlabel 'kL / 2{/Symbol p}'
set format x

set ylabel 'P(k) / P_{SN}(k)'
set yrange [0.:2.]
unset log y
set format y

plot 1 w l lt -1 noti,\
     'power_11.dat' u ($1/(2.*pi)):($2/$3):($5/$3) w e lc 1 noti,\
     'power_12.dat' u ($1/(2.*pi)):($2/$3):($5/$3) w e lc 2 noti,\
     'power_22.dat' u ($1/(2.*pi)):($2/$3):($5/$3) w e lc 3 noti

unset multiplot
