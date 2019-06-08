reset

set term aqua dashed

set log x
set xlabel 'kL / 2{/Symbol p}'
set xrange [0.9:*]

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

set key top left

f(a,x,n)=a*x**n

a3=1e-4
a5=1e-6

plot f(a3,x,3) w l lc 1 dt 2 ti '{/Symbol a} k^3',\
     f(a5,x,5) w l lc 3 dt 2 ti '{/Symbol a} k^5',\
     'power_randoms.dat'   u ($1/(2.*pi)):2 w p lw 2 lc 1 ti 'Random',\
     'power_grid.dat'      u ($1/(2.*pi)):2 w p lw 2 lc 2 ti 'Grid',\
     'power_poorglass.dat' u ($1/(2.*pi)):2 w p lw 2 lc 3 ti 'Poor glass',\
     'power_randoms.dat'   u ($1/(2.*pi)):3 w l lc -1 dt 2 ti 'Shot noise'
