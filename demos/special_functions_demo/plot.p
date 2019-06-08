reset

#set xrange [0:10]
set xlabel 'x'

#set yrange [0:1.05]
set ylabel 'f(x)'

plot 0 w l lt -1 noti,\
     'results.dat' u 1:2 w l lw 3 noti#,\
     0.5*(1.-cos(2.*pi*x/10.))
