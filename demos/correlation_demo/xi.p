reset

#rmin=0.1
#rmax=250.
set xlabel 'r / h^{-1} Mpc'
set log x
#set xrange [rmin:rmax]
set xrange [*:*]

ximin=-1.
ximax=1.
set ylabel '{/Symbol x}(r)'
#set yrange [ximin:ximax]
set yrange [*:*]

plot 0 w l lt -1 noti,\
     'xi.dat' u 1:2 w p pt 7 lc 1 noti,\
     'xi.dat' u 1:2 w l lw 2 lc 1 noti
