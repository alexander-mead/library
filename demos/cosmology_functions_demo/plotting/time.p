reset

#set term aqua
set term qt

time='data/time.dat'

set xrange [0:1]
set xlabel 'a'

set yrange [0:10]
set ylabel 't / (h^{-1} Gyr)'

set key top left

plot time u 1:2 w l lw 3 ti 'Universe age',\
     time u 1:3 w l lw 3 ti 'Look-back time'
