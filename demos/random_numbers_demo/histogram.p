reset

data='histogram.dat'

set xlabel 'x'

set ylabel 'P(x)'
#unset ytics
#set format y ''
set yrange [0:*]

plot data u 1:2 w l lw 3 ti 'Random numbers',\
     data u 1:3 w l lw 3 ti 'Underlying function'
