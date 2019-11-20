reset

if(!exists('print')){print=0}
#if(print==0){set term aqua dashed dl 1}
if(print==0){set term qt dl 1}
if(print==1){set term post enh col sol; set output 'Omegas.eps'}

density='data/density.dat'

amin=1e-5
amax=1.
set log x
set xrange [amin:amax]
set xlabel 'a'
set format x '10^{%T}'
set mxtics 10

rho_min=1e-5
rho_max=1e15
set log y
set format y '10^{%T}'
set yrange [rho_min:rho_max]
set ylabel '{/Symbol r}_i(a)'

set key outside

plot density u 1:10 w l lw 3 lc rgb 'black' ti '{/Symbol W}',\
     density u 1:2  w l lw 3 lc 1 dt 1 ti '{/Symbol W}_m',\
     density u 1:3  w l lw 3 lc 1 dt 2 ti '{/Symbol W}_c',\
     density u 1:4  w l lw 3 lc 1 dt 3 ti '{/Symbol W}_b',\
     density u 1:5  w l lw 3 lc 2 dt 1 ti '{/Symbol W}_r',\
     density u 1:6  w l lw 3 lc 2 dt 2 ti '{/Symbol W}_{/Symbol g}',\
     density u 1:7  w l lw 3 lc 5 dt 1 ti '{/Symbol W}_{/Symbol n}',\
     density u 1:8  w l lw 3 lc 3 dt 1 ti '{/Symbol W}_{/Symbol L}',\
     density u 1:9  w l lw 3 lc 4 dt 1 ti '{/Symbol W}_w'

