reset

if(!exists('print')){print=0}
#if(print==0){set term aqua dashed dl 1}
if(print==0){set term qt dashed dl 1}
if(print==1){set term post enh col sol; set output 'Omegas.eps'}

file='data/omegas.dat'

amin=1e-5
amax=1.
set log x
set xrange [amin:amax]
set mxtics 10

Omega_min=1e-3
Omega_max=1.3
set log y
set ylabel '{/Symbol W}_i(a)'

set multiplot layout 2,1 margins 0.10,0.88,0.10,0.97 spacing 0.01,0.02

do for [i=1:2]{

    if(i==1){set log y; set format y '10^{%T}'; set xlabel ''; set format x ''; Omega_min=1e-3; Omega_max=1.3; set key outside }#box}
    if(i==2){unset log y; set format y; set xlabel 'a'; set format x '10^{%T}'; Omega_min=0.; Omega_max=1.1; unset key}

    set yrange [Omega_min:Omega_max]

    

    plot file u 1:10 w l lw 3 lc rgb 'black' ti '{/Symbol W}',\
        file u 1:2  w l lw 3 lc 1 dt 1 ti '{/Symbol W}_m',\
        file u 1:3  w l lw 3 lc 1 dt 2 ti '{/Symbol W}_c',\
        file u 1:4  w l lw 3 lc 1 dt 3 ti '{/Symbol W}_b',\
        file u 1:5  w l lw 3 lc 2 dt 1 ti '{/Symbol W}_r',\
        file u 1:6  w l lw 3 lc 2 dt 2 ti '{/Symbol W}_{/Symbol g}',\
        file u 1:7  w l lw 3 lc 5 dt 1 ti '{/Symbol W}_{/Symbol n}',\
        file u 1:8  w l lw 3 lc 3 dt 1 ti '{/Symbol W}_{/Symbol L}',\
        file u 1:9  w l lw 3 lc 4 dt 1 ti '{/Symbol W}_w'

}

unset multiplot

