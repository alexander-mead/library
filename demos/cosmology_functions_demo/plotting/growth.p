reset

if(!exists('print')){print=0}
if(print==0) set term aqua dashed

set key bottom right

#if(!exists('ilog')){ilog=0}

set multiplot layout 1,2

do for [ilog=0:1]{

if(ilog==0){amin=0.; gmin=0.}
if(ilog==1){amin=1e-4; set log x; set log y}

if(ilog==0){amin=0.}
if(ilog==1){amin=1e-4; set log x; set format x '10^{%T}'}
amax=1.
set xlabel 'a'
set xrange [amin:amax]

if(ilog==0){gmin=0.}
if(ilog==1){gmin=1e-4; set log y}
gmax=1.05
set ylabel 'g(a)'
set yrange [gmin:gmax]

if(ilog==1){unset key}

plot x w l ls -1 noti,\
     'growth.dat' u 1:2 w l lc 1 dt 1 lw 2 ti 'Growth function',\
     'growth.dat' u 1:3 w l lc 2 dt 1 lw 2 ti 'Unnormalised growth',\
     'growth.dat' u 1:4 w l lc 3 dt 1 lw 2 ti 'Growth rate',\
     'growth.dat' u 1:5 w l lc 4 dt 1 lw 2 ti 'Accumulated growth',\
     'growth.dat' u 1:6 w l lc 1 dt 2 lw 3 ti 'Linder growth function approximation',\
     'growth.dat' u 1:7 w l lc 1 dt 3 lw 3 ti 'CPT growth function approximation',\
     'growth.dat' u 1:8 w l lc 3 dt 2 lw 3 ti 'Linder growth rate approximation'

}

unset multiplot