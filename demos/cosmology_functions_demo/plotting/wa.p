reset

#if(print==0) set term aqua dashed
if(print==0) set term qt dashed
if(print==1) set term post enh col; set output 'wa.eps'

wa='data/wa.dat'

set key top left

amin=1e-5
amax=1e0
set log x
set xrange [amin:amax]
set xlabel 'a'
set format x '10^{%T}'

set yrange [-1.55:1.55]
set ylabel 'w(a)'

plot 0 w l lt -1 noti,\
     1 w l dt 2 lc 'black' noti,\
     -1 w l dt 2 lc 'black' noti,\
     -1 w l lw 3 ti 'w_{vac}',\
     wa u 1:2 w l lw 3 ti 'w(a): {/Symbol W}_w only',\
     wa u 1:3 w l lw 3 ti 'w(a): {/Symbol W}_w and {/Symbol W}_v weighted',\
     wa u 1:4 w l lw 3 ti 'w_{eff}'
