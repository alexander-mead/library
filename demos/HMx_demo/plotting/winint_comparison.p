unset multiplot
reset

set term qt dashed

n=13    # Number of methods
ibase=4 # Base-line method (should be best)

# Minimum and maximum for W(k) plot
wmin=1e-5
wmax=2.

# Minimum and maximum for the residual
rmin=0.99
rmax=1.01

set log x
set mxtics 10

file(i)=sprintf('data/results_%d.dat',i)

set key top left

set lmargin 10
set rmargin 2

do for [j=1:n] {

set multiplot layout 2,1

set xlabel ''
set format x ''

set log y
set ylabel 'W(k)'
set yrange [wmin:wmax]
set format y '10^{%T}'

set key bottom left

plot for [i=1:j] file(i) u 1:(+$2) w l dt 1 lc i lw 3 ti 'Method '.i.'',\
     for [i=1:j] file(i) u 1:(-$2) w l dt 2 lc i lw 3 noti

set format x
set xlabel 'k / h Mpc^{-1}'
set format x '10^{%T}'


unset log y
set ylabel 'W_{i}(k) / W_{base}(k)'
set yrange [rmin:rmax]
set format y

error=1e-3
plot 1.+error w l dt 2 lc -1 noti,\
     1.-error w l dt 2 lc -1 noti,\
     for [i=1:j] '<paste '.file(i).' '.file(ibase).'' u 1:($2/$4) w l dt 1 lc i lw 3 noti

unset multiplot

if(j<n) {pause 0.5}

}
