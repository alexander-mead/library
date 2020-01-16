reset

if(!exists('print')){print=0}
#if(print==0){set term aqua}
if(print==0){set term qt}
if(print==1){set term post enh col; set output 'power.eps'}

power='data/power.dat'

klab='k / h Mpc^{-1}'

pmin=1e-10
pmax=1e3
plab='{/Symbol D}^2(k)'

rmin=0.95
rmax=1.05
rlab='P_{cold}(k) / P_{total}(k)'

set log x

set key top left

set multiplot layout 2, 1 margins 0.06,0.98,0.06,0.98

set xlabel ''
set format x ''

set log y
set ylabel 
set format y '10^{%T}'
set yrange [pmin:pmax]

plot power u 1:2 w l lw 3 lc 1 ti 'All matter',\
     power u 1:3 w l lw 2 lc 2 ti 'Cold matter',\
     power u 1:4 w l lw 2 lc 3 ti 'Cold matter (un-normalised)'

set xlabel klab
set format x

unset log y
set yrange [rmin:rmax]
set format y
set ylabel rlab

plot 1 w l lt -1 noti,\
   power u 1:($3/$2) w l lw 2 lc 2 noti,\
   power u 1:($4/$2) w l lw 2 lc 3 noti

unset multiplot

show output
