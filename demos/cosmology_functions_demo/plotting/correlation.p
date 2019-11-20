reset

#set term aqua dashed
set term qt dashed

print ''

corr='data/correlation.dat'

if(!exists('iplot')){iplot=2}
print 'iplot = 1: Plot xi(r)'
print 'iplot = 2: Plot r^2 x(i)'
print 'iplot = ', iplot
print ''

set log x
set xlabel 'r / h^{-1} Mpc'
set format x '10^{%T}'

if(iplot==1) {set ylabel '{/Symbol x}(r)'; set log y; set format y '10^{%T}'}
if(iplot==2) {set ylabel '4{/Symbol p}r^2{/Symbol x}(r) [Mpc/h]^2'}

unset key

if(iplot==1){
plot corr u 1:2 w l lc 1 dt 1 lw 3
}

if(iplot==2){
plot corr u 1:3 w l lc 1 dt 1 lw 3
}
