reset

if(!exists('print')) {print=0}
if(print==-1) set term x11
if(print==0) set term aqua size 1000,1000
if(print==1) set term post enh col sol font ',10'

set view map
set pm3d map

dcb=0.001
#set cbrange [1-dcb:1+dcb]
set cblabel 'Ratio of interpolation to truth'
set palette rgb 33,13,10

unset key

if(!exists('order')) {order=1}
if(order==1) file='results_linear.dat'; tit='linear interpolation'
if(order==3) file='results_cubic.dat'; tit='cubic interpolation'
print 'order = ', order
print 'file = ', file

set title tit

splot file u 1:2:5 w image


