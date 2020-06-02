unset multiplot
reset

if(!exists('print')) {print=0}
if(print==0) set term qt size 900,900
if(print==1) set term post enh col sol font ',10'

# Initial white space
print ''

# Parameters
func_min = 0.
func_max = 4.
ratio_min = 0.95
ratio_max = 1.05
difference_min = -0.05
difference_max = 0.05

# Choose plot
if(!exist('iplot')){iplot=3}
print 'Set plot via iplot'
print 'iplot = 1: Truth'
print 'iplit = 2: Interpolation'
print 'iplot = 3: Ratio'
print 'iplot = ', iplot
print ''

set view map
set pm3d map

#dcb=0.001
#set cbrange [1-dcb:1+dcb]
#set cblabel 'Ratio of interpolation to truth'
set palette rgb 33,13,10

unset key

if(!exists('order')) {order=1}
if(order==0) {c=4}
if(order==1) {c=5}
if(order==3) {c=6}
file = 'results_2D.dat'
print 'order = ', order
print 'file = ', file
print ''

set multiplot layout 2, 2

set title 'Truth'
set cbrange [func_min:func_max]
splot file u 1:2:(column(3)) w image

set title 'Interpolation'
splot file u 1:2:(column(c)) w image

set title 'Ratio'
set cbrange [ratio_min:ratio_max]
splot file u 1:2:(column(c)/column(3)) w image

set title 'Difference'
set cbrange [difference_min:difference_max]
splot file u 1:2:(column(c)-column(3)) w image

unset multiplot