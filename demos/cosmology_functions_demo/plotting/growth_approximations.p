reset

#set term aqua dashed
set term qt dashed

growth='data/growth.dat'

set xlabel 'a'

set ylabel 'Growth approximation / truth'

set multiplot layout 1,2

do for [ilog=0:1]{

if(ilog==1){set log x}

plot 1 w l lt -1 noti,\
     growth u 1:($6/$2) w l lc 1 dt 2 lw 2 ti 'Linder growth',\
     growth u 1:($7/$2) w l lc 1 dt 3 lw 2 ti 'CPT growth',\
     growth u 1:($8/$4) w l lc 3 dt 2 lw 2 ti 'Linder rate'

}

unset multiplot
