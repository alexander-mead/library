unset multiplot
reset

if(!exists('print')){print=0}
#if(print==0) set term aqua dashed
if(print==0) set term qt dashed

dc_file='data/delta_c.dat'
Dv_file='data/Delta_v.dat'

#Scale factor a axis
axis=1
if(axis==1){
amin=0.01
amax=1.
c=1
set log x
set xlabel 'a'
set xrange [amin:amax]
dcmin=1.64
dcmax=1.69
Dvmin=100.
Dvmax=400.
}

#Omega_m x axis
if(axis==2){
c=2
Om_min=0.08
Om_max=1.
set log x
set xlabel '{/Symbol W}_m(a)'
set xrange [Om_min:Om_max]
dcmin=1.64
dcmax=1.69
Dvmin=100.
Dvmax=800.
}

#Exact values for EdS cosmology
dc0=(3./20.)*(12.*pi)**(2./3.)
Dv0=18.*pi**2

set multiplot layout 1,2

#if(axis==1){set key top left}
#if(axis==2){set key bottom right}
set key bottom right

set ylabel '{/Symbol d}_c(a)'
set yrange [dcmin:dcmax]

plot dc0 w l lw 2 lc -1 ti '{/Symbol d}_c = 1.686...',\
     dc_file u 1:2 w l lw 2 ti 'Nakamura & Suto (1998)',\
     dc_file u 1:3 w l lw 2 ti 'Mead (2017)',\
     dc_file u 1:4 w l lw 2 ti 'Calculation'

set ylabel '{/Symbol D}_v(a)'
set yrange [Dvmin:Dvmax]

plot Dv0 w l lw 3 lc -1 ti '{/Symbol D}_v = 177.7...',\
     Dv_file u 1:2 w l lw 2 ti 'Bryan & Norman (1998)',\
     Dv_file u 1:3 w l lw 2 ti 'Mead (2017)',\
     Dv_file u 1:4 w l lw 2 ti 'Calculation'

unset multiplot
