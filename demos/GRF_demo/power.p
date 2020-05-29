reset

CAMB_file='LCDM_matterpower.txt'
GRF_power='data/power.dat'

kmin=0.03
kmax=1.
set log x
set xlabel 'k / h Mpc^{-1}'
set xrange [kmin:kmax]

set log y
set ylabel '{/Symbol D}^2(k)'

set key top left

plot GRF_power u 1:2:5 w e lw 1.5 ti 'Measured',\
     CAMB_file u 1:((4.*pi*($1/(2.*pi))**3)*$2) w l lw 2 lc -1 ti 'Input'
