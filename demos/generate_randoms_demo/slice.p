reset

set size square

slices="'slice_randoms.dat' 'slice_grid.dat' 'slice_poorglass.dat'"

do for [i=1:3] {
plot word(slices,i) u 1:2 w p lc i ti word(slices,i)
pause -1
}

