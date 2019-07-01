# Returns maximum of x, y
max(x, y) = (x > y ? x : y)

# Linear spacing between xmin and xmax
progression(xmin,xmax,i,n) = xmin+(xmax-xmin)*(real(i)-1)/(real(n)-1)