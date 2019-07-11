# Returns maximum of x, y
max(x, y) = (x > y ? x : y)

# Linear spacing between xmin and xmax
progression(xmin,xmax,i,n) = xmin+(xmax-xmin)*(real(i)-1)/(real(n)-1)

#color_scheme(c1,c2,c3)=sprintf('%d, %d, %d', c1, c2, c3)

# Colour schemes
#array rainbow[3]
#rainbow[1] = 33
#rainbow[2] = 13
#rainbow[3] = 10

#array AFM_hot[3]
#AFM_hot[1] = 34
#AFM_hot[2] = 35
#AFM_hot[3] = 36