MODULE interpolate

  USE fix_polynomial
  USE table_integer
  USE array_operations

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: find
  PUBLIC :: interpolate_array

  INTERFACE find
     MODULE PROCEDURE find_1D
     MODULE PROCEDURE find_2D
     MODULE PROCEDURE find_3D
  END INTERFACE find

CONTAINS

  REAL FUNCTION find_1D(x,xin,yin,n,iorder,ifind,imeth)

    ! Given two arrays x and y this routine interpolates to find the y_i value at position x_i
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN) :: xin(n)
    REAL, INTENT(IN) :: yin(n)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: iorder
    INTEGER, INTENT(IN) :: ifind
    INTEGER, INTENT(IN) :: imeth
    REAL ::  xtab(n), ytab(n)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i  

    ! This version interpolates if the value is off either end of the array
    ! Care should be chosen to insert x, xtab, ytab as log if this might give beter results from the interpolation

    ! If the value required is off the table edge the interpolation is always linear

    ! iorder = 1 => linear interpolation
    ! iorder = 2 => quadratic interpolation
    ! iorder = 3 => cubic interpolation

    ! ifind = 1 => find x in xtab quickly assuming the table is linearly spaced
    ! ifind = 2 => find x in xtab by crudely searching from x(1) to x(n)
    ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

    ! imeth = 1 => Uses standard polynomials for interpolation
    ! imeth = 2 => Uses Lagrange polynomials for interpolation

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       ! Reverse the arrays in this case
       CALL reverse_array(xtab,n)
       CALL reverse_array(ytab,n)
    END IF

    IF(x<xtab(1)) THEN

       ! Do a linear interpolation beyond the table boundary

       x1=xtab(1)
       x2=xtab(2)

       y1=ytab(1)
       y2=ytab(2)

       IF(imeth==1) THEN
          CALL fix_line(a,b,x1,y1,x2,y2)
          find_1D=a*x+b
       ELSE IF(imeth==2) THEN
          find_1D=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND_1D: Error, method not specified correctly'
       END IF

    ELSE IF(x>xtab(n)) THEN

       ! Do a linear interpolation beyond the table boundary

       x1=xtab(n-1)
       x2=xtab(n)

       y1=ytab(n-1)
       y2=ytab(n)

       IF(imeth==1) THEN
          CALL fix_line(a,b,x1,y1,x2,y2)
          find_1D=a*x+b
       ELSE IF(imeth==2) THEN
          find_1D=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND_1D: Error, method not specified correctly'
       END IF

    ELSE IF(iorder==1) THEN

       IF(n<2) STOP 'FIND_1D: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x1=xtab(1)
          x2=xtab(2)

          y1=ytab(1)
          y2=ytab(2)

       ELSE IF (x>=xtab(n-1)) THEN

          x1=xtab(n-1)
          x2=xtab(n)

          y1=ytab(n-1)
          y2=ytab(n)

       ELSE

          i=select_table_integer(x,xtab,n,ifind)

          x1=xtab(i)
          x2=xtab(i+1)

          y1=ytab(i)
          y2=ytab(i+1)

       END IF

       IF(imeth==1) THEN
          CALL fix_line(a,b,x1,y1,x2,y2)
          find_1D=a*x+b
       ELSE IF(imeth==2) THEN
          find_1D=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND_1D: Error, method not specified correctly'
       END IF

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND_1D: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x1=xtab(1)
             x2=xtab(2)
             x3=xtab(3)

             y1=ytab(1)
             y2=ytab(2)
             y3=ytab(3)

          ELSE IF (x>=xtab(n-1)) THEN

             x1=xtab(n-2)
             x2=xtab(n-1)
             x3=xtab(n)

             y1=ytab(n-2)
             y2=ytab(n-1)
             y3=ytab(n)

          END IF

          IF(imeth==1) THEN
             CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
             find_1D=a*(x**2)+b*x+c
          ELSE IF(imeth==2) THEN
             find_1D=Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))
          ELSE
             STOP 'FIND_1D: Error, method not specified correctly'
          END IF

       ELSE

          i=select_table_integer(x,xtab,n,ifind)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          IF(imeth==1) THEN
             ! In this case take the average of two separate quadratic spline values
             CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
             find_1D=(a*x**2+b*x+c)/2.
             CALL fix_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
             find_1D=find_1D+(a*x**2+b*x+c)/2.
          ELSE IF(imeth==2) THEN
             ! In this case take the average of two quadratic Lagrange polynomials
             find_1D=(Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))+Lagrange_polynomial(x,2,(/x2,x3,x4/),(/y2,y3,y4/)))/2.
          ELSE
             STOP 'FIND_1D: Error, method not specified correctly'
          END IF

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND_1D: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x1=xtab(1)
          x2=xtab(2)
          x3=xtab(3)
          x4=xtab(4)        

          y1=ytab(1)
          y2=ytab(2)
          y3=ytab(3)
          y4=ytab(4)

       ELSE IF (x>=xtab(n-2)) THEN

          x1=xtab(n-3)
          x2=xtab(n-2)
          x3=xtab(n-1)
          x4=xtab(n)

          y1=ytab(n-3)
          y2=ytab(n-2)
          y3=ytab(n-1)
          y4=ytab(n)

       ELSE

          i=select_table_integer(x,xtab,n,ifind)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       IF(imeth==1) THEN
          CALL fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
          find_1D=a*x**3+b*x**2+c*x+d
       ELSE IF(imeth==2) THEN
          find_1D=Lagrange_polynomial(x,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
       ELSE
          STOP 'FIND_1D: Error, method not specified correctly'
       END IF

    ELSE

       STOP 'FIND_1D: Error, interpolation order specified incorrectly'

    END IF

  END FUNCTION find_1D

  REAL FUNCTION find_2D(x,xin,y,yin,fin,nx,ny,iorder,ifind,imeth)

    ! A 2D interpolation routine to find value f(x,y) at position x, y
    ! TODO: Loops over coordinates to avoid repetition?
    ! TOOD: Check linear method, there is a nice thing called bilinear interpolation that has a nice geometric interpretation, see find_3D
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN) :: xin(nx)
    REAL, INTENT(IN) :: y
    REAL, INTENT(IN) :: yin(ny)
    REAL, INTENT(IN) :: fin(nx,ny)
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    INTEGER, INTENT(IN) :: iorder
    INTEGER, INTENT(IN) :: ifind
    INTEGER, INTENT(IN) :: imeth
    REAL ::  xtab(nx), ytab(ny), ftab(nx,ny)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    REAL :: f11, f12, f13, f14
    REAL :: f21, f22, f23, f24
    REAL :: f31, f32, f33, f34
    REAL :: f41, f42, f43, f44
    REAL :: f10, f20, f30, f40
    REAL :: f01, f02, f03, f04
    INTEGER :: i1, i2, i3, i4
    INTEGER :: j1, j2, j3, j4
    REAL :: findx, findy
    INTEGER :: i, j

    ! This version interpolates if the value is off either end of the array!
    ! Care should be chosen to insert x, xtab, ytab as log if this might give better!
    ! Results from the interpolation!

    ! If the value required is off the table edge the interpolation is always linear

    ! iorder = 1 => linear interpolation
    ! iorder = 2 => quadratic interpolation
    ! iorder = 3 => cubic interpolation

    ! ifind = 1 => find x in xtab by crudely searching from x(1) to x(n)
    ! ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
    ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

    ! imeth = 1 => Uses cubic polynomials for interpolation
    ! imeth = 2 => Uses Lagrange polynomials for interpolation

    IF(imeth==2) STOP 'FIND_2D: No Lagrange polynomials for you'

    xtab=xin
    ytab=yin
    ftab=fin

    IF(xtab(1)>xtab(nx)) STOP 'FIND_2D: x array in wrong order'
    IF(ytab(1)>ytab(ny)) STOP 'FIND_2D: y array in wrong order'

    IF((x<xtab(1) .OR. x>xtab(nx)) .AND. (y>ytab(ny) .OR. y<ytab(1))) THEN
       WRITE(*,*) 'FIND_2D: array xmin:', xtab(1)
       WRITE(*,*) 'FIND_2D: array xmax:', xtab(nx)
       WRITE(*,*) 'FIND_2D: requested x:', x
       WRITE(*,*) 'FIND_2D: array ymin:', ytab(1)
       WRITE(*,*) 'FIND_2D: array ymax:', ytab(ny)
       WRITE(*,*) 'FIND_2D: requested y:', y
       STOP 'FIND_2D: Desired point is outside x AND y array range'
    END IF

    IF(iorder==1) THEN

       ! TODO: This might be wasteful or even wrong, there is a nice clear geometric way of doing this (https://en.wikipedia.org/wiki/Bilinear_interpolation)

       IF(nx<2) STOP 'FIND_2D: Not enough x points in your array for linear interpolation'
       IF(ny<2) STOP 'FIND_2D: Not enough y points in your array for linear interpolation'

       IF(x<=xtab(2)) THEN

          i=1

       ELSE IF (x>=xtab(nx-1)) THEN

          i=nx-1

       ELSE

          i=select_table_integer(x,xtab,nx,ifind)

       END IF

       i1=i
       i2=i+1

       x1=xtab(i1)
       x2=xtab(i2)      

       IF(y<=ytab(2)) THEN

          j=1

       ELSE IF (y>=ytab(ny-1)) THEN

          j=ny-1

       ELSE

          j=select_table_integer(y,ytab,ny,ifind)

       END IF

       j1=j
       j2=j+1

       y1=ytab(j1)
       y2=ytab(j2)

       !!

       f11=ftab(i1,j1)
       f12=ftab(i1,j2)

       f21=ftab(i2,j1)
       f22=ftab(i2,j2)

       !! y direction interpolation

       CALL fix_line(a,b,x1,f11,x2,f21)
       f01=a*x+b

       CALL fix_line(a,b,x1,f12,x2,f22)
       f02=a*x+b

       CALL fix_line(a,b,y1,f01,y2,f02)
       findy=a*y+b

       !!

       !! x direction interpolation

       CALL fix_line(a,b,y1,f11,y2,f12)
       f10=a*y+b

       CALL fix_line(a,b,y1,f21,y2,f22)
       f20=a*y+b

       CALL fix_line(a,b,x1,f10,x2,f20)
       findx=a*x+b

       !!

       ! Final result is an average over each direction
       find_2D=(findx+findy)/2.

    ELSE IF(iorder==2) THEN

       STOP 'FIND_2D: Quadratic 2D interpolation not implemented - also probably pointless'

    ELSE IF(iorder==3) THEN

       IF(x<xtab(1) .OR. x>xtab(nx)) THEN

          IF(nx<2) STOP 'FIND_2D: Not enough x points in your array for linear interpolation'
          IF(ny<4) STOP 'FIND_2D: Not enough y points in your array for cubic interpolation'

          ! x is off the table edge

          IF(x<xtab(1)) THEN

             i1=1
             i2=2

          ELSE

             i1=nx-1
             i2=nx

          END IF

          x1=xtab(i1)
          x2=xtab(i2)

          IF(y<=ytab(4)) THEN

             j=2

          ELSE IF (y>=ytab(ny-3)) THEN

             j=ny-2

          ELSE

             j=select_table_integer(y,ytab,ny,ifind)

          END IF

          j1=j-1
          j2=j
          j3=j+1
          j4=j+2

          y1=ytab(j1)
          y2=ytab(j2)
          y3=ytab(j3)
          y4=ytab(j4)

          f11=ftab(i1,j1)
          f12=ftab(i1,j2)
          f13=ftab(i1,j3)
          f14=ftab(i1,j4)

          f21=ftab(i2,j1)
          f22=ftab(i2,j2)
          f23=ftab(i2,j3)
          f24=ftab(i2,j4)

          !! y interpolation
          
          CALL fix_cubic(a,b,c,d,y1,f11,y2,f12,y3,f13,y4,f14)
          f10=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f21,y2,f22,y3,f23,y4,f24)
          f20=a*y**3+b*y**2+c*y+d

          !!

          !! x interpolation
          
          CALL fix_line(a,b,x1,f10,x2,f20)
          find_2D=a*x+b

          !!

       ELSE IF(y<ytab(1) .OR. y>ytab(ny)) THEN

          ! y is off the table edge

          IF(nx<4) STOP 'FIND_2D: Not enough x points in your array for cubic interpolation'
          IF(ny<2) STOP 'FIND_2D: Not enough y points in your array for linear interpolation'

          IF(x<=xtab(4)) THEN

             i=2

          ELSE IF (x>=xtab(nx-3)) THEN

             i=nx-2

          ELSE

             i=select_table_integer(x,xtab,nx,ifind)

          END IF

          i1=i-1
          i2=i
          i3=i+1
          i4=i+2

          x1=xtab(i1)
          x2=xtab(i2)
          x3=xtab(i3)
          x4=xtab(i4)

          IF(y<ytab(1)) THEN

             j1=1
             j2=2

          ELSE

             j1=ny-1
             j2=ny

          END IF

          y1=ytab(j1)
          y2=ytab(j2)

          f11=ftab(i1,j1)
          f21=ftab(i2,j1)
          f31=ftab(i3,j1)
          f41=ftab(i4,j1)

          f12=ftab(i1,j2)
          f22=ftab(i2,j2)
          f32=ftab(i3,j2)
          f42=ftab(i4,j2)

          ! x interpolation

          CALL fix_cubic(a,b,c,d,x1,f11,x2,f21,x3,f31,x4,f41)
          f01=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f12,x2,f22,x3,f32,x4,f42)
          f02=a*x**3+b*x**2+c*x+d

          ! y interpolation

          CALL fix_line(a,b,y1,f01,y2,f02)
          find_2D=a*y+b

       ELSE

          ! Points exists within table boundardies (normal)

          IF(nx<4) STOP 'FIND_2D: Not enough x points in your array for cubic interpolation'
          IF(ny<4) STOP 'FIND_2D: Not enough y points in your array for cubic interpolation'

          IF(x<=xtab(4)) THEN

             i=2

          ELSE IF (x>=xtab(nx-3)) THEN

             i=nx-2

          ELSE

             i=select_table_integer(x,xtab,nx,ifind)

          END IF

          i1=i-1
          i2=i
          i3=i+1
          i4=i+2

          x1=xtab(i1)
          x2=xtab(i2)
          x3=xtab(i3)
          x4=xtab(i4)

          IF(y<=ytab(4)) THEN

             j=2

          ELSE IF (y>=ytab(ny-3)) THEN

             j=ny-2

          ELSE

             j=select_table_integer(y,ytab,ny,ifind)

          END IF

          j1=j-1
          j2=j
          j3=j+1
          j4=j+2

          y1=ytab(j1)
          y2=ytab(j2)
          y3=ytab(j3)
          y4=ytab(j4)

          !

          f11=ftab(i1,j1)
          f12=ftab(i1,j2)
          f13=ftab(i1,j3)
          f14=ftab(i1,j4)

          f21=ftab(i2,j1)
          f22=ftab(i2,j2)
          f23=ftab(i2,j3)
          f24=ftab(i2,j4)

          f31=ftab(i3,j1)
          f32=ftab(i3,j2)
          f33=ftab(i3,j3)
          f34=ftab(i3,j4)

          f41=ftab(i4,j1)
          f42=ftab(i4,j2)
          f43=ftab(i4,j3)
          f44=ftab(i4,j4)

          ! x interpolation

          CALL fix_cubic(a,b,c,d,x1,f11,x2,f21,x3,f31,x4,f41)
          f01=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f12,x2,f22,x3,f32,x4,f42)
          f02=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f13,x2,f23,x3,f33,x4,f43)
          f03=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,x1,f14,x2,f24,x3,f34,x4,f44)
          f04=a*x**3+b*x**2+c*x+d

          CALL fix_cubic(a,b,c,d,y1,f01,y2,f02,y3,f03,y4,f04)
          findy=a*y**3+b*y**2+c*y+d

          ! y interpolation

          CALL fix_cubic(a,b,c,d,y1,f11,y2,f12,y3,f13,y4,f14)
          f10=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f21,y2,f22,y3,f23,y4,f24)
          f20=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f31,y2,f32,y3,f33,y4,f34)
          f30=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,y1,f41,y2,f42,y3,f43,y4,f44)
          f40=a*y**3+b*y**2+c*y+d

          CALL fix_cubic(a,b,c,d,x1,f10,x2,f20,x3,f30,x4,f40)
          findx=a*x**3+b*x**2+c*x+d

          ! Final result is an average over each direction
          find_2D=(findx+findy)/2.

       END IF

    ELSE

       STOP 'FIND_2D: order for interpolation not specified correctly'

    END IF

  END FUNCTION find_2D

  REAL FUNCTION find_3D(x,xin,y,yin,z,zin,fin,nx,ny,nz,iorder,ifind,imeth)

    ! A 3D interpolation routine to find value f(x,y,z) given a function evalated on arrays
    ! The linear version implemented here is also know as 'trilinear interpolation'
    ! TODO: Implement loops over coordinates to avoid repetition
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN) :: xin(nx)
    REAL, INTENT(IN) :: y
    REAL, INTENT(IN) :: yin(ny)
    REAL, INTENT(IN) :: z
    REAL, INTENT(IN) :: zin(nz)
    REAL, INTENT(IN) :: fin(nx,ny,nz)
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    INTEGER, INTENT(IN) :: nz
    INTEGER, INTENT(IN) :: iorder
    INTEGER, INTENT(IN) :: ifind
    INTEGER, INTENT(IN) :: imeth
    REAL :: x1, x2, y1, y2, z1, z2
    REAL :: F(2,2,2)
    REAL :: V(2,2,2) 
    INTEGER :: ix, ix1, ix2, iy, iy1, iy2, iz, iz1, iz2

    ! If the value required is off the table edge then this will halt

    ! iorder = 1 => linear interpolation

    ! ifind = 1 => find x in xtab by crudely searching from x(1) to x(n)
    ! ifind = 2 => find x in xtab quickly assuming the table is linearly spaced
    ! ifind = 3 => find x in xtab using midpoint splitting (iterations=ceiling(log2(n)))

    ! imeth = 1 => Uses cubic polynomials for interpolation

    IF(imeth .NE. 1) STOP 'FIND_3D: No Lagrange polynomials for you, only regular polynomails, imeth=1'

    IF(xin(1)>xin(nx)) STOP 'FIND_3D: x array in wrong order'
    IF(yin(1)>yin(ny)) STOP 'FIND_3D: y array in wrong order'
    IF(zin(1)>zin(nz)) STOP 'FIND_3D: z array in wrong order'

    ! No interpolation if the desired point is outside of the array range
    IF((x<xin(1) .OR. x>xin(nx)) .OR. (y>yin(ny) .OR. y<yin(1)) .OR. (z>zin(ny) .OR. z<zin(1))) THEN
       WRITE(*,*) 'FIND_3D: array xmin:', xin(1)
       WRITE(*,*) 'FIND_3D: array xmax:', xin(nx)
       WRITE(*,*) 'FIND_3D: requested x:', x
       WRITE(*,*) 'FIND_3D: array ymin:', yin(1)
       WRITE(*,*) 'FIND_3D: array ymax:', yin(ny)
       WRITE(*,*) 'FIND_3D: requested y:', y
       WRITE(*,*) 'FIND_3D: array zmin:', zin(1)
       WRITE(*,*) 'FIND_3D: array zmax:', zin(ny)
       WRITE(*,*) 'FIND_3D: requested z:', z
       STOP 'FIND_3D: Desired point is outside array range'
    END IF

    IF(iorder==1) THEN

       IF(nx<2) STOP 'FIND_3D: Not enough x points in your array for linear interpolation'
       IF(ny<2) STOP 'FIND_3D: Not enough y points in your array for linear interpolation'
       IF(nz<2) STOP 'FIND_3D: Not enough z points in your array for linear interpolation'

       !! Get the x,y,z values !!

       ! Get the integer coordinates in the x direction
       IF(x<=xin(2)) THEN
          ix=1
       ELSE IF (x>=xin(nx-1)) THEN
          ix=nx-1
       ELSE
          ix=select_table_integer(x,xin,nx,ifind)
       END IF
       ix1=ix
       ix2=ix+1

       ! Get the x values
       x1=xin(ix1)
       x2=xin(ix2)

       ! Get the integer coordinates in the y direction
       IF(y<=yin(2)) THEN
          iy=1
       ELSE IF (y>=yin(ny-1)) THEN
          iy=ny-1
       ELSE
          iy=select_table_integer(y,yin,ny,ifind)
       END IF
       iy1=iy
       iy2=iy+1

       ! Get the y values
       y1=yin(iy1)
       y2=yin(iy2)

       ! Get the integer coordinates in the z direction
       IF(z<=zin(2)) THEN
          iz=1
       ELSE IF (z>=zin(nz-1)) THEN
          iz=nz-1
       ELSE
          iz=select_table_integer(z,zin,nz,ifind)
       END IF
       iz1=iz
       iz2=iz+1

       ! Get the z values
       z1=zin(iz1)
       z2=zin(iz2)

       !! Function values !!

       f(1,1,1)=fin(ix1,iy1,iz1)
       f(1,1,2)=fin(ix1,iy1,iz2)
       f(1,2,1)=fin(ix1,iy2,iz1)
       f(1,2,2)=fin(ix1,iy2,iz2)
       f(2,1,1)=fin(ix2,iy1,iz1)
       f(2,1,2)=fin(ix2,iy1,iz2)
       f(2,2,1)=fin(ix2,iy2,iz1)
       f(2,2,2)=fin(ix2,iy2,iz2)

       !! !!

!!$       A11=(x2-x)*(y2-y)*F22
!!$       A12=(x2-x)*(y-y1)*F12
!!$       A21=(x-x1)*(y2-y)*F21
!!$       A22=(x-x1)*(y-y1)*F22
!!$
!!$       find_3D=(A11+A12+A21+A22)/((x2-x1)(y2-y1))

       V(1,1,1)=(x2-x)*(y2-y)*(z2-z)*F(1,1,1)
       V(1,2,1)=(x2-x)*(y-y1)*(z2-z)*F(1,2,1)
       V(2,1,1)=(x-x1)*(y2-y)*(z2-z)*F(2,1,1)
       V(2,2,1)=(x-x1)*(y-y1)*(z2-z)*F(2,2,1)
       V(1,1,2)=(x2-x)*(y2-y)*(z-z1)*F(1,1,2)
       V(1,2,2)=(x2-x)*(y-y1)*(z-z1)*F(1,2,2)
       V(2,1,2)=(x-x1)*(y2-y)*(z-z1)*F(2,1,2)
       V(2,2,2)=(x-x1)*(y-y1)*(z-z1)*F(2,2,2)

       find_3D=sum(V)/((x2-x1)*(y2-y1)*(z2-z1))

!!$       !! y direction interpolation
!!$
!!$       CALL fix_line(a,b,x1,f11,x2,f21)
!!$       f01=a*x+b
!!$
!!$       CALL fix_line(a,b,x1,f12,x2,f22)
!!$       f02=a*x+b
!!$
!!$       CALL fix_line(a,b,y1,f01,y2,f02)
!!$       findy=a*y+b
!!$
!!$       !!
!!$
!!$       !! x direction interpolation
!!$
!!$       CALL fix_line(a,b,y1,f11,y2,f12)
!!$       f10=a*y+b
!!$
!!$       CALL fix_line(a,b,y1,f21,y2,f22)
!!$       f20=a*y+b
!!$
!!$       CALL fix_line(a,b,x1,f10,x2,f20)
!!$       findx=a*x+b
!!$
!!$       !!
!!$
!!$       ! Final result is an average over each direction
!!$       find_3D=(findx+findy)/2.

    ELSE

       STOP 'FIND_3D: order for interpolation not specified correctly, only linear implemented'

    END IF

  END FUNCTION find_3D

  SUBROUTINE interpolate_array(x1,y1,n1,x2,y2,n2,iorder,ifind,imeth)

    ! Interpolates array 'x1-y1' onto new 'x' values x2 and output y2
    IMPLICIT NONE
    REAL, INTENT(IN) :: x1(n1), y1(n1), x2(n2)
    REAL, INTENT(OUT) :: y2(n2)
    INTEGER, INTENT(IN) :: n1, n2
    INTEGER, INTENT(IN) :: iorder, ifind, imeth
    INTEGER :: i

    ! Could be more efficient, but probably not worth the hassle: it does 'find integer' every time

    DO i=1,n2
       y2(i)=find(x2(i),x1,y1,n1,iorder,ifind,imeth)
    END DO

  END SUBROUTINE interpolate_array

END MODULE interpolate
