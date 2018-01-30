MODULE solve_equations

  USE interpolate
  
CONTAINS

  FUNCTION find_solve(a,xtab,ytab,n)

    !Solves y(x)=a for x
    IMPLICIT NONE
    REAL :: find_solve
    REAL, INTENT(IN) :: a, xtab(n), ytab(n)
    INTEGER, INTENT(IN) :: n

    find_solve=find(a,ytab,xtab,n,3,3,2)
    
  END FUNCTION find_solve

  FUNCTION bisect_solve(xtab,ytab,n,acc)

    USE logical_operations

    !Solves y(x)=0 for x, f(x) should be monotonic and cross f=0. once only
    IMPLICIT NONE
    REAL :: bisect_solve
    REAL, INTENT(IN) :: xtab(n), ytab(n)
    INTEGER, INTENT(IN) :: n
    REAL :: x1, x2, y1, y2, x, y
    REAL, INTENT(IN) :: acc
    INTEGER :: i

    !Initial values taken from top and bottom of table
    x1=xtab(1)
    x2=xtab(n)
    y1=ytab(1)
    y2=ytab(n)

    !Now iterate until desired accuracy is reached
    i=0
    DO
       i=i+1
       x=0.5*(x1+x2)
       y=find(x,xtab,ytab,n,3,3,2)
       WRITE(*,*) i, x
       IF(ABS(y)<acc) THEN
          EXIT
       ELSE IF(positive(y1) .EQV. positive(y)) THEN
          x1=x
          y1=y
       ELSE
          x2=x
          y2=y
       END IF
    END DO

    bisect_solve=x
    
  END FUNCTION bisect_solve

!!$  SUBROUTINE fill_table(min,max,arr,n)
!!$
!!$    !Fills array 'arr' in equally spaced intervals
!!$    !I'm not sure if inputting an array like this is okay
!!$    IMPLICIT NONE
!!$    INTEGER :: i
!!$    REAL, INTENT(IN) :: min, max
!!$    REAL, ALLOCATABLE :: arr(:)
!!$    INTEGER, INTENT(IN) :: n
!!$
!!$    !Allocate the array, and deallocate it if it is full
!!$    IF(ALLOCATED(arr)) DEALLOCATE(arr)
!!$    ALLOCATE(arr(n))
!!$    arr=0.
!!$
!!$    IF(n==1) THEN
!!$       arr(1)=min
!!$    ELSE IF(n>1) THEN
!!$       DO i=1,n
!!$          arr(i)=min+(max-min)*float(i-1)/float(n-1)
!!$       END DO
!!$    END IF
!!$
!!$  END SUBROUTINE fill_table
!!$
!!$  FUNCTION Lagrange_polynomial(x,n,xv,yv)
!!$
!!$    !Computes the result of the nth order Lagrange polynomial at point x, L(x)
!!$    IMPLICIT NONE
!!$    REAL :: Lagrange_polynomial
!!$    REAL, INTENT(IN) :: x, xv(n+1), yv(n+1)
!!$    REAL :: l(n+1)
!!$    INTEGER, INTENT(IN) :: n
!!$    INTEGER :: i, j
!!$
!!$    !Initialise variables, one for sum and one for multiplication
!!$    Lagrange_polynomial=0.
!!$    l=1.
!!$
!!$    !Loops to find the polynomials, one is a sum and one is a multiple
!!$    DO i=0,n
!!$       DO j=0,n
!!$          IF(i .NE. j) l(i+1)=l(i+1)*(x-xv(j+1))/(xv(i+1)-xv(j+1))
!!$       END DO
!!$       Lagrange_polynomial=Lagrange_polynomial+l(i+1)*yv(i+1)
!!$    END DO
!!$    
!!$  END FUNCTION Lagrange_polynomial
!!$
!!$  FUNCTION find(x,xin,yin,n,iorder,ifind,imeth)
!!$
!!$    !Given two arrays x and y this routine interpolates to find the y_i value at position x_i
!!$    IMPLICIT NONE
!!$    REAL :: find
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(IN) :: x, xin(n), yin(n)
!!$    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
!!$    REAL :: a, b, c, d
!!$    REAL :: x1, x2, x3, x4
!!$    REAL :: y1, y2, y3, y4
!!$    INTEGER :: i
!!$    INTEGER, INTENT(IN) :: imeth, iorder, ifind
!!$
!!$    !This version interpolates if the value is off either end of the array!
!!$    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
!!$    !Results from the interpolation!
!!$
!!$    !If the value required is off the table edge the interpolation is always linear
!!$
!!$    !iorder = 1 => linear interpolation
!!$    !iorder = 2 => quadratic interpolation
!!$    !iorder = 3 => cubic interpolation
!!$
!!$    !ifind = 1 => find x in xtab quickly assuming the table is linearly spaced
!!$    !ifind = 2 => find x in xtab by crudely searching from x(1) to x(n)
!!$    !ifind = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))
!!$
!!$    !imeth = 1 => Uses cubic polynomials for interpolation
!!$    !imeth = 2 => Uses Lagrange polynomials for interpolation
!!$
!!$    ALLOCATE(xtab(n),ytab(n))
!!$
!!$    xtab=xin
!!$    ytab=yin
!!$
!!$    IF(xtab(1)>xtab(n)) THEN
!!$       !Reverse the arrays in this case
!!$       CALL reverse(xtab,n)
!!$       CALL reverse(ytab,n)
!!$    END IF
!!$
!!$    IF(x<xtab(1)) THEN
!!$
!!$       !Do a linear interpolation beyond the table boundary
!!$
!!$       x1=xtab(1)
!!$       x2=xtab(2)
!!$
!!$       y1=ytab(1)
!!$       y2=ytab(2)
!!$
!!$       IF(imeth==1) THEN
!!$          CALL fit_line(a,b,x1,y1,x2,y2)
!!$          find=a*x+b
!!$       ELSE IF(imeth==2) THEN
!!$          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
!!$       ELSE
!!$          STOP 'FIND: Error, method not specified correctly'
!!$       END IF
!!$       
!!$    ELSE IF(x>xtab(n)) THEN
!!$
!!$       !Do a linear interpolation beyond the table boundary
!!$       
!!$       x1=xtab(n-1)
!!$       x2=xtab(n)
!!$
!!$       y1=ytab(n-1)
!!$       y2=ytab(n)
!!$
!!$       IF(imeth==1) THEN
!!$          CALL fit_line(a,b,x1,y1,x2,y2)
!!$          find=a*x+b
!!$       ELSE IF(imeth==2) THEN
!!$          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
!!$       ELSE
!!$          STOP 'FIND: Error, method not specified correctly'
!!$       END IF
!!$
!!$    ELSE IF(iorder==1) THEN
!!$
!!$       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'
!!$
!!$       IF(x<=xtab(2)) THEN
!!$
!!$          x1=xtab(1)
!!$          x2=xtab(2)
!!$
!!$          y1=ytab(1)
!!$          y2=ytab(2)
!!$
!!$       ELSE IF (x>=xtab(n-1)) THEN
!!$
!!$          x1=xtab(n-1)
!!$          x2=xtab(n)
!!$
!!$          y1=ytab(n-1)
!!$          y2=ytab(n)
!!$
!!$       ELSE
!!$
!!$          i=table_integer(x,xtab,n,ifind)
!!$          
!!$          x1=xtab(i)
!!$          x2=xtab(i+1)
!!$
!!$          y1=ytab(i)
!!$          y2=ytab(i+1)
!!$
!!$       END IF
!!$
!!$       IF(imeth==1) THEN
!!$          CALL fit_line(a,b,x1,y1,x2,y2)
!!$          find=a*x+b
!!$       ELSE IF(imeth==2) THEN
!!$          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
!!$       ELSE
!!$          STOP 'FIND: Error, method not specified correctly'
!!$       END IF
!!$
!!$    ELSE IF(iorder==2) THEN
!!$
!!$       IF(n<3) STOP 'FIND: Not enough points in your table'
!!$
!!$       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN
!!$
!!$          IF(x<=xtab(2)) THEN
!!$
!!$             x1=xtab(1)
!!$             x2=xtab(2)
!!$             x3=xtab(3)
!!$
!!$             y1=ytab(1)
!!$             y2=ytab(2)
!!$             y3=ytab(3)
!!$
!!$          ELSE IF (x>=xtab(n-1)) THEN
!!$
!!$             x1=xtab(n-2)
!!$             x2=xtab(n-1)
!!$             x3=xtab(n)
!!$
!!$             y1=ytab(n-2)
!!$             y2=ytab(n-1)
!!$             y3=ytab(n)
!!$
!!$          END IF
!!$
!!$          IF(imeth==1) THEN
!!$             CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
!!$             find=a*(x**2.)+b*x+c
!!$          ELSE IF(imeth==2) THEN
!!$             find=Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))
!!$          ELSE
!!$             STOP 'FIND: Error, method not specified correctly'
!!$          END IF
!!$             
!!$       ELSE
!!$
!!$          i=table_integer(x,xtab,n,ifind)
!!$
!!$          x1=xtab(i-1)
!!$          x2=xtab(i)
!!$          x3=xtab(i+1)
!!$          x4=xtab(i+2)
!!$
!!$          y1=ytab(i-1)
!!$          y2=ytab(i)
!!$          y3=ytab(i+1)
!!$          y4=ytab(i+2)
!!$
!!$          IF(imeth==1) THEN
!!$             !In this case take the average of two separate quadratic spline values
!!$             CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
!!$             find=(a*(x**2.)+b*x+c)/2.
!!$             CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
!!$             find=find+(a*(x**2.)+b*x+c)/2.
!!$          ELSE IF(imeth==2) THEN
!!$             !In this case take the average of two quadratic Lagrange polynomials
!!$             find=(Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))+Lagrange_polynomial(x,2,(/x2,x3,x4/),(/y2,y3,y4/)))/2.
!!$          ELSE
!!$             STOP 'FIND: Error, method not specified correctly'
!!$          END IF
!!$
!!$       END IF
!!$
!!$    ELSE IF(iorder==3) THEN
!!$
!!$       IF(n<4) STOP 'FIND: Not enough points in your table'
!!$
!!$       IF(x<=xtab(3)) THEN
!!$
!!$          x1=xtab(1)
!!$          x2=xtab(2)
!!$          x3=xtab(3)
!!$          x4=xtab(4)        
!!$
!!$          y1=ytab(1)
!!$          y2=ytab(2)
!!$          y3=ytab(3)
!!$          y4=ytab(4)
!!$
!!$       ELSE IF (x>=xtab(n-2)) THEN
!!$
!!$          x1=xtab(n-3)
!!$          x2=xtab(n-2)
!!$          x3=xtab(n-1)
!!$          x4=xtab(n)
!!$
!!$          y1=ytab(n-3)
!!$          y2=ytab(n-2)
!!$          y3=ytab(n-1)
!!$          y4=ytab(n)
!!$
!!$       ELSE
!!$
!!$          i=table_integer(x,xtab,n,ifind)
!!$
!!$          x1=xtab(i-1)
!!$          x2=xtab(i)
!!$          x3=xtab(i+1)
!!$          x4=xtab(i+2)
!!$
!!$          y1=ytab(i-1)
!!$          y2=ytab(i)
!!$          y3=ytab(i+1)
!!$          y4=ytab(i+2)
!!$
!!$       END IF
!!$
!!$       IF(imeth==1) THEN
!!$          CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
!!$          find=a*x**3.+b*x**2.+c*x+d
!!$       ELSE IF(imeth==2) THEN
!!$          find=Lagrange_polynomial(x,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
!!$       ELSE
!!$          STOP 'FIND: Error, method not specified correctly'
!!$       END IF
!!$
!!$    ELSE
!!$
!!$       STOP 'FIND: Error, interpolation order specified incorrectly'
!!$
!!$    END IF
!!$
!!$  END FUNCTION find
!!$
!!$  FUNCTION table_integer(x,xtab,n,imeth)
!!$
!!$    !Chooses between ways to find the integer location below some value in an array
!!$    IMPLICIT NONE
!!$    INTEGER :: table_integer
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(IN) :: x, xtab(n)
!!$    INTEGER, INTENT(IN) :: imeth
!!$
!!$    IF(imeth==1) THEN
!!$       table_integer=linear_table_integer(x,xtab,n)
!!$    ELSE IF(imeth==2) THEN
!!$       table_integer=search_int(x,xtab,n)
!!$    ELSE IF(imeth==3) THEN
!!$       table_integer=int_split(x,xtab,n)
!!$    ELSE
!!$       STOP 'TABLE INTEGER: Method specified incorrectly'
!!$    END IF
!!$
!!$  END FUNCTION table_integer
!!$
!!$  FUNCTION linear_table_integer(x,xtab,n)
!!$
!!$    !Assuming the table is exactly linear this gives you the integer position
!!$    IMPLICIT NONE
!!$    INTEGER :: linear_table_integer
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(IN) :: x, xtab(n)
!!$    REAL :: x1, x2, xn
!!$    REAL :: acc
!!$
!!$    !Returns the integer (table position) below the value of x
!!$    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 6
!!$    !Assumes table is organised linearly (care for logs)
!!$
!!$    !n=SIZE(xtab)
!!$    x1=xtab(1)
!!$    x2=xtab(2)
!!$    xn=xtab(n)
!!$
!!$    !Test for linear table
!!$    acc=0.001
!!$
!!$    IF(x1>xn) STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
!!$    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'
!!$
!!$    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))
!!$
!!$  END FUNCTION linear_table_integer
!!$
!!$  FUNCTION search_int(x,xtab,n)
!!$
!!$    !Does a stupid search through the table from beginning to end to find integer
!!$    IMPLICIT NONE
!!$    INTEGER :: search_int
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(IN) :: x, xtab(n)
!!$    INTEGER :: i
!!$
!!$    IF(xtab(1)>xtab(n)) STOP 'SEARCH_INT: table in wrong order'
!!$
!!$    DO i=1,n
!!$       IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
!!$    END DO
!!$
!!$    search_int=i
!!$
!!$  END FUNCTION search_int
!!$
!!$  FUNCTION int_split(x,xtab,n)
!!$
!!$    !Finds the position of the value in the table by continually splitting it in half
!!$    IMPLICIT NONE
!!$    INTEGER :: int_split
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(IN) :: x, xtab(n)
!!$    INTEGER :: i1, i2, imid
!!$
!!$    IF(xtab(1)>xtab(n)) STOP 'INT_SPLIT: table in wrong order'
!!$
!!$    i1=1
!!$    i2=n
!!$
!!$    DO
!!$       
!!$       imid=NINT((i1+i2)/2.)
!!$
!!$       IF(x<xtab(imid)) THEN
!!$          i2=imid
!!$       ELSE
!!$          i1=imid
!!$       END IF
!!$
!!$       IF(i2==i1+1) EXIT
!!$
!!$    END DO
!!$    
!!$    int_split=i1
!!$
!!$  END FUNCTION int_split
!!$
!!$  SUBROUTINE fit_line(a1,a0,x1,y1,x2,y2)
!!$
!!$    !Given xi, yi i=1,2 fits a line between these points
!!$    IMPLICIT NONE
!!$    REAL, INTENT(OUT) :: a0, a1
!!$    REAL, INTENT(IN) :: x1, y1, x2, y2   
!!$
!!$    a1=(y2-y1)/(x2-x1)
!!$    a0=y1-a1*x1
!!$
!!$  END SUBROUTINE fit_line
!!$
!!$  SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)
!!$
!!$    !Given xi, yi i=1,2,3 fits a quadratic between these points
!!$    IMPLICIT NONE
!!$    REAL, INTENT(OUT) :: a0, a1, a2
!!$    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3   
!!$
!!$    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
!!$    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
!!$    a0=y1-a2*(x1**2.)-a1*x1
!!$
!!$  END SUBROUTINE fit_quadratic
!!$
!!$  SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
!!$
!!$    !Given xi, yi i=1,2,3,4 fits a cubic between these points
!!$    IMPLICIT NONE
!!$    REAL, INTENT(OUT) :: a, b, c, d
!!$    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
!!$    REAL :: f1, f2, f3    
!!$
!!$    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
!!$    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
!!$    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))
!!$
!!$    a=f1-f2-f3
!!$
!!$    f1=(y3-y1)/((x3-x2)*(x3-x1))
!!$    f2=(y2-y1)/((x2-x1)*(x3-x2))
!!$    f3=a*(x3+x2+x1)
!!$
!!$    b=f1-f2-f3
!!$
!!$    f1=(y4-y1)/(x4-x1)
!!$    f2=a*(x4**2.+x4*x1+x1**2.)
!!$    f3=b*(x4+x1)
!!$
!!$    c=f1-f2-f3
!!$
!!$    d=y1-a*x1**3.-b*x1**2.-c*x1
!!$
!!$  END SUBROUTINE fit_cubic
!!$
!!$  SUBROUTINE reverse(arry,n)
!!$
!!$    !This reverses the contents of arry!
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(INOUT) :: arry(n)
!!$    INTEGER :: i
!!$    REAL, ALLOCATABLE :: hold(:) 
!!$
!!$    ALLOCATE(hold(n))
!!$
!!$    hold=arry
!!$
!!$    DO i=1,n
!!$       arry(i)=hold(n-i+1)
!!$    END DO
!!$
!!$    DEALLOCATE(hold)
!!$
!!$  END SUBROUTINE reverse

END MODULE solve_equations
