PROGRAM interpolate_test

  USE interpolate
  USE constants
  USE array_operations

  IMPLICIT NONE
  REAL :: xmin, xmax, x, ymax, ymin, y
  REAL, ALLOCATABLE :: xtab(:), ytab(:), f(:,:)
  REAL :: f_int1, f_int3, f_true
  REAL :: lin, quad, cube, tru
  REAL :: rlin, rquad, rcube
  INTEGER :: i, j, n, m, nx, nxtab, ny, nytab
  INTEGER :: iexample, imeth, itest

  WRITE(*,*)
  WRITE(*,*) 'Routines for interpolating from a table of data'
  WRITE(*,*)

  WRITE(*,*) '1 - 1D'
  WRITE(*,*) '2 - 2D'
  READ(*,*) itest

  IF(itest==1) THEN

     WRITE(*,*) 'Number of points for tables (low is a good test):'
     READ(*,*) n
     WRITE(*,*)

     WRITE(*,*) 'Method to test'
     WRITE(*,*) '1 - Standard polynomial fitting'
     WRITE(*,*) '2 - Lagrange polynomial fitting'
     READ(*,*) imeth
     WRITE(*,*)

     ALLOCATE(xtab(n),ytab(n))

     WRITE(*,*) '0 - Linear (0 to 3)'
     WRITE(*,*) '1 - Quadratic (0 to 3)'
     WRITE(*,*) '2 - Cubic (0 to 3)'
     WRITE(*,*) '3 - Sin (0 to pi)'
     WRITE(*,*) '4 - Exp (0 to 3)'
     READ(*,*) iexample
     WRITE(*,*)

     IF(iexample==0) THEN
        xmin=0.
        xmax=3.
     ELSE IF(iexample==1) THEN
        xmin=0.
        xmax=3.
     ELSE IF(iexample==2) THEN
        xmin=0.
        xmax=3.
     ELSE IF(iexample==3) THEN
        xmin=0.
        xmax=pi
     ELSE IF(iexample==4) THEN
        xmin=0.
        xmax=3.
     ELSE
        STOP 'Error, example not specified correctly'
     END IF

     DO i=1,n
        xtab(i)=xmin+(xmax-xmin)*float(i-1)/float(n-1)
     END DO

     IF(iexample==0) ytab=xtab
     IF(iexample==1) ytab=xtab**2
     IF(iexample==2) ytab=xtab**3
     IF(iexample==3) ytab=sin(xtab)
     IF(iexample==4) ytab=exp(xtab)

     OPEN(7,file='table.dat')
     DO i=1,n
        WRITE(7,*) xtab(i), ytab(i)
     END DO
     CLOSE(7)

     m=10*n

     WRITE(*,*) 'Writing tables'
     OPEN(7,file='results.dat')
     OPEN(8,file='ratio.dat')
     DO i=1,m
        x=xmin+(xmax-xmin)*float(i-1)/float(m-1)
        lin=find(x,xtab,ytab,n,1,2,imeth)
        quad=find(x,xtab,ytab,n,2,2,imeth)
        cube=find(x,xtab,ytab,n,3,2,imeth)
        tru=truth(iexample,x)
        WRITE(7,*) x, lin, quad, cube, tru
        IF(tru==0.) THEN
           rlin=1.
           rquad=1.
           rcube=1.
        ELSE
           rlin=lin/tru
           rquad=quad/tru
           rcube=cube/tru
        END IF
        WRITE(8,*) x, rlin, rquad, rcube
     END DO
     CLOSE(7)
     CLOSE(8)
     WRITE(*,*) 'Done'
     WRITE(*,*)

     DO

        WRITE(*,*) 'Value to *find* for tests (-1 exits):'
        READ(*,*) x

        IF(x==-1.) EXIT

        WRITE(*,*) 'Linear interpolation:', find(x,xtab,ytab,n,1,3,imeth)
        WRITE(*,*) 'Quadratic interpolation:', find(x,xtab,ytab,n,2,3,imeth)
        WRITE(*,*) 'Cubic interpolation:', find(x,xtab,ytab,n,3,3,imeth)
        WRITE(*,*) 'Truth:', truth(iexample,x)
        WRITE(*,*)

     END DO

     WRITE(*,*)

  ELSE IF(itest==2) THEN

     nxtab=8
     xmin=0.
     xmax=1.
     CALL fill_table(xmin,xmax,xtab,nxtab)

     nytab=8
     ymin=0.
     ymax=1.
     CALL fill_table(ymin,ymax,ytab,nytab)

     ALLOCATE(f(nxtab,nytab))

     DO i=1,nxtab
        DO j=1,nytab
           f(i,j)=func(xtab(i),ytab(j))
           WRITE(*,*) i, j, f(i,j)
        END DO
     END DO

     x=0.9
     y=0.3

     WRITE(*,*) 'x:', x
     WRITE(*,*) 'y:', y
     WRITE(*,*) 'Linear Interpolation:', find2d(x,xtab,y,ytab,f,nxtab,nytab,1,3,1)
     WRITE(*,*) 'Cubic Interpolation:', find2d(x,xtab,y,ytab,f,nxtab,nytab,3,3,1)
     WRITE(*,*) 'Truth:', func(x,y)
     WRITE(*,*)

     xmin=0.
     xmax=1.
     nx=21

     ymin=0.
     ymax=1.
     ny=31

     OPEN(7,file='results_linear.dat')
     OPEN(8,file='results_cubic.dat')
     DO i=1,nx
        DO j=1,ny

           x=xmin+(xmax-xmin)*float(i-1)/float(nx-1)
           y=ymin+(ymax-ymin)*float(j-1)/float(ny-1)

           f_true=func(x,y)
           f_int1=find2d(x,xtab,y,ytab,f,nxtab,nytab,1,3,1)
           f_int3=find2d(x,xtab,y,ytab,f,nxtab,nytab,3,3,1)

           WRITE(7,*) x, y, f_true, f_int1, f_int1/f_true
           WRITE(8,*) x, y, f_true, f_int3, f_int3/f_true

        END DO
     END DO
     CLOSE(7)
     CLOSE(8)

  END IF

CONTAINS

  FUNCTION truth(iexample,x)

    IMPLICIT NONE
    REAL :: x, truth
    INTEGER :: iexample

    IF(iexample==0) THEN
       truth=x
    ELSE IF(iexample==1) THEN
       truth=x**2
    ELSE IF(iexample==2) THEN
       truth=x**3
    ELSE IF(iexample==3) THEN
       truth=sin(x)
    ELSE IF(iexample==4) THEN
       truth=exp(x)
    ELSE
       STOP 'Example not specified correctly'
    END IF

  END FUNCTION truth

  FUNCTION func(x,y)

    IMPLICIT NONE
    REAL :: func
    REAL, INTENT(IN) :: x, y

    func=sin(x*y)+1.

  END FUNCTION func

END PROGRAM interpolate_test
