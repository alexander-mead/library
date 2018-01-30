PROGRAM calculus_table_test

  USE calculus_table
  USE constants

  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:), y(:)
  REAL :: sum1, sum2, sum3
  REAL :: xmin, xmax
  REAL :: truth
  INTEGER :: i, n, m
  INTEGER :: iexample, itest
  REAL :: xv, lin, quad, cube, tru
  REAL, ALLOCATABLE :: xtab(:), ytab(:)

  WRITE(*,*) '1 - Test differentiation'
  WRITE(*,*) '2 - Test integration'
  READ(*,*) itest
  WRITE(*,*)

  IF(itest==1) THEN

     WRITE(*,*)
     WRITE(*,*) 'Routines for taking derivatives from a table of data'
     WRITE(*,*)

     WRITE(*,*) 'Number of points for tables (low is a good test):'
     READ(*,*) n
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
        xmin=0.01
        xmax=3.
     ELSE IF(iexample==1) THEN
        xmin=0.01
        xmax=3.
     ELSE IF(iexample==2) THEN
        xmin=0.01
        xmax=3.
     ELSE IF(iexample==3) THEN
        xmin=0.
        xmax=pi
     ELSE IF(iexample==4) THEN
        xmin=0.
        xmax=3.
     ELSE
        STOP 'Incorrect example chosen'
     END IF

     DO i=1,n
        xtab(i)=xmin+(xmax-xmin)*float(i-1)/float(n-1)
     END DO

     IF(iexample==0) ytab=xtab
     IF(iexample==1) ytab=xtab**2.
     IF(iexample==2) ytab=xtab**3.
     IF(iexample==3) ytab=sin(xtab)
     IF(iexample==4) ytab=exp(xtab)

     m=10*n

     OPEN(7,file='results_linear.dat')
     OPEN(8,file='results_quadratic.dat')
     OPEN(9,file='results_cubic.dat')
     DO i=1,m
        xv=xmin+(xmax-xmin)*float(i-1)/float(m-1)
        lin=derivative_table(xv,xtab,ytab,n,1,3)
        quad=derivative_table(xv,xtab,ytab,n,2,3)
        cube=derivative_table(xv,xtab,ytab,n,3,3)
        tru=truth_deriv(iexample,xv)
        WRITE(7,*) xv, tru, lin, lin/tru
        WRITE(8,*) xv, tru, quad, quad/tru
        WRITE(9,*) xv, tru, cube, cube/tru
     END DO
     CLOSE(7)
     CLOSE(8)
     CLOSE(9)

  ELSE IF(itest==2) THEN

     !This integrates a 1D table of data (x,y)

     WRITE(*,*)
     WRITE(*,*) 'Routines for integrating a table of data'
     WRITE(*,*)

     WRITE(*,*) 'Number of points for tables:'
     READ(*,*) n
     WRITE(*,*)

     ALLOCATE(x(n),y(n))

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
     END IF

     DO i=1,n
        x(i)=xmin+(xmax-xmin)*float(i-1)/float(n-1)
     END DO

     IF(iexample==0) y=x
     IF(iexample==1) y=x**2.
     IF(iexample==2) y=x**3.
     IF(iexample==3) y=sin(x)
     IF(iexample==4) y=exp(x)

     sum1=integrate_table(x,y,n,1,n,1)
     sum2=integrate_table(x,y,n,1,n,2)
     sum3=integrate_table(x,y,n,1,n,3)

     IF(iexample==0) truth=truth_x(xmin,xmax)
     IF(iexample==1) truth=truth_x2(xmin,xmax)
     IF(iexample==2) truth=truth_x3(xmin,xmax)
     IF(iexample==3) truth=truth_sin(xmin,xmax)
     IF(iexample==4) truth=truth_exp(xmin,xmax)

     WRITE(*,*) 'Integral (Trapezium):', sum1
     WRITE(*,*) 'Integral (Quadratic):', sum2
     WRITE(*,*) 'Integral (Cubic):', sum3
     WRITE(*,*) 'Truth:', truth
     WRITE(*,*)

  END IF

CONTAINS

  FUNCTION truth_deriv(iexample,x)

    IMPLICIT NONE
    REAL :: truth_deriv
    REAL, INTENT(IN) :: x
    INTEGER :: iexample

    IF(iexample==0) THEN
       truth_deriv=1.
    ELSE IF(iexample==1) THEN
       truth_deriv=2.*x
    ELSE IF(iexample==2) THEN
       truth_deriv=3.*(x**2.)
    ELSE IF(iexample==3) THEN
       truth_deriv=cos(x)
    ELSE IF(iexample==4) THEN
       truth_deriv=exp(x)
    ELSE

    END IF

  END FUNCTION truth_deriv

  FUNCTION truth_x(xmin,xmax)

    IMPLICIT NONE
    REAL :: truth_x, xmin, xmax

    truth_x=(xmax**2.-xmin**2.)/2.

  END FUNCTION truth_x

  FUNCTION truth_x2(xmin,xmax)

    IMPLICIT NONE
    REAL :: truth_x2, xmin, xmax

    truth_x2=(xmax**3.-xmin**3.)/3.

  END FUNCTION truth_x2

  FUNCTION truth_x3(xmin,xmax)

    IMPLICIT NONE
    REAL :: truth_x3, xmin, xmax

    truth_x3=(xmax**4.-xmin**4.)/4.

  END FUNCTION truth_x3

  FUNCTION truth_sin(xmin,xmax)

    IMPLICIT NONE
    REAL :: truth_sin, xmin, xmax

    truth_sin=-cos(xmax)+cos(xmin)

  END FUNCTION truth_sin

  FUNCTION truth_exp(xmin,xmax)

    IMPLICIT NONE
    REAL :: truth_exp, xmin, xmax

    truth_exp=exp(xmax)-exp(xmin)

  END FUNCTION truth_exp

END PROGRAM calculus_table_test
