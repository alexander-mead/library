PROGRAM fix_polynomial_test

  USE fix_polynomial
  
  IMPLICIT NONE
  REAL :: x1, x2, x3, x4, y1, y2, y3, y4
  REAL :: x, f, L
  INTEGER :: i, n
  REAL :: a, b, c, d
  REAL :: f1, f2, f3, f4, L1, L2, L3, L4
  INTEGER :: imode

  WRITE(*,*)
  WRITE(*,*) 'Cubic fit tester'
  WRITE(*,*) '=========='
  WRITE(*,*)
  WRITE(*,*) '1 - Line'
  WRITE(*,*) '2 - Slight curve'
  WRITE(*,*) '3 - Curve'
  WRITE(*,*) '4 - Crazy curve'
  WRITE(*,*) '5 - Crazier curve'
  WRITE(*,*) '6 - Craziest curve'
  READ(*,*) imode
  WRITE(*,*)

  x1=1.
  x2=2.
  x3=3.
  x4=4.

  IF(imode==1) THEN

     y1=1.
     y2=2.
     y3=3.
     y4=4.

  ELSE IF(imode==2) THEN

     y1=1.
     y2=1.8
     y3=3.2
     y4=4.

  ELSE IF(imode==3) THEN

     y1=1.
     y2=1.5
     y3=3.5
     y4=4.

  ELSE IF(imode==4) THEN

     y1=1.
     y2=1.
     y3=4.
     y4=4.

  ELSE IF(imode==5) THEN

     y1=1.
     y2=3.
     y3=2.
     y4=4.

  ELSE IF(imode==6) THEN

     y1=1.
     y2=1.2
     y3=0.6
     y4=4.

  ELSE

     STOP 'Curve not specified correctly'

  END IF

  CALL fix_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

  OPEN(7,file='points.dat')
  WRITE(7,*) x1, y1
  WRITE(7,*) x2, y2
  WRITE(7,*) x3, y3
  WRITE(7,*) x4, y4
  CLOSE(7)

  n=100
  OPEN(7,file='cubic.dat')
  DO i=1,n
     x=x1+(x4-x1)*float(i-1)/float(n-1)
     f=a*(x**3.)+b*(x**2.)+c*x+d
     L=Lagrange_polynomial(x,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
     WRITE(7,*) x, f, L
  END DO
  CLOSE(7)

  f1=a*(x1**3.)+b*(x1**2.)+c*x1+d
  f2=a*(x2**3.)+b*(x2**2.)+c*x2+d
  f3=a*(x3**3.)+b*(x3**2.)+c*x3+d
  f4=a*(x4**3.)+b*(x4**2.)+c*x4+d

  L1=Lagrange_polynomial(x1,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
  L2=Lagrange_polynomial(x2,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
  L3=Lagrange_polynomial(x3,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
  L4=Lagrange_polynomial(x4,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))

  WRITE(*,*) '======================================='
  WRITE(*,*) '      x_i       y_i      y(x)      L(x)'
  WRITE(*,*) '======================================='
  WRITE(*,fmt='(4F10.5)') x1, y1, f1, L1
  WRITE(*,fmt='(4F10.5)') x2, y2, f2, L2
  WRITE(*,fmt='(4F10.5)') x3, y3, f3, L3
  WRITE(*,fmt='(4F10.5)') x4, y4, f4, L4
  WRITE(*,*) '======================================='

END PROGRAM fix_polynomial_test
