PROGRAM fitting_test

  USE fitting

  IMPLICIT NONE
  REAL :: df, min, min2(2)
  REAL :: k(100), d2(100), d2mod(100)
  INTEGER :: i
  INTEGER :: itest
  REAL :: x1, x2, x3, f1, f2, f3
  REAL :: x

  !Accuracy parameter
  REAL, PARAMETER :: acc=1e-5

  !Initial guess
  REAL, PARAMETER :: x0=3.

  !One-dimensional gradient fit

!!$  IF(a<0.) STOP 'Error, quadratic function has no minimum'
!!$
!!$  WRITE(*,*)
!!$  WRITE(*,*) 'Finding minimum of quadratic'
!!$  WRITE(*,*) 'a:', a
!!$  WRITE(*,*) 'b:', b
!!$  WRITE(*,*) 'c:', c
!!$  WRITE(*,*)
!!$
!!$  WRITE(*,*) 'Analytical minimum:', quadratic_extrema()

  WRITE(*,*) '1 - Test gradient fit 1'
  WRITE(*,*) '2 - Test gradient fit 2'
  WRITE(*,*) '3 - Test grid fit 1'
  WRITE(*,*) '4 - Test grid fit 2'
  WRITE(*,*) '5 - Test grid fit 3'
  WRITE(*,*) '6 - Test quadratic fit 1'
  READ(*,*) itest
  WRITE(*,*)

  IF(itest==1) THEN

     WRITE(*,*) 'Initial guess:', x0

     df=derivative(quadratic,x0,acc)
     WRITE(*,*) 'Derivative at initial guess', df

     min=gradient_fit_1(quadratic,x0,acc)

     WRITE(*,*) 'Gradient minimum:', min
     !WRITE(*,*) 'Error:', min/quadratic_extrema()
     WRITE(*,*)

  ELSE IF(itest==2) THEN

     min2=gradient_fit_2(test,10.,3.,acc)
     WRITE(*,*) min2(1)
     WRITE(*,*) min2(2)

  ELSE IF(itest==3) THEN


     DO i=1,100
        k(i)=float(i)
        d2(i)=k(i)**3.2
     END DO

     CALL grid_fit_1(0.,5.,501,k,d2)

  ELSE IF(itest==4) THEN

!!$     DO i=1,100
!!$        k(i)=float(i)
!!$        d2(i)=k(i)**3.+k(i)**0.5
!!$     END DO
!!$
!!$     CALL fit2(0.,5.,501,0.,5.,501,k,d2)

  ELSE IF(itest==5) THEN

!!$     DO i=1,100
!!$        k(i)=float(i)
!!$        d(i)=k(i)**2.+sin(3.*k(i))+5.*k(i)
!!$     END DO
!!$
!!$     CALL fit(0.,10.,101,0.,10.,101,0.,10.,101,k,d)

  ELSE IF(itest==6) THEN

     x1=1.
     x2=2.
     x3=3.

     f1=quadratic(x1)
     f2=quadratic(x2)
     f3=quadratic(x3)

     x=quadratic_fit_1(x1,f1,x2,f2,x3,f3)

     WRITE(*,*) 'Minimum:', x
     WRITE(*,*)

  END IF

CONTAINS

  FUNCTION quadratic(x)

    IMPLICIT NONE
    REAL :: quadratic
    REAL, INTENT(IN) :: x

    REAL, PARAMETER :: a=1.
    REAL, PARAMETER :: b=0.
    REAL, PARAMETER :: c=1.

    quadratic=a*x**2+b*x+c

  END FUNCTION quadratic

!!$  FUNCTION quadratic_extrema()
!!$
!!$    IMPLICIT NONE
!!$    REAL :: quadratic_extrema
!!$
!!$    quadratic_extrema=-b/(2.*a)
!!$
!!$  END FUNCTION quadratic_extrema

  FUNCTION test(x,y)

    IMPLICIT NONE
    REAL :: test
    REAL, INTENT(IN) :: x, y
    REAL :: A, B, C, D, E, F

    A=1.
    B=1.
    F=2.

    test=A*(x**2.)+B*(y**2.)+C*(x*y)+D*x+E*y+F

  END FUNCTION test

END PROGRAM fitting_test
