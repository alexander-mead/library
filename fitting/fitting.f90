MODULE fitting

  USE calculus

  IMPLICIT NONE

CONTAINS

  SUBROUTINE grid_fit_1(Amin,Amax,Asteps,x,y)

    IMPLICIT NONE
    REAL, INTENT(IN) :: Amin, Amax
    REAL, INTENT(IN) :: x(:), y(:)
    INTEGER, INTENT(IN) :: Asteps
    REAL :: ymod, fom, fombest, A, Abest
    INTEGER :: i, j

    !This fits a one parameter model

    fombest=10000000.

    DO j=1,Asteps

       A=Amin+(Amax-Amin)*(float(j-1)/(Asteps-1))

       fom=0.

       DO i=1,SIZE(x)

          !This is the model that you wish to fit
          ymod=x(i)**A

          fom=fom+(ymod-y(i))**2.

       END DO

       IF(fom<fombest) THEN
          fombest=fom
          Abest=A
       END IF

    END DO

    WRITE(*,*) 'Best fit is:', Abest

  END SUBROUTINE grid_fit_1

!!$  SUBROUTINE grid_fit_2(Amin,Amax,Asteps,Bmin,Bmax,Bsteps,x,y)
!!$
!!$    IMPLICIT NONE
!!$
!!$    REAL, INTENT(IN) :: Amin, Amax, Bmin, Bmax
!!$    REAL, INTENT(IN) :: x(:), y(:)
!!$    INTEGER, INTENT(IN) :: Asteps, Bsteps
!!$    REAL :: fom, fombest, A, Abest, B, Bbest, ymod
!!$    INTEGER :: i, j, k
!!$
!!$    !This fits a two parameter model
!!$
!!$    fombest=10000000.
!!$
!!$    DO j=1,Asteps
!!$
!!$       A=Amin+(Amax-Amin)*(float(j-1)/(float(Asteps-1)))
!!$
!!$       DO k=1,Bsteps
!!$
!!$          B=Bmin+(Bmax-Bmin)*(float(k-1)/(float(Bsteps-1)))
!!$
!!$          fom=0.
!!$          DO i=1,SIZE(x)
!!$
!!$             !This is the model that you wish to fit
!!$             ymod=model(A,B,x(i))
!!$
!!$             !This looks at the raw difference between models!
!!$             fom=fom+(-1.+ymod/y(i))**2.
!!$
!!$             !         WRITE(*,*) y(i), A, B, ymod, fom
!!$
!!$          END DO
!!$
!!$          IF(fom<fombest) THEN
!!$             fombest=fom
!!$             Abest=A
!!$             Bbest=B
!!$          END IF
!!$
!!$       END DO
!!$    END DO
!!$
!!$    WRITE(*,*) 'Best fit A is:', Abest
!!$    WRITE(*,*) 'Best fit B is:', Bbest
!!$
!!$    OPEN(8,file='results.dat')
!!$    DO i=1,SIZE(x)
!!$       WRITE(8,*) x(i), y(i), model(Abest,Bbest,x(i))
!!$    END DO
!!$    CLOSE(8)
!!$
!!$  END SUBROUTINE grid_fit_2

!!$  SUBROUTINE grid_fit_3(Amin,Amax,Asteps,Bmin,Bmax,Bsteps,Cmin,Cmax,Csteps,x,y)
!!$
!!$    IMPLICIT NONE
!!$
!!$    REAL, INTENT(IN) :: Amin, Amax, Bmin, Bmax, Cmin, Cmax
!!$    REAL, INTENT(IN) :: x(:), y(:)
!!$    INTEGER, INTENT(IN) :: Asteps, Bsteps, Csteps
!!$    REAL :: fom, fombest, A, Abest, B, Bbest, C, Cbest, ymod
!!$    INTEGER :: i, j, k, l
!!$
!!$    !This fits a two parameter model
!!$
!!$    fombest=10000000.
!!$
!!$    DO j=1,Asteps
!!$
!!$       IF(Asteps==1) THEN
!!$          A=Amin
!!$       ELSE
!!$          A=Amin+(Amax-Amin)*(float(j-1)/(Asteps-1))
!!$       END IF
!!$
!!$       DO k=1,Bsteps
!!$
!!$          IF(Bsteps==1) THEN
!!$             B=Bmin
!!$          ELSE
!!$             B=Bmin+(Bmax-Bmin)*(float(k-1)/(Bsteps-1))
!!$          END IF
!!$
!!$          DO l=1,Csteps
!!$
!!$             IF(Csteps==1) THEN
!!$                C=Cmin
!!$             ELSE
!!$                C=Cmin+(Cmax-Cmin)*(float(l-1)/(Csteps-1))
!!$             END IF
!!$
!!$             fom=0.
!!$             DO i=1,SIZE(x)
!!$
!!$                !This is the model that you wish to fit
!!$                ymod=model(A,B,C,x(i))
!!$
!!$                !This looks at the raw difference between models!
!!$                fom=fom+(-1.+ymod/y(i))**2.
!!$
!!$                !         WRITE(*,*) y(i), A, B, ymod, fom
!!$
!!$             END DO
!!$
!!$             !         WRITE(10,*) A, B, C, fom
!!$
!!$             IF(fom<fombest) THEN
!!$                fombest=fom
!!$                Abest=A
!!$                Bbest=B
!!$                Cbest=C
!!$             END IF
!!$
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    WRITE(*,*) 'Best fit A is:', Abest
!!$    WRITE(*,*) 'Best fit B is:', Bbest
!!$    WRITE(*,*) 'Best fit C is:', Cbest
!!$
!!$    OPEN(8,file='results.dat')
!!$    DO i=1,SIZE(x)
!!$       WRITE(8,*) x(i), y(i), model(Abest,Bbest,Cbest,x(i))
!!$    END DO
!!$    CLOSE(8)
!!$
!!$  END SUBROUTINE grid_fit_3

  FUNCTION gradient_fit_1(f,x,acc)

    !Finds the minimum of the function 'f' starting at value x0 with accuracy acc!

    IMPLICIT NONE
    REAL :: gradient_fit_1
    REAL, INTENT(IN) :: x, acc
    REAL :: xold, xnew, df
    INTEGER :: i

    INTEGER, PARAMETER :: n=1000
    REAL, PARAMETER :: fac=1e3

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    STOP 'GRADIENT_FIT_1: This does not work!'

    xold=x

    DO i=1,n

       !       df=derivative(f,xold,acc/100.)
       !       xnew=xold-df*acc*100.

       WRITE(*,*) i, xold, xnew

       df=derivative(f,xold,acc/fac)
       xnew=xold-df*acc*fac

       WRITE(*,*) i, xold, xnew

       IF(i>1 .AND. ABS(xnew/xold-1.)<acc) THEN
          gradient_fit_1=xnew
          EXIT
       ELSE
          xold=xnew
       END IF

    END DO

  END FUNCTION gradient_fit_1

  FUNCTION gradient_fit_2(f,x,y,acc)

    !Finds the minimum of the function 'f' starting at value x0 with accuracy acc!

    IMPLICIT NONE
    REAL :: gradient_fit_2(2)
    REAL, INTENT(IN) :: x, y, acc
    REAL :: xold, xnew, yold, ynew, dfx, dfy
    INTEGER :: i

    INTERFACE
       FUNCTION f(x,y)
         REAL :: f
         REAL, INTENT(IN) :: x,y
       END FUNCTION f
    END INTERFACE

    STOP 'GRADIENT_FIT_2: This does not work!'

    xold=x
    yold=y

    DO i=1,1000

       dfx=derivative_x(f,xold,yold,acc/100.)
       dfy=derivative_y(f,xold,yold,acc/100.)

       xnew=xold-dfx*acc*100.
       ynew=yold-dfy*acc*100.

       IF(i>1 .AND. ABS(xnew/xold-1.)<acc .AND. ABS(ynew/yold-1.)<acc) THEN
          gradient_fit_2(1)=xnew
          gradient_fit_2(2)=ynew
          EXIT
       ELSE
          xold=xnew
          yold=ynew
       END IF

    END DO

  END FUNCTION gradient_fit_2

  FUNCTION quadratic_fit_1(x1,y1,x2,y2,x3,y3)

    !This is not really fitting, it is finding a minimum, I should change the name
    IMPLICIT NONE
    REAL :: quadratic_fit_1
    REAL, INTENT(IN) :: x1, x2, x3, y1, y2, y3
    REAL :: a, b

    !This calculates the extrema of a function under the assumption 
    !that it is quadratic. Based on taking 3 points!

    b=(y2-y1)*(x3+x1)/((x2-x1)*(x3-x2))-(y3-y1)*(x2+x1)/((x3-x1)*(x3-x2))
    a=((y2-y1)-b*(x2-x1))/((x2-x1)*(x2+x1))

    quadratic_fit_1=-b/(2.*a)

  END FUNCTION quadratic_fit_1

END MODULE fitting
