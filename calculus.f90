MODULE calculus

CONTAINS

  FUNCTION derivative(f,x,acc)

    !Calculates the derivative of a function 'f' at the point x to accuracy acc!
    IMPLICIT NONE
    REAL :: derivative
    REAL, INTENT(IN) :: x, acc
    REAL :: dnew, dold, dx 
    INTEGER :: i

    INTEGER, PARAMETER :: imin=5 !Minimum number of iterations
    INTEGER, PARAMETER :: n=100 !Maximum number of iterations

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    dold=0.
    dx=1. !Is this a good choice?

    DO i=1,n

       !dnew=(f(x+dx)-f(x))/dx !Old method
       dnew=(f(x+dx/2.)-f(x-dx/2.))/dx !New, using equal sided derivative

       IF(i>=imin .AND. ABS(dnew/dold-1.)<acc) THEN
          !derivative=dnew
          EXIT
       ELSE IF(i==n) THEN
          !derivative=dnew
          STOP 'DERIVATIVE: Error, maximum number of iterations exceeded'
       ELSE
          dold=dnew
          dx=dx/2.
       END IF

    END DO

    derivative=dnew

  END FUNCTION derivative

  FUNCTION derivative_x(f,x,y,acc)

    !Calculates the derivative of a function 'f' at the point x to accuracy acc!

    IMPLICIT NONE
    REAL :: acc, x, dnew, dold, dx, derivative_x, y
    INTEGER :: i

    INTEGER, PARAMETER :: n=100

    INTERFACE
       REAL FUNCTION f(x,y)
         REAL, INTENT(IN) :: x,y
       END FUNCTION f
    END INTERFACE

    dold=0.
    dx=4.
    
    DO i=1,n

       dnew=(f(x+dx,y)-f(x,y))/dx

       !WRITE(*,*) i, dx, dnew, dold

       IF(i>1 .AND. ABS(dnew/dold-1.)<acc) THEN
          !derivative_x=dnew
          EXIT
       ELSE IF(i==n) THEN
          STOP 'DERIVATIVE_X: Error, maximum number of iterations exceeded'
       ELSE
          dold=dnew
          dx=dx/2.
       END IF
       
    END DO

    derivative_x=dnew

  END FUNCTION derivative_x

  FUNCTION derivative_y(f,x,y,acc)

    !Calculates the derivative of a function 'f' at the point x to accuracy acc!

    IMPLICIT NONE
    REAL :: acc, x, y, dnew, dold, dy, derivative_y
    INTEGER :: i

    INTEGER, PARAMETER :: n=100

    INTERFACE
       REAL FUNCTION f(x,y)
         REAL, INTENT(IN) :: x,y
       END FUNCTION f
    END INTERFACE

    dy=4.
    dold=0.
    
    DO i=1,n

       dnew=(f(x,y+dy)-f(x,y))/dy

       !WRITE(*,*) i, dx, dnew, dold

       IF(i>1 .AND. ABS(dnew/dold-1.)<acc) THEN
          !derivative_y=dnew
          EXIT
       ELSE IF(i==n) THEN
          STOP 'DERIVATIVE_Y: Error, maximum number of iterations exceeded'
       ELSE
          dold=dnew
          dy=dy/2.
       END IF
       
    END DO

    derivative_y=dnew

  END FUNCTION derivative_y

  FUNCTION integrate_basic(a,b,f,n,iorder)

    !Integrates between a and b with nint points!
    !Not adaptive so no error control
    IMPLICIT NONE
    REAL :: integrate_basic
    REAL, INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: n, iorder
    INTEGER :: i
    REAL :: x, dx, weight
    DOUBLE PRECISION :: sum

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       integrate_basic=0.

    ELSE

       !Set the sum variable
       sum=0.d0

       DO i=1,n
          
          x=a+(b-a)*REAL(i-1)/REAL(n-1)

          IF(iorder==1) THEN
             !Composite trapezium weights
             IF(i==1 .OR. i==n) THEN
                weight=0.5
             ELSE
                weight=1.
             END IF
          ELSE IF(iorder==2) THEN
             !Composite extended formula weights
             IF(i==1 .OR. i==n) THEN
                weight=0.4166666666
             ELSE IF(i==2 .OR. i==n-1) THEN
                weight=1.0833333333
             ELSE
                weight=1.
             END IF
          ELSE IF(iorder==3) THEN
             !Composite Simpson weights
             IF(i==1 .OR. i==n) THEN
                weight=0.375
             ELSE IF(i==2 .OR. i==n-1) THEN
                weight=1.1666666666
             ELSE IF(i==3 .OR. i==n-2) THEN
                weight=0.9583333333
             ELSE
                weight=1.
             END IF
          ELSE
             STOP 'INTEGERATE_BASIC: Error, order specified incorrectly'
          END IF

          sum=sum+weight*f(x)

       END DO

       dx=(b-a)/REAL(n-1)
       integrate_basic=REAL(sum)*dx

       !WRITE(*,*) 'INTEGRATE_BASIC: Order:', iorder
       !WRITE(*,*) 'INTEGRATE_BASIC: Nint:', n

    END IF

  END FUNCTION integrate_basic

!!$  FUNCTION integrate_old(a,b,f,acc,iorder)
!!$
!!$    !Integrates between a and b until desired accuracy is reached
!!$    IMPLICIT NONE
!!$    REAL :: integrate_old
!!$    REAL, INTENT(IN) :: a, b, acc
!!$    INTEGER, INTENT(IN) :: iorder
!!$    INTEGER :: i, j
!!$    INTEGER :: n
!!$    REAL :: x, weight, dx
!!$    DOUBLE PRECISION :: sum1, sum2
!!$    INTEGER, PARAMETER :: jmax=20
!!$    INTEGER, PARAMETER :: ninit=8
!!$
!!$    INTERFACE
!!$       FUNCTION f(x)
!!$         REAL :: f
!!$         REAL, INTENT(IN) :: x
!!$       END FUNCTION f
!!$    END INTERFACE
!!$
!!$    IF(a==b) THEN
!!$
!!$       integrate_old=0.
!!$
!!$    ELSE
!!$
!!$       !Set the sum variables
!!$       sum1=0.d0
!!$       sum2=0.d0
!!$
!!$       DO j=1,jmax
!!$
!!$          n=ninit*(2**(j-1))+1
!!$
!!$          DO i=1,n
!!$
!!$             x=a+(b-a)*REAL(i-1)/REAL(n-1)
!!$
!!$             IF(iorder==1) THEN
!!$                !Composite trapezium weights
!!$                IF(i==1 .OR. i==n) THEN
!!$                   weight=0.5
!!$                ELSE
!!$                   weight=1.
!!$                END IF
!!$             ELSE IF(iorder==2) THEN
!!$                !Composite extended formula weights
!!$                IF(i==1 .OR. i==n) THEN
!!$                   weight=0.416666666666666666
!!$                ELSE IF(i==2 .OR. i==n-1) THEN
!!$                   weight=1.083333333333333333
!!$                ELSE
!!$                   weight=1.
!!$                END IF
!!$             ELSE IF(iorder==3) THEN
!!$                !Composite Simpson weights
!!$                IF(i==1 .OR. i==n) THEN
!!$                   weight=0.375
!!$                ELSE IF(i==2 .OR. i==n-1) THEN
!!$                   weight=1.166666666666666666
!!$                ELSE IF(i==3 .OR. i==n-2) THEN
!!$                   weight=0.958333333333333333
!!$                ELSE
!!$                   weight=1.
!!$                END IF
!!$             ELSE IF(iorder==4) THEN
!!$                !Mead weights, these are less good than the composite Simpson above
!!$                IF(i==1 .OR. i==n) THEN
!!$                   weight=0.458333333333333333
!!$                !ELSE IF(i==2 .OR. i==n-1) THEN
!!$                !   weight=1.
!!$                ELSE IF(i==3 .OR. i==n-2) THEN
!!$                   weight=1.041666666666666666
!!$                ELSE
!!$                   weight=1.
!!$                END IF
!!$             ELSE
!!$                STOP 'INTEGERATE: Error, order specified incorrectly'
!!$             END IF
!!$             
!!$             sum2=sum2+weight*f(x)
!!$
!!$          END DO
!!$
!!$          dx=(b-a)/REAL(n-1)
!!$          sum2=sum2*dx
!!$
!!$          IF(j .NE. 1 .AND. ABS(-1.d0+sum2/sum1)<acc) THEN
!!$             !integrate_old=REAL(sum2)
!!$             !WRITE(*,*) 'INTEGRATE_OLD: Order:', iorder
!!$             !WRITE(*,*) 'INTEGRATE_OLD: Nint:', n
!!$             EXIT
!!$          ELSE IF(j==jmax) THEN
!!$             STOP 'INTEGRATE_OLD: Integration timed out'
!!$          ELSE
!!$             sum1=sum2
!!$             sum2=0.d0
!!$          END IF
!!$
!!$       END DO
!!$
!!$       integrate_old=REAL(sum2)
!!$
!!$    END IF
!!$
!!$  END FUNCTION integrate_old

  FUNCTION integrate(a,b,f,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
    
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax
          
          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             !The first go is just the trapezium of the end points
             f1=f(a)
             f2=f(b)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             !integrate=REAL(sum_new)
             !WRITE(*,*) 'INTEGRATE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate=REAL(sum_new)

    END IF

  END FUNCTION integrate

!!$  FUNCTION integrate_pointer(a,b,f,acc,iorder)
!!$
!!$    !Integrates between a and b until desired accuracy is reached
!!$    !Stores information to reduce function calls
!!$    IMPLICIT NONE
!!$    REAL :: integrate_pointer
!!$    REAL, INTENT(IN) :: a, b, acc
!!$    INTEGER, INTENT(IN) :: iorder
!!$    INTEGER :: i, j
!!$    INTEGER :: n
!!$    REAL :: x, dx
!!$    REAL :: f1, f2, fx
!!$    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
!!$    INTEGER, PARAMETER :: jmin=5
!!$    INTEGER, PARAMETER :: jmax=30
!!$
!!$    ABSTRACT INTERFACE
!!$       FUNCTION f(x)
!!$         REAL :: f
!!$         REAL, INTENT(IN) :: x
!!$       END FUNCTION f
!!$    END INTERFACE
!!$
!!$    procedure (func), pointer :: f_ptr => null ()
!!$
!!$    IF(a==b) THEN
!!$
!!$       !Fix the answer to zero if the integration limits are identical
!!$       integrate_pointer=0.
!!$
!!$    ELSE
!!$
!!$       !Set the sum variable for the integration
!!$       sum_2n=0.d0
!!$       sum_n=0.d0
!!$       sum_old=0.d0
!!$       sum_new=0.d0
!!$
!!$       DO j=1,jmax
!!$          
!!$          !Note, you need this to be 1+2**n for some integer n
!!$          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
!!$          n=1+2**(j-1)
!!$
!!$          !Calculate the dx interval for this value of 'n'
!!$          dx=(b-a)/REAL(n-1)
!!$
!!$          IF(j==1) THEN
!!$             
!!$             !The first go is just the trapezium of the end points
!!$             f1=f(a)
!!$             f2=f(b)
!!$             sum_2n=0.5d0*(f1+f2)*dx
!!$             sum_new=sum_2n
!!$             
!!$          ELSE
!!$
!!$             !Loop over only new even points to add these to the integral
!!$             DO i=2,n,2
!!$                x=a+(b-a)*REAL(i-1)/REAL(n-1)
!!$                fx=f(x)
!!$                sum_2n=sum_2n+fx
!!$             END DO
!!$
!!$             !Now create the total using the old and new parts
!!$             sum_2n=sum_n/2.d0+sum_2n*dx
!!$
!!$             !Now calculate the new sum depending on the integration order
!!$             IF(iorder==1) THEN  
!!$                sum_new=sum_2n
!!$             ELSE IF(iorder==3) THEN         
!!$                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
!!$             ELSE
!!$                STOP 'INTEGRATE_POINTER: Error, iorder specified incorrectly'
!!$             END IF
!!$
!!$          END IF
!!$
!!$          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
!!$             EXIT
!!$          ELSE IF(j==jmax) THEN
!!$             STOP 'INTEGRATE_POINTER: Integration timed out'
!!$          ELSE
!!$             !Integral has not converged so store old sums and reset sum variables
!!$             sum_old=sum_new
!!$             sum_n=sum_2n
!!$             sum_2n=0.d0
!!$          END IF
!!$
!!$       END DO
!!$
!!$       integrate_pointer=REAL(sum_new)
!!$
!!$    END IF
!!$
!!$  END FUNCTION integrate_pointer

  FUNCTION integrate_log(a,b,f,acc,iorder,ilog)

    !Integrates between a and b until desired accuracy is reached!
    IMPLICIT NONE
    REAL :: integrate_log
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: ilog, iorder
    INTEGER :: i, j, n
    REAL :: x, weight, dx, lima, limb 
    DOUBLE PRECISION :: sum1, sum2
    
    INTEGER, PARAMETER :: jmax=20
    INTEGER, PARAMETER :: ninit=8

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       integrate_log=0.

    ELSE

       !Set the sum variables
       sum1=0.d0
       sum2=0.d0

       IF(ilog==1) THEN
          lima=log(a)
          limb=log(b)
       ELSE
          lima=a
          limb=b
       END IF

       DO j=1,jmax

          n=ninit*(2**(j-1))

          DO i=1,n

             x=lima+(limb-lima)*REAL(i-1)/REAL(n-1)

             IF(ilog==1) THEN
                x=exp(x)
             END IF

             IF(iorder==1) THEN
                !Composite trapezium weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.5
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==2) THEN
                !Composite extended formula weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.4166666666
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.0833333333
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==3) THEN
                !Composite Simpson weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.375
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.1666666666
                ELSE IF(i==3 .OR. i==n-2) THEN
                   weight=0.9583333333
                ELSE
                   weight=1.
                END IF
             ELSE
                STOP 'INTEGERATE_LOG: Error, order specified incorrectly'
             END IF

             IF(ilog==0) THEN
                sum2=sum2+weight*f(x)
             ELSE IF(ilog==1) THEN
                sum2=sum2+weight*f(x)*x
             ELSE
                STOP 'INTEGRATE_LOG: Error, ilog specified incorrectly'
             END IF

          END DO

          dx=(limb-lima)/REAL(n-1)
          sum2=sum2*dx

          IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
             !integrate_log=REAL(sum2)
             !WRITE(*,*) 'INTEGRATE_LOG: Order:', iorder
             !WRITE(*,*) 'INTEGRATE_LOG: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_LOG: Integration timed out'
          ELSE
             sum1=sum2
             sum2=0.d0
          END IF

       END DO

       integrate_log=REAL(sum2)

    END IF

  END FUNCTION integrate_log

  FUNCTION cubeint(a,b,f,acc)

    USE fix_polynomial

    !Integrates between a and b until desired accuracy is reached!
    !Fits a cubic between successive 4 points
    !Only useful if points are not eqaully spaced, thus this routine is probably redundant
    IMPLICIT NONE
    REAL :: cubeint 
    REAL, INTENT(IN) :: a, b, acc
    INTEGER :: i, j, nint, nsec
    REAL :: a3, a2, a1, a0
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    DOUBLE PRECISION :: sum1, sum2
    
    INTEGER, PARAMETER :: jmax=20
    INTEGER, PARAMETER :: ni=1

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       cubeint=0.

    ELSE

       !Set the sum variables
       sum1=0.d0
       sum2=0.d0

       DO j=1,jmax

          !This is the number of cubic sections (each of which has four function evaluations)
          nsec=ni*2**(j-1)

          !Number of function evaluation points so as to be able to fit a cubic (4,7,10 ...)
          nint=3*nsec+1

          DO i=1,nsec

             IF(i==1) THEN

                x1=a+(b-a)*float(3*(i-1)+1-1)/float(nint-1)
                y1=f(x1)

             ELSE

                x1=x4
                y1=y4

             END IF

             x2=a+(b-a)*float(3*(i-1)+2-1)/float(nint-1)
             x3=a+(b-a)*float(3*(i-1)+3-1)/float(nint-1)
             x4=a+(b-a)*float(3*(i-1)+4-1)/float(nint-1)

             y2=f(x2)
             y3=f(x3)
             y4=f(x4) 
           
             CALL fix_cubic(a3,a2,a1,a0,x1,y1,x2,y2,x3,y3,x4,y4)

             !Add the (analytical) intergal of a cubic between points x1 and x4 to the total
             sum2=sum2+(a3/4.)*(x4**4.-x1**4.)+(a2/3.)*(x4**3.-x1**3.)+(a1/2.)*(x4**2.-x1**2.)+a0*(x4-x1)

          END DO

          IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
             !cubeint=REAL(sum2)
             !WRITE(*,*) 'CUBEINT: Number of sections', nsec
             !WRITE(*,*) 'CUBEINT: Number of function points', nint
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'CUBEINT: Integration timed out'
          ELSE
             sum1=sum2
             sum2=0.d0
          END IF

       END DO

       cubeint=REAL(sum2)

    END IF

  END FUNCTION cubeint

  FUNCTION integrate_jac(a,b,f,acc,iorder,g,gi,dg)

    !Integrates between a and b until desired accuracy is reached
    !Uses a Jacobian to speed up the integration
    IMPLICIT NONE
    REAL :: integrate_jac
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    INTEGER :: i, j, n
    REAL :: dy, alim, blim
    REAL :: x, y, weight
    DOUBLE PRECISION :: sum1, sum2
    
    INTEGER, PARAMETER :: jmax=20
    INTEGER, PARAMETER :: ninit=8

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
       FUNCTION g(x)
         REAL :: g
         REAL, INTENT(IN) :: x
       END FUNCTION g
       FUNCTION gi(x)
         REAL :: gi
         REAL, INTENT(IN) :: x
       END FUNCTION gi
       FUNCTION dg(x)
         REAL :: dg
         REAL, INTENT(IN) :: x
       END FUNCTION dg
    END INTERFACE
    
    IF(a==b) THEN

       integrate_jac=0.

    ELSE

       !Set the sum variables
       sum1=0.d0
       sum2=0.d0

       alim=g(a)
       blim=g(b)

       DO j=1,jmax

          n=ninit*(2**(j-1))

          DO i=1,n

             y=alim+(blim-alim)*REAL(i-1)/REAL(n-1)

             !IF(i==1 .OR. i==n) THEN
                !multiple of 0.5 for beginning and end and multiple of 2 for middle points!
             !   weight=0.5
             !ELSE
             !   weight=1.
             !END IF

             IF(iorder==1) THEN
                !Composite trapezium weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.5
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==2) THEN
                !Composite extended formula weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.4166666666
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.0833333333
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==3) THEN
                !Composite Simpson weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.375
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.1666666666
                ELSE IF(i==3 .OR. i==n-2) THEN
                   weight=0.9583333333
                ELSE
                   weight=1.
                END IF
             ELSE
                STOP 'INTEGERATE_JAC: Error, order specified incorrectly'
             END IF

             x=gi(y)

             sum2=sum2+weight*f(x)/dg(x)

          END DO

          dy=(blim-alim)/REAL(n-1)
          sum2=sum2*dy

          IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
             !integrate_jac=REAL(sum2)
             !WRITE(*,*) 'INTEGRATE_JAC: Order:', iorder
             !WRITE(*,*) 'INTEGRATE_JAC: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_JAC: Integration timed out'
          ELSE
             sum1=sum2
             sum2=0.d0
          END IF

       END DO

       integrate_jac=REAL(sum2)

    END IF

  END FUNCTION integrate_jac

END MODULE calculus
