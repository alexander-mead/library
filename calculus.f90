MODULE calculus

  USE array_operations

CONTAINS

  REAL FUNCTION derivative(f,x,acc)

    ! Calculates the derivative of a function 'f' at the point x to accuracy acc!
    IMPLICIT NONE
    REAL, INTENT(IN) :: x, acc
    REAL :: dnew, dold, dx 
    INTEGER :: i

    INTEGER, PARAMETER :: imin=5 ! Minimum number of iterations
    INTEGER, PARAMETER :: n=100 ! Maximum number of iterations

    INTERFACE
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    dold=0.
    dx=1. ! Is this a good choice?

    DO i=1,n

       dnew=(f(x+dx/2.)-f(x-dx/2.))/dx !New, using equal sided derivative

       IF(i>=imin .AND. abs(dnew/dold-1.)<acc) THEN
          EXIT
       ELSE IF(i==n) THEN
          STOP 'DERIVATIVE: Error, maximum number of iterations exceeded'
       ELSE
          dold=dnew
          dx=dx/2.
       END IF

    END DO

    derivative=dnew

  END FUNCTION derivative

  FUNCTION derivative_x(f,x,y,acc)

    ! Calculates the derivative of a function 'f' at the point x to accuracy acc!
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

       IF(i>1 .AND. abs(dnew/dold-1.)<acc) THEN
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

  REAL FUNCTION derivative_y(f,x,y,acc)

    ! Calculates the derivative of a function 'f' at the point x to accuracy acc!

    IMPLICIT NONE
    REAL :: acc, x, y, dnew, dold, dy
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

       IF(i>1 .AND. abs(dnew/dold-1.)<acc) THEN
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

  REAL FUNCTION integrate_basic(a,b,f,n,iorder)

    ! Integrates between a and b with n points; not adaptive so no error control
    IMPLICIT NONE
    REAL, INTENT(IN) :: a, b ! Integration limits
    INTEGER, INTENT(IN) :: n ! Number of points
    INTEGER, INTENT(IN) :: iorder ! Order for integration
    INTEGER :: i
    REAL :: x, dx, weight
    DOUBLE PRECISION :: sum

    INTERFACE       
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       integrate_basic=0.

    ELSE

       ! Set the sum variable
       sum=0.

       DO i=1,n
          
          !x=a+(b-a)*real(i-1)/real(n-1)
          x=progression(a,b,i,n)

          IF(iorder==1) THEN
             ! Composite trapezium weights
             IF(i==1 .OR. i==n) THEN
                weight=0.5
             ELSE
                weight=1.
             END IF
          ELSE IF(iorder==2) THEN
             ! Composite extended formula weights
             IF(i==1 .OR. i==n) THEN
                weight=0.4166666666
             ELSE IF(i==2 .OR. i==n-1) THEN
                weight=1.0833333333
             ELSE
                weight=1.
             END IF
          ELSE IF(iorder==3) THEN
             ! Composite Simpson weights
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

       dx=(b-a)/real(n-1)
       integrate_basic=real(sum)*dx

       !WRITE(*,*) 'INTEGRATE_BASIC: Order:', iorder
       !WRITE(*,*) 'INTEGRATE_BASIC: Nint:', n

    END IF

  END FUNCTION integrate_basic

  FUNCTION integrate(a,b,f,acc,iorder)

    ! Integrates between a and b until desired accuracy is reached
    ! Stores information to reduce function calls
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
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f       
    END INTERFACE

    IF(a==b) THEN

       ! Fix the answer to zero if the integration limits are identical
       integrate=0.

    ELSE

       ! Set the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax
          
          ! Note, you need this to be 1+2**n for some integer n
          ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          ! Calculate the dx interval for this value of 'n'
          dx=(b-a)/real(n-1)

          IF(j==1) THEN
             
             ! The first go is just the trapezium of the end points
             f1=f(a)
             f2=f(b)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             ! Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(a,b,i,n)
                fx=f(x)
                sum_2n=sum_2n+fx
             END DO

             ! Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             ! Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (abs(-1.+sum_new/sum_old)<acc)) THEN
             ! Converged
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE: Integration timed out'
          ELSE
             ! Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

       integrate=real(sum_new)

    END IF

  END FUNCTION integrate

  FUNCTION integrate_log(a,b,f,acc,iorder,ilog)

    ! Integrates between a and b until desired accuracy is reached!
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
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       integrate_log=0.

    ELSE

       ! Set the sum variables
       sum1=0.
       sum2=0.

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

             x=lima+(limb-lima)*real(i-1)/real(n-1)

             IF(ilog==1) THEN
                x=exp(x)
             END IF

             IF(iorder==1) THEN
                ! Composite trapezium weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.5
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==2) THEN
                ! Composite extended formula weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.4166666666
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.0833333333
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==3) THEN
                ! Composite Simpson weights
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

          dx=(limb-lima)/real(n-1)
          sum2=sum2*dx

          IF(j .NE. 1 .AND. abs(-1.+sum2/sum1)<acc) THEN
             !integrate_log=real(sum2)
             !WRITE(*,*) 'INTEGRATE_LOG: Order:', iorder
             !WRITE(*,*) 'INTEGRATE_LOG: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_LOG: Integration timed out'
          ELSE
             sum1=sum2
             sum2=0.
          END IF

       END DO

       integrate_log=real(sum2)

    END IF

  END FUNCTION integrate_log

  FUNCTION cubeint(a,b,f,acc)

    USE fix_polynomial

    ! Integrates between a and b until desired accuracy is reached!
    ! Fits a cubic between successive 4 points
    ! Only useful if points are not eqaully spaced, thus this routine is probably redundant
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
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       cubeint=0.

    ELSE

       ! Set the sum variables
       sum1=0.
       sum2=0.

       DO j=1,jmax

          ! This is the number of cubic sections (each of which has four function evaluations)
          nsec=ni*2**(j-1)

          ! Number of function evaluation points so as to be able to fit a cubic (4,7,10 ...)
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

             ! Add the (analytical) intergal of a cubic between points x1 and x4 to the total
             sum2=sum2+(a3/4.)*(x4**4.-x1**4.)+(a2/3.)*(x4**3.-x1**3.)+(a1/2.)*(x4**2.-x1**2.)+a0*(x4-x1)

          END DO

          IF(j .NE. 1 .AND. abs(-1.+sum2/sum1)<acc) THEN
             !WRITE(*,*) 'CUBEINT: Number of sections', nsec
             !WRITE(*,*) 'CUBEINT: Number of function points', nint
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'CUBEINT: Integration timed out'
          ELSE
             sum1=sum2
             sum2=0.
          END IF

       END DO

       cubeint=real(sum2)

    END IF

  END FUNCTION cubeint

  FUNCTION integrate_jac(a,b,f,acc,iorder,g,gi,dg)

    ! Integrates between a and b until desired accuracy is reached
    ! Uses a Jacobian to speed up the integration
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
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f       
       REAL FUNCTION g(x)
         REAL, INTENT(IN) :: x
       END FUNCTION g       
       REAL FUNCTION gi(x)
         REAL, INTENT(IN) :: x
       END FUNCTION gi       
       REAL FUNCTION dg(x)
         REAL, INTENT(IN) :: x
       END FUNCTION dg       
    END INTERFACE
    
    IF(a==b) THEN

       integrate_jac=0.

    ELSE

       ! Set the sum variables
       sum1=0.
       sum2=0.

       alim=g(a)
       blim=g(b)

       DO j=1,jmax

          n=ninit*(2**(j-1))

          DO i=1,n

             y=alim+(blim-alim)*real(i-1)/real(n-1)

             IF(iorder==1) THEN
                ! Composite trapezium weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.5
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==2) THEN
                ! Composite extended formula weights
                IF(i==1 .OR. i==n) THEN
                   weight=0.4166666666
                ELSE IF(i==2 .OR. i==n-1) THEN
                   weight=1.0833333333
                ELSE
                   weight=1.
                END IF
             ELSE IF(iorder==3) THEN
                ! Composite Simpson weights
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

          dy=(blim-alim)/real(n-1)
          sum2=sum2*dy

          IF(j .NE. 1 .AND. abs(-1.+sum2/sum1)<acc) THEN
             !integrate_jac=real(sum2)
             !WRITE(*,*) 'INTEGRATE_JAC: Order:', iorder
             !WRITE(*,*) 'INTEGRATE_JAC: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE_JAC: Integration timed out'
          ELSE
             sum1=sum2
             sum2=0.
          END IF

       END DO

       integrate_jac=real(sum2)

    END IF

  END FUNCTION integrate_jac

  REAL FUNCTION integrate_monte_carlo(a,b,f,n)

    ! Integrates between a and b with n points; not adaptive so no error control
    USE random_numbers
    IMPLICIT NONE
    REAL, INTENT(IN) :: a, b ! Integration limits
    INTEGER, INTENT(IN) :: n ! Number of points
    INTEGER :: i
    REAL :: x, dx
    DOUBLE PRECISION :: sum

    INTERFACE       
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       integrate_monte_carlo=0.

    ELSE

       sum=0.d0
       DO i=1,n
          x=random_uniform(a,b)
          sum=sum+f(x)          
       END DO

       dx=(b-a)/real(n)

       sum=sum*dx

       integrate_monte_carlo=real(sum)

    END IF

  END FUNCTION integrate_monte_carlo

END MODULE calculus
