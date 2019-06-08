PROGRAM calculus_test

  USE random_numbers
  USE calculus
  USE constants
  
  IMPLICIT NONE
  REAL :: sol, true
  REAL :: lim1, lim2
  REAL :: t1, t2
  REAL :: df
  INTEGER :: i, j, iorder, n
  INTEGER :: iexample, itest
  REAL :: a, b, c, x0, xmin, xmax
  
  REAL, PARAMETER :: acc=1e-5
  INTEGER, PARAMETER :: iseed=0
  
  WRITE(*,*)
  
  WRITE(*,*) '1 - Test derivatives'
  WRITE(*,*) '2 - Test integrals'
  WRITE(*,*) '3 - Compare integration of constants'
  READ(*,*) itest
  WRITE(*,*)

  WRITE(*,*) 'Accuracy parameter:', acc
  WRITE(*,*)

  IF(itest==1) THEN

     WRITE(*,*) '1 - Derivative of quadratic (ax^2+bx+c)'
     READ(*,*) iexample
     WRITE(*,*)

     IF(iexample==1) THEN

        !Parameters of quadratics
        a=1.
        b=0.
        c=1.

        !Point at which to evaluate derivative
        x0=2.

        !Write quadratic parameters
        WRITE(*,*) 'a:', a
        WRITE(*,*) 'b:', b
        WRITE(*,*) 'c:', c

        !Compute derivative numerically and compare it to the truth
        !df=derivative(quadratic,x0,acc)
        !true=quadratic_derivative(x0)
        !WRITE(*,*) 'x0:', x0
        !WRITE(*,*) 'Derivative at x0:', df
        !WRITE(*,*) 'True solution:', true
        !WRITE(*,*) 'Error:', df/true-1.
        !WRITE(*,*)

        xmin=-10.
        xmax=10.
        n=51
        
        WRITE(*,*)
        WRITE(*,*) '  Iteration            x        numerical          analyic                error'
        WRITE(*,*) '==============================================================================='
        DO i=1,n
           x0=xmin-(xmax-xmin)*REAL(i-1)/REAL(n-1)
           df=derivative(quadratic,x0,acc)
           true=quadratic_derivative(x0)
           WRITE(*,*) i, x0, df, true, df/true-1.
        END DO
        WRITE(*,*) '==============================================================================='

     END IF

     WRITE(*,*)

  ELSE IF(itest==2) THEN

     WRITE(*,*) '1 - Integeral - x^2 (0 -> 5)'
     WRITE(*,*) '2 - Integeral - sin(x) (0 -> pi)'
     WRITE(*,*) '3 - Integeral - exp(x) (0 -> 2)'
     WRITE(*,*) '4 - Integeral - sin(x) (0 -> 151pi)'
     READ(*,*) iexample
     WRITE(*,*)

     IF(iexample==1) THEN
        lim1=1.
        lim2=5.
     ELSE IF(iexample==2) THEN
        lim1=0.
        lim2=pi
     ELSE IF(iexample==3) THEN
        lim1=0.
        lim2=3.
     ELSE IF(iexample==4) THEN
        lim1=0.
        lim2=151.*pi
     ELSE
        STOP 'INT_FUNCTION: Error, iexample not specified correctly'
     END IF

     true=solution(lim1,lim2,iexample)
     WRITE(*,*) 'True solution:', true
     WRITE(*,*)  

     WRITE(*,*) 'Integration accuracy for adaptive methods:', acc
     WRITE(*,*)

     CALL RNG_set(iseed)

     DO iorder=1,3,2

        n=16
        sol=integrate_basic(lim1,lim2,func,n,iorder)
        WRITE(*,*) 'Order:', iorder
        WRITE(*,*) 'Number of points:', n
        WRITE(*,*) 'Basic:', sol
        WRITE(*,*) 'Accuracy:', sol/true
        WRITE(*,*)

        n=1024
        sol=integrate_monte_carlo(lim1,lim2,func,n)
        WRITE(*,*) 'Number of points:', n
        WRITE(*,*) 'Monte Carlo:', sol
        WRITE(*,*) 'Accuracy:', sol/true
        WRITE(*,*)

        sol=integrate(lim1,lim2,func,acc,iorder)
        WRITE(*,*) 'Order:', iorder
        WRITE(*,*) 'Storage:', sol
        WRITE(*,*) 'Accuracy:', sol/true
        WRITE(*,*)

        !This is just to check integrate_log gives the same answer as integrate
        !when the log sampling is turned off
        !sol=integrate_log(lim1,lim2,func,acc,iorder,0)
        !WRITE(*,*) 'Log integration w/o log:', sol
        !WRITE(*,*) 'Accuracy:', sol/true
        !WRITE(*,*)

        IF(iexample==1) THEN

           sol=integrate_log(lim1,lim2,func,acc,iorder,1)
           WRITE(*,*) 'Log integration:', sol
           WRITE(*,*) 'Accuracy:', sol/true
           WRITE(*,*)

           sol=integrate_jac(lim1,lim2,func,acc,iorder,jac,jaci,djac)
           WRITE(*,*) 'Jacobian integration:', sol
           WRITE(*,*) 'Accuracy:', sol/true
           WRITE(*,*)

           sol=integrate_jac(lim1,lim2,func,acc,iorder,jac_w,jaci_w,djac_w)
           WRITE(*,*) 'Jacobian integration (wrong J):', sol
           WRITE(*,*) 'Accuracy:', sol/true
           WRITE(*,*)

        END IF

     END DO

     sol=integrate_cubic(lim1,lim2,func,acc)
     WRITE(*,*) 'Cubic integration:', sol
     WRITE(*,*) 'Accuracy:', sol/true
     WRITE(*,*)

     n=100000
     iorder=3
     WRITE(*,*) 'Timings for number of intergrations:', n
     WRITE(*,fmt='(A10,I5)') 'Order:', iorder
     WRITE(*,*)

     DO j=2,2

        CALL cpu_time(t1)
        DO i=1,n     
           !IF(j==1) sol=integrate_old(lim1,lim2,func,acc,iorder)
           IF(j==2) sol=integrate(lim1,lim2,func,acc,iorder)
        END DO
        CALL cpu_time(t2)
        IF(j==1) WRITE(*,*) 'Standard'
        IF(j==2) WRITE(*,*) 'Storage'
        WRITE(*,fmt='(A10,F20.10)') 'Solution:', sol
        WRITE(*,fmt='(A10,F20.10)') 'Time [s]:', t2-t1
        WRITE(*,*)

     END DO

  ELSE IF(itest==3) THEN

     lim1=0
     lim2=10
     iorder=3
     n=100000

     WRITE(*,*) 'Lower limit:', lim1
     WRITE(*,*) 'Upper limit:', lim2
     WRITE(*,*) 'Integration order:', iorder
     WRITE(*,*) 'Number of integrations to do:', n
     WRITE(*,*)
     
     DO i=1,6

        WRITE(*,*) 'log10(constant):', log10(const(1.))

        iexample=i
    
        CALL cpu_time(t1)
        DO j=1,n
           sol=integrate(lim1,lim2,const,acc,iorder)
        END DO
        CALL cpu_time(t2)

        WRITE(*,*) 'log10(Solution):', log10(sol)
        WRITE(*,*) 'Time [s]:', t2-t1
        WRITE(*,*)

     END DO

  END IF

CONTAINS

  FUNCTION quadratic(x)

    IMPLICIT NONE
    REAL :: quadratic
    REAL, INTENT(IN) :: x

    quadratic=a*x**2+b*x+c

  END FUNCTION quadratic

  FUNCTION quadratic_derivative(x)

    IMPLICIT NONE
    REAL :: quadratic_derivative
    REAL, INTENT(IN) :: x

    quadratic_derivative=2.*a*x+b

  END FUNCTION quadratic_derivative

  FUNCTION const(x)

    IMPLICIT NONE
    REAL :: const
    REAL, INTENT(IN) :: x
    REAL :: crap

    !Prevent compile-time warnings
    crap=x

    const=10**iexample
    
  END FUNCTION const

  FUNCTION func(x)

    IMPLICIT NONE
    REAL :: func
    REAL, INTENT(IN) :: x

    IF(iexample==1) THEN
       func=x**2.
    ELSE IF(iexample==2 .OR. iexample==4) THEN
       func=sin(x)
    ELSE IF(iexample==3) THEN
       func=exp(x)
    ELSE
       STOP 'FUNC: Error, iexample not set correctly'
    END IF

  END FUNCTION func

  FUNCTION solution(a,b,iexample)

    IMPLICIT NONE
    REAL :: solution
    REAL, INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: iexample

    IF(iexample==1) THEN
       solution=(b**3.-a**3.)/3.
    ELSE IF(iexample==2 .OR. iexample==4) THEN
       solution=-cos(b)+cos(a)
    ELSE IF(iexample==3) THEN
       solution=exp(b)-exp(a)
    ELSE
       STOP 'SOLUTION: Error, iexample not specified correctly'
    END IF

  END FUNCTION solution

  FUNCTION jac(x)

    IMPLICIT NONE
    REAL :: jac
    REAL, INTENT(IN) :: x

    jac=x**3.

  END FUNCTION jac

  FUNCTION jaci(x)

    IMPLICIT NONE
    REAL :: jaci
    REAL, INTENT(IN) :: x

    jaci=x**(1./3.)

  END FUNCTION jaci

  FUNCTION djac(x)

    IMPLICIT NONE
    REAL :: djac
    REAL, INTENT(IN) :: x

    djac=3.*(x**2.)

  END FUNCTION djac

  FUNCTION jac_w(x)

    IMPLICIT NONE
    REAL :: jac_w
    REAL, INTENT(IN) :: x

    jac_w=x**(-1.)

  END FUNCTION jac_w

  FUNCTION jaci_w(x)

    IMPLICIT NONE
    REAL :: jaci_w
    REAL, INTENT(IN) :: x

    jaci_w=x**(-1.)

  END FUNCTION jaci_w

  FUNCTION djac_w(x)

    IMPLICIT NONE
    REAL :: djac_w
    REAL, INTENT(IN) :: x

    djac_w=-x**(-2.)

  END FUNCTION djac_w

END PROGRAM calculus_test
