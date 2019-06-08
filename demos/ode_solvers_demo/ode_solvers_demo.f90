PROGRAM ODE_solvers_test

  USE ODE_solvers
  USE array_operations
  
  IMPLICIT NONE
  REAL, ALLOCATABLE :: ts(:), t(:), v(:)
  REAL, ALLOCATABLE :: x(:), xs(:)
  !DOUBLE PRECISION, ALLOCATABLE :: a(:)
  INTEGER :: i, n, iode, imeth
  LOGICAL :: ilog, iadapt
  REAL :: ti, tf, xi, vi
  
  REAL, PARAMETER :: acc=0.001 !For adaptive integrations
  INTEGER, PARAMETER :: m=200 !For reduced size arrays

  WRITE(*,*)
  WRITE(*,*) 'ODE examples'
  WRITE(*,*) '1 - dx/dt=-x'
  WRITE(*,*) '2 - dx/dt= -x*t'
  WRITE(*,*) '3 - SHM'
  WRITE(*,*) '4 - SHM damped'
  WRITE(*,*) '5 - Polynomial solution'
  READ(*,*) iode
  WRITE(*,*)

  WRITE(*,*) 'Adaptive solover:'
  WRITE(*,*) 'F - Non adaptive'
  WRITE(*,*) 'T - Adaptive'
  READ(*,*) iadapt
  WRITE(*,*)

  WRITE(*,*) 'Method:'
  WRITE(*,*) '1 - Crude'
  WRITE(*,*) '2 - Mid point'
  WRITE(*,*) '3 - RK4'
  READ(*,*) imeth
  WRITE(*,*)

  WRITE(*,*) 'Log time spacing:'
  WRITE(*,*) 'F - No'
  WRITE(*,*) 'T - Yes'
  READ(*,*) ilog
  WRITE(*,*)

  IF(iode==1 .OR. iode==2) THEN
     ti=1.
     tf=2.
     xi=1.
     vi=1.
  ELSE IF(iode==3 .OR. iode==4) THEN
     ti=0.
     tf=10.
     xi=1.
     vi=0.
  ELSE IF(iode==5) THEN
     ti=1.
     tf=5.
     xi=2.
     vi=0.
  ELSE
     STOP 'ODE_SOLVER_TEST: Error, iode not specified correctly'
  END IF

  WRITE(*,*) 'ODE_SOLVER_TEST: Initial time:', ti
  WRITE(*,*) 'ODE_SOLVER_TEST: Final time:', tf
  WRITE(*,*) 'ODE_SOLVER_TEST: Initial position:', xi
  WRITE(*,*) 'ODE_SOLVER_TEST: Initial velocity:', vi
  WRITE(*,*)

  IF(iadapt) THEN
     CALL ode_adaptive(x,v,t,ti,tf,xi,vi,fx,fv,acc,imeth,ilog)
  ELSE
     n=100
     CALL ode(x,v,t,ti,tf,xi,vi,fx,fv,n,imeth,ilog)
  END IF

  IF(iadapt) THEN
     n=SIZE(x)
     ALLOCATE(xs(m),ts(m))
     CALL reduce_array(x,n,xs,m)
     CALL reduce_array(t,n,ts,m)
  ELSE
     ALLOCATE(xs(n),ts(n))
     xs=x
     ts=t
  END IF

  OPEN(7,file='results.dat')
  DO i=1,SIZE(xs)
     WRITE(7,*) ts(i), analytic(ts(i)), xs(i)
  END DO
  CLOSE(7)

CONTAINS

  FUNCTION fx(x,v,t)

    IMPLICIT NONE
    REAL :: fx
    REAL, INTENT(IN) :: x, v, t

    IF(iode==1) THEN
       fx=decay(x,v,t)
    ELSE IF(iode==2) THEN
       fx=decay2(x,v,t)
    ELSE IF(iode==3 .OR. iode==4 .OR. iode==5) THEN
       fx=vel(x,v,t)
    ELSE
       STOP 'FX: Error, iode not specified correctly'
    END IF

  END FUNCTION fx

  FUNCTION fv(x,v,t)

    IMPLICIT NONE
    REAL :: fv
    REAL, INTENT(IN) :: x, v, t

    IF(iode==1 .OR. iode==2) THEN
       fv=one(x,v,t)
    ELSE IF(iode==3) THEN
       fv=shm(x,v,t)
    ELSE IF(iode==4) THEN
       fv=shm_damp(x,v,t)
    ELSE IF(iode==5) THEN
       fv=poly(x,v,t)
    ELSE
       STOP 'FV: Error, iode not specified correctly'
    END IF

  END FUNCTION fv

  FUNCTION poly(x,v,t)

    IMPLICIT NONE
    REAL :: poly
    REAL, INTENT(IN) :: x, v, t

    poly=(x/(t**2.))-v/t

  END FUNCTION poly

  FUNCTION vel(x,v,t)

    IMPLICIT NONE
    REAL :: vel
    REAL, INTENT(IN) :: x, v, t
    REAL :: crap

    !Avoid compile-time warnings
    crap=x
    crap=t

    vel=v

  END FUNCTION vel

  FUNCTION one(x,v,t)

    IMPLICIT NONE
    REAL :: one
    REAL, INTENT(IN) :: x, v, t
    REAL :: crap

    !Avoid compile-time warnings
    crap=x
    crap=v
    crap=t

    one=1.
   
  END FUNCTION one

  FUNCTION shm(x,v,t)

    IMPLICIT NONE
    REAL :: shm
    REAL, INTENT(IN) :: x, v, t
    REAL :: crap

    !Avoid compile-time warnings
    crap=v
    crap=t

    shm=-x

  END FUNCTION shm

  FUNCTION shm_damp(x,v,t)

    IMPLICIT NONE
    REAL :: shm_damp
    REAL, INTENT(IN) :: x, v, t
    REAL :: crap

    !Avoid compile-time warnings
    crap=t

    shm_damp=-x-0.3*v

  END FUNCTION shm_damp

  FUNCTION decay(x,v,t)

    IMPLICIT NONE
    REAL :: decay
    REAL, INTENT(IN) :: x, v, t
    REAL :: crap

    !Avoid compile-time warnings
    crap=v
    crap=t

    decay=-x

  END FUNCTION decay

  FUNCTION decay2(x,v,t)

    IMPLICIT NONE
    REAL :: decay2
    REAL, INTENT(IN) :: x, v, t
    REAL :: crap

    !Avoid compile-time warnings
    crap=v
    crap=t

    decay2=-x*t

  END FUNCTION decay2

  FUNCTION analytic(t)

    IMPLICIT NONE
    REAL :: analytic
    REAL, INTENT(IN) :: t

    IF(iode==1) THEN
       analytic=exp(1.-t)
    ELSE IF(iode==2) THEN
       analytic=exp((1.-t**2.)/2.)
    ELSE IF(iode==3) THEN
       analytic=cos(t)
    ELSE IF(iode==4) THEN
       analytic=exp(-0.3*t/2.)*(cos(1.97737*t/2.)+0.1517*sin(1.97737*t/2.))
    ELSE IF(iode==5) THEN
       analytic=t+1./t
    ELSE
       STOP 'Incorrect choice of ODE made'
    END IF

  END FUNCTION analytic

END PROGRAM ODE_solvers_test
