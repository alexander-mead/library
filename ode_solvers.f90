MODULE ODE_solvers

  USE array_operations

CONTAINS

  SUBROUTINE ODE(x,v,t,ti,tf,xi,vi,fx,fv,n,imeth,ilog)

    !Solves 2nd order ODE x''(t) from ti to tf and creates arrays of x, v, t values
    !I have sometimes called this ODE_crass
    !It has a fixed number of time steps, n
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:), v(:), t(:)
    REAL, INTENT(IN) :: xi, vi, ti, tf
    INTEGER, INTENT(IN) :: n, imeth
    LOGICAL, INTENT(IN) :: ilog
    !REAL :: dt, x4, v4, t4
    !REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4    
    DOUBLE PRECISION, ALLOCATABLE :: x8(:), v8(:), t8(:)
    INTEGER :: i

    !imeth sets ODE solving method
    !imeth = 1: Crude method
    !imeth = 2: Mid-point method
    !imeth = 3: Runge-Kutta
    
    INTERFACE

       !fx is what x' is equal to
       FUNCTION fx(x,v,t)
         REAL :: fx
         REAL, INTENT(IN) :: x, v, t
       END FUNCTION fx

       !fv is what v' is equal to
       FUNCTION fv(x,v,t)
         REAL :: fv
         REAL, INTENT(IN) :: x, v, t
       END FUNCTION fv
       
    END INTERFACE

    !Allocate arrays
    ALLOCATE(x8(n),v8(n),t8(n))

    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    x8(1)=xi
    v8(1)=vi

    !Fill time array
    IF(ilog) THEN
       CALL fill_array8(log(ti),log(tf),t8,n)
       t8=exp(t8)
    ELSE
       CALL fill_array8(ti,tf,t8,n)
    END IF

    DO i=1,n-1

       CALL ODE_advance(x8(i),x8(i+1),v8(i),v8(i+1),t8(i),t8(i+1),fx,fv,imeth)

!!$       x4=REAL(x8(i))
!!$       v4=REAL(v8(i))
!!$       t4=REAL(t8(i))
!!$
!!$       !Time steps are varible length in the log case
!!$       dt=REAL(t8(i+1)-t8(i))
!!$
!!$       IF(imeth==1) THEN
!!$
!!$          !Crude method!
!!$          kx1=dt*fx(x4,v4,t4)
!!$          kv1=dt*fv(x4,v4,t4)
!!$
!!$          x8(i+1)=x8(i)+kx1
!!$          v8(i+1)=v8(i)+kv1
!!$
!!$       ELSE IF(imeth==2) THEN
!!$
!!$          !Mid-point method!
!!$          kx1=dt*fx(x4,v4,t4)
!!$          kv1=dt*fv(x4,v4,t4)
!!$          kx2=dt*fx(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$          kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$
!!$          x8(i+1)=x8(i)+kx2
!!$          v8(i+1)=v8(i)+kv2
!!$
!!$       ELSE IF(imeth==3) THEN
!!$
!!$          !RK4 (Holy Christ, this is so fast compared to above methods)!
!!$          kx1=dt*fx(x4,v4,t4)
!!$          kv1=dt*fv(x4,v4,t4)
!!$          kx2=dt*fx(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$          kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$          kx3=dt*fx(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
!!$          kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
!!$          kx4=dt*fx(x4+kx3,v4+kv3,t4+dt)
!!$          kv4=dt*fv(x4+kx3,v4+kv3,t4+dt)
!!$
!!$          x8(i+1)=x8(i)+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.d0
!!$          v8(i+1)=v8(i)+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.d0
!!$
!!$       ELSE
!!$
!!$          STOP 'ODE: Error, imeth specified incorrectly'
!!$
!!$       END IF

    END DO

    ALLOCATE(x(n),v(n),t(n))
    x=REAL(x8)
    v=REAL(v8)
    t=REAL(t8)

    WRITE(*,*) 'ODE: Integration complete in steps:', n

  END SUBROUTINE ODE

  SUBROUTINE ODE_adaptive(x,v,t,ti,tf,xi,vi,fx,fv,acc,imeth,ilog)

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    !acc is the desired accuracy across the entire solution
    !time steps are increased until convergence is achieved
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:), t(:), v(:)
    REAL, INTENT(IN) :: xi, vi, ti, tf, acc
    INTEGER, INTENT(IN) :: imeth
    LOGICAL, INTENT(IN) :: ilog
    !REAL :: dt, x4, v4, t4
    !REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    DOUBLE PRECISION, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)    
    INTEGER :: i, j, n, k, np, ifail, kn
   
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=100

    !imeth sets ODE solving method
    !imeth = 1: Crude method
    !imeth = 2: Mid-point method
    !imeth = 3: Runge-Kutta   

    INTERFACE

       !fx is what x' is equal to
       FUNCTION fx(x,v,t)
         REAL :: fx
         REAL, INTENT(IN) :: x, v, t
       END FUNCTION fx

       !fv is what v' is equal to
       FUNCTION fv(x,v,t)
         REAL :: fv
         REAL, INTENT(IN) :: x, v, t
       END FUNCTION fv
       
    END INTERFACE

    DO j=1,jmax

       n=1+ninit*(2**(j-1))

       ALLOCATE(x8(n),v8(n),t8(n))

       x8=0.d0
       v8=0.d0
       t8=0.d0

       !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
       x8(1)=xi
       v8(1)=vi

       !Fill time array
       IF(ilog) THEN
          CALL fill_array8(log(ti),log(tf),t8,n)
          t8=exp(t8)
       ELSE
          CALL fill_array8(ti,tf,t8,n)
       END IF

       ifail=0

       DO i=1,n-1

          CALL ODE_advance(x8(i),x8(i+1),v8(i),v8(i+1),t8(i),t8(i+1),fx,fv,imeth)
          
!!$          x4=REAL(x8(i))
!!$          v4=REAL(v8(i))
!!$          t4=REAL(t8(i))
!!$
!!$          dt=REAL(t8(i+1)-t8(i))
!!$
!!$          IF(imeth==1) THEN
!!$
!!$             !Crude method!
!!$             kx1=dt*fx(x4,v4,t4)
!!$             kv1=dt*fv(x4,v4,t4)
!!$             
!!$             x8(i+1)=x8(i)+kx1
!!$             v8(i+1)=v8(i)+kv1
!!$
!!$          ELSE IF(imeth==2) THEN
!!$
!!$             !Mid-point method!
!!$             kx1=dt*fx(x4,v4,t4)
!!$             kv1=dt*fv(x4,v4,t4)
!!$             kx2=dt*fx(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$             
!!$             x8(i+1)=x8(i)+kx2
!!$             v8(i+1)=v8(i)+kv2
!!$
!!$          ELSE IF(imeth==3) THEN
!!$
!!$             !RK4 (Holy Christ, this is so fast compared to above methods)!
!!$             kx1=dt*fx(x4,v4,t4)
!!$             kv1=dt*fv(x4,v4,t4)
!!$             kx2=dt*fx(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
!!$             kx3=dt*fx(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
!!$             kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
!!$             kx4=dt*fx(x4+kx3,v4+kv3,t4+dt)
!!$             kv4=dt*fv(x4+kx3,v4+kv3,t4+dt)
!!$
!!$             x8(i+1)=x8(i)+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.d0
!!$             v8(i+1)=v8(i)+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.d0
!!$
!!$          ELSE
!!$
!!$             STOP 'ODE: Error, imeth specified incorrectly'
!!$
!!$          END IF
          
       END DO

       IF(j==1) ifail=1

       IF(j .NE. 1) THEN

          np=1+(n-1)/2

          DO k=1,1+(n-1)/2

             kn=2*k-1

             IF(ifail==0) THEN

                IF(xh(k)>acc .AND. x8(kn)>acc .AND. (ABS(xh(k)/x8(kn))-1.)>acc) ifail=1
                IF(vh(k)>acc .AND. v8(kn)>acc .AND. (ABS(vh(k)/v8(kn))-1.)>acc) ifail=1

                IF(ifail==1) THEN
                   DEALLOCATE(xh,th,vh)
                   EXIT
                END IF

             END IF
          END DO

       END IF

       IF(ifail==0) THEN
          WRITE(*,*) 'ODE: Integration complete in steps:', n-1
          WRITE(*,*)
          ALLOCATE(x(n),t(n),v(n))
          x=REAL(x8)
          v=REAL(v8)
          t=REAL(t8)
          EXIT
       END IF

       WRITE(*,*) 'ODE: Integration at:', n-1
       ALLOCATE(xh(n),th(n),vh(n))
       xh=x8
       vh=v8
       th=t8
       DEALLOCATE(x8,t8,v8)

    END DO

  END SUBROUTINE ODE_adaptive

  SUBROUTINE ODE_advance(x1,x2,v1,v2,t1,t2,fx,fv,imeth)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x1, v1, t1, t2
    DOUBLE PRECISION, INTENT(OUT) :: x2, v2
    INTEGER, INTENT(IN) :: imeth
    REAL :: x, v, t, dt
    REAL :: kx1, kx2, kx3, kx4
    REAL :: kv1, kv2, kv3, kv4

     INTERFACE

       !fx is what x' is equal to
       FUNCTION fx(x,v,t)
         REAL :: fx
         REAL, INTENT(IN) :: x, v, t
       END FUNCTION fx

       !fv is what v' is equal to
       FUNCTION fv(x,v,t)
         REAL :: fv
         REAL, INTENT(IN) :: x, v, t
       END FUNCTION fv
       
    END INTERFACE
    
    x=REAL(x1)
    v=REAL(v1)
    t=REAL(t1)

    dt=REAL(t2-t1)

    IF(imeth==1) THEN

       !Crude method!
       kx1=dt*fx(x,v,t)
       kv1=dt*fv(x,v,t)

       x2=x1+kx1
       v2=v1+kv1

    ELSE IF(imeth==2) THEN

       !Mid-point method!
       kx1=dt*fx(x,v,t)
       kv1=dt*fv(x,v,t)
       kx2=dt*fx(x+kx1/2.,v+kv1/2.,t+dt/2.)
       kv2=dt*fv(x+kx1/2.,v+kv1/2.,t+dt/2.)

       x2=x1+kx2
       v2=v1+kv2

    ELSE IF(imeth==3) THEN

       !RK4 (Holy Christ, this is so fast compared to above methods)!
       kx1=dt*fx(x,v,t)
       kv1=dt*fv(x,v,t)
       kx2=dt*fx(x+kx1/2.,v+kv1/2.,t+dt/2.)
       kv2=dt*fv(x+kx1/2.,v+kv1/2.,t+dt/2.)
       kx3=dt*fx(x+kx2/2.,v+kv2/2.,t+dt/2.)
       kv3=dt*fv(x+kx2/2.,v+kv2/2.,t+dt/2.)
       kx4=dt*fx(x+kx3,v+kv3,t+dt)
       kv4=dt*fv(x+kx3,v+kv3,t+dt)

       x2=x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.d0
       v2=v1+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.d0

    ELSE

       STOP 'ODE_ADVANCE: Error, imeth specified incorrectly'

    END IF

  END SUBROUTINE ODE_advance

END MODULE ODE_solvers
