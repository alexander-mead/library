MODULE ODE_solvers

   USE array_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ODE
   PUBLIC :: ODE_adaptive

   INTEGER, PARAMETER :: iode_crude = 1
   INTEGER, PARAMETER :: iode_mid = 2
   INTEGER, PARAMETER :: iode_RK4 = 3

CONTAINS

   SUBROUTINE ODE(x, v, t, ti, tf, xi, vi, fx, fv, n, iode, ilog)

      ! Solves 2nd order ODE d2x/dt2 from ti to tf and creates arrays of x, v, t values
      ! I have sometimes called this ODE_crass because it has a fixed number of time steps, n
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n, iode
      LOGICAL, INTENT(IN) :: ilog
      DOUBLE PRECISION :: x8(n), v8(n)
      DOUBLE PRECISION, ALLOCATABLE :: t8(:)
      INTEGER :: i

      INTERFACE

         ! fx is what dx/dt is equal to, this is almost always just v
         FUNCTION fx(x, v, t)
            REAL, INTENT(IN) :: x, v, t
         END FUNCTION fx

         ! fv is what dv/dt is equal to
         FUNCTION fv(x, v, t)
            REAL, INTENT(IN) :: x, v, t
         END FUNCTION fv

      END INTERFACE

      ! xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
      x8(1) = xi
      v8(1) = vi

      ! Fill the time array
      IF (ilog) THEN
         CALL fill_array_double(dble(log(ti)), dble(log(tf)), t8, n)
         t8 = exp(t8)
      ELSE
         CALL fill_array_double(dble(ti), dble(tf), t8, n)
      END IF

      ! Advance the system through all n-1 time steps
      DO i = 1, n-1
         CALL ODE_advance(x8(i), x8(i+1), v8(i), v8(i+1), t8(i), t8(i+1), fx, fv, iode)
      END DO

      ! Allocate arrays for final solution and copy double-precision to single-precision
      ALLOCATE (x(n), v(n), t(n))
      x = real(x8)
      v = real(v8)
      t = real(t8)

      !WRITE(*,*) 'ODE: Integration complete in steps:', n

   END SUBROUTINE ODE

   SUBROUTINE ODE_adaptive(x, v, t, ti, tf, xi, vi, fx, fv, acc, iode, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
      ! acc is the desired accuracy across the entire solution
      ! time steps are increased until convergence is achieved
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iode
      LOGICAL, INTENT(IN) :: ilog
      DOUBLE PRECISION, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
      INTEGER :: i, j, n, k, np, ifail, kn

      INTEGER, PARAMETER :: jmax = 30   ! Maximum number of goes
      INTEGER, PARAMETER :: ninit = 100 ! Initial number of points

      INTERFACE

         ! fx is what dx/dt is equal to this is almost always just v
         FUNCTION fx(x, v, t)
            REAL, INTENT(IN) :: x, v, t
         END FUNCTION fx

         ! fv is what dv/dt is equal to
         FUNCTION fv(x, v, t)
            REAL, INTENT(IN) :: x, v, t
         END FUNCTION fv

      END INTERFACE

      ! Loop over attemps for number of time steps
      DO j = 1, jmax

         ! Set the number of time steps; always 1+m*(2**n)
         n = 1+ninit*(2**(j-1))

         ! Allocate double-precision arrays
         ALLOCATE (x8(n), v8(n), t8(n))

         ! Set the array entries to zero manually
         x8 = 0.d0
         v8 = 0.d0
         t8 = 0.d0

         ! xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
         x8(1) = xi
         v8(1) = vi

         ! Fill time-step array
         IF (ilog) THEN
            CALL fill_array_double(dble(log(ti)), dble(log(tf)), t8, n)
            t8 = exp(t8)
         ELSE
            CALL fill_array_double(dble(ti), dble(tf), t8, n)
         END IF

         ! Set the fail flag
         ifail = 0

         ! Loop over all time steps
         DO i = 1, n-1
            CALL ODE_advance(x8(i), x8(i+1), v8(i), v8(i+1), t8(i), t8(i+1), fx, fv, iode)
         END DO

         ! Automatically fail on the first go
         IF (j == 1) ifail = 1

         ! Check accuracy of result compared to previous go
         IF (j .NE. 1) THEN

            ! This is the number of points in the previous try
            np = 1+(n-1)/2

            ! Loop over the number of points in the previous attempt; k
            DO k = 1, np

               ! kn is k-new
               kn = 2*k-1

               ! If still okay then check the result
               IF (ifail == 0) THEN

                  ! Fail conditions (is the first part of this not dodgy?)
                  IF (xh(k) > acc .AND. x8(kn) > acc .AND. (abs(xh(k)/x8(kn))-1.) > acc) ifail = 1
                  IF (vh(k) > acc .AND. v8(kn) > acc .AND. (abs(vh(k)/v8(kn))-1.) > acc) ifail = 1

                  ! Deallocate arrays and exit loop if a failure has occured
                  IF (ifail == 1) THEN
                     DEALLOCATE (xh, th, vh)
                     EXIT
                  END IF

               END IF

            END DO

         END IF

         ! If the integration was successful then fill single-precicion arrays with solution and exit
         IF (ifail == 0) THEN
            !WRITE(*,*) 'ODE: Integration complete in steps:', n-1
            !WRITE(*,*)
            ALLOCATE (x(n), t(n), v(n))
            x = real(x8)
            v = real(v8)
            t = real(t8)
            EXIT
         END IF

         ! Otherwise allocate arrays to store this solution and move on to the next attempt
         ALLOCATE (xh(n), th(n), vh(n))
         xh = x8
         vh = v8
         th = t8
         DEALLOCATE (x8, t8, v8)

      END DO

   END SUBROUTINE ODE_adaptive

   SUBROUTINE ODE_advance(x1, x2, v1, v2, t1, t2, fx, fv, iode)

      ! Advances the ODE system from t1 to t2, updating x1 to x2 and v1 to v2
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: x1
      DOUBLE PRECISION, INTENT(OUT) :: x2
      DOUBLE PRECISION, INTENT(IN) :: v1
      DOUBLE PRECISION, INTENT(OUT) :: v2
      DOUBLE PRECISION, INTENT(IN) :: t1
      DOUBLE PRECISION, INTENT(IN) :: t2
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: iode ! ODE solving method
      REAL :: x, v, t, dt
      REAL :: kx1, kx2, kx3, kx4
      REAL :: kv1, kv2, kv3, kv4

      INTERFACE

         ! fx is what dx/dt is equal to; this is almost always just v for 2ODE
         FUNCTION fx(x, v, t)
            REAL, INTENT(IN) :: x, v, t
         END FUNCTION fx

         ! fv is what dv/dt is equal to
         FUNCTION fv(x, v, t)
            REAL, INTENT(IN) :: x, v, t
         END FUNCTION fv

      END INTERFACE

      ! Set x, v, t to be the initial state of the system
      x = real(x1)
      v = real(v1)
      t = real(t1)

      ! Calculate dt between the final and intial state
      dt = real(t2-t1)
   
      IF (iode == iode_crude) THEN

         ! Crude method
         kx1 = dt*fx(x, v, t)
         kv1 = dt*fv(x, v, t)

         x2 = x1+kx1
         v2 = v1+kv1

         
      ELSE IF (iode == iode_mid) THEN

         ! Mid-point method
         kx1 = dt*fx(x, v, t)
         kv1 = dt*fv(x, v, t)
         kx2 = dt*fx(x+kx1/2., v+kv1/2., t+dt/2.)
         kv2 = dt*fv(x+kx1/2., v+kv1/2., t+dt/2.)

         x2 = x1+kx2
         v2 = v1+kv2

         
      ELSE IF (iode == iode_RK4) THEN

         ! 4th order Rungge-Kutta
         kx1 = dt*fx(x, v, t)
         kv1 = dt*fv(x, v, t)
         kx2 = dt*fx(x+kx1/2., v+kv1/2., t+dt/2.)
         kv2 = dt*fv(x+kx1/2., v+kv1/2., t+dt/2.)
         kx3 = dt*fx(x+kx2/2., v+kv2/2., t+dt/2.)
         kv3 = dt*fv(x+kx2/2., v+kv2/2., t+dt/2.)
         kx4 = dt*fx(x+kx3, v+kv3, t+dt)
         kv4 = dt*fv(x+kx3, v+kv3, t+dt)

         x2 = x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.d0
         v2 = v1+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.d0

      ELSE

         STOP 'ODE_ADVANCE: Error, iode specified incorrectly'

      END IF

   END SUBROUTINE ODE_advance

END MODULE ODE_solvers
