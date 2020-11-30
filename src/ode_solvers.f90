MODULE ODE_solvers

   USE basic_operations
   USE array_operations

   IMPLICIT NONE

   PRIVATE

   ! ODE integration routines
   PUBLIC :: ODE1
   PUBLIC :: ODE2
   PUBLIC :: ODEcoupled
   PUBLIC :: ODE1_adaptive
   PUBLIC :: ODE2_adaptive
   PUBLIC :: ODEcoupled_adaptive

   ! ODE integration methods
   PUBLIC :: iode_crude
   PUBLIC :: iode_mid
   PUBLIC :: iode_RK4

   ! ODE integration methods
   INTEGER, PARAMETER :: iode_crude = 1
   INTEGER, PARAMETER :: iode_mid = 2
   INTEGER, PARAMETER :: iode_RK4 = 3

   ! Adaptive ODE
   INTEGER, PARAMETER :: jmax_adapt = 20 ! Maximum number of attempts at solution
   INTEGER, PARAMETER :: nmin_adapt = 10 ! Minimum number of points in solution

CONTAINS

   SUBROUTINE ODE1(x, t, ti, tf, xi, fx, n, iode, ilog)

      ! Solves 2nd order ODE d2x/dt2 from ti to tf and creates arrays of x, v, t values
      ! I have sometimes called this ODE_crass because it has a fixed number of time steps, n
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      INTEGER :: i

      ! x' = fx
      INTERFACE       
         FUNCTION fx(x_in, t_in)
            REAL, INTENT(IN) :: x_in, t_in
         END FUNCTION fx
      END INTERFACE

      ! Allocate and set xi and vi to the initial values of x and v (i.e. x(ti), v(ti))
      ALLOCATE(x(n))
      x(1) = xi

      ! Fill the time array
      CALL fill_array(ti, tf, t, n, ilog)

      ! Advance the system through all n-1 time steps
      DO i = 1, n-1
         CALL ODE1_advance(x(i), x(i+1), t(i), t(i+1), fx, iode)
      END DO

   END SUBROUTINE ODE1

   SUBROUTINE ODE2(x, v, t, ti, tf, xi, vi, fv, n, iode, ilog)

      ! Integrates second-order ODE: x'' + ax' + bx = c from ti to tf
      ! Outputs array of x(t) and v(t)
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      INTEGER :: i

      ! x'' = v' = fv(x, x', t)
      INTERFACE        
         FUNCTION fv(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fv
      END INTERFACE
   
      ! Allocate and set xi and vi to the initial values of x and v
      ALLOCATE(x(n), v(n))
      x(1) = xi
      v(1) = vi

      ! Fill the time array
      CALL fill_array(ti, tf, t, n, ilog)

      ! Advance the system through all n-1 time steps
      DO i = 1, n-1
         CALL ODE2_advance(x(i), x(i+1), v(i), v(i+1), t(i), t(i+1), velocity, fv, iode)
      END DO

   END SUBROUTINE ODE2

   SUBROUTINE ODEcoupled(x, v, t, ti, tf, xi, vi, fx, fv, n, iode, ilog)

      ! Solves 2nd order ODE d2x/dt2 from ti to tf and creates arrays of x, v, t values
      ! I have sometimes called this ODE_crass because it has a fixed number of time steps, n
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      INTEGER :: i

      ! x' = fx; fx = v for a second-order ODE
      ! v' = fv
      INTERFACE        
         FUNCTION fx(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fx        
         FUNCTION fv(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fv
      END INTERFACE
   
      ! Allocate and set xi and vi to the initial values of x and v (i.e. x(ti), v(ti))
      ALLOCATE(x(n), v(n))
      x(1) = xi
      v(1) = vi

      ! Fill the time array
      CALL fill_array(ti, tf, t, n, ilog)

      ! Advance the system through all n-1 time steps
      DO i = 1, n-1
         CALL ODE2_advance(x(i), x(i+1), v(i), v(i+1), t(i), t(i+1), fx, fv, iode)
      END DO

   END SUBROUTINE ODEcoupled

   SUBROUTINE ODE1_adaptive(x, t, ti, tf, xi, fx, n, iode, acc, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
      ! acc is the desired accuracy across the entire solution
      ! time steps are increased until convergence is achieved
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, EXTERNAL :: fx
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      REAL, INTENT(IN) :: acc
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      REAL, ALLOCATABLE :: xh(:), th(:)
      INTEGER :: j, in, ip, nn, np
      LOGICAL :: fail
      INTEGER, PARAMETER :: jmax = jmax_adapt ! Maximum number of goes
      INTEGER, PARAMETER :: nmin = nmin_adapt ! Minimum number of points

      ! x' = fx
      INTERFACE       
         FUNCTION fx(x_in, t_in)
            REAL, INTENT(IN) :: x_in, t_in
         END FUNCTION fx
      END INTERFACE

      IF (n <= 1) STOP 'ODE1_ADAPTIVE: Error, you cannot solve an ODE with n=1 (or less)'

      ! Loop over attemps for number of time steps
      DO j = 1, jmax

         ! Set the new number of time steps
         nn = 1+(n-1)*2**(j-1)

         ! Solve the ODE with 'nn' time steps
         CALL ODE1(x, t, ti, tf, xi, fx, nn, iode, ilog)

         ! Automatically fail on the first go; or check for convergence with respect to previous solution
         ! TODO: Replace fail conditions with requal?
         IF (j == 1 .OR. nn < nmin) THEN          
            fail = .TRUE.
         ELSE
            fail = .FALSE.     
            DO ip = 1, np       
               in = 2*ip-1   
               IF (.NOT. almost_equal(xh(ip), x(in), acc)) THEN
               !IF ((.NOT. requal(xh(k), x(kn), acc))) THEN 
                  fail = .TRUE.
                  DEALLOCATE (xh, th)
                  EXIT
               END IF
            END DO

         END IF
        
         IF (fail) THEN
            ! If fail then store this solution and move on to the next attempt
            xh = x
            th = t
            np = nn
         ELSE
            EXIT ! If the integration was successful then exit
         END IF

      END DO

      IF (fail) THEN
         STOP 'ODE1_ADAPTIVE: Error, failed to converge'
      ELSE
         xh = x
         th = t
         DEALLOCATE(x, t)
         ALLOCATE(x(n), t(n))
         DO ip = 1, n
            in = 1+(ip-1)*(nn-1)/(n-1)
            x(ip) = xh(in)
            t(ip) = th(in)
         END DO
      END IF

   END SUBROUTINE ODE1_adaptive

   SUBROUTINE ODE2_adaptive(x, v, t, ti, tf, xi, vi, fv, n, iode, acc, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
      ! acc is the desired accuracy across the entire solution
      ! time steps are increased until convergence is achieved
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      REAL, INTENT(IN) :: acc
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      REAL, ALLOCATABLE :: xh(:), vh(:), th(:)
      INTEGER :: j, ip, in, np, nn
      LOGICAL :: fail
      INTEGER, PARAMETER :: jmax = jmax_adapt ! Maximum number of goes
      INTEGER, PARAMETER :: nmin = nmin_adapt ! Minimum number of points

      ! v' = fv
      INTERFACE       
         FUNCTION fv(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fv
      END INTERFACE

      IF (n <= 1) STOP 'ODE2_ADAPTIVE: Error, you cannot solve an ODE with n=1 (or less)'

      ! Loop over attemps for number of time steps
      DO j = 1, jmax

         ! Set the new number of time steps and save the old number
         nn = 1+(n-1)*2**(j-1)

         ! Solve the ODE with 'n' points
         CALL ODE2(x, v, t, ti, tf, xi, vi, fv, nn, iode, ilog)

         ! Either automatically fail or check for convergence with respect to previous solution
         ! TODO: Replace fail conditions with requal?
         IF (j == 1 .OR. nn < nmin) THEN
            fail = .TRUE.
         ELSE
            fail = .FALSE.
            DO ip = 1, np       
               in = 2*ip-1
               IF ((.NOT. almost_equal(xh(ip), x(in), acc)) .OR. (.NOT. almost_equal(vh(ip), v(in), acc))) THEN
               !IF ((.NOT. requal(xh(k), x(kn), acc)) .OR. (.NOT. requal(vh(k), v(kn), acc))) THEN 
                  fail = .TRUE.
                  DEALLOCATE (xh, vh, th)
                  EXIT
               END IF
            END DO

         END IF
        
         IF (fail) THEN
            ! If fail then store this solution and move on to the next attempt
            xh = x
            vh = v
            th = t
            np = nn
         ELSE
            EXIT ! If the integration was successful then exit
         END IF

      END DO

      IF (fail) THEN
         STOP 'ODE2_ADAPTIVE: Error, failed to converge'
      ELSE
         xh = x
         vh = v
         th = t
         DEALLOCATE(x, v, t)
         ALLOCATE(x(n), v(n), t(n))
         DO ip = 1, n
            in = 1+(ip-1)*(nn-1)/(n-1)
            x(ip) = xh(in)
            v(ip) = vh(in)
            t(ip) = th(in)
         END DO
      END IF

   END SUBROUTINE ODE2_adaptive

   SUBROUTINE ODEcoupled_adaptive(x, v, t, ti, tf, xi, vi, fx, fv, n, iode, acc, ilog)

      ! Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
      ! acc is the desired accuracy across the entire solution
      ! time steps are increased until convergence is achieved
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: t(:)
      REAL, INTENT(IN) :: ti
      REAL, INTENT(IN) :: tf
      REAL, INTENT(IN) :: xi
      REAL, INTENT(IN) :: vi
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iode
      REAL, INTENT(IN) :: acc
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog
      REAL, ALLOCATABLE :: xh(:), vh(:), th(:)
      INTEGER :: j, ip, in, np, nn
      LOGICAL :: fail
      INTEGER, PARAMETER :: jmax = jmax_adapt ! Maximum number of attempts
      INTEGER, PARAMETER :: nmin = nmin_adapt ! Minimum number of points

      ! x' = fx; v' = fv
      INTERFACE       
         FUNCTION fx(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fx
         FUNCTION fv(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fv
      END INTERFACE

      IF (n <= 1) STOP 'ODECOUPLED_ADAPTIVE: Error, you cannot solve an ODE with n=1 (or less)'

      ! Loop over attemps for number of time steps
      DO j = 1, jmax

         ! Set the number of time steps
         nn = 1+(n-1)*2**(j-1)

         ! Solve the ODE with 'nn' time steps
         CALL ODEcoupled(x, v, t, ti, tf, xi, vi, fx, fv, nn, iode, ilog)

         ! Automatically fail on the first go, otherwise check with previous attempt for convergence
         ! TODO: Replace fail conditions with requal?
         IF (j == 1 .OR. nn < nmin) THEN      
            fail = .TRUE.
         ELSE
            fail = .FALSE.
            DO ip = 1, np       
               in = 2*ip-1 ! kn is k-new      
               IF ((.NOT. almost_equal(xh(ip), x(in), acc)) .OR. (.NOT. almost_equal(vh(ip), v(in), acc))) THEN
               !IF ((.NOT. requal(xh(k), x(kn), acc)) .OR. (.NOT. requal(vh(k), v(kn), acc))) THEN 
                  fail = .TRUE.
                  DEALLOCATE (xh, vh, th)
                  EXIT
               END IF
            END DO
         END IF
        
         IF (fail) THEN
            ! If fail then store this solution and move on to the next attempt
            xh = x
            vh = v
            th = t
            np = nn
         ELSE
            EXIT ! If the integration was successful then exit
         END IF

      END DO

      IF (fail) THEN
         STOP 'ODECOUPLED_ADAPTIVE: Error, failed to converge'
      ELSE
         xh = x
         vh = v
         th = t
         DEALLOCATE(x, v, t)
         ALLOCATE(x(n), v(n), t(n))
         DO ip = 1, n
            in = 1+(ip-1)*(nn-1)/(n-1)
            x(ip) = xh(in)
            v(ip) = vh(in)
            t(ip) = th(in)
         END DO
      END IF

   END SUBROUTINE ODEcoupled_adaptive

   SUBROUTINE ODE1_advance(x1, x2, t1, t2, fx, iode)

      ! Advances the ODE system from t1 to t2, updating x1 to x2 and v1 to v2
      REAL, INTENT(IN) :: x1
      REAL, INTENT(OUT) :: x2
      REAL, INTENT(IN) :: t1
      REAL, INTENT(IN) :: t2
      REAL, EXTERNAL :: fx
      INTEGER, INTENT(IN) :: iode ! ODE solving method
      REAL :: x, t, dt
      REAL :: kx1, kx2, kx3, kx4

      ! fx is what dx/dt is equal to=
      INTERFACE     
         FUNCTION fx(x_in, t_in)
            REAL, INTENT(IN) :: x_in, t_in
         END FUNCTION fx
      END INTERFACE

      ! Set x, v, t to be the initial state of the system
      x = x1
      t = t1

      ! Calculate dt between the final and intial state
      dt = t2-t1
   
      IF (iode == iode_crude) THEN

         ! Crude method
         kx1 = dt*fx(x, t)

         x2 = x1+kx1
       
      ELSE IF (iode == iode_mid) THEN

         ! Mid-point method
         kx1 = dt*fx(x, t)
         kx2 = dt*fx(x+kx1/2., t+dt/2.)

         x2 = x1+kx2
        
      ELSE IF (iode == iode_RK4) THEN

         ! 4th order Rungge-Kutta
         kx1 = dt*fx(x, t)
         kx2 = dt*fx(x+kx1/2., t+dt/2.)
         kx3 = dt*fx(x+kx2/2., t+dt/2.)
         kx4 = dt*fx(x+kx3, t+dt)

         x2 = x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.

      ELSE

         STOP 'ODE_ADVANCE: Error, iode specified incorrectly'

      END IF

   END SUBROUTINE ODE1_advance

   SUBROUTINE ODE2_advance(x1, x2, v1, v2, t1, t2, fx, fv, iode)

      ! Advances the ODE system from t1 to t2, updating x1 to x2 and v1 to v2
      REAL, INTENT(IN) :: x1
      REAL, INTENT(OUT) :: x2
      REAL, INTENT(IN) :: v1
      REAL, INTENT(OUT) :: v2
      REAL, INTENT(IN) :: t1
      REAL, INTENT(IN) :: t2
      REAL, EXTERNAL :: fx
      REAL, EXTERNAL :: fv
      INTEGER, INTENT(IN) :: iode ! ODE solving method
      REAL :: x, v, t, dt
      REAL :: kx1, kx2, kx3, kx4
      REAL :: kv1, kv2, kv3, kv4

      ! x' = fx; v' = fv
      INTERFACE
         FUNCTION fx(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fx
         FUNCTION fv(x_in, v_in, t_in)
            REAL, INTENT(IN) :: x_in, v_in, t_in
         END FUNCTION fv
      END INTERFACE

      ! Set x, v, t to be the initial state of the system
      x = x1
      v = v1
      t = t1

      ! Calculate dt between the final and intial state
      dt = t2-t1
   
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

         x2 = x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.
         v2 = v1+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.

      ELSE

         STOP 'ODE_ADVANCE: Error, iode specified incorrectly'

      END IF

   END SUBROUTINE ODE2_advance

   REAL FUNCTION velocity(x, v, t)

      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: v
      REAL, INTENT(IN) :: t
      REAL :: crap

      ! Prevent compile-time warnings
      crap = x
      crap = t

      velocity = v

   END FUNCTION velocity

   LOGICAL FUNCTION almost_equal(x, y, eps)

      ! TODO: Surely requal can be used instead of this?
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: eps

      IF ((x > eps) .AND. (y > eps) .AND. (abs(-1.+x/y) > eps)) THEN
         almost_equal = .FALSE.
      ELSE
         almost_equal = .TRUE.
      END IF

   END FUNCTION almost_equal

END MODULE ODE_solvers
