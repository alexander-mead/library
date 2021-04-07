MODULE vectors

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: shift_angle_to_circle
   PUBLIC :: unit
   PUBLIC :: modulus
   PUBLIC :: cross_product
   PUBLIC :: rotation
   PUBLIC :: matrix_multiply
   PUBLIC :: matrix_vector
   PUBLIC :: distance
   PUBLIC :: rotate_vector
   PUBLIC :: rotate_vector_fast
   PUBLIC :: trace
   PUBLIC :: determinant

   ! Unit vectors
   PUBLIC  :: xhat
   PUBLIC  :: yhat
   PUBLIC  :: zhat

   REAL, PARAMETER :: xhat(3) = [1., 0., 0.]
   REAL, PARAMETER :: yhat(3) = [0., 1., 0.]
   REAL, PARAMETER :: zhat(3) = [0., 0., 1.]

   !INTERFACE determinant
   !   PROCEDURE determinant_2
   !   PROCEDURE determinant_3
   !END INTERFACE determinant

CONTAINS

   REAL FUNCTION trace(A)

      REAL, INTENT(IN) :: A(:, :)
      INTEGER :: i

      IF(size(A, 1) /= size(A, 2)) STOP 'TRACE: Error, trace only defined for square matrices'
      trace = 0.
      DO i = 1, size(A)
         trace = trace+A(i, i)
      END DO

   END FUNCTION trace

   REAL FUNCTION determinant(A)

      ! TODO: Write a general recursive function
      REAL, INTENT(IN) :: A(:, :)
      INTEGER :: n

      n = size(A, 1)
      IF (n /= size(A,2)) STOP 'DETERMINANT: Error, determinant only defined for square matrices'
      IF (n == 1) THEN
         determinant = A(1, 1) ! NOTE: Can be negative, should not have abs
      ELSE IF (n == 2) THEN
         determinant = determinant_2(A)
      ELSE IF (n == 3) THEN
         determinant = determinant_3(A)
      ELSE
         STOP 'DETERMINANT: Function not written for this size of matix'
      END IF

   END FUNCTION determinant

   REAL FUNCTION determinant_2(A)

      REAL, INTENT(IN) :: A(2, 2)

      determinant_2 = A(1, 1)*A(2, 2)-A(1, 2)*A(2, 1)

   END FUNCTION determinant_2

   REAL FUNCTION determinant_3(A)

      ! TODO: Surely this can be much neater
      REAL, INTENT(IN) :: A(3, 3)
      REAL :: B(2, 2), M
      INTEGER :: i

      determinant_3 = 0.
      DO i = 1, 3
         IF (i == 1) THEN
            M = A(1, 1)
            B(1, 1) = A(2, 2)
            B(1, 2) = A(2, 3)
            B(2, 1) = A(3, 2)
            B(2, 2) = A(3, 3)
         ELSE IF (i == 2) THEN
            M = -A(1, 2)
            B(1, 1) = A(2, 1)
            B(1, 2) = A(2, 3)
            B(2, 1) = A(3, 1)
            B(2, 2) = A(3, 3)
         ELSE IF (i == 3) THEN
            M = A(1, 3)
            B(1, 1) = A(2, 1)
            B(1, 2) = A(2, 2)
            B(2, 1) = A(3, 1)
            B(2, 2) = A(3, 2)
         ELSE
            STOP 'DETERMINANT_3: Error, something went very wrong'
         END IF
         determinant_3 = determinant_3+M*determinant_2(B)
      END DO

   END FUNCTION determinant_3

   REAL FUNCTION shift_angle_to_circle(theta)

      ! Changes an angle thetta [rad] to be in the interval [0:2pi]
      ! e.g., if t is -pi/2 it gets shifted to 3pi/2
      USE constants
      REAL, INTENT(IN) :: theta

      shift_angle_to_circle = theta-twopi*floor(theta/twopi)

   END FUNCTION shift_angle_to_circle

   FUNCTION unit(x)

      ! Returns the unit vector of 'x'
      REAL, INTENT(IN) :: x(:)
      REAL :: unit(size(x))

      unit = x/modulus(x)

   END FUNCTION unit

   REAL FUNCTION modulus(x)

      ! Returns the modulus of vector 'x'
      REAL, INTENT(IN) :: x(:)

      modulus = sqrt(sum(x**2))

   END FUNCTION modulus

!!$  !Note - in F95 there is an intrinsic dot_product function
!!$  !Therefore this is now obsolete
!!$  FUNCTION dot_product(x,y)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: dot_product
!!$    REAL, INTENT(IN) :: x(3), y(3)
!!$
!!$    dot_product=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
!!$
!!$  END FUNCTION dot_product

   FUNCTION cross_product(x, y)

      ! Computes the cross produce of 3-vectors x and y
      REAL :: cross_product(3)
      REAL, INTENT(IN) :: x(3)
      REAL, INTENT(IN) :: y(3)

      cross_product(1) = x(2)*y(3)-x(3)*y(2)
      cross_product(2) = x(3)*y(1)-x(1)*y(3)
      cross_product(3) = x(1)*y(2)-x(2)*y(1)

   END FUNCTION cross_product

   FUNCTION rotation(u, t)

      ! Creates a rotation matrix for rotation about vector 'u' by angle 't'
      ! u must be a unit vector
      REAL :: rotation(3, 3)
      REAL, INTENT(IN) :: u(3)
      REAL, INTENT(IN) :: t
      REAL :: s, c

      ! Pre-compute sine and cosine for speed
      s = sin(t)
      c = cos(t)

      rotation(1, 1) = c+(1.-c)*u(1)**2.
      rotation(1, 2) = u(1)*u(2)*(1.-c)-u(3)*s
      rotation(1, 3) = u(1)*u(3)*(1.-c)+u(2)*s

      rotation(2, 1) = u(2)*u(1)*(1.-c)+u(3)*s
      rotation(2, 2) = c+(1.-c)*u(2)**2.
      rotation(2, 3) = u(2)*u(3)*(1.-c)-u(1)*s

      rotation(3, 1) = u(3)*u(1)*(1.-c)-u(2)*s
      rotation(3, 2) = u(3)*u(2)*(1.-c)+u(1)*s
      rotation(3, 3) = c+(1.-c)*u(3)**2.

   END FUNCTION rotation

   FUNCTION matrix_multiply(A, B) result(C)

      ! Multiplies matrices 'A' and 'B' C=A*B
      REAL, INTENT(IN) :: A(:, :), B(:, :)
      REAL :: C(size(A, 1), size(B, 2))
      INTEGER :: i, j, k, ni, nj, nk

      ni = size(A, 1)
      nj = size(B, 2)
      nk = size(A, 2)
      IF (nk /= size(B, 1)) STOP 'MATRIX_MULTIPLY: Error, you cannot multiply these matrices'

      ! Fix this 'sum' variable to zero
      C = 0.

      DO i = 1, ni
         DO j = 1, nj
            DO k = 1, nk
               C(i, j) = C(i, j)+A(i, k)*B(k, j)
            END DO
         END DO
      END DO

      !matrix_multiply = C

   END FUNCTION matrix_multiply

   FUNCTION matrix_vector(A, b) result(c)

      ! Multiplies matrix a and vector b (c=A*b) 
      REAL, INTENT(IN) :: A(:, :)
      REAL, INTENT(IN) :: b(:)  
      REAL :: c(size(A, 1))
      INTEGER :: ix, iy, nx, ny

      nx = size(A, 1)
      ny = size(A, 2)
      IF (ny /= size(b)) STOP 'MATRIX_VECTOR: Error, array wrong sizes'

      ! Fix this 'sum' variable to zero
      c = 0.

      DO ix = 1, nx
         DO iy = 1, ny
            c(ix) = c(ix)+A(ix, iy)*b(iy)
         END DO
      END DO

   END FUNCTION matrix_vector

   REAL FUNCTION distance(x1, x2)

      ! Calculates the scalar distance between vectors x1 and x2
      ! Note that this is always a positive number
      REAL, INTENT(IN) :: x1(:)
      REAL, INTENT(IN) :: x2(:)
      INTEGER :: i, n

      n = size(x1)
      IF (n /= size(x2)) STOP 'DISTANCE: Error, x1 and x2 must be the same size'

      ! Fix this 'sum' variable to zero
      distance = 0.

      DO i = 1, n
         distance = distance+(x1(i)-x2(i))**2
      END DO

      distance = sqrt(distance)

   END FUNCTION distance

   FUNCTION rotate_vector(v, k, theta)

      ! Rotates vector 'v' about axis 'k' by angle 'theta'
      REAL :: rotate_vector(3)
      REAL, INTENT(IN) :: v(3)
      REAL, INTENT(IN) :: k(3)
      REAL, INTENT(IN) :: theta
      REAL :: R(3, 3)

      ! Get the rotation matrix
      ! Infact, I think you do not need the whole matrix to do this
      ! Look up Rodriguez Formula in case this is important
      ! Probably it is faster
      R = rotation(k, theta)

      ! Perform the rotation
      rotate_vector = matrix_vector(R, v)

   END FUNCTION rotate_vector

   FUNCTION rotate_vector_fast(v, k, theta)

      ! Rotates vector 'v' about axis 'k' by angle 'theta'
      ! Note that 'k' needs to be a unit vector
      ! See https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
      REAL :: rotate_vector_fast(3)
      REAL, INTENT(IN) :: v(3)
      REAL, INTENT(IN) :: k(3)
      REAL, INTENT(IN) :: theta
      REAL :: r1(3), r2(3), r3(3), s, c

      s = sin(theta)
      c = cos(theta)

      r1 = v*c
      r2 = cross_product(k, v)*s
      r3 = k*(dot_product(k, v))*(1.-c)

      rotate_vector_fast = r1+r2+r3

   END FUNCTION rotate_vector_fast

END MODULE vectors
