MODULE vectors

  REAL, PARAMETER :: xhat(3)=[1.,0.,0.]
  REAL, PARAMETER :: yhat(3)=[0.,1.,0.]
  REAL, PARAMETER :: zhat(3)=[0.,0.,1.]

CONTAINS

  FUNCTION unit(x,n)

    !Returns the unit vector of 'x'
    IMPLICIT NONE
    REAL :: unit(n)
    REAL, INTENT(IN) :: x(n)
    INTEGER, INTENT(IN) :: n

    unit=x/modulus(x,n)

  END FUNCTION unit

  FUNCTION modulus(x,n)

    !Returns the modulus of vector 'x'
    IMPLICIT NONE
    REAL :: modulus
    REAL, INTENT(IN) :: x(n)
    INTEGER, INTENT(IN) :: n
    
    modulus=sqrt(SUM(x**2))

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

  FUNCTION cross_product(x,y)

    !Makes the cross produce of 3-vectors x nad y
    IMPLICIT NONE
    REAL :: cross_product(3)
    REAL, INTENT(IN) :: x(3), y(3)

    cross_product(1)=x(2)*y(3)-x(3)*y(2)
    cross_product(2)=x(3)*y(1)-x(1)*y(3)
    cross_product(3)=x(1)*y(2)-x(2)*y(1)

  END FUNCTION cross_product

  FUNCTION rotation(u,t)

    !Creates a rotation matrix for rotation about vector 'u' by angle 't'
    !u must be a unit vector!!!
    IMPLICIT NONE
    REAL :: rotation(3,3)
    REAL, INTENT(IN) :: u(3), t
    REAL :: s, c

    !Pre-compute sine and cosine for speed
    s=sin(t)
    c=cos(t)
    
    rotation(1,1)=c+(1.-c)*u(1)**2.
    rotation(1,2)=u(1)*u(2)*(1.-c)-u(3)*s
    rotation(1,3)=u(1)*u(3)*(1.-c)+u(2)*s

    rotation(2,1)=u(2)*u(1)*(1.-c)+u(3)*s
    rotation(2,2)=c+(1.-c)*u(2)**2.
    rotation(2,3)=u(2)*u(3)*(1.-c)-u(1)*s

    rotation(3,1)=u(3)*u(1)*(1.-c)-u(2)*s
    rotation(3,2)=u(3)*u(2)*(1.-c)+u(1)*s
    rotation(3,3)=c+(1.-c)*u(3)**2.

  END FUNCTION rotation

  FUNCTION matrix_multiply(A,B,n)

    !Multiplies matrices 'A' and 'B' C=A*B
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL :: matrix_multiply(n,n)
    REAL, INTENT(IN) :: A(n,n), B(n,n)
    REAL :: C(n,n)
    INTEGER :: i, j, k

    !Fix this 'sum' variable to zero
    C=0.

    DO i=1,n
       DO j=1,n
          DO k=1,n
             C(i,j)=C(i,j)+A(i,k)*B(k,j)
          END DO
       END DO
    END DO

    matrix_multiply=C

  END FUNCTION matrix_multiply

  FUNCTION matrix_vector(A,b,n)

    !Multiplies matrix a and vector b (c=A*b)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL :: matrix_vector(n)
    REAL, INTENT(IN) :: A(n,n), b(n)
    REAL :: c(n)
    INTEGER :: i, j

    !Fix this 'sum' variable to zero
    c=0.

    DO i=1,n
       DO j=1,n
          c(i)=c(i)+A(i,j)*b(j)
       END DO
    END DO

    matrix_vector=c

  END FUNCTION matrix_vector

  FUNCTION distance(x1,x2,n)

    !Calculates the distance between n-vectors x1 and x2
    IMPLICIT NONE
    REAL :: distance
    REAL, INTENT(IN) :: x1(n), x2(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    !Fix this 'sum' variable to zero
    distance=0.

    DO i=1,n
       distance=distance+(x1(i)-x2(i))**2
    END DO

    distance=sqrt(distance)

  END FUNCTION distance

  FUNCTION rotate_vector(v,k,theta)

    !Rotates vector 'v' about axis 'k' by angle 'theta'
    IMPLICIT NONE
    REAL :: rotate_vector(3)
    REAL, INTENT(IN) :: v(3), k(3), theta
    REAL :: R(3,3)

    !Get the rotation matrix
    !Infact, I think you do not need the whole matrix to do this
    !Look up Rodriguez Formula in case this is important
    !Probably it is faster
    R=rotation(k,theta)

    !Perform the rotation
    rotate_vector=matrix_vector(R,v,3)
    
  END FUNCTION rotate_vector

  FUNCTION rotate_vector_fast(v,k,theta)

    !Rotates vector 'v' about axis 'k' by angle 'theta'
    !Note that 'k' needs to be a unit vector
    !See https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    IMPLICIT NONE
    REAL :: rotate_vector_fast(3)
    REAL, INTENT(IN) :: v(3), k(3), theta
    REAL :: r1(3), r2(3), r3(3), s, c

    s=sin(theta)
    c=cos(theta)

    r1=v*c
    r2=cross_product(k,v)*s
    r3=k*(dot_product(k,v))*(1.-c)
    
    rotate_vector_fast=r1+r2+r3
    
  END FUNCTION rotate_vector_fast

END MODULE vectors
