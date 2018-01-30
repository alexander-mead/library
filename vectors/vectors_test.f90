PROGRAM vectors_test

  USE vectors
  USE constants

  !Some functions to manipulate 3D Cartesian vectors
  IMPLICIT NONE
  REAL :: a(3), b(3), c(3), R(3,3), theta
  REAL :: xhat(3), yhat(3), zhat(3)
  !REAL, PARAMETER :: pi=3.141592654

  a(1)=1.
  a(2)=2.
  a(3)=3.

  b(1)=1.
  b(2)=0.
  b(3)=2.

  xhat(1)=1.
  xhat(2)=0.
  xhat(3)=0.

  yhat(1)=0.
  yhat(2)=1.
  yhat(3)=0.

  zhat(1)=0.
  zhat(2)=0.
  zhat(3)=1.

  c=0.

  WRITE(*,*)
  WRITE(*,*) 'Vector a:', a(1), a(2), a(3)
  WRITE(*,*) 'Vector b:', b(1), b(2), b(3)
  WRITE(*,*)

  c=unit(a,3)
  WRITE(*,*) 'Unit a:', c(1), c(2), c(3)
  c=unit(b,3)
  WRITE(*,*) 'Unit b:', c(1), c(2), c(3)

  WRITE(*,*) 'Mod(a):', modulus(a,3)
  WRITE(*,*) 'Mod(b):', modulus(b,3)

  WRITE(*,*) 'a.b:', dot_product(a,b)

  c=cross_product(a,b)

  WRITE(*,*) 'a x b:', c(1), c(2), c(3)
  WRITE(*,*)

  theta=pi/2.
  WRITE(*,*) 'Rotation angle:', theta
  WRITE(*,*) 'Rotation axis b:', b(1), b(2), b(3)
  R=rotation(b,theta)
  c=matrix_vector(R,a,3)

  WRITE(*,*) 'Rotation of a about axis with angle:', c(1), c(2), c(3)
  WRITE(*,*) 'Modulus after rotation:', modulus(a,3)
  WRITE(*,*)

END PROGRAM vectors_test
