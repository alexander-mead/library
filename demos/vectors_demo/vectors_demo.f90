PROGRAM vectors_test

  USE vectors
  USE constants
  USE array_operations

  !Some functions to manipulate 3D Cartesian vectors
  IMPLICIT NONE
  CHARACTER(len=256) :: test
  REAL :: a(3), b(3), c(3), R(3,3), theta
  INTEGER :: itest, i, n
  REAL :: tmin, tmax, t, ts

  CALL get_command_argument(1,test)
  IF(test=='') THEN
     itest=-1
  ELSE
     READ(test,*) itest
  END IF

  IF(itest==-1) THEN
     WRITE(*,*) 'VECTORS_TEST: Choose test'
     WRITE(*,*) 'VECTORS_TEST: 1 - Standard test'
     WRITE(*,*) 'VECTORS_TEST: 2 - Shift to circle test'
     READ(*,*) itest
     WRITE(*,*)
  END IF

  IF(itest==1) THEN

     a(1)=1.
     a(2)=2.
     a(3)=3.

     b(1)=1.
     b(2)=0.
     b(3)=2.
     
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

  ELSE IF(itest==2) THEN

     tmin=-2.*pi
     tmax=4.*pi
     n=51

     DO i=1,n
        t=progression(tmin,tmax,i,n)
        ts=shift_angle_to_circle(t)
        WRITE(*,*) t, ts, t*rad2deg, ts*rad2deg
     END DO

  ELSE

     STOP 'VECTORS_TEST: Error, test specified incorrectly'

  END IF

END PROGRAM vectors_test
