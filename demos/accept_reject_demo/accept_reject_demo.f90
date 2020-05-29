PROGRAM accept_reject_demo

  USE random_numbers
  
  IMPLICIT NONE
  INTEGER :: i
  CHARACTER(len=256) :: outfile
  REAL, ALLOCATABLE :: x(:)
  REAL, PARAMETER :: x1=0.
  REAL, PARAMETER :: x2=1.
  INTEGER, PARAMETER :: n=100000
  INTEGER, PARAMETER :: iseed=0

  ! Set random number generator
  CALL RNG_seed(iseed)

  ! Allocate array for draws
  ALLOCATE(x(n))

  ! Make draws and write results
  outfile='results.dat'
  OPEN(7,file=outfile)
  DO i=1,n
     x(i)=accept_reject(f,x1,x2,f(x2))
     WRITE(7,*) x(i)
  END DO
  CLOSE(7)
  
CONTAINS

  REAL FUNCTION f(xx)

    IMPLICIT NONE
    REAL, INTENT(IN) :: xx

    f=xx**2

  END FUNCTION f

END PROGRAM accept_reject_demo
