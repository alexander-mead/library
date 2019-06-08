PROGRAM logical_operations_demo

  USE array_operations
  USE logical_operations
  USE random_numbers

  IMPLICIT NONE
  INTEGER :: i, n, nx
  INTEGER :: test
  LOGICAL :: lresult
  REAL :: xmin, xmax, x, x0, eps, y

  ! First digit test
  INTEGER, PARAMETER :: iseed=0    ! Seed for RNG
  INTEGER, PARAMETER :: n_first=10 ! Number to test
  REAL, PARAMETER :: r1=0.         ! Lower limit for uniform random
  REAL, PARAMETER :: r2=1.         ! Upper limit for uniform random

  WRITE(*,*)

  WRITE(*,*) 'Choose test'
  WRITE(*,*) '1 - Test odd'
  WRITE(*,*) '2 - Test positive'
  WRITE(*,*) '3 - Test requal'
  WRITE(*,*) '4 - Test first digit'
  READ(*,*) test
  WRITE(*,*)

  IF(test==1 .OR. test==2) THEN

     n=10

     DO i=-n,n

        IF(test==1) THEN
           lresult=odd(i)
        ELSE IF(test==2) THEN
           lresult=positive(REAL(i))
        ELSE
           STOP 'LOGICAL_OPERATIONS_TEST: Error, test not specified correctly'
        END IF

        WRITE(*,*) i, lresult

     END DO

     WRITE(*,*)

  ELSE IF(test==3) THEN

     x0=200.
     xmin=199.
     xmax=201.
     nx=201

     eps=1e-4

     DO i=1,nx
        x=progression(xmin,xmax,i,nx)
        IF(requal(x,x0,eps)) THEN
           WRITE(*,*) 'According to requal', x, 'equals', x0
        END IF
     END DO

  ELSE IF(test==4) THEN

     CALL RNG_set(iseed)

     WRITE(*,*) '        i         x      y  first_digit'
     WRITE(*,*) '======================================='
     DO i=1,n_first
        x=random_uniform(r1,r2)
        y=log(x)
        WRITE(*,fmt='(I10,2F10.5,I10)') i, x, y, first_digit(y)
     END DO
     WRITE(*,*) '======================================='
     WRITE(*,*)
     
  ELSE

     STOP 'LOGICAL_OPERATIONS_TEST: Error, test specified incorrectly'

  END IF

END PROGRAM logical_operations_demo
