PROGRAM logical_operations_test

  USE logical_operations

  IMPLICIT NONE
  INTEGER :: i, n
  INTEGER :: test
  LOGICAL :: lresult

  WRITE(*,*) '1 - Test odd'
  WRITE(*,*) '2 - Test positive'
  READ(*,*) test
  WRITE(*,*)

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

END PROGRAM logical_operations_test
