PROGRAM numerology_test

  USE random_numbers
  USE numerology

  IMPLICIT NONE
  REAL :: x, y
  INTEGER :: i, n

  CALL RNG_set(0)

  n=100
  DO i=1,n
     x=uniform(0.,1.)
     y=log(x)
     WRITE(*,fmt='(I10,2F10.5,I10)') i, x, y, first_digit(y)
  END DO
  
END PROGRAM numerology_test
