MODULE numerology

CONTAINS

  FUNCTION first_digit(x)

    IMPLICIT NONE
    INTEGER :: first_digit
    REAL, INTENT(IN) :: x
    REAL :: y

    y=ABS(x)

    DO
       IF(y==1.) THEN
          first_digit=1
          EXIT
       ELSE IF(y>1. .AND. y<10.) THEN
          first_digit=FLOOR(y)
          EXIT
       ELSE IF(y>=10.) THEN
          y=y/10.
       ELSE IF(y<1.) THEN
          y=y*10.
       END IF
    END DO

  END FUNCTION first_digit

  SUBROUTINE swap_real(a,b)

    !Swaps the values of variables a and b
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a, b
    REAL :: c

    c=a
    a=b
    b=c

  END SUBROUTINE swap_real

  SUBROUTINE swap_int(a,b)

    !Swaps the values of variables a and b
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: a, b
    INTEGER :: c    

    c=a
    a=b
    b=c
    
  END SUBROUTINE swap_int
  
END MODULE numerology
