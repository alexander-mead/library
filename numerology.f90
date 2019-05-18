MODULE numerology

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: first_digit
  PUBLIC :: swap

  INTERFACE swap
     MODULE PROCEDURE swap_real
     MODULE PROCEDURE swap_int
  END INTERFACE swap

CONTAINS

  INTEGER FUNCTION first_digit(x)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL :: y

    y=abs(x)

    DO
       IF(y==1.) THEN
          first_digit=1
          EXIT
       ELSE IF(y>1. .AND. y<10.) THEN
          first_digit=floor(y)
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

!!$  SUBROUTINE swap_int(a,b)
!!$
!!$    !Swaps the values of integers a and b
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(INOUT) :: a, b
!!$    INTEGER :: c    
!!$
!!$    c=a
!!$    a=b
!!$    b=c
!!$    
!!$  END SUBROUTINE swap_int

  SUBROUTINE swap_int(n,m)

    ! Swap integers n and m in the most memory-efficient way possible
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: n
    INTEGER, INTENT(INOUT) :: m

    n=n+m ! n' = n+m
    m=n-m ! m' = n'-m = n+m-m = n
    n=n-m ! n'' = n'-m' = n+m-n = m
    
  END SUBROUTINE swap_int
  
END MODULE numerology
