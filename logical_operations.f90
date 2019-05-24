MODULE logical_operations

  ! TODO: Rename logical_operations -> basic_operations

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fix_min
  PUBLIC :: fix_max
  PUBLIC :: read_command_argument
  PUBLIC :: positive
  PUBLIC :: negative
  PUBLIC :: even
  PUBLIC :: odd
  PUBLIC :: requal
  PUBLIC :: present_and_correct
  PUBLIC :: first_digit
  PUBLIC :: swap

  INTERFACE swap
     MODULE PROCEDURE swap_real
     MODULE PROCEDURE swap_int
  END INTERFACE swap

  INTERFACE read_command_argument
     MODULE PROCEDURE read_command_argument_real
     MODULE PROCEDURE read_command_argument_integer
     MODULE PROCEDURE read_command_argument_character
  END INTERFACE read_command_argument

CONTAINS

  SUBROUTINE fix_min(x,xmin)

    ! If x is below xmin then set to xmin
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x ! Value to fix
    REAL, INTENT(IN) :: xmin ! Minimum value for x

    IF(x<xmin) x=xmin
    
  END SUBROUTINE fix_min

  SUBROUTINE fix_max(x,xmax)

    ! If x is above xmax then set to xmax
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x ! Value to fix
    REAL, INTENT(IN) :: xmax ! Maximum value for x

    IF(x>xmax) x=xmax
    
  END SUBROUTINE fix_max

  SUBROUTINE read_command_argument_real(i,x,desc)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i             ! Position of command-line argument
    REAL, INTENT(OUT) :: x               ! Real number to be assigned command-line argument
    CHARACTER(len=*), INTENT(IN) :: desc ! Description of command-line argument
    CHARACTER(len=256) :: word

    CALL get_command_argument(i,word)
    IF(word=='') THEN
       WRITE(*,*) 'READ_COMMAND_ARGUMENT_REAL: Missing command-line argument number:', i
       WRITE(*,*) 'READ_COMMAND_ARGUMENT_REAL: Missing command-line argument: ', TRIM(desc)
       STOP
    ELSE
       READ(word,*) x
    END IF

  END SUBROUTINE read_command_argument_real

  SUBROUTINE read_command_argument_integer(i,x,desc)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i             ! Position of command-line argument
    INTEGER, INTENT(OUT) :: x            ! Real number to be assigned command-line argument
    CHARACTER(len=*), INTENT(IN) :: desc ! Description of command-line argument
    CHARACTER(len=256) :: word

    CALL get_command_argument(i,word)
    IF(word=='') THEN
       WRITE(*,*) 'READ_COMMAND_ARGUMENT_INTEGER: Missing command-line argument:', i
       WRITE(*,*) 'READ_COMMAND_ARGUMENT_INTEGER: Missing command-line argument: ', TRIM(desc)
       STOP
    ELSE
       READ(word,*) x
    END IF

  END SUBROUTINE read_command_argument_integer

  SUBROUTINE read_command_argument_character(i,x,desc)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i             ! Position of command-line argument
    CHARACTER(len=256), INTENT(OUT) :: x ! String to be assigned command-line argument
    CHARACTER(len=*), INTENT(IN) :: desc ! Description of command-line argument
    CHARACTER(len=256) :: word

    CALL get_command_argument(i,word)
    IF(word=='') THEN
       WRITE(*,*) 'READ_COMMAND_ARGUMENT_CHARACTER: Missing command-line argument:', i
       WRITE(*,*) 'READ_COMMAND_ARGUMENT_CHARACTER: Missing command-line argument: ', TRIM(desc)
       STOP
    ELSE
       x=word
    END IF

  END SUBROUTINE read_command_argument_character

  LOGICAL FUNCTION positive(x)

    ! Logical function that returns .TRUE. if x>=0.
    ! If x=0. then it will return true
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    IF(x<0.) THEN
       positive=.FALSE.
    ELSE
       positive=.TRUE.
    END IF
    
  END FUNCTION positive

  LOGICAL FUNCTION negative(x)

    ! Returns true if argument is negative
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    IF(positive(x)) THEN
       negative=.FALSE.
    ELSE
       negative=.TRUE.
    END IF
    
  END FUNCTION negative

  LOGICAL FUNCTION even(i)

    ! Tests for i being even, returns true if even
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i

    IF(mod(i,2)==0) THEN
       even=.TRUE.
    ELSE
       even=.FALSE.
    END IF
    
  END FUNCTION even

  LOGICAL FUNCTION odd(i)

    ! Tests for i being odd, returns true if odd
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i  

    IF(even(i)) THEN
       odd=.FALSE.
    ELSE
       odd=.TRUE.
    END IF

  END FUNCTION odd
  
  LOGICAL FUNCTION requal(x,y,eps)

    ! Tests if two REAL numbers are within some tolerance of each other
    ! If they are, then they should be considered equal
    ! Adapted from https://stackoverflow.com/questions/4915462/how-should-i-do-floating-point-comparison/4915891#4915891
    IMPLICIT NONE
    REAL, INTENT(IN) :: x, y
    REAL, INTENT(IN) :: eps
    REAL :: absx, absy, diff

    absx=abs(x)
    absy=abs(y)
    diff=abs(x-y)

    IF(x==y) THEN
       requal=.TRUE.
    ELSE IF(x==0. .OR. y==0. .OR. diff<tiny(x)) THEN
       IF(diff<eps*tiny(x)) THEN
          requal=.TRUE.
       ELSE
          requal=.FALSE.
       END IF
    ELSE
       IF(diff/(absx+absy)<eps) THEN
          requal=.TRUE.
       ELSE
          requal=.FALSE.
       END IF
    END IF
    
  END FUNCTION requal

  LOGICAL FUNCTION present_and_correct(x)

    IMPLICIT NONE
    LOGICAL, OPTIONAL, INTENT(IN) :: x

    IF(present(x)) THEN
       IF(x) THEN
          present_and_correct=.TRUE.
       ELSE
          present_and_correct=.FALSE.
       END IF
    ELSE
       present_and_correct=.FALSE.
    END IF
    
  END FUNCTION present_and_correct

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

    ! Swaps the values of variables a and b
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

END MODULE logical_operations
