PROGRAM logical_operations_test

  USE basic_operations

  IMPLICIT NONE
  LOGICAL :: fails

  WRITE(*,*)
  
  CALL test_first_digit(fails)
  CALL test_swap(fails)
  CALL test_positive_negative(fails)
  CALL test_even_odd(fails)
  
CONTAINS

  SUBROUTINE test_first_digit(fail)

    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: fail
    REAL, ALLOCATABLE :: x(:)
    INTEGER, ALLOCATABLE :: f(:)
    INTEGER :: i, j
    INTEGER, PARAMETER :: n=5
    
    fail=.FALSE.    

    ALLOCATE(x(n),f(n))
    
    x(1)=0.00123
    f(1)=1
    
    x(2)=2.3
    f(2)=2
    
    x(3)=3.5e6
    f(3)=3
    
    x(4)=0.0000567
    f(4)=5
    
    x(5)=24.7
    f(5)=2   

    DO i=1,n
       j=first_digit(x(i))
       IF(j .NE. f(i)) THEN
          WRITE(*,*) 'TEST_FIRST_DIGIT: Fail:', x(i), f(i)
          fail=.TRUE.
          STOP
       END IF
    END DO

    IF(.NOT. fail) THEN
       WRITE(*,*) 'TEST_FIRST_DIGIT: Pass'
       WRITE(*,*)
    END IF
    
  END SUBROUTINE test_first_digit

  SUBROUTINE test_swap(fail)

    ! TODO: Add integer test
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: fail    
    REAL, ALLOCATABLE :: a(:), b(:)
    REAL :: aa, bb
    INTEGER :: i
    INTEGER, PARAMETER :: n=3
    
    fail=.FALSE.

    ALLOCATE(a(n),b(n))
    
    a(1)=1.
    b(1)=2.

    a(2)=44.
    b(2)=44.

    a(3)=26.
    b(3)=7.

    DO i=1,n
       aa=a(i)
       bb=b(i)
       CALL swap(aa,bb)
       IF((aa .NE. b(i)) .OR. (bb .NE. a(i))) THEN
          WRITE(*,*) 'TEST_SWAP: Fail:', a(i), b(i)
          fail=.TRUE.
          STOP
       END IF
    END DO

    IF(.NOT. fail) THEN
       WRITE(*,*) 'TEST_SWAP: Pass'
       WRITE(*,*)
    END IF

  END SUBROUTINE test_swap

  SUBROUTINE test_positive_negative(fail)

    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: fail    
    REAL, ALLOCATABLE :: test(:)
    LOGICAL, ALLOCATABLE :: ans(:)
    LOGICAL :: pos, neg
    INTEGER :: i
    INTEGER, PARAMETER :: n=5

    fail=.FALSE.

    ALLOCATE(test(n),ans(n))

    test(1)=0.
    ans(1)=.TRUE.
    
    test(2)=1.
    ans(2)=.TRUE.

    test(3)=-1.
    ans(3)=.FALSE.

    test(4)=-1000.
    ans(4)=.FALSE.

    test(5)=1e18
    ans(5)=.TRUE.

    DO i=1,n
       pos=positive(test(i))
       neg=negative(test(i))
       IF(pos .EQV. neg)     fail=.TRUE.
       IF(pos .NEQV. ans(i)) fail=.TRUE.
       IF(neg .EQV. ans(i))  fail=.TRUE.
       IF(fail) THEN
          WRITE(*,*) 'TEST_POSITIVE_NEGATIVE: Fail:', i, test(i), ans(i), pos, neg
          STOP
       END IF
    END DO
    
    IF(.NOT. fail) THEN
       WRITE(*,*) 'TEST_POSITIVE_NEGATIVE: Pass'
       WRITE(*,*)
    END IF

  END SUBROUTINE test_positive_negative

  SUBROUTINE test_even_odd(fail)

    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: fail    
    INTEGER, ALLOCATABLE :: test(:)
    LOGICAL, ALLOCATABLE :: ans(:)
    LOGICAL :: even_test, odd_test
    INTEGER :: i
    INTEGER, PARAMETER :: n=5

    fail=.FALSE.

    ALLOCATE(test(n),ans(n))

    test(1)=0
    ans(1)=.TRUE.
    
    test(2)=2
    ans(2)=.TRUE.

    test(3)=-1
    ans(3)=.FALSE.

    test(4)=-1000
    ans(4)=.TRUE.

    test(5)=10101022
    ans(5)=.TRUE.

    DO i=1,n
       even_test=even(test(i))
       odd_test=odd(test(i))
       IF(even_test .EQV.  odd_test) fail=.TRUE.
       IF(even_test .NEQV. ans(i))   fail=.TRUE.
       IF(odd_test  .EQV.  ans(i))   fail=.TRUE.
       IF(fail) THEN
          WRITE(*,*) 'TEST_EVEN_ODD: Fail:', test(i), ans(i), even_test, odd_test
          STOP
       END IF
    END DO
    
    IF(.NOT. fail) THEN
       WRITE(*,*) 'TEST_EVEN_ODD: Pass'
       WRITE(*,*)
    END IF

  END SUBROUTINE test_even_odd

END PROGRAM logical_operations_test
