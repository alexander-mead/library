MODULE solve_equations

  ! TODO: Write a wrapper for solve routines

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: find_solve
  PUBLIC :: bisect_solve
  
CONTAINS

  REAL FUNCTION find_solve(a,xtab,ytab,n)

    ! Solves y(x)=a for x
    ! TODO: Change so that a=0.
    USE interpolate
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    REAL, INTENT(IN) :: xtab(n)
    REAL, INTENT(IN) :: ytab(n)
    INTEGER, INTENT(IN) :: n

    find_solve=find(a,ytab,xtab,n,3,3,2)
    
  END FUNCTION find_solve

  REAL FUNCTION bisect_solve(xtab,ytab,n,acc)   

    ! Solves y(x)=0 for x, f(x) should be monotonic and cross f=0. once only
    USE logical_operations
    USE interpolate
    IMPLICIT NONE
    REAL, INTENT(IN) :: xtab(n)
    REAL, INTENT(IN) :: ytab(n)
    INTEGER, INTENT(IN) :: n
    REAL :: x1, x2, y1, y2, x, y
    REAL, INTENT(IN) :: acc
    INTEGER :: i

    ! Initial values taken from top and bottom of table
    x1=xtab(1)
    x2=xtab(n)
    y1=ytab(1)
    y2=ytab(n)

    ! Now iterate until desired accuracy is reached
    i=0
    DO
       i=i+1
       x=0.5*(x1+x2)
       y=find(x,xtab,ytab,n,3,3,2)
       IF(abs(y)<acc) THEN
          EXIT
       ELSE IF(positive(y1) .EQV. positive(y)) THEN
          x1=x
          y1=y
       ELSE
          x2=x
          y2=y
       END IF
    END DO

    bisect_solve=x
    
  END FUNCTION bisect_solve

END MODULE solve_equations