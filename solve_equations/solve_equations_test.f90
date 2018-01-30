PROGRAM solve_equations_test

  USE array_operations

  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:), LHS(:), RHS(:)
  REAL :: xmin, xmax, sol
  INTEGER :: i, n

  !Set xtable and range
  xmin=0.
  xmax=2.
  n=128
  CALL fill_table(xmin,xmax,x,n)

  !This example program solves the equation log(2+x)=x
  !which has solution x=1.14619 (Wolfram Alpha)
  !and has no analytic solution
  !Fill L(x)=log(2+x) and R(x)=x
  ALLOCATE(LHS(n),RHS(n))
  DO i=1,n
     LHS(i)=log(2.+x(i))
     RHS(i)=x(i)
     !WRITE(*,*) x(i), LHS(i), RHS(i)
  END DO

  sol=bisect_solve(x,RHS-LHS,n,1e-6)

  WRITE(*,*) 'Solution:', sol

  sol=find_solve(0.,x,RHS-LHS,n)

  WRITE(*,*) 'Solution:', sol

END PROGRAM solve_equations_test
