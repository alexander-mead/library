PROGRAM solve_equations_test

  USE array_operations
  USE solve_equations

  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:), LHS(:), RHS(:)
  REAL :: sol
  INTEGER :: i

  REAL, PARAMETER :: xmin=0.
  REAL, PARAMETER :: xmax=2.
  INTEGER, PARAMETER :: n=128

  REAL, PARAMETER :: alpha=1.14619322062058

  WRITE(*,*)
  WRITE(*,*) 'Testing routines to solve equations'
  WRITE(*,*)

  WRITE(*,*) 'Solving equation x = ln(2+x)'
  WRITE(*,*) 'Wolfram alpha solution: x =', alpha
  WRITE(*,*)
  
  !Fill array of x value
  CALL fill_array(xmin,xmax,x,n)

  !This example program solves the equation log(2+x) = x
  !This equation has solution x=1.14619 (Wolfram Alpha) but no analytic solution
  !Fill L(x)=log(2+x) and R(x)=x
  ALLOCATE(LHS(n),RHS(n))
  DO i=1,n
     LHS(i)=log(2.+x(i))
     RHS(i)=x(i)
  END DO

  WRITE(*,*) 'Solving using bisect solve'
  sol=bisect_solve(x,RHS-LHS,n,1e-6)
  WRITE(*,*) 'Solution:', sol
  WRITE(*,*) 'Accuracy:', -1.+sol/alpha
  WRITE(*,*)

  WRITE(*,*) 'Solving using find solve'
  sol=find_solve(0.,x,RHS-LHS,n)
  WRITE(*,*) 'Solution:', sol
  WRITE(*,*) 'Accuracy:', -1.+sol/alpha
  WRITE(*,*)

END PROGRAM solve_equations_test
