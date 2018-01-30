PROGRAM table_integer_test

  USE table_integer
  
  IMPLICIT NONE
  REAL :: xmin, xmax, x
  REAL, ALLOCATABLE :: xtab(:)
  INTEGER :: i, n, m, int1, int2, int3

  n=6

  ALLOCATE(xtab(n))

  xmin=0.
  xmax=1.

  WRITE(*,*)
  WRITE(*,*) 'Original table'
  WRITE(*,*) '=============='
  DO i=1,n
     xtab(i)=xmin+(xmax-xmin)*float(i-1)/float(n-1)
     WRITE(*,*) i, xtab(i)
  END DO
  WRITE(*,*)

  m=10*n
  
  WRITE(*,*) 'Table positions'
  WRITE(*,*) '==============='
  DO i=1,m
     x=xmin+(xmax-xmin)*float(i-1)/float(m-1)
     int1=select_table_integer(x,xtab,n,1)
     int2=select_table_integer(x,xtab,n,2)
     int3=select_table_integer(x,xtab,n,3)
     WRITE(*,*) x, int1, int2, int3
  END DO
  WRITE(*,*)

END PROGRAM table_integer_test
