PROGRAM correlation_test

  USE random_numbers
  USE simulations

  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:,:), w(:)
  REAL, ALLOCATABLE :: r(:), xi(:)
  INTEGER, ALLOCATABLE :: nbin(:)
  INTEGER :: i
  
  INTEGER, PARAMETER :: n=10000
  REAL, PARAMETER :: L=500.
  REAL, PARAMETER :: rmin=L/real(n)**(1./3.)
  REAL, PARAMETER :: rmax=L/2.
  INTEGER, PARAMETER :: nr=100
  INTEGER, PARAMETER :: iseed=0
  REAL, PARAMETER :: nbar=n/L**3

  CALL RNG_set(iseed)

  ALLOCATE(x(3,n),w(n))
  CALL generate_randoms(x,n,L)
  w=1.

  ALLOCATE(r(nr),xi(nr),nbin(nr))
  CALL correlation_function(rmin,rmax,r,xi,nbin,nr,x,x,w,w,n,n,L)

  OPEN(7,file='xi.dat')
  DO i=1,nr
     IF(nbin(i) .NE. 0) THEN
        WRITE(7,*) r(i), xi(i), nbin(i)
     END IF
  END DO
  CLOSE(7)
  
END PROGRAM correlation_test
