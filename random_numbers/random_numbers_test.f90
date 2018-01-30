PROGRAM random_numbers_test

  USE random_numbers

  IMPLICIT NONE
  INTEGER :: i, nbins, j, n, iran, ri, iseed
  REAL, ALLOCATABLE :: binval(:), func(:), bin(:), x(:)
  REAL :: rn, binmin, binmax, mean, sig, pn
  REAL, PARAMETER :: pi=3.141592654

  WRITE(*,*)
  WRITE(*,*) 'Random number generation'
  WRITE(*,*)

  !This sets the 0 point for the random generator!
  !Not optional unless you want the same numbers on each use!
  WRITE(*,*) 'Seed for RNG set (0 sets it by the clock):'
  READ(*,*) iseed
  WRITE(*,*)
  CALL RNG_set(iseed)

  WRITE(*,*) '1 - Uniform distribution'
  WRITE(*,*) '2 - Rayleigh distribution'
  WRITE(*,*) '3 - Poisson distribution'
  WRITE(*,*) '4 - Gaussian distribution'
  WRITE(*,*) '5 - Polynomial distribution'
  WRITE(*,*) '6 - Lorentzian distribution'
  WRITE(*,*) '7 - Uniform integer distribution'
  READ(*,*) iran
  WRITE(*,*)

  WRITE(*,*) 'Number of random numbers'
  READ(*,*) n
  WRITE(*,*) 
  !n=100000

  nbins=50

  IF(iran==1) THEN
     binmin=0.
     binmax=1.
     mean=0.5
     sig=1.
  ELSE IF(iran==2) THEN
     binmin=0.
     binmax=10.
     sig=2.
  ELSE IF(iran==3) THEN
     binmin=0.
     binmax=10.
     mean=2.
  ELSE IF(iran==4) THEN
     binmin=-10.
     binmax=10.
     sig=2.
     mean=2.
  ELSE IF(iran==5) THEN
     binmin=0.
     binmax=1.
     WRITE(*,*) 'What order polynomial?'
     READ(*,*) pn
  ELSE IF(iran==6) THEN
     binmin=0.
     binmax=10.
  ELSE IF(iran==7) THEN
     binmin=0.
     binmax=10.
  ELSE
     STOP 'Error, demonstration not chosen correctly'
  END IF

  ALLOCATE(binval(nbins+1),bin(nbins),func(nbins),x(nbins))

  bin=0

  DO i=1,nbins+1
     binval(i)=binmin+(binmax-binmin)*float(i-1)/float(nbins)
  END DO

  DO i=1,nbins
     x(i)=(binval(i)+binval(i+1))/2.
  END DO

  IF(iran==1) func=1.
  IF(iran==2) func=x*exp(-(x**2.)/(2.*(sig**2.)))/(sig**2.)
  IF(iran==3) func=exp(-x/mean)/mean
  IF(iran==4) func=exp(-((x-mean)**2.)/(2.*(sig**2.)))/(sig*sqrt(2.*pi))
  IF(iran==5) func=(pn+1.)*(x**(pn))
  IF(iran==6) func=2./(pi*(1.+x**2.))
  IF(iran==7) func=1.

  WRITE(*,*) 'Creating random numbers:', n

  OPEN(8,file='random_numbers.dat')
  DO i=1,n

     IF(iran==1) rn=uniform(0.,1.)
     IF(iran==2) rn=Rayleigh(sig)
     IF(iran==3) rn=Poisson(mean)
     IF(iran==4) rn=Gaussian(mean,sig)
     IF(iran==5) rn=polynomial(pn)
     IF(iran==6) rn=Lorentzian()
     IF(iran==7) ri=random_integer(1,9)

     WRITE(*,fmt='(I10,F14.7)') i, rn

     IF(iran .NE. 7) THEN

        WRITE(8,*) i, rn

        DO j=1,nbins
           IF(rn>=binval(j) .AND. rn<binval(j+1)) THEN
              bin(j)=bin(j)+1.
              EXIT
           END IF
        END DO

     ELSE

        WRITE(8,*) i, ri

        DO j=1,nbins
           IF(float(ri)>=binval(j) .AND. float(ri)<binval(j+1)) THEN
              bin(j)=bin(j)+1.
              EXIT
           END IF
        END DO

     END IF

  END DO
  CLOSE(8)

  DO i=1,nbins
     bin(i)=bin(i)/(n*(binval(i+1)-binval(i)))
  END DO

  OPEN(7,file='rng.dat')
  DO i=1,nbins
     WRITE(7,*) x(i), bin(i), func(i)
  END DO
  CLOSE(7)

  WRITE(*,*) 'Done'
  WRITE(*,*)

 END PROGRAM random_numbers_test
