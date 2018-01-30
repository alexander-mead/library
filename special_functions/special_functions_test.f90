PROGRAM special_functions_test

  USE special_functions

  IMPLICIT NONE
  REAL :: x, f
  REAL :: xmin, xmax
  INTEGER :: i, n, m
  REAL :: x1, x2, pow
  INTEGER :: test

  WRITE(*,*) 'Choose a special function to test'
  WRITE(*,*) '1 - Apodise'
  WRITE(*,*) '2 - Blob'
  WRITE(*,*) '3 - Smooth apodise'
  WRITE(*,*) '4 - Smooth blob'
  WRITE(*,*) '5 - Sinc'
  WRITE(*,*) '6 - Top hat Fourier transform'
  WRITE(*,*) '7 - J0'
  WRITE(*,*) '8 - J1'
  WRITE(*,*) '9 - J2'
  READ(*,*) test
  WRITE(*,*)
  
  xmin=0.
  xmax=10.
  n=500

  IF(test==1 .OR. test==3) THEN
     x1=4.
     x2=6.
     pow=0.5
  ELSE IF(test==2 .OR. test==4) THEN
     x1=2.
     x2=8.
     pow=0.1
  END IF

  OPEN(7,file='results.dat')
  DO i=1,n
     x=xmin+(xmax-xmin)*float(i-1)/float(n-1)
     IF(test==1) THEN
        f=apodise(x,x1,x2,pow)
     ELSE IF(test==2) THEN
        f=blob(x,x1,x2,pow)
     ELSE IF(test==3) THEN
        f=smoothapodise(x,x1,x2,pow)
     ELSE IF(test==4) THEN
        f=smoothblob(x,x1,x2,pow)
     ELSE IF(test==5) THEN
        f=sinc(x)
     ELSE IF(test==6) THEN
        f=wk_tophat(x)
     ELSE IF(test==7 .OR. test==8 .OR. test==9) THEN
        IF(test==7) m=0
        IF(test==8) m=1
        IF(test==9) m=2
        f=Bessel(m,x)
     ELSE
        STOP 'SPECIAL_FUNCTION_TEST: Error, function not specified correctly'
     END IF
     WRITE(7,*) x, f
  END DO
  CLOSE(7)
  
END PROGRAM special_functions_test
