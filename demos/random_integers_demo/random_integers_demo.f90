PROGRAM random_integers_test

  USE random_numbers
  IMPLICIT NONE
  INTEGER :: i
  INTEGER :: r
  INTEGER :: itest
  INTEGER :: n
  INTEGER, PARAMETER :: iseed=0
  INTEGER, PARAMETER :: imin=1
  INTEGER, PARAMETER :: imax=10

  WRITE(*,*)
  
  WRITE(*,*) 'Choose test'
  WRITE(*,*) '1 - Write out 50 random integers between 1 and 10'
  WRITE(*,*) '2 - Test random sign'
  WRITE(*,*) '3 - Check lots of random integers to make sure they are inside bounds'
  READ(*,*) itest
  WRITE(*,*)

  IF(itest==1 .OR. itest==2) THEN
     n=50
  ELSE IF(itest==3) THEN
     n=10000000
  ELSE
     STOP 'Error, something went wrong setting n'
  END IF

  CALL RNG_set(iseed)
  
  DO i=1,n
     
     IF(itest==1 .OR. itest==3) THEN
        r=random_integer(imin,imax)
     ELSE IF(itest==2) THEN
        r=random_sign()
     END IF

     IF(itest==1 .OR. itest==2) THEN
        WRITE(*,*) i, r
     ELSE IF(itest==3 .AND. (r<imin .OR. r>imax)) THEN
        WRITE(*,*) i, r
        WRITE(*,*) 'A giant fucking disaster has occurred'
     END IF
     
  END DO
  
END PROGRAM random_integers_test
