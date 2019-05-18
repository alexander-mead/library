MODULE camb_stuff
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_CAMB_Pk

CONTAINS

  SUBROUTINE read_CAMB_Pk(k,p,n,infile)

    USE file_info
    USE constants
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:), p(:)
    INTEGER, INTENT(OUT) :: n
    INTEGER :: i
    
    n=file_length(infile,verbose=.FALSE.)
    n=n-1
    WRITE(*,*) 'READ_CAMB_PK: CAMB file: ', trim(infile)
    WRITE(*,*) 'READ_CAMB_PK: Number of points:', n
    WRITE(*,*)

    ALLOCATE(k(n),p(n))

    OPEN(7,file=infile)
    DO i=0,n
       IF(i==0) THEN
          READ(7,*)
       ELSE
          READ(7,*) k(i), p(i)
       END IF
    END DO
    CLOSE(7)

    p=4.*pi*p*(k**3)/twopi**3

  END SUBROUTINE read_CAMB_Pk
  
END MODULE camb_stuff
