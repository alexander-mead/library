MODULE sorting

CONTAINS

  SUBROUTINE bubble(a,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: a(n) 
    REAL :: hold
    INTEGER :: i, isort

    !This is a dodgy module!
    !I think it overwrites some of the things
    !in the original array somehow!

    DO
       isort=0
       DO i=1,n-1
          IF(a(i)>a(i+1)) THEN
             hold=a(i+1)
             a(i+1)=a(i)
             a(i)=hold
             isort=1
          END IF
       END DO
       IF(isort==0) EXIT
    END DO

  END SUBROUTINE bubble

  SUBROUTINE bubble_index(n,a,ind)

    IMPLICIT NONE
    INTEGER :: n, i, isort, hold
    REAL, INTENT(IN) :: a(:)
    INTEGER, INTENT(OUT) :: ind(:)

    WRITE(*,*) 'Starting bubble index'

    DO i=1,n
       ind(i)=i
    END DO

    DO
       isort=0
       DO i=1,n-1
          IF(a(ind(i))>a(ind(i+1))) THEN
             hold=ind(i+1)
             ind(i+1)=ind(i)
             ind(i)=hold
             isort=1
          END IF
       END DO
       IF(isort==0) EXIT
    END DO

    WRITE(*,*) 'Bubble index finished'
    WRITE(*,*)

  END SUBROUTINE bubble_index

  SUBROUTINE selection_sort(a,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: a(n)
    REAL :: hold, min
    INTEGER :: i, j, minl

    DO i=1,n-1
       min=a(i)
       minl=i
       DO j=i+1,n
          IF(a(j)<min) THEN
             min=a(j)
             minl=j
          END IF
       END DO
       hold=a(i)
       a(i)=min
       a(minl)=hold
    END DO

  END SUBROUTINE selection_sort
  
END MODULE sorting
