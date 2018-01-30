PROGRAM sorting_test

  USE sorting
  
  IMPLICIT NONE
  REAL :: a(16), b(16)
  INTEGER :: i, isort, idx(16)

  a(1)=1
  a(2)=6
  a(3)=1
  a(4)=2
  a(5)=10
  a(6)=8
  a(7)=1
  a(8)=4
  a(9)=1
  a(10)=5
  a(11)=3
  a(12)=13
  a(13)=5
  a(14)=16
  a(15)=15
  a(16)=5

  b=a

  WRITE(*,*) '1 - Bubble sort'
  WRITE(*,*) '2 - Selection sort'
  WRITE(*,*) '3 - Bubble index'
  !WRITE(*,*) '11 - NR Indexing'
  READ(*,*) isort
  WRITE(*,*)

  IF(isort==1)  CALL bubble(b,16)
  IF(isort==2)  CALL selection_sort(b,16)
  !IF(isort==11) CALL indexx(16,a,idx)
  IF(isort==3) CALL bubble_index(16,a,idx)

  WRITE(*,*) 'List'
  WRITE(*,*) '****'

  IF(isort==1 .OR. isort==2) THEN
     !Not index methods
     DO i=1,SIZE(b)
        WRITE(*,*) i, a(i), b(i)
     END DO
  ELSE IF(isort==3) THEN
     !Index methods
     DO i=1,SIZE(a)
        WRITE(*,*) i, a(i), a(idx(i))
     END DO     
  ELSE
     STOP 'SORTING_TEST: Error, sorting method not specified correctly'
  END IF

END PROGRAM sorting_test
