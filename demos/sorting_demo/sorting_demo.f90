PROGRAM sorting_test

  USE sorting
  USE random_numbers

  IMPLICIT NONE
  REAL, ALLOCATABLE :: a(:), b(:)
  INTEGER, ALLOCATABLE :: idx(:)
  INTEGER :: i, isort, itest, imode, n
  REAL :: t1, t2
  LOGICAL :: sorted

  REAL, PARAMETER :: rmin=0.    ! Minimum random number
  REAL, PARAMETER :: rmax=100.  ! Maximum random number
  INTEGER, PARAMETER :: iseed=0 ! Random number seed
  INTEGER, PARAMETER :: imeth=1 ! Sorting method for timing tests

  WRITE(*,*)

  CALL RNG_set(iseed)

  WRITE(*,*) 'SORTING_TESTS: Choose test'
  WRITE(*,*) '1 - Veracity tests'
  WRITE(*,*) '2 - Timing tests'
  READ(*,*) imode
  WRITE(*,*)
  
  IF(imode==1) THEN

     ! Small tests
     n=16

     ALLOCATE(a(n),b(n),idx(n))

     ! Initial white space
     WRITE(*,*)

     ! Array elements
     DO i=1,n-1
        a(i)=CEILING(random_uniform(rmin,rmax))
     END DO
     a(n)=a(1) ! Having a repeated element is good for tests

     ! Write original array to screen
     WRITE(*,*) 'SORTING_TEST: Original array'
     WRITE(*,*) '=============='
     DO i=1,16
        WRITE(*,fmt='(I5,F10.5)') i, a(i)
     END DO
     WRITE(*,*) '=============='
     WRITE(*,*)

     DO itest=1,2

        IF(itest==1) THEN

           DO isort=1,3

              ! Needed because these methods overwrite the input array
              b=a

              CALL sort(b,n,isort)

              ! Direct array re-arranging
              WRITE(*,*) 'SORTING_TEST: Sorting method:', isort
              WRITE(*,*) 'SORTING_TEST: Sorted array'
              WRITE(*,*) '=============='
              DO i=1,SIZE(b)
                 WRITE(*,fmt='(I5,F10.5)') i, b(i)
              END DO
              WRITE(*,*) '=============='

              sorted=check_sorted(b,n)
              WRITE(*,*) 'SORTING_TEST: Is array sorted?:', sorted
              IF(.NOT. sorted) STOP 'SORTING_TEST: Error, something has fucked up'
              WRITE(*,*)

           END DO

           DEALLOCATE(b)

        ELSE IF(itest==2) THEN

           ! Indexing methods
           DO isort=1,2

              CALL index(a,idx,16,isort)

              WRITE(*,*) 'SORTING_TEST: Index method:', isort
              WRITE(*,*) 'SORTING_TEST: Sorted array'
              WRITE(*,*) '=============='
              DO i=1,n
                 WRITE(*,fmt='(I5,F10.5)') i, a(idx(i))
              END DO
              WRITE(*,*) '=============='

              sorted=check_sorted_index(a,idx,n)
              WRITE(*,*) 'SORTING_TEST: Is array sorted?:', sorted
              IF(.NOT. sorted) STOP 'SORTING_TEST: Error, something has fucked up'
              WRITE(*,*)

           END DO

        ELSE
           STOP 'SORTING_TEST: Error, test not specified correctly'
        END IF

     END DO

     WRITE(*,*) 'SORTING_TEST: Tests complete'
     WRITE(*,*)

  ELSE IF(imode==2) THEN

     ! Sorting timing tests

     WRITE(*,*) 'SORTING_TEST: Method', imeth
     WRITE(*,*) '================================='
     WRITE(*,*) '        n         T [s]    sorted'
     WRITE(*,*) '================================='
     DO itest=1,5
        
        n=10**itest
        ALLOCATE(a(n))

        DO i=1,n
           a(i)=random_uniform(rmin,rmax)
        END DO

        CALL cpu_time(t1)
        CALL sort(a,n,imeth)
        CALL cpu_time(t2)

        WRITE(*,fmt='(I10,F14.7,L10)') n, t2-t1, check_sorted(a,n)       

        DEALLOCATE(a)
        
     END DO
     WRITE(*,*) '================================='

  ELSE

  END IF

END PROGRAM sorting_test
