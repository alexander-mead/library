MODULE array_operations

CONTAINS

  FUNCTION sum_double(a,n)

    IMPLICIT NONE
    REAL :: sum_double
    REAL, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION :: sum
    INTEGER :: i
    
    sum=0.d0

    DO i=1,n
       sum=sum+a(i)
    END DO

    sum_double=REAL(sum)

  END FUNCTION sum_double

  SUBROUTINE amputate(arr,n_old,n_new)

    !Chop an array down to a smaller size
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
    REAL, ALLOCATABLE :: hold(:)
    INTEGER, INTENT(IN) :: n_new
    INTEGER, INTENT(IN) :: n_old
    INTEGER :: i    

    IF(n_old<n_new) STOP 'AMPUTATE: Error, new array should be smaller than the old one'

    ALLOCATE(hold(n_old))
    hold=arr
    DEALLOCATE(arr)
    ALLOCATE(arr(n_new))
    
    DO i=1,n_new
       arr(i)=hold(i)
    END DO
    
    DEALLOCATE(hold)

  END SUBROUTINE amputate

  SUBROUTINE reduce(arr1,n1,arr2,n2)

    !Reduces the size of array1 to the size of array2
    !This will not preserve the spacing of entries in array1
    !Thus it might be a terrible idea in many cases
    IMPLICIT NONE
    REAL, INTENT(IN) :: arr1(n1)
    REAL, INTENT(OUT) :: arr2(n2)
    INTEGER, INTENT(IN) :: n1, n2
    INTEGER :: i, j

    DO i=1,n2
       j=1+CEILING(float((n1-1)*(i-1))/float(n2-1))
       arr2(i)=arr1(j)
    END DO

  END SUBROUTINE reduce

  SUBROUTINE reduceto(arr1,n)

    !Reduces the array from whatever size to size 'n'
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: arr1(:)
    INTEGER, INTENT(IN) :: n
    REAL, ALLOCATABLE :: hold(:)
    INTEGER :: i, j

    ALLOCATE(hold(n))

    !n1=SIZE(arr1)

    DO i=1,n
       j=1+CEILING(REAL((n-1)*(i-1))/REAL(n-1))
       hold(i)=arr1(j)
    END DO

    DEALLOCATE(arr1)
    ALLOCATE(arr1(n))

    arr1=hold

    DEALLOCATE(hold)

  END SUBROUTINE reduceto

  SUBROUTINE reverse(arry,n)

    !This reverses the contents of arry!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: arry(n)
    INTEGER :: i
    REAL, ALLOCATABLE :: hold(:) 

    ALLOCATE(hold(n))

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

    DEALLOCATE(hold)

  END SUBROUTINE reverse

!!$  SUBROUTINE reverse(arry)
!!$
!!$    IMPLICIT NONE
!!$    INTEGER :: n, i
!!$    REAL, ALLOCATABLE :: hold(:)
!!$    REAL :: arry(:)
!!$
!!$    !This reverses the contents of arry!
!!$
!!$    n=SIZE(arry)
!!$
!!$    ALLOCATE(hold(n))
!!$
!!$    hold=arry
!!$
!!$    DO i=1,n
!!$       arry(i)=hold(n-i+1)
!!$    END DO
!!$
!!$    DEALLOCATE(hold)
!!$
!!$  END SUBROUTINE reverse

  FUNCTION splay(a,n1,n2,n3)

    IMPLICIT NONE
    REAL :: splay(n1*n2*n3)
    REAL, INTENT(IN) :: a(n1,n2,n3)
    INTEGER, INTENT(IN) :: n1, n2, n3
    INTEGER :: i, j, k, ii

    !This splays out a 3d array 'a' into a 1d array 'b' of the same size (n1*n2*n3)

    !n1=SIZE(a,1)
    !n2=SIZE(a,2)
    !n3=SIZE(a,3)

    !ALLOCATE(b(n1*n2*n3))

    !b=0.
    ii=0

    DO i=1,n1
       DO j=1,n2
          DO k=1,n3             
             ii=ii+1
             splay(ii)=a(i,j,k)
          END DO
       END DO
    END DO

  END FUNCTION splay

  SUBROUTINE binning(a,a1,a2,n,b,c,ilog)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, ilog
    REAL, INTENT(IN) :: a(:)
    REAL, ALLOCATABLE :: b(:), c(:)
    REAL :: a1, a2, min, max
    REAL, ALLOCATABLE :: binlim(:)
    INTEGER :: i, j

    min=a1
    max=a2

    WRITE(*,*) 'Binning'
    WRITE(*,*) 'Min:', min
    WRITE(*,*) 'Max:', max

    ALLOCATE(binlim(n+1),b(n),c(n))

    IF(ilog==1) THEN
       min=log10(min)
       max=log10(max)
    END IF

    !This sets the limits for the bins!
    DO i=1,n+1
       binlim(i)=min+(max-min)*float(i-1)/float(n)
    END DO

    !This sets the centre value for each bin!
    DO i=1,n
       b(i)=(binlim(i)+binlim(i+1))/2.
    END DO

    IF(ilog==1) THEN
       binlim=10.**binlim
       b=10.**b
    END IF

    c=0.

    DO i=1,SIZE(a)
       DO j=1,n
          IF(a(i)>=binlim(j) .AND. a(i)<=binlim(j+1)) THEN
             c(j)=c(j)+1.
          END IF
       END DO
    END DO

    DEALLOCATE(binlim)

    WRITE(*,*) 'Binning complete'
    WRITE(*,*)

  END SUBROUTINE binning

!!$  SUBROUTINE fill_array(array,min,max,n,ilog)
!!$
!!$    IMPLICIT NONE
!!$    REAL, ALLOCATABLE :: array(:)
!!$    REAL :: min, max
!!$    REAL :: a, b
!!$    INTEGER :: n, ilog
!!$    INTEGER :: i
!!$
!!$    IF(ilog==1) THEN
!!$       a=log10(min)
!!$       b=log10(max)
!!$    ELSE
!!$       a=min
!!$       b=max
!!$    END IF
!!$
!!$    ALLOCATE(array(n))
!!$
!!$    DO i=1,n
!!$       array(i)=a+(b-a)*float(i-1)/float(n-1)
!!$    END DO
!!$
!!$    IF(ilog==1) array=10.**array
!!$
!!$  END SUBROUTINE fill_array
!!$
!!$  SUBROUTINE fill_linear(a,b,arr)
!!$
!!$    IMPLICIT NONE
!!$    INTEGER :: i, n
!!$    REAL, INTENT(IN) :: a, b
!!$    REAL, INTENT(OUT) :: arr(:)
!!$
!!$    n=SIZE(arr)
!!$      
!!$    arr=0.
!!$
!!$    DO i=1,n
!!$       arr(i)=a+(b-a)*float(i-1)/float(n-1)
!!$    END DO
!!$
!!$  END SUBROUTINE fill_linear

  SUBROUTINE merge_arrays(a,na,b,nb,c,nc)

    !Takes arrays a and b and merges them together to make c
    !with length SIZE(a)+SIZE(b)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a(na), b(nb)
    REAL, ALLOCATABLE, INTENT(OUT) :: c(:)
    INTEGER, INTENT(IN) :: na, nb
    INTEGER, INTENT(OUT) :: nc
    INTEGER :: i
    
    nc=na+nb

    IF(ALLOCATED(c)) DEALLOCATE(c)
    ALLOCATE(c(nc))

    DO i=1,na
       c(i)=a(i)
    END DO

    DO i=1,nb
       c(i+na)=b(i)
    END DO

  END SUBROUTINE merge_arrays

  FUNCTION concatenate_arrays(a,na,b,nb)

    !Concatenate arrays a and b to form new array with length SIZE(a)+SIZE(b)
    IMPLICIT NONE
    REAL :: concatenate_arrays(na+nb)
    REAL, INTENT(IN) :: a(na), b(nb)
    INTEGER, INTENT(IN) :: na, nb
    INTEGER :: i

    DO i=1,na
       concatenate_arrays(i)=a(i)
    END DO

    DO i=1,nb
       concatenate_arrays(i+na)=b(i)
    END DO

  END FUNCTION concatenate_arrays

  SUBROUTINE fill_array(min,max,arr,n)

    !Fills array 'arr' in equally spaced intervals
    !I'm not sure if inputting an array like this is okay
    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
    INTEGER, INTENT(IN) :: n

    !Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          !arr(i)=min+(max-min)*REAL(i-1)/REAL(n-1)
          arr(i)=progression(min,max,i,n)
       END DO
    END IF

  END SUBROUTINE fill_array

  SUBROUTINE fill_array8(min,max,arr,n)

    !Fills array 'arr' in equally spaced intervals
    !I'm not sure if inputting an array like this is okay
    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: arr(:)
    INTEGER, INTENT(IN) :: n

    !Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.d0

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          !arr(i)=min+(max-min)*DBLE(i-1)/DBLE(n-1)
          arr(i)=progression8(min,max,i,n)
       END DO
    END IF

  END SUBROUTINE fill_array8

  FUNCTION progression(xmin,xmax,i,n)

    IMPLICIT NONE
    REAL :: progression
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression=xmin+(xmax-xmin)*REAL(i-1)/REAL(n-1)
    
  END FUNCTION progression

  FUNCTION progression_log(xmin,xmax,i,n)

    IMPLICIT NONE
    REAL :: progression_log
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression_log=exp(progression(log(xmin),log(xmax),i,n))
    
  END FUNCTION progression_log

  FUNCTION progression8(xmin,xmax,i,n)

    IMPLICIT NONE
    DOUBLE PRECISION :: progression8
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression8=xmin+(xmax-xmin)*DBLE(i-1)/DBLE(n-1)
    
  END FUNCTION progression8

  FUNCTION maximum(x,y,n)

    USE fix_polynomial
    !From an array y(x) finds the x location of the first maximum
    IMPLICIT NONE
    REAL :: maximum
    REAL, INTENT(IN) :: x(n), y(n)
    INTEGER, INTENT(IN) :: n
    REAL :: x1, x2, x3, y1, y2, y3, a, b, c
    INTEGER :: i

    !Need this to stop a compile-time warning
    maximum=0.

    DO i=1,n-1
       IF(y(i+1)<y(i)) THEN

          !Get the x positions
          x1=x(i-1)
          x2=x(i)
          x3=x(i+1)

          !Get the y values
          y1=y(i-1)
          y2=y(i)
          y3=y(i+1)

          !Fix a quadratic around the maximum
          CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          !Read off the maximum x from the parabola
          maximum=-b/(2.*a)

          !Exit the loop
          EXIT
       ELSE IF(i<n-1) THEN
          CYCLE
       ELSE
          STOP 'MAXIMUM: Error, array does not have a maximum'
       END IF
    END DO

  END FUNCTION maximum

END MODULE array_operations
