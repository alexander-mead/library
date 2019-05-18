MODULE array_operations

  ! These interfaces did not work if -fdefault-real-8 and -fdefault-real-8 were set
!!$  INTERFACE fill_array
!!$     MODULE PROCEDURE fill_array_single
!!$     MODULE PROCEDURE fill_array_double
!!$  END INTERFACE fill_array
!!$
!!$  INTERFACE progression
!!$     MODULE PROCEDURE progression_single
!!$     MODULE PROCEDURE progression_double
!!$  END INTERFACE progression
!!$
!!$  INTERFACE progression_log
!!$     MODULE PROCEDURE progression_log_single
!!$     MODULE PROCEDURE progression_log_double
!!$  END INTERFACE progression_log

  PRIVATE

  PUBLIC :: within_array
  PUBLIC :: swap_arrays
  PUBLIC :: add_to_array  
  PUBLIC :: splay
  PUBLIC :: write_array_list
  PUBLIC :: maximum
  PUBLIC :: sum_double
  PUBLIC :: repeated_entries
  PUBLIC :: array_position  
  PUBLIC :: reverse_array
  PUBLIC :: remove_array_element
  PUBLIC :: array_positions
  PUBLIC :: amputate_array
  PUBLIC :: unique_index  
  PUBLIC :: reduce_array
  PUBLIC :: fill_pixels
  PUBLIC :: binning
  PUBLIC :: merge_arrays
  PUBLIC :: mask

  PUBLIC :: unique_entries
  PUBLIC :: number_of_appearances
  PUBLIC :: remove_repeated_array_elements
  PUBLIC :: remove_repeated_two_array_elements

  PUBLIC :: progression
  PUBLIC :: progression_double
  PUBLIC :: progression_log
  PUBLIC :: progression_log_double
  
  PUBLIC :: fill_array
  PUBLIC :: fill_array_double

  INTERFACE add_to_array
     MODULE PROCEDURE add_to_array_2D
     MODULE PROCEDURE add_to_array_3D
  END INTERFACE add_to_array

  INTERFACE splay
     PROCEDURE splay_2D
     PROCEDURE splay_3D
  END INTERFACE splay

  INTERFACE write_array_list
     PROCEDURE write_array_list_real
     PROCEDURE write_array_list_int
  END INTERFACE write_array_list

CONTAINS

  LOGICAL FUNCTION within_array(x,a,n)

    ! Checks to see if x falls within the range of values in array a
    IMPLICIT NONE
    REAL, INTENT(IN) :: x    ! Value to check
    REAL, INTENT(IN) :: a(n) ! Array (of x values, presumably)
    INTEGER, INTENT(IN) :: n ! Size of array

    IF(x>=minval(a) .AND. x<=maxval(a)) THEN
       within_array=.TRUE.
    ELSE
       within_array=.FALSE.
    END IF
    
  END FUNCTION within_array

  SUBROUTINE swap_arrays(x,y,n)

    ! Swap arrays x<->y in a memory-efficient way
    ! Only one excess real number is ever stored
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(n) ! Array 1
    REAL, INTENT(INOUT) :: y(n) ! Array 2
    INTEGER, INTENT(IN) :: n    ! Size of arrays
    INTEGER :: i
    REAL :: h

    ! Loop over array and swap element by element
    DO i=1,n
       h=x(i)
       x(i)=y(i)
       y(i)=h
    END DO
    
  END SUBROUTINE swap_arrays

  SUBROUTINE add_to_array_2D(a,m,v,i)

    ! Add value 'v' to array 'a' at location 'i' in array
    ! If 'i' is outside the array range then this routine does nothing
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a(m,m)
    INTEGER, INTENT(IN) :: m
    REAL, INTENT(IN) :: v
    INTEGER, INTENT(IN) :: i(2)
    INTEGER :: j
    LOGICAL :: bin
    INTEGER, PARAMETER :: dim=2

    bin=.TRUE.
    DO j=1,dim
       IF(i(j)<1 .OR. i(j)>m) THEN
          bin=.FALSE.
          EXIT
       END IF
    END DO

    IF(bin) a(i(1),i(2))=a(i(1),i(2))+v
    
  END SUBROUTINE add_to_array_2D

   SUBROUTINE add_to_array_3D(a,m,v,i)

    ! Add value 'v' to array 'a' at location 'i' in array
    ! If 'i' is outside the array range then this routine does nothing
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a(m,m,m)
    INTEGER, INTENT(IN) :: m
    REAL, INTENT(IN) :: v
    INTEGER, INTENT(IN) :: i(3)
    INTEGER :: j
    LOGICAL :: bin
    INTEGER, PARAMETER :: dim=3

    bin=.TRUE.
    DO j=1,dim
       IF(i(j)<1 .OR. i(j)>m) THEN
          bin=.FALSE.
          EXIT
       END IF
    END DO

    IF(bin) a(i(1),i(2),i(3))=a(i(1),i(2),i(3))+v
    
  END SUBROUTINE add_to_array_3D

  INTEGER FUNCTION array_position(x,a,n)

    ! Returns the location in a(n) of value x
    ! If x is not in array then returns zero
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: x    ! Value to check if it is in array
    INTEGER, INTENT(IN) :: a(n) ! Array to check 
    INTEGER, INTENT(IN) :: n    ! Size of array
    INTEGER :: i

    array_position=0
    
    DO i=1,n
       IF(a(i)==x) THEN
          array_position=i
          EXIT
       END IF
    END DO
    
  END FUNCTION array_position

  INTEGER FUNCTION number_of_appearances(x,a,n)

    ! Returns the number of appearances in a(n) of value x
    ! If x is not in a(n) then returns zero
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: x    ! Value to check if it is in array
    INTEGER, INTENT(IN) :: a(n) ! Array to check 
    INTEGER, INTENT(IN) :: n    ! Size of array
    INTEGER :: i

    number_of_appearances=0
    
    DO i=1,n
       IF(a(i)==x) THEN
          number_of_appearances=number_of_appearances+1
       END IF
    END DO
    
  END FUNCTION number_of_appearances

  SUBROUTINE array_positions(x,a,n,b,m)

    ! Returns the location in the array of value x
    ! If x is not in array then returns zero
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: x    ! Value to check if it is in array
    INTEGER, INTENT(IN) :: a(n) ! Array to check 
    INTEGER, INTENT(IN) :: n    ! Size of array
    INTEGER, ALLOCATABLE, INTENT(OUT) :: b(:)
    INTEGER, INTENT(OUT) :: m
    INTEGER :: i, p

    m=number_of_appearances(x,a,n)
    ALLOCATE(b(m))

    p=0
    DO i=1,n
       IF(a(i)==x) THEN
          p=p+1
          b(p)=i
       END IF
    END DO

  END SUBROUTINE array_positions
  
  FUNCTION sum_double(a,n)

    ! Sum using double precision, which is necessary for many array elements
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

    sum_double=real(sum)

  END FUNCTION sum_double

!!$  SUBROUTINE amputate(arr,n_old,n_new)
!!$
!!$    ! Chop an array down to a smaller size
!!$    ! TODO: Retire
!!$    IMPLICIT NONE
!!$    REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
!!$    REAL, ALLOCATABLE :: hold(:)
!!$    INTEGER, INTENT(IN) :: n_new
!!$    INTEGER, INTENT(IN) :: n_old
!!$    INTEGER :: i    
!!$
!!$    IF(n_old<n_new) STOP 'AMPUTATE: Error, new array should be smaller than the old one'
!!$
!!$    ALLOCATE(hold(n_old))
!!$    hold=arr
!!$    DEALLOCATE(arr)
!!$    ALLOCATE(arr(n_new))
!!$    
!!$    DO i=1,n_new
!!$       arr(i)=hold(i)
!!$    END DO
!!$    
!!$    DEALLOCATE(hold)
!!$
!!$  END SUBROUTINE amputate

  SUBROUTINE amputate_array(a,n,i1,i2)

    ! Chop an array of size a(n) down to a smaller size demarked by indices i1, i2
    ! If i1=1 and i2=n then this does nothing
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: i1, i2
    REAL, ALLOCATABLE :: b(:)
    INTEGER :: i, m

    IF(i2<i1) THEN
       STOP 'AMPUTATE: Error, i2 should be greater than i1'
    END IF

    m=i2-i1+1
    IF(n<m) THEN
       STOP 'AMPUTATE: Error, new array should be smaller than the old one'
    END IF

    ! Store input array and then deallocate
    ALLOCATE(b(n))
    b=a
    DEALLOCATE(a)

    ! Allocate new output array
    ALLOCATE(a(m))

    ! Fill new output array
    DO i=1,m
       a(i)=b(i+i1-1)
    END DO

    ! Deallocate holding array
    DEALLOCATE(b)

  END SUBROUTINE amputate_array

  SUBROUTINE reduce_array(arr1,n1,arr2,n2)

    ! Reduces the size of array1 to the size of array2
    ! This will not preserve the spacing of entries in array1 and might be a terrible idea in many cases
    IMPLICIT NONE
    REAL, INTENT(IN) :: arr1(n1)
    REAL, INTENT(OUT) :: arr2(n2)
    INTEGER, INTENT(IN) :: n1, n2
    INTEGER :: i, j

    DO i=1,n2
       j=1+ceiling(real((n1-1)*(i-1))/real(n2-1))
       arr2(i)=arr1(j)
    END DO

  END SUBROUTINE reduce_array

!!$  SUBROUTINE reduceto(arr1,n)
!!$
!!$    ! Reduces the array from whatever size to size 'n'
!!$    ! TODO: Remove
!!$    IMPLICIT NONE
!!$    REAL, ALLOCATABLE, INTENT(INOUT) :: arr1(:)
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, ALLOCATABLE :: hold(:)
!!$    INTEGER :: i, j
!!$
!!$    ALLOCATE(hold(n))
!!$
!!$    DO i=1,n
!!$       j=1+ceiling(real((n-1)*(i-1))/real(n-1))
!!$       hold(i)=arr1(j)
!!$    END DO
!!$
!!$    DEALLOCATE(arr1)
!!$    ALLOCATE(arr1(n))
!!$
!!$    arr1=hold
!!$
!!$    DEALLOCATE(hold)
!!$
!!$  END SUBROUTINE reduceto

  SUBROUTINE reverse_array(arry,n)

    ! Reverses the contents of arry
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: arry(n)
    INTEGER :: i
    REAL :: hold(n) 

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

  END SUBROUTINE reverse_array

  SUBROUTINE remove_array_element(a,n,i)

    ! Remove element i from array a(n)
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Input array
    INTEGER, INTENT(IN) :: n                 ! Original length of array, it will change to n-1
    INTEGER, INTENT(IN) :: i                 ! Element to remove
    REAL :: b(n-1)
    INTEGER :: j, jj

    IF(i<1 .OR. i>n) THEN
       WRITE(*,*) 'Array size:', n
       WRITE(*,*) 'Element to remove:', i
       STOP 'REMOVE_ARRAY_ELEMENT: Error, element to remove is outside array bounds'
    END IF

    jj=0
    DO j=1,n
       IF(j==i) THEN
          CYCLE
       ELSE
          jj=jj+1
          b(jj)=a(j)
       END IF
    END DO

    DEALLOCATE(a)
    ALLOCATE(a(n-1))
    a=b
    
  END SUBROUTINE remove_array_element

  SUBROUTINE remove_repeated_array_elements(a,n,m)

    ! Remove any repeated entries from the array
    ! Assumes the array is sorted
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: a(:) ! Array to consider
    INTEGER, INTENT(IN) :: n                 ! Orignal size of array
    INTEGER, INTENT(OUT) :: m                ! Final size of array
    INTEGER :: i
    LOGICAL :: remove(n)

    remove=.FALSE.

    m=n
    DO i=1,n-1
       IF(a(i)==a(i+1)) THEN
          remove(i+1)=.TRUE.
          m=m-1
       END IF
    END DO

    a=PACK(a,.NOT. remove)
    
  END SUBROUTINE remove_repeated_array_elements

  SUBROUTINE remove_repeated_two_array_elements(a,b,n,m)

    ! Remove any repeated entries in the array a from both a and b
    ! Assumes the array a is sorted
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
    REAL, ALLOCATABLE, INTENT(INOUT) :: b(:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: m
    INTEGER :: i
    LOGICAL :: remove(n)

    remove=.FALSE.

    m=n
    DO i=1,n-1
       IF(a(i)==a(i+1)) THEN
          remove(i+1)=.TRUE.
          m=m-1
       END IF
    END DO

    a=PACK(a,.NOT. remove)
    b=PACK(b,.NOT. remove)
    
  END SUBROUTINE remove_repeated_two_array_elements

  SUBROUTINE write_array_list_real(a,n)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    WRITE(*,*) 'WRITE_ARRAY_LIST: Writing array'
    DO i=1,n
       WRITE(*,*) 'WRITE_ARRAY_LIST:', i, a(i)
    END DO
    WRITE(*,*) 'WRITE_ARRAY_LIST: Done'
    WRITE(*,*)

  END SUBROUTINE write_array_list_real

  SUBROUTINE write_array_list_int(a,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    WRITE(*,*) 'WRITE_ARRAY_LIST: Writing array'
    DO i=1,n
       WRITE(*,*) 'WRITE_ARRAY_LIST:', i, a(i)
    END DO
    WRITE(*,*) 'WRITE_ARRAY_LIST: Done'
    WRITE(*,*)

  END SUBROUTINE write_array_list_int

  FUNCTION splay_2D(a,n1,n2)

    ! This splays out a 3d array 'a' into a 1d array 'b' of the same size (n1*n2*n3)
    IMPLICIT NONE
    REAL :: splay_2D(n1*n2)
    REAL, INTENT(IN) :: a(n1,n2)
    INTEGER, INTENT(IN) :: n1, n2
    INTEGER :: i, j, ii

    ! Set sum integer to zero
    ii=0

    DO j=1,n2
       DO i=1,n1
          ii=ii+1
          splay_2D(ii)=a(i,j)
       END DO
    END DO

  END FUNCTION splay_2D

  FUNCTION splay_3D(a,n1,n2,n3)

    ! This splays out a 3d array 'a' into a 1d array 'b' of the same size (n1*n2*n3)
    ! TODO: Should i, j, k order of loops be reversed?
    IMPLICIT NONE
    REAL :: splay_3D(n1*n2*n3)
    REAL, INTENT(IN) :: a(n1,n2,n3)
    INTEGER, INTENT(IN) :: n1, n2, n3
    INTEGER :: i, j, k, ii

    ! Set sum integer to zero
    ii=0

    DO i=1,n1
       DO j=1,n2
          DO k=1,n3             
             ii=ii+1
             splay_3D(ii)=a(i,j,k)
          END DO
       END DO
    END DO

  END FUNCTION splay_3D

  SUBROUTINE binning(a,a1,a2,n,b,c,ilog)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, ilog
    REAL, INTENT(IN) :: a(:)
    REAL, INTENT(OUT) :: b(n), c(n)
    REAL :: a1, a2, min, max
    REAL, ALLOCATABLE :: binlim(:)
    INTEGER :: i, j

    min=a1
    max=a2

    WRITE(*,*) 'Binning'
    WRITE(*,*) 'Min:', min
    WRITE(*,*) 'Max:', max

    IF(ilog==1) THEN
       min=log10(min)
       max=log10(max)
    END IF

    ! This sets the limits for the bins!
    CALL fill_array(min,max,binlim,n+1)

    ! This sets the centre value for each bin!
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

    WRITE(*,*) 'Binning complete'
    WRITE(*,*)

  END SUBROUTINE binning

  SUBROUTINE merge_arrays(a,na,b,nb,c,nc)

    ! Takes arrays a and b and merges them together to make c with length SIZE(a)+SIZE(b)
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

    ! Concatenate arrays a and b to form new array with length SIZE(a)+SIZE(b)
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

    ! Fills array 'arr' in linearly spaced intervals
    ! e.g., 4 values between 0 and 1 would be 0, 1/3, 2/3, and 1
    ! This means that min and max are included in the array
    IMPLICIT NONE
    REAL, INTENT(IN) :: min ! Minimum value for array
    REAL, INTENT(IN) :: max ! Maximum value for array
    REAL, ALLOCATABLE, INTENT(OUT) :: arr(:) ! Output array
    INTEGER, INTENT(IN) :: n ! Number of entries in array
    INTEGER :: i

    ! Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=progression(min,max,i,n)
       END DO
    END IF

  END SUBROUTINE fill_array

  SUBROUTINE fill_pixels(min,max,arr,n)

    ! Fill an array between min and max as if the values should correspond to central pixel values
    ! e.g., 3 pixel values between 0 and 1 would be 1/6, 1/2 and 5/6
    ! This means that min and max are not included in the array
    IMPLICIT NONE
    REAL, INTENT(IN) :: min ! Minimum value bordering array
    REAL, INTENT(IN) :: max ! Maximum value bordering array
    REAL, ALLOCATABLE, INTENT(OUT) :: arr(:) ! Output array
    INTEGER, INTENT(IN) :: n ! Number of pixels in array
    REAL :: mmin, mmax, pixel_size

    pixel_size=(max-min)/real(n) ! Calculate the pixel size
    mmin=min+pixel_size/2. ! Offset minimum value by half a pixel into the region
    mmax=max-pixel_size/2. ! Offset maximum value by half a pixel into the region

    CALL fill_array(mmin,mmax,arr,n) ! Linearly fill array
    
  END SUBROUTINE fill_pixels

  SUBROUTINE fill_array_double(min,max,arr,n)

    ! Fills array 'arr' in equally spaced intervals
    ! TODO: Not sure if inputting an array like this is okay
    IMPLICIT NONE   
    DOUBLE PRECISION, INTENT(IN) :: min, max
    DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: arr(:)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    ! Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.d0

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=progression_double(min,max,i,n)
       END DO
    END IF

  END SUBROUTINE fill_array_double

  REAL FUNCTION progression(xmin,xmax,i,n)

    IMPLICIT NONE
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    IF(n==1) THEN
       progression=xmin
    ELSE
       progression=xmin+(xmax-xmin)*real(i-1)/real(n-1)
    END IF
    
  END FUNCTION progression

  DOUBLE PRECISION FUNCTION progression_double(xmin,xmax,i,n)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression_double=xmin+(xmax-xmin)*DBLE(i-1)/DBLE(n-1)
    
  END FUNCTION progression_double

  REAL FUNCTION progression_log(xmin,xmax,i,n)

    IMPLICIT NONE
    REAL, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression_log=exp(progression(log(xmin),log(xmax),i,n))
    
  END FUNCTION progression_log

  DOUBLE PRECISION FUNCTION progression_log_double(xmin,xmax,i,n)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: i, n

    progression_log_double=exp(progression_double(log(xmin),log(xmax),i,n))
    
  END FUNCTION progression_log_double

  REAL FUNCTION maximum(x,y,n)

    USE fix_polynomial
    
    ! From an array y(x) finds the x location of the first maximum
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(n), y(n)
    INTEGER, INTENT(IN) :: n
    REAL :: x1, x2, x3, y1, y2, y3, a, b, c
    INTEGER :: i

    ! Need this to stop a compile-time warning
    maximum=0.

    DO i=1,n-1
       
       IF(y(i+1)<y(i)) THEN

          ! Get the x positions
          x1=x(i-1)
          x2=x(i)
          x3=x(i+1)

          ! Get the y values
          y1=y(i-1)
          y2=y(i)
          y3=y(i+1)

          ! Fix a quadratic around the maximum
          CALL fix_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          ! Read off the maximum x from the parabola
          maximum=-b/(2.*a)

          ! Exit the loop
          EXIT
          
       ELSE IF(i<n-1) THEN
          
          CYCLE
          
       ELSE
          
          STOP 'MAXIMUM: Error, array does not have a maximum'
          
       END IF
       
    END DO

  END FUNCTION maximum

  SUBROUTINE mask(okay,m,n,min,max)

    ! Flags objects that make the cut as 'okay'
    ! Can be applied to any scalar array, not just mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m(n), min, max
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(OUT) :: okay(n)
    INTEGER :: i, o

    WRITE(*,*) 'MASK: Imposing property cut'    
    WRITE(*,*) 'MASK: Minimum value:', min
    WRITE(*,*) 'MASK: Maximum value:', max
    WRITE(*,*) 'MASK: Original number of objects', n

    okay=.FALSE.

    DO i=1,n
       IF(m(i)>=min .AND. m(i)<=max) okay(i)=.TRUE.
    END DO

    o=COUNT(okay)

    WRITE(*,*) 'MASK: Final number of objects:', o
    WRITE(*,*) 'MASK: Fraction remaining:', real(o)/real(n)
    WRITE(*,*) 'MASK: Fraction culled:', 1.-real(o)/real(n)
    WRITE(*,*) 'MASK: Done'
    WRITE(*,*)

  END SUBROUTINE mask

  SUBROUTINE add_to_array(a,n,b)

    ! Append b to array a(n) to make a new array a(n+1)
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: b
    INTEGER :: hold(n)
    INTEGER :: i

    hold=a
    DEALLOCATE(a)
    ALLOCATE(a(n+1))

    DO i=1,n
       a(i)=hold(i)
    END DO
    a(n+1)=b

  END SUBROUTINE add_to_array

  INTEGER FUNCTION unique_entries(a,n)

    ! Counts the total number of unique entries in array 'a'
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: a(n) ! Input array
    INTEGER, INTENT(IN) :: n    ! Size of input array
    INTEGER :: i, j

    ! Initially assume all entries are unique
    unique_entries=n

    ! Double loop to each each pair against each other
    DO i=1,n
       DO j=i+1,n
          IF(a(j)==a(i)) THEN
             unique_entries=unique_entries-1 ! A non-unique entry has been discovered, so subtract one
             EXIT
          END IF
       END DO
    END DO

  END FUNCTION unique_entries

  SUBROUTINE unique_index(array,n,unique,m,match)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: array(n) ! Array to find the unique indices of
    INTEGER, INTENT(IN) :: n        ! Number of entries in input array
    INTEGER, ALLOCATABLE, INTENT(OUT) :: unique(:) ! Output array of unique indices
    INTEGER, INTENT(OUT) :: m                      ! Number of unique indices
    INTEGER, INTENT(OUT) :: match(n)               ! Array for matching input and unique arrays
    INTEGER :: i, j, p
    LOGICAL :: increment
    
    ! First count the number of unique entries
    m=unique_entries(array,n)
    ALLOCATE(unique(m))

    ! Set the first unique entry to be the first entry
    unique(1)=array(1)

    p=1
    DO i=1,n
       increment=.FALSE.
       DO j=1,p
          IF(array(i) /= unique(j)) THEN
             unique(p+1)=array(i)
             increment=.TRUE.
             EXIT
          END IF
       END DO
       IF(increment) p=p+1
       IF(p==m) EXIT ! Added this in haste
    END DO

    ! Now fill the matching array
    DO j=1,m
       DO i=1,n
          IF(unique(j)==array(i)) THEN
             match(i)=j
          END IF
       END DO
    END DO
  
  END SUBROUTINE unique_index

!!$  LOGICAL FUNCTION in_array(b,a,n)
!!$
!!$    ! Test to see if b is in a(n)
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: b
!!$    INTEGER, INTENT(IN) :: a(n)
!!$    INTEGER, INTENT(IN) :: n
!!$    INTEGER :: i
!!$
!!$    DO i=1,n
!!$       IF(b==a(i)) THEN
!!$          in_array=.TRUE.
!!$       END IF
!!$       EXIT
!!$    END DO
!!$
!!$  END FUNCTION in_array

  LOGICAL FUNCTION repeated_entries(a,n)

    ! Checks for repeated entries in a(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j

    repeated_entries=.FALSE.
    
    DO i=1,n
       DO j=i+1,n
          IF(a(i)==a(j)) THEN
             repeated_entries=.TRUE.
             EXIT
          END IF
       END DO
       IF(repeated_entries) EXIT
    END DO
    
  END FUNCTION repeated_entries

END MODULE array_operations
