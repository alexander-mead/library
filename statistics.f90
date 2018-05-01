MODULE statistics

  USE table_integer
  USE array_operations
  
CONTAINS

  FUNCTION mean(x,n)

    IMPLICIT NONE
    REAL :: mean
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x(n)
    DOUBLE PRECISION :: sum
    INTEGER :: i

    sum=0.d0
    DO i=1,n
       sum=sum+x(i)
    END DO

    mean=REAL(sum)/REAL(n)

  END FUNCTION mean

  FUNCTION variance(x,n)

    IMPLICIT NONE
    REAL :: variance
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x(n)
    DOUBLE PRECISION :: sum
    REAL :: avg
    INTEGER :: i

    avg=mean(x,n)

    sum=0.d0
    DO i=1,n
       sum=sum+(x(i)-avg)**2
    END DO

    variance=REAL(sum)/REAL(n)

  END FUNCTION variance

  SUBROUTINE histogram(xmin,xmax,x,hist,n,data,m)

    IMPLICIT NONE
    REAL, INTENT(IN) :: xmin, xmax, data(m)
    INTEGER, INTENT(IN) :: n, m
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: hist(:)
    INTEGER :: i, j

    WRITE(*,*) 'HISTOGRAM: Assiging arrays'

    !Fill the table for the xrange and allocate the histogram array
    CALL fill_array(xmin,xmax,x,n+1)
    ALLOCATE(hist(n))
    hist=0

    WRITE(*,*) 'HISTOGRAM: Constructing histogram'

    !Make the histogram from the data
    DO i=1,m
       IF(data(i)<xmin) THEN
          CYCLE
       ELSE IF(data(i)>xmax) THEN
          CYCLE
       ELSE
          j=select_table_integer(data(i),x,n,1)
          hist(j)=hist(j)+1
       END IF
    END DO

    WRITE(*,*) 'HISTOGRAM: Fraction of data assigned to histogram:', REAL(SUM(hist))/REAL(m)
    WRITE(*,*) 'HISTOGRAM: Done'
    WRITE(*,*)

  END SUBROUTINE histogram

END MODULE statistics
