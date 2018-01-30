PROGRAM adaptive_grid

  IMPLICIT NONE
  REAL :: xmin(3), lo(3), hi(3), ref(3)
  REAL :: a, b, c
  INTEGER :: ngrid(3), nref

  a=3.
  b=4.
  c=5.

  xmin=1.
  lo=-10.
  hi=10.
  ref=5.
  ngrid=10

  WRITE(*,*)
  WRITE(*,*) 'Adaptive grid search algorithm'
  WRITE(*,*)

  CALL adapt3(xmin,quadratic,lo,hi,ref,ngrid,5)

  WRITE(*,*) 'Best fit minimum at:', xmin
  WRITE(*,*) 'True miniumum at:', a, b, c
  WRITE(*,*)

CONTAINS

  FUNCTION quadratic(x)

    IMPLICIT NONE
    REAL :: x(3), quadratic

    quadratic=(x(1)-a)**2.+(x(2)-b)**2.+(x(3)-c)**2.

  END FUNCTION quadratic

  SUBROUTINE adapt3(xmin,func,min,max,ref,ngrid,nref)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: xmin(3)
    REAL, INTENT(IN) :: min(3), max(3), ref(3)
    INTEGER, INTENT(IN) :: ngrid(3), nref
    REAL, EXTERNAL :: func
    REAL :: fom, fommin, range(3), x(3)
    REAL :: x1(3), x2(3)
    INTEGER :: i, j, k, l, imin, jmin, kmin

    !xmin is output
    !min/max give you the initial bounds for the search
    !ref is the refinement level for each separate grid search
    !n is the number of points per grid search
    !m is the number of refinement levels

    x1=min
    x2=max

    DO l=1,nref

       DO i=1,ngrid(1)

          x(1)=x1(1)+(x2(1)-x1(1))*float(i-1)/float(ngrid(1)-1)

          DO j=1,ngrid(2)

             x(2)=x1(2)+(x2(2)-x1(2))*float(j-1)/float(ngrid(2)-1)
            
             DO k=1,ngrid(3)
                
                x(3)=x1(3)+(x2(3)-x1(3))*float(k-1)/float(ngrid(3)-1)

                fom=func(x)
                
                IF((i==1 .AND. j==1 .AND. k==1) .OR. fom<fommin) THEN
                   fommin=fom
                   xmin=x
                   imin=i
                   jmin=j
                   kmin=k
                END IF

             END DO
             
          END DO

       END DO

       IF(imin==1 .OR. imin==ngrid(1)) STOP 'Error: best fit x(1) at edge of range'
       IF(jmin==1 .OR. jmin==ngrid(2)) STOP 'Error: best fit x(2) at edge of range'
       IF(kmin==1 .OR. kmin==ngrid(3)) STOP 'Error: best fit x(3) at edge of range'

       IF(l .NE. nref) THEN
          range=x2-x1
          x1=xmin-range/ref
          x2=xmin+range/ref
       END IF

    END DO

  END SUBROUTINE adapt3

END PROGRAM
