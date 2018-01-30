PROGRAM adaptive_grid

  IMPLICIT NONE
  REAL :: xmin
  REAL :: a, b, c

  a=1.
  b=-4.
  c=2.

  WRITE(*,*)
  WRITE(*,*) 'Adaptive grid search algorithm'
  WRITE(*,*)

  CALL adapt(xmin,quadratic,-10.,10.,10.,10,5)

  WRITE(*,*) 'Best fit minimum at:', xmin
  WRITE(*,*) 'True miniumum at:', quadratic_minimum()
  WRITE(*,*)

CONTAINS

  FUNCTION quadratic(x)

    IMPLICIT NONE
    REAL :: x, quadratic

    quadratic=a*(x**2.)+b*x+c

  END FUNCTION quadratic

  FUNCTION quadratic_minimum()

    IMPLICIT NONE
    REAL :: quadratic_minimum
    
    quadratic_minimum=-b/(2.*a)

  END FUNCTION quadratic_minimum

  SUBROUTINE adapt(xmin,func,min,max,ref,ngrid,nref)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: xmin
    REAL, INTENT(IN) :: min, max, ref
    INTEGER, INTENT(IN) :: ngrid, nref
    REAL, EXTERNAL :: func
    REAL :: fom, fommin, range, x
    REAL :: x1, x2
    INTEGER :: i, j, imin

    !xmin is output
    !min/max give you the initial bounds for the search
    !ref is the refinement level for each separate grid search
    !n is the number of points per grid search
    !m is the number of refinement levels

    x1=min
    x2=max

    DO j=1,nref

       DO i=1,ngrid

          x=x1+(x2-x1)*float(i-1)/float(ngrid-1)
          fom=func(x)

!          WRITE(*,*) j, x

          IF(i==1 .OR. fom<fommin) THEN
             fommin=fom
             xmin=x
             imin=i
          END IF

       END DO

       IF(imin==1 .OR. imin==ngrid) STOP 'Error: best fit at edge of range'

       IF(j .NE. nref) THEN
          range=x2-x1
          x1=xmin-range/ref
          x2=xmin+range/ref
       END IF

    END DO

  END SUBROUTINE adapt

END PROGRAM
