PROGRAM down_sample_test

  USE logical_operations

  IMPLICIT NONE
  REAL :: x, y, z
  INTEGER :: j, k, n
  INTEGER :: ix, iy, iz

  WRITE(*,*) 'Down sampling test'
  WRITE(*,*)

  WRITE(*,*) 'Cube root number of points in grid (e.g. 8)'
  READ(*,*) n
  WRITE(*,*)

  k=0

  ix=0
  iy=0
  iz=0

  OPEN(7,file='cube.dat')
  OPEN(8,file='accepted.dat')
  DO j=1,n**3

     IF(mod(j-1,n**2)==0) THEN
        ix=1
        iy=1
        iz=iz+1
     ELSE IF(mod(j-1,n)==0) THEN
        ix=1
        iy=iy+1
     ELSE
        ix=ix+1
     END IF

     x=(float(ix)-0.5)/float(n)
     y=(float(iy)-0.5)/float(n)
     z=(float(iz)-0.5)/float(n)

     WRITE(7,*) x, y, z, j 

     IF(odd(CEILING(float(j)/float(n**2)))) THEN
        IF(odd(CEILING(float(j)/float(n)))) THEN
           IF(odd(CEILING(float(j)))) THEN

              k=k+1

              WRITE(*,*) 'Accepted:', j

              WRITE(8,*) x, y, z, j

           END IF
        END IF
     END IF

  END DO
  CLOSE(7)
  CLOSE(8)

  WRITE(*,*)
  WRITE(*,*) 'Total points:', n**3
  WRITE(*,*) 'Total accepted:', k
  WRITE(*,*) 'Fractional acceptance:', float(k)/float(n**3)
  WRITE(*,*)

END PROGRAM down_sample_test
