MODULE tipsy_stuff

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_tipsy

   CONTAINS

   SUBROUTINE read_tipsy(filename, x, v, m, n, a)

      IMPLICIT NONE
      CHARACTER(len=256), INTENT(IN) :: filename
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :), v(:, :)
      REAL, INTENT(OUT) :: m, a
      INTEGER, INTENT(OUT) :: n
      REAL :: spam
      INTEGER :: i, j

      ! Read header
      OPEN(7, file=filename)
      READ(7, *) n
      READ(7, *) spam
      READ(7, *) a
      WRITE(*, *) 'READ_TIPSY: Scale factor:', a
  
      ! Read positions and velocities
      ALLOCATE(x(3, n), v(3, n))
      x = 0.; v = 0.
      DO j = 1, 7
         DO i = 1, n
            IF(j==1) READ(7, *) m
            IF(j==2) READ(7, *) x(1, i)
            IF(j==3) READ(7, *) x(2, i)
            IF(j==4) READ(7, *) x(3, i)
            IF(j==5) READ(7, *) v(1, i)
            IF(j==6) READ(7, *) v(2, i)
            IF(j==7) READ(7, *) v(3, i)
         END DO
      END DO
      CLOSE(7)

      ! Convert units
      x = x/1000.
      v = v*sqrt(a)
      m = m*1e10
      WRITE(*, *) 'READ_TIPSY: Particle mass [log10(M_sun/h)]:', log10(m)
      WRITE(*, *) 'READ_TIPSY: Done'
      WRITE(*, *)

   END SUBROUTINE read_tipsy

END MODULE tipsy_stuff