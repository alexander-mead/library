MODULE gadget

   USE precision
   USE array_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_gadget
   PUBLIC :: write_gadget
   PUBLIC :: read_catalogue
   PUBLIC :: write_catalogue

   REAL, PARAMETER :: Lunit = 1e-3  ! Convert from Gadget2 kpc/h to Mpc/h
   REAL, PARAMETER :: Munit = 1e10  ! Convert from Gadget2 10^10 Msun/h to Msun/h

CONTAINS

   SUBROUTINE read_gadget(x, v, id, L, Om_m, Om_v, h, m, a, z, n, infile)

      IMPLICIT NONE
      REAL*4, ALLOCATABLE, INTENT(OUT) :: x(:, :)
      REAL*4, ALLOCATABLE, INTENT(OUT) :: v(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
      REAL*8, INTENT(OUT) :: L
      REAL*8, INTENT(OUT) :: Om_m
      REAL*8, INTENT(OUT) :: Om_v
      REAL*8, INTENT(OUT) :: h
      REAL*8, INTENT(OUT) :: m
      REAL*8, INTENT(OUT) :: a
      REAL*8, INTENT(OUT) :: z
      INTEGER, INTENT(OUT) :: n
      CHARACTER(len=*), INTENT(IN) :: infile
      REAL*8 :: mass(6)!, z8, a8, L8, om_m8, om_v8, h8
      INTEGER :: np(6), np2(6), crap
      LOGICAL :: lexist

      WRITE (*, *) 'READ_GADGET: Reading in Gadget-2 file: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) STOP 'READ_GADGET: Error, input file does not exist'

      OPEN (7, file=infile, form='unformatted', status='old')
      READ (7) np, mass, a, z, crap, crap, np2, crap, crap, L, Om_m, Om_v, h
      CLOSE (7)

      ! Convert Gadget doubles to my reals
      !a = real(a8)
      !z = real(z8)
      !om_m = real(om_m8)
      !om_v = real(om_v8)
      !h = real(h8)
      !L = real(L8)/Lunit

      ! Multiply the masses by 1e10 to get in units of M_sun/h
      m = mass(2)*Munit
      L = L*Lunit
      WRITE (*, *) 'READ_GADGET: Particle number:', np(2)
      WRITE (*, *) 'READ_GADGET: Which is:', nint(np(2)**(1./3.)), 'cubed.'
      WRITE (*, *) 'READ_GADGET: Particle mass log10([M_sun/h]):', log10(m)
      WRITE (*, *) 'READ_GADGET: Box size [Mpc/h]:', L
      WRITE (*, *) 'READ_GADGET: a:', a
      WRITE (*, *) 'READ_GADGET: z:', z
      WRITE (*, *) 'READ_GADGET: Om_m:', Om_m
      WRITE (*, *) 'READ_GADGET: Om_v:', Om_v

      ! Fix the total number of simulation particles
      n = np(2)

      ! Allocate arrays for position, velocity and particle ID number
      ALLOCATE (x(3, n), v(3, n), id(n))

      ! Read in the binary data, skip the header line
      OPEN (7, file=infile, form='unformatted', status='old')
      READ (7)
      READ (7) x
      READ (7) v
      READ (7) id
      CLOSE (7)

      ! kpc -> Mpc conversion
      x = x*real(Lunit)

      ! Change from weird Gadget units to peculiar velocities
      v = v*sqrt(real(a))

      WRITE (*, *) 'READ_GADGET: Finished reading in file'
      WRITE (*, *)

   END SUBROUTINE read_gadget

   SUBROUTINE write_gadget(x, v, id, L, Om_m, Om_v, h, m, a, z, n, outfile)

      IMPLICIT NONE
      REAL*4, INTENT(IN) :: x(3, n)
      REAL*4, INTENT(IN) :: v(3, n)
      INTEGER, INTENT(IN) :: id(n)
      REAL*8, INTENT(IN) :: L
      REAL*8, INTENT(IN) :: Om_m
      REAL*8, INTENT(IN) :: Om_v
      REAL*8, INTENT(IN) :: h
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: a
      REAL*8, INTENT(IN) :: z
      INTEGER, INTENT(IN) :: n
      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL*8 :: mass(6), crap(12)
      INTEGER :: np(6), crapi

      WRITE (*, *) 'WRITE_GADGET: Outputting particle data in Gadget2 format: ', trim(outfile)

      np = 0
      np(2) = n
      mass = 0.
      mass(2) = m/Munit ! Convert mass to Gadget units
      !a8 = a
      !z8 = z
      !crapi = 0
      crap = 0.
      !om_m8 = om_m
      !om_v8 = om_v
      !h8 = h
      !L8 = L*Lunit
      !Lout = L/Lunit

      WRITE (*, *) 'WRITE_GADGET: Particle number:', n
      WRITE (*, *) 'WRITE_GADGET: Which is:', nint(n**(1./3.)), 'cubed'
      WRITE (*, *) 'WRITE_GADGET: Box size [Mpc/h]:', L
      WRITE (*, *) 'WRITE_GADGET: a:', a
      WRITE (*, *) 'WRITE_GADGET: z:', z
      WRITE (*, *) 'WRITE_GADGET: Particle mass [Msun/h]:', m
      WRITE (*, *) 'WRITE_GADGET: Om_m:', Om_m
      WRITE (*, *) 'WRITE_GADGET: Om_v:', Om_v

      OPEN (7, file=outfile, form='unformatted', status='replace')
      WRITE (7) np, mass, a, z, crapi, crapi, np, crapi, crapi, L/Lunit, Om_m, Om_v, h, crap
      WRITE (7) x/Lunit
      WRITE (7) v/sqrt(a)
      WRITE (7) id
      CLOSE (7)

      WRITE (*, *) 'WRITE_GADGET: Finished writing file'
      WRITE (*, *)

   END SUBROUTINE write_gadget

   SUBROUTINE read_catalogue(x, v, m, npart, disp, c, env, Dv, rmax, avg_r, rms_r, n, infile)

      USE file_info
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :), v(:, :), m(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: disp(:), c(:), env(:), Dv(:), rmax(:), avg_r(:), rms_r(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: npart(:)
      INTEGER, INTENT(OUT) :: n
      CHARACTER(len=*), INTENT(IN) :: infile
      INTEGER :: i
      LOGICAL :: lexist

      STOP 'READ_CATALOGUE: need to change this for single/double precision'

      WRITE (*, *) 'READ_CATALOGUE: Reading in catalogue: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) STOP 'READ_CATALOGUE: Error, input file does not exist'

      n = file_length(infile, verbose=.FALSE.)
      ALLOCATE (x(3, n), v(3, n), m(n), npart(n))
      ALLOCATE (disp(n), c(n), env(n), Dv(n), rmax(n), avg_r(n), rms_r(n))

      OPEN (7, file=infile)
      DO i = 1, n
         READ (7, *) x(1, i), x(2, i), x(3, i), v(1, i), v(2, i), v(3, i), &
            m(i), npart(i), disp(i), c(i), env(i), Dv(i), rmax(i), avg_r(i), rms_r(i)
      END DO
      CLOSE (7)

      WRITE (*, *) 'READ_CATALOGUE: Min x [Mpc/h]    :', minval(x)
      WRITE (*, *) 'READ_CATALOGUE: Max x [Mpc/h]    :', maxval(x)
      WRITE (*, *) 'READ_CATALOGUE: Min v [km/s]     :', minval(v)
      WRITE (*, *) 'READ_CATALOGUE: Max v [km/s]     :', maxval(v)
      WRITE (*, *) 'READ_CATALOGUE: Min mass [Msun/h]:', minval(m)
      WRITE (*, *) 'READ_CATALOGUE: Max mass [Msun/h]:', maxval(m)
      WRITE (*, *) 'READ_CATALOGUE: Finished reading catalogue file'
      WRITE (*, *)

   END SUBROUTINE read_catalogue

   SUBROUTINE write_catalogue(x, v, m, npart, disp, c, env, Dv, rmax, avg_r, rms_r, n, outfile)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL, INTENT(IN) :: x(3, n), v(3, n), m(n)
      REAL, INTENT(IN) :: disp(n), c(n), env(n), Dv(n), rmax(n), avg_r(n), rms_r(n)
      INTEGER, INTENT(IN) :: npart(n)
      INTEGER :: i

      STOP 'WRITE_CATALOGUE: need to change this for single/double precision'

      WRITE (*, *) 'WRITE_CATALOGUE: Outputting catalogue'

      OPEN (7, file=outfile)
      DO i = 1, n
         WRITE (7, *) x(1, i), x(2, i), x(3, i), v(1, i), v(2, i), v(3, i), &
            m(i), npart(i), disp(i), c(i), env(i), Dv(i), rmax(i), avg_r(i), rms_r(i)
      END DO
      CLOSE (7)

      WRITE (*, *) 'WRITE_CATALOGUE: Finished writing catalogue file'
      WRITE (*, *)

   END SUBROUTINE write_catalogue

END MODULE gadget
