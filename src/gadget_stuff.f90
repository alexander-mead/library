MODULE gadget_stuff

   USE precision
   USE array_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: write_gadget_simulation_format_Pk
   PUBLIC :: read_gadget
   PUBLIC :: write_gadget
   PUBLIC :: read_catalogue
   PUBLIC :: write_catalogue

   REAL, PARAMETER :: Lunit = 1e-3  ! Convert from Gadget2 kpc/h to Mpc/h
   REAL, PARAMETER :: Munit = 1e10  ! Convert from Gadget2 10^10 Msun/h to Msun/h

   ! Gadget file format precisions
   INTEGER, PARAMETER :: gadget_sp = 4
   INTEGER, PARAMETER :: gadget_dp = 8

CONTAINS

   SUBROUTINE write_gadget_simulation_format_Pk(k, Pk, outfile, verbose)

      ! Converts a CAMB P(k) file to a Gadget format P(k) file 
      REAL, INTENT(IN) :: k(:)  ! Input k [h/Mpc]
      REAL, INTENT(IN) :: Pk(:) ! Input dimensionless Delta^2(k)
      CHARACTER(len=*), INTENT(IN) :: outfile ! Output file
      LOGICAL, INTENT(IN) :: verbose ! Verbose or not
      INTEGER :: i, nk
      REAL, ALLOCATABLE :: Gadget_k(:), Gadget_Pk(:)

      nk = size(k)
      IF (nk /= size(Pk)) STOP 'WRITE_GADGET_SIMULATION_FORMAT: Error, k and Pk must be the same size'
      ALLOCATE(Gadget_k(nk), Gadget_Pk(nk))

      IF(verbose) WRITE(*, *) 'WRITE_GADGET_INPUT_POWER: Converting units to Gadget format'

      Gadget_k = log10(k)-3. ! Convert Mpc/h to kpc/h and take log10
      Gadget_Pk = log10(Pk)  ! Take log10

      IF(verbose) WRITE(*, *) 'WRITE_GADGET_INPUT_POWER: Output file: ', trim(outfile)
      OPEN(7, file=outfile)
      DO i = 1, nk
         WRITE(7,*) Gadget_k(i), Gadget_Pk(i)
      END DO
      CLOSE(7)
      IF(verbose) THEN
         WRITE(*, *) 'WRITE_GADGET_INPUT_POWER: Done'
         WRITE(*, *)
      END IF

   END SUBROUTINE write_gadget_simulation_format_Pk

   SUBROUTINE read_gadget(x, v, id, L, Om_m, Om_v, h, m, a, z, n, infile)

      ! Read in particle data from a Gadget formatted binary file
      USE precision
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
      REAL, INTENT(OUT) :: L
      REAL, INTENT(OUT) :: Om_m
      REAL, INTENT(OUT) :: Om_v
      REAL, INTENT(OUT) :: h
      REAL, INTENT(OUT) :: m
      REAL, INTENT(OUT) :: a
      REAL, INTENT(OUT) :: z
      INTEGER, INTENT(OUT) :: n
      CHARACTER(len=*), INTENT(IN) :: infile
      REAL(gadget_dp) :: mass_in(6), L_in, Om_m_in, Om_v_in, h_in, a_in, z_in
      REAL(gadget_sp), ALLOCATABLE :: x_in(:, :), v_in(:, :)
      INTEGER :: np(6), craps(6), crap
      LOGICAL :: lexist

      WRITE (*, *) 'READ_GADGET: Reading in Gadget-2 file: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) STOP 'READ_GADGET: Error, input file does not exist'

      OPEN (7, file=infile, form='unformatted', status='old')
      READ (7) np, mass_in, a_in, z_in, crap, crap, craps, crap, crap, L_in, Om_m_in, Om_v_in, h_in
      CLOSE (7)

      ! Convert Gadget header numbers to my numbers
      ! TODO: Possible precision change
      m = mass_in(2)
      a = a_in
      z = z_in
      Om_m = Om_m_in
      Om_v = Om_v_in
      h = h_in
      L = L_in
      n = np(2)

      ! Convert units
      m = m*Munit
      L = L*Lunit

      WRITE (*, *) 'READ_GADGET: Particle number:', n
      WRITE (*, *) 'READ_GADGET: Which is:', nint(n**(1./3.)), 'cubed.'
      WRITE (*, *) 'READ_GADGET: Particle mass log10([M_sun/h]):', log10(m)
      WRITE (*, *) 'READ_GADGET: Box size [Mpc/h]:', L
      WRITE (*, *) 'READ_GADGET: a:', a
      WRITE (*, *) 'READ_GADGET: z:', z
      WRITE (*, *) 'READ_GADGET: Om_m:', Om_m
      WRITE (*, *) 'READ_GADGET: Om_v:', Om_v
      WRITE (*, *) 'READ_GADGET: h:', h

      ! Allocate arrays for position, velocity and particle ID number
      ALLOCATE (x_in(3, n), v_in(3, n), id(n))

      ! Read in the binary data, skip the header line
      OPEN (7, file=infile, form='unformatted', status='old')
      READ (7)
      READ (7) x_in
      READ (7) v_in
      READ (7) id
      CLOSE (7)

      ! Convert precision
      ALLOCATE(x(3, n))
      x = x_in
      DEALLOCATE(x_in)
      
      ! Convert precision
      ALLOCATE(v(3, n))
      v = v_in
      DEALLOCATE(v_in)

      ! Convert units
      x = x*Lunit
      v = v*sqrt(a)

      WRITE (*, *) 'READ_GADGET: Finished reading file'
      WRITE (*, *)

   END SUBROUTINE read_gadget

   SUBROUTINE write_gadget(x, v, id, L, Om_m, Om_v, h, m, a, z, outfile)

      ! Write particle data to a Gadget formatted particle file
      REAL, INTENT(IN) :: x(:, :)
      REAL, INTENT(IN) :: v(:, :)
      INTEGER, INTENT(IN) :: id(:)
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: Om_m
      REAL, INTENT(IN) :: Om_v
      REAL, INTENT(IN) :: h
      REAL, INTENT(IN) :: m
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: z
      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL(dp) :: mass(6), crap(12), Om_m_out, Om_v_out, h_out, L_out, a_out, z_out
      INTEGER :: np(6), crapi, n

      IF(size(x, 1) /= 3 .OR. size(v, 1) /= 3 ) STOP 'WRITE_GADGET: Error, input x, v arrays should be 3D'
      n = size(x, 2)
      IF (n /= size(x, 2) .OR. n /= size(v, 2) .OR. n /= size(id)) THEN
         STOP 'WRITE_GADGET: Error, x, v, id do not have the same number of entries'
      END IF

      WRITE (*, *) 'WRITE_GADGET: Outputting particle data in Gadget2 format: ', trim(outfile)

      np = 0
      np(2) = n
      mass = 0.
      mass(2) = real(m/Munit, gadget_dp) ! Convert mass to Gadget units
      a_out = real(a, gadget_dp)
      z_out = real(z, gadget_dp)
      crapi = 0
      crap = 0.
      Om_m_out = real(Om_m, gadget_dp)
      Om_v_out = real(Om_v, gadget_dp)
      h_out = real(h, gadget_dp)
      L_out = real(L/Lunit, gadget_dp)

      WRITE (*, *) 'WRITE_GADGET: Particle number:', n
      WRITE (*, *) 'WRITE_GADGET: Which is:', nint(n**(1./3.)), 'cubed'
      WRITE (*, *) 'WRITE_GADGET: Box size [Mpc/h]:', L
      WRITE (*, *) 'WRITE_GADGET: a:', a
      WRITE (*, *) 'WRITE_GADGET: z:', z
      WRITE (*, *) 'WRITE_GADGET: Particle mass log10([Msun/h]):', log10(m)
      WRITE (*, *) 'WRITE_GADGET: Om_m:', Om_m
      WRITE (*, *) 'WRITE_GADGET: Om_v:', Om_v
      WRITE (*, *) 'WRITE_GADGET: h:', h

      OPEN (7, file=outfile, form='unformatted', status='replace')
      WRITE (7) np, mass, a_out, z_out, crapi, crapi, np, crapi, crapi, L_out, Om_m_out, Om_v_out, h_out, crap
      WRITE (7) real(x/Lunit, gadget_sp)
      WRITE (7) real(v/sqrt(a), gadget_sp)
      WRITE (7) id
      CLOSE (7)

      WRITE (*, *) 'WRITE_GADGET: Finished writing file'
      WRITE (*, *)

   END SUBROUTINE write_gadget

   SUBROUTINE read_catalogue(x, v, m, npart, disp, c, env, Dv, rmax, avg_r, rms_r, n, infile)

      USE file_info
      REAL(sp), ALLOCATABLE, INTENT(OUT) :: x(:, :)
      REAL(sp), ALLOCATABLE, INTENT(OUT) :: v(:, :)
      REAL(sp), ALLOCATABLE, INTENT(OUT) :: m(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: npart(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: disp(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: c(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: env(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Dv(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: rmax(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: avg_r(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: rms_r(:) 
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

   SUBROUTINE write_catalogue(x, v, m, npart, disp, c, env, Dv, rmax, avg_r, rms_r, outfile)

      REAL, INTENT(IN) :: x(:, :)
      REAL, INTENT(IN) :: v(:, :)
      REAL, INTENT(IN) :: m(:)
      INTEGER, INTENT(IN) :: npart(:)
      REAL, INTENT(IN) :: disp(:)
      REAL, INTENT(IN) :: c(:)
      REAL, INTENT(IN) :: env(:)
      REAL, INTENT(IN) :: Dv(:)
      REAL, INTENT(IN) :: rmax(:)
      REAL, INTENT(IN) :: avg_r(:)
      REAL, INTENT(IN) :: rms_r(:)
      CHARACTER(len=*), INTENT(IN) :: outfile
      INTEGER :: i, n

      STOP 'WRITE_CATALOGUE: Need to change this for single/double precision'
      STOP 'WRITE_CATALOGUE: Add checks that sizes of catalogue entries are all the same size'

      IF(size(x,1) /= 3 .OR. size(v,1) /= 3 ) STOP 'WRITE_CATALOGUE: Error, input x, v arrays should be 3D'
      n = size(x,2)
      IF (n /= size(x,2) .OR. n /= size(v,2)) THEN
         STOP 'WRITE_CATALOGUE: Error, x, v, id do not have the same number of entries'
      END IF

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

END MODULE gadget_stuff
