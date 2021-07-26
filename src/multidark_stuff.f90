MODULE multidark_stuff

   USE array_operations
   USE string_operations
   USE file_info
   USE table_integer

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_multidark_halo_catalogue
   PUBLIC :: read_multidark_particles
   PUBLIC :: write_multidark_halo_catalogue
   PUBLIC :: multidark_snapshot
   PUBLIC :: multidark_scale_factor
   PUBLIC :: nearest_multidark_snapshot
   PUBLIC :: nearest_multidark_snapshot_sig8

   PUBLIC :: column_Mvir
   PUBLIC :: column_Mtot
   PUBLIC :: column_M200
   PUBLIC :: column_M200c
   PUBLIC :: column_M500c
   PUBLIC :: column_M2500c

   ! Available Multidark scale factors
   REAL, PARAMETER :: as(35) = [0.257, 0.287, 0.318, 0.348, 0.378, 0.409, 0.439, &
                                 0.470, 0.500, 0.530, 0.561, 0.591, 0.621, 0.652, &
                                 0.682, 0.713, 0.728, 0.743, 0.758, 0.773, 0.788, &
                                 0.804, 0.819, 0.834, 0.849, 0.864, 0.880, 0.895, &
                                 0.910, 0.925, 0.940, 0.956, 0.971, 0.986, 1.001]

   ! Snapshots corresponding to scale factors
   INTEGER, PARAMETER :: snaps(35) = [36, 38, 40, 42, 44, 46, 48, &
                                       50, 52, 54, 56, 58, 60, 62, &
                                       64, 66, 67, 68, 69, 70, 71, &
                                       72, 73, 74, 75, 76, 77, 78, &
                                       79, 80, 81, 82, 83, 84, 85]

   ! sigma8 values corresponding to scale factors
   REAL, PARAMETER :: sig8s(35) = [0.275, 0.306, 0.338, 0.368, 0.398, 0.428, 0.456, &
                                    0.484, 0.511, 0.536, 0.562, 0.586, 0.609, 0.631, &
                                    0.652, 0.673, 0.682, 0.692, 0.701, 0.710, 0.718, &
                                    0.727, 0.736, 0.744, 0.752, 0.759, 0.767, 0.774, &
                                    0.781, 0.788, 0.795, 0.802, 0.808, 0.814, 0.820]

   ! Halo catalogue files
   ! TODO: Note that Rockstar catalogues have 6 mass entries, while BDMV only have 2 mass entries
   ! NOTE: Only have BDMV catalogues for Bolshoi
   INTEGER, PARAMETER :: column_Mvir = 1
   INTEGER, PARAMETER :: column_Mtot = 2
   INTEGER, PARAMETER :: column_M200 = 3
   INTEGER, PARAMETER :: column_M200c = 4
   INTEGER, PARAMETER :: column_M500c = 5
   INTEGER, PARAMETER :: column_M2500c = 6 

   ! Snapshot finding
   INTEGER, PARAMETER :: ifind_snap = ifind_split
   REAL, PARAMETER :: eps_snap = 0.

CONTAINS

   ! SUBROUTINE read_multidark_haloes(infile, mmin, x, m, n)

   !    ! TODO: Remove, this is unnecessary if you use the SQL interface
   !    IMPLICIT NONE
   !    CHARACTER(len=*), INTENT(IN) :: infile      ! File to read in
   !    REAL, INTENT(IN) :: mmin                    ! Minimum halo virial mass [Msun/h]
   !    REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :)    ! Position array [Mpc/h]
   !    REAL, ALLOCATABLE, INTENT(OUT) :: m(:, :)    ! Virial halo mass [Msun/h]
   !    INTEGER, INTENT(OUT) :: n                   ! Total number of haloes
   !    LOGICAL :: lexist
   !    INTEGER :: i, j
   !    INTEGER :: p, pid
   !    REAL :: mm!, c
   !    REAL, ALLOCATABLE :: data(:)

   !    INTEGER, PARAMETER :: hash_lines = 58 ! Number of lines beginning with # (Both Multidark and Bolshoi have 58)
   !    !REAL, PARAMETER :: mmin=1.74e12     ! Minimum halo mass [Msun/h] (corresponds to N>200 for MDR1; 1e12 would mean n~115)
   !    !INTEGER, PARAMETER :: columns=73       ! Total number of columns in file
   !    INTEGER, PARAMETER :: columns = 41       ! Total number of columns to read from file (weird errors if this is set to 73)
   !    INTEGER, PARAMETER :: column_pid = 6     ! Column for PID (-1 if unique halo)
   !    INTEGER, PARAMETER :: column_mv = 11     ! Column for virial mass (unbinding done) [Msun/h]
   !    INTEGER, PARAMETER :: column_x = 18      ! Column for x position [Mpc/h]
   !    INTEGER, PARAMETER :: column_y = 19      ! Column for y position [Mpc/h]
   !    INTEGER, PARAMETER :: column_z = 20      ! Column for z position [Mpc/h]
   !    INTEGER, PARAMETER :: column_vx = 21     ! Column for vx position [km/s]
   !    INTEGER, PARAMETER :: column_vy = 22     ! Column for vy position [km/s]
   !    INTEGER, PARAMETER :: column_vz = 23     ! Column for vz position [km/s]
   !    INTEGER, PARAMETER :: column_mvu = 37    ! Column for total virial mass (no unbinding) [Msun/h]
   !    INTEGER, PARAMETER :: column_m200 = 38   ! Column for M200 [Msun/h]
   !    INTEGER, PARAMETER :: column_m200c = 39  ! Column for M200 critical [Msun/h]
   !    INTEGER, PARAMETER :: column_m500c = 40  ! Column for M500 critical [Msun/h]
   !    INTEGER, PARAMETER :: column_m2500c = 41 ! Column for M2500 critical [Msun/h]

   !    ! Check file exists
   !    INQUIRE (file=infile, exist=lexist)
   !    IF (.NOT. lexist) THEN
   !       WRITE (*, *) 'READ_MULTIDARK_HALOES: Error, catalogue file does not exist: ', trim(infile)
   !       STOP
   !    END IF

   !    ! Welcome message
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Reading in halo catalogue: ', trim(infile)

   !    ! Find file length
   !    n = file_length(infile, verbose=.FALSE.)
   !    n = n-hash_lines

   !    ! Allocate the data array with the total number of columns so that you can read in a full line
   !    ALLOCATE (data(columns))

   !    ! Count unique haloes
   !    p = 0 ! Set sum variable to zero
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Total number of haloes:', n
   !    OPEN (7, file=infile)
   !    DO i = 1, hash_lines
   !       READ (7, *)
   !    END DO
   !    DO i = 1, n
   !       READ (7, *) (data(j), j=1, columns)
   !       mm = data(column_mv) ! Read virial mass
   !       pid = nint(data(column_pid))
   !       IF (mm > mmin .AND. pid == -1) p = p+1
   !    END DO
   !    CLOSE (7)
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Minimum allowed halo mass [Msun/h]:', mmin
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Total number of distinct haloes:', p
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Fraction of distinct haloes:', REAL(p)/REAL(n)

   !    ALLOCATE (x(3, p), m(6, p))
   !    x = 0.
   !    m = 0.

   !    ! Fill arrays with unique halo properties
   !    p = 0 ! Reset sum variable to zero
   !    OPEN (7, file=infile)
   !    DO i = 1, hash_lines
   !       READ (7, *)
   !    END DO
   !    DO i = 1, n
   !       ! Virial mass is column 11, x,y,z positions are columns 18,19,20
   !       ! Halo is unique iff pid=-1 (pid is column 6)
   !       READ (7, *) (data(j), j=1, columns)
   !       mm = data(column_mv)
   !       pid = nint(data(column_pid))
   !       IF (mm > mmin .AND. pid == -1) THEN
   !          p = p+1
   !          !WRITE(*,*) p
   !          x(1, p) = data(column_x)
   !          x(2, p) = data(column_y)
   !          x(3, p) = data(column_z)
   !          m(1, p) = data(column_mv)
   !          m(2, p) = data(column_mvu)
   !          m(3, p) = data(column_m200)
   !          m(4, p) = data(column_m200c)
   !          m(5, p) = data(column_m500c)
   !          m(6, p) = data(column_m2500c)
   !       END IF
   !    END DO
   !    CLOSE (7)

   !    ! Set the total output number
   !    n = p

   !    ! Write min/max values to screen
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Minimum x [Mpc/h]:', MINVAL(x(1, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Maximum x [Mpc/h]:', MAXVAL(x(1, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Minimum y [Mpc/h]:', MINVAL(x(2, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Maximum y [Mpc/h]:', MAXVAL(x(2, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Minimum z [Mpc/h]:', MINVAL(x(3, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Maximum z [Mpc/h]:', MAXVAL(x(3, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Minimum virial halo mass [Msun/h]:', MINVAL(m(1, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Maximum virial halo mass [Msun/h]:', MAXVAL(m(1, :))
   !    WRITE (*, *) 'READ_MULTIDARK_HALOES: Done'
   !    WRITE (*, *)

   ! END SUBROUTINE read_multidark_haloes

   SUBROUTINE read_multidark_halo_catalogue(infile, x, m, n)

      ! New version for halo catalogues downloaded via https://www.cosmosim.org/query
      CHARACTER(len=*), INTENT(IN) :: infile
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: m(:, :)
      INTEGER, INTENT(OUT) :: n
      INTEGER :: i, j, crap, nmass

      IF (snippet_in_string('rockstar', infile)) THEN
         nmass = 6
      ELSE
         nmass = 2
      END IF

      n = file_length(infile, verbose=.FALSE.)
      n = n-1 ! Remove comment line
      ALLOCATE (x(3, n), m(nmass, n))

      WRITE (*, *) 'READ_MULTIDARK_HALO_CATALOGUE: ', trim(infile)
      WRITE (*, *) 'READ_MULTIDARK_HALO_CATALOGUE: Number of haloes:', n
      OPEN (7, file=infile, status='old')
      READ (7, *) ! Comment line
      DO i = 1, n
         READ (7, *) crap, (x(j, i), j=1,3), (m(j, i), j=1,nmass)
      END DO
      CLOSE (7)
      WRITE (*, *) 'READ_MULTIDARK_HALO_CATALOGUE: Done'
      WRITE (*, *)

   END SUBROUTINE read_multidark_halo_catalogue

   SUBROUTINE write_multidark_halo_catalogue(outfile, x, m, idx)

      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL, INTENT(IN) :: x(:, :)
      REAL, INTENT(IN) :: m(:, :)
      INTEGER, INTENT(IN) :: idx(:)
      INTEGER :: i, j, k
      INTEGER :: nd, nm, n

      nd = size(x, 1)
      nm = size(m, 1)
      n = size(x, 2)

      ! Write out the little catalogue
      WRITE (*, *) 'WRITE_MULTIDARK_HALO_CATALOGUE: Writing outfile: ', trim(outfile)
      OPEN (7, file=outfile)
      DO i = 1, n
         j = idx(n+1-i)
         WRITE (7, *) (x(k, j), k=1, nd), (m(k, j), k=1,nm)
      END DO
      CLOSE (7)
      WRITE (*, *) 'WRITE_MULTIDARK_HALO_CATALOGUE: Done'
      WRITE (*, *)

   END SUBROUTINE write_multidark_halo_catalogue

   SUBROUTINE read_multidark_particles(infile, x, n)

      CHARACTER(len=*), INTENT(IN) :: infile
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :)
      INTEGER, INTENT(OUT) :: n
      INTEGER :: i
      REAL :: crap

      WRITE (*, *) 'READ_MULTIDARK_PARTICLES: ', trim(infile)

      n = file_length(infile, verbose=.FALSE.)
      n = n-1 ! First line is comment
      ALLOCATE (x(3, n))

      WRITE (*, *) 'READ_MULTIDARK_PARTICLES: Number of particles:', n
      OPEN (7, file=infile, status='old')
      READ (7, *) ! First line is comment
      DO i = 1, n
         READ (7, *) crap, x(1, i), x(2, i), x(3, i)
      END DO
      CLOSE (7)
      WRITE (*, *) 'READ_MULTIDARK_PARTICLES: Done'
      WRITE (*, *)

   END SUBROUTINE read_multidark_particles

   INTEGER FUNCTION multidark_snapshot(a)

      ! Exact multidark snapshot corresponding to a given scale factor
      REAL, INTENT(IN) :: a
      INTEGER :: i
      REAL, PARAMETER :: eps = eps_snap

      i = array_position(a, as, eps)
      IF (i == 0) THEN
         WRITE(*, *) 'MULTIDARK_SNAPSHOT: a:', a
         STOP 'MULTIDARK_SNAPSHOT: Error, no snapshot corresponding to the scale factor'
      ELSE
         multidark_snapshot = i
      END IF

   END FUNCTION multidark_snapshot

   REAL FUNCTION multidark_scale_factor(s)

      ! Multidark scale factor corresponding to a given snapshot
      ! Inverse function of multidark_snapshot
      INTEGER, INTENT(IN) :: s
      INTEGER ::  i

      IF (s == 416) THEN
         multidark_scale_factor = 1.000 ! This is actually Bolshoi, this is lazy
      ELSE
         i = array_position(s, snaps)
         IF (i == 0) THEN
            WRITE(*, *) 'MULTIDARK_SCALE_FACTOR: snapshot:', s
            STOP 'MULTIDARK_SNAPSHOT: Error, no scale factor corresponding to the snapshot'
         ELSE
            multidark_scale_factor = as(i)
         END IF
      ENDIF

   END FUNCTION multidark_scale_factor

   INTEGER FUNCTION nearest_multidark_snapshot(a)

      ! Nearest multidark snapshot to scale factor 'a' in terms of 'a' distance
      REAL, INTENT(IN) :: a
      INTEGER :: i
      INTEGER, PARAMETER :: ifind = ifind_snap

      i = nearest_table_integer(a, as, ifind)
      nearest_multidark_snapshot = snaps(i)

   END FUNCTION nearest_multidark_snapshot

   INTEGER FUNCTION nearest_multidark_snapshot_sig8(sig8)

      ! Nearest multidark snapshot to scale factor in terms of 'sigma_8' distance
      REAL, INTENT(IN) :: sig8
      INTEGER :: i
      INTEGER, PARAMETER :: ifind = ifind_snap

      i = nearest_table_integer(sig8, sig8s, ifind)
      nearest_multidark_snapshot_sig8 = snaps(i)

   END FUNCTION nearest_multidark_snapshot_sig8

END MODULE multidark_stuff
