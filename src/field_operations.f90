MODULE field_operations

   USE precision
   USE constants
   USE basic_operations
   USE array_operations
   USE fft

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: NGP_cell
   PUBLIC :: cell_position
   PUBLIC :: make_Gaussian_random_field
   PUBLIC :: generate_displacement_fields
   PUBLIC :: write_3D_field_projection_ascii
   PUBLIC :: compress_field
   PUBLIC :: add_to_stack_3D
   PUBLIC :: project_3D_to_2D
   PUBLIC :: clip
   PUBLIC :: anticlip
   PUBLIC :: periodic_distance
   PUBLIC :: check_Hermiticity

   PUBLIC :: box_mode_power
   PUBLIC :: compute_power_spectrum_pole ! Remove?
   !PUBLIC :: compute_power_spectrum_rsd
   !PUBLIC :: compute_power_spectrum_rsd2

   PUBLIC :: field_correlation_function

   PUBLIC :: read_field_binary
   PUBLIC :: write_field_binary
   PUBLIC :: write_field_ascii
   PUBLIC :: compute_power_spectrum
   PUBLIC :: compute_multipole_power_spectrum
   PUBLIC :: smooth
   PUBLIC :: sharpen
   PUBLIC :: sharpen_k
   PUBLIC :: empty_cells

   INTERFACE read_field_binary
      MODULE PROCEDURE read_field_binary_2D
      MODULE PROCEDURE read_field_binary_3D
   END INTERFACE read_field_binary

   INTERFACE write_field_binary
      MODULE PROCEDURE write_field_binary_2D
      MODULE PROCEDURE write_field_binary_3D
   END INTERFACE write_field_binary

   INTERFACE write_field_ascii
      MODULE PROCEDURE write_field_ascii_2D
   END INTERFACE write_field_ascii

   INTERFACE compute_power_spectrum
      MODULE PROCEDURE compute_power_spectrum_2D
      MODULE PROCEDURE compute_power_spectrum_3D
   END INTERFACE compute_power_spectrum

   INTERFACE compute_multipole_power_spectrum
      MODULE PROCEDURE compute_multipole_power_spectrum_3D
   END INTERFACE compute_multipole_power_spectrum

   INTERFACE smooth
      MODULE PROCEDURE smooth_2D
      MODULE PROCEDURE smooth_3D
   END INTERFACE smooth

   INTERFACE sharpen
      MODULE PROCEDURE sharpen_2D
      MODULE PROCEDURE sharpen_3D
   END INTERFACE sharpen

   INTERFACE sharpen_k
      MODULE PROCEDURE sharpen_k_2D
      MODULE PROCEDURE sharpen_k_3D
   END INTERFACE sharpen_k

   INTERFACE empty_cells
      MODULE PROCEDURE empty_cells_3D
   END INTERFACE empty_cells

   ! Full complex Fourier transforms, or assume real and therefore Hermitian arrays
   LOGICAL, PARAMETER :: complex = .TRUE.

   ! Gaussian random mode generation
   ! 1 - Assign mode amplitude and phase
   ! 2 - Assign mode real and imaginary part
   ! 3 - Mean variance to all modes
   INTEGER, PARAMETER :: imode_GRF = 1

CONTAINS

   INTEGER FUNCTION NGP_cell(x, L, m)

      ! Find the integer coordinates of the cell that coordinate x is in
      REAL, INTENT(IN) :: x    ! Particle position
      REAL, INTENT(IN) :: L    ! Box size (could be Mpc/h or angle or something else)
      INTEGER, INTENT(IN) :: m ! Number of mesh cells in grid

      IF (x == 0.) THEN
         ! Catch this edge case
         NGP_cell = 1
      ELSE
         NGP_cell = ceiling(x*real(m)/L)
      END IF

   END FUNCTION NGP_cell

   REAL FUNCTION cell_position(i, L, m)

      ! Gets the coordinates of cell centre i in box of length L with m cells
      INTEGER, INTENT(IN) :: i ! Integer label for cell
      REAL, INTENT(IN) :: L    ! Actual size corresponding to volume
      INTEGER, INTENT(IN) :: m ! Number of mesh cells

      cell_position = L*(i-0.5)/real(m)

   END FUNCTION cell_position

   ! REAL FUNCTION random_mode_amplitude(Pk, use_average)

   !    ! This calculates the Fourier amplitudes of the density field
   !    USE interpolate
   !    USE constants
   !    USE random_numbers
   !    REAL, INTENT(IN) :: Pk             ! NOTE: Not Delta^2(k)
   !    LOGICAL, INTENT(IN) :: use_average
   !    REAL :: sigma

   !    ! Sigma parameter in the Rayleigh distribution
   !    sigma = sqrt(Pk/2.) ! Half the variance to each of the real and imaginary parts

   !    IF (use_average) THEN
   !       ! Fixed mode amplitudes
   !       ! TODO: The mean of Rayleigh is sqrt(pi/2)*sigma, shouldn't this be here instead?
   !       random_mode_amplitude = sigma*sqrt(2.)
   !    ELSE
   !       ! Correctly assigned random mode amplitudes
   !       random_mode_amplitude = random_Rayleigh(sigma)
   !    END IF

   ! END FUNCTION random_mode_amplitude

   SUBROUTINE make_Gaussian_random_modes(dk, m, L, k_tab, Pk_tab)

      ! Uses a tablulated P(k) to make a Gaussian Random Field realisation
      ! Assumes that P(k) is real and corresponds to a real array, so the output modes are Hermitian
      ! TODO: Upgrade so that reduced arrays (mn, m, m) can be filled as well as full (m, m, m) arrays
      USE constants
      USE FFT
      USE random_numbers
      USE interpolate

      COMPLEX, ALLOCATABLE, INTENT(OUT) :: dk(:, :, :)
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: k_tab(:)
      REAL, INTENT(IN) :: Pk_tab(:)
      INTEGER :: ix, iy, iz
      REAL :: kx, ky, kz, k
      REAL :: amp, D2, Pk, sig, modes(2)
      COMPLEX :: rot
      TYPE(interpolator1D) :: Pk_interp
      INTEGER, PARAMETER :: iorder = 3
      INTEGER, PARAMETER :: iextrap = iextrap_lin
      LOGICAL, PARAMETER :: store = .TRUE.
      LOGICAL, PARAMETER :: logk = .TRUE.
      LOGICAL, PARAMETER :: logPk = .TRUE.

      ! Write useful things to screen
      WRITE (*, *) 'MAKE_GAUSSIAN_RANDOM_MODES: Creating Fourier realisation of Gaussian field'
      WRITE (*, *) 'MAKE_GAUSSIAN_RANDOM_MODES: Mesh size:', m
      WRITE (*, *) 'MAKE_GAUSSIAN_RANDOM_MODES: Box size [Mpc/h]:', L

      ! Create an interpolator for P(k)
      CALL init_interpolator(k_tab, Pk_tab, Pk_interp, iorder, iextrap, store, logk, logPk)

      ! Allocate the Fourier array
      ALLOCATE(dk(m, m, m))

      ! This fills up all of k space, no Hermiticity imposed
      DO iz = 1, m
         DO iy = 1, m
            DO ix = 1, m

               ! Wavenumber
               CALL k_fft(ix, iy, iz, m, kx, ky, kz, k, L)

               IF (ix == 1 .AND. iy == 1 .AND. iz == 1) THEN

                  ! Set the zero mode to zero
                  dk(ix, iy, iz) = (0.d0, 0.d0)

               ELSE

                  ! Get mode standard error
                  D2 = evaluate_interpolator(k, Pk_interp)
                  Pk = D2/(4.*pi*(k*L/twopi)**3)
                  sig = sqrt(Pk/2.) ! Half of the variance to each of the real and imaginary part

                  IF (imode_GRF == 1) THEN
                     amp = random_Rayleigh(sig)
                     rot = random_complex_unit()
                     dk(ix, iy, iz) = amp*rot
                  ELSE IF (imode_GRF == 2)  THEN                 
                     modes = random_Gaussian_pair(0., sig) 
                     dk(ix, iy, iz) = cmplx(modes(1), modes(2))
                  ELSE IF (imode_GRF == 3) THEN
                     !amp = sig*sqrt(2.)
                     amp = sqrt(pi/2.)*sig
                     rot = random_complex_unit()
                     dk(ix, iy, iz) = amp*rot
                  ELSE
                     STOP 'MAKE_GAUSSIAN_RANDOM_MODES: Error, mode creation scheme not specified'
                  END IF

               END IF

            END DO
         END DO
      END DO

      WRITE (*, *) 'MAKE_GAUSSIAN_RANDOM_MODES: Initial modes generated'

      ! Make the array Hermitian
      CALL enforce_Hermiticity(dk, verbose=.FALSE.)  

      WRITE (*, *) 'MAKE_GAUSSIAN_RANDOM_MODES: Hermitian field generated'
      WRITE (*, *)

   END SUBROUTINE make_Gaussian_random_modes

   SUBROUTINE enforce_Hermiticity(dk, verbose)

      ! Take an input complex array and force it in to a Hermitian state
      USE FFT
      COMPLEX, INTENT(INOUT) :: dk(:, :, :)
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: ix, iy, iz
      INTEGER :: jx, jy, jz
      INTEGER :: mx, my, mz

      IF (verbose) WRITE (*, *) 'ENFORCE_HERMITICITY: Enforcing Hermiticity'

      mx = size(dk, 1)
      my = size(dk, 2)
      mz = size(dk, 3)

      ! Enforce Hermiticity - probably could save a load of operations above
      DO iz = 1, mz
         jz = conjugate_mode_location(iz, mz)
         DO iy = 1, my
            jy = conjugate_mode_location(iy, my)
            DO ix = 1, mx
               jx = conjugate_mode_location(ix, mx)
               IF ((ix == jx) .AND. (iy == jy) .AND. (iz == jz)) THEN
                  dk(ix, iy, iz) = real(dk(ix, iy, iz)) ! Nyquist modes should be real
               ELSE
                  dk(ix, iy, iz) = conjg(dk(jx, jy, jz)) ! Other modes should be conjugate
               END IF
            END DO
         END DO
      END DO

      ! Because all modes have been conjugated above
      dk = conjg(dk)

      IF (verbose) THEN
         WRITE (*, *) 'ENFORCE_HERMITICITY: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE enforce_Hermiticity

   LOGICAL FUNCTION check_Hermiticity(dk)

      ! TODO: Should there be a tolerance parameter here?
      USE FFT
      COMPLEX, INTENT(IN) :: dk(:, :, :)
      INTEGER :: ix, iy, iz
      INTEGER :: mx, my, mz
      INTEGER :: jx, jy, jz

      mx = size(dk, 1)
      my = size(dk, 2)
      mz = size(dk, 3)

      check_Hermiticity = .TRUE.
      DO iz = 1, mz
         jz = conjugate_mode_location(iz, mz)
         DO iy = 1, my
            jy = conjugate_mode_location(iy, my)
            DO ix = 1, mx
               jx = conjugate_mode_location(ix, mx)
               IF ((ix == jx) .AND. (iy == jy) .AND. (iz == jz)) THEN
                  IF (aimag(dk(ix, iy, iz)) .NE. 0.) THEN
                     check_Hermiticity = .FALSE. ! Nyquist modes should be real
                     RETURN
                  END IF
               ELSE
                  IF (dk(ix, iy, iz) .NE. conjg(dk(jx, jy, jz))) THEN
                     check_Hermiticity = .FALSE. ! Other modes should be conjugate
                     RETURN
                  END IF
               END IF
            END DO
         END DO
      END DO

   END FUNCTION check_Hermiticity

   SUBROUTINE apodise_field(dk, L, ka)

      ! Remove Fourier modes that exist beyond ka
      ! These modes are not fully sampled in k space because the shell is cut off
      ! TODO: Different choices for apodising function
      USE constants
      USE FFT

      COMPLEX, INTENT(INOUT) :: dk(:, :, :)
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: ka
      REAL :: kx, ky, kz, k
      INTEGER :: ix, iy, iz
      INTEGER :: m

      m = size(dk, 1)
      IF ((m .NE. size(dk, 2)) .OR. (m .NE. size(dk, 3))) STOP

      DO iz = 1, m
         DO iy = 1, m
            DO ix = 1, m
               CALL k_FFT(ix, iy, iz, m, kx, ky, kz, k, L)
               IF (k > ka) dk(ix, iy, iz) = (0., 0.)
            END DO
         END DO
      END DO

   END SUBROUTINE apodise_field

   SUBROUTINE make_Gaussian_random_field(d, m, L, k_tab, Pk_tab)

      ! Uses a tablulated P(k) to make a Gaussian Random Field realisation
      ! TODO: Should modes be removed beyond the Nyquist frequency?

      REAL, ALLOCATABLE, INTENT(OUT) :: d(:, :, :)
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: k_tab(:)
      REAL, INTENT(IN) :: Pk_tab(:)
      COMPLEX, ALLOCATABLE :: dk(:, :, :)
      COMPLEX :: dk_new(m, m, m)

      CALL make_Gaussian_random_modes(dk, m, L, k_tab, Pk_tab)

      WRITE (*, *) 'MAKE_GAUSSIAN_RANDOM_FIELD: Transform to real space'

      ! FT the displacement field from k-space to real space!
      CALL fft3(dk, dk_new, m, m, m, 1)
      dk = dk_new
      d = real(dk)

      WRITE (*, *) 'MAKE_GAUSSIAN_RANDOM_FIELD: Done'
      WRITE (*, *)

   END SUBROUTINE make_Gaussian_random_field

   SUBROUTINE generate_displacement_fields(f, m, L, k_tab, Pk_tab)

      REAL, ALLOCATABLE, INTENT(OUT) :: f(:, :, :, :)
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: k_tab(:)
      REAL, INTENT(IN) :: Pk_tab(:)
      COMPLEX :: d(m, m, m), fk(3, m, m, m)
      COMPLEX, ALLOCATABLE :: dk(:, :, :)
      INTEGER :: i, ix, iy, iz
      REAL :: kx, ky, kz, k, kNy

      ! Make a Gaussian realisation of overdensity modes
      CALL make_Gaussian_random_modes(dk, m, L, k_tab, Pk_tab)

      ! Remove modes above the Nyquist frequency
      ! TODO: Think about this
      kNy = Nyquist_frequency(L, m)
      CALL apodise_field(dk, L, kNy)

      WRITE (*, *) 'GENERATE_DISPLACEMENT_FIELDS: Creating realisation of displacement field'

      ! Allocate displacement field array
      ALLOCATE(f(3, m, m, m))

      ! This fills up displacement array in all of k space!
      DO iz = 1, m
         DO iy = 1, m
            DO ix = 1, m

               ! Get the wave vectors
               CALL k_fft(ix, iy, iz, m, kx, ky, kz, k, L)

               IF (ix == 1 .AND. iy == 1 .AND. iz == 1) THEN

                  ! Fix the zero mode to zero, avoiding a division by zero
                  fk(1, ix, iy, iz) = (0.d0, 0.d0)
                  fk(2, ix, iy, iz) = (0.d0, 0.d0)
                  fk(3, ix, iy, iz) = (0.d0, 0.d0)

               ! ELSE IF (ix == 1+m/2 .OR. iy == 1+m/2 .OR. iz == 1+m/2) THEN

               !    ! Sets Nyquist modes to 0.!
               !    ! Maybe all modes with mod(k)>k_ny should be set to 0.?!
               !    ! Bridit Falck wrote a 'corner modes' paper about this
               !    ! https://arxiv.org/abs/1610.04862
               !    fk(1, ix, iy, iz) = (0.d0, 0.d0)
               !    fk(2, ix, iy, iz) = (0.d0, 0.d0)
               !    fk(3, ix, iy, iz) = (0.d0, 0.d0)

               ELSE

                  ! Assign values to the displacement field
                  fk(1, ix, iy, iz) = (0.d0, -1.d0)*dk(ix, iy, iz)*kx/k**2
                  fk(2, ix, iy, iz) = (0.d0, -1.d0)*dk(ix, iy, iz)*ky/k**2
                  fk(3, ix, iy, iz) = (0.d0, -1.d0)*dk(ix, iy, iz)*kz/k**2

               END IF

            END DO
         END DO
      END DO
      WRITE (*, *) 'GENERATE_DISPLACEMENT_FIELDS: Fourier displacement fields done'

      ! Reverse Fourier transform back to configuration space
      DO i = 1, 3
         dk = fk(i, :, :, :)
         CALL fft3(dk, d, m, m, m, 1)
         f(i, :, :, :) = real(d(:, :, :))
      END DO

      WRITE (*, *) 'GENERATE_DISPLACEMENT_FIELDS: Minimum 1D displacement [Mpc/h]:', minval(f)
      WRITE (*, *) 'GENERATE_DISPLACEMENT_FIELDS: Maximum 1D displacement [Mpc/h]:', maxval(f)
      WRITE (*, *) 'GENERATE_DISPLACEMENT_FIELDS: Real-space displacement fields generated'
      WRITE (*, *)

   END SUBROUTINE generate_displacement_fields

   SUBROUTINE read_field_binary_2D(d, m, L, infile)

      ! Read in a binary 'field' file
      USE statistics

      REAL, ALLOCATABLE, INTENT(OUT) :: d(:, :)
      INTEGER, INTENT(OUT) :: m
      REAL, INTENT(OUT) :: L
      CHARACTER(len=*), INTENT(IN) :: infile
      LOGICAL :: lexist

      ! Write file name to screen and check it exists
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Binary input: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) STOP 'READ_SIMULATION_POWER_SPECTRUM: File does not exist'

      ! First get mesh size and physical size
      OPEN (7, file=infile, form='unformatted', access='stream', status='old')
      READ (7) m
      READ (7) L
      CLOSE (7)
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Mesh size:', m
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Field size:', L

      ALLOCATE (d(m, m))

      ! Read unformatted data
      OPEN (7, file=infile, form='unformatted', access='stream', status='old')
      READ (7) m
      READ (7) L
      READ (7) d
      CLOSE (7)
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Minimum field value:', minval(d)
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Maximum field value:', maxval(d)
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Average:', mean(splay(d, m, m))
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Variance:', variance(splay(d, m, m))
      WRITE (*, *) 'READ_FIELD_BINARY_2D: Done'
      WRITE (*, *)

   END SUBROUTINE read_field_binary_2D

   SUBROUTINE read_field_binary_3D(d, m, L, infile)

      ! Read in a binary 'field' file
      USE statistics

      REAL, ALLOCATABLE, INTENT(OUT) :: d(:, :, :)
      INTEGER, INTENT(OUT) :: m
      REAL, INTENT(OUT) :: L
      CHARACTER(len=*), INTENT(IN) :: infile
      LOGICAL :: lexist

      ! Write file name to screen and check it exists
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Binary input: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) STOP 'READ_SIMULATION_POWER_SPECTRUM: File does not exist'

      ! First get mesh size and physical size
      OPEN (7, file=infile, form='unformatted', access='stream', status='old')
      READ (7) m
      READ (7) L
      CLOSE (7)
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Mesh size:', m
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Field size:', L

      ALLOCATE (d(m, m, m))

      ! Read unformatted data
      OPEN (7, file=infile, form='unformatted', access='stream', status='old')
      READ (7) m
      READ (7) L
      READ (7) d
      CLOSE (7)
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Minval:', minval(d)
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Maxval:', maxval(d)
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Average:', mean(splay(d, m, m, m))
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Variance:', variance(splay(d, m, m, m))
      WRITE (*, *) 'READ_FIELD_BINARY_3D: Done'
      WRITE (*, *)

   END SUBROUTINE read_field_binary_3D

   SUBROUTINE write_field_binary_2D(d, m, L, outfile)

      ! Write out a binary 'field' file
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d(m, m)
      REAL, INTENT(IN) :: L
      CHARACTER(len=*), INTENT(IN) :: outfile

      WRITE (*, *) 'WRITE_FIELD_BINARY_2D: Binary output: ', trim(outfile)
      WRITE (*, *) 'WRITE_FIELD_BINARY_2D: Mesh size:', m
      WRITE (*, *) 'WRITE_FIELD_BINARY_2D: Field size:', L
      WRITE (*, *) 'WRITE_FIELD_BINARY_2D: Minimum field value:', minval(d)
      WRITE (*, *) 'WRITE_FIELD_BINARY_2D: Maximum field value:', maxval(d)
      WRITE (*, *) 'WRITE_FIELD_BINARY_2D: Using new version with access=stream'
      OPEN (7, file=outfile, form='unformatted', access='stream', status='replace')
      WRITE (7) m
      WRITE (7) L
      WRITE (7) d
      CLOSE (7)
      WRITE (*, *) 'WRITE_FIELD_BINARY_2D: Done'
      WRITE (*, *)

   END SUBROUTINE write_field_binary_2D

   SUBROUTINE write_field_binary_3D(d, m, L, outfile)

      ! Write out a binary 'field' file
      ! Used to be called write_field
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d(m, m, m) 
      REAL, INTENT(IN) :: L
      CHARACTER(len=*), INTENT(IN) :: outfile

      WRITE (*, *) 'WRITE_FIELD_BINARY_3D: Binary output: ', trim(outfile)
      WRITE (*, *) 'WRITE_FIELD_BINARY_3D: Mesh size:', m
      WRITE (*, *) 'WRITE_FIELD_BINARY_3D: Field size:', L
      WRITE (*, *) 'WRITE_FIELD_BINARY_3D: Minimum field value:', minval(d)
      WRITE (*, *) 'WRITE_FIELD_BINARY_3D: Maximum field value:', maxval(d)
      WRITE (*, *) 'WRITE_FIELD_BINARY_3D: Using new version with access=stream'
      OPEN (7, file=outfile, form='unformatted', access='stream', status='replace')
      WRITE (7) m
      WRITE (7) L
      WRITE (7) d
      CLOSE (7)
      WRITE (*, *) 'WRITE_FIELD_BINARY_3D: Done'
      WRITE (*, *)

   END SUBROUTINE write_field_binary_3D

   SUBROUTINE write_field_ascii_2D(d, m, L, outfile)

      ! Used to be called print_2D_field
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d(m, m), L   
      CHARACTER(len=*), INTENT(IN) :: outfile
      INTEGER :: i, j
      REAL :: x, y

      WRITE (*, *) 'WRITE_FIELD_ASCII_2D: Writing to: ', trim(outfile)
      WRITE (*, *) 'WRITE_FIELD_ASCII_2D: Field mesh size:', m
      WRITE (*, *) 'WRITE_FIELD_ASCII_2D: Field physical size:', L
      WRITE (*, *) 'WRITE_FIELD_ASCII_2D: Minimum field value:', minval(d)
      WRITE (*, *) 'WRITE_FIELD_ASCII_2D: Maximum field value:', maxval(d)

      OPEN (8, file=outfile)
      DO j = 1, m
         DO i = 1, m

            x = cell_position(i, L, m)
            y = cell_position(j, L, m)

            WRITE (8, *) x, y, d(i, j)

         END DO
      END DO
      CLOSE (8)

      WRITE (*, *) 'WRITE_FIELD_ASCII_2D: Done'
      WRITE (*, *)

   END SUBROUTINE write_field_ascii_2D

   SUBROUTINE write_3D_field_projection_ascii(d, m, L, nz, outfile)

      ! Used to be called print_projected_field
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d(m, m, m)
      REAL, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: nz
      CHARACTER(len=*), INTENT(IN) :: outfile
      INTEGER :: i, j, k
      REAL :: x, y
      REAL :: sum

      WRITE (*, *) 'WRITE_3D_FIELD_PROJECTION_ASCII: Writing to: ', trim(outfile)
      WRITE (*, *) 'WRITE_3D_FIELD_PROJECTION_ASCII: Cells projecting:', nz

      OPEN (8, file=outfile)
      DO j = 1, m
         DO i = 1, m

            x = cell_position(i, L, m)
            y = cell_position(j, L, m)

            sum = 0.
            DO k = 1, nz
               sum = sum+d(i, j, k)
            END DO
            sum = sum/real(nz)

            WRITE (8, *) x, y, sum

         END DO
      END DO
      CLOSE (8)

      WRITE (*, *) 'WRITE_3D_FIELD_PROJECTION_ASCII: Done'
      WRITE (*, *)

   END SUBROUTINE write_3D_field_projection_ascii

   SUBROUTINE compress_field(d, ds, m)

      ! Shrinks a 3D field size by a factor of 2
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d(m, m, m)
      REAL, ALLOCATABLE, INTENT(OUT) :: ds(:, :, :)
      INTEGER :: i, j, k

      ! Allocate the small array
      ALLOCATE (ds(m/2, m/2, m/2))
      ds = 0.

      ! Fill up the small array by summing blocks of 8 from the larger array
      DO k = 1, m/2
         DO j = 1, m/2
            DO i = 1, m/2
               ds(i, j, k) = ds(i, j, k)+d(2*i-1, 2*j-1, 2*k-1)
               ds(i, j, k) = ds(i, j, k)+d(2*i, 2*j-1, 2*k-1)
               ds(i, j, k) = ds(i, j, k)+d(2*i-1, 2*j, 2*k-1)
               ds(i, j, k) = ds(i, j, k)+d(2*i-1, 2*j-1, 2*k)
               ds(i, j, k) = ds(i, j, k)+d(2*i, 2*j, 2*k-1)
               ds(i, j, k) = ds(i, j, k)+d(2*i, 2*j-1, 2*k)
               ds(i, j, k) = ds(i, j, k)+d(2*i-1, 2*j, 2*k)
               ds(i, j, k) = ds(i, j, k)+d(2*i, 2*j, 2*k)
            END DO
         END DO
      END DO

      ! Divide by the number of blocks that are being averaged over
      ds = ds/8.

   END SUBROUTINE compress_field

   SUBROUTINE sharpen_2D(d, m, ibin)

      ! Sharpen a 3D configuration-space array to account for the binning
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(INOUT) :: d(m, m)  
      INTEGER, INTENT(IN) :: ibin
      COMPLEX, ALLOCATABLE :: dk(:, :)
      COMPLEX :: dkout(m, m)
      REAL :: dc(m, m)
      INTEGER :: mn

      WRITE (*, *) 'SHARPEN_2D: Correcting for binning by sharpening field'
      WRITE (*, *) 'SHARPEN_2D: Mesh size:', m

      IF (complex) THEN
         mn = m
         WRITE (*, *) 'SHARPEN_2D: Doing complex->complex FFT'
      ELSE
         mn = m/2+1
         WRITE (*, *) 'SHARPEN_2D: Doing real->complex FFT'
      END IF
      WRITE (*, *) 'SHARPEN_2D: Mesh size general:', m
      WRITE (*, *) 'SHARPEN_2D: Mesh size first-dimension:', mn

      ALLOCATE (dk(mn, m))

      IF (complex) THEN
         dk = d
         CALL fft2(dk, dkout, m, m, -1)
         dk = dkout
      ELSE
         dc = d
         CALL fft2(dc, dk, m, m, -1)
      END IF

      CALL sharpen_k_2D(dk, ibin)

      IF (complex) THEN
         CALL fft2(dk, dkout, m, m, 1)
         dk = dkout
         d = real(real(dk))/real(m**2)
      ELSE
         CALL fft2(dc, dk, m, m, 1)
         d = real(dc)
      END IF

      WRITE (*, *) 'SHARPEN_2D: Sharpening complete'
      WRITE (*, *)

   END SUBROUTINE sharpen_2D

   SUBROUTINE sharpen_3D(d, m, ibin)

      ! Sharpen a 3D configuration-space array to account for the binning
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(INOUT) :: d(m, m, m)
      INTEGER, INTENT(IN) :: ibin
      COMPLEX, ALLOCATABLE :: dk(:, :, :)
      COMPLEX :: dkout(m, m, m)
      REAL :: dc(m, m, m)
      INTEGER :: mn

      WRITE (*, *) 'SHARPEN_3D: Correcting for binning by sharpening field'
      WRITE (*, *) 'SHARPEN_3D: Mesh size:', m

      IF (complex) THEN
         mn = m
         WRITE (*, *) 'SHARPEN_3D: Doing complex->complex FFT'
      ELSE
         mn = m/2+1
         WRITE (*, *) 'SHARPEN_3D: Doing real->complex FFT'
      END IF
      WRITE (*, *) 'SHARPEN_3D: Mesh size general:', m
      WRITE (*, *) 'SHARPEN_3D: Mesh size first-dimension:', mn

      ALLOCATE (dk(mn, m, m))

      IF (complex) THEN
         dk = d
         CALL fft3(dk, dkout, m, m, m, -1)
         dk = dkout
      ELSE
         dc = d
         CALL fft3(dc, dk, m, m, m, -1)
      END IF

      CALL sharpen_k_3D(dk, ibin)

      IF (complex) THEN
         CALL fft3(dk, dkout, m, m, m, 1)
         dk = dkout
         d = real(real(dk))/real(m**3)
      ELSE
         CALL fft3(dc, dk, m, m, m, 1)
         d = real(dc)
      END IF

      WRITE (*, *) 'SHARPEN_3D: Sharpening complete'
      WRITE (*, *)

   END SUBROUTINE sharpen_3D

   SUBROUTINE sharpen_k_2D(dk, ibin)

      USE special_functions

      ! Sharpens a 3D Fourier array to account for the binning
      COMPLEX, INTENT(INOUT) :: dk(:, :)
      INTEGER, INTENT(IN) :: ibin
      INTEGER :: i, j, mn, m
      REAL :: kx, ky, kmod
      REAL :: kxh, kyh
      REAL :: fcx, fcy, fcorr
      REAL :: crap
      REAL, PARAMETER :: L = 1. ! This does not matter for this routine

      mn = size(dk, 1)
      m = size(dk, 2)

      ! Check that the array is sensible
      IF (mn == m .OR. mn == m/2+1) THEN
         ! Do nothing
      ELSE
         WRITE (*, *) 'SHARPEN_K_2D: Array first-dimension size:', mn
         WRITE (*, *) 'SHARPEN_K_2D: Array general size:', m
         STOP 'SHARPEN_K_2D: Error, the array is not sensible'
      END IF

      ! Now correct for binning
      DO j = 1, m
         DO i = 1, mn

            CALL k_fft(i, j, 1, m, kx, ky, crap, kmod, L)

            kxh = L*kx/(2.*real(m))
            kyh = L*ky/(2.*real(m))

            fcx = sinc(kxh)
            fcy = sinc(kyh)

            IF (ibin == 1) THEN
               fcorr = fcx*fcy
            ELSE IF (ibin == 2) THEN
               fcorr = (fcx*fcy)**2
            ELSE
               STOP 'SHARPEN_K_2D: Error, ibin specified incorrectly'
            END IF

            dk(i, j) = dk(i, j)/fcorr

         END DO
      END DO

   END SUBROUTINE sharpen_k_2D

   SUBROUTINE sharpen_k_3D(dk, ibin)

      USE special_functions

      ! Sharpens a 3D array to account for the binning
      COMPLEX, INTENT(INOUT) :: dk(:, :, :)    
      INTEGER, INTENT(IN) :: ibin
      INTEGER :: i, j, k, m, mn
      REAL :: kx, ky, kz, kmod
      REAL :: kxh, kyh, kzh
      REAL :: fcx, fcy, fcz, fcorr
      REAL, PARAMETER :: L = 1. ! This does not matter for this routine

      mn = size(dk, 1)
      m = size(dk, 2)

      ! Check that the array is sensible
      IF (mn == m .OR. mn == m/2+1) THEN
         ! Do nothing; array is sensible
      ELSE
         WRITE (*, *) 'SHARPEN_K_3D: Array first-dimension size:', mn
         WRITE (*, *) 'SHARPEN_K_3D: Array general size:', m
         STOP 'SHARPEN_K_3D: Error, the array is not sensible'
      END IF

      ! Now correct for binning
      DO k = 1, m
         DO j = 1, m
            DO i = 1, mn

               CALL k_fft(i, j, k, m, kx, ky, kz, kmod, L)

               kxh = L*kx/(2.*real(m))
               kyh = L*ky/(2.*real(m))
               kzh = L*kz/(2.*real(m))

               fcx = sinc(kxh)
               fcy = sinc(kyh)
               fcz = sinc(kzh)

               IF (ibin == 1) THEN
                  fcorr = fcx*fcy*fcz
               ELSE IF (ibin == 2) THEN
                  fcorr = (fcx*fcy*fcz)**2
               ELSE
                  STOP 'SHARPEN_K_3D: Error, ibin specified incorrectly'
               END IF

               dk(i, j, k) = dk(i, j, k)/fcorr

            END DO
         END DO
      END DO

   END SUBROUTINE sharpen_k_3D

   SUBROUTINE smooth_2D(d, m, r, L)

      !arr(n,n): input array of size n x n
      !r: smoothing scale in Mpc/h
      !L: box size in Mpc/h
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(INOUT) :: d(m, m)
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: L
      REAL :: kx, ky, kz, k
      COMPLEX, ALLOCATABLE :: dk(:, :)
      COMPLEX :: dkout(m, m)
      REAL :: dc(m, m)
      INTEGER :: i, j, mn

      WRITE (*, *) 'SMOOTH_2D: Smoothing array'
      WRITE (*, *) 'SMOOTH_2D: Assuming array is periodic'
      WRITE (*, *) 'SMOOTH_2D: Smoothing scale [Mpc/h]:', r

      IF (complex) THEN
         mn = m
         WRITE (*, *) 'SMOOTH_2D: Doing complex->complex FFT'
      ELSE
         mn = m/2+1
         WRITE (*, *) 'SMOOTH_2D: Doing real->complex FFT'
      END IF
      WRITE (*, *) 'SMOOTH_2D: Mesh size:', m
      WRITE (*, *) 'SMOOTH_2D: Mesh size x:', mn
      ALLOCATE (dk(mn, m))

      ! Fourier transform
      IF (complex) THEN
         dk = d
         CALL fft2(dk, dkout, m, m, -1)
         dk = dkout
      ELSE
         dc = d
         CALL fft2(dc, dk, m, m, -1)
      END IF

      DO j = 1, m
         DO i = 1, mn
            CALL k_fft(i, j, 1, m, kx, ky, kz, k, L)
            dk(i, j) = dk(i, j)*exp(-((k*r)**2)/2.)
         END DO
      END DO

      ! Normalise post Fourier transform
      IF (complex) THEN
         CALL fft2(dk, dkout, m, m, 1)
         dk = dkout
         dk = dk/real(m**2)
         d = real(real(dk))
      ELSE
         CALL fft2(dc, dk, m, m, 1)
         dc = dc/real(m**2)
         d = real(dc)
      END IF
      DEALLOCATE (dk)

      WRITE (*, *) 'SMOOTH_2D: Done'
      WRITE (*, *)

   END SUBROUTINE smooth_2D

!!$  SUBROUTINE smooth2D_nonperiodic(arr,n,r,L)
!!$
!!$    !arr(n,n): input array of size n x n
!!$    !r: smoothing scale in Mpc/h
!!$    !L: box size in Mpc/h
!!$    USE fft
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(INOUT) :: arr(n,n)
!!$    REAL, INTENT(IN) :: r, L
!!$    REAL :: kx, ky, kz, k
!!$    DOUBLE COMPLEX, ALLOCATABLE :: ac(:,:), ac_out(:,:)
!!$    INTEGER :: i, j, m
!!$
!!$    INTEGER, PARAMETER :: pad=2 !Padding because we generally will not be continuous
!!$
!!$    WRITE(*,*) 'SMOOTH2D: Smoothing array'
!!$    WRITE(*,*) 'SMOOTH2D: Smoothing scale [Mpc/h]:', r
!!$
!!$    !For padding, I cant imagine that x2 would ever be insufficient!
!!$    m=pad*n
!!$
!!$    ALLOCATE(ac(m,m),ac_out(m,m))
!!$
!!$    !Not sure if this is necessary
!!$    ac=(0.d0,0.d0)
!!$    ac_out=(0.d0,0.d0)
!!$
!!$    !Put image into complex array, padded with zeroes where image is not!
!!$    DO j=1,n
!!$       DO i=1,n
!!$          ac(i,j)=arr(i,j)
!!$       END DO
!!$    END DO
!!$
!!$    CALL fft2(ac,ac_out,m,m,-1)
!!$    ac=ac_out
!!$
!!$    !Smoothing length in terms of image(m x m) size!
!!$    !r=pix/float(m)
!!$    DO j=1,m
!!$       DO i=1,m
!!$          CALL k_fft(i,j,1,m,kx,ky,kz,k,pad*L)
!!$          ac(i,j)=ac(i,j)*exp(-((k*r)**2)/2.)
!!$       END DO
!!$    END DO
!!$
!!$    CALL fft2(ac,ac_out,m,m,1)
!!$    ac=ac_out
!!$    DEALLOCATE(ac_out)
!!$
!!$    !Normalise post Fourier transform!
!!$    ac=ac/real(m**2)
!!$
!!$    !Retrieve smooth image from complex array!
!!$    !Need a loop because arr and ac will have different sizes
!!$    DO j=1,n
!!$       DO i=1,n
!!$          arr(i,j)=real(real(ac(i,j)))
!!$       END DO
!!$    END DO
!!$
!!$    WRITE(*,*) 'SMOOTH2D: Done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE smooth2D_nonperiodic

   SUBROUTINE smooth_3D(d, m, r, L)

      INTEGER, INTENT(IN) :: m
      REAL, INTENT(INOUT) :: d(m, m, m)
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: L
      REAL :: kx, ky, kz, kmod
      COMPLEX, ALLOCATABLE :: dk(:, :, :)
      COMPLEX :: dkout(m, m, m)
      REAL :: dc(m/2+1, m, m)
      INTEGER :: i, j, k, mn

      WRITE (*, *) 'SMOOTH_3D: Smoothing array'

      IF (complex) THEN
         mn = m
         WRITE (*, *) 'SMOOTH_3D: Doing complex->complex FFT'
      ELSE
         mn = m/2+1
         WRITE (*, *) 'SMOOTH_3D: Doing real->complex FFT'
      END IF
      WRITE (*, *) 'SMOOTH_3D: Mesh size:', m
      WRITE (*, *) 'SMOOTH_3D: Mesh size x:', mn
      ALLOCATE (dk(mn, m, m))

      ! Move to Fourier space
      IF (complex) THEN
         dk = d
         CALL fft3(dk, dkout, m, m, m, -1)
         dk = dkout
      ELSE
         dc = d
         CALL fft3(dc, dk, m, m, m, -1)
      END IF

      DO k = 1, m
         DO j = 1, m
            DO i = 1, mn
               CALL k_fft(i, j, k, m, kx, ky, kz, kmod, L)
               !dk(i,j,k)=dk(i,j,k)*sinc(kx*r/2.)*sinc(ky*r/2.)*sinc(kz*r/2.) ! Surely this should be a Gaussian?
               dk(i, j, k) = dk(i, j, k)*exp(-((kmod*r)**2)/2.)
            END DO
         END DO
      END DO

      ! Move back to real space and normalise
      IF (complex) THEN
         CALL fft3(dk, dkout, m, m, m, 1)
         dk = dkout
         dk = dk/(real(m)**3)
         d = real(real(dk))
      ELSE
         CALL fft3(dc, dk, m, m, m, 1)
         dc = dc/real(m)**3
         d = real(dc)
      END IF
      DEALLOCATE (dk)

      WRITE (*, *) 'SMOOTH_3D: Done'
      WRITE (*, *)

   END SUBROUTINE smooth_3D

   SUBROUTINE add_to_stack_3D(x, stack, Ls, ms, back, Lb, mb)

      ! Adds some points in a density field to a stack
      ! Assumes the background field is periodic
      ! 'stack' should have been previously allocated
      ! 'stack' should be set to zero before using this subroutine
      ! '*s' variables refer to the stacked field
      ! '*_back' variables refer to the background field
      INTEGER :: i, j, k, is(3), ib(3), d
      INTEGER, INTENT(IN) :: ms, mb
      REAL, INTENT(IN) :: x(3), Ls, Lb
      REAL, INTENT(INOUT) :: stack(ms, ms, ms)
      REAL, INTENT(IN) :: back(mb, mb, mb)
      REAL :: xb(3)

      ! Loop over cells on stacking mesh
      DO i = 1, ms
         DO j = 1, ms
            DO k = 1, ms

               ! Set the stack integer array
               is(1) = i
               is(2) = j
               is(3) = k

               DO d = 1, 3

                  ! Get coordinates of position on the stack
                  ! This changes coordiantes from stack to simulation coordinates
                  xb(d) = x(d)+Ls*(0.5+real(is(d)-1))/real(ms)-Ls/2.

                  ! Bring the coordinates back into the simulation box if they are outside
                  IF (xb(d) <= 0.) THEN
                     xb(d) = xb(d)+Lb
                  ELSE IF (xb(d) > Lb) THEN
                     xb(d) = xb(d)-Lb
                  END IF

                  ! Find the integer coordinates of mesh cell in the background mesh
                  ! This is just an NGP-type scheme. Could/should be improved?
                  ib(d) = ceiling(real(mb)*xb(d)/Lb)

               END DO

               ! Add the value to the stack
               ! Should there be a volume factor here?
               stack(is(1), is(2), is(3)) = stack(is(1), is(2), is(3))+back(ib(1), ib(2), ib(3))

            END DO
         END DO
      END DO

   END SUBROUTINE add_to_stack_3D

   SUBROUTINE project_3D_to_2D(d3d, d2d, m)

      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d3d(m, m, m) 
      REAL, INTENT(OUT) :: d2d(m, m)
      INTEGER :: i, j, k

      WRITE (*, *) 'PROJECT_3D_TO_2D: Projecting 3D stack into 2D'
      d2d = 0.
      DO i = 1, m
         DO j = 1, m
            DO k = 1, m
               d2d(i, j) = d2d(i, j)+d3d(i, j, k)
            END DO
         END DO
      END DO
      d2d = d2d/real(m)
      WRITE (*, *) 'PROJECT_3D_TO_2D: Minimum value of 2D stack:', minval(d2d)
      WRITE (*, *) 'PROJECT_3D_TO_2D: Maximum value of 2D stack:', maxval(d2d)
      WRITE (*, *) 'PROJECT_3D_TO_2D: Done'
      WRITE (*, *)

   END SUBROUTINE project_3D_to_2D

   SUBROUTINE clip(d, m1, m2, m3, d0, verbose)

      USE statistics

      REAL, INTENT(INOUT) :: d(:, :, :)
      REAL, INTENT(IN) :: d0
      INTEGER, INTENT(IN) :: m1, m2, m3
      LOGICAL, INTENT(IN) :: verbose
      REAL :: var1, av1, max1, var2, av2, max2
      INTEGER :: i, j, k

      IF (verbose) THEN
         WRITE (*, *) 'CLIP: Clipping density field'
         WRITE (*, *) 'CLIP: Threshold:', d0
         WRITE (*, *) 'CLIP: Mesh:', m1, m2, m3
      END IF

      av1 = mean(splay(d, m1, m2, m3))
      var1 = variance(splay(d, m1, m2, m3))
      max1 = maxval(d)

      IF (verbose) THEN
         WRITE (*, *) 'CLIP: Average over-density pre-clipping:', av1
         WRITE (*, *) 'CLIP: Variance in over-density pre-clipping:', var1
         WRITE (*, *) 'CLIP: Maximum density pre-clipping:', max1
      END IF

      !    dep=0.25*(1.+erf(d0/(sqrt(2.*var1))))**2.
      !    IF(verbose==1) WRITE(*,*) 'Expected large-scale power depletion factor:', dep

      ! Now do the clipping
      DO k = 1, m3
         DO j = 1, m2
            DO i = 1, m1
               IF (d(i, j, k) > d0) d(i, j, k) = d0
            END DO
         END DO
      END DO

      IF (verbose) WRITE (*, *) 'CLIP: Density field clipped'

      av2 = mean(splay(d, m1, m2, m3))
      var2 = variance(splay(d, m1, m2, m3))
      max2 = maxval(d)

      IF (verbose) THEN
         WRITE (*, *) 'CLIP: Average over-density post-clipping:', av2
         WRITE (*, *) 'CLIP: Variance in over-density post-clipping:', var2
         WRITE (*, *) 'CLIP: Maximum density post-clipping:', max2
         WRITE (*, *)
      END IF

   END SUBROUTINE clip

   SUBROUTINE anticlip(d, m1, m2, m3, d0, verbose)

      USE statistics

      INTEGER, INTENT(IN) :: m1, m2, m3
      REAL, INTENT(INOUT) :: d(m1, m2, m3) 
      REAL, INTENT(IN) :: d0
      LOGICAL, INTENT(IN) :: verbose
      REAL :: var1, av1, min1, var2, av2, min2
      INTEGER :: i, j, k, m

      IF (verbose) THEN
         WRITE (*, *) 'Anti-clipping over-density field'
         WRITE (*, *) 'Threshold:', d0
         WRITE (*, *) 'Mesh:', m
      END IF

      av1 = mean(splay(d, m1, m2, m3))
      var1 = variance(splay(d, m1, m2, m3))
      min1 = minval(d)

      IF (verbose) THEN
         WRITE (*, *) 'Average over-density pre-clipping:', av1
         WRITE (*, *) 'Variance in over-density pre-clipping:', var1
         WRITE (*, *) 'Minimum over-density pre-clipping:', min1
      END IF

      !    dep=0.25*(1.+erf(d0/(sqrt(2.*var1))))**2.
      !    IF(verbose==1) WRITE(*,*) 'Expected large-scale power depletion factor:', dep

      ! Now do the clipping
      DO k = 1, m
         DO j = 1, m
            DO i = 1, m
               IF (d(i, j, k) < d0) d(i, j, k) = d0
            END DO
         END DO
      END DO

      IF (verbose) WRITE (*, *) 'Over-density field clipped'

      av2 = mean(splay(d, m1, m2, m3))
      var2 = variance(splay(d, m1, m2, m3))
      min2 = minval(d)

      IF (verbose) THEN
         WRITE (*, *) 'Average over-density post-clipping:', av2
         WRITE (*, *) 'Variance in over-density post-clipping:', var2
         WRITE (*, *) 'Minimum over-density post-clipping:', min2
         WRITE (*, *)
      END IF

   END SUBROUTINE anticlip

   INTEGER FUNCTION count_empty_cells(d, m)

      USE precision
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d(m, m, m)
      INTEGER(int8) :: sum
      INTEGER :: i, j, k

      sum = 0
      DO k = 1, m
         DO j = 1, m
            DO i = 1, m
               IF (d(i, j, k) == 0.) THEN
                  sum = sum+1
               END IF
            END DO
         END DO
      END DO

      count_empty_cells = INT(sum)

   END FUNCTION count_empty_cells

   SUBROUTINE compute_power_spectrum_2D(dk1, dk2, m, L, kmin, kmax, nk, k, pow, nmodes, sigma, linear_k_range)

      USE table_integer 

      ! Takes in a dk(m,m) array and computes the power spectrum
      ! NOTE: Leave the double complex as it allows the running to determine complex vs real
      COMPLEX, INTENT(IN) :: dk1(:, :) ! Fourier components of field 1
      COMPLEX, INTENT(IN) :: dk2(:, :) ! Fourier components of field 2
      INTEGER, INTENT(IN) :: m         ! mesh size for fields
      REAL, INTENT(IN) :: L            ! box size [Mpc/h]
      REAL, INTENT(IN) :: kmin         ! minimum and maximum wavenumber [h/Mpc]
      REAL, INTENT(IN) :: kmax         ! minimum and maximum wavenumber [h/Mpc]
      INTEGER, INTENT(IN) :: nk        ! number of k bins
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)         ! Output k values
      REAL, ALLOCATABLE, INTENT(OUT) :: pow(:)       ! Output Delta^2(k) values
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nmodes(:) ! Output Number of modes contributing to the k bin
      REAL, ALLOCATABLE, INTENT(OUT) :: sigma(:)     ! Output varaicnce in bin
      LOGICAL, OPTIONAL, INTENT(IN) :: linear_k_range
      INTEGER :: i, ix, iy, n, mn
      REAL :: kx, ky, kmod, Dk, f, crap
      REAL, ALLOCATABLE :: kbin(:)
      !INTEGER(int8) :: nmodes8(nk)

      REAL, PARAMETER :: dbin = 1e-3 ! Bin slop parameter for first and last bin edges
      LOGICAL, PARAMETER :: logmeank = .FALSE. ! Enable this to assign k to the log-mean of the bin (foolish)

      ERROR STOP 'COMPUTE_POWER_SPECTRUM_2D: Check this subroutine very carefully'
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_2D: Computing isotropic power spectrum'

      ! Set summation variables to 0.d0
      !k8 = 0.d0
      !pow8 = 0.d0
      !nmodes8 = 0
      !sigma8 = 0.d0

      ! Allocate arrays
      ALLOCATE (k(nk), pow(nk), nmodes(nk), sigma(nk))
      k = 0.; pow = 0.; nmodes = 0; sigma = 0.

      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_2D: Binning power'
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_2D: Mesh:', m
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_2D: Bins:', nk
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_2D: k_min [h/Mpc]:', kmin
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_2D: k_max [h/Mpc]:', kmax

      ! Fill array of k bins with linear/log spacing
      IF (present_and_correct(linear_k_range)) THEN
         CALL fill_array(kmin, kmax, kbin, nk+1)
      ELSE
         CALL fill_array_log(kmin, kmax, kbin, nk+1)
      END IF

      ! Explicitly extend the first and last bins to be sure to include *all* modes
      ! This is necessary due to rounding errors!
      kbin(1) = kbin(1)*(1.-dbin)
      kbin(nk+1) = kbin(nk+1)*(1.+dbin)

      ! Cell location of Nyquist
      mn = m/2+1

      ! Loop over all independent elements of dk
      DO iy = 1, m
         DO ix = 1, mn

            ! Cycle for the zero mode (k=0)
            IF (ix == 1 .AND. iy == 1) CYCLE

            ! Cycle for the repeated zero modes and Nyquist modes
            ! I *think* this is correct to avoid double counting zero modes and Nyquist modes
            ! For example 0,1 is the same as 0,-1
            IF ((ix == 1 .OR. ix == mn) .AND. iy > mn) CYCLE

            CALL k_fft(ix, iy, 1, m, kx, ky, crap, kmod, L)

            ! Find integer 'n' in bins from place in table
            IF (kmod >= kbin(1) .AND. kmod <= kbin(nk+1)) THEN
               n = find_table_integer(kmod, kbin, ifind_split)
               IF (n < 1 .OR. n > nk) THEN
                  CYCLE
               ELSE
                  k(n) = k(n)+kmod
                  f = real(conjg(dk1(ix, iy))*dk2(ix, iy))
                  f = f/real(m)**2 ! Note the division by m^4 here
                  pow(n) = pow(n)+f
                  sigma(n) = sigma(n)+f**2
                  nmodes(n) = nmodes(n)+1
               END IF
            END IF

         END DO
      END DO

      ! Deallocate and reallocate arrays
      ! TODO: Do I need to bother deallocating these arrays?
      ! TODO: Should I pass in allocatable arrays or should they already be allocated?
      !IF (ALLOCATED(k))      DEALLOCATE (k)
      !IF (ALLOCATED(pow))    DEALLOCATE (pow)
      !IF (ALLOCATED(nmodes)) DEALLOCATE (nmodes)
      !IF (ALLOCATED(sigma))  DEALLOCATE (sigma)

      ! Now create the power spectrum and k array
      DO i = 1, nk
         IF (nmodes(i) == 0) THEN
            k(i) = sqrt(kbin(i+1)*kbin(i))
            pow(i) = 0.d0
            sigma(i) = 0.d0
         ELSE
            IF (logmeank) THEN
               k(i) = sqrt(kbin(i+1)*kbin(i))
            ELSE
               k(i) = real(k(i))/real(nmodes(i))
            END IF
            pow(i) = pow(i)/real(nmodes(i))
            IF (nmodes(i) == 1) THEN
               sigma(i) = 0
            ELSE
               sigma(i) = sigma(i)/real(nmodes(i)) ! Create <P(k)^2>
               sigma(i) = sqrt(sigma(i)-pow(i)**2) ! Create biased estimate of sigma
               sigma(i) = sigma(i)*real(nmodes(i))/real(nmodes(i)-1) ! Correct for bias
               sigma(i) = sigma(i)/sqrt(real(nmodes(i))) ! Convert to error on the mean
            END IF
            Dk = twopi*(k(i)*L/twopi)**2
            pow(i) = pow(i)*Dk
            sigma(i) = sigma(i)*Dk
         END IF
      END DO

      ! Convert from double to single
      !pow = real(pow8)
      !sigma = real(sigma8)
      !nmodes = INT(nmodes8)

      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_2D: Power computed'
      WRITE (*, *)

   END SUBROUTINE compute_power_spectrum_2D

   !SUBROUTINE pk(dk1,dk2,m,L,kmin,kmax,bins,k,pow,nbin)
   SUBROUTINE compute_power_spectrum_3D(dk1, dk2, m, L, kmin, kmax, nk, k, pow, nmodes, sigma, linear_k_range)

      ! Takes in a dk(m,m,m) array and computes the power spectrum
      ! TODO: Care with large sums
      USE table_integer

      COMPLEX, INTENT(IN) :: dk1(:, :, :) ! Fourier components of field 1
      COMPLEX, INTENT(IN) :: dk2(:, :, :) ! Fourier components of field 2
      INTEGER, INTENT(IN) :: m            ! Mesh size for fields
      REAL, INTENT(IN) :: L               ! Box size [Mpc/h]
      REAL, INTENT(IN) :: kmin            ! Minimum wavenumber [h/Mpc]
      REAL, INTENT(IN) :: kmax            ! Maximum wavenumber [h/Mpc]
      INTEGER, INTENT(IN) :: nk           ! Number of k bins
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)         ! Output k values
      REAL, ALLOCATABLE, INTENT(OUT) :: pow(:)       ! Output Delta^2(k) values
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nmodes(:) ! Output Number of modes contributing to the k bin
      REAL, ALLOCATABLE, INTENT(OUT) :: sigma(:)     ! Output varaicnce in bin
      LOGICAL, OPTIONAL, INTENT(IN) :: linear_k_range
      INTEGER :: i, ix, iy, iz, n, mn
      REAL :: kx, ky, kz, kmod, Dk, f
      REAL, ALLOCATABLE :: kbin(:)

      REAL, PARAMETER :: dbin = 1e-3           ! Bin slop parameter for first and last bin edges
      LOGICAL, PARAMETER :: logmeank = .FALSE. ! Enable this to assign k to the log-mean of the bin (foolish)

      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: Computing isotropic power spectrum'
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: Binning power'
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: Mesh:', m
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: Bins:', nk
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: k_min [h/Mpc]:', kmin
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: k_max [h/Mpc]:', kmax

      ! Fill array of k bins with linear/log spacing
      IF (present_and_correct(linear_k_range)) THEN
         CALL fill_array(kmin, kmax, kbin, nk+1)
      ELSE
         CALL fill_array_log(kmin, kmax, kbin, nk+1)
      END IF

      ! Explicitly extend the first and last bins to be sure to include *all* modes
      ! This is necessary due to rounding errors!
      kbin(1) = kbin(1)*(1.-dbin)
      kbin(nk+1) = kbin(nk+1)*(1.+dbin)

      ! Allocate arrays and set to zero because they are used for large sums
      ALLOCATE (k(nk), pow(nk), nmodes(nk), sigma(nk))
      k = 0.
      pow = 0.
      nmodes = 0
      sigma = 0.

      ! Cell location of Nyquist
      !IF (.NOT. even(m)) STOP 'COMPUTE_POWER_SPECTRUM_3D: Error, m must be even'
      mn = m/2+1

      ! Loop over all independent elements of dk
      DO iz = 1, m
         DO iy = 1, m
            DO ix = 1, mn

               ! Cycle for the zero mode (k=0)
               IF (ix == 1 .AND. iy == 1 .AND. iz == 1) CYCLE

               ! Cycle for the repeated Nyquist modes
               ! I *think* this is correct to avoid double counting zero modes and Nyquist modes
               ! For example (0, 1, 0) is the same as (0, -1, 0)
               IF ((ix == 1 .OR. ix == mn) .AND. (iy > mn .OR. iz > mn)) CYCLE

               ! Get the wavenumbers corresponding to the array position
               CALL k_FFT(ix, iy, iz, m, kx, ky, kz, kmod, L)

               ! Find integer 'n' in bins from place in table
               IF (kmod >= kbin(1) .AND. kmod <= kbin(nk+1)) THEN
                  n = find_table_integer(kmod, kbin, ifind_split)
                  IF (n < 1 .OR. n > nk) THEN
                     CYCLE
                  ELSE
                     k(n) = k(n)+kmod
                     f = real(conjg(dk1(ix, iy, iz))*dk2(ix, iy, iz))
                     f = f/real(m)**6
                     pow(n) = pow(n)+f
                     sigma(n) = sigma(n)+f**2
                     nmodes(n) = nmodes(n)+1
                  END IF
               END IF

            END DO
         END DO
      END DO

      ! Now create the power spectrum and k array
      DO i = 1, nk
         IF (nmodes(i) == 0) THEN
            k(i) = sqrt(kbin(i+1)*kbin(i))
            pow(i) = 0.d0
            sigma(i) = 0.d0
         ELSE
            IF (logmeank) THEN
               k(i) = sqrt(kbin(i+1)*kbin(i))
            ELSE
               k(i) = k(i)/nmodes(i)
            END IF
            pow(i) = pow(i)/nmodes(i) ! Create <P(k)>
            IF (nmodes(i) == 1) THEN
               sigma(i) = 0.
            ELSE
               sigma(i) = sigma(i)/nmodes(i) ! Create <P(k)^2>
               sigma(i) = sigma(i)-pow(i)**2 ! Create biased estimate of variance
               sigma(i) = sigma(i)*nmodes(i)/(nmodes(i)-1) ! Correct for bias
               sigma(i) = sqrt(sigma(i)) ! Create estimate of standard deviation
               sigma(i) = sigma(i)/sqrt(real(nmodes(i))) ! Convert to error-on-the-mean
            END IF
            Dk = 4.*pi*(k(i)*L/twopi)**3
            pow(i) = pow(i)*Dk
            sigma(i) = sigma(i)*Dk
         END IF
      END DO

      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: Delta^2 at min k:', pow(1)
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: Delta^2 at max k:', pow(nk)
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_3D: Power computed'
      WRITE (*, *)

   END SUBROUTINE compute_power_spectrum_3D

   SUBROUTINE compute_multipole_power_spectrum_3D(dk1, dk2, m, L, ipole, izz, kmin, kmax, nk, k, pow, nmodes, sigma, linear_k_range)

      ! Takes in a dk(m,m,m) array and computes the power spectrum
      USE table_integer
      USE special_functions
      COMPLEX, INTENT(IN) :: dk1(:, :, :) ! Fourier components of field 1
      COMPLEX, INTENT(IN) :: dk2(:, :, :) ! Fourier components of field 2
      INTEGER, INTENT(IN) :: m     ! Mesh size for fields
      REAL, INTENT(IN) :: L        ! Box size [Mpc/h]
      INTEGER, INTENT(IN) :: ipole ! Multipole to compute
      INTEGER, INTENT(IN) :: izz   ! Redshift direction
      REAL, INTENT(IN) :: kmin     ! Minimum wavenumber [h/Mpc]
      REAL, INTENT(IN) :: kmax     ! Maximum wavenumber [h/Mpc]
      INTEGER, INTENT(IN) :: nk    ! Number of k bins
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)         ! Output k values
      REAL, ALLOCATABLE, INTENT(OUT) :: pow(:)       ! Output Delta^2(k) values
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nmodes(:) ! Output Number of modes contributing to the k bin
      REAL, ALLOCATABLE, INTENT(OUT) :: sigma(:)     ! Output varaicnce in bin
      LOGICAL, OPTIONAL, INTENT(IN) :: linear_k_range
      INTEGER :: i, ix, iy, iz, n, mn
      REAL :: kx, ky, kz, kmod, Dk, mu, f
      REAL, ALLOCATABLE :: kbin(:)
      REAL, PARAMETER :: dbin = 1e-3           ! Bin slop parameter for first and last bin edges
      LOGICAL, PARAMETER :: logmeank = .FALSE. ! Enable this to assign k to the log-mean of the bin (foolish)

      ERROR STOP 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: Check this subroutine very carefully'
      WRITE (*, *) 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: Computing isotropic power spectrum'

      ! Set summation variables to 0.d0
      !k8 = 0.d0
      !pow8 = 0.d0
      !nmodes8 = 0
      !sigma8 = 0.d0
      ALLOCATE (k(nk), pow(nk), nmodes(nk), sigma(nk))
      k = 0.; pow = 0.; nmodes = 0; sigma = 0.

      WRITE (*, *) 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: Binning power'
      WRITE (*, *) 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: Mesh:', m
      WRITE (*, *) 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: Bins:', nk
      WRITE (*, *) 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: k_min [h/Mpc]:', kmin
      WRITE (*, *) 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: k_max [h/Mpc]:', kmax

      ! Fill array of k bins with linear/log spacing
      IF (present_and_correct(linear_k_range)) THEN
         CALL fill_array(kmin, kmax, kbin, nk+1)
      ELSE
         CALL fill_array_log(kmin, kmax, kbin, nk+1)
      END IF

      ! Explicitly extend the first and last bins to be sure to include *all* modes
      ! This is necessary due to rounding errors!
      kbin(1) = kbin(1)*(1.-dbin)
      kbin(nk+1) = kbin(nk+1)*(1.+dbin)

      ! Cell location of Nyquist
      mn = m/2+1

      ! Loop over all independent elements of dk
      DO iz = 1, m
         DO iy = 1, m
            DO ix = 1, mn

               ! Cycle for the zero mode (k=0)
               IF (ix == 1 .AND. iy == 1 .AND. iz == 1) CYCLE

               ! Cycle for the repeated zero modes and Nyquist modes
               ! I *think* this is correct to avoid double counting zero modes and Nyquist modes
               ! For example 0,1,0 is the same as 0,-1,0
               IF ((ix == 1 .OR. ix == mn) .AND. (iy > mn .OR. iz > mn)) CYCLE

               CALL k_fft(ix, iy, iz, m, kx, ky, kz, kmod, L)

               IF (izz == 1) THEN
                  mu = kx/kmod
               ELSE IF (izz == 2) THEN
                  mu = ky/kmod
               ELSE IF (izz == 3) THEN
                  mu = kz/kmod
               ELSE
                  STOP 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: Error, iz not specified correctly'
               END IF

               ! Find integer 'n' in bins from place in table
               IF (kmod >= kbin(1) .AND. kmod <= kbin(nk+1)) THEN
                  n = find_table_integer(kmod, kbin, ifind_split)
                  IF (n < 1 .OR. n > nk) THEN
                     CYCLE
                  ELSE
                     k(n) = k(n)+kmod
                     f = real(conjg(dk1(ix, iy, iz))*dk2(ix, iy, iz))
                     f = f/real(m)**6
                     f = f*Legendre_polynomial(ipole, mu)*(2.*ipole+1.)
                     pow(n) = pow(n)+f
                     sigma(n) = sigma(n)+f**2
                     nmodes(n) = nmodes(n)+1
                  END IF
               END IF

            END DO
         END DO
      END DO

      ! Deallocate and reallocate arrays
      ! TODO: Do I need to bother deallocating these arrays?
      ! TODO: Should I pass in allocatable arrays or should they already be allocated?
      !IF (ALLOCATED(k))      DEALLOCATE (k)
      !IF (ALLOCATED(pow))    DEALLOCATE (pow)
      !IF (ALLOCATED(nmodes)) DEALLOCATE (nmodes)
      !IF (ALLOCATED(sigma))  DEALLOCATE (sigma)
      !ALLOCATE (k(nk), pow(nk), nmodes(nk), sigma(nk))

      ! Now create the power spectrum and k array
      DO i = 1, nk
         IF (nmodes(i) == 0) THEN
            k(i) = sqrt(kbin(i+1)*kbin(i))
            pow(i) = 0.d0
            sigma(i) = 0.d0
         ELSE
            IF (logmeank) THEN
               k(i) = sqrt(kbin(i+1)*kbin(i))
            ELSE
               k(i) = real(k(i))/real(nmodes(i))
            END IF
            pow(i) = pow(i)/real(nmodes(i)) ! Create <P(k)>
            IF (nmodes(i) == 1) THEN
               sigma(i) = 0.
            ELSE
               sigma(i) = sigma(i)/real(nmodes(i)) ! Create <P(k)^2>
               sigma(i) = sigma(i)-pow(i)**2 ! Create biased estimate of variance
               sigma(i) = sigma(i)*real(nmodes(i))/real(nmodes(i)-1) ! Correct for bias
               sigma(i) = sqrt(sigma(i)) ! Create estimate of standard deviation
               sigma(i) = sigma(i)/sqrt(real(nmodes(i))) ! Convert to error-on-the-mean
            END IF
            Dk = 4.*pi*(k(i)*L/twopi)**3
            pow(i) = pow(i)*Dk
            sigma(i) = sigma(i)*Dk
         END IF
      END DO

      ! Convert from double to single
      !pow = real(pow8)
      !sigma = real(sigma8)
      !nmodes = INT(nmodes8)

      WRITE (*, *) 'COMPUTE_MULTIPOLE_POWER_SPECTRUM_3D: Power computed'
      WRITE (*, *)

   END SUBROUTINE compute_multipole_power_spectrum_3D

   SUBROUTINE compute_power_spectrum_pole(d, m, L, ipole, iz, kmin, kmax, nk, kval, pow, nmodes)

      USE special_functions
      INTEGER, INTENT(IN) :: m
      COMPLEX, INTENT(IN) :: d(m, m, m)
      REAL, INTENT(IN) :: kmin, kmax, L
      INTEGER, INTENT(IN) :: iz, ipole, nk
      REAL, ALLOCATABLE, INTENT(OUT) :: pow(:), kval(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nmodes(:)
      INTEGER :: i, j, k, n
      REAL :: kx, ky, kz, kmod, mu
      REAL :: kbin(nk+1)

      STOP 'COMPUTE_POWER_SPECTRUM_POLE: Check this subroutine very carefully, see notes in code'
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_POLE: Computing isotropic power spectrum'

      ALLOCATE(pow(nk), kval(nk), nmodes(nk))
      kval = 0.; pow = 0.; nmodes = 0

      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_POLE: Binning power'
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_POLE: Bins:', nk
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_POLE: k_min:', kmin
      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_POLE: k_max:', kmax

      ! Log-spaced bins
      ! TODO: Surely this needs to be in a loop? Use allocate function?
      kbin(i) = progression(kmin, kmax, i, nk+1, ilog=.TRUE.)

      !Explicitly extend the bins to be sure to include all modes
      !This is necessary due to rounding errors!
      !    kbin(1)=kbin(1)*0.999
      !    kbin(bins+1)=kbin(bins+1)*1.001
      !m=SIZE(d(:,1,1))

      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_POLE: Mesh:', m

      DO k = 1, m
         DO j = 1, m
            DO i = 1, m

               IF (i == 1 .AND. j == 1 .AND. k == 1) CYCLE

               CALL k_fft(i, j, k, m, kx, ky, kz, kmod, L)

               IF (iz == 1) THEN
                  mu = kx/kmod
               ELSE IF (iz == 2) THEN
                  mu = ky/kmod
               ELSE IF (iz == 3) THEN
                  mu = kz/kmod
               END IF

               !             DO o=1,bins
               !                IF(kmod>=kbin(o) .AND. kmod<=kbin(o+1)) THEN
               !                   pow8(o)=pow8(o)+(abs(d(i,j,k))**2.)*legendre(ipole,mu)*(2.*float(ipole)+1.)!/2.
               !                   kval(o)=kval(o)+kmod
               !                   nbin8(o)=nbin8(o)+1
               !                   EXIT
               !                END IF
               !             END DO

               ! Find integer automatically from place in table. Assumes log-spaced bins
               ! Implemented (27/08/15) so could be a source of bugs
               ! Differences will appear due to k modes that are on the boundary
               n = 1+floor(real(nk)*log(kmod/kmin)/log(kmax/kmin))
               IF (n < 1 .OR. n > nk) THEN
                  CYCLE
               ELSE
                  pow(n) = pow(n)+(abs(d(i, j, k))**2)*Legendre_polynomial(ipole, mu)*(2.*ipole+1.)
                  kval(n) = kval(n)+kmod
                  nmodes(n) = nmodes(n)+1
               END IF

            END DO
         END DO
      END DO

      !Deallocate and re-allocate arrays
      !IF (ALLOCATED(kval)) DEALLOCATE (kval)
      !IF (ALLOCATED(pow)) DEALLOCATE (pow)
      !IF (ALLOCATED(nmodes)) DEALLOCATE (nmodes)
      !ALLOCATE (kval(nk), pow(nk), nmodes(nk))

      DO i = 1, nk
         kval(i) = (kbin(i+1)+kbin(i))/2. ! Should be sqrt(kmax*kmin) for log
         IF (nmodes(i) == 0) THEN
            pow(i) = 0.
         ELSE
            pow(i) = pow(i)/real(nmodes(i))
            pow(i) = pow(i)*((L*kval(i))**3)/(2.*pi**2) ! Not sure this is good; create P(k) and then Delta(k)?
         END IF
      END DO

      pow = pow/real(m)**6

      !Divide by 2 because double count Hermitian conjugates
      nmodes = nmodes/2

      WRITE (*, *) 'COMPUTE_POWER_SPECTRUM_POLE: Power computed'
      WRITE (*, *)

   END SUBROUTINE compute_power_spectrum_pole

   ! SUBROUTINE compute_power_spectrum_rsd(d, L, kmin, kmax, nk, kv, mu, pow, nmodes, iz)

   !    INTEGER :: i, j, k, m, ii, jj, nk, iz
   !    REAL :: kx, ky, kz, kmod, L, kmin, kmax, a, b, mus
   !    REAL :: pow(nk, nk), kv(nk), kbin(nk+1), mu(nk), mubin(nk+1)
   !    !DOUBLE PRECISION :: pow8(nk, nk)
   !    INTEGER :: nmodes(nk, nk)
   !    !INTEGER*8 :: nmodes8(nk, nk)
   !    INTEGER(int8) :: nmodes8(nk, nk)
   !    !DOUBLE COMPLEX :: d(:, :, :)
   !    COMPLEX :: d(:, :, :)

   !    STOP 'COMPUTE_POWER_SPECTRUM_RSD: Check this very carefully'

   !    WRITE (*, *) 'Computing RSD power spectrum'

   !    kbin = 0.
   !    mubin = 0.
   !    kv = 0.
   !    mu = 0.
   !    pow = 0.
   !    nmodes = 0

   !    pow8 = 0.d0
   !    nmodes8 = 0

   !    WRITE (*, *) 'Binning power'
   !    WRITE (*, *) 'Bins:', nk
   !    WRITE (*, *) 'k_min:', kmin
   !    WRITE (*, *) 'k_max:', kmax

   !    a = kmin
   !    b = kmax

   !    a = log10(a)
   !    b = log10(b)

   !    !DO i=1,bins+1
   !    !   kbin(i)=a+(b-a)*float(i-1)/float(bins)
   !    !END DO
   !    kbin(i) = progression(a, b, i, nk+1)

   !    !DO i=1,bins+1
   !    !   mubin(i)=float(i-1)/float(bins)
   !    !END DO
   !    mubin(i) = progression(0., 1., i, nk+1)

   !    DO i = 1, nk
   !       kv(i) = (kbin(i)+kbin(i+1))/2.
   !       mu(i) = (mubin(i)+mubin(i+1))/2.
   !    END DO

   !    kbin = 10.**kbin
   !    kv = 10.**kv

   !    !Explicitly extend the bins to be sure to include all modes
   !    !This is necessary due to rounding errors!
   !    kbin(1) = kbin(1)*0.999
   !    kbin(nk+1) = kbin(nk+1)*1.001
   !    mubin(1) = -0.001
   !    mubin(nk+1) = 1.001

   !    m = SIZE(d(:, 1, 1))

   !    WRITE (*, *) 'Mesh:', m

   !    DO k = 1, m
   !       DO j = 1, m
   !          DO i = 1, m

   !             IF (i == 1 .AND. j == 1 .AND. k == 1) CYCLE

   !             CALL k_fft(i, j, k, m, kx, ky, kz, kmod, L)

   !             IF (iz == 1) THEN
   !                mus = kx/kmod
   !             ELSE IF (iz == 2) THEN
   !                mus = ky/kmod
   !             ELSE IF (iz == 3) THEN
   !                mus = kz/kmod
   !             ELSE
   !                STOP 'COMPUTE_POWER_SPECTRUM_RSD: Error, iz not specified correctly'
   !             END IF

   !             mus = abs(mus)

   !             !             WRITE(*,*) mus, kbin(1), kbin(2), kmod
   !             !             IF(i==10) STOP

   !             DO jj = 1, nk
   !                IF (kmod >= kbin(jj) .AND. kmod <= kbin(jj+1)) THEN
   !                   DO ii = 1, nk
   !                      IF (mus >= mubin(ii) .AND. mus <= mubin(ii+1)) THEN
   !                         pow8(ii, jj) = pow8(ii, jj)+abs(d(i, j, k))**2.
   !                         nmodes8(ii, jj) = nmodes8(ii, jj)+1
   !                         EXIT
   !                      END IF
   !                   END DO
   !                   EXIT
   !                END IF
   !             END DO

   !          END DO
   !       END DO
   !    END DO

   !    DO jj = 1, nk
   !       DO ii = 1, nk
   !          IF (nmodes8(ii, jj) == 0) THEN
   !             pow8(ii, jj) = 0.
   !          ELSE
   !             pow8(ii, jj) = pow8(ii, jj)/real(nmodes8(ii, jj))
   !             pow8(ii, jj) = pow8(ii, jj)*((L*kv(jj))**3.)/(2.*pi**2.)
   !          END IF
   !       END DO
   !    END DO

   !    pow = real(pow8)/(real(m)**6)

   !    !Divide by 2 because double count Hermitian conjugates
   !    nmodes = INT(nmodes8)/2

   !    WRITE (*, *) 'Power computed'
   !    WRITE (*, *)

   ! END SUBROUTINE compute_power_spectrum_rsd

   ! SUBROUTINE compute_power_spectrum_rsd2(d, L, kmin, kmax, nk, kpar, kper, pow, nmodes, iz)

   !    INTEGER :: i, j, k, m, ii, jj, nk, iz
   !    REAL :: kx, ky, kz, kmod, L, kmin, kmax, a, b, kpers, kpars
   !    REAL :: pow(nk, nk), kpar(nk), kparbin(nk+1), kper(nk), kperbin(nk+1)
   !    !DOUBLE PRECISION :: pow8(nk, nk)
   !    REAL :: pow8(nk, nk)
   !    INTEGER :: nmodes(nk, nk)
   !    !INTEGER*8 :: nmodes8(nk, nk)
   !    INTEGER(int8) :: nmodes8(nk, nk)
   !    !DOUBLE COMPLEX :: d(:, :, :)
   !    COMPLEX :: d(:, :, :)

   !    STOP 'COMPUTE_POWER_SPECTRUM_RSD2: Check this very carefully'

   !    WRITE (*, *) 'Computing rsd power spectrum'

   !    kparbin = 0.
   !    kperbin = 0.
   !    kpar = 0.
   !    kper = 0.
   !    pow = 0.
   !    nmodes = 0

   !    pow8 = 0.d0
   !    nmodes8 = 0

   !    WRITE (*, *) 'Binning power'
   !    WRITE (*, *) 'Bins:', nk
   !    WRITE (*, *) 'k_min:', kmin
   !    WRITE (*, *) 'k_max:', kmax

   !    a = kmin
   !    b = kmax

   !    a = log10(a)
   !    b = log10(b)

   !    !DO i=1,bins+1
   !    !   kparbin(i)=a+(b-a)*float(i-1)/float(bins)
   !    !END DO
   !    kparbin(i) = progression(a, b, i, nk+1)

   !    DO i = 1, nk
   !       kpar(i) = (kparbin(i)+kparbin(i+1))/2.
   !    END DO

   !    kparbin = 10.**kparbin
   !    kpar = 10.**kpar

   !    !Explicitly extend the bins to be sure to include all modes
   !    !This is necessary due to rounding errors!
   !    kparbin(1) = kparbin(1)*0.999
   !    kparbin(nk+1) = kparbin(nk+1)*1.001

   !    kperbin = kparbin
   !    kper = kpar

   !    m = SIZE(d(:, 1, 1))

   !    WRITE (*, *) 'Mesh:', m

   !    DO k = 1, m
   !       DO j = 1, m
   !          DO i = 1, m

   !             IF (i == 1 .AND. j == 1 .AND. k == 1) CYCLE

   !             CALL k_fft(i, j, k, m, kx, ky, kz, kmod, L)

   !             IF (iz == 1) THEN
   !                kpars = abs(kx)
   !                kpers = sqrt(ky**2.+kz**2.)
   !             ELSE IF (iz == 2) THEN
   !                kpars = abs(ky)
   !                kpers = sqrt(kz**2.+kx**2.)
   !             ELSE IF (iz == 3) THEN
   !                kpars = abs(kz)
   !                kpers = sqrt(kx**2.+ky**2.)
   !             ELSE
   !                STOP 'COMPUTE_POWER_SPECTRUM_RSD2: Error, iz not specified correctly'
   !             END IF

   !             DO jj = 1, nk
   !                IF (kpars >= kparbin(jj) .AND. kpars <= kparbin(jj+1)) THEN
   !                   DO ii = 1, nk
   !                      IF (kpers >= kperbin(ii) .AND. kpers <= kperbin(ii+1)) THEN
   !                         pow8(ii, jj) = pow8(ii, jj)+abs(d(i, j, k))**2.
   !                         nmodes8(ii, jj) = nmodes8(ii, jj)+1
   !                         EXIT
   !                      END IF
   !                   END DO
   !                   EXIT
   !                END IF
   !             END DO

   !          END DO
   !       END DO
   !    END DO

   !    !    DO jj=1,bins
   !    !       DO ii=1,bins
   !    !          pow(ii,jj)=pow(ii,jj)/log(kperbin(ii+1)/kperbin(ii))
   !    !          pow(ii,jj)=pow(ii,jj)/log(kparbin(jj+1)/kparbin(jj))
   !    !       END DO
   !    !    END DO

   !    DO jj = 1, nk
   !       DO ii = 1, nk
   !          IF (nmodes8(ii, jj) == 0) THEN
   !             pow8(ii, jj) = 0.
   !          ELSE
   !             pow8(ii, jj) = pow8(ii, jj)/real(nmodes8(ii, jj))
   !             pow8(ii, jj) = pow8(ii, jj)*(L**3.*kpar(jj)*kper(ii)**2.)/(2.*pi**2.)
   !          END IF
   !       END DO
   !    END DO

   !    pow = real(pow8)/(real(m)**6)

   !    !Divide by 2 because double count Hermitian conjugates
   !    nmodes = INT(nmodes8)/2

   !    WRITE (*, *) 'Power computed'
   !    WRITE (*, *)

   ! END SUBROUTINE compute_power_spectrum_rsd2

   REAL FUNCTION box_mode_power(dk)

      COMPLEX, INTENT(IN) :: dk(:, :, :)    

      box_mode_power = real(abs(dk(2, 1, 1))**2+abs(dk(1, 2, 1))**2+abs(dk(1, 1, 2))**2)/3.

   END FUNCTION box_mode_power

   FUNCTION periodic_distance(x1, x2, L)

      ! Calculates the distance between x1 and x2 assuming that they are coordinates in a periodic box
      ! This is in field_operations because it needs coordinates, *not* necessarily particles
      REAL :: periodic_distance
      REAL, INTENT(IN) :: x1(3), x2(3), L
      REAL :: dx(3)
      INTEGER :: i

      ! Initially dx is just the absolute vector difference
      dx = abs(x2-x1)

      ! Now check if any legs are greater than half-box size
      ! Note the Cartesian distance *cannot* be larger than L/2
      DO i = 1, 3
         IF (dx(i) > L/2.) THEN
            dx(i) = L-dx(i)
         END IF
      END DO

      periodic_distance = sqrt(dx(1)**2+dx(2)**2+dx(3)**2)

   END FUNCTION periodic_distance

   FUNCTION periodic_mean(x1, x2, L)

      ! Calculates the periodic mean of two coordinates in a box
      ! This is in field_operations because it needs coordinates, *not* necessarily particles
      REAL :: periodic_mean(3)
      REAL, INTENT(IN) :: x1(3), x2(3), L
      REAL :: dx(3)
      INTEGER :: i

      ! Initially dx is just the absolute vector difference
      dx = abs(x2-x1)

      DO i = 1, 3
         periodic_mean(i) = 0.5*(x1(i)+x2(i))
         IF (dx(i) > L/2.) THEN
            periodic_mean(i) = periodic_mean(i)+L/2.
         END IF
      END DO

   END FUNCTION periodic_mean

   SUBROUTINE field_correlation_function(r_array, xi_array, n_array, n, d, m, L)

      ! This double counts, so time could be at least halved
      ! Also could be parrallelised
      ! Also could just not be complete shit, but it should get the job done
      USE precision
      USE table_integer

      INTEGER, INTENT(IN) :: n, m
      REAL, INTENT(OUT) :: xi_array(n)
      REAL, INTENT(IN) :: L, d(m, m, m), r_array(n)
      !INTEGER*8, INTENT(OUT) :: n_array(n)
      INTEGER(int8), INTENT(OUT) :: n_array(n)
      REAL:: rmin, rmax
      !DOUBLE PRECISION :: xi8_array(n)
      INTEGER :: i1, i2, i3, j1, j2, j3, i(3), j(3), k, dim
      REAL :: r, x1(3), x2(3)

      ERROR STOP 'FIELD_CORRELATION_FUNCTION: Check this very carefully'

      rmin = r_array(1)
      rmax = r_array(n)

      WRITE (*, *) 'FIELD_CORRELATION_FUNCTION: rmin [Mpc/h]:', rmin
      WRITE (*, *) 'FIELD_CORRELATION_FUNCTION: rmax [Mpc/h]:', rmax
      WRITE (*, *) 'FIELD_CORRELATION_FUNCTION: number of r bins:', n

      xi_array = 0.d0
      n_array = 0

      DO i3 = 1, m
         DO i2 = 1, m
            DO i1 = 1, m

               i(1) = i1
               i(2) = i2
               i(3) = i3

               DO dim = 1, 3
                  x1(dim) = cell_position(i(dim), L, m)
               END DO

               DO j3 = 1, m
                  DO j2 = 1, m
                     DO j1 = 1, m

                        j(1) = j1
                        j(2) = j2
                        j(3) = j3

                        DO dim = 1, 3
                           x2(dim) = cell_position(j(dim), L, m)
                        END DO

                        r = periodic_distance(x1, x2, L)

                        IF (r < rmin .OR. r > rmax) THEN
                           CYCLE
                        ELSE
                           k = find_table_integer(r, r_array, ifind_split)
                           IF (k < 1 .OR. k > n) ERROR STOP 'FIELD_CORRELATION_FUNCTION: Integer finding has failed'
                           xi_array(k) = xi_array(k)+d(i(1), i(2), i(3))*d(j(1), j(2), j(3))
                           n_array(k) = n_array(k)+1
                        END IF

                     END DO
                  END DO
               END DO

            END DO
         END DO
      END DO
      xi_array = xi_array/n_array

      WRITE (*, *) 'FIELD_CORRELATION_FUNCTION: done'
      WRITE (*, *)

   END SUBROUTINE field_correlation_function

   INTEGER FUNCTION empty_cells_3D(d, m)

      USE precision

      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: d(m, m, m)   
      !INTEGER*8 :: sum
      INTEGER(int8) :: sum
      INTEGER :: i, j, k

      sum = 0
      DO k = 1, m
         DO j = 1, m
            DO i = 1, m
               IF (d(i, j, k) == 0.) THEN
                  sum = sum+1
               END IF
            END DO
         END DO
      END DO

      empty_cells_3D = int(sum)

   END FUNCTION empty_cells_3D

END MODULE field_operations
