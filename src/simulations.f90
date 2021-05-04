MODULE simulations

   ! This module should contain only routines that pertain to simualtion particles specifically
   ! Anything that involves only the fields should go in field_operations.f90
   ! Each routine should take particle properties (e.g., positions) as an argument

   USE constants
   USE basic_operations
   USE array_operations
   USE table_integer
   USE interpolate
   USE field_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: particle_bin
   PUBLIC :: particle_bin_average
   PUBLIC :: correlation_function
   PUBLIC :: create_mass_function
   PUBLIC :: halo_mass_cut
   PUBLIC :: halo_mass_weights
   PUBLIC :: compute_and_write_power_spectrum
   PUBLIC :: write_power_spectrum
   PUBLIC :: write_density_slice_ascii
   PUBLIC :: Zeldovich_ICs
   PUBLIC :: replace
   PUBLIC :: sharp_Fourier_density_contrast
   PUBLIC :: power_spectrum_particles
   PUBLIC :: cross_spectrum_particles
   PUBLIC :: folded_power_spectrum_particles
   PUBLIC :: combine_folded_power
   PUBLIC :: fold_particles

   PUBLIC :: generate_randoms
   PUBLIC :: generate_poor_glass
   PUBLIC :: generate_grid

   PUBLIC :: sparse_sample

   PUBLIC :: random_translation
   PUBLIC :: random_inversion
   PUBLIC :: random_rotation

   PUBLIC :: shot_noise
   PUBLIC :: shot_noise_simple
   PUBLIC :: shot_noise_mass
   PUBLIC :: shot_noise_k

   PUBLIC :: SOD
   PUBLIC :: find_pairs
   PUBLIC :: zshift
   PUBLIC :: write_slice_ascii
   PUBLIC :: write_adaptive_field

   PUBLIC :: random_halo_particles
   PUBLIC :: make_HOD

   ! Random NFW profile generation
   INTEGER, PARAMETER :: nr_random_NFW = 65
   INTEGER, PARAMETER :: iorder_random_NFW = 3
   INTEGER, PARAMETER :: ifind_random_NFW = ifind_split
   INTEGER, PARAMETER :: iinterp_random_NFW = iinterp_Lagrange

   INTERFACE particle_bin
      MODULE PROCEDURE particle_bin_2D
      MODULE PROCEDURE particle_bin_3D
   END INTERFACE particle_bin

   INTERFACE particle_bin_average
      MODULE PROCEDURE particle_bin_average_3D
   END INTERFACE particle_bin_average

   INTERFACE NGP
      MODULE PROCEDURE NGP_2D
      MODULE PROCEDURE NGP_3D
   END INTERFACE NGP

   INTERFACE CIC
      MODULE PROCEDURE CIC_2D
      MODULE PROCEDURE CIC_3D
   END INTERFACE CIC

CONTAINS

   SUBROUTINE correlation_function(rmin, rmax, r, xi, n, nr, x1, x2, w1, w2, n1, n2, L)

      ! Calculate the correlation function for the sample with positions 'x' and weights 'w'
      ! CARE: Removed double-precision sums
      USE table_integer
      INTEGER, INTENT(IN) :: n1, n2            ! Total numbers of particles
      INTEGER, INTENT(IN) :: nr                ! Required number of log-spaced bins
      REAL, INTENT(IN) :: rmin, rmax           ! Maximum and maximum distances to calculate xi for [Mpc/h]
      REAL, INTENT(OUT) :: r(nr)               ! Output array of r [Mpc/h]
      REAL, INTENT(OUT) :: xi(nr)              ! Output array of xi
      INTEGER, INTENT(OUT) :: n(nr)            ! Output array of number of pairs in bin     
      REAL, INTENT(IN) :: x1(3, n1), x2(3, n2) ! Particle position arrays [Mpc/h]
      REAL, INTENT(IN) :: w1(n1), w2(n2)       ! Particle weight arrays     
      REAL, INTENT(IN) :: L                    ! Periodic box size [Mpc/h]
      INTEGER :: i, i1, i2
      REAL, ALLOCATABLE :: rbin(:)
      REAL :: rpair, V, nbar1, nbar2, sum1, sum2
      REAL :: r8(nr), xi8(nr) ! TODO: Remove this crap with the 8 suffix

      INTEGER, PARAMETER :: ifind = ifind_split

      ! Allocate arrays and fill the bin edges
      CALL fill_array(log(rmin), log(rmax), rbin, nr+1)
      rbin = exp(rbin)

      ! Set these to zero as they will be used for sums
      r8 = 0.
      xi8 = 0.
      n = 0

      WRITE (*, *) 'CORRELATION_FUNCTION: Calculating correlation function'
      WRITE (*, *) 'CORRELATION_FUNCTION: Summing directly over pairs (slow)'

      ! Loop over first particle
      DO i1 = 1, n1

         ! Do not waste time with zero weights
         IF (w1(i1) == 0.) CYCLE

         ! Loop over second particle
         DO i2 = 1, n2

            ! Do not waste time with zero weights
            IF (w2(i2) == 0.) CYCLE

            ! Calculate the distance between the particle pair
            rpair = periodic_distance(x1(:, i1), x2(:, i2), L)

            ! Bin depending on the distance
            IF (rpair < rmin .OR. rpair > rmax) THEN
               CYCLE
            ELSE
               i = find_table_integer(rpair, rbin, ifind)
               r8(i) = r8(i)+rpair
               xi8(i) = xi8(i)+w1(i1)*w2(i2)
               n(i) = n(i)+1
            END IF

         END DO
      END DO

      WRITE (*, *) 'CORRELATION_FUNCTION: Done sum'

      sum1 = sum(w1)
      sum2 = sum(w2)
      nbar1 = sum1/L**3
      nbar2 = sum2/L**3

      WRITE (*, *) 'CORRELATION_FUNCTION: sum 1:', sum1
      WRITE (*, *) 'CORRELATION_FUNCTION: sum 2:', sum2
      WRITE (*, *) 'CORRELATION_FUNCTION: nbar 1:', nbar1
      WRITE (*, *) 'CORRELATION_FUNCTION: nbar 2:', nbar2
      WRITE (*, *) 'CORRELATION_FUNCTION: Doing final bit (not sure if this is correct)'

      ! Calculate the actual correlation functions from the sums above
      DO i = 1, nr
         IF (n(i) == 0) THEN
            r(i) = sqrt(rbin(i+1)*rbin(i))
            xi(i) = -1.
         ELSE
            r(i) = real(r8(i))/real(n(i))
            V = 4.*pi*(rbin(i+1)-rbin(i))*r(i)**2
            xi(i) = -1.+real(xi8(i))/(sum1*sum2*V/L**3)
         END IF
      END DO

      WRITE (*, *) 'CORRELATION_FUNCTION: Done'
      WRITE (*, *)

   END SUBROUTINE correlation_function

   SUBROUTINE create_mass_function(mmin, mmax, halo_masses, n_haloes, mass_bins, mass_function, n_bins, L)

      USE statistics

      INTEGER, INTENT(IN) :: n_haloes                    ! Total number of haloes
      REAL, INTENT(IN) :: mmin                           ! Minimum halo mass for the mass function [Msun/h]
      REAL, INTENT(IN) :: mmax                           ! Maximum halo mass for the mass function [Msun/h]
      REAL, INTENT(IN) :: halo_masses(n_haloes)          ! Array of halo masses [Msun/h] 
      REAL, ALLOCATABLE, INTENT(OUT) :: mass_bins(:)     ! Mass function bin values [Msun/h]
      REAL, ALLOCATABLE, INTENT(OUT) :: mass_function(:) ! Mass function
      INTEGER, INTENT(IN) :: n_bins                      ! Number of bins in mass function
      REAL, INTENT(IN) :: L                              ! Simulation volume [Mpc/h]
      INTEGER :: i
      REAL, ALLOCATABLE :: mass_bin_edges(:)
      INTEGER, ALLOCATABLE :: counts(:)

      CALL histogram(log(mmin), log(mmax), mass_bin_edges, counts, n_bins, log(halo_masses))
      mass_bin_edges = exp(mass_bin_edges)

      ALLOCATE (mass_bins(n_bins))
      ALlOCATE (mass_function(n_bins))

      DO i = 1, n_bins
         mass_bins(i) = sqrt(mass_bin_edges(i)*mass_bin_edges(i+1))
         mass_function(i) = real(counts(i))/((mass_bin_edges(i+1)-mass_bin_edges(i))*L**3)
      END DO

   END SUBROUTINE create_mass_function

   SUBROUTINE halo_mass_cut(mmin, mmax, x, m, n)

      REAL, INTENT(IN) :: mmin
      REAL, INTENT(IN) :: mmax
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:, :)
      REAL, ALLOCATABLE, INTENT(INOUT) :: m(:)
      INTEGER, INTENT(INOUT) :: n
      REAL :: x_store(3, n), m_store(n)
      INTEGER :: i, j, n_store

      WRITE (*, *) 'HALO_MASS_CUT: Number of haloes before cut:', n
      IF (minval(m) == 0.) THEN
         WRITE (*, *) 'HALO_MASS_CUT: Minimum halo mass before cut [Msun/h]:', minval(m)
      ELSE
         WRITE (*, *) 'HALO_MASS_CUT: Minimum halo mass before cut log10([Msun/h]):', log10(minval(m))
      END IF    
      WRITE (*, *) 'HALO_MASS_CUT: Maximum halo mass before cut log10([Msun/h]):', log10(maxval(m))
      WRITE (*, *) 'HALO_MASS_CUT: Minimum mass for cut log10([Msun/h]):', log10(mmin)
      WRITE (*, *) 'HALO_MASS_CUT: Maximum mass for cut log10([Msun/h]):', log10(mmax)

      ! Initially store the initial values of the inputs
      x_store = x
      m_store = m
      n_store = n

      ! Find how many haloes satisfy this cut
      n = 0
      DO i = 1, n_store
         IF (m_store(i) > mmin .AND. m_store(i) < mmax) THEN
            n = n+1
         END IF
      END DO

      WRITE (*, *) 'HALO_MASS_CUT: Number of haloes after cut:', n

      ! Deallocate and reallocate input arrays
      DEALLOCATE (x, m)
      ALLOCATE (x(3, n), m(n))

      ! Now fill up arrays for output
      j = 0
      DO i = 1, n_store
         IF (m_store(i) > mmin .AND. m_store(i) < mmax) THEN
            j = j+1
            x(:, j) = x_store(:, i)
            m(j) = m_store(i)
         END IF
      END DO

      WRITE (*, *) 'HALO_MASS_CUT: Minimum halo mass after cut log10([Msun/h]):', log10(minval(m))
      WRITE (*, *) 'HALO_MASS_CUT: Maximum halo mass after cut log10([Msun/h]):', log10(maxval(m))
      WRITE (*, *) 'HALO_MASS_CUT: Done'
      WRITE (*, *)

   END SUBROUTINE halo_mass_cut

   SUBROUTINE halo_mass_weights(mmin, mmax, m, w)

      ! Set the weight, w(i), to one if the halo is in the mass range, zero otherwise
      REAL, INTENT(IN) :: mmin
      REAL, INTENT(IN) :: mmax
      REAL, INTENT(IN) :: m(:)
      REAL, INTENT(OUT) :: w(:)
      INTEGER :: i, n

      n = size(m)
      IF (n .NE. size(w)) STOP 'HALO_MASS_WEIGHTS: Error, m and w must be the same size'

      w = 0.
      DO i = 1, n
         IF (m(i) >= mmin .AND. m(i) < mmax) w(i) = 1.
      END DO

   END SUBROUTINE halo_mass_weights

   SUBROUTINE compute_and_write_power_spectrum(x, L, m, nk, outfile)

      ! Compute and write the power spectrum out in some standard format
      REAL, INTENT(IN) :: x(:, :)    
      REAL, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(IN) :: nk
      CHARACTER(len=*), INTENT(IN) :: outfile
      REAL, ALLOCATABLE :: k(:), Pk(:), sig(:), shotk(:)
      INTEGER, ALLOCATABLE :: nbin(:)
      REAL :: shot
      INTEGER :: n

      IF (size(x, 1) .NE. 3) STOP 'COMPUTE_AND_WRITE_POWER_SPECTRUM: Error, needs to be 3D'
      n = size(x, 2)

      ! Compute the power spectrum from the particle positions
      ! This is only correct if all particles have the same mass
      CALL power_spectrum_particles(x, L, m, nk, k, Pk, nbin, sig)

      ! Compute the shot noise assuming all particles have equal mass
      shot = shot_noise_simple(L, n)
      shotk = shot_noise_k(k, shot)

      CALL write_power_spectrum(k, Pk, shotk, nbin, sig, outfile)

   END SUBROUTINE compute_and_write_power_spectrum

   SUBROUTINE write_power_spectrum(k, Pk, shot, nbin, sig, outfile)

      REAL, INTENT(IN) :: k(:)
      REAL, INTENT(IN) :: Pk(:)
      REAL, INTENT(IN) :: shot(:)
      INTEGER, INTENT(IN) :: nbin(:)
      REAL, INTENT(IN) :: sig(:)
      CHARACTER(len=*), INTENT(IN) :: outfile
      INTEGER :: nk
      INTEGER :: u
      INTEGER :: ik

      ! Number of k elements
      nk = size(k)

      ! Write to screen
      WRITE (*, *) 'WRITE_POWER_SPECTRUM: Outfile: ', trim(outfile)
      WRITE (*, *) 'WRITE_POWER_SPECTRUM: Number of k points:', nk

      ! Write to file in standard format
      OPEN (newunit=u, file=outfile)
      DO ik = 1, nk
         IF (nbin(ik) == 0) CYCLE
         WRITE (u, *) k(ik), Pk(ik), shot(ik), nbin(ik), sig(ik)
      END DO
      CLOSE (u)

      ! Write to screen
      WRITE (*, *) 'WRITE_POWER_SPECTRUM: Done'
      WRITE (*, *)

   END SUBROUTINE write_power_spectrum

   SUBROUTINE power_spectrum_particles(x, L, m, nk, k, Pk, nbin, sig, kmin_opt, kmax_opt, linear_k)

      REAL, INTENT(IN) :: x(:, :)
      REAL, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(IN) :: nk
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nbin(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: sig(:)
      REAL, OPTIONAL, INTENT(IN) :: kmin_opt
      REAL, OPTIONAL, INTENT(IN) :: kmax_opt
      LOGICAL, OPTIONAL, INTENT(IN) :: linear_k
      COMPLEX :: dk(m, m, m)
      REAL :: kmin, kmax, kmin_def, kmax_def
      INTEGER :: n

      ! Default minimum and maximum wavenumbers
      kmin_def = twopi/L
      kmax_def = pi*m/L
      kmin = default_or_optional(kmin_def, kmin_opt)
      kmax = default_or_optional(kmax_def, kmax_opt)

      IF (size(x, 1) .NE. 3) STOP 'POWER_SPECTRUM_PARTICLES: Error, particles must be in 3D'
      n = size(x, 2)

      ! Convert the particle positions into a density field
      ! This assumes all particles have the same mass
      CALL sharp_Fourier_density_contrast(x, n, L, dk, m)

      ! Compute the power spectrum from the density field  
      CALL compute_power_spectrum(dk, dk, m, L, kmin, kmax, nk, k, Pk, nbin, sig, linear_k)

   END SUBROUTINE power_spectrum_particles

   SUBROUTINE folded_power_spectrum_particles(x, L, m, nk, k, Pk, nbin, sig, nfold_opt)

      ! TODO: Can it be ensured that k overlap between folds somehow?
      ! TODO: Combine measurements from different folds to get a reasonable total
      REAL, INTENT(IN) :: x(:, :)
      REAL, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(IN) :: nk
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nbin(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: sig(:, :)
      INTEGER, OPTIONAL, INTENT(IN) :: nfold_opt
      REAL, ALLOCATABLE :: xfold(:, :)
      REAL, ALLOCATABLE :: kfold(:), Pkfold(:), sigfold(:)
      REAL :: Lfold
      INTEGER, ALLOCATABLE :: nbinfold(:)
      INTEGER :: ifold
      INTEGER :: nfold
      INTEGER, PARAMETER :: nfold_def = 0
      LOGICAL, PARAMETER :: linear_k = .TRUE.

      ! Number of foldings
      nfold = default_or_optional(nfold_def, nfold_opt)
      IF (nfold < 0) STOP 'FOLDED_POWER_SPECTRUM_PARTICLES: Error, number of foldings must be zero or positive'

      ! Make a copy of the particle array to fold
      ! TODO: Expensive in memory
      xfold = x
      Lfold = L

      ! Allocate arrays
      ALLOCATE(k(nk, nfold+1), Pk(nk, nfold+1), nbin(nk, nfold+1), sig(nk, nfold+1))

      ! Loop over the folds
      DO ifold = 0, nfold

         ! Make a fold
         IF (ifold > 0) THEN
            CALL fold_particles(xfold, Lfold)
         END IF

         ! Compute power at this level of fold and correct for folding
         CALL power_spectrum_particles(xfold, Lfold, m, nk, kfold, Pkfold, nbinfold, sigfold)
         Pkfold = Pkfold*(L/Lfold)**3
         sigfold = sigfold*(L/Lfold)**3

         ! Fill arrays
         k(:, ifold+1) = kfold
         Pk(:, ifold+1) = Pkfold
         nbin(:, ifold+1) = nbinfold
         sig(:, ifold+1) = sigfold

      END DO

   END SUBROUTINE folded_power_spectrum_particles

   SUBROUTINE cross_spectrum_particles(x1, x2, L, m, nk, k, Pk, nbin, sig, kmin_opt, kmax_opt)

      USE FFT
      REAL, INTENT(IN) :: x1(:, :)
      REAL, INTENT(IN) :: x2(:, :)
      REAL, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(IN) :: nk
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nbin(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: sig(:)
      REAL, OPTIONAL, INTENT(IN) :: kmin_opt
      REAL, OPTIONAL, INTENT(IN) :: kmax_opt
      COMPLEX :: dk1(m, m, m), dk2(m, m, m)
      REAL :: kmin, kmax, kmin_def, kmax_def
      INTEGER :: n1, n2
      
      ! Default minimum and maximum wavenumbers
      kmin_def = twopi/L
      kmax_def = pi*m/L

      ! Check that array sizes make sense
      IF ((size(x1, 1) .NE. 3) .OR. (size(x2, 1) .NE. 3)) STOP 'CROSS_SPECTRUM_PARTICLES: Error, must be 3D'
      n1 = size(x1, 2)
      n2 = size(x2, 2)

      ! Convert the particle positions into a density field
      ! This assumes all particles have the same weight
      CALL sharp_Fourier_density_contrast(x1, n1, L, dk1, m)
      CALL sharp_Fourier_density_contrast(x2, n2, L, dk2, m)

      ! Compute the power spectrum from the density field     
      kmin = default_or_optional(kmin_def, kmin_opt)
      kmax = default_or_optional(kmax_def, kmax_opt)
      CALL compute_power_spectrum(dk1, dk2, m, L, kmin, kmax, nk, k, Pk, nbin, sig)

   END SUBROUTINE cross_spectrum_particles

   SUBROUTINE sharp_Fourier_density_contrast(x, n, L, dk, m)

      ! Bin particles and create the Fourier modes with appropriate sharpening for the binning strategy
      ! TODO: Upgrade to support both m and mn, m arrays
      USE FFT
      INTEGER, INTENT(IN) :: n, m
      REAL, INTENT(IN) :: x(3, n), L
      COMPLEX, INTENT(OUT) :: dk(m, m, m)
      REAL :: w(n), dbar
      REAL :: d(m, m, m)
      COMPLEX :: dk_out(m, m, m)
      INTEGER, PARAMETER :: ibin = 2 ! 2 - CIC binning

      ! Bin the particles with equal weight to create the particle-number field in cells [dimensionless]
      w = 1.
      CALL particle_bin(x, n, L, w, d, m, ibin, all=.TRUE., periodic=.TRUE., verbose=.TRUE.)

      ! Now compute the mean particle number density in a cell [dimensionless]
      dbar = real(n)/real(m)**3

      ! Convert the particle-number field to the density-contrast field [dimensionless]
      d = d/dbar

      ! Do the Fourier Transform, no normalisation
      dk = d
      CALL fft3(dk, dk_out, m, m, m, -1)
      dk = dk_out

      ! Sharpen the Fourier Transform for the binning
      CALL sharpen_k(dk, ibin)

   END SUBROUTINE sharp_Fourier_density_contrast

   SUBROUTINE write_density_slice_ascii(x, n, z1, z2, L, m, outfile)

      ! Write out a slice of density field
      INTEGER, INTENT(IN) :: n   ! Number of particles
      REAL, INTENT(IN) :: x(3, n) ! Particle positions  
      REAL, INTENT(IN) :: z1     ! Starting point for slab
      REAL, INTENT(IN) :: z2     ! Starting point for slab
      REAL, INTENT(IN) :: L      ! Simulation box size
      INTEGER, INTENT(IN) :: m   ! Mesh size for output
      CHARACTER(len=*), INTENT(IN) :: outfile ! Output file name
      REAL :: d(m, m)
      REAL :: xp, yp, smoothing
      INTEGER :: i, j

      ! Make the projected 2D density
      CALL make_projected_density(x, n, z1, z2, L, d, m)

      ! Smooth the density field
      smoothing = L/real(m) !Set the smoothing scale to be the mesh size
      CALL smooth(d, m, smoothing, L)

      ! Write out to file
      WRITE (*, *) 'WRITE_DENSITY_SLICE_ASCII: Writing density map'
      WRITE (*, *) 'WRITE_DENSITY_SLICE_ASCII: File: ', trim(outfile)
      OPEN (9, file=outfile)
      DO i = 1, m
         DO j = 1, m
            xp = cell_position(i, L, m)
            yp = cell_position(j, L, m)
            WRITE (9, *) xp, yp, d(i, j)
         END DO
      END DO
      CLOSE (9)

      WRITE (*, *) 'WRITE_DENSITY_SLICE_ASCII: Done'
      WRITE (*, *)

   END SUBROUTINE write_density_slice_ascii

   SUBROUTINE make_projected_density(x, n, z1, z2, L, d, m)

      ! Write out a slice of density field
      INTEGER, INTENT(IN) :: m ! Mesh
      INTEGER, INTENT(IN) :: n ! Number of particles
      REAL, INTENT(IN) :: x(3, n) ! Particle positions   
      REAL, INTENT(IN) :: z1 ! Starting z for slab
      REAL, INTENT(IN) :: z2 ! Finishing z for slab
      REAL, INTENT(IN) :: L  ! Length of simulation
      REAL, INTENT(OUT) :: d(m, m) ! Projected density field   
      REAL :: vfac, dbar
      REAL, ALLOCATABLE :: x2D(:, :), w2D(:)
      INTEGER :: i, i2D, n2D
      INTEGER, PARAMETER :: ibin = 2 ! 2 - CIC binning

      ! Count the number of particles falling into the slice
      n2D = 0
      DO i = 1, n
         IF (x(3, i) >= z1 .AND. x(3, i) <= z2) n2D = n2D+1
      END DO
      ALLOCATE (x2D(2, n2D))

      ! Make the 2D position array
      i2D = 0
      DO i = 1, n
         IF (x(3, i) >= z1 .AND. x(3, i) <= z2) THEN
            i2D = i2D+1
            x2D(1, i2D) = x(1, i)
            x2D(2, i2D) = x(2, i)
         END IF
      END DO

      ! Bin for the 2D density field and convert to relative density
      ALLOCATE (w2D(n2D))
      w2D = 1.
      CALL particle_bin(x2D, n2D, L, w2D, d, m, ibin, all=.TRUE., periodic=.TRUE., verbose=.TRUE.)
      DEALLOCATE (x2D)
      vfac = (z2-z1)/L
      dbar = (real(n)*vfac)/real(m**2)
      d = d/dbar

      ! Write useful things to screen
      WRITE (*, *) 'MAKE_PROJECTED_DENSITY: Thickness in z [Mpc/h]:', z2-z1
      WRITE (*, *) 'MAKE_PROJECTED_DENSITY: Volume factor:', vfac
      WRITE (*, *) 'MAKE_PROJECTED_DENSITY: Mean cell particle density:', dbar
      WRITE (*, *) 'MAKE_PROJECTED_DENSITY: Done'
      WRITE (*, *)

   END SUBROUTINE make_projected_density

   SUBROUTINE Zeldovich_ICs(x, v, L, k_tab, Pk_tab, vfac, m)

      ! Generate Zeldovich displacement fields and move particles
      REAL, INTENT(INOUT) :: x(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:, :)
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: k_tab(:)
      REAL, INTENT(IN) :: Pk_tab(:)
      REAL, INTENT(IN) :: vfac
      INTEGER, INTENT(IN) :: m 
      REAL :: ips, maxf
      INTEGER :: n
      REAL, ALLOCATABLE :: f(:, :, :, :)
      INTEGER, PARAMETER :: dim = 3

      ! Check array dimensions
      IF ((size(x, 1) .NE. dim)) STOP 'ZELDOVICH_ICS: Error, x must be three dimensional'

      ! Make the displacement field
      CALL generate_displacement_fields(f, m, L, k_tab, Pk_tab)

      ! Calculate some useful things
      n = size(x, 2)
      ips = L/real(n**(1./3.)) ! Mean ID inter-particle spacing
      maxf = maxval(f)         ! Maximum value of 1D displacement

      ! Calculate the particle velocities first
      ALLOCATE(v(3, n))
      CALL Zeldovich_velocity(x, v, vfac*f, L)

      ! Then do the particle positions
      CALL Zeldovich_displacement(x, f, L)

      ! Write some useful things to the screen
      WRITE (*, *) 'ZELDOVICH_ICS: Max 1D displacement [Mpc/h]:', maxf
      WRITE (*, *) 'ZELDOVICH_ICS: Inter-particle spacing [Mpc/h]:', ips
      WRITE (*, *) 'ZELDOVICH_ICS: Max 1D displacement in units of the IPS:', maxf/ips
      WRITE (*, *) 'ZELDOVICH_ICS: Done'
      WRITE (*, *)

   END SUBROUTINE Zeldovich_ICs

   SUBROUTINE Zeldovich_displacement(x, s, L)

      ! Displace particles using the Zeldovich approximation given a displacement field
      ! TODO: Could do CIC reverse interpolation here
      REAL, INTENT(INOUT) :: x(:, :)    ! Initial and final particle positions [Mpc/h]
      REAL, INTENT(IN) :: s(:, :, :, :) ! Dispalcement field [Mpc/h]
      REAL, INTENT(IN) :: L             ! Box size [Mpc/h]
      INTEGER :: i, j, ix(3)
      INTEGER :: n, m
      INTEGER, PARAMETER :: dim = 3

      ! Check fields are 3D
      IF (size(x, 1) .NE. dim) STOP 'ZELDOVICH_DISPLACEMENTS: Error, only works in three dimensions'
      IF (size(s, 1) .NE. dim) STOP 'ZELDOVICH_DISPLACEMENTS: Error, only works in three dimensions'

      ! Total number of particles
      n = size(x, 2)
      m = size(s, 2)
      IF ((m .NE. size(s, 3)) .OR. (m .NE. size(s, 4))) STOP 'ZELDOVICH_DISPLACEMENTS: Error, displacement field must be square'

      WRITE (*, *) 'ZELDOVICH_DISPLACEMENT: Displacing particles'

      ! Loop over all particles
      DO i = 1, n
         DO j = 1, 3
            ix(j) = NGP_cell(x(j, i), L, m) ! Find the integer-coordinates for which cell you are in
         END DO
         x(:, i) = x(:, i)+s(:, ix(1), ix(2), ix(3)) ! Do the displacement
      END DO

      WRITE (*, *) 'ZELDOVICH_DISPLACEMENT: Done'
      WRITE (*, *)

      ! Replace particles that may have strayed
      CALL replace(x, L)

   END SUBROUTINE Zeldovich_displacement

   SUBROUTINE Zeldovich_velocity(x, v, s, L)

      ! Give particles velocities from a velocity field
      ! TODO: Could do CIC reverse interpolation here
      REAL, INTENT(IN) :: x(:, :)               ! Particle position array [Mpc/h]
      REAL, ALLOCATABLE, INTENT(OUT) :: v(:, :) ! Particle velocity array [km/s]
      REAL, INTENT(IN) :: s(:, :, :, :)         ! Velocity field [km/s]
      REAL, INTENT(IN) :: L                     ! Box size [Mpc/h]
      INTEGER :: i, j, ix(3)
      INTEGER :: n, m
      INTEGER, PARAMETER :: dim = 3

      ! Check array dimensions
      IF ((dim .NE. size(x, 1)) .OR. (dim .NE. size(s, 1))) STOP 'ZELDOVICH_VELOCITY: This only works in three dimensions'
      IF ((size(s, 2) .NE. size(s, 3)) .OR. (size(s, 2) .NE. size(s, 4))) STOP 'ZELDOVICH_DISPLACEMENTS: Error, displacement field must be square'

      ! Array sizes
      n = size(x, 2)
      m = size(s, 2)

      ! Write to screen
      WRITE (*, *) 'ZELDOVICH_VELOCITY: Assigining particle velocties'
      WRITE (*, *) 'ZELDOVICH_VELOCITY: Number of particles:', n
      WRITE (*, *) 'ZELDOVICH_VELOCITY: Mesh size:', m
      WRITE (*, *) 'ZELDOVICH_VELOCITY: Box size [Mpc/h]:', L

      ! Allocate array for particle velocities
      ALLOCATE(v(dim, n))

      ! Loop over all particles
      DO i = 1, n
         DO j = 1, 3
            ix(j) = NGP_cell(x(j, i), L, m) ! Find the integer-coordinates for the cell
         END DO
         v(:, i) = s(:, ix(1), ix(2), ix(3)) ! Assign the velocity to the particles from the field
      END DO

      WRITE (*, *) 'ZELDOVICH_VELOCITY: Done'
      WRITE (*, *)

   END SUBROUTINE Zeldovich_velocity

   SUBROUTINE generate_randoms(x, L)

      ! Generate random x, y, z positions in a cube of size L^3
      USE random_numbers
      REAL, INTENT(OUT) :: x(:, :)
      REAL, INTENT(IN) :: L   
      INTEGER :: i, j
      INTEGER :: n, d
      REAL :: dx
      REAL, PARAMETER :: eps = 1e-6

      ! To prevent the particle being exactly at L
      dx = eps*L

      ! Write to screen
      n = size(x, 2)
      d = size(x, 1)
      WRITE (*, *) 'GENERATE_RANDOMS: Generating a uniform-random particle distribution'
      WRITE (*, *) 'GENERATE_RANDOMS: Number of particles:', n
      WRITE (*, *) 'GENERATE_RANDOMS: Number of dimensions:', d

      ! Loop over all particles and coordinates and assign randomly
      DO i = 1, n
         DO j = 1, d
            x(j, i) = random_uniform(0., L-dx)
         END DO
      END DO

      ! Write to screen
      WRITE (*, *) 'GENERATE_RANDOMS: Done'
      WRITE (*, *)

   END SUBROUTINE generate_randoms

   SUBROUTINE generate_grid(x, L)

      ! Generate a grid of positions in a cube of size L^3
      ! These positions will be the locations of cube centres that tessalate the volume
      REAL, INTENT(OUT) :: x(:, :)
      REAL, INTENT(IN) :: L    
      INTEGER :: ix, iy, iz, i, m, n
      INTEGER, PARAMETER :: dim = 3

      IF (size(x, 1) .NE. dim) STOP 'GENERATE_GRID: Error, this only works for a three-dimensional grid'

      ! Check that the particle number is cubic
      n = size(x, 2)
      m = nint(n**(1./3.))
      IF (m**3 .NE. n) STOP 'GENERATE_GRID: Error, you need a cubic number of particles for a grid'

      WRITE (*, *) 'GENERATE_GRID: Generating a grid particle distribution'
      WRITE (*, *) 'GENERATE_GRID: Number of particles (should be a cube number):', n
      WRITE (*, *) 'GENERATE_GRID: Mesh size:', m
      WRITE (*, *) 'GENERATE_GRID: Mesh size cubed (should eqaul number of particles):', m**3
      WRITE (*, *) 'GENERATE_GRID: Box size [Mpc/h]:', L

      ! Loop over all particles
      i = 0 ! Set the particle counting variable to zero
      DO iz = 1, m
         DO iy = 1, m
            DO ix = 1, m
               i = i+1 !Increment the particle counter
               x(1, i) = cell_position(ix, L, m)
               x(2, i) = cell_position(iy, L, m)
               x(3, i) = cell_position(iz, L, m)
            END DO
         END DO
      END DO

      ! Replace in case any end up on boundary
      CALL replace(x, L)

      WRITE (*, *) 'GENERATE_GRID: Done'
      WRITE (*, *)

   END SUBROUTINE generate_grid

   SUBROUTINE generate_poor_glass(x, L)

      ! Generate a poor man's glass in a cube of size L^3
      USE random_numbers
      REAL, INTENT(OUT) :: x(:, :) 
      REAL, INTENT(IN) :: L
      INTEGER :: i, j, m, n
      REAL :: dx
      INTEGER, PARAMETER :: dim = 3

      IF (size(x, 1) .NE. dim) STOP 'GENERATE_POOR_GLASS: Error, this only works in three dimensions'

      ! First genrate a grid
      CALL generate_grid(x, L)

      ! The cube-root of the number of particles
      n = size(x, 2)
      m = nint(n**(1./3.))

      ! How far can the particles be shifted in x, y, z
      ! They need to stay in their initial cube region
      dx = L/real(m)
      dx = dx/2.

      ! Write to screen
      WRITE (*, *) 'GENERATE_POOR_GLASS: Generating a poor man glass from the grid'
      WRITE (*, *) 'GENERATE_POOR_GLASS: Number of particles:', n
      WRITE (*, *) 'GENERATE_POOR_GLASS: Cube root of number of particles:', m
      WRITE (*, *) 'GENERATE_POOR_GLASS: Box size [Mpc/h]:', L
      WRITE (*, *) 'GENERATE_POOR_GLASS: Cell size [Mpc/h]', 2.*dx
      WRITE (*, *) 'GENERATE_POOR_GLASS: Maximum displacement [Mpc/h]', dx

      ! Loop over the particles and do the displacement
      DO i = 1, n
         DO j = 1, dim
            x(j, i) = x(j, i)+random_uniform(-dx, dx)
         END DO
      END DO

      ! Write to screen
      WRITE (*, *) 'GENERATE_POOR_GLASS: Done'
      WRITE (*, *)

   END SUBROUTINE generate_poor_glass

   SUBROUTINE sparse_sample(x, v, n, f)

      USE random_numbers
      REAL, ALLOCATABLE, INTENT(INOUT) :: x(:, :)
      REAL, ALLOCATABLE, INTENT(INOUT) :: v(:, :)
      REAL, INTENT(IN) :: f
      INTEGER, INTENT(INOUT) :: n
      REAL :: x_old(3, n), v_old(3, n)
      INTEGER :: keep(n), n_old, j, i

      n_old = n

      keep = 0
      n = 0

      STOP 'SPARSE_SAMPLE: Test this'

      WRITE (*, *) 'SPARSE_SAMPLE: Sparse sampling'

      DO i = 1, n_old

         IF (random_unit() < f) THEN
            n = n+1
            keep(i) = 1
         END IF

      END DO

      x_old = x
      v_old = v

      DEALLOCATE (x, v)
      ALLOCATE (x(3, n), v(3, n))

      j = 0

      DO i = 1, n_old

         IF (keep(i) == 1) THEN
            j = j+1
            x(1, j) = x_old(1, i)
            x(2, j) = x_old(2, i)
            x(3, j) = x_old(3, i)
            v(1, j) = v_old(1, i)
            v(2, j) = v_old(2, i)
            v(3, j) = v_old(3, i)
         END IF

      END DO

      WRITE (*, *) 'SPARSE_SAMPLE: Complete'
      WRITE (*, *) 'SPARSE_SAMPLE: Before:', n_old
      WRITE (*, *) 'SPARSE_SAMPLE: After:', n
      WRITE (*, *) 'SPARSE_SAMPLE: Ratio:', real(n)/real(n_old)
      WRITE (*, *)

   END SUBROUTINE sparse_sample

   SUBROUTINE replace(x, L, verbose)

      ! Ensures/enforces periodicity by cycling particles round that may have strayed
      ! This forces all particles to be 0<=x<L, so they cannot be exactly at x=L
      REAL, INTENT(INOUT) :: x(:, :)    
      REAL, INTENT(IN) :: L
      LOGICAL, OPTIONAL :: verbose
      INTEGER :: i, j
      INTEGER :: m

      IF (present_and_correct(verbose)) WRITE (*, *) 'REPLACE: Replacing particles that may have strayed'

      ! Loop over all particles and coordinates
      DO i = 1, size(x, 2)
         DO j = 1, size(x, 1)
            m = FLOOR(x(j, i)/L)
            x(j, i) = x(j, i)-m*L
            IF (x(j, i) == -0.) x(j, i) = 0.
         END DO
      END DO

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'REPLACE: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE replace

   SUBROUTINE random_translation(x, n, L, verbose)

      USE random_numbers

      INTEGER, INTENT(IN) :: n
      REAL, INTENT(INOUT) :: x(3, n)    
      REAL, INTENT(IN) :: L
      LOGICAL, OPTIONAL :: verbose
      INTEGER :: i
      REAL :: T(3)

      IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_TRANSLATION: Applying random translation'

      ! Calculate and apply the translations
      DO i = 1, 3
         T(i) = random_uniform(0., L)
         IF (present_and_correct(verbose)) THEN
            WRITE (*, *) 'RANDOM_TRANSLATION: Direction', i
            WRITE (*, *) 'RANDOM_TRANSLATION: Translation [Mpc/h]:', T(i)
         END IF
         x(i, :) = x(i, :)+T(i)
      END DO

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'RANDOM_TRANSLATION: Done'
         WRITE (*, *)
      END IF

      ! Replace particles within the cube
      CALL replace(x, L, verbose)

   END SUBROUTINE random_translation

   SUBROUTINE random_inversion(x, n, L, verbose)

      USE random_numbers

      INTEGER, INTENT(IN) :: n
      REAL, INTENT(INOUT) :: x(3, n)
      REAL, INTENT(IN) :: L
      LOGICAL, OPTIONAL :: verbose
      INTEGER :: i
      INTEGER :: pm(3)

      IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_INVERSION: Applying random inversion'

      ! Calculate and apply the translations
      DO i = 1, 3
         pm(i) = random_sign()
         IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_INVERSION: Direction', i, 'Sign:', pm(i)
         x(i, :) = x(i, :)*pm(i)
      END DO

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'RANDOM_INVERSION: Done'
         WRITE (*, *)
      END IF

      ! Replace particles within the cube
      CALL replace(x, L, verbose)

   END SUBROUTINE random_inversion

   SUBROUTINE random_rotation(x, n, verbose)

      ! Do a random rotation of the box by 90 degrees, I think this also captures inversions
      ! This way is memory efficient
      ! TODO: This raises an 'array temporary' warning because non-contiguous sections of x are passed to swap_arrays
      USE random_numbers

      INTEGER, INTENT(IN) :: n
      REAL, INTENT(INOUT) :: x(3, n)  
      LOGICAL, OPTIONAL :: verbose
      INTEGER :: type

      IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_ROTATION: Applying random rotation'

      ! Choose random rotation
      type = random_integer(1, 6)

      ! Apply random rotation
      IF (type == 1) THEN
         IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_ROTATION: Doing 123'
         ! Do nothing here
      ELSE IF (type == 2) THEN
         IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_ROTATION: Doing 132'
         CALL swap_arrays(x(2, :), x(3, :)) ! Swap 2<->3 (123 -> 132)
      ELSE IF (type == 3) THEN
         IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_ROTATION: Doing 213'
         CALL swap_arrays(x(1, :), x(2, :)) ! Swap 1<->2 (123 -> 213)
      ELSE IF (type == 4) THEN
         IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_ROTATION: Doing 231'
         CALL swap_arrays(x(1, :), x(2, :)) ! Swap 1<->2 (123 -> 213)
         CALL swap_arrays(x(2, :), x(3, :)) ! Swap 2<->3 (213 -> 231)
      ELSE IF (type == 5) THEN
         IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_ROTATION: Doing 312'
         CALL swap_arrays(x(2, :), x(3, :)) ! Swap 1<->2 (123 -> 132)
         CALL swap_arrays(x(1, :), x(2, :)) ! Swap 2<->3 (132 -> 312)
      ELSE IF (type == 6) THEN
         IF (present_and_correct(verbose)) WRITE (*, *) 'RANDOM_ROTATION: Doing 321'
         CALL swap_arrays(x(1, :), x(3, :)) ! Swap 1<->3 (123 -> 321)
      ELSE
         STOP 'RANDOM_ROTATION: Error, something went very wrong'
      END IF

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'RANDOM_ROTATION: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE random_rotation

   SUBROUTINE particle_bin_2D(x, n, L, w, d, m, ibin, all, periodic, verbose)

      ! Bin particle properties onto a mesh, summing as you go
      ! TODO: Adapt for different lengths and different meshes in x,y
      INTEGER, INTENT(IN) :: n         ! Total number of particles
      INTEGER, INTENT(IN) :: m         ! Mesh size for density field
      REAL, INTENT(IN) :: x(2, n)       ! 2D particle positions    
      REAL, INTENT(IN) :: L            ! Area side length
      REAL, INTENT(IN) :: w(n)         ! Weight array for each particle
      REAL, INTENT(OUT) :: d(m, m)      ! Output of eventual 2D density field
      INTEGER, INTENT(IN) :: ibin      ! Binning strategy
      LOGICAL, INTENT(IN) :: all       ! Should all particles be contributing to the binning?
      LOGICAL, INTENT(IN) :: periodic  ! Is the volume periodic?
      LOGICAL, INTENT(IN) :: verbose   ! Verbose?

      IF (periodic .AND. (.NOT. all)) STOP 'PARTICLE_BIN_2D: Very strange to have periodic and not all particles contribute'

      IF (ibin == 1) THEN
         CALL NGP_2D(x, n, L, w, d, m, all, verbose)
      ELSE IF (ibin == 2) THEN
         CALL CIC_2D(x, n, L, w, d, m, all, periodic, verbose)
      ELSE
         STOP 'PARTICLE_BIN_2D: Error, ibin not specified correctly'
      END IF

   END SUBROUTINE particle_bin_2D

   SUBROUTINE particle_bin_3D(x, n, L, w, d, m, ibin, all, periodic, verbose)

      ! Bin particle properties onto a mesh, summing as you go
      INTEGER, INTENT(IN) :: n        ! Total number of particles
      INTEGER, INTENT(IN) :: m        ! Mesh size for density field
      REAL, INTENT(IN) :: x(3, n)      ! 3D particle positions  
      REAL, INTENT(IN) :: L           ! Volume side length
      REAL, INTENT(IN) :: w(n)        ! Weight array for each particle
      REAL, INTENT(INOUT) :: d(m, m, m) ! Output of eventual 3D density field
      INTEGER, INTENT(IN) :: ibin     ! Binning strategy
      LOGICAL, INTENT(IN) :: all      ! Should all particles be contributing to the binning?
      LOGICAL, INTENT(IN) :: periodic ! Is the volume periodic?
      LOGICAL, INTENT(IN) :: verbose  ! Verbose?

      IF (periodic .AND. (.NOT. all)) STOP 'PARTICLE_BIN_3D: Very strange to have periodic and not all particles contribute'

      IF (ibin == 1) THEN
         CALL NGP_3D(x, n, L, w, d, m, all, verbose)
      ELSE IF (ibin == 2) THEN
         CALL CIC_3D(x, n, L, w, d, m, all, periodic, verbose)
      ELSE
         STOP 'PARTICLE_BIN_3D: Error, ibin not specified correctly'
      END IF

   END SUBROUTINE particle_bin_3D

   SUBROUTINE particle_bin_average_3D(x, n, L, w, d, m, ibin, all, periodic)

      ! Bin particle properties onto a mesh, averaging properties over cells
      ! TODO: This should probably not be used if there are any empty cells
      INTEGER, INTENT(IN) :: n, m
      REAL, INTENT(IN) :: x(3, n)    
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: w(n)
      REAL, INTENT(OUT) :: d(m, m, m)  
      INTEGER, INTENT(IN) :: ibin
      LOGICAL, INTENT(IN) :: all
      LOGICAL, INTENT(IN) :: periodic
      REAL :: number(m, m, m), ones(n)
      INTEGER :: i, j, k, sum

      IF (periodic .AND. (.NOT. all)) STOP 'PARTICLE_BIN_3D: Very strange to have periodic and not all particles contribute'

      ! Need an array of ones
      ! TODO: wasteful for memory, could use pointer maybe?
      ones = 1.

      !Call the binning twice, first to bin particle property and second to count
      IF (ibin == 1) THEN
         CALL NGP_3D(x, n, L, w, d, m, all, verbose=.TRUE.)
         CALL NGP_3D(x, n, L, ones, number, m, all, verbose=.TRUE.)
      ELSE IF (ibin == 2) THEN
         CALL CIC_3D(x, n, L, w, d, m, all, periodic, verbose=.TRUE.)
         CALL CIC_3D(x, n, L, ones, number, m, all, periodic, verbose=.TRUE.)
      ELSE
         STOP 'PARTICLE_BIN_AVERAGE: Error, ibin not specified correctly'
      END IF

      ! Now want to average over all elements in each cell, but need to be wary of cells with zero particles
      ! Loop over all elements of the field and average over the number of contributions
      ! Not sure how this will work for CIC binning if cell only gets a small contribution from one particle
      sum = 0
      WRITE (*, *) 'PARTICLE_BIN_AVERAGE: Averaging over particles in cells'
      DO k = 1, m
         DO j = 1, m
            DO i = 1, m
               IF (number(i, j, k) == 0.) THEN
                  sum = sum+1
               ELSE
                  d(i, j, k) = d(i, j, k)/number(i, j, k)
               END IF
            END DO
         END DO
      END DO
      WRITE (*, *) 'PARTICLE_BIN_AVERAGE: Numer of empty cells:', n
      WRITE (*, *) 'PARTICLE_BIN_AVERAGE: Fraction of empty cells', real(sum)/real(m**3)
      WRITE (*, *) 'PARTICLE_BIN_AVERAGE: Done'
      WRITE (*, *)

   END SUBROUTINE particle_bin_average_3D

   SUBROUTINE NGP_2D(x, n, L, w, d, m, all, verbose)

      ! Nearest-grid-point binning routine
      INTEGER, INTENT(IN) :: n       ! Total number of particles in area
      INTEGER, INTENT(IN) :: m       ! Mesh size for density field
      REAL, INTENT(IN) :: x(2, n)    ! particle positions     
      REAL, INTENT(IN) :: L          ! Area side length
      REAL, INTENT(IN) :: w(n)       ! Weight array
      REAL, INTENT(OUT) :: d(m, m)   ! Output of eventual 2D density field  
      LOGICAL, INTENT(IN) :: all     ! Should all particles be contributing to the binning?
      LOGICAL, INTENT(IN) :: verbose ! Verbose
      INTEGER :: i, j, ix(2)
      LOGICAL :: outside
      INTEGER, PARAMETER :: dim = 2

      IF (verbose) THEN
         WRITE (*, *) 'NGP_2D: Binning particles and creating field'
         WRITE (*, *) 'NGP_2D: Binning region size:', L
         WRITE (*, *) 'NGP_2D: Cells:', m
      END IF

      ! Fix to zero
      d = 0.

      DO i = 1, n

         ! Get the integer coordinates of the cell
         DO j = 1, dim
            ix(j) = NGP_cell(x(j, i), L, m)
         END DO

         ! Check if particles are outside the mesh region
         outside = .FALSE.
         DO j = 1, dim
            IF (ix(j) < 1 .OR. ix(j) > m) outside = .TRUE.
         END DO

         ! If the particle is outside the mesh region then either make an error or cycle
         IF (outside) THEN
            IF (all) THEN
               DO j = 1, dim
                  WRITE (*, *) 'NGP_2D: Coordinate:', j
                  WRITE (*, *) 'NGP_2D: position:', x(j, i)
                  WRITE (*, *) 'NGP_2D: cell: ', ix(j)
               END DO
               STOP 'NGP_2D: Error, particle outside boundary'
            ELSE
               CYCLE
            END IF
         END IF

         ! Do the binning
         d(ix(1), ix(2)) = d(ix(1), ix(2))+w(i)

      END DO

      IF (verbose) THEN
         WRITE (*, *) 'NGP_2D: Binning complete'
         WRITE (*, *)
      END IF

   END SUBROUTINE NGP_2D

   SUBROUTINE NGP_3D(x, n, L, w, d, m, all, verbose)

      ! Nearest-grid-point binning routine
      INTEGER, INTENT(IN) :: n, m     ! Total number of particles in area
      REAL, INTENT(IN) :: x(3, n)     ! particle positions    
      REAL, INTENT(IN) :: L           ! Area side length
      REAL, INTENT(IN) :: w(n)        ! Weight array
      REAL, INTENT(OUT) :: d(m, m, m) ! Output of eventual 3D density field    
      LOGICAL, INTENT(IN) :: all      ! Should all particles be contributing to the binning?
      LOGICAL, INTENT(IN) :: verbose  ! Verbose
      INTEGER :: i, j, ix(3)
      LOGICAL :: outside
      INTEGER, PARAMETER :: dim = 3

      IF (verbose) THEN
         WRITE (*, *) 'NGP_3D: Binning particles and creating field'
         WRITE (*, *) 'NGP_3D: Binning region size:', L
         WRITE (*, *) 'NGP_3D: Cells:', m
      END IF

      ! Fix to zero
      d = 0.

      DO i = 1, n

         ! Get the integer coordinates of the cell
         DO j = 1, dim
            ix(j) = NGP_cell(x(j, i), L, m)
         END DO

         ! Check if particles are outside the mesh region
         outside = .FALSE.
         DO j = 1, dim
            IF (ix(j) < 1 .OR. ix(j) > m) outside = .TRUE.
         END DO

         ! If the particle is outside the mesh region then either make an error or cycle
         IF (outside) THEN
            IF (all) THEN
               DO j = 1, dim
                  WRITE (*, *) 'NGP_3D: Coordinate:', j
                  WRITE (*, *) 'NGP_3D: position:', x(j, i)
                  WRITE (*, *) 'NGP_3D: cell: ', ix(j)
               END DO
               STOP 'NGP_3D: Error, particle outside boundary'
            ELSE
               CYCLE
            END IF
         END IF

         ! Do the binning
         d(ix(1), ix(2), ix(3)) = d(ix(1), ix(2), ix(3))+w(i)

      END DO

      IF (verbose) THEN
         WRITE (*, *) 'NGP_3D: Binning complete'
         WRITE (*, *)
      END IF

   END SUBROUTINE NGP_3D

   SUBROUTINE CIC_2D(x, n, L, w, d, m, all, periodic, verbose)

      ! Cloud-in-cell binning routine
      INTEGER, INTENT(IN) :: n, m     ! Total number of particles in area
      REAL, INTENT(IN) :: x(2, n)     ! 2D particle positions     
      REAL, INTENT(IN) :: L           ! Area side length
      REAL, INTENT(IN) :: w(n)        ! Weight array
      REAL, INTENT(OUT) :: d(m, m)    ! Output of eventual 2D density field
      LOGICAL, INTENT(IN) :: all      ! Should all particles be contributing to the binning?
      LOGICAL, INTENT(IN) :: periodic ! Is the volume periodic?
      LOGICAL, INTENT(IN) :: verbose  ! Verbose
      INTEGER :: i, j, k
      INTEGER :: ix(2), iy(2), ic(2)
      REAL :: dx(2), eps
      LOGICAL :: outside
      INTEGER, PARAMETER :: dim = 2

      IF (verbose) THEN
         WRITE (*, *) 'CIC_2D: Binning particles and creating density field'
         WRITE (*, *) 'CIC_2D: Binning region size:', L
         WRITE (*, *) 'CIC_2D: Cells:', m
      END IF

      ! Needs to be set to zero
      eps = 0.
      d = 0.

      DO i = 1, n

         ! Get the cell interger coordinates
         DO j = 1, dim
            ix(j) = NGP_cell(x(j, i), L, m)
         END DO

         ! Check to make sure particles are all within the binning region
         IF (all) THEN

            ! See if particles are outside region
            outside = .FALSE.
            DO j = 1, dim
               IF (ix(j) < 1 .OR. ix(j) > m) outside = .TRUE.
            END DO

            ! Write to screen if they are
            IF (outside) THEN
               DO j = 1, dim
                  WRITE (*, *) 'CIC_2D: Coordinate', j
                  WRITE (*, *) 'CIC_2D: Position:', x(j, i)
                  WRITE (*, *) 'CIC_2D: Cell: ', ix(j)
               END DO
               STOP 'CIC_2D: Error, particle outside boundary'
            END IF

         END IF

         ! dx, dy in cell units, away from cell centre
         DO j = 1, dim
            dx(j) = (x(j, i)/L)*real(m)-(real(ix(j))-0.5)
         END DO

         ! Find which other cell needs a contribution
         DO j = 1, dim
            IF (dx(j) >= 0.) THEN
               iy(j) = ix(j)+1
            ELSE
               iy(j) = ix(j)-1
               dx(j) = -dx(j) ! So that dx is always positive
            END IF
         END DO

         ! Deal with periodicity
         IF (periodic) THEN
            DO j = 1, dim
               IF (iy(j) == m+1) THEN
                  iy(j) = 1
               ELSE IF (iy(j) == 0) THEN
                  iy(j) = m
               END IF
            END DO
         END IF

         ! Carry out CIC binning
         DO k = 1, 4
            IF (k == 1) THEN
               ! Main cell
               eps = (1.-dx(1))*(1.-dx(2))
               ic(1) = ix(1)
               ic(2) = ix(2)
            ELSE IF (k == 2) THEN
               ! Offset in x
               eps = dx(1)*(1.-dx(2))
               ic(1) = iy(1)
               ic(2) = ix(2)
            ELSE IF (k == 3) THEN
               ! Offset in y
               eps = (1.-dx(1))*dx(2)
               ic(1) = ix(1)
               ic(2) = iy(2)
            ELSE IF (k == 4) THEN
               ! Offset in xy
               eps = dx(1)*dx(2)
               ic(1) = iy(1)
               ic(2) = iy(2)
            ELSE
               STOP 'CIC_2D: Error, something went terribly wrong'
            END IF
            CALL add_to_array(d, eps*w(i), ic)
         END DO

      END DO

      IF (verbose) THEN
         WRITE (*, *) 'CIC_2D: Binning complete'
         WRITE (*, *)
      END IF

   END SUBROUTINE CIC_2D

   SUBROUTINE CIC_3D(x, n, L, w, d, m, all, periodic, verbose)

      ! Cloud-in-cell binning routine
      INTEGER, INTENT(IN) :: n, m     ! Total number of particles in area
      REAL, INTENT(IN) :: x(3, n)     ! 3D particle positions    
      REAL, INTENT(IN) :: L           ! Area side length
      REAL, INTENT(IN) :: w(n)        ! Weight array
      REAL, INTENT(OUT) :: d(m, m, m) ! Output of eventual 2D density field
      LOGICAL, INTENT(IN) :: all      ! Should all particles be contributing to the binning?
      LOGICAL, INTENT(IN) :: periodic ! Is the volume periodic?
      LOGICAL, INTENT(IN) :: verbose  ! Verbose
      INTEGER :: i, j, k
      INTEGER :: ix(3), iy(3), ic(3)
      REAL :: dx(3), eps
      LOGICAL :: outside
      INTEGER, PARAMETER :: dim = 3

      IF (verbose) THEN
         WRITE (*, *) 'CIC_3D: Binning particles and creating density field'
         WRITE (*, *) 'CIC_3D: Number of particles:', n
         WRITE (*, *) 'CIC_3D: Binning region size:', L
         WRITE (*, *) 'CIC_3D: Cells:', m
      END IF

      ! Needs to be set to zero
      eps = 0.
      d = 0.

      DO i = 1, n

         ! Get the cell interger coordinates
         DO j = 1, dim
            ix(j) = NGP_cell(x(j, i), L, m)
         END DO

         ! Check to make sure particles are all within the binning region
         IF (all) THEN

            ! See if particles are outside region
            outside = .FALSE.
            DO j = 1, dim
               IF (ix(j) < 1 .OR. ix(j) > m) outside = .TRUE.
            END DO

            ! Write to screen if they are
            IF (outside) THEN
               DO j = 1, dim
                  WRITE (*, *) 'CIC_3D: Coordinate', j
                  WRITE (*, *) 'CIC_3D: Position:', x(j, i)
                  WRITE (*, *) 'CIC_3D: Cell: ', ix(j)
               END DO
               STOP 'CIC_3D: Error, particle outside boundary'
            END IF

         END IF

         ! dx, dy in cell units, away from cell centre
         DO j = 1, dim
            dx(j) = (x(j, i)/L)*real(m)-(real(ix(j))-0.5)
         END DO

         ! Find which other cell needs a contribution
         DO j = 1, dim
            IF (dx(j) >= 0.) THEN
               iy(j) = ix(j)+1
            ELSE
               iy(j) = ix(j)-1
               dx(j) = -dx(j) ! So that dx is always positive
            END IF
         END DO

         ! Deal with periodicity
         IF (periodic) THEN
            DO j = 1, dim
               IF (iy(j) == m+1) THEN
                  iy(j) = 1
               ELSE IF (iy(j) == 0) THEN
                  iy(j) = m
               END IF
            END DO
         END IF

         ! Carry out CIC binning
         DO k = 1, 8
            IF (k == 1) THEN
               ! Main cell
               eps = (1.-dx(1))*(1.-dx(2))*(1-dx(3))
               ic(1) = ix(1)
               ic(2) = ix(2)
               ic(3) = ix(3)
            ELSE IF (k == 2) THEN
               ! Offset in x
               eps = dx(1)*(1.-dx(2))*(1-dx(3))
               ic(1) = iy(1)
               ic(2) = ix(2)
               ic(3) = ix(3)
            ELSE IF (k == 3) THEN
               ! Offset in y
               eps = (1.-dx(1))*dx(2)*(1-dx(3))
               ic(1) = ix(1)
               ic(2) = iy(2)
               ic(3) = ix(3)
            ELSE IF (k == 4) THEN
               ! Offset in z
               eps = (1.-dx(1))*(1.-dx(2))*dx(3)
               ic(1) = ix(1)
               ic(2) = ix(2)
               ic(3) = iy(3)
            ELSE IF (k == 5) THEN
               ! Offset in xy
               eps = dx(1)*dx(2)*(1-dx(3))
               ic(1) = iy(1)
               ic(2) = iy(2)
               ic(3) = ix(3)
            ELSE IF (k == 6) THEN
               ! Offset in yz
               eps = (1.-dx(1))*dx(2)*dx(3)
               ic(1) = ix(1)
               ic(2) = iy(2)
               ic(3) = iy(3)
            ELSE IF (k == 7) THEN
               ! Offset in xz
               eps = dx(1)*(1.-dx(2))*dx(3)
               ic(1) = iy(1)
               ic(2) = ix(2)
               ic(3) = iy(3)
            ELSE IF (k == 8) THEN
               ! Offset in xyz
               eps = dx(1)*dx(2)*dx(3)
               ic(1) = iy(1)
               ic(2) = iy(2)
               ic(3) = iy(3)
            ELSE
               STOP 'CIC_3D: Error, something went terribly wrong'
            END IF
            CALL add_to_array(d, eps*w(i), ic)
         END DO

      END DO

      IF (verbose) THEN
         WRITE (*, *) 'CIC_3D: Binning complete'
         WRITE (*, *)
      END IF

   END SUBROUTINE CIC_3D

   SUBROUTINE SOD(x, n, L, w, d, m, Rs)

      ! Spherical density routine
      INTEGER, INTENT(IN) :: n, m
      REAL, INTENT(IN) :: x(3, n)   
      REAL, INTENT(IN) :: L
      REAL, INTENT(IN) :: w(n)
      REAL, INTENT(OUT) :: d(m, m, m)
      REAL, INTENT(IN) :: Rs
      !LOGICAL, INTENT(IN) :: all
      !LOGICAL, INTENT(IN) :: periodic
      REAL :: xc(3, m, m, m)
      REAL :: r, dx, dy, dz
      INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz, kx, ky, kz

      WRITE (*, *) 'SOD: Binning particles and creating density field'
      WRITE (*, *) 'SOD: Cells:', m

      d = 0.

      ! Fill array with sphere centre positions
      DO k = 1, m
         DO j = 1, m
            DO i = 1, m
               xc(1, i, j, k) = cell_position(i, L, m)
               xc(2, i, j, k) = cell_position(j, L, m)
               xc(3, i, j, k) = cell_position(k, L, m)
            END DO
         END DO
      END DO

      ! Loop over all particles and assign to spheres
      DO i = 1, n

         ! Coordinate of nearest mesh cell
         ix = NGP_cell(x(1, i), L, m)
         iy = NGP_cell(x(2, i), L, m)
         iz = NGP_cell(x(3, i), L, m)

         ! See if particle is within spheres of over 27 neighbouring spheres
         DO jx = -1, 1
            DO jy = -1, 1
               DO jz = -1, 1

                  ! Find the x integer for the sphere
                  kx = ix+jx
                  IF (kx < 1) THEN
                     kx = kx+m
                  ELSE IF (kx > m) THEN
                     kx = kx-m
                  END IF

                  ! Find the y integer for the sphere
                  ky = iy+jy
                  IF (ky < 1) THEN
                     ky = ky+m
                  ELSE IF (ky > m) THEN
                     ky = ky-m
                  END IF

                  ! Find the z integer for the sphere
                  kz = iz+jz
                  IF (kz < 1) THEN
                     kz = kz+m
                  ELSE IF (kz > m) THEN
                     kz = kz-m
                  END IF

                  !WRITE(*,*) 'Mesh:', jx, jy, jz, kx, ky, kz

                  ! Calculate the displacement between the particle and the sphere centre
                  dx = x(1, i)-xc(1, kx, ky, kz)
                  dy = x(2, i)-xc(2, kx, ky, kz)
                  dz = x(3, i)-xc(3, kx, ky, kz)
                  !!

                  !WRITE(*,*) 'Displacements:', dx, dy, dz

                  !!
                  !Correct if particle is on the other side of the box - x
                  IF (dx > L/2.) THEN
                     dx = dx-L
                  ELSE IF (dx < -L/2.) THEN
                     dx = dx+L
                  END IF

                  !Correct if particle is on the other side of the box - y
                  IF (dy > L/2.) THEN
                     dy = dy-L
                  ELSE IF (dy < -L/2.) THEN
                     dy = dy+L
                  END IF

                  !Correct if particle is on the other side of the box - z
                  IF (dz > L/2.) THEN
                     dz = dz-L
                  ELSE IF (dz < -L/2.) THEN
                     dz = dz+L
                  END IF
                  !!

                  !WRITE(*,*) 'Displacements:', dx, dy, dz

                  !!
                  !Calculate the radial distance between the sphere centre and the particle
                  r = sqrt(dx**2.+dy**2.+dz**2.)

                  !Add to density if within sphere radius
                  IF (r < Rs) THEN
                     d(kx, ky, kz) = d(kx, ky, kz)+w(i)
                     !WRITE(*,*) 'Added:', kx, ky, kz
                  END IF
                  !!

               END DO
            END DO
         END DO
         !STOP

      END DO

      WRITE (*, *) 'SOD: Binning complete'
      WRITE (*, *)

   END SUBROUTINE SOD

   SUBROUTINE find_pairs(x, okay, n, rmin, rmax, L, outfile)

      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: x(3, n)
      LOGICAL, INTENT(IN) :: okay(n)
      REAL, INTENT(IN) :: rmin, rmax, L
      CHARACTER(len=*), INTENT(IN) :: outfile
      INTEGER :: i, j, np
      REAL :: r

      WRITE (*, *) 'FIND_PAIRS: Minimum separation [Mpc/h]:', rmin
      WRITE (*, *) 'FIND_PAIRS: Maximum separation [Mpc/h]:', rmax
      WRITE (*, *) 'FIND_PAIRS: Doing finding'

      !Fix the number counter to be zero
      np = 0

      OPEN (7, file=outfile)
      DO i = 1, n
         DO j = i+1, n
            IF (okay(i) .AND. okay(j)) THEN
               r = periodic_distance(x(:, i), x(:, j), L)
               IF (r >= rmin .AND. r <= rmax) THEN
                  !CALL add_to_array(x(:,i),x(:,j),pairs,np)
                  np = np+1
                  WRITE (7, *) x(1, i), x(2, i), x(3, i), x(1, j), x(2, j), x(3, j)
               END IF
            END IF
         END DO
      END DO
      CLOSE (7)

      WRITE (*, *) 'FIND_PAIRS: Number of pairs:', np
      WRITE (*, *) 'FIND_PAIRS: Done'
      WRITE (*, *)

   END SUBROUTINE find_pairs

   SUBROUTINE zshift(x, v, n, Om_m, Om_v, z, iz)

      ! Shift particles to redshift space
      INTEGER, INTENT(IN) :: n       ! Total number of particles
      REAL, INTENT(INOUT) :: x(3, n) ! Particle positions [Mpc/h]
      REAL, INTENT(IN) :: v(3, n)    ! Particle velocities [km/s] 
      REAL, INTENT(IN) :: Om_m, Om_v ! Cosmological parameters
      REAL, INTENT(IN) :: z          ! Redshift
      INTEGER, INTENT(IN) :: iz      ! Dimension along which to shift (1,2,3 for x,y,z)
      INTEGER :: i
      REAL :: H, a

      WRITE (*, *) 'ZSHIFT: Shifting particles into redshift space'

      H = 100.*sqrt(Hubble2_simple(z, Om_m, Om_v))

      a = 1./(1.+z)

      DO i = 1, n
         x(iz, i) = x(iz, i)+v(iz, i)/(a*H)
      END DO

      WRITE (*, *) 'ZSHIFT: Done'
      WRITE (*, *)

   END SUBROUTINE zshift

   SUBROUTINE fold_particles(x, L)

      ! Fold (really translate) particles by half in each dimension
      REAL, INTENT(INOUT) :: x(:, :) ! Particles
      REAL, INTENT(INOUT) :: L       ! Length of periodic box
      INTEGER :: n, dim
      INTEGER :: i, d

      ! Pre-calculations
      n = size(x, 2)   ! Number of particles
      dim = size(x, 1) ! Number of dimensions
      L = L/2.         ! Half box size
      
      ! Loop and do 'folding'
      DO i = 1, n
         DO d = 1, dim
            IF (x(d, i) >= L) x(d, i) = x(d, i)-L
         END DO
      END DO

   END SUBROUTINE fold_particles

   REAL FUNCTION Hubble2_simple(z, Om_m, Om_v)

      ! This calculates the dimensionless squared hubble parameter squared at redshift z!
      ! It ignores contributions from radiation (not accurate at very high z)!
      ! It also ignores anything other than vacuum and matter
      ! TODO: Include dark energy?
      ! TODO: Remove in favour of cosmology module?
      REAL, INTENT(IN) :: z ! Redshift
      REAL, INTENT(IN) :: Om_m, Om_v ! Cosmological parameters

      Hubble2_simple = Om_m*(1.+z)**3+Om_v+(1.-Om_m-Om_v)*(1.+z)**2

   END FUNCTION Hubble2_simple

   ! This subroutine was previously called slice
   SUBROUTINE write_slice_ascii(x, n, x1, x2, y1, y2, z1, z2, outfile)

      INTEGER, INTENT(IN) :: n                   ! Total number of particles
      REAL, INTENT(IN) :: x(3, n)                ! Particle positions [Mpc/h]
      REAL, INTENT(IN) :: x1, x2, y1, y2, z1, z2 ! Limits of the slice [Mpc/h]
      CHARACTER(len=*), INTENT(IN) :: outfile    ! Output file
      INTEGER :: i

      WRITE (*, *) 'WRITE_SLICE_ASCII: Writing slice'
      OPEN (10, file=outfile)
      WRITE (*, *) 'WRITE_SLICE_ASCII: Thickness in z [Mpc/h]:', (z2-z1)
      DO i = 1, n
         IF (x1 < x(1, i) .AND. x(1, i) <= x2 .AND. y1 < x(2, i) .AND. x(2, i) <= y2 .AND. z1 < x(3, i) .AND. x(3, i) <= z2) THEN
            WRITE (10, *) x(1, i), x(2, i), x(3, i)
         END IF
      END DO
      CLOSE (10)
      WRITE (*, *) 'WRITE_SLICE_ASCII: Slice written: ', trim(outfile)
      WRITE (*, *)

   END SUBROUTINE write_slice_ascii

   REAL FUNCTION shot_noise_simple(L, n)

      ! Calculate simulation shot noise for equal mass matter particlers [(Mpc/h)^3]
      USE precision
      REAL, INTENT(IN) :: L ! Box size [Mpc/h]
      !INTEGER*8, INTENT(IN) :: n ! Total number of particles
      INTEGER(int8), INTENT(IN) :: n ! Total number of particles

      ! Calculate number density
      shot_noise_simple = L**3/real(n)

   END FUNCTION shot_noise_simple

   REAL FUNCTION shot_noise(u, v, n, L)

      ! Calculates shot noise in P_uv(k) [(Mpc/h)^3]
      ! See appendix in HMx paper
      ! Note if this is for mass then u = v = particle_mass/total_mass
      ! CARE: Removed double-precision sums
      INTEGER, INTENT(IN) :: n ! Number of particles
      REAL, INTENT(IN) :: u(n) ! Contributions to the total field u and v per particle
      REAL, INTENT(IN) :: v(n) ! Contributions to the total field u and v per particle
      REAL, INTENT(IN) :: L    ! Box size [Mpc/h]

      ! Do the sum of the two fields, making sure to use doubles
      shot_noise = sum(u*v)

      ! Multiply through by factors of volume and mesh
      shot_noise = (L**3)*shot_noise

   END FUNCTION shot_noise

   REAL FUNCTION shot_noise_mass(L, m, n)

      ! Calculate simulation shot noise constant P(k) for different-mass tracers [(Mpc/h)^3]
      ! CARE: Removed double-precision sums
      INTEGER, INTENT(IN) :: n ! Number of particles
      REAL, INTENT(IN) :: L    ! Box size [Mpc/h]
      REAL, INTENT(IN) :: m(n) ! Array of particle masses
      REAL :: Nbar

      ! Calculate the effective mean number of tracers
      Nbar = sum(m)**2/sum(m**2)

      ! Calculate number density, this makes units of P(k) [(Mpc/h)^3]
      shot_noise_mass = L**3/Nbar

   END FUNCTION shot_noise_mass

   ELEMENTAL REAL FUNCTION shot_noise_k(k, shot)

      ! Calculates shot noise as Delta^2(k) from a constant-P(k) thing with units [(Mpc/h)^3]
      ! This shot noise is then dimensionless, exactly like Delta^2(k)
      USE cosmology_functions
      REAL, INTENT(IN) :: k    ! Wave vector [h/Mpc]
      REAL, INTENT(IN) :: shot ! The constant shot-noise term: P(k) [(Mpc/h)^3]

      shot_noise_k = Delta_Pk(shot, k)

   END FUNCTION shot_noise_k

   !SUBROUTINE adaptive_density(xc,yc,Lsub,z1,z2,m,r,x,w,n,L,nbar,outfile)
   SUBROUTINE write_adaptive_field(xc, yc, Lsub, z1, z2, m, r, x, w, n, L, outfile)

      ! Makes a pretty picture of a density field but uses adaptive meshes to make it nice
      USE string_operations
      REAL, INTENT(IN) :: xc, yc ! Coordinates of image centre [Mpc/h]
      REAL, INTENT(IN) :: Lsub ! Size of image [Mpc/h]
      REAL, INTENT(IN) :: z1, z2 ! Front and back of image [Mpc/h]
      INTEGER, INTENT(IN) :: m ! Image size in pixels (e.g., 2048)
      INTEGER, INTENT(IN) :: r ! Integer number of levels of refinement
      INTEGER, INTENT(IN) :: n ! Number of particles
      REAL, INTENT(IN) :: x(3, n) ! Particle position array [Mpc/h]
      REAL, INTENT(IN) :: w(n) ! Weight array (e.g., mass, or just an array of ones)
      REAL, INTENT(IN) :: L ! Box size [Mpc/h]
      CHARACTER(len=*), INTENT(IN) :: outfile ! Output file

      REAL :: x1, x2, y1, y2, Lx, Ly, Lz, delta, npexp
      INTEGER :: np
      INTEGER :: m1, m2, m3, m4, m5
      INTEGER :: i, j
      REAL, ALLOCATABLE :: y(:, :), d(:, :), u(:), ones(:)
      REAL, ALLOCATABLE :: count1(:, :), count2(:, :), count3(:, :), count4(:, :), count5(:, :)
      REAL, ALLOCATABLE :: field1(:, :), field2(:, :), field3(:, :), field4(:, :), field5(:, :)
      CHARACTER(len=256) :: base, ext, output
      LOGICAL :: all, periodic

      REAL, PARAMETER :: dc = 4. ! Refinement conditions (particles-per-cell)
      REAL, PARAMETER :: fcell = 1. ! Smoothing factor over cell sizes (maybe should be set to 1; 0.75 looks okay)
      LOGICAL, PARAMETER :: test = .FALSE. ! Activate test mode
      INTEGER, PARAMETER :: ibin = 2 ! 2- CIC binnings

      IF (r > 4) STOP 'WRITE_ADAPTIVE_FIELD: Error, too many refinement leves requested. Maximum is 4'

      ! Write information to screen
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Adaptive density field generator'
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Refinement conditions (particles-per-cell):', dc
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Smoothing factor in cell size:', fcell
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Number of refinements:', r
      WRITE (*, *)

      ! Set the region boundaries in Mpc/h
      x1 = xc-Lsub/2.
      x2 = xc+Lsub/2.
      y1 = yc-Lsub/2.
      y2 = yc+Lsub/2.

      ! Set the region thicknesses
      Lx = x2-x1
      Ly = y2-y1
      Lz = z2-z1

      ! Write information to screen
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Subvolume:'
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: x1:', x1
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: x2:', x2
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Lx:', Lx
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: y1:', y1
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: y2:', y2
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Ly:', Ly
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: z1:', z1
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: z2:', z2
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Lz:', Lz
      WRITE (*, *)

      ! First pass to count the number of particles in the region
      np = 0
      DO i = 1, n
         IF (x(1, i) > x1 .AND. x(1, i) < x2 .AND. x(2, i) > y1 .AND. x(2, i) < y2 .AND. x(3, i) > z1 .AND. x(3, i) < z2) THEN
            np = np+1
         END IF
      END DO

      ! Calculate particle-density statistics
      npexp = real(n)*(Lx*Ly*Lz/L**3.)
      delta = -1.+real(np)/real(npexp)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Particle-density statistics'
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Particles in subvolume:', np
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Expected number of partilces in subvolume:', npexp
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Particle over-density in subvolume:', delta
      WRITE (*, *)

      ! Allocate arrays for 2D position and 2D weight
      ALLOCATE (y(2, np), u(np), ones(np))
      ones = 1. ! Need an array of ones for particle binning

      ! Second pass to add particles in the region to 2D array y
      j = 0
      DO i = 1, n
         IF (x(1, i) > x1 .AND. x(1, i) < x2 .AND. x(2, i) > y1 .AND. x(2, i) < y2 .AND. x(3, i) > z1 .AND. x(3, i) < z2) THEN
            j = j+1
            y(1, j) = x(1, i)-xc+Lsub/2. ! 2D position array
            y(2, j) = x(2, i)-yc+Lsub/2. ! 2D position array
            u(j) = w(i) ! u is the 2D weight array
         END IF
      END DO

      ! Boost the weight array to account for the smaller volume
      ! Note that this is a bit complicated
      ! The weight array for overdensity is m_i/M, where M is the total box mass
      ! i.e., the weight is the contribution of each particle to total box mass
      ! For pressure it is the pressure contribution to the total box pressure
      ! Need to account for the fact that we are actually working in a subvolume now
      u = u*(L**3/(Lx*Ly*Lz))

      ! Set the sizes of all of the adaptive meshes, which all differ by a factor of 2
      ! m1=m, m2=m/2, m3=m/4, m4=m/8, m5=m/16
      m1 = m
      m2 = m1/2
      m3 = m2/2
      m4 = m3/2
      m5 = m4/2

      ! Write out sizes to screen
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh size 1:', m1
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh size 2:', m2
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh size 3:', m3
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh size 4:', m4
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh size 5:', m5
      WRITE (*, *)

      ! Allocate arrays for particle counts and fields
      ! This is ugly, but I do not think it can be avoided with e.g., field(5,m,m) because m are all different
      ALLOCATE (field1(m1, m1), count1(m1, m1))
      ALLOCATE (field2(m2, m2), count2(m2, m2))
      ALLOCATE (field3(m3, m3), count3(m3, m3))
      ALLOCATE (field4(m4, m4), count4(m4, m4))
      ALLOCATE (field5(m5, m5), count5(m5, m5))

      all = .TRUE. ! All particles in array y contribute to the binning
      IF (Lsub == L) THEN
         periodic = .TRUE. ! Only possibly periodic if the subvolume size is the same as the actual size
      ELSE
         periodic = .FALSE.
      END IF

      ! Do binning of particles on each mesh resolution
      CALL particle_bin(y, np, Lsub, ones, count1, m1, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, ones, count2, m2, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, ones, count3, m3, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, ones, count4, m4, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, ones, count5, m5, ibin, all, periodic, verbose=.TRUE.)

      ! Now bin field values, rather than particle numbers
      ! Need to multiply the weights through by factors of mesh to give the contribution to the mesh cell
      ! For example, for overdensity u is m_i/M, where M is now the subvolume mass
      ! This changes it to be the contribution to the density per mesh cell
      ! Same for pressure contributions
      CALL particle_bin(y, np, Lsub, u*m1**2, field1, m1, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, u*m2**2, field2, m2, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, u*m3**2, field3, m3, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, u*m4**2, field4, m4, ibin, all, periodic, verbose=.TRUE.)
      CALL particle_bin(y, np, Lsub, u*m5**2, field5, m5, ibin, all, periodic, verbose=.TRUE.)

      ! Smooth fields
      CALL smooth(field1, m1, fcell*Lsub/real(m1), Lsub)
      CALL smooth(field2, m2, fcell*Lsub/real(m2), Lsub)
      CALL smooth(field3, m3, fcell*Lsub/real(m3), Lsub)
      CALL smooth(field4, m4, fcell*Lsub/real(m4), Lsub)
      CALL smooth(field5, m5, fcell*Lsub/real(m5), Lsub)

      ! Write out density statistics on each mesh
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh:', m1
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max count:', maxval(count1)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min count:', minval(count1)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max overdensity:', maxval(field1)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min overdensity:', minval(field1)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh:', m2
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max count:', maxval(count2)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min count:', minval(count2)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max overdensity:', maxval(field2)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min overdensity:', minval(field2)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh:', m3
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max count:', maxval(count3)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min count:', minval(count3)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max overdensity:', maxval(field3)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min overdensity:', minval(field3)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh:', m4
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max count:', maxval(count4)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min count:', minval(count4)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max overdensity:', maxval(field4)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min overdensity:', minval(field4)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Mesh:', m5
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max count:', maxval(count5)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min count:', minval(count5)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max overdensity:', maxval(field5)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min overdensity:', minval(field5)
      WRITE (*, *)

      ! Write out each mesh if testing
      IF (test) THEN

         base = 'density_'
         ext = '.dat'

         output = number_file_zeroes(base, m1, 4, ext)
         WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Writing: ', trim(output)
         CALL write_field_ascii(field1, m1, Lsub, output)

         output = number_file_zeroes(base, m2, 4, ext)
         WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Writing: ', trim(output)
         CALL write_field_ascii(field2, m2, Lsub, output)

         output = number_file_zeroes(base, m3, 4, ext)
         WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Writing: ', trim(output)
         CALL write_field_ascii(field3, m3, Lsub, output)

         output = number_file_zeroes(base, m4, 4, ext)
         WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Writing: ', trim(output)
         CALL write_field_ascii(field4, m4, Lsub, output)

         output = number_file_zeroes(base, m5, 4, ext)
         WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Writing: ', trim(output)
         CALL write_field_ascii(field5, m5, Lsub, output)

      END IF

      ! No refinement required
      IF (r == 0) THEN

         ALLOCATE (d(m, m))
         d = field1

         ! Start the refinements
      ELSE

         ! First refinement
         IF (r >= 4) THEN

            ! Allocate d with the size of the second-corsest image
            ALLOCATE (d(m4, m4))
            DO i = 1, m5
               DO j = 1, m5
                  IF (count5(i, j) < dc) THEN
                     ! Use coarser density
                     d(2*i-1, 2*j-1) = field5(i, j)
                     d(2*i, 2*j-1) = field5(i, j)
                     d(2*i-1, 2*j) = field5(i, j)
                     d(2*i, 2*j) = field5(i, j)
                  ELSE
                     ! Use finer density
                     d(2*i-1, 2*j-1) = field4(2*i-1, 2*j)
                     d(2*i, 2*j-1) = field4(2*i, 2*j-1)
                     d(2*i-1, 2*j) = field4(2*i-1, 2*j)
                     d(2*i, 2*j) = field4(2*i, 2*j)
                  END IF
               END DO
            END DO

         END IF

         ! Second refinement
         IF (r >= 3) THEN

            IF (ALLOCATED(d)) THEN
               field4 = d
               DEALLOCATE (d)
            END IF
            ALLOCATE (d(m3, m3))
            DO i = 1, m4
               DO j = 1, m4
                  IF (count4(i, j) < dc) THEN
                     ! Use coarser density
                     d(2*i-1, 2*j-1) = field4(i, j)
                     d(2*i, 2*j-1) = field4(i, j)
                     d(2*i-1, 2*j) = field4(i, j)
                     d(2*i, 2*j) = field4(i, j)
                  ELSE
                     ! Use finer density
                     d(2*i-1, 2*j-1) = field3(2*i-1, 2*j)
                     d(2*i, 2*j-1) = field3(2*i, 2*j-1)
                     d(2*i-1, 2*j) = field3(2*i-1, 2*j)
                     d(2*i, 2*j) = field3(2*i, 2*j)
                  END IF
               END DO
            END DO

         END IF

         ! Third refinement
         IF (r >= 2) THEN

            IF (ALLOCATED(d)) THEN
               field3 = d
               DEALLOCATE (d)
            END IF
            ALLOCATE (d(m2, m2))
            DO i = 1, m3
               DO j = 1, m3
                  IF (count3(i, j) < dc) THEN
                     ! Use coarser density
                     d(2*i-1, 2*j-1) = field3(i, j)
                     d(2*i, 2*j-1) = field3(i, j)
                     d(2*i-1, 2*j) = field3(i, j)
                     d(2*i, 2*j) = field3(i, j)
                  ELSE
                     ! Use finer density
                     d(2*i-1, 2*j-1) = field2(2*i-1, 2*j)
                     d(2*i, 2*j-1) = field2(2*i, 2*j-1)
                     d(2*i-1, 2*j) = field2(2*i-1, 2*j)
                     d(2*i, 2*j) = field2(2*i, 2*j)
                  END IF
               END DO
            END DO

         END IF

         ! Fourth refinement
         IF (r >= 1) THEN

            IF (ALLOCATED(d)) THEN
               field2 = d
               DEALLOCATE (d)
            END IF
            ALLOCATE (d(m1, m1))
            DO i = 1, m2
               DO j = 1, m2
                  IF (count2(i, j) < dc) THEN
                     ! Use coarser density
                     d(2*i-1, 2*j-1) = field2(i, j)
                     d(2*i, 2*j-1) = field2(i, j)
                     d(2*i-1, 2*j) = field2(i, j)
                     d(2*i, 2*j) = field2(i, j)
                  ELSE
                     ! Use finer density
                     d(2*i-1, 2*j-1) = field1(2*i-1, 2*j)
                     d(2*i, 2*j-1) = field1(2*i, 2*j-1)
                     d(2*i-1, 2*j) = field1(2*i-1, 2*j)
                     d(2*i, 2*j) = field1(2*i, 2*j)
                  END IF
               END DO
            END DO

         END IF

      END IF

      ! Deallocate arrays
      DEALLOCATE (field1, field2, field3, field4, field5)
      DEALLOCATE (count1, count2, count3, count4, count5)

      ! Write info to screen
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Adaptive density field'
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max density:', maxval(d)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min density:', minval(d)
      WRITE (*, *)

      ! Ensure all cells are positive (note this is to make pretty pictures, not for science)
      DO i = 1, m1
         DO j = 1, m1
            IF (d(i, j) < 0.) d(i, j) = 0.
         END DO
      END DO

      ! Write field statistics
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Smoothed adaptive density field'
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Max density:', maxval(d)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Min density:', minval(d)
      WRITE (*, *)

      ! Finally write out an ascii file for plotting
      CALL write_field_ascii(d, m1, Lsub, outfile)
      WRITE (*, *) 'WRITE_ADAPTIVE_FIELD: Done'
      WRITE (*, *)

   END SUBROUTINE write_adaptive_field

   SUBROUTINE make_HOD(xh, nph, rv, rs, irho, L, x, np)

      ! Make an HOD realisation
      ! TODO: Could speed up for NFW by not regenerating the iversion/interpolation halo unless rv, rs change
      REAL, ALLOCATABLE, INTENT(IN) :: xh(:, :) ! Halo position array [Mpc/h]
      INTEGER, INTENT(IN) :: nph(:)             ! Number of particles in each halo
      REAL, INTENT(IN) :: rv(:)                 ! Halo virial radii [Mpc/h]
      REAL, INTENT(IN) :: rs(:)                 ! Halo scale radii [Mpc/h]
      INTEGER, INTENT(IN) :: irho               ! Halo profile specifier
      REAL, INTENT(IN) :: L                     ! Simulation box size [Mpc/h]
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :) ! HOD particle position array [Mpc/h]
      INTEGER, INTENT(OUT) :: np                ! Total number of HOD particles generated  
      INTEGER :: i, j, k
      INTEGER :: nh
      REAL, ALLOCATABLE :: xi(:, :)

      nh = size(xh, 2)
      IF((nh .NE. size(nph)) .OR. (nh .NE. size(rv))) THEN
         STOP 'MAKE_HOD: Error, halo-position array, number of particles and virial radii must all be the same size'
      END IF

      ! Calculate total number of particles to generate
      np = SUM(nph)
      ALLOCATE (x(3, np))

      ! Generate the particles
      k = 0
      DO i = 1, nh
         IF (nph(i) > 0) THEN
            ALLOCATE(xi(3, nph(i)))
            CALL random_halo_particles(xi, rv(i), rs(i), irho)
            DO j = 1, nph(i)
               k = k + 1
               x(:, k) = xh(:, i)+xi(:, j)
            END DO
            DEALLOCATE(xi)
         END IF
      END DO

      ! Replace particles that may have strayed outside the volume
      CALL replace(x, L)

   END SUBROUTINE make_HOD

   SUBROUTINE random_halo_particles(rr, rv, rs, irho)

      ! Get a set of x,y,z coordiantes for a random point in an artificial spherical halo
      ! The halo centre is fixed to be the coordinate zero point
      USE random_numbers
      USE HMx
      REAL, INTENT(OUT) :: rr(:, :) ! Halo coordinates (should be shape 3,n) [Mpc/h]
      REAL, INTENT(IN) :: rv        ! Halo virial radius [Mpc/h]
      REAL, INTENT(IN) :: rs        ! Halo scale radius [Mpc/h]
      INTEGER, INTENT(IN) :: irho   ! Switch for halo profile (from HMx)
      REAL, ALLOCATABLE :: r(:)
      REAL :: theta, phi
      INTEGER :: i, n

      IF (size(rr, 1) /= 3) STOP 'RANDOM_HALO_PARTICLES: Error, output array must be shape (3, n)'
      n = size(rr, 2)
      ALLOCATE(r(n))

      IF (irho == irho_NFW) THEN
         CALL random_rs_NFW(r, rv, rs)
      ELSE
         DO i = 1, n
            ! Get radial coordinate
            IF (irho == irho_tophat) THEN
               r(i) = random_r_constant(rv)
            ELSE IF (irho == irho_iso) THEN
               r(i) = random_r_isothermal(rv)
            ELSE IF (irho == irho_shell) THEN
               r(i) = rv
            ELSE IF (irho == irho_delta) THEN
               r(i) = 0.
            ELSE
               STOP 'RANDOM_SPHERE: Error, irho specified incorrectly or not supported'
            END IF
         END DO
      END IF

      ! Generate the full 3D coordinates assuming a spherical profile (isotropic angles)
      DO i = 1, n
         theta = random_spherical_theta()    ! Angle from z: Between 0 and pi
         phi = random_uniform(0., twopi)     ! x-y plane angle
         rr(1, i) = r(i)*sin(theta)*cos(phi) ! x
         rr(2, i) = r(i)*sin(theta)*sin(phi) ! y
         rr(3, i) = r(i)*cos(theta)          ! z
      END DO

   END SUBROUTINE random_halo_particles

   REAL FUNCTION random_r_constant(rv)

      ! The radial random weighting for a constant-density halo profile
      USE random_numbers
      REAL, INTENT(IN) :: rv ! halo virial radius [Mpc/h]

      random_r_constant = rv*random_polynomial(2.)

   END FUNCTION random_r_constant

   REAL FUNCTION random_r_isothermal(rv)

      ! The radial random weighting for an isothermal halo profile
      USE random_numbers
      REAL, INTENT(IN) :: rv ! halo virial radius [Mpc/h]

      random_r_isothermal = random_uniform(0., rv)

   END FUNCTION random_r_isothermal

   SUBROUTINE random_rs_NFW(rr, rv, rs)

      ! Get a set of random radii drawn from an NFW distribution
      ! Has to perform a numerical inversion/interpolation
      ! If F(x) = ln(1+x) - x/(1+x) then f = F(r/rs) / F(c) with f in [0->1] must be inverted for r
      USE random_numbers
      USE interpolate
      USE HMx
      REAL, INTENT(OUT) :: rr(:)
      REAL, INTENT(IN) :: rv
      REAL, INTENT(IN) :: rs
      REAL, ALLOCATABLE :: r(:), f(:)
      REAL :: u
      INTEGER :: i
      INTEGER, PARAMETER :: nr = nr_random_NFW
      INTEGER, PARAMETER :: iorder = iorder_random_NFW
      INTEGER, PARAMETER :: ifind = ifind_random_NFW
      INTEGER, PARAMETER :: iinterp = iinterp_random_NFW

      ! Fill array of r vs. f(r) for interpolation/inversion
      CALL fill_array(0., rv, r, nr)
      ALLOCATE(f(nr))
      DO i = 1, nr
         f(i) = NFW_factor(r(i)/rs)/NFW_factor(rv/rs)
      END DO

      ! Get n different values of the radius using numerical inversion
      DO i = 1, size(rr)
         u = random_unit()
         rr(i) = find(u, f, r, nr, iorder, ifind, iinterp)
      END DO

   END SUBROUTINE random_rs_NFW
   
   SUBROUTINE combine_folded_power(k, Pk, sig, k_bin, k_tot, Pk_tot, sig_tot)

      ! TODO: Should weight enter the calculation of mean k? Or should this just be an array of ones?
      ! TODO: Remove power above half-Nyquist at each fold
      USE statistics
      REAL, INTENT(IN) :: k(:, :)
      REAL, INTENT(IN) :: Pk(:, :)
      REAL, INTENT(IN) :: sig(:, :)
      REAL, INTENT(IN) :: k_bin(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: k_tot(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: Pk_tot(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: sig_tot(:)
      REAL, ALLOCATABLE :: weight(:, :), ones(:)
      REAL, ALLOCATABLE :: w_k(:), w_Pk(:), w_sig(:)
      INTEGER :: nk, nf
      INTEGER :: ifold, ik

      ! Check k dimension
      nk = size(k, 1)
      IF (size(Pk, 1) /= nk) STOP 'COMBINE_FOLDED_POWER_DATA: Error Pk_fold should have the same first dimension as k_fold'
      IF (size(sig, 1) /= nk) STOP 'COMBINE_FOLDED_POWER_DATA: Error sig_fold should have the same first dimension as k_fold'

      ! Check fold dimension
      nf = size(k, 2)
      IF (size(Pk, 2) /= nf) STOP 'COMBINE_FOLDED_POWER_DATA: Error Pk_fold should have the same second dimension as k_fold'
      IF (size(sig, 2) /= nf) STOP 'COMBINE_FOLDED_POWER_DATA: Error sig_fold should have the same second dimension as k_fold'
      nf = nf-1

      ! Create weight array from errors
      ALLOCATE(weight(nk, nf+1))
      DO ifold = 1, nf+1
         DO ik = 1, nk
            IF (sig(ik, ifold) == 0.) THEN
               weight(ik, ifold) = 0.
            ELSE
               weight(ik, ifold) = 1./sig(ik, ifold)
            END IF
         END DO
      END DO
      ALLOCATE(ones(nk*(nf+1)))
      ones = 1.

      ! Histogram k, P(k) and sig(k) 
      CALL advanced_histogram(reshape(k, [nk*(nf+1)]), reshape(k, [nk*(nf+1)]), ones, k_bin, k_tot, w_k)
      CALL advanced_histogram(reshape(k, [nk*(nf+1)]), reshape(Pk, [nk*(nf+1)]), reshape(weight, [nk*(nf+1)]), k_bin, Pk_tot, w_Pk)
      CALL advanced_histogram(reshape(k, [nk*(nf+1)]), reshape(weight**2, [nk*(nf+1)]), ones, k_bin, sig_tot, w_sig)

      ! Normalise results
      DO ik = 1, size(k_tot)
         k_tot(ik) = k_tot(ik)/w_k(ik)
         IF (w_Pk(ik) == 0.) THEN
            Pk_tot(ik) = 0.
            sig_tot(ik) = 0.
         ELSE         
            Pk_tot(ik) = Pk_tot(ik)/w_Pk(ik)
            sig_tot(ik) = sqrt(1./sig_tot(ik))
         END IF
      END DO

   END SUBROUTINE combine_folded_power

END MODULE simulations
