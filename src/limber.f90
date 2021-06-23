MODULE Limber

   USE constants
   USE basic_operations
   USE array_operations
   USE interpolate
   USE cosmology_functions

   IMPLICIT NONE

   PRIVATE

   ! Types
   PUBLIC :: projection

   ! Functions
   PUBLIC :: maxdist
   PUBLIC :: calculate_Cl
   PUBLIC :: xcorr_type
   PUBLIC :: Cl_contribution_ell
   PUBLIC :: Cl_contribution
   PUBLIC :: fill_projection_kernels
   PUBLIC :: read_nz
   PUBLIC :: write_projection_kernel
   PUBLIC :: write_angular_xi
   PUBLIC :: write_nz
   PUBLIC :: write_efficiency
   PUBLIC :: calculate_angular_xi
   PUBLIC :: write_Cl
   PUBLIC :: k_ell
   PUBLIC :: xpow_pka
   PUBLIC :: xpow_halomod
   PUBLIC :: set_field_for_xpow

   ! Tracers
   PUBLIC :: n_tracers
   PUBLIC :: tracer_RCSLenS
   PUBLIC :: tracer_Compton_y
   PUBLIC :: tracer_CMB_lensing
   PUBLIC :: tracer_CFHTLenS_vanWaerbeke2013
   PUBLIC :: tracer_KiDS
   PUBLIC :: tracer_KiDS_bin1
   PUBLIC :: tracer_KiDS_bin2
   PUBLIC :: tracer_KiDS_bin3
   PUBLIC :: tracer_KiDS_bin4
   PUBLIC :: tracer_gravity_wave
   PUBLIC :: tracer_KiDS_450
   PUBLIC :: tracer_KiDS_450_fat_bin1
   PUBLIC :: tracer_KiDS_450_fat_bin2
   PUBLIC :: tracer_KiDS_450_highz
   PUBLIC :: tracer_CIB_353
   PUBLIC :: tracer_CIB_545
   PUBLIC :: tracer_CIB_857
   PUBLIC :: tracer_galaxies
   PUBLIC :: tracer_lensing_z1p00
   PUBLIC :: tracer_lensing_z0p75
   PUBLIC :: tracer_lensing_z0p50
   PUBLIC :: tracer_lensing_z0p25
   PUBLIC :: tracer_KiDS_450_bin1
   PUBLIC :: tracer_KiDS_450_bin2
   PUBLIC :: tracer_KiDS_450_bin3
   PUBLIC :: tracer_KiDS_450_bin4
   PUBLIC :: tracer_CFHTLenS_Kilbinger2013

   ! Projection quantities that need to be calculated only once; these relate to the Limber integrals
   TYPE projection
      REAL, ALLOCATABLE :: X(:), r_X(:)
      REAL, ALLOCATABLE :: q(:), r_q(:)
      REAL, ALLOCATABLE :: nz(:), z_nz(:)
      INTEGER :: nX, nq, nnz, order_nz
   END TYPE projection

   ! P(k,a) look-up table parameters
   REAL, PARAMETER :: kmin_pka = 1e-4   ! k' value for P(k,a) table; P(k<k',a)=0
   REAL, PARAMETER :: kmax_pka = 1e2    ! k' value for P(k,a) table; P(k>k',a)=0
   REAL, PARAMETER :: amin_pka = 1e-3   ! a' value for P(k,a) table; P(k,a<a')=0
   REAL, PARAMETER :: amax_pka = 1.     ! a' value for P(k,a) table; P(k,a>a')=0

   ! Maximum k that corresponding to a given ell (prevents infinities for low z)
   REAL, PARAMETER :: k_ell_max = 1e3

   ! xcorr - C(l) calculation
   LOGICAL, PARAMETER :: verbose_Limber = .FALSE. ! Verbosity

   ! Maxdist
   REAL, PARAMETER :: dr_max = 0.01 ! Small subtraction from maxdist to prevent numerical issues

   ! Limber integration parameters
   INTEGER, PARAMETER  :: n_cont = 1025 ! Number of samples to take in r for Cl ell contribution calculation
   REAL, PARAMETER :: lcorr = 0.5       ! 1/2 in LoVerde (2008) Limber correction k(r)=(l+1/2)/f_k(r)
   REAL, PARAMETER :: acc_Limber = 1e-4 ! Accuracy parameter for Limber integration

   ! n(z)
   REAL, PARAMETER :: zmin_nz = 0.  ! Minimum redshift for the analytic n(z) tables
   REAL, PARAMETER :: zmax_nz = 2.5 ! Maximum redshift for the analytic n(z) tables
   INTEGER, PARAMETER :: n_nz = 129 ! Number of entries in the analytic n(z) tables
   LOGICAL, PARAMETER :: norm_nz = .TRUE. ! Explicitly normalise the input n(z)

   ! Projection kernel options
   REAL, PARAMETER :: rmin_kernel = 0.         ! Minimum r for table
   INTEGER, PARAMETER :: nx_kernel = 129       ! Number of entires in look-up table
   REAL, PARAMETER :: rmin_lensing = 0.        ! Minimum distance in integral
   INTEGER, PARAMETER :: nX_lensing = 129      ! Number of entries in X(r) table
   INTEGER, PARAMETER :: nq_efficiency = 129   ! Number of entries in q(r) table
   !INTEGER, PARAMETER :: iorder_efficiency = 1 ! Order for lensing-efficiency integration over n(z)

   ! Gravitational waves
   REAL, PARAMETER :: A_gwave = 1.
   REAL, PARAMETER :: rmin_gwave = 10.

   ! Method for angular correlation function
   ! 1 - Summation with accuracy parameter (converges very slowly, accuracy condition seems to not work)
   ! 2 - Standard integration (does not converge well, crap)
   ! 3 - Standard summation
   ! 4 - Summation assuming C(l) precomputed (fastest)
   INTEGER, PARAMETER :: method_angular_xi = 4

   ! Angular xi verbosity
   LOGICAL, PARAMETER :: verbose_xi = .FALSE. ! Verbosity

   ! Tracer types
   INTEGER, PARAMETER :: tracer_RCSLenS = 1
   INTEGER, PARAMETER :: tracer_Compton_y = 2
   INTEGER, PARAMETER :: tracer_CMB_lensing = 3
   INTEGER, PARAMETER :: tracer_CFHTLenS_vanWaerbeke2013 = 4
   INTEGER, PARAMETER :: tracer_KiDS = 5
   INTEGER, PARAMETER :: tracer_KiDS_bin1 = 6
   INTEGER, PARAMETER :: tracer_KiDS_bin2 = 7
   INTEGER, PARAMETER :: tracer_KiDS_bin3 = 8
   INTEGER, PARAMETER :: tracer_KiDS_bin4 = 9
   INTEGER, PARAMETER :: tracer_gravity_wave = 10
   INTEGER, PARAMETER :: tracer_KiDS_450 = 11
   INTEGER, PARAMETER :: tracer_KiDS_450_fat_bin1 = 12
   INTEGER, PARAMETER :: tracer_KiDS_450_fat_bin2 = 13
   INTEGER, PARAMETER :: tracer_KiDS_450_highz = 14
   INTEGER, PARAMETER :: tracer_CIB_353 = 15
   INTEGER, PARAMETER :: tracer_CIB_545 = 16
   INTEGER, PARAMETER :: tracer_CIB_857 = 17
   INTEGER, PARAMETER :: tracer_galaxies = 18
   INTEGER, PARAMETER :: tracer_lensing_z1p00 = 19
   INTEGER, PARAMETER :: tracer_lensing_z0p75 = 20
   INTEGER, PARAMETER :: tracer_lensing_z0p50 = 21
   INTEGER, PARAMETER :: tracer_lensing_z0p25 = 22
   INTEGER, PARAMETER :: tracer_KiDS_450_bin1 = 23
   INTEGER, PARAMETER :: tracer_KiDS_450_bin2 = 24
   INTEGER, PARAMETER :: tracer_KiDS_450_bin3 = 25
   INTEGER, PARAMETER :: tracer_KiDS_450_bin4 = 26
   INTEGER, PARAMETER :: tracer_CFHTLenS_Kilbinger2013 = 27
   INTEGER, PARAMETER :: n_tracers = 27

CONTAINS

   CHARACTER(len=256) FUNCTION xcorr_type(ix)

      ! Names for cross-correlation field types
      INTEGER, INTENT(IN) :: ix

      xcorr_type = ''
      IF (ix == tracer_RCSLenS) xcorr_type = 'RCSLenS lensing'
      IF (ix == tracer_Compton_y) xcorr_type = 'Compton y'
      IF (ix == tracer_CMB_lensing) xcorr_type = 'CMB lensing'
      IF (ix == tracer_CFHTLenS_vanWaerbeke2013) xcorr_type = 'CFHTLenS lensing (van Waerbeke 2013)'
      IF (ix == tracer_KiDS) xcorr_type = 'KiDS lensing (z = 0.1 -> 0.9)'
      IF (ix == tracer_KiDS_bin1) xcorr_type = 'KiDS lensing (z = 0.1 -> 0.3)'
      IF (ix == tracer_KiDS_bin2) xcorr_type = 'KiDS lensing (z = 0.3 -> 0.5)'
      IF (ix == tracer_KiDS_bin3) xcorr_type = 'KiDS lensing (z = 0.5 -> 0.7)'
      IF (ix == tracer_KiDS_bin4) xcorr_type = 'KiDS lensing (z = 0.7 -> 0.9)'
      IF (ix == tracer_gravity_wave) xcorr_type = 'Gravitational waves'
      IF (ix == tracer_KiDS_450) xcorr_type = 'KiDS 450 (z = 0.1 -> 0.9)'
      IF (ix == tracer_KiDS_450_fat_bin1) xcorr_type = 'KiDS 450 (z = 0.1 -> 0.5)'
      IF (ix == tracer_KiDS_450_fat_bin2) xcorr_type = 'KiDS 450 (z = 0.5 -> 0.9)'
      IF (ix == tracer_KiDS_450_highz) xcorr_type = 'KiDS 450 (z = 0.9 -> 3.5)'
      IF (ix == tracer_CIB_353) xcorr_type = 'CIB 353 GHz'
      IF (ix == tracer_CIB_545) xcorr_type = 'CIB 545 GHz'
      IF (ix == tracer_CIB_857) xcorr_type = 'CIB 857 GHz'
      IF (ix == tracer_galaxies) xcorr_type = 'Galaxies'
      IF (ix == tracer_lensing_z1p00) xcorr_type = 'Lensing with fixed z=1.00 source plane'
      IF (ix == tracer_lensing_z0p75) xcorr_type = 'Lensing with fixed z=0.75 source plane'
      IF (ix == tracer_lensing_z0p50) xcorr_type = 'Lensing with fixed z=0.50 source plane'
      IF (ix == tracer_lensing_z0p25) xcorr_type = 'Lensing with fixed z=0.25 source plane'
      IF (ix == tracer_KiDS_450_bin1) xcorr_type = 'KiDS 450 (z = 0.1 -> 0.3)'
      IF (ix == tracer_KiDS_450_bin2) xcorr_type = 'KiDS 450 (z = 0.3 -> 0.5)'
      IF (ix == tracer_KiDS_450_bin3) xcorr_type = 'KiDS 450 (z = 0.5 -> 0.7)'
      IF (ix == tracer_KiDS_450_bin4) xcorr_type = 'KiDS 450 (z = 0.7 -> 0.9)'
      IF (ix == tracer_CFHTLenS_Kilbinger2013) xcorr_type = 'CFHTLenS lensing (Kilbinger 2013)'
      IF (xcorr_type == '') STOP 'XCORR_TYPE: Error, ix not specified correctly'

   END FUNCTION xcorr_type

   REAL FUNCTION k_ell(ell, a, cosm)

      ! Finds the k that corresponds to ell at the given a
      REAL, INTENT(IN) :: ell ! angular wave number
      REAL, INTENT(IN) :: a   ! scale factor
      TYPE(cosmology), INTENT(INOUT) :: cosm ! cosmology
      REAL :: r

      IF (a == 1.) THEN
         ! This should really be infinite, but this stops a division by infinity
         k_ell = k_ell_max
      ELSE
         r = comoving_distance(a, cosm)
         k_ell = (ell+lcorr)/f_k(r, cosm)
      END IF

   END FUNCTION k_ell

   SUBROUTINE write_nz(proj, output)

      ! Write out the n(z) to a file
      TYPE(projection), INTENT(IN) :: proj
      CHARACTER(len=*), INTENT(IN) :: output
      INTEGER :: i

      OPEN (7, file=output)
      DO i = 1, proj%nnz
         WRITE (7, *) proj%z_nz(i), proj%nz(i)
      END DO
      CLOSE (7)

   END SUBROUTINE write_nz

   SUBROUTINE fill_projection_kernels(ix, proj, cosm)

      ! Fill look-up tables for the two projection kerels X_ij
      ! TODO: Change to take ix(n)?
      INTEGER, INTENT(IN) :: ix(2) ! Label for the type of projection kernel needed
      TYPE(projection), INTENT(OUT) :: proj(2) ! Output the projection kernel
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmological model
      INTEGER :: nk, i

      IF (ix(1) == ix(2)) THEN
         nk = 1
      ELSE
         nk = 2
      END IF

      ! Loop over the two kernels
      DO i = 1, nk
         CALL fill_projection_kernel(ix(i), proj(i), cosm)
      END DO

      ! In case the autospectrum is being considered
      IF (nk == 1) THEN
         proj(2) = proj(1)
      END IF

   END SUBROUTINE fill_projection_kernels

   SUBROUTINE fill_projection_kernel(ix, proj, cosm)

      ! Fill look-up table for a projection kernel X_i
      INTEGER, INTENT(IN) :: ix
      TYPE(projection), INTENT(OUT) :: proj
      TYPE(cosmology), INTENT(INOUT) :: cosm
      !TYPE(lensing) :: lens

      IF (ix == tracer_Compton_y .OR. &
          ix == tracer_gravity_wave .OR. &
          ix == tracer_CIB_353 .OR. &
          ix == tracer_CIB_545 .OR. &
          ix == tracer_CIB_857 .OR. &
          ix == tracer_galaxies) THEN
         CALL fill_kernel(ix, proj, cosm)
      ELSE IF (ix == tracer_RCSLenS .OR. &
               ix == tracer_CFHTLenS_vanWaerbeke2013 .OR. &
               ix == tracer_CMB_lensing .OR. &
               ix == tracer_KiDS .OR. &
               ix == tracer_KiDS_bin1 .OR. &
               ix == tracer_KiDS_bin2 .OR. &
               ix == tracer_KiDS_bin3 .OR. &
               ix == tracer_KiDS_bin4 .OR. &
               ix == tracer_KiDS_450 .OR. &
               ix == tracer_KiDS_450_fat_bin1 .OR. &
               ix == tracer_KiDS_450_fat_bin2 .OR. &
               ix == tracer_KiDS_450_highz .OR. &
               ix == tracer_lensing_z1p00 .OR. &
               ix == tracer_lensing_z0p75 .OR. &
               ix == tracer_lensing_z0p50 .OR. &
               ix == tracer_lensing_z0p25 .OR. &
               ix == tracer_KiDS_450_bin1 .OR. &
               ix == tracer_KiDS_450_bin2 .OR. &
               ix == tracer_KiDS_450_bin3 .OR. &
               ix == tracer_KiDS_450_bin4 .OR. &
               ix == tracer_CFHTLenS_Kilbinger2013) THEN
         CALL fill_lensing_kernel(ix, proj, cosm)
      ELSE
         STOP 'FILL_PROJECTION_KERNEL: Error, tracer type specified incorrectly'
      END IF

   END SUBROUTINE fill_projection_kernel

   SUBROUTINE calculate_Cl(r1, r2, ell, Cell, nl, k, a, pow, nk, na, proj, cosm)

      ! Calculates C(l) using the Limber approximation
      ! Note that using Limber and flat-sky for sensible results limits lmin to ell~10
      INTEGER, INTENT(IN) :: nl, nk, na ! Number of ell values
      REAL, INTENT(IN) :: r1, r2 ! Maximum and minimum comoving distances for integration
      REAL, INTENT(IN) :: ell(nl) ! Input array of desired ell values for C(ell)
      REAL, INTENT(OUT) :: Cell(nl) ! Output array for C(ell) values  
      REAL, INTENT(IN) :: k(nk), a(na), pow(nk, na) ! Input k, a and P(k,a) arrays    
      TYPE(projection), INTENT(IN) :: proj(2) ! Projection kernels for the Limber integration
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      REAL :: logk(nk), loga(na), logpow(nk, na)
      INTEGER :: i, j

      ! Create log tables to speed up 2D find routine in find_pkz
      logk = log(k)
      loga = log(a)
      DO j = 1, na
         logpow(:, j) = log((2.*pi**2)*pow(:, j)/k**3)
      END DO

      ! Write some useful things to the screen
      IF (verbose_Limber) THEN
         WRITE (*, *) 'CALCULATE_CL: ell min:', real(ell(1))
         WRITE (*, *) 'CALCULATE_CL: ell max:', real(ell(nl))
         WRITE (*, *) 'CALCULATE_CL: number of ell:', nl
         WRITE (*, *) 'CALCULATE_CL: kmin [h/Mpc]:', k(1)
         WRITE (*, *) 'CALCULATE_CL: kmax [h/Mpc]:', k(nk)
         WRITE (*, *) 'CALCULATE_CL: number of k:', nk
         WRITE (*, *) 'CALCULATE_CL: amin [h/Mpc]:', a(1)
         WRITE (*, *) 'CALCULATE_CL: amax [h/Mpc]:', a(na)
         WRITE (*, *) 'CALCULATE_CL: number of a:', na
         WRITE (*, *) 'CALCULATE_CL: Minimum distance [Mpc/h]:', real(r1)
         WRITE (*, *) 'CALCULATE_CL: Maximum distance [Mpc/h]:', real(r2)
      END IF

      ! Finally do the integration
      IF (verbose_Limber) WRITE (*, *) 'CALCULATE_CL: Doing calculation'
      DO i = 1, nl
         Cell(i) = integrate_Limber(ell(i), r1, r2, logk, loga, logpow, nk, na, acc_Limber, 3, proj, cosm)
      END DO
      IF (verbose_Limber) THEN
         WRITE (*, *) 'CALCULATE_CL: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE calculate_Cl

   SUBROUTINE Cl_contribution_ell(r1, r2, k, a, pow, nk, na, proj, cosm, fbase, fext)

      USE string_operations

      ! Calculates the contribution to each ell of C(l) as a function of z, k, r
      ! Note that using Limber and flat-sky for sensible results limits lmin to ~10
      INTEGER, INTENT(IN) :: nk, na ! Number of entries in k and a
      REAL, INTENT(IN) :: r1, r2 ! Integration range
      REAL, INTENT(IN) :: k(nk), a(na), pow(nk, na) ! Input arrays of k, a and P(k,a) 
      TYPE(projection), INTENT(IN) :: proj(2) ! Projection kernels
      TYPE(cosmology), INTENT(INOUT) :: cosm ! Cosmology
      CHARACTER(len=*), INTENT(IN) :: fbase
      CHARACTER(len=*), INTENT(IN) :: fext
      REAL :: logk(nk), loga(na), logpow(nk, na)
      REAL, ALLOCATABLE :: ell(:)
      INTEGER :: i, j
      CHARACTER(len=256) :: outfile

      REAL, PARAMETER :: lmin = 1e1
      REAL, PARAMETER :: lmax = 1e4
      INTEGER, PARAMETER :: nl = 17 ! Number of ell values to take, from l=1 to l=2**(n-1); 2^15 ~ 32,000

      ! Create log tables to speed up 2D find routine in find_pka
      logk = log(k)
      loga = log(a)
      DO j = 1, na
         logpow(:, j) = log((2.*pi**2)*pow(:, j)/k**3)
      END DO

      CALL fill_array(log(lmin), log(lmax), ell, nl)
      ell = exp(ell)

      ! Now call the contribution subroutine
      !fbase='data/Cl_contribution_ell_'
      !fext='.dat'
      DO i = 1, nl
         outfile = number_file(fbase, i, fext)
         CALL Limber_contribution(ell(i), r1, r2, logk, loga, logpow, nk, na, proj, cosm, outfile)
      END DO

   END SUBROUTINE Cl_contribution_ell

   SUBROUTINE write_Cl(l, Cl, nl, outfile, verbose)

      ! Write C(l) to a file; writes l, C(l), l(l+1)C(l)/2pi
      INTEGER, INTENT(IN) :: nl
      REAL, INTENT(IN) :: l(nl)
      REAL, INTENT(IN) :: Cl(nl)
      CHARACTER(len=256) :: outfile
      LOGICAL, INTENT(IN) :: verbose
      INTEGER :: i

      IF (verbose) WRITE (*, *) 'WRITE_CL: Writing data to: ', TRIM(outfile)
      OPEN (7, file=outfile)
      DO i = 1, nl
         WRITE (7, *) l(i), Cl(i), l(i)*(l(i)+1.)*Cl(i)/twopi
      END DO
      CLOSE (7)
      IF (verbose) THEN
         WRITE (*, *) 'WRITE_CL: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE write_Cl

   SUBROUTINE calculate_angular_xi(ibessel, m, th_tab, xi_tab, n, l_tab, cl_tab, nl, lmax)

      ! Calcuate the two-dimensional correlation functions given a C(ell) array
      ! This uses the direct sum over discrete ell do the transformation
      ! TODO: Can I bunch together ell in the sum when ell is high?
      INTEGER, INTENT(IN) :: m, n, nl   ! Number of Bessel function transforms
      INTEGER, INTENT(IN) :: ibessel(m) ! Bessel functions to use   
      REAL, INTENT(IN) :: th_tab(n)     ! Angles for correlation function [degrees]
      REAL, INTENT(OUT) :: xi_tab(m, n) ! Correlation function    
      REAL, INTENT(IN) :: l_tab(nl)     ! Array of ell
      REAL, INTENT(IN) :: cl_tab(nl)    ! Array of C(ell)
      INTEGER, INTENT(IN) :: lmax       ! Maximum ell to go to in summation
      INTEGER :: i, j, l
      INTEGER :: ifind
      REAL :: logl(nl), logCl(nl), Cl(lmax)
      REAL :: theta, xi(m)
      INTEGER, PARAMETER :: iorder=3
      INTEGER, PARAMETER :: imeth=2

      ! Method
      ! 1 - Summation with accuracy parameter (converges very slowly, accuracy condition seems to not work)
      ! 2 - Standard integration (does not converge well, crap)
      ! 3 - Standard summation
      ! 4 - Summation assuming C(l) precomputed (fastest)
      INTEGER, PARAMETER :: method = method_angular_xi

      ! Speed up find routine by doing logarithms in advance
      logl = log(l_tab)
      logCl = log(cl_tab)

      IF (method == 4) THEN
         ! Need to precompute C(l) at each integer l for this method
         IF (regular_spacing(logl)) THEN
            ifind = 1 ! Array is linearly spaced
         ELSE
            ifind = 3 ! Mid-point method
         END IF
         DO l = 1, lmax
            Cl(l) = exp(find(log(real(l)), logl, logCl, nl, iorder, ifind, imeth))
         END DO
      END IF

      IF (verbose_xi) WRITE (*, *) 'CALCULATE_ANGULAR_XI: Computing correlation functions'

      DO i = 1, n

         ! Get theta value and convert from degrees to radians
         theta = th_tab(i)/rad2deg

         ! Calculate the correlation function
         DO j = 1, m

            IF (method == 1) THEN
               xi(j) = angular_xi_summation_adaptive(theta, ibessel(j), logl, logCl, nl, lmax, acc=acc_Limber)
            ELSE IF (method == 2) THEN
               xi(j) = integrate_angular_xi(theta, ibessel(j), logl, logCl, nl, acc_Limber)
            ELSE IF (method == 3) THEN
               IF (j == 2) EXIT
               xi = angular_xi_summation(theta, ibessel(j), m, logl, logCl, nl, lmax)
            ELSE IF (method == 4) THEN
               IF (j == 2) EXIT
               xi = angular_xi_summation_precompute(theta, ibessel, m, Cl, lmax)
            ELSE
               STOP 'CALCULATE_ANGULAR_XI: Error, method specified incorrectly'
            END IF

         END DO

         ! Populate tables
         xi_tab(:, i) = xi

      END DO

      IF (verbose_xi) THEN
         WRITE (*, *) 'CALCULATE_ANGULAR_XI: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE calculate_angular_xi

   REAL FUNCTION angular_xi_summation_adaptive(theta, nbessel, logl, logCl, nl, lmax, acc)

      USE special_functions

      ! This routine converges very slowly. The adaptive part seems to be terrible
      ! TODO: Probably remove this, although I am not sure why it is so bad
      INTEGER, INTENT(IN) :: nl       ! Number of entries in l-Cl arrays
      REAL, INTENT(IN) :: theta       ! Angle [rad]
      INTEGER, INTENT(IN) :: nbessel  ! Order form Bessel functions to use
      REAL, INTENT(IN) :: logl(nl)    ! Array of log l
      REAL, INTENT(IN) :: logCl(nl)   ! Array of log Cl  
      INTEGER, INTENT(IN) :: lmax     ! Maximum number of l to sum
      REAL, INTENT(IN) :: acc         ! Accuracy
      INTEGER :: l
      REAL :: Cl, xi
      REAL :: sum
      INTEGER, PARAMETER :: iorder = 3    ! Order for interpolation to find C(l) (3 - cubic)
      INTEGER, PARAMETER :: ifind = 3     ! Finding scheme for interpolation of C(l) (3 - mid-point)
      INTEGER, PARAMETER :: imeth = 2     ! Interpolation polynomial (2 - Lagrange polynomial)
      INTEGER, PARAMETER :: lmin = 100    ! Minimum number of l to sum

      ! Set values to zero before summing
      sum = 0.

      ! Do the conversion from Cl to xi as a summation over integer l
      DO l = 1, lmax

         ! Get the C(l) value associated with integer l
         Cl = exp(find(log(real(l)), logl, logCl, nl, iorder, ifind, imeth))

         ! Add to the running todays for the correlation functions
         xi = (2.*real(l)+1.)*Cl*Bessel(nbessel, real(l)*theta)
         sum = sum+xi

         ! Exit loop if sufficiently accurate
         IF (l > lmin .AND. ABS(xi/sum) < acc) THEN
            EXIT
         END IF

      END DO

      ! Divide by correct pre-factor
      angular_xi_summation_adaptive = real(sum)/(4.*pi)

   END FUNCTION angular_xi_summation_adaptive

   FUNCTION angular_xi_summation(theta, ibessel, n, logl, logCl, nl, lmax)

      USE special_functions

      INTEGER, INTENT(IN) :: n, nl          ! Number of Bessel functions to use
      REAL :: angular_xi_summation(n)
      REAL, INTENT(IN) :: theta         ! Angle [rad]
      INTEGER, INTENT(IN) :: ibessel(n) ! Order form Bessel functions to use    
      REAL, INTENT(IN) :: logl(nl)      ! Array of log l
      REAL, INTENT(IN) :: logCl(nl)     ! Array of log Cl  
      INTEGER, INTENT(IN) :: lmax       ! Maximum number of l to sum
      INTEGER :: l, i
      REAL :: Cl, xi(n)
      INTEGER :: ifind
      REAL :: sum(n)
      INTEGER, PARAMETER :: iorder = 3  ! Order for interpolation to find C(l) (3 - cubic)
      INTEGER, PARAMETER :: imeth = 2   ! Interpolation polynomial (2 - Lagrange polynomial)

      IF (regular_spacing(logl)) THEN
         ifind = 1 ! Array is linearly spaced
      ELSE
         ifind = 3 ! Mid-point method
      END IF

      ! Set values to zero before summing
      sum = 0.

      ! Do the conversion from Cl to xi as a summation over integer l
      DO l = 1, lmax

         ! Get the C(l) value associated with integer l
         Cl = exp(find(log(real(l)), logl, logCl, nl, iorder, ifind, imeth))

         ! Add to the running totals for the correlation functions
         DO i = 1, n
            xi(i) = (2.*real(l)+1.)*Cl*Bessel(ibessel(i), real(l)*theta)
            sum(i) = sum(i)+xi(i)
         END DO

      END DO

      ! Divide by correct pre-factor
      angular_xi_summation = real(sum)/(4.*pi)

   END FUNCTION angular_xi_summation

   FUNCTION angular_xi_summation_precompute(theta, ibessel, n, Cl, nl)

      USE special_functions

      INTEGER, INTENT(IN) :: n, nl          ! Number of Bessel functions to use
      REAL :: angular_xi_summation_precompute(n)
      REAL, INTENT(IN) :: theta         ! Angle [rad]
      INTEGER, INTENT(IN) :: ibessel(n) ! Order form Bessel functions to use 
      REAL, INTENT(IN) :: Cl(nl)     ! Array of log Cl 
      INTEGER :: l, i
      REAL :: xi(n)
      REAL :: sum(n)

      ! Set values to zero before summing
      sum = 0.

      ! Do the conversion from Cl to xi as a summation over integer l
      DO l = 1, nl

         ! Add to the running totals for the correlation functions
         DO i = 1, n
            xi(i) = (2.*real(l)+1.)*Cl(l)*Bessel(ibessel(i), real(l)*theta)
            sum(i) = sum(i)+xi(i)
         END DO

      END DO

      ! Divide by correct pre-factor
      angular_xi_summation_precompute = real(sum)/(4.*pi)

   END FUNCTION angular_xi_summation_precompute

   REAL FUNCTION integrate_angular_xi(theta, nbessel, logl, logCl, nl, acc)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      INTEGER, INTENT(IN) :: nl
      REAL, INTENT(IN) :: theta
      INTEGER, INTENT(IN) :: nbessel
      REAL, INTENT(IN) :: logl(nl)
      REAL, INTENT(IN) :: logCl(nl)
      REAL, INTENT(IN) :: acc
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass
      REAL, PARAMETER :: a = 0.
      REAL, PARAMETER :: b = 1.
      INTEGER, PARAMETER :: jmin = 5   ! Standard integration parameters
      INTEGER, PARAMETER :: jmax = 30  ! Standard integration parameters
      INTEGER, PARAMETER :: iorder = 3 ! Order for integration

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_angular_xi = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.
         sum_n = 0.
         sum_old = 0.
         sum_new = 0.

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = angular_xi_integrand(a, theta, nbessel, logl, logCl, nl)
               f2 = angular_xi_integrand(b, theta, nbessel, logl, logCl, nl)
               sum_2n = 0.5*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = progression(a, b, i, n)
                  fx = angular_xi_integrand(x, theta, nbessel, logl, logCl, nl)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_ANGULAR_XI: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (j == jmax) THEN
               STOP 'INTEGRATE_COSM_1: Integration timed out'
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO

         integrate_angular_xi = sum_new

      END IF

   END FUNCTION integrate_angular_xi

   REAL FUNCTION angular_xi_integrand(t, theta, nbessel, logl, logCl, nl)

      USE special_functions

      ! TODO: Write integration routine for this
      INTEGER, INTENT(IN) :: nl       ! Number of entries in l-Cl arrays
      REAL, INTENT(IN) :: t           ! Integration variable: 1+l = 1/t
      REAL, INTENT(IN) :: theta       ! Angle [rad]
      INTEGER, INTENT(IN) :: nbessel  ! Order of Bessel functions to use
      REAL, INTENT(IN) :: logl(nl)    ! Array of log l
      REAL, INTENT(IN) :: logCl(nl)   ! Array of log Cl  
      REAL :: l, Cl
      INTEGER, PARAMETER :: iorder = 3 ! Order for interpolation to find C(l) (3 - cubic)
      INTEGER, PARAMETER :: ifind = 3  ! Finding scheme for interpolation of C(l) (3 - mid-point)
      INTEGER, PARAMETER :: imeth = 2  ! Interpolation polynomial (2 - Lagrange polynomial)

      IF (t == 0.) THEN
         ! Corresponds to l=infinity, where Cl is zero (hopefully)
         angular_xi_integrand = 0.
      ELSE IF (t == 1.) THEN
         ! Corresponds to l=0, where Cl is zero (hopefully)
         angular_xi_integrand = 0.
      ELSE IF (t < 0. .OR. t > 1.) THEN
         STOP 'ANGULAR_XI_INTEGRAND: Error, t is out of bounds'
      ELSE
         l = -1.+1./t
         Cl = exp(find(log(l), logl, logCl, nl, iorder, ifind, imeth))
         angular_xi_integrand = l*Cl*Bessel(nbessel, l*theta)/t**2
      END IF

   END FUNCTION

   SUBROUTINE write_angular_xi(th_tab, xi_tab, nb, nth, outfile)

      ! Write angular correlation functions to disk
      ! TODO: Is this a bit absurd?
      INTEGER, INTENT(IN) :: nb, nth          ! Number of theta values
      REAL, INTENT(IN) :: th_tab(nth)     ! Array of theta values [probably degress]
      REAL, INTENT(IN) :: xi_tab(nb, nth) ! Correlation functions
      CHARACTER(len=256), INTENT(IN) :: outfile ! Output file
      INTEGER :: i, j

      OPEN (7, file=outfile)
      DO i = 1, nth
         WRITE (7, *) th_tab(i), (xi_tab(j, i), j=1, nb)
      END DO
      CLOSE (7)

   END SUBROUTINE write_angular_xi

   FUNCTION maxdist(proj)

      ! Calculates the maximum distance necessary for the lensing integration
      REAL :: maxdist
      TYPE(projection), INTENT(IN) :: proj(2)
      REAL :: rmax1, rmax2

      !Fix the maximum redshift and distance (which may fixed at source plane)
      rmax1 = maxval(proj(1)%r_x)
      rmax2 = maxval(proj(2)%r_x)

      !Subtract a small distance here because of rounding errors in recalculating zmax
      maxdist = min(rmax1, rmax2)-dr_max

   END FUNCTION maxdist

   SUBROUTINE write_projection_kernel(proj, cosm, outfile)

      TYPE(projection), INTENT(IN) :: proj
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=256), INTENT(IN) :: outfile
      INTEGER :: i
      REAL :: r, z

      ! Kernel
      IF (verbose_Limber) WRITE (*, *) 'WRITE_PROJECTION_KERNEL: Writing out kernel: ', trim(outfile)
      OPEN (7, file=outfile)
      DO i = 1, proj%nX
         r = proj%r_X(i)
         z = redshift_r(r, cosm)
         WRITE (7, *) r, z, proj%X(i)
      END DO
      CLOSE (7)
      IF (verbose_Limber) THEN
         WRITE (*, *) 'WRITE_PROJECTION_KERNEL: Writing done'
         WRITE (*, *)
      END IF

   END SUBROUTINE write_projection_kernel

   SUBROUTINE fill_lensing_kernel(ix, proj, cosm)

      ! Fill the lensing projection kernel
      INTEGER, INTENT(IN) :: ix
      !TYPE(lensing) :: lens
      TYPE(cosmology), INTENT(INOUT) :: cosm
      TYPE(projection), INTENT(OUT) :: proj
      REAL :: zmin, zmax, rmax, amax, r, dz
      INTEGER :: i, nX

      ! Choose either n(z) or fixed z_s
      IF (ix == tracer_CMB_lensing .OR. &
          ix == tracer_lensing_z1p00 .OR. &
          ix == tracer_lensing_z0p75 .OR. &
          ix == tracer_lensing_z0p50 .OR. &
          ix == tracer_lensing_z0p25) THEN
         zmin = 0.
         IF (ix == tracer_CMB_lensing) THEN
            zmax = cosm%z_cmb
         ELSE IF (ix == tracer_lensing_z1p00) THEN
            zmax = 1.00
         ELSE IF (ix == tracer_lensing_z0p75) THEN
            zmax = 0.75
         ELSE IF (ix == tracer_lensing_z0p50) THEN
            zmax = 0.50
         ELSE IF (ix == tracer_lensing_z0p25) THEN
            zmax = 0.25
         ELSE
            STOP 'FILL_LENSING_KERNEL: Error, tracer is specified incorrectly'
         END IF
         IF (verbose_Limber) WRITE (*, *) 'FILL_LENSING_KERNEL: Source plane redshift:', real(zmax)
      ELSE IF (ix == tracer_RCSLenS .OR. &
               ix == tracer_CFHTLenS_vanWaerbeke2013 .OR. &
               ix == tracer_KiDS .OR. &
               ix == tracer_KiDS_bin1 .OR. &
               ix == tracer_KiDS_bin2 .OR. &
               ix == tracer_KiDS_bin3 .OR. &
               ix == tracer_KiDS_bin4 .OR. &
               ix == tracer_KiDS_450 .OR. &
               ix == tracer_KiDS_450_fat_bin1 .OR. &
               ix == tracer_KiDS_450_fat_bin2 .OR. &
               ix == tracer_KiDS_450_highz .OR. &
               ix == tracer_KiDS_450_bin1 .OR. &
               ix == tracer_KiDS_450_bin2 .OR. &
               ix == tracer_KiDS_450_bin3 .OR. &
               ix == tracer_KiDS_450_bin4 .OR. &
               ix == tracer_CFHTLenS_Kilbinger2013) THEN
         CALL read_nz(ix, proj)
         IF (proj%order_nz == 0) THEN
            dz = proj%z_nz(2)-proj%z_nz(1)
            zmin = proj%z_nz(1)-dz/2.
            zmax = proj%z_nz(proj%nnz)+dz/2.
         ELSE
            zmin = proj%z_nz(1)
            zmax = proj%z_nz(proj%nnz)
         END IF
      ELSE
         STOP 'FILL_LENSING_KERNEL: Error, tracer is specified incorrectly'
      END IF

      ! Get the distance range for the lensing kernel
      amax = scale_factor_z(zmax)
      rmax = comoving_distance(amax, cosm)
      IF (verbose_Limber) THEN
         WRITE (*, *) 'FILL_LENSING_KERNEL: minimum r [Mpc/h]:', real(rmin_lensing)
         WRITE (*, *) 'FILL_LENSING_KERNEL: maximum r [Mpc/h]:', real(rmax)
         WRITE (*, *) 'FILL_LENSING_KERNEL: minimum z:', real(zmin)
         WRITE (*, *) 'FILL_LENSING_KERNEL: maximum z:', real(zmax)
         WRITE (*, *)
      END IF

      ! Fill the q(r) table
      CALL fill_efficiency(ix, rmin_lensing, rmax, zmax, proj, cosm)

      ! Assign arrays for the projection function
      nX = nX_lensing
      proj%nx = nx
      CALL fill_array(rmin_lensing, rmax, proj%r_x, nx)
      IF (verbose_Limber) WRITE (*, *) 'FILL_LENSING_KERNEL: number of points:', nx
      IF (ALLOCATED(proj%x)) DEALLOCATE (proj%x)
      ALLOCATE (proj%x(proj%nx))

      DO i = 1, nx
         ! Get r and fill X(r)
         r = proj%r_x(i)
         proj%x(i) = lensing_kernel(r, proj, cosm)
         ! Enforce that the kernel must not be negative (is this necessary?)
         IF (proj%x(i) < 0.) proj%x(i) = 0.
      END DO
      IF (verbose_Limber) THEN
         WRITE (*, *) 'FILL_LENSING_KERNEL: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE fill_lensing_kernel

   SUBROUTINE fill_efficiency(ix, rmin, rmax, zmax, proj, cosm)

      ! Fill the table for lensing efficiency
      INTEGER, INTENT(IN) :: ix
      REAL, INTENT(IN) :: rmin, rmax, zmax
      TYPE(projection), INTENT(INOUT) :: proj
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: r, z
      INTEGER :: i
      INTEGER, PARAMETER :: iorder=3

      ! Fill the r vs. q(r) tables
      proj%nq = nq_efficiency
      IF (verbose_Limber) WRITE (*, *) 'FILL_EFFICIENCY: number of points:', proj%nq
      CALL fill_array(rmin, rmax, proj%r_q, proj%nq)
      IF (ALLOCATED(proj%q)) DEALLOCATE (proj%q)
      ALLOCATE (proj%q(proj%nq))

      DO i = 1, proj%nq
         r = proj%r_q(i)
         z = redshift_r(r, cosm)
         !WRITE(*,*) i, z, zmax
         IF (r == 0.) THEN
            ! To avoid division by zero
            proj%q(i) = 1.
         ELSE IF (i == proj%nq) THEN
            proj%q(i) = 0.
         ELSE
            IF (ix == tracer_CMB_lensing .OR. &
                ix == tracer_lensing_z1p00 .OR. &
                ix == tracer_lensing_z0p75 .OR. &
                ix == tracer_lensing_z0p50 .OR. &
                ix == tracer_lensing_z0p25) THEN
               ! q(r) for a fixed source plane
               proj%q(i) = f_k(rmax-r, cosm)/f_k(rmax, cosm)
            ELSE IF (ix == tracer_RCSLenS .OR. &
                     ix == tracer_CFHTLenS_vanWaerbeke2013 .OR. &
                     ix == tracer_KiDS .OR. &
                     ix == tracer_KiDS_bin1 .OR. &
                     ix == tracer_KiDS_bin2 .OR. &
                     ix == tracer_KiDS_bin3 .OR. &
                     ix == tracer_KiDS_bin4 .OR. &
                     ix == tracer_KiDS_450 .OR. &
                     ix == tracer_KiDS_450_fat_bin1 .OR. &
                     ix == tracer_KiDS_450_fat_bin2 .OR. &
                     ix == tracer_KiDS_450_highz .OR. &
                     ix == tracer_KiDS_450_bin1 .OR. &
                     ix == tracer_KiDS_450_bin2 .OR. &
                     ix == tracer_KiDS_450_bin3 .OR. &
                     ix == tracer_KiDS_450_bin4 .OR. &
                     ix == tracer_CFHTLenS_Kilbinger2013) THEN
               ! q(r) for a n(z) distribution
               proj%q(i) = integrate_q(r, z, zmax, acc_Limber, iorder, proj, cosm)
            ELSE
               STOP 'FILL_EFFICIENCY: Error, tracer specified incorrectly'
            END IF
         END IF
      END DO

      IF (verbose_Limber) THEN
         WRITE (*, *) 'FILL_EFFICIENCY: Done writing'
         WRITE (*, *)
      END IF

   END SUBROUTINE fill_efficiency

   FUNCTION q_r(r, proj)

      ! Interpolation function for q(r)
      REAL :: q_r
      REAL, INTENT(IN) :: r
      TYPE(projection), INTENT(IN) :: proj

      q_r = find(r, proj%r_q, proj%q, proj%nq, 3, 3, 2)

   END FUNCTION q_r

   SUBROUTINE write_efficiency(proj, cosm, outfile)

      ! Write lensing efficiency q(r) function to a file
      TYPE(projection), INTENT(INOUT) :: proj
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=256), INTENT(IN) :: outfile
      REAL :: r, z, q
      INTEGER :: i

      WRITE (*, *) 'WRITE_EFFICIENCY: Writing q(r): ', trim(outfile)
      OPEN (7, file=outfile)
      DO i = 1, proj%nq
         r = proj%r_q(i)
         z = redshift_r(r, cosm)
         q = proj%q(i)
         WRITE (7, *) r, z, q
      END DO
      CLOSE (7)
      WRITE (*, *) 'WRITE_EFFICIENCY: Done'
      WRITE (*, *)

   END SUBROUTINE write_efficiency

   FUNCTION lensing_kernel(r, proj, cosm)

      ! The lensing projection kernel
      REAL :: lensing_kernel
      REAL, INTENT(IN) :: r
      TYPE(projection) :: proj
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: z, q

      ! Get z(r)
      z = redshift_r(r, cosm)

      ! Get the lensing efficiency
      q = q_r(r, proj)

      ! This is then the projection kernel (X_kappa)
      lensing_kernel = (1.+z)*f_k(r, cosm)*q
      lensing_kernel = lensing_kernel*1.5*cosm%om_m/(Hdist**2)

   END FUNCTION lensing_kernel

   SUBROUTINE fill_kernel(ix, proj, cosm)

      ! Fill a table of kernel values
      INTEGER, INTENT(IN) :: ix
      TYPE(projection), INTENT(OUT) :: proj
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i, nX
      REAL :: rmax, r

      nX = nx_kernel
      proj%nx = nx

      ! Get the distance range for the projection function
      ! Use the same as that for the distance calculation
      ! Assign arrays for the kernel function
      !rmax=comoving_distance(0.,cosm) ! Cheat to ensure that init_distance has been run
      !rmax=maxval(cosm%r)
      rmax = comoving_distance(scale_factor_z(cosm%z_CMB), cosm) ! TODO: Worry
      CALL fill_array(rmin_kernel, rmax, proj%r_X, proj%nX)
      IF (verbose_Limber) THEN
         WRITE (*, *) 'FILL_KERNEL: minimum r [Mpc/h]:', real(rmin_kernel)
         WRITE (*, *) 'FILL_KERNEL: maximum r [Mpc/h]:', real(rmax)
         WRITE (*, *) 'FILL_KERNEL: number of points:', nX
      END IF

      IF (ALLOCATED(proj%X)) DEALLOCATE (proj%X)
      ALLOCATE (proj%X(nX))

      ! Now fill the kernels
      DO i = 1, nX
         r = proj%r_x(i)
         IF (ix == tracer_Compton_y) THEN
            proj%x(i) = y_kernel(r, cosm)
         ELSE IF (ix == tracer_gravity_wave) THEN
            proj%x(i) = gwave_kernel(r, cosm)
         ELSE IF (ix == tracer_CIB_353 .OR. ix == tracer_CIB_545 .OR. ix == tracer_CIB_857) THEN
            proj%x(i) = CIB_kernel(r, cosm)
         ELSE IF (ix == tracer_galaxies) THEN
            proj%x(i) = galaxy_kernel(r, cosm)
         ELSE
            STOP 'FILL_KERNEL: Error, tracer specified incorrectly'
         END IF
      END DO

      IF (verbose_Limber) THEN
         WRITE (*, *) 'FILL_KERNEL: Done:'
         WRITE (*, *)
      END IF

   END SUBROUTINE fill_kernel

   FUNCTION y_kernel(r, cosm)

      ! The Compton-y projection kernel
      ! TODO: (re)check the power of 'a'
      REAL :: y_kernel
      REAL, INTENT(IN) :: r
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: z, a

      ! Get the scale factor
      z = redshift_r(r, cosm)
      a = scale_factor_z(z)

      ! Make the kernel and do some unit conversions
      y_kernel = yfac                     ! yfac = sigma_T / m_e c^2 [kg^-1 s^2]
      y_kernel = y_kernel*Mpc/cosm%h      ! Add Mpc/h units from the dr in the integral (h is new)
      y_kernel = y_kernel/a**2            ! These come from 'a^-3' for pressure multiplied by 'a' for comoving distance
      y_kernel = y_kernel*eV*(0.01)**(-3) ! Convert units of pressure spectrum from [eV/cm^3] to [J/m^3]

   END FUNCTION y_kernel

   REAL FUNCTION gwave_kernel(r, cosm)

      ! Projection kernel for gravitational waves (~ 1/r)
      REAL, INTENT(IN) :: r
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: crap

      ! To stop compile-time warnings
      crap = cosm%Om_m

      IF (r < rmin_gwave) THEN
         gwave_kernel = A_gwave/rmin_gwave
      ELSE
         gwave_kernel = A_gwave/r
      END IF

   END FUNCTION gwave_kernel

   REAL FUNCTION CIB_kernel(r, cosm)

      ! Projection kernel CIB emission
      REAL, INTENT(IN) :: r
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: z, a
      REAL :: rmin = 1. ! Minimum radius [Mpc/h]

      IF (r < rmin) THEN
         CIB_kernel = 0.
      ELSE
         z = redshift_r(r, cosm)
         a = scale_factor_z(z)
         CIB_kernel = 1./((1.+z)*luminosity_distance(a, cosm)**2)
      END IF

   END FUNCTION CIB_kernel

   REAL FUNCTION galaxy_kernel(r, cosm)

      USE special_functions

      REAL, INTENT(IN) :: r
      TYPE(cosmology), INTENT(INOUT) :: cosm
      REAL :: z, a
      REAL, PARAMETER :: sig = 0.5

      z = redshift_r(r, cosm)
      a = scale_factor_z(z)

      !nz=Rayleigh(z,sig) ! No, because this is included in the window function for galaxiex

      galaxy_kernel = sqrt(Hubble2(a, cosm))/Hdist

   END FUNCTION galaxy_kernel

   SUBROUTINE read_nz(ix, proj)
    
      USE calculus_table

      ! The the n(z) function for lensing
      INTEGER, INTENT(IN) :: ix
      TYPE(projection), INTENT(INOUT) :: proj
      REAL :: norm
      LOGICAL, PARAMETER :: normalise = norm_nz

      IF (ix == tracer_RCSLenS .OR. ix == tracer_CFHTLenS_vanWaerbeke2013) THEN
         CALL fill_analytic_nz_table(ix, proj)
      ELSE IF (ix == tracer_KiDS .OR. &
               ix == tracer_KiDS_bin1 .OR. &
               ix == tracer_KiDS_bin2 .OR. &
               ix == tracer_KiDS_bin3 .OR. &
               ix == tracer_KiDS_bin4 .OR. &
               ix == tracer_KiDS_450 .OR. &
               ix == tracer_KiDS_450_fat_bin1 .OR. &
               ix == tracer_KiDS_450_fat_bin2 .OR. &
               ix == tracer_KiDS_450_highz .OR. &
               ix == tracer_KiDS_450_bin1 .OR. &
               ix == tracer_KiDS_450_bin2 .OR. &
               ix == tracer_KiDS_450_bin3 .OR. &
               ix == tracer_KiDS_450_bin4 .OR. &
               ix == tracer_CFHTLenS_Kilbinger2013) THEN
         CALL fill_nz_table(ix, proj)
      ELSE
         STOP 'READ_NZ: Error, tracer specified incorrectly'
      END IF

      IF (ix == tracer_RCSLenS .OR. ix == tracer_CFHTLenS_vanWaerbeke2013) THEN
         proj%order_nz = 3 ! Continuous functions - treat as cubics
      ELSE IF (ix == tracer_KiDS .OR. &
               ix == tracer_KiDS_bin1 .OR. &
               ix == tracer_KiDS_bin2 .OR. &
               ix == tracer_KiDS_bin3 .OR. &
               ix == tracer_KiDS_bin4) THEN
         proj%order_nz = 3 ! These seem to integrate more accurately to unity this way
      ELSE IF (ix == tracer_KiDS_450 .OR. &
               ix == tracer_KiDS_450_fat_bin1 .OR. &
               ix == tracer_KiDS_450_fat_bin2 .OR. &
               ix == tracer_KiDS_450_highz .OR. &
               ix == tracer_KiDS_450_bin1 .OR. &
               ix == tracer_KiDS_450_bin2 .OR. &
               ix == tracer_KiDS_450_bin3 .OR. &
               ix == tracer_KiDS_450_bin4 .OR. &
               ix == tracer_CFHTLenS_Kilbinger2013) THEN
         proj%order_nz = 0 ! These are supposed to be histograms
      ELSE
         STOP 'READ_NZ: Error, order for integration not specified correctly'
      END IF

      IF (normalise) CALL normalise_nz(proj)

      IF (verbose_Limber) THEN
         norm = integrate_nz(proj)
         WRITE (*, *) 'READ_NZ: Normalisation (should be unity):', norm
         WRITE (*, *) 'READ_NZ: Order for integration:', proj%order_nz
         WRITE (*, *) 'READ_NZ: zmin:', proj%z_nz(1)
         WRITE (*, *) 'READ_NZ: zmax:', proj%z_nz(proj%nnz)
         WRITE (*, *) 'READ_NZ: nz:', proj%nnz
         WRITE (*, *)
      END IF

   END SUBROUTINE read_nz

   SUBROUTINE normalise_nz(proj)

      ! Explicitly normalise the input n(z)
      TYPE(projection), INTENT(INOUT) :: proj
      REAL :: norm

      ! Explicitly normalise
      norm = integrate_nz(proj)
      proj%nz = proj%nz/norm

   END SUBROUTINE normalise_nz

   REAL FUNCTION integrate_nz(proj)

      USE calculus_table

      ! Integrate the n(z) function over z (which should give unity)
      TYPE(projection), INTENT(INOUT) :: proj

      integrate_nz = integrate_table(proj%z_nz, proj%nz, 1, proj%nnz, proj%order_nz)

   END FUNCTION integrate_nz

   SUBROUTINE fill_analytic_nz_table(ix, proj)

      INTEGER, INTENT(IN) :: ix
      TYPE(projection), INTENT(INOUT) :: proj
      INTEGER :: i

      ! From analytical function
      proj%nnz = n_nz
      IF (ALLOCATED(proj%z_nz)) DEALLOCATE (proj%z_nz)
      IF (ALLOCATED(proj%nz)) DEALLOCATE (proj%nz)
      ALLOCATE (proj%z_nz(proj%nnz), proj%nz(proj%nnz))

      ! Fill the look-up tables
      CALL fill_array(zmin_nz, zmax_nz, proj%z_nz, proj%nnz)
      DO i = 1, n_nz
         proj%nz(i) = nz_distribution(proj%z_nz(i), ix)
      END DO

   END SUBROUTINE fill_analytic_nz_table

   SUBROUTINE fill_nz_table(ix, proj)

      USE file_info

      INTEGER, INTENT(IN) :: ix
      TYPE(projection), INTENT(INOUT) :: proj
      INTEGER :: i
      REAL :: spam
      CHARACTER(len=256) :: input

      ! Get file name
      IF (ix == tracer_KiDS) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/KiDS_z0.1-0.9.txt'
      ELSE IF (ix == tracer_KiDS_bin1) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/KiDS_z0.1-0.3.txt'
      ELSE IF (ix == tracer_KiDS_bin2) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/KiDS_z0.3-0.5.txt'
      ELSE IF (ix == tracer_KiDS_bin3) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/KiDS_z0.5-0.7.txt'
      ELSE IF (ix == tracer_KiDS_bin4) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/KiDS_z0.7-0.9.txt'
      ELSE IF (ix == tracer_KiDS_450 .OR. &
               ix == tracer_KiDS_450_fat_bin1 .OR. &
               ix == tracer_KiDS_450_fat_bin2 .OR. &
               ix == tracer_KiDS_450_highz) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/KiDS-450_fat_bin_nofz.txt'
      ELSE IF (ix == tracer_KiDS_450_bin1) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/Nz_DIR_z0.1t0.3.asc'
      ELSE IF (ix == tracer_KiDS_450_bin2) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/Nz_DIR_z0.3t0.5.asc'
      ELSE IF (ix == tracer_KiDS_450_bin3) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/Nz_DIR_z0.5t0.7.asc'
      ELSE IF (ix == tracer_KiDS_450_bin4) THEN
         input = '/Users/Mead/Physics/data/KiDS/nz/Nz_DIR_z0.7t0.9.asc'
      ELSE IF (ix == tracer_CFHTLenS_Kilbinger2013) THEN
         input = '/Users/Mead/Physics/data/CFHTLenS/nz.txt'
      ELSE
         STOP 'FILL_NZ_TABLE: tracer not specified correctly'
      END IF
      IF (verbose_Limber) WRITE (*, *) 'FILL_NZ_TABLE: Input file: ', trim(input)

      ! Allocate arrays
      proj%nnz = count_number_of_lines(input)
      IF (ALLOCATED(proj%z_nz)) DEALLOCATE (proj%z_nz)
      IF (ALLOCATED(proj%nz))   DEALLOCATE (proj%nz)
      ALLOCATE (proj%z_nz(proj%nnz), proj%nz(proj%nnz))

      ! Read in n(z) table
      OPEN (7, file=input)
      DO i = 1, proj%nnz
         IF (ix == tracer_KiDS .OR. &
             ix == tracer_KiDS_bin1 .OR. &
             ix == tracer_KiDS_bin2 .OR. &
             ix == tracer_KiDS_bin3 .OR. &
             ix == tracer_KiDS_bin4 .OR. &
             ix == tracer_KiDS_bin4 .OR. &
             ix == tracer_KiDS_450_bin1 .OR. &
             ix == tracer_KiDS_450_bin2 .OR. &
             ix == tracer_KiDS_450_bin3 .OR. &
             ix == tracer_KiDS_450_bin4 .OR. &
             ix == tracer_CFHTLenS_Kilbinger2013) THEN
            READ (7, *) proj%z_nz(i), proj%nz(i) ! Second column
         ELSE IF (ix == tracer_KiDS_450) THEN
            READ (7, *) proj%z_nz(i), proj%nz(i) ! Second column (z = 0.1 -> 0.9)
         ELSE IF (ix == tracer_KiDS_450_fat_bin1) THEN
            READ (7, *) proj%z_nz(i), spam, proj%nz(i) ! Third column (z = 0.1 -> 0.5)
         ELSE IF (ix == tracer_KiDS_450_fat_bin2) THEN
            READ (7, *) proj%z_nz(i), spam, spam, proj%nz(i) ! Fourth column (z = 0.5 -> 0.9)
         ELSE IF (ix == tracer_KiDS_450_highz) THEN
            READ (7, *) proj%z_nz(i), spam, spam, spam, proj%nz(i) ! Fifth column (z = 0.9 -> 3.5)
         ELSE
            STOP 'FILL_NZ_TABLE: tracer not specified correctly'
         END IF
      END DO
      CLOSE (7)

      ! Do this because the KiDS-450 files contain the lower left edge of histograms
      ! The bin sizes are 0.05 in z, so need to add 0.05/2 = 0.025
      ! This is also true for the CFHTLenS Kilbinger (2013) data
      IF (ix == tracer_KiDS_450 .OR. &
          ix == tracer_KiDS_450_fat_bin1 .OR. &
          ix == tracer_KiDS_450_fat_bin2 .OR. &
          ix == tracer_KiDS_450_highz .OR. &
          ix == tracer_KiDS_450_bin1 .OR. &
          ix == tracer_KiDS_450_bin2 .OR. &
          ix == tracer_KiDS_450_bin3 .OR. &
          ix == tracer_KiDS_450_bin4 .OR. &
          ix == tracer_CFHTLenS_Kilbinger2013) THEN
         proj%z_nz = proj%z_nz+0.025
      END IF

      IF (verbose_Limber) THEN
         WRITE (*, *) 'FILL_NZ_TABLE: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE fill_nz_table

   REAL FUNCTION nz_distribution(z, ix)

      ! Analytical n(z) for different surveys
      REAL, INTENT(IN) :: z
      INTEGER, INTENT(IN) :: ix
      REAL :: a, b, c, d, e, f, g, h, i
      REAL :: norm
      REAL :: n1, n2, n3
      REAL :: z1, z2

      IF (ix == tracer_RCSLenS) THEN
         ! RCSLenS
         a = 2.94
         b = -0.44
         c = 1.03
         d = 1.58
         e = 0.40
         f = 0.25
         g = 0.38
         h = 0.81
         i = 0.12
         n1 = a*z*exp(-(z-b)**2/c**2)
         n2 = d*z*exp(-(z-e)**2/f**2)
         n3 = g*z*exp(-(z-h)**2/i**2)
         nz_distribution = n1+n2+n3
      ELSE IF (ix == tracer_CFHTLenS_vanWaerbeke2013) THEN
         ! CFHTLenS
         z1 = 0.7 ! Not a free parameter in Van Waerbeke et al. fit (2013)
         z2 = 1.2 ! Not a free parameter in Van Waerbeke et al. fit (2013)
         a = 1.50
         b = 0.32
         c = 0.20
         d = 0.46
         norm = 1.0129840620118542 ! This is to ensure normalisation; without this integrates to ~1.013 according to python
         nz_distribution = (a*exp(-((z-z1)/b)**2)+c*exp(-((z-z2)/d)**2))/norm
      ELSE
         STOP 'NZ_DISTRIBUTION: Error, tracer specified incorrectly'
      END IF

   END FUNCTION nz_distribution

   REAL FUNCTION integrate_q(r, a, b, acc, iorder, proj, cosm)

      ! Integrates between a and b until desired accuracy is reached
      ! Stores information to reduce function calls
      ! TODO: This could be changed to a general integrate_proj_cosm routine
      REAL, INTENT(IN) :: r
      REAL, INTENT(IN) :: a
      REAL, INTENT(IN) :: b
      REAL, INTENT(IN) :: acc
      INTEGER, INTENT(IN) :: iorder
      TYPE(cosmology), INTENT(INOUT) :: cosm
      TYPE(projection), INTENT(IN) :: proj
      INTEGER :: i, j, n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = 5  ! Standard integration parameters
      INTEGER, PARAMETER :: jmax = 30 ! Standard integration parameters

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_q = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.
         sum_n = 0.
         sum_old = 0.
         sum_new = 0.

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = q_integrand(a, r, proj, cosm)
               f2 = q_integrand(b, r, proj, cosm)
               sum_2n = 0.5*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = progression(a, b, i, n)
                  fx = q_integrand(x, r, proj, cosm)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_Q: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (j == jmax) THEN
               STOP 'INTEGRATE_COSM_1: Integration timed out'
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO

         integrate_q = sum_new

      END IF

   END FUNCTION integrate_q

   FUNCTION q_integrand(z, r, proj, cosm)

      ! The lensing efficiency integrand, which is a function of z
      ! z is integrated over while r is just a parameter
      ! This is only called for n(z)
      REAL :: q_integrand
      REAL, INTENT(IN) :: r, z
      TYPE(cosmology), INTENT(INOUT) :: cosm
      TYPE(projection), INTENT(IN) :: proj
      REAL :: rdash, nzz, a
      INTEGER, PARAMETER :: ifind = 3
      INTEGER, PARAMETER :: imeth = 2

      a = scale_factor_z(z)

      IF (z == 0.) THEN
         q_integrand = 0.
      ELSE
         ! Find the r'(z) variable that is integrated over
         rdash = comoving_distance(a, cosm)
         ! Find the n(z)
         nzz = find(z, proj%z_nz, proj%nz, proj%nnz, proj%order_nz, ifind, imeth)
         ! This is then the integrand
         q_integrand = nzz*f_k(rdash-r, cosm)/f_k(rdash, cosm)
      END IF

   END FUNCTION q_integrand

   REAL FUNCTION integrate_Limber(l, a, b, logktab, logatab, logptab, nk, na, acc, iorder, proj, cosm)

      ! Integrates between a and b until desired accuracy is reached
      INTEGER, INTENT(IN) :: nk, na
      REAL, INTENT(IN) :: a, b, acc
      INTEGER, INTENT(IN) :: iorder
      REAL, INTENT(IN) :: l
      REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk, na)  
      TYPE(projection), INTENT(IN) :: proj(2)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER :: i, j
      INTEGER :: n
      REAL :: x, dx
      REAL :: f1, f2, fx
      REAL :: sum_n, sum_2n, sum_new, sum_old
      LOGICAL :: pass

      INTEGER, PARAMETER :: jmin = 5  ! Standard integration parameters
      INTEGER, PARAMETER :: jmax = 25 ! Standard integration parameters

      IF (a == b) THEN

         ! Fix the answer to zero if the integration limits are identical
         integrate_Limber = 0.

      ELSE

         ! Set the sum variable for the integration
         sum_2n = 0.
         sum_n = 0.
         sum_old = 0.
         sum_new = 0.

         DO j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n = 1+2**(j-1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b-a)/real(n-1)

            IF (j == 1) THEN

               ! The first go is just the trapezium of the end points
               f1 = Limber_integrand(a, l, logktab, logatab, logptab, nk, na, proj, cosm)
               f2 = Limber_integrand(b, l, logktab, logatab, logptab, nk, na, proj, cosm)
               sum_2n = 0.5*(f1+f2)*dx
               sum_new = sum_2n

            ELSE

               ! Loop over only new even points to add these to the integral
               DO i = 2, n, 2
                  x = progression(a, b, i, n)
                  fx = Limber_integrand(x, l, logktab, logatab, logptab, nk, na, proj, cosm)
                  sum_2n = sum_2n+fx
               END DO

               ! Now create the total using the old and new parts
               sum_2n = sum_n/2.+sum_2n*dx

               ! Now calculate the new sum depending on the integration order
               IF (iorder == 1) THEN
                  sum_new = sum_2n
               ELSE IF (iorder == 3) THEN
                  sum_new = (4.*sum_2n-sum_n)/3. ! This is Simpson's rule and cancels error
               ELSE
                  STOP 'INTEGRATE_LIMBER: Error, iorder specified incorrectly'
               END IF

            END IF

            IF (sum_old == 0.d0 .OR. j<jmin) THEN
               pass = .FALSE.
            ELSE IF (j == jmax) THEN
               STOP 'INTEGRATE_COSM_1: Integration timed out'
            ELSE IF (abs(-1.d0+sum_new/sum_old) < acc) THEN
               pass = .TRUE.
            ELSE
               pass = .FALSE.
            END IF

            IF (pass) THEN
               EXIT
            ELSE
               ! Integral has not converged so store old sums and reset sum variables
               sum_old = sum_new
               sum_n = sum_2n
               sum_2n = 0.d0
            END IF

         END DO

         integrate_Limber = sum_new

      END IF

   END FUNCTION integrate_Limber

   REAL FUNCTION Limber_integrand(r, l, logktab, logatab, logptab, nk, na, proj, cosm)

      ! The integrand for the Limber integral
      INTEGER, INTENT(IN) :: nk, na
      REAL, INTENT(IN) :: r, l
      REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk, na)
      
      TYPE(cosmology), INTENT(INOUT) :: cosm
      TYPE(projection), INTENT(IN) :: proj(2)
      REAL :: z, a, k, X(2)
      INTEGER :: i

      IF (r == 0.) THEN

         ! Must be set to zero here to prevent division by zero
         Limber_integrand = 0.

      ELSE

         ! Get the two kernels
         DO i = 1, 2
            X(i) = find(r, proj(i)%r_X, proj(i)%X, proj(i)%nX, 3, 3, 2)
            !X(i)=exp(find(log(r),log(proj(i)%r_X),log(proj(i)%X),proj(i)%nX,3,3,2)) ! Barfed with this
         END DO

         ! Get variables r, z(r) and k(r) for P(k,z)
         z = redshift_r(r, cosm)
         a = scale_factor_z(z)
         k = (l+lcorr)/f_k(r, cosm) ! LoVerde et al. (2008) Limber correction
         !k=k_ell(l,a,cosm) ! Does not make sense since this does an internal a->r conversion

         ! Construct the integrand
         Limber_integrand = X(1)*X(2)*find_pka(k, a, logktab, logatab, logptab, nk, na)/f_k(r, cosm)**2

      END IF

   END FUNCTION Limber_integrand

   SUBROUTINE Limber_contribution(l, r1, r2, logktab, logatab, logptab, nk, na, proj, cosm, outfile)

      INTEGER, INTENT(IN) :: nk, na
      REAL, INTENT(IN) :: l, r1, r2
      REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk, na)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      TYPE(projection), INTENT(IN) :: proj(2)
      CHARACTER(len=256), INTENT(IN) :: outfile
      REAL :: k, r, z, total, int, a
      INTEGER :: i

      ! Calculate the integral for this value of ell
      total = integrate_Limber(l, r1, r2, logktab, logatab, logptab, nk, na, acc_Limber, 3, proj, cosm)

      ! Now split up the contributions
      ! You need the Jacobian and to remember that the contribution is split in ln(k), ln(z) and ln(R)
      ! This means factors of:
      ! k - f_k(r)/f'_k(r) (=r if flat)
      ! z - z/H(z)
      ! r - r
      OPEN (7, file=trim(outfile))
      DO i = 1, n_cont
         r = progression(r1, r2, i, n_cont)
         IF (r == 0.) THEN
            ! Avoid trouble for r = 0 exactly
            ! Maybe should do some sort of awful Taylor expansion here
            CYCLE
         ELSE
            k = (l+lcorr)/f_k(r, cosm)
            !k=k_ell(l,a,cosm) ! Does not make sense since this does an internal a->r conversion
            z = redshift_r(r, cosm)
            a = scale_factor_z(z)
            int = Limber_integrand(r, l, logktab, logatab, logptab, nk, na, proj, cosm)
            WRITE (7, *) k, int*f_k(r, cosm)/(fdash_k(r, cosm)*total), z, int*z/(sqrt(Hubble2(a, cosm))*total), r, int*r/total
         END IF
      END DO
      CLOSE (7)

   END SUBROUTINE Limber_contribution

   REAL FUNCTION find_pka(k, a, logktab, logatab, logptab, nk, na)

      ! Looks up the power as a 2D function of k and a
      ! Note that it cuts P(k,a) off above and below certain wavenumbers defined in the header (kmin_pka, kmax_pka)
      ! It will interpolate in log(a) outside range
      ! It will interpolate in log(k) outside range of ktab until kmin_pka/kmax_pka
      REAL, INTENT(IN) :: k, a ! Input desired values of k and a
      INTEGER, INTENT(IN) :: nk, na ! Number of entries of k and a in arrays
      REAL, INTENT(IN) :: logktab(nk), logatab(na), logptab(nk, na) ! Arrays of log(k), log(a) and log(P(k,a))
      REAL :: logk, loga
      INTEGER, PARAMETER :: iorder = 3 ! 3 - Cubic interpolation
      INTEGER, PARAMETER :: ifind = 3  ! 3 - Midpoint finding scheme
      INTEGER, PARAMETER :: imeth = 1  ! 1 - Polynomial method

      logk = log(k)
      loga = log(a)

      ! Get the power
      IF ((logk < logktab(1) .OR. logk > logktab(nk)) .AND. (loga < logatab(1) .OR. loga > logatab(na))) THEN
         find_pka = 0.
      ELSE IF (k < kmin_pka .OR. k > kmax_pka) THEN
         find_pka = 0.
      ELSE IF (a < amin_pka .OR. a > amax_pka) THEN
         find_pka = 0.
      ELSE
         find_pka = exp(find(logk, logktab, loga, logatab, logptab, nk, na, iorder, ifind, ifind, imeth))
      END IF

   END FUNCTION find_pka

   SUBROUTINE xpow_pka(ix, ell, Cl, nl, k, a, pow, nk, na, cosm)

      ! Calculates the C(l) for the cross correlation of fields ix(1) and ix(2) given P(k,a)
      ! TODO: Change to take in ix(n)?
      INTEGER, INTENT(IN) :: nk, na, nl
      INTEGER, INTENT(INOUT) :: ix(2)
      REAL, INTENT(IN) :: ell(nl)
      REAL, INTENT(OUT) :: Cl(nl)    
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(IN) :: a(na)
      REAL, INTENT(IN) :: pow(nk, na)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      TYPE(projection) :: proj(2)
      REAL :: r1, r2

      ! Fill out the projection kernel
      CALL fill_projection_kernels(ix, proj, cosm)

      ! Set the range in comoving distance for the Limber integral [Mpc]
      r1 = 0.
      r2 = maxdist(proj)

      ! Actually calculate the C(ell), but only for the full halo model part
      ! TODO: Array temporary
      CALL calculate_Cl(r1, r2, ell, Cl, nl, k, a, pow, nk, na, proj, cosm)

   END SUBROUTINE xpow_pka

   SUBROUTINE Cl_contribution(ix, k, a, pow, nk, na, cosm, fbase, fext)

      ! Calculates the C(l) for the cross correlation of fields ix(1) and ix(2) given P(k,a)
      ! TODO: Change to take in ix(n)?
      INTEGER, INTENT(IN) :: nk, na
      INTEGER, INTENT(INOUT) :: ix(2)
      REAL, INTENT(IN) :: k(nk)
      REAL, INTENT(IN) :: a(na)
      REAL, INTENT(IN) :: pow(nk, na)
      TYPE(cosmology), INTENT(INOUT) :: cosm
      CHARACTER(len=*), INTENT(IN) :: fbase
      CHARACTER(len=*), INTENT(IN) :: fext
      TYPE(projection) :: proj(2)
      REAL :: r1, r2

      ! Fill out the projection kernel
      CALL fill_projection_kernels(ix, proj, cosm)

      ! Set the range in comoving distance for the Limber integral [Mpc]
      r1 = 0.
      r2 = maxdist(proj)

      CALL Cl_contribution_ell(r1, r2, k, a, pow, nk, na, proj, cosm, fbase, fext)

   END SUBROUTINE Cl_contribution

   SUBROUTINE xpow_halomod(ix, nx, ell, Cl, nl, cosm, ihm, verbose)

      USE HMx

      ! Calculates the C(l) for the cross correlation of fields ix(1) and ix(2)
      ! TODO: Speed up if there are repeated fields in ix(n) (should this ever happen?)
      ! TODO: Add bin theory option
      ! TODO: Move this because it requires HMx?
      INTEGER, INTENT(IN) :: nx, nl
      INTEGER, INTENT(INOUT) :: ix(nx)
      REAL, INTENT(IN) :: ell(nl)
      REAL, INTENT(OUT) :: Cl(nl, nx, nx)
      !TYPE(halomod), INTENT(INOUT) :: hmod
      TYPE(cosmology), INTENT(INOUT) :: cosm
      INTEGER, INTENT(INOUT) :: ihm
      LOGICAL, INTENT(IN) :: verbose
      REAL, ALLOCATABLE :: a(:), k(:), pow_li(:, :), pow_2h(:, :, :, :), pow_1h(:, :, :, :), pow_hm(:, :, :, :)
      REAL :: lmin, lmax
      INTEGER :: ixx(2), ip(nx), nk, na, i, j, ii, jj, nnx, match(nx)
      INTEGER, ALLOCATABLE :: iix(:)
      REAL, ALLOCATABLE :: uCl(:, :, :)

      REAL, PARAMETER :: kmin_xpow = 1e-3 ! Minimum k to calculate P(k); it will be extrapolated below this by Limber
      REAL, PARAMETER :: kmax_xpow = 1e2  ! Maximum k to calculate P(k); it will be extrapolated above this by Limber
      INTEGER, PARAMETER :: nk_xpow = 257 ! Number of log-spaced k values (used to be 32)
      REAL, PARAMETER :: amin_xpow = 0.1  ! Minimum scale factor (problems with one-halo term if amin is less than 0.1)
      REAL, PARAMETER :: amax_xpow = 1.0  ! Maximum scale factor
      INTEGER, PARAMETER :: na_xpow = 17  ! Number of linearly-spaced scale factores

      !IF(repeated_entries(ix,n)) STOP 'XPOW: Error, repeated tracers'
      CALL unique_index(ix, nx, iix, nnx, match)
      ALLOCATE (uCl(nl, nnx, nnx))

      ! Set the k range
      nk = nk_xpow
      CALL fill_array(log(kmin_xpow), log(kmax_xpow), k, nk)
      k = exp(k)

      ! Set the a range
      na = na_xpow
      CALL fill_array(amin_xpow, amax_xpow, a, na)

      ! Set the ell range
      lmin = ell(1)
      lmax = ell(nl)

      ! Write to screen
      IF (verbose) THEN
         WRITE (*, *) 'XPOW: Cross-correlation information'
         WRITE (*, *) 'XPOW: Number of tracers:', nx
         WRITE (*, *) 'XPOW: P(k) minimum k [h/Mpc]:', REAL(kmin_xpow)
         WRITE (*, *) 'XPOW: P(k) maximum k [h/Mpc]:', REAL(kmax_xpow)
         WRITE (*, *) 'XPOW: Number of k:', nk
         WRITE (*, *) 'XPOW: minimum a:', REAL(amin_xpow)
         WRITE (*, *) 'XPOW: maximum a:', REAL(amax_xpow)
         WRITE (*, *) 'XPOW: number of a:', na
         WRITE (*, *) 'XPOW: minimum ell:', REAL(lmin)
         WRITE (*, *) 'XPOW: maximum ell:', REAL(lmax)
         WRITE (*, *) 'XPOW: number of ell:', nl
         WRITE (*, *)
      END IF

      ! Use the xpowlation type to set the necessary halo profiles
      DO i = 1, nnx
         CALL set_field_for_xpow(ix(i), ip(i))
      END DO

      ! Do the halo model power spectrum calculation
      !CALL calculate_HMx(ip, nnx, k, nk, a, na, pow_li, pow_2h, pow_1h, pow_hm, hmod, cosm, verbose)
      CALL calculate_HMx_full(ip, k, a, pow_li, pow_2h, pow_1h, pow_hm, nnx, nk, na, cosm, ihm)

      ! Set the Cl to zero initially
      uCl = 0.

      ! Loop over triangle combinations
      DO i = 1, nnx
         ixx(1) = ix(i)
         DO j = i, nnx
            ixx(2) = ix(j)
            CALL xpow_pka(ixx, ell, uCl(:, i, j), nl, k, a, pow_hm(i, j, :, :), nk, na, cosm)
         END DO
      END DO

      ! Fill the symmetric cross terms
      DO i = 1, nnx
         DO j = i+1, nnx
            uCl(:, j, i) = uCl(:, i, j)
         END DO
      END DO

      ! Now fill the full arrays from the unique arrays
      DO i = 1, nx
         DO j = 1, nx
            ii = match(i)
            jj = match(j)
            Cl(:, i, j) = uCl(:, ii, jj)
         END DO
      END DO

      ! Write to screen
      IF (verbose) THEN
         WRITE (*, *) 'XPOW: Done'
         WRITE (*, *)
      END IF

   END SUBROUTINE xpow_halomod

   SUBROUTINE set_field_for_xpow(ix, ip)

      USE HMx

      ! Set the cross-correlation type
      ! TODO: Move this because it requires HMx?
      INTEGER, INTENT(INOUT) :: ix ! Tracer for cross correlation
      INTEGER, INTENT(OUT) :: ip   ! Corresponding field for power spectrum
      INTEGER :: j

      IF (ix == -1) THEN
         WRITE (*, *) 'SET_FIELDS_FOR_XPOW: Choose field'
         WRITE (*, *) '================================='
         DO j = 1, n_tracers
            WRITE (*, fmt='(I3,A3,A30)') j, '- ', TRIM(xcorr_type(j))
         END DO
         READ (*, *) ix
         WRITE (*, *) '==================================='
         WRITE (*, *)
      END IF

      IF (ix == tracer_Compton_y) THEN
         ! Compton y
         ip = field_electron_pressure
      ELSE IF (ix == tracer_gravity_wave) THEN
         ! Gravitational waves
         ip = field_dmonly
      ELSE IF (ix == tracer_RCSLenS .OR. &
               ix == tracer_CFHTLenS_vanWaerbeke2013 .OR. &
               ix == tracer_CMB_lensing .OR. &
               ix == tracer_KiDS .OR. &
               ix == tracer_KiDS_bin1 .OR. &
               ix == tracer_KiDS_bin2 .OR. &
               ix == tracer_KiDS_bin3 .OR. &
               ix == tracer_KiDS_bin4 .OR. &
               ix == tracer_KiDS_450 .OR. &
               ix == tracer_KiDS_450_fat_bin1 .OR. &
               ix == tracer_KiDS_450_fat_bin2 .OR. &
               ix == tracer_KiDS_450_highz .OR. &
               ix == tracer_lensing_z1p00 .OR. &
               ix == tracer_lensing_z0p75 .OR. &
               ix == tracer_lensing_z0p50 .OR. &
               ix == tracer_lensing_z0p25 .OR. &
               ix == tracer_KiDS_450_bin1 .OR. &
               ix == tracer_KiDS_450_bin2 .OR. &
               ix == tracer_KiDS_450_bin3 .OR. &
               ix == tracer_KiDS_450_bin4 .OR. &
               ix == tracer_CFHTLenS_Kilbinger2013) THEN
         ! Lensing
         ip = field_matter
      ELSE IF (ix == tracer_CIB_353) THEN
         ip = field_CIB_353
      ELSE IF (ix == tracer_CIB_545) THEN
         ip = field_CIB_545
      ELSE IF (ix == tracer_CIB_857) THEN
         ip = field_CIB_857
      ELSE IF (ix == tracer_galaxies) THEN
         ip = field_galaxies
      ELSE
         STOP 'SET_FIELD_FOR_XPOW: Error, tracer specified incorrectly'
      END IF

   END SUBROUTINE set_field_for_xpow

END MODULE Limber
