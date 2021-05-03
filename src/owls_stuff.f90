MODULE owls_stuff

   USE constants
   USE basic_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_mccarthy
   PUBLIC :: read_mccarthy_gas
   PUBLIC :: convert_kT_to_comoving_electron_pressure
   PUBLIC :: write_mccarthy
   PUBLIC :: BAHAMAS_snapshot
   PUBLIC :: BAHAMAS_scale_factor
   PUBLIC :: BAHAMAS_scale_factors

   PUBLIC :: nh_cut
   PUBLIC :: mue

   ! BAHAMAS simulation parameters
   REAL, PARAMETER :: fh = 0.752      ! Hydrogen mass fraction
   REAL, PARAMETER :: mup = 0.61      ! Mean particle mass relative to proton
   REAL, PARAMETER :: Xeh = 1.17      ! Number of electrons per hydrogen (X_{e/H} in my notation)
   REAL, PARAMETER :: Xih = 1.08      ! Number of ions per hydrogen (X_{i/H}; all gas particles are either electrons or ions)
   REAL, PARAMETER :: mfac = 1e10     ! Mass conversion factor to get Msun/h
   REAL, PARAMETER :: eV_erg = eV*1e7 ! eV in ergs

   ! BAHAMAS derived parameters
   REAL, PARAMETER :: mue = mup*(Xeh+Xih)/Xeh ! Mean particle mass per electron relative to proton

   ! Cuts to remove neutral gas
   LOGICAL, PARAMETER :: apply_nh_cut = .TRUE. ! Apply a cut in hydrogen density?
   REAL, PARAMETER :: nh_cut = 0.1             ! Cut in the hydrogen number density [#/cm^3] gas denser than this is not ionised

CONTAINS

   SUBROUTINE read_mccarthy(x, m, n, infile)

      ! Read in a McCarthy format particle data file
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: m(:)
      INTEGER, INTENT(OUT) :: n
      CHARACTER(len=*), INTENT(IN) :: infile
      LOGICAL :: lexist

      ! Open the file using stream
      WRITE (*, *) 'READ_MCCARTHY: Reading in binary file: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) STOP 'READ_MCCARTHY: Error, input file does not exist'
      OPEN (7, file=infile, form='unformatted', access='stream', status='old')
      READ (7) n
      CLOSE (7)

      ! In case the array is empty, but actually Ian has n=1 set (e.g., UVB_stars)
      IF (n == 1) THEN
         n = 0
      END IF

      ! Write information to screen
      WRITE (*, *) 'READ_MCCARTHY: Particle number:', n
      WRITE (*, *) 'READ_MCCARTHY: Which is ~', nint(n**(1./3.)), 'cubed.'

      ! Allocate arrays
      ALLOCATE (x(3, n), m(n))

      IF (n .NE. 0) THEN

         ! Need to read in 'n' again with stream access
         OPEN (7, file=infile, form='unformatted', access='stream', status='old')
         READ (7) n
         READ (7) m
         READ (7) x
         CLOSE (7)

         ! Multiply by mass factor
         m = m*mfac

         ! Write information to screen
         WRITE (*, *) 'READ_MCCARTHY: Minimum particle mass [Msun/h]:', minval(m)
         WRITE (*, *) 'READ_MCCARTHY: Maximum particle mass [Msun/h]:', maxval(m)
         WRITE (*, *) 'READ_MCCARTHY: Total particle mass [Msun/h]:', sum(m)
         WRITE (*, *) 'READ_MCCARTHY: Minimum particle x coordinate [Mpc/h]:', minval(x(1, :))
         WRITE (*, *) 'READ_MCCARTHY: Maximum particle x coordinate [Mpc/h]:', maxval(x(1, :))
         WRITE (*, *) 'READ_MCCARTHY: Minimum particle y coordinate [Mpc/h]:', minval(x(2, :))
         WRITE (*, *) 'READ_MCCARTHY: Maximum particle y coordinate [Mpc/h]:', maxval(x(2, :))
         WRITE (*, *) 'READ_MCCARTHY: Minimum particle z coordinate [Mpc/h]:', minval(x(3, :))
         WRITE (*, *) 'READ_MCCARTHY: Maximum particle z coordinate [Mpc/h]:', maxval(x(3, :))
         WRITE (*, *) 'READ_MCCARTHY: Finished reading in file'

      END IF

      ! Final white space
      WRITE (*, *)

   END SUBROUTINE read_mccarthy

   SUBROUTINE write_mccarthy(x, m, n, outfile)

      ! Write a particle data file using McCarthy format
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: x(3, n)
      REAL, INTENT(IN) :: m(n)
      CHARACTER(len=*), INTENT(IN) :: outfile

      WRITE (*, *) 'WRITE_MCCARTHY: Outputting binary file: ', trim(outfile)

      WRITE (*, *) 'WRITE_MCCARTHY: Particle number:', n
      WRITE (*, *) 'WRITE_MCCARTHY: Which is ~', nint(n**(1./3.)), 'cubed.'

      OPEN (7, file=outfile, form='unformatted', access='stream', status='replace')
      WRITE (7) n
      WRITE (7) m/mfac
      WRITE (7) x
      CLOSE (7)

      WRITE (*, *) 'WRITE_MCCARTHY: Finished writing file'
      WRITE (*, *)

   END SUBROUTINE write_mccarthy

   SUBROUTINE read_mccarthy_gas(x, m, kT, rho, n, infile)

      ! Read in a McCarthy format gas file
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:, :) ! Particle positions [Mpc/h]
      REAL, ALLOCATABLE, INTENT(OUT) :: m(:)   ! Particle mass [Msun/h]
      REAL, ALLOCATABLE, INTENT(OUT) :: kT(:)  ! Particle internal energy [eV]
      REAL, ALLOCATABLE, INTENT(OUT) :: rho(:) ! Particle SPH physical density [mp/cm^3]
      INTEGER, INTENT(OUT) :: n                ! Total number of gas particles
      CHARACTER(len=*), INTENT(IN) :: infile   ! Input file name
      REAL, ALLOCATABLE :: nh(:), ep(:)
      LOGICAL :: lexist

      ! Read in the binary file
      WRITE (*, *) 'READ_MCCARTHY_GAS: Reading in binary file: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) STOP 'READ_MCCARTHY: Error, input file does not exist'
      OPEN (7, file=infile, form='unformatted', access='stream', status='old')
      READ (7) n
      CLOSE (7)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Particle number:', n
      WRITE (*, *) 'READ_MCCARTHY_GAS: Which is ~', nint(n**(1./3.)), 'cubed.'

      ! Allocate arrays for quantities in the file
      ALLOCATE (x(3, n), m(n), ep(n), nh(n))

      ! Need to read in 'n' again with stream access
      OPEN (7, file=infile, form='unformatted', access='stream', status='old')
      READ (7) n
      READ (7) m
      READ (7) x
      READ (7) ep ! physical electron pressure for the particle [erg/cm^3]
      READ (7) nh ! physical hydrogen number density for the partcle in [/cm^3]
      CLOSE (7)

      ! Convert particle masses
      m = m*mfac ! [Msun/h]

      ! Write information to the screen
      WRITE (*, *) 'READ_MCCARTHY_GAS: Calculating kT from physical electron pressure'
      WRITE (*, *) 'READ_MCCARTHY_GAS: Note that the electron pressure calcuated here is *not* comoving'
      WRITE (*, *) 'READ_MCCARTHY_GAS: Using numbers appropriate for BAHAMAS'
      WRITE (*, *) 'READ_MCCARTHY_GAS: Hydrogen mass fraction: f_H, Y_H:', fh
      WRITE (*, *) 'READ_MCCARTHY_GAS: Mean particle mass: mu_p [m_p]:', mup
      WRITE (*, *) 'READ_MCCARTHY_GAS: Number of electrons per hydrogen: X_e/X_H:', Xeh
      WRITE (*, *) 'READ_MCCARTHY_GAS: Number of ions per hydrogen: X_i/X_H:', Xih
      WRITE (*, *) 'READ_MCCARTHY_GAS: Mean particle mass per electron: mu_e [m_p]:', mue

      ! Convert the physical hydrogen number density into a physical particle mass density [mp/cm^3]
      ! Note that these densities are physical *not* comoving
      ALLOCATE (rho(n))
      rho = nh/fh ! Convert physical hydrogen number density [#/cm^3] to physsical particle SPH density [mp/cm^3]
      DEALLOCATE (nh) ! Deallocate the physical electron pressure array

      ! Convert the physical electron pressure [erg/cm^3] and hydrogen density [#/cm^3] into kT [erg]
      ! This is the temperature of gas particles (equal for all species)
      ! Temperature is neither comoving nor physical
      ALLOCATE (kT(n))
      kT = (ep/rho)*mue
      kT = kT/eV_erg ! Convert internal energy from erg to eV
      DEALLOCATE (ep) ! Deallocate the physical electron pressure array

      ! Write information to the screen
      WRITE (*, *) 'READ_MCCARTHY_GAS: Minimum particle mass [Msun/h]:', minval(m)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Maximum particle mass [Msun/h]:', maxval(m)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Total mass of all particles [Msun/h]:', sum(m)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Minimum particle x coordinate [Mpc/h]:', minval(x(1, :))
      WRITE (*, *) 'READ_MCCARTHY_GAS: Maximum particle x coordinate [Mpc/h]:', maxval(x(1, :))
      WRITE (*, *) 'READ_MCCARTHY_GAS: Minimum particle y coordinate [Mpc/h]:', minval(x(2, :))
      WRITE (*, *) 'READ_MCCARTHY_GAS: Maximum particle y coordinate [Mpc/h]:', maxval(x(2, :))
      WRITE (*, *) 'READ_MCCARTHY_GAS: Minimum particle z coordinate [Mpc/h]:', minval(x(3, :))
      WRITE (*, *) 'READ_MCCARTHY_GAS: Maximum particle z coordinate [Mpc/h]:', maxval(x(3, :))
      WRITE (*, *) 'READ_MCCARTHY_GAS: Minimum particle internal energy [eV]:', minval(kT)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Maximum particle internal energy [eV]:', maxval(kT)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Minimum particle SPH density [mp/cm^3]:', minval(rho)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Maximum particle SPH density [mp/cm^3]:', maxval(rho)
      WRITE (*, *) 'READ_MCCARTHY_GAS: Finished reading in file'
      WRITE (*, *)

   END SUBROUTINE read_mccarthy_gas

   SUBROUTINE convert_kT_to_comoving_electron_pressure(kT, rho, m, n, L, h)

      ! This routine converts the input particle internal energy kT [eV] to comoving electron pressure, Pe [eV/cm^3]
      ! Note very well that Pe will be the contribution to the total pressure in the volume per particle
      ! CARE: I removed factors of 'm', pressure is now contribution to entire volume, rather than the contribution per mesh cell
      INTEGER, INTENT(IN) :: n     ! total number of particles
      REAL, INTENT(INOUT) :: kT(n) ! particle internal energy [eV], output as electron pressure [eV/cm^3]
      REAL, INTENT(IN) :: rho(n)   ! physical gas particle SPH density [mp/cm^3]
      REAL, INTENT(IN) :: m(n)     ! hydrodynamic particle mass [Msun/h] 
      REAL, INTENT(IN) :: L        ! box size [Mpc/h]
      REAL, INTENT(IN) :: h        ! Hubble parameter (necessary because pressure will be in eV/cm^3 without h factors)
      REAL :: V
      REAL :: units, kT_dble(n) ! CARE: Previously this was double precision
      INTEGER, PARAMETER :: scheme = 3

      ! Exclude gas that is sufficiently dense to not be ionised and be forming stars
      IF (apply_nh_cut) CALL exclude_nh(nh_cut, kT, rho, n)

      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Converting kT to comoving electron pressure'
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Using numbers appropriate for BAHAMAS'
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Note that this is COMOVING'
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Y_H:', fh
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: mu_p [mp]:', mup
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: X_e/X_H:', Xeh
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: X_i/X_H:', Xih
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: mu_e [mp]:', mue

      ! Use double precision because all the constants are dreadful
      kT_dble = kT ! [eV]

      IF (scheme == 1 .OR. scheme == 2) THEN

         ! Convert to particle internal energy that can be mapped to grid
         IF (scheme == 1) THEN
            kT_dble = kT_dble*(m/mue) ! [eV*Msun/h] BUG: I think there should be a factor of h here due to Msun/h units
         ELSE IF (scheme == 2) THEN
            kT_dble = kT_dble*(m/(h*mue)) ! [eV*Msun]
         ELSE
            STOP 'CONVERT_KT_TO_ELECTRON_PRESSURE: Error, scheme specified incorrectly'
         END IF

         ! Comoving volume
         V = (L/h)**3 ! [(Mpc)^3]

         ! This is now comoving electron pressure
         kT_dble = kT_dble/V ! [eV*Msun/Mpc^3]

         ! Convert units of comoving electron pressure
         ! Note that there are no h factors here
         units = msun             ! [kg]
         units = units/mp         ! Divide out proton mass here [dimensionless]
         units = units/(Mpc/0.01) ! Necessary to do this in stages due to overflow [kg/cm^1]
         units = units/(Mpc/0.01) ! Necessary to do this in stages due to overflow [kg/cm^2]
         units = units/(Mpc/0.01) ! Necessary to do this in stages due to overflow [kg/cm^3]
         kT_dble = kT_dble*units  ! [eV/cm^3]

      ELSE IF (scheme == 3) THEN

         ! Comoving volume
         V = L**3 ! [(Mpc)^3]

         ! Big blob of units (see notes)
         units = (Msun*(h**2)*(0.01)**3)
         units = units/Mpc
         units = units/Mpc
         units = units/Mpc

         ! Final quantity
         kT_dble = units*(m/(mue*mp))*kT_dble/V ! [eV/cm^3]

      ELSE
         STOP 'CONVERT_KT_TO_ELECTRON_PRESSURE: Error, scheme specified incorrectly'
      END IF

      ! Go back to single precision
      kT = real(kT_dble) ! [eV/cm^3]

      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Minimum electron pressure [eV/cm^3]: ', minval(kT)
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Maximum electron pressure [eV/cm^3]: ', maxval(kT)
      WRITE (*, *) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Done'
      WRITE (*, *)

   END SUBROUTINE convert_kT_to_comoving_electron_pressure

   SUBROUTINE exclude_nh(nhcut, kT, rho, n)

      ! Set the gas particle internal energy to zero for high-density particles that have nh > nhcut
      INTEGER, INTENT(IN) :: n      ! Total number of particles
      REAL, INTENT(IN) :: nhcut     ! Cut to impose on hydrogen number density [#/cm^3]
      REAL, INTENT(INOUT) :: kT(n)  ! Gas particle internal energy [eV]
      REAL, INTENT(IN) :: rho(n)    ! Gap particle SPH density [mp/cm^3]
      INTEGER :: i
      REAL :: rhocut

      ! Convert the hydrogen number density cut into a cut on particle SPH density
      rhocut = nhcut/fh

      DO i = 1, n
         IF (rho(i) > rhocut) kT(i) = 0.
      END DO

   END SUBROUTINE exclude_nh

   CHARACTER(len=32) FUNCTION BAHAMAS_snapshot(z)

      ! Set the redshift corresponding to the snapsho
      REAL, INTENT(IN) :: z
      REAL, PARAMETER :: eps = 1e-3
   
      IF (requal(z, 0.0, eps)) THEN
         BAHAMAS_snapshot = 'snap32'
      ELSE IF (requal(z, 0.125, eps)) THEN
         BAHAMAS_snapshot = 'snap31'
      ELSE IF (requal(z, 0.25, eps)) THEN
         BAHAMAS_snapshot = 'snap30'
      ELSE IF (requal(z, 0.375, eps)) THEN
         BAHAMAS_snapshot = 'snap29'
      ELSE IF (requal(z, 0.5, eps)) THEN
         BAHAMAS_snapshot = 'snap28'
      ELSE IF (requal(z, 0.75, eps)) THEN
         BAHAMAS_snapshot = 'snap27'
      ELSE IF (requal(z, 1.0, eps)) THEN
         BAHAMAS_snapshot = 'snap26'
      ELSE IF (requal(z, 1.25, eps)) THEN
         BAHAMAS_snapshot = 'snap25'
      ELSE IF (requal(z, 1.5, eps)) THEN
         BAHAMAS_snapshot = 'snap24'
      ELSE IF (requal(z, 1.75, eps)) THEN
         BAHAMAS_snapshot = 'snap23'
      ELSE IF (requal(z, 2.0, eps)) THEN
         BAHAMAS_snapshot = 'snap22'
      ELSE IF (requal(z, 2.25, eps)) THEN
         BAHAMAS_snapshot = 'snap21'
      ELSE IF (requal(z, 2.5, eps)) THEN
         BAHAMAS_snapshot = 'snap20'
      ELSE IF (requal(z, 2.75, eps)) THEN
         BAHAMAS_snapshot = 'snap19'
      ELSE IF (requal(z, 3.0, eps)) THEN
         BAHAMAS_snapshot = 'snap18'
      ELSE
         WRITE (*, *) 'BAHAMAS_SNAPSHOT: z', z
         STOP 'BAHAMAS_SNAPSHOT: Error, redshift specified incorrectly'
      END IF

   END FUNCTION BAHAMAS_snapshot

   ELEMENTAL REAL FUNCTION BAHAMAS_scale_factor(isnap)

      USE cosmology_functions
      INTEGER, INTENT(IN) :: isnap
      REAL :: z

      IF (isnap == 18) THEN
         z = 3.000
      ELSE IF (isnap == 19) THEN
         z = 2.750
      ELSE IF (isnap == 20) THEN
         z = 2.500
      ELSE IF (isnap == 21) THEN
         z = 2.250
      ELSE IF (isnap == 22) THEN
         z = 2.000
      ELSE IF (isnap == 23) THEN
         z = 1.750
      ELSE IF (isnap == 24) THEN
         z = 1.500
      ELSE IF (isnap == 25) THEN
         z = 1.250
      ELSE IF (isnap == 26) THEN
         z = 1.000
      ELSE IF (isnap == 27) THEN
         z = 0.750
      ELSE IF (isnap == 28) THEN
         z = 0.500
      ELSE IF (isnap == 29) THEN
         z = 0.375
      ELSE IF (isnap == 30) THEN
         z = 0.250
      ELSE IF (isnap == 31) THEN
         z = 0.125
      ELSE IF (isnap == 32) THEN
         z = 0.000
      !ELSE
      !   STOP 'BAHAMAS_scale_factors'
      END IF

      BAHAMAS_scale_factor = scale_factor_z(z)

   END FUNCTION BAHAMAS_scale_factor

   SUBROUTINE BAHAMAS_scale_factors(a, n)

      ! Returns an array of all the BAHAMAS scale factors for a given set of n
      USE cosmology_functions    
      REAL, ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: n
      INTEGER :: snaps(n)

      IF (ALLOCATED(a)) DEALLOCATE(a)
      ALLOCATE(a(n))

      IF (n == 4) THEN
         snaps = [22, 26, 28, 32]
      ELSE IF (n == 11) THEN
         snaps = [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
      ELSE IF (n == 15) THEN
         snaps = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
      ELSE
         STOP 'BAHAMAS_SCALE_FACTORS: Error, n must be either 4, 11 or 15'
      END IF

      a = BAHAMAS_scale_factor(snaps)

   END SUBROUTINE BAHAMAS_scale_factors

END MODULE owls_stuff
