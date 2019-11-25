MODULE cosmic_emu_stuff

   USE file_info
   USE array_operations
   USE cosmology_functions
   USE interpolate

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: get_Cosmic_Emu_power
   PUBLIC :: get_Franken_Emu_power
   PUBLIC :: get_Mira_Titan_power

CONTAINS

   SUBROUTINE get_Cosmic_Emu_power(k, P, n, z, cosm, rebin)

      USE logical_operations
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:), P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, INTENT(IN) :: rebin
      INTEGER :: i
      REAL, ALLOCATABLE :: k2(:), P2(:)
      CHARACTER(len=256) :: crap
      REAL :: h
      REAL, PARAMETER :: kmin = 1e-2 ! Minimum k if rebinnning
      REAL, PARAMETER :: kmax = 7.   ! Minimum k if rebinnning
      INTEGER, PARAMETER :: nk = 128 ! Number of k values if rebinning
      INTEGER, PARAMETER :: nh = 10  ! Length of header
      INTEGER, PARAMETER :: h_li = 8 ! Line that h in on
      CHARACTER(len=256), PARAMETER :: params = 'emu_params.txt'
      CHARACTER(len=256), PARAMETER :: output = 'emu_power.dat'
      CHARACTER(len=256), PARAMETER :: exe = '/Users/Mead/Physics/CosmicEmu/emu.exe'
      REAL, PARAMETER :: eps_h = 3e-3

      ! Remove previous parameter and power file
      CALL SYSTEM('rm '//trim(params))
      CALL SYSTEM('rm '//trim(output))

      ! Write a new parameter file
      OPEN (7, file=params)
      WRITE (7, fmt='(A20,7F10.5)') trim(output), cosm%om_m*cosm%h**2, cosm%om_b*cosm%h**2, cosm%n, cosm%sig8, cosm%w, z
      CLOSE (7)

      ! Run emu
      CALL SYSTEM(trim(exe)//' '//trim(params))! > /dev/null')

      ! Get length of emu file
      n = file_length(output)
      n = n-nh
      IF (ALLOCATED(k)) DEALLOCATE (k)
      IF (ALLOCATED(P)) DEALLOCATE (P)
      ALLOCATE (k(n), P(n))

      ! Write useful things to screen
      WRITE (*, *) 'GET_COSMIC_EMU_POWER: z:', z
      WRITE (*, *) 'GET_COSMIC_EMU_POWER: P(k) file length:', n
      WRITE (*, *) 'GET_COSMIC_EMU_POWER: Reading in P(k): ', trim(output)

      ! Read in data file
      OPEN (7, file=output)
      DO i = 1-nh, n
         IF (i == h_li-nh) THEN
            READ (7, *) crap, crap, crap, crap, crap, crap, crap, crap, h
            WRITE (*, *) 'GET_COSMIC_EMU_POWER: CMB derived h:', h
            WRITE (*, *) 'GET_COSMIC_EMU_POWER: Cosmology h:', cosm%h
            WRITE (*, *) 'GET_COSMIC_EMU_POWER: h ratio:', cosm%h/h
            IF (.NOT. requal(h, cosm%h, eps_h)) STOP 'GET_COSMIC_EMU_POWER: Error, h values differ'
         ELSE IF (i < 1) THEN
            READ (7, *)
         ELSE
            READ (7, *) k(i), P(i)
         END IF
      END DO
      CLOSE (7)

      !Convert k to k/h
      k = k/cosm%h

      ! Rebin on a log-linear axis if desired
      IF (rebin) THEN
         CALL fill_array(log(kmin), log(kmax), k2, nk)
         k2 = exp(k2)
         ALLOCATE (P2(nk))
         CALL interpolate_array(log(k), log(P), n, log(k2), P2, nk, 3, 3, 2)
         P2 = exp(P2)
         DEALLOCATE (k, P)
         ALLOCATE (k(nk), P(nk))
         k = k2
         P = P2
         n = nk
         DEALLOCATE (k2, P2)
      END IF

      ! Done
      WRITE (*, *) 'GET_COSMIC_EMU_POWER: Done'
      WRITE (*, *)

   END SUBROUTINE get_Cosmic_Emu_power

   SUBROUTINE get_Franken_Emu_power(k, P, n, z, cosm, rebin)

      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:), P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, INTENT(IN) :: rebin
      INTEGER :: i
      REAL, ALLOCATABLE :: k2(:), P2(:)
      REAL, PARAMETER :: kmin = 1e-2 ! Minimum k if rebinnning
      REAL, PARAMETER :: kmax = 7.   ! Minimum k if rebinnning
      INTEGER, PARAMETER :: nk = 128 ! Number of k values if rebinning
      INTEGER, PARAMETER :: nh = 5   ! Length of header
      CHARACTER(len=256), PARAMETER :: params = 'emu_params.txt'
      CHARACTER(len=256), PARAMETER :: output = 'emu_power.dat'
      CHARACTER(len=256), PARAMETER :: exe = '/Users/Mead/Physics/FrankenEmu/emu.exe'

      ! Remove previous parameter and power file
      !CALL SYSTEM('rm emu_params.txt')
      !CALL SYSTEM('rm emu_power.dat')
      CALL SYSTEM('rm '//trim(params))
      CALL SYSTEM('rm '//trim(output))

      ! Write a new parameter file
      OPEN (7, file=params)
      WRITE (7, fmt='(A20,7F10.5)') trim(output), cosm%om_m*cosm%h**2, cosm%om_b*cosm%h**2, cosm%n, -cosm%w, cosm%sig8, cosm%h, z
      CLOSE (7)

      ! Run emu
      !CALL SYSTEM('/Users/Mead/Physics/FrankenEmu/emu.exe emu_params.txt')! > /dev/null')
      CALL SYSTEM(trim(exe)//' '//trim(params))! > /dev/null')

      ! Get length of emu file
      n = file_length(output, verbose=.FALSE.)
      n = n-nh
      IF (ALLOCATED(k)) DEALLOCATE (k)
      IF (ALLOCATED(P)) DEALLOCATE (P)
      ALLOCATE (k(n), P(n))

      ! Write useful things to screen
      WRITE (*, *) 'GET_FRANKENEMU_POWER: z:', z
      WRITE (*, *) 'GET_FRANKENEMU_POWER: P(k) file length:', n
      WRITE (*, *) 'GET_FRANKENEMU_POWER: Reading in P(k): ', trim(output)

      ! Read in data file
      OPEN (7, file=output)
      DO i = 1-nh, n
         IF (i < 1) THEN
            READ (7, *)
         ELSE
            READ (7, *) k(i), P(i)
         END IF
      END DO
      CLOSE (7)

      !Convert k to k/h
      k = k/cosm%h

      ! Rebin on a log-linear axis if desired
      IF (rebin) THEN
         CALL fill_array(log(kmin), log(kmax), k2, nk)
         k2 = exp(k2)
         ALLOCATE (P2(nk))
         CALL interpolate_array(log(k), log(P), n, log(k2), P2, nk, 3, 3, 2)
         P2 = exp(P2)
         DEALLOCATE (k, P)
         ALLOCATE (k(nk), P(nk))
         k = k2
         P = P2
         n = nk
         DEALLOCATE (k2, P2)
      END IF

      ! Done
      WRITE (*, *) 'GET_FRANKENEMU_POWER: Done'
      WRITE (*, *)

   END SUBROUTINE get_Franken_Emu_power

   SUBROUTINE get_Mira_Titan_power(k, P, n, z, cosm, rebin)

      USE constants
      IMPLICIT NONE
      REAL, ALLOCATABLE, INTENT(OUT) :: k(:)
      REAL, ALLOCATABLE, INTENT(OUT) :: P(:)
      INTEGER, INTENT(OUT) :: n
      REAL, INTENT(IN) :: z
      TYPE(cosmology), INTENT(IN) :: cosm
      LOGICAL, INTENT(IN) :: rebin
      CHARACTER(len=256) :: output
      INTEGER :: i, nk
      REAL, ALLOCATABLE :: k2(:), P2(:)
      REAL :: kmin, kmax

      CALL SYSTEM('rm xstar.dat')
      CALL SYSTEM('rm EMU0.txt')

      OPEN (7, file='xstar.dat')
    WRITE (7, *) (cosm%Om_m*cosm%h**2), (cosm%Om_b*cosm%h**2), cosm%sig8, cosm%h, cosm%n, cosm%w, cosm%wa, (cosm%om_nu*cosm%h**2), z
      CLOSE (7)

      CALL SYSTEM('/Users/Mead/Physics/MiraTitan/P_tot/emu.exe')
      output = 'EMU0.txt'

      n = file_length(output, verbose=.FALSE.)
      IF (ALLOCATED(k)) DEALLOCATE (k)
      IF (ALLOCATED(P)) DEALLOCATE (P)
      ALLOCATE (k(n), P(n))

      WRITE (*, *) 'GET_MIRA_TITAN_POWER: z:', z
      WRITE (*, *) 'GET_MIRA_TITAN_POWER: P(k) file length:', n
      WRITE (*, *) 'GET_MIRA_TITAN_POWER: Reading in P(k)'

      OPEN (7, file=output)
      DO i = 1, n
         READ (7, *) k(i), P(i)
      END DO
      CLOSE (7)

      !Convert P(k) to Delta^2(k)
      P = P*(k**3)*4.*pi/(2.*pi)**3

      !Convert k to k/h
      k = k/cosm%h

      ! Rebin on a log-linear axis
      IF (rebin) THEN
         kmin = 1e-2
         kmax = 7.
         nk = 128
         CALL fill_array(log(kmin), log(kmax), k2, nk)
         k2 = exp(k2)
         ALLOCATE (P2(nk))
         CALL interpolate_array(log(k), log(P), n, log(k2), P2, nk, 3, 3, 2)
         P2 = exp(P2)
         DEALLOCATE (k, P)
         ALLOCATE (k(nk), P(nk))
         k = k2
         P = P2
         n = nk
         DEALLOCATE (k2, P2)
      END IF

      WRITE (*, *) 'GET_MIRA_TITAN_POWER: Done'
      WRITE (*, *)

   END SUBROUTINE get_Mira_Titan_power

END MODULE cosmic_emu_stuff
