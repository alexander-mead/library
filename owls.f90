MODULE owls

  IMPLICIT NONE

  ! BAHAMAS simulation parameters
  REAL, PARAMETER :: fh=0.752 ! Hydrogen mass fraction
  REAL, PARAMETER :: mu=0.61 ! Mean molecular weight relative to proton
  REAL, PARAMETER :: Xe=1.17 ! Electron fraction (number of electrons per hydrogen)
  REAL, PARAMETER :: Xi=1.08 ! Ion fraction (number of ionisations per hydrogen)
  
CONTAINS

   SUBROUTINE read_mccarthy(x,m,n,infile)

    ! Read in a McCarthy format particle data file
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), m(:)
    INTEGER, INTENT(OUT) :: n
    REAL, PARAMETER :: mfac=1e10

    ! Write to screen
    WRITE(*,*) 'READ_MCCARTHY: Reading in binary file: ', TRIM(infile)

    ! Open the file using stream
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    CLOSE(7)

    ! In case the array is empty, but actually Ian has n=1 set (e.g., UVB_stars)
    IF(n==1) THEN
       n=0
    END IF

    ! Write information to screen
    WRITE(*,*) 'READ_MCCARTHY: Particle number:', n
    WRITE(*,*) 'READ_MCCARTHY: Which is ~', NINT(n**(1./3.)), 'cubed.'

    ! Allocate arrays
    ALLOCATE(x(3,n),m(n))

    IF(n .NE. 0) THEN

       ! Need to read in 'n' again with stream access
       OPEN(7,file=infile,form='unformatted',access='stream',status='old')
       READ(7) n
       READ(7) m
       READ(7) x
       CLOSE(7)

       ! Multiply by mass factor
       m=m*mfac

       ! Write information to screen
       WRITE(*,*) 'READ_MCCARTHY: Minimum particle mass [Msun/h]:', MINVAL(m)
       WRITE(*,*) 'READ_MCCARTHY: Maximum particle mass [Msun/h]:', MAXVAL(m)
       WRITE(*,*) 'READ_MCCARTHY: Total particle mass [Msun/h]:', SUM(m)
       WRITE(*,*) 'READ_MCCARTHY: Minimum x coordinate [Mpc/h]:', MINVAL(x(1,:))
       WRITE(*,*) 'READ_MCCARTHY: Minimum x coordinate [Mpc/h]:', MAXVAL(x(1,:))
       WRITE(*,*) 'READ_MCCARTHY: Minimum y coordinate [Mpc/h]:', MINVAL(x(2,:))
       WRITE(*,*) 'READ_MCCARTHY: Minimum y coordinate [Mpc/h]:', MAXVAL(x(2,:))
       WRITE(*,*) 'READ_MCCARTHY: Minimum z coordinate [Mpc/h]:', MINVAL(x(3,:))
       WRITE(*,*) 'READ_MCCARTHY: Minimum z coordinate [Mpc/h]:', MAXVAL(x(3,:))
       WRITE(*,*) 'READ_MCCARTHY: Finished reading in file'

    END IF

    ! Final white space
    WRITE(*,*)

  END SUBROUTINE read_mccarthy

  SUBROUTINE read_mccarthy_gas(x,m,kT,nh,n,infile)

    ! Read in a McCarthy format gas file
    USE constants
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), m(:), nh(:), kT(:)
    REAL, ALLOCATABLE :: ep(:)
    INTEGER, INTENT(OUT) :: n
    
    REAL, PARAMETER :: mfac=1e10 ! Convert mass to Solar masses
    REAL, PARAMETER :: eV_erg=eV*1e7 ! eV in ergs

    ! Read in the binary file
    WRITE(*,*) 'READ_MCCARTHY_GAS: Reading in binary file: ', TRIM(infile)
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    CLOSE(7)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Particle number:', n
    WRITE(*,*) 'READ_MCCARTHY_GAS: Which is ~', NINT(n**(1./3.)), 'cubed.'

    ! Allocate arrays for quantities in the file
    ALLOCATE(x(3,n),m(n),ep(n),nh(n))
    
    ! Need to read in 'n' again with stream access
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    READ(7) m
    READ(7) x
    READ(7) ep ! physical electron pressure for the particle in erg/cm^3
    READ(7) nh ! hydrogen number density for the partcle in /cm^3
    CLOSE(7)

    ! Convert masses into Solar masses
    m=m*mfac

    ! Write information to the screen
    WRITE(*,*) 'READ_MCCARTHY_GAS: Calculating kT from physical electron pressure'
    WRITE(*,*) 'READ_MCCARTHY_GAS: Note that the electron pressure is *not* comoving'
    WRITE(*,*) 'READ_MCCARTHY_GAS: Using numbers appropriate for BAHAMAS'
    WRITE(*,*) 'READ_MCCARTHY_GAS: YH:', fh
    WRITE(*,*) 'READ_MCCARTHY_GAS: mu_H:', mu
    WRITE(*,*) 'READ_MCCARTHY_GAS: Xe:', Xe
    WRITE(*,*) 'READ_MCCARTHY_GAS: Xi:', Xi
    
    ! Convert the physical electron pressure [erg/cm^3] and hydrogen density [#/cm^3] into kT [erg]
    ! This is the temperature of gas particles (equal for all species)
    ! Temperature is neither comoving nor physical
    ALLOCATE(kT(n))
    kT=((Xe+Xi)/Xe)*(ep/nh)*mu*fh

    ! Convert internal energy from erg to eV
    kT=kT/eV_erg

    ! Deallocate the physical electron pressure array
    DEALLOCATE(ep)

    ! Write information to the screen
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum particle mass [Msun/h]:', MINVAL(m)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum particle mass [Msun/h]:', MAXVAL(m)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum x coordinate [Mpc/h]:', MINVAL(x(1,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum x coordinate [Mpc/h]:', MAXVAL(x(1,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum y coordinate [Mpc/h]:', MINVAL(x(2,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum y coordinate [Mpc/h]:', MAXVAL(x(2,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum z coordinate [Mpc/h]:', MINVAL(x(3,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum z coordinate [Mpc/h]:', MAXVAL(x(3,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum internal energy [eV]:', MINVAL(kT)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum internal energy [eV]:', MAXVAL(kT)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum hydrogen number density [cm^-3]:', MINVAL(nh)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum hydrogen number density [cm^-3]:', MAXVAL(nh)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Finished reading in file'
    WRITE(*,*)

  END SUBROUTINE read_mccarthy_gas

  SUBROUTINE convert_kT_to_comoving_electron_pressure(kT,nh,mass,n,L,h,m)

    ! kT is particle internal energy input in units of [eV]
    ! nh is hydrogen number density [/cm^3]
    ! mass is particle mass in units of [Msun]
    ! n is the total number of particles
    ! L is the box size in units of [Mpc/h]
    ! h is the dimensionless hubble parameter
    ! m is the mesh size onto which the pressure will be binned
    USE constants
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: kT(n)
    REAL, INTENT(IN) :: mass(n), nh(n), L, h
    INTEGER, INTENT(IN) :: n, m
    REAL :: V
    DOUBLE PRECISION :: units, kT_dble(n)
    
    LOGICAL, PARAMETER :: apply_nh_cut=.TRUE. ! Apply a cut in hydrogen density
    REAL, PARAMETER :: nh_cut=0.1 ! Cut in the hydrogen number density [cm^-3] gas denser than this is not ionised

    ! Exclude gas that is sufficiently dense to not be ionised and be forming stars
    IF(apply_nh_cut) CALL exclude_nh(nh_cut,kT,nh,n)

    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Converting kT to comoving electron pressure'
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Using numbers appropriate for BAHAMAS'
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Note that this is COMOVING'
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Y_H:', fh
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: mu_H:', mu
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Xe:', Xe
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Xi:', Xi

    ! Use double precision because all the constants are dreadful 
    kT_dble=kT ! [eV]
    
    ! Convert to particle internal energy that needs to be mapped to grid
    kT_dble=kT_dble*(mass/mu)*Xe/(Xe+Xi) ! [eV*Msun]

    ! Comoving cell volume
    V=(L/REAL(m))**3 ! [(Mpc/h)^3]
    V=V/h**3 ! remove h factors [Mpc^3]

    ! This is now comoving electron pressure
    kT_dble=kT_dble/V ! [Msun*eV/Mpc^3]

    ! Convert units of comoving electron pressure
    ! Note that there are no h factors here
    units=msun
    units=units/mp ! Divide out proton mass here [eV/Mpc^3]
    units=units/(Mpc/cm)
    units=units/(Mpc/cm)
    units=units/(Mpc/cm)
    kT_dble=kT_dble*units! [eV/cm^3]

    ! Go back to single precision
    kT=REAL(kT_dble)! [eV/cm^3]

    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Done'
    WRITE(*,*)

  END SUBROUTINE convert_kT_to_comoving_electron_pressure

  SUBROUTINE write_mccarthy(x,m,n,outfile)

    ! Write a particle data file using McCarthy format
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: outfile
    REAL, INTENT(IN) :: x(3,n), m(n)
    INTEGER, INTENT(IN) :: n
    REAL, PARAMETER :: mfac=1e10

    WRITE(*,*) 'WRITE_MCCARTHY: Outputting binary file: ', TRIM(outfile)

    WRITE(*,*) 'WRITE_MCCARTHY: Particle number:', n
    WRITE(*,*) 'WRITE_MCCARTHY: Which is ~', NINT(n**(1./3.)), 'cubed.'

    OPEN(7,file=outfile,form='unformatted',access='stream',status='replace')
    WRITE(7) n
    WRITE(7) m/mfac
    WRITE(7) x
    CLOSE(7)
    
    WRITE(*,*) 'WRITE_MCCARTHY: Finished writing file'
    WRITE(*,*)

  END SUBROUTINE write_mccarthy

  SUBROUTINE exclude_nh(nhcut,ep,nh,n)

    ! Set the electron pressure to zero of any high-density particle that has nh > nhcut
    IMPLICIT NONE
    REAL, INTENT(IN) :: nhcut, nh(n)
    REAL, INTENT(INOUT) :: ep(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    DO i=1,n
       IF(nh(i)>nhcut) ep(i)=0.
    END DO
    
  END SUBROUTINE exclude_nh

END MODULE owls
