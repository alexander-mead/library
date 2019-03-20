MODULE owls

  USE constants
  
  IMPLICIT NONE

  ! BAHAMAS simulation parameters
  REAL, PARAMETER :: fh=0.752 ! Hydrogen mass fraction
  REAL, PARAMETER :: mup=0.61 ! Mean particle mass relative to proton
  REAL, PARAMETER :: Xe=1.17 ! Number of electrons per hydrogen (X_{e/H} in my notation)
  REAL, PARAMETER :: Xi=1.08 ! Number of ions per hydrogen (X_{i/H}; note that all gas particles are either electrons or ions)
  REAL, PARAMETER :: mfac=1e10 ! Mass conversion factor to get Msun/h
  REAL, PARAMETER :: eV_erg=eV*1e7 ! eV in ergs

  LOGICAL, PARAMETER :: apply_nh_cut=.TRUE. ! Apply a cut in hydrogen density
  REAL, PARAMETER :: nh_cut=0.1 ! Cut in the hydrogen number density [cm^-3] gas denser than this is not ionised
  
CONTAINS

   SUBROUTINE read_mccarthy(x,m,n,infile)

    ! Read in a McCarthy format particle data file
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), m(:)
    INTEGER, INTENT(OUT) :: n
    LOGICAL :: lexist

    ! Write to screen
    WRITE(*,*) 'READ_MCCARTHY: Reading in binary file: ', trim(infile)

    ! Open the file using stream
    INQUIRE(file=infile, exist=lexist)
    IF(.NOT. lexist) STOP 'READ_MCCARTHY: Error, input file does not exist'
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    CLOSE(7)

    ! In case the array is empty, but actually Ian has n=1 set (e.g., UVB_stars)
    IF(n==1) THEN
       n=0
    END IF

    ! Write information to screen
    WRITE(*,*) 'READ_MCCARTHY: Particle number:', n
    WRITE(*,*) 'READ_MCCARTHY: Which is ~', nint(n**(1./3.)), 'cubed.'

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
       WRITE(*,*) 'READ_MCCARTHY: Minimum particle mass [Msun/h]:', minval(m)
       WRITE(*,*) 'READ_MCCARTHY: Maximum particle mass [Msun/h]:', maxval(m)
       WRITE(*,*) 'READ_MCCARTHY: Total particle mass [Msun/h]:', sum(m)
       WRITE(*,*) 'READ_MCCARTHY: Minimum x coordinate [Mpc/h]:', minval(x(1,:))
       WRITE(*,*) 'READ_MCCARTHY: Maximum x coordinate [Mpc/h]:', maxval(x(1,:))
       WRITE(*,*) 'READ_MCCARTHY: Minimum y coordinate [Mpc/h]:', minval(x(2,:))
       WRITE(*,*) 'READ_MCCARTHY: Maximum y coordinate [Mpc/h]:', maxval(x(2,:))
       WRITE(*,*) 'READ_MCCARTHY: Minimum z coordinate [Mpc/h]:', minval(x(3,:))
       WRITE(*,*) 'READ_MCCARTHY: Maximum z coordinate [Mpc/h]:', maxval(x(3,:))
       WRITE(*,*) 'READ_MCCARTHY: Finished reading in file'

    END IF

    ! Final white space
    WRITE(*,*)

  END SUBROUTINE read_mccarthy

  SUBROUTINE read_mccarthy_gas(x,m,kT,nh,n,infile)

    ! Read in a McCarthy format gas file
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), m(:), nh(:), kT(:)
    REAL, ALLOCATABLE :: ep(:)
    INTEGER, INTENT(OUT) :: n
    LOGICAL :: lexist   
    REAL :: mue
    
    ! Read in the binary file
    INQUIRE(file=infile, exist=lexist)
    IF(.NOT. lexist) STOP 'READ_MCCARTHY: Error, input file does not exist'
    WRITE(*,*) 'READ_MCCARTHY_GAS: Reading in binary file: ', trim(infile)
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    CLOSE(7)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Particle number:', n
    WRITE(*,*) 'READ_MCCARTHY_GAS: Which is ~', nint(n**(1./3.)), 'cubed.'

    ! Allocate arrays for quantities in the file
    ALLOCATE(x(3,n),m(n),ep(n),nh(n))
    
    ! Need to read in 'n' again with stream access
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    READ(7) m
    READ(7) x
    READ(7) ep ! physical electron pressure for the particle [erg/cm^3]
    READ(7) nh ! hydrogen number density for the partcle in [/cm^3]
    CLOSE(7)

    ! Convert masses into Solar masses
    m=m*mfac

    ! Write information to the screen
    WRITE(*,*) 'READ_MCCARTHY_GAS: Calculating kT from physical electron pressure'
    WRITE(*,*) 'READ_MCCARTHY_GAS: Note that the electron pressure is *not* comoving'
    WRITE(*,*) 'READ_MCCARTHY_GAS: Using numbers appropriate for BAHAMAS'
    WRITE(*,*) 'READ_MCCARTHY_GAS: YH:', fh
    WRITE(*,*) 'READ_MCCARTHY_GAS: mu_p:', mup
    WRITE(*,*) 'READ_MCCARTHY_GAS: Xe:', Xe
    WRITE(*,*) 'READ_MCCARTHY_GAS: Xi:', Xi

    ! Calculate and write the 'particle mass per free electron: mu_e'
    mue=mup*(Xe+Xi)/Xe
    WRITE(*,*) 'READ_MCCARTHY_GAS: mu_e:', mue
    
    ! Convert the physical electron pressure [erg/cm^3] and hydrogen density [#/cm^3] into kT [erg]
    ! This is the temperature of gas particles (equal for all species)
    ! Temperature is neither comoving nor physical
    ALLOCATE(kT(n))
    !kT=((Xe+Xi)/Xe)*(ep/nh)*mu*fh
    kT=(ep/nh)*mue*fh

    ! Convert internal energy from erg to eV
    kT=kT/eV_erg

    ! Deallocate the physical electron pressure array
    DEALLOCATE(ep)

    ! Write information to the screen
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum particle mass [Msun/h]:', minval(m)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum particle mass [Msun/h]:', maxval(m)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum x coordinate [Mpc/h]:', minval(x(1,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum x coordinate [Mpc/h]:', maxval(x(1,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum y coordinate [Mpc/h]:', minval(x(2,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum y coordinate [Mpc/h]:', maxval(x(2,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum z coordinate [Mpc/h]:', minval(x(3,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum z coordinate [Mpc/h]:', maxval(x(3,:))
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum internal energy [eV]:', minval(kT)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum internal energy [eV]:', maxval(kT)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Minimum hydrogen number density [cm^-3]:', minval(nh)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Maximum hydrogen number density [cm^-3]:', maxval(nh)
    WRITE(*,*) 'READ_MCCARTHY_GAS: Finished reading in file'
    WRITE(*,*)

  END SUBROUTINE read_mccarthy_gas

  SUBROUTINE convert_kT_to_comoving_electron_pressure(kT,nh,m,n,L,h)

    ! This routine converts the input particle internal energy kT [eV] to electron pressure, Pe [eV/cm^3]
    ! Note very well that Pe will be the contribution to the total pressure in the volume per particle
    ! CARE: I removed factors of 'm', pressure is now contribution to entire volume, rather than the contribution per mesh cell
    USE constants
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: kT(n) ! particle internal energy [eV]
    REAL, INTENT(IN) :: nh(n)    ! hydrogen number density [/cm^3]
    REAL, INTENT(IN) :: m(n)     ! hydrodynamic particle mass [Msun/h]
    INTEGER, INTENT(IN) :: n     ! total number of particles
    REAL, INTENT(IN) :: L        ! Box size [Mpc/h]
    REAL, INTENT(IN) :: h        ! Hubble parameter (necessary because pressure will be in eV/cm^3 without h factors) 
    REAL :: mue, V
    DOUBLE PRECISION :: units, kT_dble(n)

    ! Exclude gas that is sufficiently dense to not be ionised and be forming stars
    IF(apply_nh_cut) CALL exclude_nh(nh_cut,kT,nh,n)

    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Converting kT to comoving electron pressure'
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Using numbers appropriate for BAHAMAS'
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Note that this is COMOVING'
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Y_H:', fh
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: mu_p:', mup
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Xe:', Xe
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Xi:', Xi

    mue=mup*(Xe+Xi)/Xe
    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: mu_e:', mue

    ! Use double precision because all the constants are dreadful 
    kT_dble=kT ! [eV]
    
    ! Convert to particle internal energy that can be mapped to grid
    !kT_dble=kT_dble*(m/mu)*Xe/(Xe+Xi) ! [eV*Msun]
    kT_dble=kT_dble*(m/mue) ! [eV*Msun]

    ! Comoving volume
    V=(L/h)**3 ! [(Mpc)^3]

    ! This is now comoving electron pressure
    kT_dble=kT_dble/V ! [Msun*eV/Mpc^3]

    ! Convert units of comoving electron pressure
    ! Note that there are no h factors here
    units=msun
    units=units/mp ! Divide out proton mass here [eV/Mpc^3]
    units=units/(Mpc/0.01)
    units=units/(Mpc/0.01)
    units=units/(Mpc/0.01)
    kT_dble=kT_dble*units ! [eV/cm^3]

    ! Go back to single precision
    kT=real(kT_dble) ! [eV/cm^3]

    WRITE(*,*) 'CONVERT_KT_TO_ELECTRON_PRESSURE: Done'
    WRITE(*,*)

  END SUBROUTINE convert_kT_to_comoving_electron_pressure

  SUBROUTINE write_mccarthy(x,m,n,outfile)

    ! Write a particle data file using McCarthy format
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: outfile
    REAL, INTENT(IN) :: x(3,n), m(n)
    INTEGER, INTENT(IN) :: n

    WRITE(*,*) 'WRITE_MCCARTHY: Outputting binary file: ', trim(outfile)

    WRITE(*,*) 'WRITE_MCCARTHY: Particle number:', n
    WRITE(*,*) 'WRITE_MCCARTHY: Which is ~', nint(n**(1./3.)), 'cubed.'

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
