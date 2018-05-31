MODULE owls

  IMPLICIT NONE

  !Simulation parameters
  REAL, PARAMETER :: fh=0.752 !Hydrogen mass fraction
  REAL, PARAMETER :: mu=0.61 !Mean molecular weight
  REAL, PARAMETER :: Xe=1.17 !Electron fraction
  REAL, PARAMETER :: Xi=1.08 !Ionisation fraction

  !Physical constants
  REAL, PARAMETER :: msun=1.989e30 !Solar mass in kg
  REAL, PARAMETER :: mp=1.6726e-27 !Proton mass in kg
  REAL, PARAMETER :: Mpc=3.0857e22 !Mpc in m
  REAL, PARAMETER :: cm=0.01 !cm in m
  REAL, PARAMETER :: eV=1.60218e-12 !eV in erg
  
CONTAINS

   SUBROUTINE read_mccarthy(x,m,n,infile)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), m(:)
    INTEGER, INTENT(OUT) :: n
    REAL, PARAMETER :: mfac=1e10

    WRITE(*,*) 'READ_MCCARTHY: Reading in binary file: ', TRIM(infile)

    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    CLOSE(7)

    !In case the array is empty, but actually Ian has n=1 set (e.g., UVB_stars)
    IF(n==1) THEN
       n=0
    END IF

    WRITE(*,*) 'READ_MCCARTHY: Particle number:', n
    WRITE(*,*) 'READ_MCCARTHY: Which is ~', NINT(n**(1./3.)), 'cubed.'

    ALLOCATE(x(3,n),m(n))

    IF(n .NE. 0) THEN

       !Need to read in 'n' again with stream access
       OPEN(7,file=infile,form='unformatted',access='stream',status='old')
       READ(7) n
       READ(7) m
       READ(7) x
       CLOSE(7)

       m=m*mfac
       
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

    WRITE(*,*)

  END SUBROUTINE read_mccarthy

  SUBROUTINE read_mccarthy_gas(x,m,kT,nh,n,infile)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), m(:), nh(:), kT(:)
    REAL, ALLOCATABLE :: ep(:)
    INTEGER, INTENT(OUT) :: n
    
    REAL, PARAMETER :: mfac=1e10 !Convert mass to Solar masses

    !Read in the binary file
    WRITE(*,*) 'READ_MCCARTHY: Reading in binary file: ', TRIM(infile)
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    CLOSE(7)
    WRITE(*,*) 'READ_MCCARTHY: Particle number:', n
    WRITE(*,*) 'READ_MCCARTHY: Which is ~', NINT(n**(1./3.)), 'cubed.'

    !Allocate arrays for quantities in the file
    ALLOCATE(x(3,n),m(n),ep(n),nh(n))
    
    !Need to read in 'n' again with stream access
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) n
    READ(7) m
    READ(7) x
    READ(7) ep !electron pressure in erg/cm^3
    READ(7) nh !hydrogen number density in /cm^3
    CLOSE(7)

    !Convert masses into Solar masses
    m=m*mfac

    !Convert the electron pressure [erg/cm^3] and hydrogen density [#/cm^3] into kT
    !Units of kT will be [erg]
    ALLOCATE(kT(n))
    kT=((Xe+Xi)/Xe)*(ep/nh)*mu*fh

    !Convert internal energy from erg to eV
    kT=kT/eV

    !Deallocate the electron pressure array
    DEALLOCATE(ep)

    !Write information to the screen
    WRITE(*,*) 'READ_MCCARTHY: Minimum particle mass [Msun/h]:', MINVAL(m)
    WRITE(*,*) 'READ_MCCARTHY: Maximum particle mass [Msun/h]:', MAXVAL(m)
    WRITE(*,*) 'READ_MCCARTHY: Minimum x coordinate [Mpc/h]:', MINVAL(x(1,:))
    WRITE(*,*) 'READ_MCCARTHY: Minimum x coordinate [Mpc/h]:', MAXVAL(x(1,:))
    WRITE(*,*) 'READ_MCCARTHY: Minimum y coordinate [Mpc/h]:', MINVAL(x(2,:))
    WRITE(*,*) 'READ_MCCARTHY: Minimum y coordinate [Mpc/h]:', MAXVAL(x(2,:))
    WRITE(*,*) 'READ_MCCARTHY: Minimum z coordinate [Mpc/h]:', MINVAL(x(3,:))
    WRITE(*,*) 'READ_MCCARTHY: Minimum z coordinate [Mpc/h]:', MAXVAL(x(3,:))
    WRITE(*,*) 'READ_MCCARTHY: Minimum internal energy [eV]:', MINVAL(kT)
    WRITE(*,*) 'READ_MCCARTHY: Maximum internal energy [eV]:', MAXVAL(kT)
    WRITE(*,*) 'READ_MCCARTHY: Minimum hydrogen number density [cm^-3]:', MINVAL(nh)
    WRITE(*,*) 'READ_MCCARTHY: Maximum hydrogen number density [cm^-3]:', MAXVAL(nh)
    WRITE(*,*) 'READ_MCCARTHY: Finished reading in file'
    WRITE(*,*)

  END SUBROUTINE read_mccarthy_gas

  SUBROUTINE convert_kT_to_electron_pressure(kT,nh,mass,n,L,h,m)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: kT(n)
    REAL, INTENT(IN) :: mass(n), nh(n), L, h
    INTEGER, INTENT(IN) :: n, m
    REAL :: V
    DOUBLE PRECISION :: units, kT_dble(n)
    
    LOGICAL, PARAMETER :: apply_nh_cut=.TRUE. !Apply a cut in hydrogen density
    REAL, PARAMETER :: nh_cut=0.1 !Cut in the hydrogen number density [cm^-3]

    !Exclude gas that is sufficiently dense to not be ionised and be forming stars
    IF(apply_nh_cut) CALL exclude_nh(nh_cut,kT,nh,n)

    !Use double precision because all the constants are dreadful
    kT_dble=kT  

    !Convert to particle internal energy that needs to be mapped to grid [eV*Msun]
    kT_dble=kT_dble*(mass/mu)*Xe/(Xe+Xi)
    !kT=kT*(mass/mu)*Xe/(Xe+Xi)

    !Cell volume in [(Mpc/h)^3]
    V=(L/REAL(m))**3

    !Cell volume in [Mpc^3]
    V=V/h**3

    !This is now pressure [eV/Mpc^3 * Msun]
    kT_dble=kT_dble/V

    !Convert units of pressure [eV/cm^3]
    units=msun
    units=units/mp
    units=units/(Mpc/cm)
    units=units/(Mpc/cm)
    units=units/(Mpc/cm)
    kT_dble=kT_dble*units

    !Go back to single precision
    kT=REAL(kT_dble)

  END SUBROUTINE convert_kT_to_electron_pressure
  
!!$  SUBROUTINE read_mccarthy_gas_old(x,m,ep,nh,n,infile)
!!$
!!$    IMPLICIT NONE
!!$    CHARACTER(len=256), INTENT(IN) :: infile
!!$    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), m(:), ep(:), nh(:)
!!$    INTEGER, INTENT(OUT) :: n
!!$    REAL, PARAMETER :: mfac=1e10
!!$
!!$    WRITE(*,*) 'READ_MCCARTHY: Reading in binary file: ', TRIM(infile)
!!$
!!$    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
!!$    READ(7) n
!!$    CLOSE(7)
!!$
!!$    WRITE(*,*) 'READ_MCCARTHY: Particle number:', n
!!$    WRITE(*,*) 'READ_MCCARTHY: Which is ~', NINT(n**(1./3.)), 'cubed.'
!!$
!!$    ALLOCATE(x(3,n),m(n),ep(n),nh(n))
!!$    
!!$    !Need to read in 'n' again with stream access
!!$    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
!!$    READ(7) n
!!$    READ(7) m
!!$    READ(7) x
!!$    READ(7) ep
!!$    READ(7) nh
!!$    CLOSE(7)
!!$
!!$    m=m*mfac
!!$
!!$    WRITE(*,*) 'READ_MCCARTHY: Minimum and maximum mass [Msun/h]:', MINVAL(m), MAXVAL(m)
!!$    WRITE(*,*) 'READ_MCCARTHY: Minimum and maximum x coordinate [Mpc/h]:', MINVAL(x(1,:)), MAXVAL(x(1,:))
!!$    WRITE(*,*) 'READ_MCCARTHY: Minimum and maximum y coordinate [Mpc/h]:', MINVAL(x(2,:)), MAXVAL(x(2,:))
!!$    WRITE(*,*) 'READ_MCCARTHY: Minimum and maximum z coordinate [Mpc/h]:', MINVAL(x(3,:)), MAXVAL(x(3,:))
!!$    WRITE(*,*) 'READ_MCCARTHY: Minimum and maximum electron pressure [erg/cm^3]:', MINVAL(ep), MAXVAL(ep)
!!$    WRITE(*,*) 'READ_MCCARTHY: Minimum and maximum hydrogen number density [cm^-3]:', MINVAL(nh), MAXVAL(nh)
!!$    WRITE(*,*) 'READ_MCCARTHY: Finished reading in file'
!!$
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE read_mccarthy_gas_old

  SUBROUTINE write_mccarthy(x,m,n,outfile)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: outfile
    REAL, INTENT(IN) :: x(3,n), m(n)
    INTEGER, INTENT(IN) :: n
    REAL, PARAMETER :: mfac=1e10

    WRITE(*,*) 'WRITE_MCCARTHY: Outputting binary file: ', TRIM(outfile)

    WRITE(*,*) 'WRITE_MCCARTHY: Particle number:', n
    WRITE(*,*) 'WRITE_MCCARTHY: Which is ~', NINT(n**(1./3.)), 'cubed.'

    !Need to read in 'n' again with stream access
    OPEN(7,file=outfile,form='unformatted',access='stream',status='replace')
    WRITE(7) n
    WRITE(7) m/mfac
    WRITE(7) x
    CLOSE(7)
    
    WRITE(*,*) 'WRITE_MCCARTHY: Finished writing file'
    WRITE(*,*)

  END SUBROUTINE write_mccarthy

  SUBROUTINE exclude_nh(nhcut,ep,nh,n)

    !Set the electron pressure to zero of any particle that has nh > nhcut
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
