MODULE gadget

  USE array_operations
  
CONTAINS

  SUBROUTINE read_gadget(x,v,id,L,om_m,om_v,h,m,a,z,n,infile)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), v(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
    REAL, INTENT(OUT) :: L, z, a, om_m, om_v, h, m
    INTEGER, INTENT(OUT) :: n
    DOUBLE PRECISION :: massarr(6), z8, a8, L8, om_m8, om_v8, h8
    INTEGER :: np(6), np2(6), crap

    !Parameters
    REAL, PARAMETER :: Lunit=1000. !Convert from Gadget2 kpc to Mpc
    REAL, PARAMETER :: Munit=1e10 !Convert from Gadget2 10^10 Msun to Msun

    WRITE(*,*) 'READ_GADGET: Reading in Gadget-2 file: ', trim(infile)

    OPEN(7,file=infile,form='unformatted',status='old')
    READ(7) np, massarr, a8, z8, crap, crap, np2, crap, crap, L8, om_m8, om_v8, h8
    CLOSE(7)

    !Convert Gadget doubles to my reals
    a=real(a8)
    z=real(z8)
    om_m=real(om_m8)
    om_v=real(om_v8)
    h=real(h8)
    L=real(L8)/Lunit

    !Multiply the masses by 1e10 to get in units of M_sun/h
    m=real(massarr(2))*Munit
    WRITE(*,*) 'READ_GADGET: Particle number:', np(2)
    WRITE(*,*) 'READ_GADGET: Which is:', nint(np(2)**(1./3.)), 'cubed.'
    WRITE(*,*) 'READ_GADGET: Particle mass [M_sun/h]:', m
    WRITE(*,*) 'READ_GADGET: Box size [Mpc/h]:', L
    WRITE(*,*) 'READ_GADGET: a:', a
    WRITE(*,*) 'READ_GADGET: z:', z
    WRITE(*,*) 'READ_GADGET: Om_m:', om_m
    WRITE(*,*) 'READ_GADGET: Om_v:', om_v

    !Fix the total number of simulation particles
    n=np(2)

    !Allocate arrays for position, velocity and particle ID number
    ALLOCATE(x(3,n),v(3,n),id(n))

    !Read in the binary data, skip the header line
    OPEN(7,file=infile,form='unformatted',status='old') !CAREFUL - I added status='old' without checking
    READ(7)
    READ(7) x
    READ(7) v
    READ(7) id
    CLOSE(7)

    !kpc -> Mpc conversion!
    x=x/Lunit

    !Change from weird Gadget units to peculiar velocities
    v=v*sqrt(real(a8))

    WRITE(*,*) 'READ_GADGET: Finished reading in file'
    WRITE(*,*)

  END SUBROUTINE read_gadget

  SUBROUTINE write_gadget(x,v,id,L,om_m,om_v,h,m,a,z,n,outfile)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, id(n)
    REAL, INTENT(IN) :: x(3,n), v(3,n)
    REAL, INTENT(IN) :: L, a, z, om_m, om_v, h, m
    CHARACTER(len=*), INTENT(IN) :: outfile
    DOUBLE PRECISION :: massarr(6), z8, a8, L8, om_m8, om_v8, h8, crap8(12)
    INTEGER :: np(6), crapi

    !Parameters
    REAL, PARAMETER :: Lunit=1000. !Convert from Gadget2 kpc to Mpc
    REAL, PARAMETER :: Munit=1e10 !Convert from Gadget2 10^10 Msun to Msun

    WRITE(*,*) 'WRITE_GADGET: Outputting particle data in Gadget2 format: ', trim(outfile)

    np=0
    np(2)=n
    massarr=0.d0
    massarr(2)=m/Munit
    a8=a
    z8=z
    crapi=0
    crap8=0.d0
    om_m8=om_m
    om_v8=om_v
    h8=h
    L8=L*Lunit
    
    WRITE(*,*) 'WRITE_GADGET: Particle number:', n
    WRITE(*,*) 'WRITE_GADGET: Which is:', nint(n**(1./3.)), 'cubed'
    WRITE(*,*) 'WRITE_GADGET: Box size [Mpc/h]:', L
    WRITE(*,*) 'WRITE_GADGET: a:', a
    WRITE(*,*) 'WRITE_GADGET: z:', z
    WRITE(*,*) 'WRITE_GADGET: Particle mass [Msun/h]:', m
    WRITE(*,*) 'WRITE_GADGET: Om_m:', om_m
    WRITE(*,*) 'WRITE_GADGET: Om_v:', om_v

    !x=x*1000.
    !v=v/sqrt(a)

    OPEN(7,file=outfile,form='unformatted',status='replace')
    WRITE(7) np, massarr, a8, z8, crapi, crapi, np, crapi, crapi, L8, om_m8, om_v8, h8, crap8
    WRITE(7) x*Lunit
    WRITE(7) v/sqrt(a)
    WRITE(7) id
    CLOSE(7)

    !Incase these are to be used again outwith the subroutine
    !x=x/1000.
    !v=v*sqrt(a)

    WRITE(*,*) 'WRITE_GADGET: Finished writing file'
    WRITE(*,*)

  END SUBROUTINE write_gadget

  SUBROUTINE read_catalogue(x,v,m,npart,disp,c,env,Dv,rmax,avg_r,rms_r,n,infile)

    USE file_info
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), v(:,:), m(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: disp(:), c(:), env(:), Dv(:), rmax(:), avg_r(:), rms_r(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: npart(:)
    INTEGER, INTENT(OUT) :: n
    INTEGER :: i

    WRITE(*,*) 'READ_CATALOGUE: Reading in catalogue: ', trim(infile)

    n=file_length(infile,verbose=.FALSE.)
    ALLOCATE(x(3,n),v(3,n),m(n),npart(n))
    ALLOCATE(disp(n),c(n),env(n),Dv(n),rmax(n),avg_r(n),rms_r(n))
    
    OPEN(7,file=infile)
    DO i=1,n
       READ(7,*) x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i), m(i), npart(i), disp(i), c(i), env(i), Dv(i), rmax(i), avg_r(i), rms_r(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'READ_CATALOGUE: Min x [Mpc/h]    :', minval(x)
    WRITE(*,*) 'READ_CATALOGUE: Max x [Mpc/h]    :', maxval(x)
    WRITE(*,*) 'READ_CATALOGUE: Min v [km/s]     :', minval(v)
    WRITE(*,*) 'READ_CATALOGUE: Max v [km/s]     :', maxval(v)
    WRITE(*,*) 'READ_CATALOGUE: Min mass [Msun/h]:', minval(m)
    WRITE(*,*) 'READ_CATALOGUE: Max mass [Msun/h]:', maxval(m)
    WRITE(*,*) 'READ_CATALOGUE: Finished reading catalogue file'
    WRITE(*,*)

  END SUBROUTINE read_catalogue

  SUBROUTINE write_catalogue(x,v,m,npart,disp,c,env,Dv,rmax,avg_r,rms_r,n,outfile)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=*), INTENT(IN) :: outfile
    REAL, INTENT(IN) :: x(3,n), v(3,n), m(n)
    REAL, INTENT(IN) :: disp(n), c(n), env(n), Dv(n), rmax(n), avg_r(n), rms_r(n)
    INTEGER, INTENT(IN) :: npart(n)
    INTEGER :: i

    WRITE(*,*) 'WRITE_CATALOGUE: Outputting catalogue'

    OPEN(7,file=outfile)
    DO i=1,n
       WRITE(7,*) x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i), m(i), npart(i), disp(i), c(i), env(i), Dv(i), rmax(i), avg_r(i), rms_r(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'WRITE_CATALOGUE: Finished writing catalogue file'
    WRITE(*,*)

  END SUBROUTINE write_catalogue

!!$  SUBROUTINE write_power(k,D2,nbin,nk,np,L,outfile)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: k(nk), D2(nk), L
!!$    INTEGER, INTENT(IN) :: nbin(nk), nk, np
!!$    CHARACTER(len=256), INTENT(IN) :: outfile
!!$    INTEGER :: i
!!$
!!$    WRITE(*,*) 'WRITE_POWER: Output file: ', trim(outfile)
!!$    OPEN(7,file=outfile)
!!$    DO i=1,nk
!!$       !WRITE(*,*) i
!!$       IF(nbin(i)==0) THEN
!!$          CYCLE
!!$       ELSE
!!$          WRITE(7,*) k(i), D2(i), shot_noise(k(i),L,np), nbin(i)
!!$       END IF
!!$    END DO
!!$    CLOSE(7)
!!$    WRITE(*,*) 'WRITE_POWER: Done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE write_power
  
END MODULE gadget
