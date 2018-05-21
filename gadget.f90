MODULE gadget

  USE array_operations
  
CONTAINS

  SUBROUTINE read_gadget(x,v,id,L,om_m,om_v,h,m,a,z,n,infile)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), v(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
    REAL, INTENT(OUT) :: L, z, a, om_m, om_v, h, m
    INTEGER, INTENT(OUT) :: n
    DOUBLE PRECISION :: massarr(6), z8, a8, L8, om_m8, om_v8, h8
    INTEGER :: np(6), np2(6), crap

    !Parameters
    REAL, PARAMETER :: Lunit=1000. !Convert from Gadget2 kpc to Mpc
    REAL, PARAMETER :: Munit=1e10 !Convert from Gadget2 10^10 Msun to Msun

    WRITE(*,*) 'READ_GADGET: Reading in Gadget-2 file: ', TRIM(infile)

    OPEN(7,file=infile,form='unformatted',status='old')
    READ(7) np, massarr, a8, z8, crap, crap, np2, crap, crap, L8, om_m8, om_v8, h8
    CLOSE(7)

    !Convert Gadget doubles to my reals
    a=REAL(a8)
    z=REAL(z8)
    om_m=REAL(om_m8)
    om_v=REAL(om_v8)
    h=REAL(h8)

    !Multiply the masses by 1e10 to get in units of M_sun/h
    m=REAL(massarr(2))*Munit
    WRITE(*,*) 'READ_GADGET: Particle number:', np(2)
    WRITE(*,*) 'READ_GADGET: Which is:', nint(np(2)**(1./3.)), 'cubed.'
    WRITE(*,*) 'READ_GADGET: Particle mass (M_sun/h):', m
    WRITE(*,*) 'READ_GADGET: Box size (Mpc/h):', L
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
    L=REAL(L8)/Lunit

    !Change from weird Gadget units to peculiar velocities
    v=v*sqrt(REAL(a8))

    WRITE(*,*) 'READ_GADGET: Finished reading in file'
    WRITE(*,*)

  END SUBROUTINE read_gadget

  SUBROUTINE write_gadget(x,v,id,L,om_m,om_v,h,m,a,z,n,outfile)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, id(n)
    REAL, INTENT(IN) :: x(3,n), v(3,n)
    REAL, INTENT(IN) :: L, a, z, om_m, om_v, h, m
    CHARACTER(len=256) :: outfile
    DOUBLE PRECISION :: massarr(6), z8, a8, L8, om_m8, om_v8, h8, crap8(12)
    INTEGER :: np(6), crapi

    !Parameters
    REAL, PARAMETER :: Lunit=1000. !Convert from Gadget2 kpc to Mpc
    REAL, PARAMETER :: Munit=1e10 !Convert from Gadget2 10^10 Msun to Msun

    STOP 'BUG CHECK. I changed the WRITE statements to have the unit conversions in them (e.g., x*Lunit)'

    WRITE(*,*) 'WRITE_GADGET: Outputting particle data in Gadget2 format: ', TRIM(outfile)

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
    WRITE(*,*) 'WRITE_GADGET: Box size (Mpc/h):', L
    WRITE(*,*) 'WRITE_GADGET: a:', a
    WRITE(*,*) 'WRITE_GADGET: z:', z
    WRITE(*,*) 'WRITE_GADGET: Particle mass:', m
    WRITE(*,*) 'WRITE_GADGET: Om_m:', om_m
    WRITE(*,*) 'WRITE_GADGET: Om_v:', om_v

    !x=x*1000.
    !v=v/sqrt(a)

    OPEN(7,file=outfile,form='unformatted')
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
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:), v(:,:), m(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: disp(:), c(:), env(:), Dv(:), rmax(:), avg_r(:), rms_r(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: npart(:)
    INTEGER, INTENT(OUT) :: n
    INTEGER :: i

    WRITE(*,*) 'READ_CATALOGUE: Reading in catalogue: ', TRIM(infile)

    n=file_length(infile)
    ALLOCATE(x(3,n),v(3,n),m(n),npart(n))
    ALLOCATE(disp(n),c(n),env(n),Dv(n),rmax(n),avg_r(n),rms_r(n))
    
    OPEN(7,file=infile)
    DO i=1,n
       READ(7,*) x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i), m(i), npart(i), disp(i), c(i), env(i), Dv(i), rmax(i), avg_r(i), rms_r(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'READ_CATALOGUE: Min x [Mpc/h]    :', MINVAL(x)
    WRITE(*,*) 'READ_CATALOGUE: Max x [Mpc/h]    :', MAXVAL(x)
    WRITE(*,*) 'READ_CATALOGUE: Min v [km/s]     :', MINVAL(v)
    WRITE(*,*) 'READ_CATALOGUE: Max v [km/s]     :', MAXVAL(v)
    WRITE(*,*) 'READ_CATALOGUE: Min mass [Msun/h]:', MINVAL(m)
    WRITE(*,*) 'READ_CATALOGUE: Max mass [Msun/h]:', MAXVAL(m)
    WRITE(*,*) 'READ_CATALOGUE: Finished reading catalogue file'
    WRITE(*,*)

  END SUBROUTINE read_catalogue

  SUBROUTINE write_catalogue(x,v,m,npart,disp,c,env,Dv,rmax,avg_r,rms_r,n,outfile)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=256), INTENT(IN) :: outfile
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

  SUBROUTINE replace(x,n,L)

    !Ensures/enforces periodicity by cycling particles round that may have strayed
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(3,n)
    REAL, INTENT(IN) :: L
    INTEGER :: i, j
    INTEGER, INTENT(IN) :: n

    DO i=1,n
       DO j=1,3
          IF(x(j,i)>L)   x(j,i)=x(j,i)-L
          IF(x(j,i)<=0.) x(j,i)=x(j,i)+L
       END DO
    END DO

  END SUBROUTINE replace

  SUBROUTINE particle_bin(x,n,L,w,d,m,ibin)

    !Bin particle properties onto a mesh, summing as you go
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, m, ibin
    REAL, INTENT(INOUT) :: d(m,m,m)
    REAL, INTENT(IN) :: x(3,n), L, w(n)

    IF(ibin==1) THEN
       CALL NGP(x,n,L,w,d,m)
    ELSE IF(ibin==2) THEN
       CALL CIC(x,n,L,w,d,m)
    ELSE
       STOP 'PARTICLE_BIN: Error, ibin not specified correctly'
    END IF

  END SUBROUTINE particle_bin

  SUBROUTINE particle_bin_average(x,n,L,w,d,m,ibin)

    !Bin particle properties onto a mesh, averaging properties over cells
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, m, ibin
    REAL, INTENT(INOUT) :: d(m,m,m)
    REAL, INTENT(IN) :: x(3,n), L, w(n)
    REAL :: number(m,m,m), one(n)
    INTEGER :: i, j, k, sum

    !Need an array of ones (wasteful for memory, could use pointer maybe?)
    one=1.

    !Call the binning twice, first to bin particle property and second to count
    IF(ibin==1) THEN
       CALL NGP(x,n,L,w,d,m)
       CALL NGP(x,n,L,one,number,m)
    ELSE IF(ibin==2) THEN
       CALL CIC(x,n,L,w,d,m)
       CALL CIC(x,n,L,one,number,m)
    ELSE
       STOP 'PARTICLE_BIN_AVERAGE: Error, ibin not specified correctly'
    END IF

    !Now loop over all elements of the field and average over the number of contributions
    !Not sure how this will work for CIC binning if cell only gets a small contribution from one particle
    sum=0
    WRITE(*,*) 'PARTICLE_BIN_AVERAGE: Averaging over particles in cells'
    DO k=1,m
       DO j=1,m
          DO i=1,m
             IF(number(i,j,k)==0.) THEN
                sum=sum+1
             ELSE
                d(i,j,k)=d(i,j,k)/number(i,j,k)
             END IF
          END DO
       END DO
    END DO
    WRITE(*,*) 'PARTICLE_BIN_AVERAGE: Numer of empty cells:', n
    WRITE(*,*) 'PARTICLE_BIN_AVERAGE: Fraction of empty cells', REAL(sum)/REAL(m**3)
    WRITE(*,*) 'PARTICLE_BIN_AVERAGE: Done'
    WRITE(*,*)

  END SUBROUTINE particle_bin_average

  SUBROUTINE NGP(x,n,L,w,d,m)

    !Nearest-grid-point binning routine
    USE statistics
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, m
    INTEGER :: i
    REAL, INTENT(INOUT) :: d(m,m,m)
    REAL, INTENT(IN) :: x(3,n), L, w(n)
    INTEGER :: ix, iy, iz

    WRITE(*,*) 'NGP: Binning particles and creating field'
    WRITE(*,*) 'NGP: Cells:', m

    d=0.

    DO i=1,n

       ix=CEILING(x(1,i)*float(m)/L)
       iy=CEILING(x(2,i)*float(m)/L)
       iz=CEILING(x(3,i)*float(m)/L)

       IF(ix>m .OR. ix<1) THEN
          WRITE(*,*) 'x:', i, x(1,i)
          STOP 'NGP: Warning, point outside box'
       END IF

       IF(iy>m .OR. iy<1) THEN
          WRITE(*,*) 'y:', i, x(2,i)
          STOP 'NGP: Warning, point outside box'
       END IF

       IF(iz>m .OR. iz<1) THEN
          WRITE(*,*) 'z:', i, x(3,i)
          STOP 'NGP: Warning, point outside box'
       END IF

       d(ix,iy,iz)=d(ix,iy,iz)+w(i)

    END DO

    WRITE(*,*) 'NGP: Average:', mean(d,m)
    WRITE(*,*) 'NGP: RMS:', sqrt(variance(d,m))
    WRITE(*,*) 'NGP: Minimum:', MINVAL(REAL(d))
    WRITE(*,*) 'NGP: Maximum:', MAXVAL(REAL(d))
    WRITE(*,*) 'NGP: Binning complete'
    WRITE(*,*)

  END SUBROUTINE NGP

  SUBROUTINE CIC(x,n,L,w,d,m)

    !Cloud-in-cell binning routine
    USE statistics
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, m
    REAL, INTENT(IN) :: x(3,n), L, w(n)
    REAL, INTENT(INOUT) :: d(m,m,m)
    INTEGER :: ix, iy, iz, ixn, iyn, izn
    INTEGER :: i
    REAL :: dx, dy, dz

    WRITE(*,*) 'CIC: Binning particles and creating field'
    WRITE(*,*) 'CIC: Cells:', m

    d=0.

    DO i=1,n

       ix=CEILING(x(1,i)*float(m)/L)
       iy=CEILING(x(2,i)*float(m)/L)
       iz=CEILING(x(3,i)*float(m)/L)

       IF(ix>m .OR. ix<1) THEN
          WRITE(*,*) 'x:', i, x(1,i)
          STOP 'CIC: Warning, point outside box'
       END IF

       IF(iy>m .OR. iy<1) THEN
          WRITE(*,*) 'y:', i, x(2,i)
          STOP 'CIC: Warning, point outside box'
       END IF

       IF(iz>m .OR. iz<1) THEN
          WRITE(*,*) 'z:', i, x(3,i)
          STOP 'CIC: Warning, point outside box'
       END IF

       !dx, dy, dz in box units
       dx=(x(1,i)/L)*float(m)-(float(ix)-0.5)
       dy=(x(2,i)/L)*float(m)-(float(iy)-0.5)
       dz=(x(3,i)/L)*float(m)-(float(iz)-0.5)

       IF(dx>=0.) THEN
          ixn=ix+1
          IF(ixn>m) ixn=1
       ELSE
          ixn=ix-1
          dx=-dx  
          IF(ixn<1) ixn=m    
       END IF

       IF(dy>=0.) THEN
          iyn=iy+1
          IF(iyn>m) iyn=1
       ELSE
          iyn=iy-1
          dy=-dy
          IF(iyn<1) iyn=m
       END IF

       IF(dz>=0.) THEN
          izn=iz+1
          IF(izn>m) izn=1
       ELSE
          izn=iz-1
          dz=-dz
          IF(izn<1) izn=m
       END IF

       d(ix,iy,iz)=d(ix,iy,iz)+(1.-dx)*(1.-dy)*(1.-dz)*w(i)

       d(ix,iy,izn)=d(ix,iy,izn)+(1.-dx)*(1.-dy)*dz*w(i)
       d(ix,iyn,iz)=d(ix,iyn,iz)+(1.-dx)*dy*(1.-dz)*w(i)
       d(ixn,iy,iz)=d(ixn,iy,iz)+dx*(1.-dy)*(1.-dz)*w(i)

       d(ix,iyn,izn)=d(ix,iyn,izn)+(1.-dx)*dy*dz*w(i)
       d(ixn,iyn,iz)=d(ixn,iyn,iz)+dx*dy*(1.-dz)*w(i)
       d(ixn,iy,izn)=d(ixn,iy,izn)+dx*(1.-dy)*dz*w(i)

       d(ixn,iyn,izn)=d(ixn,iyn,izn)+dx*dy*dz*w(i)

    END DO

    WRITE(*,*) 'CIC: Average:', mean(d,m)
    WRITE(*,*) 'CIC: RMS:', sqrt(variance(d,m))
    WRITE(*,*) 'CIC: Minimum:', MINVAL(REAL(d))
    WRITE(*,*) 'CIC: Maximum:', MAXVAL(REAL(d))
    WRITE(*,*) 'CIC: Binning complete'
    WRITE(*,*)

  END SUBROUTINE CIC

  SUBROUTINE CIC2D(x,y,L,d,m)

    !This could probably be usefully combined with CIC somehow
    IMPLICIT NONE
    INTEGER :: ix, iy, ixn, iyn
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: m
    REAL, INTENT(IN) :: x(:), y(:), L
    REAL, ALLOCATABLE, INTENT(INOUT) :: d(:,:)
    REAL :: dx, dy

    WRITE(*,*) 'CIC: Binning particles and creating density field'
    WRITE(*,*) 'CIC: Cells:', m

    IF(ALLOCATED(d)) DEALLOCATE(d)
    ALLOCATE(d(m,m))
    d=0.

    n=SIZE(x)

    DO i=1,n

       ix=CEILING(x(i)*float(m)/L)
       iy=CEILING(y(i)*float(m)/L)

       IF(ix>m .OR. ix<1) THEN
          WRITE(*,*) 'x:', i, x(i), ix
          STOP 'CIC: Warning, point outside box'
       END IF

       IF(iy>m .OR. iy<1) THEN
          WRITE(*,*) 'y:', i, y(i), iy
          STOP 'CIC: Warning, point outside box'
       END IF

       !dx, dy in cell units, away from cell centre
       dx=(x(i)/L)*float(m)-(float(ix)-0.5)
       dy=(y(i)/L)*float(m)-(float(iy)-0.5)

       !Find CIC weights in x
       IF(dx>0.) THEN
          ixn=ix+1
          IF(ixn>m) ixn=1
       ELSE
          ixn=ix-1
          dx=-dx  
          IF(ixn<1) ixn=m    
       END IF

       !Find CIC weights in y
       IF(dy>=0.) THEN
          iyn=iy+1
          IF(iyn>m) iyn=1
       ELSE
          iyn=iy-1
          dy=-dy
          IF(iyn<1) iyn=m
       END IF

       !Carry out CIC binning
       d(ix,iy)=d(ix,iy)+(1.-dx)*(1.-dy)
       d(ix,iyn)=d(ix,iyn)+(1.-dx)*dy
       d(ixn,iy)=d(ixn,iy)+dx*(1.-dy)
       d(ixn,iyn)=d(ixn,iyn)+dx*dy

    END DO

    WRITE(*,*) 'CIC: Binning complete'
    WRITE(*,*)

  END SUBROUTINE CIC2D
  
  SUBROUTINE SOD(x,n,L,w,d,m,Rs)

    !Spherical density routine
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, ALLOCATABLE, INTENT(INOUT) :: d(:,:,:)
    REAL, ALLOCATABLE :: xc(:,:,:,:)
    REAL, INTENT(IN) :: x(3,n), L, w(n), Rs
    REAL :: r, dx, dy, dz
    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz, kx, ky, kz
    INTEGER, INTENT(IN) :: m

    WRITE(*,*) 'SOD: Binning particles and creating density field'
    WRITE(*,*) 'SOD: Cells:', m

    ALLOCATE(d(m,m,m),xc(3,m,m,m))
    d=0.
    
    !Fill array with sphere centre positions
    DO i=1,m
       DO j=1,m
          DO k=1,m
             xc(1,i,j,k)=L*(float(i)-0.5)/float(m)
             xc(2,i,j,k)=L*(float(j)-0.5)/float(m)
             xc(3,i,j,k)=L*(float(k)-0.5)/float(m)
          END DO
       END DO
    END DO

    !Loop over all particles and assign to spheres
    DO i=1,n

       !WRITE(*,*) 'Particle coordinates:', i, x(1,i), x(2,i), x(3,i)

       !Coordinate of nearest mesh cell
       ix=CEILING(x(1,i)*float(m)/L)
       iy=CEILING(x(2,i)*float(m)/L)
       iz=CEILING(x(3,i)*float(m)/L)

       !WRITE(*,*) 'Mesh cell coordinates:', ix, iy, iz

       !Check that particles are within box
       IF(ix>m .OR. ix<1) THEN
          WRITE(*,*) 'x:', i, x(1,i)
          STOP 'SOD: Error, particle outside box'
       END IF

       !Check that particles are within box
       IF(iy>m .OR. iy<1) THEN
          WRITE(*,*) 'y:', i, x(2,i)
          STOP 'SOD: Error, particle outside box'
       END IF

       !Check that particles are within box
       IF(iz>m .OR. iz<1) THEN
          WRITE(*,*) 'z:', i, x(3,i)
          STOP 'SOD: Error, particle outside box'
       END IF

       !See if particle is within spheres of over 27 neighbouring spheres
       DO jx=-1,1
          DO jy=-1,1
             DO jz=-1,1

                !Find the x integer for the sphere
                kx=ix+jx
                IF(kx<1) THEN
                   kx=kx+m
                ELSE IF(kx>m) THEN
                   kx=kx-m
                END IF

                !Find the y integer for the sphere
                ky=iy+jy
                IF(ky<1) THEN
                   ky=ky+m
                ELSE IF(ky>m) THEN
                   ky=ky-m
                END IF

                !Find the z integer for the sphere
                kz=iz+jz
                IF(kz<1) THEN
                   kz=kz+m
                ELSE IF(kz>m) THEN
                   kz=kz-m
                END IF
                !Done

                !WRITE(*,*) 'Mesh:', jx, jy, jz, kx, ky, kz

                !Calculate the displacement between the particle and the sphere centre
                dx=x(1,i)-xc(1,kx,ky,kz)
                dy=x(2,i)-xc(2,kx,ky,kz)
                dz=x(3,i)-xc(3,kx,ky,kz)
                !!!

                !WRITE(*,*) 'Displacements:', dx, dy, dz
                
                !!!
                !Correct if particle is on the other side of the box - x
                IF(dx>L/2.) THEN
                   dx=dx-L
                ELSE IF(dx<-L/2.) THEN
                   dx=dx+L
                END IF

                !Correct if particle is on the other side of the box - y
                IF(dy>L/2.) THEN
                   dy=dy-L
                ELSE IF(dy<-L/2.) THEN
                   dy=dy+L
                END IF

                !Correct if particle is on the other side of the box - z
                IF(dz>L/2.) THEN
                   dz=dz-L
                ELSE IF(dz<-L/2.) THEN
                   dz=dz+L
                END IF
                !!!

                !WRITE(*,*) 'Displacements:', dx, dy, dz
                
                !!!
                !Calculate the radial distance between the sphere centre and the particle
                r=sqrt(dx**2.+dy**2.+dz**2.)

                !Add to density if within sphere radius
                IF(r<Rs) THEN
                   d(kx,ky,kz)=d(kx,ky,kz)+w(i)
                   !WRITE(*,*) 'Added:', kx, ky, kz
                END IF
                !!!
                
             END DO
          END DO
       END DO
       !STOP

    END DO

    WRITE(*,*) 'SOD: Binning complete'
    WRITE(*,*)

  END SUBROUTINE SOD

  FUNCTION periodic_distance(x1,x2,L)

    !Calculates the distance between x1 and x2 assuming
    !that they are coordinates in a periodic box
    IMPLICIT NONE
    REAL :: periodic_distance
    REAL, INTENT(IN) :: x1(3), x2(3), L
    REAL :: dx(3)
    INTEGER :: i

    !Initially dx is just the absolute vector difference
    dx=ABS(x2-x1)

    !Now check if any legs are greater than half-box size
    !Note the Cartesian distance *cannot* be larger than L/2
    DO i=1,3
       IF(dx(i)>L/2.) THEN
          dx(i)=L-dx(i)
       !ELSE IF(dx(j)<-L/2.) THEN
       !   dx(j)=dx(j)+L
       END IF
    END DO

    periodic_distance=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
    
  END FUNCTION periodic_distance

  FUNCTION periodic_mean(x1,x2,L)

    !Calculates the periodic mean of two coordinates in a box
    IMPLICIT NONE
    REAL :: periodic_mean(3)
    REAL, INTENT(IN) :: x1(3), x2(3), L
    REAL :: dx(3)
    INTEGER :: i

    !Initially dx is just the absolute vector difference
    dx=ABS(x2-x1)
    
    DO i=1,3
       periodic_mean(i)=0.5*(x1(i)+x2(i))
       IF(dx(i)>L/2.) THEN
          periodic_mean(i)=periodic_mean(i)+L/2.
       END IF
    END DO
    
  END FUNCTION periodic_mean

  SUBROUTINE find_pairs(x,okay,n,rmin,rmax,L,outfile)!pairs,np)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x(3,n)
    LOGICAL, INTENT(IN) :: okay(n)
    REAL, INTENT(IN) :: rmin, rmax, L   
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=256), INTENT(IN) :: outfile
    !REAL, ALLOCATABLE, INTENT(OUT) :: pairs(:,:)
    !INTEGER, INTENT(OUT) :: np
    INTEGER :: i, j, np
    REAL :: r  
    
    !Only scan haloes that fall within the mass range
    !okay=.FALSE.
    !DO i=1,n
    !   IF(m(i)>=mmin .AND. m(i)<=mmax) THEN
    !      okay(i)=.TRUE.
    !   END IF
    !END DO

    WRITE(*,*) 'FIND_PAIRS: Minimum separation [Mpc/h]:', rmin
    WRITE(*,*) 'FIND_PAIRS: Maximum separation [Mpc/h]:', rmax
    WRITE(*,*) 'FIND_PAIRS: Doing finding'

    !Fix the number counter to be zero
    np=0

    OPEN(7,file=outfile)
    DO i=1,n
       DO j=i+1,n
          IF(okay(i) .AND. okay(j)) THEN
             r=periodic_distance(x(:,i),x(:,j),L)
             IF(r>=rmin .AND. r<=rmax) THEN
                !CALL add_to_array(x(:,i),x(:,j),pairs,np)
                np=np+1
                !WRITE(*,*) i, j, r
                WRITE(7,*) x(1,i), x(2,i), x(3,i), x(1,j), x(2,j), x(3,j)
             END IF
          END IF
       END DO
    END DO
    CLOSE(7)

    WRITE(*,*) 'FIND_PAIRS: Number of pairs:', np
    WRITE(*,*) 'FIND_PAIRS: Done'
    WRITE(*,*) 
    
  END SUBROUTINE find_pairs

  SUBROUTINE cut(okay,m,n,min,max)

    !Flags objects that make the cut as 'okay'
    !Can be applied to any scalar array, not just mass
    IMPLICIT NONE
    REAL, INTENT(IN) :: m(n), min, max
    INTEGER, INTENT(IN) :: n
    LOGICAL, ALLOCATABLE, INTENT(OUT) :: okay(:)
    INTEGER :: i, o

    WRITE(*,*) 'CUT: Imposing property cut'    
    WRITE(*,*) 'CUT: Minimum value:', min
    WRITE(*,*) 'CUT: Maximum value:', max
    WRITE(*,*) 'CUT: Original number of objects', n

    ALLOCATE(okay(n))
    okay=.FALSE.
    
    DO i=1,n
       IF(m(i)>=min .AND. m(i)<=max) okay(i)=.TRUE.
    END DO

    o=COUNT(okay)

    WRITE(*,*) 'CUT: Final number of objects:', o
    WRITE(*,*) 'CUT: Fraction remaining:', REAL(o)/REAL(n)
    WRITE(*,*) 'CUT: Fraction culled:', 1.-REAL(o)/REAL(n)
    WRITE(*,*) 'CUT: Done'
    WRITE(*,*)
    
  END SUBROUTINE cut

  SUBROUTINE zshift(x,v,n,Om_m,Om_v,z,iz)

    !Shift particles to redshift space
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, iz
    REAL, INTENT(IN) :: v(3,n), Om_m, Om_v, z
    REAL, INTENT(INOUT) :: x(3,n)
    INTEGER :: i
    REAL :: H, a

    WRITE(*,*) 'ZSHIFT: Shifting particles into redshift space'

    H=100.*sqrt(Hubble2_simple(z,Om_m,Om_v))

    a=1./(1.+z)

    DO i=1,n
       x(iz,i)=x(iz,i)+v(iz,i)/(a*H)
    END DO

    WRITE(*,*) 'ZSHIFT: Done'
    WRITE(*,*)

  END SUBROUTINE zshift

  FUNCTION Hubble2_simple(z,Om_m,Om_v)

    !This calculates the dimensionless squared hubble parameter squared at redshift z!
    !It ignores contributions from radiation (not accurate at very high z)!
    !It also ignores anything other than vacuum and matter
    IMPLICIT NONE
    REAL :: Hubble2_simple, z
    REAL :: Om_m, Om_v

    Hubble2_simple=(Om_m*(1.+z)**3.)+Om_v+((1.-Om_m-Om_v)*(1.+z)**2.)

  END FUNCTION Hubble2_simple

  SUBROUTINE slice(x,x1,x2,y,y1,y2,z,z1,z2,filename)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x1, x2, y1, y2, z1, z2
    CHARACTER(len=64), INTENT(IN) :: filename
    REAL, INTENT(IN) :: x(:), y(:), z(:)
    INTEGER :: i

    WRITE(*,*) 'Writing slice'
    OPEN(10,file=filename)
    WRITE(*,*) 'Thickness in z/(Mpc/h):', (z2-z1)
    DO i=1,SIZE(x)
       IF(x1<x(i) .AND. x(i)<=x2 .AND. y1<y(i) .AND. y(i)<=y2 .AND. z1<z(i) .AND. z(i)<=z2) THEN
          WRITE(10,*) x(i), y(i), z(i)
       END IF
    END DO
    CLOSE(10)
    WRITE(*,*) 'Slice written: ', filename
    WRITE(*,*)

  END SUBROUTINE slice
  
  FUNCTION shot_noise_simple(L,n)

    !Calculate simulation shot noise
    USE constants
    IMPLICIT NONE
    REAL :: shot_noise_simple
    REAL, INTENT(IN) :: L
    INTEGER*8, INTENT(IN) :: n

    !Calculate number density
    shot_noise_simple=L**3/REAL(n)

  END FUNCTION shot_noise_simple

  FUNCTION shot_noise(L,m,n)

    !Calculate simulation shot noise
    USE constants
    IMPLICIT NONE
    REAL :: shot_noise
    REAL, INTENT(IN) :: m(n), L
    INTEGER, INTENT(IN) :: n
    REAL :: Nbar

    !Calculate the effective mean number of tracers
    Nbar=sum_double(m,n)**2/sum_double(m**2,n)
    
    !Calculate number density
    shot_noise=L**3/Nbar
    
  END FUNCTION shot_noise

!!$  SUBROUTINE write_power(k,D2,nbin,nk,np,L,outfile)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: k(nk), D2(nk), L
!!$    INTEGER, INTENT(IN) :: nbin(nk), nk, np
!!$    CHARACTER(len=256), INTENT(IN) :: outfile
!!$    INTEGER :: i
!!$
!!$    WRITE(*,*) 'WRITE_POWER: Output file: ', TRIM(outfile)
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
