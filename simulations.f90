MODULE simulations

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Zeldovich_displacement(x,n,L,s,m)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(3,n)
    REAL, INTENT(IN) :: s(3,m,m,m), L
    INTEGER, INTENT(IN) :: n,m
    INTEGER :: i, ix(3)

    !Loop over all particles
    DO i=1,n
       ix=NGP_cell(x(:,i),L,m) !Find the integer-coordinates for which cell you are in
       x(:,i)=x(:,i)+s(:,ix(1),ix(2),ix(3)) !Do the displacement
    END DO
    
  END SUBROUTINE Zeldovich_displacement

  FUNCTION NGP_cell(x,L,m)

    !Find the integer coordinates of the cell the particle x is in
    IMPLICIT NONE
    INTEGER :: NGP_cell(3)
    REAL, INTENT(IN) :: x(3), L
    INTEGER, INTENT(IN) :: m
    INTEGER :: i

    DO i=1,3
       NGP_cell(i)=NINT(0.5+m*x(i)/L)
    END DO
    
  END FUNCTION NGP_cell
  
  SUBROUTINE generate_randoms(x,n,L)

    !Generate random x,y,z positions in a cube of size L^3
    USE random_numbers
    IMPLICIT NONE
    REAL, INTENT(OUT) :: x(3,n)
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j

    DO i=1,n
       DO j=1,3
          x(j,i)=random_uniform(0.,L)
       END DO
    END DO

  END SUBROUTINE generate_randoms

  SUBROUTINE generate_grid(x,n,L)

    !Generate a grid of positions in a cube of size L^3
    IMPLICIT NONE
    REAL, INTENT(OUT) :: x(3,n)
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: n
    INTEGER :: ix, iy, iz, i, m

    !Check that the particle number is cubic
    m=NINT(n**(1./3.))
    IF(m**3 .NE. n) STOP 'GENERATE_GRID: Error, you need a cubic number of particles for a grid'

    !Set the particle counting variable to zero
    i=0

    !Loop over all particles
    DO iz=1,m
       DO iy=1,m
          DO ix=1,m
             i=i+1 !Increment the particle counter
             x(1,i)=L*(ix-0.5)/REAL(m) !Assign x
             x(2,i)=L*(iy-0.5)/REAL(m) !Assign y
             x(3,i)=L*(iz-0.5)/REAL(m) !Assign z
          END DO
       END DO
    END DO
    
  END SUBROUTINE generate_grid

  SUBROUTINE generate_poor_glass(x,n,L)

    !Generate a poor man's glass in a cube of size L^3
    USE random_numbers
    IMPLICIT NONE
    REAL, INTENT(OUT) :: x(3,n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: L
    INTEGER :: i, j, m
    REAL :: dx

    !First genrate a grid
    CALL generate_grid(x,n,L)

    !The cube-root of the number of particles
    m=NINT(n**(1./3.))

    !How far can the particles be shifted in x,y,z
    !They need to stay in their initial cube region
    dx=L/REAL(m)
    dx=dx/2.

    !Loop over the particles and do the displacement
    DO i=1,n
       DO j=1,3
          x(j,i)=x(j,i)+random_uniform(-dx,dx)
       END DO
    END DO

  END SUBROUTINE generate_poor_glass

  SUBROUTINE sparse_sample(x,v,n,f)

    USE random_numbers
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: x(:,:), v(:,:)
    REAL, ALLOCATABLE :: x_old(:,:), v_old(:,:)
    REAL, INTENT(IN) :: f
    INTEGER, INTENT(INOUT) :: n
    INTEGER :: keep(n), n_old, j, i
    REAL :: rand

    CALL RNG_set(0)
    
    n_old=n

    keep=0
    n=0

    WRITE(*,*) 'Sparse sampling'

    DO i=1,n_old
       
       IF(random_uniform(0.,1.)<f) THEN
          n=n+1
          keep(i)=1
       END IF

    END DO
    
    ALLOCATE(x_old(3,n_old),v_old(3,n_old))

    x_old=x
    v_old=v

    DEALLOCATE(x,v)
    ALLOCATE(x(3,n),v(3,n))

    j=0

    DO i=1,n_old

       IF(keep(i)==1) THEN
          j=j+1
          x(1,j)=x_old(1,i)
          x(2,j)=x_old(2,i)
          x(3,j)=x_old(3,i)
          v(1,j)=v_old(1,i)
          v(2,j)=v_old(2,i)
          v(3,j)=v_old(3,i)
       END IF

    END DO

    DEALLOCATE(x_old,v_old)

    WRITE(*,*) 'Complete'
    WRITE(*,*) 'Before:', n_old
    WRITE(*,*) 'After:', n
    WRITE(*,*) 'Ratio:', float(n)/float(n_old)
    WRITE(*,*)    

  END SUBROUTINE sparse_sample

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
    INTEGER, INTENT(IN) :: n, m
    INTEGER, INTENT(INOUT) :: ibin
    REAL, INTENT(INOUT) :: d(m,m,m)
    REAL, INTENT(IN) :: x(3,n), L, w(n)

    IF(ibin==-1) THEN
       WRITE(*,*) 'Choose binning strategy'
       WRITE(*,*) '1 - NGP'
       WRITE(*,*) '2 - CIC'
       READ(*,*) ibin
       WRITE(*,*)
    END IF

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
                !!

                !WRITE(*,*) 'Displacements:', dx, dy, dz

                !!
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
                !!

                !WRITE(*,*) 'Displacements:', dx, dy, dz

                !!
                !Calculate the radial distance between the sphere centre and the particle
                r=sqrt(dx**2.+dy**2.+dz**2.)

                !Add to density if within sphere radius
                IF(r<Rs) THEN
                   d(kx,ky,kz)=d(kx,ky,kz)+w(i)
                   !WRITE(*,*) 'Added:', kx, ky, kz
                END IF
                !!

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

  REAL FUNCTION Hubble2_simple(z,Om_m,Om_v)

    !This calculates the dimensionless squared hubble parameter squared at redshift z!
    !It ignores contributions from radiation (not accurate at very high z)!
    !It also ignores anything other than vacuum and matter
    IMPLICIT NONE
    REAL, INTENT(IN) :: z, Om_m, Om_v

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

  REAL FUNCTION shot_noise_simple(L,n)

    !Calculate simulation shot noise
    IMPLICIT NONE
    REAL, INTENT(IN) :: L
    INTEGER*8, INTENT(IN) :: n

    !Calculate number density
    shot_noise_simple=L**3/REAL(n)

  END FUNCTION shot_noise_simple

  REAL FUNCTION shot_noise(L,m,n)

    !Calculate simulation shot noise
    USE array_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: m(n), L
    INTEGER, INTENT(IN) :: n
    REAL :: Nbar

    !Calculate the effective mean number of tracers
    Nbar=sum_double(m,n)**2/sum_double(m**2,n)

    !Calculate number density
    shot_noise=L**3/Nbar

  END FUNCTION shot_noise

  REAL FUNCTION shot_noise_k(k,shot)

    USE constants
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, shot

    shot_noise_k=shot*4.*pi*(k/(2.*pi))**3

  END FUNCTION shot_noise_k

END MODULE simulations
