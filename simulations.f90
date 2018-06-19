MODULE simulations

  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_power_spectrum(x,n,L,m,nk,outfile)

    USE constants
    USE field_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(3,n), L
    INTEGER, INTENT(IN) :: n, m, nk
    CHARACTER(len=*), INTENT(IN) :: outfile
    INTEGER :: i
    DOUBLE COMPLEX :: dk(m,m,m)
    REAL, ALLOCATABLE :: k(:), Pk(:)
    INTEGER, ALLOCATABLE :: nbin(:)
    REAL :: kmin, kmax, shot

    CALL sharp_Fourier_density_contrast(x,n,L,dk,m)

    kmin=twopi/L
    kmax=REAL(m)*pi/L
    CALL compute_power_spectrum(dk,dk,m,L,kmin,kmax,nk,k,Pk,nbin)

    WRITE(*,*) 'WRITE_POWER_SPECTRUM: Outfile: ', TRIM(outfile)

    shot=shot_noise_simple(L,INT8(n))
    OPEN(7,file=outfile)
    DO i=1,nk
       IF(nbin(i)==0) CYCLE
       WRITE(7,*) k(i), Pk(i), shot_noise_k(k(i),shot), nbin(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'WRITE_POWER_SPECTRUM: Done'
    WRITE(*,*)

  END SUBROUTINE write_power_spectrum

  SUBROUTINE sharp_Fourier_density_contrast(x,n,L,dk,m)

    USE fft
    USE field_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(3,n), L
    INTEGER, INTENT(IN) :: n, m
    DOUBLE COMPLEX, INTENT(OUT) :: dk(m,m,m)
    REAL :: w(n), dbar
    REAL :: d(m,m,m)
    DOUBLE COMPLEX :: dk_out(m,m,m)

    !Things that would like to be PARAMETERS
    INTEGER :: ibin=2 !Set the binning strategy

    !Assign weight=1 and bin the particles
    w=1.
    CALL particle_bin(x,n,L,w,d,m,ibin)
    dbar=REAL(n)/REAL(m)**3
    d=d/dbar

    dk=d
    CALL fft3(dk,dk_out,m,m,m,-1)
    dk=dk_out

    CALL sharpen_k(dk,m,L,ibin)    

  END SUBROUTINE sharp_Fourier_density_contrast

  SUBROUTINE write_density_slice_ascii(x,n,z1,z2,L,m,outfile)

    !Write out a slice of density field
    !x(n), y(n), z(n): particle positions
    !x1->x2, y1->y2, z1->z2: range for the slice
    !L: box size [Mpc/h]
    !s: Smoothing length [Mpc/h]
    !m: mesh size for the density field
    !outfile: output file
    USE field_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(3,n), L
    REAL, INTENT(IN) :: z1, z2    
    INTEGER, INTENT(IN) :: n, m
    REAL :: d(m,m)
    CHARACTER(len=*), INTENT(IN) :: outfile
    REAL :: xp, yp, smoothing
    INTEGER :: i, j    

    !Make the projected 2D density  
    CALL make_projected_density(x,n,z1,z2,L,d,m)

    !Smooth the density field
    smoothing=L/REAL(m) !Set the smoothing scale to be the mesh size
    CALL smooth2D(d,m,smoothing,L)

    !Write out to file
    WRITE(*,*) 'WRITE_DENSITY_SLICE_ASCII: Writing density map'
    WRITE(*,*) 'WRITE_DENSITY_SLICE_ASCII: File: ', TRIM(outfile)
    OPEN(9,file=outfile)
    DO i=1,m
       DO j=1,m
          !xp=(REAL(i)-0.5)/REAL(m)
          !yp=(REAL(j)-0.5)/REAL(m)
          xp=cell_position(i,L,m)
          yp=cell_position(j,L,m)
          WRITE(9,*) xp, yp, d(i,j)
       END DO
    END DO
    CLOSE(9)

    WRITE(*,*) 'WRITE_DENSITY_SLICE_ASCII: Done'
    WRITE(*,*)  

  END SUBROUTINE write_density_slice_ascii

  SUBROUTINE make_projected_density(x,n,z1,z2,L,d,m)

    !Write out a slice of density field
    !x(n), y(n), z(n): particle positions
    !x1->x2, y1->y2, z1->z2: range for the slice
    !L: box size [Mpc/h]
    !m: mesh size for the density field
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(3,n), L
    REAL, INTENT(IN) :: z1, z2
    REAL, INTENT(OUT) :: d(m,m)
    INTEGER, INTENT(IN) :: n, m
    REAL :: vfac, dbar
    REAL, ALLOCATABLE :: x2D(:,:)
    INTEGER :: i, i2D, n2D
    INTEGER :: ibin

    !Count the number of particles falling into the slice
    n2D=0
    DO i=1,n
       IF(x(3,i)>=z1 .AND. x(3,i)<=z2) n2D=n2D+1
    END DO
    ALLOCATE(x2D(2,n2D))

    !Make the 2D position array
    i2D=0
    DO i=1,n
       IF(x(3,i)>=z1 .AND. x(3,i)<=z2) THEN
          i2D=i2D+1
          x2D(1,i2D)=x(1,i)
          x2D(2,i2D)=x(2,i)
       END IF
    END DO

    !Bin for the 2D density field and convert to relative density
    ibin=2
    CALL particle_bin_2D(x2D,n2D,L,d,m,ibin)
    DEALLOCATE(x2D)
    vfac=(z2-z1)/L
    dbar=(REAL(n)*vfac)/REAL(m**2)
    d=d/dbar

    WRITE(*,*) 'MAKE_PROJECTED_DENSITY: Thickness in z [Mpc/h]:', z2-z1
    WRITE(*,*) 'MAKE_PROJECTED_DENSITY: Volume factor:', vfac
    WRITE(*,*) 'MAKE_PROJECTED_DENSITY: Mean cell particle density:', dbar
    WRITE(*,*) 'MAKE_PROJECTED_DENSITY: Done'
    WRITE(*,*)

  END SUBROUTINE make_projected_density

  SUBROUTINE Zeldovich_ICs(x,v,n,L,logk_tab,logPk_tab,nk,vfac,m,use_average)

    USE field_operations
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(3,n), v(3,n)
    REAL, INTENT(IN) :: logk_tab(nk), logPk_tab(nk), L, vfac
    INTEGER, INTENT(IN) :: n, m, nk
    LOGICAL, INTENT(IN) :: use_average
    REAL :: f(3,m,m,m), ips, maxf

    !Make the displacement field
    CALL generate_displacement_fields(f,m,L,logk_tab,logPk_tab,nk,use_average)

    !Calculate some useful things
    ips=L/REAL(n**(1./3.)) !Mean ID inter-particle spacing
    maxf=MAXVAL(f) !Maximum value of 1D displacement    

    !Calculate the particle velocities first
    CALL Zeldovich_velocity(x,v,n,L,vfac*f,m)

    !Then do the particle positions
    CALL Zeldovich_displacement(x,n,L,f,m)

    !Write some useful things to the screen
    WRITE(*,*) 'ZELDOVICH_ICS: Max 1D displacement [Mpc/h]:', maxf
    WRITE(*,*) 'ZELDOVICH_ICS: Inter-particle spacing [Mpc/h]:', ips
    WRITE(*,*) 'ZELDOVICH_ICS: Max 1D displacement in units of the IPS:', maxf/ips
    WRITE(*,*) 'ZELDOVICH_ICS: Done'
    WRITE(*,*)

  END SUBROUTINE Zeldovich_ICs

  SUBROUTINE Zeldovich_displacement(x,n,L,s,m)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(3,n)
    REAL, INTENT(IN) :: s(3,m,m,m), L
    INTEGER, INTENT(IN) :: n, m
    INTEGER :: i, j, ix(3)

    WRITE(*,*) 'ZELDOVICH_DISPLACEMENT: Displacing particles'

    !Loop over all particles
    DO i=1,n
       DO j=1,3
          ix(j)=NGP_cell(x(j,i),L,m) !Find the integer-coordinates for which cell you are in
       END DO
       x(:,i)=x(:,i)+s(:,ix(1),ix(2),ix(3)) !Do the displacement
    END DO

    CALL replace(x,n,L)

    WRITE(*,*) 'ZELDOVICH_DISPLACEMENT: Done'
    WRITE(*,*)

  END SUBROUTINE Zeldovich_displacement

  SUBROUTINE Zeldovich_velocity(x,v,n,L,s,m)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: v(3,n)
    REAL, INTENT(IN) :: x(3,n), s(3,m,m,m), L
    INTEGER, INTENT(IN) :: n, m
    INTEGER :: i, j, ix(3)

    WRITE(*,*) 'ZELDOVICH_VELOCITY: Assigining particle velocties'
    WRITE(*,*) 'ZELDOVICH_VELOCITY: Any previous velocity set to zero'

    !Loop over all particles
    DO i=1,n
       DO j=1,3
          ix(j)=NGP_cell(x(j,i),L,m) !Find the integer-coordinates for which cell you are in
       END DO
       v(:,i)=s(:,ix(1),ix(2),ix(3)) !Assign the velocity
    END DO

    WRITE(*,*) 'ZELDOVICH_VELOCITY: Done'
    WRITE(*,*)

  END SUBROUTINE Zeldovich_velocity

  INTEGER FUNCTION NGP_cell(x,L,m)

    !Find the integer coordinates of the cell the particle x is in
    IMPLICIT NONE
    REAL, INTENT(IN) :: x, L
    INTEGER, INTENT(IN) :: m

    NGP_cell=NINT(0.5+m*x/L)

    IF(NGP_cell<1 .OR. NGP_cell>m) THEN
       WRITE(*,*) 'NGP_CELL: Particle position [Mpc/h]:', x
       WRITE(*,*) 'NGP_CELL: Box size [Mpc/h]:', L
       WRITE(*,*) 'NGP_CELL: Mesh size:', m 
       WRITE(*,*) 'NGP_CELL: Assigned cell:', NGP_cell
       STOP 'NGP_CELL: Error, the assigned cell position is outside the mesh'
    END IF

  END FUNCTION NGP_cell

  REAL FUNCTION cell_position(i,L,m)

    IMPLICIT NONE
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: i, m

    cell_position=L*(i-0.5)/REAL(m)

  END FUNCTION cell_position

  SUBROUTINE generate_randoms(x,n,L)

    !Generate random x,y,z positions in a cube of size L^3
    USE random_numbers
    IMPLICIT NONE
    REAL, INTENT(OUT) :: x(3,n)
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j
    REAL :: dx

    !To prevent the particle being exactly at zero
    dx=L/1e8

    !Set the random-number generator
    !CALL RNG_set(0)

    WRITE(*,*) 'GENERATE_RANDOMS: Generating a uniform-random particle distribution'

    !Loop over all particles and coordinates and assign randomly
    DO i=1,n
       DO j=1,3
          x(j,i)=random_uniform(dx,L)
       END DO
    END DO

    WRITE(*,*) 'GENERATE_RANDOMS: Done'
    WRITE(*,*)

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

    WRITE(*,*) 'GENERATE_GRID: Generating a grid particle distribution'

    !Loop over all particles
    i=0 !Set the particle counting variable to zero
    DO iz=1,m
       DO iy=1,m
          DO ix=1,m
             i=i+1 !Increment the particle counter
             !x(1,i)=L*(ix-0.5)/REAL(m) !Assign x
             !x(2,i)=L*(iy-0.5)/REAL(m) !Assign y
             !x(3,i)=L*(iz-0.5)/REAL(m) !Assign z
             x(1,i)=cell_position(ix,L,m)
             x(2,i)=cell_position(iy,L,m)
             x(3,i)=cell_position(iz,L,m)
          END DO
       END DO
    END DO

    WRITE(*,*) 'GENERATE_GRID: Done'
    WRITE(*,*)

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

    WRITE(*,*) 'GENERATE_POOR_GLASS: Generating a poor man glass from the grid'

    !Loop over the particles and do the displacement
    DO i=1,n
       DO j=1,3
          x(j,i)=x(j,i)+random_uniform(-dx,dx)
       END DO
    END DO

    WRITE(*,*) 'GENERATE_POOR_GLASS: Done'
    WRITE(*,*)

  END SUBROUTINE generate_poor_glass

  SUBROUTINE sparse_sample(x,v,n,f)

    !USE random_numbers
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: x(:,:), v(:,:)
    REAL :: x_old(3,n), v_old(3,n)
    REAL, INTENT(IN) :: f
    INTEGER, INTENT(INOUT) :: n
    INTEGER :: keep(n), n_old, j, i
    REAL :: rand

    !CALL RNG_set(0)

    n_old=n

    keep=0
    n=0

    WRITE(*,*) 'Sparse sampling'

    DO i=1,n_old

       IF(rand(0)<f) THEN
          n=n+1
          keep(i)=1
       END IF

    END DO

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

    WRITE(*,*) 'Complete'
    WRITE(*,*) 'Before:', n_old
    WRITE(*,*) 'After:', n
    WRITE(*,*) 'Ratio:', REAL(n)/REAL(n_old)
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

  SUBROUTINE particle_bin_2D(x,n,L,d,m,ibin)

    !Bin particle properties onto a mesh, summing as you go
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, m
    INTEGER, INTENT(INOUT) :: ibin
    REAL, INTENT(OUT) :: d(m,m)
    REAL, INTENT(IN) :: x(2,n), L

    IF(ibin==-1) THEN
       WRITE(*,*) 'Choose binning strategy'
       WRITE(*,*) '1 - NGP'
       WRITE(*,*) '2 - CIC'
       READ(*,*) ibin
       WRITE(*,*)
    END IF

    IF(ibin==1) THEN
       CALL NGP2D(x,n,L,d,m)
    ELSE IF(ibin==2) THEN
       CALL CIC2D(x,n,L,d,m)
    ELSE
       STOP 'PARTICLE_BIN: Error, ibin not specified correctly'
    END IF

  END SUBROUTINE particle_bin_2D

  SUBROUTINE NGP2D(x,n,L,d,m)

    !Nearest-grid-point binning routine
    USE statistics
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(2,n), L
    INTEGER, INTENT(IN) :: n, m
    REAL, INTENT(OUT) :: d(m,m)
    INTEGER :: i, ix, iy

    WRITE(*,*) 'NGP2D: Binning particles and creating field'
    WRITE(*,*) 'NGP2D: Cells:', m

    !Set array to zero explicitly
    d=0.

    DO i=1,n

       !ix=CEILING(x(1,i)*REAL(m)/L)
       !iy=CEILING(x(2,i)*REAL(m)/L)
       ix=NGP_cell(x(1,i),L,m)
       iy=NGP_cell(x(2,i),L,m)

!!$       IF(ix>m .OR. ix<1) THEN
!!$          WRITE(*,*) 'x:', i, x(1,i)
!!$          STOP 'NGP2D: Warning, point outside box'
!!$       END IF
!!$
!!$       IF(iy>m .OR. iy<1) THEN
!!$          WRITE(*,*) 'y:', i, x(2,i)
!!$          STOP 'NGP: Warning, point outside box'
!!$       END IF

       d(ix,iy)=d(ix,iy)+1.

    END DO

!!$    WRITE(*,*) 'NGP2D: Average:', mean(d,m)
!!$    WRITE(*,*) 'NGP2D: RMS:', sqrt(variance(d,m))
!!$    WRITE(*,*) 'NGP2D: Minimum:', MINVAL(REAL(d))
!!$    WRITE(*,*) 'NGP2D: Maximum:', MAXVAL(REAL(d))
    WRITE(*,*) 'NGP2D: Binning complete'
    WRITE(*,*)

  END SUBROUTINE NGP2D

  SUBROUTINE CIC2D(x,n,L,d,m)

    !This could probably be usefully combined with CIC somehow
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(2,n), L
    INTEGER, INTENT(IN) :: n, m
    REAL, INTENT(OUT) :: d(m,m)
    INTEGER :: i, ix, iy, ixn, iyn
    REAL :: dx, dy

    WRITE(*,*) 'CIC2D: Binning particles and creating density field'
    WRITE(*,*) 'CIC2D: Cells:', m

    !Set array to zero explicitly
    d=0.

    DO i=1,n

       ! Get the cell interger coordinates
       ix=NGP_cell(x(1,i),L,m)
       iy=NGP_cell(x(2,i),L,m)

       !dx, dy in cell units, away from cell centre
       dx=(x(1,i)/L)*REAL(m)-(REAL(ix)-0.5)
       dy=(x(2,i)/L)*REAL(m)-(REAL(iy)-0.5)

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

    WRITE(*,*) 'CIC2D: Binning complete'
    WRITE(*,*)

  END SUBROUTINE CIC2D

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

    !Set array to zero explicitly
    d=0.

    DO i=1,n

       ix=CEILING(x(1,i)*REAL(m)/L)
       iy=CEILING(x(2,i)*REAL(m)/L)
       iz=CEILING(x(3,i)*REAL(m)/L)

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

    !Set array to zero explicitly
    d=0.

    DO i=1,n

       ix=CEILING(x(1,i)*REAL(m)/L)
       iy=CEILING(x(2,i)*REAL(m)/L)
       iz=CEILING(x(3,i)*REAL(m)/L)

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
       dx=(x(1,i)/L)*REAL(m)-(REAL(ix)-0.5)
       dy=(x(2,i)/L)*REAL(m)-(REAL(iy)-0.5)
       dz=(x(3,i)/L)*REAL(m)-(REAL(iz)-0.5)

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

  SUBROUTINE SOD(x,n,L,w,d,m,Rs)

    !Spherical density routine
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(OUT) :: d(m,m,m)
    REAL :: xc(3,m,m,m)
    REAL, INTENT(IN) :: x(3,n), L, w(n), Rs
    REAL :: r, dx, dy, dz
    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz, kx, ky, kz
    INTEGER, INTENT(IN) :: m

    WRITE(*,*) 'SOD: Binning particles and creating density field'
    WRITE(*,*) 'SOD: Cells:', m

    d=0.

    !Fill array with sphere centre positions
    DO k=1,m
       DO j=1,m
          DO i=1,m
             !xc(1,i,j,k)=L*(REAL(i)-0.5)/REAL(m)
             !xc(2,i,j,k)=L*(REAL(j)-0.5)/REAL(m)
             !xc(3,i,j,k)=L*(REAL(k)-0.5)/REAL(m)
             xc(1,i,j,k)=cell_position(i,L,m)
             xc(2,i,j,k)=cell_position(j,L,m)
             xc(3,i,j,k)=cell_position(k,L,m)
          END DO
       END DO
    END DO

    !Loop over all particles and assign to spheres
    DO i=1,n

       !WRITE(*,*) 'Particle coordinates:', i, x(1,i), x(2,i), x(3,i)

       !Coordinate of nearest mesh cell
       ix=CEILING(x(1,i)*REAL(m)/L)
       iy=CEILING(x(2,i)*REAL(m)/L)
       iz=CEILING(x(3,i)*REAL(m)/L)

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

  SUBROUTINE find_pairs(x,okay,n,rmin,rmax,L,outfile)!pairs,np)

    USE field_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(3,n)
    LOGICAL, INTENT(IN) :: okay(n)
    REAL, INTENT(IN) :: rmin, rmax, L   
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=*), INTENT(IN) :: outfile
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
    LOGICAL, INTENT(OUT) :: okay(n)
    INTEGER :: i, o

    WRITE(*,*) 'CUT: Imposing property cut'    
    WRITE(*,*) 'CUT: Minimum value:', min
    WRITE(*,*) 'CUT: Maximum value:', max
    WRITE(*,*) 'CUT: Original number of objects', n

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

    Hubble2_simple=Om_m*(1.+z)**3+Om_v+(1.-Om_m-Om_v)*(1.+z)**2

  END FUNCTION Hubble2_simple

  !Used to be called slice
  SUBROUTINE write_slice_ascii(x,x1,x2,y,y1,y2,z,z1,z2,filename)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x1, x2, y1, y2, z1, z2
    CHARACTER(len=*), INTENT(IN) :: filename
    REAL, INTENT(IN) :: x(:), y(:), z(:)
    INTEGER :: i

    WRITE(*,*) 'WRITE_SLICE_ASCII: Writing slice'
    OPEN(10,file=filename)
    WRITE(*,*) 'WRITE_SLICE_ASCII: Thickness in z [Mpc/h]:', (z2-z1)
    DO i=1,SIZE(x)
       IF(x1<x(i) .AND. x(i)<=x2 .AND. y1<y(i) .AND. y(i)<=y2 .AND. z1<z(i) .AND. z(i)<=z2) THEN
          WRITE(10,*) x(i), y(i), z(i)
       END IF
    END DO
    CLOSE(10)
    WRITE(*,*) 'WRITE_SLICE_ASCII: Slice written: ', TRIM(filename)
    WRITE(*,*)

  END SUBROUTINE write_slice_ascii

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

    shot_noise_k=shot*4.*pi*(k/twopi)**3

  END FUNCTION shot_noise_k

  SUBROUTINE adaptive_density(xc,yc,Lsub,z1,z2,m,r,x,n,L,outfile)    

    USE string_operations
    USE field_operations
    IMPLICIT NONE
    REAL, INTENT(IN) :: xc, yc, Lsub, z1, z2
    INTEGER, INTENT(IN) :: m, n, r
    REAL, INTENT(IN) :: x(3,n), L
    CHARACTER(len=*), INTENT(IN) :: outfile

    REAL :: x1, x2, y1, y2, Lx, Ly, Lz, delta, nbar, npexp
    INTEGER :: np
    INTEGER :: m1, m2, m3, m4, m5
    INTEGER :: i, j
    REAL, ALLOCATABLE :: y(:,:), d(:,:)
    REAL, ALLOCATABLE :: c1(:,:), c2(:,:), c3(:,:), c4(:,:), c5(:,:)
    REAL, ALLOCATABLE :: d1(:,:), d2(:,:), d3(:,:), d4(:,:), d5(:,:)
    REAL, ALLOCATABLE :: di(:,:,:), ci(:,:,:)
    CHARACTER(len=256) :: base, ext, output

    !INTEGER, PARAMETER :: r=4 ! Number of refinements
    REAL, PARAMETER :: dc=4. ! Refinement conditions (particles-per-cell)
    REAL, PARAMETER :: fcell=1. ! Smoothing factor over cell sizes (maybe should be set to 1; 0.75 looks okay)
    LOGICAL, PARAMETER :: test=.FALSE. ! Activate test mode

    IF(r>4) STOP 'ADAPTIVE_DENSITY: Error, too many refinement leves requested. Maximum is 4'

    ! Write information to screen
    WRITE(*,*) 'ADAPTIVE_DENSITY: Adaptive density field generator'
    WRITE(*,*) 'ADAPTIVE_DENISTY: Refinement conditions (particles-per-cell):', dc
    WRITE(*,*) 'ADAPTIVE_DENSITY: Smoothing factor in cell size:', fcell
    WRITE(*,*) 'ADAPTIVE_DENSITY: Number of refinements:', r
    WRITE(*,*)

    ! Calculate the 2D average particle number density
    nbar=REAL(n)/L**2

    ! Set the region boundaries in Mpc/h
    x1=xc-Lsub/2.
    x2=xc+Lsub/2.
    y1=yc-Lsub/2.
    y2=yc+Lsub/2.

    ! Set the region thicknesses
    Lx=x2-x1
    Ly=y2-y1
    Lz=z2-z1

    ! Write information to screen
    WRITE(*,*) 'ADAPTIVE_DENSITY: Subvolume:'
    WRITE(*,*) 'ADAPTIVE_DENSITY: x1:', x1
    WRITE(*,*) 'ADAPTIVE_DENSITY: x2:', x2
    WRITE(*,*) 'ADAPTIVE_DENSITY: Lx:', Lx
    WRITE(*,*) 'ADAPTIVE_DENSITY: y1:', y1
    WRITE(*,*) 'ADAPTIVE_DENSITY: y2:', y2
    WRITE(*,*) 'ADAPTIVE_DENSITY: Ly:', Ly
    WRITE(*,*) 'ADAPTIVE_DENSITY: z1:', z1
    WRITE(*,*) 'ADAPTIVE_DENSITY: z2:', z2
    WRITE(*,*) 'ADAPTIVE_DENSITY: Lz:', Lz
    WRITE(*,*)

    ! First pass to count the number of particles in the region
    np=0
    DO i=1,n
       IF(x(1,i)>x1 .AND. x(1,i)<x2 .AND. x(2,i)>y1 .AND. x(2,i)<y2 .AND. x(3,i)>z1 .AND. x(3,i)<z2) THEN
          np=np+1
       END IF
    END DO

    ! Caclulate density statistics
    npexp=REAL(n)*(Lx*Ly*Lz/L**3.)
    delta=-1.+REAL(np)/REAL(npexp)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Density statistics'
    WRITE(*,*) 'ADAPTIVE_DENSITY: Particles in region:', np
    WRITE(*,*) 'ADAPTIVE_DENSITY: Expectation of partilces in region:', npexp
    WRITE(*,*) 'ADAPTIVE_DENSITY: Region over-density:', delta
    WRITE(*,*)

    ALLOCATE(y(2,np))

    ! Second pass to add particles in the region to 2D array y
    j=0
    DO i=1,n
       IF(x(1,i)>x1 .AND. x(1,i)<x2 .AND. x(2,i)>y1 .AND. x(2,i)<y2 .AND. x(3,i)>z1 .AND. x(3,i)<z2) THEN
          j=j+1
          y(1,j)=x(1,i)-xc+Lsub/2.
          y(2,j)=x(2,i)-yc+Lsub/2.
       END IF
    END DO

    ! Set the sizes of all of the adaptive meshes, which all differ by a factor of 2
    ! m1=m, m2=m/2, m3=m/4, m4=m/8, m5=m/16
    m1=m
    m2=m1/2
    m3=m2/2
    m4=m3/2
    m5=m4/2

    ! Write out sizes to screen
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 1:', m1
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 2:', m2
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 3:', m3
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 4:', m4
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 5:', m5
    WRITE(*,*)

    ! Allocate arrays for particle counts and overdensity
    ! This is very ugly, but I do not think it can be avoided with d(5,m,m) because m are all different
    ! Could remove count array if necessary, but including it makes the calculation more obvious
    ALLOCATE(d1(m1,m1),c1(m1,m1))
    ALLOCATE(d2(m2,m2),c2(m2,m2))
    ALLOCATE(d3(m3,m3),c3(m3,m3))
    ALLOCATE(d4(m4,m4),c4(m4,m4))
    ALLOCATE(d5(m5,m5),c5(m5,m5))    

    ! Do CIC binning on each mesh resolution
    ! These CIC routines assume the volume is periodic, so you get some weird edge effects
    CALL CIC2D(y,np,Lsub,c1,m1)
    CALL CIC2D(y,np,Lsub,c2,m2)
    CALL CIC2D(y,np,Lsub,c3,m3)
    CALL CIC2D(y,np,Lsub,c4,m4)
    CALL CIC2D(y,np,Lsub,c5,m5)

    ! Convert to overdensities (1+delta = rho/rhobar)
    ! Could certainly save memory here by having d1=c1, but having the c arrays makes the calculation more obvious
    d1=(c1/((Lsub/REAL(m1))**2))/nbar
    d2=(c2/((Lsub/REAL(m2))**2))/nbar
    d3=(c3/((Lsub/REAL(m3))**2))/nbar
    d4=(c4/((Lsub/REAL(m4))**2))/nbar
    d5=(c5/((Lsub/REAL(m5))**2))/nbar

    ! Smooth density fields
    CALL smooth2D(d1,m1,fcell*Lsub/REAL(m1),Lsub)
    CALL smooth2D(d2,m2,fcell*Lsub/REAL(m2),Lsub)
    CALL smooth2D(d3,m3,fcell*Lsub/REAL(m3),Lsub)
    CALL smooth2D(d4,m4,fcell*Lsub/REAL(m4),Lsub)
    CALL smooth2D(d5,m5,fcell*Lsub/REAL(m5),Lsub)

    ! Write out density statistics on each mesh
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m1
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max count:', MAXVAL(c1)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min count:', MINVAL(c1)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max overdensity:', MAXVAL(d1)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min overdensity:', MINVAL(d1)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m2
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max count:', MAXVAL(c2)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min count:', MINVAL(c2)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max overdensity:', MAXVAL(d2)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min overdensity:', MINVAL(d2)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m3
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max count:', MAXVAL(c3)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min count:', MINVAL(c3)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max overdensity:', MAXVAL(d3)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min overdensity:', MINVAL(d3)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m4
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max count:', MAXVAL(c4)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min count:', MINVAL(c4)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max overdensity:', MAXVAL(d4)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min overdensity:', MINVAL(d4)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m5
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max count:', MAXVAL(c5)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min count:', MINVAL(c5)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max overdensity:', MAXVAL(d5)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min overdensity:', MINVAL(d5)
    WRITE(*,*)

    ! Write out each mesh if testing
    IF(test) THEN       
       
       base='density_'
       ext='.dat'
       
       output=number_file_zeroes(base,m1,4,ext)
       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
       CALL write_2D_field_ascii(d1,m1,Lsub,output)
       
       output=number_file_zeroes(base,m2,4,ext)
       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
       CALL write_2D_field_ascii(d2,m2,Lsub,output)
       
       output=number_file_zeroes(base,m3,4,ext)       
       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
       CALL write_2D_field_ascii(d3,m3,Lsub,output)
       
       output=number_file_zeroes(base,m4,4,ext)
       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
       CALL write_2D_field_ascii(d4,m4,Lsub,output)
       
       output=number_file_zeroes(base,m5,4,ext)
       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
       CALL write_2D_field_ascii(d5,m5,Lsub,output)
       
    END IF

    ! No refinement required
    IF(r==0) THEN
       
       ALLOCATE(d(m,m))
       d=d1

    ! Start the refinements
    ELSE
       
       ! First refinement
       IF(r>=4) THEN     

          ! Allocate d with the size of the second-corsest image
          ALLOCATE(d(m4,m4))
          DO i=1,m5
             DO j=1,m5
                IF(c5(i,j)<dc) THEN
                   ! Use coarser density
                   d(2*i-1,2*j-1)=d5(i,j)
                   d(2*i,2*j-1)=d5(i,j)
                   d(2*i-1,2*j)=d5(i,j)
                   d(2*i,2*j)=d5(i,j)
                ELSE
                   ! Use finer density
                   d(2*i-1,2*j-1)=d4(2*i-1,2*j)
                   d(2*i,2*j-1)=d4(2*i,2*j-1)
                   d(2*i-1,2*j)=d4(2*i-1,2*j)
                   d(2*i,2*j)=d4(2*i,2*j)
                END IF
             END DO
          END DO

       END IF

       ! Second refinement
       IF(r>=3) THEN

          IF(ALLOCATED(d)) THEN
             d4=d
             DEALLOCATE(d)
          END IF
          ALLOCATE(d(m3,m3))
          DO i=1,m4
             DO j=1,m4
                IF(c4(i,j)<dc) THEN
                   ! Use coarser density
                   d(2*i-1,2*j-1)=d4(i,j)
                   d(2*i,2*j-1)=d4(i,j)
                   d(2*i-1,2*j)=d4(i,j)
                   d(2*i,2*j)=d4(i,j)
                ELSE
                   ! Use finer density
                   d(2*i-1,2*j-1)=d3(2*i-1,2*j)
                   d(2*i,2*j-1)=d3(2*i,2*j-1)
                   d(2*i-1,2*j)=d3(2*i-1,2*j)
                   d(2*i,2*j)=d3(2*i,2*j)
                END IF
             END DO
          END DO

       END IF

       ! Third refinement
       IF(r>=2) THEN

          IF(ALLOCATED(d)) THEN
             d3=d
             DEALLOCATE(d)
          END IF
          ALLOCATE(d(m2,m2))
          DO i=1,m3
             DO j=1,m3
                IF(c3(i,j)<dc) THEN
                   ! Use coarser density
                   d(2*i-1,2*j-1)=d3(i,j)
                   d(2*i,2*j-1)=d3(i,j)
                   d(2*i-1,2*j)=d3(i,j)
                   d(2*i,2*j)=d3(i,j)
                ELSE
                   ! Use finer density
                   d(2*i-1,2*j-1)=d2(2*i-1,2*j)
                   d(2*i,2*j-1)=d2(2*i,2*j-1)
                   d(2*i-1,2*j)=d2(2*i-1,2*j)
                   d(2*i,2*j)=d2(2*i,2*j)
                END IF
             END DO
          END DO

       END IF

       ! Fourth refinement
       IF(r>=1) THEN

          IF(ALLOCATED(d)) THEN
             d2=d
             DEALLOCATE(d)
          END IF
          ALLOCATE(d(m1,m1))
          DO i=1,m2
             DO j=1,m2
                IF(c2(i,j)<dc) THEN
                   ! Use coarser density
                   d(2*i-1,2*j-1)=d2(i,j)
                   d(2*i,2*j-1)=d2(i,j)
                   d(2*i-1,2*j)=d2(i,j)
                   d(2*i,2*j)=d2(i,j)
                ELSE
                   ! Use finer density
                   d(2*i-1,2*j-1)=d1(2*i-1,2*j)
                   d(2*i,2*j-1)=d1(2*i,2*j-1)
                   d(2*i-1,2*j)=d1(2*i-1,2*j)
                   d(2*i,2*j)=d1(2*i,2*j)
                END IF
             END DO
          END DO

       END IF

    END IF

    ! Write info to screen
    WRITE(*,*) 'ADAPTIVE_DENSITY: Adaptive density field'
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max density:', MAXVAL(d)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min density:', MINVAL(d)
    WRITE(*,*)

    ! Ensure all cells are positive (note this is to make pretty pictures, not for science)
    DO i=1,m1
       DO j=1,m1
          IF(d(i,j)<0.) d(i,j)=0.
       END DO
    END DO

    ! Write field statistics
    WRITE(*,*) 'ADAPTIVE_DENSITY: Smoothed adaptive density field'
    WRITE(*,*) 'ADAPTIVE_DENSITY: Max density:', MAXVAL(d)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Min density:', MINVAL(d)
    WRITE(*,*)

    ! Finally write out an ascii file for plotting
    CALL write_2D_field_ascii(d,m1,Lsub,outfile)
    WRITE(*,*) 'ADAPTIVE_DENSITY: Done'
    WRITE(*,*)

  END SUBROUTINE adaptive_density

!!$  SUBROUTINE adaptive_density_original(xc,yc,Lsub,z1,z2,m,dcc,fac,x,n,L,outfile)    
!!$
!!$    USE string_operations
!!$    USE field_operations
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: xc, yc, Lsub, z1, z2, dcc, fac
!!$    INTEGER, INTENT(IN) :: m, n
!!$    REAL, INTENT(IN) :: x(3,n), L
!!$    CHARACTER(len=*), INTENT(IN) :: outfile
!!$
!!$    REAL :: x1, x2, y1, y2, Lx, Ly, Lz, delta, dc, nbar, npexp
!!$    INTEGER :: np
!!$    INTEGER :: m1, m2, m3, m4, m5
!!$    INTEGER :: i, j
!!$    REAL, ALLOCATABLE :: y(:,:), d(:,:)
!!$    REAL, ALLOCATABLE :: d1(:,:), d2(:,:), d3(:,:), d4(:,:), d5(:,:)
!!$    CHARACTER(len=256) :: base, ext, output
!!$
!!$    REAL :: fcell=1. ! Smoothing factor over cell sizes (maybe should be set to 1; 0.75 looks okay)
!!$    LOGICAL, PARAMETER :: test=.TRUE. ! Activate test mode
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Adaptive density field generator'
!!$
!!$    ! Calculate the 2D average particle number density
!!$    nbar=REAL(n)/L**2
!!$
!!$    ! Set the region boundaries in Mpc/h
!!$    x1=xc-Lsub/2.
!!$    x2=xc+Lsub/2.
!!$    y1=yc-Lsub/2.
!!$    y2=yc+Lsub/2.
!!$
!!$    ! Set the region thicknesses
!!$    Lx=x2-x1
!!$    Ly=y2-y1
!!$    Lz=z2-z1
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Subvolume:'
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: x1:', x1
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: x2:', x2
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Lx:', Lx
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: y1:', y1
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: y2:', y2
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Ly:', Ly
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: z1:', z1
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: z2:', z2
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Lz:', Lz
!!$    WRITE(*,*)
!!$
!!$    ! First pass to count the number of particles in the region
!!$    np=0
!!$    DO i=1,n
!!$       IF(x(1,i)>x1 .AND. x(1,i)<x2 .AND. x(2,i)>y1 .AND. x(2,i)<y2 .AND. x(3,i)>z1 .AND. x(3,i)<z2) THEN
!!$          np=np+1
!!$       END IF
!!$    END DO
!!$
!!$    npexp=REAL(n)*(Lx*Ly*Lz/L**3.)
!!$    delta=-1.+REAL(np)/REAL(npexp)
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Density statistics'
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Particles in region:', np
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Expectation of partilces in region:', npexp
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Region over-density:', delta
!!$    WRITE(*,*)
!!$
!!$    ALLOCATE(y(2,np))
!!$
!!$    ! Second pass to add particles in the region to 2D array y
!!$    j=0
!!$    DO i=1,n
!!$       IF(x(1,i)>x1 .AND. x(1,i)<x2 .AND. x(2,i)>y1 .AND. x(2,i)<y2 .AND. x(3,i)>z1 .AND. x(3,i)<z2) THEN
!!$          j=j+1
!!$          y(1,j)=x(1,i)
!!$          y(2,j)=x(2,i)
!!$       END IF
!!$    END DO
!!$
!!$    ! Set the sizes of all of the adaptive meshes
!!$    m1=m
!!$    m2=m1/2
!!$    m3=m2/2
!!$    m4=m3/2
!!$    m5=m4/2
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 1:', m1
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 2:', m2
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 3:', m3
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 4:', m4
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh size 5:', m5
!!$    WRITE(*,*)
!!$
!!$    ALLOCATE(d1(m1,m1))
!!$    ALLOCATE(d2(m2,m2))
!!$    ALLOCATE(d3(m3,m3))
!!$    ALLOCATE(d4(m4,m4))
!!$    ALLOCATE(d5(m5,m5))
!!$
!!$    ! Do CIC binning on each mesh resolution
!!$    ! Note well. These CIC routines assume the volume is periodic
!!$    ! This means you may get some weird edge effects
!!$    CALL CIC2D(y,np,Lsub,d1,m1)
!!$    CALL CIC2D(y,np,Lsub,d2,m2)
!!$    CALL CIC2D(y,np,Lsub,d3,m3)
!!$    CALL CIC2D(y,np,Lsub,d4,m4)
!!$    CALL CIC2D(y,np,Lsub,d5,m5)
!!$
!!$    ! Convert to over-densities (1+delta = rho/rhobar)
!!$    d1=(d1/((Lsub/REAL(m1))**2))/nbar
!!$    d2=(d2/((Lsub/REAL(m2))**2))/nbar
!!$    d3=(d3/((Lsub/REAL(m3))**2))/nbar
!!$    d4=(d4/((Lsub/REAL(m4))**2))/nbar
!!$    d5=(d5/((Lsub/REAL(m5))**2))/nbar
!!$
!!$    ! Smooth density fields
!!$    CALL smooth2D(d1,m1,fcell*Lsub/REAL(m1),Lsub)
!!$    CALL smooth2D(d2,m2,fcell*Lsub/REAL(m2),Lsub)
!!$    CALL smooth2D(d3,m3,fcell*Lsub/REAL(m3),Lsub)
!!$    CALL smooth2D(d4,m4,fcell*Lsub/REAL(m4),Lsub)
!!$    CALL smooth2D(d5,m5,fcell*Lsub/REAL(m5),Lsub)
!!$
!!$    ! Write out density statistics on each mesh
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m1
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Max dens:', MAXVAL(d1)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Min dens:', MINVAL(d1)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m2
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Max dens:', MAXVAL(d2)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Min dens:', MINVAL(d2)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m3
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Max dens:', MAXVAL(d3)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Min dens:', MINVAL(d3)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m4
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Max dens:', MAXVAL(d4)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Min dens:', MINVAL(d4)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Mesh:', m5
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Max dens:', MAXVAL(d5)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Min dens:', MINVAL(d5)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Recommended refinement level:', 1+CEILING(log10(MAXVAL(d1)))
!!$    WRITE(*,*)
!!$
!!$    IF(test) THEN
!!$       ! Write out each mesh if testing
!!$       
!!$       base='density_'
!!$       ext='.dat'
!!$       
!!$       output=number_file_zeroes(base,m1,4,ext)
!!$       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
!!$       CALL write_2D_field_ascii(d1,m1,Lsub,output)
!!$       
!!$       output=number_file_zeroes(base,m2,4,ext)
!!$       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
!!$       CALL write_2D_field_ascii(d2,m2,Lsub,output)
!!$       
!!$       output=number_file_zeroes(base,m3,4,ext)       
!!$       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
!!$       CALL write_2D_field_ascii(d3,m3,Lsub,output)
!!$       
!!$       output=number_file_zeroes(base,m4,4,ext)
!!$       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
!!$       CALL write_2D_field_ascii(d4,m4,Lsub,output)
!!$       
!!$       output=number_file_zeroes(base,m5,4,ext)
!!$       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
!!$       CALL write_2D_field_ascii(d5,m5,Lsub,output)
!!$       
!!$    END IF
!!$
!!$    ! The exact details of the refinement here could certainly be improved
!!$    ! This would involve fiddling with dc and fac
!!$    ! Also all the refinements could be put in a loop
!!$
!!$    ! First refinement
!!$
!!$    ! Allocate d with the size of the second-corsest image
!!$    ALLOCATE(d(m4,m4))
!!$    dc=dcc
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Initial refinement 1+delta:', dc
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Refinement factor:', fac
!!$
!!$    DO i=1,m5
!!$       DO j=1,m5
!!$          IF(d5(i,j)<dc) THEN
!!$             ! Use coarser density
!!$             d(2*i-1,2*j-1)=d5(i,j)
!!$             d(2*i,2*j-1)=d5(i,j)
!!$             d(2*i-1,2*j)=d5(i,j)
!!$             d(2*i,2*j)=d5(i,j)
!!$          ELSE
!!$             ! Use finer density
!!$             d(2*i-1,2*j-1)=d4(2*i-1,2*j)
!!$             d(2*i,2*j-1)=d4(2*i,2*j-1)
!!$             d(2*i-1,2*j)=d4(2*i-1,2*j)
!!$             d(2*i,2*j)=d4(2*i,2*j)
!!$             ! Stops it being below the threshold (looks grainy)
!!$             IF(d(2*i-1,2*j-1)<dc) d(2*i-1,2*j-1)=dc
!!$             IF(d(2*i,2*j-1)<dc)   d(2*i,2*j-1)=dc
!!$             IF(d(2*i-1,2*j)<dc)   d(2*i-1,2*j)=dc
!!$             IF(d(2*i,2*j)<dc)     d(2*i,2*j)=dc
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    ! Second refinement
!!$    d4=d
!!$    DEALLOCATE(d)
!!$    ALLOCATE(d(m3,m3))
!!$    dc=fac*dc
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Second refinement 1+delta:', dc
!!$
!!$    DO i=1,m4
!!$       DO j=1,m4
!!$          IF(d4(i,j)<dc) THEN
!!$             ! Use coarser density
!!$             d(2*i-1,2*j-1)=d4(i,j)
!!$             d(2*i,2*j-1)=d4(i,j)
!!$             d(2*i-1,2*j)=d4(i,j)
!!$             d(2*i,2*j)=d4(i,j)
!!$          ELSE
!!$             ! Use finer density
!!$             d(2*i-1,2*j-1)=d3(2*i-1,2*j)
!!$             d(2*i,2*j-1)=d3(2*i,2*j-1)
!!$             d(2*i-1,2*j)=d3(2*i-1,2*j)
!!$             d(2*i,2*j)=d3(2*i,2*j)
!!$             IF(d(2*i-1,2*j-1)<dc) d(2*i-1,2*j-1)=dc
!!$             IF(d(2*i,2*j-1)<dc)   d(2*i,2*j-1)=dc
!!$             IF(d(2*i-1,2*j)<dc)   d(2*i-1,2*j)=dc
!!$             IF(d(2*i,2*j)<dc)     d(2*i,2*j)=dc
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    ! Third refinement
!!$    d3=d
!!$    DEALLOCATE(d)
!!$    ALLOCATE(d(m2,m2))
!!$    dc=fac*dc
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Third refinement 1+delta:', dc
!!$
!!$    DO i=1,m3
!!$       DO j=1,m3
!!$          IF(d3(i,j)<dc) THEN
!!$             ! Use coarser density
!!$             d(2*i-1,2*j-1)=d3(i,j)
!!$             d(2*i,2*j-1)=d3(i,j)
!!$             d(2*i-1,2*j)=d3(i,j)
!!$             d(2*i,2*j)=d3(i,j)
!!$          ELSE
!!$             ! Use finer density
!!$             d(2*i-1,2*j-1)=d2(2*i-1,2*j)
!!$             d(2*i,2*j-1)=d2(2*i,2*j-1)
!!$             d(2*i-1,2*j)=d2(2*i-1,2*j)
!!$             d(2*i,2*j)=d2(2*i,2*j)
!!$             IF(d(2*i-1,2*j-1)<dc) d(2*i-1,2*j-1)=dc
!!$             IF(d(2*i,2*j-1)<dc)   d(2*i,2*j-1)=dc
!!$             IF(d(2*i-1,2*j)<dc)   d(2*i-1,2*j)=dc
!!$             IF(d(2*i,2*j)<dc)     d(2*i,2*j)=dc
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    ! Fourth refinement
!!$    d2=d
!!$    DEALLOCATE(d)
!!$    ALLOCATE(d(m1,m1))
!!$    dc=fac*dc
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Fourth refinement 1+delta:', dc
!!$
!!$    DO i=1,m2
!!$       DO j=1,m2
!!$          IF(d2(i,j)<dc) THEN
!!$             ! Use coarser density
!!$             d(2*i-1,2*j-1)=d2(i,j)
!!$             d(2*i,2*j-1)=d2(i,j)
!!$             d(2*i-1,2*j)=d2(i,j)
!!$             d(2*i,2*j)=d2(i,j)
!!$          ELSE
!!$             ! Use finer density
!!$             d(2*i-1,2*j-1)=d1(2*i-1,2*j)
!!$             d(2*i,2*j-1)=d1(2*i,2*j-1)
!!$             d(2*i-1,2*j)=d1(2*i-1,2*j)
!!$             d(2*i,2*j)=d1(2*i,2*j)
!!$             IF(d(2*i-1,2*j-1)<dc) d(2*i-1,2*j-1)=dc
!!$             IF(d(2*i,2*j-1)<dc)   d(2*i,2*j-1)=dc
!!$             IF(d(2*i-1,2*j)<dc)   d(2*i-1,2*j)=dc
!!$             IF(d(2*i,2*j)<dc)     d(2*i,2*j)=dc
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Adaptive density field'
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Max density:', MAXVAL(d)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Min density:', MINVAL(d)
!!$    WRITE(*,*)
!!$
!!$    IF(test) THEN
!!$       output='density.dat'
!!$       WRITE(*,*) 'ADAPTIVE_DENSITY: Writing: ', TRIM(output)
!!$       CALL write_2D_field_ascii(d,m1,Lsub,output)
!!$    END IF
!!$
!!$    ! Smooth density field on scales of cells
!!$    CALL smooth2D(d,m1,fcell*Lsub/REAL(m1),Lsub)
!!$
!!$    ! Ensure all cells are positive (note this is to make pretty pictures, not for science)
!!$    DO i=1,m1
!!$       DO j=1,m1
!!$          IF(d(i,j)<0.) d(i,j)=0.
!!$       END DO
!!$    END DO
!!$
!!$    ! Write field statistics
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Smoothed adaptive density field'
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Max density:', MAXVAL(d)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Min density:', MINVAL(d)
!!$    WRITE(*,*)
!!$
!!$    ! Finally write out an ascii file for plotting
!!$    CALL write_2D_field_ascii(d,m1,Lsub,outfile)
!!$    WRITE(*,*) 'ADAPTIVE_DENSITY: Done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE adaptive_density_original

END MODULE simulations
