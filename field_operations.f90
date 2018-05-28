MODULE field_operations

CONTAINS

  REAL FUNCTION random_mode_amplitude(k,L,logk_tab,logPk_tab,nk,use_average)

    !This calculates the Fourier amplitudes of the density field
    USE interpolate
    USE constants
    USE random_numbers
    IMPLICIT NONE
    REAL, INTENT(IN) :: k, L, logk_tab(nk), logPk_tab(nk)
    LOGICAL, INTENT(IN) :: use_average
    INTEGER, INTENT(IN) :: nk
    REAL :: sigma

    !! EXTREME CAUTION: FUDGE FACTOR IN RAYLEIGH !!

    !Sigma parameter in the Rayleigh distribution
    sigma=sqrt(exp(find(log(k),logk_tab,logPk_tab,nk,3,3,2))/(4.*pi*(L*k/(2.*pi))**3))

    IF(use_average) THEN
       !Fixed mode amplitudes       
       random_mode_amplitude=sigma
    ELSE
       !Correctly assigned random mode amplitudes
       random_mode_amplitude=random_Rayleigh(sigma)/sqrt(2.) !sqrt(2) is a FUDGE (something to do with average of Rayleigh?)
    END IF

    !! EXTREME CAUTION: FUDGE FACTOR IN RAYLEIGH !!

  END FUNCTION random_mode_amplitude

  COMPLEX FUNCTION random_complex_phase()

    !Get a complex phase with theta between 0 and 2pi
    USE constants
    USE random_numbers
    IMPLICIT NONE
    REAL :: theta

    theta=random_uniform(0.,twopi)
    random_complex_phase=CMPLX(cos(theta),sin(theta))

  END FUNCTION random_complex_phase

  SUBROUTINE make_Gaussian_random_modes(dk,m,L,logk_tab,logPk_tab,nk,use_average)

    !Uses a tablulated P(k) to make a Gaussian Random Field realisation
    USE fft
    USE random_numbers
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(OUT) :: dk(m,m,m)
    REAL, INTENT(IN) :: logk_tab(nk), logPk_tab(nk), L
    INTEGER, INTENT(IN) :: m, nk
    LOGICAL, INTENT(IN) :: use_average
    INTEGER :: ix, iy, iz, ixx, iyy, izz
    REAL :: kx, ky, kz, k
    REAL :: amp
    COMPLEX :: rot

    dk=(0.d0,0.d0)

    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_MODES: Creating Fourier realisation of Gaussian field'
    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_MODES: Using average mode power:', use_average
    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_MODES: Mesh size:', m
    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_MODES: Box size [Mpc/h]:', L

    !This fills up displacement array in all of k space!
    DO iz=1,m
       DO iy=1,m
          DO ix=1,m

             CALL k_fft(ix,iy,iz,m,kx,ky,kz,k,L)

             IF(ix==1 .AND. iy==1 .AND. iz==1) THEN

                !Set the zero mode to zero
                dk(ix,iy,iz)=0.d0

             ELSE IF(ix==1+m/2 .OR. iy==1+m/2 .OR. iz==1+m/2) THEN 

                !Sets Nyquist modes to 0.!
                !Maybe all modes with mod(k)>k_ny should be set to 0.?!
                !Bridit Falck wrote a 'corner modes' paper about this
                !https://arxiv.org/abs/1610.04862
                dk(ix,iy,iz)=0.d0

             ELSE

                !Get mode amplitudes and phases
                amp=random_mode_amplitude(k,L,logk_tab,logPk_tab,nk,use_average)
                rot=random_complex_phase()

                !Assign values to the density field
                dk(ix,iy,iz)=amp*rot

             END IF

          END DO
       END DO
    END DO

    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_MODES: Done'
    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_MODES: Enforcing Hermiticity'

    !Enforce Hermiticity - probably could save a load of operations above
    DO iz=1,m
       DO iy=1,m
          DO ix=1,m

             ixx=m-ix+2
             iyy=m-iy+2
             izz=m-iz+2

             IF(ix==1) ixx=1
             IF(iy==1) iyy=1
             IF(iz==1) izz=1

             !Do the enforcing
             dk(ix,iy,iz)=CONJG(dk(ixx,iyy,izz))

          END DO
       END DO
    END DO

    !Is this a good idea?
    dk=CONJG(dk)

    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_MODES: Hermitian field generated'
    WRITE(*,*)

  END SUBROUTINE make_Gaussian_random_modes

  SUBROUTINE make_Gaussian_random_field(d,m,L,logk_tab,logPk_tab,nk,use_average)

    !Uses a tablulated P(k) to make a Gaussian Random Field realisation
    USE fft
    USE random_numbers
    IMPLICIT NONE
    REAL, INTENT(OUT) :: d(m,m,m)
    DOUBLE COMPLEX :: dk(m,m,m), dk_new(m,m,m)
    REAL, INTENT(IN) :: logk_tab(nk), logPk_tab(nk), L
    INTEGER, INTENT(IN) :: m, nk
    LOGICAL, INTENT(IN) :: use_average

    CALL make_Gaussian_random_modes(dk,m,L,logk_tab,logPk_tab,nk,use_average)

    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_FIELD: Transform to real space'

    !FT the displacement field from k-space to real space!
    CALL fft3(dk,dk_new,m,m,m,1)
    dk=dk_new
    d=REAL(REAL(dk))

    WRITE(*,*) 'MAKE_GAUSSIAN_RANDOM_FIELD: Done'
    WRITE(*,*)

  END SUBROUTINE make_Gaussian_random_field

  SUBROUTINE generate_displacement_fields(f,m,L,logk_tab,logPk_tab,nk,use_average)

    USE fft
    USE array_operations
    !USE statistics
    IMPLICIT NONE
    REAL, INTENT(OUT) :: f(3,m,m,m)
    INTEGER, INTENT(IN) :: m, nk
    REAL, INTENT(IN) :: L, logk_tab(nk), logPk_tab(nk)
    LOGICAL, INTENT(IN) :: use_average
    DOUBLE COMPLEX :: d(m,m,m), dk(m,m,m), fk(3,m,m,m)
    INTEGER :: i, ix, iy, iz
    REAL :: kx, ky, kz, k

    CALL make_Gaussian_random_modes(dk,m,L,logk_tab,logPk_tab,nk,use_average)

    WRITE(*,*) 'GENERATE_DISPLACEMENT_FIELDS: Creating realisation of displacement field'

    !This fills up displacement array in all of k space!
    DO iz=1,m
       DO iy=1,m
          DO ix=1,m

             !Get the wave vectors
             CALL k_fft(ix,iy,iz,m,kx,ky,kz,k,L)

             IF(ix==1 .AND. iy==1 .AND. iz==1) THEN

                !Set the DC mode to zero
                fk(1,ix,iy,iz)=0.d0
                fk(2,ix,iy,iz)=0.d0
                fk(3,ix,iy,iz)=0.d0

             ELSE IF(ix==1+m/2 .OR. iy==1+m/2 .OR. iz==1+m/2) THEN 

                !Sets Nyquist modes to 0.!
                !Maybe all modes with mod(k)>k_ny should be set to 0.?!
                !Bridit Falck wrote a 'corner modes' paper about this
                !https://arxiv.org/abs/1610.04862
                fk(1,ix,iy,iz)=0.d0
                fk(2,ix,iy,iz)=0.d0
                fk(3,ix,iy,iz)=0.d0

             ELSE

                !Assign values to the displacement field
                fk(1,ix,iy,iz)=(0.d0,-1.d0)*dk(ix,iy,iz)*kx/k**2
                fk(2,ix,iy,iz)=(0.d0,-1.d0)*dk(ix,iy,iz)*ky/k**2
                fk(3,ix,iy,iz)=(0.d0,-1.d0)*dk(ix,iy,iz)*kz/k**2

             END IF

          END DO
       END DO
    END DO
    WRITE(*,*) 'GENERATE_DISPLACEMENT_FIELDS: Fourier displacement fields done'

    DO i=1,3
       dk=fk(i,:,:,:)
       CALL fft3(dk,d,m,m,m,1)
       f(i,:,:,:)=REAL(REAL(d(:,:,:)))
    END DO

    WRITE(*,*) 'GENERATE_DISPLACEMENT_FIELDS: Minimum 1D displacement [Mpc/h]:', MINVAL(f)
    WRITE(*,*) 'GENERATE_DISPLACEMENT_FIELDS: Maximum 1D displacement [Mpc/h]:', MAXVAL(f)
    WRITE(*,*) 'GENERATE_DISPLACEMENT_FIELDS: Real-space displacement fields generated'
    WRITE(*,*)

  END SUBROUTINE generate_displacement_fields

  SUBROUTINE read_field(d,m,infile)

    !Read in a binary 'field' file
    USE array_operations
    USE statistics
    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL, INTENT(OUT) :: d(m,m,m)
    INTEGER, INTENT(IN) :: m

    !Output unformatted data
    WRITE(*,*) 'READ_FIELD: Binary input: ', TRIM(infile)
    WRITE(*,*) 'READ_FIELD: Mesh size:', m
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) d
    CLOSE(7)
    WRITE(*,*) 'READ_FIELD: Minval:', MINVAL(d)
    WRITE(*,*) 'READ_FIELD: Maxval:', MAXVAL(d)
    WRITE(*,*) 'READ_FIELD: Average:', mean(splay(d,m,m,m),m**3)
    WRITE(*,*) 'READ_FIELD: Variance:', variance(splay(d,m,m,m),m**3)
    WRITE(*,*) 'READ_FIELD: Done'
    WRITE(*,*)

  END SUBROUTINE read_field

  SUBROUTINE read_field8(d,m,infile)

    USE array_operations
    USE statistics
    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL, INTENT(OUT) :: d(m,m,m)
    DOUBLE PRECISION :: d8(m,m,m)
    INTEGER, INTENT(IN) :: m

    !Input unformatted data
    WRITE(*,*) 'READ_FIELD: Binary input: ', TRIM(infile)
    WRITE(*,*) 'READ_FIELD: Mesh size:', m
    OPEN(7,file=infile,form='unformatted',access='stream',status='old')
    READ(7) d8
    CLOSE(7)
    d=REAL(d8)
    WRITE(*,*) 'READ_FIELD: Minval:', MINVAL(d)
    WRITE(*,*) 'READ_FIELD: Maxval:', MAXVAL(d)
    WRITE(*,*) 'READ_FIELD: Average:', REAL(mean(splay(d,m,m,m),m**3))
    WRITE(*,*) 'READ_FIELD: Variance:', REAL(variance(splay(d,m,m,m),m**3))
    WRITE(*,*) 'READ_FIELD: Done'
    WRITE(*,*)

  END SUBROUTINE read_field8

  !Used to be called write_field
  SUBROUTINE write_field_binary(d,m,outfile)

    !Write out a binary 'field' file
    IMPLICIT NONE    
    REAL, INTENT(IN) :: d(m,m,m)
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=256), INTENT(IN) :: outfile

    WRITE(*,*) 'WRITE_FIELD_BINARY: Binary output: ', TRIM(outfile)
    WRITE(*,*) 'WRITE_FIELD_BINARY: Mesh size:', m
    WRITE(*,*) 'WRITE_FIELD_BINARY: Minval:', MINVAL(d)
    WRITE(*,*) 'WRITE_FIELD_BINARY: Maxval:', MAXVAL(d)
    WRITE(*,*) 'WRITE_FIELD_BINARY: Using new version with access=stream'
    OPEN(7,file=outfile,form='unformatted',access='stream',status='replace')
    WRITE(7) d
    CLOSE(7)
    WRITE(*,*) 'WRITE_3D_FIELD_BINARY: Done'
    WRITE(*,*)

  END SUBROUTINE write_field_binary

  !Used to be called print_2D_field
  SUBROUTINE write_2D_field_ascii(d,m,L,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d(m,m), L
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=256), INTENT(IN) :: outfile
    INTEGER :: i, j
    REAL :: x, y

    WRITE(*,*) 'WRITE_2D_FIELD_ASCII: Writing to: ', TRIM(outfile)

    OPEN(8,file=outfile)
    DO j=1,m
       DO i=1,m

          x=L*(REAL(i)-0.5)/REAL(m)
          y=L*(REAL(j)-0.5)/REAL(m)

          WRITE(8,*) x, y, d(i,j)

       END DO
    END DO
    CLOSE(8)

    WRITE(*,*) 'WRITE_2D_FIELD_ASCII: Done'
    WRITE(*,*)

  END SUBROUTINE write_2D_field_ascii

  !Used to be called print_projected_field
  SUBROUTINE write_3D_field_projection_ascii(d,m,L,nz,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d(m,m,m), L
    INTEGER, INTENT(IN) :: m, nz
    CHARACTER(len=256), INTENT(IN) :: outfile
    INTEGER :: i, j, k
    REAL :: x, y
    REAL :: sum

    WRITE(*,*) 'WRITE_3D_FIELD_PROJECTION_ASCII: Writing to: ', TRIM(outfile)
    WRITE(*,*) 'WRITE_3D_FIELD_PROJECTION_ASCII: Cells projecting:', nz

    OPEN(8,file=outfile)
    DO j=1,m
       DO i=1,m

          x=L*(REAL(i)-0.5)/REAL(m)
          y=L*(REAL(j)-0.5)/REAL(m)

          sum=0.
          DO k=1,nz
             sum=sum+d(i,j,k)
          END DO
          sum=sum/float(nz)

          WRITE(8,*) x, y, sum

       END DO
    END DO
    CLOSE(8)

    WRITE(*,*) 'WRITE_3D_FIELD_PROJECTION_ASCII: Done'
    WRITE(*,*)

  END SUBROUTINE write_3D_field_projection_ascii

  SUBROUTINE compress_field(d,ds,m)

    !Shrinks a 3D field size by a factor of 2
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m
    REAL, INTENT(IN) :: d(m,m,m)
    REAL, ALLOCATABLE, INTENT(OUT) :: ds(:,:,:)
    INTEGER :: i, j, k

    !Allocate the small array
    ALLOCATE(ds(m/2,m/2,m/2))
    ds=0.

    !Fill up the small array by summing blocks of 8 from the larger array
    DO k=1,m/2
       DO j=1,m/2
          DO i=1,m/2
             ds(i,j,k)=ds(i,j,k)+d(2*i-1,2*j-1,2*k-1)
             ds(i,j,k)=ds(i,j,k)+d(2*i,2*j-1,2*k-1)
             ds(i,j,k)=ds(i,j,k)+d(2*i-1,2*j,2*k-1)
             ds(i,j,k)=ds(i,j,k)+d(2*i-1,2*j-1,2*k)
             ds(i,j,k)=ds(i,j,k)+d(2*i,2*j,2*k-1)
             ds(i,j,k)=ds(i,j,k)+d(2*i,2*j-1,2*k)
             ds(i,j,k)=ds(i,j,k)+d(2*i-1,2*j,2*k)
             ds(i,j,k)=ds(i,j,k)+d(2*i,2*j,2*k)
          END DO
       END DO
    END DO

    !Divide by the number of blocks that are being averaged over
    ds=ds/8.

  END SUBROUTINE compress_field

  SUBROUTINE sharpen(d,m,L,ibin)

    USE fft
    IMPLICIT NONE
    INTEGER :: m
    REAL, INTENT(INOUT) :: d(m,m,m)
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: ibin
    DOUBLE COMPLEX :: dk(m,m,m), dkout(m,m,m)

    !ibin = 1 NGP
    !ibin = 2 CIC

    WRITE(*,*) 'SHARPEN: Correcting for binning by sharpening field'
    WRITE(*,*) 'SHARPEN: Mesh size:', m

    dk=d

    CALL fft3(dk,dkout,m,m,m,-1)
    dk=dkout

    CALL sharpen_k(dk,m,L,ibin)

    CALL fft3(dk,dkout,m,m,m,1)
    dk=dkout

    d=REAL(REAL(dk))/REAL(m**3)

    WRITE(*,*) 'SHARPEN: Sharpening complete'
    WRITE(*,*)

  END SUBROUTINE sharpen

  SUBROUTINE sharpen_k(dk,m,L,ibin)

    USE special_functions
    USE fft
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(INOUT) :: dk(m,m,m)
    INTEGER, INTENT(IN) :: m, ibin
    REAL, INTENT(IN) :: L
    INTEGER :: i, j, k
    REAL :: kx, ky, kz, kmod
    REAL :: kxh, kyh, kzh
    REAL :: fcx, fcy, fcz, fcorr

    !Now correct for binning!
    DO k=1,m
       DO j=1,m
          DO i=1,m

             CALL k_fft(i,j,k,m,kx,ky,kz,kmod,L)

             kxh=L*kx/(2.*REAL(m))
             kyh=L*ky/(2.*REAL(m))
             kzh=L*kz/(2.*REAL(m))

             fcx=sinc(kxh)
             fcy=sinc(kyh)
             fcz=sinc(kzh)

             IF(ibin==1) THEN
                fcorr=fcx*fcy*fcz
             ELSE IF(ibin==2) THEN
                fcorr=(fcx*fcy*fcz)**2
             ELSE
                STOP 'SHARPEN_K: Error, ibin specified incorrectly'
             END IF

             dk(i,j,k)=dk(i,j,k)/fcorr

          END DO
       END DO
    END DO

  END SUBROUTINE sharpen_k

  SUBROUTINE smooth2D(d,m,r,L)

    !arr(n,n): input array of size n x n
    !r: smoothing scale in Mpc/h
    !L: box size in Mpc/h
    USE fft
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: d(m,m)
    REAL, INTENT(IN) :: r, L
    INTEGER, INTENT(IN) :: m
    REAL :: kx, ky, kz, k
    DOUBLE COMPLEX :: dc(m,m), dk(m,m)
    INTEGER :: i, j

    WRITE(*,*) 'SMOOTH2D: Smoothing array'
    WRITE(*,*) 'SMOOTH2D: Assuming array is periodic'
    WRITE(*,*) 'SMOOTH2D: Smoothing scale [Mpc/h]:', r

    dc=d
    CALL fft2(dc,dk,m,m,-1)

    DO j=1,m
       DO i=1,m
          CALL k_fft(i,j,1,m,kx,ky,kz,k,L)
          dk(i,j)=dk(i,j)*exp(-((k*r)**2)/2.)
       END DO
    END DO

    !Normalise post Fourier transform
    CALL fft2(dk,dc,m,m,1)
    dc=dc/REAL(m**2)
    d=REAL(REAL(dc))

    WRITE(*,*) 'SMOOTH2D: Done'
    WRITE(*,*)

  END SUBROUTINE smooth2D

!!$  SUBROUTINE smooth2D_nonperiodic(arr,n,r,L)
!!$
!!$    !arr(n,n): input array of size n x n
!!$    !r: smoothing scale in Mpc/h
!!$    !L: box size in Mpc/h
!!$    USE fft
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(INOUT) :: arr(n,n)
!!$    REAL, INTENT(IN) :: r, L
!!$    REAL :: kx, ky, kz, k
!!$    DOUBLE COMPLEX, ALLOCATABLE :: ac(:,:), ac_out(:,:)
!!$    INTEGER :: i, j, m
!!$
!!$    INTEGER, PARAMETER :: pad=2 !Padding because we generally will not be continuous
!!$
!!$    WRITE(*,*) 'SMOOTH2D: Smoothing array'
!!$    WRITE(*,*) 'SMOOTH2D: Smoothing scale [Mpc/h]:', r
!!$
!!$    !For padding, I cant imagine that x2 would ever be insufficient!
!!$    m=pad*n
!!$
!!$    ALLOCATE(ac(m,m),ac_out(m,m))
!!$
!!$    !Not sure if this is necessary
!!$    ac=(0.d0,0.d0)
!!$    ac_out=(0.d0,0.d0)
!!$
!!$    !Put image into complex array, padded with zeroes where image is not!
!!$    DO j=1,n
!!$       DO i=1,n
!!$          ac(i,j)=arr(i,j)
!!$       END DO
!!$    END DO
!!$
!!$    CALL fft2(ac,ac_out,m,m,-1)
!!$    ac=ac_out
!!$
!!$    !Smoothing length in terms of image(m x m) size!
!!$    !r=pix/float(m)
!!$    DO j=1,m
!!$       DO i=1,m
!!$          CALL k_fft(i,j,1,m,kx,ky,kz,k,pad*L)
!!$          ac(i,j)=ac(i,j)*exp(-((k*r)**2)/2.)
!!$       END DO
!!$    END DO
!!$
!!$    CALL fft2(ac,ac_out,m,m,1)
!!$    ac=ac_out
!!$    DEALLOCATE(ac_out)
!!$
!!$    !Normalise post Fourier transform!
!!$    ac=ac/REAL(m**2)
!!$
!!$    !Retrieve smooth image from complex array!
!!$    !Need a loop because arr and ac will have different sizes
!!$    DO j=1,n
!!$       DO i=1,n
!!$          arr(i,j)=REAL(REAL(ac(i,j)))
!!$       END DO
!!$    END DO
!!$
!!$    WRITE(*,*) 'SMOOTH2D: Done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE smooth2D_nonperiodic

  SUBROUTINE smooth3D(arr,m,r,L)

    USE fft
    USE special_functions
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: arr(m,m,m)
    REAL, INTENT(IN) :: r, L
    INTEGER, INTENT(IN) :: m
    REAL :: kx, ky, kz, kmod
    DOUBLE COMPLEX :: ac(m,m,m), ac_out(m,m,m)
    INTEGER :: i, j, k

    WRITE(*,*) 'Smoothing array'

    ac=arr

    !Move to Fourier space
    CALL fft3(ac,ac_out,m,m,m,-1)
    ac=ac_out

    DO k=1,m
       DO j=1,m
          DO i=1,m
             CALL k_fft(i,j,k,m,kx,ky,kz,kmod,L)
             ac(i,j,k)=ac(i,j,k)*sinc(kx*r/2.)*sinc(ky*r/2.)*sinc(kz*r/2.)
          END DO
       END DO
    END DO

    !Move back to real space
    CALL fft3(ac,ac_out,m,m,m,1)
    ac=ac_out

    !Normalise post Fourier transform!
    ac=ac/(REAL(m)**3)

    arr=REAL(REAL(ac))

    WRITE(*,*) 'Done'
    WRITE(*,*)

  END SUBROUTINE smooth3D

  SUBROUTINE add_to_stack_3D(x,stack,Ls,ms,back,Lb,mb)

    !Adds some points in a density field to a stack
    IMPLICIT NONE
    INTEGER :: i, j, k, is(3), ib(3), d
    INTEGER, INTENT(IN) :: ms, mb
    REAL, INTENT(IN) :: x(3), Ls, Lb
    REAL, INTENT(INOUT) :: stack(ms,ms,ms)
    REAL, INTENT(IN) :: back(mb,mb,mb)
    REAL :: xb(3)

    !Assumes the background field is periodic
    !'stack' should have been previously allocated
    !'stack' should be set to zero before using this subroutine
    !'*s' variables refer to the stacked field
    !'*_back' variables refer to the background field

    !Loop over cells on stacking mesh
    DO i=1,ms
       DO j=1,ms
          DO k=1,ms

             !Set the stack integer array
             is(1)=i
             is(2)=j
             is(3)=k

             DO d=1,3

                !Get coordinates of position on the stack
                !This changes coordiantes from stack to simulation coordinates
                xb(d)=x(d)+Ls*(0.5+float(is(d)-1))/float(ms)-Ls/2.

                !Bring the coordinates back into the simulation box if they are outside
                IF(xb(d)<=0.) THEN
                   xb(d)=xb(d)+Lb
                ELSE IF(xb(d)>Lb) THEN
                   xb(d)=xb(d)-Lb
                END IF

                !Find the integer coordinates of mesh cell in the background mesh
                !This is just an NGP-type scheme. Could/should be improved?
                ib(d)=CEILING(float(mb)*xb(d)/Lb)

             END DO

             !Add the value to the stack
             !Should there be a volume factor here?
             stack(is(1),is(2),is(3))=stack(is(1),is(2),is(3))+back(ib(1),ib(2),ib(3))

          END DO
       END DO
    END DO

  END SUBROUTINE add_to_stack_3D

  SUBROUTINE project_3D_to_2D(d3d,d2d,m)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d3d(m,m,m)
    INTEGER, INTENT(IN) :: m
    REAL, INTENT(OUT) :: d2d(m,m)
    INTEGER :: i, j, k

    WRITE(*,*) 'PROJECT_3D_TO_2D: Projecting 3D stack into 2D'
    d2d=0.
    DO i=1,m
       DO j=1,m
          DO k=1,m
             d2d(i,j)=d2d(i,j)+d3d(i,j,k)
          END DO
       END DO
    END DO
    d2d=d2d/REAL(m)
    WRITE(*,*) 'PROJECT_3D_TO_2D: Minimum value of 2D stack:', MINVAL(d2d)
    WRITE(*,*) 'PROJECT_3D_TO_2D: Maximum value of 2D stack:', MAXVAL(d2d)
    WRITE(*,*) 'PROJECT_3D_TO_2D: Done'
    WRITE(*,*)

  END SUBROUTINE project_3D_to_2D

  SUBROUTINE clip(d,m1,m2,m3,d0,talk)

    USE statistics
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: d(:,:,:)
    REAL, INTENT(IN) :: d0
    INTEGER, INTENT(IN) :: m1, m2, m3
    LOGICAL, INTENT(IN) :: talk
    REAL :: var1, av1, max1, var2, av2, max2
    INTEGER :: i, j, k

    IF(talk) THEN
       WRITE(*,*) 'CLIP: Clipping density field'
       WRITE(*,*) 'CLIP: Threshold:', d0
       WRITE(*,*) 'CLIP: Mesh:', m1, m2, m3
    END IF

    av1=mean(splay(d,m1,m2,m3),m1*m2*m3)
    var1=variance(splay(d,m1,m2,m3),m1*m2*m3)
    max1=MAXVAL(d)

    IF(talk) THEN
       WRITE(*,*) 'CLIP: Average over-density pre-clipping:', av1
       WRITE(*,*) 'CLIP: Variance in over-density pre-clipping:', var1
       WRITE(*,*) 'CLIP: Maximum density pre-clipping:', max1
    END IF

    !    dep=0.25*(1.+erf(d0/(sqrt(2.*var1))))**2.
    !    IF(talk==1) WRITE(*,*) 'Expected large-scale power depletion factor:', dep

    !Now do the clipping
    DO k=1,m3
       DO j=1,m2
          DO i=1,m1
             IF(d(i,j,k)>d0) d(i,j,k)=d0
          END DO
       END DO
    END DO

    IF(talk) WRITE(*,*) 'CLIP: Density field clipped'

    av2=mean(splay(d,m1,m2,m3),m1*m2*m3)
    var2=variance(splay(d,m1,m2,m3),m1*m2*m3)
    max2=MAXVAL(d)

    IF(talk) THEN
       WRITE(*,*) 'CLIP: Average over-density post-clipping:', av2
       WRITE(*,*) 'CLIP: Variance in over-density post-clipping:', var2
       WRITE(*,*) 'CLIP: Maximum density post-clipping:', max2
       WRITE(*,*)
    END IF

  END SUBROUTINE clip

  SUBROUTINE anticlip(d,m1,m2,m3,d0,talk)

    USE statistics
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: d(m1,m2,m3)
    INTEGER, INTENT(IN) :: m1, m2, m3
    REAL, INTENT(IN) :: d0
    LOGICAL, INTENT(IN) :: talk
    REAL :: var1, av1, min1, var2, av2, min2
    INTEGER :: i, j, k, m

    IF(talk) THEN
       WRITE(*,*) 'Anti-clipping over-density field'
       WRITE(*,*) 'Threshold:', d0
       WRITE(*,*) 'Mesh:', m
    END IF

    av1=mean(splay(d,m1,m2,m3),m1*m2*m3)
    var1=variance(splay(d,m1,m2,m3),m1*m2*m3)
    min1=MINVAL(d)

    IF(talk) THEN
       WRITE(*,*) 'Average over-density pre-clipping:', av1
       WRITE(*,*) 'Variance in over-density pre-clipping:', var1
       WRITE(*,*) 'Minimum over-density pre-clipping:', min1
    END IF

    !    dep=0.25*(1.+erf(d0/(sqrt(2.*var1))))**2.
    !    IF(talk==1) WRITE(*,*) 'Expected large-scale power depletion factor:', dep

    !Now do the clipping
    DO k=1,m
       DO j=1,m
          DO i=1,m
             IF(d(i,j,k)<d0) d(i,j,k)=d0
          END DO
       END DO
    END DO

    IF(talk) WRITE(*,*) 'Over-density field clipped'

    av2=mean(splay(d,m1,m2,m3),m1*m2*m3)
    var2=variance(splay(d,m1,m2,m3),m1*m2*m3)
    min2=MINVAL(d)

    IF(talk) THEN
       WRITE(*,*) 'Average over-density post-clipping:', av2
       WRITE(*,*) 'Variance in over-density post-clipping:', var2
       WRITE(*,*) 'Minimum over-density post-clipping:', min2
       WRITE(*,*)
    END IF

  END SUBROUTINE anticlip

  FUNCTION empty_cells(d,m)

    IMPLICIT NONE
    INTEGER :: empty_cells
    REAL, INTENT(IN) :: d(m,m,m)
    INTEGER, INTENT(IN) :: m
    INTEGER*8 :: sum
    INTEGER :: i, j, k

    sum=0
    DO k=1,m
       DO j=1,m
          DO i=1,m
             IF(d(i,j,k)==0.) THEN
                sum=sum+1
             END IF
          END DO
       END DO
    END DO

    empty_cells=INT(sum)

  END FUNCTION empty_cells

  !SUBROUTINE pk(dk1,dk2,m,L,kmin,kmax,bins,k,pow,nbin)
  SUBROUTINE compute_power_spectrum(dk1,dk2,m,L,kmin,kmax,bins,k,pow,nbin)

    !Takes in a dk(m,m,m) array and computes the power spectrum
    !dk1 - Fourier components of field 1
    !dk2 - field 2
    !m - mesh size for fields
    !L - box size in Mpc/h
    !kmin - minimum wavenumber
    !kmax - maximum wavenumber
    !bins - number of k bins
    USE table_integer
    USE constants
    USE array_operations
    USE fft
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(IN) :: dk1(m,m,m), dk2(m,m,m)
    REAL, INTENT(OUT) :: pow(bins), k(bins)
    INTEGER, INTENT(OUT) :: nbin(bins)
    INTEGER, INTENT(IN) :: m, bins
    REAL, INTENT(IN) :: L, kmin, kmax
    INTEGER :: i, ix, iy, iz, n
    REAL :: kx, ky, kz, kmod  
    REAL, ALLOCATABLE :: kbin(:)  
    DOUBLE PRECISION :: pow8(bins), k8(bins)    
    INTEGER*8 :: nbin8(bins)

    WRITE(*,*) 'PK: Computing isotropic power spectrum'

    !Set summation variables to 0.d0
    k8=0.d0
    pow8=0.d0
    nbin8=0

    WRITE(*,*) 'PK: Binning power'
    WRITE(*,*) 'PK: Mesh:', m
    WRITE(*,*) 'PK: Bins:', bins
    WRITE(*,*) 'PK: k_min [h/Mpc]:', kmin
    WRITE(*,*) 'PK: k_max [h/Mpc]:', kmax

    !Fill array of k bins
    CALL fill_array(log(kmin),log(kmax),kbin,bins+1)
    kbin=exp(kbin)

    !Explicitly extend the first and last bins to be sure to include *all* modes
    !This is necessary due to rounding errors!
    kbin(1)=kbin(1)*0.999
    kbin(bins+1)=kbin(bins+1)*1.001    

    !Loop over all elements of dk
    DO iz=1,m
       DO iy=1,m
          DO ix=1,m

             !Cycle for the zero mode (k=0)
             IF(ix==1 .AND. iy==1 .AND. iz==1) CYCLE

             CALL k_fft(ix,iy,iz,m,kx,ky,kz,kmod,L)

             !Find integer 'n' in bins from place in table
             IF(kmod>=kbin(1) .AND. kmod<=kbin(bins+1)) THEN
                n=select_table_integer(kmod,kbin,bins+1,3)
                IF(n<1 .OR. n>bins) THEN
                   CYCLE
                ELSE
                   k8(n)=k8(n)+kmod
                   pow8(n)=pow8(n)+REAL(dk1(ix,iy,iz)*CONJG(dk2(ix,iy,iz)))               
                   nbin8(n)=nbin8(n)+1
                END IF
             END IF

          END DO
       END DO
    END DO

    !Now create the power spectrum and k array
    DO i=1,bins
       k(i)=sqrt(kbin(i+1)*kbin(i))
       IF(nbin8(i)==0) THEN
          !k(i)=sqrt(kbin(i+1)*kbin(i))       
          pow8(i)=0.
       ELSE
          !k(i)=k8(i)/float(nbin8(i))
          pow8(i)=pow8(i)/float(nbin8(i))
          pow8(i)=pow8(i)*((L*k(i))**3.)/(2.*pi**2.)
       END IF
    END DO

    !Do this to account for m^3 factors in P(k)
    pow=REAL(pow8/(DBLE(m)**6))

    !Divide by 2 because up to now we have double count Hermitian conjugates
    nbin=INT(nbin8)/2

    WRITE(*,*) 'PK: Power computed'
    WRITE(*,*) 

  END SUBROUTINE compute_power_spectrum

  SUBROUTINE compute_power_spectrum_pole(d,L,ipole,iz,kmin,kmax,bins,kval,pow,nbin)

    USE constants
    USE special_functions
    USE fft
    IMPLICIT NONE
    INTEGER :: i, j, k, m, n
    REAL :: kx, ky, kz, kmod, mu
    REAL, INTENT(IN) :: kmin, kmax
    REAL, INTENT(IN) :: L
    REAL :: kbin(bins+1)
    REAL, INTENT(OUT) :: pow(bins), kval(bins)
    DOUBLE PRECISION :: pow8(bins), kval8(bins)
    INTEGER, INTENT(OUT) :: nbin(bins)
    INTEGER*8 :: nbin8(bins)
    INTEGER, INTENT(IN) :: iz, ipole, bins
    DOUBLE COMPLEX, INTENT(IN) :: d(:,:,:)

    WRITE(*,*) 'Computing isotropic power spectrum'

    kval=0.
    pow=0.
    nbin=0

    kval8=0.d0
    pow8=0.d0
    nbin8=0

    WRITE(*,*) 'Binning power'
    WRITE(*,*) 'Bins:', bins
    WRITE(*,*) 'k_min:', kmin
    WRITE(*,*) 'k_max:', kmax

    !Log-spaced bins
    DO i=1,bins+1
       kbin(i)=exp(log(kmin)+log(kmax/kmin)*float(i-1)/float(bins))
    END DO

    !Explicitly extend the bins to be sure to include all modes
    !This is necessary due to rounding errors!
    !    kbin(1)=kbin(1)*0.999
    !    kbin(bins+1)=kbin(bins+1)*1.001

    m=SIZE(d(:,1,1))

    WRITE(*,*) 'Mesh:', m

    DO k=1,m
       DO j=1,m
          DO i=1,m

             IF(i==1 .AND. j==1 .AND. k==1) CYCLE

             CALL k_fft(i,j,k,m,kx,ky,kz,kmod,L)

             IF(iz==1) THEN
                mu=kx/kmod
             ELSE IF(iz==2) THEN
                mu=ky/kmod
             ELSE IF(iz==3) THEN
                mu=kz/kmod
             END IF

             !             DO o=1,bins
             !                IF(kmod>=kbin(o) .AND. kmod<=kbin(o+1)) THEN
             !                   pow8(o)=pow8(o)+(ABS(d(i,j,k))**2.)*legendre(ipole,mu)*(2.*float(ipole)+1.)!/2.
             !                   kval(o)=kval(o)+kmod
             !                   nbin8(o)=nbin8(o)+1
             !                   EXIT
             !                END IF
             !             END DO

             !Find integer automatically from place in table. Assumes log-spaced bins
             !Recently implemented (27/08/15) so could be a source of bugs
             !Differences will appear due to k modes that are on the boundary
             n=1+FLOOR(float(bins)*log(kmod/kmin)/log(kmax/kmin))
             IF(n<1 .OR. n>bins) THEN
                CYCLE
             ELSE
                pow8(n)=pow8(n)+(ABS(d(i,j,k))**2.)*Legendre_polynomial(ipole,mu)*(2.*float(ipole)+1.)
                kval8(n)=kval8(n)+kmod
                nbin8(n)=nbin8(n)+1
             END IF

          END DO
       END DO
    END DO

    DO i=1,bins
       !       kval(i)=(kbin(i+1)+kbin(i))/2.
       !       kval(i)=sqrt(kbin(i+1)*kbin(i))
       IF(nbin8(i)==0) THEN
          kval(i)=(kbin(i+1)+kbin(i))/2.
          !       kval(i)=sqrt(kbin(i+1)*kbin(i))
          pow8(i)=0.
       ELSE
          !          kval(i)=(kbin(i+1)+kbin(i))/2.
          kval(i)=REAL(kval8(i))/REAL(nbin8(i))
          pow8(i)=pow8(i)/float(nbin8(i))
          pow8(i)=pow8(i)*((L*kval(i))**3.)/(2.*pi**2.)
       END IF
    END DO

    pow=REAL(pow8)/(REAL(m)**6)

    !Divide by 2 because double count Hermitian conjugates
    nbin=INT(nbin8)/2

    WRITE(*,*) 'Power computed'
    WRITE(*,*) 

  END SUBROUTINE compute_power_spectrum_pole

  SUBROUTINE compute_power_spectrum_rsd(d,L,kmin,kmax,bins,kv,mu,pow,nbin,iz)

    USE constants
    USE fft
    IMPLICIT NONE
    INTEGER :: i, j, k, m, ii, jj, bins, iz
    REAL :: kx, ky, kz, kmod, L, kmin, kmax, a, b, mus
    REAL :: pow(bins,bins), kv(bins), kbin(bins+1), mu(bins), mubin(bins+1)
    DOUBLE PRECISION :: pow8(bins,bins)
    INTEGER :: nbin(bins,bins)
    INTEGER*8 :: nbin8(bins,bins)
    DOUBLE COMPLEX :: d(:,:,:)

    WRITE(*,*) 'Computing RSD power spectrum'

    kbin=0.
    mubin=0.
    kv=0.
    mu=0.
    pow=0.
    nbin=0

    pow8=0.d0
    nbin8=0

    WRITE(*,*) 'Binning power'
    WRITE(*,*) 'Bins:', bins
    WRITE(*,*) 'k_min:', kmin
    WRITE(*,*) 'k_max:', kmax

    a=kmin
    b=kmax

    a=log10(a)
    b=log10(b)

    DO i=1,bins+1
       kbin(i)=a+(b-a)*float(i-1)/float(bins)
    END DO

    DO i=1,bins+1
       mubin(i)=float(i-1)/float(bins)
    END DO

    DO i=1,bins
       kv(i)=(kbin(i)+kbin(i+1))/2.
       mu(i)=(mubin(i)+mubin(i+1))/2.
    END DO

    kbin=10.**kbin
    kv=10.**kv

    !Explicitly extend the bins to be sure to include all modes
    !This is necessary due to rounding errors!
    kbin(1)=kbin(1)*0.999
    kbin(bins+1)=kbin(bins+1)*1.001
    mubin(1)=-0.001
    mubin(bins+1)=1.001

    m=SIZE(d(:,1,1))

    WRITE(*,*) 'Mesh:', m

    DO k=1,m
       DO j=1,m
          DO i=1,m

             IF(i==1 .AND. j==1 .AND. k==1) CYCLE

             CALL k_fft(i,j,k,m,kx,ky,kz,kmod,L)

             IF(iz==1) THEN
                mus=kx/kmod
             ELSE IF(iz==2) THEN
                mus=ky/kmod
             ELSE IF(iz==3) THEN
                mus=kz/kmod
             ELSE
                STOP 'COMPUTE_POWER_SPECTRUM_RSD: Error, iz not specified correctly'
             END IF

             mus=ABS(mus)

             !             WRITE(*,*) mus, kbin(1), kbin(2), kmod
             !             IF(i==10) STOP

             DO jj=1,bins
                IF(kmod>=kbin(jj) .AND. kmod<=kbin(jj+1)) THEN                
                   DO ii=1,bins
                      IF(mus>=mubin(ii) .AND. mus<=mubin(ii+1)) THEN
                         pow8(ii,jj)=pow8(ii,jj)+ABS(d(i,j,k))**2.
                         nbin8(ii,jj)=nbin8(ii,jj)+1
                         EXIT
                      END IF
                   END DO
                   EXIT
                END IF
             END DO

          END DO
       END DO
    END DO

    DO jj=1,bins
       DO ii=1,bins
          IF(nbin8(ii,jj)==0) THEN
             pow8(ii,jj)=0.
          ELSE
             pow8(ii,jj)=pow8(ii,jj)/float(nbin8(ii,jj))
             pow8(ii,jj)=pow8(ii,jj)*((L*kv(jj))**3.)/(2.*pi**2.)
          END IF
       END DO
    END DO

    pow=REAL(pow8)/(REAL(m)**6)

    !Divide by 2 because double count Hermitian conjugates
    nbin=INT(nbin8)/2

    WRITE(*,*) 'Power computed'
    WRITE(*,*) 

  END SUBROUTINE compute_power_spectrum_rsd

  SUBROUTINE compute_power_spectrum_rsd2(d,L,kmin,kmax,bins,kpar,kper,pow,nbin,iz)

    USE constants
    USE fft
    IMPLICIT NONE
    INTEGER :: i, j, k, m, ii, jj, bins, iz
    REAL :: kx, ky, kz, kmod, L, kmin, kmax, a, b, kpers, kpars
    REAL :: pow(bins,bins), kpar(bins), kparbin(bins+1), kper(bins), kperbin(bins+1)
    DOUBLE PRECISION :: pow8(bins,bins)
    INTEGER :: nbin(bins,bins)
    INTEGER*8 :: nbin8(bins,bins)
    DOUBLE COMPLEX :: d(:,:,:)

    !    STOP 'Not updated according to JAP prescription'

    WRITE(*,*) 'Computing rsd power spectrum'

    kparbin=0.
    kperbin=0.
    kpar=0.
    kper=0.
    pow=0.
    nbin=0

    pow8=0.d0
    nbin8=0

    WRITE(*,*) 'Binning power'
    WRITE(*,*) 'Bins:', bins
    WRITE(*,*) 'k_min:', kmin
    WRITE(*,*) 'k_max:', kmax

    a=kmin
    b=kmax

    a=log10(a)
    b=log10(b)

    DO i=1,bins+1
       kparbin(i)=a+(b-a)*float(i-1)/float(bins)
    END DO

    DO i=1,bins
       kpar(i)=(kparbin(i)+kparbin(i+1))/2.
    END DO

    kparbin=10.**kparbin
    kpar=10.**kpar

    !Explicitly extend the bins to be sure to include all modes
    !This is necessary due to rounding errors!
    kparbin(1)=kparbin(1)*0.999
    kparbin(bins+1)=kparbin(bins+1)*1.001

    kperbin=kparbin
    kper=kpar

    m=SIZE(d(:,1,1))

    WRITE(*,*) 'Mesh:', m

    DO k=1,m
       DO j=1,m
          DO i=1,m

             IF(i==1 .AND. j==1 .AND. k==1) CYCLE

             CALL k_fft(i,j,k,m,kx,ky,kz,kmod,L)

             IF(iz==1) THEN
                kpars=ABS(kx)
                kpers=sqrt(ky**2.+kz**2.)
             ELSE IF(iz==2) THEN
                kpars=ABS(ky)
                kpers=sqrt(kz**2.+kx**2.)
             ELSE IF(iz==3) THEN
                kpars=ABS(kz)
                kpers=sqrt(kx**2.+ky**2.)
             ELSE
                STOP 'COMPUTE_POWER_SPECTRUM_RSD2: Error, iz not specified correctly'
             END IF

             DO jj=1,bins
                IF(kpars>=kparbin(jj) .AND. kpars<=kparbin(jj+1)) THEN                
                   DO ii=1,bins
                      IF(kpers>=kperbin(ii) .AND. kpers<=kperbin(ii+1)) THEN
                         pow8(ii,jj)=pow8(ii,jj)+ABS(d(i,j,k))**2.
                         nbin8(ii,jj)=nbin8(ii,jj)+1
                         EXIT
                      END IF
                   END DO
                   EXIT
                END IF
             END DO

          END DO
       END DO
    END DO

    !    DO jj=1,bins
    !       DO ii=1,bins
    !          pow(ii,jj)=pow(ii,jj)/log(kperbin(ii+1)/kperbin(ii))
    !          pow(ii,jj)=pow(ii,jj)/log(kparbin(jj+1)/kparbin(jj))
    !       END DO
    !    END DO

    DO jj=1,bins
       DO ii=1,bins
          IF(nbin8(ii,jj)==0) THEN
             pow8(ii,jj)=0.
          ELSE
             pow8(ii,jj)=pow8(ii,jj)/float(nbin8(ii,jj))
             pow8(ii,jj)=pow8(ii,jj)*(L**3.*kpar(jj)*kper(ii)**2.)/(2.*pi**2.)
          END IF
       END DO
    END DO

    pow=REAL(pow8)/(REAL(m)**6)

    !Divide by 2 because double count Hermitian conjugates
    nbin=INT(nbin8)/2

    WRITE(*,*) 'Power computed'
    WRITE(*,*) 

  END SUBROUTINE compute_power_spectrum_rsd2

  FUNCTION box_mode_power(dk,m)

    IMPLICIT NONE
    REAL :: box_mode_power
    DOUBLE COMPLEX, INTENT(IN) :: dk(m,m,m)
    INTEGER, INTENT(IN) :: m    

    box_mode_power=REAL(ABS(dk(2,1,1))**2.+ABS(dk(1,2,1))**2.+ABS(dk(1,1,2))**2.)/3.

  END FUNCTION box_mode_power

END MODULE field_operations
