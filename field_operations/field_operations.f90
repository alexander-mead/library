MODULE field_operations
  
CONTAINS

  SUBROUTINE read_field(d,m,infile)

    !Read in a binary 'field' file
    USE array_operations
    USE statistics
    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: d(:,:,:)
    INTEGER, INTENT(IN) :: m

    ALLOCATE(d(m,m,m))

    !Output unformatted data
    WRITE(*,*) 'READ_FIELD: Binary input: ', TRIM(infile)
    WRITE(*,*) 'READ_FIELD: Mesh size:', m
    OPEN(7,file=infile,form='unformatted')
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
    REAL, ALLOCATABLE, INTENT(OUT) :: d(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: d8(:,:,:)
    INTEGER, INTENT(IN) :: m

    ALLOCATE(d(m,m,m),d8(m,m,m))

    !Input unformatted data
    WRITE(*,*) 'READ_FIELD: Binary input: ', TRIM(infile)
    WRITE(*,*) 'READ_FIELD: Mesh size:', m
    OPEN(7,file=infile,form='unformatted',access='stream')
    READ(7) d8
    CLOSE(7)
    d=REAL(d8)
    DEALLOCATE(d8)
    WRITE(*,*) 'READ_FIELD: Minval:', MINVAL(d)
    WRITE(*,*) 'READ_FIELD: Maxval:', MAXVAL(d)
    WRITE(*,*) 'READ_FIELD: Average:', REAL(mean(splay(d,m,m,m),m**3))
    WRITE(*,*) 'READ_FIELD: Variance:', REAL(variance(splay(d,m,m,m),m**3))
    WRITE(*,*) 'READ_FIELD: Done'
    WRITE(*,*)

  END SUBROUTINE read_field8

  SUBROUTINE write_field(d,m,outfile)

    !Write out a binary 'field' file
    IMPLICIT NONE    
    REAL, INTENT(IN) :: d(m,m,m)
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=256), INTENT(IN) :: outfile

    WRITE(*,*) 'WRITE_FIELD: Binary output: ', TRIM(outfile)
    WRITE(*,*) 'WRITE_FIELD: Mesh size:', m
    WRITE(*,*) 'WRITE_FIELD: Minval:', MINVAL(d)
    WRITE(*,*) 'WRITE_FIELD: Maxval:', MAXVAL(d)
    OPEN(7,file=outfile,form='unformatted')
    WRITE(7) d
    CLOSE(7)
    WRITE(*,*) 'WRITE_FIELD: Done'
    WRITE(*,*)

  END SUBROUTINE write_field

  SUBROUTINE print_2D_field(d,m,L,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d(m,m), L
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=256), INTENT(IN) :: outfile
    INTEGER :: i, j
    REAL :: x, y

    WRITE(*,*) 'PRINT_2D_FIELD: Writing to: ', TRIM(outfile)
    
    OPEN(8,file=outfile)
    DO j=1,m
       DO i=1,m

          x=L*(REAL(i)-0.5)/REAL(m)
          y=L*(REAL(j)-0.5)/REAL(m)

          !sum=0.
          !DO k=1,nz
          !   sum=sum+d(i,j,k)
          !END DO
          !sum=sum/float(nz)

          WRITE(8,*) x, y, d(i,j)

       END DO
    END DO
    CLOSE(8)

    WRITE(*,*) 'PRINT_2D_FIELD: Done'
    WRITE(*,*)

  END SUBROUTINE print_2D_field

  SUBROUTINE print_projected_field(d,m,L,nz,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d(m,m,m), L
    INTEGER, INTENT(IN) :: m, nz
    CHARACTER(len=256), INTENT(IN) :: outfile
    INTEGER :: i, j, k
    REAL :: x, y
    REAL :: sum

    WRITE(*,*) 'PRINT_PROJECTED_FIELD: Writing to: ', TRIM(outfile)
    WRITE(*,*) 'PRINT_PROJECTED_FIELD: Cells projecting:', nz
    
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

    WRITE(*,*) 'PRINT_PROJECTED_FIELD: Done'
    WRITE(*,*)
    
  END SUBROUTINE print_projected_field

  SUBROUTINE compress_field(d,ds,m)

    !Shrinks the field size by a factor of 2
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
    DOUBLE COMPLEX :: dc(m,m,m), dcout(m,m,m)

    !ibin = 1 NGP
    !ibin = 2 CIC

    WRITE(*,*) 'SHARPEN: Correcting for binning by sharpening field'
    WRITE(*,*) 'SHARPEN: Mesh size:', m

    dc=d

    CALL fft3(dc,dcout,m,m,m,-1)
    dc=dcout

    CALL sharpen_k(dc,m,L,ibin)

    CALL fft3(dc,dcout,m,m,m,1)
    dc=dcout

    d=REAL(REAL(dc))/REAL(m**3)

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

  SUBROUTINE smooth2D(arr,n,r,L)

    USE fft
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: arr(n,n)
    REAL, INTENT(IN) :: r, L
    REAL :: kx, ky, kz, kmod
    DOUBLE COMPLEX, ALLOCATABLE :: ac(:,:), acout(:,:)
    INTEGER :: i, j, m

    WRITE(*,*) 'Smoothing array'
    WRITE(*,*) 'Smoothing scale [Mpc/h]:', r

    !For padding, I cant imagine that x2 would ever be insufficient!
    m=2*n

    ALLOCATE(ac(m,m),acout(m,m))

    !Not sure if this is necessary
    ac=(0.d0,0.d0)
    acout=(0.d0,0.d0)

    !Put image into complex array, padded with 0s where image is not!
    DO j=1,n
       DO i=1,n
          ac(i,j)=arr(i,j)
       END DO
    END DO

    CALL fft2(ac,acout,n,n,-1)
    ac=acout

    !Smoothing length in terms of image(m x m) size!
    !r=pix/float(m)

    DO j=1,m
       DO i=1,m
          CALL k_fft(i,j,1,m,kx,ky,kz,kmod,L*2)
          ac(i,j)=ac(i,j)*exp(-((kmod*r)**2.)/2.)          
       END DO
    END DO

    CALL fft2(ac,acout,n,n,1)
    ac=acout

    !Normalise post Fourier transform!
    ac=ac/(REAL(m)**2)

    !Retrieve smooth image from complex array!
    DO j=1,n
       DO i=1,n
          arr(i,j)=REAL(REAL(ac(i,j)))
       END DO
    END DO

    WRITE(*,*) 'Done'
    WRITE(*,*)

  END SUBROUTINE smooth2D

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

    !Allocate complex array
    !ALLOCATE(ac(m,m,m))
    !ac=(0.d0,0.d0)

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
    REAL, ALLOCATABLE, INTENT(OUT) :: d2d(:,:)
    INTEGER :: i, j, k

    WRITE(*,*) 'PROJECT_3D_TO_2D: Projecting 3D stack into 2D'
    ALLOCATE(d2d(m,m))
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

  SUBROUTINE field_correlation_function(r_array,xi_array,n_array,n,d,m,L)

    USE gadget
    USE table_integer
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, m
    REAL, INTENT(OUT) :: xi_array(n)
    REAL, INTENT(IN) :: L, d(m,m,m), r_array(n)
    INTEGER*8, INTENT(OUT) :: n_array(n)
    REAL:: rmin, rmax
    DOUBLE PRECISION, ALLOCATABLE :: xi8_array(:)
    INTEGER :: i1, i2, i3, j1, j2, j3, i(3), j(3), k, dim
    REAL :: r, x1(3), x2(3)

    !This double counts, so time could be at least halved
    !Also could be parrallelised
    !Also could just not be complete shit, but it should get the job done

    rmin=r_array(1)
    rmax=r_array(n)

    WRITE(*,*) 'CORRELATION_FUNCTION: rmin [Mpc/h]:', rmin
    WRITE(*,*) 'CORRELATION_FUNCTION: rmax [Mpc/h]:', rmax
    WRITE(*,*) 'CORRELATION_FUNCTION: number of r bins:', n

    ALLOCATE(xi8_array(n))
    xi8_array=0.d0
    n_array=0

    DO i3=1,m
       DO i2=1,m
          DO i1=1,m

             i(1)=i1
             i(2)=i2
             i(3)=i3
             !x1(1)=L*(i1-0.5)/float(m)
             !x1(2)=L*(j1-0.5)/float(m)
             !x1(3)=L*(k1-0.5)/float(m)
             DO dim=1,3
                x1(dim)=L*(i(dim)-0.5)/float(m)
             END DO

             DO j3=1,m
                DO j2=1,m
                   DO j1=1,m

                      j(1)=j1
                      j(2)=j2
                      j(3)=j3
                      !x2(1)=L*(i2-0.5)/float(m)
                      !x2(2)=L*(j2-0.5)/float(m)
                      !x2(3)=L*(k2-0.5)/float(m)
                      DO dim=1,3
                         x2(dim)=L*(j(dim)-0.5)/float(m)
                      END DO

                      r=periodic_distance(x1,x2,L)

                      IF(r<rmin .OR. r>rmax) THEN
                         CYCLE
                      ELSE
                         k=select_table_integer(r,r_array,n,3)
                         IF(k<1 .OR. k>n) STOP 'Integer finding has fucked up'
                         xi8_array(k)=xi8_array(k)+d(i(1),i(2),i(3))*d(j(1),j(2),j(3))
                         n_array(k)=n_array(k)+1
                      END IF

                   END DO
                END DO
             END DO

          END DO
       END DO
    END DO

    xi_array=REAL(xi8_array/float(n_array))

    DEALLOCATE(xi8_array)

    WRITE(*,*) 'CORRELATION_FUNCTION: done'
    WRITE(*,*)

  END SUBROUTINE field_correlation_function

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

  SUBROUTINE pk(dk1,dk2,m,L,kmin,kmax,bins,k,pow,nbin)

    !Takes in a dk(m,m,m) array and computes the power spectrum
    USE table_integer
    USE constants
    USE array_operations
    USE fft
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(IN) :: dk1(m,m,m), dk2(m,m,m)
    REAL, ALLOCATABLE, INTENT(OUT) :: pow(:), k(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: nbin(:)
    INTEGER, INTENT(IN) :: m, bins
    REAL, INTENT(IN) :: L, kmin, kmax
    INTEGER :: i, ix, iy, iz, n
    REAL :: kx, ky, kz, kmod  
    REAL, ALLOCATABLE :: kbin(:)  
    REAL*8, ALLOCATABLE :: pow8(:), k8(:)    
    INTEGER*8, ALLOCATABLE :: nbin8(:)

    WRITE(*,*) 'PK: Computing isotropic power spectrum'

    !Allocate arrays used in this calculation
    ALLOCATE(kbin(bins+1),k(bins))
    ALLOCATE(pow(bins),nbin(bins))
    ALLOCATE(k8(bins),pow8(bins),nbin8(bins))

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
    nbin=INT(nbin8/2)

    !Deallocate arrays
    DEALLOCATE(kbin,pow8,nbin8,k8)

    WRITE(*,*) 'PK: Power computed'
    WRITE(*,*) 

  END SUBROUTINE pk

  FUNCTION box_mode_power(dk,m)

    IMPLICIT NONE
    REAL :: box_mode_power
    DOUBLE COMPLEX, INTENT(IN) :: dk(m,m,m)
    INTEGER, INTENT(IN) :: m    

    box_mode_power=REAL(ABS(dk(2,1,1))**2.+ABS(dk(1,2,1))**2.+ABS(dk(1,1,2))**2.)/3.

  END FUNCTION box_mode_power

END MODULE field_operations
