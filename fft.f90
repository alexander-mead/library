MODULE fft

  !This use statement seems not to be necessary any more
  !USE iso_c_binding

  !Include the FFTW libraray locations
  INCLUDE '/usr/local/include/fftw3.f' !Mac
  !INCLUDE '/usr/include/fftw3.f' !Linux'

CONTAINS

  SUBROUTINE k_fft(ix,iy,iz,m,kx,ky,kz,kmod,L)

    USE constants
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix, iy, iz, m
    REAL, INTENT(IN) :: L
    REAL, INTENT(OUT) :: kx, ky, kz, kmod
    !REAL, PARAMETER :: pi=3.141592654

    IF(MOD(m,2) .NE. 0) STOP 'K_FFT: Fourier transform does not have an even mesh'

    kx=float(ix-1)
    ky=float(iy-1)
    kz=float(iz-1)

    IF(ix>m/2+1) kx=-float(m-ix+1)
    IF(iy>m/2+1) ky=-float(m-iy+1)
    IF(iz>m/2+1) kz=-float(m-iz+1)

    kx=kx*2.*pi/L
    ky=ky*2.*pi/L
    kz=kz*2.*pi/L

    kmod=sqrt((kx**2)+(ky**2)+(kz**2))

  END SUBROUTINE k_fft

  SUBROUTINE FFT1(in,out,nx,ifb)

    !Wraps the 1D FFT
    IMPLICIT NONE
    !INCLUDE '/usr/local/include/fftw3.f' !Mac
    DOUBLE COMPLEX, INTENT(IN) :: in(nx)
    DOUBLE COMPLEX, INTENT(OUT) :: out(nx)
    INTEGER, INTENT(IN) :: ifb
    INTEGER, INTENT(IN) :: nx
    INTEGER*8 :: plan

    IF(ifb .NE. 1 .AND. ifb .NE. -1) THEN
       WRITE(*,*) 'FFT1: Error - need to specify forwards or backwards'
    END IF

    WRITE(*,*) 'FFT1: Starting FFT - creating plan'
    IF(ifb==-1) THEN
       call dfftw_plan_dft_1d(plan,nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
    ELSE IF(ifb==1) THEN
       call dfftw_plan_dft_1d(plan,nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
    END IF

    !This computes the FFT!
    WRITE(*,*) 'FFT1: Executing FFTW'
    call dfftw_execute(plan)
    WRITE(*,*) 'FFT1: FFTW complete'

    !And this destroys the plan!
    call dfftw_destroy_plan(plan)
    WRITE(*,*) 'FFT1: Plan destroyed'
    WRITE(*,*)

  END SUBROUTINE FFT1

  SUBROUTINE FFT2(in,out,nx,ny,ifb)

    !Wraps the 2D FFT
    !INCLUDE '/usr/local/include/fftw3.f' !Mac
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(IN) :: in(nx,ny)
    DOUBLE COMPLEX, INTENT(OUT) :: out(nx,ny)
    INTEGER, INTENT(IN) :: ifb
    INTEGER, INTENT(IN) :: nx, ny 
    INTEGER*8 :: plan

    IF(ifb .NE. 1 .AND. ifb .NE. -1) THEN
       WRITE(*,*) 'FFT2: Error - need to specify forwards or backwards'
    END IF

    !WRITE(*,*) 'FFT2: Starting FFT - creating plan'
    IF(ifb==-1) THEN
       CALL dfftw_plan_dft_2d(plan,nx,ny,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
    ELSE IF(ifb==1) THEN
       CALL dfftw_plan_dft_2d(plan,nx,ny,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
    END IF

    !This computes the FFT!
    !WRITE(*,*) 'FFT2: Executing FFTW'
    CALL dfftw_execute(plan)
    !WRITE(*,*) 'FFT2: FFTW complete'

    !And this destroys the plan!
    CALL dfftw_destroy_plan(plan)
    !WRITE(*,*) 'FFT2: Plan destroyed'
    !WRITE(*,*)

  END SUBROUTINE FFT2

  SUBROUTINE FFT3(in,out,nx,ny,nz,ifb)

    IMPLICIT NONE
    !INCLUDE '/usr/local/include/fftw3.f' !Mac
    DOUBLE COMPLEX, INTENT(IN)  :: in(nx,ny,nz)
    DOUBLE COMPLEX, INTENT(OUT) :: out(nx,ny,nz)
    INTEGER, INTENT(IN) :: nx, ny, nz
    INTEGER, INTENT(IN) :: ifb
    INTEGER*8 :: plan

    IF(ifb .NE. 1 .AND. ifb .NE. -1) THEN
       WRITE(*,*) 'Error - need to specify forwards or backwards'
    END IF

    !WRITE(*,*) 'Starting FFT - creating plan'
    IF(ifb==-1) THEN
       CALL dfftw_plan_dft_3d(plan,nx,ny,nz,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
    ELSE IF(ifb==1) THEN
       CALL dfftw_plan_dft_3d(plan,nx,ny,nz,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
    END IF

    !This computes the FFT!
    !WRITE(*,*) 'Executing FFTW'
    CALL dfftw_execute(plan)
    !WRITE(*,*) 'FFTW complete'

    !And this destroys the plan!
    CALL dfftw_destroy_plan(plan)
    !WRITE(*,*) 'Plan destroyed'
    !WRITE(*,*) ''

  END SUBROUTINE FFT3

END MODULE fft
