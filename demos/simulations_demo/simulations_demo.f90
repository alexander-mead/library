PROGRAM simulations_test

  USE random_numbers
  USE gadget
  USE simulations
  
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:,:), v(:,:), w(:)
  INTEGER, ALLOCATABLE :: id(:)
  REAL, ALLOCATABLE :: dNGP_2D(:,:), dCIC_2D(:,:)
  REAL, ALLOCATABLE :: dNGP_3D(:,:,:), dCIC_3D(:,:,:)
  REAL :: L, Om_m, Om_v, h, a, z
  REAL :: crap
  INTEGER :: n, m
  INTEGER :: itest
  INTEGER :: i, j, k
  INTEGER :: dim
  REAL :: dx(3)

  REAL, PARAMETER :: ddx=0.5
  REAL, PARAMETER :: ddy=0.25
  REAL, PARAMETER :: ddz=2.
  INTEGER, PARAMETER :: ntest=10
  INTEGER, PARAMETER :: iseed=0
  CHARACTER(len=256), PARAMETER :: infile='/Users/Mead/Physics/Gadget_test/N128/L256/Data_006'
  INTEGER, PARAMETER :: ingp=1
  INTEGER, PARAMETER :: icic=2
  LOGICAL, PARAMETER :: verbose=.TRUE.
  
  WRITE(*,*)
  WRITE(*,*) 'Simulations test'
  WRITE(*,*) '1 - Test replace'
  WRITE(*,*) '2 - Test random displacements'
  WRITE(*,*) '3 - Test random inversion'
  WRITE(*,*) '4 - Test random rotation'
  WRITE(*,*) '5 - Test 2D binning'
  WRITE(*,*) '6 - Test 3D binning'
  READ(*,*) itest
  WRITE(*,*)

  IF(itest==1 .OR. itest==2 .OR. itest==3 .OR. itest==4) THEN
     
     CALL read_gadget(x,v,id,L,Om_m,Om_v,h,crap,a,z,n,infile)

     ! Force some particles to be at boundary
     x(1,1)=0.
     x(1,2)=L

  END IF

  IF(itest==1) THEN

     WRITE(*,*) 'Original positions'
     CALL write_particles(x,ntest)

     ! Set displacements
     dx(1)=L*ddx
     dx(2)=L*ddy
     dx(3)=L*ddz

     WRITE(*,*) 'x displacement [Mpc/h]:', dx(1)
     WRITE(*,*) 'y displacement [Mpc/h]:', dx(2)
     WRITE(*,*) 'z displacement [Mpc/h]:', dx(3)
     WRITE(*,*)

     ! Displace particles
     DO i=1,3
        x(i,:)=x(i,:)+dx(i)
     END DO

     WRITE(*,*) 'Displaced positions'
     CALL write_particles(x,ntest)

     CALL replace(x,n,L,verbose)
     
     WRITE(*,*) 'Replaced positions'
     CALL write_particles(x,ntest)

  ELSE IF(itest==2 .OR. itest==3 .OR. itest==4) THEN

     CALL RNG_set(iseed)

     WRITE(*,*) 'Original positions'
     CALL write_particles(x,ntest)

     IF(itest==2) THEN
        CALL random_translation(x,n,L)
     ELSE IF(itest==3) THEN
        CALL random_inversion(x,n,L)
     ELSE IF(itest==4) THEN
        CALL random_rotation(x,n)
     ELSE
        STOP 'SIMULATIONS_TEST: Error, test not specified correctly'
     END IF

     WRITE(*,*) 'New positions'
     CALL write_particles(x,ntest)

  ELSE IF(itest==5 .OR. itest==6) THEN

     IF(itest==5) dim=2
     IF(itest==6) dim=3

     n=1
     IF(itest==5) m=3
     IF(itest==6) m=3
     
     ALLOCATE(x(dim,n),w(n))

     IF(itest==5) THEN
        ALLOCATE(dNGP_2D(m,m))
        ALLOCATE(dCIC_2D(m,m))
     ELSE IF(itest==6) THEN
        ALLOCATE(dNGP_3D(m,m,m))
        ALLOCATE(dCIC_3D(m,m,m))
     END IF

     L=1.
     x=0.33333/2. ! Place in centre of first cell
     x(1,1)=x(1,1)+0.01 ! Displace in x slightly
     w=1.

     WRITE(*,*) 'Particle position:', (x(j,1), j=1,dim)

     IF(itest==5) THEN
        CALL particle_bin(x,n,L,w,dNGP_2D,m,ingp,all=.TRUE.,periodic=.TRUE.,verbose=.TRUE.)
        CALL particle_bin(x,n,L,w,dCIC_2D,m,icic,all=.TRUE.,periodic=.TRUE.,verbose=.TRUE.)
     ELSE IF(itest==6) THEN
        CALL particle_bin(x,n,L,w,dNGP_3D,m,ingp,all=.TRUE.,periodic=.TRUE.,verbose=.TRUE.)
        CALL particle_bin(x,n,L,w,dCIC_3D,m,icic,all=.TRUE.,periodic=.TRUE.,verbose=.TRUE.)
     END IF

     IF(itest==5) THEN

        DO i=1,m
           DO j=1,m
              WRITE(*,*) i, j, dNGP_2D(i,j), dCIC_2D(i,j)
           END DO
        END DO      

     ELSE IF(itest==6) THEN       

        DO i=1,m
           DO j=1,m
              DO k=1,m
                 WRITE(*,*) i, j, k, dNGP_3D(i,j,k), dCIC_3D(i,j,k)
              END DO
           END DO
        END DO

     END IF
     
  ELSE

     STOP 'SIMULATIONS_TEST: Error, test not specified correctly'

  END IF

CONTAINS

  SUBROUTINE write_particles(x,n)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x(3,n)
    INTEGER, INTENT(IN) :: n

    DO i=1,n
       WRITE(*,*) i, x(1,i), x(2,i), x(3,i)
    END DO
    WRITE(*,*)

   END SUBROUTINE write_particles
  
END PROGRAM simulations_test
