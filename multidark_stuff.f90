MODULE multidark_stuff

  USE file_info

CONTAINS

  SUBROUTINE read_multidark_haloes(infile,mmin,x,m,n)

    ! TODO: Make minimum mass (or minimum halo-particle number) an input
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, INTENT(IN) :: mmin
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:)
    REAL, ALLOCATABLE, INTENT(OUT) :: m(:)
    INTEGER, INTENT(OUT) :: n
    LOGICAL :: lexist
    !REAL :: c
    INTEGER :: i, j
    INTEGER :: p, pid
    REAL :: mm!, xx, yy, zz
    REAL, ALLOCATABLE :: data(:)
    
    INTEGER, PARAMETER :: hash_lines=58 ! Number of lines beginning with # (Both Multidark and Bolshoi have 58)
    !REAL, PARAMETER :: mmin=1.74e12     ! Minimum halo mass [Msun/h] (corresponds to N>200 for MDR1; 1e12 would mean n~115)
    INTEGER, PARAMETER :: columns=73       ! Total number of columns in file
    INTEGER, PARAMETER :: column_pid=6     ! Column for PID (-1 if unique halo)
    INTEGER, PARAMETER :: column_mv=11     ! Column for virial mass (unbinding done) [Msun/h]
    INTEGER, PARAMETER :: column_x=18      ! Column for x position [Mpc/h]
    INTEGER, PARAMETER :: column_y=19      ! Column for y position [Mpc/h]
    INTEGER, PARAMETER :: column_z=20      ! Column for z position [Mpc/h]
    INTEGER, PARAMETER :: column_vx=21     ! Column for vx position [km/s]
    INTEGER, PARAMETER :: column_vy=22     ! Column for vy position [km/s]
    INTEGER, PARAMETER :: column_vz=23     ! Column for vz position [km/s]
    INTEGER, PARAMETER :: column_mvu=37    ! Column for total virial mass (no particle unbinding; for distinct haloes this difference is only 1-2%; Msun/h)
    INTEGER, PARAMETER :: column_m200=38   ! Column for M200 [Msun/h]
    INTEGER, PARAMETER :: column_m200c=39  ! Column for M200 critical [Msun/h]
    INTEGER, PARAMETER :: column_m500c=40  ! Column for M500 critical [Msun/h]
    INTEGER, PARAMETER :: column_m2500c=41 ! Column for M2500 critical [Msun/h]

    ! Check file exists
    INQUIRE(file=infile, exist=lexist)
    IF(.NOT. lexist) STOP 'READ_MULTIDARK_HALOES: Error, catalogue file does not exist'

    ! Welcome message
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Reading in halo catalogue: ', trim(infile)

    ! Find file length
    n=file_length(infile,verbose=.FALSE.)
    n=n-hash_lines

    ! Allocate the data array with the total number of columns so that you can read in a full line
    ALLOCATE(data(columns))

    ! Count unique haloes
    p=0 ! Set sum variable to zero
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Total number of haloes:', n    
    OPEN(7,file=infile)
    DO i=1,hash_lines
       READ(7,*)
    END DO
    DO i=1,n
       !READ(7,*) c, c, c, c, c, pid, c, c, c, c, mm
       READ(7,*) (data(j), j=1,columns)
       mm=data(column_mvu)
       pid=NINT(data(column_pid))
       IF(mm>mmin .AND. pid==-1) p=p+1
    END DO
    CLOSE(7)
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Minimum allowed halo mass [Msun/h]:', mmin
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Total number of distinct haloes:', p   
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Fraction of distinct haloes:', REAL(p)/REAL(n)

    ALLOCATE(x(3,p),m(p))

    ! Fill arrays with unique halo properties
    p=0
    OPEN(7,file=infile)
    DO i=1,hash_lines
       READ(7,*)
    END DO
    DO i=1,n
       ! Virial mass is column 11, x,y,z positions are columns 18,19,20
       ! Halo is unique iff pid=-1 (pid is column 6)
       !READ(7,*) c, c, c, c, c, pid, c, c, c, c, mm, c, c, c, c, c, c, xx, yy, zz
       READ(7,*) (data(j), j=1,columns)
       mm=data(column_mvu)
       pid=NINT(data(column_pid))
       IF(mm>=mmin .AND. pid==-1) THEN
          p=p+1
          !m(p)=mm
          !x(1,p)=xx
          !x(2,p)=yy
          !x(3,p)=zz
          m(p)=mm
          x(1,p)=data(column_x)
          x(2,p)=data(column_y)
          x(3,p)=data(column_z)
       END IF
    END DO
    CLOSE(7)

    ! Set the total output number
    n=p

    ! Write min/max values to screen
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Minimum x [Mpc/h]:', MINVAL(x(1,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Maximum x [Mpc/h]:', MAXVAL(x(1,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Minimum y [Mpc/h]:', MINVAL(x(2,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Maximum y [Mpc/h]:', MAXVAL(x(2,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Minimum z [Mpc/h]:', MINVAL(x(3,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Maximum z [Mpc/h]:', MAXVAL(x(3,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Minimum halo mass [Msun/h]:', MINVAL(m)
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Maximum halo mass [Msun/h]:', MAXVAL(m)
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Done'
    WRITE(*,*)

  END SUBROUTINE read_multidark_haloes

  SUBROUTINE read_multidark_halo_catalogue(infile,m,x,n)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: m(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:)
    INTEGER, INTENT(OUT) :: n
    INTEGER :: i

    n=file_length(infile,verbose=.FALSE.)
    ALLOCATE(m(n),x(3,n))

    WRITE(*,*) 'READ_MULTIDARK_HALO_CATALOGUE: ', trim(infile)
    WRITE(*,*) 'READ_MULTIDARK_HALO_CATALOGUE: Number of haloes:', n
    OPEN(7,file=infile,status='old')
    DO i=1,n
       READ(7,*) m(i), x(1,i), x(2,i), x(3,i)
    END DO
    CLOSE(7)
    WRITE(*,*) 'READ_MULTIDARK_HALO_CATALOGUE: Done'
    WRITE(*,*)
    
  END SUBROUTINE read_multidark_halo_catalogue

  SUBROUTINE read_multidark_particles(infile,x,n)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:)
    INTEGER, INTENT(OUT) :: n
    INTEGER :: i
    REAL :: crap

    WRITE(*,*) 'READ_MULTIDARK_PARTICLES: ', trim(infile)

    n=file_length(infile,verbose=.FALSE.)
    n=n-1 ! First line is comment
    ALLOCATE(x(3,n))
   
    WRITE(*,*) 'READ_MULTIDARK_PARTICLES: Number of particles:', n
    OPEN(7,file=infile,status='old')
    READ(7,*) ! First line is comment
    DO i=1,n
       READ(7,*) crap, x(1,i), x(2,i), x(3,i)
    END DO
    CLOSE(7)
    WRITE(*,*) 'READ_MULTIDARK_PARTICLES: Done'
    WRITE(*,*)
    
  END SUBROUTINE read_multidark_particles
  
END MODULE multidark_stuff
