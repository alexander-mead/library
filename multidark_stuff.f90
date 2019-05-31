MODULE multidark_stuff

  USE file_info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_multidark_haloes
  PUBLIC :: read_multidark_halo_catalogue
  PUBLIC :: read_multidark_particles
  PUBLIC :: write_multidark_halo_catalogue

CONTAINS

  SUBROUTINE read_multidark_haloes(infile,mmin,x,m,n)

    ! TODO: Make minimum mass (or minimum halo-particle number) an input
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile      ! File to read in
    REAL, INTENT(IN) :: mmin                    ! Minimum halo virial mass [Msun/h]
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:)    ! Position array [Mpc/h]
    REAL, ALLOCATABLE, INTENT(OUT) :: m(:,:)    ! Virial halo mass [Msun/h]
    INTEGER, INTENT(OUT) :: n                   ! Total number of haloes
    LOGICAL :: lexist
    INTEGER :: i, j
    INTEGER :: p, pid
    REAL :: mm!, c
    REAL, ALLOCATABLE :: data(:)
    
    INTEGER, PARAMETER :: hash_lines=58 ! Number of lines beginning with # (Both Multidark and Bolshoi have 58)
    !REAL, PARAMETER :: mmin=1.74e12     ! Minimum halo mass [Msun/h] (corresponds to N>200 for MDR1; 1e12 would mean n~115)
    !INTEGER, PARAMETER :: columns=73       ! Total number of columns in file
    INTEGER, PARAMETER :: columns=41       ! Total number of columns to read from file
    INTEGER, PARAMETER :: column_pid=6     ! Column for PID (-1 if unique halo)
    INTEGER, PARAMETER :: column_mv=11     ! Column for virial mass (unbinding done) [Msun/h]
    INTEGER, PARAMETER :: column_x=18      ! Column for x position [Mpc/h]
    INTEGER, PARAMETER :: column_y=19      ! Column for y position [Mpc/h]
    INTEGER, PARAMETER :: column_z=20      ! Column for z position [Mpc/h]
    INTEGER, PARAMETER :: column_vx=21     ! Column for vx position [km/s]
    INTEGER, PARAMETER :: column_vy=22     ! Column for vy position [km/s]
    INTEGER, PARAMETER :: column_vz=23     ! Column for vz position [km/s]
    INTEGER, PARAMETER :: column_mvu=37    ! Column for total virial mass (no unbinding) [Msun/h]
    INTEGER, PARAMETER :: column_m200=38   ! Column for M200 [Msun/h]
    INTEGER, PARAMETER :: column_m200c=39  ! Column for M200 critical [Msun/h]
    INTEGER, PARAMETER :: column_m500c=40  ! Column for M500 critical [Msun/h]
    INTEGER, PARAMETER :: column_m2500c=41 ! Column for M2500 critical [Msun/h]

    ! Check file exists
    INQUIRE(file=infile, exist=lexist)
    IF(.NOT. lexist) THEN
       WRITE(*,*) 'READ_MULTIDARK_HALOES: Error, catalogue file does not exist: ', trim(infile)
       STOP
    END IF
       
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
       mm=data(column_mv) ! Read virial mass
       pid=nint(data(column_pid))
       IF(mm>mmin .AND. pid==-1) p=p+1
    END DO
    CLOSE(7)
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Minimum allowed halo mass [Msun/h]:', mmin
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Total number of distinct haloes:', p   
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Fraction of distinct haloes:', REAL(p)/REAL(n)

    ALLOCATE(x(3,p),m(6,p))
    x=0.
    m=0.

    ! Fill arrays with unique halo properties
    p=0 ! Reset sum variable to zero
    OPEN(7,file=infile)
    DO i=1,hash_lines
       READ(7,*)
    END DO
    DO i=1,n
       ! Virial mass is column 11, x,y,z positions are columns 18,19,20
       ! Halo is unique iff pid=-1 (pid is column 6)
       !READ(7,*) c, c, c, c, c, pid, c, c, c, c, mm, c, c, c, c, c, c, xx, yy, zz
       READ(7,*) (data(j), j=1,columns)
       mm=data(column_mv)
       pid=nint(data(column_pid))
       IF(mm>mmin .AND. pid==-1) THEN
          p=p+1
          !WRITE(*,*) p
          x(1,p)=data(column_x)
          x(2,p)=data(column_y)
          x(3,p)=data(column_z)
          !WRITE(*,*) x(1,p), x(2,p), x(3,p)
          m(1,p)=data(column_mv)
          m(2,p)=data(column_mvu)
          m(3,p)=data(column_m200)
          m(4,p)=data(column_m200c)
          m(5,p)=data(column_m500c)
          m(6,p)=data(column_m2500c)
          !WRITE(*,*) m(1,p)
          !WRITE(*,*)
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
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Minimum virial halo mass [Msun/h]:', MINVAL(m(1,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Maximum virial halo mass [Msun/h]:', MAXVAL(m(1,:))
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Done'
    WRITE(*,*)

  END SUBROUTINE read_multidark_haloes

  SUBROUTINE read_multidark_halo_catalogue(infile,x,m,n)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:)
    REAL, ALLOCATABLE, INTENT(OUT) :: m(:,:)   
    INTEGER, INTENT(OUT) :: n
    INTEGER :: i

    n=file_length(infile,verbose=.FALSE.)
    ALLOCATE(x(3,n),m(6,n))

    WRITE(*,*) 'READ_MULTIDARK_HALO_CATALOGUE: ', trim(infile)
    WRITE(*,*) 'READ_MULTIDARK_HALO_CATALOGUE: Number of haloes:', n
    OPEN(7,file=infile,status='old')
    DO i=1,n
       READ(7,*) x(1,i), x(2,i), x(3,i), m(1,i), m(2,i), m(3,i), m(4,i), m(5,i), m(6,i)
    END DO
    CLOSE(7)
    WRITE(*,*) 'READ_MULTIDARK_HALO_CATALOGUE: Done'
    WRITE(*,*)
    
  END SUBROUTINE read_multidark_halo_catalogue

  SUBROUTINE write_multidark_halo_catalogue(outfile,x,m,idx,n)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: outfile
    REAL, INTENT(IN) :: x(3,n)
    REAL, INTENT(IN) :: m(6,n)
    INTEGER, INTENT(IN) :: idx(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j

    ! Write out the little catalogue
    WRITE(*,*) 'WRITE_MULTIDARK_HALO_CATALOGUE: Writing outfile: ', trim(outfile)
    OPEN(7,file=outfile)
    DO i=1,n
       j=idx(n+1-i)
       WRITE(7,*) x(1,j), x(2,j), x(3,j), m(1,j), m(2,j), m(3,j), m(4,j), m(5,j), m(6,j)
    END DO
    CLOSE(7)
    WRITE(*,*) 'WRITE_MULTIDARK_HALO_CATALOGUE: Done'
    WRITE(*,*)

  END SUBROUTINE write_multidark_halo_catalogue

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
