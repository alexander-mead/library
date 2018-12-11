MODULE multidark_stuff

CONTAINS

  SUBROUTINE read_multidark_haloes(infile,x,m,n)

    USE file_info
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:,:)
    REAL, ALLOCATABLE, INTENT(OUT) :: m(:)
    INTEGER, INTENT(OUT) :: n
    LOGICAL :: lexist
    REAL :: c
    INTEGER :: i, p, pid
    REAL :: mm, xx, yy, zz
    INTEGER, PARAMETER :: hash_lines=58 ! Lines beginning with #
    REAL, PARAMETER :: mmin=1.75e12 ! Minimum halo mass [Msun/h]

    ! Check file exists
    INQUIRE(file=infile, exist=lexist)
    IF(.NOT. lexist) STOP 'READ_MULTIDARK_HALOES: Error, catalogue file does not exist'

    ! Welcome message
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Reading in halo catalogue'

    ! Find file length
    n=file_length(infile,verbose=.FALSE.)
    n=n-hash_lines

    ! Count unique haloes
    p=0 ! Set sum variable to zero
    WRITE(*,*) 'READ_MULTIDARK_HALOES: Total number of haloes:', n    
    OPEN(7,file=infile)
    DO i=1,hash_lines
       READ(7,*)
    END DO
    DO i=1,n
       READ(7,*) c, c, c, c, c, pid, c, c, c, c, mm
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
       ! Mass is 11, x,y,z are 18,19,20
       READ(7,*) c, c, c, c, c, pid, c, c, c, c, mm, c, c, c, c, c, c, xx, yy, zz
       IF(mm>mmin .AND. pid==-1) THEN
          p=p+1
          m(p)=mm
          x(1,p)=xx
          x(2,p)=yy
          x(3,p)=zz
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

  END SUBROUTINE
  
END MODULE multidark_stuff
