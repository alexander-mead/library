PROGRAM multidark_demo

  USE multidark_stuff
  IMPLICIT NONE
  CHARACTER(len=256) :: infile
  REAL, ALLOCATABLE :: x(:,:), m(:)
  INTEGER :: i, n
  REAL, PARAMETER :: mmin=1e12 ! Minimum halo mass to read in

  ! Initial white space
  WRITE(*,*)

  ! Halo catalogue to read in
  !infile='/Volumes/Storage/Multidark/hlist_1.00109.list'
  infile='/Users/Mead/Physics/data/Multidark/hlist_1.00109.list'

  ! Read in the haloes
  CALL read_multidark_haloes(infile,mmin,x,m,n)

  OPEN(7,file='haloes_1.00109.cat')
  DO i=1,n
     WRITE(7,*) m(i), x(1,i), x(2,i), x(3,i)
  END DO
  CLOSE(7)
  
END PROGRAM multidark_demo
