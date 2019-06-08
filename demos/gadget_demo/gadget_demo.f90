PROGRAM gadget_test

  USE gadget
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:,:), v(:,:)
  INTEGER, ALLOCATABLE :: id(:)
  REAL :: om_m, om_v, h, a, z, L, pm
  INTEGER :: n
  CHARACTER(len=256) :: infile, outfile
  LOGICAL :: lexist

  !!!FROM READ_WRITE_GADGET!!!

  CALL get_command_argument(1,infile)
  IF(infile=='') STOP 'Specify input Gadget file'
  INQUIRE(file=infile,exist=lexist)
  IF(lexist .EQV. .FALSE.) STOP 'This input file does not exist'

  CALL get_command_argument(2,outfile)
  IF(outfile=='') STOP 'Specify output Gadget file'
  
  CALL read_gadget(x,v,id,L,om_m,om_v,h,pm,a,z,n,infile)

  CALL write_gadget(x,v,id,L,om_m,om_v,h,pm,a,z,n,outfile)

  !!! !!!

!!$  !!!from READ_WRITE_CATALOGUE!!!
!!$
!!$  CALL get_command_argument(1,infile)
!!$  IF(infile=='') STOP 'Specify input catalogue file'
!!$  INQUIRE(file=infile,exist=lexist)
!!$  IF(lexist .EQV. .FALSE.) STOP 'This input catalogue file does not exist'
!!$
!!$  CALL get_command_argument(2,outfile)
!!$  IF(outfile=='') STOP 'Specify output catalogue file'
!!$
!!$  CALL read_catalogue(x,v,m,npart,disp,c,env,Dv,rmax,avg_r,rms_r,n,infile)
!!$
!!$  CALL write_catalogue(x,v,m,npart,disp,c,env,Dv,rmax,avg_r,rms_r,n,outfile)
!!$
!!$  !!! !!!

END PROGRAM gadget_test
