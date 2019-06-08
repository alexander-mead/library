PROGRAM field_operations_test

  USE field_operations

  IMPLICIT NONE
  CHARACTER(len=256) :: infile, outfile, mesh
  INTEGER :: m
  REAL :: L
  REAL, ALLOCATABLE :: d(:,:,:)
  LOGICAL :: lexist

  CALL get_command_argument(1,infile)
  IF(infile=='') STOP 'Specifify input file'
  INQUIRE(file=infile,exist=lexist)
  IF(lexist .EQV. .FALSE.) STOP 'This input file does not exist'  

  CALL get_command_argument(2,mesh)
  IF(mesh=='') STOP 'Specify mesh size'
  READ(mesh,*) m

  CALL get_command_argument(3,outfile)
  IF(outfile=='') STOP 'Specifify output file'

  CALL read_field(d,m,infile)
  L=1.
  CALL write_field_binary(d,m,L,outfile)
  
END PROGRAM field_operations_test
