PROGRAM file_info_test

  USE file_info
  
  IMPLICIT NONE
  CHARACTER(len=256) :: tod
  INTEGER :: n

  tod='cunt.dat'

  n=file_length(tod,verbose=.TRUE.)

CONTAINS

END PROGRAM file_info_test

