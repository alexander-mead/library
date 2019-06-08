PROGRAM HMx_test

  USE cosmology_functions
  USE HMx
  IMPLICIT NONE
  CHARACTER(len=256) :: cosmo, halomodel
  INTEGER :: icosmo, ihm
  TYPE(cosmology) :: cosm
  TYPE(halomod) :: hmod
  LOGICAL, PARAMETER :: verbose=.TRUE.
  REAL, PARAMETER :: z=0.0

  CALL get_command_argument(1,cosmo)
  IF(cosmo=='') THEN
     icosmo=-1
  ELSE
     READ(cosmo,*) icosmo
  END IF

  CALL get_command_argument(2,halomodel)
  IF(halomodel=='') THEN
     ihm=-1
  ELSE
     READ(halomodel,*) ihm
  END IF
  
  ! Assigns the cosmological model
  CALL assign_cosmology(icosmo,cosm,verbose)
  CALL init_cosmology(cosm)
  CALL print_cosmology(cosm)

  ! Initiliasation for the halomodel calcualtion
  CALL assign_halomod(ihm,hmod,verbose)
  CALL init_halomod(mmin_HMx,mmax_HMx,scale_factor_z(z),hmod,cosm,verbose)
  CALL print_halomod(hmod,cosm,verbose)
  
END PROGRAM HMx_test

