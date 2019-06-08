PROGRAM generate_randoms_tests

  USE constants
  USE simulations
  USE field_operations
  USE fft
  USE random_numbers
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x(:,:)
  INTEGER :: itest
  CHARACTER(len=256) :: outfile_power, outfile_slice
  
  ! PARAMETERS
  REAL, PARAMETER :: L=1. ! Box size
  INTEGER, PARAMETER :: nroot=64 ! Cube root of number of particles
  INTEGER, PARAMETER :: n=nroot**3 ! Number of particles
  INTEGER, PARAMETER :: m=256 ! Mesh size  
  INTEGER, PARAMETER :: nk=50 ! Number of bins for power spectrum
  INTEGER, PARAMETER :: iseed=1 ! Random number seed

  ! Initial white space
  WRITE(*,*)
  WRITE(*,*) 'GENERATE_RANDOMS_TEST: Generating random particle distributions'
  WRITE(*,*) 'GENERATE_RANDOMS_TEST: Box size:', L
  WRITE(*,*) 'GENERATE_RANDOMS_TEST: Number of particles:', n
  WRITE(*,*) 'GENERATE_RANDOMS_TEST: Mesh size for power spectrum:', m
  WRITE(*,*) 'GENERATE_RANDOMS_TEST: Number of k bins for power spectrum:', nk
  WRITE(*,*) 'GENERATE_RANDOMS_TEST: Random seed:', iseed
  WRITE(*,*)

  CALL RNG_set(iseed)

  ! Allocate position array
  ALLOCATE(x(3,n))

  DO itest=1,3

     ! Generate different particle loads
     IF(itest==1) THEN
        CALL generate_randoms(x,n,L)
        outfile_power='power_randoms.dat'
        outfile_slice='slice_randoms.dat'
     ELSE IF(itest==2) THEN
        CALL generate_grid(x,n,L)
        outfile_power='power_grid.dat'
        outfile_slice='slice_grid.dat'
     ELSE IF(itest==3) THEN
        CALL generate_poor_glass(x,n,L)
        outfile_power='power_poorglass.dat'
        outfile_slice='slice_poorglass.dat'
     ELSE
        STOP 'SIMULATIONS_TEST: Error, itest not specified correctly'
     END IF

     ! Write a slice file to look at
     CALL write_slice_ascii(x,n,0.,L/4.,0.,L/4.,0.,L/REAL(nroot),outfile_slice)

     ! Calculate the power spectrim
     CALL write_power_spectrum(x,n,L,m,nk,outfile_power)

  END DO
 
END PROGRAM generate_randoms_tests
