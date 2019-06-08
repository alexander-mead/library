PROGRAM dice_demo

  USE random_numbers
  USE logical_operations

  ! Calculate the score you get rolling N D3 with values 0,1,2
  IMPLICIT NONE
  INTEGER :: roll, nroll, dmin, dmax, ndice, i
  INTEGER, PARAMETER :: iseed=0
  CHARACTER(len=256), PARAMETER :: outfile='results.dat'

  ! Initial white space
  WRITE(*,*)

  ! Read in number of dice from command line
  CALL read_command_argument(1,dmin,'Specify the minimum number on the die (e.g., 1)')
  CALL read_command_argument(2,dmax,'Specify the maximum number on the die (e.g., 6)')
  CALL read_command_argument(3,ndice,'Specify the number of dice to roll (e.g., 2)')
  CALL read_command_argument(4,nroll,'Specify the number of rolls to make (e.g., 1)')

  IF(nroll<=0) STOP 'DICE_DEMO: Error, number of rolls must be greater than zero'

  ! Set the random number generator
  ! TODO: Should add verbose argument to RNG set to make it shut the fuck up
  CALL RNG_set(iseed,verbose=.FALSE.)

  IF(nroll>1) OPEN(7,file=outfile)

  ! Loop over dice rolls
  DO i=1,nroll

     ! Roll the dice and sum the score
     roll=dice(dmin,dmax,ndice)
     
     IF(nroll==1) THEN
        ! Write results to screen if just one roll
        WRITE(*,*) 'You rolled this many dice:', ndice
        WRITE(*,*) 'You rolled a score of:', roll
        WRITE(*,*)
     ELSE
        ! Otherwise write to file
        WRITE(7,*) roll
     END IF

  END DO

  IF(nroll>1) THEN
     WRITE(*,*) 'Data written to: ', trim(outfile)
     WRITE(*,*)
     CLOSE(7)
  END IF

END PROGRAM dice_demo
