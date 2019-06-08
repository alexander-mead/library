PROGRAM orbits_test

  USE constants
  USE orbits
  
  IMPLICIT NONE
  REAL :: M
  REAL :: a1, e1, w1, n1
  REAL :: a2, e2, w2, n2
  REAL :: x(3), v(3)

  ! Initial white space
  WRITE(*,*)

  ! Star parameters
  M=1.00 ! Central mass

  ! Planet parameters
  a1=1.20 ! Semi-major axis [au]
  e1=0.80 ! Eccentricity
  n1=pi*1.5 ! True anomaly [rad]
  w1=pi*1.5 ! Argument of periapsis [rad]

  ! Write the input orbit to screen
  WRITE(*,*) 'Input orbit:'
  WRITE(*,*) 'Semi-major axis [au]:', a1
  WRITE(*,*) 'Eccentricity:', e1
  WRITE(*,*) 'True anomaly [rad,deg]:', n1, n1*rad2deg
  WRITE(*,*) 'Argument of periapsis [rad,deg]:', w1, w1*rad2deg 
  WRITE(*,*)

  ! Calculate the orbit positions
  x=orbit_position(a1,e1,n1,w1)
  v=orbit_velocity(a1,e1,n1,w1,M)

  ! Write the orbit coordinates to screen
  WRITE(*,*) 'Orbit coordinates:'
  WRITE(*,*) 'x [au]:', x(1)
  WRITE(*,*) 'y [au]:', x(2)
  WRITE(*,*) 'z [au]:', x(3)
  WRITE(*,*) 'vx [2piau/yr]:', v(1)
  WRITE(*,*) 'vy [2piau/yr]:', v(2)
  WRITE(*,*) 'vz [2piau/yr]:', v(3)
  WRITE(*,*)

  ! Calculate orbital parameters again
  a2=semi_major_axis(x,v,M)
  e2=eccentricity(x,v,M)
  n2=true_anomaly(x,v,M)
  w2=argument_of_periapsis(x,v,M)  
  
  ! Write the output orbit to screen
  WRITE(*,*) 'Output orbit:'
  WRITE(*,*) 'Semi-major axis [au]:', a2
  WRITE(*,*) 'Eccentricity:', e2
  WRITE(*,*) 'True anomaly [rad,deg]:', n2, n2*rad2deg  
  WRITE(*,*) 'Argument of periapsis [rad,deg]:', w2, w2*rad2deg
  WRITE(*,*)

  ! Write ratios to screen to check parameters are recovered correctly
  WRITE(*,*) 'Ratio of input to output parameters:'
  WRITE(*,*) 'Semi-major axis:', a2/a1
  WRITE(*,*) 'Eccentricity:', e2/e1
  WRITE(*,*) 'True anomaly:', n2/n1
  WRITE(*,*) 'Argument of periapsis:', w2/w1
  WRITE(*,*)
  
END PROGRAM
