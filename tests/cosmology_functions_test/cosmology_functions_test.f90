PROGRAM cosmology_functions_test

  USE constants
  USE array_operations
  USE basic_operations
  USE cosmology_functions

  IMPLICIT NONE
  LOGICAL :: fail_main
  LOGICAL, PARAMETER :: verbose_main = .FALSE.

  WRITE(*,*)

  CALL test_EdS(fail_main, verbose_main)

CONTAINS

  SUBROUTINE test_EdS(fail, verbose)

    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: fail
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: icosmo, i
    TYPE(cosmology) :: cosm
    REAL :: a, p, t, g, h
    REAL :: p_test, t_test, g_test, h_test
    REAL, ALLOCATABLE :: as(:)
    
    REAL, PARAMETER :: amin=1e-3
    REAL, PARAMETER :: amax=1e0
    INTEGER, PARAMETER :: na=64
    REAL, PARAMETER :: eps=1e-3

    fail=.FALSE.

    ! 15 - TCDM model
    icosmo=15 

    ! Assign the cosmological model
    CALL assign_cosmology(icosmo,cosm,verbose)
    CALL init_cosmology(cosm)
    CALL print_cosmology(cosm)

    CALL fill_array(log(amin),log(amax),as,na)
    as=exp(as)
    as(na)=amax ! To stop going above 1.

    DO i=1,na

       ! Scale factor
       a=as(i)

       ! Values from cosmology routines
       g=grow(a,cosm)
       h=acc_growth(a,cosm)
       p=comoving_particle_horizon(a,cosm)
       t=cosmic_time(a,cosm)       
       
       ! Tests
       g_test=a ! Growth in EdS
       h_test=a ! Accumulated growth in EdS
       p_test=p_EdS(a,cosm%Om_m,cosm%Om_r) ! Particle horizon in EdS
       t_test=t_EdS(a,cosm%Om_m,cosm%Om_r) ! Cosmic time in EdS
       
       !IF(verbose) WRITE(*,*) i, a, r/r_test, t/t_test, g/g_test, h/h_test
       IF(verbose)THEN
          WRITE(*,*) 'Scale factor:', a
          WRITE(*,*) 'Growth:', g, g_test, g/g_test
          WRITE(*,*) 'Accumulated growth:', h, h_test, h/h_test
          WRITE(*,*) 'Particle horizon [Mpc/h]:', p, p_test, p/p_test
          WRITE(*,*) 'Cosmic time [Gyr/h]:', t, t_test, t/t_test          
          WRITE(*,*)
       END IF
       
       IF(.NOT. requal(g,g_test,eps)) fail=.TRUE.
       IF(.NOT. requal(p,p_test,eps)) fail=.TRUE.
       IF(.NOT. requal(t,t_test,eps)) fail=.TRUE.
       IF(.NOT. requal(h,h_test,eps)) fail=.TRUE.
       
    END DO

    IF(fail) THEN
       STOP 'TEST_EDS: Fail'
    ELSE
       WRITE(*,*) 'TEST_EDS: Pass'
       WRITE(*,*)
    END IF

  END SUBROUTINE test_EdS

  REAL FUNCTION p_EdS(a,Om_m,Om_r)

    ! Comoving particle horizon at scale factor 'a' in Einstein-de Sitter (with radiation) [Mpc/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    REAL, INTENT(IN) :: Om_m
    REAL, INTENT(IN) :: Om_r

    p_EdS=2.*Hdist*(sqrt(Om_m*a+Om_r)/Om_m-sqrt(Om_r)/Om_m)

  END FUNCTION p_EdS

  REAL FUNCTION t_EdS(a,Om_m,Om_r)

    ! Comoving time at scale factor 'a' in Einstein-de Sitter (with radiation) [Gyr/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    REAL, INTENT(IN) :: Om_m
    REAL, INTENT(IN) :: Om_r

    t_EdS=2.*Htime*((Om_m*a-2.*Om_r)*sqrt(Om_m*a+Om_r)+2.*Om_r**(3./2.))/(3.*Om_m**2)

  END FUNCTION t_EdS

  REAL FUNCTION t_LCDM(a,Om_m,Om_v)

    ! Comoving time at scale factor 'a' in LCDM (no radiation) [Gyr/h]
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    REAL, INTENT(IN) :: Om_m
    REAL, INTENT(IN) :: Om_v

    t_LCDM=Htime*(2./(3.*sqrt(Om_v)))*asinh(sqrt(Om_v*a**3/Om_m))

  END FUNCTION t_LCDM

END PROGRAM
