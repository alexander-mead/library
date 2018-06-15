MODULE cosmology_functions

  USE interpolate
  USE constants

  IMPLICIT NONE

  !Contains cosmological parameters that only need to be calculated once
  TYPE cosmology     
     REAL :: Om_m, Om_b, Om_v, Om_w, Om_nu, h, n, sig8, w, wa, m_wdm, YHe !Primary parameters
     REAL :: z_CMB, T_CMB, neff, Om_r, age, horizon !Secondary parameters
     REAL :: Om, k, Om_k, Om_c, A !Derived parameters
     REAL :: YH, mue, mup!, epfac !Derived thermal parameters
     REAL :: a1, a2, ns, ws, am, dm, wm !DE parameters     
     REAL :: Om_ws, as, a1n, a2n !Derived DE parameters
     REAL :: alpha, eps, Gamma, M0, Astar, whim !Baryon parameters
     REAL :: mgal, HImin, HImax !HOD parameters
     INTEGER :: iw !Switches
     REAL, ALLOCATABLE :: log_sigma(:), log_r_sigma(:) !Arrays for sigma(R)
     REAL, ALLOCATABLE :: log_a_growth(:), log_growth(:), growth_rate(:), log_acc_growth(:) !Arrays for growth
     REAL, ALLOCATABLE :: r(:), a_r(:) !Arrays for distance
     REAL, ALLOCATABLE :: log_plin(:), log_k_plin(:) !Arrays for input linear P(k)
     REAL, ALLOCATABLE :: log_a_dcDv(:), dc(:), Dv(:) !Arrays for spherical-collapse parameters
     INTEGER :: n_sigma, n_growth, n_r, n_plin, n_dcDv !Array entries
     REAL :: gnorm
     CHARACTER(len=256) :: name = ""
     LOGICAL :: wdm
     LOGICAL :: has_distance, has_growth, has_sigma, has_spherical, has_power
     LOGICAL :: is_init, is_normalised
     LOGICAL :: external_plin     
  END TYPE cosmology

  !Global parameters
  REAL, PARAMETER :: acc_cosm=1e-4 !Accuacy for the integrations
  INTEGER, PARAMETER :: ibox=0 !Consider the simulation volume
  REAL, PARAMETER :: Lbox=400. !Simulation box size
  LOGICAL, PARAMETER :: verbose_cosmology=.TRUE. !Global verbose parameter

CONTAINS

  SUBROUTINE assign_cosmology(icosmo,cosm)

    !Assigns the cosmological parameters
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(INOUT) :: icosmo    
    INTEGER :: i
    REAL :: Om_c, Om_g_h2, rho_g

    !Names of pre-defined cosmologies    
    INTEGER, PARAMETER :: ncosmo=14
    CHARACTER(len=256) :: names(0:ncosmo)
    names(0)='User defined'
    names(1)='Boring'
    names(2)='WMAP7 (cosmo-OWLS version; 1312.5462)'
    names(3)='Planck 2013 (cosmo-OWLS version; 1312.5462)'
    names(4)='WMAP9 (BAHAMAS version: 1712.02411)'
    names(5)='Boring open model'
    names(6)='Einstein de-Sitter'
    names(7)='IDE I (user)'
    names(8)='IDE II (user)'
    names(9)='IDE III (user)'
    names(10)='IDE3'
    names(11)='IDE10'
    names(12)='LCDM (user)'
    names(13)='w(a)CDM (user)'
    names(14)='Boring WDM'

    IF(icosmo==-1) THEN
       WRITE(*,*) 'ASSIGN_COSMOLOGY: Choose cosmological model'
       WRITE(*,*) '==========================================='
       DO i=0,SIZE(names)-1
          WRITE(*,*) i, '- ', TRIM(names(i))
       END DO
       READ(*,*) icosmo
       WRITE(*,*) '==========================================='
    END IF

    !Set the name of the cosmological model
    cosm%name=names(icosmo)

    !Boring default cosmology
    cosm%Om_m=0.3
    cosm%Om_b=0.05
    cosm%Om_v=1.-cosm%Om_m
    cosm%Om_w=0.
    cosm%Om_nu=0.
    cosm%h=0.7
    cosm%sig8=0.8
    cosm%n=0.96
    cosm%w=-1.
    cosm%wa=0.
    cosm%T_CMB=2.725 !CMB temperature [K]
    cosm%z_CMB=1100. !Redshift of the last-scatting surface
    cosm%neff=3.046 !Effective number of relativistic neutrinos
    cosm%YHe=0.24 !Helium mass fraction

    !Default dark energy is Lambda
    cosm%iw=1

    !Default to have no WDM
    cosm%wdm=.FALSE.
    cosm%m_wdm=1. !WDM mass [keV]   

    !Initially set the normalisation to 1
    cosm%A=1.

    !Default values of baryon parameters
    !TODO: These should eventually be HMx parameters
    cosm%alpha=1.
    cosm%eps=1.
    cosm%Gamma=1.17
    cosm%M0=1e14 !Halo mass that has lost half gas
    cosm%Astar=0.02 !Maximum star-formation efficiency
    cosm%whim=1e6 !WHIM temperature [K]

    !Default values of the HOD parameters
    !TODO: These should eventually be HMx parameters
    cosm%mgal=1e13

    !Default values for the HI parameters
    !TODO: These should eventually be HMx parameters
    cosm%HImin=1e9
    cosm%HImax=1e12

    !Set all 'has' logicals to FALSE
    cosm%is_init=.FALSE.
    cosm%has_distance=.FALSE.
    cosm%has_growth=.FALSE.
    cosm%has_sigma=.FALSE.
    cosm%has_power=.FALSE.
    cosm%has_spherical=.FALSE.
    cosm%is_normalised=.FALSE.
    cosm%external_plin=.FALSE.    

    IF(icosmo==0) THEN
       STOP 'Need to implement user decision here'
    ELSE IF(icosmo==1) THEN
       !Boring - do nothing
    ELSE IF(icosmo==2) THEN
       !cosmo-OWLS - WMAP7 (1312.5462)
       cosm%Om_m=0.272
       cosm%Om_b=0.0455
       cosm%Om_v=1.-cosm%Om_m
       cosm%Om_nu=0.
       cosm%h=0.704
       cosm%sig8=0.81
       cosm%n=0.967
    ELSE IF(icosmo==3) THEN
       !cosmo-OWLS - Planck 2013 (1312.5462)
       cosm%Om_m=0.3175
       cosm%Om_b=0.0490
       cosm%Om_v=1.-cosm%Om_m
       cosm%h=0.6711
       cosm%n=0.9624
       cosm%sig8=0.834
    ELSE IF(icosmo==4) THEN
       !BAHAMAS - WMAP9 (1712.02411)
       cosm%h=0.7
       Om_c=0.2330
       cosm%Om_b=0.0463
       cosm%Om_m=Om_c+cosm%Om_b
       cosm%Om_v=1.-cosm%Om_m       
       cosm%Om_nu=0.
       cosm%n=0.9720
       cosm%sig8=0.8211
    ELSE IF(icosmo==5) THEN
       !Boring open model
       cosm%Om_v=0.
    ELSE IF(icosmo==6) THEN
       !Einstein-de Sitter
       cosm%Om_m=1.
       cosm%Om_v=0.
    ELSE IF(icosmo==7) THEN
       !IDE I
       cosm%iw=5
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%as
       WRITE(*,*) 'n*:'
       READ(*,*) cosm%ns
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       cosm%Om_m=0.3
       cosm%Om_v=0.7
    ELSE IF(icosmo==8) THEN
       !IDE II model
       cosm%iw=6      
       WRITE(*,*) 'n*:'
       READ(*,*) cosm%ns
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%as
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       cosm%Om_m=0.3
       cosm%Om_w=0.7
       cosm%Om_v=0. !No vacuum necessary here
    ELSE IF(icosmo==9) THEN
       !IDE III model
       cosm%iw=7
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%as
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       WRITE(*,*) 'w*:'
       READ(*,*) cosm%ws
       cosm%Om_m=0.3
       cosm%Om_w=0.7
       cosm%Om_v=0.
    ELSE IF(icosmo==10 .OR. icosmo==11) THEN
       cosm%iw=6
       cosm%Om_w=cosm%Om_v
       cosm%Om_v=0.
       IF(icosmo==10) cosm%ns=3
       IF(icosmo==11) cosm%ns=10
       cosm%as=0.01
       cosm%Om_ws=0.2
    ELSE IF(icosmo==12) THEN
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Om_v:'
       READ(*,*) cosm%Om_v
    ELSE IF(icosmo==13) THEN
       cosm%iw=3
       WRITE(*,*) 'w0:'
       READ(*,*) cosm%w
       WRITE(*,*) 'wa:'
       READ(*,*) cosm%wa
    ELSE IF(icosmo==14) THEN
       cosm%wdm=.TRUE.
    ELSE
       STOP 'ASSIGN_COSMOLOGY: Error, icosmo not specified correctly'
    END IF

    !Radiation density
    rho_g=(4.*SBconst*cosm%T_CMB**4/c_light**3)
    Om_g_h2=rho_g*(8.*pi*bigG/3.)/H0**2
    cosm%Om_r=Om_g_h2*(1.+0.227*cosm%neff)/cosm%h**2

    !Correction to vacuum density in order for radiation to maintain flatness
    cosm%Om_v=cosm%Om_v-cosm%Om_r

    WRITE(*,*) 'ASSIGN_COSMOLOGY: Cosmology assigned'
    WRITE(*,*) 'ASSIGN_COSMOLOGY: Cosmology: ', TRIM(cosm%name)
    WRITE(*,*)

    CALL init_cosmology(cosm)

  END SUBROUTINE assign_cosmology

  SUBROUTINE print_cosmology(cosm)

    !Prints the cosmological parameters to the screen
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm

    WRITE(*,*) 'COSMOLOGY: ', TRIM(cosm%name)
    WRITE(*,*) '===================================='
    WRITE(*,*) 'COSMOLOGY: Standard parameters'
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_m:', REAL(cosm%Om_m)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_b:', REAL(cosm%Om_b)    
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_v:', REAL(cosm%Om_v)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_w:', REAL(cosm%Om_w)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_r:', REAL(cosm%Om_r)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'h:', REAL(cosm%h)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w_0:', REAL(cosm%w)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w_a:', REAL(cosm%wa)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'sigma_8:', REAL(cosm%sig8)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n_s:', REAL(cosm%n)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'T_CMB [K]:', REAL(cosm%T_CMB)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'z_CMB:', REAL(cosm%z_CMB)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n_eff:', REAL(cosm%neff)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Y_He:', REAL(cosm%YHe)
    IF(cosm%wdm) THEN
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'm_wdm [keV]:', REAL(cosm%m_wdm)
    END IF
    WRITE(*,*) '===================================='
    WRITE(*,*) 'COSMOLOGY: Derived parameters'
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega:', cosm%Om
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_c:', cosm%Om_c
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Omega_k:', cosm%Om_k
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'k [Mpc/h]^-2:', REAL(cosm%k)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Y_H:', REAL(cosm%YH)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'mu_p:', REAL(cosm%mup)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'mu_e:', REAL(cosm%mue)
    !WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'ep_fac:', REAL(cosm%epfac)
    WRITE(*,*) '===================================='
    IF(cosm%iw==1) THEN
       WRITE(*,*) 'COSMOLOGY: Vacuum energy'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w:', -1.
    ELSE IF(cosm%iw==2) THEN
       WRITE(*,*) 'COSMOLOGY: QUICC dark energy prescription'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w0:', cosm%w
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'wm:', cosm%wm
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'am:', cosm%am
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'dm:', cosm%dm
    ELSE IF(cosm%iw==3) THEN
       WRITE(*,*) 'COSMOLOGY: w(a) = w0+wa(1.-a)'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w0:', cosm%w
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'wa:', cosm%wa
    ELSE IF(cosm%iw==4) THEN
       WRITE(*,*) 'COSMOLOGY: Constant w'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w:', cosm%w
    ELSE IF(cosm%iw==5) THEN
       WRITE(*,*) 'COSMOLOGY: IDE I'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a*:', cosm%as
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n*:', cosm%ns
    ELSE IF(cosm%iw==6) THEN
       WRITE(*,*) 'COSMOLOGY: IDE II'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a*:', cosm%as
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'n*:', cosm%ns
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a1^n (derived):', cosm%a1n
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a2^n (derived):', cosm%a2n
    ELSE IF(cosm%iw==7) THEN
       WRITE(*,*) 'COSMOLOGY: IDE III'
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'a*:', cosm%a1
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Om_w(a*):', cosm%Om_ws
       WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'w*:', cosm%ws
    END IF
    WRITE(*,*) '===================================='
    WRITE(*,*) 'COSMOLOGY: HOD'
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'log10(M_gal):', log10(cosm%mgal)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'log10(M_HI_-):', log10(cosm%HImin)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'log10(M_HI_+):', log10(cosm%HImax)
    WRITE(*,*) '===================================='
    WRITE(*,*) 'COSMOLOGY: Baryon model'
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'alpha:', cosm%alpha
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'epsilon:', cosm%eps
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Gamma:', cosm%Gamma
    IF(cosm%M0 .NE. 0.) WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'log10(M0):', log10(cosm%M0)
    WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'Astar:', cosm%Astar
    IF(cosm%whim .NE. 0.) WRITE(*,fmt='(A11,A15,F11.5)') 'COSMOLOGY:', 'log10(WHIM):', log10(cosm%whim)
    WRITE(*,*) '===================================='
    WRITE(*,*)

  END SUBROUTINE print_cosmology

  SUBROUTINE init_cosmology(cosm)

    !Calcualtes derived parameters
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: Xs, f1, f2

    !Derived cosmological parameters    
    cosm%Om_c=cosm%Om_m-cosm%Om_b-cosm%Om_nu
    cosm%Om=cosm%Om_m+cosm%Om_v+cosm%Om_r+cosm%Om_w
    cosm%Om_k=1.-cosm%Om
    cosm%k=(cosm%Om-1.)/(Hdist**2)

    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_c:', REAL(cosm%Om_c)
       WRITE(*,*) 'INIT_COSMOLOGY: Omega:', REAL(cosm%Om)
       WRITE(*,*) 'INIT_COSMOLOGY: Omega_k:', REAL(cosm%Om_k)
       WRITE(*,*) 'INIT_COSMOLOGY: k [Mpc/h]^-2:', REAL(cosm%k)
    END IF

    !Gas parameters
    cosm%YH=1.-cosm%YHe !Hydrogen mass density
    cosm%mup=4./(5.*cosm%YH+3.) !Nuclear mass per particle (~0.588 if fH=0.76)
    cosm%mue=2./(1.+cosm%YH) !Nuclear mass per electron (~1.136 if fH=0.76)
    !cosm%epfac=2.*(1.+cosm%YH)/(5.*cosm%YH+3.) !Thermal -> electron pressure conversion; ratio of electron density to total particle density; (~0.518 if fH=0.76)

    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Y_H:', REAL(cosm%YH)
       WRITE(*,*) 'INIT_COSMOLOGY: mu_p:', REAL(cosm%mup)
       WRITE(*,*) 'INIT_COSMOLOGY: mu_e:', REAL(cosm%mue)
       !WRITE(*,*) 'INIT_COSMOLOGY: ep_fac:', REAL(cosm%epfac)
    END IF
    
    cosm%is_init=.TRUE.

    !Dark energy models
    IF(cosm%iw==5) THEN
       !STOP 'INIT_COSMOLOGY: IDE I does not work yet'
       !WRITE(*,*) H2(astar), X(astar)
       !Om_w=Om_w*(Om_m*astar**(-3)+Om_v)/(X(astar)*(1.-Om_w))
       cosm%Om_w=cosm%Om_ws*(Hubble2(cosm%a,cosm)-cosm%Om_ws*X_de(cosm%as,cosm)+cosm%Om_ws*cosm%as**(-2))/(X_de(cosm%as,cosm)*(1.-cosm%Om_ws)+cosm%Om_ws*cosm%as**(-2))
    ELSE IF(cosm%iw==6) THEN
       !STOP 'INIT_COSMOLOGY: Check this carefully: IDE II'
       !Define a1^n
       cosm%a1n=cosm%as**cosm%ns

       !Necessary for first step below
       cosm%a2n=cosm%a1n

       !All neccessary to convert parameters to a1,a2
       f1=cosm%Om_ws*(Hubble2(cosm%as,cosm)-cosm%Om_w*X_de(cosm%as,cosm))
       f2=cosm%Om_w*(1.-cosm%Om_ws)
       Xs=f1/f2
       !a1=astar**nstar
       !Xstar=Xstar**(-nstar/6.)
       Xs=Xs**(cosm%ns/6.)

       !Top and bottom of fraction
       !f1=2.-a1*Xstar*(1.+1./a1)
       !f2=Xstar*(1.+1./a1)-2.
       f1=cosm%a1n*(2.*Xs-(1.+cosm%a1n))
       f2=(1.+cosm%a1n)-2.*Xs*cosm%a1n

       !Finally! a2
       cosm%a2n=f1/f2
       !IF(a2<a1) a2=a1
    ELSE IF(cosm%iw==7) THEN
       !STOP 'INIT_COSMOLOGY: Check this carefully: IDE III'
       cosm%a1=cosm%as !Scale-factor at which Om_w(a*) is most important
       cosm%a2=cosm%as !Needs to be set for X(a*) and H2(a*) below (which cancel each other)
       f1=cosm%Om_ws*(Hubble2(cosm%as,cosm)-cosm%Om_w*X_de(cosm%as,cosm))
       f2=cosm%Om_w*(1.-cosm%Om_ws)
       cosm%a2=cosm%as*(f1/f2)**(1./(3.*(1.+cosm%ws)))
    END IF

    !Not sure this is necessary
    !Cacluate the age of the universe
    !cosm%age=age_of_universe(cosm)
    !IF(verbose_cosmology) THEN
    !   WRITE(*,*) 'INIT_COSMOLOGY: age of Universe [Gyrs/h]:', REAL(cosm%age)
    !   WRITE(*,*)
    !END IF

    !Ensure deallocate distances
    cosm%has_distance=.FALSE.
    IF(ALLOCATED(cosm%r)) DEALLOCATE(cosm%r)
    IF(ALLOCATED(cosm%a_r)) DEALLOCATE(cosm%a_r)

    !Ensure deallocate growth
    cosm%has_growth=.FALSE.
    IF(ALLOCATED(cosm%log_a_growth))    DEALLOCATE(cosm%log_a_growth)
    IF(ALLOCATED(cosm%log_growth))      DEALLOCATE(cosm%log_growth)
    IF(ALLOCATED(cosm%growth_rate))     DEALLOCATE(cosm%growth_rate)
    IF(ALLOCATED(cosm%log_acc_growth))  DEALLOCATE(cosm%log_acc_growth)

    !Ensure deallocate sigma
    cosm%has_sigma=.FALSE.
    IF(ALLOCATED(cosm%log_r_sigma)) DEALLOCATE(cosm%log_r_sigma)
    IF(ALLOCATED(cosm%log_sigma))   DEALLOCATE(cosm%log_sigma)

    !Ensure deallocate plin
    cosm%has_power=.FALSE.
    IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
    IF(ALLOCATED(cosm%log_plin))   DEALLOCATE(cosm%log_plin)

    cosm%has_spherical=.FALSE.
    IF(ALLOCATED(cosm%log_a_dcDv)) DEALLOCATE(cosm%log_a_dcDv)
    IF(ALLOCATED(cosm%dc))         DEALLOCATE(cosm%dc)
    IF(ALLOCATED(cosm%Dv))         DEALLOCATE(cosm%Dv)

    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_COSMOLOGY: Done'
       WRITE(*,*)
    END IF
    
!!$    !Set the power normalisation to unity initially
!!$    !cosm%is_normalised=.FALSE.
!!$    cosm%A=1.
!!$
!!$    !Calculate the initial sigma_8 value (will not be correct)
!!$    sigi=sqrt(sigma_integral0(8.,1.,cosm,acc_cosm))
!!$
!!$    IF(verbose_cosmology) WRITE(*,*) 'INIT_COSMOLOGY Initial sigma_8:', REAL(sigi)
!!$
!!$    !Reset the normalisation to give the correct sigma8
!!$    cosm%A=cosm%sig8/sigi
!!$    !cosm%A=391.0112 !Appropriate for sig8=0.8 in the boring model (for tests)
!!$
!!$    !Recalculate sigma8, should be correct this time
!!$    sigi=sqrt(sigma_integral0(8.,1.,cosm,acc_cosm))
!!$
!!$    !Write to screen
!!$    IF(verbose_cosmology) THEN
!!$       WRITE(*,*) 'INIT_COSMOLOGY Normalisation factor:', REAL(cosm%A)
!!$       WRITE(*,*) 'INIT_COSMOLOGY Target sigma_8:', REAL(cosm%sig8)
!!$       WRITE(*,*) 'INIT_COSMOLOGY Final sigma_8 (calculated):', REAL(sigi)
!!$       WRITE(*,*) 'INIT_COSMOLOGY Complete'
!!$       WRITE(*,*)
!!$    END IF
!!$
!!$    !Change flag
!!$    cosm%is_normalised=.TRUE.

    !Normalise the power spectrum
    CALL normalise_power(cosm)

  END SUBROUTINE init_cosmology

  SUBROUTINE normalise_power(cosm)

    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: sigi

    !Set the power normalisation to unity initially
    cosm%A=1.

    !Calculate the initial sigma_8 value (will not be correct)
    sigi=sqrt(sigma_integral0(8.,1.,cosm,acc_cosm))

    IF(verbose_cosmology) WRITE(*,*) 'NORMALISE_POWER: Initial sigma_8:', REAL(sigi)

    !Reset the normalisation to give the correct sigma8
    cosm%A=cosm%sig8/sigi
    !cosm%A=391.0112 !Appropriate for sig8=0.8 in the boring model (for tests)

    !Recalculate sigma8, should be correct this time
    sigi=sqrt(sigma_integral0(8.,1.,cosm,acc_cosm))

    !Write to screen
    IF(verbose_cosmology) THEN
       WRITE(*,*) 'NORMALISE_POWER: Normalisation factor:', REAL(cosm%A)
       WRITE(*,*) 'NORMALISE_POWER: Target sigma_8:', REAL(cosm%sig8)
       WRITE(*,*) 'NORMALISE_POWER: Final sigma_8 (calculated):', REAL(sigi)
       WRITE(*,*) 'NORMALISE_POWER: Done'
       WRITE(*,*)
    END IF

    !Change the flag
    cosm%is_normalised=.TRUE.

  END SUBROUTINE normalise_power

  FUNCTION comoving_critical_density(a,cosm)

    !Comoving critical density in (Msun/h) / (Mpc/h)^3
    IMPLICIT NONE
    REAL :: comoving_critical_density
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    comoving_critical_density=physical_critical_density(a,cosm)*a**3

  END FUNCTION comoving_critical_density

  FUNCTION physical_critical_density(a,cosm)

    !Physical critical density in (Msun/h) / (Mpc/h)^3
    IMPLICIT NONE
    REAL :: physical_critical_density
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_critical_density=critical_density*Hubble2(a,cosm)

  END FUNCTION physical_critical_density

  FUNCTION comoving_matter_density(cosm)

    !Comoving matter density in (Msun/h) / (Mpc/h)^3
    !Not a function of redshift!
    IMPLICIT NONE
    REAL :: comoving_matter_density
    TYPE(cosmology), INTENT(INOUT) :: cosm

    comoving_matter_density=critical_density*cosm%Om_m

  END FUNCTION comoving_matter_density

  FUNCTION physical_matter_density(a,cosm)

    !Physical matter density in (Msun/h) / (Mpc/h)^3
    IMPLICIT NONE
    REAL :: physical_matter_density
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_matter_density=comoving_matter_density(cosm)*a**(-3)

  END FUNCTION physical_matter_density

  FUNCTION Hubble2(a,cosm)

    !Calculates Hubble^2 in units such that H^2(z=0)=1.
    IMPLICIT NONE
    REAL :: Hubble2
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'HUBBLE2: Error, cosmology is not initialised'
    Hubble2=cosm%Om_m*a**(-3)+cosm%Om_nu*a**(-3)+cosm%Om_r*a**(-4)+cosm%Om_v+cosm%Om_w*X_de(a,cosm)+(1.-cosm%om)*a**(-2)

  END FUNCTION Hubble2

  FUNCTION Hubble2_norad(a,cosm)

    IMPLICIT NONE
    REAL :: Hubble2_norad
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Hubble2_norad=Hubble2(a,cosm)-cosm%Om_r*a**(-4)
    
  END FUNCTION Hubble2_norad

  FUNCTION Hubble2a4_highz(cosm)

    !Calculates Hubble^2a^4 in units such that H^2(z=0)=1.
    !This is only valid at high z, when only radiation is important
    !Makes some assumptions that DE is *not* important at high z
    !Need to worry if Omega_de is scaling anything like a^-4 (e.g., kinetic dominated a^-6)
    IMPLICIT NONE
    REAL :: Hubble2a4_highz
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'HUBBLE2A4_HIGHZ: Error, cosmology is not initialised'
    Hubble2a4_highz=cosm%Om_r

  END FUNCTION Hubble2a4_highz

  FUNCTION AH(a,cosm)

    !\ddot{a}/a
    IMPLICIT NONE
    REAL :: AH
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'AH: Error, cosmology is not initialised'
    AH=cosm%Om_m*a**(-3)+cosm%Om_nu*a**(-3)+2.*cosm%Om_r*a**(-4)-2.*cosm%Om_v+cosm%Om_w*(1.+3.*w_de(a,cosm))*X_de(a,cosm)
    AH=-AH/2.

  END FUNCTION AH

  FUNCTION AH_norad(a,cosm)

    IMPLICIT NONE
    REAL :: AH_norad
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    AH_norad=AH(a,cosm)+cosm%Om_r*a**(-4)
    
  END FUNCTION AH_norad

  REAL FUNCTION Omega_m(a,cosm)

    !This calculates Omega_m variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_m=cosm%Om_m*a**(-3)/Hubble2(a,cosm)

  END FUNCTION Omega_m

  REAL FUNCTION Omega_m_norad(a,cosm)

    !This calculates Omega_m variations with a, but ignoring any radiation component
    !This ensures that Omega_m_norad(a->0) -> 1
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_m_norad=cosm%Om_m*a**(-3)/Hubble2_norad(a,cosm)

  END FUNCTION Omega_m_norad

  REAL FUNCTION Omega_b(a,cosm)

    !This calculates Omega_m variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_b=cosm%Om_b*a**(-3)/Hubble2(a,cosm)

  END FUNCTION Omega_b

  REAL FUNCTION Omega_c(a,cosm)

    !This calculates Omega_c variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_c=cosm%Om_c*a**(-3)/Hubble2(a,cosm)

  END FUNCTION Omega_c

  REAL FUNCTION Omega_nu(a,cosm)

    !This calculates Omega_nu variations with a, assuming that nu scales like matter
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_nu=cosm%Om_nu*a**(-3)/Hubble2(a,cosm)

  END FUNCTION Omega_nu

  REAL FUNCTION Omega_r(a,cosm)

    !This calculates Omega_r variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'Omega_r: Error, cosmology is not initialised'
    Omega_r=cosm%Om_r*a**(-4)/Hubble2(a,cosm)

  END FUNCTION Omega_r

  REAL FUNCTION Omega_v(a,cosm)

    !This calculates Omega_v variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_v=cosm%Om_v/Hubble2(a,cosm)

  END FUNCTION Omega_v

  REAL FUNCTION Omega_w(a,cosm)

    !This calculates Omega_w variations with a
    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Omega_w=cosm%Om_w*X_de(a,cosm)/Hubble2(a,cosm)

  END FUNCTION Omega_w

  REAL FUNCTION Omega(a,cosm)

    !This calculates total Omega variations with a
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%is_init .EQV. .FALSE.) STOP 'OMEGA: Error, cosmology is not initialised'
    Omega=Omega_m(a,cosm)+Omega_nu(a,cosm)+Omega_r(a,cosm)+Omega_v(a,cosm)+Omega_w(a,cosm)

  END FUNCTION Omega

!!$  FUNCTION w_de(a,cosm)
!!$
!!$    !w(a) for DE models
!!$    IMPLICIT NONE
!!$    REAL :: w_de
!!$    REAL, INTENT(IN) :: a
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    w_de=cosm%w+(1.-a)*cosm%wa
!!$
!!$  END FUNCTION w_de
!!$
!!$  FUNCTION X_de(a,cosm)
!!$
!!$    !The time evolution for Om_w for w(a) DE models
!!$    !X_de = rho_w(z)/rho_w(z=0)
!!$    IMPLICIT NONE
!!$    REAL :: X_de
!!$    REAL, INTENT(IN) :: a
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$
!!$    X_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))
!!$
!!$  END FUNCTION X_de

  FUNCTION w_de(a,cosm)

    !Variations of the dark energy parameter w(a)
    IMPLICIT NONE
    REAL :: w_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: p1, p2, p3, p4
    DOUBLE PRECISION :: f1, f2, f3, f4

    IF(cosm%iw==1) THEN
       !LCDM
       w_de=-1.
    ELSE IF(cosm%iw==2) THEN
       !QUICC parameterisation
       p1=1.+exp(cosm%am/cosm%dm)
       p2=1.-exp(-(a-1.)/cosm%dm)
       p3=1.+exp(-(a-cosm%am)/cosm%dm)
       p4=1.-exp(1./cosm%dm)
       w_de=cosm%w+(cosm%wm-cosm%w)*p1*p2/(p3*p4)
    ELSE IF(cosm%iw==3) THEN
       !w(a)CDM
       w_de=cosm%w+(1.-a)*cosm%wa
    ELSE IF(cosm%iw==4) THEN
       !wCDM
       w_de=cosm%w
    ELSE IF(cosm%iw==5) THEN
       !IDE I
       w_de=((a/cosm%as)**cosm%ns-1.)/((a/cosm%as)**cosm%ns+1.)
    ELSE IF(cosm%iw==6) THEN
       !IDE II
       f1=a**cosm%ns-cosm%a1n
       f2=a**cosm%ns+cosm%a1n
       f3=a**cosm%ns-cosm%a2n
       f4=a**cosm%ns+cosm%a2n
       w_de=-1.+REAL(f1/f2-f3/f4)
    ELSE IF(cosm%iw==7) THEN
       !IDE III
       IF(a<cosm%a1) THEN
          w_de=-1.
       ELSE IF(cosm%a1<=a .AND. a<cosm%a2) THEN
          w_de=cosm%ws
       ELSE IF(a>=cosm%a2) THEN
          w_de=-1.
       ELSE
          STOP 'W_DE: Error, something went wrong'
       END IF
    ELSE
       STOP 'W_DE: Error, value of iw set incorrectly'
    END IF

  END FUNCTION w_de

  FUNCTION w_de_total(a,cosm)

    !Do an average over the DE components
    IMPLICIT NONE
    REAL :: w_de_total
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%Om_nu .NE. 0.) STOP 'W_DE_TOTAL: Error, does not support massive neutrinos'

    IF(cosm%Om_v==0. .AND. cosm%Om_w==0.) THEN
       w_de_total=-1.
    ELSE
       w_de_total=w_de(a,cosm)*Omega_w(a,cosm)-Omega_v(a,cosm)
       w_de_total=w_de_total/(Omega_w(a,cosm)+Omega_v(a,cosm))
    END IF

  END FUNCTION w_de_total

  FUNCTION w_eff(a,cosm)

    !Do an average over the DE components
    IMPLICIT NONE
    REAL :: w_eff
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%Om_nu .NE. 0.) STOP 'W_EFF: Error, does not support massive neutrinos'

    w_eff=w_de(a,cosm)*Omega_w(a,cosm)-Omega_v(a,cosm)+Omega_r(a,cosm)/3.
    w_eff=w_eff/Omega(a,cosm)

  END FUNCTION w_eff

  FUNCTION X_de(a,cosm)

    !Redshift scaling for dark energy (i.e., if w=0 x(a)=a^-3, if w=-1 x(a)=const etc.)
    IMPLICIT NONE
    REAL :: X_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    DOUBLE PRECISION :: f1, f2, f3, f4
    !REAL, PARAMETER :: acc=1e-3

    IF(cosm%iw==1) THEN
       !LCDM
       X_de=1.
    ELSE IF(cosm%iw==2) THEN
       STOP 'X_DE: Error, need to support explicit integration here'
    ELSE IF(cosm%iw==3) THEN
       !w(a)CDM
       X_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))
    ELSE IF(cosm%iw==4) THEN
       !wCDM
       X_de=a**(-3.*(1.+cosm%w))
    ELSE IF(cosm%iw==5) THEN
       !IDE I
       X_de=((1.+(a/cosm%as)**cosm%ns)/(1.+(1./cosm%as)**cosm%ns))**(-6./cosm%ns)
    ELSE IF(cosm%iw==6) THEN
       !IDE II
       f1=a**cosm%ns+cosm%a1n
       f2=1.+cosm%a1n
       f3=1.+cosm%a2n
       f4=a**cosm%ns+cosm%a2n
       X_de=REAL(f1*f3/(f2*f4))**(-6./cosm%ns)
    ELSE IF(cosm%iw==7) THEN
       !IDE III
       IF(a<cosm%a1) THEN
          X_de=(cosm%a1/cosm%a2)**(-3.*(1.+cosm%ws))
       ELSE IF(cosm%a1<=a .AND. a<cosm%a2) THEN
          X_de=(a/cosm%a2)**(-3.*(1.+cosm%ws))
       ELSE IF(a>=cosm%a2) THEN
          X_de=1.
       ELSE
          STOP 'X_DE: Error, something went wrong'
       END IF
    ELSE
       !Generally true, doing this integration can make calculations very slow
       !Difficult to implement into library because of cosm dependence of w_de
       !See the commented out bit below
       !X_de=(a**(-3))*exp(3.*integrate_log(a,1.,integrand_de,acc,3,1))
       STOP 'X_DE: Error, this integration routine has not been tested'
       X_de=(a**(-3))*exp(3.*integrate_cosm(a,1.,integrand_de,cosm,acc_cosm,3))
!!$    ELSE
!!$       STOP 'X_DE: Error, iw not specified correctly'
    END IF

  END FUNCTION X_de

  FUNCTION integrand_de(a,cosm)

    !The integrand for the X_de(a) integral
    IMPLICIT NONE
    REAL :: integrand_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    integrand_de=w_de(a,cosm)/a

  END FUNCTION integrand_de

  FUNCTION redshift_a(a)

    !The redshift corresponding to scale-factor a
    IMPLICIT NONE
    REAL :: redshift_a
    REAL, INTENT(IN) :: a

    IF(a==0. .OR. a>1.) THEN
       WRITE(*,*) 'REDSHIFT_A: a', a
       STOP 'REDSHIFT_A: Error, routine called with weird a'
    END IF

    !IF(a==0.) THEN
    !   WRITE(*,*) 'REDSHIFT_A: a', a
    !   STOP 'REDSHIFT_A: Error, routine called with a=0'
    !END IF

    redshift_a=-1.+1./a

  END FUNCTION redshift_a

  FUNCTION scale_factor_z(z)

    !The scale factor corresponding to redshift z
    IMPLICIT NONE
    REAL :: scale_factor_z
    REAL, INTENT(IN) :: z

    IF(z<0.) THEN
       WRITE(*,*) 'SCALE_FACTOR_Z: z', z
       STOP 'SCALE_FACTOR_Z: Error, routine called for z<0'
    END IF

    scale_factor_z=1./(1.+z)

  END FUNCTION scale_factor_z

  FUNCTION redshift_r(r,cosm)

    !The redshift corresponding to comoving distance r
    IMPLICIT NONE
    REAL :: redshift_r
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_distance .EQV. .FALSE.) CALL init_distances(cosm)
    redshift_r=redshift_a(find(r,cosm%r,cosm%a_r,cosm%n_r,3,3,2))

  END FUNCTION redshift_r

  FUNCTION f_k(r,cosm)

    !Curvature function
    IMPLICIT NONE
    REAL :: f_k
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%k==0.) THEN
       f_k=r
    ELSE IF(cosm%k<0.) THEN
       f_k=sinh(sqrt(-cosm%k)*r)/sqrt(-cosm%k)
    ELSE IF(cosm%k>0.) THEN
       f_k=sin(sqrt(cosm%k)*r)/sqrt(cosm%k)
    ELSE
       STOP 'F_K: Something went wrong'
    END IF

  END FUNCTION f_k

  FUNCTION fdash_k(r,cosm)

    !Derivative of curvature function
    IMPLICIT NONE
    REAL :: fdash_k
    REAL, INTENT(IN) :: r
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%k==0.) THEN
       fdash_k=1.
    ELSE IF(cosm%k<0.) THEN
       fdash_k=cosh(sqrt(-cosm%k)*r)
    ELSE IF(cosm%k>0.) THEN
       fdash_k=cos(sqrt(cosm%k)*r)
    ELSE
       STOP 'FDASH_K: Something went wrong'
    END IF

  END FUNCTION fdash_k

  FUNCTION distance_integrand(a,cosm)

    !The integrand for the cosmic-distance calculation
    IMPLICIT NONE
    REAL :: distance_integrand
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: amin=1e-5

    IF(a<amin) THEN
       distance_integrand=Hdist/sqrt(Hubble2a4_highz(cosm))
    ELSE
       distance_integrand=Hdist/(sqrt(Hubble2(a,cosm))*a**2)
    END IF

  END FUNCTION distance_integrand

  FUNCTION comoving_distance(a,cosm)

    !The comoving distance to a galaxy at scale-factor a
    IMPLICIT NONE
    REAL :: comoving_distance
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_distance .EQV. .FALSE.) CALL init_distances(cosm)
    comoving_distance=find(a,cosm%a_r,cosm%r,cosm%n_r,3,3,2)

  END FUNCTION comoving_distance

  FUNCTION physical_distance(a,cosm)

    !The physical distance to a galaxy at scale-factor a
    IMPLICIT NONE
    REAL :: physical_distance
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_distance=comoving_distance(a,cosm)*a

  END FUNCTION physical_distance

  FUNCTION physical_angular_distance(a,cosm)

    !The physical angular-diameter distance to a galaxy at scale-factor a
    IMPLICIT NONE
    REAL :: physical_angular_distance
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    physical_angular_distance=f_k(physical_distance(a,cosm),cosm)

  END FUNCTION physical_angular_distance

  FUNCTION comoving_angular_distance(a,cosm)

    !The physical angular-diameter distance to a galaxy at scale-factor a
    IMPLICIT NONE
    REAL :: comoving_angular_distance
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    comoving_angular_distance=a*physical_angular_distance(a,cosm)

  END FUNCTION comoving_angular_distance

  FUNCTION luminosity_distance(a,cosm)

    !The luminosity distance to a galaxy at scale-factor a
    IMPLICIT NONE
    REAL :: luminosity_distance
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    luminosity_distance=f_k(comoving_distance(a,cosm),cosm)/a

  END FUNCTION luminosity_distance

  SUBROUTINE init_distances(cosm)

    !Fill up tables of a vs. r(a) (comoving distance)
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: zmin, zmax, amin, amax
    INTEGER :: i

    INTEGER, PARAMETER :: nr=128

    zmin=0.
    zmax=cosm%z_CMB
    amin=scale_factor_z(zmax)
    amax=scale_factor_z(zmin)
    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_DISTANCE: Redshift range for r(z) tables'
       WRITE(*,*) 'INIT_DISTANCE: minimum z:', REAL(zmin)
       WRITE(*,*) 'INIT_DISTANCE: maximum z:', REAL(zmax)
       WRITE(*,*) 'INIT_DISTANCE: minimum a:', REAL(amin)
       WRITE(*,*) 'INIT_DISTANCE: maximum a:', REAL(amax)
    END IF
    cosm%n_r=nr
    CALL fill_array(amin,amax,cosm%a_r,cosm%n_r)
    IF(ALLOCATED(cosm%r)) DEALLOCATE(cosm%r)
    ALLOCATE(cosm%r(cosm%n_r))

    !Now do the r(z) calculation
    DO i=1,cosm%n_r
       cosm%r(i)=integrate_cosm(cosm%a_r(i),1.,distance_integrand,cosm,acc_cosm,3)
    END DO
    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_DISTANCE: minimum r [Mpc/h]:', REAL(cosm%r(cosm%n_r))
       WRITE(*,*) 'INIT_DISTANCE: maximum r [Mpc/h]:', REAL(cosm%r(1))
    END IF

    !Find the horizon distance in your cosmology
    cosm%horizon=integrate_cosm(0.,1.,distance_integrand,cosm,acc_cosm,3)
    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_DISTANCE: Horizon distance [Mpc/h]:', REAL(cosm%horizon)
       WRITE(*,*) 'INIT_DISTANCE: Done'
       WRITE(*,*)
    END IF

    cosm%has_distance=.TRUE.

  END SUBROUTINE init_distances

  SUBROUTINE write_distances(cosm)

    !Write file of z vs. r(z)
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    CHARACTER(len=256) :: output
    INTEGER :: i
    REAL :: z

    !Now write the results of r(z) calculation
    output='projection/distance.dat'
    WRITE(*,*) 'WRITE_DISTANCE: Writing r(a): ', TRIM(output)
    OPEN(7,file=output)
    DO i=1,cosm%n_r
       z=redshift_a(cosm%a_r(i))
       WRITE(7,*) z, cosm%r(i), f_k(cosm%r(i),cosm)
    END DO
    CLOSE(7)
    WRITE(*,*) 'WRITE_DISTANCE: Done'
    WRITE(*,*)

  END SUBROUTINE write_distances

  FUNCTION age_of_universe(cosm)

    IMPLICIT NONE
    REAL :: age_of_universe
    TYPE(cosmology), INTENT(INOUT) :: cosm

    age_of_universe=cosmic_time(1.,cosm)

  END FUNCTION age_of_universe

  FUNCTION cosmic_time(a,cosm)

    !The age of the universe at scale-factor 'a'
    IMPLICIT NONE
    REAL :: cosmic_time
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    cosmic_time=integrate_cosm(0.,a,time_integrand,cosm,acc_cosm,3)

  END FUNCTION cosmic_time

  FUNCTION look_back_time(a,cosm)

    !The time in the past that photons at scale-factor 'a' were emitted
    IMPLICIT NONE
    REAL :: look_back_time
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    look_back_time=integrate_cosm(a,1.,time_integrand,cosm,acc_cosm,3)

  END FUNCTION look_back_time

  FUNCTION time_integrand(a,cosm)

    !The integrand for the cosmic-distance calculation
    IMPLICIT NONE
    REAL :: time_integrand
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: amin=1e-5

    IF(a<amin) THEN
       time_integrand=a*Htime/sqrt(Hubble2a4_highz(cosm))
    ELSE
       time_integrand=Htime/(a*sqrt(Hubble2(a,cosm)))
    END IF

  END FUNCTION time_integrand

!!$  SUBROUTINE random_cosmology(cosm)
!!$
!!$    !Generate some random cosmological parameter
!!$    USE random_numbers
!!$    IMPLICIT NONE
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    REAL :: Om_m_min, Om_m_max, Om_b_min, Om_b_max, n_min, n_max
!!$    REAL :: w_min, w_max, h_min, h_max, sig8_min, sig8_max, wa_min, wa_max
!!$
!!$    !Needs to be set to normalise P_lin
!!$    cosm%A=1.
!!$
!!$    Om_m_min=0.1
!!$    Om_m_max=1.
!!$    cosm%Om_m=random_uniform(Om_m_min,Om_m_max)
!!$
!!$    cosm%Om_v=1.-cosm%Om_m
!!$
!!$    Om_b_min=0.005
!!$    Om_b_max=MIN(0.095,cosm%Om_m)
!!$    cosm%Om_b=random_uniform(Om_b_min,Om_b_max)
!!$
!!$    cosm%Om_c=cosm%Om_m-cosm%Om_b
!!$
!!$    n_min=0.5
!!$    n_max=1.5
!!$    cosm%n=random_uniform(n_min,n_max)
!!$
!!$    h_min=0.4
!!$    h_max=1.2
!!$    cosm%h=random_uniform(h_min,h_max)
!!$
!!$    w_min=-1.5
!!$    w_max=-0.5
!!$    cosm%w=random_uniform(w_min,w_max)
!!$
!!$    wa_min=-1.
!!$    wa_max=-cosm%w*0.8
!!$    cosm%wa=random_uniform(wa_min,wa_max)
!!$
!!$    sig8_min=0.2
!!$    sig8_max=1.5
!!$    cosm%sig8=random_uniform(sig8_min,sig8_max)
!!$
!!$  END SUBROUTINE random_cosmology

  FUNCTION Tk(k,cosm)

    !Transfer function selection
    IMPLICIT NONE
    REAL :: Tk, k
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Tk=Tk_eh(k,cosm)

    !Damp if considering WDM
    IF(cosm%wdm) Tk=Tk*Tk_wdm(k,cosm)

  END FUNCTION Tk

  REAL FUNCTION Tk_defw(k,cosm)

    !This function was written by John Peacock
    !The DEFW transfer function approximation
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: keff, q, tk
    DOUBLE PRECISION :: q8, tk8

    keff=0.172+0.011*log(cosm%gamma/0.36)*log(cosm%gamma/0.36)
    q=1.e-20 + k/cosm%gamma
    q8=1.e-20 + keff/cosm%gamma
    tk=1./(1.+(6.4*q+(3.0*q)**1.5+(1.7*q)**2)**1.13)**(1./1.13)
    tk8=1./(1.+(6.4*q8+(3.0*q8)**1.5+(1.7*q8)**2)**1.13)**(1./1.13)

    tk_defw=tk/REAL(tk8)

  END FUNCTION Tk_defw

  FUNCTION Tk_eh(yy,cosm)

    USE special_functions

    !Eisenstein & Hu fitting function (arXiv: 9709112)
    !JP - the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
    !JP - remember I use k/h, whereas they use pure k, Om_m is cdm + baryons
    IMPLICIT NONE
    REAL :: Tk_eh
    REAL, INTENT(IN) :: yy
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: rk, e, thet, b1, b2, zd, ze, rd, re, rke, s, rks
    REAL :: q
    REAL :: y, g, ab
    REAL :: a1, a2, ac
    REAL :: bc
    REAL :: f, fac
    REAL :: c1, c2, tc
    REAL :: bb, bn, ss, tb
    REAL :: Om_m, Om_b, h, Omh2, Obh2

    !Define some useful variables
    Om_m=cosm%Om_m
    Om_b=cosm%Om_b
    h=cosm%h
    Omh2=Om_m*h**2
    Obh2=Om_b*h**2

    !Wave-number
    rk=yy*h

    !2.718...
    e=exp(1.)

    !CMB temperature (Section 2)
    thet=cosm%T_CMB/2.7
    
    b1=0.313*(Omh2)**(-0.419)*(1.+0.607*(Omh2)**0.674) !Equation (4)
    b2=0.238*(Omh2)**0.223 !Equation (4)
    zd=1291.*(1.+b1*(Obh2)**b2)*(Omh2)**0.251/(1.+0.659*(Omh2)**0.828) !Drag redshift; equation (4)
    ze=2.50e4*Omh2/thet**4 !z_eq; equation (2)
    rd=31500.*Obh2/thet**4/(1.+zd) !Equation (5); changed from /zd -> /(1+zd) because it lines up with http://background.uchicago.edu/~whu/transfer/tf_fit.c (thanks Steven Murray)
    re=31500.*Obh2/thet**4/ze !Equation (5)
    rke=7.46e-2*Omh2/thet**2 !k_eq; equation (3)
    s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1.+sqrt(re))) !Sound horizon at drag; equation (6)
    rks=1.6*( (Obh2)**0.52 ) * ( (Omh2)**0.73 ) * (1.+(10.4*Omh2)**(-0.95)) !Silk k; equation(7)

    q=rk/13.41/rke !Equation (10)

    y=(1.+ze)/(1.+zd) !y that enters in equation G(y) in equations (14) and (15)
    g=y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.))) !Equation (15)
    ab=g*2.07*rke*s/(1.+rd)**(0.75) !Equation (14)

    a1=(46.9*Omh2)**0.670*(1+(32.1*Omh2)**(-0.532)) !Equation (11)
    a2=(12.0*Omh2)**0.424*(1+(45.0*Omh2)**(-0.582)) !Equation (11)
    ac=(a1**(-Om_b/Om_m)) * (a2**(-(Om_b/Om_m)**3)) !Equation (11)

    b1=0.944/(1.+(458.*Omh2)**(-0.708)) !Equation (12)
    b2=(0.395*Omh2)**(-0.0266) !Equation (12)
    bc=1./(1.+b1*((1.-Om_b/Om_m)**b2-1.)) !Equation (12)

    f=1./(1.+(rk*s/5.4)**4) !Equation (18)

    c1=14.2 + 386./(1.+69.9*q**1.08) !Equation (20) without alpha_c in the 14.2/alpha_c first bit
    c2=14.2/ac + 386./(1.+69.9*q**1.08) !Equation (20) (C function should have explicity alpha_c dependence in paper)
    tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q**2) +(1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q**2) !Equation (17)

    bb=0.5+(Om_b/Om_m) + (3.-2.*Om_b/Om_m)*sqrt((17.2*Omh2)**2+1.) !Equation (24)
    bn=8.41*(Omh2)**0.435 !Equation (23)
    ss=s/(1.+(bn/rk/s)**3)**(1./3.) !Equation (22)
    tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q**2)/(1.+(rk*s/5.2)**2) !First term in equation (21)
    !IF((rk/rks**1.4)>7.) THEN
    !   fac=0.
    !ELSE
    !Removed this IF statement as it produced a discontinuity in P_lin(k) as cosmology
    !was varied - thanks David Copeland for pointing this out
    fac=exp(-(rk/rks)**1.4) !Silk-damping factor from equation (21)
    !END IF
    
    !tb=(tb+ab*fac/(1.+(bb/rk/s)**3))*sin(rk*ss)/rk/ss !Equation (21)
    tb=(tb+ab*fac/(1.+(bb/rk/s)**3))*sinc(rk*ss) !Equation (21)

    tk_eh=real((Om_b/Om_m)*tb+(1.-Om_b/Om_m)*tc) !The weighted mean of baryon and CDM transfer functions

  END FUNCTION TK_EH

  REAL FUNCTION Tk_wdm(k,cosm)

    !Warm dark matter 'correction' to the standard transfer function
    !This version and equation references were taken from arxiv:1605.05973
    !Originally from Bode et al. (2001; arixv:0010389)
    IMPLICIT NONE
    REAL, INTENT(IN) :: k
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: alpha, mu

    alpha=0.074*0.7*cosm%m_wdm**(-1.15) !alpha from equation (5), units Mpc/h
    mu=1.12 !mu from equation (4), dimensionless

    Tk_wdm=(1.+(alpha*k)**(2.*mu))**(-5./mu) !Equation (2)
    
  END FUNCTION Tk_wdm

  FUNCTION p_lin(k,a,cosm)

    !Linear matter power spectrum
    !P(k) should have been previously normalised so as to get the amplitude 'A' correct
    IMPLICIT NONE
    REAL :: p_lin
    REAL, INTENT (IN) :: k, a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    REAL, PARAMETER :: kmax=1e8

    !Using init_power seems to provide no significant speed improvements to HMx
    !IF((cosm%has_power .EQV. .FALSE.) .AND. (cosm%external_plin .EQV. .FALSE.)) CALL init_power(cosm)

    IF(k==0.) THEN
       !If p_lin happens to be foolishly called for 0 mode
       !This call should never happen, but may in integrals
       p_lin=0.
    ELSE IF(k>kmax) THEN
       !Avoids some issues if p_lin is called for very (absurdly) high k values
       !For some reason crashes can occur if this is the case
       p_lin=0.
    ELSE IF(ibox==1 .AND. k<2.*pi/Lbox) THEN
       !If investigating effects caused by a finite box size
       p_lin=0.
    ELSE
       IF(cosm%external_plin .OR. cosm%has_power) THEN
          p_lin=cosm%A**2*exp(find(log(k),cosm%log_k_plin,cosm%log_plin,cosm%n_plin,3,3,2))
       ELSE
          !In this case get the power from the transfer function
          p_lin=(cosm%A**2)*(Tk(k,cosm)**2)*(k**(cosm%n+3.))
       END IF
       !'Grow' the power from z=0 to the redshift of interest
       p_lin=p_lin*grow(a,cosm)**2
    END IF

  END FUNCTION p_lin

  SUBROUTINE init_power(cosm)

    !Fill a look-up table for the linear power spectrum
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    REAL :: k

    !PARAMETERS
    REAL, PARAMETER :: kmin=1e-3
    REAL, PARAMETER :: kmax=1e2
    INTEGER, PARAMETER :: n_plin=256

    !Set the numer of points in the look-up tables
    !Note you need to have enough to resolve the BAO well
    !Probably some non-log/linear spacing would be best (CAMB does this)
    cosm%n_plin=n_plin

    !Allocate arrays
    IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
    IF(ALLOCATED(cosm%log_plin)) DEALLOCATE(cosm%log_plin)
    ALLOCATE(cosm%log_k_plin(cosm%n_plin),cosm%log_plin(cosm%n_plin))

    !Get values for the linear power spectrum
    DO i=1,cosm%n_plin
       k=progression_log(kmin,kmax,i,cosm%n_plin)
       cosm%log_k_plin(i)=k
       cosm%log_plin(i)=(Tk(k,cosm)**2)*k**(cosm%n+3.)
    END DO

    !Take logarithms
    cosm%log_k_plin=log(cosm%log_k_plin)
    cosm%log_plin=log(cosm%log_plin)

    !Change the flag to true
    cosm%has_power=.TRUE.
    
  END SUBROUTINE init_power

  SUBROUTINE init_sigma(cosm)

    !This fills up tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigma(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: r, sigma
    INTEGER :: i

    !These values of 'r' work fine for any power spectrum of cosmological importance
    !Having nsig as a 2** number is most efficient for the look-up routines
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is a power-law below and above these values of R   
    INTEGER, PARAMETER :: nsig=128 !Number of entries for sigma(R) tables
    REAL, PARAMETER :: rmin=1e-4 !Minimum r value (NB. sigma(R) needs to be power-law below)
    REAL, PARAMETER :: rmax=1e3 !Maximum r value (NB. sigma(R) needs to be power-law above)
    REAL, PARAMETER :: rsplit=1e-2 !Scale split between integration methods

    !IF(cosm%is_normalised .EQV. .FALSE.) CALL normalise_power(cosm)

    IF(cosm%wdm) STOP 'INIT_SIGMA: This will crash with WDM'

    !Deallocate tables if they are already allocated
    IF(ALLOCATED(cosm%log_r_sigma)) DEALLOCATE(cosm%log_r_sigma)
    IF(ALLOCATED(cosm%log_sigma))   DEALLOCATE(cosm%log_sigma)

    !Allocate arrays
    cosm%n_sigma=nsig
    ALLOCATE(cosm%log_r_sigma(nsig),cosm%log_sigma(nsig))

    !Write crud
    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_SIGMA: Filling sigma interpolation table'
       WRITE(*,*) 'INIT_SIGMA: R minimum [Mpc/h]:', REAL(rmin) !If I put REAL() here I get an error with goslow for some reason!?
       WRITE(*,*) 'INIT_SIGMA: R maximum [Mpc/h]:', REAL(rmax)
       WRITE(*,*) 'INIT_SIGMA: number of points:', nsig
    END IF

    DO i=1,nsig

       !Equally spaced r in log
       r=progression_log(rmin,rmax,i,nsig)

       !Integration method changes depending on r to make this as fast as possible
       IF(r>=rsplit) THEN
          sigma=sqrt(sigma_integral0(r,1.,cosm,acc_cosm))
       ELSE IF(r<rsplit) THEN
          sigma=sqrt(sigma_integral1(r,1.,cosm,acc_cosm)+sigma_integral2(r,1.,cosm,acc_cosm))
       ELSE
          STOP 'INIT_SIGMA: Error, something went wrong'
       END IF

       !Fill look-up tables
       cosm%log_r_sigma(i)=r
       cosm%log_sigma(i)=sigma

    END DO

    !Make the tables logarithmic
    cosm%log_r_sigma=log(cosm%log_r_sigma)
    cosm%log_sigma=log(cosm%log_sigma)

    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_SIGMA: Done'
       WRITE(*,*)
    END IF

    !Change flag so that it is known that the look-up tables are filled
    cosm%has_sigma=.TRUE.

  END SUBROUTINE init_sigma

  FUNCTION sigma(r,a,cosm)

    !Finds sigma_cold from look-up tables
    IMPLICIT NONE
    REAL :: sigma
    REAL, INTENT(IN) :: r, a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !IF(cosm%is_normalised .EQV. .FALSE.) CALL normalise_power(cosm)
    IF(cosm%has_sigma .EQV. .FALSE.) CALL init_sigma(cosm)
    sigma=grow(a,cosm)*exp(find(log(r),cosm%log_r_sigma,cosm%log_sigma,cosm%n_sigma,3,3,2))

  END FUNCTION sigma

  FUNCTION sigma_integrand(k,R,a,cosm)

    !The integrand for the sigma(R) integrals
    USE special_functions
    IMPLICIT NONE
    REAL :: sigma_integrand
    REAL, INTENT(IN) :: k, R, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: y, w_hat

    IF(k==0.) THEN
       sigma_integrand=0.
    ELSE
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand=p_lin(k,a,cosm)*(w_hat**2)/k
    END IF

  END FUNCTION sigma_integrand

  FUNCTION sigma_integrand_transformed(t,R,f,a,cosm)

    !The integrand for the sigma(R) integrals
    USE special_functions
    IMPLICIT NONE
    REAL :: sigma_integrand_transformed
    REAL, INTENT(IN) :: t, R, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: k, y, w_hat

    INTERFACE
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    !Integrand to the sigma integral in terms of t. Defined by k=(1/t-1)/f(R) where f(R) is *any* function

    IF(t==0.) THEN
       !t=0 corresponds to k=infintiy when W(kR)=0
       sigma_integrand_transformed=0.
    ELSE IF(t==1.) THEN
       !t=1 corresponds to k=0 when P(k)=0
       sigma_integrand_transformed=0.
    ELSE
       !f(R) can be *any* function of R here to improve integration speed
       k=(-1.+1./t)/f(R)
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand_transformed=p_lin(k,a,cosm)*(w_hat**2)/(t*(1.-t))
    END IF

  END FUNCTION sigma_integrand_transformed

  FUNCTION sigma_integral0(r,a,cosm,acc)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigma_integral0
    REAL, INTENT(IN) :: r, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, INTENT(IN) :: acc    
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    REAL, PARAMETER :: b=0. !Integration lower limit (corresponts to k=inf)
    REAL, PARAMETER :: c=1. !Integration upper limit (corresponds to k=0)
    INTEGER, PARAMETER :: iorder=3

    IF(b==c) THEN

       !Fix the answer to zero if the integration limits are identical
       sigma_integral0=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(c-b)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=sigma_integrand_transformed(b,r,f0_rapid,a,cosm)
             f2=sigma_integrand_transformed(c,r,f0_rapid,a,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(b,c,i,n)
                fx=sigma_integrand_transformed(x,r,f0_rapid,a,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGMA_INTEGRAL0: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             sigma_integral0=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             sigma_integral0=0.d0
             STOP 'SIGMA_INTEGRAL0: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sigma_integral0=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION sigma_integral0

  FUNCTION f0_rapid(r)

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error
    IMPLICIT NONE
    REAL :: f0_rapid
    REAL, INTENT(IN) :: r
    REAL :: alpha

    REAL, PARAMETER :: rsplit=1e-2

    IF(r>rsplit) THEN
       !alpha 0.3-0.5 works well
       alpha=0.5
    ELSE
       !If alpha=1 this goes tits up
       !alpha 0.7-0.9 works well
       alpha=0.8
    END IF

    f0_rapid=r**alpha

  END FUNCTION f0_rapid

  FUNCTION sigma_integral1(r,a,cosm,acc)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigma_integral1
    REAL, INTENT(IN) :: r, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, INTENT(IN) :: acc
    REAL :: b, c
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: iorder=3

    b=r/(r+r**.5)
    c=1.

    IF(b==c) THEN

       !Fix the answer to zero if the integration limits are identical
       sigma_integral1=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(c-b)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=sigma_integrand_transformed(b,r,f1_rapid,a,cosm)
             f2=sigma_integrand_transformed(c,r,f1_rapid,a,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(b,c,i,n)
                fx=sigma_integrand_transformed(x,r,f1_rapid,a,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=REAL(sum_2n)
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGMA_INTEGRAL1: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             sigma_integral1=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             sigma_integral1=0.d0
             STOP 'SIGMA_INTEGRAL1: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sigma_integral1=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION sigma_integral1

  FUNCTION f1_rapid(r)

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error
    IMPLICIT NONE
    REAL :: f1_rapid
    REAL, INTENT(IN) :: r

    REAL, PARAMETER :: alpha=0.5

    f1_rapid=r**alpha

  END FUNCTION f1_rapid

  FUNCTION sigma_integral2(r,a,cosm,acc)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigma_integral2
    REAL, INTENT(IN) :: r, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL, INTENT(IN) :: acc
    REAL :: b, c
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    REAL, PARAMETER :: CC=10. !How far to go out in 1/r units for integral
    INTEGER, PARAMETER :: iorder=3

    b=1./r
    c=CC/r

    IF(b==c) THEN

       !Fix the answer to zero if the integration limits are identical
       sigma_integral2=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(c-b)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=sigma_integrand(b,r,a,cosm)
             f2=sigma_integrand(c,r,a,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(b,c,i,n)
                fx=sigma_integrand(x,r,a,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGMA_INTEGRAL2: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             sigma_integral2=REAL(sum_new)
             !WRITE(*,*) 'INTEGRATE_STORE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             sigma_integral2=0.d0
             STOP 'SIGMA_INTEGRAL2: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sigma_integral2=0.d0
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION sigma_integral2

  FUNCTION sigmaV(R,a,cosm)

    IMPLICIT NONE
    REAL :: sigmaV
    REAL, INTENT(IN) :: R, a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !IF(cosm%is_normalised .EQV. .FALSE.) CALL normalise_power(cosm)
    sigmaV=sigmaV_integral(R,a,cosm,acc_cosm)

  END FUNCTION sigmaV

  FUNCTION sigmaV_integral(R,a,cosm,acc)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigmaV_integral
    REAL, INTENT(IN) :: R, a, acc
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: b, c
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: iorder=3

    !Integration range for integration parameter
    !Note 0 -> infinity in k has changed to 0 -> 1 in x
    b=0.
    c=1.

    IF(b==c) THEN

       !Fix the answer to zero if the integration limits are identical
       sigmaV_integral=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.
       sum_n=0.
       sum_old=0.
       sum_new=0.

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(c-b)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=sigmaV_integrand(b,R,a,cosm)
             f2=sigmaV_integrand(c,R,a,cosm)
             sum_2n=0.5*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=progression(b,c,i,n)
                fx=sigmaV_integrand(x,R,a,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGMAV: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             sigmaV_integral=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             sigmaV_integral=0.
             STOP 'SIGMAV: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sigmaV_integral=0.
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.
          END IF

       END DO

    END IF

  END FUNCTION sigmaV_integral

  FUNCTION sigmaV_integrand(theta,R,a,cosm)

    !This is the integrand for the velocity dispersion integral
    USE special_functions
    IMPLICIT NONE
    REAL :: sigmaV_integrand
    REAL, INTENT(IN) :: theta, a, R
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: k

    REAL, PARAMETER :: alpha=1.65 !Speeds up integral for large 'R'
    REAL, PARAMETER :: Rsplit=10. !Value to impliment speed up

    !Note that I have not included the speed up alpha and Rsplit
    !The choice of alpha=1.65 seemed to work well for R=100.
    !Rsplit=10 is thoughlessly chosen (only because 100.>10.)
    !Including this seems to make things slower (faster integration but slower IF statements?)

    IF(theta==0. .OR. theta==1.) THEN
       sigmaV_integrand=0.
    ELSE
       !IF(r>Rsplit) THEN
       !   k=(-1.+1./theta)/r**alpha
       !ELSE
       k=(-1.+1./theta)
       sigmaV_integrand=(p_lin(k,a,cosm)/k**2)*(wk_tophat(k*r)**2)/(theta*(1.-theta))
    END IF

  END FUNCTION sigmaV_integrand

  FUNCTION grow(a,cosm)

    !Scale-independent growth function | normalised g(z=0)=1
    IMPLICIT NONE
    REAL :: grow
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)
    IF(a==1.) THEN
       grow=1.
    ELSE       
       grow=exp(find(log(a),cosm%log_a_growth,cosm%log_growth,cosm%n_growth,3,3,2))
    END IF

  END FUNCTION grow

  FUNCTION ungrow(a,cosm)

    !Unnormalised growth function
    IMPLICIT NONE
    REAL :: ungrow
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
   
    ungrow=cosm%gnorm*grow(a,cosm)

  END FUNCTION ungrow

  FUNCTION growth_rate(a,cosm)

    !Growth rate: dln(g) / dln(a)
    IMPLICIT NONE
    REAL :: growth_rate
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)    
    growth_rate=find(log(a),cosm%log_a_growth,cosm%growth_rate,cosm%n_growth,3,3,2)

  END FUNCTION growth_rate

  FUNCTION acc_growth(a,cosm)

    !Accumulated growth function: int_0^a g(a)/a da
    IMPLICIT NONE
    REAL :: acc_growth
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_growth .EQV. .FALSE.) CALL init_growth(cosm)    
    acc_growth=exp(find(log(a),cosm%log_a_growth,cosm%log_acc_growth,cosm%n_growth,3,3,2))

  END FUNCTION acc_growth

  FUNCTION growth_rate_Linder(a,cosm)

    !Approximation for the growth rate from Linder xxxx.xxxx
    IMPLICIT NONE
    REAL :: growth_rate_Linder
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: gam

    IF(cosm%w<-1.) THEN
       gam=0.55+0.02*(1.+cosm%w)
    ELSE IF(cosm%w>-1) THEN
       gam=0.55+0.05*(1.+cosm%w)
    ELSE
       gam=0.55
    END IF

    growth_rate_Linder=Omega_m_norad(a,cosm)**gam
    
  END FUNCTION growth_rate_Linder

  FUNCTION growth_Linder_integrand(a,cosm)

    !Integrand for the approximate growth integral using Linder approximate growth rate
    IMPLICIT NONE
    REAL :: growth_Linder_integrand
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    growth_Linder_integrand=growth_rate_Linder(a,cosm)/a

  END FUNCTION growth_Linder_integrand

  FUNCTION growth_Linder(a,cosm)

    !Calculate the growth function from the Linder growth rate via integration
    IMPLICIT NONE
    REAL :: growth_Linder
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    growth_Linder=exp(-integrate_cosm(a,1.,growth_Linder_integrand,cosm,acc_cosm,3))
    
  END FUNCTION growth_Linder

  FUNCTION growth_CPT(a,cosm)

    !Carrol, Press & Turner (1992) approximation to growth function (good to 5%)
    IMPLICIT NONE
    REAL :: growth_CPT
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: Om_mz, Om_vz, Om_m, Om_v

    !Get all necessary Omega values
    Om_mz=Omega_m_norad(a,cosm)
    Om_vz=Omega_v(a,cosm)+Omega_w(a,cosm)
    Om_m=cosm%Om_m
    Om_v=cosm%Om_v+cosm%Om_w

    !Now call CPT twice, second time to normalise it
    growth_CPT=CPT(a,Om_mz,Om_vz)/CPT(1.,Om_m,Om_v)

  END FUNCTION growth_CPT

  FUNCTION CPT(a,Om_m,Om_v)

    !The main CPT approximation from 1992
    IMPLICIT NONE
    REAL :: CPT
    REAL, INTENT(IN) :: a, Om_m, Om_v

    CPT=a*Om_m/((Om_m**(4./7.))-Om_v+(1.+Om_m/2.)*(1.+Om_v/70.))
    
  END FUNCTION CPT

  SUBROUTINE init_growth(cosm)

    USE calculus_table

    !Fills a table of the growth function vs. a
    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i, na
    REAL :: a
    REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL :: ainit, amax, dinit, vinit
    REAL :: g0, f0, bigG0

    INTEGER, PARAMETER :: n=128 !Number of entries for growth tables

    !TODO: Figure out why if I set amax=10, rather than amax=1, I start getting weird f(a) around a=0.001

    !The calculation should start at a z when Omega_m(a)=1, so that the assumption
    !of starting in the g\propto a growing mode is valid (this will not work for early DE)
    !Maximum a should be a=1. unless considering models in the future
    ainit=1e-3
    amax=1.

    !! First do growth factor and growth rate !!
    
    !These set the initial conditions to be the Omega_m(a)=1 growing mode
    dinit=ainit
    vinit=1.

    !Write some useful information to the screen
    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_GROWTH: Solving growth equation'
       WRITE(*,*) 'INIT_GROWTH: Minimum scale factor:', ainit
       WRITE(*,*) 'INIT_GROWTH: Maximum scale factor:', amax
       WRITE(*,*) 'INIT_GROWTH: Number of points for look-up tables:', n
    END IF

    !Solve the ODE
    CALL ODE_adaptive_cosmology(d_tab,v_tab,0.,a_tab,cosm,ainit,amax,dinit,vinit,fd,fv,acc_cosm,3,.FALSE.)
    IF(verbose_cosmology) WRITE(*,*) 'INIT_GROWTH: ODE done'
    na=SIZE(a_tab)

    !Convert dv/da to f = dlng/dlna for later, so v_tab should really be f_tab from now on
    v_tab=v_tab*a_tab/d_tab

    !Normalise so that g(z=0)=1
    cosm%gnorm=find(1.,a_tab,d_tab,na,3,3,2)
    IF(verbose_cosmology) WRITE(*,*) 'INIT_GROWTH: unnormalised growth at z=0:', REAL(cosm%gnorm)
    d_tab=d_tab/cosm%gnorm   

    !Allocate arrays
    IF(ALLOCATED(cosm%log_a_growth))    DEALLOCATE(cosm%log_a_growth)
    IF(ALLOCATED(cosm%log_growth))      DEALLOCATE(cosm%log_growth)
    IF(ALLOCATED(cosm%growth_rate))     DEALLOCATE(cosm%growth_rate)
    IF(ALLOCATED(cosm%log_acc_growth))  DEALLOCATE(cosm%log_acc_growth)
    cosm%n_growth=n

    !This downsamples the tables that come out of the ODE solver (which can be a bit long)
    !Could use some table-interpolation routine here to save time
    ALLOCATE(cosm%log_a_growth(n),cosm%log_growth(n),cosm%growth_rate(n))
    DO i=1,n
       a=progression(ainit,amax,i,n)
       cosm%log_a_growth(i)=a
       cosm%log_growth(i)=exp(find(log(a),log(a_tab),log(d_tab),na,3,3,2))
       cosm%growth_rate(i)=find(log(a),log(a_tab),v_tab,na,3,3,2)
    END DO

    !! !!
 
    !! Table integration to calculate G(a)=int_0^a g(a')/a' da' !!

    !Allocate array
    ALLOCATE(cosm%log_acc_growth(n))

    !Set to zero, because I have an x=x+y thing later on
    cosm%log_acc_growth=0.
    
    !Do the integral up to table position i, which fills the accumulated growth table
    DO i=1,n

       !Do the integral using the arrays
       IF(i>1) THEN
          cosm%log_acc_growth(i)=integrate_table(cosm%log_a_growth,cosm%gnorm*cosm%log_growth/cosm%log_a_growth,n,1,i,3)
       END IF
       
       !Then add on the section that is missing from the beginning
       !NB. g(a=0)/0 = 1, so you just add on a rectangle of height g*a/a=g
       cosm%log_acc_growth(i)=cosm%log_acc_growth(i)+cosm%gnorm*cosm%log_growth(1)
       
    END DO

    !! !!

    !Write stuff about growth parameter at a=1 to the screen    
    IF(verbose_cosmology) THEN
       f0=find(1.,cosm%log_a_growth,cosm%growth_rate,n,3,3,2)
       g0=find(1.,cosm%log_a_growth,cosm%log_growth,n,3,3,2)
       bigG0=find(1.,cosm%log_a_growth,cosm%log_acc_growth,n,3,3,2)      
       WRITE(*,*) 'INIT_GROWTH: normalised growth at z=0:', g0
       WRITE(*,*) 'INIT_GROWTH: growth rate at z=0:', f0
       WRITE(*,*) 'INIT_GROWTH: integrated growth at z=0:', bigG0
       WRITE(*,*)
    END IF

    !Make the some of the tables log for easier interpolation
    cosm%log_a_growth=log(cosm%log_a_growth)
    cosm%log_growth=log(cosm%log_growth)
    cosm%log_acc_growth=log(cosm%log_acc_growth)

    !Set the flag to true so that this subroutine is only called once
    cosm%has_growth=.TRUE.

    IF(verbose_cosmology) THEN
       WRITE(*,*) 'INIT_GROWTH: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE init_growth

  FUNCTION fd(d,v,k,a,cosm)

    !Needed for growth function solution
    !This is the fd in \dot{\delta}=fd
    IMPLICIT NONE
    REAL :: fd
    REAL, INTENT(IN) :: d, v, k, a
    REAL :: crap
    TYPE(cosmology), INTENT(INOUT) :: cosm

    !To prevent compile-time warnings
    crap=d
    crap=k
    crap=cosm%A
    crap=a    

    fd=v

  END FUNCTION fd

  FUNCTION fv(d,v,k,a,cosm)

    !Needed for growth function solution
    !This is the fv in \ddot{\delta}=fv
    IMPLICIT NONE
    REAL :: fv
    REAL, INTENT(IN) :: d, v, k, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: f1, f2
    REAL :: crap

    !To prevent compile-time warning
    crap=k

    f1=3.*Omega_m_norad(a,cosm)*d/(2.*(a**2))
    f2=(2.+AH_norad(a,cosm)/Hubble2_norad(a,cosm))*(v/a)
    fv=f1-f2

  END FUNCTION fv

  FUNCTION fvnl(d,v,k,a,cosm)

    !Function used for ODE solver in non-linear growth calculation
    IMPLICIT NONE
    REAL :: fvnl
    REAL, INTENT(IN) :: d, v, k, a
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: f1, f2, f3
    REAL :: crap

    !To prevent compile-time warning
    crap=k

    f1=3.*Omega_m_norad(a,cosm)*d*(1.+d)/(2.*(a**2))
    f2=-(2.+AH_norad(a,cosm)/Hubble2_norad(a,cosm))*(v/a)
    f3=4.*(v**2)/(3.*(1.+d))

    fvnl=f1+f2+f3

  END FUNCTION fvnl

  FUNCTION dc_NakamuraSuto(a,cosm)

    !Nakamura & Suto (1997) fitting formula for LCDM
    IMPLICIT NONE
    REAL :: dc_NakamuraSuto
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    dc_NakamuraSuto=dc0*(1.+0.0123*log10(Omega_m_norad(a,cosm)))

  END FUNCTION dc_NakamuraSuto

  FUNCTION Dv_BryanNorman(a,cosm)

    !Bryan & Norman (1998) spherical over-density fitting function
    !Here overdensity is defined relative to the background matter density, rather than the critical density
    IMPLICIT NONE
    REAL :: Dv_BryanNorman
    REAL :: x, Om_m
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    Om_m=Omega_m_norad(a,cosm)
    x=Om_m-1.

    IF(cosm%Om_v==0. .AND. cosm%Om_w==0.) THEN
       !Open model results
       Dv_BryanNorman=Dv0+60.*x-32.*x**2
       Dv_BryanNorman=Dv_BryanNorman/Om_m
    ELSE
       !LCDM results
       Dv_BryanNorman=Dv0+82.*x-39.*x**2
       Dv_BryanNorman=Dv_BryanNorman/Om_m
    END IF

  END FUNCTION Dv_BryanNorman

  FUNCTION dc_Mead(a,cosm)

    !delta_c fitting function from Mead (2017)
    IMPLICIT NONE
    REAL :: dc_Mead
    REAL, INTENT(IN) :: a !scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: lg, bG, Om_m

    !See Appendix A of Mead (2016) for naming convention
    REAL, PARAMETER :: p10=-0.0069
    REAL, PARAMETER :: p11=-0.0208
    REAL, PARAMETER :: p12=0.0312
    REAL, PARAMETER :: p13=0.0021
    INTEGER, PARAMETER :: a1=1
    REAL, PARAMETER :: p20=0.0001
    REAL, PARAMETER :: p21=-0.0647
    REAL, PARAMETER :: p22=-0.0417
    REAL, PARAMETER :: p23=0.0646
    INTEGER, PARAMETER :: a2=0

    lg=ungrow(a,cosm)
    bG=acc_growth(a,cosm)
    Om_m=Omega_m_norad(a,cosm)

    !WRITE(*,*) 'DC_STUFF:', a, lg, bG, Om_m

    dc_Mead=1.
    dc_Mead=dc_Mead+f_Mead(lg/a,bG/a,p10,p11,p12,p13)*log10(Om_m)**a1
    dc_Mead=dc_Mead+f_Mead(lg/a,bG/a,p20,p21,p22,p23)
    dc_Mead=dc_Mead*dc0

  END FUNCTION dc_Mead

  FUNCTION Dv_Mead(a,cosm)

    !Delta_v fitting function from Mead (2017)
    IMPLICIT NONE
    REAL :: Dv_Mead
    REAL, INTENT(IN) :: a !scale factor
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: lg, bG, Om_m

    !See Appendix A of Mead (2017) for naming convention
    REAL, PARAMETER :: p30=-0.79
    REAL, PARAMETER :: p31=-10.17
    REAL, PARAMETER :: p32=2.51
    REAL, PARAMETER :: p33=6.51
    INTEGER, PARAMETER :: a3=1
    REAL, PARAMETER :: p40=-1.89
    REAL, PARAMETER :: p41=0.38
    REAL, PARAMETER :: p42=18.8
    REAL, PARAMETER :: p43=-15.87
    INTEGER, PARAMETER :: a4=2

    lg=ungrow(a,cosm)
    bG=acc_growth(a,cosm)
    Om_m=Omega_m_norad(a,cosm)

    Dv_Mead=1.
    Dv_Mead=Dv_Mead+f_Mead(lg/a,bG/a,p30,p31,p32,p33)*log10(Om_m)**a3
    Dv_Mead=Dv_Mead+f_Mead(lg/a,bG/a,p40,p41,p42,p43)*log10(Om_m)**a4
    Dv_Mead=Dv_Mead*Dv0

  END FUNCTION Dv_Mead

  PURE FUNCTION f_Mead(x,y,p0,p1,p2,p3)

    !Equation A3 in Mead (2017)
    IMPLICIT NONE
    REAL :: f_Mead
    REAL, INTENT(IN) :: x, y
    REAL, INTENT(IN) :: p0, p1, p2, p3

    f_Mead=p0+p1*(1.-x)+p2*(1.-x)**2+p3*(1.-y)

  END FUNCTION f_Mead

  REAL FUNCTION dc_spherical(a,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)
    
    IF(log(a)<cosm%log_a_dcDv(1)) THEN
       dc_spherical=dc0
    ELSE
       dc_spherical=find(log(a),cosm%log_a_dcDv,cosm%dc,cosm%n_dcDv,3,3,2)
    END IF

  END FUNCTION dc_spherical

  REAL FUNCTION Dv_spherical(a,cosm)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm

    IF(cosm%has_spherical .EQV. .FALSE.) CALL init_spherical_collapse(cosm)
    
    IF(log(a)<cosm%log_a_dcDv(1)) THEN
       Dv_spherical=Dv0
    ELSE
       Dv_spherical=find(log(a),cosm%log_a_dcDv,cosm%Dv,cosm%n_dcDv,3,3,2)
    END IF

  END FUNCTION Dv_spherical

  SUBROUTINE init_spherical_collapse(cosm)

    IMPLICIT NONE
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: dinit, ainit, vinit, ac, dc
    REAL :: av, a_rmax, d_rmax, Dv, rmax, rv
    REAL, ALLOCATABLE :: d(:), a(:), v(:)
    REAL, ALLOCATABLE :: dnl(:), vnl(:), rnl(:)
    REAL, ALLOCATABLE :: a_coll(:), r_coll(:)
    INTEGER :: i, j, k, k2   

    REAL, PARAMETER :: amax=2. !Maximum scale factor to consider
    REAL, PARAMETER :: dmin=1e-7 !Minimum starting value for perturbation
    REAL, PARAMETER :: dmax=1e-3 !Maximum starting value for perturbation
    INTEGER, PARAMETER :: m=128 !Number of collapse scale-factors to try to calculate (you usually get fewer)
    INTEGER, PARAMETER :: n=100000 !Number of points for ODE calculations (needs to be large (~1e5) to capture final stages of collapse

    IF(verbose_cosmology) WRITE(*,*) 'SPHERICAL_COLLAPSE: Doing integration'

    !Allocate arrays
    IF(ALLOCATED(cosm%log_a_dcDv)) DEALLOCATE(cosm%log_a_dcDv)
    IF(ALLOCATED(cosm%dc))         DEALLOCATE(cosm%dc)
    IF(ALLOCATED(cosm%Dv))         DEALLOCATE(cosm%Dv)
    ALLOCATE(cosm%log_a_dcDv(m),cosm%dc(m),cosm%Dv(m))
    cosm%log_a_dcDv=0.
    cosm%dc=0.
    cosm%Dv=0.

    IF(verbose_cosmology) THEN
       WRITE(*,*) 'SPHERICAL_COLLAPSE: delta min', dmin
       WRITE(*,*) 'SPHERICAL_COLLAPSE: delta max', dmax
       WRITE(*,*) 'SPHERICAL_COLLAPSE: number of collapse points attempted', m
    END IF

    !BCs for integration. Note ainit=dinit means that collapse should occur around a=1 for dmin
    !amax should be slightly greater than 1 to ensure at least a few points for a>0.9 (i.e not to miss out a=1)    
    ainit=dmin
    vinit=1.*(dmin/ainit) !vinit=1 is EdS growing mode solution

    !Now loop over all initial density fluctuations
    DO j=1,m       

       !log range of initial delta
       dinit=progression_log(dmin,dmax,j,m)

       !Do both with the same a1 and a2 and using the same number of time steps
       !This means that arrays a, and anl will be identical, which simplifies calculation
       CALL ODE_spherical(dnl,vnl,0.,a,cosm,ainit,amax,dinit,vinit,fd,fvnl,n,3,.TRUE.)
       DEALLOCATE(a)
       CALL ODE_spherical(d,v,0.,a,cosm,ainit,amax,dinit,vinit,fd,fv,n,3,.TRUE.)

       !If this condition is met then collapse occured some time a<amax
       IF(dnl(n)==0.) THEN

          !! delta_c calcualtion !!

          ALLOCATE(rnl(n))

          rnl=a*(1.+dnl)**(-1./3.)

          !Find the collapse point (very crude)
          !More accurate calculations seem to be worse
          !I think this is due to the fact that delta spikes very quickly
          DO i=1,n
             IF(dnl(i)==0.) THEN
                !k is the new maxium size of the arrays
                k=i-1
                EXIT
             END IF
          END DO

          !Cut away parts of the arrays for a>ac
          CALL amputate(a,n,k)
          CALL amputate(d,n,k)
          CALL amputate(dnl,n,k)
          CALL amputate(rnl,n,k)

          !Collapse has occured so use previous a as ac and d as dc
          ac=a(k)
          dc=d(k)

          !! !!

          !! Now to Delta_v calculation !!

          !Find the a values when the perturbation is maximum size
          a_rmax=maximum(a,rnl,k)

          !Find the over-density at this point
          d_rmax=exp(find(log(a_rmax),log(a),log(dnl),SIZE(a),1,3,2))

          !Find the maximum radius
          rmax=find(log(a_rmax),log(a),rnl,SIZE(a),1,3,2)

          !The radius of the perturbation when it is virialised is half maximum
          !This might not be appropriate for LCDM models (or anything with DE)
          rv=rmax/2.

          !Need to assign new arrays for the collapse branch of r such that it is monotonic
          k2=int_split(d_rmax,dnl,k)

          !Allocate collapse branch arrays
          ALLOCATE(a_coll(k-k2+1),r_coll(k-k2+1))

          !Fill collapse branch arrays
          DO i=k2,k
             a_coll(i-k2+1)=a(i)
             r_coll(i-k2+1)=rnl(i)
          END DO

          !Find the scale factor when the perturbation has reached virial radius
          av=exp(find(rv,r_coll,log(a_coll),SIZE(r_coll),3,3,2))

          !Deallocate collapse branch arrays
          DEALLOCATE(a_coll,r_coll)

          !Spherical model approximation is that perturbation is at virial radius when
          !'collapse' is considered to have occured, which has already been calculated
          Dv=exp(find(log(av),log(a),log(dnl),SIZE(a),1,3,2))*(ac/av)**3.
          Dv=Dv+1.

          !!

          cosm%log_a_dcDv(j)=ac
          cosm%dc(j)=dc
          cosm%Dv(j)=Dv

          DEALLOCATE(rnl)

       END IF

       !Deallocate arrays ready for next calculation
       DEALLOCATE(d,v,a)
       DEALLOCATE(dnl,vnl)

    END DO

    IF(verbose_cosmology) WRITE(*,*) 'SPHERICAL COLLAPSE: calculation complete'

    !Reverse the arrays so that they run lowest a to highest a
    CALL reverse(cosm%log_a_dcDv,m)
    CALL reverse(cosm%dc,m)
    CALL reverse(cosm%Dv,m)

    !IF(verbose_cosmology) THEN
    !   WRITE(*,*) '===================================='
    !   WRITE(*,*) 'Point  scalefactor  delta_c  Delta_v'
    !   WRITE(*,*) '===================================='
    !   DO i=1,m
    !      IF(cosm%a_dcDv(i)==0.) EXIT
    !      WRITE(*,fmt='(I5,F13.4,F9.4,F9.1)') i, cosm%a_dcDv(i), cosm%dc(i), cosm%Dv(i)
    !   END DO
    !   WRITE(*,*) '===================================='
    !END IF

    !Calculate the maximum sizes for these new arrays
    DO i=1,m
       IF(cosm%log_a_dcDv(i)==0.) EXIT
    END DO
    cosm%n_dcDv=i-1

    IF(verbose_cosmology) THEN
       WRITE(*,*) 'SPHERICAL_COLLAPSE: number of collapse points:', cosm%n_dcDv
       WRITE(*,*)
    END IF

    !Remove bits of the array that are unnecessary
    CALL amputate(cosm%log_a_dcDv,m,cosm%n_dcDv)
    CALL amputate(cosm%dc,m,cosm%n_dcDv)
    CALL amputate(cosm%Dv,m,cosm%n_dcDv)

    !Take a logarithm
    cosm%log_a_dcDv=log(cosm%log_a_dcDv)

    !Set the flag
    cosm%has_spherical=.TRUE.

  END SUBROUTINE init_spherical_collapse

  SUBROUTINE ODE_spherical(x,v,kk,t,cosm,ti,tf,xi,vi,fx,fv,n,imeth,ilog)

    !Solves 2nd order ODE x''(t) from ti to tf and creates arrays of x, v, t values
    !I have sometimes called this ODE_crass
    !It has a fixed number of time steps, n
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:), v(:), t(:)
    REAL, INTENT(IN) :: kk, xi, vi, ti, tf
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(IN) :: imeth, n
    LOGICAL, INTENT(IN) :: ilog     
    DOUBLE PRECISION, ALLOCATABLE :: x8(:), v8(:), t8(:)
    INTEGER :: i

    !imeth sets ODE solving method
    !imeth = 1: Crude method
    !imeth = 2: Mid-point method
    !imeth = 3: Runge-Kutta

    INTERFACE

       !fx is what x' is equal to
       FUNCTION fx(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL :: fx
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fx

       !fv is what v' is equal to
       FUNCTION fv(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL :: fv
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fv

    END INTERFACE

    !Allocate arrays
    ALLOCATE(x8(n),v8(n),t8(n))

    !Need to be set to zero for this to work in the spherical-collapse case
    x8=0.d0
    v8=0.d0
    t8=0.d0

    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    x8(1)=xi
    v8(1)=vi

    !Fill time array
    IF(ilog) THEN
       CALL fill_array_double(log(ti),log(tf),t8,n)
       t8=exp(t8)
    ELSE
       CALL fill_array_double(ti,tf,t8,n)
    END IF

    DO i=1,n-1
       
       CALL ODE_advance_cosmology(x8(i),x8(i+1),v8(i),v8(i+1),t8(i),t8(i+1),fx,fv,imeth,kk,cosm)
       
       !Needed to escape from the ODE solver when the perturbation is ~collapsed
       IF(x8(i+1)>1e8) EXIT
       
    END DO

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(v)) DEALLOCATE(v)
    IF(ALLOCATED(t)) DEALLOCATE(t)
    ALLOCATE(x(n),v(n),t(n))
    x=REAL(x8)
    v=REAL(v8)
    t=REAL(t8)

    !WRITE(*,*) 'ODE: Integration complete in steps:', n

  END SUBROUTINE ODE_spherical

  SUBROUTINE ODE_adaptive_cosmology(x,v,kk,t,cosm,ti,tf,xi,vi,fx,fv,acc,imeth,ilog)

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    !acc is the desired accuracy across the entire solution
    !time steps are increased until convergence is achieved
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: x(:), t(:), v(:)
    REAL, INTENT(IN) :: kk, xi, vi, ti, tf, acc
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(IN) :: imeth
    LOGICAL, INTENT(IN) :: ilog
    DOUBLE PRECISION, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)      
    INTEGER :: i, j, n, k, np, ifail, kn

    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=100

    !imeth sets ODE solving method
    !imeth = 1: Crude method
    !imeth = 2: Mid-point method
    !imeth = 3: Runge-Kutta   

    INTERFACE

       !fx is what x' is equal to
       FUNCTION fx(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL :: fx
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fx

       !fv is what v' is equal to
       FUNCTION fv(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL :: fv
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fv

    END INTERFACE

    DO j=1,jmax

       n=1+ninit*(2**(j-1))

       ALLOCATE(x8(n),v8(n),t8(n))

       x8=0.d0
       v8=0.d0
       t8=0.d0

       !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
       x8(1)=xi
       v8(1)=vi

       !Fill time array
       IF(ilog) THEN
          CALL fill_array_double(log(ti),log(tf),t8,n)
          t8=exp(t8)
       ELSE
          CALL fill_array_double(ti,tf,t8,n)
       END IF

       ifail=0

       DO i=1,n-1
          CALL ODE_advance_cosmology(x8(i),x8(i+1),v8(i),v8(i+1),t8(i),t8(i+1),fx,fv,imeth,kk,cosm)
       END DO

       IF(j==1) ifail=1

       IF(j .NE. 1) THEN

          np=1+(n-1)/2

          DO k=1,1+(n-1)/2

             kn=2*k-1

             IF(ifail==0) THEN

                IF(xh(k)>acc .AND. x8(kn)>acc .AND. (ABS(xh(k)/x8(kn))-1.)>acc) ifail=1
                IF(vh(k)>acc .AND. v8(kn)>acc .AND. (ABS(vh(k)/v8(kn))-1.)>acc) ifail=1

                IF(ifail==1) THEN
                   DEALLOCATE(xh,th,vh)
                   EXIT
                END IF

             END IF
          END DO

       END IF

       IF(ifail==0) THEN
          IF(ALLOCATED(x)) DEALLOCATE(x)
          IF(ALLOCATED(v)) DEALLOCATE(v)
          IF(ALLOCATED(t)) DEALLOCATE(t)
          ALLOCATE(x(n),v(n),t(n))
          x=REAL(x8)
          v=REAL(v8)
          t=REAL(t8)
          EXIT
       END IF

       ALLOCATE(xh(n),th(n),vh(n))
       xh=x8
       vh=v8
       th=t8
       DEALLOCATE(x8,t8,v8)

    END DO

  END SUBROUTINE ODE_adaptive_cosmology

  SUBROUTINE ODE_advance_cosmology(x1,x2,v1,v2,t1,t2,fx,fv,imeth,k,cosm)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x1, v1, t1, t2
    DOUBLE PRECISION, INTENT(OUT) :: x2, v2
    INTEGER, INTENT(IN) :: imeth
    REAL, INTENT(IN) :: k
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: x, v, t, dt
    REAL :: kx1, kx2, kx3, kx4
    REAL :: kv1, kv2, kv3, kv4

    INTERFACE

       !fx is what x' is equal to
       FUNCTION fx(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL :: fx
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fx

       !fv is what v' is equal to
       FUNCTION fv(x,v,k,t,cosm)
         IMPORT :: cosmology
         REAL :: fv
         REAL, INTENT(IN) :: x, v, k, t
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION fv

    END INTERFACE

    x=REAL(x1)
    v=REAL(v1)
    t=REAL(t1)

    dt=REAL(t2-t1)

    IF(imeth==1) THEN

       !Crude method!
       kx1=dt*fx(x,v,k,t,cosm)
       kv1=dt*fv(x,v,k,t,cosm)

       x2=x1+kx1
       v2=v1+kv1

    ELSE IF(imeth==2) THEN

       !Mid-point method!
       kx1=dt*fx(x,v,k,t,cosm)
       kv1=dt*fv(x,v,k,t,cosm)
       kx2=dt*fx(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)
       kv2=dt*fv(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)

       x2=x1+kx2
       v2=v1+kv2

    ELSE IF(imeth==3) THEN

       !RK4 (Holy Christ, this is so fast compared to above methods)!
       kx1=dt*fx(x,v,k,t,cosm)
       kv1=dt*fv(x,v,k,t,cosm)
       kx2=dt*fx(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)
       kv2=dt*fv(x+kx1/2.,v+kv1/2.,k,t+dt/2.,cosm)
       kx3=dt*fx(x+kx2/2.,v+kv2/2.,k,t+dt/2.,cosm)
       kv3=dt*fv(x+kx2/2.,v+kv2/2.,k,t+dt/2.,cosm)
       kx4=dt*fx(x+kx3,v+kv3,k,t+dt,cosm)
       kv4=dt*fv(x+kx3,v+kv3,k,t+dt,cosm)

       x2=x1+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.d0
       v2=v1+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.d0

    ELSE

       STOP 'ODE_ADVANCE: Error, imeth specified incorrectly'

    END IF

  END SUBROUTINE ODE_advance_cosmology

   FUNCTION integrate_cosm(a,b,f,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate_cosm
    REAL, INTENT(IN) :: a, b, acc
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER, INTENT(IN) :: iorder
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old

    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       REAL FUNCTION f(a,cosm)
         IMPORT :: cosmology
         REAL, INTENT(IN) :: a
         TYPE(cosmology), INTENT(INOUT) :: cosm
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate_cosm=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax

          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN

             !The first go is just the trapezium of the end points
             f1=f(a,cosm)
             f2=f(b,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n

          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             !integrate=REAL(sum_new)
             !WRITE(*,*) 'INTEGRATE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate_cosm=REAL(sum_new)

    END IF

  END FUNCTION integrate_cosm

!!$  FUNCTION integrate_time(a,b,cosm,acc)
!!$
!!$    !Integrates between a and b until desired accuracy is reached
!!$    !Stores information to reduce function calls
!!$    IMPLICIT NONE
!!$    REAL :: integrate_time
!!$    REAL, INTENT(IN) :: a, b, acc
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    INTEGER :: i, j
!!$    INTEGER :: n
!!$    REAL :: x, dx
!!$    REAL :: f1, f2, fx
!!$    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
!!$
!!$    INTEGER, PARAMETER :: jmin=5
!!$    INTEGER, PARAMETER :: jmax=30
!!$    INTEGER, PARAMETER :: iorder=3
!!$
!!$    IF(a==b) THEN
!!$
!!$       !Fix the answer to zero if the integration limits are identical
!!$       integrate_time=0.
!!$
!!$    ELSE
!!$
!!$       !Set the sum variable for the integration
!!$       sum_2n=0.
!!$       sum_n=0.
!!$       sum_old=0.
!!$       sum_new=0.
!!$
!!$       DO j=1,jmax
!!$
!!$          !Note, you need this to be 1+2**n for some integer n
!!$          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
!!$          n=1+2**(j-1)
!!$
!!$          !Calculate the dx interval for this value of 'n'
!!$          dx=(b-a)/REAL(n-1)
!!$
!!$          IF(j==1) THEN
!!$
!!$             !The first go is just the trapezium of the end points
!!$             f1=time_integrand(a,cosm)
!!$             f2=time_integrand(b,cosm)
!!$             sum_2n=0.5*(f1+f2)*dx
!!$             sum_new=sum_2n
!!$
!!$          ELSE
!!$
!!$             !Loop over only new even points to add these to the integral
!!$             DO i=2,n,2
!!$                x=progression(a,b,i,n)
!!$                fx=time_integrand(x,cosm)
!!$                sum_2n=sum_2n+fx
!!$             END DO
!!$
!!$             !Now create the total using the old and new parts
!!$             sum_2n=sum_n/2.+sum_2n*dx
!!$
!!$             !Now calculate the new sum depending on the integration order
!!$             IF(iorder==1) THEN  
!!$                sum_new=sum_2n
!!$             ELSE IF(iorder==3) THEN         
!!$                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
!!$             ELSE
!!$                STOP 'INTEGRATE_TIME: Error, iorder specified incorrectly'
!!$             END IF
!!$
!!$          END IF
!!$
!!$          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
!!$             !jmin avoids spurious early convergence
!!$             integrate_time=REAL(sum_new)
!!$             EXIT
!!$          ELSE IF(j==jmax) THEN
!!$             integrate_time=0.d0
!!$             STOP 'INTEGRATE_TIME: Integration timed out'
!!$          ELSE
!!$             !Integral has not converged so store old sums and reset sum variables
!!$             integrate_time=0.d0
!!$             sum_old=sum_new
!!$             sum_n=sum_2n
!!$             sum_2n=0.
!!$          END IF
!!$
!!$       END DO
!!$
!!$    END IF
!!$
!!$  END FUNCTION integrate_time

!!$  FUNCTION integrate_distance(a,b,cosm,acc)
!!$
!!$    !Integrates between a and b until desired accuracy is reached
!!$    !Stores information to reduce function calls
!!$    IMPLICIT NONE
!!$    REAL :: integrate_distance
!!$    REAL, INTENT(IN) :: a, b, acc
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    INTEGER :: i, j
!!$    INTEGER :: n
!!$    REAL :: x, dx
!!$    REAL :: f1, f2, fx
!!$    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
!!$
!!$    INTEGER, PARAMETER :: jmin=5
!!$    INTEGER, PARAMETER :: jmax=30
!!$    INTEGER, PARAMETER :: iorder=3
!!$
!!$    IF(a==b) THEN
!!$
!!$       !Fix the answer to zero if the integration limits are identical
!!$       integrate_distance=0.
!!$
!!$    ELSE
!!$
!!$       !Set the sum variable for the integration
!!$       sum_2n=0.
!!$       sum_n=0.
!!$       sum_old=0.
!!$       sum_new=0.
!!$
!!$       DO j=1,jmax
!!$
!!$          !Note, you need this to be 1+2**n for some integer n
!!$          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
!!$          n=1+2**(j-1)
!!$
!!$          !Calculate the dx interval for this value of 'n'
!!$          dx=(b-a)/REAL(n-1)
!!$
!!$          IF(j==1) THEN
!!$
!!$             !The first go is just the trapezium of the end points
!!$             f1=distance_integrand(a,cosm)
!!$             f2=distance_integrand(b,cosm)
!!$             sum_2n=0.5*(f1+f2)*dx
!!$             sum_new=sum_2n
!!$
!!$          ELSE
!!$
!!$             !Loop over only new even points to add these to the integral
!!$             DO i=2,n,2
!!$                x=progression(a,b,i,n)
!!$                fx=distance_integrand(x,cosm)
!!$                sum_2n=sum_2n+fx
!!$             END DO
!!$
!!$             !Now create the total using the old and new parts
!!$             sum_2n=sum_n/2.+sum_2n*dx
!!$
!!$             !Now calculate the new sum depending on the integration order
!!$             IF(iorder==1) THEN  
!!$                sum_new=sum_2n
!!$             ELSE IF(iorder==3) THEN         
!!$                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
!!$             ELSE
!!$                STOP 'INTEGRATE_DISTANCE: Error, iorder specified incorrectly'
!!$             END IF
!!$
!!$          END IF
!!$
!!$          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
!!$             !jmin avoids spurious early convergence
!!$             integrate_distance=REAL(sum_new)
!!$             EXIT
!!$          ELSE IF(j==jmax) THEN
!!$             integrate_distance=0.d0
!!$             STOP 'INTEGRATE_DISTANCE: Integration timed out'
!!$          ELSE
!!$             !Integral has not converged so store old sums and reset sum variables
!!$             integrate_distance=0.d0
!!$             sum_old=sum_new
!!$             sum_n=sum_2n
!!$             sum_2n=0.
!!$          END IF
!!$
!!$       END DO
!!$
!!$    END IF
!!$
!!$  END FUNCTION integrate_distance

!!$  FUNCTION growint(a,cosm,acc)
!!$
!!$    !Integrates between a and b until desired accuracy is reached
!!$    !Stores information to reduce function calls
!!$    IMPLICIT NONE
!!$    REAL :: growint
!!$    REAL, INTENT(IN) :: a
!!$    REAL, INTENT(IN) :: acc
!!$    TYPE(cosmology), INTENT(INOUT) :: cosm
!!$    REAL :: b
!!$    INTEGER :: i, j
!!$    INTEGER :: n
!!$    REAL :: x, dx
!!$    REAL :: f1, f2, fx
!!$    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
!!$
!!$    INTEGER, PARAMETER :: jmin=5
!!$    INTEGER, PARAMETER :: jmax=30
!!$    INTEGER, PARAMETER :: iorder=3   
!!$
!!$    !Integration range for integration parameter
!!$    !Note a -> 1
!!$    b=1.
!!$
!!$    IF(a==b) THEN
!!$
!!$       !Fix the answer to zero if the integration limits are identical
!!$       growint=exp(0.)
!!$
!!$    ELSE
!!$
!!$       !Reset the sum variable for the integration
!!$       sum_2n=0.
!!$       sum_n=0.
!!$       sum_old=0.
!!$       sum_new=0.
!!$
!!$       DO j=1,jmax
!!$
!!$          !Note, you need this to be 1+2**n for some integer n
!!$          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
!!$          n=1+2**(j-1)
!!$
!!$          !Calculate the dx interval for this value of 'n'
!!$          dx=(b-a)/REAL(n-1)
!!$
!!$          IF(j==1) THEN
!!$
!!$             !The first go is just the trapezium of the end points
!!$             f1=growint_integrand(a,cosm)
!!$             f2=growint_integrand(b,cosm)
!!$             sum_2n=0.5*(f1+f2)*dx
!!$             sum_new=sum_2n
!!$
!!$          ELSE
!!$
!!$             !Loop over only new even points to add these to the integral
!!$             DO i=2,n,2
!!$                !x=a+(b-a)*REAL(i-1)/REAL(n-1)
!!$                x=progression(a,b,i,n)
!!$                fx=growint_integrand(x,cosm)
!!$                sum_2n=sum_2n+fx
!!$             END DO
!!$
!!$             !Now create the total using the old and new parts
!!$             sum_2n=sum_n/2.+sum_2n*dx
!!$
!!$             !Now calculate the new sum depending on the integration order
!!$             IF(iorder==1) THEN  
!!$                sum_new=sum_2n
!!$             ELSE IF(iorder==3) THEN         
!!$                sum_new=(4.*sum_2n-sum_n)/3. !This is Simpson's rule and cancels error
!!$             ELSE
!!$                STOP 'GROWINT: Error, iorder specified incorrectly'
!!$             END IF
!!$
!!$          END IF
!!$
!!$          IF((j>=jmin) .AND. (ABS(-1.+sum_new/sum_old)<acc)) THEN
!!$             !jmin avoids spurious early convergence
!!$             growint=exp(REAL(sum_new))
!!$             EXIT
!!$          ELSE IF(j==jmax) THEN
!!$             growint=0.d0
!!$             STOP 'GROWINT: Integration timed out'
!!$          ELSE
!!$             !Integral has not converged so store old sums and reset sum variables
!!$             growint=0.d0
!!$             sum_old=sum_new
!!$             sum_n=sum_2n
!!$             sum_2n=0.
!!$          END IF
!!$
!!$       END DO
!!$
!!$    END IF
!!$
!!$  END FUNCTION growint

END MODULE cosmology_functions
