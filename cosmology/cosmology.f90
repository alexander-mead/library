MODULE cosmology

  TYPE cosmo
     !Contains only things that do not need to be recalculated with each new z
     REAL :: Om_m, Om_v, Om_w, Om_r, Om !Standard Omegas
     REAL :: Om_mp, Om_vp !Matter and vacuum perturbation variables     
     REAL :: a1, a2, a1n, a2n, astar, nstar, wstar, Om_ws !IDE models
     REAL :: w0, wa, wm, am, dm !QUICC models
     REAL :: mua
     REAL :: H0rc, c2, c3
     REAL :: omeg, phi0
     INTEGER :: iw, ihub, img !Switches
     !DOUBLE PRECISION :: h, n, sig8, w, wa, om_nu, Om_c, Om_b
     !DOUBLE PRECISION :: om, k, z_cmb, om_r, T_cmb
     !DOUBLE PRECISION :: A
     !DOUBLE PRECISION, ALLOCATABLE :: sigma(:), r_sigma(:)
     !DOUBLE PRECISION, ALLOCATABLE :: logsigma(:), logr_sigma(:)
     !DOUBLE PRECISION, ALLOCATABLE :: growth(:), a_growth(:)
     !DOUBLE PRECISION, ALLOCATABLE :: r(:), a_r(:)
     !INTEGER :: nsig, ng, nr
     !CHARACTER(len=256) :: name
  END TYPE cosmo

CONTAINS

  FUNCTION Hubble2(a,cosm)

    !This calculates the dimensionless squared hubble parameter at redshift z!
    !and it ignores contributions from radiation (not accurate at high z)!
    IMPLICIT NONE
    REAL :: Hubble2
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    IF(cosm%ihub==1) THEN
       Hubble2=cosm%Om_r*(a**(-4))+cosm%Om_m*(a**(-3))+cosm%Om_w*X_de(a,cosm)+cosm%Om_v+(1.-cosm%Om)*(a**(-2))
    ELSE IF(cosm%ihub==2) THEN
       !Assumes flatness I think, and no radiation or dark energy
       IF(a<1.e-5) THEN
          !To avoid floating exception with -6 power in full expression
          Hubble2=cosm%Om_m*a**(-3)
       ELSE
          Hubble2=0.5*(cosm%Om_m*(a**(-3))+sqrt((cosm%Om_m**2)*a**(-6)+4.*(1.-cosm%Om_m)))
       END IF
    ELSE
       STOP 'HUBBLE2: Error, Hubble parameter paramter specified incorrectly'
    END IF

  END FUNCTION Hubble2

  FUNCTION AH(a,cosm)

    !\ddot a/a
    IMPLICIT NONE
    REAL :: AH
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm
    REAL :: f1, f2, f3

    IF(cosm%ihub==1) THEN
       AH=2.*cosm%Om_r*a**(-4)+cosm%Om_m*a**(-3)+cosm%Om_w*(1.+3.*w_de(a,cosm))*X_de(a,cosm)-2.*cosm%Om_v
       AH=-AH/2.
    ELSE IF(cosm%ihub==2) THEN
       IF(a<1e-5) THEN
          !To avoid floating expception with -6 power in full expression
          AH=-0.5*cosm%Om_m*a**(-3)
       ELSE
          f1=cosm%Om_m*a**(-3)
          f2=(cosm%Om_m**2)*a**(-6)-8.*(1.-cosm%Om_m)
          f3=sqrt((cosm%Om_m**2)*a**(-6)+4.*(1.-cosm%Om_m))
          AH=-0.25*(f1+f2/f3)
       END IF
    ELSE
       STOP 'AH: Error, AH parameter parameter specified incorrectly'
    END IF

  END FUNCTION AH

  FUNCTION Omega_m(a,cosm)

    !This calculates Omega_m variations with z!
    IMPLICIT NONE
    REAL :: Omega_m
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    Omega_m=cosm%Om_m*(a**(-3))/Hubble2(a,cosm)

  END FUNCTION Omega_m

  FUNCTION Omega_r(a,cosm)

    !This calculates Omega_r variations with z!
    IMPLICIT NONE
    REAL :: Omega_r
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    Omega_r=cosm%Om_r*(a**(-4))/Hubble2(a,cosm)

  END FUNCTION Omega_r

  FUNCTION Omega_v(a,cosm)

    !This calculates Omega_v variations with z!
    IMPLICIT NONE
    REAL :: Omega_v
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    Omega_v=cosm%Om_v/Hubble2(a,cosm)

  END FUNCTION Omega_v

  FUNCTION Omega_w(a,cosm)

    !This calculates Omega_w variations with z!
    IMPLICIT NONE
    REAL :: Omega_w
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    Omega_w=cosm%Om_w*X_de(a,cosm)/Hubble2(a,cosm)

  END FUNCTION Omega_w

  FUNCTION Omega(a,cosm)

    !This calculates total Omega variations with z!
    IMPLICIT NONE
    REAL :: Omega
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    Omega=Omega_r(a,cosm)+Omega_m(a,cosm)+Omega_v(a,cosm)+Omega_w(a,cosm)

  END FUNCTION Omega

  FUNCTION w_de(a,cosm)

    !Variations of the dark energy parameter w(a)
    IMPLICIT NONE
    REAL :: w_de
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm
    REAL :: p1, p2, p3, p4
    DOUBLE PRECISION :: f1, f2, f3, f4

    IF(cosm%iw==0) THEN
       !LCDM
       w_de=-1.
    ELSE IF(cosm%iw==1) THEN
       !QUICC parameterisation
       p1=1.+exp(cosm%am/cosm%dm)
       p2=1.-exp(-(a-1.)/cosm%dm)
       p3=1.+exp(-(a-cosm%am)/cosm%dm)
       p4=1.-exp(1./cosm%dm)
       w_de=cosm%w0+(cosm%wm-cosm%w0)*p1*p2/(p3*p4)
    ELSE IF(cosm%iw==2) THEN
       !w(a)CDM
       w_de=cosm%w0+(1.-a)*cosm%wa
    ELSE IF(cosm%iw==3) THEN
       !wCDM
       w_de=cosm%w0
    ELSE IF(cosm%iw==4) THEN
       !IDE I
       w_de=((a/cosm%astar)**cosm%nstar-1.)/((a/cosm%astar)**cosm%nstar+1.)
    ELSE IF(cosm%iw==5) THEN
       !IDE II
       f1=a**cosm%nstar-cosm%a1n
       f2=a**cosm%nstar+cosm%a1n
       f3=a**cosm%nstar-cosm%a2n
       f4=a**cosm%nstar+cosm%a2n
       w_de=-1.+REAL(f1/f2-f3/f4)
    ELSE IF(cosm%iw==6) THEN
       !IDE III
       IF(a<cosm%a1) THEN
          w_de=-1.
       ELSE IF(cosm%a1<=a .AND. a<cosm%a2) THEN
          w_de=cosm%wstar
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
    TYPE(cosmo), INTENT(IN) :: cosm

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
    TYPE(cosmo), INTENT(IN) :: cosm

    w_eff=w_de(a,cosm)*Omega_w(a,cosm)-Omega_v(a,cosm)+Omega_r(a,cosm)/3.
    w_eff=w_eff/Omega(a,cosm)

  END FUNCTION w_eff

  FUNCTION X_de(a,cosm)

    !USE calculus
    
    !Redshift scaling for dark energy (i.e., if w=0 x(a)=a^-3, if w=-1 x(a)=const etc.)
    IMPLICIT NONE
    REAL :: X_de
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm
    DOUBLE PRECISION :: f1, f2, f3, f4
    REAL, PARAMETER :: acc=1e-3

    IF(cosm%iw==0) THEN
       !LCDM
       X_de=1.   
    ELSE IF(cosm%iw==2) THEN
       !w(a)CDM
       X_de=(a**(-3.*(1.+cosm%w0+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))
    ELSE IF(cosm%iw==3) THEN
       !wCDM
       X_de=a**(-3.*(1.+cosm%w0))
    ELSE IF(cosm%iw==4) THEN
       !IDE I
       X_de=((1.+(a/cosm%astar)**cosm%nstar)/(1.+(1./cosm%astar)**cosm%nstar))**(-6./cosm%nstar)
    ELSE IF(cosm%iw==5) THEN
       !IDE II
       f1=a**cosm%nstar+cosm%a1n
       f2=1.+cosm%a1n
       f3=1.+cosm%a2n
       f4=a**cosm%nstar+cosm%a2n
       X_de=REAL(f1*f3/(f2*f4))**(-6./cosm%nstar)
    ELSE IF(cosm%iw==6) THEN
       !IDE III
       IF(a<cosm%a1) THEN
          X_de=(cosm%a1/cosm%a2)**(-3.*(1.+cosm%wstar))
       ELSE IF(cosm%a1<=a .AND. a<cosm%a2) THEN
          X_de=(a/cosm%a2)**(-3.*(1.+cosm%wstar))
       ELSE IF(a>=cosm%a2) THEN
          X_de=1.
       ELSE
          STOP 'X_DE: Error, something went wrong'
       END IF
    ELSE
       STOP 'X_DE: Error, iw not specified correctly'
!!$    ELSE
!!$       !Generally true, doing this integration can make calculations very slow
!!$       !Difficult to implement into library because of cosm dependence of w_de !!!
!!$       !X_de=(a**(-3))*exp(3.*integrate_log(a,1.,integrand_de,acc,3,1))
    END IF

  END FUNCTION X_de

!!$  FUNCTION integrand_de(a)!,cosm)
!!$
!!$    !The integrand for the X_de(a) integral
!!$    IMPLICIT NONE
!!$    REAL :: integrand_de
!!$    REAL, INTENT(IN) :: a
!!$    !TYPE(cosmo), INTENT(IN) :: cosm
!!$
!!$    integrand_de=w_de(a,cosm)/a
!!$
!!$  END FUNCTION integrand_de

  FUNCTION phibar(a,cosm)

    !JBD phibar approximation
    IMPLICIT NONE
    REAL :: phibar
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    phibar=cosm%phi0*a**(1./(1.+cosm%omeg))
    
  END FUNCTION phibar

  FUNCTION beta_DGP(a,cosm)

    !DGP beta function
    IMPLICIT NONE
    REAL :: beta_DGP
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm
    
    !H0rc updated to be dimensionless parameter
    
    beta_DGP=1.+(2./3.)*sqrt(Hubble2(a,cosm))*cosm%H0rc*(2.+AH(a,cosm)/Hubble2(a,cosm))

  END FUNCTION beta_DGP

  FUNCTION beta1(a,cosm)

    IMPLICIT NONE
    REAL :: beta1
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm
    REAL :: f1, f2, f3
    
    f1=cosm%c2
    f2=4.*cosm%c3*(phiddot(a,cosm)+2.*sqrt(Hubble2(a,cosm))*phidot(a,cosm))
    f3=2.*(cosm%c3**2)*phidot(a,cosm)**4

    beta1=(-f1-f2+f3)/(6.*cosm%c3)
    
  END FUNCTION beta1

  FUNCTION beta2(a,cosm)

    IMPLICIT NONE
    REAL :: beta2
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    beta2=2.*beta1(a,cosm)/(phidot(a,cosm)**2.)

    !beta2=2.*H2(a)*beta1(a)/(6.*(1.-om_m))
    
  END FUNCTION beta2

  FUNCTION phidot(a,cosm)

    IMPLICIT NONE
    REAL :: phidot
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm

    phidot=sqrt(6.*(1.-cosm%Om_m))/sqrt(Hubble2(a,cosm))
    
  END FUNCTION phidot

  FUNCTION phiddot(a,cosm)

    IMPLICIT NONE
    REAL :: phiddot
    REAL, INTENT(IN) :: a
    TYPE(cosmo), INTENT(IN) :: cosm
    
    phiddot=sqrt(6.*(1.-cosm%Om_m))*(1.-AH(a,cosm)/Hubble2(a,cosm))
    
  END FUNCTION phiddot

  SUBROUTINE assign_cosmology(cosm)

    !Set cosmology here
    IMPLICIT NONE
    TYPE(cosmo), INTENT(OUT) :: cosm
    INTEGER :: imod
    !REAL :: f1, f2, Om_ws, Xstar

    cosm%Om_m=1.0
    cosm%Om_w=0.0
    cosm%Om_v=0.0
    cosm%w0=-1.
    cosm%wa=0.
    cosm%wm=-1.
    cosm%am=0.
    cosm%dm=1.
    cosm%mua=0.
    cosm%iw=0
    cosm%img=0
    cosm%ihub=1
    
    WRITE(*,*)
    WRITE(*,*) 'Welcome to perturbation integrator'
    WRITE(*,*)
    WRITE(*,*) ' 1 - EdS'
    WRITE(*,*) ' 2 - Flat LCDM'
    WRITE(*,*) ' 3 - w(a)CDM'
    WRITE(*,*) ' 4 - INV1'
    WRITE(*,*) ' 5 - INV2'
    WRITE(*,*) ' 6 - SUGRA'
    WRITE(*,*) ' 7 - 2EXP'
    WRITE(*,*) ' 8 - AS'
    WRITE(*,*) ' 9 - CNR'
    WRITE(*,*) '10 - General LCDM (can be non-flat)'
    WRITE(*,*) '11 - Millennium (WMAP1ish)'
    WRITE(*,*) '12 - MG, constant mu'
    WRITE(*,*) '13 - Simpson mu parametrisation'
    WRITE(*,*) '14 - wCDM (constant w)'
    WRITE(*,*) '15 - SLICS-LE (WMAP9ish)'
    WRITE(*,*) '16 - Flat DGP'
    WRITE(*,*) '17 - DGP (Om_m=1.)'
    WRITE(*,*) '18 - IDE I'
    WRITE(*,*) '19 - IDE II'
    WRITE(*,*) '20 - IDE III'
    WRITE(*,*) '21 - QCDM - Cubic Galileon'
    WRITE(*,*) '22 - Cubic Galileon'
    WRITE(*,*) '23 - Background/perturbation Omegas'
    WRITE(*,*) '24 - Jordan-Brans-Dicke cosmology'
    WRITE(*,*) '25 - Thawing quintessence'
    READ(*,*) imod
    WRITE(*,*)

    IF(imod==1) THEN
       !EdS
    ELSE IF(imod==2) THEN
       !Flat LCDM
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       cosm%Om_v=1.-cosm%Om_m
    ELSE IF(imod==3) THEN
       !w(a)CDM
       cosm%iw=2
       WRITE(*,*) 'w(a)=w0+wa*(1.-a)'
       WRITE(*,*) 'w0:'
       READ(*,*) cosm%w0
       WRITE(*,*) 'wa:'
       READ(*,*) cosm%wa
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Om_w:'
       READ(*,*) cosm%Om_w
    ELSE IF(imod==4) THEN
       !QUICC - INV1
       cosm%Om_m=0.26
       cosm%Om_w=0.74
       cosm%iw=1
       cosm%w0=-0.4
       cosm%wm=-0.27
       cosm%am=0.18
       cosm%dm=0.5
    ELSE IF(imod==5) THEN
       !QUICC - INV2
       cosm%Om_m=0.26
       cosm%Om_w=0.74
       cosm%iw=1
       cosm%w0=-0.79
       cosm%wm=-0.67
       cosm%am=0.29
       cosm%dm=0.4
    ELSE IF(imod==6) THEN
       !QUICC - SUGRA
       cosm%Om_m=0.26
       cosm%Om_w=0.74
       cosm%iw=1
       cosm%w0=-0.82
       cosm%wm=-0.18
       cosm%am=0.1
       cosm%dm=0.7
    ELSE IF(imod==7) THEN
       !QUICC - 2EXP
       cosm%Om_m=0.26
       cosm%Om_w=0.74
       cosm%iw=1
       cosm%w0=-1.
       cosm%wm=0.01
       cosm%am=0.19
       cosm%dm=0.043
    ELSE IF(imod==8) THEN
       !QUICC - AS
       cosm%Om_m=0.26
       cosm%Om_w=0.74
       cosm%iw=1
       cosm%w0=-0.96
       cosm%wm=-0.01
       cosm%am=0.53
       cosm%dm=0.13
    ELSE IF(imod==9) THEN
       !QUICC - CNR
       cosm%Om_m=0.26
       cosm%Om_w=0.74
       cosm%iw=1
       cosm%w0=-1.
       cosm%wm=0.1
       cosm%am=0.15
       cosm%dm=0.016
    ELSE IF(imod==10) THEN
       !General LCDM (non-flat)
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Om_v:'
       READ(*,*) cosm%Om_v
    ELSE IF(imod==11) THEN
       !Millennium (WMAP1ish)
       cosm%Om_m=0.25
       cosm%Om_v=0.75
    ELSE IF(imod==12) THEN
       !Silly MG parameterisation
       cosm%img=1
       WRITE(*,*) 'Om_m=1.'
       WRITE(*,*) 'MG: G -> (1.+mu*a)G, mu is constant'
       WRITE(*,*) 'mu:'
       READ(*,*) cosm%mua   
    ELSE IF(imod==13) THEN
       !Another silly MG parameterisation
       cosm%img=2
       WRITE(*,*) 'MG: G -> (1.+mu)G, mu=mu0*Om_v(a)/Om_v'
       WRITE(*,*) 'mu0:'
       READ(*,*) cosm%mua
    ELSE IF(imod==14) THEN
       !wCDM (constant w)
       cosm%iw=3
       WRITE(*,*) 'Omega_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Omega_w:'
       READ(*,*) cosm%Om_w
       WRITE(*,*) 'w (constant):'
       READ(*,*) cosm%w0
    ELSE IF(imod==15) THEN
       !SLICS-LE (WMAP9ish)
       cosm%Om_m=0.2905
       cosm%Om_v=0.7095
    ELSE IF(imod==16) THEN
       !Flat DGP
       cosm%img=3
       WRITE(*,*) 'DGP'
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*)
       cosm%Om_v=1.-cosm%Om_m
       WRITE(*,*) 'H0*rc [dimensionless]:'
       READ(*,*) cosm%H0rc
    ELSE IF(imod==17) THEN
       !DGP with EdS background
       cosm%img=3
       cosm%Om_m=1.
       cosm%Om_w=0.
       WRITE(*,*) 'DGP (with Om_m=1.)'
       WRITE(*,*) 'H0*rc [dimensionless]:'
       READ(*,*) cosm%H0rc
    ELSE IF(imod==18) THEN
       !IDE I
       cosm%iw=4
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Om_v:'
       READ(*,*) cosm%Om_v
       !WRITE(*,*) 'Om_w:'
       !READ(*,*) Om_w
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%astar
       WRITE(*,*) 'n:'
       READ(*,*) cosm%nstar
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_w
    ELSE IF(imod==19) THEN
       !IDE II model
       cosm%iw=5
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Om_w:'
       READ(*,*) cosm%Om_w
       cosm%Om_v=0. !No vacuum necessary herea
       WRITE(*,*) 'n:'
       READ(*,*) cosm%nstar
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%astar
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       !nstar=5
       !astar=0.05
       !Om_ws=0.2
    ELSE IF(imod==20) THEN
       !IDE III model
       cosm%iw=6
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Om_w:'
       READ(*,*) cosm%Om_w
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%astar
       WRITE(*,*) 'Om_w(a*):'
       READ(*,*) cosm%Om_ws
       WRITE(*,*) 'w*:'
       READ(*,*) cosm%wstar
       cosm%Om_v=0. !No vacuum necessary here
       !a1=0.1
       !a2=0.2
       !Om_m=0.3
       !Om_w=0.7
       !Om_v=0.
       !wstar=1
       !Om_ws=0.2
       !astar=0.1
    ELSE IF(imod==21 .OR. imod==22) THEN
       !QCDM and Galileon
       cosm%ihub=2
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       cosm%Om_v=1.-cosm%Om_m
       IF(imod==22) THEN
          cosm%img=4
          cosm%c2=-1.
          cosm%c3=1./(6.*sqrt(6.*(1.-cosm%Om_m)))
       END IF
    ELSE IF(imod==23) THEN
       !Different background/perturbation Omegas
       cosm%img=5
       WRITE(*,*) 'Background Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'Perturbation Om_m:'
       READ(*,*) cosm%Om_mp
       cosm%Om_v=1.-cosm%Om_m
       cosm%Om_vp=1.-cosm%Om_mp
    ELSE IF(imod==24) THEN
       !JBD
       cosm%img=6
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       WRITE(*,*) 'JBD omega:'
       READ(*,*) cosm%omeg
       WRITE(*,*)
       cosm%Om_v=1.-cosm%Om_m
       cosm%phi0=1.
    ELSE IF(imod==25) THEN
       !Thawing quintessence (this is the same as IDE I, just different parameters user specified)
       cosm%iw=4
       WRITE(*,*) 'Om_m:'
       READ(*,*) cosm%Om_m
       !WRITE(*,*) 'Om_v:'
       !READ(*,*) Om_v
       WRITE(*,*) 'Om_w:'
       READ(*,*) cosm%Om_w
       WRITE(*,*) 'a*:'
       READ(*,*) cosm%astar
       WRITE(*,*) 'n:'
       READ(*,*) cosm%nstar
    ELSE
       STOP 'Model not specified correctly'
    END IF
    WRITE(*,*)
    
  END SUBROUTINE assign_cosmology

  SUBROUTINE derive_cosmology(cosm)

    !Calculate derived cosmological parameters
    IMPLICIT NONE
    TYPE(cosmo), INTENT(INOUT) :: cosm
    REAL :: Xstar, f1, f2

    !Curvature Omega
    cosm%Om=cosm%Om_m+cosm%Om_w+cosm%Om_v+cosm%Om_r

    !Dark energy models
    IF(cosm%iw==4) THEN
       !WRITE(*,*) H2(astar), X(astar)
       !Om_w=Om_w*(Om_m*astar**(-3)+Om_v)/(X(astar)*(1.-Om_w))
       cosm%Om_w=cosm%Om_w*(Hubble2(cosm%astar,cosm)-cosm%Om_w*X_de(cosm%astar,cosm)+cosm%Om_w*cosm%astar**(-2))/(X_de(cosm%astar,cosm)*(1.-cosm%Om_w)+cosm%Om_w*cosm%astar**(-2))
    ELSE IF(cosm%iw==5) THEN
       !Define a1^n
       cosm%a1n=cosm%astar**cosm%nstar

       !Necessary for first step below
       cosm%a2n=cosm%a1n

       !All neccessary to convert parameters to a1,a2
       f1=cosm%Om_ws*(Hubble2(cosm%astar,cosm)-cosm%Om_w*X_de(cosm%astar,cosm))
       f2=cosm%Om_w*(1.-cosm%Om_ws)
       Xstar=f1/f2
       !a1=astar**nstar
       !Xstar=Xstar**(-nstar/6.)
       Xstar=Xstar**(cosm%nstar/6.)
            
       !Top and bottom of fraction
       !f1=2.-a1*Xstar*(1.+1./a1)
       !f2=Xstar*(1.+1./a1)-2.
       f1=cosm%a1n*(2.*Xstar-(1.+cosm%a1n))
       f2=(1.+cosm%a1n)-2.*Xstar*cosm%a1n
       
       !Finally! a2
       cosm%a2n=f1/f2
       !IF(a2<a1) a2=a1
    ELSE IF(cosm%iw==6) THEN
       cosm%a1=cosm%astar !Scale-factor at which Om_w(a*) is most important
       cosm%a2=cosm%astar !Needs to be set for X(a*) and H2(a*) below (which cancel each other)
       f1=cosm%Om_ws*(Hubble2(cosm%astar,cosm)-cosm%Om_w*X_de(cosm%astar,cosm))
       f2=cosm%Om_w*(1.-cosm%Om_ws)
       cosm%a2=cosm%astar*(f1/f2)**(1./(3.*(1.+cosm%wstar)))
    END IF
    
  END SUBROUTINE derive_cosmology

  SUBROUTINE write_cosmology(cosm)

    IMPLICIT NONE
    TYPE(cosmo), INTENT(IN) :: cosm

    !Standard Omegas
    WRITE(*,*) 'Parameters:'
    WRITE(*,*) 'Om_m:', cosm%Om_m
    WRITE(*,*) 'Om_w:', cosm%Om_w
    WRITE(*,*) 'Om_v:', cosm%Om_v
    WRITE(*,*) 'Om_r:', cosm%Om_r    
    WRITE(*,*) 'Om:', cosm%Om
    WRITE(*,*) 'Om_k:', 1.-cosm%Om

    !Dark energy
    IF(cosm%iw==0) THEN
       WRITE(*,*) 'Vacuum energy'
       WRITE(*,*) 'w:', -1
    ELSE IF(cosm%iw==1) THEN
       WRITE(*,*) 'QUICC dark energy prescription'
       WRITE(*,*) 'w0:', cosm%w0
       WRITE(*,*) 'wm:', cosm%wm
       WRITE(*,*) 'am:', cosm%am
       WRITE(*,*) 'dm:', cosm%dm
    ELSE IF(cosm%iw==2) THEN
       WRITE(*,*) 'w(a) = w0+wa(1.-a)'
       WRITE(*,*) 'w0:', cosm%w0
       WRITE(*,*) 'wa:', cosm%wa
    ELSE IF(cosm%iw==3) THEN
       WRITE(*,*) 'Constant w'
       WRITE(*,*) 'w:', cosm%w0
    ELSE IF(cosm%iw==4) THEN
       WRITE(*,*) 'IDE I'
       WRITE(*,*) 'a*:', cosm%astar
       WRITE(*,*) 'n:', cosm%nstar
    ELSE IF(cosm%iw==5) THEN
       WRITE(*,*) 'IDE II'
       WRITE(*,*) 'a*:', cosm%astar
       WRITE(*,*) 'Om_w(a*):', cosm%Om_ws
       WRITE(*,*) 'n*:', cosm%nstar
       WRITE(*,*) 'a1^n (derived):', cosm%a1n
       WRITE(*,*) 'a2^n (derived):', cosm%a2n
    ELSE IF(cosm%iw==6) THEN
       WRITE(*,*) 'IDE III'
       WRITE(*,*) 'a*:', cosm%a1
       WRITE(*,*) 'Om_w(a*):', cosm%Om_ws
       WRITE(*,*) 'w*:', cosm%wstar
    END IF
    
    !Modified gravity
    IF(cosm%img==3) THEN
       WRITE(*,*) 'nDGP gravity'
       WRITE(*,*) 'Beta(z=0):', beta_DGP(1.,cosm)
    ELSE IF(cosm%img==5) THEN
       WRITE(*,*) 'Different Omegas gravity'
       WRITE(*,*) 'Perturbation Om_m:', cosm%Om_mp
       WRITE(*,*) 'Perturbation Om_v:', cosm%Om_vp
    END IF
    WRITE(*,*)
    
  END SUBROUTINE write_cosmology
  
END MODULE cosmology
