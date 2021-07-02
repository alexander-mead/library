MODULE HOD_functions

   USE basic_operations
   USE array_operations
   USE special_functions
   USE random_numbers

   IMPLICIT NONE

   PRIVATE

   ! Type
   PUBLIC :: hodmod

   ! Theory
   PUBLIC :: assign_HOD
   PUBLIC :: mean_centrals
   PUBLIC :: mean_satellites
   PUBLIC :: mean_galaxies
   PUBLIC :: variance_centrals
   PUBLIC :: variance_satellites
   PUBLIC :: covariance_centrals_satellites
   PUBLIC :: variance_galaxies

   ! Realisation
   !PUBLIC :: random_number_of_centrals
   !PUBLIC :: random_number_of_satellites
   PUBLIC :: random_number_of_galaxies

   ! HOD models
   PUBLIC :: ihod_toy
   PUBLIC :: ihod_toy_int
   PUBLIC :: ihod_toy_noscatter
   PUBLIC :: ihod_Zehavi
   PUBLIC :: ihod_Zheng
   PUBLIC :: ihod_Zhai

   ! Type for HOD
   TYPE hodmod
      INTEGER :: ihod
      REAL :: Mmin, Mmax
      REAL :: M0, M1, Mcen, Mcut, sigma, alpha
      INTEGER :: stats_cen, stats_sat
   END TYPE hodmod

   ! Statistical models
   INTEGER, PARAMETER :: stats_delta = 0           ! Delta function statistics
   INTEGER, PARAMETER :: stats_Bernoulli = 1       ! Bernoulli statistics (central galaxies if only 0 or 1 possible)
   INTEGER, PARAMETER :: stats_Poisson = 2         ! Poisson statistics (usual for satellite galaxies)
   INTEGER, PARAMETER :: stats_Poisson_central = 3 ! Modified Poisson via central condition

   ! Defaults
   REAL, PARAMETER :: Mmin_def = 1e7
   REAL, PARAMETER :: Mmax_def = 1e17
   REAL, PARAMETER :: Mcen_def = 1e13
   REAL, PARAMETER :: sigma_def = 0.5
   REAL, PARAMETER :: M0_def = Mcen_def
   REAL, PARAMETER :: M1_def = Mcen_def
   REAL, PARAMETER :: alpha_def = 1.
   REAL, PARAMETER :: Mcut_def = 1e14

   ! Satellite occupation
   REAL, PARAMETER :: eps_sat = 1e-6

   ! HOD models
   INTEGER, PARAMETER :: ihod_Zehavi = 1        ! Zehavi et al. (2005; https://arxiv.org/abs/astro-ph/0408569)
   INTEGER, PARAMETER :: ihod_Zheng = 2         ! Zheng et al. (2005; https://arxiv.org/abs/astro-ph/0408564)
   INTEGER, PARAMETER :: ihod_toy = 3           ! Toy model
   INTEGER, PARAMETER :: ihod_toy_int = 4       ! Toy model with integer numbers of galaxies
   INTEGER, PARAMETER :: ihod_toy_noscatter = 5 ! Toy model with integer numbers of galaxies and no scatter
   INTEGER, PARAMETER :: ihod_Zhai = 6          ! Zhai et al. (2017; https://arxiv.org/abs/1607.05383)

   CONTAINS

   SUBROUTINE assign_HOD(ihod, hod)

      ! Assign a halo-occupation model
      INTEGER, INTENT(IN) :: ihod      ! Integer to choose model
      TYPE(hodmod), INTENT(OUT) :: hod ! HOD model

      ! Set HOD model
      IF (ihod == -1) THEN
         WRITE(*, *) '1 - Zehavi et al. (2005)'
         WRITE(*, *) '2 - Zheng et al. (2005)'
         WRITE(*, *) '3 - Toy HOD'
         WRITE(*, *) '4 - Toy HOD with no scatter'
         WRITE(*, *) '5 - Toy HOD'
         WRITE(*, *) '6 - Zhai et al. (2017)'
         READ(*, *) hod%ihod
      ELSE
         hod%ihod = ihod
      END IF

      ! Set the statistical models for the distribution of central and satellite galaxies
      hod%stats_cen = stats_Bernoulli
      hod%stats_sat = stats_Poisson_central

      ! Minimum and maximum halo masses to host a galaxy
      hod%Mmin = Mmin_def
      hod%Mmax = Mmax_def

      ! Default HOD parameters
      hod%Mcen = Mcen_def
      hod%sigma = sigma_def
      hod%M0 = M0_def
      hod%M1 = M1_def
      hod%alpha = alpha_def
      hod%Mcut = Mcut_def

      ! Set specific HOD model parameters
      IF (hod%ihod == ihod_toy_noscatter) hod%stats_sat = stats_delta

   END SUBROUTINE assign_HOD

   REAL FUNCTION mean_centrals(M, hod)

      ! Mean number of central galaxies in a halo
      ! Note that this is the mean, not necessarily 0 or 1, but should (probably) be between 0 and 1
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model

      IF (between(M, hod%Mmin, hod%Mmax)) THEN
         IF (is_in_array(hod%ihod, [ihod_Zheng, ihod_Zhai])) THEN
            IF (hod%sigma == 0.) THEN
               mean_centrals = Heaviside(M-hod%Mcen, 1.)
            ELSE
               mean_centrals = 0.5*(1.+erf(log10(M/hod%Mcen)/hod%sigma)) ! Weird combination of erf with log10
            END IF
         ELSE IF (is_in_array(hod%ihod, [ihod_Zehavi, ihod_toy, ihod_toy_int, ihod_toy_noscatter])) THEN
            mean_centrals = Heaviside(M-hod%Mcen, 1.) ! mean_centrals is only ever 0 or 1 in these HOD models
         ELSE
            STOP 'MEAN_CENTRALS: Error, HOD not recognised'
         END IF
      ELSE
         mean_centrals = 0.
      END IF

   END FUNCTION mean_centrals

   REAL FUNCTION mean_satellites(M, hod)

      ! Mean number of satellite galaxies in a halo
      ! Note that this is not an integer in general
      ! Note that this may be reduced from the expectation from the HOD if the central condition is imposed
      REAL, INTENT(IN) :: M            ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod  ! HOD model
      REAL, PARAMETER :: eps = eps_sat ! To stop negative values: TODO: Necessary?

      IF (between(M, hod%Mmin, hod%Mmax)) THEN
         IF (hod%ihod == ihod_Zehavi) THEN
            mean_satellites = (M/hod%M1)**hod%alpha
         ELSE IF (hod%ihod == ihod_Zheng) THEN
            ! No satellite galaxies if M < M0
            ! It is annoying that the denominator is not M1-M0
            mean_satellites = Heaviside(M-hod%M0)*((M-hod%M0)/hod%M1)**hod%alpha
         ELSE IF (hod%ihod == ihod_Zhai) THEN
            mean_satellites = ((M/hod%M1)**hod%alpha)*exp(-hod%Mcut/M)
         ELSE IF (is_in_array(hod%ihod, [ihod_toy, ihod_toy_int, ihod_toy_noscatter])) THEN
            mean_satellites = M/hod%M1
            IF (is_in_array(hod%ihod, [ihod_toy_int, ihod_toy_noscatter])) mean_satellites = floor(mean_satellites)
         ELSE
            STOP 'MEAN_SATELLITES: Error, HOD not recognised'
         END IF
      ELSE
         mean_satellites = 0.
      END IF

      ! Note very carefully this in this case the mean you get from the standard HOD is *not* the actual mean
      ! The actual mean is always less; this is very important
      IF (hod%stats_sat == stats_Poisson_central) mean_satellites = mean_satellites*mean_centrals(M, hod)

   END FUNCTION mean_satellites

   REAL FUNCTION mean_galaxies(M, hod)

      ! Mean number of galaxies in a halo
      ! Calculated as the sum of the means of central and satellite galaxies
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model

      mean_galaxies = mean_centrals(M, hod)+mean_satellites(M, hod)

   END FUNCTION mean_galaxies

   REAL FUNCTION variance_centrals(M, hod)

      ! Variance in the number of central galaxies in a halo
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      REAL :: p

      p = mean_centrals(M, hod)
      IF (p .NE. 0.) THEN
         IF (hod%stats_cen == stats_delta) THEN
            variance_centrals = 0.
         ELSE IF (hod%stats_cen == stats_Bernoulli) THEN
            ! Variance for a Bernoulli distribution is p-p^2
            ! If p = 1 or p = 0 then the variance is automatically 0, as expected  
            variance_centrals = p*(1.-p)
         ELSE
            STOP 'VARIANCE_CENTRALS: Error, central statistics not recognised'
         END IF
      ELSE
         variance_centrals = 0.
      END IF

   END FUNCTION variance_centrals

   REAL FUNCTION variance_satellites(M, hod)

      ! Variance in the number of satellite galaxies in a halo
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      REAL :: p, lam, lam_old

      lam = mean_satellites(M, hod)
      IF (lam .NE. 0.) THEN
         IF (hod%stats_sat == stats_delta) THEN
            variance_satellites = 0.
         ELSE IF (hod%stats_sat == stats_Poisson) THEN
            variance_satellites = lam ! Variance for a Poisson distribution is equal to the mean
         ELSE IF (hod%stats_sat == stats_Poisson_central) THEN
            p = mean_centrals(M, hod)
            IF (p == 0.) THEN
               variance_satellites = 0.
            ELSE
               lam_old = lam/p ! Mean of underlying Poisson distribution
               variance_satellites = p*lam_old*(1.+lam_old*(1.-p)) ! See notes; if Nc=1 then this is Poisson result
            END IF
         ELSE
            STOP 'VARIANCE_SATELLITES: Error, satellite statistics not recognised'
         END IF
      ELSE
         variance_satellites = 0.
      END IF

   END FUNCTION variance_satellites

   REAL FUNCTION covariance_centrals_satellites(M, hod)

      ! Covariance between the number of central and satellite galaxies in a halo
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      REAL :: p, lam, lam_old

      IF (hod%stats_sat == stats_Poisson_central) THEN
         p = mean_centrals(M, hod)
         lam = mean_satellites(M, hod)
         IF (p /= 0.) THEN
            lam_old = lam/p ! Mean of underlying Poisson distribution
            covariance_centrals_satellites = lam_old*(1.-p)
         ELSE
            covariance_centrals_satellites = 0.
         END IF
      ELSE
         covariance_centrals_satellites = 0.
      END IF

   END FUNCTION covariance_centrals_satellites

   REAL FUNCTION variance_galaxies(M, hod)

      ! Variance in the total number of galaxies in a halo
      ! Ignores any possible covariance between the number of centrals and number of satellites
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model

      variance_galaxies = variance_centrals(M, hod)+variance_satellites(M, hod)
      variance_galaxies = variance_galaxies+2.*covariance_centrals_satellites(M, hod)

   END FUNCTION variance_galaxies

   SUBROUTINE random_number_of_galaxies(Nc, Ns, M, hod)

      ! Draw random numbers of central and satellite galaxies
      ! This must be done jointly if the central condition is to be imposed
      INTEGER, INTENT(OUT) :: Nc, Ns  ! Numbers of central and satellite galaxies
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      REAL :: p, lam, lam_old

      p = mean_centrals(M, hod)
      IF (hod%stats_cen == stats_delta) THEN
         Nc = nint(p)
      ELSE IF (hod%stats_cen == stats_Bernoulli) THEN
         Nc = random_Bernoulli(p)
      ELSE
         STOP 'RANDOM_NUMBER_OF_GALAXIES: Error, central statistics not recognised'
      END IF

      lam = mean_satellites(M, hod)
      IF (hod%stats_sat == stats_delta) THEN
         Ns = nint(lam)
      ELSE IF (hod%stats_sat == stats_Poisson) THEN
         Ns = random_Poisson(lam)
      ELSE IF (hod%stats_sat == stats_Poisson_central) THEN
         IF (p == 0) THEN
            Ns = 0 ! Central condition
         ELSE
            lam_old = lam/p ! Mean of the underlying Poisson distribution
            Ns = random_Poisson(lam_old) ! Use underlying mean here, which is *not* the actual mean of the distribution
         END IF
      ELSE
         STOP 'RANDOM_NUMBER_OF_GALAXIES: Error, satellite statistics not recognised'
      END IF

   END SUBROUTINE random_number_of_galaxies

   ! INTEGER FUNCTION random_number_of_centrals(M, hod)

   !    ! Draw from the random distribution to get a number of central galaxies in a halo
   !    ! This should probably be either 0 or 1
   !    REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
   !    TYPE(hodmod), INTENT(IN) :: hod ! HOD model
   !    REAL :: Nc, Ns

   !    !mean = mean_centrals(M, hod)
   !    !IF (hod%stats_cen == stats_delta) THEN
   !    !   random_number_of_centrals = nint(mean)
   !    !ELSE IF (hod%stats_cen == stats_Bernoulli) THEN
   !    !   random_number_of_centrals = random_Bernoulli(mean)
   !    !ELSE
   !    !   STOP 'RANDOM_NUMBER_OF_CENTRALS: Error, central statistics not recognised'
   !    !END IF
   !    CALL draw_random_number_of_galaxies(Nc, Ns, M, hod)
   !    random_number_of_satellites = Nc

   ! END FUNCTION random_number_of_centrals

   ! INTEGER FUNCTION random_number_of_satellites(M, hod)

   !    ! Draw from the random distribution to get a number of satellite galaxies in a halo
   !    REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
   !    TYPE(hodmod), INTENT(IN) :: hod ! HOD model
   !    REAL :: Nc, Ns

   !    !mean = mean_satellites(M, hod)
   !    !IF (hod%stats_sat == stats_delta) THEN
   !    !   random_number_of_satellites = nint(mean)
   !    !ELSE IF (hod%stats_sat == stats_Poisson) THEN
   !    !   random_number_of_satellites = random_Poisson(mean)
   !    !ELSE IF (hod%stats_sat == stats_Poisson_central) THEN
   !    !   STOP 'RANDOM_NUMBER_OF_SATELLITES: Error, if the central condition is imposed you must jointly draw central and satellite galaxies'
   !    !ELSE
   !    !   STOP 'RANDOM_NUMBER_OF_SATELLITES: Error, satellite statistics not recognised'
   !    !END IF
   !    CALL draw_random_number_of_galaxies(Nc, Ns, M, hod)
   !    random_number_of_satellites = Ns

   ! END FUNCTION random_number_of_satellites

   ! INTEGER FUNCTION random_number_of_galaxies(M, hod)

   !    ! Draw from the random distribution to get a total number of galaxies in a halo
   !    REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
   !    TYPE(hodmod), INTENT(IN) :: hod ! HOD model
   !    INTEGER :: Nc, Ns

   !    !Nc = random_number_of_centrals(M, hod)
   !    !Ns = random_number_of_satellites(M, hod)
   !    CALL draw_random_number_of_galaxies(Nc, Ns, M, hod)
   !    random_number_of_galaxies = Nc+Ns

   ! END FUNCTION random_number_of_galaxies

END MODULE HOD_functions