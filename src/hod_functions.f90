MODULE HOD_functions

   USE basic_operations
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
   PUBLIC :: variance_galaxies

   ! Realisation
   PUBLIC :: random_number_of_centrals
   PUBLIC :: random_number_of_satellites
   PUBLIC :: random_number_of_galaxies

   ! Type for HOD
   TYPE hodmod
      INTEGER :: ihod
      REAL :: Mmin, Mmax
      REAL :: M0, M1, sigma, alpha
      INTEGER :: stats_cen, stats_sat
   END TYPE hodmod

   ! Statistical models
   INTEGER, PARAMETER :: stats_delta = 0
   INTEGER, PARAMETER :: stats_Bernoulli = 1
   INTEGER, PARAMETER :: stats_Poisson = 2

   REAL, PARAMETER :: Mmin_def = 1e12
   REAL, PARAMETER :: Mmax_def = 1e18

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
         READ(*, *) hod%ihod
      ELSE
         hod%ihod = ihod
      END IF

      ! Set the statistical models for the distribution of central and satellite galaxies
      hod%stats_cen = stats_Bernoulli
      hod%stats_sat = stats_Poisson

      ! Minimum and maximum halo masses to host a galaxy
      hod%Mmin = Mmin_def
      hod%Mmax = Mmax_def

      ! Set HOD model parameters
      IF (hod%ihod == 1) THEN
         ! Zehavi et al. (2005; https://arxiv.org/abs/astro-ph/0408569)
         hod%M1 = hod%Mmin
         hod%alpha = 1.
      ELSE IF (hod%ihod == 2) THEN
         ! Zheng et al. (2005; https://arxiv.org/abs/astro-ph/0408564)
         hod%sigma = 0.1
         hod%M0 = hod%Mmin
         hod%M1 = hod%Mmin
         hod%alpha = 1.
      ELSE IF (hod%ihod == 3) THEN
         ! Toy model
      ELSE IF (hod%ihod == 4) THEN
         ! Toy model with integer mean numbers of galaxies
      ELSE IF (hod%ihod == 5) THEN
         ! Toy model with integer mean numbers of galaxies and no variance
         hod%stats_sat = stats_delta
      END IF

   END SUBROUTINE assign_HOD

   REAL FUNCTION mean_centrals(M, hod)

      ! Mean number of central galaxies in a halo
      ! Note that this is the mean, not necessarily 0 or 1 but should be between 0 and 1
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model

      IF (between(M, hod%Mmin, hod%Mmax)) THEN
         IF (hod%ihod == 1) THEN
            mean_centrals = Heaviside(M-hod%Mmin)
         ELSE IF (hod%ihod == 2) THEN
            mean_centrals = 0.5*(1.+erf(log(M/hod%Mmin)/hod%sigma))
         ELSE IF (hod%ihod == 3 .OR. hod%ihod == 4 .OR. hod%ihod == 5) THEN
            mean_centrals = 1.
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
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      REAL, PARAMETER :: eps = 1e-6

      IF (between(M, hod%Mmin, hod%Mmax)) THEN
         IF (hod%ihod == 1) THEN
            mean_satellites = (M/hod%M1)**hod%alpha
         ELSE IF (hod%ihod == 2) THEN
            mean_satellites = ((M-hod%M0)/hod%M1)**hod%alpha
         ELSE IF (hod%ihod == 3 .OR. hod%ihod == 4 .OR. hod%ihod == 5) THEN
            mean_satellites = (M/hod%Mmin)-1.
            IF (hod%ihod == 4 .OR. hod%ihod == 5) mean_satellites = floor(mean_satellites+eps)
         ELSE
            STOP 'MEAN_SATELLITES: Error, HOD not recognised'
         END IF
      ELSE
         mean_satellites = 0.
      END IF

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
      REAL :: mean

      mean = mean_centrals(M, hod)
      IF (mean .NE. 0.) THEN
         IF (hod%stats_cen == stats_delta) THEN
            variance_centrals = 0.
         ELSE IF (hod%stats_cen == stats_Bernoulli) THEN
            ! Variance for a Bernoulli distribution is p-p^2
            ! If p = 1 or p = 0 then the variance is automatically 0, as expected  
            variance_centrals = mean-mean**2
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
      REAL :: mean

      mean = mean_satellites(M, hod)
      IF (mean .NE. 0.) THEN
         IF (hod%stats_sat == stats_delta) THEN
            variance_satellites = 0.
         ELSE IF (hod%stats_sat == stats_Poisson) THEN
            ! Variance for a Poisson distribution is equal to the mean
            variance_satellites = mean
         ELSE
            STOP 'VARIANCE_SATELLITES: Error, satellite statistics not recognised'
         END IF
      ELSE
         variance_satellites = 0.
      END IF

   END FUNCTION variance_satellites

   REAL FUNCTION variance_galaxies(M, hod)

      ! Variance in the total number of galaxies in a halo
      ! Ignores possible covariance between the number of centrals and number of satellites
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model

      variance_galaxies = variance_centrals(M, hod)+variance_satellites(M, hod)

   END FUNCTION variance_galaxies

   INTEGER FUNCTION random_number_of_centrals(M, hod)

      ! Draw from the random distribution to get a number of central galaxies in a halo
      ! This must be either 0 or 1
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      REAL :: mean

      mean = mean_centrals(M, hod)
      IF (hod%stats_cen == stats_delta) THEN
         random_number_of_centrals = nint(mean)
      ELSE IF (hod%stats_cen == stats_Bernoulli) THEN
         random_number_of_centrals = random_Bernoulli(mean)
      ELSE
         STOP 'RANDOM_NUMBER_OF_CENTRALS: Error, central statistics not recognised'
      END IF

   END FUNCTION random_number_of_centrals

   INTEGER FUNCTION random_number_of_satellites(M, hod)

      ! Draw from the random distribution to get a number of satellite galaxies in a halo
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      REAL :: mean

      mean = mean_satellites(M, hod)
      IF (hod%stats_sat == stats_delta) THEN
         random_number_of_satellites = nint(mean)
      ELSE IF (hod%stats_sat == stats_Poisson) THEN
         random_number_of_satellites = random_Poisson(mean)
      ELSE
         STOP 'RANDOM_NUMBER_OF_SATELLITES: Error, satellite statistics not recognised'
      END IF

   END FUNCTION random_number_of_satellites

   INTEGER FUNCTION random_number_of_galaxies(M, hod)

      ! Draw from the random distribution to get a total number of galaxies in a halo
      REAL, INTENT(IN) :: M           ! Halo mass [Msun/h]
      TYPE(hodmod), INTENT(IN) :: hod ! HOD model
      INTEGER :: Nc, Ns

      Nc = random_number_of_centrals(M, hod)
      Ns = random_number_of_satellites(M, hod)
      random_number_of_galaxies = Nc+Ns

   END FUNCTION random_number_of_galaxies

END MODULE HOD_functions