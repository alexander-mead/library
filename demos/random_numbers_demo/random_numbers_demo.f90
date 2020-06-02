PROGRAM random_numbers_demo

   USE constants
   USE random_numbers
   USE special_functions

   IMPLICIT NONE
   INTEGER :: i, j, iran, ri
   REAL, ALLOCATABLE :: binval(:), func(:), bin(:), x(:)
   REAL :: rn, binmin, binmax, mean, sig, pn
   LOGICAL :: continuous

   INTEGER, PARAMETER :: n = 100000 !Number of random numbers to generate
   INTEGER, PARAMETER :: nbins = 200 !Number of bins for histogram
   INTEGER, PARAMETER :: iseed = 0 !Random number seed

   WRITE (*, *)
   WRITE (*, *) 'Random number generation'
   WRITE (*, *)

   !Set the random number generator
   !CALL RNG_set(iseed)
   CALL random_generator_seed(iseed)

   WRITE (*, *) 'Select distribution'
   WRITE (*, *) '1 - Uniform distribution'
   WRITE (*, *) '2 - Rayleigh distribution'
   WRITE (*, *) '3 - exponential distribution'
   WRITE (*, *) '4 - Gaussian distribution'
   WRITE (*, *) '5 - Polynomial distribution'
   WRITE (*, *) '6 - Lorentzian distribution'
   WRITE (*, *) '7 - Uniform integer distribution'
   WRITE (*, *) '8 - Log-normal distribution'
   READ (*, *) iran
   WRITE (*, *)

   IF (iran == 1) THEN
      binmin = 0.
      binmax = 1.
      mean = 0.5
      sig = 1.
      continuous = .TRUE.
   ELSE IF (iran == 2) THEN
      binmin = 0.
      binmax = 10.
      sig = 2.
      continuous = .TRUE.
   ELSE IF (iran == 3) THEN
      binmin = 0.
      binmax = 10.
      mean = 2.
      continuous = .TRUE.
   ELSE IF (iran == 4) THEN
      binmin = -10.
      binmax = 10.
      sig = 2.
      mean = 2.
      continuous = .TRUE.
   ELSE IF (iran == 5) THEN
      binmin = 0.
      binmax = 1.
      !WRITE(*,*) 'What order polynomial?'
      !READ(*,*) pn
      pn = 2.
      continuous = .TRUE.
   ELSE IF (iran == 6) THEN
      binmin = 0.
      binmax = 10.
      continuous = .TRUE.
   ELSE IF (iran == 7) THEN
      binmin = 0.
      binmax = 10.
      continuous = .FALSE.
   ELSE IF (iran == 8) THEN
      binmin = 0.
      binmax = 10.
      mean = 2.
      sig = 1.
      continuous = .TRUE.
   ELSE
      STOP 'RANDOM_NUMBERS_TEST: Error, demonstration not chosen correctly'
   END IF

   !Allocate arrays for the bins
   ALLOCATE (binval(nbins+1), bin(nbins), func(nbins), x(nbins))
   bin = 0

   !Make the bin-edges and store values in array
   DO i = 1, nbins+1
      binval(i) = binmin+(binmax-binmin)*float(i-1)/float(nbins)
   END DO

   !Loop over the bins
   DO i = 1, nbins

      !Set x value to be the bin mean
      x(i) = (binval(i)+binval(i+1))/2.

      !Draw the function values corresponding to the bin mean
      IF (iran == 1) func(i) = uniform(x(i), 0., 1.)
      IF (iran == 2) func(i) = Rayleigh(x(i), sig)
      IF (iran == 3) func(i) = exponential(x(i), mean)
      IF (iran == 4) func(i) = Gaussian(x(i), mean, sig)
      IF (iran == 5) func(i) = polynomial(x(i), pn)
      IF (iran == 6) func(i) = Lorentzian(x(i))
      IF (iran == 7) func(i) = 1.
      IF (iran == 8) func(i) = lognormal(x(i), mean, sig)

   END DO

   !Generate the random numbers and write to disk
   WRITE (*, *) 'Creating random numbers:', n
   OPEN (8, file='random_numbers.dat')
   DO i = 1, n

      IF (iran == 1) rn = random_uniform(0., 1.)
      IF (iran == 2) rn = random_Rayleigh(sig)
      IF (iran == 3) rn = random_exponential(mean)
      IF (iran == 4) rn = random_Gaussian(mean, sig)
      IF (iran == 5) rn = random_polynomial(pn)
      IF (iran == 6) rn = random_Lorentzian()
      IF (iran == 7) ri = random_integer(1, 9)
      IF (iran == 8) rn = random_lognormal(mean, sig)

      IF (i <= 20) WRITE (*, fmt='(I10,F14.7)') i, rn

      IF (continuous) THEN

         WRITE (8, *) i, rn

         DO j = 1, nbins
            IF (rn >= binval(j) .AND. rn < binval(j+1)) THEN
               bin(j) = bin(j)+1.
               EXIT
            END IF
         END DO

      ELSE

         WRITE (8, *) i, ri

         DO j = 1, nbins
            IF (REAL(ri) >= binval(j) .AND. REAL(ri) < binval(j+1)) THEN
               bin(j) = bin(j)+1.
               EXIT
            END IF
         END DO

      END IF

   END DO
   CLOSE (8)

   !Normalise the histogram
   DO i = 1, nbins
      bin(i) = bin(i)/(n*(binval(i+1)-binval(i)))
   END DO

   !Write the histogram to disk
   OPEN (7, file='histogram.dat')
   DO i = 1, nbins
      WRITE (7, *) x(i), bin(i), func(i)
   END DO
   CLOSE (7)

   !Finish
   WRITE (*, *) 'Done'
   WRITE (*, *)

END PROGRAM random_numbers_demo
