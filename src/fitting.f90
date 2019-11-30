MODULE fitting

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: fit_constant

CONTAINS

   REAL FUNCTION fit_constant(data, weight, n)

      IMPLICIT NONE
      REAL :: data(n)
      REAL :: weight(n)
      INTEGER :: n
      REAL :: total_weight

      total_weight = sum(weight)
      IF(total_weight==0.) THEN
         STOP 'FIT_CONSTANT: Error, weights sum to zero'
      ELSE
         fit_constant = sum(data*weight)/total_weight
      END IF

   END FUNCTION fit_constant

END MODULE fitting
