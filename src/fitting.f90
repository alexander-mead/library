MODULE fitting

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: fit_constant

CONTAINS

   REAL FUNCTION fit_constant(data, weight)

      REAL, INTENT(IN) :: data(:)
      REAL, INTENT(IN) :: weight(:)
      INTEGER :: n
      REAL :: total_weight

      n = size(data)
      IF (n /= size(weight)) STOP 'FIT_CONSTANT: Error, data and weight must be the same size'

      total_weight = sum(weight)
      IF (total_weight == 0.) THEN
         STOP 'FIT_CONSTANT: Error, weights sum to zero'
      ELSE
         fit_constant = sum(data*weight)/total_weight
      END IF

   END FUNCTION fit_constant

END MODULE fitting
