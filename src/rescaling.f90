MODULE rescaling

   USE cosmology_functions

   IMPLICIT NONE

CONTAINS

   SUBROUTINE calculate_rescale_coefficients_power(kmin, kmax, cosm, cosm_target)

      REAL, INTENT(IN) :: kmin
      REAL, INTENT(OUT) :: kmax
      TYPE(cosmology), INTENT(IN) :: cosm
      TYPE(cosmology), INTENT(IN) :: cosm_target

   END SUBROUTINE calculate_rescale_coefficients_power

END MODULE rescaling