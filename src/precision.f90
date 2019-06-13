MODULE precision

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: sp
   PUBLIC :: dl

   INTEGER, PARAMETER :: sp = KIND(1.0)
   INTEGER, PARAMETER :: dl = KIND(1.d0)

END MODULE precision
