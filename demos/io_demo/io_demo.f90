PROGRAM owls_demo

   USE io
   IMPLICIT NONE

   CALL read_VD20_demo()

   CONTAINS

      SUBROUTINE read_VD20_demo()

         IMPLICIT NONE
         REAL, ALLOCATABLE :: k(:), z(:), Pk(:, :)
         INTEGER :: ik, iz, nk, nz
         CHARACTER(len=128), PARAMETER :: name = 'BAHAMAS_Theat7.6_nu0_WMAP9'

         CALL read_VD20_power(k, z, Pk, nk, nz, name)

         WRITE(*, *) 'READ_VD20_DEMO: k array'
         DO ik = 1, nk
            WRITE(*, *) ik, k(ik)
         END DO
         WRITE(*, *)

         WRITE(*, *) 'READ_VD20_DEMO: z array'
         DO iz = 1, nz
            WRITE(*, *) iz, z(iz)
         END DO
         WRITE(*, *)

         WRITE(*, *) 'READ_VD20_DEMO: Power array at z=0'
         DO ik = 1, nk
            WRITE(*, *) ik, k(ik), Pk(ik, nz)
         END DO
         WRITE(*, *)

      END SUBROUTINE read_VD20_demo

END PROGRAM owls_demo