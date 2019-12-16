PROGRAM limber_demo

   USE cosmology_functions
   USE Limber
   IMPLICIT NONE
   INTEGER :: idemo
   INTEGER :: icosmo

   WRITE(*, *) 
   WRITE(*, *) 'Choose demo'
   WRITE(*, *) '1 - n(z) normalisation check'
   WRITE(*, *) '2 - Lensing diagnostics'
   READ(*, *) idemo
   WRITE(*, *)  

   IF(idemo == 1) THEN
      CALL nz_normalisation()
   ELSE IF (idemo == 2) THEN
      icosmo = -1
      CALL lensing_diagnostics(icosmo)
   ELSE
      STOP 'LIMBER_DEMO: Error, idemo specified incorrectly'
   END IF

CONTAINS

SUBROUTINE nz_normalisation()

   ! n(z) normalisation check
   USE calculus_table
   IMPLICIT NONE
   TYPE(projection) :: proj

   INTEGER :: i
   INTEGER :: nz
   INTEGER, PARAMETER :: nnz = 16

   WRITE (*, *) 'HMx_DRIVER: Checking n(z) functions'
   WRITE (*, *)

   ! Number of n(z) to check
   DO i = 1, nnz
      IF (i == 1) nz = tracer_RCSLenS
      IF (i == 2) nz = tracer_CFHTLenS_vanWaerbeke2013
      IF (i == 3) nz = tracer_KiDS
      IF (i == 4) nz = tracer_KiDS_bin1
      IF (i == 5) nz = tracer_KiDS_bin2
      IF (i == 6) nz = tracer_KiDS_bin3
      IF (i == 7) nz = tracer_KiDS_bin4
      IF (i == 8) nz = tracer_KiDS_450
      IF (i == 9) nz = tracer_KiDS_450_fat_bin1
      IF (i == 10) nz = tracer_KiDS_450_fat_bin2
      IF (i == 11) nz = tracer_KiDS_450_highz
      IF (i == 12) nz = tracer_KiDS_450_bin1
      IF (i == 13) nz = tracer_KiDS_450_bin2
      IF (i == 14) nz = tracer_KiDS_450_bin3
      IF (i == 15) nz = tracer_KiDS_450_bin4
      IF (i == 16) nz = tracer_CFHTLenS_Kilbinger2013
      WRITE (*, *) 'HMx_DRIVER: n(z) number:', nz
      WRITE (*, *) 'HMx_DRIVER: n(z) name: ', trim(xcorr_type(nz))
      CALL read_nz(nz, proj)
      WRITE (*, *) 'HMx_DRIVER: integration order: ', proj%order_nz
      WRITE (*, *) 'HMx_DRIVER: n(z) integral (hist):', integrate_table(proj%z_nz, proj%nz, proj%nnz, 1, proj%nnz, 0)
      WRITE (*, *) 'HMx_DRIVER: n(z) integral (linear):', integrate_table(proj%z_nz, proj%nz, proj%nnz, 1, proj%nnz, 1)
      WRITE (*, *) 'HMx_DRIVER: n(z) integral (quadratic):', integrate_table(proj%z_nz, proj%nz, proj%nnz, 1, proj%nnz, 2)
      WRITE (*, *) 'HMx_DRIVER: n(z) integral (cubic):', integrate_table(proj%z_nz, proj%nz, proj%nnz, 1, proj%nnz, 3)
      WRITE (*, *)
   END DO

END SUBROUTINE nz_normalisation

SUBROUTINE lensing_diagnostics(icosmo)

   ! Projection diagnostics
   IMPLICIT NONE
   INTEGER, INTENT(INOUT) :: icosmo
   INTEGER :: i, j, ix(2), ip(2)
   CHARACTER(len=256) :: outfile
   TYPE(cosmology) :: cosm
   TYPE(projection) :: proj(2)

   LOGICAL, PARAMETER :: verbose = .TRUE.

   ! Assigns the cosmological model
   CALL assign_cosmology(icosmo, cosm, verbose)
   CALL init_cosmology(cosm)
   CALL print_cosmology(cosm)

   ! Set the field types
   ix = -1
   DO i = 1, 2
      CALL set_field_for_xpow(ix(i), ip(i))
   END DO

   ! Fill the projection kernels (plural)
   CALL fill_projection_kernels(ix, proj, cosm)

   DO j = 1, 2

      IF (j == 1) outfile = 'data/nz1.dat'
      IF (j == 2) outfile = 'data/nz2.dat'
      OPEN (7, file=outfile)
      IF (ALLOCATED(proj(j)%nz)) THEN
         CALL write_nz(proj(j), outfile)
      ELSE
         WRITE (7, *) 0., 0.
      END IF
      CLOSE (7)

      IF (j == 1) outfile = 'data/efficiency1.dat'
      IF (j == 2) outfile = 'data/efficiency2.dat'
      OPEN (7, file=outfile)
      IF (ALLOCATED(proj(j)%q)) THEN
         CALL write_efficiency(proj(j), cosm, outfile)
      ELSE
         WRITE (7, *) 0., 0.
      END IF
      CLOSE (7)

      IF (j == 1) outfile = 'data/kernel1.dat'
      IF (j == 2) outfile = 'data/kernel2.dat'
      CALL write_projection_kernel(proj(j), cosm, outfile)

   END DO

END SUBROUTINE lensing_diagnostics

END PROGRAM