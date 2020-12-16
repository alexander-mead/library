MODULE fft

   ! TODO: Re-enable single and double precision versions of functions
   ! TODO: Remove unnecessary array-size arguments from functions
   ! TODO: Clean up

   !USE iso_c_binding ! This use statement seems not to be necessary any more, not sure why
   USE basic_operations
   USE precision

   IMPLICIT NONE

   !Include the FFTW libraray locations
   !INCLUDE '/usr/local/include/fftw3.f' !Mac
   !INCLUDE '/usr/include/fftw3.f' !Linux'
   INCLUDE 'fftw3.f'

   PRIVATE

   PUBLIC :: FFT1
   PUBLIC :: FFT2
   PUBLIC :: FFT3

   PUBLIC :: k_FFT

   PUBLIC :: fundamental_frequency
   PUBLIC :: Nyquist_frequency
   PUBLIC :: Nyquist_cell
   PUBLIC :: conjugate_mode_location

   INTERFACE FFT1
      MODULE PROCEDURE FFT1_complex_double
      MODULE PROCEDURE FFT1_complex_single
      MODULE PROCEDURE FFT1_real_double
      MODULE PROCEDURE FFT1_real_single
   END INTERFACE FFT1

   INTERFACE FFT2
      MODULE PROCEDURE FFT2_complex_double
      MODULE PROCEDURE FFT2_complex_single
      MODULE PROCEDURE FFT2_real_double
      MODULE PROCEDURE FFT2_real_single
   END INTERFACE FFT2

   INTERFACE FFT3
      MODULE PROCEDURE FFT3_complex_double
      MODULE PROCEDURE FFT3_complex_single
      MODULE PROCEDURE FFT3_real_double
      MODULE PROCEDURE FFT3_real_single
   END INTERFACE FFT3

! This works, but cannot be called FFT as this is the same name as fft.f90
!!$  INTERFACE FFTs
!!$     MODULE PROCEDURE FFT1_complex_double
!!$     MODULE PROCEDURE FFT1_complex_single
!!$     MODULE PROCEDURE FFT1_real_double
!!$     MODULE PROCEDURE FFT1_real_single
!!$     MODULE PROCEDURE FFT2_complex_double
!!$     MODULE PROCEDURE FFT2_complex_single
!!$     MODULE PROCEDURE FFT2_real_double
!!$     MODULE PROCEDURE FFT2_real_single
!!$     MODULE PROCEDURE FFT3_complex_double
!!$     MODULE PROCEDURE FFT3_complex_single
!!$     MODULE PROCEDURE FFT3_real_double
!!$     MODULE PROCEDURE FFT3_real_single
!!$  END INTERFACE FFTs

CONTAINS

   INTEGER FUNCTION Nyquist_cell(m)

      ! Location of Nyquist cell for Fourier transform array of size m
      INTEGER, INTENT(IN) :: m

      IF(odd(m)) STOP 'NYQUIST_CELL: Error, need to think about Nyquist for odd FFT'
      Nyquist_cell = 1+m/2

   END FUNCTION Nyquist_cell

   INTEGER FUNCTION conjugate_mode_location(i, m)

      ! Assumes that m is even
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN) :: m

      IF (i == 1) THEN
         conjugate_mode_location = 1
      ELSE
         conjugate_mode_location = m+2-i
      END IF

   END FUNCTION conjugate_mode_location

   REAL FUNCTION fundamental_frequency(L)

      USE constants
      REAL, INTENT(IN) :: L

      fundamental_frequency =  twopi/L

   END FUNCTION fundamental_frequency

   REAL FUNCTION Nyquist_frequency(L, m)

      USE constants
      REAL, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: m

      Nyquist_frequency =  pi*m/L

   END FUNCTION Nyquist_frequency

   SUBROUTINE k_FFT(ix, iy, iz, m, kx, ky, kz, kmod, L)

      ! Finds the wavevector associated with the FFT cell
      USE constants
      INTEGER, INTENT(IN) :: ix, iy, iz, m
      REAL, INTENT(IN) :: L
      REAL, INTENT(OUT) :: kx, ky, kz, kmod

      IF (odd(m)) STOP 'K_FFT: Fourier transform does not have an even mesh'

      kx = real(ix-1)
      ky = real(iy-1)
      kz = real(iz-1)

      IF (ix > m/2+1) kx = -real(m-ix+1)
      IF (iy > m/2+1) ky = -real(m-iy+1)
      IF (iz > m/2+1) kz = -real(m-iz+1)

      kx = kx*twopi/L
      ky = ky*twopi/L
      kz = kz*twopi/L

      kmod = sqrt(kx**2+ky**2+kz**2)

   END SUBROUTINE k_FFT

   SUBROUTINE FFT1_complex_double(in, out, n, ifb)

      ! Wrapper for the 1D FFTW
      INTEGER, INTENT(IN) :: n
      COMPLEX(dp), INTENT(IN) :: in(n) 
      COMPLEX(dp), INTENT(OUT) :: out(n)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(n)) STOP 'FFT1_COMPLEX_DOUBLE Error, the array should be even'

      ! Create the plan for the FFT
      IF (verbose) WRITE (*, *) 'FFT1_COMPLEX_DOUBLE Starting FFT - creating plan'
      IF (ifb == -1) THEN
         CALL dfftw_plan_dft_1d(plan, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
      ELSE IF (ifb == 1) THEN
         CALL dfftw_plan_dft_1d(plan, n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
      ELSE
         STOP 'FFT1_COMPLEX_DOUBLE Error - need to specify forwards or backwards'
      END IF

      ! This computes the FFT...
      IF (verbose) WRITE (*, *) 'FFT1_COMPLEX_DOUBLE Executing FFTW'
      CALL dfftw_execute(plan)
      IF (verbose) WRITE (*, *) 'FFT1_COMPLEX_DOUBLE FFTW complete'

      ! ...and this destroys the plan
      CALL dfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT1_COMPLEX_DOUBLE Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT1_complex_double

   SUBROUTINE FFT1_complex_single(in, out, n, ifb)

      ! Wrapper for the 1D FFTW
      INTEGER, INTENT(IN) :: n
      COMPLEX(sp), INTENT(IN) :: in(n)
      COMPLEX(sp), INTENT(OUT) :: out(n)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(n)) STOP 'FFT1_COMPELX_SINGLE: Error, the array should be even'

      ! Create the plan for the FFT
      IF (verbose) WRITE (*, *) 'FFT1_COMPELX_SINGLE: Starting FFT - creating plan'
      IF (ifb == -1) THEN
         CALL sfftw_plan_dft_1d(plan, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
      ELSE IF (ifb == 1) THEN
         CALL sfftw_plan_dft_1d(plan, n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
      ELSE
         STOP 'FFT1_COMPELX_SINGLE: Error - need to specify forwards or backwards'
      END IF

      ! This computes the FFT...
      IF (verbose) WRITE (*, *) 'FFT1_COMPELX_SINGLE: Executing FFTW'
      CALL sfftw_execute(plan)
      IF (verbose) WRITE (*, *) 'FFT1_COMPELX_SINGLE: FFTW complete'

      ! ...and this destroys the plan
      CALL sfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT1_COMPELX_SINGLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT1_complex_single

   SUBROUTINE FFT1_real_double(rspace, fspace, n, ifb)

      ! Wrapper for the real 1D FFTW
      INTEGER, INTENT(IN) :: n
      REAL(dp), INTENT(INOUT) :: rspace(n)
      COMPLEX(dp), INTENT(INOUT) :: fspace(n/2+1)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan

      IF (odd(n)) STOP 'FFT1_REAL_DOUBLE: Error, the array should be even'

      IF (ifb == -1) THEN
         CALL dfftw_plan_dft_r2c_1d(plan, N, rspace, fspace, FFTW_ESTIMATE)
         CALL dfftw_execute_dft_r2c(plan, rspace, fspace)
      ELSE IF (ifb == 1) THEN
         CALL dfftw_plan_dft_c2r_1d(plan, N, fspace, rspace, FFTW_ESTIMATE)
         CALL dfftw_execute_dft_c2r(plan, fspace, rspace)
      ELSE
         STOP 'FFT1_REAL_DOUBLE: Error - need to specify forwards or backwards'
      END IF

      CALL dfftw_destroy_plan(plan)

   END SUBROUTINE FFT1_real_double

   SUBROUTINE FFT1_real_single(rspace, fspace, n, ifb)

      ! Wrapper for the real 1D FFTW
      INTEGER, INTENT(IN) :: n
      REAL(sp), INTENT(INOUT) :: rspace(n)   
      COMPLEX(sp), INTENT(INOUT) :: fspace(n/2+1)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan

      IF (odd(n)) STOP 'FFT1_REAL_SINGLE: Error, the array should be even'

      IF (ifb == -1) THEN
         CALL sfftw_plan_dft_r2c_1d(plan, N, rspace, fspace, FFTW_ESTIMATE)
         CALL sfftw_execute_dft_r2c(plan, rspace, fspace)
      ELSE IF (ifb == 1) THEN
         CALL sfftw_plan_dft_c2r_1d(plan, N, fspace, rspace, FFTW_ESTIMATE)
         CALL sfftw_execute_dft_c2r(plan, fspace, rspace)
      ELSE
         STOP 'FFT1_REAL_SINGLE: Error - need to specify forwards or backwards'
      END IF

      CALL sfftw_destroy_plan(plan)

   END SUBROUTINE FFT1_real_single

   SUBROUTINE FFT2_complex_double(in, out, nx, ny, ifb)

      ! Wrapper for the 2D FFTW
      INTEGER, INTENT(IN) :: nx, ny
      COMPLEX(dp), INTENT(IN) :: in(nx, ny)
      COMPLEX(dp), INTENT(OUT) :: out(nx, ny)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT2_COMPLEX_DOUBLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT2_COMPLEX_DOUBLE: Error, the array should be even'

      ! Create the plans for the FFT
      IF (verbose) WRITE (*, *) 'FFT2_COMPLEX_DOUBLE: Starting FFT - creating plan'
      IF (ifb == -1) THEN
         CALL dfftw_plan_dft_2d(plan, nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
      ELSE IF (ifb == 1) THEN
         CALL dfftw_plan_dft_2d(plan, nx, ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
      ELSE
         WRITE (*, *) 'FFT2_COMPLEX_DOUBLE: Error - need to specify forwards or backwards'
      END IF

      ! This computes the FFT...
      IF (verbose) WRITE (*, *) 'FFT2_COMPLEX_DOUBLE: Executing FFTW'
      CALL dfftw_execute(plan)
      IF (verbose) WRITE (*, *) 'FFT2_COMPLEX_DOUBLE: FFTW complete'

      ! ...and this destroys the plan!
      CALL dfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT2_COMPLEX_DOUBLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT2_complex_double

   SUBROUTINE FFT2_complex_single(in, out, nx, ny, ifb)

      ! Wrapper for the 2D FFTW
      INTEGER, INTENT(IN) :: nx, ny
      COMPLEX(sp), INTENT(IN) :: in(nx, ny)
      COMPLEX(sp), INTENT(OUT) :: out(nx, ny)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT2_COMPLEX_SINGLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT2_COMPLEX_SINGLE: Error, the array should be even'

      ! Create the plans for the FFT
      IF (verbose) WRITE (*, *) 'FFT2_COMPLEX_SINGLE: Starting FFT - creating plan'
      IF (ifb == -1) THEN
         CALL sfftw_plan_dft_2d(plan, nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
      ELSE IF (ifb == 1) THEN
         CALL sfftw_plan_dft_2d(plan, nx, ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
      ELSE
         WRITE (*, *) 'FFT2_COMPLEX_SINGLE: Error - need to specify forwards or backwards'
      END IF

      ! This computes the FFT...
      IF (verbose) WRITE (*, *) 'FFT2_COMPLEX_SINGLE: Executing FFTW'
      CALL sfftw_execute(plan)
      IF (verbose) WRITE (*, *) 'FFT2_COMPLEX_SINGLE: FFTW complete'

      ! ...and this destroys the plan!
      CALL sfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT2_COMPLEX_SINGLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT2_complex_single

   SUBROUTINE FFT2_real_double(rspace, fspace, nx, ny, ifb)

      ! Wrapper for the real 3D FFTW
      INTEGER, INTENT(IN) :: nx, ny
      REAL(dp), INTENT(INOUT)  :: rspace(nx, ny)
      COMPLEX(dp), INTENT(INOUT) :: fspace(nx/2+1, ny)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT2_REAL_DOUBLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT2_REAL_DOUBLE: Error, the array should be even'

      ! Create the plans for the FFT and execute
      IF (verbose) WRITE (*, *) 'FFT2_REAL_DOUBLE: Starting FFT - creating plan and executing simulatanesouly'
      IF (ifb == -1) THEN
         CALL dfftw_plan_dft_r2c_2d(plan, nx, ny, rspace, fspace, FFTW_ESTIMATE)
         CALL dfftw_execute_dft_r2c(plan, rspace, fspace)
      ELSE IF (ifb == 1) THEN
         CALL dfftw_plan_dft_c2r_2d(plan, nx, ny, fspace, rspace, FFTW_ESTIMATE)
         CALL dfftw_execute_dft_c2r(plan, fspace, rspace)
      ELSE
         WRITE (*, *) 'FFT2_REAL_DOUBLE: Error - need to specify forwards or backwards'
      END IF

      ! ..and this destroys the plan!
      CALL dfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT2_REAL_DOUBLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT2_real_double

   SUBROUTINE FFT2_real_single(rspace, fspace, nx, ny, ifb)

      ! Wrapper for the real 3D FFTW
      INTEGER, INTENT(IN) :: nx, ny
      REAL(sp), INTENT(INOUT)  :: rspace(nx, ny)
      COMPLEX(sp), INTENT(INOUT) :: fspace(nx/2+1, ny)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT2_REAL_SINGLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT2_REAL_SINGLE: Error, the array should be even'

      ! Create the plans for the FFT and execute
      IF (verbose) WRITE (*, *) 'FFT2_REAL_SINGLE: Starting FFT - creating plan and executing simulatanesouly'
      IF (ifb == -1) THEN
         CALL sfftw_plan_dft_r2c_2d(plan, nx, ny, rspace, fspace, FFTW_ESTIMATE)
         CALL sfftw_execute_dft_r2c(plan, rspace, fspace)
      ELSE IF (ifb == 1) THEN
         CALL sfftw_plan_dft_c2r_2d(plan, nx, ny, fspace, rspace, FFTW_ESTIMATE)
         CALL sfftw_execute_dft_c2r(plan, fspace, rspace)
      ELSE
         WRITE (*, *) 'FFT2_REAL_SINGLE: Error - need to specify forwards or backwards'
      END IF

      ! ..and this destroys the plan!
      CALL sfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT2_REAL_SINGLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT2_real_single

   SUBROUTINE FFT3_complex_double(in, out, nx, ny, nz, ifb)

      ! Wrapper for the 3D FFTW
      INTEGER, INTENT(IN) :: nx, ny, nz
      COMPLEX(dp), INTENT(IN)  :: in(nx, ny, nz)
      COMPLEX(dp), INTENT(OUT) :: out(nx, ny, nz)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT3_COMPLEX_DOUBLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT3_COMPLEX_DOUBLE: Error, the array should be even'
      IF (odd(nz)) STOP 'FFT3_COMPLEX_DOUBLE: Error, the array should be even'

      ! Create the plans for the FFT
      IF (verbose) WRITE (*, *) 'FFT3_COMPLEX_DOUBLE: Starting FFT - creating plan'
      IF (ifb == -1) THEN
         CALL dfftw_plan_dft_3d(plan, nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
      ELSE IF (ifb == 1) THEN
         CALL dfftw_plan_dft_3d(plan, nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
      ELSE
         WRITE (*, *) 'FFT3_COMPLEX_DOUBLE: Error - need to specify forwards or backwards'
      END IF

      ! This computes the FFT...
      IF (verbose) WRITE (*, *) 'FFT3_COMPLEX_DOUBLE: Executing FFTW'
      CALL dfftw_execute(plan)
      IF (verbose) WRITE (*, *) 'FFT3_COMPLEX_DOUBLE: FFTW complete'

      ! ...and this destroys the plan!
      CALL dfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT3_COMPLEX_DOUBLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT3_complex_double

   SUBROUTINE FFT3_complex_single(in, out, nx, ny, nz, ifb)

      ! Wrapper for the 3D FFTW
      INTEGER, INTENT(IN) :: nx, ny, nz
      COMPLEX(sp), INTENT(IN)  :: in(nx, ny, nz)
      COMPLEX(sp), INTENT(OUT) :: out(nx, ny, nz)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT3_COMPLEX_SINGLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT3_COMPLEX_SINGLE: Error, the array should be even'
      IF (odd(nz)) STOP 'FFT3_COMPLEX_SINGLE: Error, the array should be even'

      ! Create the plans for the FFT
      IF (verbose) WRITE (*, *) 'FFT3_COMPLEX_SINGLE: Starting FFT - creating plan'
      IF (ifb == -1) THEN
         CALL sfftw_plan_dft_3d(plan, nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
      ELSE IF (ifb == 1) THEN
         CALL sfftw_plan_dft_3d(plan, nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
      ELSE
         WRITE (*, *) 'FFT3_COMPLEX_SINGLE: Error - need to specify forwards or backwards'
      END IF

      ! This computes the FFT...
      IF (verbose) WRITE (*, *) 'FFT3_COMPLEX_SINGLE: Executing FFTW'
      CALL sfftw_execute(plan)
      IF (verbose) WRITE (*, *) 'FFT3_COMPLEX_SINGLE: FFTW complete'

      ! ...and this destroys the plan!
      CALL sfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT3_COMPLEX_SINGLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT3_complex_single

   SUBROUTINE FFT3_real_double(rspace, fspace, nx, ny, nz, ifb)

      ! Wrapper for the real 3D FFTW
      INTEGER, INTENT(IN) :: nx, ny, nz
      REAL(dp), INTENT(INOUT)  :: rspace(nx, ny, nz)
      COMPLEX(dp), INTENT(INOUT) :: fspace(nx/2+1, ny, nz)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT3_REAL_DOUBLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT3_REAL_DOUBLE: Error, the array should be even'
      IF (odd(nz)) STOP 'FFT3_REAL_DOUBLE: Error, the array should be even'

      ! Create the plans for the FFT and execute
      IF (verbose) WRITE (*, *) 'FFT3_REAL_DOUBLE: Starting FFT - creating plan and executing simultanesouly'
      IF (ifb == -1) THEN
         CALL dfftw_plan_dft_r2c_3d(plan, nx, ny, nz, rspace, fspace, FFTW_ESTIMATE)
         CALL dfftw_execute_dft_r2c(plan, rspace, fspace)
      ELSE IF (ifb == 1) THEN
         CALL dfftw_plan_dft_c2r_3d(plan, nx, ny, nz, fspace, rspace, FFTW_ESTIMATE)
         CALL dfftw_execute_dft_c2r(plan, fspace, rspace)
      ELSE
         WRITE (*, *) 'FFT3_REAL_DOUBLE: Error - need to specify forwards or backwards'
      END IF

      ! ..and this destroys the plan!
      CALL dfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT3_REAL_DOUBLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT3_real_double

   SUBROUTINE FFT3_real_single(rspace, fspace, nx, ny, nz, ifb)

      ! Wrapper for the real 3D FFTW
      INTEGER, INTENT(IN) :: nx, ny, nz
      REAL(sp), INTENT(INOUT)  :: rspace(nx, ny, nz)
      COMPLEX(sp), INTENT(INOUT) :: fspace(nx/2+1, ny, nz)
      INTEGER, INTENT(IN) :: ifb
      INTEGER(int8) :: plan
      LOGICAL, PARAMETER :: verbose = .FALSE.

      IF (odd(nx)) STOP 'FFT3_REAL_SINGLE: Error, the array should be even'
      IF (odd(ny)) STOP 'FFT3_REAL_SINGLE: Error, the array should be even'
      IF (odd(nz)) STOP 'FFT3_REAL_SINGLE: Error, the array should be even'

      ! Create the plans for the FFT and execute
      IF (verbose) WRITE (*, *) 'FFT3_REAL_SINGLE: Starting FFT - creating plan and executing simulatanesouly'
      IF (ifb == -1) THEN
         CALL sfftw_plan_dft_r2c_3d(plan, nx, ny, nz, rspace, fspace, FFTW_ESTIMATE)
         CALL sfftw_execute_dft_r2c(plan, rspace, fspace)
      ELSE IF (ifb == 1) THEN
         CALL sfftw_plan_dft_c2r_3d(plan, nx, ny, nz, fspace, rspace, FFTW_ESTIMATE)
         CALL sfftw_execute_dft_c2r(plan, fspace, rspace)
      ELSE
         WRITE (*, *) 'FFT3_REAL_SINGLE: Error - need to specify forwards or backwards'
      END IF

      ! ..and this destroys the plan!
      CALL sfftw_destroy_plan(plan)
      IF (verbose) THEN
         WRITE (*, *) 'FFT3_REAL_SINGLE: Plan destroyed'
         WRITE (*, *)
      END IF

   END SUBROUTINE FFT3_real_single

END MODULE fft
