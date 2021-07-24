MODULE basic_operations

   IMPLICIT NONE
 
   PRIVATE

   PUBLIC :: progression
   PUBLIC :: fix_minimum
   PUBLIC :: fix_maximum
   PUBLIC :: read_command_argument
   PUBLIC :: positive
   PUBLIC :: negative
   PUBLIC :: even
   PUBLIC :: odd
   PUBLIC :: requal
   PUBLIC :: present_and_correct
   PUBLIC :: first_digit
   PUBLIC :: swap
   PUBLIC :: increment
   PUBLIC :: default_or_optional
   PUBLIC :: between
   PUBLIC :: equals_or

   INTERFACE swap
      MODULE PROCEDURE swap_real
      MODULE PROCEDURE swap_integer
   END INTERFACE swap

   INTERFACE read_command_argument
      MODULE PROCEDURE read_command_argument_real
      MODULE PROCEDURE read_command_argument_integer
      MODULE PROCEDURE read_command_argument_logical
      MODULE PROCEDURE read_command_argument_character
   END INTERFACE read_command_argument

   INTERFACE increment
      MODULE PROCEDURE increment_real
      MODULE PROCEDURE increment_integer
   END INTERFACE increment

   INTERFACE default_or_optional
      MODULE PROCEDURE default_or_optional_real
      MODULE PROCEDURE default_or_optional_integer
   END INTERFACE default_or_optional

   INTERFACE between
      MODULE PROCEDURE between_real
      MODULE PROCEDURE between_integer
   END INTERFACE between

CONTAINS

   LOGICAL FUNCTION equals_or(a, i)

      ! Returns true if a=i and b=j *or* a=j and b=i
      ! Useful for symmetric things
      INTEGER, INTENT(IN) :: a(:)
      INTEGER, INTENT(IN) :: i(:)

      equals_or = (a(1) == i(1) .AND. a(2) == i(2)) .OR. (a(1) == i(2) .AND. a(2) == i(1))

   END FUNCTION equals_or

   LOGICAL FUNCTION between_real(x, xmin, xmax)

      ! True if xmin <= x <= xmax
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: xmin
      REAL, INTENT(IN) :: xmax

      between_real = (x >= xmin .AND. x <= xmax)

   END FUNCTION between_real

   LOGICAL FUNCTION between_integer(x, xmin, xmax)

      ! True if xmin <= x <= xmax
      INTEGER, INTENT(IN) :: x
      INTEGER, INTENT(IN) :: xmin
      INTEGER, INTENT(IN) :: xmax

      between_integer = (x >= xmin .AND. x <= xmax)

   END FUNCTION between_integer

   REAL FUNCTION default_or_optional_real(x_default, x_optional)

      REAL, INTENT(IN) :: x_default
      REAL, OPTIONAL, INTENT(IN) :: x_optional

      IF (present(x_optional)) THEN
         default_or_optional_real = x_optional
      ELSE
         default_or_optional_real = x_default
      END IF

   END FUNCTION default_or_optional_real

   INTEGER FUNCTION default_or_optional_integer(x_default, x_optional)

      INTEGER, INTENT(IN) :: x_default
      INTEGER, OPTIONAL, INTENT(IN) :: x_optional

      IF (present(x_optional)) THEN
         default_or_optional_integer = x_optional
      ELSE
         default_or_optional_integer = x_default
      END IF

   END FUNCTION default_or_optional_integer

   RECURSIVE REAL FUNCTION progression(xmin, xmax, i, n, ilog)

      ! Split the region xmin -> xmax into n linearly spaced increments (including both xmin and xmax)
      ! Returns the value at the i-th point
      REAL, INTENT(IN) :: xmin
      REAL, INTENT(IN) :: xmax
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN) :: n
      LOGICAL, OPTIONAL, INTENT(IN) :: ilog

      IF (n == 1) THEN
         IF (xmin /= xmax) STOP 'PROGRESSION: Error, if n=1 then xmin and xmax must be identical'
         progression = xmin
      ELSE
         IF (present_and_correct(ilog)) THEN
            progression = exp(log(xmin)+log(xmax/xmin)*real(i-1)/real(n-1))
         ELSE
            progression = xmin+(xmax-xmin)*real(i-1)/real(n-1)
         END IF
      END IF

   END FUNCTION progression

   SUBROUTINE increment_real(x, y)

      ! Adds the value of y to x: x -> x+y
      REAL, INTENT(INOUT) :: x ! Value to be added to
      REAL, INTENT(IN) :: y    ! Value to add

      x = x+y

   END SUBROUTINE increment_real

   SUBROUTINE increment_integer(x, y)

      ! Adds the value of y to x: x -> x+y
      INTEGER, INTENT(INOUT) :: x ! Value to be added to
      INTEGER, INTENT(IN) :: y    ! Value to add

      x = x+y

   END SUBROUTINE increment_integer

   SUBROUTINE fix_minimum(x, xmin)

      ! If x is below xmin then set to xmin
      REAL, INTENT(INOUT) :: x ! Value to fix
      REAL, INTENT(IN) :: xmin ! Minimum value for x

      IF (x < xmin) x = xmin

   END SUBROUTINE fix_minimum

   SUBROUTINE fix_maximum(x, xmax)

      ! If x is above xmax then set to xmax
      REAL, INTENT(INOUT) :: x ! Value to fix
      REAL, INTENT(IN) :: xmax ! Maximum value for x

      IF (x > xmax) x = xmax

   END SUBROUTINE fix_maximum

   SUBROUTINE read_command_argument_real(i, x, desc, def)

      ! Read command-line argument i and use to to fill variable x
      INTEGER, INTENT(IN) :: i             ! Position of command-line argument
      REAL, INTENT(OUT) :: x               ! Real number to be assigned command-line argument
      CHARACTER(len=*), INTENT(IN) :: desc ! Description of command-line argument
      REAL, INTENT(IN), OPTIONAL :: def    ! Default value
      CHARACTER(len=256) :: word

      CALL get_command_argument(i, word)
      IF (word == '') THEN
         IF (present(def)) THEN
            x = def
         ELSE
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_REAL: Missing command-line argument number:', i
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_REAL: ', TRIM(desc)
            STOP
         END IF
      ELSE
         READ (word, *) x
      END IF

   END SUBROUTINE read_command_argument_real

   SUBROUTINE read_command_argument_integer(i, x, desc, def)

      ! Read command-line argument i and use to to fill variable x
      INTEGER, INTENT(IN) :: i             ! Position of command-line argument
      INTEGER, INTENT(OUT) :: x            ! Real number to be assigned command-line argument
      CHARACTER(len=*), INTENT(IN) :: desc ! Description of command-line argument
      INTEGER, INTENT(IN), OPTIONAL :: def ! Default value
      CHARACTER(len=256) :: word

      CALL get_command_argument(i, word)
      IF (word == '') THEN
         IF (present(def)) THEN
            x = def
         ELSE
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_INTEGER: Missing command-line argument:', i
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_INTEGER: ', TRIM(desc)
            STOP
         END IF
      ELSE
         READ (word, *) x
      END IF

   END SUBROUTINE read_command_argument_integer

   SUBROUTINE read_command_argument_logical(i, x, desc, def)

      ! Read command-line argument i and use to to fill variable x
      INTEGER, INTENT(IN) :: i             ! Position of command-line argument
      LOGICAL, INTENT(OUT) :: x            ! Real number to be assigned command-line argument
      CHARACTER(len=*), INTENT(IN) :: desc ! Description of command-line argument
      LOGICAL, INTENT(IN), OPTIONAL :: def ! Default value
      CHARACTER(len=256) :: word

      CALL get_command_argument(i, word)
      IF (word == '') THEN
         IF (present(def)) THEN
            x = def
         ELSE
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_LOGICAL: Missing command-line argument:', i
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_LOGICAL: ', TRIM(desc)
            STOP
         END IF
      ELSE
         IF (word == 'TRUE') THEN
            x = .TRUE.
         ELSE IF (word == 'FALSE') THEN
            x = .FALSE.
         ELSE
            STOP 'READ_COMMAND_ARGUMENT_LOGICAL: Error, should be either TRUE or FALSE'
         END IF
      END IF

   END SUBROUTINE read_command_argument_logical

   SUBROUTINE read_command_argument_character(i, x, desc, def)

      ! Read command-line argument i and use to to fill variable x
      INTEGER, INTENT(IN) :: i                      ! Position of command-line argument
      CHARACTER(len=*), INTENT(OUT) :: x            ! String to be assigned command-line argument
      CHARACTER(len=*), INTENT(IN) :: desc          ! Description of command-line argument
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: def ! Default value
      CHARACTER(len=256) :: word

      CALL get_command_argument(i, word)
      IF (word == '') THEN
         IF (present(def)) THEN
            x = def
         ELSE
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_CHARACTER: Missing command-line argument:', i
            WRITE (*, *) 'READ_COMMAND_ARGUMENT_CHARACTER: ', TRIM(desc)
            STOP
         END IF
      ELSE
         x = trim(word)
      END IF

   END SUBROUTINE read_command_argument_character

   LOGICAL FUNCTION positive(x)

      ! Returns .TRUE. if x>=0., if x=0. will return true
      REAL, INTENT(IN) :: x

      positive = (.NOT. negative(x))

   END FUNCTION positive

   LOGICAL FUNCTION negative(x)

      ! Returns true if argument is negative
      REAL, INTENT(IN) :: x

      negative = (x < 0.)

   END FUNCTION negative

   LOGICAL FUNCTION even(i)

      ! Tests for i being even, returns true if even
      INTEGER, INTENT(IN) :: i

      even = (mod(i, 2) == 0)

   END FUNCTION even

   LOGICAL FUNCTION odd(i)

      ! Tests for i being odd, returns true if odd
      INTEGER, INTENT(IN) :: i

      odd = (.NOT. even(i))

   END FUNCTION odd

   LOGICAL FUNCTION requal(x, y, eps)

      ! Tests if two real numbers are within some tolerance of each other
      ! If they are, then they should be considered equal
      ! Adapted from https://stackoverflow.com/questions/4915462/how-should-i-do-floating-point-comparison/4915891#4915891
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: y
      REAL, INTENT(IN) :: eps
      REAL :: absx, absy, diff

      absx = abs(x)
      absy = abs(y)
      diff = abs(x-y)
      IF (x == y) THEN
         requal = .TRUE.
      ELSE IF (x == 0. .OR. y == 0. .OR. diff < tiny(x)) THEN
         requal = (diff < eps*tiny(x))
      ELSE
         requal = (diff/(absx+absy) < 0.5*eps) ! June 2020: Added 0.5
      END IF

   END FUNCTION requal

   LOGICAL FUNCTION present_and_correct(x)

      ! Returns true if an optional argument is both present and true
      LOGICAL, OPTIONAL, INTENT(IN) :: x

      IF (present(x)) THEN
         present_and_correct = x
      ELSE
         present_and_correct = .FALSE.
      END IF

   END FUNCTION present_and_correct

   INTEGER FUNCTION first_digit(x)

      ! Returns the first non-zero digit of a real number
      REAL, INTENT(IN) :: x
      REAL :: y

      y = abs(x)
      DO
         IF (y == 1.) THEN
            first_digit = 1
            EXIT
         ELSE IF (y > 1. .AND. y < 10.) THEN
            first_digit = floor(y)
            EXIT
         ELSE IF (y >= 10.) THEN
            y = y/10.
         ELSE IF (y < 1.) THEN
            y = y*10.
         END IF
      END DO

   END FUNCTION first_digit

   SUBROUTINE swap_real(a, b)

      ! Swaps the values of variables a and b
      REAL, INTENT(INOUT) :: a
      REAL, INTENT(INOUT) :: b
      REAL :: c

      c = a
      a = b
      b = c

   END SUBROUTINE swap_real

   SUBROUTINE swap_integer(n, m)

      ! Swap integers n and m in the most memory-efficient way possible
      INTEGER, INTENT(INOUT) :: n
      INTEGER, INTENT(INOUT) :: m

      n = n+m ! n' = n+m
      m = n-m ! m' = n'-m = n+m-m = n
      n = n-m ! n'' = n'-m' = n+m-n = m

   END SUBROUTINE swap_integer

END MODULE basic_operations
