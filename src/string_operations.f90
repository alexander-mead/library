MODULE string_operations

   USE basic_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: number_file
   PUBLIC :: number_file2
   PUBLIC :: number_file_zeroes
   PUBLIC :: integer_to_string
   PUBLIC :: real_to_string
   PUBLIC :: exp_to_string
   PUBLIC :: snippet_in_string

CONTAINS

!!$  ! I cannot remember where (or why) I stole this from
!!$  elemental subroutine str2int(str,int,stat)
!!$
!!$    implicit none
!!$    character(len=*),intent(in) :: str
!!$    integer,intent(out)         :: int
!!$    integer,intent(out)         :: stat
!!$
!!$    read(str,*,iostat=stat) int
!!$
!!$  end subroutine str2int

   LOGICAL FUNCTION snippet_in_string(snip, string)

      ! Returns true if 'str1' exists within 'str2'
      CHARACTER(len=*), INTENT(IN) :: snip
      CHARACTER(len=*), INTENT(IN) :: string
      INTEGER :: i
      INTEGER :: l1, l2

      l1 = len(trim(snip))
      l2 = len(trim(string))

      snippet_in_string = .FALSE.

      DO i = 1, l2-l1+1
         IF (trim(snip) == string(i:i+l1-1)) THEN
            snippet_in_string = .TRUE.
            EXIT
         END IF
      END DO

   END FUNCTION snippet_in_string

   CHARACTER(len=8) FUNCTION integer_to_string(i, n)

      INTEGER, INTENT(IN) :: i
      INTEGER, OPTIONAL, INTENT(IN) :: n
      CHARACTER(len=8) :: fmt 

      IF (present(n)) THEN

         IF (i < 0) STOP 'INTEGER_TO_STRING: Error, cannot pad with negative number'
         IF (log10(real(i)) >= n) STOP 'INTEGER_TO_STRING: Error, number greater than padding size'
      
         IF (n == 1) THEN
            fmt = '(I0.1)'
         ELSE IF (n == 2) THEN
            fmt = '(I0.2)'
         ELSE IF (n == 3) THEN
            fmt = '(I0.3)'
         ELSE IF (n == 4) THEN
            fmt = '(I0.4)'
         ELSE
            STOP 'INTEGER_TO_STRING: Error, need to support padding > 4'
         END IF

      ELSE

         IF (i < 0) THEN
            IF (between(i, -9, -1)) THEN
               fmt = '(I2)'
            ELSE
               STOP 'INTEGER_TO_STRING: Error, your integer is not supported'
            END IF
         ELSE
            IF (between(i, 0, 9)) THEN
               fmt = '(I1)'
            ELSE IF (between(i, 10, 99)) THEN
               fmt = '(I2)'
            ELSE IF (between(i, 100, 999)) THEN
               fmt = '(I3)'
            ELSE IF (between(i, 1000, 9999)) THEN
               fmt = '(I4)'
            ELSE
               STOP 'INTEGER_TO_STRING: Error, your integer is not supported'
            END IF
         END IF

      END IF

      WRITE (integer_to_string, fmt=fmt) i

   END FUNCTION integer_to_string

   CHARACTER(len=8) FUNCTION real_to_string(x, pre_decimal, post_decimal)

      REAL, INTENT(IN) :: x               ! Real number to convert to string
      INTEGER, INTENT(IN) :: pre_decimal  ! Number of characters pre decimal point
      INTEGER, INTENT(IN) :: post_decimal ! Number of characters post decimal point
      CHARACTER(len=8) :: fmt

      IF (pre_decimal < 0 .OR. post_decimal < 0) THEN
         STOP 'REAL_TO_STRING: Error, neither pre or post decimal should be negative'
      ELSE IF (post_decimal == 0) THEN
         real_to_string = integer_to_string(int(x))
      ELSE
         IF (x < 0) THEN
            fmt = '(F'//trim(integer_to_string(pre_decimal+post_decimal+2))//'.'//trim(integer_to_string(post_decimal))//')'
         ELSE
            fmt = '(F'//trim(integer_to_string(pre_decimal+post_decimal+1))//'.'//trim(integer_to_string(post_decimal))//')'
         END IF
         WRITE (real_to_string, fmt=fmt) x
      END IF

   END FUNCTION real_to_string

   CHARACTER(len=16) FUNCTION exp_to_string(x, pre_decimal, post_decimal, exponent)

      REAL, INTENT(IN) :: x               ! Real number to convert to string
      INTEGER, INTENT(IN) :: pre_decimal  ! Number of characters pre decimal point
      INTEGER, INTENT(IN) :: post_decimal ! Number of characters post decimal point
      INTEGER, INTENT(IN) :: exponent     ! Exponent to write with

      exp_to_string = trim(real_to_string(x/10.**exponent, pre_decimal, post_decimal))//'e'//trim(integer_to_string(exponent))

   END FUNCTION exp_to_string

   FUNCTION number_file(fbase, i, fext)

      CHARACTER(len=256) :: number_file
      CHARACTER(len=*), INTENT(IN) :: fbase, fext
      INTEGER, INTENT(IN) :: i
      CHARACTER(len=8) num

      num = integer_to_string(i)
      number_file = trim(fbase)//trim(num)//trim(fext)

   END FUNCTION number_file

   FUNCTION number_file2(fbase, i1, mid, i2, fext)

      CHARACTER(len=256) ::number_file2
      CHARACTER(len=*), INTENT(IN) :: fbase, fext, mid
      INTEGER, INTENT(IN) :: i1, i2
      CHARACTER(len=8) :: num1, num2

      num1 = integer_to_string(i1)
      num2 = integer_to_string(i2)
      number_file2 = trim(fbase)//trim(num1)//trim(mid)//trim(num2)//trim(fext)

   END FUNCTION number_file2

   FUNCTION number_file_zeroes(fbase, i, num, fext)

      !Number a file with zero padding
      !Num specifies the number of digits
      CHARACTER(len=256) :: number_file_zeroes
      CHARACTER(len=*), INTENT(IN) :: fbase, fext
      CHARACTER(len=4) :: num4
      CHARACTER(len=3) :: num3
      CHARACTER(len=2) :: num2
      CHARACTER(len=1) :: num1
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN) :: num
      INTEGER, PARAMETER :: maxnum = 4

      IF (i < 0) STOP 'NUMBER_FILE_ZEROES: Error: cannot write negative number file names'

      IF (num > maxnum) STOP 'NUMBER_FILE_ZEROES: Error: need to add extra number capacity'

      IF (num == 1) THEN
         IF (i >= 10) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
         WRITE (num1, fmt='(I0.1)') i
         number_file_zeroes = trim(fbase)//num1//trim(fext)
      ELSE IF (num == 2) THEN
         IF (i >= 100) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
         WRITE (num2, fmt='(I0.2)') i
         number_file_zeroes = trim(fbase)//num2//trim(fext)
      ELSE IF (num == 3) THEN
         IF (i >= 1000) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
         WRITE (num3, fmt='(I0.3)') i
         number_file_zeroes = trim(fbase)//num3//trim(fext)
      ELSE IF (num == 4) THEN
         IF (i >= 10000) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
         WRITE (num4, fmt='(I0.4)') i
         number_file_zeroes = trim(fbase)//trim(num4)//trim(fext)
      END IF

   END FUNCTION number_file_zeroes

END MODULE string_operations
