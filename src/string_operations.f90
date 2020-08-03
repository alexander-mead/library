MODULE string_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: number_file
   PUBLIC :: number_file2
   PUBLIC :: number_file_zeroes
   PUBLIC :: integer_to_string
   PUBLIC :: real_to_string
   PUBLIC :: string_in_string

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

   LOGICAL FUNCTION string_in_string(str1, str2)

      ! Returns true if 'str1' exists within 'str2'
      CHARACTER(len=*), INTENT(IN) :: str1
      CHARACTER(len=*), INTENT(IN) :: str2
      INTEGER :: i
      INTEGER :: l1, l2

      l1 = len(trim(str1))
      l2 = len(trim(str2))

      string_in_string = .FALSE.

      DO i = 1, l2-l1+1
         IF (trim(str1) == str2(i:i+l1-1)) THEN
            string_in_string = .TRUE.
            EXIT
         END IF
      END DO

   END FUNCTION string_in_string

   CHARACTER(len=8) FUNCTION integer_to_string(i)

      INTEGER, INTENT(IN) :: i
      CHARACTER(len=8) :: fmt   

      IF (i >= 0 .AND. i < 10) THEN
         fmt = '(I1)'
      ELSE IF (i >= 10 .AND. i < 100) THEN
         fmt = '(I2)'
      ELSE IF (i >= 100 .AND. i < 1000) THEN
         fmt = '(I3)'
      ELSE IF (i >= 1000 .AND. i < 10000) THEN
         fmt = '(I4)'
      ELSE
         STOP 'INTEGER_TO_STRING: Error, your integer is not supported'
      END IF

      WRITE (integer_to_string, fmt=fmt) i

   END FUNCTION integer_to_string

   CHARACTER(len=8) FUNCTION real_to_string(x, pre_decimal, post_decimal)

      REAL, INTENT(IN) :: x
      INTEGER, INTENT(IN) :: pre_decimal
      INTEGER, INTENT(IN) :: post_decimal
      CHARACTER(len=8) :: fmt

      IF (pre_decimal < 0 .OR. post_decimal < 0) THEN
         STOP 'REAL_TO_STRING: Error, neither pre or post decimal should be negative'
      END IF

      IF (post_decimal == 0) THEN
         real_to_string = integer_to_string(int(x))
      ELSE
         fmt = '(F'//trim(integer_to_string(pre_decimal+post_decimal+1))//'.'//trim(integer_to_string(post_decimal))//')'
         WRITE (real_to_string, fmt=fmt) x
      END IF

   END FUNCTION real_to_string

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
