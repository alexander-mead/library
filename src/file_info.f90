MODULE file_info

   USE basic_operations

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: file_length
   PUBLIC :: file_columns
   PUBLIC :: count_number_of_lines
   PUBLIC :: check_file_exists
   PUBLIC :: file_exists

CONTAINS

   INTEGER FUNCTION file_length(file_name, verbose)

      ! Get the number of lines in the file 
      ! TODO: Remove GOTO-ish thing   
      CHARACTER(len=*), INTENT(IN) :: file_name
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: n, u
      LOGICAL :: lexist

      IF (present_and_correct(verbose)) WRITE (*, *) 'FILE_LENGTH: File: ', trim(file_name)
      INQUIRE (file=file_name, exist=lexist)
      IF (.NOT. lexist) THEN
         WRITE (*, *) 'FILE_LENGTH: File: ', trim(file_name)
         STOP 'FILE_LENGTH: Error, file does not exist'
      END IF
   
      n = 0
      OPEN (newunit=u, file=file_name, status='old')
      DO
         n = n+1
         READ (u, *, end=301) ! Horrid go to
      END DO    
301   CLOSE (u) ! 301 is just the label to jump to when the end of the file is reached
      file_length = n-1

      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'FILE_LENGTH: Length:', file_length
         WRITE (*, *)
      END IF

   END FUNCTION file_length

   INTEGER FUNCTION file_columns(infile, verbose)

      ! Count the number of columns in a file
      ! Based on: https://stackoverflow.com/questions/7314216/reading-data-file-in-fortran-with-known-number-of-lines-but-unknown-number-of-en
      CHARACTER(len=*) :: infile
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose
      INTEGER :: i, error, u
      !CHARACTER(len=16) :: crap
      CHARACTER(len=1024) :: line
      CHARACTER(len=16), ALLOCATABLE :: array(:)
      LOGICAL :: lexist
      !INTEGER :: file_rows
      INTEGER, PARAMETER :: max_columns = 1000 ! Maximum possible number of entries

      ! Check the file exists
      IF (present_and_correct(verbose)) WRITE (*, *) 'FILE_COLUMNS: File: ', trim(infile)
      INQUIRE (file=infile, exist=lexist)
      IF (.NOT. lexist) THEN
         WRITE (*, *) 'FILE_COLUMNS: File: ', trim(infile)
         STOP 'FILE_COLUMNS: Error, file does not exist'
      END IF
      
      ! Read in first line as string
      OPEN(newunit=u, file=infile, status='old')
      READ(u, '(A)') line
      CLOSE(u)

      ! Analyse string for how many bits it can be split in to
      ALLOCATE(array(max_columns))
      DO i = 1, max_columns       
         READ(line, *, iostat=error) array(1:i)
         IF (error .NE. 0) THEN
            file_columns = i-1
            EXIT
         END IF        
      END DO

      ! Write to screen
      IF (present_and_correct(verbose)) THEN
         WRITE (*, *) 'FILE_COLUMNS: Columns:', file_columns
         WRITE (*, *)
      END IF

   END FUNCTION file_columns

   FUNCTION count_number_of_lines(filename) result(n)

      ! Tilman's version of file_length (nice because no GOTO)
      CHARACTER(len=*), INTENT(IN) :: filename
      INTEGER :: n
      INTEGER :: file_unit, iostat
      CHARACTER :: c

      OPEN (newunit=file_unit, file=filename, status='old', iostat=iostat)
      IF (iostat > 0) THEN
         PRINT *, "Failed to open file ", filename
         n = -1
         RETURN
      END IF

      n = 0
      DO
         READ (unit=file_unit, fmt=*, iostat=iostat) c
         IF (iostat == 0) THEN
            IF (c == '#') THEN
               print *, "Detected comment lines in file ", filename, ". Check your data file."
               CYCLE
            END IF

            n = n+1
         ELSE IF (iostat < 0) THEN
            !print *, "Reached end of file."
            EXIT
         ELSE
            PRINT *, "Error reading file ", filename
            n = -1
            RETURN
         END IF
      END DO

      CLOSE (file_unit)

   END FUNCTION count_number_of_lines

   SUBROUTINE check_file_exists(file)

      CHARACTER(len=*), INTENT(IN) :: file

      ! Check file exists
      IF (.NOT. file_exists(file)) THEN
         WRITE (*, *) 'CHECK_FILE_EXISTS: File: ', trim(file)
         STOP 'CHECK_FILE_EXISTS: File does not exist'
      END IF

   END SUBROUTINE check_file_exists

   LOGICAL FUNCTION file_exists(file)

      ! Returns true if a file exists, or false otherwise
      CHARACTER(len=*), INTENT(IN) :: file

      INQUIRE (file=file, exist=file_exists)

   END FUNCTION file_exists

END MODULE file_info

