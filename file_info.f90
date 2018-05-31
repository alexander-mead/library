MODULE file_info

CONTAINS

  INTEGER FUNCTION file_length(file_name)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: file_name
    INTEGER :: n

    WRITE(*,*) 'FILE_LENGTH: File: ', TRIM(file_name)
    OPEN(7,file=file_name)

    !Newer version that lacks 'data' seems okay
    n=0
    DO
       n=n+1
       READ(7,*, end=301)
    END DO

    !301 is just the label to jump to when the end of the file is reached

301 CLOSE(7)

    file_length=n-1

    WRITE(*,*) 'FILE_LENGTH: Length:', file_length
    WRITE(*,*)

  END FUNCTION file_length

  !Tilman's version of file_length (nice because no GOTO)
  FUNCTION count_number_of_lines(filename)result(n)
    
    IMPLICIT NONE
    CHARACTER(len=*), intent(in) :: filename
    INTEGER :: n, file_unit, iostat
    CHARACTER :: c

    OPEN(newunit=file_unit, file=filename, status='old', iostat=iostat)
    IF(iostat > 0) THEN
        PRINT *, "Failed to open file ", filename
        n = -1
        RETURN
    END IF

    n = 0
    DO
        READ(unit=file_unit, fmt=*, iostat=iostat) c
        IF(iostat == 0) THEN
            IF(c == '#') THEN
                print *, "Detected comment lines in file ", filename, ". Check your data file."
                CYCLE
            END IF

            n = n + 1
        ELSE IF(iostat < 0) THEN
            !print *, "Reached end of file."
            EXIT
        ELSE
            PRINT *, "Error reading file ", filename
            n = -1
            RETURN
        END IF
    END DO

    CLOSE(file_unit)
    
  END FUNCTION count_number_of_lines

END MODULE file_info

