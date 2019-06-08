PROGRAM string_operations_demo

  USE string_operations

  IMPLICIT NONE
  CHARACTER(len=256) :: part1, part2
  INTEGER :: i
!!$  CHARACTER(len=20) :: str(3)
!!$  INTEGER :: int(3), stat(3)

  part1='hello_'
  part2='.turd'

  DO i=1,100
     WRITE(*,*) i, TRIM(number_file(part1,i,part2))
     WRITE(*,*) i, TRIM(number_file_zeroes(part1,i,3,part2))
     WRITE(*,*)
  END DO

!!$  !!!String-to-integer test
!!$
!!$  str(1) = '123' ! Valid integer
!!$  str(2) = '-1'  ! Also valid
!!$  str(3) = 'one' ! invalid
!!$
!!$  call str2int(str,int,stat)
!!$
!!$  do i=1,3
!!$     write(*,*) int(i)
!!$    if ( stat(i) == 0 ) then
!!$      print *,i,int(i)
!!$    else
!!$      print *,'Conversion of string ',i,' failed!'
!!$    endif
!!$  enddo
!!$
!!$!!!

END PROGRAM string_operations_demo
