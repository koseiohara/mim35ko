!
! Function
!   check range of the variables
!
! Arguements (in)
!   im   : number of x-direction grid
!   jm   : number of y-direction grid
!   km   : number of z-direction grid
!   var  : variable to be checked
!   min  : possible minimum value of var
!   max  : possible maximum value of var
!   from : subroutine name (for message)
!   name : variable name (for message)
!
subroutine check_range( im, jm, km, var, min, max, from, name )
  use myglobal
  implicit none

  integer,intent(in) :: im, jm, km
  real(4),intent(in) :: var(im, jm, km)
  real(4),intent(in) :: min, max
  character(*),intent(in) :: from, name
  integer,save :: warn_count = 1
  integer,parameter :: max_warn_count = 10000
  
  integer :: i, j, k
  
  do k=1, km
     do j=1, jm
        do i=1, im

           if( var(i,j,k) < min .or. var(i,j,k) > max ) then
              write(warn_unit,*) 'warning in check_range() : ' 
              write(warn_unit,*) trim(name) // ' should be between ', min, 'and', max
              write(warn_unit,*) '(i,j,k) = ', i, j, k
              write(warn_unit,*) 'var = ', var(i,j,k)
              write(warn_unit,*) 'called from ' // trim(from)
              write(warn_unit,*) 'warning count = ', warn_count
              
              if( warn_count > max_warn_count ) then
                 write(warn_unit,*) 'warning count exceeds maximum number'
                 stop
              end if
              warn_count = warn_count + 1

           end if

        end do
     end do
  end do
  
end subroutine check_range
