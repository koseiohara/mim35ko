!
! Function
!   interpolate/extrapolate undef data (pressure level only)
!
! Arguements (in)
!   im    : number of x-direction grid
!   jm    : number of y-direction grid
!   km    : number of z-direction grid
!   undef : undef value
!   pin   : pressure levels
!
! Arguements (inout)
!   var   : variable
!
subroutine undef_fill( im, jm, km, undef, pin, &
     &                 var )
  implicit none
  integer,intent(in)    :: im, jm, km
  real(4),intent(in)    :: undef
  real(4),intent(in)    :: pin(km)
  real(4),intent(inout) :: var(im, jm, km)

  integer :: i, j, k
  real(4) :: temp

  do k=3, km
     do j=1, jm
        do i=1, im

           if( var(i,j,k) == undef ) then
              call undef_hokan( pin(k), &
                   &            pin(k-1), var(i,j,k-1), &
                   &            pin(k-2), var(i,j,k-2), &
                   &            temp )
              var(i,j,k) = temp
           end if

        end do
     end do
  end do

end subroutine undef_fill



!
! Function
!   interpolate/extrapolate undef data
!
! Arguements (in)
!   p    : interpolated/extrapolated level
!   p1   : level (1) to be used in interpolation/extrapolation
!   A1   : value at p=p1
!   p2   : level (2) to be used in interpolation/extrapolation
!   A2   : value at p=p2
!
! Arguements (inout)
!   ret  : value at p
!
! Note
!   -log(p) linear interpolation/extrapolation is used.
!
subroutine undef_hokan( p, p1, A1, p2, A2, ret )
  implicit none
  real(4),intent(in)  :: p, p1, A1, p2, A2
  real(4),intent(out) :: ret

  ret = ( (A2-A1)*log(p) + A1*log(p2) - A2*log(p1) ) / (log(p2)-log(p1))
!  ret = ( (A2-A1)*p + A1*p2 - A2*p1 ) / (p2-p1)
end subroutine undef_hokan
