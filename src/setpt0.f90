!
! Function
!   get potential temperature at the pressure levels
!
! Arguements (in)
!   t      : temperature [K]
!   p_sfcs : surface pressure [hPa]
!
! Arguements (out)
!   pt     : potential temperature [K]
!   pt_sfc : surface potential temperature [K]
!
! Note
!   -If the atmosphere is unstable, potential temperature profile is adjusted
!    (not so appropriate until now)
!   -Potential temperature under the ground is set to its surface value
!
subroutine setpt0( t, p_sfc, pt, pt_sfc )
  use parameter, only : rkappa, pt_max, pt_min
  use com_var, only : im, jm, km, pin
  implicit none
  real(4),intent(in)  :: t(im, jm, km)
  real(4),intent(in)  :: p_sfc(im, jm)
  real(4),intent(out) :: pt(im, jm, km)
  real(4),intent(out) :: pt_sfc(im, jm)

  integer :: i, j, k, l
  real(4) :: d

  ! temperature -> potential temperature
  do k=1, km
     do j=1, jm
        do i=1, im
           pt(i,j,k) = t(i,j,k) * ( 1000.0 / pin(k) )**rkappa
        end do
     end do
  end do
  
  ! if atmosphere is instable, modify pt by extrapolation
  ! Note: energy conservation is NOT satisfied until now
  !       (it should be modified later...)
  do k=2, km
     do j=1, jm
        do i=1, im
           if( k > 1 .and. pt(i,j,k) >= pt(i,j,k-1) ) then
              pt(i,j,k) = pt(i,j,k-1) - 1.0
           end if
        end do
     end do
  end do
  
  ! get pt_sfc
  do i=1, im
     do j=1, jm

        ! search for the level which is nearest to the ground
        l = 1
        do while( pin(l+1) < p_sfc(i,j) .and. l <= km-2 )
           l = l + 1
        end do
        
        ! interpolate to ground level
        d = ( p_sfc(i,j) - pin(l) ) / ( pin(l+1) - pin(l) )
        pt_sfc(i,j) = pt(i,j,l) * (1.0-d) + pt(i,j,l+1) * d
        
        if( pt_sfc(i,j) < 0 )then
           write(6,*) "warning(setpt0): pt_sfc<0", &
                &    i, j, l, pt_sfc(i,j), pt(i,j,l), pt(i,j,l+1)
        end if

        ! pt = pt_sfc under the ground
        do k=1, km
           if( pin(k) > p_sfc(i,j) ) pt(i,j,k) = pt_sfc(i,j)        
        end do
        
     end do
  end do
  
  ! check value
  call check_range( im, jm, km, pt, pt_min, pt_max, 'setpt0()', 'pt' )
  call check_range( im, jm, 1, pt_sfc, pt_min, pt_max, 'setpt0()', 'pt_sfc' )

  return
end subroutine setpt0
