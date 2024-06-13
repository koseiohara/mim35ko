!
! EP flux divergence (meridional component)
!
subroutine epflux_div_y( epy, depy )
  use parameter, only : radius, divf_min, divf_max
  use com_var, only : jm, ko, costbl, alat, rho
  implicit none
  real(4),intent(in)  :: epy(jm, ko)
  real(4),intent(out) :: depy(jm, ko)

  integer :: j, k

  do k=1, ko

     depy(1,k) = 0.0

     do j=2, jm-1
        depy(j,k) = 1.0 / ( rho(k) * (radius*costbl(j))**2 ) &
             &              * ( ( epy(j+1,k) * costbl(j+1) ) &
             &                - ( epy(j-1,k) * costbl(j-1) ) ) &
             &    / ( alat(j+1) - alat(j-1) )
     end do

     depy(jm,k) = 0.0

  end do

  call check_range( 1, jm, ko, depy, divf_min, divf_max, &
       &            'epflux_div_y()', 'depy' )

  return
end subroutine epflux_div_y



!
! EP flux divergence (vertical  component)
!
subroutine epflux_div_z( epz, depz )
  use parameter, only : radius, divf_min, divf_max
  use com_var, only : jm, ko, costbl, rho, zd
  implicit none
  real(4),intent(in) :: epz(jm, ko)
  real(4),intent(out) :: depz(jm, ko)

  call derivative_z( 1, jm, ko, zd, epz, &
       &             depz )
  depz(:,:) = depz(:,:) &
       &    / ( spread(rho,1,jm) * radius * spread(costbl,2,ko) )

  depz(1,:) = 0.0
  depz(jm,:) = 0.0

  call check_range( 1, jm, ko, depz, divf_min, divf_max, &
       &            'epflux_div_z()', '(one of the) depz' )

  return
end subroutine epflux_div_z
