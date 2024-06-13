!
! meridional component of G flux 
! (G flux: eddy transport in meridional momentum equation)
!
subroutine gflux_y( v_zm, vv_zm, gy )
  use parameter, only : radius
  use com_var, only : jm, ko, rho, costbl
  implicit none
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: vv_zm(jm, ko)
  real(4),intent(out) :: gy(jm, ko)
  integer :: j, k

  do k=1, ko
     do j=1, jm
!        gy(j,k) = -rho(k) * ( vv_zm(j,k) - v_zm(j,k)**2 )
        gy(j,k) = -rho(k) * radius * ( vv_zm(j,k) - v_zm(j,k)**2 )
     end do
  end do

end subroutine gflux_y



!
! G flux divergence (meridional component)
!
subroutine gflux_div_y( gy, dgy )
  use parameter, only : radius
  use com_var, only : jm, ko, rho, alat, costbl
  implicit none
  real(4),intent(in)  :: gy(jm, ko)
  real(4),intent(out) :: dgy(jm, ko)

  real(4) :: temp(jm, ko)

  temp(:,:) = gy(:,:) * spread( costbl, 2, ko )

  call derivative_y( 1, jm, ko, alat, temp, dgy )

!  dgy(:,:) = dgy(:,:) / ( radius * spread(rho,1,jm) * spread(costbl,2,ko) )
  dgy(:,:) = dgy(:,:) / ( radius * radius * spread(rho,1,jm) * spread(costbl,2,ko) )
  dgy(1,:) = 0
  dgy(jm,:) = 0
  
end subroutine gflux_div_y



!
! vertical component of G flux 
! (G flux: eddy transport in meridional momentum equation)
!
subroutine gflux_z( pt_zm, v_zm, vv_zm, v_pt_dot_x_zm, &
     &              gz )
  use parameter, only : radius
  use com_var, only : im, jm, ko, rho
  implicit none
  real(4),intent(in)  :: pt_zm(jm, ko) 
  real(4),intent(in)  :: v_zm(jm, ko) 
  real(4),intent(in)  :: vv_zm(jm, ko) 
  real(4),intent(in)  :: v_pt_dot_x_zm(jm, ko)
  real(4),intent(out) :: gz(jm, ko)

  real(4) :: v_v_x_zm(jm, ko)

  real(4) :: gz_uv(jm, ko), gz_ut(jm, ko)

  v_v_x_zm(:,:) = vv_zm(:,:) - v_zm(:,:)**2

  call get_var_w_x_zm( pt_zm, v_v_x_zm, v_pt_dot_x_zm, &
       &               gz_uv, gz_ut, gz )

  gz(:,:) = -spread( rho, 1, jm ) * radius &
       &  * gz(:,:)

end subroutine gflux_z



!
! G flux divergence (vertical component)
!
subroutine gflux_div_z( gz, dgz )
  use parameter, only : radius
  use com_var, only : jm, ko, zd, rho, costbl
  implicit none
  real(4),intent(in)  :: gz(jm, ko)
  real(4),intent(out) :: dgz(jm, ko)

  call derivative_z( 1, jm, ko, zd, gz, &
       &             dgz )
  dgz(:,:) = dgz(:,:) &
       &   / ( spread( rho, 1, jm ) * radius )

  return
end subroutine gflux_div_z
