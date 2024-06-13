subroutine energy_conv_az_kz( alt, p_pds, v_zm, pt_zm, pt_sfc, phi_dagger, &
     &                        c_az_kz )
  use parameter, only : radius, grav, econv_min, econv_max
  use com_var, only : im, jm, ko, alat, pout
  implicit none
  real(4),intent(in)  :: alt(im, jm)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(in)  :: pt_sfc(im,jm)
  real(4),intent(in)  :: phi_dagger(jm, ko)
  real(4),intent(out) :: c_az_kz(jm, ko)

  real(4) :: phi_dagger_y(jm, ko)
  real(4) :: c_az_kz_modify(jm, ko)

  call derivative_y( 1, jm, ko, alat, phi_dagger, &
       &             phi_dagger_y )
  c_az_kz(:,:) = -v_zm(:,:) / radius * phi_dagger_y(:,:)
  
  ! Mountain modify (if necessary)
  call mount_modify( pt_sfc, alt, &
       &             v_zm, p_pds, pt_zm, phi_dagger, &
       &             c_az_kz_modify )
  c_az_kz(:,:) = c_az_kz_modify(:,:)

  call check_range( 1, jm, ko, c_az_kz, econv_min, econv_max, &
       &           'energy_conv_az_kz()', 'c_az_kz' )

end subroutine energy_conv_az_kz



subroutine energy_conv_kz_ae( u_zm, v_zm, depz_form, dz_dlat_zm, c_az_kz, &
     &                        c_kz_ae_u, c_kz_ae_v, c_kz_ae )
  use parameter, only : radius, grav, econv_min, econv_max
  use com_var, only : jm, ko
  implicit none
  real(4),intent(in)  :: u_zm(jm, ko)
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: depz_form(jm, ko)
  real(4),intent(in)  :: dz_dlat_zm(jm, ko)
  real(4),intent(in)  :: c_az_kz(jm, ko)
  real(4),intent(out) :: c_kz_ae_u(jm, ko)
  real(4),intent(out) :: c_kz_ae_v(jm, ko)
  real(4),intent(out) :: c_kz_ae(jm, ko)

  c_kz_ae_u(:,:) = -u_zm(:,:) * depz_form(:,:)

  c_kz_ae_v(:,:) = v_zm(:,:) / radius * grav * dz_dlat_zm(:,:) &
       &           + c_az_kz(:,:)

  c_kz_ae(:,:) = c_kz_ae_u(:,:) + c_kz_ae_v(:,:)

  call check_range( 1, jm, ko, c_kz_ae, econv_min, econv_max, &
       &            'energy_conv_kz_ae()', 'c_kz_ae' )
end subroutine energy_conv_kz_ae


subroutine energy_conv_ae_ke( v_zm, c_kz_ae_u, &
     &                        dz_dlat_zm, v_dz_dlat_zm, u_dz_dlon_zm, &
     &                        c_ae_ke_u, c_ae_ke_v, c_ae_ke )
  use parameter, only : radius, grav, econv_min, econv_max
  use com_var, only : jm, ko, costbl
  implicit none
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: c_kz_ae_u(jm, ko)
  real(4),intent(in)  :: dz_dlat_zm(jm, ko)
  real(4),intent(in)  :: v_dz_dlat_zm(jm, ko)
  real(4),intent(in)  :: u_dz_dlon_zm(jm, ko)
  real(4),intent(out) :: c_ae_ke_u(jm, ko)
  real(4),intent(out) :: c_ae_ke_v(jm, ko)
  real(4),intent(out) :: c_ae_ke(jm, ko)
  integer :: j, k

  do k=1, ko
     c_ae_ke_u(1,k) = 0
     do j=2, jm-1
        c_ae_ke_u(j,k) = -grav / radius / costbl(j) * u_dz_dlon_zm(j,k) &
             &           + c_kz_ae_u(j,k)
     end do
     c_ae_ke_u(jm,k) = 0
  end do
  
  c_ae_ke_v(:,:) = -grav / radius &
       &           * ( v_dz_dlat_zm(:,:) - v_zm(:,:) * dz_dlat_zm(:,:) )
  c_ae_ke(:,:) = c_ae_ke_u(:,:) + c_ae_ke_v(:,:)

  call check_range( 1, jm, ko, c_ae_ke, econv_min, econv_max, &
       &            'energy_conv_ae_ke()', 'c_ae_ke' )
end subroutine energy_conv_ae_ke




subroutine energy_conv_kz_ke( u_zm, v_zm, u_u_x_zm,  &
     &                        depy, depz_uw, dgy, dgz, &
     &                        c_kz_ke_uy, c_kz_ke_uz, &
     &                        c_kz_ke_vy, c_kz_ke_vz, c_kz_ke_tan, &
     &                        c_kz_ke )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, alat, costbl, tantbl
  implicit none
  real(4),intent(in)  :: u_zm(jm, ko)
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: u_u_x_zm(jm, ko)
  real(4),intent(in)  :: depy(jm, ko)
  real(4),intent(in)  :: depz_uw(jm, ko)
  real(4),intent(in)  :: dgy(jm, ko)
  real(4),intent(in)  :: dgz(jm, ko)
  real(4),intent(out) :: c_kz_ke_uy(jm, ko)
  real(4),intent(out) :: c_kz_ke_uz(jm, ko)
  real(4),intent(out) :: c_kz_ke_vy(jm, ko)
  real(4),intent(out) :: c_kz_ke_vz(jm, ko)
  real(4),intent(out) :: c_kz_ke_tan(jm, ko)
  real(4),intent(out) :: c_kz_ke(jm, ko)

  c_kz_ke_uy(:,:)  = -u_zm(:,:) * depy(:,:)

  c_kz_ke_uz(:,:)  = -u_zm(:,:) * depz_uw(:,:)

  c_kz_ke_vy(:,:)  = -v_zm(:,:) * dgy(:,:)

  c_kz_ke_vz(:,:)  = -v_zm(:,:) * dgz(:,:)

  c_kz_ke_tan(:,:) = v_zm(:,:) * u_u_x_zm(:,:) &
       &           / radius * spread( tantbl, 2, ko )
  c_kz_ke_tan(1,:) = 0.0
  c_kz_ke_tan(jm,:) = 0.0

  c_kz_ke(:,:) = c_kz_ke_uy(:,:) + c_kz_ke_uz(:,:) &
       &       + c_kz_ke_vy(:,:) + c_kz_ke_vz(:,:) &
       &       + c_kz_ke_tan(:,:)

  call check_range( 1, jm, ko, c_kz_ke, econv_min, econv_max, &
       &            'energy_conv_kz_ke()', 'c_kz_ke' )
end subroutine energy_conv_kz_ke
