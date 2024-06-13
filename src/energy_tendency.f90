! advection by v
! -1 / [ a cos(phi) ] * d/d(phi) [ Kz v cos(phi) ]
subroutine energy_tendency_dkzdt_vkz( v_zm, kz_zm, dkzdt_vke )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, alat, costbl
  implicit none
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: kz_zm(jm, ko)
  real(4),intent(out) :: dkzdt_vke(jm, ko)
  real(4) :: temp(jm, ko)

  temp(:,:) = v_zm(:,:) * kz_zm(:,:) * spread( costbl, 2, ko )

  call derivative_y( 1, jm, ko, alat, temp, &
       &             dkzdt_vke )
  dkzdt_vke(:,:) = -dkzdt_vke(:,:) / ( radius * spread(costbl,2,ko) )
  dkzdt_vke(1,:) = 0   ! suppress sivergence
  dkzdt_vke(jm,:) = 0  ! suppress sivergence

  call check_range( 1, jm, ko, dkzdt_vke, econv_min, econv_max, &
       &            'energy_tendency_dkzdt_vkz()', 'dkzdt_vke' )

end subroutine energy_tendency_dkzdt_vkz


! advection by w+
! -1/rho * d/d(z+) ( Kz w )
subroutine energy_tendency_dkzdt_wkz( w_zm, kz_zm, dkzdt_wkz )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, alat, costbl, rho, zd
  implicit none
  real(4),intent(in)  :: w_zm(jm, ko)
  real(4),intent(in)  :: kz_zm(jm, ko)
  real(4),intent(out) :: dkzdt_wkz(jm, ko)

  real(4) :: temp(jm, ko)

  temp(:,:) = spread( rho, 1, jm ) * w_zm(:,:) * kz_zm(:,:)

  call derivative_z( 1, jm, ko, zd, temp, &
       &             dkzdt_wkz )

  dkzdt_wkz(:,:) = -dkzdt_wkz(:,:) / spread( rho, 1, jm )

  call check_range( 1, jm, ko, dkzdt_wkz, econv_min, econv_max, &
       &            'energy_tendency_dkzdt_wkz()', 'dkzdt_wkz' )

end subroutine energy_tendency_dkzdt_wkz






! advection by v
! -1/ [ a cos(phi) ] * d/d(phi) [ Ke v cos(phi) ]
subroutine energy_tendency_dkedt_vke( v_ke_zm, dkedt_vke )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, alat, costbl
  implicit none
  real(4),intent(in)  :: v_ke_zm(jm, ko)
  real(4),intent(out) :: dkedt_vke(jm, ko)
  real(4) :: temp(jm, ko)

  temp(:,:) = v_ke_zm(:,:) * spread( costbl, 2, ko )

  call derivative_y( 1, jm, ko, alat, temp, &
       &             dkedt_vke )
  dkedt_vke(:,:) = -dkedt_vke(:,:) / ( radius * spread(costbl,2,ko) )
  dkedt_vke(1,:) = 0   ! suppress sivergence
  dkedt_vke(jm,:) = 0  ! suppress sivergence

  call check_range( 1, jm, ko, dkedt_vke, econv_min, econv_max, &
       &            'energy_tendency_dkedt_vke()', 'dkedt_vke' )

end subroutine energy_tendency_dkedt_vke



! advection by w+
! -1/rho * d/d(z+) ( Ke w )
subroutine energy_tendency_dkedt_wke( pt_zm, v_ke_zm, pt_dot_ke_zm, &
     &                                dkedt_wke )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, alat, costbl, rho, zd
  implicit none
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(in)  :: v_ke_zm(jm, ko)
  real(4),intent(in)  :: pt_dot_ke_zm(jm, ko)
  real(4),intent(out) :: dkedt_wke(jm, ko)

  real(4) :: w_ke_zm(jm, ko)
  real(4) :: dummy1(jm, ko), dummy2(jm, ko)

  call get_var_w_x_zm( pt_zm, v_ke_zm, pt_dot_ke_zm, &
       &               dummy1, dummy2, w_ke_zm )
  w_ke_zm(:,:) = w_ke_zm(:,:) * spread( rho, 1, jm )

  call derivative_z( 1, jm, ko, zd, w_ke_zm, &
       &             dkedt_wke )

  dkedt_wke(:,:) = -dkedt_wke(:,:) / spread( rho, 1, jm )

  dkedt_wke(1,:) = 0.0
  dkedt_wke(jm,:) = 0.0

  call check_range( 1, jm, ko, dkedt_wke, econv_min, econv_max, &
       &            'energy_tendency_dkedt_wke()', 'dkedt_wke' )

end subroutine energy_tendency_dkedt_wke



! 1 / [ rho a^2 cos(phi) ] * d/d(phi) ( u * Fy )
subroutine energy_tendency_dkedt_uy( u_zm, epy, dkedt_uy )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, alat, costbl, rho
  implicit none
  real(4),intent(in)  :: u_zm(jm, ko)
  real(4),intent(in)  :: epy(jm, ko)
  real(4),intent(out) :: dkedt_uy(jm, ko)

  real(4) :: temp(jm, ko)

  temp = u_zm(:,:) * epy(:,:)

  call derivative_y( 1, jm, ko, alat, temp, &
       &             dkedt_uy )

  dkedt_uy(:,:) = dkedt_uy(:,:) &
       &        / ( spread(rho,1,jm) * radius * radius * spread(costbl,2,ko) )
  dkedt_uy(1,:) = 0   ! suppress sivergence
  dkedt_uy(jm,:) = 0  ! suppress sivergence

  call check_range( 1, jm, ko, dkedt_uy, econv_min, econv_max, &
       &            'energy_tendency_dkedt_uy()', 'dkedt_uy' )

end subroutine energy_tendency_dkedt_uy


subroutine energy_tendency_dkedt_vy( v_zm, gy, dkedt_vy )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, alat, costbl, rho
  implicit none
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: gy(jm, ko)
  real(4),intent(out) :: dkedt_vy(jm, ko)

  real(4) :: temp(jm, ko)

  temp(:,:) = v_zm(:,:) * gy(:,:) * spread( costbl, 2, ko )

  call derivative_y( 1, jm, ko, alat, temp, &
       &             dkedt_vy )

  dkedt_vy(:,:) = dkedt_vy(:,:) &
       &        / ( spread(rho,1,jm) * radius * radius * spread(costbl,2,ko) )
!       &        / ( spread(rho,1,jm) * radius * spread(costbl,2,ko) )
  dkedt_vy(1,:) = 0   ! suppress sivergence
  dkedt_vy(jm,:) = 0  ! suppress sivergence

  call check_range( 1, jm, ko, dkedt_vy, econv_min, econv_max, &
       &            'energy_tendency_dkedt_vy()', 'dkedt_vy' )

end subroutine energy_tendency_dkedt_vy



subroutine energy_tendency_dkedt_uz( u_zm, epz_uw, dkedt_uz )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, costbl, rho, zd
  implicit none
  real(4),intent(in)  :: u_zm(jm, ko)
  real(4),intent(in)  :: epz_uw(jm, ko)
  real(4),intent(out) :: dkedt_uz(jm, ko)

  real(4) :: temp(jm, ko)

  temp = u_zm(:,:) * epz_uw(:,:)

  call derivative_z( 1, jm, ko, zd, temp, &
       &             dkedt_uz )

  dkedt_uz(:,:) = dkedt_uz(:,:) &
       &        / ( spread(rho,1,jm) * radius * spread(costbl,2,ko) )

  dkedt_uz(1,:) = 0.0
  dkedt_uz(jm,:) = 0.0

  call check_range( 1, jm, ko, dkedt_uz, econv_min, econv_max, &
       &            'energy_tendency_dkedt_uz()', 'dkedt_uz' )

end subroutine energy_tendency_dkedt_uz



subroutine energy_tendency_dkedt_vz( v_zm, gz, dkedt_vz )
  use parameter, only : radius, econv_min, econv_max
  use com_var, only : jm, ko, rho, zd
  implicit none
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: gz(jm, ko)
  real(4),intent(out) :: dkedt_vz(jm, ko)

  real(4) :: temp(jm, ko)

  temp = v_zm(:,:) * gz(:,:)

  call derivative_z( 1, jm, ko, zd, temp, &
       &             dkedt_vz )

  dkedt_vz(:,:) = dkedt_vz(:,:) &
       &        / ( spread(rho,1,jm) * radius )

  dkedt_vz(1,:) = 0.0
  dkedt_vz(jm,:) = 0.0

  call check_range( 1, jm, ko, dkedt_vz, econv_min, econv_max, &
       &            'energy_tendency_dkedt_vz()', 'dkedt_vz' )

end subroutine energy_tendency_dkedt_vz



!subroutine energy_tendency_daedt_vae()
!  use mim_var  ! temporal
!  use parameter, only : rkappa, grav, cp, radius
!  use com_var, only : im, jm, km, ko, pout, costbl, alat
!  implicit none
!  integer :: i, j, k, h
!  real(4) :: w
!  integer :: kl, ku
!  real(4) :: v_pt(im, jm, ko)
!
!  real(4) :: pdiff_zm(jm, ko)
!  real(4) :: temp_zm(jm, ko)
!  real(4) :: v_ae_vint(jm)
!  real(4) :: daedt_vae_vint(jm)
!  real(4) :: const
!
!  const = cp * (100000.0)**(-rkappa) / (1+rkappa) / grav
!
!
!  ! v(p) -> v(pd) interpolate with pt
!  do i=1, im
!     do j=1, jm
!        do k=1, ko
!
!           if( pt_zm(j,k) < pt(i,j,ko) ) then
!              !write(*,*) i, j, k
!              v_pt(i,j,k) = v_pt(i,j,k-1)  ! approx. surface value
!              cycle
!           end if
!
!           if( pt_zm(j,k) > pt(i,j,1) ) then
!              kl = 1
!              ku = 2
!           else
!              do h=1, ko-1
!                 if( pt(i,j,h+1) <= pt_zm(j,k) &
!                      &       .and. pt_zm(j,k) <= pt(i,j,h) ) then
!                    kl = h
!                    ku = h + 1
!                    exit
!                 end if
!              end do
!           end if
!
!
!           w = ( pt_zm(j,k) - pt(i,j,kl) ) / ( pt(i,j,ku) - pt(i,j,kl) )
!           v_pt(i,j,k) = w * v(i,j,ku) + (1-w) * v(i,j,kl)
!
!!           write(*,*) i, j, k, pt(i,j,kl), pt_zm(j,k), pt(i,j,ku)
!!           write(*,*) '       ->', v(i,j,kl), v_pt(i,j,k), v(i,j,ku)
!
!        end do
!     end do
!  end do
!
!
!
!  ! integrand
!  do k=1, ko
!     do j=1, jm
!
!        pdiff_zm(j,k) = 0.0
!        do i=1, im
!           pdiff_zm(j,k) = pdiff_zm(j,k) &
!                &        + ( (p_pd(i,j,k)*100.0)**(1+rkappa) &
!                &            - (p_zm(j,k)*100.0)**(1+rkappa) ) &
!                &        * v_pt(i,j,k)
!        end do
!        pdiff_zm(j,k) = pdiff_zm(j,k) / real(im) * const
!
!!        write(*,*) j, k, pdiff_zm(j,k)
!     end do
!  end do
!
!!  call derivative_y(1, jm, ko, alat, pdiff_zm, &
!!       &            temp_zm)
!
!
!  call integral_pt(im, jm, ko, pout, p_pds, pt_sfc, pt_zm, pdiff_zm, &
!       &           v_ae_vint)
!!  call integral_pt(im, jm, ko, pout, p_pds, pt_sfc, pt_zm, temp_zm, &
!!       &           daedt_vae_vint)
!
!!  do j=1, jm
!!     write(*,*) j, ae_zm_vint(j), v_ae_vint(j)
!!  end do
!  v_ae_vint(:) = v_ae_vint(:) * costbl(:)
!
!
!
!  call derivative_y(1, jm, 1, alat, v_ae_vint, &
!       &            daedt_vae_vint)
!  daedt_vae_vint(:) = -daedt_vae_vint(:) / ( radius * costbl(:) )
!  daedt_vae_vint(1) = 0   ! suppress sivergence
!  daedt_vae_vint(jm) = 0  ! suppress sivergence
!
!
!  do j=1, jm
!     write(*,*) j, daedt_vae_vint(j)
!  end do
!
!
!end subroutine energy_tendency_daedt_vae






!
!
!
!subroutine energy_tendency_dpedt_vt(t, t_dagger, pt, v, v_zm, z, phi_dagger, dpedt_vt)
!  use parameter, only : radius, gasr, cp, grav, econv_min, econv_max
!  use com_var, only : im, jm, km, ko, costbl, alat
!  use biseki, only : biseki_biseki
!  implicit none
!  real(4),intent(in)  :: t(im, jm, km)
!  real(4),intent(in)  :: t_dagger(jm, ko)
!  real(4),intent(in)  :: pt(im, jm, km)
!  real(4),intent(in)  :: v(im, jm, km)
!  real(4),intent(in)  :: v_zm(jm, ko)
!  real(4),intent(in)  :: z(im, jm, km)
!  real(4),intent(in)  :: phi_dagger(jm, ko)
!  real(4),intent(out) :: dpedt_vt(jm, ko)
!
!  real(4) :: t_v(im, jm, km)
!  real(4) :: t_v_zm(jm, ko)
!  real(4) :: temp(jm, ko)
!  integer :: j, k
!
!  real(4) :: p_energy(im, jm, km)
!  real(4) :: p_energy_zm(jm, ko)
!  real(4) :: pe(im, jm, km)
!  real(4) :: pe_zm(jm, ko)







!  t_v(:,:,:) = t(:,:,:) * v(:,:,:)
!  call biseki_biseki(t_v, t_v_zm)
!
!  temp(:,:) = gasr * ( t_v_zm(:,:) - t_dagger(:,:) * v_zm(:,:) ) &
!       &    * spread(costbl,2,ko)
!
!  call derivative_y(1, jm, ko, alat, temp, &
!       &            dpedt_vt)
!
!  dpedt_vt(:,:) = -dpedt_vt(:,:) / ( radius * spread(costbl,2,ko) )
!
!  dpedt_vt(1,:) = 0.0
!  dpedt_vt(jm,:) = 0.0

!  call check_range(1, jm, ko, dpedt_vt, econv_min, econv_max, &
!       &           'energy_tendency_dpedt_vt()', 'dpedt_vt')


!  pe(:,:,:) = pt(:,:,:) / ( 100000.0**rkappa ) &
!       &    * ( rkappa *   )
!
!  p_energy(:,:,:) = ( cp - gasr ) * t(:,:,:) + z(:,:,:) * grav
!
!  call biseki_biseki(p_energy, p_energy_zm)
!
!  pe_zm(:,:) = p_energy_zm(:,:) - ( cp - gasr ) * t_dagger(:,:) - phi_dagger(:,:)
!
!  dpedt_vt(:,:) = pe_zm(:,:) ! dummy

!  do k=1, ko
!     write(*,*) k, p_energy_zm(10,k), pe_zm(10,k)
!  end do
!  dpedt_vt(:,:) = gasr * ( t_v_zm(:,:) - t_dagger(:,:) * v_zm(:,:) )
!  dpedt_vt(:,:) = gasr * t_v_zm(:,:)

!end subroutine energy_tendency_dpedt_vt
