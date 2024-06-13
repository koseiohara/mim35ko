!
! Vertically integrated A_Z
!
! Note:
!   It is not used for global mean A_Z (because of the accuracy)
!
subroutine energy_az_vint( p_pds, p_pdds, pd_pdd, pt_pdds, pt_ym, &
     &                     az_zm_vint )
  use parameter, only : rkappa, grav, cp, pai
  use com_var, only : jm, ko, alat, pout
  implicit none
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: p_pdds(1)
  real(4),intent(in)  :: pd_pdd(jm, ko)
  real(4),intent(in)  :: pt_pdds(1)
  real(4),intent(in)  :: pt_ym(ko)
  real(4),intent(out) :: az_zm_vint(jm)

  real(4) :: pd_ym(ko)
  real(4) :: integ(jm, ko)
  integer :: j, k
  real(4) :: const
  real(4) :: az_modify(jm)

  ! pd_ym : global mean p+ at the p++ levels
  !         pd_ym must be almost equal to standard p++ levels 
  !         except under the ground.
  call integral_meridional( 1, jm, ko, alat, pd_pdd, &
       &                    pd_ym )
  
  ! get integrand
  const = cp * (1.0e+5)**(-rkappa) / (1+rkappa) / grav
  do k=1, ko
     do j=1, jm
        integ(j,k) = const * ( ( pd_pdd(j,k)*100 )**(1+rkappa) &
             &                - ( pd_ym(k)*100)**(1+rkappa) )
     end do
  end do

  ! integrate with pt
  call integral_pt_ym( jm, ko, pout, p_pdds, pt_ym, pt_pdds, integ, &
       &               az_zm_vint )

  ! lower boundary modification
  !   it is proportional to pt_ymin * ( p+s^{kappa+1} - p++s^{kappa+1} ).
  az_modify(:) = const * pt_pdds(1) &
       &       * ( p_pds(:)**(rkappa+1) - p_pdds(1)**(rkappa+1) )

  az_zm_vint(:) = az_zm_vint(:) + az_modify(:)

end subroutine energy_az_vint



!
! Global mean A_Z
!
! Note:
!   order of the integration is different from energy_az_vint()
!
subroutine energy_az_gmean( p_pds, p_pdds, pd_pdd, pd_ym, &
     &                      pt_pdds, pt_ym, &
     &                      az_gmean )
  use parameter, only : rkappa, grav, cp, pai
  use com_var, only : jm, ko, alat, costbl, pout
  implicit none
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: p_pdds(1)
  real(4),intent(in)  :: pd_pdd(jm, ko)
  real(4),intent(in)  :: pd_ym(ko)
  real(4),intent(in)  :: pt_pdds(1)
  real(4),intent(in)  :: pt_ym(ko)
  real(4),intent(out) :: az_gmean(1)

  real(4) :: integ(jm, ko), integ_temp(ko)
  integer :: j, k
  real(4) :: temp1(jm), az_modify(1)
  
  real(4) :: const

  ! get integrand
  const = cp * (1.0e+5)**(-rkappa) / (1+rkappa) / grav
  do k=1, ko
     do j=1, jm
        integ(j,k) = const * ( ( pd_pdd(j,k)*100 )**( 1+rkappa ) &
             &               - ( pd_ym(k)*100 )**( 1+rkappa ) )
     end do
  end do
  
  ! meridional mean
  call integral_meridional( 1, jm, ko, alat, integ, &
       &                    integ_temp )


  ! integrate with pt
  call integral_pt_ym( 1, ko, pout, p_pdds, pt_ym, pt_pdds, integ_temp, &
       &               az_gmean )

  ! lower boundary modification
  !   it is proportional to pt_ymin * ( p+s^{kappa+1} - p++s^{kappa+1} ).
  temp1(:) = const * pt_pdds(1)  &
       &   * ( p_pds(:)**(rkappa+1) - p_pdds(1)**(rkappa+1) )
  call integral_meridional( 1, jm, 1, alat, temp1, &
       &                    az_modify )
     
  az_gmean(1) = az_gmean(1) + az_modify(1)
  
end subroutine energy_az_gmean



!
! 2-dimensional A_E
!   = P - Pz
!   = const * pt d/d(p+) [ p^(kappa+1) - p+^(kappa+1) ]
!
! Note:
!   It is not used for estimate of vertically integrated A_E 
!     since mountain effects have not been neglected here.
!
subroutine energy_ae_total( p_pt, p_zm, p_pds, pt_zm, &
     &                      ae_total_zm )
  use parameter, only : rkappa, grav, cp
  use com_var, only : im, jm, km, ko, pout, pin
  use biseki, only : biseki_biseki
  implicit none
  real(4),intent(in)  :: p_pt(im, jm, ko)
  real(4),intent(in)  :: p_zm(jm, ko)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(out) :: ae_total_zm(jm, ko)

  real(4) :: p_diff(jm, ko)
  real(4) :: const

  real(4) :: p_kappa(im, jm, ko)
  real(4) :: p_kappa_zm(jm, ko)
  integer :: i, j, k

  ! pending
!  do k=1, km
!     do j=1, jm
!        do i=1, im
!           p_kappa(i,j,k) = (100*pin(k))**rkappa
!        end do
!     end do
!  end do
!  call biseki_biseki( p_kappa, p_kappa_zm )
!
!  do k=1, ko
!     do j=1, jm
!        ae_total_zm(j,k) = ( p_kappa_zm(j,k) - (100*pout(k))**rkappa ) &
!             &           * cp * pt_zm(j,k) * ( 1.0e+5**(-rkappa) )
!     end do
!  end do
!
!  return


  const = cp * (1.0e+5)**(-rkappa) / (1+rkappa)
  
  p_diff(:,:) = sum( (p_pt*100)**(rkappa+1), dim=1 ) / real(im) &
       &      - ( p_zm(:,:)*100 )**(rkappa+1)

  call derivative_p( 1, jm, ko, pout*100, p_pds*100, p_diff, &
       &             ae_total_zm )

  ae_total_zm(:,:) = ae_total_zm(:,:) * pt_zm(:,:) * const

end subroutine energy_ae_total



!
! Vettically integrated A_E
!
! ae_zm_vint : proportional to int[ (p^(rkappa+1))_zm - pd^(rkappa+1) ] d(pt)
!
subroutine energy_ae_vint( p_pt, p_zm, p_sfc, p_pds, pt_zm, pt_pds, &
     &                     ae_zm_vint )
  use parameter, only : rkappa, grav, cp
  use com_var, only : im, jm, ko, pout
  implicit none
  real(4),intent(in)  :: p_pt(im, jm, ko)
  real(4),intent(in)  :: p_zm(jm, ko)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(in)  :: pt_pds(jm)
  real(4),intent(in)  :: p_sfc(im, jm)
  real(4),intent(out) :: ae_zm_vint(jm)

  real(4) :: p_kappa_p1_zm(jm, ko)
  real(4) :: p_diff(jm, ko)
  real(4) :: const
  integer :: i, j
  real(4) :: temp1(im, jm), ae_modify(jm)

  const = cp * (1.0e+5)**(-rkappa) / (1+rkappa) / grav


  ! Note: even if Taylor expansion is used, results are similar to the below
  do j=1, jm

     ! zonal mean p^(1+rkappa)
     p_kappa_p1_zm(j,1:ko) = 0
     do i=1, im
        p_kappa_p1_zm(j,1:ko) = p_kappa_p1_zm(j,1:ko) &
             &                + ( p_pt(i,j,1:ko)*100 )**(1+rkappa)
     end do
     p_kappa_p1_zm(j,1:ko) = p_kappa_p1_zm(j,1:ko) / real(im)

     ! proportional to { p^(1+rkappa) }_zm - pd^(1+rkappa)
     p_diff(j,1:ko) = const &
          &         * ( p_kappa_p1_zm(j,1:ko) &
          &           - ( p_zm(j,1:ko)*100 )**(1+rkappa) )
  end do

  ! integrate with pt
  call integral_pt( jm, ko, pout, p_pds, pt_zm, pt_pds, p_diff, &
       &            ae_zm_vint )

  ! lower boundary modification
  !   it is proportional to pt_ymin * ( ps^{kappa+1} - p+s^{kappa+1} ).
  temp1(:,:) = const * spread(pt_pds,1,im) &
       &     *  ( p_sfc(:,:)**(rkappa+1) - spread(p_pds,1,im)**(rkappa+1) )
  ae_modify(:) = sum( temp1, dim=1 ) / real(im)
  ae_zm_vint = ae_zm_vint + ae_modify

end subroutine energy_ae_vint
