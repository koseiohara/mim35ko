subroutine mount_modify( pt_sfc, alt, v_zm, p_pds, pt_zm, phi_dagger, &
     &                   c_az_kz_modify )
  use parameter, only : h0, grav, radius, cp, rkappa, gasr
  use com_var, only : im, jm, ko, alat, pout
  implicit none
  real(4),intent(in)  :: pt_sfc(im, jm)
  real(4),intent(in)  :: alt(im, jm)
  real(4),intent(in)  :: v_zm(jm, ko)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(in)  :: phi_dagger(jm, ko)
  real(4),intent(out) :: c_az_kz_modify(jm, ko)

  real(4) :: phi_dagger_modify(jm, ko)
  real(4) :: pout_modify(jm, ko)
  real(4) :: v_zm_modify(jm, ko)
  real(4) :: pt_zm_modify(jm, ko)

  real(4) :: p_sfc_temp1(im, jm)
  real(4) :: p_sfc_temp2(im, jm)
  real(4) :: p_sfc_temp2_zm(jm)
  real(4) :: p_sfc_modify(im, jm)
  real(4) :: p_sfc_modify_max(jm)
  
  real(4) :: pt_pds(jm)
  real(4) :: alt_min(jm)
  real(4) :: grad
  real(4) :: tmp(jm)
  
  real(4) :: d
  integer :: i, j, k, l
  integer :: ju, jl, ku, kl


  ! p_sfc_temp1 : 2-dimensional surface pressure normalized by p+s
  !               estimated from altitude
  tmp = sum( exp(-alt/h0), dim=1 ) / im
  do i=1, im
     do j=1, jm
        if( tmp(j) /= 0 ) then
           p_sfc_temp1(i,j) = p_pds(j) * exp(-alt(i,j)/h0) / tmp(j)
        else
           p_sfc_temp1(i,j) = 0
        end if
     end do
  end do

  
  ! p_sfc_modify : 2-dimensional surface pressure normalized by p+s
  !                estimated from phi_dagger
  do i=1, im
     do j=1, jm

        do k=2, ko

           if( p_sfc_temp1(i,j) <= pout(k) ) then
              p_sfc_temp2(i,j) = pout(k-1) &
                   &  + ( grav * alt(i,j) - phi_dagger(j,k-1) ) &
                   &    * ( pout(k) - pout(k-1) ) &
                   &    / ( phi_dagger(j,k) - phi_dagger(j,k-1) )
              exit

           else if( k == ko ) then
              p_sfc_temp2(i,j) = pout(ko) &
                   &  + ( grav * alt(i,j) - phi_dagger(j,ko) ) &
                   &    * ( pout(ko) - pout(ko-1) ) &
                   &    / ( phi_dagger(j,ko) - phi_dagger(j,ko-1) )
              exit

           end if

        end do

     end do
  end do
  
  p_sfc_temp2_zm = sum( p_sfc_temp2, dim=1 ) / real(im)
  
  do i=1, im
     do j=1, jm
        p_sfc_modify(i,j) = p_pds(j) * p_sfc_temp2(i,j) / p_sfc_temp2_zm(j)
     end do
  end do
  

  ! pout_modify : actual p+
  !               pout_modify != pout only near the surface
  do k=1, ko
     do j=1, jm

        pout_modify(j,k) = 0
        do i=1, im
           pout_modify(j,k) = pout_modify(j,k) + min( pout(k), p_sfc_modify(i,j) )
        end do
        pout_modify(j,k) = pout_modify(j,k) / im

     end do
  end do
  

  ! interpolate (pout levels -> pout_modify levels)

  !naiso theta & v_zm  
  do j=1, jm

     ! k=1
     v_zm_modify(j,1) = v_zm(j,1)     
     pt_zm_modify(j,1) = pt_zm(j,1)       

     do k=2, ko
        do l=ko, 2, -1

           if( pout_modify(j,k) == pout(l) ) then
              v_zm_modify(j,k) = v_zm(j,l)
              pt_zm_modify(j,k) = pt_zm(j,l)

           else if( pout_modify(j,k) < pout(l) &
                &   .and. pout_modify(j,k) > pout(l-1) ) then
              d = ( pout_modify(j,k) - pout(l-1) ) &
                   & / ( pout(l) - pout(l-1) )
              v_zm_modify(j,k) = v_zm(j,l-1) * (1.0-d) + v_zm(j,l) * d
              pt_zm_modify(j,k) = pt_zm(j,l-1) * (1.0-d) + pt_zm(j,l) * d

           end if

        end do
     end do

  end do
  
  alt_min(:) = minval( alt, dim=1 )
  p_sfc_modify_max(:) = maxval( p_sfc_modify, dim=1 )
  pt_pds(:) = minval( pt_sfc, dim=1 )

  !end naiso
  

  ! modify phi_dagger

  !pt_new -> pt_zm  &  pa_sea -->  pasmax   
  do k=ko, 1, -1
     do j=1, jm

        if( pout(k) >= p_sfc_modify_max(j) ) then
           phi_dagger_modify(j,k) = alt_min(j) * grav

        else if( k == ko ) then
           phi_dagger_modify(j,k) = alt_min(j) * grav &
                &   + gasr &
                &   * ( pt_pds(j) * (p_sfc_modify_max(j)/1000.0)**rkappa &
                &     + pt_zm_modify(j,k) * (pout(k)/1000.0)**rkappa ) / 2.0 &
                &   * ( log(p_sfc_modify_max(j)*100) - log(pout(k)*100) )

        else if( pout(k+1) >= p_sfc_modify_max(j) ) then
           phi_dagger_modify(j,k) = alt_min(j) * grav &
                &   + gasr &
                &   * ( pt_pds(j) * (p_sfc_modify_max(j)/1000.0)**rkappa &
                &     + pt_zm_modify(j,k) * (pout(k)/1000.0)**rkappa ) / 2.0 &
                &   * ( log(p_sfc_modify_max(j)*100) - log(pout(k)*100) )

        else 
           phi_dagger_modify(j,k) = phi_dagger_modify(j,k+1) &
                &   + gasr &
                &   * ( pt_zm_modify(j,k+1) * (pout(k+1)/1000.0)**rkappa &
                &     + pt_zm_modify(j,k) * (pout(k)/1000.0)**rkappa ) / 2.0 &
                &   * ( log(pout(k+1)*100) - log(pout(k)*100) )
        end if

     end do
  end do
  
  
  do j=1, jm
     grad = 0.0
     do k=3, ko

        if( pout(k) >= p_sfc_modify_max(j) .and. pout(k-1) < p_sfc_modify_max(j) ) then
           grad = ( phi_dagger_modify(j,k-1) - phi_dagger_modify(j,k-2) ) &
                & / ( log(pout(k-1)*100) - log(pout(k-2)*100) )
        end if

        if( pout(k) > p_sfc_modify_max(j) ) then
           phi_dagger_modify(j,k) = phi_dagger_modify(j,k-1) &
                &       - grad * ( log(pout(k-1)*100) - log(pout(k)*100) )
        end if

     end do
  end do


  ! modify C(Az,Kz)
  ! Probably, d(pout_modify)/d(pout) is multipled 
  ! for the accurate vertical integration.
  ! c_az_kz_modify may be on the pout_modify levels, not on the pout levels
  
  !v_new->v_zm
  do k=1, ko
     kl = k - 1
     ku = k + 1
     if( k == 1 )  kl = k
     if( k == ko ) ku = k
     
     do j=1, jm
        jl = j - 1
        ju = j + 1
        if( j == 1 )  jl = j
        if( j == jm ) ju = j
        
        c_az_kz_modify(j,k) &
             &  = -v_zm_modify(j,k) / radius  &
             &  * ( phi_dagger_modify(ju,k) - phi_dagger_modify(jl,k) ) &
             &  / ( alat(ju) - alat(jl) ) &
             &  * ( pout_modify(j,kl) - pout_modify(j,ku) ) &
             &  / ( pout(kl) - pout(ku) )
        
     end do
  end do

  return
end subroutine mount_modify
