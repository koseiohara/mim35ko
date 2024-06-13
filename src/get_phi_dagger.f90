!
! Function
!   get geopotential by integrating temperature
!
! Arguements (in)
!   alt      : altitude
!   p_pds    : p+s
!   pt_pds   : pt at the surface in p+ system
!   t_dagger : T+ (temperature in p+ system)
!
! Arguements (out)
!   phi_dagger : Phi+ (geopotential in p+ system)
!
! Note
!   -phi_dagger and z_zm is different:
!      phi_dagger is obtained from temperature.
!      z_zm is obtained directly from input geopotential data.
!
subroutine get_phi_dagger( alt, p_pds, pt_pds, t_dagger, &
     &                     phi_dagger )
  use parameter, only : grav, rkappa, gasr
  use com_var, only : im, jm, ko, pout
  implicit none
  real(4),intent(in)  :: alt(im, jm)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: pt_pds(jm)
  real(4),intent(in)  :: t_dagger(jm, ko)
  real(4),intent(out) :: phi_dagger(jm, ko)

  real(4) :: alt_zm(jm)
  real(4) :: t_dagger_pds(jm)
  real(4) :: grad
  integer :: j, k

  alt_zm = sum( alt, dim=1 ) / real(im)
  t_dagger_pds(:) = pt_pds(:) * ( p_pds(:) / 1000.0 )**rkappa

  ! integrate : R int[ T_dagger ] dlog(p+)
  do j=1, jm
     do k=ko, 1, -1
        
        if( pout(k) >= p_pds(j) )then
           phi_dagger(j,k) = alt_zm(j) * grav
        else if( k == ko ) then
           phi_dagger(j,k) = alt_zm(j) * grav &
                &     + gasr * ( t_dagger_pds(j) + t_dagger(j,k) ) / 2.0 &
                &            * ( log( p_pds(j)*100 ) - log( pout(k)*100 ) )
        else if( pout(k+1) > p_pds(j) ) then
           phi_dagger(j,k) = alt_zm(j) * grav &
                &     + gasr * ( t_dagger_pds(j) + t_dagger(j,k) ) / 2.0 &
                &            * ( log( p_pds(j)*100 ) - log( pout(k)*100 ) )
        else 
           phi_dagger(j,k) = phi_dagger(j,k+1) &
                &     + gasr * ( t_dagger(j,k+1) + t_dagger(j,k) ) / 2.0 &
                &            * ( log( pout(k+1)*100 ) - log( pout(k)*100 ) )
        endif
        
     end do
  end do


  ! modify
  do j=1, jm
     grad = 0.0
     do k=3, ko

        if( pout(k) >= p_pds(j) .and. pout(k-1) < p_pds(j) ) then
           grad = ( phi_dagger(j,k-1) - phi_dagger(j,k-2) ) &
                & / ( log( pout(k-1)*100 ) - log( pout(k-2)*100 ) )
        end if

        if( pout(k) > p_pds(j) ) then
           phi_dagger(j,k) = phi_dagger(j,k-1) &
                &          - grad * ( log(pout(k-1)*100) - log(pout(k)*100) )
        end if

     end do
  end do

end subroutine get_phi_dagger
