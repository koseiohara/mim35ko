!
! get_z_pt() - get geopotential height on p+ levels
!
!   t      : 3D temperature on p levels [K]
!   z      : 3D geopotential height on p levels [m]
!   p_pt   : 3D pressure on p+ levels [hPa]
!   p_pds  : p+ at the surface [hPa]
!   alt    : altitude [m]
!   z_pt   : 3D geopotential height on p+ levels [m]
!   z_zm   : zonal mean geopotential height z(y,p+) [m]
!            for EP Flux calculation
!
subroutine get_z_pt( alt, t, z, p_pt, p_pds, &
     &               z_pt, z_zm )
  use parameter, only : rkappa, grav, cp, h0, z_min, z_max
  use com_var, only : im, jm, km, ko, pin, pout
  implicit none
  real(4),intent(in)  :: alt(im, jm)
  real(4),intent(in)  :: t(im, jm, km)
  real(4),intent(in)  :: z(im, jm, km)
  real(4),intent(in)  :: p_pt(im, jm, ko)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(out) :: z_pt(im, jm, ko)
  real(4),intent(out) :: z_zm(jm, ko)

  real(4) :: ganma(im, jm, km)
  real(4) :: hwork1, hwork2
  integer :: l, i, k, j, kl, ku
  integer :: sw  ! interpolation method
  real(4) :: grad
  real(4) :: const

  const = rkappa * cp / grav

  do j=1, jm
     do i=1, im

        do k=1, km
           if( k == 1 ) then
              ganma(i,j,k) = ( t(i,j,k+1) - t(i,j,k) ) / &
                   &         ( z(i,j,k+1) - z(i,j,k) )
           else
              ganma(i,j,k) = ( t(i,j,k) - t(i,j,k-1) ) / &
                   &         ( z(i,j,k) - z(i,j,k-1) )
           end if
        end do


        do k=1, ko

           ! get kl & ku
           do l=2, km
              if( pin(l-1) < p_pt(i,j,k) .and. p_pt(i,j,k) <= pin(l) ) then
                 kl = l - 1
                 ku = l
                 sw = 1
                 exit
              end if
           end do

           if( p_pt(i,j,k) <= pin(1) ) then
              kl = 1
              sw = 2
           else if( p_pt(i,j,k) > pin(km) ) then
              kl = km
              sw = 2
           end if
           
           ! interpolation
           if( sw == 1 ) then
              hwork1 = -const &
                   &   * ( log( p_pt(i,j,k) / pin(ku) ) ) * t(i,j,ku) &
                   & + 0.5 * const * h0 * ganma(i,j,ku) &
                   &   * ( log( p_pt(i,j,k) / pin(ku) ) )**2
              
              hwork2 = -const &
                   &   * ( log( pin(kl) / pin(ku) ) ) * t(i,j,ku) &
                   & + 0.5 * const * h0 * ganma(i,j,ku) &
                   &   * ( log( pin(kl) / pin(ku) ) )**2 
              
              z_pt(i,j,k) = z(i,j,ku) &
                   &      + ( z(i,j,kl) - z(i,j,ku) ) * hwork1 / hwork2
              
           else if( sw == 2 ) then
              z_pt(i,j,k) = z(i,j,kl) &
                   &      - const &
                   &        * ( log( p_pt(i,j,k) / pin(kl) ) ) &
                   &        * t(i,j,kl) &
                   &      + 0.5 * const * h0 * ganma(i,j,kl) &
                   &        * ( log( p_pt(i,j,k) / pin(kl) ) )**2 
              
           end if

           ! modify
           if( z_pt(i,j,k) < 0.0 ) z_pt(i,j,k) = 0.0
           if( z_pt(i,j,k) <= alt(i,j) ) z_pt(i,j,k) = alt(i,j)
           if( k /= 1 ) then
              if( z_pt(i,j,k) > z_pt(i,j,k-1) ) z_pt(i,j,k) = z_pt(i,j,k-1)
           end if
           
        end do
        
     end do
  end do

  call check_range( im, jm, ko, z_pt, z_min, z_max, 'get_z_pt()', 'z_pt' )

  
  ! z_zm
  z_zm(:,:) = sum( z_pt, dim=1 ) / im

  ! modify
  do j=1, jm

     grad = 0.0
     do k=2, ko

        if( pout(k-1) < p_pds(j) .and. p_pds(j) <= pout(k) ) then
           grad = ( z_zm(j,k-1) - z_zm(j,k-2) ) &
                & / log( pout(k-1) / pout(k-2) )
        end if

        if( p_pds(j) < pout(k) ) then
           z_zm(j,k) = z_zm(j,k-1) &
                &    - grad * log( pout(k-1) / pout(k) )
        end if

     end do

  end do

  call check_range( 1, jm, ko, z_zm, z_min, z_max, 'get_z_pt()', 'z_zm' )

end subroutine get_z_pt
