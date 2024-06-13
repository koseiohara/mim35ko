!
! 1/2 int[ var cos(phi) ] d(phi)
!
subroutine integral_meridional( im, jm, ko, alat, var, &
     &                          var_ym )
  use parameter, only: pai
  implicit none
  integer,intent(in)  :: im, jm, ko
  real(4),intent(in)  :: alat(jm)
  real(4),intent(in)  :: var(im, jm, ko)
  real(4),intent(out) :: var_ym(im, ko)

  integer :: i, j, k
  real(4) :: pai2
  pai2 = 0.5 * pai
  
  do i=1, im
     do k=1, ko
        
        ! trapezoidal integration
        var_ym(i,k) = var(i,1,k) * cos(alat(1)) * &
             &        ( pai2 - alat(1) ) / 4

        do j=2, jm
           var_ym(i,k) = var_ym(i,k) + &
                &        ( var(i,j-1,k) * cos(alat(j-1)) + &
                &          var(i,j,k)   * cos(alat(j))     ) &
                &      * ( alat(j-1) - alat(j) ) / 4
        end do

        var_ym(i,k) = var_ym(i,k) + &
             &        var(i,jm,k) * cos(alat(jm)) * &
             &        ( alat(jm) + pai2 ) / 4

     end do
  end do
     
end subroutine integral_meridional



!
! xint = 1/g int[ x_zm ] dp
!
!   p: p+ or p++
!
!  if p=p+  -> ps_zm = p_pds
!  if p=p++ -> ps_zm = p_pdds
!
subroutine integral_p( jm, ko, pout, ps_zm, x, &
     &                 xint )
  use parameter, only : grav
  implicit none
  integer,intent(in)  :: jm, ko
  real(4),intent(in)  :: x(jm, ko)
  real(4),intent(in)  :: ps_zm(jm)
  real(4),intent(in)  :: pout(ko)
  real(4),intent(out) :: xint(jm)

  integer :: j, k
  real(4) :: x_tmp(jm)
  
  ! trapeziodal integration
  do j=1, jm
     
     xint(j) = x(j,1) * pout(1)
     
     do k=2, ko
        
        if( pout(k) > ps_zm(j) .and. pout(k-1) < ps_zm(j) ) then
           
           x_tmp(j) = ( x(j,k) + x(j,k-1) ) &
                &   * ( ps_zm(j) - pout(k-1) ) &
                &   / ( pout(k)-pout(k-1) ) / 2
           
           xint(j) = 0.5 * ( x_tmp(j) + x(j,k-1) ) &
                &  * ( ps_zm(j) - pout(k-1) ) &
                &  + xint(j)
           
        else if( pout(k) <= ps_zm(j) ) then
           
           xint(j) = 0.5 * ( x(j,k) + x(j,k-1) ) &
                &  * ( pout(k) - pout(k-1) ) &
                &  + xint(j)
           
        else
           xint(j) = xint(j)
        end if
        
     end do
     
     xint(j) = 100.0 / grav * xint(j)
     
  end do
  
  return
end subroutine integral_p



!
! x_vint = int[ x_zm ] d(pt)
!   [ pt_min : pt_top ]
!
subroutine integral_pt( jm, ko, pout, p_pds, pt_zm, pt_pds, x_zm, &
     &                  x_vint )
  implicit none
  integer,intent(in)  :: jm, ko
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(in)  :: pt_pds(jm)
  real(4),intent(in)  :: x_zm(jm, ko)
  real(4),intent(out) :: x_vint(jm)
  integer :: j, k
  
  ! trapeziodal integration
  do j=1, jm

     !x_vint(j) = pt_zm(j,1) * x_zm(j,1)
     x_vint(j) = 0  ! above top -> 0

     do k=2, ko

        if( pout(k) <= p_pds(j) ) then

           x_vint(j) = x_vint(j) + ( pt_zm(j,k-1) - pt_zm(j,k) ) &
                &                * (  x_zm(j,k-1) +  x_zm(j,k) ) / 2

           if( k == ko .and. pt_zm(j,ko) >= pt_pds(j) ) then
              x_vint(j) = x_vint(j) &
                   &    + ( pt_zm(j,ko) - pt_pds(j) ) * x_zm(j,ko)
           end if

        else if( ( pout(k) > p_pds(j) ) &
             &   .and. ( pout(k-1) <= p_pds(j) ) &
             &   .and. ( pt_zm(j,k-1) >= pt_pds(j) ) ) then

           x_vint(j) = x_vint(j) + ( pt_zm(j,k-1) - pt_pds(j) ) &
                &                * (  x_zm(j,k-1) + x_zm(j,k) ) / 2

        end if

     end do
  end do
  
end subroutine integral_pt




!
! x_vint = int[ x_pdd ] d(pt)
!   [ pt_ymin : pt_top ]
!
subroutine integral_pt_ym( jm, ko, pout, p_pdds, pt_ym, pt_pdds, x_pdd, &
     &                     x_vint )
  implicit none
  integer,intent(in)  :: jm, ko
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: p_pdds(1)
  real(4),intent(in)  :: pt_ym(ko)
  real(4),intent(in)  :: pt_pdds(1)
  real(4),intent(in)  :: x_pdd(jm, ko)
  real(4),intent(out) :: x_vint(jm)
  integer :: j, k
  
  ! trapeziodal integration
  do j=1, jm

     x_vint(j) = 0  ! above top -> 0

     do k=2, ko

        if( pout(k) <= p_pdds(1) ) then

           x_vint(j) = x_vint(j) + (   pt_ym(k-1) -    pt_ym(k) ) &
                &                * ( x_pdd(j,k-1) +  x_pdd(j,k) ) / 2

           if( k == ko .and. pt_ym(ko) >= pt_pdds(1) ) then
              x_vint(j) = x_vint(j) &
                   &    + ( pt_ym(ko) - pt_pdds(1) ) * x_pdd(j,ko)
           end if

        else if( ( pout(k) > p_pdds(1) ) &
             &   .and. ( pout(k-1) <= p_pdds(1) ) &
             &   .and. ( pt_ym(k-1) >= pt_pdds(1) ) ) then

           x_vint(j) = x_vint(j) + (    pt_ym(k-1) - pt_pdds(1) ) &
                &                * (  x_pdd(j,k-1) + x_pdd(j,k) ) / 2

        end if

     end do
  end do
  
end subroutine integral_pt_ym
