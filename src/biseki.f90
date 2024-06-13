module biseki
  use com_var, only : im, jm, km, ko, pout
  implicit none
  private
  public :: biseki_ini, biseki_biseki, biseki_sekibun, biseki_bibun

  real(4),pointer :: p_pd(:,:,:)  ! pressure at the p+ levels [hPa]
  real(4),pointer :: p_sfc(:,:)   ! pressure at the surface [hPa]
  real(4),pointer :: p_pds(:)     ! p+s (p-dagger at the surface) [hPa]
  integer,pointer :: nlev(:,:,:)  ! level number for interpolation
  real(4),pointer :: dlev(:,:,:)  ! weight for interpolation

contains

  !
  ! initialize
  !
  !   -create symbolic link to the actual data
  ! 
  subroutine biseki_ini( p_pd_in, p_sfc_in, p_pds_in, nlev_in, dlev_in )
    real(4),intent(in),target :: p_pd_in(im, jm, km)
    real(4),intent(in),target :: p_sfc_in(im, jm)
    real(4),intent(in),target :: p_pds_in(jm)
    integer,intent(in),target :: nlev_in(im, jm, km)
    real(4),intent(in),target :: dlev_in(im, jm, km)

    p_pd  => p_pd_in
    p_sfc => p_sfc_in
    p_pds => p_pds_in
    nlev  => nlev_in
    dlev  => dlev_in

  end subroutine biseki_ini


  !
  ! Function
  !   call biseki_sekibun() and biseki_bibun()
  !
  ! Npte
  !   -Normally, it is not necessary to call biseki_sekibun()
  !    and biseki_bibun(), independently.
  ! 
  subroutine biseki_biseki( x, x_zm )
    real(4),intent(in)  :: x(im, jm, km)
    real(4),intent(out) :: x_zm(jm, ko)
    real(4) :: x_pt(im, jm, ko)
    real(4) :: xint_zm(jm, ko)

    call biseki_sekibun( x, x_pt, xint_zm )
    call biseki_bibun( xint_zm, x_zm )
  end subroutine biseki_biseki



  !
  ! Function
  !   integrate vertically for coordinate transformation from p to p+
  !
  ! Arguements (in)
  !   x       : 3-dimensional pressure coordinate date
  !
  ! Arguements (out)
  !   x_pd    : 3-dimensional p+ coordinate data
  !   xint_zm : vertically integrated data for the next step (biseki_bibun)
  !
  ! Note
  !   -Normally, biseki_bibun() is called just after this subroutine.
  !
  subroutine biseki_sekibun( x, x_pd, xint_zm )
    real(4),intent(in)  :: x(im, jm, km)
    real(4),intent(out) :: x_pd(im, jm, ko)
    real(4),intent(out) :: xint_zm(jm, ko)
    
    real(4) :: xint(im, jm, ko)
    integer :: i, j, k, n
    real(4) :: d
    
    do j=1, jm
       do i=1, im
          
          ! 3D pressure coordinate -> 3D p+ coordinate
          do k=1, ko
             n = nlev(i,j,k)
             d = dlev(i,j,k)
             x_pd(i,j,k) = x(i,j,n) * (1-d) + x(i,j,n+1) * d
!             if( d < 0 ) then              
!                x_pd(i,j,k) = x(i,j,k)
!             else 
!                x_pd(i,j,k) = x(i,j,l) * (1-d) + x(i,j,l+1) * d
!             endif
          end do

          
          ! integrate with respect to p
          xint(i,j,1) = x_pd(i,j,1) * p_pd(i,j,1)
          
          do k=2, ko  ! upper -> lower

             if( p_pd(i,j,k) > p_sfc(i,j) ) then  ! near (or under) the ground
                xint(i,j,k) = 0.5 * ( x_pd(i,j,k) + x_pd(i,j,k-1) ) * &
                     &              ( p_sfc(i,j)  - p_pd(i,j,k-1) ) &
                     &      + xint(i,j,k-1)
             else
                xint(i,j,k) = 0.5 * ( x_pd(i,j,k) + x_pd(i,j,k-1) ) * &
                     &              ( p_pd(i,j,k) - p_pd(i,j,k-1) ) &
                     &      + xint(i,j,k-1)
             end if

          end do
          
       end do
    end do
    
    ! zonal mean
    xint_zm = sum( xint, dim=1 ) / im
    
    return
  end subroutine biseki_sekibun



  !
  ! Function
  !   differentiate with respect to p+ in order to transform coordinate
  !
  ! Arguement (in)
  !   xint_zm : result of biseki_sekibun()
  !
  ! Arguement (out)
  !   x_zm    : 2-dimensional p+ coordinate data
  !
  ! Note
  !   -Normally, biseki_sekibun() is run before this subroutine.
  !
  subroutine biseki_bibun( xint_zm, x_zm )
    implicit none
    real(4),intent(in)  :: xint_zm(jm, ko)
    real(4),intent(out) :: x_zm(jm, ko)
    
    integer :: k, j
    real(4) :: x_hlf(jm, ko+1), p_hlf(jm, ko+1), d
    integer :: ncalc

    ! switch for derivation
    ncalc = 1   ! nacalc=1 ---> standard
    !ncalc = 2   ! for comparision with mochi-version
    
    if( ncalc == 1 ) then
       call derivative_p( 1, jm, ko, pout, p_pds, xint_zm, &
            &             x_zm )

!       do j=1, jm
!          do k=1, ko
!             if( pout(k) > p_pds(j) ) then
!                x_zm(j,k) = 0
!             end if
!
!          end do
!       end do
       
!       do k=1, ko
!          write(*,*) k, x_zm(60,k)
!       end do

       
    else ! ncalc ne 1 
       
       !========intpl ---->> bibun=============================
       !   if vertical grid is uniform ,calc.2 will be better than calc.1
       !-----------interpolate xint_zm(k) on half level-------------
       
       !set new half level
       do j=1, jm
          p_hlf(j,1) = max( 0.5*pout(1), -0.5*pout(2)+1.5*pout(1) )
          
          do k=2, ko+1
             
             p_hlf(j,k) = 2 * pout(k-1) - p_hlf(j,k-1)
             
             if( p_hlf(j,k) >= pout(k) .and. k /= ko+1 )then
                
                write(6,*) &
                     &  "warning: vertical grid is not appropriate. ", &
                     &  "p must be pout(k+1)>p_hlf(k)>pout(k))", &
                     &  k, p_hlf(j,k), pout(k)
                
             end if
             
          end do
          
          d = ( p_hlf(j,1) - pout(1) ) / ( pout(2) - pout(1) )
          x_hlf(j,1) = xint_zm(j,1) * (1.0-d) + xint_zm(j,2) * d !k=1 extpl d>0
          
          do k=2, ko
             
             d = ( p_hlf(j,k) - pout(k-1) ) / ( pout(k) - pout(k-1) ) !intpl
             x_hlf(j,k) = xint_zm(j,k-1) * (1.0-d) + xint_zm(j,k) * d
             
          end do
          
          d = ( p_hlf(j,ko+1) - pout(ko-1) ) / ( pout(ko) - pout(ko-1) )
          !k=ko+1 extpl d>1
          x_hlf(j,ko+1) = xint_zm(j,ko-1) * (1.0-d) + xint_zm(j,ko) * d 
          
          !----------bibun -------------------------
          
          do k=1, ko
             
             x_zm(j,k) = ( x_hlf(j,k+1) - x_hlf(j,k) ) / &
                  &      ( p_hlf(j,k+1) - p_hlf(j,k) )
             
          end do
       end do
       
    end if
  


    return
  
  end subroutine biseki_bibun



end module biseki
