module biseki_y
  use com_var, only : jm, ko, pout, alat
  implicit none
  private
  public :: biseki_y_ini,  biseki_y_biseki_y, &
       &    biseki_y_sekibun_y, biseki_y_bibun_y

  real(4),pointer :: pd_pdd(:,:)  ! p+ at the p++ levels [hPa]
  real(4),pointer :: p_pds(:)     ! p+s (p+ at the surface) [hPa]
  real(4),pointer :: p_pdds(:)    ! p++s (p++ at the surface) [hPa]
  integer,pointer :: nlev_y(:,:)  ! level number for interpolation
  real(4),pointer :: dlev_y(:,:)  ! weight for interpolation

contains

  !
  ! initialize
  ! 
  !   -create symbolic link to the actual data
  !
  subroutine biseki_y_ini( pd_pdd_in, p_pds_in, p_pdds_in, &
       &                   nlev_y_in, dlev_y_in )
    real(4),intent(in),target :: pd_pdd_in(jm, ko)
    real(4),intent(in),target :: p_pds_in(jm)
    real(4),intent(in),target :: p_pdds_in(1)
    integer,intent(in),target :: nlev_y_in(jm, ko)
    real(4),intent(in),target :: dlev_y_in(jm, ko)

    pd_pdd  => pd_pdd_in
    p_pds  => p_pds_in
    p_pdds => p_pdds_in
    nlev_y => nlev_y_in
    dlev_y => dlev_y_in

  end subroutine biseki_y_ini


  !
  ! Function
  !   call biseki_y_sekibun_y() and biseki_y_bibun_y()
  !
  ! Npte
  !   -Normally, it is not necessary to call biseki_y_sekibun_y()
  !    and biseki_y_bibun_y(), independently.
  ! 
  subroutine biseki_y_biseki_y( x_zm, x_ym )
    real(4),intent(in)  :: x_zm(jm, ko)
    real(4),intent(out) :: x_ym(ko)
    real(4) :: x_pt(jm, ko)
    real(4) :: xint_ym(ko)
    
    call biseki_y_sekibun_y( x_zm, x_pt, xint_ym )
    call biseki_y_bibun_y( xint_ym, x_ym )
  end subroutine biseki_y_biseki_y


  !
  ! Function
  !   integrate vertically for coordinate transformation from p+ to p++
  !
  ! Arguements (in)
  !   x_zm    : 2-dimensional p+ coordinate date
  !
  ! Arguements (out)
  !   x_pdd   : 2-dimensional p++ coordinate data
  !   xint_ym : vertically integrated data for the next step (biseki_y_bibun_y)
  !
  ! Note
  !   -Normally, biseki_y_bibun_y() is called just after this subroutine.
  !
  subroutine biseki_y_sekibun_y( x_zm, x_pdd, xint_ym )
    real(4),intent(in)  :: x_zm(jm, ko)
    real(4),intent(out) :: x_pdd(jm, ko)
    real(4),intent(out) :: xint_ym(ko)
    
    real(4) :: xint(jm, ko)
    integer :: j, k, n
    real(4) :: d
    
    do j=1, jm

       ! 2D p+ coordinate -> 2D p++ coordinate
       do k=1, ko
          n = nlev_y(j,k)
          d = dlev_y(j,k)
          x_pdd(j,k) = x_zm(j,n) * (1-d) + x_zm(j,n+1) * d             
!          if( d < 0 ) then              
!             x_pdd(j,k) = x_zm(j,k)
!          else 
!             x_pdd(j,k) = x_zm(j,n) * (1-d) + x_zm(j,n+1) * d
!          endif
       end do
          

       ! integrate with respect to p+
       xint(j,1) = x_pdd(j,1) * pd_pdd(j,1)
          
       do k=2, ko  ! upper -> lower
             
          if( pd_pdd(j,k) > p_pds(j) ) then
             xint(j,k) = 0.5 * ( x_pdd(j,k) + x_pdd(j,k-1) ) * &
                  &            ( p_pds(j) - pd_pdd(j,k-1) )  &
                  &    + xint(j,k-1)
          else
             xint(j,k) = 0.5 * ( x_pdd(j,k) + x_pdd(j,k-1) ) * &
                  &            ( pd_pdd(j,k) - pd_pdd(j,k-1) ) &
                  &    + xint(j,k-1)
          end if
       end do
       
    end do

    ! meridional mean
    call integral_meridional( 1, jm, ko, alat, xint, xint_ym )
      
    return
  end subroutine biseki_y_sekibun_y


  !
  ! Function
  !   differentiate with respect to p++ in order to transform coordinate
  !
  ! Arguement (in)
  !   xint_ym : result of biseki_y_sekibun_y()
  !
  ! Arguement (out)
  !   x_ym    : 1-dimensional p++ coordinate data
  !
  ! Note
  !   -Normally, biseki_y_sekibun_y() is run before this subroutine.
  !
  subroutine biseki_y_bibun_y(xint_ym, x_ym)
    implicit none
    real(4),intent(in)  :: xint_ym(ko)
    real(4),intent(out) :: x_ym(ko)
    
    call derivative_p( 1, 1, ko, pout, p_pdds, xint_ym, &
         &             x_ym )
  
    return
  end subroutine biseki_y_bibun_y

end module biseki_y
