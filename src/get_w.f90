!
! Function
!   estimate vertical velocity using continuity equation
!
! Arguments (in)
!   u        : x-wind
!   v        : y-wind
!   ps       : surface pressure [hPa]
!
! Arguments (out)
!   omega    : vertical velocity
!
! Note
!   - This subroutine is used, if input does not have vertical velocity
!
subroutine get_omega(u, v, ps, &
     &                omega)
  use com_var, only: im, jm, km, pin
  implicit none
  real(4),intent(in) :: u(im, jm, km)
  real(4),intent(in) :: v(im, jm, km)
  real(4),intent(in) :: ps(im, jm)
  !
  real(4),intent(out) :: omega(im, jm, km)
  !
  real(4) :: div(im, jm, km)
  integer :: i, j, k
  real(4) :: p_d, p_u, var_d, var_u, dp

  ! calculate horizontal divergence
  call get_div(u, v, ps, div)

  ! calculate omega from vertical integration of horizonal divergence
  omega = 0.0

  do k = 2, km
     do j = 1, jm
        do i = 1, im

           !---- lower
           if( pin(k) < ps(i,j) ) then
              p_d = pin(k)
              var_d = div(i,j,k)
           else
              p_d = ps(i,j)
              var_d = div(i,j,k)
           end if

           !---- upper
           if( pin(k-1) < ps(i,j) ) then
              p_u = pin(k-1)
              var_u = div(i,j,k-1)
           else
              p_u = ps(i,j)
              var_u = div(i,j,k-1)
           end if

           dp = max(p_d-p_u, 0.0)

           ! vertical integration
           omega(i,j,k) = 0.5 * (var_d+var_u) * dp * 100.0 &
                &       + omega(i,j,k-1)

        end do
     end do
  end do

  omega(:,:,:) = -omega(:,:,:)

  return
end subroutine get_omega


!===============================================================
subroutine get_div(u,v, ps, &
     &                div)
  use parameter, only: radius, pai, radian
  use com_var, only: im, jm, km, alat, pin
  implicit none
  real(4),intent(in) :: u(im, jm, km)
  real(4),intent(in) :: v(im, jm, km)
  real(4),intent(in) :: ps(im, jm)
  !
  real(4),intent(out) :: div(im, jm, km)
  !
  integer :: i, j, k
  real(4) :: dlon, dlat
  real(4) :: fn, fs, fe, fw
  real(4) :: int, ds
  real(4) :: dsn, dss
  real(4) :: stg_costbl(jm-1), stg_sintbl(jm-1)
  real(4) :: rjm, phi(jm-1)
  real(4),allocatable :: uft(:,:,:), vft(:,:,:)
  !
  dlon = 360.0 / real(im)  * radian     ! [radian]

  !--- get staggard sin&cos
  do j = 1, jm-1
     phi(j) = 0.5 * (alat(j) + alat(j+1))
  end do
  stg_sintbl(:) = sin( phi(:) )
  stg_costbl(:) = cos( phi(:) )

  !--- check lower boundary & substitute 0 below surface
  allocate( uft(0:im+1, jm, km), vft(im, jm, km) )
  !
  uft(1:im,:,:) = u(:,:,:)
  vft(:,:,:) = v(:,:,:)

  do k = 1, km
     do j = 1, jm
        do i = 1, im

           if( pin(k) > ps(i,j) ) then
              uft(i,j,k) = 0.0
              vft(i,j,k) = 0.0
           end if

        end do
     end do
  end do
  !
  uft(0, 1:jm, 1:km)    = uft(im, 1:jm, 1:km)
  uft(im+1, 1:jm, 1:km) = uft(1, 1:jm, 1:km)

  !---- calculate divergence
  div = 0.0

  do k = 1, km
     do j = 2, jm-1
        do i = 1, im

           !---- calculate wind at boundaries
           ! north
           fn = 0.5 * (vft(i,j,k) + vft(i,j-1,k))
           ! south
           fs = 0.5 * (vft(i,j+1,k) + vft(i,j,k))
           ! east
           fe = 0.5 * (uft(i,j,k) + uft(i+1,j,k))
           ! west
           fw = 0.5 * (uft(i-1,j,k) + uft(i,j,k))

           dlat =   0.5 * ( alat(j) + alat(j-1) ) &
                & - 0.5 * ( alat(j) + alat(j+1) )              ! [radian]


           !--- calculate integration
           int =   radius * dlat * (fe-fw) &
                & -radius * stg_costbl(j)   * dlon * fs &
                & +radius * stg_costbl(j-1) * dlon * fn

           !---- ds: area element
           ds = radius**2 * (stg_sintbl(j-1) - stg_sintbl(j)) * dlon

           !--- divergence
           if( pin(k) < ps(i,j) ) then
              div(i,j,k) = int / ds
           else
              div(i,j,k) = 0.0
           end if
        end do
     end do


     !--- pole
     fn = 0.0
     fs = 0.0
     do i = 1, im
        fn = -0.5 * (vft(i, 1,k) + vft(i,   2,k)) * radius &
             &  * stg_costbl(1)    * dlon + fn

        fs =  0.5 * (vft(i, jm,k) + vft(i,jm-1,k)) * radius &
             &  * stg_costbl(jm-1) * dlon + fs
     end do
     dsn = 2 * pai * radius**2 * ( 1-stg_sintbl(1))
     dss = 2 * pai * radius**2 * (-1-stg_sintbl(jm-1))
     !
     if( pin(k) < ps(1,1) ) then
        div(1,1,k) = fn / dsn
     else
        div(1,1,k) = 0.0
     end if

     if( pin(k) < ps(1,jm) ) then
        div(1,jm,k) = fs / dss
     else
        div(1,jm,k) = 0.0
     end if
     !
     do i = 2, im
        div(i, 1,k) = div(1, 1,k)
        div(i,jm,k) = div(1,jm,k)
     end do
  end do


  deallocate(uft, vft)

  return
end subroutine get_div
