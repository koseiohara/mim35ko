!
! Function
!   estimate D(pt)/Dt (i.e. diabatic heating) from total differential
!
! Arguements (in)
!   dt        : timestep of the input date
!   u         : x-wind
!   v         : y-wind
!   omega     : omega velocity
!   pt        : potential temperature at the pressure levels
!   pt_before : pt 1-step before 
!
! Arguements (out)
!   pt_dot    : D(pt)/Dt
!
! Note
!   -If omega=0 (no input), vertical advection of pt will be neglected.
!   -If pt=pt_before (e.g. first time step), local pt change will be neglected.
!
subroutine get_pt_dot_omega( dt, u, u_before, v, v_before, &
     &                       omega, omega_before, pt, pt_before, &
     &                       pt_dot )
  use parameter, only: radius, pai
  use com_var, only : im, jm, km, alat, costbl, pin
  implicit none
  real(4),intent(in)  :: dt
  real(4),intent(in)  :: u(im, jm, km)
  real(4),intent(in)  :: u_before (im, jm, km)
  real(4),intent(in)  :: v(im, jm, km)
  real(4),intent(in)  :: v_before (im, jm, km)
  real(4),intent(in)  :: omega(im, jm, km)
  real(4),intent(in)  :: omega_before (im, jm, km)
  real(4),intent(in)  :: pt(im, jm, km)
  real(4),intent(in)  :: pt_before(im, jm, km)
  real(4),intent(out) :: pt_dot(im, jm, km)

  integer :: i, j, k, nt

!Å@Splitting analysis step into ntm short time steps
  data ntm/6/


! tendency of pt
  do k=1, km                    
     do j=1, jm
        do i=1, im

           ! d(pt)/dt
           pt_dot(i,j,k) = ( pt(i,j,k) - pt_before(i,j,k) ) / dt 

        end do
     end do
  end do

! averaging advection terms computed at short time steps

  do nt=1, ntm
  do k=1, km                    
     do j=1, jm
        do i=1, im

!Å@tint  ( 0 < tint < 1 ) time interpolation parameter

           tint = (nt - 0.5) / ntm

           ! u d(pt)/dx
           if( i == 1 ) then
              pt_east = pt_before(2,j,k) *(1-tint) + pt(2,j,k) *tint 
              pt_west = pt_before(im,j,k)*(1-tint) + pt(im,j,k)*tint 

            !  pt_dot(i,j,k) = pt_dot(i,j,k)     &
            !       &        + u(i,j,k) * ( pt(i+1,j,k) - pt(im,j,k) )&
            !       &          / ( 4 * pai * radius * costbl(j) / im ) 
           else if( i == im ) then
              pt_east = pt_before(1,j,k)  *(1-tint) + pt(1,j,k)  *tint 
              pt_west = pt_before(i-1,j,k)*(1-tint) + pt(i-1,j,k)*tint

            ! pt_dot(i,j,k) = pt_dot(i,j,k) &
            !     &        + u(i,j,k) * ( pt(1,j,k) - pt(i-1,j,k) ) &
            !     &          / ( 4 * pai * radius * costbl(j) / im )
           else
              pt_east = pt_before(i+1,j,k)*(1-tint)+pt(i+1,j,k)*tint
              pt_west = pt_before(i-1,j,k)*(1-tint)+pt(i-1,j,k)*tint

            ! pt_dot(i,j,k) = pt_dot(i,j,k) &
            !    &        + u(i,j,k) * ( pt(i+1,j,k) - pt(i-1,j,k) ) &
            !    &          / ( 4 * pai * radius * costbl(j) / im )
           endif

           u_int  = u_before(i,j,k) * (1-tint) + u(i,j,k) * tint

           pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &      + u_int * ( pt_east - pt_west ) &
                   &      / ( 4 * pai * radius * costbl(j) / (im*ntm) )

           ! v d(pt)/dy
           if( j == 1 ) then
              pt_nth = pt_before(i,j+1,k)*(1-tint) + pt(i,j+1,k)*tint
              pt_sth = pt_before(i,j,k)  *(1-tint) + pt(i,j,k)  *tint
              dy = ( alat(j+1) - alat(j) ) * radius

              ! pt_dot(i,j,k) = pt_dot(i,j,k) &
              !     &        + v(i,j,k) * ( pt(i,j+1,k) - pt(i,1,k) ) &
              !     &          / ( alat(j+1) - alat(1) ) / radius 
           else if( j == jm ) then
              pt_nth = pt_before (i,j,k) *(1-tint) + pt(i,j,k)  *tint
              pt_sth = pt_before(i,j-1,k)*(1-tint) + pt(i,j-1,k)*tint
              dy = ( alat(j) - alat(j-1) ) * radius

              ! pt_dot(i,j,k) = pt_dot(i,j,k) &
              !     &        + v(i,j,k) * ( pt(i,jm,k) - pt(i,j-1,k) ) &
              !     &          / ( alat(jm) - alat(j-1) ) / radius 
           else
              pt_nth = pt_before(i,j+1,k) *(1-tint) + pt(i,j+1,k)*tint
              pt_sth = pt_before(i,j-1,k) *(1-tint) + pt(i,j-1,k)*tint
              dy = ( alat(j+1) - alat(j-1) ) * radius

              ! pt_dot(i,j,k) = pt_dot(i,j,k) &
              !    &        + v(i,j,k) * ( pt(i,j+1,k) - pt(i,j-1,k) ) &
              !    &          / ( alat(j+1) - alat(j-1) ) / radius 
           endif

           v_int  = v_before(i,j,k) * (1-tint) + v(i,j,k) * tint
           pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &        + v_int * ( pt_nth - pt_sth ) / dy / ntm

           ! omega d(pt)/dp
           if( k == 1 ) then
              continue
           else if( k == km ) then
              pt_d = pt_before(i,j,k) *(1-tint) + pt(i,j,k)*tint
              pt_u = pt_before(i,j,k-1) *(1-tint) + pt(i,j,k-1)*tint
              omega_int  = omega_before(i,j,k)*(1-tint)+omega(i,j,k)*tint

              pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &           + omega_int * (pt_d - pt_u ) &
                   ! &           + omega(i,j,k) * (pt(i,j,k) - pt(i,j,k-1) ) &
                   &             / ( ( pin(k) - pin(k-1) ) * 100 ) 
           else
                   pt_d = pt_before(i,j,k+1) *(1-tint) + pt(i,j,k+1)*tint
                   pt_u = pt_before(i,j,k-1) *(1-tint) + pt(i,j,k-1)*tint
                   omega_int  = omega_before(i,j,k)*(1-tint)+omega(i,j,k)*tint

                   pt_dot(i,j,k) = pt_dot(i,j,k) &
                   &           + omega_int * (pt_d - pt_u ) &
                   ! &          + omega(i,j,k) * ( pt(i,j,k+1) - pt(i,j,k-1) ) &
                   &            / ( ( pin(k+1) - pin(k-1) ) * 100 )
           endif

        end do
     end do
  end do
  end do
end subroutine get_pt_dot_omega



!
! Function
!   estimate D(pt)/Dt (i.e. diabatic heating) from diabatic heating input
!
! Arguements (in)
!   q_3d      : input diabatic heating data
!
! Arguements (out)
!   pt_dot    : D(pt)/Dt
!
subroutine get_pt_dot_q( q_3d, pt_dot )
  use parameter, only: cp, rkappa
  use com_var, only : im, jm, km, pin
  implicit none
  real(4),intent(in)  :: q_3d(im, jm, km)  ! [J/(kg s)]
  real(4),intent(out) :: pt_dot(im, jm, km)
  
  integer :: i, j, k

  do k=1, km
     do j=1, jm
        do i=1, im
           pt_dot(i,j,k) = q_3d(i,j,k) / ( cp * (pin(k)/1000.0)**rkappa )
        end do
     end do
  end do

end subroutine get_pt_dot_q


subroutine get_pt_dot_q_inv( pt_dot, q_3d )
  use parameter, only: cp, rkappa
  use com_var, only : im, jm, km, pin
  implicit none
  real(4),intent(in)  :: pt_dot(im, jm, km)
  real(4),intent(out) :: q_3d(im, jm, km)  ! [J/(kg s)]
  
  integer :: i, j, k

  do k=1, km
     do j=1, jm
        do i=1, im
           q_3d(i,j,k) =  pt_dot(i,j,k) * ( cp * (pin(k)/1000.0)**rkappa )
        end do
     end do
  end do

end subroutine get_pt_dot_q_inv


