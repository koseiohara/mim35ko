!
! meridional component of EP flux
!
subroutine epflux_y( u_d_v_d_zm, epy )
  use parameter, only: radius
  use com_var, only : jm, ko, costbl, pout, rho
  implicit none
  real(4),intent(in)  :: u_d_v_d_zm(jm, ko)
  real(4),intent(out) :: epy(jm, ko)
  
  integer :: j, k
  
  do j=1, jm
     do k=1, ko
        epy(j,k) = -radius * costbl(j) * rho(k) &
             &   * u_d_v_d_zm(j,k)
     end do
  end do
  
  return
end subroutine epflux_y



!
! one of the vertical component of EP flux
!
subroutine epflux_z_uw( pt_zm, u_v_x_zm, u_pt_dot_x_zm, &
     &                  epz_uv, epz_ut, epz_uw )
  use parameter, only : grav, radius
  use com_var, only : jm, ko, pout, alat, costbl, rho
  implicit none
  real(4),intent(in)  :: pt_zm(jm, ko) 
  real(4),intent(in)  :: u_v_x_zm(jm, ko)
  real(4),intent(in)  :: u_pt_dot_x_zm(jm, ko)
  real(4),intent(out) :: epz_uv(jm, ko)
  real(4),intent(out) :: epz_ut(jm, ko)
  real(4),intent(out) :: epz_uw(jm, ko)

  call get_var_w_x_zm( pt_zm, u_v_x_zm, u_pt_dot_x_zm, &
       &              epz_uv, epz_ut, epz_uw )
  epz_uv(:,:) = -spread(rho,1,jm) * radius * spread(costbl,2,ko) &
       &      * epz_uv(:,:)
  epz_ut(:,:) = -spread(rho,1,jm) * radius * spread(costbl,2,ko) &
       &      * epz_ut(:,:)
  epz_uw(:,:) = -spread(rho,1,jm) * radius * spread(costbl,2,ko) &
       &      * epz_uw(:,:)

end subroutine epflux_z_uw





!
! Function
!   get form drag (one of the vertical component of EP flux)
!   
! Arguements (in)
!   p_pd   : pressure at the p+ surface (if j,k=fixed -> pt surface) [hPa]
!   z_pt   : geopotential height z(x,y,p+) [m]
!   p_sfc  : surface pressure [hPa]
!
! Arguements (out)
!   epz_form : form drag Fz(y,p+) [kg/s^2]
!
subroutine epflux_z_form( p_pd, z_pd, p_sfc, &
     &                    epz_form )
  use parameter, only : pai, radius
  use com_var, only : im, jm, ko, costbl
  implicit none
  real,intent(in)  :: p_pd(im, jm, ko)
  real,intent(in)  :: z_pd(im, jm, ko)
  real,intent(in)  :: p_sfc(im, jm)
  real,intent(out) :: epz_form(jm, ko)
  
  integer i, j, k
  real(8) :: dz, pdz
  real(4) :: p_pd2(im, jm, ko)

  p_pd2(:,:,:) = p_pd(:,:,:)
  do i=1, im
     do j=1, jm
        do k=1, ko
           if( p_pd2(i,j,k) >= p_sfc(i,j) ) then
              p_pd2(i,j,k) = p_sfc(i,j)
           end if
        end do
     end do
  end do
!  where( p_pd2 >= spread( p_sfc, 3, ko ) )
!     p_pd2 = spread( p_sfc, 3, ko )
!  end where
    
  do j=1, jm
     do k=1, ko
        
        pdz = 0 
        
        !*** calculate p * 2dz  (center difference)

        ! i=1
        dz = z_pd(1,j,k) - z_pd(im-1,j,k)
        pdz = pdz + ( p_pd2(im,j,k) * 100 ) * dz
        
        do i=2, im-1
           dz = z_pd(i+1,j,k) - z_pd(i-1,j,k) 
           pdz = pdz + ( p_pd2(i,j,k) * 100 ) * dz
        end do
        
        ! i=im
        dz = z_pd(2,j,k) - z_pd(im,j,k)
        pdz = pdz + ( p_pd2(1,j,k) * 100 ) * dz

        ! Iwasaki(1989) e.q.(D.3) M' -> Phi, p' -> p
        ! F_z = a cos(phi) (1/N) sum_{i=1}^N p_i (dz_i/dx)
        !     = a cos(phi) 1/(Ndx) sum_{i=1}^N p_i dz_i
        !     = 1/(4 pai) sum_{i=1}^N p_i 2dz_i
        epz_form(j,k) = radius * costbl(j) * &
             &          pdz / ( 4 * pai * radius * costbl(j) )
        
     end do
  end do

  return
end subroutine epflux_z_form



!
! Function
!   wave number deconposition of the form drag
!   
! Arguements (in)
!   p_pd   : pressure at the p+ surface (if j,k=fixed -> pt surface) [hPa]
!   z_pt   : geopotential height z(x,y,p+) [m]
!   p_sfc  : surface pressure [hPa]
!   wmax   : maximum wave number to decompose
!
! Arguements (out)
!   p_pd_wave     : #w component of p
!   z_pd_wave     : #w component of z
!   epz_form_wave : #w component of form drag
!
subroutine epflux_z_form_wave( p_pd, z_pd, p_sfc, wmax, &
     &                         p_pd_wave, z_pd_wave, &
     &                         epz_form_wave )
  use com_var, only : im, jm, km, ko, costbl
  implicit none
  real(4),intent(in)  :: p_pd(im, jm, ko)
  real(4),intent(in)  :: z_pd(im, jm, ko)
  real(4),intent(in)  :: p_sfc(im, jm)
  integer,intent(in)  :: wmax
  real(4),intent(out) :: p_pd_wave(wmax, im, jm, ko)
  real(4),intent(out) :: z_pd_wave(wmax,im, jm, ko)
  real(4),intent(out) :: epz_form_wave(wmax, jm, ko)

  integer :: w
  integer :: j, k

  ! for FFT
  real(8),allocatable :: wave_r(:), wave_i(:)   ! complex number
  real(8),allocatable :: kwave_r(:), kwave_i(:) ! complex number
  complex(8),allocatable :: wave(:)
  complex(8),allocatable :: kwave(:)
  real(8),parameter :: pi=3.1415926535898
  allocate( wave_r(0:im-1) )
  allocate( wave_i(0:im-1) )
  allocate( wave(0:im-1) )
  allocate( kwave_r(0:im-1) )
  allocate( kwave_i(0:im-1) )
  allocate( kwave(0:im-1) )

  !***** wavenumber decomposition *****!
  do k=1, ko
     do j=1, jm

!        wave_r(0:im-1) = p_pd(1:im,j,k)
!        wave_i(0:im-1) = 0
!        call fft_quick(im, 2*pi, wave_r, wave_i)

        wave(0:im-1) = cmplx( p_pd(1:im,j,k), 0.0 )
        call fft( im, +1, wave )


        do w=1, wmax

!           kwave_r = 0
!           kwave_i = 0
!           kwave_r(w) = wave_r(w)
!           kwave_i(w) = wave_i(w)
!           kwave_r(im-w) = wave_r(im-w)
!           kwave_i(im-w) = wave_i(im-w)

           kwave(:) = cmplx( 0.0, 0.0 )
           kwave(w) = wave(w)
           kwave(im-w) = wave(im-w)


!           call fft_quick(im, -2*pi, kwave_r, kwave_i)
!           kwave_r = kwave_r / im
!           kwave_i = kwave_i / im

           call fft( im, -1, kwave )
           kwave = kwave / im
           

!           p_pd_wave(w,1:im,j,k) = kwave_r(0:im-1)

           p_pd_wave(w,1:im,j,k) = real( kwave(0:im-1) )
           

        end do


!        wave_r(0:im-1) = z_pd(1:im,j,k)
!        wave_i(0:im-1) = 0
!        call fft_quick(im, 2*pi, wave_r, wave_i)

        wave(0:im-1) = cmplx( z_pd(1:im,j,k), 0.0 )
        call fft( im, +1, wave )

        do w=1, wmax

!           kwave_r = 0
!           kwave_i = 0
!           kwave_r(w) = wave_r(w)
!           kwave_i(w) = wave_i(w)
!           kwave_r(im-w) = wave_r(im-w)
!           kwave_i(im-w) = wave_i(im-w)

           kwave(:) = cmplx( 0.0, 0.0 )
           kwave(w) = wave(w)
           kwave(im-w) = wave(im-w)


!           call fft_quick(im, -2*pi, kwave_r, kwave_i)
!           kwave_r = kwave_r / im
!           kwave_i = kwave_i / im

           call fft( im, -1, kwave )
           kwave = kwave / im


!           z_pd_wave(w,1:im,j,k) = kwave_r(0:im-1)

           z_pd_wave(w,1:im,j,k) = real( kwave(0:im-1) )


        end do


     end do
  end do


  !***** epz_form *****!
  do w=1, wmax
     
     call epflux_z_form( p_pd_wave(w,:,:,:), z_pd_wave(w,:,:,:), p_sfc, &
          &              epz_form_wave(w,:,:) )
  end do

 
  deallocate( wave_r )
  deallocate( wave_i )
  deallocate( kwave_r )
  deallocate( kwave_i )
  deallocate( wave )
  deallocate( kwave )

end subroutine epflux_z_form_wave


