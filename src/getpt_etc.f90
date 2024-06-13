!
! Function
!   get 1st approximation of potential temperature
!
! Arguements (in)
!   im     : number of input data grid point in x-direction
!   jm     : number of input data grid point in y-direction
!   km     : number of input data grid point in z-direction
!   ko     : number of output data grid point in z-direction
!   pin    : input pressure levels
!   pout   : output p+ levels
!   pt     : potential temperature at the pressure levels [K]
!   p_sfc  : surface pressure [hPa]
!
! Arguements (out)
!   pt_pd  : potential temperature approximately at the p+ levels
!   nlev   : interpolation parameter
!   dlev   : interpolation parameter
!
! Note
!   -Exactly speaking, pt_zm must be potential temperature at the p+ levels.
!    Instead, potential temperature at the pressure levels is used
!    as a first approximation. Not only pt_zm but also nlev and dlev
!    will be modified in the other subroutines.
!
subroutine getpt_pt1( im, jm, km, ko, pin, pout, pt, p_sfc, &   ! Input
     &                pt_pd, nlev, dlev )                       ! Output
  use parameter
  implicit none
  integer,intent(in)  :: im, jm, km, ko
  real(4),intent(in)  :: pin(km)
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: pt(im, jm, km)
  real(4),intent(in)  :: p_sfc(im, jm)
  real(4),intent(out) :: pt_pd(im, jm, ko)  ! pt on p+ levels
  integer,intent(out) :: nlev(im, jm, ko)
  real(4),intent(out) :: dlev(im, jm, ko)

  integer :: i, j, k, n
  real(4) :: d

  nlev(:,:,:) = 0
  dlev(:,:,:) = 0.0

  do k=1, ko
     do j=1, jm
        do i=1, im

           ! determine n which is one level lower than or equal to k-level
           ! if pin == pout -> k=n
           n = 1
           if( k > 1 ) n = max( 1, nlev(i,j,k-1) )

           do while( pin(n+1) < pout(k) .and. n <= km-2 )
              n = n + 1
           end do

           ! Prepare for interpolation
           !   pin(n) < pout(k) < pin(n+1)
           if( p_sfc(i,j) > pin(n+1) ) then  ! in the atmosphere

              d = min(1.0, ( pout(k)-pin(n) ) / ( pin(n+1)-pin(n) ) )
              nlev(i,j,k) = n
              dlev(i,j,k) = d
              pt_pd(i,j,k) = pt(i,j,n) * (1.0-d) + pt(i,j,n+1) * d

           else  ! underground

              n = nlev(i,j,k-1)
              d = ( p_sfc(i,j)-pin(n) ) / ( pin(n+1)-pin(n) )
              nlev(i,j,k) = n
              dlev(i,j,k) = d
              pt_pd(i,j,k) = pt(i,j,n) * (1.0-d) + pt(i,j,n+1) * d

           end if

        end do
     end do
  end do

  return
end subroutine getpt_pt1




!
! getpt_lev() - get interpolation parameter
!
!   pt      pressure  variable
!
!  pt(n)     pin(n)    A(n)         <- standard pressure levels
!     |        |
!  pt_zm     (ppp)      A           <- standard p+ levels (interpolation point)
!     |        |
!  pt(n+1)   pin(n+1)  A(n+1)       <- standard pressure levels
!
! interpolation eq. : A = d A(n+1) + (1-d) A(n)
!   nlev <-> n
!   dlev <-> d (interpolation weight)
!
! if p+ -> p++ interpolation, i is replaced by j
! (e.g. im -> jm, pt_zm -> pt_ym)
!
subroutine getpt_lev( im, km, ko, pin, pt, pt_sfc, pt_zm, p_sfc,  &
     &                nlev, dlev )
  use parameter, only : rkappa
  implicit none
  integer,intent(in)  :: im, km, ko
  real(4),intent(in)  :: pin(km)
  real(4),intent(in)  :: pt(im,km)
  real(4),intent(in)  :: pt_sfc(im)
  real(4),intent(in)  :: pt_zm(ko)
  real(4),intent(in)  :: p_sfc(im)
  integer,intent(out) :: nlev(im, ko)
  real(4),intent(out) :: dlev(im, ko)

  real(8) :: ppp, dppp
  integer :: i, k, n

  do k=1, ko
     do i=1, im

        !********** get pressure ppp corresponding to pt_zm **********!

        if( pt_zm(k) > pt_sfc(i) ) then  ! if pt_zm(k) is in the atmosphere

           !***** find n which satisfies pt(n) > pt_zm(k) > pt(n+1)
           n = 1
           do while( pt(i,n+1) > pt_zm(k) .and. n <= ko-2 )
              n = n + 1
           end do

           !***** interpolate to ppp *****!

           ! if pt_zm(k) is below the lowermost grid point,
           ! linear interpolation using pt(i,km) and pt_sfc(i)
           if( pt_zm(k) <= pt(i,km) ) then
              ppp = pin(km) + ( pt_zm(k) - pt(i,km) ) &
                   &        / ( pt_sfc(i) - pt(i,km) ) &
                   &        * ( p_sfc(i) - pin(km) )

           ! if pin(n+1) is under the ground,
           ! linear interpolation using pt(i,n) and pt_sfc(i)
           else if( pin(n+1) > p_sfc(i) ) then
              ppp = pin(n) + ( pt_zm(k) - pt(i,n) ) &
                   &       / ( pt_sfc(i) - pt(i,n) ) &
                   &       * ( p_sfc(i) - pin(n) )

           !  if pt_zm(k) is above the uppermost grid point,
           !  use definition of potential temperature and assume T=const.
           else if( pt_zm(k) > pt(i,1) ) then
              ppp = ( pt(i,1) / pt_zm(k) )**(1/rkappa) * pin(1)

           else ! log(p) interpolation
              dppp = log(pin(n)) + ( pt_zm(k) - pt(i,n) ) &
                   &             / ( pt(i,n+1) - pt(i,n) ) &
                   &             * ( log(pin(n+1)) - log(pin(n)) )
              ppp = exp(dppp)
!              ppp = pin(n) + ( pt_zm(k) - pt(i,n) ) / &
!                   &         ( pt(i,n+1) - pt(i,n) ) &
!                   &       * ( pin(n+1) - pin(n) )

           end if

        else ! if pt_zm(k) is under the ground

           !***** find n which satisfies pin(n) < p_sfc(k) < pin(n+1)
           n = 1
           do while( pin(n+1) < p_sfc(i) .and. n <= ko-2 )
              n = n + 1
           end do

           ! use surface value (from the definition)
           ppp = p_sfc(i)

        end if


        nlev(i,k) = n

        !********** calculate interpolation weight **********!
        ! log(p) interpolation
        dlev(i,k) = ( log(ppp) - log(pin(n)) ) &
             &    / ( log(pin(n+1)) - log(pin(n)) )
!        dlev(i,k) = ( ppp - pin(n) ) / ( pin(n+1) - pin(n) )
!        write(*,*) i, k, ppp, pin(n), pin(n+1), dlev(i,k)

        ! if pt_zm(k) is above the uppermost grid point,
        ! use linear interpolation to avoid too small dlev.
        if( dlev(i,k) < 0 .and. ppp < pin(1) ) then
           dlev(i,k) = ( ppp - pin(n) ) / ( pin(n+1) - pin(n) )
        end if

     end do
  end do

  return
end subroutine getpt_lev




!
! Function
!   get pressure at the p+ levels
!     or get p+ at the p++ levels
!
! if im != 1,
!   pin   : pressure levels
!   pt    : potential temperature at the pressure levels
!   pt_zm : potential temperature at the p+ levels
!   p_pd  : pressure at the p+ levels
!
! if im == 1, be careful...
!   pin   : p+ levels
!   pt    : potential temperature at the p+ levels
!   pt_zm : potential temperature at the p++ levels
!   p_pd  : p+ at the p++ levels
!
subroutine getpt_p( im, jm, km, ko, nlev, dlev, pin, pt, pt_zm, &  ! Input
     &              p_pd )                                         ! Output
  implicit none
  integer,intent(in)  :: im, jm, ko, km
  integer,intent(in)  :: nlev(im, jm, ko)
  real(4),intent(in)  :: dlev(im, jm, ko)
  real(4),intent(in)  :: pin(km)
  real(4),intent(in)  :: pt(im, jm, km)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(out) :: p_pd(im, jm, ko)
  real(4) :: d
  integer :: i, j, k, l

  do i=1, im
     do j=1, jm
        do k=1, ko

           ! linear interpolation/extrapolation
           l = nlev(i,j,k)
           d = dlev(i,j,k)
!           p_pd(i,j,k) = (1.0-d) * pin(l) + d * pin(l+1)
           ! nlev and dlev is evaluated with log(p)-interpolation,
           ! so below is better than above.
           p_pd(i,j,k) = exp( (1.0-d) * log( pin(l) ) + d * log( pin(l+1) ) )

           ! if p<=0 -> extrapolate again with log(p) interpolation
           if( p_pd(i,j,k) <= 0 ) then
              l = nlev(i,j,k)
              d = ( pt_zm(j,k) - pt(i,j,l) ) / ( pt(i,j,l+1) - pt(i,j,l) )
              p_pd(i,j,k) = pin(l) * ( pin(l+1) / pin(l) )**d
           end if

           ! if longitudinal variation is too large
!           if( p_pd(i,j,k) <= pin(1)*1.0e-5 ) then
!              if( i > 2 ) then
!                 p_pd(i,j,k) = p_pd(i-1,j,k)
!              else
!                 p_pd(i,j,k) = pt_zm(j,k)
!              end if
!           end if

        end do
     end do
  end do

  return
end subroutine getpt_p



!
! Function
!   interpolate pt_zm using pt_zm_old in order to make p_zm close to pout
!
! pt_zm_old -> pt_zm
!
subroutine getpt_ptiter( ko, pout, pt_zm_old, p_zm, p_pds, pt_pds, &
     &                   pt_zm )
  implicit none
  integer,intent(in)  :: ko
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: pt_zm_old(ko)
  real(4),intent(in)  :: p_zm(ko)
  real(4),intent(in)  :: p_pds, pt_pds
  real(4),intent(out) :: pt_zm(ko)
  integer :: k, kk
  !
  real(4),allocatable :: p_zm_itr(:), pt_zm_itr(:)
  real(4) :: pl, pu, ptl, ptu
  real(8) :: a

  !-- prepare pressure/potential temperature arrays including surface values
  allocate( p_zm_itr(1:ko+1) )
  allocate( pt_zm_itr(1:ko+1) )
  !
  p_zm_itr(ko+1) = p_pds
  p_zm_itr(1:ko) = p_zm(1:ko)
  pt_zm_itr(ko+1) = pt_pds
  pt_zm_itr(1:ko) = pt_zm_old(1:ko)
  do k = 1, ko+1
     if( p_zm_itr(k) > p_pds ) then
        p_zm_itr(k) = p_pds
     end if
     if( pt_zm_itr(k) < pt_pds ) then
        pt_zm_itr(k) = pt_pds
     end if
  end do


  pt_zm(:) = -999.0             ! substitute unrealistic values for checking

  kk = 1
  kloop: do k=2, ko+1                  ! loop for p_zm_itr
11   continue

     !-- find layers that sandwitch pout(kk)
     a = ( p_zm_itr(k-1) - pout(kk) )*( p_zm_itr(k) - pout(kk) )
     if( pout(kk) < p_zm_itr(1)  .or. a <= 0.0 ) then

        if( pout(kk) < p_zm_itr(1) ) then ! use k=1 and 2
           pl = p_zm_itr(2)
           pu = p_zm_itr(1)
           ptl = pt_zm_itr(2)
           ptu = pt_zm_itr(1)
        else                    ! a <= 0.0
           pl = p_zm_itr(k)
           pu = p_zm_itr(k-1)
           ptl = pt_zm_itr(k)
           ptu = pt_zm_itr(k-1)
        end if

        pt_zm(kk) = ptl + ( ptl - ptu ) &
             &          / (  pl - pu ) &
             &          * ( pout(kk) - pl )

        kk = kk + 1

        if(kk == ko+1) then     ! check whether all pt_zm are updated or not
           exit kloop
        end if
        goto 11

     end if

  end do kloop

  !-- surface check
  do kk = 1, ko
     if( pt_zm(kk) < 0 ) then
        if( pout(kk) < p_pds ) then
           pt_zm(kk) = pt_zm_old(kk) + ( pt_zm_old(kk) - pt_zm_old(kk-1) ) &
                &                    / ( p_zm(kk) - p_zm(kk-1) ) &
                &                    * ( pout(kk) - p_zm(kk) )
        else
           pt_zm(kk) = pt_pds
        end if
     end if
  end do

  deallocate( p_zm_itr, pt_zm_itr )

  return
end subroutine getpt_ptiter


!
! getpt_stable() - check whether the atmosphere is stable or not
!
subroutine getpt_stable( jm, ko, pt, p_sfc, pt_sfc )
  implicit none
  integer,intent(in) :: jm, ko
  real(4),intent(in) :: pt(jm, ko), p_sfc(jm), pt_sfc(jm)

  integer :: j, k

  do j=1, jm
     do k=2, ko

        if( pt(j,k) > pt(j,k-1) ) then
           write(6,*) "error in getpt_stable() : pt is instable j=", j, &
                &             "k=", k, k-1, "pt=", pt(j,k), pt(j,k-1), &
                &             "p_sfc=", p_sfc(j), "pt_sfc=", pt_sfc(j)
           write(6,*) pt(j,:), pt_sfc(j)
           stop
        end if

     end do
  end do

  return
end subroutine getpt_stable


!
! getp_stable() - check whether the atmosphere is stable or not
!
subroutine getp_stable( jm, ko, p, p_sfc)
  implicit none
  integer,intent(in) :: jm, ko
  real(4),intent(in) :: p(jm, ko), p_sfc(jm)

  integer :: j, k

  do j=1, jm
     do k=2, ko

        if( p(j,k) < p(j,k-1) ) then
           write(6,*) "error in getp_stable() : p is instable j=", j, &
                &             "k=", k, k-1, "p=", p(j,k), p(j,k-1), &
                &             "p_sfc=", p_sfc(j)
           stop
        end if

     end do
  end do

  return
end subroutine getp_stable
