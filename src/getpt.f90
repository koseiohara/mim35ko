!
! Function
!   prepare for p -> p+ coordinate transformation
!
! Arguements (in)
!   im     : number of input data grid point in x-direction
!   km     : number of input data grid point in z-direction
!   ko     : number of output data grid point in z-direction
!   icount : current time step
!   pin    : pressure levels
!   pout   : p+ levels
!   p_sfc  : surface pressure
!   p_pds  : p+s
!   pt     : potential temperature
!   pt_sfc : potential temperature at the surface
!   pt_pds : potential temperature at the surface in the p+ coordinate
!
! Arguements (out)
!   dlev   : interpolation parameter (weight)
!   nlev   : interpolation parameter (grid number)
!   p_pd   : pressure at the p+ levels
!   p_zm   : zonal mean pressure at the p+ levels (it should approximate pout)
!   pt_zm  : (zonal mean) potential temperature at the p+ levels
!
! Note
!   p_zm will be close to pout after running this subroutine.
!
subroutine getpt( im, km, ko, icount, pin, pout, &
     &            p_sfc, p_pds, pt, pt_sfc, pt_pds, &  ! Input
     &            dlev, nlev, p_pd, p_zm, pt_zm )      ! Output
  use parameter, only : p_min, p_max, pt_min, pt_max
  implicit none
  integer,intent(in)  :: im, km, ko, icount
  real(4),intent(in)  :: pin(km)
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: p_sfc(im)
  real(4),intent(in)  :: p_pds
  real(4),intent(in)  :: pt_pds
  real(4),intent(in)  :: pt(im, km)
  real(4),intent(in)  :: pt_sfc(im)
  real(4),intent(out) :: dlev(im, ko)
  integer,intent(out) :: nlev(im, ko)
  real(4),intent(out) :: p_pd(im, ko)
  real(4),intent(out) :: p_zm(ko)
  real(4),intent(out) :: pt_zm(ko)

  real(4) :: pt_pd(im, ko)  ! potential temperature at the p+ levels
  real(4) :: pt_zm_old(ko)
  real(4) :: dr
  integer :: itmax = 10  ! maximum iteration number
  integer :: k, it
  !
  real(4) :: dlev_bst(im, ko)
  real(4) :: p_pd_bst(im, ko), p_zm_bst(ko), pt_zm_bst(ko)
  real(4) :: dr_all(ko), dr_old = 0.0
  integer :: nlev_bst(im, ko)

  ! get 1st approximation of pt_zm
  call getpt_pt1( im, 1, km, ko, pin, pout, pt, p_sfc, &
       &          pt_pd, nlev, dlev )
  pt_zm = sum( pt_pd, dim=1 ) / im

  ! get p_pd ( pressure at the standard p+ levels ) and p_zm ( = p+ )
  call getpt_p( im, 1, km, ko, nlev, dlev, pin, pt, pt_zm, &
       &        p_pd )
  p_zm = sum( p_pd, dim=1 ) / im         ! zonal mean
  where( pout > spread( p_pds, 1, ko ) )
     pt_zm = spread( pt_pds, 1, ko )     ! if underground, pt_zm = pt_pds
  end where

  !***** iteration (itmax or less times) *****!
  do it=1, itmax
!     write(*,*) it

     ! unstable -> stable
     do k=2, ko-1
        if( pt_zm(k) < pt_zm(ko) ) then
           if( k == ko) then
              pt_zm(k) = ( pt_zm(k-1) + pt_zm(ko) ) / 2
           else if( k /= ko) then
              pt_zm(k) = ( pt_zm(k-1) + pt_zm(k+1) ) / 2
           end if

        end if
     end do

     ! previous value
     pt_zm_old(:) = pt_zm(:)

     ! interpolate pt_zm using pt_zm_old in order to make p_zm close to pout
     call getpt_ptiter( ko, pout, pt_zm_old, p_zm, p_pds, pt_pds, &
          &             pt_zm )

     ! get interpolation parameter nlev & dlev
     call getpt_lev( im, km, ko, pin, pt, pt_sfc, pt_zm, p_sfc, &
          &          nlev, dlev )

     ! get p_pd ( pressure on standard p+ levels ) and p_zm ( = p+ )
     call getpt_p( im, 1, km, ko, nlev, dlev, pin, pt, pt_zm, &
          &        p_pd )
     p_zm = sum(p_pd, dim=1) / im         ! zonal mean
     where( pout > spread(p_pds,1,ko) )
        pt_zm = spread(pt_pds, 1, ko)     ! if underground, pt_zm = pt_pds
     end where

     ! check convergence condition and finish if appropriate
     do k=1, ko
        dr = abs( p_zm(k) / pout(k) - 1.0 )
        dr_all(k) = dr
!        if( dr > 0.001 ) exit
     end do

     ! If the results are betther than previous one, then store it
     if( it == 1 .or. sum(dr_all(:)) < dr_old ) then
        dlev_bst(:,:) = dlev(:,:)
        nlev_bst(:,:) = nlev(:,:)
        p_zm_bst(:)   = p_zm(:)
        pt_zm_bst(:)  = pt_zm(:)
        p_pd_bst(:,:) = p_pd(:,:)
        dr_old = sum(dr_all(:))
     end if

!     if( dr <= 0.001 ) exit
     if( maxval(dr_all(:)) <= 0.001 ) exit

  end do

  ! If above iteration is not converged, then use best values
  if( it == itmax + 1 ) then
     dlev(:,:) = dlev_bst(:,:)
     nlev(:,:) = nlev_bst(:,:)
     p_zm(:)   = p_zm_bst(:)
     pt_zm(:)  = pt_zm_bst(:)
     p_pd(:,:) = p_pd_bst(:,:)
  end if

  ! check p_zm and warn
  do k=1, ko
     dr = abs( p_zm(k) / pout(k) - 1.0 )
     if( p_pds > pout(k) .and. dr > 0.001 ) then
        write(65,*) "time=", icount, "k=", k, "pout=", pout(k), &
             &      "p_zm=", p_zm(k), "err=", dr, "p_pds=", p_pds
     end if
  end do

  ! check whether the atmosphere is stable or not
  call getpt_stable( 1, ko, pt_zm, p_pds, pt_sfc )
  call getp_stable( 1, ko, p_zm, p_pds)

  call check_range( im, 1, ko, p_pd, p_min, p_max, 'getpt()', 'p_pd' )
  call check_range( 1, 1, ko, pt_zm, pt_min, pt_max, 'getpt()', 'pt_zm' )

  return
end subroutine getpt



!
! Function
!   prepare for p+ -> p++ coordinate transformation
!
! Arguements (in)
!   jm     : number of input data grid point in y-direction
!   km     : number of input (p+) data grid point in z-direction
!   ko     : number of output (p++) data grid point in z-direction
!   icount : current time step
!   pin    : p+ levels (NOT equals to pin in com_var.f90)
!   pout   : p++ levels
!   alat   : latitude in radian
!   p_pds  : p+s
!   p_pdds : p++s
!   pt_zm  : (zonal mean) potential temperature at the p+ levels
!   pt_pds : potential temperature at the surface in the p+ coordinate
!   pt_pdds: potential temperature at the surface in the p++ coordinate
!
! Arguements (out)
!   dlev_y : interpolation parameter (weight)
!   nlev_y : interpolation parameter (grid number)
!   pd_pdd : p+ at the p++ levels
!   pd_ym  : global mean pressure at the p++ levels
!            (it should approximate pout)
!   pt_ym  : (global mean) potential temperature at the p++ levels
!
! Note
!   pd_ym will be close to pout after running this subroutine.
!
subroutine getpt_y( jm, km, ko, icount, pin, pout, alat, &
     &              p_pds, p_pdds, pt_zm, pt_pds, pt_pdds, &
     &              dlev_y, nlev_y, pd_pdd, pd_ym, pt_ym )
  use parameter, only : p_min, p_max, pt_min, pt_max
  implicit none
  integer,intent(in)  :: jm, km, ko, icount
  real(4),intent(in)  :: pin(km)
  real(4),intent(in)  :: alat(jm)
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: p_pdds
  real(4),intent(in)  :: pt_zm(jm, km)
  real(4),intent(in)  :: pt_pds(jm)
  real(4),intent(in)  :: pt_pdds
  real(4),intent(out) :: dlev_y(jm, ko)
  integer,intent(out) :: nlev_y(jm, ko)
  real(4),intent(out) :: pd_pdd(jm, ko)
  real(4),intent(out) :: pd_ym(ko)
  real(4),intent(out) :: pt_ym(ko)

  real(4) :: pt_pdd(jm, ko)  ! potential temperature at the p++ levels
  real(4) :: pt_ym_old(ko)
  real(4) :: dr
  integer :: itmax = 10  ! maximum iteration number
  integer :: k, it
  !
  real(4) :: dlev_y_bst(jm, ko)
  integer :: nlev_y_bst(jm, ko)
  real(4) :: pd_pdd_bst(jm, ko), pd_ym_bst(ko), pt_ym_bst(ko)
  real(4) :: dr_all(ko), dr_old = 0.0

  ! get 1st approximation of pt_ym
  call getpt_pt1( 1, jm, km, ko, pin, pout, pt_zm, p_pds, &
       &          pt_pdd, nlev_y, dlev_y )
  call integral_meridional( 1, jm, ko, alat, pt_pdd, &
       &                    pt_ym )

  ! get pd_pdd ( p+ at the standard p++ levels ) and pd_ym ( = p++ )
  call getpt_p( 1, jm, km, ko, nlev_y, dlev_y, pin, pt_zm, pt_ym, &
       &        pd_pdd )
  call integral_meridional( 1, jm, ko, alat, pd_pdd, &
       &                    pd_ym )
  where( pout > spread( p_pdds, 1, ko ) )
     pt_ym = spread( pt_pdds, 1, ko )    ! if underground, pt_ym = pt_pdds
  end where


  !***** iteration (itmax times or less) *****!
  do it=1, itmax
!     write(*,*) it

     ! unstable -> stable
     do k=2, ko
        if( pt_ym(k) < pt_ym(ko) ) then
           if( k == ko ) then
              pt_ym(k) = ( pt_ym(k-1) + pt_ym(ko) ) / 2
           else if( k /= ko ) then
              pt_ym(k) = ( pt_ym(k-1) + pt_ym(k+1) ) / 2
           end if

        end if
     end do

     ! previous value
     pt_ym_old(:) = pt_ym(:)

     ! interpolate pt_ym using pt_ym_old in order to make pd_ym close to pout
     call getpt_ptiter( ko, pout, pt_ym_old, pd_ym, p_pdds, pt_pdds, &
          &             pt_ym )

     ! get interpolation parameter nlev_y & dlev_y
     call getpt_lev( jm, km, ko, pin, pt_zm, pt_pds, pt_ym, p_pds, &
                     nlev_y, dlev_y )

     ! get pd_pdd ( p+ on standard p++ levels ) and p_ym ( = p++ )
     call getpt_p( 1, jm, km, ko, nlev_y, dlev_y, pin, pt_zm, pt_ym, &
          &        pd_pdd )
     call integral_meridional( 1, jm, ko, alat, pd_pdd, &
          &                    pd_ym )
     where( pout > spread( p_pdds, 1, ko ) )
        pt_ym = spread( pt_pdds, 1, ko )     ! if underground, pt_ym = pt_pdds
     end where

     ! check convergence condition and finish if appropriate
     do k=1, ko
        dr = abs( pd_ym(k) / pout(k) - 1.0 )
        dr_all(k) = dr
        ! if( dr > 0.001 ) exit
     end do

     ! If the results are betther than previous one, then store it
     if( it == 1 .or. sum(dr_all(:)) < dr_old ) then
        dlev_y_bst(:,:) = dlev_y(:,:)
        nlev_y_bst(:,:) = nlev_y(:,:)
        pd_ym_bst(:)  = pd_ym(:)
        pt_ym_bst(:) = pt_ym(:)
        pd_pdd_bst(:,:)   = pd_pdd(:,:)
        dr_old = sum(dr_all(:))
     end if

!     if( dr <= 0.001 ) exit
     if( maxval(dr_all(:)) <= 0.001 ) exit

  end do

  ! If above iteration is not converged, then use best values
  if( it == itmax + 1 ) then
     dlev_y(:,:) = dlev_y_bst(:,:)
     nlev_y(:,:) = nlev_y_bst(:,:)
     pd_ym(:)  = pd_ym_bst(:)
     pt_ym(:) = pt_ym_bst(:)
     pd_pdd(:,:)   = pd_pdd_bst(:,:)
  end if

  ! check pd_ym and warn
  do k=1, ko
     dr = abs( pd_ym(k) / pout(k) - 1.0 )
     if( p_pdds > pout(k) .and. dr > 0.001 ) then
        write(65,*) "time=", icount, "k=", k, "pout=", pout(k), &
             &      "pd_ym=", pd_ym(k), "err=", dr, "p_pdds=", p_pdds
     end if
  end do

  ! check whether the atmosphere is stable or not
  call getpt_stable( 1, ko, pt_ym, p_pdds, pt_pdds )
  call getp_stable( 1, ko, pd_ym, p_pdds)

  call check_range( 1, jm, ko, pd_pdd, p_min, p_max, 'getpt_y()', 'pd_pdd' )
  call check_range( 1, 1, ko, pt_ym, pt_min, pt_max, 'getpt()', 'pt_ym' )

  return
end subroutine getpt_y
