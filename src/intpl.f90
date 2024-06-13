!
! p_zm : p+ (zonal mean pressure) at the p+ levels 
! pd_p : 3-dimensional p+ at the pressure levels
!
subroutine intpl_pd_p( im, jm, km, ko, p_zm, pt, pt_zm, &
     &                 pd_p )
  implicit none
  integer,intent(in)  :: im, jm, km, ko
  real(4),intent(in)  :: p_zm(jm, ko)
  real(4),intent(in)  :: pt(im, jm, km)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(out) :: pd_p(im, jm, km)

  integer :: i, j, k, h
  integer :: kl, ku
  real(4) :: w

  do i=1, im
     do j=1, jm
        do k=1, km

           if( pt(i,j,k) > pt_zm(j,1) ) then
              kl = 1
              ku = 2
           else if( pt(i,j,k) < pt_zm(j,ko) ) then
              kl = ko - 1
              ku = ko
           else
              do h=1, ko-1
                 if( pt_zm(j,h+1) <= pt(i,j,k) &
                      &        .and. pt(i,j,k) <= pt_zm(j,h) ) then
                    kl = h
                    ku = h + 1
                    exit
                 end if
              end do
           end if
           if( i == 1 .and. j == 1)then

           endif
           ! linear
           w = ( pt(i,j,k) - pt_zm(j,kl) ) / ( pt_zm(j,ku) - pt_zm(j,kl) )
           pd_p(i,j,k) = w * p_zm(j,ku) + (1-w) * p_zm(j,kl)

           ! log(p)
           if( pd_p(i,j,k) <= 0 ) then
              w = ( pt(i,j,k) - pt_zm(j,kl) ) &
                   &  / ( pt_zm(j,ku) - pt_zm(j,kl) )
              pd_p(i,j,k) = p_zm(j,kl) * ( p_zm(j,ku) / p_zm(j,kl) )**w
           end if

        end do
     end do
  end do

end subroutine intpl_pd_p



!
! pd_ym : p++ (global mean pressure) at the p++ levels 
! pdd_pd : 2-dimensional p++ at the p+ levels
!
subroutine intpl_pdd_pd( jm, ko, pd_ym, pt_zm, pt_ym, &
     &                   pdd_pd )
  implicit none
  integer,intent(in)  :: jm, ko
  real(4),intent(in)  :: pd_ym(ko)
  real(4),intent(in)  :: pt_zm(jm, ko)
  real(4),intent(in)  :: pt_ym(ko)
  real(4),intent(out) :: pdd_pd(jm, ko)

  integer :: j, k, h
  integer :: kl, ku
  real(4) :: w

  do j=1, jm
     do k=1, ko
        
        if( pt_zm(j,k) > pt_ym(1) ) then
           kl = 1
           ku = 2
        else if( pt_zm(j,k) < pt_ym(ko) ) then
           kl = ko - 1
           ku = ko
        else
           do h=1, ko-1
              if( pt_ym(h+1) <= pt_zm(j,k) &
                   &      .and. pt_zm(j,k) <= pt_ym(h) ) then
                 kl = h
                 ku = h + 1
                 exit
              end if
           end do
        end if

        ! linear interpolation
        w = ( pt_zm(j,k) - pt_ym(kl) ) / ( pt_ym(ku) - pt_ym(kl) )
        pdd_pd(j,k) = w * pd_ym(ku) + (1-w) * pd_ym(kl)

        if( pdd_pd(j,k) <= 0 ) then
           ! log(p) interpolation
           w = ( pt_zm(j,k) - pt_ym(kl) ) &
             & / ( pt_ym(ku) - pt_ym(kl) )
           pdd_pd(j,k) = pd_ym(kl) * ( pd_ym(ku) / pd_ym(kl) )**w
        end if

     end do
  end do

end subroutine intpl_pdd_pd
