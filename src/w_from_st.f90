!
! get vertical velocity from mass streamfunction
!
subroutine w_from_st( jm, ko, sintbl, pout, p_pds, st_zm, &
     &                w_zm )
  use parameter, only : h0, radius, pai, grav
  implicit none
  integer,intent(in)  :: jm, ko
  real(4),intent(in)  :: sintbl(jm)
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: p_pds(jm)
  real(4),intent(in)  :: st_zm(jm, ko)
  real(4),intent(out) :: w_zm(jm, ko) 

  real(4) :: stzm
  integer :: j, k
  
  do k=1, ko
     
     stzm = 2.0 * st_zm(1,k) - st_zm(2,k)
     w_zm(1,k) = - grav / ( 2.0 * pai * radius ) &
          &    * ( stzm - st_zm(2,k) ) &
          &    / radius / ( ( sintbl(1) - sintbl(2) ) * 2 ) &
          &    * ( -h0 / (100.0*pout(k)) )
     
     do j=2, jm-1
        if( p_pds(j) > pout(k) ) then
           w_zm(j,k) = - grav / ( 2.0 * pai * radius ) &
                &    * ( st_zm(j-1,k) - st_zm(j+1,k) ) &
                &    / radius / ( sintbl(j-1) - sintbl(j+1) ) &
                &    * ( -h0 / (100.0*pout(k)) )
        else
           w_zm(j,k)=0.0
        end if
     end do
     
     stzm = 2.0 * st_zm(jm,k) - st_zm(jm-1,k)
     w_zm(jm,k) = - grav / ( 2.0 * pai * radius ) &
          &     * ( st_zm(jm-1,k) - stzm ) &
          &     / radius / ( ( sintbl(jm-1) - sintbl(jm) ) * 2 ) &
          &             * ( -h0 / (100.*pout(k)) )
     
  end do
  return
end subroutine w_from_st
