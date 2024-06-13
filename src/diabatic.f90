!
! Qz(diabatic A_Z generation rate
!   (global mean, without vertically integrated value)
!
subroutine diabatic_qz( jm, ko, q_ex_zm, pout, pdd_pd, &
     &                  qz_pdd )
  use parameter, only : cp, rkappa
  use biseki_y, only : biseki_y_biseki_y
  implicit none
  integer,intent(in)  :: jm, ko
  real(4),intent(in)  :: q_ex_zm(jm, ko)
  real(4),intent(in)  :: pout(ko)
  real(4),intent(in)  :: pdd_pd(jm, ko)
  real(4),intent(out) :: qz_pdd(ko)

  real(4) :: temp(jm, ko)
  integer :: j, k

  do j=1, jm
     do k=1, ko
        temp(j,k) = q_ex_zm(j,k) * cp * (1000.0**(-rkappa)) &
             &    * ( pout(k)**rkappa - pdd_pd(j,k)**rkappa )
     end do
  end do
 
  call biseki_y_biseki_y( temp, qz_pdd )

end subroutine diabatic_qz





subroutine diabatic_qe( im, jm, km, ko, q_3d, pin, pd_p, &
     &                  qe_zm )
  use parameter, only : rkappa
  use biseki, only : biseki_biseki
  implicit none
  integer,intent(in)  :: im, jm, km, ko
  real(4),intent(in)  :: q_3d(im, jm, km)
  real(4),intent(in)  :: pin(km)
  real(4),intent(in)  :: pd_p(im, jm, km)
  real(4),intent(out) :: qe_zm(jm, ko)

  real(4) :: temp(im, jm, km)
  integer :: i, j, k

  do i=1, im
     do j=1, jm
        do k=1, km
           temp(i,j,k) = q_3d(i,j,k) * ( 1 - (pd_p(i,j,k)/pin(k))**(rkappa) )
        end do
     end do
  end do

  call biseki_biseki( temp, qe_zm )

end subroutine diabatic_qe
