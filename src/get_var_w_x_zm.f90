!
! ret_var_w : (var' w_dagger')_zm
! ret_var_v : var' v' component
! ret_var_t : var' pt_dot' component
!
subroutine get_var_w_x_zm( pt_zm, var_v_x_zm, var_pt_dot_x_zm, &
     &                     ret_var_v, ret_var_t, ret_var_w )
  use parameter, only : grav, radius
  use com_var, only : jm, ko, pout, alat, costbl, rho
  implicit none
  real(4),intent(in)  :: pt_zm(jm, ko) 
  real(4),intent(in)  :: var_v_x_zm(jm, ko)
  real(4),intent(in)  :: var_pt_dot_x_zm(jm, ko)
  real(4),intent(out) :: ret_var_v(jm, ko)
  real(4),intent(out) :: ret_var_t(jm, ko)
  real(4),intent(out) :: ret_var_w(jm, ko)
  
  integer j, k
  integer :: kl, ku


  !***** ret_var_v *****!
  do k=1, ko
     ret_var_v(1,k)  = 0.0

     do j=2, jm-1

        kl = k - 1
        ku = k + 1
        if( k == 1 )  kl = 1
        if( k == ko ) ku = ko
        
        if( abs(pt_zm(j,ku)-pt_zm(j,kl)) <= 0.01 ) then
           ret_var_v(j,k) = 0.0
        else
           ret_var_v(j,k) = var_v_x_zm(j,k) / ( rho(k) * grav * radius ) &
                &      * ( pt_zm(j+1,k) - pt_zm(j-1,k) ) &
                &      / ( alat(j+1) - alat(j-1) ) &
                &      * ( pout(ku) - pout(kl) ) * 100 &
                &      / ( pt_zm(j,ku) - pt_zm(j,kl) )
        end if
        
     end do

     ret_var_v(jm,k) = 0.0
  end do
  
  
  !***** ret_var_t *****!
  do j=1, jm
     
     ret_var_t(j,1) = 0.0
     
     do k=2, ko
        kl = k - 1
        ku = k + 1
        if( k == 1 )  kl = 1
        if( k == ko ) ku = ko

        if( abs(pt_zm(j,ku)-pt_zm(j,kl)) < 0.01 ) then
           ret_var_t(j,k) = 0.0
        else
           ret_var_t(j,k) = -var_pt_dot_x_zm(j,k) / ( rho(k) * grav ) &
                &      * ( pout(ku) - pout(kl) ) * 100 &
                &      / ( pt_zm(j,ku) - pt_zm(j,kl) )
        end if
        
     end do

  end do


  !***** ret_var_w *****!
  ret_var_w(:,:) = ret_var_v(:,:) + ret_var_t(:,:)

  return 
end subroutine get_var_w_x_zm
