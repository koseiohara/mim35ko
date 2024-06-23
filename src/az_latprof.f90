module az_latprof

    use parameter, only : rkappa, grav, cp
    use parameter, only : pai -> pi
    use com_var  , only : jm, ko, alat, pout

    implicit none

    private
    public ::

    contains

    !
    ! Vertically integrated A_Z
    !
    ! Note:
    !   It is not used for global mean A_Z (because of the accuracy)
    !
    subroutine energy_az_vint( p_pds, p_pdds, pd_pdd, pt_pdds, pt_ym, &
         &                     az_zm_vint )
      use parameter, only : rkappa, grav, cp, pai
      use com_var, only : jm, ko, alat, pout
      implicit none
      real(4),intent(in)  :: p_pds(jm)
      real(4),intent(in)  :: p_pdds(1)
      real(4),intent(in)  :: pd_pdd(jm, ko)         ! p+ at 
      real(4),intent(in)  :: pt_pdds(1)             ! Potential Temperature at the surface p++ level
      real(4),intent(in)  :: pt_ym(ko)              ! Potential Temperature at each p++ level
      real(4),intent(out) :: az_zm_vint(jm)

      real(4) :: pd_ym(ko)
      real(4) :: integ(jm, ko)
      integer :: j, k
      real(4) :: const
      real(4) :: az_modify(jm)

      ! pd_ym : global mean p+ at the p++ levels
      !         pd_ym must be almost equal to standard p++ levels 
      !         except under the ground.
      call integral_meridional( 1, jm, ko, alat, pd_pdd, &
           &                    pd_ym )
      
      ! get integrand
      const = cp * (1.0e+5)**(-rkappa) / (1+rkappa) / grav
      do k=1, ko
         do j=1, jm
            integ(j,k) = const * ( ( pd_pdd(j,k)*100 )**(1+rkappa) &
                 &                - ( pd_ym(k)*100)**(1+rkappa) )
         end do
      end do

      ! integrate with pt
      call integral_pt_ym( jm, ko, pout, p_pdds, pt_ym, pt_pdds, integ, &
           &               az_zm_vint )

      ! lower boundary modification
      !   it is proportional to pt_ymin * ( p+s^{kappa+1} - p++s^{kappa+1} ).
      az_modify(:) = const * pt_pdds(1) &
           &       * ( p_pds(:)**(rkappa+1) - p_pdds(1)**(rkappa+1) )

      az_zm_vint(:) = az_zm_vint(:) + az_modify(:)

    end subroutine energy_az_vint


end module az_latprof

