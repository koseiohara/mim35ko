module az_latprof

    use parameter, only : rkappa, grav, cp
    !use parameter, only : pai
    use com_var  , only : jm, ko, alat, pout
    use mim_var  , only : pd_pdd

    implicit none

    private
    public :: az_latprof_vint, energy_az_highPrecision

    contains

    !
    ! Vertically integrated A_Z
    !
    ! Note:
    !   It is not used for global mean A_Z (because of the accuracy)
    !
    subroutine az_latprof_vint( p_pds, p_pdds, pd_pdd, pt_pdds, pt_ym, &
         &                     az_zm_vint )
      use parameter, only : rkappa, grav, cp, pai
      use com_var, only : jm, ko, alat, pout
      implicit none
      real(4),intent(in)  :: p_pds(jm)
      real(4),intent(in)  :: p_pdds(1)
      real(4),intent(in)  :: pd_pdd(jm, ko)         ! p+ at p++ surfaces
      real(4),intent(in)  :: pt_pdds(1)             ! Potential Temperature at the surface p++ level
      real(4),intent(in)  :: pt_ym(ko)              ! Potential Temperature at each p++ level
      real(4),intent(out) :: az_zm_vint(jm)

      real(4) :: pd_ym(ko)
      real(4) :: integ(jm, ko)
      integer :: j, k
      real(4) :: const
      real(4) :: az_modify(jm)

      write(*,*) 'Marker1'

      ! pd_ym : global mean p+ at the p++ levels
      !         pd_ym must be almost equal to standard p++ levels 
      !         except under the ground.
      call integral_meridional( 1, jm, ko, alat, pd_pdd, &
           &                    pd_ym )
      
      ! get integrand
      call energy_az_highPrecision(pd_ym(1:ko)   , &  !! IN
                                 & integ(1:jm,1:ko))  !! OUT

      ! integrate with pt
      call integral_pt_ym( jm, ko, pout, p_pdds, pt_ym, pt_pdds, integ, &
           &               az_zm_vint )

      ! lower boundary modification
      !   it is proportional to pt_ymin * ( p+s^{kappa+1} - p++s^{kappa+1} ).
      az_modify(:) = const * pt_pdds(1) &
           &       * ( p_pds(:)**(rkappa+1) - p_pdds(1)**(rkappa+1) )

      az_zm_vint(:) = az_zm_vint(:) + az_modify(:)

    end subroutine az_latprof_vint


    subroutine energy_az_simple(pd_ym, integrand)
        real(4), intent(in)  :: pd_ym(ko)
        real(4), intent(out) :: integrand(jm,ko)

        real(4) :: expansion(jm)                      !! p_dagger - p_dagger_dagger

        integer, parameter :: coeff = cp / (100000._4**rkappa * (1._4 + rkappa) * grav)
        integer :: k

        do k = 1, ko
            expansion(1:jm) = (pd_pdd(1:jm,k) - pd_ym(k)) / pd_ym(k)
            integrand(1:jm,k)  &
            & = coeff*(pd_ym(k)*100._4)**(1+rkappa) * (                                                              &  !! coefficient
            &   (rkappa+1._4)*rkappa/2._4                                                                            &  !! 2st term
            &                       *expansion(1:jm)*expansion(1:jm)                                                 &
            & + (rkappa+1._4)*rkappa*(rkappa-1._4)/6._4                                                              &  !! 3nd term
            &                       *expansion(1:jm)*expansion(1:jm)*expansion(1:jm)                                 &
            & + (rkappa+1._4)*rkappa*(rkappa-1._4)*(rkappa-2._4)/24._4                                               &  !! 4rd term
            &                       *expansion(1:jm)*expansion(1:jm)*expansion(1:jm)*expansion(1:jm)                 &
            & + (rkappa+1._4)*rkappa+(rkappa-1._4)*(rkappa-2._4)*(rkappa-3._4)/120._4                                &  !! 5th term
            &                       *expansion(1:jm)*expansion(1:jm)*expansion(1:jm)*expansion(1:jm)*expansion(1:jm) &
            & )
        enddo

    end subroutine energy_az_simple


    subroutine energy_az_highPrecision(pd_ym, integrand)
        real(4), intent(in)  :: pd_ym(ko)
        real(4), intent(out) :: integrand(jm,ko)

        real(4) :: ratio(jm)

        integer, parameter :: coeff = cp / (100000._4**rkappa * (1._4 + rkappa) * grav)
        integer, parameter :: C2 = (rkappa+1)*rkappa/2._4
        integer, parameter :: C3 = (rkappa+1)*rkappa*(rkappa-1)/6._4
        integer, parameter :: C4 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)/24._4
        integer, parameter :: C5 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)*(rkappa-3)/120._4
        integer :: k

        do k = 1, ko

            ratio(1:jm) = pd_pdd(1:jm,k) / pd_ym(k)

            integrand(1:jm,k) = C5 * ratio(1:jm)**5
            integrand(1:jm,k) = integrand(1:jm,k) + ( C4 - 5*C5                ) * ratio(1:jm)**4
            integrand(1:jm,k) = integrand(1:jm,k) + ( C3 - 4*C4 + 10*C5        ) * ratio(1:jm)**3
            integrand(1:jm,k) = integrand(1:jm,k) + ( C2 - 3*C3 +  6*C4 - 10*C5) * ratio(1:jm)**2
            integrand(1:jm,k) = integrand(1:jm,k) + (-C2 + 2*C3 -  3*C4 +  4*C5) * ratio(1:jm)
            
            integrand(1:jm,k) = integrand(1:jm,k) * coeff * (pd_ym(k)*100._4)**(1+rkappa)

        enddo

    end subroutine energy_az_highPrecision


end module az_latprof

