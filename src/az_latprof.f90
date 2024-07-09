module az_latprof

    use parameter, only : rkappa, grav, cp
    !use parameter, only : pai
    use com_var  , only : jm, ko, alat, pout
    use mim_var  , only : pd_pdd

    implicit none

    private
    public :: az_latprof_vint

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
      real(4), parameter :: const = cp * (1.0e+5)**(-rkappa) / (1+rkappa) / grav
      real(4) :: az_modify(jm)

      ! pd_ym : global mean p+ at the p++ levels
      !         pd_ym must be almost equal to standard p++ levels 
      !         except under the ground.
      call integral_meridional( 1, jm, ko, alat, pd_pdd, &
           &                    pd_ym )

      !call energy_az_simple(pd_ym(1:ko)   , &  !! IN
      !                    & integ(1:jm,1:ko))  !! OUT
     
      ! get integrand
      !call energy_az_decompose_tight(pd_ym(1:ko)   , &  !! IN
      !                             & integ(1:jm,1:ko))  !! OUT

      call energy_az_decompose_loose(pd_ym(1:ko)   , &  !! IN
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

        integer, parameter :: coeff = cp / (1.0e5**rkappa * (1._4 + rkappa) * grav)
        real(4), parameter :: C2 = (rkappa+1)*rkappa/2._4
        real(4), parameter :: C3 = (rkappa+1)*rkappa*(rkappa-1)/6._4
        real(4), parameter :: C4 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)/24._4
        real(4), parameter :: C5 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)*(rkappa-3)/120._4
        real(4), parameter :: C6 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)*(rkappa-3)*(rkappa-4)/720._4
        integer :: k

        do k = 1, ko
            expansion(1:jm) = (pd_pdd(1:jm,k) - pd_ym(k)) / pd_ym(k)

            integrand(1:jm,k) = C6*expansion(1:jm)**6
            integrand(1:jm,k) = integrand(1:jm,k) + C5*expansion(1:jm)**5
            integrand(1:jm,k) = integrand(1:jm,k) + C4*expansion(1:jm)**4
            integrand(1:jm,k) = integrand(1:jm,k) + C3*expansion(1:jm)**3
            integrand(1:jm,k) = integrand(1:jm,k) + C2*expansion(1:jm)**2

            integrand(1:jm,k) = integrand(1:jm,k) * coeff * (pd_ym(k)*100._4)**(1._4+rkappa)

        enddo

    end subroutine energy_az_simple


    subroutine energy_az_decompose_tight(pd_ym, integrand)
        real(4), intent(in)  :: pd_ym(ko)
        real(4), intent(out) :: integrand(jm,ko)

        real(4) :: ratio(jm)

        real(4), parameter :: coeff = cp / (1.0e5**rkappa * (1._4 + rkappa) * grav)
        real(4), parameter :: C2 = (rkappa+1)*rkappa/2._4
        real(4), parameter :: C3 = (rkappa+1)*rkappa*(rkappa-1)/6._4
        real(4), parameter :: C4 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)/24._4
        real(4), parameter :: C5 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)*(rkappa-3)/120._4
        real(4), parameter :: C6 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)*(rkappa-3)*(rkappa-4)/720._4
        integer :: k

        do k = 1, ko

            ratio(1:jm) = pd_pdd(1:jm,k) / pd_ym(k)

            integrand(1:jm,k) = C6 * ratio(1:jm)**6
            integrand(1:jm,k) = integrand(1:jm,k) + (                                C5 -  6._4*C6) * ratio(1:jm)**5
            integrand(1:jm,k) = integrand(1:jm,k) + (                     C4 -  5._4*C5 + 15._4*C6) * ratio(1:jm)**4
            integrand(1:jm,k) = integrand(1:jm,k) + (           C3 - 4._4*C4 + 10._4*C5 - 20._4*C6) * ratio(1:jm)**3
            integrand(1:jm,k) = integrand(1:jm,k) + ( C2 - 3._4*C3 + 6._4*C4 - 10._4*C5 + 15._4*C6) * ratio(1:jm)**2
            integrand(1:jm,k) = integrand(1:jm,k) + (-C2 + 2._4*C3 - 3._4*C4 +  4._4*C5 -  5._4*C6) * ratio(1:jm)
            
            integrand(1:jm,k) = integrand(1:jm,k) * coeff * (pd_ym(k)*100._4)**(1+rkappa)

        enddo

    end subroutine energy_az_decompose_tight


    subroutine energy_az_decompose_loose(pd_ym, integrand)
        real(4), intent(in)  :: pd_ym(ko)
        real(4), intent(out) :: integrand(jm,ko)

        real(4) :: ratio(jm)

        real(4), parameter :: coeff = cp / (1.0e5**rkappa * (1._4 + rkappa) * grav)
        real(4), parameter :: C2 = (rkappa+1)*rkappa/2._4
        real(4), parameter :: C3 = (rkappa+1)*rkappa*(rkappa-1)/6._4
        real(4), parameter :: C4 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)/24._4
        real(4), parameter :: C5 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)*(rkappa-3)/120._4
        real(4), parameter :: C6 = (rkappa+1)*rkappa*(rkappa-1)*(rkappa-2)*(rkappa-3)*(rkappa-4)/720._4
        integer :: k

        do k = 1, ko

            ratio(1:jm) = pd_pdd(1:jm,k) / pd_ym(k)

            integrand(1:jm,k) = C6 * ratio(1:jm)**6
            integrand(1:jm,k) = integrand(1:jm,k) + (                                      C5 -  6._4*C6) * ratio(1:jm)**5
            integrand(1:jm,k) = integrand(1:jm,k) + (                           C4 -  5._4*C5 + 15._4*C6) * ratio(1:jm)**4
            integrand(1:jm,k) = integrand(1:jm,k) + (                 C3 - 4._4*C4 + 10._4*C5 - 20._4*C6) * ratio(1:jm)**3
            integrand(1:jm,k) = integrand(1:jm,k) + (       C2 - 3._4*C3 + 6._4*C4 - 10._4*C5 + 15._4*C6) * ratio(1:jm)**2
            integrand(1:jm,k) = integrand(1:jm,k) + (- 2._4*C2 + 3._4*C3 - 4._4*C4 +  5._4*C5 -  6._4*C6) * ratio(1:jm)
            integrand(1:jm,k) = integrand(1:jm,k) + (       C2 -      C3 +      C4 -       C5 +       C6)
            
            integrand(1:jm,k) = integrand(1:jm,k) * coeff * (pd_ym(k)*100._4)**(1+rkappa)

        enddo

    end subroutine energy_az_decompose_loose


end module az_latprof

