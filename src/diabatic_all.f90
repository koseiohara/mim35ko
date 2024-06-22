module diabatic_all

    use parameter, only : rkappa, cp
    use com_var  , only : im, jm, km, ko, pin, pout
    use mim_var  , only : pd_p, pdd_pd, p_pdds
    use biseki, only : biseki_biseki
    use biseki_y, only : biseki_y_biseki_y

    implicit none

    private
    public :: diabaticHeating

    contains


    ! Arguments :
    !     1. 3d distribution of diabatic heating
    !     2. zonal mean diabatic heating
    !     3. diabatic heating to the zonal mean state
    !     4. Generation rate of eddy available potential energy
    !     5. Generation rate of zonal available potential energy
    subroutine diabaticHeating(heating_3d, heating_zm, heating_gz_zm, eddy_generation, zonal_generation_gmean)
        real(4), intent(in)  :: heating_3d(im,jm,km)
        real(4), intent(out) :: heating_zm(jm,ko)
        real(4), intent(out) :: heating_gz_zm(jm,ko)
        real(4), intent(out) :: eddy_generation(jm,ko)
        real(4), intent(out) :: zonal_generation_gmean(1)

        real(4) :: heating_exner_zm(jm,ko)
        real(4) :: heating_pdd(ko)
        real(4) :: zonal_generation(ko)

        !call heating_ZonalMean(heating_3d(1:im,1:jm,1:km), &  !! IN
        !                     & heating_zm(1:jm,1:ko)       )  !! OUT


        ! Zonal Mean Diabatic Heanting
        call biseki_biseki(heating_3d(1:im,1:jm,1:km), &  !! IN
                         & heating_zm(1:jm,1:ko)       )  !! OUT

        call heating_per_Exner(heating_3d(1:im,1:jm,1:km), &  !! IN
                             & heating_exner_zm(1:jm,1:ko) )  !! OUT

        call heating_ZonalMeanState(heating_exner_zm(1:jm,1:ko), &  !! IN
                                  & heating_gz_zm(1:jm,1:ko)     )  !! OUT

        !call heating_Eddy_simple(heating_zm(1:jm,1:ko)   , &  !! IN
        !                       & heating_gz_zm(1:jm,1:ko), &  !! IN
        !                       & eddy_generation(1:jm,1:ko))  !! OUT

        call heating_Eddy_highPrecision(heating_3d(1:im,1:jm,1:km), &  !! IN
                                      & eddy_generation(1:jm,1:ko)  )  !! OUT

        !call heating_GroundState(heating_exner_zm(1:jm,1:ko), &  !! IN
        !                       & heating_pdd(1:ko)            )  !! OUT

        !call heating_Zonal_simple(heating_gz_zm(1:jm,1:ko), &  !! IN
        !                        & heating_pdd(1:ko)       , &  !! IN
        !                        & zonal_generation(1:ko)    )  !! OUT

        call heating_Zonal_highPrecision(heating_exner_zm(1:jm,1:ko), &  !! IN
                                       & zonal_generation(1:ko)       )  !! OUT

        call integral_p(1                      , &  !! IN   size in lat-direction
                      & ko                     , &  !! IN   size in p-direction
                      & pout(1:ko)             , &  !! IN   output levels
                      & p_pdds(1)              , &  !! IN   p_dagger_dagger at the surface
                      & zonal_generation(1:ko) , &  !! IN   vertical profile of parameter
                      & zonal_generation_gmean(1))  !! OUT  vertically integrated parameter

    end subroutine diabaticHeating


    ! Zonal mean diabatic heating bar{Q^*}
    ! This subroutine is not used because this is just a wrapper routine for biseki_biseki()
    subroutine heating_ZonalMean(heating_3d, heating_zm)
        real(4), intent(in)  :: heating_3d(im,jm,km)
        real(4), intent(out) :: heating_zm(jm,ko)

        call biseki_biseki(heating_3d(1:im,1:jm,1:km), &  !! IN
                         & heating_zm(1:jm,1:ko)       )  !! OUT

    end subroutine heating_ZonalMean


    ! Zonal mean diabatic heating divided by Exner function bar{[Q/Exner(p)]^*}
    subroutine heating_per_Exner(heating_3d, heating_exner_zm)
        real(4), intent(in)  :: heating_3d(im,jm,km)
        !real(4), intent(out) :: heating_gz_zm(jm,ko)
        real(4), intent(out) :: heating_exner_zm(jm,ko)

        real(4) :: heating_exner_3d(im,jm,km)

        integer :: k

        do k = 1, km
            heating_exner_3d(1:im,1:jm,k) = heating_3d(1:im,1:jm,k) / (cp * (pin(k)*0.001_4)**rkappa)
        enddo

        ! Zonal Mean
        call biseki_biseki(heating_exner_3d(1:im,1:jm,1:km), &  !! IN
                         & heating_exner_zm(1:jm,1:ko)       )  !! OUT

    end subroutine heating_per_Exner


    ! Diabatic heating to the zonal mean state bar{[Q/Exner(p)]^*}Exner(p_dagger)
    subroutine heating_ZonalMeanState(heating_exner_zm, heating_gz_zm)
        real(4), intent(in)  :: heating_exner_zm(jm,ko)
        real(4), intent(out) :: heating_gz_zm(jm,ko)

        integer :: k

        do k = 1, ko
            heating_gz_zm(1:jm,k) = heating_exner_zm(1:jm,k) * cp * (pout(k)*0.001_4)**rkappa
        enddo

    end subroutine heating_ZonalMeanState


    ! Generation rate of eddy available potential energy bar{Q^*} - bar{[Q/Exner(p)]^*}Exner(p_dagger)
    ! This should not be used beacause of bad precision
    subroutine heating_Eddy_simple(heating_zm, heating_gz_zm, eddy_generation)
        real(4), intent(in)  :: heating_zm(jm,ko)
        real(4), intent(in)  :: heating_gz_zm(jm,ko)
        real(4), intent(out) :: eddy_generation(jm,ko)

        eddy_generation(1:jm,1:ko) = heating_zm(1:jm,1:ko) - heating_gz_zm(1:jm,1:ko)

    end subroutine heating_Eddy_simple


    ! Generation rate of eddy available potential energy bar{Q^*} - bar{[Q/Exner(p)]^*}Exner(p_dagger)
    ! This should not be used beacause of bad precision
    subroutine heating_Eddy_highPrecision(heating_3d, eddy_generation)
        real(4), intent(in)  :: heating_3d(im,jm,km)
        real(4), intent(out) :: eddy_generation(jm,ko)

        real(4) :: work_eddy_generation(im,jm,km)

        integer :: k

        do k = 1, km
            work_eddy_generation(1:im,1:jm,k) = heating_3d(1:im,1:jm,k) * (1._4 - (pd_p(1:im,1:jm,k) / pin(k))**rkappa)
        enddo

        call biseki_biseki(work_eddy_generation(1:im,1:jm,1:km), &  !! IN
                         & eddy_generation(1:jm,1:ko)            )  !! OUT

    end subroutine heating_Eddy_highPrecision


    ! Diabatic heating to the ground state bar{bar{[Q/Exner(o)]partial p / partial p_dagger_dagger}}Exner(p_dagger_dagger)
    subroutine heating_GroundState(heating_exner_zm, heating_pdd)
        real(4), intent(in)  :: heating_exner_zm(jm,ko)
        real(4), intent(out) :: heating_pdd(ko)

        real(4) :: work_heating_pdd(jm,ko)

        work_heating_pdd(1:jm,1:ko) = heating_exner_zm(1:jm,1:ko) * cp * (pdd_pd(1:jm,1:ko)*0.001_4)**rkappa

        call biseki_y_biseki_y(work_heating_pdd(1:jm,1:ko), &  !! IN
                             & heating_pdd(1:ko)            )  !! OUT

    end subroutine heating_GroundState


    ! Generation rate of zonal mean available potential energy
    ! bar{[Q/Exner(p)]^*}Exner(p_dagger) - bar{bar{[Q/Exner(o)]partial p / partial p_dagger_dagger}}Exner(p_dagger_dagger)
    ! This should not be used beacause bad precision
    subroutine heating_Zonal_simple(heating_gz_zm, heating_pdd, zonal_generation)
        real(4), intent(in)  :: heating_gz_zm(jm,ko)
        real(4), intent(in)  :: heating_pdd(ko)
        real(4), intent(out) :: zonal_generation(ko)

        call biseki_y_biseki_y(heating_gz_zm(1:jm,1:ko), &  !! IN
                             & zonal_generation(1:ko)    )  !! OUT

        zonal_generation(1:ko) = zonal_generation(1:ko) - heating_pdd(1:ko)

    end subroutine heating_Zonal_simple


    subroutine heating_Zonal_highPrecision(heating_exner_zm, zonal_generation)
        real(4), intent(in)  :: heating_exner_zm(jm,ko)
        real(4), intent(out) :: zonal_generation(ko)

        real(4) :: work_zonal_generation(jm,ko)

        integer :: k

        do k = 1, ko
            work_zonal_generation(1:jm,k) = heating_exner_zm(1:jm,k) * &
                                          & cp * ((pout(k)*0.001_4)**rkappa - (pdd_pd(1:jm,k)*0.001_4)**rkappa)
        enddo

        call biseki_y_biseki_y(work_zonal_generation(1:jm,1:ko), &  !! IN
                             & zonal_generation(1:ko)            )  !! OUT

    end subroutine heating_Zonal_highPrecision


end module diabatic_all

