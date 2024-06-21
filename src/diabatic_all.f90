module diabatic_all

    use parameter, only : rkappa, cp
    use com_var  , only : im, jm, km, ko, pin, pout
    use mim_var  , only : pdd_pd
    use biseki, only : biseki_biseki
    use biseki_y, only : biseki_y_biseki_y

    implicit none

    private
    public :: 

    contains

    ! Zonal mean diabatic heating bar{Q^*}
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
            heating_exner_3d(1:im,1:jm,k) = heating_3d(1:im,1:jm,k) / (cp * (pin(k)*0.001)**rkappa)
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
            heating_gz_zm(1:jm,k) = heating_exner_zm(1:jm,k) * cp * (pout(k)*0.001)**rkappa
        enddo

    end subroutine heating_ZonalMeanState


    ! Generation rate of eddy available potential energy bar{Q^*} - bar{[Q/Exner(p)]^*}Exner(p_dagger)
    subroutine heating_Eddy(heating_zm, heating_gz_zm, eddy_generation)
        real(4), intent(in)  :: heating_zm(jm,ko)
        real(4), intent(in)  :: heating_gz_zm(jm,ko)
        real(4), intent(out) :: eddy_generation(jm,ko)

        eddy_generation(1:jm,1:ko) = heating_zm(jm,ko) - heating_gz_zm(jm,ko)

    end subroutine heating_Eddy


    ! Diabatic heating to the ground state bar{bar{[Q/Exner(o)]partial p / partial p_dagger_dagger}}Exner(p_dagger_dagger)
    subroutine heating_GroundState(heating_exner_zm, heating_pdd)
        real(4), intent(in)  :: heating_exner_zm(jm,ko)
        real(4), intent(out) :: heating_pdd(ko)

        real(4) :: work_heating_pdd(jm,ko)

        work_heating_pdd(1:jm,1:ko) = heating_exner_zm(1:jm,1:ko) * cp * (pdd_pd(1:jm,1:ko)*0.001)**rkappa

        call biseki_y_biseki_y(work_heating_pdd(1:jm,1:ko), &  !! IN
                             & heating_pdd(1:ko)            )  !! OUT

    end subroutine heating_GroundState


    ! Generation rate of zonal mean available potential energy
    ! bar{[Q/Exner(p)]^*}Exner(p_dagger) - bar{bar{[Q/Exner(o)]partial p / partial p_dagger_dagger}}Exner(p_dagger_dagger)
    subroutine heating_Zonal(heating_gz_zm, heating_pdd, zonal_generation)
        real(4), intent(in)  :: heating_gz_zm(jm,ko)
        real(4), intent(in)  :: heating_pdd(ko)
        real(4), intent(out) :: zonal_generation(ko)

        call biseki_y_biseki_y(heating_gz_zm(1:jm,1:ko), &  !! IN
                             & zonal_generation(1:ko)    )  !! OUT

        zonal_generation(1:ko) = zonal_generation(1:ko) + heating_pdd(1:ko)

    end subroutine heating_Zonal


end module diabatic_all

