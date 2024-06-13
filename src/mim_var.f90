!#########################################################
!
!  Variable list for mim.f90
!
!    "use variable" should be specified only in mim.f90
!
!#########################################################
module mim_var
  implicit none

  !********** input **********!

  !***** pressure-surface *****!
  real(4),allocatable :: u(:,:,:)     ! zonal wind [m/s]
  real(4),allocatable :: v(:,:,:)     ! meridional Wind [m/s]
  real(4),allocatable :: t(:,:,:)     ! temperature [K]
  real(4),allocatable :: z(:,:,:)     ! geopotential height [m]
  real(4),allocatable :: omega(:,:,:) ! p-velocity [Pa/s]

  !***** surface *****!
  real(4),allocatable :: p_sfc(:,:)   ! surface pressure [hPa]
  real(4),allocatable :: msl(:,:)     ! mean sea level pressure [hPa]
  real(4),allocatable :: t_sfc(:,:)   ! surface temperature [K]
  real(4),allocatable :: alt(:,:)     ! surface altitude [m]


  !********** derived **********!

  !*** 3D basic variables on the standard pressure levels ***!
  real(4),allocatable :: pt(:,:,:)     ! potential temperature
  real(4),allocatable :: pt_dot(:,:,:) ! D(pt)/Dt

  !*** 2D zonal mean basic variables on the standard p+ levels ***!
  real(4),allocatable :: p_zm (:,:)    ! exact p+ [hPa] (approx. equal pout)
  real(4),allocatable :: u_zm (:,:)    ! zonal wind [m/s]
  real(4),allocatable :: v_zm(:,:)     ! meridional wind [m/s]
  real(4),allocatable :: pt_zm(:,:)    ! potential temperature [K]
  real(4),allocatable :: t_dagger(:,:) ! temperature dagger (from pt_zm) [K]
!  real(4),allocatable :: t_zm(:,:)     ! temperature (from t) [K]
  real(4),allocatable :: st_zm(:,:)    ! mass streamfunction [kg/s]
  real(4),allocatable :: w_zm (:,:)    ! vertical velocity (from st_zm) [m/s]
  real(4),allocatable :: z_zm(:,:)     ! geopotential height [m]
  real(4),allocatable :: pt_dot_zm(:,:)! D(pt)/Dt [K/s]

  !*** correlarion ***!
  real(4),allocatable :: u_u_zm(:,:)      ! (u^2)_zm [m^2/s^2]
  real(4),allocatable :: u_u_x_zm(:,:)    ! (u'u')_zm [m^2/s^2]
  real(4),allocatable :: v_v_zm(:,:)
  real(4),allocatable :: v_v_x_zm(:,:)
  real(4),allocatable :: u_v_zm(:,:)
  real(4),allocatable :: u_v_x_zm(:,:)
  real(4),allocatable :: u_u_v_zm(:,:)     ! (u^2 v)_zm
  real(4),allocatable :: v_v_v_zm(:,:)     ! (v^3)_zm
  real(4),allocatable :: u_pt_dot_zm(:,:)
  real(4),allocatable :: u_pt_dot_x_zm(:,:)
  real(4),allocatable :: v_pt_dot_zm(:,:)
  real(4),allocatable :: v_pt_dot_x_zm(:,:)
  real(4),allocatable :: u_u_pt_dot_zm(:,:)
  real(4),allocatable :: v_v_pt_dot_zm(:,:)
  real(4),allocatable :: v_ke_zm(:,:)      ! (v Ke)_zm
  real(4),allocatable :: pt_dot_ke_zm(:,:) ! (pt_dot Ke)_zm

  !*** EP flux & G Flux ***!
  real(4),allocatable :: epy(:,:)      ! Fy [kg/s^2]
  real(4),allocatable :: depy(:,:)     ! divF due to epy [m/s^2]
  real(4),allocatable :: epz_form(:,:) ! Fz (Form Drag) [kg/s^2]
  real(4),allocatable :: depz_form(:,:)! divF due to epz_form [m/s^2]
  real(4),allocatable :: epz_uv(:,:)   ! part of Fz (u'w') [kg/s^2]
  real(4),allocatable :: depz_uv(:,:)  ! divF due to epz_uv [m/s^2]
  real(4),allocatable :: epz_ut(:,:)   ! part of Fz (u'w') [kg/s^2]
  real(4),allocatable :: depz_ut(:,:)  ! divF due to epz_ut [m/s^2]
  real(4),allocatable :: epz_uw(:,:)   ! Fz (u'w') = epz_uv + epz_ut [kg/s^2]
  real(4),allocatable :: depz_uw(:,:)  ! divF due to epz_uw [m/s^2]
  real(4),allocatable :: epz(:,:)      ! Fz (Total) [kg/s^2]
  real(4),allocatable :: depz(:,:)     ! divF due to epz [m/s^2]
  real(4),allocatable :: divf(:,:)     ! EP Flux Divergence [m/s^2]
  real(4),allocatable :: gy(:,:)       ! Gy [kg/s^2]
  real(4),allocatable :: dgy(:,:)      ! dGy/dy [m/s^2]
  real(4),allocatable :: gz(:,:)       ! Gz [kg/s^2]
  real(4),allocatable :: dgz(:,:)      ! dGz/dy [m/s^2]

  !*** energy ***!
  real(4),allocatable :: kz_zm(:,:)    ! zonal kinetic energy [m^2/s^2]
  real(4),allocatable :: ke_zm(:,:)    ! eddy kinetic energy  [m^2/s^2]
  real(4),allocatable :: pz_zm(:,:)    ! potential energy [m^2/s^2]
  real(4),allocatable :: ae_zm_vint(:) ! eddy available potential energy [J/m^2]
  real(4),allocatable :: ae_total_zm(:,:) ! (not recommended to use)
  real(4),allocatable :: ae_s_zm(:)    ! surface correction
  real(4),allocatable :: az_zm_vint(:) ! (not recommended to use)
  real(4)             :: az_gmean(1)   ! zonal available potential energy [J/m^2]

  !*** energy conversion (for global mean) ***!
  real(4),allocatable :: c_az_kz(:,:)   ! C(Az,Kz) [W/m^2]
  real(4),allocatable :: c_az_kz_modify(:,:)
  real(4),allocatable :: c_kz_ae(:,:)   ! C(Kz,Ae) [W/m^2]
  real(4),allocatable :: c_kz_ae_u(:,:)
  real(4),allocatable :: c_kz_ae_v(:,:)
  real(4),allocatable :: c_ae_ke(:,:)   ! C(Ae,Ke) [W/m^2]
  real(4),allocatable :: c_ae_ke_u(:,:)
  real(4),allocatable :: c_ae_ke_v(:,:)
  real(4),allocatable :: c_kz_ke(:,:)   ! C(Kz,ke) [W/m^2]
  real(4),allocatable :: c_kz_ke_uy(:,:)
  real(4),allocatable :: c_kz_ke_uz(:,:)
  real(4),allocatable :: c_kz_ke_vy(:,:)
  real(4),allocatable :: c_kz_ke_vz(:,:)
  real(4),allocatable :: c_kz_ke_tan(:,:)
  real(4),allocatable :: c_kz_w(:,:)   ! C(Kz,W) [W/m^2]

  !*** for p -> pd ***!
  integer,allocatable :: nlev(:,:,:)
  real(4),allocatable :: dlev(:,:,:)
  real(4),allocatable :: p_pd(:,:,:)
  real(4),allocatable :: x_pd(:,:,:)
  real(4),allocatable :: xint_zm(:,:)
  real(4),allocatable :: pt_sfc(:,:)
  real(4),allocatable :: pt_pds(:)    ! pt_{xmin}
  real(4),allocatable :: pd_p(:,:,:)  ! pd(x,y,p)

  !*** for pd -> pdd ***!
  integer,allocatable :: nlev_y(:,:)
  real(4),allocatable :: dlev_y(:,:)
  real(4),allocatable :: pd_pdd(:,:)
  real(4),allocatable :: pd_ym(:)
  real(4),allocatable :: pt_ym(:)
  real(4)             :: pt_pdds(1)   ! pt_pds_{ymin}
  real(4),allocatable :: pdd_pd(:,:)  ! pdd(y,pd)


  real(4),allocatable :: dz_dlat_zm(:,:)
  real(4),allocatable :: v_dz_dlat_zm(:,:)
  real(4),allocatable :: u_dz_dlon_zm(:,:)
  real(4),allocatable :: temp_vint(:) !xxx
  real(4),allocatable :: phi_dagger(:,:)
  real(4),allocatable :: phi_dagger_y(:,:)
  real(4),allocatable :: p_dphi_dt(:,:)
  real(4),allocatable :: p_dz_dt_zm(:,:)
  real(4),allocatable :: p_dz_dt(:,:,:)
  real(4),allocatable :: divz_tzm(:,:)
  real(4),allocatable :: divphi_t(:,:)
  real(4),allocatable :: dwdt(:,:)
  real(4),allocatable :: uuv_tmp(:,:)
!  real(4),allocatable :: d_u_epy(:,:)
  real(4),allocatable :: d_u_epz(:,:) !local
!  real(4),allocatable :: dkedt_uvu_y(:,:)

  !*** energy conversion (local) ***!
  real(4),allocatable :: dkzdt_vkz(:,:)
  real(4),allocatable :: dkzdt_wkz(:,:)
  real(4),allocatable :: dkedt_uy(:,:)
  real(4),allocatable :: dkedt_vy(:,:)
  real(4),allocatable :: dkedt_uz(:,:)
  real(4),allocatable :: dkedt_vz(:,:)
  real(4),allocatable :: dkedt_vke(:,:)
  real(4),allocatable :: dkedt_wke(:,:)
  real(4),allocatable :: dpedt_vt(:,:)
  real(4),allocatable :: dpedt_wt(:,:)

  real(4),allocatable :: z_pd(:,:,:)
  real(4),allocatable :: pt_past(:,:,:)
  real(4),allocatable :: u_past(:,:,:)
  real(4),allocatable :: v_past(:,:,:)
  real(4),allocatable :: omega_past(:,:,:)
  real(4),allocatable :: phi_dagger_past(:,:)
  real(4),allocatable :: z_pd_past(:,:,:)
  real(4),allocatable :: p_pd_past(:,:,:)

  !*** wavenumber decomposition ***!
  real(4),allocatable :: epz_wave(:,:,:)
  real(4),allocatable :: z_pt_wave(:,:,:,:)
  real(4),allocatable :: p_pt_wave(:,:,:,:)

  !*** diabatic heating ***!
  real(4),allocatable :: q_3d(:,:,:)    ! 3D diabatic heating
  real(4),allocatable :: q_zm(:,:)      ! zonal mean diabatic heating
  real(4),allocatable :: q_ex_3d(:,:,:)
  real(4),allocatable :: q_ex_zm(:,:)
  real(4),allocatable :: qgz_zm(:,:)
  real(4),allocatable :: qz_pdd(:)
  real(4)             :: qz_gmean(1)  ! Az generation by diabatic heating
  real(4),allocatable :: qe_zm(:,:)   ! Ae generation by diabatic heating

  !*** others ***!
  real(4),allocatable :: p_pds(:)    ! p+s: zonal mean surface pressure
  real(4)             :: p_pdds(1)   ! p++s: global mean surface pressure [hPa]
  real(4),allocatable :: work(:,:,:) ! temporal work space (im,jm,km)


contains
  subroutine mim_var_ini(im, jm, km, ko, wmax)
    integer,intent(in) :: im   ! x( East -> West ) direction grid
    integer,intent(in) :: jm   ! y( North -> South ) direction grid
    integer,intent(in) :: km   ! input data z direction levels
    integer,intent(in) :: ko   ! output data z direction levels
    integer,intent(in) :: wmax ! maximum wave number

    allocate( u(im,jm,km) )
    allocate( v(im,jm,km) )
    allocate( t(im,jm,km) )
    allocate( z(im,jm,km) )
    allocate( omega(im,jm,km) )

    allocate( p_sfc(im,jm) )
    allocate( msl(im,jm) )
    allocate( t_sfc(im,jm) )
    allocate( alt(im,jm) )

    allocate( pt(im,jm,km) )
    allocate( pt_dot(im,jm,km) )

    allocate( p_zm(jm,ko) )
    allocate( u_zm(jm,ko) )
    allocate( v_zm(jm,ko) )
    allocate( pt_zm(jm,ko) )
    allocate( t_dagger(jm,ko) )
!    allocate( t_zm(jm,ko) )
    allocate( st_zm(jm,ko) )
    allocate( w_zm(jm,ko) )
    allocate( z_zm(jm,ko) )
    allocate( pt_dot_zm(jm,ko) )

    allocate( u_u_zm(jm,ko) )
    allocate( u_u_x_zm(jm,ko) )
    allocate( v_v_zm(jm,ko) )
    allocate( v_v_x_zm(jm,ko) )
    allocate( u_v_zm(jm,ko) )
    allocate( u_v_x_zm(jm,ko) )
    allocate( u_u_v_zm(jm,ko) )
    allocate( v_v_v_zm(jm,ko) )
    allocate( u_pt_dot_zm(jm,ko) )
    allocate( u_pt_dot_x_zm(jm,ko) )
    allocate( v_pt_dot_zm(jm,ko) )
    allocate( v_pt_dot_x_zm(jm,ko) )
    allocate( u_u_pt_dot_zm(jm,ko) )
    allocate( v_v_pt_dot_zm(jm,ko) )
    allocate( v_ke_zm(jm,ko) )
    allocate( pt_dot_ke_zm(jm,ko) )

    allocate( epy(jm,ko) )
    allocate( depy(jm,ko) )
    allocate( epz_form(jm,ko) )
    allocate( depz_form(jm,ko) )
    allocate( epz_uv(jm,ko) )
    allocate( depz_uv(jm,ko) )
    allocate( epz_ut(jm,ko) )
    allocate( depz_ut(jm,ko) )
    allocate( epz_uw(jm,ko) )
    allocate( depz_uw(jm,ko) )
    allocate( epz(jm,ko) )
    allocate( depz(jm,ko) )
    allocate( divf(jm,ko) )
    allocate( gy(jm,ko) )
    allocate( dgy(jm,ko) )
    allocate( gz(jm,ko) )
    allocate( dgz(jm,ko) )

    allocate( kz_zm(jm,ko) )
    allocate( ke_zm(jm,ko) )
    allocate( pz_zm(jm,ko) )
    allocate( ae_zm_vint(jm) )
    allocate( ae_total_zm(jm,ko) )
    allocate( ae_s_zm(jm) )
    allocate( az_zm_vint(jm) )

    allocate( c_az_kz(jm,ko) )
    allocate( c_az_kz_modify(jm,ko) )
    allocate( c_kz_ae(jm,ko) )
    allocate( c_kz_ae_u(jm,ko) )
    allocate( c_kz_ae_v(jm,ko) )
    allocate( c_ae_ke(jm,ko) )
    allocate( c_ae_ke_u(jm,ko) )
    allocate( c_ae_ke_v(jm,ko) )
    allocate( c_kz_ke(jm,ko) )
    allocate( c_kz_ke_uy(jm,ko) )
    allocate( c_kz_ke_uz(jm,ko) )
    allocate( c_kz_ke_vy(jm,ko) )
    allocate( c_kz_ke_vz(jm,ko) )
    allocate( c_kz_ke_tan(jm,ko) )
    allocate( c_kz_w(jm,ko) )


    allocate( nlev(im,jm,ko) )
    allocate( dlev(im,jm,ko) )
    allocate( p_pd(im,jm,ko) )
    allocate( x_pd(im,jm,km) )
    allocate( xint_zm(jm,ko) )
    allocate( pt_sfc(im,jm) )
    allocate( pt_pds(jm) )
    allocate( pd_p(im,jm,km) )

    allocate( dlev_y(jm,ko) )
    allocate( nlev_y(jm,ko) )
    allocate( pd_pdd(jm,ko) )
    allocate( pd_ym(ko) )
    allocate( pt_ym(ko) )
    allocate( pdd_pd(jm,ko) )




    allocate( dz_dlat_zm(jm,ko) )
    allocate( v_dz_dlat_zm(jm,ko) )
    allocate( u_dz_dlon_zm(jm,ko) )
    allocate( temp_vint(jm) )
    allocate( phi_dagger(jm,ko) )
    allocate( phi_dagger_y(jm,ko) )
    allocate( p_dphi_dt(jm,ko) )
    allocate( p_dz_dt_zm(jm,ko) )
    allocate( p_dz_dt(im,jm,km) )
    allocate( divz_tzm(jm,ko) )
    allocate( divphi_t(jm,ko) )
    allocate( dwdt(jm,ko) )
    allocate( uuv_tmp(jm,ko) )
!    allocate(d_u_epy(jm,ko))
    allocate( d_u_epz(jm,ko) )
!    allocate(dkedt_uvu_y(jm,ko))

    allocate( dkzdt_vkz(jm,ko) )
    allocate( dkzdt_wkz(jm,ko) )
    allocate( dkedt_uy(jm,ko) )
    allocate( dkedt_vy(jm,ko) )
    allocate( dkedt_uz(jm,ko) )
    allocate( dkedt_vz(jm,ko) )
    allocate( dkedt_vke(jm,ko) )
    allocate( dkedt_wke(jm,ko) )
    allocate( dpedt_vt(jm,ko) )
    allocate( dpedt_wt(jm,ko) )


    allocate( z_pd(im,jm,ko) )
    allocate( pt_past(im,jm,km) )
    allocate( u_past(im,jm,km) )
    allocate( v_past(im,jm,km) )
    allocate( omega_past(im,jm,km) )
    allocate( phi_dagger_past(jm,ko) )
    allocate( z_pd_past(im,jm,ko) )
    allocate( p_pd_past(im,jm,ko) )

    if( wmax >= 1 ) then
       allocate( epz_wave(wmax,jm,ko) )
       allocate( z_pt_wave(wmax,im,jm,ko) )
       allocate( p_pt_wave(wmax,im,jm,ko) )
    end if

    allocate( q_3d(im,jm,km) )
    allocate( q_zm(jm,ko) )
    allocate( q_ex_3d(im,jm,km) )
    allocate( q_ex_zm(jm,ko) )
    allocate( qgz_zm(jm,ko) )
    allocate( qz_pdd(ko) )
    allocate( qe_zm(jm,ko) )

    allocate( p_pds(jm) )
    allocate( work(im,jm,km) )


  end subroutine mim_var_ini



  subroutine mim_var_end()
    deallocate( u )
    deallocate( v )
    deallocate( t )
    deallocate( z )
    deallocate( omega )

    deallocate( p_sfc )
    deallocate( msl )
    deallocate( t_sfc )
    deallocate( alt )

    deallocate( pt )
    deallocate( pt_dot )

    deallocate( p_zm )
    deallocate( u_zm )
    deallocate( v_zm )
    deallocate( pt_zm )
    deallocate( t_dagger )
!    deallocate( t_zm )
    deallocate( st_zm )
    deallocate( w_zm )
    deallocate( z_zm )
    deallocate( pt_dot_zm )

    deallocate( u_u_zm )
    deallocate( u_u_x_zm )
    deallocate( v_v_zm )
    deallocate( v_v_x_zm )
    deallocate( u_v_zm )
    deallocate( u_v_x_zm )
    deallocate( u_u_v_zm )
    deallocate( v_v_v_zm )
    deallocate( u_pt_dot_zm )
    deallocate( u_pt_dot_x_zm )
    deallocate( v_pt_dot_zm )
    deallocate( v_pt_dot_x_zm )
    deallocate( u_u_pt_dot_zm )
    deallocate( v_v_pt_dot_zm )
    deallocate( v_ke_zm )
    deallocate( pt_dot_ke_zm )

    deallocate( epy )
    deallocate( depy )
    deallocate( epz_form )
    deallocate( depz_form )
    deallocate( epz_uv )
    deallocate( depz_uv )
    deallocate( epz_ut )
    deallocate( depz_ut )
    deallocate( epz_uw )
    deallocate( depz_uw )
    deallocate( epz )
    deallocate( depz )
    deallocate( divf )
    deallocate( gy )
    deallocate( dgy )
    deallocate( gz )
    deallocate( dgz )

    deallocate( kz_zm )
    deallocate( ke_zm )
    deallocate( pz_zm )
    deallocate( ae_zm_vint )
    deallocate( ae_total_zm )
    deallocate( ae_s_zm )
    deallocate( az_zm_vint )

    deallocate( c_az_kz )
    deallocate( c_az_kz_modify )
    deallocate( c_kz_ae )
    deallocate( c_kz_ae_u )
    deallocate( c_kz_ae_v )
    deallocate( c_ae_ke )
    deallocate( c_ae_ke_u )
    deallocate( c_ae_ke_v )
    deallocate( c_kz_ke )
    deallocate( c_kz_ke_uy )
    deallocate( c_kz_ke_uz )
    deallocate( c_kz_ke_vy )
    deallocate( c_kz_ke_vz )
    deallocate( c_kz_ke_tan )
    deallocate( c_kz_w )


    deallocate( nlev )
    deallocate( dlev )
    deallocate( p_pd )
    deallocate( x_pd )
    deallocate( xint_zm )
    deallocate( pt_sfc )
    deallocate( pt_pds )
    deallocate( pd_p )

    deallocate( dlev_y )
    deallocate( nlev_y )
    deallocate( pd_pdd )
    deallocate( pd_ym )
    deallocate( pt_ym )
    deallocate( pdd_pd )


    deallocate( dz_dlat_zm )
    deallocate( v_dz_dlat_zm )
    deallocate( u_dz_dlon_zm )
    deallocate( temp_vint )
    deallocate( phi_dagger )
    deallocate( phi_dagger_y )
    deallocate( p_dphi_dt )
    deallocate( p_dz_dt_zm )
    deallocate( p_dz_dt )
    deallocate( divz_tzm )
    deallocate( divphi_t )
    deallocate( dwdt )
    deallocate( uuv_tmp )
!    deallocate(d_u_epy)
    deallocate( d_u_epz )
!    deallocate(dkedt_uvu_y)

    deallocate( dkzdt_vkz )
    deallocate( dkzdt_wkz )
    deallocate( dkedt_uy )
    deallocate( dkedt_vy )
    deallocate( dkedt_uz )
    deallocate( dkedt_vz )
    deallocate( dkedt_vke )
    deallocate( dkedt_wke )
    deallocate( dpedt_vt )
    deallocate( dpedt_wt )


    deallocate( z_pd )
    deallocate( pt_past )
    deallocate( u_past )
    deallocate( v_past )
    deallocate( omega_past )
    deallocate( phi_dagger_past )
    deallocate( z_pd_past )
    deallocate( p_pd_past )

    if( allocated( epz_wave ) ) deallocate( epz_wave )
    if( allocated( z_pt_wave ) ) deallocate( z_pt_wave )
    if( allocated( p_pt_wave ) ) deallocate( p_pt_wave )

    deallocate( q_3d )
    deallocate( q_zm )
    deallocate( q_ex_3d )
    deallocate( q_ex_zm )
    deallocate( qgz_zm )
    deallocate( qz_pdd )
    deallocate( qe_zm )

    deallocate( p_pds )
    deallocate( work )


  end subroutine mim_var_end


end module mim_var
