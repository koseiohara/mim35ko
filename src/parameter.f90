!#########################################################
!
!  Parameter  
!
!#########################################################
module parameter
  implicit none

  real(4),parameter :: rkappa  = 0.286        ! r/cp
  real(4),parameter :: grav    = 9.81         ! m/sec**2
  real(4),parameter :: pai     = 3.14159      !
  real(4),parameter :: radian  = pai / 180.0  ! const. for angle
  real(4),parameter :: radius  = 6.371e+6     ! earth radius
  real(4),parameter :: h0      = 7000.0       ! scale height
  real(4),parameter :: cp      = 1003.0       ! heat capacity (j/(kg*k))
  real(4),parameter :: gasr    = rkappa * cp  !
  real(4),parameter :: twomg   = 1.45849e-4   !
  real(4),parameter :: gamma   = 6.5e-3       ! lapse rate

  ! possible range of each variable ( e.g. used in check_range() )
!  real(4),parameter :: t_min     = 100.0      ! Temperature [K]
  real(4),parameter :: t_min     = 50.0      ! Temperature [K]
  real(4),parameter :: t_max     = 10000.0    ! Temperature [K]
  real(4),parameter :: pt_min    = 100.0      ! Potential Temperature [K]
  real(4),parameter :: pt_max    = 10000.0    ! Potential Temperature [K]
  real(4),parameter :: wind_min  = -500.0     ! Horizontal Wind [m/s]
  real(4),parameter :: wind_max  = 500.0      ! Horizontal Wind [m/s]
  real(4),parameter :: omega_min = -1.0e+10   ! Omega-velocity [Pa/s]
  real(4),parameter :: omega_max = 1.0e+10    ! Omega-velocity [Pa/s]
  real(4),parameter :: w_min     = -100.0     ! Vertical Wind [m/s]
  real(4),parameter :: w_max     = 100.0      ! Vertical Wind [m/s]
  real(4),parameter :: p_min     = 1.0e-10    ! Pressure [hPa]
  real(4),parameter :: p_max     = 1200.0     ! Pressure [hPa]
  real(4),parameter :: z_min     = -3000.0    ! Height [m] (-1000 is better)
  real(4),parameter :: z_max     = 1.0e+6     ! Height [m]
  real(4),parameter :: alt_min   = -1000.0    ! Altitude [m]
  real(4),parameter :: alt_max   = 10000.0    ! Altitude [m]
  real(4),parameter :: st_min    = -1.0e+12   ! Streamfunction [kg/s]
  real(4),parameter :: st_max    = 1.0e+12    ! Streamfunction [kg/s]
  real(4),parameter :: divf_min  = -0.1       ! EP Flux divergence [m/s^2]
  real(4),parameter :: divf_max  = 0.1        ! EP Flux divergence [m/s^2]
  real(4),parameter :: econv_min = -10.0      ! Energy Conversion [W/kg]
  real(4),parameter :: econv_max = 10.0       ! Energy Conversion [W/kg]

end module parameter
