!############################################################
!
!  Common variables throughout the MIM program
!
!    Don't change variables out of com_var.f90
!
!############################################################
module com_var
  implicit none

  integer :: im      ! x (East -> West) direction grid
  integer :: jm      ! y (North -> South) direction grid
  integer :: km      ! input data z direction levels (p levels)
  integer :: ko      ! output data z direction levels (p+ levels, p++ levels)
  integer :: wmax    ! maximum wave number to analyze

  real(4),allocatable :: pin(:)     ! standard pressure levels [hPa]
  real(4),allocatable :: pout(:)    ! standard p+/p++ levels [hPa]
  real(4),allocatable :: zd(:)      ! z_dagger corresponding to pout [m]
  real(4),allocatable :: rho(:)     ! standard density corresponding to pout [kg/m^3]
  real(4),allocatable :: alat(:)    ! latitude [rad]
  real(4),allocatable :: sintbl(:)  ! sin(latitude)
  real(4),allocatable :: costbl(:)  ! cos(latitude)
  real(4),allocatable :: tantbl(:)  ! tan(latitude)


contains
  subroutine com_var_ini( im_in, jm_in, km_in, ko_in, wmax_in )
    integer,intent(in) :: im_in
    integer,intent(in) :: jm_in
    integer,intent(in) :: km_in
    integer,intent(in) :: ko_in
    integer,intent(in) :: wmax_in

    im = im_in
    jm = jm_in
    km = km_in
    ko = ko_in
    wmax = wmax_in
    
    allocate( pin(km) )
    allocate( pout(ko) )
    allocate( zd(ko) )
    allocate( rho(ko) )
    allocate( alat(jm) )
    allocate( sintbl(jm) )  
    allocate( costbl(jm) )
    allocate( tantbl(jm) )

  end subroutine com_var_ini


  subroutine com_var_end()
    deallocate( pin )
    deallocate( pout )
    deallocate( zd )
    deallocate( rho )
    deallocate( alat )
    deallocate( sintbl )  
    deallocate( costbl )
    deallocate( tantbl )
  end subroutine com_var_end


  !*********************************!
  !                                 !
  !          Set variables          !
  !                                 !
  !*********************************!

  ! pin [hPa]
  subroutine com_var_pin()
    use namelist, only : INPUT_ZDEF_LEVEL
    if( INPUT_ZDEF_LEVEL(1) < INPUT_ZDEF_LEVEL(2) ) then  ! Upper -> Lower
       pin(1:km) = INPUT_ZDEF_LEVEL(1:km)
    else                                                  ! Lower -> Upper
       pin(1:km) = INPUT_ZDEF_LEVEL(km:1:-1)
    end if
  end subroutine com_var_pin


  ! pout [hPa]
  subroutine com_var_pout()
    use namelist, only : OUTPUT_ZDEF_LEVEL
    if( OUTPUT_ZDEF_LEVEL(1) < OUTPUT_ZDEF_LEVEL(2) ) then  ! Upper -> Lower
       pout(1:ko) = OUTPUT_ZDEF_LEVEL(1:ko)
    else                                                  ! Lower -> Upper
       pout(1:ko) = OUTPUT_ZDEF_LEVEL(ko:1:-1)
    end if
  end subroutine com_var_pout


  ! zd [m]
  subroutine com_var_zd()
    use parameter, only : h0
    zd(:) = -h0 * log( pout(:) / 1000.0 )
  end subroutine com_var_zd


  ! rho [kg/m^3]
  subroutine com_var_rho()
    use parameter, only : h0, grav
    rho(:) = pout(:) * 100.0 / ( h0 * grav )
  end subroutine com_var_rho


  ! alat [rad]
  subroutine com_var_alat()
    use namelist, only : INPUT_YDEF_TYPE, INPUT_YDEF_LEVEL, &
         &               INPUT_YDEF_SOUTH, INPUT_YDEF_NORTH
    use parameter, only : pai
    real(4) :: rjm
    integer :: j
    
    ! - all variables should be YREV (north->south) 
    !   and ZREV (upper->lower) in mim.f90.
    if( INPUT_YDEF_TYPE == 'lat_radian' ) then

       if( INPUT_YDEF_LEVEL(1) < INPUT_YDEF_LEVEL(2) ) then
          alat(1:jm) = INPUT_YDEF_LEVEL(jm:1:-1)
       else
          alat(1:jm) = INPUT_YDEF_LEVEL(1:jm)
       end if
       
    else if( INPUT_YDEF_TYPE == 'lat_degree' ) then

       if( INPUT_YDEF_LEVEL(1) < INPUT_YDEF_LEVEL(2) ) then
          alat(1:jm) = INPUT_YDEF_LEVEL(jm:1:-1) * pai / 180.0
       else
          alat(1:jm) = INPUT_YDEF_LEVEL(1:jm) * pai / 180.0
       end if
       
    else if( INPUT_YDEF_TYPE == 'linear' ) then

       rjm = ( INPUT_YDEF_NORTH - INPUT_YDEF_SOUTH ) / ( jm - 1 )
       do j=1, jm
          alat(j) = ( INPUT_YDEF_NORTH - (j-1) * rjm ) * pai / 180
       end do

    end if

  end subroutine com_var_alat


  ! sine table
  subroutine com_var_sintbl()
    sintbl(:) = sin( alat(:) )
  end subroutine com_var_sintbl


  ! cosine table
  subroutine com_var_costbl()
    costbl(:) = cos( alat(:) )
  end subroutine com_var_costbl


  ! tangent table
  subroutine com_var_tantbl()
    tantbl(:) = tan( alat(:) )
  end subroutine com_var_tantbl


end module com_var
