!#########################################################
!
!  Namelistlist  
!
!#########################################################
module namelist
  implicit none

  integer,parameter :: km_max=512, ko_max=512  ! maximum km & ko
  integer,parameter :: jm_max=2048             ! maximum jm

  character(1024) :: INPUT_TYPE
  character(1024) :: INPUT_UVT_FILENAME
  character(1024) :: INPUT_U_FILENAME
  character(1024) :: INPUT_V_FILENAME
  character(1024) :: INPUT_T_FILENAME
  character(1024) :: INPUT_PS_FILENAME
  character(1024) :: INPUT_MSL_FILENAME
  character(1024) :: INPUT_TS_FILENAME
  character(1024) :: INPUT_Z_FILENAME
  character(1024) :: INPUT_OMEGA_FILENAME
  character(1024) :: INPUT_TOPO_FILENAME
  character(1024) :: INPUT_Q_FILENAME

  character(1024) :: INPUT_UNIT_Z
  character(1024) :: INPUT_UNIT_PS
  character(1024) :: INPUT_UNIT_MSL
  character(1024) :: INPUT_UNIT_TOPO
  
  character(1024) :: INPUT_UNDEF_DEFAULT
  character(1024) :: INPUT_UNDEF_UVT
  character(1024) :: INPUT_UNDEF_U
  character(1024) :: INPUT_UNDEF_V
  character(1024) :: INPUT_UNDEF_T
  character(1024) :: INPUT_UNDEF_PS
  character(1024) :: INPUT_UNDEF_MSL
  character(1024) :: INPUT_UNDEF_TS
  character(1024) :: INPUT_UNDEF_Z
  character(1024) :: INPUT_UNDEF_OMEGA
  character(1024) :: INPUT_UNDEF_Q
  real(4)         :: INPUT_UNDEF_DEFAULT_REAL  ! below derivatives
  real(4)         :: INPUT_UNDEF_UVT_REAL
  real(4)         :: INPUT_UNDEF_U_REAL
  real(4)         :: INPUT_UNDEF_V_REAL
  real(4)         :: INPUT_UNDEF_T_REAL
  real(4)         :: INPUT_UNDEF_PS_REAL
  real(4)         :: INPUT_UNDEF_MSL_REAL
  real(4)         :: INPUT_UNDEF_TS_REAL
  real(4)         :: INPUT_UNDEF_Z_REAL
  real(4)         :: INPUT_UNDEF_OMEGA_REAL
  real(4)         :: INPUT_UNDEF_Q_REAL

  character(1024) :: INPUT_ENDIAN_DEFAULT
  character(1024) :: INPUT_ENDIAN_UVT
  character(1024) :: INPUT_ENDIAN_U
  character(1024) :: INPUT_ENDIAN_V
  character(1024) :: INPUT_ENDIAN_T
  character(1024) :: INPUT_ENDIAN_PS
  character(1024) :: INPUT_ENDIAN_MSL
  character(1024) :: INPUT_ENDIAN_TS
  character(1024) :: INPUT_ENDIAN_Z
  character(1024) :: INPUT_ENDIAN_OMEGA
  character(1024) :: INPUT_ENDIAN_Q
  character(1024) :: INPUT_ENDIAN_TOPO
  integer         :: INPUT_ENDIAN_DEFAULT_INT  ! below derivatives
  integer         :: INPUT_ENDIAN_UVT_INT
  integer         :: INPUT_ENDIAN_U_INT
  integer         :: INPUT_ENDIAN_V_INT
  integer         :: INPUT_ENDIAN_T_INT
  integer         :: INPUT_ENDIAN_PS_INT
  integer         :: INPUT_ENDIAN_MSL_INT
  integer         :: INPUT_ENDIAN_TS_INT
  integer         :: INPUT_ENDIAN_Z_INT
  integer         :: INPUT_ENDIAN_OMEGA_INT
  integer         :: INPUT_ENDIAN_Q_INT
  integer         :: INPUT_ENDIAN_TOPO_INT

  integer         :: INPUT_XDEF_NUM

  character(1024) :: INPUT_YDEF_TYPE
  integer         :: INPUT_YDEF_NUM
  real(4)         :: INPUT_YDEF_LEVEL(jm_max)
  real(4)         :: INPUT_YDEF_SOUTH
  real(4)         :: INPUT_YDEF_NORTH
  integer         :: INPUT_YDEF_YREV_DEFAULT
  integer         :: INPUT_YDEF_YREV_TOPO

  integer         :: INPUT_ZDEF_NUM, INPUT_ZDEF_NUM_OMEGA
  real(4)         :: INPUT_ZDEF_LEVEL(km_max)
  integer         :: INPUT_ZDEF_ZREV

  character(1024) :: INPUT_TDEF_TYPE
  integer         :: INPUT_TDEF_DAYNUM
  integer         :: INPUT_TDEF_TSTEP
  integer         :: INPUT_TDEF_YEAR, INPUT_TDEF_MONTH
  integer         :: INPUT_TDEF_365DAY
  real(4)         :: INPUT_TDEF_DT  ! derivative

  integer         :: WAVE_MAX_NUMBER

  integer         :: OUTPUT_ZDEF_NUM
  real(4)         :: OUTPUT_ZDEF_LEVEL(ko_max)

  character(1024) :: OUTPUT_ZONAL_FILENAME
  character(1024) :: OUTPUT_VINT_FILENAME
  character(1024) :: OUTPUT_GMEAN_FILENAME
  character(1024) :: OUTPUT_WAVE_FILENAME
  character(1024) :: OUTPUT_ERROR_FILENAME

contains
  !
  ! read Namelist from standard input
  !
  subroutine namelist_init()
    
    integer, parameter :: unit_input=100
    character(128), parameter :: fname_input='../nml/input.nml'
    integer, parameter :: unit_unit=101
    character(128), parameter :: fname_unit='../nml/input_unit.nml'
    integer, parameter :: unit_undef=102
    character(128), parameter :: fname_undef='../nml/input_undef.nml'
    integer, parameter :: unit_endian=103
    character(128), parameter :: fname_endian='../nml/input_endian.nml'
    integer, parameter :: unit_xdef=104
    character(128), parameter :: fname_xdef='../nml/input_xdef.nml'
    integer, parameter :: unit_ydef=105
    character(128), parameter :: fname_ydef='../nml/input_ydef.nml'
    integer, parameter :: unit_zdef=106
    character(128), parameter :: fname_zdef='../nml/input_zdef.nml'
    integer, parameter :: unit_tdef=107
    character(128), parameter :: fname_tdef='../nml/input_tdef.nml'
    integer, parameter :: unit_wave=108
    character(128), parameter :: fname_wave='../nml/wave.nml'
    integer, parameter :: unit_output=109
    character(128), parameter :: fname_output='../nml/output.nml'
    integer, parameter :: unit_output_zdef=110
    character(128), parameter :: fname_output_zdef='../nml/output_zdef.nml'

    integer :: is

    !***** declare *****!
    namelist / INPUT  / INPUT_TYPE, &
         &              INPUT_UVT_FILENAME, &
         &              INPUT_U_FILENAME, &
         &              INPUT_V_FILENAME, &
         &              INPUT_T_FILENAME, &
         &              INPUT_PS_FILENAME, &
         &              INPUT_MSL_FILENAME, &
         &              INPUT_TS_FILENAME, &
         &              INPUT_Z_FILENAME, &
         &              INPUT_OMEGA_FILENAME, &
         &              INPUT_TOPO_FILENAME, &
         &              INPUT_Q_FILENAME

    namelist / INPUT_UNIT / INPUT_UNIT_Z, INPUT_UNIT_PS, &
         &                  INPUT_UNIT_MSL, INPUT_UNIT_TOPO

    namelist / INPUT_UNDEF / INPUT_UNDEF_DEFAULT, &
         &                   INPUT_UNDEF_UVT, INPUT_UNDEF_U, &
         &                   INPUT_UNDEF_V, INPUT_UNDEF_T, &
         &                   INPUT_UNDEF_PS, INPUT_UNDEF_MSL, &
         &                   INPUT_UNDEF_TS, INPUT_UNDEF_Z, &
         &                   INPUT_UNDEF_OMEGA, INPUT_UNDEF_Q

    namelist / INPUT_ENDIAN / INPUT_ENDIAN_DEFAULT, &
         &                    INPUT_ENDIAN_UVT, INPUT_ENDIAN_U, &
         &                    INPUT_ENDIAN_V, INPUT_ENDIAN_T, &
         &                    INPUT_ENDIAN_PS, INPUT_ENDIAN_MSL, &
         &                    INPUT_ENDIAN_TS, INPUT_ENDIAN_Z, &
         &                    INPUT_ENDIAN_OMEGA, INPUT_ENDIAN_Q, &
         &                    INPUT_ENDIAN_TOPO


    namelist / INPUT_XDEF / INPUT_XDEF_NUM

    namelist / INPUT_YDEF / INPUT_YDEF_TYPE, INPUT_YDEF_NUM, &
         &                  INPUT_YDEF_LEVEL, &
         &                  INPUT_YDEF_SOUTH, INPUT_YDEF_NORTH, &
         &                  INPUT_YDEF_YREV_DEFAULT, INPUT_YDEF_YREV_TOPO
    namelist / INPUT_ZDEF / INPUT_ZDEF_NUM, INPUT_ZDEF_NUM_OMEGA, &
         &                  INPUT_ZDEF_LEVEL, &
         &                  INPUT_ZDEF_ZREV

    namelist / INPUT_TDEF / INPUT_TDEF_TYPE, &
         &                  INPUT_TDEF_DAYNUM, &
         &                  INPUT_TDEF_TSTEP, &
         &                  INPUT_TDEF_YEAR, INPUT_TDEF_MONTH, &
         &                  INPUT_TDEF_365DAY

    namelist / WAVE / WAVE_MAX_NUMBER

    namelist / OUTPUT / OUTPUT_ZONAL_FILENAME, &
         &              OUTPUT_VINT_FILENAME, &
         &              OUTPUT_GMEAN_FILENAME, &
         &              OUTPUT_WAVE_FILENAME, &
         &              OUTPUT_ERROR_FILENAME

    namelist / OUTPUT_ZDEF / OUTPUT_ZDEF_NUM, OUTPUT_ZDEF_LEVEL


    !***** default values *****!
    INPUT_UVT_FILENAME   = ''
    INPUT_PS_FILENAME    = ''
    INPUT_MSL_FILENAME   = ''
    INPUT_TS_FILENAME    = ''
    INPUT_OMEGA_FILENAME = ''
    INPUT_Q_FILENAME     = ''

    INPUT_UNIT_Z    = 'm'
    INPUT_UNIT_PS   = 'hPa'
    INPUT_UNIT_MSL  = 'hPa'
    INPUT_UNIT_TOPO = 'm'

    INPUT_UNDEF_DEFAULT = '9.999e+20'
    INPUT_UNDEF_UVT     = ''
    INPUT_UNDEF_U       = ''
    INPUT_UNDEF_V       = ''
    INPUT_UNDEF_T       = ''
    INPUT_UNDEF_PS      = ''
    INPUT_UNDEF_MSL     = ''
    INPUT_UNDEF_TS      = ''
    INPUT_UNDEF_Z       = ''
    INPUT_UNDEF_OMEGA   = ''
    INPUT_UNDEF_Q       = ''

    INPUT_ENDIAN_DEFAULT = 'little'
    INPUT_ENDIAN_UVT     = ''
    INPUT_ENDIAN_U       = ''
    INPUT_ENDIAN_V       = ''
    INPUT_ENDIAN_T       = ''
    INPUT_ENDIAN_PS      = ''
    INPUT_ENDIAN_MSL     = ''
    INPUT_ENDIAN_TS      = ''
    INPUT_ENDIAN_Z       = ''
    INPUT_ENDIAN_OMEGA   = ''
    INPUT_ENDIAN_Q       = ''
    INPUT_ENDIAN_TOPO    = ''

    INPUT_YDEF_SOUTH        = -90.0
    INPUT_YDEF_NORTH        = 90.0
    INPUT_YDEF_YREV_DEFAULT = 0
    INPUT_YDEF_YREV_TOPO    = -1

    INPUT_ZDEF_NUM_OMEGA = 0
    INPUT_ZDEF_ZREV      = 0

    INPUT_TDEF_365DAY = 0

    WAVE_MAX_NUMBER = 0

    OUTPUT_ZONAL_FILENAME  = ''
    OUTPUT_VINT_FILENAME   = ''
    OUTPUT_GMEAN_FILENAME  = ''
    OUTPUT_WAVE_FILENAME   = ''


    open(unit_input, file=fname_input, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_input)
        error stop
    endif

    open(unit_unit, file=fname_unit, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_unit)
        error stop
    endif

    open(unit_undef, file=fname_undef, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_undef)
        error stop
    endif

    open(unit_endian, file=fname_endian, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_endian)
        error stop
    endif

    open(unit_xdef, file=fname_xdef, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_xdef)
        error stop
    endif

    open(unit_ydef, file=fname_ydef, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_ydef)
        error stop
    endif

    open(unit_zdef, file=fname_zdef, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_zdef)
        error stop
    endif

    open(unit_tdef, file=fname_tdef, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_tdef)
        error stop
    endif

    open(unit_wave, file=fname_wave, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_wave)
        error stop
    endif

    open(unit_output, file=fname_output, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_output)
        error stop
    endif

    open(unit_output_zdef, file=fname_output_zdef, action='read', iostat=is)
    if (is /= 0) then
        write(*,*) 'Failed to open ' // trim(fname_output_zdef)
        error stop
    endif


    !***** read *****!
    read(unit_input, nml=INPUT)
    read(unit_unit, nml=INPUT_UNIT)
    read(unit_undef, nml=INPUT_UNDEF)
    read(unit_endian, nml=INPUT_ENDIAN)
    read(unit_xdef, nml=INPUT_XDEF)
    read(unit_ydef, nml=INPUT_YDEF)
    read(unit_zdef, nml=INPUT_ZDEF)
    read(unit_tdef, nml=INPUT_TDEF)
    read(unit_wave, nml=WAVE)
    read(unit_output, nml=OUTPUT)
    read(unit_output_zdef, nml=OUTPUT_ZDEF)


    close(unit_input)
    close(unit_unit)
    close(unit_undef)
    close(unit_endian)
    close(unit_xdef)
    close(unit_ydef)
    close(unit_zdef)
    close(unit_tdef)
    close(unit_wave)
    close(unit_output)
    close(unit_output_zdef)


    !***** adjust *****!
    read(INPUT_UNDEF_DEFAULT,*) INPUT_UNDEF_DEFAULT_REAL

    call undef( INPUT_UNDEF_UVT, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_UVT_REAL )
    call undef( INPUT_UNDEF_U, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_U_REAL )
    call undef( INPUT_UNDEF_V, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_V_REAL )
    call undef( INPUT_UNDEF_T, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_T_REAL )
    call undef( INPUT_UNDEF_PS, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_PS_REAL )
    call undef( INPUT_UNDEF_MSL, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_MSL_REAL )
    call undef( INPUT_UNDEF_TS, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_TS_REAL )
    call undef( INPUT_UNDEF_Z, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_Z_REAL )
    call undef( INPUT_UNDEF_OMEGA, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_OMEGA_REAL )
    call undef( INPUT_UNDEF_Q, INPUT_UNDEF_DEFAULT_REAL, &
         &      INPUT_UNDEF_Q_REAL )


    if( INPUT_ENDIAN_DEFAULT == 'little' ) then
       INPUT_ENDIAN_DEFAULT_INT = 1
    else if( INPUT_ENDIAN_DEFAULT == 'big' ) then
       INPUT_ENDIAN_DEFAULT_INT = -1
    else
       write(0,*) 'error in namelist_init() : INPUT_ENDIAN_DEFAULT = ' &
            &     // INPUT_ENDIAN_DEFAULT
    end if

    call endian( INPUT_ENDIAN_UVT, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_UVT_INT )
    call endian( INPUT_ENDIAN_U, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_U_INT )
    call endian( INPUT_ENDIAN_V, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_V_INT )
    call endian( INPUT_ENDIAN_T, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_T_INT )
    call endian( INPUT_ENDIAN_PS, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_PS_INT )
    call endian( INPUT_ENDIAN_MSL, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_MSL_INT )
    call endian( INPUT_ENDIAN_TS, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_TS_INT )
    call endian( INPUT_ENDIAN_Z, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_Z_INT )
    call endian( INPUT_ENDIAN_OMEGA, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_OMEGA_INT )
    call endian( INPUT_ENDIAN_Q, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_Q_INT )
    call endian( INPUT_ENDIAN_TOPO, INPUT_ENDIAN_DEFAULT_INT, &
         &       INPUT_ENDIAN_TOPO_INT )

    if( INPUT_YDEF_YREV_TOPO == -1 ) &
         &  INPUT_YDEF_YREV_TOPO = INPUT_YDEF_YREV_DEFAULT

    if( INPUT_ZDEF_NUM_OMEGA == 0 ) then
       INPUT_ZDEF_NUM_OMEGA = INPUT_ZDEF_NUM
    end if
    
    INPUT_TDEF_DT = 24.0 * 60.0 * 60.0 / INPUT_TDEF_DAYNUM


    !***** check *****!
    call namelist_check()
    
  end subroutine namelist_init


  subroutine undef( undef_char, undef_default, undef_out )
    character(*),intent(in) :: undef_char
    real(4),intent(in)      :: undef_default
    real(4),intent(out)     :: undef_out

    if( undef_char == '' ) then
       undef_out = undef_default
    else
       read(undef_char,*) undef_out
    end if
  end subroutine undef


  subroutine endian( endian_char, endian_default, endian_out )
    character(*),intent(in) :: endian_char
    integer,intent(in)      :: endian_default
    integer,intent(out)     :: endian_out

    if( endian_char == 'little' ) then
       endian_out = 1
    else if( endian_char == 'big' ) then
       endian_out = -1
    else
       endian_out = endian_default
    end if
  end subroutine endian


  !
  ! check namelist
  !
  subroutine namelist_check()

    if( INPUT_UNIT_Z /= 'm' .and. INPUT_UNIT_Z /= 'm^2/s^2' ) then
       write(*,*) 'error: INPUT_UNIT_Z=' // TRIM(INPUT_UNIT_Z)
       stop
    end if

    if( INPUT_UNIT_PS /= 'Pa' .and. INPUT_UNIT_PS /= 'hPa' ) then
       write(*,*) 'error: INPUT_UNIT_PS=' // TRIM(INPUT_UNIT_PS)
       stop
    end if

    if( INPUT_UNIT_MSL /= 'Pa' .and. INPUT_UNIT_MSL /= 'hPa' ) then
       write(*,*) 'error: INPUT_UNIT_MSK=' // TRIM(INPUT_UNIT_MSL)
       stop
    end if

    if( INPUT_UNIT_TOPO /= 'm' .and. INPUT_UNIT_TOPO /= 'm^2/s^2' ) then
       write(*,*) 'error: INPUT_UNIT_TOPO=' // TRIM(INPUT_UNIT_TOPO)
       stop
    end if

    if( INPUT_XDEF_NUM < 1 ) then
       write(*,*) 'error: INPUT_XDEF_NUM=', INPUT_XDEF_NUM
       stop
    end if

    if(    INPUT_YDEF_TYPE /= 'lat_degree'   .and. &
         & INPUT_YDEF_TYPE /= 'lat_radian' .and. &
         & INPUT_YDEF_TYPE /= 'linear' ) then
       write(*,*) 'error: INPUT_YDEF_TYPE=' // TRIM(INPUT_YDEF_TYPE)
       stop
    end if

    if( INPUT_YDEF_NUM < 1 .or. INPUT_YDEF_NUM > jm_max ) then
       write(*,*) 'error: INPUT_YDEF_NUM=', INPUT_YDEF_NUM
       stop
    end if

    if( INPUT_ZDEF_NUM < 1 .or. INPUT_ZDEF_NUM > km_max ) then
       write(*,*) 'error: INPUT_ZDEF_NUM=', INPUT_ZDEF_NUM
       stop
    end if

    if( INPUT_ZDEF_NUM_OMEGA > INPUT_ZDEF_NUM ) then
       write(*,*) 'error: INPUT_ZDEF_NUM_OMEGA(', INPUT_ZDEF_NUM_OMEGA, &
            &     ') > INPUT_ZDEF_NUM(', INPUT_ZDEF_NUM, ')'
    end if

    if(    INPUT_TDEF_TYPE /= 'tstep'   .and. &
         & INPUT_TDEF_TYPE /= 'monthly' .and. &
         & INPUT_TDEF_TYPE /= 'annual' ) then
       write(*,*) 'error: INPUT_TDEF_TYPE=', INPUT_TDEF_TYPE
       stop
    end if

    if( OUTPUT_ZDEF_NUM < 1 .or. OUTPUT_ZDEF_NUM > km_max ) then
       write(*,*) 'error: OUTPUT_ZDEF_NUM=', OUTPUT_ZDEF_NUM
       stop
    end if

    if( OUTPUT_WAVE_FILENAME /= '' ) then
       if( WAVE_MAX_NUMBER > INPUT_XDEF_NUM ) then
          write(0,*) 'error: WAVE_MAX_NUMBER(', WAVE_MAX_NUMBER, &
               &     ') > INPUT_XDEF_NUM(', INPUT_XDEF_NUM, ')'
          stop
       else if( WAVE_MAX_NUMBER < 1 ) then
          write(*,*) 'warning: OUTPUT_WAVE_FILENAME=' &
               &  // TRIM(OUTPUT_WAVE_FILENAME) // ' is invalid'
       end if
    end if

  end subroutine namelist_check
end module namelist

