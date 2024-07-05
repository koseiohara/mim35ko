!############################################################
!
!  Read/Write GrADS data
!
!  How to use
!    1. grads_open
!    2. grads_read / grads_write
!    3. grads_close
!
!############################################################
module grads
  
  use, intrinsic :: IEEE_ARITHMETIC, only : IEEE_IS_NAN
  use myglobal, only : nan_detector
  implicit none

  type grads_info
     integer :: ionum             ! file number
     integer :: imax, jmax, kmax  ! number of grid point
     integer :: irev, jrev, krev  ! 0: nothing, 1: reverse y-direction (YREV)
     integer :: record            ! location of record
     integer :: endian            ! 1: little  -1: big
  end type grads_info

contains

  !
  ! open GrADS file
  !
  subroutine grads_open( ionum, fname, imax, jmax, kmax, irev, jrev, krev, endian, ginfo )
    integer,intent(in)      :: ionum
    character(*),intent(in) :: fname ! file name
    integer,intent(in)      :: imax, jmax, kmax
    integer,intent(in)      :: irev, jrev, krev
    integer,intent(in)      :: endian
    type(grads_info),intent(out) :: ginfo

    ginfo%ionum  = ionum
    ginfo%imax   = imax
    ginfo%jmax   = jmax
    ginfo%kmax   = kmax
    ginfo%irev   = irev
    ginfo%jrev   = jrev
    ginfo%krev   = krev
    ginfo%endian = endian
    ginfo%record = 1

    open(ionum, file=fname, form='unformatted', access='direct', recl=imax*jmax*kmax*4) 

  end subroutine grads_open


  !
  ! close GrADS file
  !
  subroutine grads_close( ginfo )
    type(grads_info),intent(inout) :: ginfo
    close( ginfo%ionum )
    ginfo%ionum = -1
  end subroutine grads_close


  !
  ! read GrADS data
  !
  subroutine grads_read( ginfo, gdata )
    type(grads_info),intent(inout) :: ginfo
    real(4),intent(out) :: gdata(ginfo%imax, ginfo%jmax, ginfo%kmax)

    call grads_read_core( ginfo%ionum, ginfo%imax, ginfo%jmax, ginfo%kmax, &
         &                ginfo%irev, ginfo%jrev, ginfo%krev, ginfo%endian, &
         &                ginfo%record, gdata )
  end subroutine grads_read


  !
  ! write GrADS data
  !
  subroutine grads_write( icount, var, ginfo, gdata )
    integer, intent(in) :: icount
    character(*), intent(in) :: var
    type(grads_info),intent(inout) :: ginfo
    real(4),intent(in) :: gdata(ginfo%imax, ginfo%jmax, ginfo%kmax)

    call grads_write_core( ginfo%ionum, ginfo%imax, ginfo%jmax, ginfo%kmax, &
         &                 ginfo%irev, ginfo%jrev, ginfo%krev, ginfo%endian, &
         &                 ginfo%record, gdata )

    if (any(IEEE_IS_NAN(gdata))) then
        write(nan_detector,'(a,i0,a)') 'icount = ', icount, '  : NaN is found in ' // trim(var)
    endif

  end subroutine grads_write


  !
  ! read GrADS data
  !
  ! Note
  !   -This subroutine works independently, but grads_read() is easier to run.
  !
  subroutine grads_read_core( ionum, imax, jmax, kmax, &
       &                      irev, jrev, krev, endian, &
       &                      record, data )
    implicit none
    
    integer,intent(in) :: ionum            ! file number
    integer,intent(in) :: imax, jmax, kmax ! number of grid point
    integer,intent(in) :: irev, jrev, krev ! 0: nothing, 
                                           ! 1: reverse y-direction (YREV)
    integer,intent(in) :: endian           ! 1: little  -1: big
    integer,intent(inout) :: record        ! location of record
    real(4),intent(out) :: data(imax, jmax, kmax) ! data
    
    integer :: i, j, k
    integer :: istart, iend, istep      ! loop parameters
    integer :: jstart, jend, jstep      ! loop parameters
    integer :: kstart, kend, kstep      ! loop parameters
    
    ! set loop parameters (i)
    if( irev == 0 ) then
       istart = 1
       iend   = imax
       istep  = 1
    else if( irev == 1 ) then
       istart = imax
       iend   = 1
       istep  = -1
    else
       write(*,*) 'error in grads_read_core : invalid irev'
       stop 999
    end if
    
    ! set loop parameters (j)
    if( jrev == 0 ) then
       jstart = 1
       jend   = jmax
       jstep  = 1
    else if( jrev == 1 ) then
       jstart = jmax
       jend   = 1
       jstep  = -1
    else
       write(*,*) 'error in grads_read_core : invalid jrev'
       stop 999
    end if
    
    ! set loop parameters (k)
    if( krev == 0 ) then
       kstart = 1
       kend   = kmax
       kstep  = 1
    else if( krev == 1 ) then
       kstart = kmax
       kend   = 1
       kstep  = -1
    else
       write(*,*) 'error in grads_read_core : invalid krev'
       stop 999
    end if
        
    ! read
    read(ionum,rec=record) &
         &  (( (data(i,j,k), &
         &   i=istart,iend,istep), j=jstart,jend,jstep ), k=kstart,kend,kstep )
    
    !if( endian == -1 ) then
    !   do k=1, kmax
    !      do j=1, jmax
    !         do i=1, imax
    !            call endian_filter( data(i,j,k) )
    !         end do
    !      end do
    !   end do
    !end if
    
    record = record + 1

    return
  end subroutine grads_read_core
  
  
  !
  ! write GrADS data
  !
  ! Note
  !   -This subroutine works independently, but grads_write() is easier to run.
  !
  subroutine grads_write_core( ionum, imax, jmax, kmax, &
       &                       irev, jrev, krev, endian, &
       &                       record, data )
    implicit none
    
    integer,intent(in) :: ionum            ! file number
    integer,intent(in) :: imax, jmax, kmax ! number of grid point
    integer,intent(in) :: irev, jrev, krev ! 0: nothing, 
                                           ! 1: reverse y-direction (YREV)
    integer,intent(in) :: endian           ! 1: little  -1: big
    integer,intent(inout) :: record        ! location of record
    real(4),intent(in) :: data(imax, jmax, kmax) ! data

    integer :: i, j, k
    integer :: istart, iend, istep      ! loop parameters
    integer :: jstart, jend, jstep      ! loop parameters
    integer :: kstart, kend, kstep      ! loop parameters
    
    ! set loop parameters (i)
    if( irev == 0 ) then
       istart = 1
       iend   = imax
       istep  = 1
    else if( irev == 1 ) then
       istart = imax
       iend   = 1
       istep  = -1
    else
       write(*,*) 'error in grads_write_core : invalid irev'
       stop 999
    end if
    
    ! set loop parameters (j)
    if( jrev == 0 ) then
       jstart = 1
       jend   = jmax
       jstep  = 1
    else if( jrev == 1 ) then
       jstart = jmax
       jend   = 1
       jstep  = -1
    else
       write(*,*) 'error in grads_write_core : invalid jrev'
       stop 999
    end if
    
    ! set loop parameters (k)
    if( krev == 0 ) then
       kstart = 1
       kend   = kmax
       kstep  = 1
    else if( krev == 1 ) then
       kstart = kmax
       kend   = 1
       kstep  = -1
    else
       write(*,*) 'error in grads_write_core : invalid krev'
       stop 999
    end if
    
    ! write
    write(ionum,rec=record) &
         &  (( (data(i,j,k), &
         &   i=istart,iend,istep), j=jstart,jend,jstep ), k=kstart,kend,kstep )
    
    !if( endian == -1 ) then
    !   do k=1, kmax
    !      do j=1, jmax
    !         do i=1, imax
    !            call endian_filter( data(i,j,k) )
    !         end do
    !      end do
    !   end do
    !end if
    
    record = record + 1

    return
  end subroutine grads_write_core



    !***** Big endian <-> Little endian *****
    subroutine endian_filter(c)
      character,intent(inout) ::  c(4)
      c(1:4) = c(4:1:-1)
    end subroutine endian_filter
  
end module grads


!***** Big endian <-> Little endian *****
subroutine endian_filter(c)
  character,intent(inout) ::  c(4)
  c(1:4) = c(4:1:-1)
end subroutine endian_filter
