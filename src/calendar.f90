!
! get timestep
!
subroutine calendar_tdef( tstep )
  use namelist, only : INPUT_TDEF_TYPE, &
       &               INPUT_TDEF_DAYNUM, &
       &               INPUT_TDEF_YEAR, INPUT_TDEF_MONTH, &
       &               INPUT_TDEF_365DAY, &
       &               INPUT_TDEF_TSTEP
  implicit none
  integer,intent(out) :: tstep
  integer :: m, tstep_temp

  if( INPUT_TDEF_TYPE == "tstep" ) then
     tstep = INPUT_TDEF_TSTEP
     
  else if( INPUT_TDEF_TYPE == "monthly" ) then
     call calendar( INPUT_TDEF_YEAR, INPUT_TDEF_MONTH, INPUT_TDEF_365DAY, tstep )
     tstep = tstep * INPUT_TDEF_DAYNUM
     
  else if( INPUT_TDEF_TYPE == "annual" ) then
     tstep = 0
     do m=1, 12
        call calendar( INPUT_TDEF_YEAR, m, INPUT_TDEF_365DAY, tstep_temp )
        tstep = tstep + tstep_temp
     end do
     tstep = tstep * INPUT_TDEF_DAYNUM
  end if

end subroutine calendar_tdef

!
! return days per month specified
! day365=0 : with leap year  day365=1 : without leap year
!
! Note: month <= 0 or month >= 13 is also available ( e.g. month=13 -> 1 )
!
subroutine calendar( year, month, day365, res )
  implicit none
  integer,intent(in) :: year, month, day365
  integer,intent(out) :: res
  integer :: month12

  ! days per month
  !                          1   2   3   4   5   6   7   8   9  10  11  12
  integer :: mdays(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

  month12 = month
  do while(month12 < 1)
     month12 = month12 + 12
  end do

  month12 = mod(month12-1, 12) + 1  
  res = mdays(month12)


  if( day365 == 1 .or. month12 /= 2 ) then
     return
  end if

  if( (year / 400) * 400 == year ) then       ! leap year
     res = res + 1

  else if( (year / 100) * 100 == year ) then  ! no leap year
     res = res

  else if( (year / 4) * 4 == year ) then      ! leap year
     res = res + 1
     
  end if
     
  return
end subroutine calendar
