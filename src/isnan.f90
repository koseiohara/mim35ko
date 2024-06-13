subroutine isnan(array)
    use, intrinsic :: IEEE_ARITHMETIC, only : IEEE_IS_NAN
    implicit none
    
    real(4), intent(in) :: array(:,:)

    if (any(IEEE_IS_NAN(array(:,:)))) then
        write(*,*)
        write(*,'(a)') 'WARNING --------------------------------------------'
        write(*,'(a)') '|   NaN is found'
        write(*,'(a)') '----------------------------------------------------'
        write(*,*)
    endif

end subroutine isnan

