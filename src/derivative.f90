!
! d(var)/dx
!
subroutine derivative_x( xdef, ydef, zdef, var, &
     &                   var_d )
  use parameter, only : pai
  implicit none
  integer,intent(in)  :: xdef, ydef, zdef
  real(4),intent(in)  :: var(xdef, ydef, zdef)
  real(4),intent(out) :: var_d(xdef, ydef, zdef)

  real(4) :: dlat
  integer :: x, y, z
  integer :: xl, xu
  
  dlat = 2 * pai / xdef

  do z=1, zdef
     do y=1, ydef
        do x=1, xdef

           xl = x - 1
           xu = x + 1
           if( x == 1 )    xl = xdef
           if( x == xdef ) xu = 1

           var_d(x,y,z) = ( var(xu,y,z) - var(xl,y,z) ) / ( 2.0 * dlat )

        end do
     end do
  end do

end subroutine derivative_x



!
! d(var)/dy
!
subroutine derivative_y( xdef, ydef, zdef, ylev, var, &
     &                   var_d )
  implicit none
  integer,intent(in)  :: xdef, ydef, zdef
  real(4),intent(in)  :: ylev(ydef)
  real(4),intent(in)  :: var(xdef, ydef, zdef)
  real(4),intent(out) :: var_d(xdef, ydef, zdef)

  integer :: x, y, z
  integer :: yl, yu

  do x=1, xdef
     do y=1, ydef
        do z=1, zdef

           yl = y - 1
           yu = y + 1
           if( y == 1 ) then 
              yl = 1
           else if( y == ydef ) then
              yu = ydef
           end if

           var_d(x,y,z) = ( var(x,yu,z) - var(x,yl,z) ) / &
                &         ( ylev(yu)    - ylev(yl) )

        end do
     end do
  end do

end subroutine derivative_y



!
! d(var)/dp
!
! Note: if the unit of plev is [hPa], then the result will be [*/hPa]
!
subroutine derivative_p( xdef, ydef, zdef, plev, ps, var, &
     &                   var_d )
  implicit none
  integer,intent(in)  :: xdef, ydef, zdef
  real(4),intent(in)  :: plev(zdef)
  real(4),intent(in)  :: ps(xdef, ydef)
  real(4),intent(in)  :: var(xdef, ydef, zdef)
  real(4),intent(out) :: var_d(xdef, ydef, zdef)

  integer :: x, y, z
  real(4) :: x_hlf(xdef, ydef, zdef+1), p_hlf(xdef, ydef, zdef+1), d

  do x=1, xdef
     do y=1, ydef
        do z=2, zdef
           
           if( plev(z) <= ps(x,y) ) then
              
              x_hlf(x,y,z) = ( var(x,y,z) - var(x,y,z-1) ) / &
                   &         ( plev(z) - plev(z-1) )
              p_hlf(x,y,z) = 0.5 * ( plev(z-1) + plev(z) )      
              
           else if( plev(z-1) < ps(x,y) ) then

              if( abs( ps(x,y) - plev(z-1)) / plev(z-1) < 0.01 ) then
                 x_hlf(x,y,z) = ( var(x,y,z-1) - var(x,y,z-2) ) / &
                      &         ( plev(z-1) - plev(z-2) )
              else
                 x_hlf(x,y,z) = ( var(x,y,z) - var(x,y,z-1) ) / &
                      &         ( ps(x,y) - plev(z-1) )
              end if

              p_hlf(x,y,z) = 0.5 * ( plev(z-1) + ps(x,y) )

           else
              
              x_hlf(x,y,z) = ( var(x,y,z) - var(x,y,z-1) ) / &
                   &         ( plev(z) - plev(z-1) )
              
              p_hlf(x,y,z) = 0.5 * ( plev(z-1) + plev(z) )     
              
           end if
           
        end do
        
        x_hlf(x,y,1) = x_hlf(x,y,2)
        x_hlf(x,y,zdef+1) = x_hlf(x,y,zdef)
        
        p_hlf(x,y,1) = max( 0.5*plev(1), -0.5*plev(2)+1.5*plev(1) )
        p_hlf(x,y,zdef+1) = -0.5 * plev(zdef-1) + 1.5 * plev(zdef)
        
        !intpl   
        do z=2, zdef-1
           d = ( plev(z) - p_hlf(x,y,z) ) / ( p_hlf(x,y,z+1) - p_hlf(x,y,z) )
           var_d(x,y,z) = x_hlf(x,y,z) * (1.0-d) + x_hlf(x,y,z+1) * d
        end do
        
        var_d(x,y,zdef) = x_hlf(x,y,zdef)
        var_d(x,y,1) = x_hlf(x,y,1)
        
     end do
  end do

  return
end subroutine derivative_p



!
! d(var)/dp
!
subroutine derivative_p_nops( xdef, ydef, zdef, plev, var, &
     &                        var_d )
  implicit none
  integer,intent(in)  :: xdef, ydef, zdef
  real(4),intent(in)  :: plev(zdef)
  real(4),intent(in)  :: var(xdef, ydef, zdef)
  real(4),intent(out) :: var_d(xdef, ydef, zdef)

  real(4) :: ps(xdef, ydef)

  ps(:,:) = 1e+10  ! dummy

  call derivative_p( xdef, ydef, zdef, plev, ps, var, &
       &             var_d )

end subroutine derivative_p_nops



!
! d(var)/dz
!
subroutine derivative_z( xdef, ydef, zdef, zlev, var, &
     &                   var_d )
  implicit none
  integer,intent(in)  :: xdef, ydef, zdef
  real(4),intent(in)  :: zlev(zdef)
  real(4),intent(in)  :: var(xdef, ydef, zdef)
  real(4),intent(out) :: var_d(xdef, ydef, zdef)

  integer :: x, y, z
  real(4) :: x_hlf(xdef, ydef, zdef+1), z_hlf(xdef, ydef, zdef+1), d

  do x=1, xdef
     do y=1, ydef

        do z=2, zdef
           x_hlf(x,y,z) = ( var(x,y,z) - var(x,y,z-1) ) / &
                &         ( zlev(z) - zlev(z-1) )
           z_hlf(x,y,z) = 0.5 * ( zlev(z-1) + zlev(z) )      
        end do
        
        x_hlf(x,y,1) = x_hlf(x,y,2)
        x_hlf(x,y,zdef+1) = x_hlf(x,y,zdef)
        
        z_hlf(x,y,1) = max( 0.5*zlev(1), -0.5*zlev(2)+1.5*zlev(1) )
        z_hlf(x,y,zdef+1) = -0.5 * zlev(zdef-1) + 1.5 * zlev(zdef)
        
        ! interpolate
        do z=2, zdef-1
           d = ( zlev(z) - z_hlf(x,y,z) ) / ( z_hlf(x,y,z+1) - z_hlf(x,y,z) )
           var_d(x,y,z) = x_hlf(x,y,z) * (1.0-d) + x_hlf(x,y,z+1) * d
        end do
        
        var_d(x,y,zdef) = x_hlf(x,y,zdef)
        var_d(x,y,1) = x_hlf(x,y,1)
        
     end do
  end do

  return
end subroutine derivative_z
