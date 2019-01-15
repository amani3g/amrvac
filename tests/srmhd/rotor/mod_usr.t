module mod_usr
  use mod_global_parameters
  use mod_srmhd_parameters
  use mod_srmhd  

  implicit none
contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian")
    call srmhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ixO^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ixO^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: InitRadius, Middle^D, omega
    logical,save     :: first=.true.

   if(mype==0)then
      if (first) then
         print*,'Setting up initial conditions for 2D rotor problem'
      end if
   end if

   {Middle^D = xprobmin^D + half*(xprobmax^D-xprobmin^D)\}
   InitRadius = 0.1d0
   omega = 9.95d0
   where({^D&(x(ixO^S,^D)-Middle^D)**2+}<=InitRadius**2)
      w(ixO^S,rho_)   = 1.0d1
      w(ixO^S,mom(1)) = -omega*(x(ixO^S,2)-Middle2)
      w(ixO^S,mom(2)) =  omega*(x(ixO^S,1)-Middle1)
      w(ixO^S,p_)     = 1.0d0
      w(ixO^S,mag(1)) = 1.0d0
      w(ixO^S,mag(2)) = 0.0d0
   elsewhere
      w(ixO^S,rho_)   = 1.0d0
      w(ixO^S,mom(1)) = 0.0d0
      w(ixO^S,mom(2)) = 0.0d0
      w(ixO^S,p_)     = 1.0d0
      w(ixO^S,mag(1)) = 1.0d0
      w(ixO^S,mag(2)) = 0.0d0
   endwhere

    call srmhd_get_4u_from_3v(ixG^L,ixO^L,w(ixG^S,mom(1):mom(ndir)))
    call srmhd_to_conserved(ixG^L,ixO^L,w,x)

    if (first) then
       if(mype==0)then
          print*,'Rotor problem with omega=',omega
          print *,'gamma=',srmhd_gamma
          print*,'initial rotor radius=',InitRadius
       end if
       first=.false.
    end if


  end subroutine initonegrid_usr

end module mod_usr
