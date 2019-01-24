module mod_usr
  use mod_global_parameters
  use mod_srmhd_parameters
  use mod_srmhd  

  implicit none
contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_set_B0 => specialset_B0

    call set_coordinate_system("Cartesian")
    call srmhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ixO^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ixO^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: Rblast,Rhoblast,Rhoout,Pblast,Pout
    logical,save     :: first=.true.

   if(mype==0)then
      if (first) then
         print*,'A simple SRMHD blastwave problem (see Mignone et al. 2007)'
	 print*,'Setup for 2D Cartesian geometry'
      end if
   end if

   Rblast = 0.1d0
   Rhoblast = 1.0d-2
   Rhoout = 1.0d-4
   Pblast = 1.0
   Pout = 3.0d-5
   where(x(ixO^S,1)<=Rblast)
      w(ixO^S,rho_)   = Rhoblast
      w(ixO^S,mom(1)) = 0.0d0
      w(ixO^S,mom(2)) = 0.0d0
      w(ixO^S,p_)     = Pblast
   elsewhere
      w(ixO^S,rho_)   = Rhoout
      w(ixO^S,mom(1)) = 0.0d0
      w(ixO^S,mom(2)) = 0.0d0
      w(ixO^S,p_)     = Pout
   endwhere

      w(ixO^S,mag(1)) = 0.0d0
      w(ixO^S,mag(2)) = 0.0d0

    call srmhd_get_4u_from_3v(ixG^L,ixO^L,w(ixG^S,mom(1):mom(ndir)))
    call srmhd_to_conserved(ixG^L,ixO^L,w,x)

    if (first) then
       if(mype==0)then
          print*,'Rblast=',Rblast
          print*,'Rhoblast=',Rhoblast
 	  print*,'Rhoout=',Rhoout
       end if
       first=.false.
    end if


  end subroutine initonegrid_usr

end module mod_usr
