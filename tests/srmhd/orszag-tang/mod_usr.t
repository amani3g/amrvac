module mod_usr
  use mod_global_parameters
  use mod_srmhd_parameters
  use mod_srmhd  

  implicit none
  double precision:: vmax, rho0, p0
contains

  !> Read this module parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ vmax, rho0, p0

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_variables

    call usr_params_read(par_files) 
    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian")
    call srmhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ixO^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ixO^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    logical,save     :: first=.true.
    double precision :: vmaxlocal

    if(first)then
      if(mype==0) then
          write(*,*)'Doing 2D ideal SRMHD, Orszag Tang problem'
          write(*,*)'rho - p - gamma:',rho0,p0,srmhd_gamma
      endif
      first=.false.
    endif

    vmaxlocal=vmax/dsqrt(two)
    w(ixO^S,rho_)=rho0
    w(ixO^S,mom(1))=-vmaxlocal*sin(x(ixO^S,2))
    w(ixO^S,mom(2))= vmaxlocal*sin(x(ixO^S,1))
    w(ixO^S,p_)=p0
    w(ixO^S,mag(1))=-sin(x(ixO^S,2))
    w(ixO^S,mag(2))= sin(two*x(ixO^S,1))
    if (srmhd_glm) w(ixO^S,psi_)     = 0.0d0
  
    call srmhd_get_4u_from_3v(ixG^L,ixO^L,w(ixG^S,mom(1):mom(ndir)))
    call srmhd_to_conserved(ixG^L,ixO^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
