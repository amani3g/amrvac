!> Special Relativistic Hydrodynamics module
! Adapted for RHD, D. Millas January 2019 (original SRMHD by Z. Meliani, modified by R. Keppens)
module mod_srhd_phys
  use mod_global_parameters, only: std_len
  use mod_srhd_parameters
  use mod_srhd_eos
  
  implicit none
  private

  ! Public methods
  public :: srhd_phys_init
  public :: srhd_kin_en_primitive
  public :: srhd_get_pthermal
  public :: srhd_get_p_total
  public :: srhd_get_v
  public :: srhd_get_v_idim
  public :: srhd_to_conserved
  public :: srhd_to_primitive
  public :: srhd_get_csound2
  public :: srhd_get_csound_prim
  public :: srhd_get_4u_from_3v
contains

  !> Read this module's parameters from a file
  subroutine srhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/ srhd_energy, srhd_eos, srhd_n_tracer, srhd_gamma, &
                          srhd_adiab, He_abundance, SI_unit, srhd_particles, &
                          srhd_glm_alpha, srhd_4th_order, &
                          typedivbfix, source_split_divb, &
                          divbdiff, divbwave, srhd_glm,&
                          boundary_divbfix, boundary_divbfix_skip, &
                          srhd_maxiterationNR, srhd_absaccNR, srhd_tolerNR,&
                          srhd_checkNR, srhd_maxdspeed, small_vec2

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine srhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine srhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = srhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine srhd_write_info

  subroutine srhd_phys_init()
    use mod_global_parameters
    use mod_particles, only: particles_init
    use mod_physics

    integer :: itr, idir


    call srhd_read_params(par_files)

    physics_type = "srhd"
    phys_energy=srhd_energy

    ! set default gamma for polytropic/isothermal process
    if(.not.srhd_energy) then
       srhd_gamma=1.d0
    else
       ! set derived quantities
       gamma_1=srhd_gamma-1.0d0
       inv_gamma_1=1.d0/gamma_1
       gamma_to_gamma_1=srhd_gamma/gamma_1
    endif
    
    use_particles=srhd_particles
    srhd_maxspeed=1.0d0-srhd_maxdspeed
  
    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (srhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(tracer(srhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! Set index for auxiliary variables
    xi_  = var_set_auxvar('xi')
    lfac_= var_set_auxvar('lfac')
    ! auxiliary variables need boundary conditions
    nwfluxbc=nwfluxbc+2

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    do idir=1,ndir
      if(ndim>1) flux_type(idir,mag(idir))=flux_tvdlf
    end do
 
    phys_get_dt              => srhd_get_dt
    phys_get_cmax            => srhd_get_cmax
    phys_get_cbounds         => srhd_get_cbounds
    phys_get_flux            => srhd_get_flux
    phys_get_v_idim          => srhd_get_v_idim
    phys_add_source_geom     => srhd_add_source_geom
    phys_add_source          => srhd_add_source
    phys_to_conserved        => srhd_to_conserved
    phys_to_primitive        => srhd_to_primitive
    phys_get_aux             => srhd_get_auxiliary
    phys_get_aux_prim        => srhd_get_auxiliary_prim
    phys_check_params        => srhd_check_params
    phys_check_w             => srhd_check_w
    phys_get_pthermal        => srhd_get_pthermal
    phys_boundary_adjust     => srhd_boundary_adjust
    phys_write_info          => srhd_write_info
    phys_handle_small_values => srhd_handle_small_values

    ! derive units from basic units
    call srhd_physical_units()

    ! Initialize particles module
    if(srhd_particles) then
      call particles_init()
      phys_req_diagonal = .true.
    end if

  end subroutine srhd_phys_init

  subroutine srhd_check_params
    use mod_global_parameters

    ! after user parameter setting
    if (.not. srhd_energy) then
       if (srhd_gamma <= 0.0d0) call mpistop ("Error: srhd_gamma <= 0")
       if (srhd_adiab <  0.0d0) call mpistop ("Error: srhd_adiab < 0")
       small_pressure = srhd_adiab*small_density**srhd_gamma
    else
       if (srhd_gamma <= 0.0d0 .or. srhd_gamma == 1.0d0) &
            call mpistop ("Error: srhd_gamma <= 0 or srhd_gamma == 1")
       small_e = small_pressure * inv_gamma_1
    end if
    small_xi=small_density+gamma_to_gamma_1*small_pressure

  end subroutine srhd_check_params

  subroutine srhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0

    ! note const_c is in cgs 
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      unit_velocity=const_c/100.0d0
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
      unit_velocity=const_c
    end if
    unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
    unit_time=unit_length/unit_velocity

  end subroutine srhd_physical_units

  subroutine srhd_check_w(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in)            :: primitive
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: w(ixI^S,nw)
    integer, intent(inout)         :: flag(ixI^S)
    double precision :: tmp(ixI^S)

    flag(ixO^S)=0
    if (primitive) then
       where(w(ixO^S, rho_) < small_density) flag(ixO^S) = rho_
    else
       where(w(ixO^S, d_)/w(ixO^S,lfac_) < small_density) flag(ixO^S) = rho_
    endif

    if (srhd_energy) then
       if (block%e_is_internal) then
          where(w(ixO^S, e_) < small_e) flag(ixO^S) = e_
       else
          if (primitive)then
            where(w(ixO^S, p_) < small_pressure) flag(ixO^S) = p_ ! p_=e_
          else 
            where(w(ixO^S, e_) < small_e) flag(ixO^S) = e_
         end if 
       end if 
    end if 

  end subroutine srhd_check_w

  !> Set auxiliary variables from a primitive state
  subroutine srhd_get_auxiliary_prim(ixI^L,ixO^L,w)
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(inout)    :: w(ixI^S, nw)
    double precision, dimension(ixO^S) :: sqrU,B2,VdotB,rhoh,sqrV

    sqrU(ixO^S)    = sum(w(ixO^S, mom(:))**2, dim=ndim+1)
    w(ixO^S,lfac_) = sqrt(1.0d0+sqrU(ixO^S))
    sqrV           = sqrU/w(ixO^S,lfac_)

    ! fill the auxiliary variable xi and density D
    call srhd_get_enthalpy(ixO^L,w(ixO^S,rho_),w(ixO^S,p_),rhoh)

    ! with enthalpy w: xi= lfac^2 rhoh
    w(ixO^S,xi_) = w(ixO^S,lfac_)**2.0d0*rhoh(ixO^S)

  end subroutine srhd_get_auxiliary_prim

 !> Calculate thermal pressure for enthalpy and density
  subroutine srhd_get_pressure_fromprimitive(ixI^L,ixO^L,w,pth)
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: w(ixI^S,nw)
    double precision, intent(out)  :: pth(ixI^S)

    double precision               :: rhoh(ixO^S)

    rhoh  = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_pressure_primitive_eos(ixI^L,ixO^L,w(ixO^S,rho_),rhoh,pth) 
  end subroutine srhd_get_pressure_fromprimitive

  !> Transform primitive variables into conservative ones
  subroutine srhd_to_conserved(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(inout)    :: w(ixI^S, nw)
    double precision, intent(in)       :: x(ixI^S, 1:ndim)
    integer                            :: idir, itr
    double precision, dimension(ixO^S) :: sqrU,B2,VdotB,rhoh,sqrV
    integer, dimension(ixO^S)          :: flag_error    
    character(len=30)                  :: subname_loc
 
    flag_error(ixO^S)=0
    where(w(ixO^S,rho_)<small_density)flag_error(ixO^S)=rho_

    sqrU(ixO^S)    = sum(w(ixO^S, mom(:))**2, dim=ndim+1)
    w(ixO^S,lfac_) = dsqrt(1.0d0+sqrU(ixO^S))
    sqrV=sqrU/w(ixO^S,lfac_)

    ! fill the auxiliary variable xi and density D
    call srhd_get_enthalpy(ixO^L,w(ixO^S,rho_),w(ixO^S,p_),rhoh)

    ! with enthalpy w: xi= lfac^2 rhoh
    w(ixO^S,xi_) = w(ixO^S,lfac_)**2.0d0*rhoh(ixO^S)

    ! density: d = lfac * rho
    w(ixO^S,rho_)=w(ixO^S,lfac_)*w(ixO^S,rho_) 
    
    !!! DM this must change !!!
    call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.false.,B2=B2,VdotB=VdotB)
    ! Convert velocity to momentum
    ! s= (xi + B^2) * v - (v.B) * B
    ! re-use rhoh as rhoh=(xi+B^2)/lfac
    rhoh(ixO^S)= (w(ixO^S,xi_)+B2(ixO^S))/w(ixO^S,lfac_)
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = rhoh(ixO^S) * w(ixO^S,mom(idir))&
                             -VdotB(ixO^S)*w(ixO^S,mag(idir))
    end do 

    if (srhd_energy) then
       ! re-use sqrU=v^2 B^2 - (v.B)^2
       sqrU(ixO^S) = B2(ixO^S)*sqrV(ixO^S)-VdotB(ixO^S)**2.0d0
       ! sqrU should positive
       where(sqrU(ixO^S)<0.0d0) flag_error(ixO^S)=e_
       ! re-use sqrV=xi - p -D
       sqrV(ixO^S) = w(ixO^S,xi_) - w(ixO^S,p_) - w(ixO^S,d_)
       where(sqrV(ixO^S)<0.0d0) flag_error(ixO^S)=e_
       ! E = xi - p +(B^2+v^2 B^2 - (v.B)^2)/2- D
       w(ixO^S,e_)=sqrV(ixO^S) +0.5d0*(B2(ixO^S) + sqrU(ixO^S))
       if(type_divb==divb_glm2) w(ixO^S,e_)=w(ixO^S,e_) &
                                  + 0.5d0*w(ixO^S,psi_)**2
    end if 
    subname_loc= 'srhd_to_conserved' 
    if (check_small_values) call srhd_handle_small_values(.false., &
                                  w, x, ixI^L, ixO^L,trim(subname_loc),&
                                  flag_error=flag_error)
  end subroutine srhd_to_conserved

    !!! DM Same as above

  !> Transform conservative variables into primitive ones
  subroutine srhd_to_primitive(ixI^L,ixO^L,w,x)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: plocal(ixI^S),B2(ixO^S),VdotB(ixO^S)
    integer                         :: itr, idir
    character(len=30)               :: subname_loc
 
    subname_loc='srhd_to_primitive'
    ! get auxiliary variables
    call srhd_get_auxiliary(.true.,w,x,ixI^L,ixO^L,subname_loc)

    w(ixO^S,rho_) = w(ixO^S,d_)/w(ixO^S,lfac_)
    if (srhd_energy) then
      if(.not.block%e_is_internal) then
          call srhd_get_pressure_fromprimitive(ixI^L,ixI^L,w,plocal)
          w(ixO^S,p_)   = plocal(ixO^S)
      else
          w(ixO^S,p_)   = gamma_1*w(ixO^S, e_)
      end if
    end if

    call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2,VdotB)
    do idir=1,ndir
       w(ixO^S,mom(idir)) = w(ixO^S,lfac_)*(w(ixO^S,mom(idir))&
                      +VdotB*w(ixO^S,mag(idir)))&
                      /(w(ixO^S,xi_)+B2)
    end do 

    if (check_small_values) call srhd_handle_small_values(.true., w, x&
                                 , ixI^L, ixO^L,trim(subname_loc))
  end subroutine srhd_to_primitive

  subroutine srhd_handle_small_values(primitive, w, x, ixI^L, ixO^L, &
                                       subname,flag_error)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    integer, optional, intent(in)   :: flag_error(ixO^S)               

    double precision :: smallone
    integer :: idir, flag(ixI^S)

    if (small_values_method == "ignore") return
    if(present(flag_error)) then
     flag(ixO^S) = flag_error(ixO^S)
    else
     call srhd_check_w(primitive, ixI^L, ixO^L, w, flag)
    end if
    if (any(flag(ixO^S) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixO^S) /= 0) w(ixO^S,rho_) = small_density

          do idir = 1, ndir
             where(flag(ixO^S) /= 0) w(ixO^S, mom(idir)) = 0.0d0
          end do

          if (srhd_energy) then
             if(primitive) then
               smallone = small_pressure
             else
               smallone = small_e
             end if
             where(flag(ixO^S) /= 0) w(ixO^S,e_) = smallone
          end if
       case ("average")
          call small_values_average(ixI^L, ixO^L, w, x, flag)
       case default
          call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       end select
    end if
  end subroutine srhd_handle_small_values

  !> Calculate thermal pressure within ixO^L
  subroutine srhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    implicit none

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)

    double precision             :: rho(ixO^S),rhoh(ixO^S)

    rho  = w(ixO^S,d_)/w(ixO^S,lfac_)
    rhoh =w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_pthermal_eos(ixI^L,ixO^L,x,rho,rhoh,w(ixO^S,e_),pth)
  end subroutine srhd_get_pthermal

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  subroutine srhd_get_csound2_prim(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    implicit none
    integer, intent(in)               :: ixI^L, ixO^L
    double precision, intent(in)      :: w(ixI^S,nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    double precision, intent(out)     :: csound2(ixO^S)

    double precision,dimension(ixO^S) :: rhoh

    
    rhoh=w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_csound2_prim_eos(ixI^L,ixO^L,x,w(ixO^S,rho_),&
                                    rhoh,w(ixO^S,p_),csound2)
  end subroutine srhd_get_csound2_prim

  !> Convert energy to entropy
  subroutine e_to_rhos(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision,intent(inout)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

!!! DM fix the energy (delete the magnetic)
    if (srhd_energy) then
      if(.not.block%e_is_internal) &
        w(ixO^S, e_) = w(ixO^S, e_) - srhd_kin_en_primitive(w, ixI^L, ixO^L) &
                                    - srhd_mag_en_primitive(w, ixI^L, ixO^L)
      w(ixO^S, e_) = gamma_1* w(ixO^S, rho_)**(1.0d0 - srhd_gamma) * &
            w(ixO^S, e_)
    else
      call mpistop("e_to_rhos can not be used without energy equation!")
    end if
  end subroutine e_to_rhos

  !> Convert entropy to energy
  subroutine rhos_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (srhd_energy) then
       w(ixO^S, e_) = w(ixO^S, rho_)**gamma_1 * w(ixO^S, e_) &
            * inv_gamma_1
       if(.not.block%e_is_internal) &
         w(ixO^S, e_) =w(ixO^S, e_) + srhd_kin_en_primitive(w, ixI^L, ixO^L) + &
            srhd_mag_en_primitive(w, ixI^L, ixO^L)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

!!! DM Modify...
  !> Calculate v vector
  subroutine srhd_get_v(w,x,ixI^L,ixO^L,v,B2,VdotB)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,1:ndir)

    double precision, optional, intent(in)  :: VdotB(ixO^S),B2(ixO^S)
    double precision  :: sub_VdotB(ixO^S),sub_B2(ixO^S)
    integer :: idir

    if(present(B2)) then
       do idir=1,ndir
          v(ixO^S,idir) = (w(ixO^S, mom(idir)) + VdotB*w(ixO^S,mag(idir)))&
                      /(w(ixO^S,xi_)+B2)
       end do
    else  
     call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     do idir=1,ndir
         v(ixO^S,idir) = (w(ixO^S, mom(idir))+sub_VdotB*w(ixO^S,mag(idir)))&
                      /(w(ixO^S,xi_)+sub_B2)
     end do 
    end if 
  end subroutine srhd_get_v

!!! DM Modify...
  !> Calculate v component
  subroutine srhd_get_v_idim(w,ixI^L,ixO^L,idim,v_idim)
    use mod_global_parameters

    integer, intent(in)                    :: ixI^L, ixO^L, idim
    double precision, intent(in)           :: w(ixI^S,nw)
    double precision, intent(out)          :: v_idim(ixI^S)

    double precision                       :: sub_VdotB(ixO^S),sub_B2(ixO^S)

    call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
    v_idim(ixO^S) = (w(ixO^S, mom(idim)) + sub_VdotB*w(ixO^S,mag(idim)))&
                      /(w(ixO^S,xi_)+sub_B2)
  end subroutine srhd_get_v_idim

!!! DM Modify...
  !> Calculate v component (different interface, includes location vector x)
  subroutine srhd_get_v_idim_loc(w,x,ixI^L,ixO^L,idim,v_idim,B2,VdotB)
    use mod_global_parameters

    integer, intent(in)                    :: ixI^L, ixO^L, idim
    double precision, intent(in)           :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out)          :: v_idim(ixO^S)
    double precision, optional, intent(in) :: VdotB(ixO^S),B2(ixO^S)

    double precision                       :: sub_VdotB(ixO^S),sub_B2(ixO^S)

    if(present(B2) .and. present(VdotB))then
     v_idim(ixO^S) = (w(ixO^S, mom(idim)) + VdotB*w(ixO^S,mag(idim)))&
                      /(w(ixO^S,xi_)+B2)
    else
     call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     v_idim(ixO^S) = (w(ixO^S, mom(idim)) + sub_VdotB*w(ixO^S,mag(idim)))&
                      /(w(ixO^S,xi_)+sub_B2)
    end if
  end subroutine srhd_get_v_idim_loc

!!! DM Modify...
  !> Calculate v^2
  subroutine srhd_get_v2(w,x,ixI^L,ixO^L,v2,B2,VdotB)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v2(ixO^S)

    double precision, optional, intent(in)  :: VdotB(ixO^S),B2(ixO^S)
    double precision                        :: sub_VdotB(ixO^S),sub_B2(ixO^S)
    double precision                        :: v_num(ixO^S,1:ndir)
    integer :: idir

    if(present(B2)) then
     do idir=1,ndir
       v_num(ixO^S,idir)=w(ixO^S, mom(idir)) + VdotB*w(ixO^S,mag(idir))
     end do 
     v2(ixO^S) = sum(v_num(ixO^S, 1:ndir)**2.0 ,dim=ndim+1)&
                      /(w(ixO^S,xi_)+B2)**2.0
    else 
     call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     do idir=1,ndir
       v_num(ixO^S,idir)=w(ixO^S, mom(idir)) + sub_VdotB*w(ixO^S,mag(idir))
     end do 
     v2(ixO^S) = sum(v_num(ixO^S, :)**2.0 ,dim=ndim+1)&
                      /(w(ixO^S,xi_)+sub_B2)**2.0
    end if
  end subroutine srhd_get_v2

!!! DM Modify...
  !> maximal speed using gammie recipe
  subroutine srhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,w,rhoh,B2,&
                                   VdotB,calfven,vidim,csound2,cmax&
                                   ,patch_gammie,from_cbound,cmin)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(in), dimension(ixO^S) :: rhoh,B2,VdotB&
                                                     ,calfven, csound2
    double precision, intent(inout)          :: cmax(ixI^S)
    double precision, optional, intent(inout):: cmin(ixI^S)
    logical, optional, intent(in)            :: from_cbound
    logical, optional, intent(in)            :: patch_gammie(ixO^S)

    double precision, dimension(ixO^S):: A,B
    double precision, dimension(ixO^S):: vidim,v2
   is_patch : if(present(patch_gammie))then
    where(patch_gammie(ixO^S))
     A(ixO^S) = csound2(ixO^S)+calfven(ixO^S)-csound2(ixO^S)*calfven(ixO^S)
     B(ixO^S)=vidim(ixO^S)**2.0d0
    end where
    if(.not.present(from_cbound))&
           where(patch_gammie(ixO^S))vidim(ixO^S)=dabs(vidim(ixO^S))
    cond_onedir_patch: if(ndir==1)then
     where(patch_gammie(ixO^S))cmax(ixO^S)=(vidim(ixO^S)+dsqrt(A))&
                                           /(1.0d0+dsqrt(A*B))
     if(present(cmin))&
        where(patch_gammie(ixO^S))cmin(ixO^S)=(vidim(ixO^S)-dsqrt(A))&
                                               /(1.0d0+dsqrt(A*B))
    else cond_onedir_patch

     call srhd_get_v2(w,x,ixI^L,ixO^L,v2,B2=B2,VdotB=VdotB)
     where(patch_gammie(ixO^S))
      cmax(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  +dsqrt(A*(1.0d0-v2(ixO^S))*((1.0d0-v2(ixO^S)*A)&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
     end where
     if(present(cmin))&
        where(patch_gammie(ixO^S))cmin(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  -dsqrt(A*(1.0d0-v2(ixO^S))*(1.0d0-v2(ixO^S)*A&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
    end if cond_onedir_patch    
   else is_patch 
    A = csound2(ixO^S)+calfven(ixO^S)-csound2(ixO^S)*calfven(ixO^S)
    ! use B to save vidim**2
    if(.not.present(from_cbound).or..not.present(cmin))&
                                        vidim(ixO^S)=dabs(vidim(ixO^S))
    B(ixO^S)=vidim(ixO^S)**2.0d0
    cond_onedir: if(ndir==1)then
     cmax(ixO^S)=(vidim(ixO^S)+dsqrt(A))/(1.0d0+dsqrt(A*B))
     if(present(cmin))cmin(ixO^S)=(vidim(ixO^S)-dsqrt(A))/(1.0d0+dsqrt(A*B))
    else cond_onedir

     call srhd_get_v2(w,x,ixI^L,ixO^L,v2,B2=B2,VdotB=VdotB)
     cmax(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  +dsqrt(A*(1.0d0-v2(ixO^S))*((1.0d0-v2(ixO^S)*A)&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
     if(present(cmin))cmin(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  -dsqrt(A*(1.0d0-v2(ixO^S))*(1.0d0-v2(ixO^S)*A&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
    end if cond_onedir

   end if is_patch 

  end subroutine srhd_get_cmax_gammie

!!! DM Modify...
  !> Calculate cmax_idim using gammie method within ixO^L
  subroutine srhd_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    
    double precision, dimension(ixO^S):: rhoh,B2,VdotB,calfven,vidim
    double precision, dimension(ixO^S):: v2,csound2

    rhoh=w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_csound2(w,x,ixI^L,ixO^L,rhoh,csound2)
    call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2,VdotB)
    call srhd_get_calfven2(ixI^L,ixO^L,x,w,.true.,calfven,rhoh=rhoh,B2=B2)
    call srhd_get_v_idim_loc(w,x,ixI^L,ixO^L,idim,vidim,B2=B2,VdotB=VdotB)

    call srhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,w,rhoh=rhoh,B2=B2,&
                               VdotB=VdotB,calfven=calfven,vidim=vidim,&
                               csound2=csound2,cmax=cmax)

  end subroutine srhd_get_cmax

!!! DM Modify...
  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine srhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,cmax,cmin)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)              :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)              :: x(ixI^S,1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision                   :: wmean(ixI^S,nw)
    double precision, dimension(ixO^S) :: umean, dmean, csound2Lp, &
                                          csound2Rp, tmp1,tmp2,tmp3, &
                                          B2Lp,B2Rp,VdotBLp,VdotBRp, &
                                          rhohLp,rhohRp, &
                                          cmaxL,cmaxR,csound2,rhoh,&
                                          B2,VdotB, vidim,vidimLp,vidimRp,&
                                          calfvenLp,calfvenRp,calfven
    character(len=30)                  :: subname_loc

    subname_loc='srhd_get_cbounds'
    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      rhohLp=wLp(ixO^S,xi_)/wLp(ixO^S,lfac_)**2.0d0
      rhohRp=wRp(ixO^S,xi_)/wRp(ixO^S,lfac_)**2.0d0

      call srhd_get_csound2(wLp,x,ixI^L,ixO^L,rhohLp,csound2Lp)
      call srhd_get_csound2(wRp,x,ixI^L,ixO^L,rhohRp,csound2Rp)
      call srhd_get_B2andVdotB(ixI^L,ixO^L,wLp,.false.,B2=B2Lp,VdotB=VdotBLp)
      call srhd_get_B2andVdotB(ixI^L,ixO^L,wRp,.false.,B2=B2Rp,VdotB=VdotBRp)

      call srhd_get_calfven2(ixI^L,ixO^L,x,wLp,.false.,calfvenLp,&
                              rhoh=rhohLp,B2=B2Lp)
      call srhd_get_calfven2(ixI^L,ixO^L,x,wRp,.false.,calfvenRp,&
                              rhoh=rhohRp,B2=B2Rp)

      call srhd_get_v_idim_loc(wLp,x,ixI^L,ixO^L,idim,vidimLp,B2=B2Lp,&
                                VdotB=VdotBLp)
      call srhd_get_v_idim_loc(wRp,x,ixI^L,ixO^L,idim,vidimRp,B2=B2Rp,&
                                VdotB=VdotBRp)

      tmp1(ixO^S)=sqrt(wLp(ixO^S,xi_)+B2Lp(ixO^S))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,xi_)+B2Rp(ixO^S))
      tmp3(ixO^S)=1.d0/(sqrt(wLp(ixO^S,xi_)+B2Lp(ixO^S))&
                  +sqrt(wRp(ixO^S,rho_)+B2Rp(ixO^S)))

      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)&
                    +wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

     call srhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,wLp,rhoh=rhohLp,B2=B2Lp,&
                                VdotB=VdotBLp,calfven=calfvenLp,vidim=vidimLp,&
                                csound2=csound2Lp,cmax=cmaxL)

     call srhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,wRp,rhoh=rhohRp,B2=B2Rp,&
                                VdotB=VdotBRp,calfven=calfvenRp,vidim=vidimRp,&
                                csound2=csound2Rp,cmax=cmaxR)


      dmean(ixO^S)=(tmp1(ixO^S)*cmaxL(ixO^S)**2&
                    +tmp2(ixO^S)*cmaxR(ixO^S)**2)*tmp3(ixO^S)+&
                    0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
                    (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2
      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S)=abs(umean(ixO^S))+dmean(ixO^S)
      end if
    else
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      ! get auxiliary variables
      call srhd_get_auxiliary(.true.,wmean,x,ixI^L,ixO^L,subname_loc)
      rhoh=wmean(ixO^S,xi_)/wmean(ixO^S,lfac_)**2.0d0

      call srhd_get_csound2(wmean,x,ixI^L,ixO^L,rhoh,csound2)
      call srhd_get_B2andVdotB(ixI^L,ixO^L,wmean,.true.,B2=B2,VdotB=VdotB)
      call srhd_get_calfven2(ixI^L,ixO^L,x,wmean,.true.,&
                               calfven,rhoh=rhoh,B2=B2) 
      call srhd_get_v_idim_loc(wmean,x,ixI^L,ixO^L,idim,vidim,B2=B2,&
                                VdotB=VdotB)
      call srhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,wmean,rhoh=rhoh,B2=B2,&
                                 VdotB=VdotB,calfven=calfven,&
                                 vidim=vidim,csound2=csound2,cmax=cmax&
                                 ,from_cbound=.true.,cmin=cmin)
      if(present(cmin)) then
        cmax(ixO^S)=min(max(cmax(ixO^S),0.0D0),1.0d0)
        cmin(ixO^S)=max(min(cmin(ixO^S),0.0D0),-1.0d0)
      else
        cmax(ixO^S)=min(cmax(ixO^S),1.0d0)
      end if
    end if
  end subroutine srhd_get_cbounds

!!! DM Modify and delete if not needed
  !> Calculate the Alfven speed
  subroutine srhd_get_calfven2(ixI^L,ixO^L,x,w,conserve,calfven,rhoh,B2)
   use mod_global_parameters
   implicit none
   integer, intent(in)           :: ixI^L, ixO^L
   double precision, intent(in)  :: w(ixI^S, 1/nw), x(ixI^S,1:ndim)
   logical,          intent(in)  :: conserve
   double precision, intent(out) :: calfven(ixO^S)
   double precision, optional    :: rhoh(ixO^S),B2(ixO^S)

   double precision              :: sub_B2(ixO^S),sub_rhoh(ixO^S)

   if(present(B2))then
    calfven=B2/(B2+rhoh)
   else
    sub_rhoh=w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_B2andVdotB(ixI^L,ixO^L,w,conserve,B2=sub_B2)
    calfven=sub_B2/(sub_B2+sub_rhoh)
   end if
  end subroutine srhd_get_calfven2

!!! DM Modify and probably delete
  !> Calculate magnetic pressure within ixO^L 
  subroutine srhd_get_p_mag(w,ixI^L,ixO^L,pmag)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(out)   :: pmag(ixI^S)

    pmag(ixO^S) = 0.5d0 * (sum(w(ixO^S, mag(:))**2, dim=ndim+1)&
            +sum(w(ixO^S, mag(:))*w(ixO^S, mom(:)), dim=ndim+1)&
               /w(ixO^S,lfac_))&
               / w(ixO^S,lfac_)

  end subroutine srhd_get_p_mag

!!! DM This is probably ok
  !> Calculate sound speed
  subroutine srhd_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: inv_rho(ixO^S)

    inv_rho=1.d0/w(ixO^S,rho_)

    if(srhd_energy) then
      csound(ixO^S)=dsqrt(srhd_gamma*w(ixO^S,p_)*inv_rho)
    else
      csound(ixO^S)=dsqrt(srhd_gamma*srhd_adiab*w(ixO^S,rho_)**gamma_1)
    end if

  end subroutine srhd_get_csound_prim

!!! DM This is probably ok
  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  subroutine srhd_get_csound2(w,x,ixI^L,ixO^L,rhoh,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw),rhoh(ixO^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixO^S)
  
    double precision                :: rho(ixO^S)

    if(srhd_energy) then
      rho=w(ixO^S,d_)/w(ixO^S,lfac_)
      call srhd_get_csound2_eos(ixI^L,ixO^L,x,rho,rhoh,csound2)
    else
      csound2(ixO^S)=srhd_gamma*srhd_adiab*rho(ixO^S)**gamma_1
    end if

  end subroutine srhd_get_csound2

!!! DM Not needed probably
  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine srhd_get_p_total(w,x,ixI^L,ixO^L,p,B2,VdotB)
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: p(ixI^S)
    double precision, optional, intent(in) :: B2(ixO^S),VdotB(ixO^S)

    double precision                :: pmag(ixI^S)

    call srhd_get_pthermal(w,x,ixI^L,ixO^L,p)
    call srhd_get_pmag(w,x,ixI^L,ixO^L,pmag,B2=B2,VdotB=VdotB)
    p(ixO^S) = p(ixO^S) + pmag(ixO^S)

  end subroutine srhd_get_p_total

!!! DM Delete 
  !> compute magnetic pressure (different interface as get_p_mag)
  subroutine srhd_get_pmag(w,x,ixI^L,ixO^L,pmag,B2,VdotB)
    use mod_global_parameters
    implicit none

    integer, intent(in)                    :: ixI^L, ixO^L
    double precision, intent(in)           :: w(ixI^S,nw)
    double precision, intent(in)           :: x(ixI^S,1:ndim)
    double precision, intent(out)          :: pmag(ixI^S)
    double precision, optional, intent(in) :: B2(ixO^S),VdotB(ixO^S)

    double precision, dimension(ixO^S)     :: sub_B2,sub_VdotB

    if(present(B2)) then
       pmag(ixO^S) = 0.5d0*(VdotB(ixO^S)**2.0d0  &
                       + B2(ixO^S)/w(ixO^S,lfac_)**2.0d0)
    else 
       call srhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
       pmag(ixO^S) = 0.5d0*(sub_VdotB(ixO^S)**2.0d0  &
                       + sub_B2(ixO^S)/w(ixO^S,lfac_)**2.0d0)
    end if 
  end subroutine srhd_get_pmag

!!! DM Modify...
  !> Calculate fluxes within ixO^L.
  subroutine srhd_get_flux(wC,wP,x,ixI^L,ixO^L,idim,f)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: wP(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: ptotal(ixI^S)
    double precision             :: B2(ixO^S),VdotB(ixO^S)
    double precision             :: v(ixI^S,1:ndir)
    integer                      :: idirmin, iw, idir


    call srhd_get_B2andVdotB(ixI^L,ixO^L,wP,.false.,B2=B2,VdotB=VdotB)
    call srhd_get_p_total(wC,x,ixI^L,ixO^L,ptotal,B2=B2,VdotB=VdotB)
    call srhd_get_v(wC,x,ixI^L,ixO^L,v,B2=B2,VdotB=VdotB)
    ! Get flux of density
    f(ixO^S,rho_)=wP(ixO^S,mom(idim))*wP(ixO^S,rho_)

    ! Get flux of tracer
    do iw=1,srhd_n_tracer
      f(ixO^S,tracer(iw))=wP(ixO^S,mom(idim))*wP(ixO^S,tracer(iw))
    end do

!!! DM Modify...
    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-B_k*b_i [+ptotal if i==k]
    Loop_idir_mom : do idir=1,ndir
      f(ixO^S,mom(idir))= v(ixO^S,idir)*wC(ixO^S,mom(idim))&
                         -wP(ixO^S,mag(idir))*(VdotB*v(ixO^S,idim)&
                         +wP(ixO^S,mag(idim))/wP(ixO^S,lfac_)**2.0d0)

    end do Loop_idir_mom
    f(ixO^S,mom(idim))=ptotal(ixO^S)+f(ixO^S,mom(idim))

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    is_notiso : if (srhd_energy) then
       is_internal : if (block%e_is_internal) then
          f(ixO^S,e_)=v(ixO^S,idim)*wP(ixO^S,p_)
       else is_internal
          f(ixO^S,e_)=v(ixO^S,idim)*(wC(ixO^S,e_) + ptotal(ixO^S))- &
             wP(ixO^S,mag(idim))*VdotB(ixO^S)
          if(type_divb==divb_glm2) then
            f(ixO^S,e_) = f(ixO^S,e_) + vmax_global*wP(ixO^S,psi_)&
                                         *wP(ixO^S,mag(idim))
          end if

       end if is_internal
    end if is_notiso

!!! DM delete
    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    Loop_idir_Bflux : do idir=1,ndir
      is_idim_Bflux : if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        is_glm_Bflux : if (srhd_glm) then
           if(type_divb==divb_glm1) then
             f(ixO^S,mag(idir))=wP(ixO^S,psi_)
           else
             f(ixO^S,mag(idir))=vmax_global*wP(ixO^S,psi_)
           end if
        else is_glm_Bflux
           f(ixO^S,mag(idir))=zero
        end if is_glm_Bflux
      else is_idim_Bflux
        f(ixO^S,mag(idir))=v(ixO^S,idim)*wP(ixO^S,mag(idir))&
                           -wP(ixO^S,mag(idim))*v(ixO^S,idir)

      end if is_idim_Bflux
    end do Loop_idir_Bflux

  end subroutine srhd_get_flux

!!! DM Deleted the divB related sources 
 !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine srhd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if (srhd_energy .and. block%e_is_internal) then
        active = .true.
        call internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
      endif
    endif
 
    call srhd_get_auxiliary(.true.,w,x,ixI^L,ixO^L,'srhd_add_source')

  end subroutine srhd_add_source

!!! DM Modify...
  subroutine internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call srhd_get_v(wCT,x,ixI^L,ixI^L,v)
    call divvector(v,ixI^L,ixO^L,divv)
    call srhd_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
    w(ixO^S,e_)=w(ixO^S,e_)-qdt*pth(ixO^S)*divv(ixO^S)

  end subroutine internal_energy_add_source

!!! DM All divB related routines deleted...

  subroutine srhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    dtnew = bigdouble

  end subroutine srhd_get_dt

!!! DM Modify...
  ! Add geometrical source terms to w
  subroutine srhd_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),ptot(ixI^S),tmp2(ixI^S)
    double precision :: B2(ixO^S),VdotB(ixO^S),VdotB2(ixO^S),xiplusB2(ixO^S)
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_
    character(len=30)  :: subname_loc

    if(typeaxial /= 'slab')then
     subname_loc='srhd_add_geom'
     ! get auxiliary variables
     call srhd_get_auxiliary(.true.,wCT,x,ixI^L,ixO^L,subname_loc)
     call srhd_get_B2andVdotB(ixI^L,ixO^L,wCT,.true.,B2=B2,VdotB=VdotB)
     call srhd_get_p_total(wCT,x,ixI^L,ixO^L,ptot)
     VdotB2=VdotB**2.0
     xiplusB2(ixO^S) = wCT(ixO^S,xi_)**2.0+B2(ixO^S)
    
     mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
     br_=mag(1); bphi_=mag(1)-1+phi_
    end if

    select case (typeaxial)
    case ('cylindrical')
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
       if(phi_>0) then
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*(ptot(ixO^S)+&
                   (wCT(ixO^S,mphi_)**2-VdotB2(ixO^S)*wCT(ixO^S,bphi_)**2)&
                   /xiplusB2(ixO^S)-wCT(ixO^S,bphi_)**2/wCT(ixO^S,lfac_)**2.0)

         w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt/x(ixO^S,1)*(&
                  -(wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)&
                  -VdotB2(ixO^S)*wCT(ixO^S,br_)*wCT(ixO^S,bphi_))&
                  /xiplusB2(ixO^S) &
                  +wCT(ixO^S,bphi_)*wCT(ixO^S,br_)/wCT(ixO^S,lfac_)**2)

         w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt/x(ixO^S,1)*&
                  (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                  -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                  /xiplusB2(ixO^S)
       else
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*ptot(ixO^S)
       end if
       if(srhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
    case ('spherical')
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       
       ! m1
       tmp(ixO^S)=ptot(ixO^S)*x(ixO^S,1) &
         *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1))/block%dvolume(ixO^S)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(idir))**2.0&
                         -(VdotB2(ixO^S)*B2(ixO^S)))&
                        /(xiplusB2(ixO^S))&
                      -wCT(ixO^S,mag(idir))**2/wCT(ixO^S,lfac_)**2.0
         end do
       end if
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b1
       if(srhd_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt/x(ixO^S,1)*2.0d0*wCT(ixO^S,psi_)
       end if

       {^NOONED
       ! m2
       !tmp(ixO^S)=tmp1(ixO^S)
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*ptot(ixO^S) &
            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
            /block%dvolume(ixO^S)

       tmp(ixO^S)=-((wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))&
                   -VdotB2(ixO^S)*wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))&
                    /xiplusB2(ixO^S) &
            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2))/wCT(ixO^S,lfac_)**2.0)

       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+((wCT(ixO^S,mom(3))**2&
              -VdotB2(ixO^S)*wCT(ixO^S,mag(3))**2.0)/xiplusB2(ixO^S) &
              -(wCT(ixO^S,mag(3))/wCT(ixO^S,lfac_))**2.0)&
              *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
       end if
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b2
       tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
            -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))/xiplusB2(ixO^S)

       if(srhd_glm) then
         tmp(ixO^S)=tmp(ixO^S) &
              + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
       end if
       w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       }

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixO^S)=-((wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))&
                       -VdotB2(ixO^S)*wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1)))&
                       /xiplusB2(ixO^S) &
                -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))&
                /wCT(ixO^S,lfac_)**2.0d0) {^NOONED &
                -((wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))&
                -VdotB2(ixO^S)*wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3)))&
                /xiplusB2(ixO^S) &
                -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))/wCT(ixO^S,lfac_)**2.0d0) &
                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
              -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))/xiplusB2(ixO^S) {^NOONED &
              -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
              /(xiplusB2(ixO^S)*dsin(x(ixO^S,2))) }
         w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
       end if
    end select
  end subroutine srhd_add_source_geom
  
!!! DM Bfield & Benergy deleted
   
  !> compute kinetic energy from primitive variables
  function srhd_kin_en_primitive(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
       ke = (w(ixO^S,xi_)-1.0d0) * inv_rho
    else
       ke = (w(ixO^S,xi_)-1.0d0) / w(ixO^S, rho_)
    end if
  end function srhd_kin_en_primitive

!!! DM Modify...
  subroutine srhd_boundary_adjust
    use mod_global_parameters
    integer :: iB, idim, iside, iigrid, igrid
    integer :: ixG^L, ixO^L, i^D

    ixG^L=ixG^LL;
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>ps(igrid)
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        do idim=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
           do iside=1,2
              i^D=kr(^D,idim)*(2*iside-3);
              if (neighbor_type(i^D,igrid)/=1) cycle
              iB=(idim-1)*2+iside
              if(.not.boundary_divbfix(iB)) cycle
              if(any(typeboundary(:,iB)=="special")) then
                ! MF nonlinear force-free B field extrapolation and data driven
                ! require normal B of the first ghost cell layer to be untouched by
                ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
                select case (idim)
                {case (^D)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin^DD=ixGmax^D+1-nghostcells+boundary_divbfix_skip(2*^D)^D%ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmax^DD;
                   else
                      ! minimal boundary
                      ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmin^D-1+nghostcells-boundary_divbfix_skip(2*^D-1)^D%ixOmax^DD=ixGmax^DD;
                   end if \}
                end select
                call fixdivB_boundary(ixG^L,ixO^L,ps(igrid)%w,ps(igrid)%x,iB)
              end if
           end do
        end do
     end do

  end subroutine srhd_boundary_adjust

!!! DM Deleted divB boundary
   
!!! DM Modify...
   !> subroutine using srhd_con2prim to calculate the enthalpy and the lorentz factor
   subroutine srhd_get_auxiliary(clipping,w,x,ixI^L,ixO^L,subname)
    use mod_global_parameters
    use mod_srhd_con2prim
    implicit none

    logical, intent(in)             :: clipping
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    character(len=*), intent(in)    :: subname
  
    integer                        :: ix^D,ierror
    integer                        :: flag_error(ixO^S)
    character(len=len(subname)+30) :: subname_loc
    double precision               :: mom_vec(ndir), mag_vec(ndir)

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    mom_vec = w(ix^D,mom(1):mom(ndir))
    mag_vec = w(ix^D,mag(1):mag(ndir))
    call srhd_con2prim(w(ix^D,d_),mom_vec,w(ix^D,e_),&
             mag_vec,w(ix^D,lfac_),w(ix^D,xi_),ierror)
    if(check_small_values)then
     if(ierror/=0) then
       flag_error(ix^D) = ierror
     else
       flag_error(ix^D) = 0
     end if 
     
    end if
   {enddo^D&\} 
   is_check_small : if(check_small_values)then
    is_flag_on : if(any(flag_error(ixO^S)/=0))then
     subname_loc='srhd_get_auxiliary from -> '//trim(subname)
     call srhd_handle_small_values(.false., &
                                    w, x, &
                                    ixI^L, ixO^L,subname_loc,&
                                    flag_error=flag_error)
    end if is_flag_on
   end if is_check_small
   end subroutine srhd_get_auxiliary

   subroutine srhd_get_4u_from_3v(ixI^L,ixO^L,vtou)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: vtou(ixI^S,1:ndir)

    double precision                :: lfac(ixO^S)
    integer                         :: idir

    lfac= 1.0d0/dsqrt(1.0d0-sum(vtou(ixO^S,1:ndir)**2.0,dim=ndim+1))
    do idir=1,ndir
       vtou(ixO^S,idir)=lfac(ixO^S)*vtou(ixO^S,idir)
    end do   
   end subroutine srhd_get_4u_from_3v
end module mod_srhd_phys
