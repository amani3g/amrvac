!> Parameter module for Special Relativistic Hydrodynamics 
module mod_srmhd_parameters
  ! made by Z. MELIANI 14/02/2018
  ! simplified by R. Keppens, adapted by D. MIllas (January 2019)
  use mod_global_parameters, only: std_len
  implicit none

   !> Whether sole energy equation
  logical, public                  :: srhd_energy = .true.

  !> Whether synge eos (T) is used, or simple polytrope (F)
  logical, public                  :: srhd_eos = .false.

  !> Whether particles module is added
  logical, public                  :: srhd_particles = .false.

  !> Number of tracer species
  integer, public                  :: srhd_n_tracer = 0

  !> Index of the moving frame density (in the w array)
  integer, public                  :: rho_

   !> Index of the lab frame density (in the w array)
  integer, public                  :: d_

  !> Indices of the momentum density
  integer, allocatable, public     :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public                  :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                  :: p_

  !> Indices of the Lorentz factor (auxiliary variable)
  integer, public     :: lfac_

  !> Indices of the inertia (auxiliary variable)
  integer, public     :: xi_

  !> Indices of the tracers
  integer, allocatable, public     :: tracer(:)

  !> Helium abundance over Hydrogen
  double precision, public      :: He_abundance=0.0d0

  !> The adiabatic index and derived values
  double precision, public                :: srhd_gamma = 5.d0/3.0d0
  double precision, public                :: gamma_1, inv_gamma_1, gamma_to_gamma_1

  !> The adiabatic constant
  double precision, public                :: srhd_adiab = 1.0d0

  !> Newton Raphson related variables for conservative to primitive
  logical,          public         :: srhd_checkNR   = .true.
  integer,          public         :: srhd_maxiterationNR=100
  double precision, public         :: srhd_absaccNR  = 1.0d-8
  double precision, public         :: srhd_tolerNR   = 1.0d-9
  double precision, public         :: srhd_maxdspeed = 1.0d-7 ! max difference from 1
  double precision, public         :: srhd_maxspeed  = 0.9999d0 ! i.e. Lfac about 100

  !> The smallest allowed energy
  double precision                 :: small_e

  !> The smallest allowed inertia
  double precision                 :: small_xi

  !> small value for speed
  double precision, public         :: small_vec2  = 0.0

  !> The number of waves
  integer :: nwwave=8

end module mod_srmhd_parameters
