!> Parameter module for Special Relativistic Magnetohydrodynamics 
! Updated and simplified on january 2019, R. Keppens
module mod_srmhd_parameters
  ! made by Z. MELIANI 14/02/2018
  use mod_global_parameters, only: std_len
  implicit none

   !> Whether sole energy equation
  logical, public                  :: srmhd_energy = .true.

  !> Whether synge eos (T) is used, or simple polytrope (F)
  logical, public                  :: srmhd_eos = .false.

  !> Whether particles module is added
  logical, public                  :: srmhd_particles = .false.

  !> Number of tracer species
  integer, public                  :: srmhd_n_tracer = 0

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

  !> Indices of the magnetic field
  integer, allocatable, public     :: mag(:)

  !> Indices of the GLM psi
  integer, public     :: psi_

  !> Indices of the Lorentz factor (auxiliary variable)
  integer, public     :: lfac_

  !> Indices of the inertia (auxiliary variable)
  integer, public     :: xi_

  !> Indices of the tracers
  integer, allocatable, public     :: tracer(:)

  !> Helium abundance over Hydrogen
  double precision, public      :: He_abundance=0.0d0

  !> The adiabatic index and derived values
  double precision, public                :: srmhd_gamma = 5.d0/3.0d0
  double precision, public                :: gamma_1, inv_gamma_1, gamma_to_gamma_1

  !> The adiabatic constant
  double precision, public                :: srmhd_adiab = 1.0d0

  !> Newton Raphson related variables for conservative to primitive
  logical,          public         :: srmhd_checkNR   = .true.
  integer,          public         :: srmhd_maxiterationNR=100
  double precision, public         :: srmhd_absaccNR  = 1.0d-8
  double precision, public         :: srmhd_tolerNR   = 1.0d-9
  double precision, public         :: srmhd_maxdspeed = 1.0d-7 ! max difference from 1
  double precision, public         :: srmhd_maxspeed  = 0.9999d0 ! i.e. Lfac about 100

  !> The smallest allowed energy
  double precision                 :: small_e

  !> The smallest allowed inertia
  double precision                 :: small_xi

  !> small value for speed
  double precision, public         :: small_vec2  = 0.0

  !> The number of waves
  integer :: nwwave=8

  !> Whether GLM is used for DivB
  logical, public                  :: srmhd_glm = .false.

  !> Whether divB cleaning sources are added split from fluid solver
  logical, public                  :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public         :: srmhd_glm_alpha = 0.5d0

  !> MHD fourth order evaluation of divB
  logical, public                  :: srmhd_4th_order = .false.

  !> Method type to clean divergence of B
  character(len=std_len), public     :: typedivbfix  = 'linde'

  !> DivB Method type in integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> To control divB=0 fix for boundary
  logical, public     :: boundary_divbfix(2*^ND)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public     :: boundary_divbfix_skip(2*^ND)=0

  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_glm1          = 1
  integer, parameter :: divb_glm2          = 2
  integer, parameter :: divb_janhunen      = 3
  integer, parameter :: divb_linde         = 4
  integer, parameter :: divb_lindejanhunen = 6
  integer, parameter :: divb_lindeglm      = 7

end module mod_srmhd_parameters
