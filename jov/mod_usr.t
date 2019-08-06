module FFTW3
  !-------------------------------------------------------------------------------
  !
  !   This must be included in any routine that uses FFTW (version 3.3.3).
  !
  !   Copyright (c) 2016, SHTOOLS
  !   All rights reserved.
  !
  !-------------------------------------------------------------------------------
        INTEGER FFTW_R2HC
        PARAMETER (FFTW_R2HC=0)
        INTEGER FFTW_HC2R
        PARAMETER (FFTW_HC2R=1)
        INTEGER FFTW_DHT
        PARAMETER (FFTW_DHT=2)
        INTEGER FFTW_REDFT00
        PARAMETER (FFTW_REDFT00=3)
        INTEGER FFTW_REDFT01
        PARAMETER (FFTW_REDFT01=4)
        INTEGER FFTW_REDFT10
        PARAMETER (FFTW_REDFT10=5)
        INTEGER FFTW_REDFT11
        PARAMETER (FFTW_REDFT11=6)
        INTEGER FFTW_RODFT00
        PARAMETER (FFTW_RODFT00=7)
        INTEGER FFTW_RODFT01
        PARAMETER (FFTW_RODFT01=8)
        INTEGER FFTW_RODFT10
        PARAMETER (FFTW_RODFT10=9)
        INTEGER FFTW_RODFT11
        PARAMETER (FFTW_RODFT11=10)
        INTEGER FFTW_FORWARD
        PARAMETER (FFTW_FORWARD=-1)
        INTEGER FFTW_BACKWARD
        PARAMETER (FFTW_BACKWARD=+1)
        INTEGER FFTW_MEASURE
        PARAMETER (FFTW_MEASURE=0)
        INTEGER FFTW_DESTROY_INPUT
        PARAMETER (FFTW_DESTROY_INPUT=1)
        INTEGER FFTW_UNALIGNED
        PARAMETER (FFTW_UNALIGNED=2)
        INTEGER FFTW_CONSERVE_MEMORY
        PARAMETER (FFTW_CONSERVE_MEMORY=4)
        INTEGER FFTW_EXHAUSTIVE
        PARAMETER (FFTW_EXHAUSTIVE=8)
        INTEGER FFTW_PRESERVE_INPUT
        PARAMETER (FFTW_PRESERVE_INPUT=16)
        INTEGER FFTW_PATIENT
        PARAMETER (FFTW_PATIENT=32)
        INTEGER FFTW_ESTIMATE
        PARAMETER (FFTW_ESTIMATE=64)
        INTEGER FFTW_WISDOM_ONLY
        PARAMETER (FFTW_WISDOM_ONLY=2097152)
        INTEGER FFTW_ESTIMATE_PATIENT
        PARAMETER (FFTW_ESTIMATE_PATIENT=128)
        INTEGER FFTW_BELIEVE_PCOST
        PARAMETER (FFTW_BELIEVE_PCOST=256)
        INTEGER FFTW_NO_DFT_R2HC
        PARAMETER (FFTW_NO_DFT_R2HC=512)
        INTEGER FFTW_NO_NONTHREADED
        PARAMETER (FFTW_NO_NONTHREADED=1024)
        INTEGER FFTW_NO_BUFFERING
        PARAMETER (FFTW_NO_BUFFERING=2048)
        INTEGER FFTW_NO_INDIRECT_OP
        PARAMETER (FFTW_NO_INDIRECT_OP=4096)
        INTEGER FFTW_ALLOW_LARGE_GENERIC
        PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
        INTEGER FFTW_NO_RANK_SPLITS
        PARAMETER (FFTW_NO_RANK_SPLITS=16384)
        INTEGER FFTW_NO_VRANK_SPLITS
        PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
        INTEGER FFTW_NO_VRECURSE
        PARAMETER (FFTW_NO_VRECURSE=65536)
        INTEGER FFTW_NO_SIMD
        PARAMETER (FFTW_NO_SIMD=131072)
        INTEGER FFTW_NO_SLOW
        PARAMETER (FFTW_NO_SLOW=262144)
        INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
        PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
        INTEGER FFTW_ALLOW_PRUNING
        PARAMETER (FFTW_ALLOW_PRUNING=1048576)

  end module FFTW3

module mod_usr
  use mod_mhd
  use mod_particles

  implicit none

  double precision :: charge = 1.60218d-19 ![C]
  double precision :: mass = 9.10938d-31 ![kg]

  ! Initial position (in m)
  double precision, parameter :: rj = 4.192d7![m]
  double precision :: x0(3) = [5*rj, 0.0d0, 0.0d0]

  ! Initial velocity (in m/s) 
  ! the speed of light const_c is defined in cgs in amrvac !
  double precision :: v0(3) = [0.0d0, 0.0d0, 0.0d0]

  ! Maxwellian velocity (per component vx, vy, vz)
  double precision :: maxwellian_velocity = 3.42527d6

  double precision, parameter :: not_used_value = -1.0d20
  double precision :: force_E0(3) = [not_used_value, 0.0d0, 0.0d0]
  double precision :: force_B0(3) = [not_used_value, 0.0d0, 0.0d0]

  ! Use an analytic field instead of an interpolated one
  logical :: use_analytic_field = .true. !GIVES A SEGMENTATION ERROR, INTERPOLATING ONLY WITH GCA
  ! If false, then in each cell of the grid the value of the grid will be interpolated
  ! If true, it will just read the x value that it's given

  ! Dipole Vectors csv file name 
  character*(30), parameter :: dipole_coord_bcsv = 'dipole_grid/dipole_vectors.csv'

  ! JRM09 Vectors csv file name 
  character*(29), parameter :: jrm09_coord_bcsv = 'jrm09_grid/jrm09_vectors.csv'

  ! Coordinate vectors of 3d bfield grid solution
  double precision :: r_vec(22), theta_vec(22), phi_vec(22), mult_vec(22)


  contains

    subroutine usr_init()
      use mod_initialize

      !declare integer i
      integer :: i

      unit_length        = 1.d0
      unit_numberdensity = 1.d0
      unit_velocity      = 1.0d0

      usr_init_one_grid => initonegrid_usr
      usr_create_particles => generate_particles
      usr_particle_fields => set_custom_field

      call set_coordinate_system("Cartesian_3D")
      call mhd_activate()
      call params_read(par_files)

      call initialize_amrvac()    ! So that we have settings available

      if (use_analytic_field) then
        if (physics_type_particles /= 'Lorentz' .and. physics_type_particles /= 'Vay' .and. physics_type_particles /= 'HC') &
             call mpistop('Analytic fields only supported with Boris, HC or Vay schemes')
        usr_particle_analytic => get_analytic_field
      end if

      ! initialize jrm09 coordinate vectors
      ! open coordinate csv
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! If this is the dipole case this needs to be changed!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open(1, file = dipole_coord_bcsv)
      do i = 1, 22
        read(1,*) mult_vec(i), r_vec(i), theta_vec(i), phi_vec(i)
      enddo
      print*, mult_vec(i), r_vec, theta_vec, phi_vec

    end subroutine usr_init

    !> Read parameters from a file
    subroutine params_read(files)
      character(len=*), intent(in) :: files(:)
      integer                      :: n

      namelist /my_list/ charge, mass, x0, v0, use_analytic_field, force_E0, &
           force_B0, maxwellian_velocity

      do n = 1, size(files)
         open(unitpar, file=trim(files(n)), status="old")
         read(unitpar, my_list, end=111)
  111    close(unitpar)
      end do

    end subroutine params_read

    subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
      ! initialize one grid
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w(ixI^S,1:nw)

      double precision:: rho0,p0,b0
      logical, save :: first=.true.

      rho0 = 1.0d0
      p0 = 1.0d0
      b0 = 1.0d0

      w(ixO^S,rho_)= rho0
      w(ixO^S,mom(1))= 0.0d0
      w(ixO^S,mom(2))= 0.0d0
      w(ixO^S,mom(3))= 0.0d0
      w(ixO^S,p_)=p0
      w(ixO^S,mag(1))= 0.0d0
      w(ixO^S,mag(2))= 0.0d0
      w(ixO^S,mag(3))= b0

      call mhd_to_conserved(ixI^L,ixO^L,w,x)

      if (first .and. mype==0 )then
        write(*,*) 'Test particles in 3D simple B-field'
        write(*,*) 'rho - p - b:',rho0,p0,b0
        first=.false.
        print *, "Running test case", iprob
      endif

    end subroutine initonegrid_usr

    subroutine generate_particles(n_particles, x, v, q, m, follow)
      use mod_particles
      integer, intent(in)           :: n_particles
      double precision, intent(out) :: x(3, n_particles)
      double precision, intent(out) :: v(3, n_particles)
      double precision, intent(out) :: q(n_particles)
      double precision, intent(out) :: m(n_particles)
      logical, intent(out)          :: follow(n_particles)
      integer                       :: n

      do n = 1, n_particles
        call get_particle(x(:, n), v(:, n), q(n), m(n), n, n_particles, iprob)
      end do

      follow(:) = .true.

      ! Scale to CGS units
      x = x * 1d2             ! m to cm
      v = v * 1d2             ! to cm/s
      q = q * const_c * 0.1d0 ! A unit charge converted to CGS units
      m = m * 1d3             ! kg to gram
      !print*,'mass in cgs:',m,'charge in cgs:',q,'q/m ratio',q/m
    end subroutine generate_particles

    ! Return field at location x (in SI units: Tesla, V/m)
    subroutine get_field(x, E, B)
      double precision, intent(in)  :: x(3)
      double precision, intent(out) :: E(3), B(3)

      ! declare: i, j, k, r_index, theta_index, phi_index as integers
      integer :: i, j, k, r_index, theta_index, phi_index
      ! declare: r_delta, theta_delta, phi_delta as 22long arrays, 
      double precision :: r_delta(22), theta_delta(22), phi_delta(22)
      ! declare: x_sphere, B_sphere as 3long array
      double precision :: x_sphere(3), B_sphere(3)
      ! declare r_min, theta_min, phi_min as double precision
      double precision :: r_min, theta_min, phi_min    
      ! character formater for bfield_grid/##_surface
      character(8) :: fmt 
      integer :: r_mult
      character(2) :: r_mult_str
      ! file names 
      character(36) :: Br_file3, Btheta_file3, Bphi_file3 !for iprob 3, dipole
      character(35) :: Br_file4, Btheta_file4, Bphi_file4 !for iprob 4, jrm09
      ! tmp arrays for reading rows (corresponding to a theta value) from files
      double precision :: Br_tmp(22), Btheta_tmp(22), Bphi_tmp(22)
      ! integers for reading through csv files
      integer :: i_, j_, k_


      select case (iprob)
      case (1)
        ! Magnetic dipole (run up to t = 100)
        E = [0.0d0, 0.0d0, 0.0d0]

        ! x is in cm, this corresponds to B = 10 T at 1 m
        ! M = 10 G * 10^4 T/m * 10^2 cm/m
        ! Set M = 4.176 G * 10^4 T/m * 10^2 cm/m
        B = 4.176 * 1d6 * [3d0 * x(1) * x(3), &
             3d0 * x(2) * x(3), &
             2d0 * x(3)**2 - x(1)**2 - x(2)**2] / &
             (x(1)**2 + x(2)**2 + x(3)**2)**2.5d0

      case (2)
        ! Magnetic dipole (run up to t = 100)
        ! Checkin to make sure there are no issues with 
        ! translating between cartesian to spherical coordinates
        E = [0.0d0, 0.0d0, 0.0d0]

        ! Convert x to x_spherical
        ! Point Transformation:
        ! r = sqrt(x^2 + y^2 + z^2), theta = cos^-1(z/r), phi = tan^-1(y/x)
        x_sphere = [(x(1)**2 + x(2)**2 + x(3)**2)**.5, acos(x(3)/((x(1)**2 + x(2)**2 + x(3)**2)**.5)), atan2(x(2),x(1))]

        ! Calculate B in spherical coordinates

        B_sphere = 4.176 * 1d6 * [2*cos(x_sphere(2)), sin(x_sphere(2)), 0.0d0] / (x_sphere(1)**3.0d0) ! x_sphere ~10^9 x_sphere^3 ~10^27

        ! Convert B(r,theta,phi) to B(x,y,z)
        ! Vector Transformation:
        ! Ax = Ar*sin(theta)*cos(phi) + Atheta*cos(theta)*cos(phi) - Aphi*sin(phi)
        ! Ay = Ar*sin(theta)*sin(phi) + Atheta*cos(theta)*sin(phi) + Aphi*cos(phi)
        ! Az = Ar*cos(theta)          - Atheta*sin(theta)          +      0
        B = [B_sphere(1)*sin(x_sphere(2))*cos(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*cos(x_sphere(3)) - B_sphere(3)*sin(x_sphere(3)), &
             B_sphere(1)*sin(x_sphere(2))*sin(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*sin(x_sphere(3)) + B_sphere(3)*cos(x_sphere(3)), &
             B_sphere(1)*cos(x_sphere(2)) - B_sphere(2)*sin(x_sphere(2))]


      case(3) ! read from dipole csv
        E = [0.0d0, 0.0d0, 0.0d0]

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! convert x -> x_spherical
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        x_sphere = [(x(1)**2 + x(2)**2 + x(3)**2)**.5, acos(x(3)/((x(1)**2 + x(2)**2 + x(3)**2)**.5)), atan2(x(2),x(1))]

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get index of the closest position value
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! loop through the r_vec

        do i = 1, 22
          ! calculate the absolute value of the difference between r_vec[i] and x_spherical(1)
          ! save the difference in r_delta array
          r_delta(i) = ABS(r_vec(i)-x_sphere(1))
        end do

        ! find the minimun value of r_delta array
        r_min = minval(r_delta)

        ! loop through the values of r_delta, and get the index
        do i = 1, 22
          ! if the minimum value is equal to r_delta[i]
          if (r_delta(i) == r_min) then
            ! set r_index = i
            r_index = i
          end if
        end do
              

        ! loop through the theta_vec
        do j = 1,22
          ! calculate the absolute difference between theta_vec[j] and x_spherical(2)
          ! save the difference in theta_delta array
          theta_delta(j) = ABS(theta_vec(j)-x_sphere(2))
        end do

        ! find the minimun value of theta_delta array
        theta_min = minval(theta_delta)

        ! loop through the values of theta_delta, and get the index
        do j= 1,22
          ! if the minimum value is equal to theta_delta[j]
          if (theta_delta(j) == theta_min) then
            ! set theta_index = j
            theta_index = j
          end if
        end do
   
        ! loop through the phi_vec
        do k = 1,22
          ! calculate the difference between phi_vec[k] and x_spherical(3)
          ! save the difference in phi_delta array
          phi_delta(k) = ABS(phi_vec(k)-x_sphere(3))
        end do

        ! find the minimun value of phi_delta array
        phi_min = minval(phi_delta)

        ! loop through the values of phi_delta, and get the index
        do k = 1,22
          ! if the minimum value is equal to phi_delta[k]
          if (phi_delta(k) == phi_min) then
            ! set phi_index = k
            phi_index = k
          end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get B value at that position from csv!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! first find the ith directory in bfield_grid_directory
        ! ith_dir = 'jrm09_spherical_grid/'+str(i)
        ! format ith num
        fmt = '(I2.2)' ! integer of width 2 with zeros at the left
        r_mult = mult_vec(r_index)
        write(r_mult_str, fmt) r_mult

        ! open ith_dir/Br.csv
        Br_file3 = 'dipole_grid/'//trim(r_mult_str)//'_rj_surface/Br.csv'
        open(2, file = Br_file3)
        ! Read up to the j-1th row of the csv
        ! Read the jth row of the Br_file3 into Br_tmp array
        do i_ = 1, (theta_index-1)
          read(2, *)
        end do 
        read(2, *) Br_tmp
        ! Set B_spherical(1) equal to the kth term in the Br_tmp array
        B_sphere(1) = Br_tmp(phi_index)
        close(2)

        ! open ith_dir/Btheta.csv
        Btheta_file3 = 'dipole_grid/'//trim(r_mult_str)//'_rj_surface/Btheta.csv'
        open(3, file = Btheta_file3)
        ! Read up to the j-1th row of the csv
        ! Read the jth row into the Btheta_tmp array 
        do j_ = 1, (theta_index-1)
          read(3, *)
        end do
        read(3, *) Btheta_tmp
        ! Set B_spherical(2) equal to the kth term in the Btheta_tmp array
        B_sphere(2) = Btheta_tmp(phi_index)
        close(3)


        ! open ith_dir/Bphi.csv
        Bphi_file3 = 'dipole_grid/'//trim(r_mult_str)//'_rj_surface/Bphi.csv'
        open(4, file = Bphi_file3)
        ! Read up to the j-1th row of the csv
        ! Read the jth row in the Bphi_tmp array
        do k_ = 1, (theta_index-1)
          read(4, *)
        end do
        ! Set B_spherical(3) equal to the kth term in the Bphi_tmp array
        B_sphere(3) = Bphi_tmp(phi_index)
        close(4)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Convert B_spherical to B 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Return B in cartesian coordinates!
        ! Convert B(r,theta,phi) to B(x,y,z)
        ! Vector Transformation:
        ! Ax = Ar*sin(theta)*cos(phi) + Atheta*cos(theta)*cos(phi) - Aphi*sin(phi)
        ! Ay = Ar*sin(theta)*sin(phi) + Atheta*cos(theta)*sin(phi) + Aphi*cos(phi)
        ! Az = Ar*cos(theta)          - Atheta*sin(theta)          +      0
        ! B[G] => B[G] * 10^4 T/m * 10^2 cm/m
        B = 1d6 * [B_sphere(1)*sin(x_sphere(2))*cos(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*cos(x_sphere(3)) - B_sphere(3)*sin(x_sphere(3)), &
             B_sphere(1)*sin(x_sphere(2))*sin(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*sin(x_sphere(3)) + B_sphere(3)*cos(x_sphere(3)), &
             B_sphere(1)*cos(x_sphere(2)) - B_sphere(2)*sin(x_sphere(2))]

      case(4) !read jrm09 from csv

        E = [0.0d0, 0.0d0, 0.0d0]

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! convert x -> x_spherical
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        x_sphere = [(x(1)**2 + x(2)**2 + x(3)**2)**.5, acos(x(3)/((x(1)**2 + x(2)**2 + x(3)**2)**.5)), atan2(x(2),x(1))]

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get index of the closest position value
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! loop through the r_vec

        do i = 1, 22
          ! calculate the absolute value of the difference between r_vec[i] and x_spherical(1)
          ! save the difference in r_delta array
          r_delta(i) = ABS(r_vec(i)-x_sphere(1))
        end do

        ! find the minimun value of r_delta array
        r_min = minval(r_delta)

        ! loop through the values of r_delta, and get the index
        do i = 1, 22
          ! if the minimum value is equal to r_delta[i]
          if (r_delta(i) == r_min) then
            ! set r_index = i
            r_index = i
          end if
        end do
              

        ! loop through the theta_vec
        do j = 1,22
          ! calculate the absolute difference between theta_vec[j] and x_spherical(2)
          ! save the difference in theta_delta array
          theta_delta(j) = ABS(theta_vec(j)-x_sphere(2))
        end do

        ! find the minimun value of theta_delta array
        theta_min = minval(theta_delta)

        ! loop through the values of theta_delta, and get the index
        do j= 1,22
          ! if the minimum value is equal to theta_delta[j]
          if (theta_delta(j) == theta_min) then
            ! set theta_index = j
            theta_index = j
          end if
        end do
   
        ! loop through the phi_vec
        do k = 1,22
          ! calculate the difference between phi_vec[k] and x_spherical(3)
          ! save the difference in phi_delta array
          phi_delta(k) = ABS(phi_vec(k)-x_sphere(3))
        end do

        ! find the minimun value of phi_delta array
        phi_min = minval(phi_delta)

        ! loop through the values of phi_delta, and get the index
        do k = 1,22
          ! if the minimum value is equal to phi_delta[k]
          if (phi_delta(k) == phi_min) then
            ! set phi_index = k
            phi_index = k
          end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get B value at that position from csv!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! first find the ith directory in bfield_grid_directory
        ! ith_dir = 'jrm09_spherical_grid/'+str(i)
        ! format ith num
        fmt = '(I2.2)' ! integer of width 2 with zeros at the left
        r_mult = mult_vec(r_index)
        write(r_mult_str, fmt) r_mult

        ! open ith_dir/Br.csv
        Br_file4 = 'jrm09_grid/'//trim(r_mult_str)//'_rj_surface/Br.csv'
        open(2, file = Br_file4)
        ! Read up to the j-1th row of the csv
        ! Read the jth row of the Br_file into Br_tmp array
        do i_ = 1, (theta_index-1)
          read(2, *)
        end do 
        read(2, *) Br_tmp
        ! Set B_spherical(1) equal to the kth term in the Br_tmp array
        B_sphere(1) = Br_tmp(phi_index)
        close(2)

        ! open ith_dir/Btheta.csv
        Btheta_file4 = 'jrm09_grid/'//trim(r_mult_str)//'_rj_surface/Btheta.csv'
        open(3, file = Btheta_file4)
        ! Read up to the j-1th row of the csv
        ! Read the jth row into the Btheta_tmp array 
        do j_ = 1, (theta_index-1)
          read(3, *)
        end do
        read(3, *) Btheta_tmp
        ! Set B_spherical(2) equal to the kth term in the Btheta_tmp array
        B_sphere(2) = Btheta_tmp(phi_index)
        close(3)


        ! open ith_dir/Bphi.csv
        Bphi_file4 = 'jrm09_grid/'//trim(r_mult_str)//'_rj_surface/Bphi.csv'
        open(4, file = Bphi_file4)
        ! Read up to the j-1th row of the csv
        ! Read the jth row in the Bphi_tmp array
        do k_ = 1, (theta_index-1)
          read(4, *)
        end do
        ! Set B_spherical(3) equal to the kth term in the Bphi_tmp array
        B_sphere(3) = Bphi_tmp(phi_index)
        close(4)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Convert B_spherical to B 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Return B in cartesian coordinates!
        ! Convert B(r,theta,phi) to B(x,y,z)
        ! Vector Transformation:
        ! Ax = Ar*sin(theta)*cos(phi) + Atheta*cos(theta)*cos(phi) - Aphi*sin(phi)
        ! Ay = Ar*sin(theta)*sin(phi) + Atheta*cos(theta)*sin(phi) + Aphi*cos(phi)
        ! Az = Ar*cos(theta)          - Atheta*sin(theta)          +      0
        ! B[G] => B[G] * 10^4 T/m * 10^2 cm/m
        B = 1d6 * [B_sphere(1)*sin(x_sphere(2))*cos(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*cos(x_sphere(3)) - B_sphere(3)*sin(x_sphere(3)), &
             B_sphere(1)*sin(x_sphere(2))*sin(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*sin(x_sphere(3)) + B_sphere(3)*cos(x_sphere(3)), &
             B_sphere(1)*cos(x_sphere(2)) - B_sphere(2)*sin(x_sphere(2))]

      case default
        call mpistop("Unknown value for iprob")
      end select

      ! The user can override these fields
      if (force_E0(1) > not_used_value) E = force_E0
      if (force_B0(1) > not_used_value) B = force_B0

    end subroutine get_field

    ! Set particle properties (in SI units)
    subroutine get_particle(x, v, q, m, ipart, n_particles, iprob)
      double precision, intent(out) :: x(3), v(3), q, m
      integer, intent(in)           :: ipart, iprob, n_particles
      double precision              :: tmp_vec(4), phi !ipart is the particle number n_particle of n_particles

      q = charge
      m = mass

      x = x0

      select case (iprob)
      
      !case (1) ! Dipole case
      !   v = v0
      !   q = (charge * ipart) / n_particles
      !
      !   ! Distribute over circle, velocity inwards. Avoid pi/4.
      !   phi = ((ipart+0.125d0) * 2 * acos(-1.0d0)) / n_particles
      !   x = norm2(x0) * [cos(phi), sin(phi), 0.0d0]

         ! Add Maxwellian velocity. Random numbers come in pairs of two
      !   tmp_vec(1:2) = rng%two_normals()
      !   tmp_vec(3:4) = rng%two_normals()
      !   v = v0 + tmp_vec(1:3) * maxwellian_velocity
      !case (2) ! Dipole case
      !   v = v0
      !   q = (charge * ipart) / n_particles
      !   if (physics_type_particles /= 'gca') then
      !      ! Assume B = 10 T, and v_x = 0 initially
      !      ! x = vx + |vy|* m/ M * q -> offset probably to put on orbit
      !      x(1) = x(1) + abs(v(2)) * m / (q * 10.0d0)
      !
      !   end if
      !case (3) ! jrm09
      !   v = v0
      !   q = (charge * ipart) / n_particles
      !   ! Assume B = 10 T, and v_x = 0 initially
      !   ! x = vx + |vy|*m/10q
      !   !x(1) = x(1) + abs(v(2)) * m / (q * 10.0d0)
      !
      !   ! Distribute over circle, velocity inwards. Avoid pi/4.
      !   phi = ((ipart+0.125d0) * 2 * acos(-1.0d0)) / n_particles
      !   x = norm2(x0) * [cos(phi), sin(phi), 0.0d0]
      !
      !   ! Add Maxwellian velocity. Random numbers come in pairs of two
      !   tmp_vec(1:2) = rng%two_normals()
      !   tmp_vec(3:4) = rng%two_normals()
      !   v = v0 + tmp_vec(1:3) * maxwellian_velocity
      !case (4)
      !   ! Add Maxwellian velocity. Random numbers come in pairs of two
      !   tmp_vec(1:2) = rng%two_normals()
      !   tmp_vec(3:4) = rng%two_normals()
      !   v = v0 + tmp_vec(1:3) * maxwellian_velocity
      !   print *, ipart, v
      !case (5)
      !   v = (v0 * ipart) / n_particles
      !case (8)
      !   ! Distribute over circle, velocity inwards. Avoid pi/4.
      !   phi = ((ipart+0.125d0) * 2 * acos(-1.0d0)) / n_particles
      !   x = norm2(x0) * [cos(phi), sin(phi), 0.0d0]
      !   v = -x * norm2(v0)/norm2(x0)
      case default
         !v = v0
         !q = (charge * ipart) / n_particles

         ! Distribute over circle, velocity inwards. Avoid pi/4.
         phi = ((ipart+0.125d0) * 2 * acos(-1.0d0)) / n_particles
         x = norm2(x0) * [cos(phi), sin(phi), 0.0d0]

         ! Add Maxwellian velocity. Random numbers come in pairs of two
         tmp_vec(1:2) = rng%two_normals()
         tmp_vec(3:4) = rng%two_normals()
         v = v0 + tmp_vec(1:3) * maxwellian_velocity
      end select
    end subroutine get_particle

    subroutine set_custom_field(w, x, E_field, B_field)
      double precision, intent(in)  :: w(ixG^T,nw)
      double precision, intent(in)  :: x(ixG^T,ndim)
      double precision, intent(out) :: E_field(ixG^T, ndir)
      double precision, intent(out) :: B_field(ixG^T, ndir)

      integer          :: i^D
      double precision :: E(3), B(3), xtmp(ndim)


      {do i^D = ixGlo^D, ixGhi^D\}
      xtmp = x(i^D, :)
      call get_field(xtmp, E, B)

      ! Convert to CGS units, 1 T -> 1e4 Gauss
      B_field(i^D, :) = B * 1.0d4

      ! Convert to CGS units
      E_field(i^D, :) = E * 1.0d6/const_c
      {end do\}
    !print*,'E in cgs:',E(1) * 1.0d6/const_c, 'B in cgs',B(3) * 1.0d4
    !print*,'E/B ratio < 1',(E(1)*1.0d6/const_c)/(B(3)*1.0d4)
    end subroutine set_custom_field

    subroutine get_analytic_field(ix, x, tloc, vec)
      integer, intent(in)           :: ix(ndir) !< Indices in gridvars
      double precision, intent(in)  :: x(ndir)
      double precision, intent(in)  :: tloc
      double precision, intent(out) :: vec(ndir)
      double precision              :: E(3), B(3)

      ! Get the field analytically and put it on the grid

      call get_field(x, E, B)

      if (ix(1) == bp(1)) then
        vec = B * 1.0d4
      else if (ix(1) == ep(1)) then
        vec = E * 1.0d6/const_c
      else
        call mpistop("get_analytic_field: unknown variable index")
      end if
    end subroutine get_analytic_field

    subroutine MakeMagGridDH(cilm, lmax, r0, a, f, rad_grid, theta_grid, &
                             phi_grid, total_grid, n, sampling, lmax_calc, &
                             pot_grid, exitstatus)
    !------------------------------------------------------------------------------
    !
    !   Given the Schmidt semi-normalized magnetic potential spherical harmonic
    !   coefficients CILM, this subroutine will compute a 2D Driscoll and Healy
    !   sampled grid of the three components and magnitude field. The output grids
    !   are cacluated on a flattened ellipsoid with semi-major axis A and flattening
    !   F. The output grids contain N samples in latitude and longitude by default,
    !   but if the optional parameter SAMPLING is set to 2, the grids will contain
    !   N samples in latitude and 2N samples in longitude. In order to calculate
    !   the entire gravitational acceleration, it is necessary that the degree-0
    !   term be set equal to 1. The radial magnetic field component is assumed to
    !   be positive when directed UPWARDS. The magnetic potential is defined
    !   according to
    !
    !       V = R0 Sum_{l=1}^LMAX (R0/r)^{l+1} Sum_{m=-l}^l C_{lm} Y_{lm}.
    !
    !   and
    !
    !   band for 90 S is not calculated, and the latitudinal sampling interval is
    !       B = - Grad V
    !
    !   The first latitudinal band of the grid corresponds to 90 N, the latitudinal
    !   180/N degrees. The first longitudinal band is 0 E, the longitudinal band for
    !   360 E is not calculated, and the longitudinal sampling interval is 360/N for
    !   equally sampled and 180/N for equally spaced grids, respectively.
    !
    !   Calling Parameters
    !
    !       IN
    !           cilm        The Schmidt semi-normalized magnetic potential
    !                       spherical harmonic coefficients.
    !           lmax        The maximum spherical harmonic degree of the function,
    !                       used to determine the number of samples N.
    !           r0          Reference radius of potential coefficients.
    !           a           The semimajor axis of the flattened ellipsoid.
    !           f           Flattening of the planet.
    !
    !       IN, OPTIONAL
    !           sampling    (1) Grid is N latitudes by N longitudes (default).
    !                       (2) Grid is N by 2N. The higher frequencies resulting
    !                       from this oversampling in longitude are discarded, and
    !                       hence not aliased into lower frequencies.
    !           lmax_calc   The maximum spherical harmonic degree to evaluate
    !                       the coefficients up to.
    !
    !       OUT
    !           rad_grid    Gridded expansion of the radial component of the
    !                       magnetic field.
    !           theta_grid  Gridded expansaion of the theta component of the
    !                       gravitational field.
    !           phi_grid    Gridded expansaion of the phi component of the
    !                       gravitational field.
    !           total_grid  Gridded expansaion of the the magnitude of the
    !                       gravitational field.
    !           N           Number of samples in latitude. Number of samples in
    !                       longitude is N when sampling is 1 (default), and is 2N
    !                       when sampling is 2.
    !
    !       OUT, OPTIONAL
    !           pot_grid    Magnetic potential on the ellipsoid in SI units.
    !           exitstatus  If present, instead of executing a STOP when an error
    !                       is encountered, the variable exitstatus will be
    !                       returned describing the error.
    !                       0 = No errors;
    !                       1 = Improper dimensions of input array;
    !                       2 = Improper bounds for input variable;
    !                       3 = Error allocating memory;
    !                       4 = File IO error.
    !
    !   Notes:
    !       1.  If lmax is greater than the the maximum spherical harmonic
    !           degree of the input file, Cilm will be ZERO PADDED!
    !           (i.e., those degrees after lmax are assumed to be zero).
    !       2.  Latitude is geocentric latitude.
    !
    !   Dependencies:   FFTW3
    !
    !   Copyright (c) 2016, SHTOOLS
    !   All rights reserved.
    !
    !-------------------------------------------------------------------------------
        use FFTW3 !module FFTW3.f95 also needs to be included, but since it just names parameters I can put it at the beginning of mod_usr.t
#ifdef FFTW3_UNDERSCORE
#define dfftw_plan_dft_c2r_1d dfftw_plan_dft_c2r_1d_
#define dfftw_execute dfftw_execute_
#define dfftw_destroy_plan dfftw_destroy_plan_
#endif
        implicit none

        real*8, intent(in) :: cilm(:,:,:), r0, a, f 
        !cilm is a 3d array (gh, l, m) of dimension (2, *, *), l[1, lmax], m[0,l]
        real*8, intent(out) :: rad_grid(:,:), theta_grid(:,:), phi_grid(:,:), &
                               total_grid(:,:)
        real*8, intent(out), optional :: pot_grid(:,:)
        integer, intent(in) ::  lmax !highest order the coefficients will be calculated to
        integer, intent(out) :: n
        integer, intent(in), optional :: sampling, lmax_calc
        integer, intent(out), optional :: exitstatus
        integer ::      l, m, i, l1, m1, lmax_comp, i_eq, i_s, astat(4), nlong
        real*8 ::       grid(4*lmax+4), pi, theta, scalef, rescalem, u, p, dp, &
                        pmm, pm1, pm2, z, tempr, r_ex, lat, prefactor(lmax), &
                        coefr0, coefu0, coefrs0, coeft0, coefts0, coefp0, coefps0, &
                        coefus0
        complex*16 :: coef(2*lmax+3), coefr(2*lmax+3), coefrs(2*lmax+3), &
                      coeft(2*lmax+3), coefts(2*lmax+3), coefp(2*lmax+3), &
                      coefps(2*lmax+3), coefu(2*lmax+3), coefus(2*lmax+3), tempc
        integer*8 :: plan
        real*8, save, allocatable :: ff1(:,:), ff2(:,:), sqr(:)
        integer*1, save, allocatable :: fsymsign(:,:)
        integer, save ::  lmax_old = 0
        logical :: calcu
        external :: dfftw_plan_dft_c2r_1d, dfftw_execute, dfftw_destroy_plan

    !$OMP   threadprivate(ff1, ff2, sqr, fsymsign, lmax_old)

        if (present(exitstatus)) exitstatus = 0

        n = 2 * lmax + 2

        if (present(sampling)) then
            if (sampling /= 1 .and. sampling /=2) then
                print*, "Error --- MakeMagGridDH"
                print*, "Optional parameter SAMPLING must be 1 (N by N) " // &
                        "or 2 (N by 2N)."
                print*, "Input value is ", sampling
                if (present(exitstatus)) then
                    exitstatus = 2
                    return
                else
                    stop
                end if

            endif
        end if

        if (size(cilm(:,1,1)) < 2) then
            print*, "Error --- MakeMagGridDH"
            print*, "CILM must be dimensioned as (2, *, *)."
            print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
                                           size(cilm(1,1,:))
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

        if (present(sampling)) then
            if (sampling == 1) then
                nlong = n
            else
                nlong = 2 * n
            end if

        else
            nlong = n

        end if

        if (size(rad_grid(:,1)) < n .or. size(rad_grid(1,:)) < nlong .or. &
                size(theta_grid(:,1)) < n .or. size(theta_grid(1,:)) < nlong &
                .or. size(phi_grid(:,1)) < n .or. size(phi_grid(1,:)) < nlong .or. &
                size(total_grid(:,1)) < n .or. size(total_grid(1,:)) < nlong) then
            print*, "Error --- MakeMagGridDH"
            if (present(sampling)) then
                if (sampling == 1) then
                    print*, "RAD_GRID, THETA_GRID, PHI_GRID, and TOTAL_GRID " // &
                            "must be dimensioned as (N, N) where N is ", n
                else if (sampling == 2) then
                    print*, "RAD_GRID, THETA_GRID, PHI_GRID, and TOTAL_GRID " // &
                            "must be dimensioned as (N, 2N) where N is ", n
                end if
            else
                print*, "RAD_GRID, THETA_GRID, PHI_GRID, and TOTAL_GRID " // &
                        "must be dimensioned as (N, N) where N is ", n
            end if

            print*, "Input dimensions are ", size(rad_grid(:,1)), &
                    size(rad_grid(1,:)), size(theta_grid(:,1)), &
                    size(theta_grid(1,:)), size(phi_grid(:,1)), &
                    size(phi_grid(1,:)), size(total_grid(:,1)), &
                    size(total_grid(1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

        if (present(pot_grid)) then
            calcu = .true.
            if (size(pot_grid(:,1)) < n .or. size(pot_grid(1,:)) < nlong) then
                print*, "Error --- MakeMagGridDH"
                if (present(sampling)) then
                    if (sampling == 1) then
                        print*, "POT_GRID must be dimensioned as (N, N) " // &
                                "where N is ", n
                    else if (sampling == 2) then
                        print*, "POT_GRID must be dimensioned as (N, 2N) " // & 
                                "where N is ", n
                    end if
                else
                    print*, "POT_GRID must be dimensioned as (N, N) where N is ", n
                end if

                print*, "Input dimensions are ", size(pot_grid(:,1)), &
                                                 size(pot_grid(1,:))
                if (present(exitstatus)) then
                    exitstatus = 1
                    return
                else
                    stop
                end if

            end if

        else 
            calcu = .false.

        end if

        pi = acos(-1.0d0)

        scalef = 1.0d-280

        if (present(lmax_calc)) then
            if (lmax_calc > lmax) then
                print*, "Error --- MakeMagGridDH"
                print*, "LMAX_CALC must be less than or equal to LMAX."
                print*, "LMAX = ", lmax
                print*, "LMAX_CALC = ", lmax_calc
                if (present(exitstatus)) then
                    exitstatus = 2
                    return
                else
                    stop
                end if

            else
                lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1, &
                                lmax_calc)

            end if

        else
            lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1)

        end if

        if (lmax_comp == 0) then
            rad_grid(1:n,1:nlong) = 0.0d0
            theta_grid(1:n,1:nlong) = 0.0d0
            phi_grid(1:n,1:nlong) = 0.0d0
            total_grid(1:n,1:nlong) = 0.0d0
            if (calcu) pot_grid(1:n,1:nlong) = 0.0d0
            return
        endif

        !--------------------------------------------------------------------------
        !
        !   Calculate recursion constants used in computing Legendre polynomials
        !
        !--------------------------------------------------------------------------
        if (lmax_comp /= lmax_old) then

            if (allocated(sqr)) deallocate(sqr)
            if (allocated(ff1)) deallocate(ff1)
            if (allocated(ff2)) deallocate(ff2)
            if (allocated(fsymsign)) deallocate(fsymsign)

            allocate (sqr(2*lmax_comp+1), stat=astat(1))
            allocate (ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
            allocate (ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
            allocate (fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))

            if (sum(astat(1:4)) /= 0) then
                print*, "Error --- MakeMagGridDH"
                print*, "Problem allocating arrays SQR, FF1, FF2, or FSYMSIGN", &
                        astat(1), astat(2), astat(3), astat(4)
                if (present(exitstatus)) then
                    exitstatus = 3
                    return
                else
                    stop
                end if

            endif

            !----------------------------------------------------------------------
            !
            !   Calculate signs used for symmetry of Legendre functions about
            !   equator. For the first derivative in theta, these signs are
            !   reversed.
            !
            !----------------------------------------------------------------------
            do l = 0, lmax_comp, 1
                do m = 0, l, 1
                    if (mod(l-m,2) == 0) then
                        fsymsign(l+1,m+1) = 1

                    else
                        fsymsign(l+1,m+1) = -1

                    end if

                end do

            end do

            !----------------------------------------------------------------------
            !
            !   Precompute square roots of integers that are used several times.
            !
            !----------------------------------------------------------------------
            do l = 1, 2 * lmax_comp+1
                sqr(l) = sqrt(dble(l))
            end do

            !----------------------------------------------------------------------
            !
            !   Precompute multiplicative factors used in recursion relationships
            !       P(l,m) = x*f1(l,m)*P(l-1,m) - P(l-2,m)*f2(l,m)
            !       k = l*(l+1)/2 + m + 1
            !   Note that prefactors are not used for the case when m=l as a
            !   different recursion is used. Furthermore, for m=l-1, Plmbar(l-2,m)
            !   is assumed to be zero.
            !
            !----------------------------------------------------------------------
            if (lmax_comp /= 0) then
                ff1(2,1) = 1.0d0
                ff2(2,1) = 0.0d0

            endif

            do l = 2, lmax_comp, 1
                ff1(l+1,1) = dble(2*l-1) / dble(l)
                ff2(l+1,1) = dble(l-1) / dble(l)

                do m = 1, l-2, 1
                    ff1(l+1,m+1) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                    ff2(l+1,m+1) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) / sqr(l-m)
                end do

                ff1(l+1,l)= dble(2*l-1) / sqr(l+m) / sqr(l-m)
                ff2(l+1,l) = 0.0d0

            end do

            lmax_old = lmax_comp

        end if

        !--------------------------------------------------------------------------
        !
        !   Determine Clms one l at a time by intergrating over latitude.
        !
        !--------------------------------------------------------------------------
        call dfftw_plan_dft_c2r_1d(plan, nlong, coef(1:nlong/2+1), grid(1:nlong), &
                                   FFTW_MEASURE)

        i_eq = n / 2 + 1  ! Index correspondong to zero latitude

        do i = 1, i_eq - 1, 1

            i_s = 2 * i_eq - i

            theta = pi * dble(i-1) / dble(n)
            z = cos(theta)
            u = sqrt( (1.0d0-z) * (1.0d0+z) )

            lat = pi / 2.0d0 - theta

            if (i==1) then      ! Reference ellipsoid radius
                r_ex = a * (1.0d0 - f)

            else
                r_ex = (1.0d0 + tan(lat)**2) &
                       / (1.0d0  + tan(lat)**2 / (1.0d0 - f)**2)
                r_ex = a * sqrt(r_ex)

            end if

            coefr(1:lmax+2) = dcmplx(0.0d0,0.0d0)
            coefr0 = 0.0d0
            coefrs(1:lmax+2) = dcmplx(0.0d0,0.0d0)
            coefrs0 = 0.0d0

            coeft(1:lmax+2) = dcmplx(0.0d0,0.0d0)
            coeft0 = 0.0d0
            coefts(1:lmax+2) = dcmplx(0.0d0,0.0d0)
            coefts0 = 0.0d0

            coefp(1:lmax+2) = dcmplx(0.0d0,0.0d0)
            coefp0 = 0.0d0
            coefps(1:lmax+2) = dcmplx(0.0d0,0.0d0)
            coefps0 = 0.0d0

            if (calcu) then
                coefu(1:lmax+2) = dcmplx(0.0d0,0.0d0)
                coefu0 = 0.0d0
                coefus(1:lmax+2) = dcmplx(0.0d0,0.0d0)
                coefus0 = 0.0d0
            end if

            pm2 = 1.0d0

            ! l = 0 terms are zero
            prefactor(1) = (r0 / r_ex)**2

            do l = 2, lmax_comp, 1
                prefactor(l) = prefactor(l-1) * r0 / r_ex
            end do

            pm1 = ff1(2,1) * z
            tempr = cilm(1,2,1) * pm1 * (-2) * prefactor(1) ! -2 = (l+1) prefactor
            coefr0 = coefr0 + tempr
            coefrs0 = coefrs0 - tempr   ! fsymsign = -1

            if (calcu) then
                tempr = cilm(1,2,1) * pm1 * prefactor(1)
                coefu0 = coefu0 + tempr
                coefus0 = coefus0 - tempr   ! fsymsign = -1
            end if

            dp = ff1(2,1)
            tempr = cilm(1,2,1) * dp * prefactor(1)
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 + tempr   ! reverse fsymsign

            do l = 2, lmax_comp, 1
                l1 = l + 1
                p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
                tempr = cilm(1,l1,1) * p * (-l1) * prefactor(l)
                coefr0 = coefr0 + tempr
                coefrs0 = coefrs0 + tempr * fsymsign(l1,1)

                if (calcu) then
                    tempr = cilm(1,l1,1) * p * prefactor(l)
                    coefu0 = coefu0 + tempr
                    coefus0 = coefus0 + tempr * fsymsign(l1,1)
                end if

                dp = l * (pm1 - z * p) / u**2

                tempr = cilm(1,l1,1) * dp * prefactor(l)
                coeft0 = coeft0 + tempr
                coefts0 = coefts0 - tempr * fsymsign(l1,1)  ! reverse fsymsign

                pm2 = pm1
                pm1 = p

            end do

            pmm = sqr(2) * scalef

            rescalem = 1.0d0 / scalef

            do m = 1, lmax_comp-1, 1

                m1 = m + 1
                rescalem = rescalem * u

                pmm = pmm * sqr(2*m+1) / sqr(2*m)
                pm2 = pmm / sqr(2*m+1)

                tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * (-m-1) &
                        * prefactor(m)    ! (m,m)
                coefr(m1) = coefr(m1) + tempc
                coefrs(m1) = coefrs(m1) + tempc
                ! fsymsign = 1

                if (calcu) then
                    tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 &
                            * prefactor(m) ! (m,m)
                    coefu(m1) = coefu(m1) + tempc
                    coefus(m1) = coefus(m1) + tempc
                    ! fsymsign = 1
                end if

                tempc = dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 &
                        * prefactor(m) * m ! (m,m)
                coefp(m1) = coefp(m1) + tempc
                coefps(m1) = coefps(m1) + tempc
                ! fsymsign = 1

                dp = -m * z * pm2 / u**2
                tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * dp &
                        * prefactor(m)  ! (m,m)
                coeft(m1) = coeft(m1) + tempc
                coefts(m1) = coefts(m1) - tempc ! reverse fsymsign

                pm1 = z * ff1(m1+1,m1) * pm2
                tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * pm1 * (-m-2) &
                        * prefactor(m+1)  ! (m+1,m)
                coefr(m1) = coefr(m1) + tempc
                coefrs(m1) = coefrs(m1) - tempc
                ! fsymsign = -1

                if (calcu) then
                    tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * pm1 &
                            * prefactor(m+1)   ! (m+1,m)
                    coefu(m1) = coefu(m1) + tempc
                    coefus(m1) = coefus(m1) - tempc
                end if

                tempc = dcmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1)) * pm1 &
                        * prefactor(m+1) * m ! (m+1,m)
                coefp(m1) = coefp(m1) + tempc
                coefps(m1) = coefps(m1) - tempc
                ! fsymsign = -1

                dp = (pm2 * sqr(2*m+1) - z * (m+1) * pm1) / u**2

                tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * dp &
                        * prefactor(m+1)    ! (m+1,m)
                coeft(m1) = coeft(m1) + tempc
                coefts(m1) = coefts(m1) + tempc ! reverse fsymsign

                do l = m+2, lmax_comp, 1
                    l1 = l + 1
                    p = z * ff1(l1,m1) * pm1 - ff2(l1,m1) * pm2
                    tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * (-l1) &
                            * prefactor(l)
                    coefr(m1) = coefr(m1) + tempc
                    coefrs(m1) = coefrs(m1) + tempc * fsymsign(l1,m1)

                    if (calcu) then
                        tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p &
                                * prefactor(l)
                        coefu(m1) = coefu(m1) + tempc
                        coefus(m1) = coefus(m1) + tempc * fsymsign(l1,m1)
                    end if

                    tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * p &
                            * prefactor(l) * m
                    coefp(m1) = coefp(m1) + tempc
                    coefps(m1) = coefps(m1) + tempc * fsymsign(l1,m1)

                    dp = ( sqr(l+m) * sqr(l-m) * pm1 - l * z * p ) / u**2
                    tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp &
                            * prefactor(l)
                    coeft(m1) = coeft(m1) + tempc

                    ! reverse fsymsign
                    coefts(m1) = coefts(m1) - tempc * fsymsign(l1,m1)

                    pm2 = pm1
                    pm1 = p

                end do

                coefr(m1) = coefr(m1) * rescalem
                coefrs(m1) = coefrs(m1) * rescalem

                if (calcu) then
                    coefu(m1) = coefu(m1) * rescalem
                    coefus(m1) = coefus(m1) * rescalem
                end if

                coeft(m1) = coeft(m1) * rescalem
                coefts(m1) = coefts(m1) * rescalem

                coefp(m1) = coefp(m1) * rescalem
                coefps(m1) = coefps(m1) * rescalem

            end do

            rescalem = rescalem * u

            pmm = pmm / sqr(2*lmax_comp) * rescalem
            tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - &
                           cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                           * (-lmax_comp-1) * prefactor(lmax_comp)
            coefr(lmax_comp+1) = coefr(lmax_comp+1) + tempc
            coefrs(lmax_comp+1) = coefrs(lmax_comp+1) + tempc
            ! fsymsign = 1

            if (calcu) then
                tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - &
                               cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                               * prefactor(lmax_comp)
                coefu(lmax_comp+1) = coefu(lmax_comp+1) + tempc
                coefus(lmax_comp+1) = coefus(lmax_comp+1) + tempc
            end if

            tempc = dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                           cilm(1,lmax_comp+1,lmax_comp+1)) * pmm &
                           * prefactor(lmax_comp) * lmax_comp
            coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc
            coefps(lmax_comp+1) = coefps(lmax_comp+1) + tempc
            ! fsymsign = 1

            dp = -lmax_comp * z * pmm / u**2
            tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - &
                           cilm(2,lmax_comp+1,lmax_comp+1)) * dp &
                           * prefactor(lmax_comp)
            coeft(lmax_comp+1) = coeft(lmax_comp+1) + tempc
            coefts(lmax_comp+1) = coefts(lmax_comp+1) - tempc   ! reverse fsymsign

            coef(1) = dcmplx(coefr0,0.0d0)
            coef(2:lmax+1) = coefr(2:lmax+1) / 2.0d0

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                end if
            end if

            call dfftw_execute(plan)   ! take fourier transform
            rad_grid(i,1:nlong) = - grid(1:nlong) * r0 / r_ex

            if (calcu) then
                coef(1) = dcmplx(coefu0,0.0d0)
                coef(2:lmax+1) = coefu(2:lmax+1) / 2.0d0

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                    end if
                end if

                call dfftw_execute(plan)   ! take fourier transform
                pot_grid(i,1:nlong) = grid(1:nlong) * r0

            end if

            if (i==1) then
                ! These two derivatives are undefined at the pole
                theta_grid(1,1:nlong) = 0.0d0
                phi_grid(1,1:nlong) = 0.0d0

            else
                coef(1) = dcmplx(coeft0,0.0d0)
                coef(2:lmax+1) = coeft(2:lmax+1) / 2.0d0

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                    end if
                end if

                call dfftw_execute(plan)   ! take fourier transform
                theta_grid(i,1:nlong) = sin(theta) * grid(1:nlong) * r0 / r_ex

                coef(1) = dcmplx(coefp0,0.0d0)
                coef(2:lmax+1) = coefp(2:lmax+1) / 2.0d0

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                    end if
                end if

                call dfftw_execute(plan)   ! take fourier transform
                    phi_grid(i,1:nlong) = - grid(1:nlong) * (r0/r_ex) / sin(theta)
                endif

                if (i /= 1) then    ! don't compute value for south pole.
                    coef(1) = dcmplx(coefrs0,0.0d0)
                    coef(2:lmax+1) = coefrs(2:lmax+1) / 2.0d0

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                    end if
                end if

                call dfftw_execute(plan)   ! take fourier transform
                rad_grid(i_s,1:nlong) = - grid(1:nlong) * r0 / r_ex

                if (calcu) then
                    coef(1) = dcmplx(coefus0,0.0d0)
                    coef(2:lmax+1) = coefus(2:lmax+1) / 2.0d0

                    if (present(sampling)) then
                        if (sampling == 2) then
                            coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                        end if
                    end if

                    call dfftw_execute(plan)   ! take fourier transform
                    pot_grid(i_s,1:nlong) = grid(1:nlong) * r0
                end if

                coef(1) = dcmplx(coefts0,0.0d0)
                coef(2:lmax+1) = coefts(2:lmax+1) / 2.0d0

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                    end if
                end if

                call dfftw_execute(plan)   ! take fourier transform
                theta_grid(i_s,1:nlong) = sin(theta)*grid(1:nlong) * r0 / r_ex

                coef(1) = dcmplx(coefps0,0.0d0)
                coef(2:lmax+1) = coefps(2:lmax+1) / 2.0d0

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                    end if
                end if

                call dfftw_execute(plan)   ! take fourier transform
                phi_grid(i_s,1:nlong) = - grid(1:nlong) * (r0 / r_ex) / sin(theta)

            end if

        end do

        ! Finally, do equator
        r_ex = a
        theta = pi / 2.0d0
        z = 0.0d0
        u = 1.0d0

        coefr(1:lmax+2) = dcmplx(0.0d0,0.0d0)
        coefr0 = 0.0d0

        if (calcu) then
            coefu(1:lmax+2) = dcmplx(0.0d0,0.0d0)
            coefu0 = 0.0d0
        end if

        coeft(1:lmax+2) = dcmplx(0.0d0,0.0d0)
        coeft0 = 0.0d0

        coefp(1:lmax+2) = dcmplx(0.0d0,0.0d0)
        coefp0 = 0.0d0

        pm2 = 1.0d0

        prefactor(1) = (r0 / r_ex)**2
        do l=2, lmax_comp,1
            prefactor(l) = prefactor(l-1) * r0 / r_ex
        end do

        pm1 =  0.0d0

        if (calcu) coefu0 = coefu0 + cilm(1,2,1) * pm1 * prefactor(1)

        dp = ff1(2,1)
        coeft0 = coeft0 + cilm(1,2,1) * dp * prefactor(1)

        do l = 2, lmax_comp, 1
            l1 = l + 1
            p = - ff2(l1,1) * pm2
            coefr0 = coefr0 + cilm(1,l1,1) * p * (-l1) * prefactor(l)

            if (calcu) coefu0 = coefu0 + cilm(1,l1,1) * p * prefactor(l)

            dp = l * (pm1)
            coeft0 = coeft0 + cilm(1,l1,1) * dp * prefactor(l)

            pm2 = pm1
            pm1 = p
        end do

        pmm = sqr(2) * scalef

        rescalem = 1.0d0 / scalef

        do m = 1, lmax_comp-1, 1
            m1 = m + 1

            pmm = pmm * sqr(2*m+1) / sqr(2*m)
            pm2 = pmm / sqr(2*m+1)

            coefr(m1) = coefr(m1) + dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 &
                        * (-m-1) * prefactor(m)
        
            if (calcu) coefu(m1) = coefu(m1) + dcmplx(cilm(1,m1,m1), - &
                                               cilm(2,m1,m1)) * pm2 * prefactor(m)

            coefp(m1) = coefp(m1) + dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 &
                                    * prefactor(m) * m

            pm1 = 0.0d0

            dp = (pm2 * sqr(2*m+1))
            coeft(m1) = coeft(m1) + dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) &
                                    * dp * prefactor(m+1)

            do l = m+2, lmax_comp, 1
                l1 = l + 1
                p = - ff2(l1,m1) * pm2
                coefr(m1) = coefr(m1) + dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) &
                                        * p * (-l1) * prefactor(l)

                if (calcu) coefu(m1) = coefu(m1) + dcmplx(cilm(1,l1,m1), - &
                                       cilm(2,l1,m1)) * p * prefactor(l)

                coefp(m1) = coefp(m1) + dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) &
                                        * p * prefactor(l) * m

                dp = ( sqr(l+m) * sqr(l-m) * pm1)
                coeft(m1) = coeft(m1) + dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) &
                                        * dp * prefactor(l)

                pm2 = pm1
                pm1 = p

            end do

            coefr(m1) = coefr(m1) * rescalem

            if (calcu) coefu(m1) = coefu(m1) * rescalem

            coefp(m1) = coefp(m1) * rescalem

            coeft(m1) = coeft(m1) * rescalem

        end do

        pmm = pmm / sqr(2*lmax_comp) * rescalem
        coefr(lmax_comp+1) = coefr(lmax_comp+1) + &
                             dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                             - cilm(2,lmax_comp+1,lmax_comp+1)) &
                             * pmm * (-lmax_comp-1) * prefactor(lmax_comp)

        if (calcu) coefu(lmax_comp+1) = coefu(lmax_comp+1) + &
                                        dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                                        - cilm(2,lmax_comp+1,lmax_comp+1)) &
                                        * pmm * prefactor(lmax_comp)

        coefp(lmax_comp+1) = coefp(lmax_comp+1) + &
                            dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1)) &
                            * pmm * prefactor(lmax_comp) * lmax_comp

        dp = -lmax_comp * z * pmm / u**2
        coeft(lmax_comp+1) = coeft(lmax_comp+1) &
                            + dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1)) &
                            * dp * prefactor(lmax_comp)

        coef(1) = dcmplx(coefr0,0.0d0)
        coef(2:lmax+1) = coefr(2:lmax+1) / 2.0d0

        if (present(sampling)) then
            if (sampling == 2) then
                coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
            end if
        end if

        call dfftw_execute(plan)   ! take fourier transform
        rad_grid(i_eq,1:nlong) = - grid(1:nlong) * r0 / r_ex

        if (calcu) then
            coef(1) = dcmplx(coefu0,0.0d0)
            coef(2:lmax+1) = coefu(2:lmax+1) / 2.0d0

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
                end if
            end if

            call dfftw_execute(plan)   ! take fourier transform
            pot_grid(i_eq,1:nlong) = grid(1:nlong) * r0

        end if

        coef(1) = dcmplx(coeft0,0.0d0)
        coef(2:lmax+1) = coeft(2:lmax+1) / 2.0d0

        if (present(sampling)) then
            if (sampling == 2) then
                coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
            end if
        end if

        call dfftw_execute(plan)   ! take fourier transform
        theta_grid(i_eq,1:nlong) = sin(theta) * grid(1:nlong) * r0 / r_ex

        coef(1) = dcmplx(coefp0,0.0d0)
        coef(2:lmax+1) = coefp(2:lmax+1) / 2.0d0

        if (present(sampling)) then
            if (sampling == 2) then
                coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)
            end if
        end if

        call dfftw_execute(plan)   ! take fourier transform
        phi_grid(i_eq,1:nlong) = - grid(1:nlong) * (r0/r_ex) / sin(theta)

        call dfftw_destroy_plan(plan)

        total_grid(1:n, 1:nlong) = sqrt(rad_grid(1:n,1:nlong)**2 &
                                        + phi_grid(1:n,1:nlong)**2 &
                                        + theta_grid(1:n,1:nlong)**2)

    end subroutine MakeMagGridDH
  end module mod_usr
