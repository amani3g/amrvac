module mod_usr
  use mod_mhd
  use mod_particles

  implicit none

  double precision :: charge = 1.60218d-19 ![C]
  double precision :: mass = 9.10938d-31 ![kg]

  ! Initial position (in m)
  double precision, parameter :: rj = 4.192d7
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

      B_sphere = 10 * 1d6 * [2*cos(x_sphere(3)), sin(x_sphere(3)), 0.0d0] / (x_sphere(1)**3.0d0) ! x_sphere ~10^9 x_sphere^3 ~10^27

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
      B = [B_sphere(1)*sin(x_sphere(2))*cos(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*cos(x_sphere(3)) - B_sphere(3)*sin(x_sphere(3)), &
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
      B = [B_sphere(1)*sin(x_sphere(2))*cos(x_sphere(3)) + B_sphere(2)*cos(x_sphere(2))*cos(x_sphere(3)) - B_sphere(3)*sin(x_sphere(3)), &
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
    double precision              :: tmp_vec(4), phi

    q = charge
    m = mass

    x = x0

    select case (iprob)
    
    case (1) ! Dipole case
       v = v0
       q = (charge * ipart) / n_particles
       ! Assume B = 10 T, and v_x = 0 initially
       ! x = vx + |vy|*m/10q
       !x(1) = x(1) + abs(v(2)) * m / (q * 10.0d0)

       ! Distribute over circle, velocity inwards. Avoid pi/4.
       phi = ((ipart+0.125d0) * 2 * acos(-1.0d0)) / n_particles
       x = norm2(x0) * [cos(phi), sin(phi), 0.0d0]

       ! Add Maxwellian velocity. Random numbers come in pairs of two
       tmp_vec(1:2) = rng%two_normals()
       tmp_vec(3:4) = rng%two_normals()
       v = v0 + tmp_vec(1:3) * maxwellian_velocity
    case (2) ! Dipole case
       v = v0
       q = (charge * ipart) / n_particles
       if (physics_type_particles /= 'gca') then
          ! Assume B = 10 T, and v_x = 0 initially
          ! x = vx + |vy|* m/ M * q -> offset probably to put on orbit
          x(1) = x(1) + abs(v(2)) * m / (q * 10.0d0)

       end if
    case (3) ! jrm09
       v = v0
       q = (charge * ipart) / n_particles
       ! Assume B = 10 T, and v_x = 0 initially
       ! x = vx + |vy|*m/10q
       !x(1) = x(1) + abs(v(2)) * m / (q * 10.0d0)

       ! Distribute over circle, velocity inwards. Avoid pi/4.
       phi = ((ipart+0.125d0) * 2 * acos(-1.0d0)) / n_particles
       x = norm2(x0) * [cos(phi), sin(phi), 0.0d0]

       ! Add Maxwellian velocity. Random numbers come in pairs of two
       tmp_vec(1:2) = rng%two_normals()
       tmp_vec(3:4) = rng%two_normals()
       v = v0 + tmp_vec(1:3) * maxwellian_velocity
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
       v = v0
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
end module mod_usr
