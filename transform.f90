program transform

    use mpi

    implicit none

    integer :: simpson_n(3), i, k, num_channel_funcs, num_basis_funcs, nc1
    real :: mixed_lims(6), single_lims_a(4), single_lims_b(4), integral, mass_a, mass_b, mass_c, &
        mass_total, mu_a, mu_b, m_a, m_b, values(20), R_a_bound, R_b_bound, mixed_int_value
    logical :: save_inputs, mixed_int

    ! MPI variables
    integer :: comm, rank, comm_size, ierr, group, sub_group, sub_ranks(1), sub_comm
    double precision :: time(4), t_diff

    ! Initialise MPI
    comm = MPI_COMM_WORLD
    call MPI_Init(ierr)

    ! Start the timer
    call MPI_Barrier(comm, ierr)
    time(1) = MPI_Wtime()

    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, comm_size, ierr)

    ! Create new group with only process 0, and create corresponding new comm
    sub_ranks(1) = 0
    call MPI_Comm_group(comm, group, ierr);
    call MPI_Group_incl(group, 1, sub_ranks, sub_group, ierr)
    call MPI_Comm_create(MPI_COMM_WORLD, sub_group, sub_comm, ierr)

    save_inputs = .true.

    ! Number of Simpson cells to split Jacobi ingegral into
    ! (R_b_n, gamma_ab_n) at the R_a boundary or (R_a_n, gamma_ab_n) at the R_b boundary
    simpson_n = (/ 250, 250, 250 /)

    ! Boundary values for R_a and R_b
    R_a_bound = 1.
    R_b_bound = 2.

    ! Limits of integration for mixed Jacobi coordinates
    ! (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
    ! At the R_a boundary, R_a = R_a_bound and remaining four limits used
    ! At the R_b boundary, R_b = R_b_bound and remaining four limits used
    ! Note: 4.*atan(1.) = pi
    mixed_lims = (/ 0., 50., 0., 50., 0., 4.*atan(1.) /)

    ! Calculate mass relations (three masses defined within)
    call calc_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

    ! Limits of integration for single Jacobi coordinates
    ! (r_a_min, r_a_max, gamma_a_min, gamma_a_max)
    single_lims_a = (/ 0., 0., 0., 4.*atan(1.) /)
    single_lims_a(2) = mu_a * ((R_a_bound / mass_c) + (mixed_lims(4) / m_b))
    ! (r_b_min, r_b_max, gamma_b_min, gamma_b_max)
    single_lims_b = (/ 0., 0., 0., -4.*atan(1.) /)
    single_lims_b(2) = mu_b * ((R_b_bound / mass_c) + (mixed_lims(2) / m_a))

    if (rank == 0 .and. save_inputs) then

            if (mixed_int) then
                mixed_int_value = 1.
            else
                mixed_int_value = 0.
            end if

            values = (/ mass_a, mass_b, mass_c, mixed_lims(1), mixed_lims(2), mixed_lims(3), &
                mixed_lims(4), mixed_lims(5), mixed_lims(6), single_lims_a(1), single_lims_a(2), &
                single_lims_a(3), single_lims_a(4), single_lims_b(1), single_lims_b(2), &
                single_lims_b(3), single_lims_b(4), R_a_bound, R_b_bound, mixed_int_value /)
            open (unit=41, file='outputs/values', form='unformatted')
            write(41) values
            close (unit=41)
    end if

    ! Check the timer before integration
    call MPI_Barrier(comm, ierr)
    time(2) = MPI_Wtime()

    num_channel_funcs = 5
    num_basis_funcs = 5
    nc1 = 3

    if (rank == 0) then
        print *, "Transforming in mixed basis..."
    end if
    mixed_int = .true.
    call transform_all_amps(num_channel_funcs, num_basis_funcs, nc1, simpson_n, mixed_lims, &
        single_lims_a, single_lims_b, R_a_bound, R_b_bound, mu_a, mu_b, m_a, m_b, mass_c, &
        mixed_int, comm)

    ! Check the timer after transformation
    call MPI_Barrier(comm, ierr)
    time(3) = MPI_Wtime()
    if (rank == 0) then
            t_diff = time(3) - time(2)
            print *, "Integration time: ", t_diff
    end if

    if (rank == 0) then
        print *, "Transforming in single basis..."
    end if
    mixed_int = .false.
    call transform_all_amps(num_channel_funcs, num_basis_funcs, nc1, simpson_n, mixed_lims, &
    single_lims_a, single_lims_b, R_a_bound, R_b_bound, mu_a, mu_b, m_a, m_b, mass_c, &
    mixed_int, comm)

    ! Check the timer after transformation, at end of program
    call MPI_Barrier(comm, ierr)
    time(4) = MPI_Wtime()
    if (rank == 0) then
            t_diff = time(4) - time(3)
            print *, "Integration time: ", t_diff

            t_diff = time(4) - time(1)
            print *, "Program time: ", t_diff
    end if

    call MPI_Finalize(ierr)

contains

    subroutine calc_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

            implicit none

            real, intent(inout) :: mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b

            ! Atom masses
            mass_a = 1.
            mass_b = 1.
            mass_c = 1.

            ! Sum of masses
            mass_total = 0.
            mass_total = mass_a + mass_b + mass_c

            ! Internal reduced masses
            m_a = 0.
            m_b = 0.
            m_a = mass_b * mass_c / (mass_b + mass_c)
            m_b = mass_a * mass_c / (mass_a + mass_c)

            ! Reduced channel masses
            mu_a = 0.
            mu_b = 0.
            mu_a = mass_a * mass_b * mass_c / (mass_total * m_a)
            mu_b = mass_a * mass_b * mass_c / (mass_total * m_b)

    end subroutine calc_masses

    ! Channel functions (phi or theta) in single Jacobi coordinates
    ! If config_a, x = R_a, r = r_a and gamma = gamma_a
    ! Else x = R_b, r = r_b and gamma = gamma_b
    function get_single_phi(i, x, r, gamma, config_a, conj) result(func)

        implicit none

        integer, intent(in) :: i
        real, intent(in) :: x, r, gamma
        logical, intent(in) :: config_a, conj

        real :: func

        if (config_a) then
            ! phi defined for 1 to nc1 (=3)
            func = x**2. * r**2. * cos(gamma)
            func = func * exp(-1. * x * r)

        else
            ! theta defined for 1 to nc2 (=nc-nc1=2)
            func = x**2. * r**2. * cos(gamma)
            func = func * exp(-1. * x * r)
        end if

    end function get_single_phi


    ! Channel functions (phi or theta) in mixed Jacobi coordinates
    ! If config_a, x = R_a, y = R_b and z = gamma_ab
    ! Else x = R_b, y = R_a and z = gamma_ab
    function get_mixed_phi(i, x, y, z, config_a, conj) result(func)

        implicit none

        integer, intent(in) :: i
        real, intent(in) :: x, y, z
        logical, intent(in) :: config_a, conj

        real :: func

        if (config_a) then
            ! phi defined for 1 to nc1 (=3)
            func = x**2. * y**2. * cos(z)
            func = func * exp(-1. * x * y)
        else
            ! theta defined for 1 to nc2 (=nc-nc1=2)
            func = x**2. * y**2. * cos(z)
            func = func * exp(-1. * x * y)
        end if

    end function get_mixed_phi


    ! Basis functions in single Jacobi coordinates
    ! If config_a, x = R_a, r = r_a and gamma = gamma_a
    ! Else x = R_b, r = r_b and gamma = gamma_b
    function get_single_psi(k, x, r, gamma, config_a, conj) result(func)

        implicit none

        integer, intent(in) :: k
        real, intent(in) :: x, r, gamma
        logical, intent(in) :: config_a, conj

        real :: func

        if (config_a) then
            ! psi defined for 1 to nc1 (=3)
            func = x**2. * r**2. * cos(gamma)
            func = func * exp(-1. * x * r)
        else
            ! psi defined for 1 to nc2 (=nc-nc1=2)
            func = x**2. * r**2. * cos(gamma)
            func = func * exp(-1. * x * r)
        end if

    end function get_single_psi


    ! Basis functions in mixed Jacobi coordinates
    ! If config_a, x = R_a, y = R_b and z = gamma_ab
    ! Else x = R_b, y = R_a and z = gamma_ab
    function get_mixed_psi(k, x, y, z, config_a, conj) result(func)

        implicit none

        integer, intent(in) :: k
        real, intent(in) :: x, y, z
        logical, intent(in) :: config_a, conj

        real :: func

        if (config_a) then
            ! psi defined for 1 to nc1 (=3)
            func = x**2. * y**2. * cos(z)
            func = func * exp(-1. * x * y)
        else
            ! psi defined for 1 to nc2 (=nc-nc1=2)
            func = x**2. * y**2. * cos(z)
            func = func * exp(-1. * x * y)
        end if

    end function get_mixed_psi


    ! Calculates total integrand
    ! For config_a x = R_a, y = R_b and z = gamma_ab
    ! For config_a x = R_b, y = R_a and z = gamma_ab
    ! x is not integrated over, but the integrand may be a function of x
    function integrand_func(x, y, z, mu_a, mass_c, m_b, trans_coords, i, j, config_a, mixed_int, &
        channel_func_1, single_func_1, channel_func_2, single_func_2) result(total_integrand)

            implicit none

            integer, intent(in) :: i, j
            real, intent(in) :: x, y, z, mu_a, mass_c, m_b, trans_coords(2)
            logical, intent(in) :: config_a, mixed_int, channel_func_1, single_func_1, &
                channel_func_2, single_func_2

            real :: R_1, R_2, gamma_ab, small_r_1, gamma_1, total_integrand, integrand_volume, integrand_1, integrand_2
            logical :: conj

            ! Volume component to be integrated
            integrand_volume = y**2. * sin(z)

            if (mixed_int) then
                ! x, y, z coords are mixed Jacobi
                ! Transformed coords are single Jacobi
                R_1 = x
                R_2 = y
                gamma_ab = z
                small_r_1 = trans_coords(1)
                gamma_1 = trans_coords(2)
            else
                ! x, y, z coords are single Jacobi
                ! Transformed coords are mixed Jacobi
                R_1 = x
                small_r_1 = y
                gamma_1 = z
                R_2 = trans_coords(1)
                gamma_ab = trans_coords(2)
            end if

            ! Function to be integrated:

            ! To do: check if i, j are correct in all cases or need to add nc1 etc.

            ! Complex conjugate - not currently used
            conj = .true.
            if (channel_func_1) then
                if (single_func_1) then
                    integrand_1 = get_single_phi(i, R_1, small_r_1, gamma_1, config_a, conj)
                else
                    integrand_1 = get_mixed_phi(i, R_1, R_2, gamma_ab, config_a, conj)
                end if
            else
                if (single_func_1) then
                    integrand_1 = get_single_psi(i, R_1, small_r_1, gamma_1, config_a, conj)
                else
                    integrand_1 = get_mixed_psi(i, R_1, R_2, gamma_ab, config_a, conj)
                end if
            end if

            conj = .false.
            if (channel_func_2) then
                if (single_func_2) then
                    integrand_2 = get_single_phi(j, R_1, small_r_1, gamma_1, config_a, conj)
                else
                    integrand_2 = get_mixed_phi(j, R_1, R_2, gamma_ab, config_a, conj)
                end if
            else
                if (single_func_2) then
                    integrand_2 = get_single_psi(j, R_1, small_r_1, gamma_1, config_a, conj)
                else
                    integrand_2 = get_mixed_psi(j, R_1, R_2, gamma_ab, config_a, conj)
                end if
            end if

            total_integrand = integrand_volume * integrand_1 * integrand_2

    end function integrand_func


    ! Use Simpson's rule to integrate f(x, y, z) over two variables: y and z
    ! Values of x, y and z defined by mixed_int and config_a flags
    ! config_a and mixed_int: x = R_a, y = R_b and z = gamma_ab
    ! not config_a and mixed_int: x = R_b, y = R_a and z = gamma_ab
    ! config_a and not mixed_int: x = R_a, y = r_a, z = gamma_a
    ! not config_a and not mixed_int: x = R_b, y = r_b, z = gamma_b
    function integrate_double_Simpson(n, lims, x, mu_1, mu_2, m_1, m_2, mass_c, &
        trans_coords, config_a, mixed_int, channel_func_1, single_func_1, channel_func_2, &
        single_func_2) result(total_integral)

            implicit none

            integer, intent(in) :: n(2)
            real, intent(in) :: lims(4), x, mu_1, mu_2, m_1, m_2, mass_c, &
                trans_coords(2, n(1)+1, n(2)+1)
            logical, intent(in) :: config_a, channel_func_1, single_func_1, channel_func_2, &
                single_func_2, mixed_int

            real :: width(2), y, z, temp_integral, z_integral, total_integral
            integer :: i, j
            logical :: verbose

            verbose = .true.

            total_integral = 0.
            temp_integral = 0.

            width(1) = abs(lims(2) - lims(1)) / real(n(1))

            ! Integrate each variable in turn, covering full limits
            ! Variables labelled x, y and z for simplicity
            do i = 0, n(1)

                    ! Set value for y at fixed x
                    y = lims(1) + real(i) * width(1)

                    ! Total z integral given x and y
                    z_integral = 0.

                    do j = 0, n(2)

                            ! Set value for z
                            width(2) = abs(lims(4) - lims(3)) / real(n(2))
                            z = lims(3) + real(j) * width(2)

                            ! Evaluate integral at this point
                            temp_integral = integrand_func(x, y, z, mu_a, mass_c, m_b, &
                                trans_coords(:, i+1, j+1), i+1, j+1, config_a, mixed_int, &
                                channel_func_1, single_func_1, channel_func_2, single_func_2)

                            if (j == 0 .or. j == n(2)) then
                                    z_integral = z_integral + temp_integral
                            else if (mod(j, 2) == 0) then
                                    z_integral = z_integral + 2. * temp_integral
                            else
                                    z_integral = z_integral + 4. * temp_integral
                            end if

                    end do

                    ! Total gamma_ab intergral for given R_a and R_b
                    z_integral = width(2) * z_integral / 3.

                    ! Use Simpon's rule to add contributions for this subinterval
                    temp_integral = z_integral

                    if (i == 0 .or. i == n(1)) then
                        total_integral = total_integral + temp_integral
                    else if (mod(i, 2) == 0) then
                        total_integral = total_integral + 2. * temp_integral
                    else
                        total_integral = total_integral + 4. * temp_integral
                    end if

            end do

            ! Total integral
            total_integral = width(1) * total_integral / 3.

            ! Scale for correct volume element in mixed basis
            if (mixed_int) then
                total_integral = ( (mu_1 * mu_2) / (m_1 * m_2) )**(3./2.) * total_integral
            end if

    end function integrate_double_Simpson


    ! Calculate r_a or r_b from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! For r_a, R_1 = R_a, R_2 = R_b, mu_1 = mu_a and m_1 = m_b
    ! For r_b, R_1 = R_b and R_2 = R_a, mu_1 = mu_b and m_1 = m_a
    function calc_r_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c) result(r)

        implicit none

        real, intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c

        real :: r


        r = (R_1 / mass_c)**2.  + (R_2 / m_1)**2. + 2. * &
            (R_1 * R_2 * cos(gamma_ab) / (mass_c * m_1))

        ! Set r to 0 for invalid coordinates, but should not occur
        if (r >= 0) then
            r = mu_1 * sqrt(r)
        else
            print *, "Warning: invalid r"
            r = 0.
        end if

    end function calc_r_single


    ! Calculate R_b or R_a from single Jacobi coordinates (R_a, r_a and gamma_a) or
    ! (R_b, r_b and gamma_b)
    ! For R_b:
        ! R_1 = R_a, small_r_1 = r_a, gamma_1 = gamma_a, mu_1 = mu_a, m_1 = m_b, config_a = .true.
    ! For R_a:
        ! R_1 = R_b, small_r_1 = r_b, gamma_1 = gamma_b, mu_1 = mu_b, m_1 = m_a, config_a = .false.
    function calc_R_mixed(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a) result(R_2)

        implicit none

        real, intent(in) :: R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real :: R_2

        R_2 = (R_1 / mass_c)**2. + (small_r_1 / mu_1)**2.

        if (config_a) then
            R_2 = R_2 + (2. * R_1 * small_r_1 * cos(gamma_1) / (mu_1 * mass_c))
        else
            R_2 = R_2 - (2. * R_1 * small_r_1 * cos(gamma_1) / (mu_1 * mass_c))
        end if

        R_2 = m_1 * sqrt(R_2)

    end function calc_R_mixed


    ! Calculate gamma_a or gamma_b from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! For gamma_a, R_1 = R_a, R_2 = R_b, mu_1 = mu_a and m = m_b
    ! For gamma_b, R_1 = R_b and R_2 = R_a, mu_1 = mu_b and m = m_a
    function calc_gamma_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a) result(gamma_1)

        implicit none

        real, intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real :: r, gamma_1, cos_numerator, err

        err = 0.000001

        ! Angle undefined for R_1 = 0
        if (R_1 == 0.) then
            gamma_1 = 0.
            return
        end if

        r = calc_r_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c)

        ! Angle undefined for r = 0
        if (r == 0.) then
            gamma_1 = 0.
            return
        end if

        ! Multiply by sgn(gamma_ba) for gamma_a, sgn(gamma_ab) fo gamma_b
        if (config_a) then
            cos_numerator = - mu_1 * ( (R_1 / mass_c) + (R_2 / m_1) * cos(gamma_ab) )
        else
            cos_numerator = mu_1 * ( (R_1 / mass_c) + (R_2 / m_1) * cos(gamma_ab) )
        end if

        gamma_1 = cos_numerator / r

        ! Temporary to prevent NaN when taking acos - there is probably a better solution!
        if (gamma_1 > 1. .and. gamma_1 <= (1.+err)) then
            gamma_1 = 1.
        end if
        if (gamma_1 < -1. .and. gamma_1 >= -(1.+err)) then
            gamma_1 = -1.
        end if

        gamma_1 = acos(gamma_1)

    end function calc_gamma_single


    ! Calculate gamma_ab from single Jacobi coordinates (R_a, r_a and gamma_a) or
    ! (R_b, r_b and gamma_b)
    function calc_gamma_mixed(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a) &
        result(gamma_ab)

        implicit none

        real, intent(in) :: R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real :: gamma_ab, err, R_2

        ! Angle undefined for R_1 = 0
        if (R_1 == 0.) then
            gamma_ab = 0.
            return
        end if

        err = 0.001

        R_2 = calc_R_mixed(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a)

        ! Angle undefined for R_2 = 0
        if (R_2 == 0.) then
            gamma_ab = 0.
            return
        end if

        if (config_a) then
            gamma_ab = (-1. * small_r_1 * cos(gamma_1) / mu_1) - (R_1 / mass_c )
        else
            gamma_ab = (small_r_1 * cos(gamma_1) / mu_1) - (R_1 / mass_c )
        end if

        gamma_ab = (m_1 / R_2) * gamma_ab

        ! Temporary to prevent NaN when taking acos - there is probably a better solution!
        if (gamma_ab > 1.) then
            if (gamma_ab <= (1.+err)) then
                gamma_ab = 1.
            else
                print *, "Warning: invalid gamma_ab"
                gamma_ab = 1.
                stop
            end if
        end if
        if (gamma_ab < -1.) then
            if (gamma_ab >= -(1.+err)) then
                gamma_ab = -1.
            else
                print *, "Warning: invalid gamma_ab"
                gamma_ab = -1.
                stop
            end if
        end if

        gamma_ab = acos(gamma_ab)


    end function calc_gamma_mixed


    ! Calculate single Jacobi coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b)
    ! Choice of single coordinates is determined by the config_a flag
    ! from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    function transform_mixed_to_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a) &
        result(single_coords)

        implicit none

        real, intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real :: single_coords(2)

        ! r_a or r_b
        single_coords(1) = calc_r_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c)

        ! gamma_a or gamma_b
        single_coords(2) = calc_gamma_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a)

    end function transform_mixed_to_single


    ! Calculate mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! from single Jacobi coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b)
    ! Choice of single coordinates is determined by the config_a flag
    function transform_single_to_mixed(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a) &
        result(mixed_coords)

        implicit none

        real, intent(in) :: R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real :: mixed_coords(2)

        ! R_b or R_a
        mixed_coords(1) = calc_R_mixed(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a)

        ! gamma_b or gamma_a
        mixed_coords(2) = calc_gamma_mixed(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a)

    end function transform_single_to_mixed


    ! Transform a grid of mixed Jacobi coordinates (R_a, R_b and gamma_ab) to single Jacobi
    ! coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b) or vice versa
    ! Choice of single coordinates is determined by the config_a flag
    ! Direction of transformation is determined by the mixed_int flag
    ! (mixed_int == .true. means input limits and x are written in mixed coordinates)
    ! Grid defined by limits [R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max]
    ! and number of points in each coordinate [n_R_a, n_R_b, n_gamma_ab]
    ! or equivalent for each set of single Jacobi coordinates
    function transform_grid(x, n, lims, mu_1, m_1, mass_c, config_a, mixed_int) &
        result(trans_coords)

        implicit none

        integer, intent(in) :: n(2)
        real, intent(in) :: x, lims(4), mu_1, m_1, mass_c
        logical, intent(in) :: config_a, mixed_int

        integer :: i, j, k
        real :: y, z, width(2), trans_coords(2, n(1)+1, n(2)+1)

        ! Spacing between points for each coordinate
        width(1) = abs(lims(2) - lims(1)) / real(n(1))
        width(2) = abs(lims(4) - lims(3)) / real(n(2))

        ! Loop over R_b (config_a) or R_a (not config_a) values
        do i = 0, n(1)
            y = lims(1) + real(i) * width(1)
            do j = 0, n(2)
                z = lims(3) + real(j) * width(2)
                if (mixed_int) then
                    trans_coords(:, i+1, j+1) = transform_mixed_to_single(x, y, z, mu_1, &
                        m_1, mass_c, config_a)
                else
                    trans_coords(:, i+1, j+1) = transform_single_to_mixed(x, y, z, mu_1, &
                        m_1, mass_c, config_a)
                end if
            end do
        end do

    end function transform_grid


    ! Transform surface amplitudes w_ik from mixed to single Jacobi coordinates
    function transform_amps(old_amps, num_channel_funcs, num_basis_funcs, simpson_n, lims, mu_1, &
        mu_2, m_1, m_2, mass_c, trans_coords, boundary_val, config_a, mixed_int, sub_comm) &
        result(amps)

        integer, intent(in) :: num_channel_funcs, num_basis_funcs, simpson_n(2)
        real, intent(in) :: old_amps(num_channel_funcs, num_basis_funcs), lims(4), mu_1, mu_2, &
            m_1, m_2, mass_c, trans_coords(2, simpson_n(1)+1, simpson_n(2)+1), boundary_val
        logical, intent(in) :: config_a, mixed_int

        ! MPI variables
        integer, intent(in) :: sub_comm

        integer :: i, k, n, m, imin, imax, progress, num_values
        real :: integral_1(num_channel_funcs), integral_2(num_basis_funcs), &
            amps(num_channel_funcs, num_basis_funcs)
        logical :: verbose, channel_func_1, single_func_1, channel_func_2, single_func_2

        ! MPI variables
        integer :: sub_comm_size, ierr, request, recv_status(MPI_STATUS_SIZE)

        call MPI_Comm_size(sub_comm, sub_comm_size, ierr)

        verbose = .true.

        ! Initialise amps array
        do i = 1, num_channel_funcs
            do k = 1, num_basis_funcs
                amps(i, k) = 0.
            end do
        end do

        ! Divide calculation of amps (i runs between 1 and num_channel_funcs overall)
        imin = 1 + rank * num_channel_funcs / sub_comm_size
        if (rank == sub_comm_size - 1) then
                imax = num_channel_funcs
        else
                imax = (rank + 1) * num_channel_funcs / sub_comm_size
        end if

        ! print *, "On rank ", rank, "Loop min = ", imin, "Loop max = ", imax
        do i = imin, imax

            ! print progress every ~20% for each process
            if ((imax - imin >= 4) .and. verbose) then
                progress = 100 * (i - imin) / (imax - imin)
                if (mod((i - imin), (imax - imin + 1) / 5) == 0 &
                .or. i == imax) then
                        print *, "On rank", rank, "...Progress... ", progress, "%"
                end if
            end if

            do k = 1, num_basis_funcs
                ! Loop over channel functions (phi or theta)
                ! Take inner product of channel function i in single basis
                ! with channel function n in mixed basis
                channel_func_1 = .true.
                single_func_1 = .true.
                channel_func_2 = .true.
                single_func_1 = .false.
                do n = 1, num_channel_funcs
                    integral_1(n) = integrate_double_Simpson(simpson_n, lims, boundary_val, &
                        mu_1, mu_2, m_1, m_2, mass_c, trans_coords, config_a, mixed_int, &
                        channel_func_1, single_func_1, channel_func_2, single_func_2)
                end do
                ! Loop over basis functions (psi)
                ! Take inner product of basis function m in mixed basis
                ! with basis function k in single basis
                channel_func_1 = .false.
                single_func_1 = .false.
                channel_func_2 = .false.
                single_func_1 = .true.
                do m = 1, num_basis_funcs
                    integral_2(m) = integrate_double_Simpson(simpson_n, lims, boundary_val, &
                        mu_1, mu_2, m_1, m_2, mass_c, trans_coords, config_a, mixed_int, &
                        channel_func_1, single_func_1, channel_func_2, single_func_2)
            end do

                do n = 1, num_channel_funcs
                    do m = 1, num_basis_funcs
                        amps(i, k) = amps(i, k) + integral_1(n) * integral_2(m) * old_amps(n, m)
                    end do
                end do
            end do
        end do

        ! Send all amps to rank 0
        if (rank == 0) then
                ! Receive from all other ranks
                do i = 1, sub_comm_size - 1
                        imin = 1 + i * num_channel_funcs / sub_comm_size
                        if (i == sub_comm_size - 1) then
                                imax = num_channel_funcs
                        else
                                imax = (i + 1) * num_channel_funcs / sub_comm_size
                        end if

                        num_values = (imax - imin + 1) * num_basis_funcs
                        call MPI_Recv(amps(imin:imax, :), num_values, &
                                MPI_REAL, i, 0, sub_comm, recv_status, ierr)
                end do
        else
                ! All non-zero ranks send to rank 0
                num_values = (imax - imin + 1) * num_basis_funcs
                call MPI_Ssend(amps(imin:imax, :), num_values, MPI_REAL, 0, 0, sub_comm, &
                    request, ierr)
        end if

        ! All sends/receives complete
        call MPI_Barrier(sub_comm, ierr)

    end function transform_amps


    ! Calculate surface amplitudes at boundary
    function calc_old_amps(num_channel_funcs, num_basis_funcs, simpson_n, lims, trans_coords, &
        boundary_val, mu_1, mu_2, m_1, m_2, mass_c, config_a, mixed_int) result(amps)

        implicit none

        integer, intent(in) :: num_channel_funcs, num_basis_funcs, simpson_n(2)
        real, intent(in) :: lims(4), trans_coords(2, simpson_n(1)+1, simpson_n(2)+1), &
            boundary_val, mu_1, mu_2, m_1, m_2, mass_c
        logical, intent(in) :: config_a, mixed_int

        real :: amps(num_channel_funcs, num_basis_funcs)
        integer :: i, k
        logical :: channel_func_1, single_func_1, channel_func_2, single_func_2

        ! Take inner product of channel function i in the mixed basis
        ! with basis function k in mixed basis
        channel_func_1 = .true.
        single_func_1 = .false.
        channel_func_2 = .false.
        single_func_1 = .false.

        do i = 1, num_channel_funcs
            do k = 1, num_basis_funcs
                ! amps(i, k) = 1.
                amps(i, k) = integrate_double_Simpson(simpson_n, lims, boundary_val, mu_1, mu_2, &
                    m_1, m_2, mass_c, trans_coords, config_a, mixed_int, channel_func_1, &
                    single_func_1, channel_func_2, single_func_2)
                amps(i, k) = amps(i, k) / boundary_val
            end do
        end do

    end function calc_old_amps


    ! Transform surface amplitudes w_ik from mixed to single Jacobi coordinates
    subroutine transform_all_amps(num_channel_funcs, num_basis_funcs, nc1, simpson_n, mixed_lims, &
        single_lims_a, single_lims_b, R_a_bound, R_b_bound, mu_a, mu_b, m_a, m_b, mass_c, &
        mixed_int, sub_comm)

        integer, intent(in) :: num_channel_funcs, num_basis_funcs, simpson_n(3)
        real, intent(in) :: mixed_lims(6), single_lims_a(4), single_lims_b(4), R_a_bound, &
            R_b_bound, mu_a, mu_b, m_a, m_b, mass_c
        logical :: mixed_int

        ! MPI variables
        integer, intent(in) :: sub_comm

        integer :: i, j, nc1, nc2, simpson_n_1(2)
        real :: old_amps(num_channel_funcs, num_basis_funcs), &
            new_amps(num_channel_funcs, num_basis_funcs), &
            trans_coords_a(2, simpson_n(2)+1, simpson_n(3)+1), &
            trans_coords_b(2, simpson_n(1)+1, simpson_n(3)+1), lims_1(4), boundary_val
        logical :: config_a

        ! Initialise old amps
        do i = 1, num_channel_funcs
            do j = 1, num_basis_funcs
                old_amps(i, j) = 0.
            end do
        end do

        ! Initialise transformed amps
        do i = 1, num_channel_funcs
            do j = 1, num_basis_funcs
                new_amps(i, j) = 0.
            end do
        end do

        nc2 = num_channel_funcs - nc1

        ! At the R_a boundary
        config_a = .true.
        simpson_n_1(:) = simpson_n(2:3)
        boundary_val = R_a_bound
        if (mixed_int) then
            lims_1(:) = mixed_lims(3:6)
        else
            lims_1(:) = single_lims_a
        end if

        ! If mixed_int, calculate single Jacobi coordinates (R_a, r_a, gamma_a)
        ! for all mixed coordinates, else calculate mixed Jacobi coordinates (R_a, R_b, gamma_ab)
        ! for all single coordinates
        trans_coords_a = transform_grid(boundary_val, simpson_n_1, lims_1, mu_a, m_b, mass_c, &
            config_a, mixed_int)

        ! Calculate amplitudes in mixed Jacobi coordinates at R_a boundary
        old_amps(:nc1, :) = calc_old_amps(nc1, num_basis_funcs, simpson_n_1, lims_1, &
            trans_coords_a, boundary_val, mu_a, mu_b, m_b, m_a, mass_c, config_a, mixed_int)

        ! Transform surface amplitudes at R_a boundary
        new_amps(:nc1, :) = transform_amps(old_amps(:nc1, :), nc1, num_basis_funcs, simpson_n_1, &
            lims_1, mu_a, mu_b, m_b, m_a, mass_c, trans_coords_a, boundary_val, config_a, &
            mixed_int, sub_comm)

        ! At the R_b boundary
        config_a = .false.
        simpson_n_1(1) = simpson_n(1)
        simpson_n_1(2) = simpson_n(3)
        boundary_val = R_b_bound
        if (mixed_int) then
            lims_1(1:2) = mixed_lims(1:2)
            lims_1(3:4) = mixed_lims(5:6)
        else
            lims_1(:) = single_lims_b
        end if

        ! If mixed_int, calculate single Jacobi coordinates (R_b, r_b, gamma_b)
        ! for all mixed coordinates, else calculate mixed Jacobi coordinates (R_a, R_b, gamma_ab)
        ! for all single coordinates
        trans_coords_b = transform_grid(boundary_val, simpson_n_1, lims_1, mu_b, m_a, mass_c, &
            config_a, mixed_int)

        old_amps(nc1+1:, :) = calc_old_amps(nc2, num_basis_funcs, simpson_n_1, lims_1, &
            trans_coords_b, boundary_val, mu_b, mu_a, m_a, m_b, mass_c, config_a, mixed_int)

        ! Transform surface amplitudes at R_b boundary
        new_amps(nc1+1:, :) = transform_amps(old_amps(nc1+1:, :), nc2, num_basis_funcs, &
            simpson_n_1, lims_1, mu_b, mu_a, m_a, m_b, mass_c, trans_coords_b, boundary_val, &
            config_a, mixed_int, sub_comm)

        ! Print transformed amplitudes - save in future
        if (rank == 0) then

            print *, "___________________________"
            print *, "Original amplitudes"
            print *, old_amps
            print *, "___________________________"

            print *, "___________________________"
            print *, "Transformed amplitudes"
            print *, new_amps
            print *, "___________________________"
        end if

    end subroutine transform_all_amps

end program transform