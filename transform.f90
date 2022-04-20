module precisn
    implicit none
    private
    integer, parameter::   wp = selected_real_kind(12)
    public wp
end module precisn

program transform

    use mpi
    use precisn, only: wp

    implicit none

    call main()

contains

    subroutine main()

        implicit none

        integer :: simpson_n(3), i, k, n_1, nc_1, n_2, nc_2, b, nt
        real(wp) :: mixed_lims(6), integral, mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, &
            m_b, values(10), mixed_int_value
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
        call MPI_Comm_create(comm, sub_group, sub_comm, ierr)

        save_inputs = .true.

        ! Number of Simpson cells to split Jacobi ingegral into
        ! (R_b_n, gamma_ab_n) at the R_a boundary or (R_a_n, gamma_ab_n) at the R_b boundary
        simpson_n = (/ 100, 100, 100 /)

        ! Limits of integration for mixed Jacobi coordinates
        ! Also defines R_a and R_b boundaries used with single Jacobi coordinates
        ! (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
        ! At the R_a boundary, R_a = R_a_max and remaining four limits used
        ! At the R_b boundary, R_b = R_b_max and remaining four limits used
        ! Note: 4._wp*atan(1._wp) = pi
        mixed_lims = (/ 0._wp, 3._wp, 0._wp, 5._wp, 0._wp, 4._wp*atan(1._wp) /)

        ! Calculate mass relations (three masses defined within)
        call calc_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

        if (rank == 0 .and. save_inputs) then

            if (mixed_int) then
                mixed_int_value = 1._wp
            else
                mixed_int_value = 0._wp
            end if

            values = (/ mass_a, mass_b, mass_c, mixed_lims(1), mixed_lims(2), mixed_lims(3), &
                mixed_lims(4), mixed_lims(5), mixed_lims(6), mixed_int_value /)
            open (unit=41, file='outputs/values', form='unformatted')
            write(41) values
            close (unit=41)
        end if

        n_1 = 3 ! Number of channel functions for configuration A
        nc_1 = 2 ! Number of radial continuum basis orbitals for configuration A
        n_2 = 2 ! Number of channel functions for configuration B
        nc_2 = 1 ! Number of radial continuum basis orbitals for configuration B
        b = 0 ! Number of quadratically integrably functions - not used
        nt = n_1 * nc_1 + n_2 * nc_2 + b ! Number of linearly indep basis functions

        if (rank == 0) then
            print *, "Transforming in mixed basis..."
        end if

        ! Check the timer before first transformation
        call MPI_Barrier(comm, ierr)
        time(2) = MPI_Wtime()

        mixed_int = .true.
        call transform_all_amps(n_1, nc_1, n_2, nc_2, nt, simpson_n, mixed_lims, mu_a, mu_b, m_a, &
            m_b, mass_c, mixed_int, comm)

        ! Check the timer after first transformation, before second transformation
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
        call transform_all_amps(n_1, nc_1, n_2, nc_2, nt, simpson_n, mixed_lims, mu_a, mu_b, m_a, &
            m_b, mass_c, mixed_int, comm)

        ! Check the timer after second transformation, at end of program
        call MPI_Barrier(comm, ierr)
        time(4) = MPI_Wtime()
        if (rank == 0) then
            t_diff = time(4) - time(3)
            print *, "Integration time: ", t_diff

            t_diff = time(4) - time(1)
            print *, "Program time: ", t_diff
        end if

        call MPI_Finalize(ierr)

    end subroutine main


    subroutine calc_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

        implicit none

        real(wp), intent(inout) :: mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b

        ! Atom masses
        mass_a = 1._wp
        mass_b = 1._wp
        mass_c = 1._wp

        ! Sum of masses
        mass_total = 0._wp
        mass_total = mass_a + mass_b + mass_c

        ! Internal reduced masses
        m_a = 0._wp
        m_b = 0._wp
        m_a = mass_b * mass_c / (mass_b + mass_c)
        m_b = mass_a * mass_c / (mass_a + mass_c)

        ! Reduced channel masses
        mu_a = 0._wp
        mu_b = 0._wp
        mu_a = mass_a * mass_b * mass_c / (mass_total * m_a)
        mu_b = mass_a * mass_b * mass_c / (mass_total * m_b)

    end subroutine calc_masses

    ! Channel functions (phi or theta) in single Jacobi coordinates
    ! If config_a, R_1 = R_a, small_r_1 = r_a and gamma_1 = gamma_a
    ! Else R_1 = R_b, small_r_1 = r_b and gamma_1 = gamma_b
    function get_single_phi(i, boundary_val_2, R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, &
        config_a, conj) result(func)

        implicit none

        integer, intent(in) :: i
        real(wp), intent(in) :: boundary_val_2, R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a, conj

        real(wp) :: R_2, func

        ! phi defined for 1 to n_1, theta defined for 1 to n_2 (passed as n_1)
        R_2 = calc_R_from_single(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a)

        func = real(i, kind=wp) * cos(gamma_1) * (R_2 - boundary_val_2)**2._wp
        func = func * exp(-0.5_wp * R_1 * small_r_1)

    end function get_single_phi


    ! Channel functions (phi or theta) in mixed Jacobi coordinates
    ! If config_a, x = R_a, y = R_b and z = gamma_ab
    ! Else x = R_b, y = R_a and z = gamma_ab
    function get_mixed_phi(i, boundary_val_2, R_1, R_2, gamma_ab, config_a, conj) result(func)

        implicit none

        integer, intent(in) :: i
        real(wp), intent(in) :: boundary_val_2, R_1, R_2, gamma_ab
        logical, intent(in) :: config_a, conj

        real(wp) :: func

        ! phi defined for 1 to n_1, theta defined for 1 to n_2 (passed as n_1)
        func = real(i, kind=wp) * cos(gamma_ab) * (R_2 - boundary_val_2)**2._wp
        func = func * exp(-0.5_wp * R_1 * R_2)

    end function get_mixed_phi


    ! Basis functions in single Jacobi coordinates
    ! If config_a, R_1 = R_a, small_r_1 = r_a and gamma_1 = gamma_a
    ! Else R_1 = R_b, small_r_1 = r_b and gamma_1 = gamma_b
    ! n_1 is the number of channel functions, nc_1 is number of radial continuum basis orbitals
    function get_single_psi(k, n_1, nc_1, boundary_val_2, R_1, small_r_1, gamma_1, mu_1, m_1, &
        mass_c, config_a, conj) result(func)

        implicit none

        integer, intent(in) :: k, n_1, nc_1
        real(wp), intent(in) :: boundary_val_2, R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a, conj

        real(wp) :: func
        integer :: i, j

        ! psi defined for 1 to n_1*nc_1 or 1 to n_2*nc_2 (passed in as n_1 and nc_2 in either case)
        func = 0._wp
        do i = 1, n_1
            do j = 1, nc_1
                func = func + (get_single_phi(i, boundary_val_2, R_1, small_r_1, gamma_1, &
                    mu_1, m_1, mass_c, config_a, conj) * (1._wp / R_1) * &
                    get_single_radial_func(R_1, i, j, config_a) * &
                    get_single_coeff(i, j, k, config_a))
            end do
        end do

    end function get_single_psi


    ! Basis functions in mixed Jacobi coordinates
    ! If config_a, R_1 = R_a, R_2 = R_b
    ! Else R_1 = R_b, R_2 = R_a
    ! n_1 is the number of channel functions, nc_1 is number of radial continuum basis orbitals
    function get_mixed_psi(k, n_1, nc_1, boundary_val_2, R_1, R_2, gamma_ab, config_a, conj) &
        result(func)

        implicit none

        integer, intent(in) :: k, n_1, nc_1
        real(wp), intent(in) :: boundary_val_2, R_1, R_2, gamma_ab
        logical, intent(in) :: config_a, conj

        real(wp) :: func
        integer :: i, j

        ! psi defined for 1 to n_1*nc_1 or 1 to n_2*nc_2 (passed in as n_1 and nc_2 in either case)
        func = 0._wp
        do i = 1, n_1
            do j = 1, nc_1
                func = func + (get_mixed_phi(i, boundary_val_2, R_1, R_2, gamma_ab, config_a, &
                    conj) * (1._wp / R_1) * get_mixed_radial_func(R_1, i, j, config_a) * &
                    get_mixed_coeff(i, j, k, config_a))
            end do
        end do

    end function get_mixed_psi


    ! Radial continuum basis orbitals in single Jacobi coordinates
    ! If config_a, R_1 = R_a, else R_1 = R_b
    function get_single_radial_func(R_1, i, j, config_a) result(radial_func)

        implicit none

        integer, intent(in) :: i, j
        real(wp), intent(in) :: R_1
        logical, intent(in) :: config_a

        real(wp) :: radial_func

        if (config_a) then
            radial_func = real(i, kind=wp) * real(j, kind=wp) * exp(-1._wp * R_1)
        else
            radial_func = real(i, kind=wp) * real(j, kind=wp) * exp(-1._wp * R_1)
        end if

    end function get_single_radial_func


    ! Radial continuum basis orbitals in mixed Jacobi coordinates
    ! If config_a, R_1 = R_a, else R_1 = R_b
    function get_mixed_radial_func(R_1, i, j, config_a) result(radial_func)

        implicit none

        integer, intent(in) :: i, j
        real(wp), intent(in) :: R_1
        logical, intent(in) :: config_a

        real(wp) :: radial_func

        if (config_a) then
            radial_func = real(i, kind=wp) * real(j, kind=wp) * exp(-1._wp * R_1)
        else
            radial_func = real(i, kind=wp) * real(j, kind=wp) * exp(-1._wp * R_1)
        end if

    end function get_mixed_radial_func


    ! Coefficients from basis funcion expansion in single Jacobi coordinates
    function get_single_coeff(i, j, k, config_a) result(coeff)

        implicit none

        integer, intent(in) :: i, j, k
        logical, intent(in) :: config_a

        real(wp) :: coeff

        if (config_a) then
            coeff = real(i, kind=wp) + real(j, kind=wp) + real(k, kind=wp)
            coeff = 1._wp / coeff
        else
            coeff = real(i, kind=wp) + real(j, kind=wp) + real(k, kind=wp)
            coeff = 1._wp / coeff
        end if


    end function get_single_coeff


    ! Coefficients from basis funcion expansion in mixed Jacobi coordinates
    function get_mixed_coeff(i, j, k, config_a) result(coeff)

        implicit none

        integer, intent(in) :: i, j, k
        logical, intent(in) :: config_a

        real(wp) :: coeff

        if (config_a) then
            coeff = real(i, kind=wp) * real(j, kind=wp) * real(k, kind=wp)**2._wp
            coeff = 1._wp / coeff
        else
            coeff = (real(i, kind=wp) * real(j, kind=wp) * real(k, kind=wp))**2._wp
            coeff = 1._wp / coeff
        end if


    end function get_mixed_coeff


    ! Calculates total integrand
    ! For config_a x = R_a, y = R_b and z = gamma_ab
    ! For config_a x = R_b, y = R_a and z = gamma_ab
    ! x is not integrated over, but the integrand may be a function of x
    function integrand_func(x, y, z, boundary_val_2, mu_1, m_1, mass_c, trans_coords, idx_1, &
        idx_2, n_1, nc_1, config_a, mixed_int, channel_func_1, single_func_1, channel_func_2, &
        single_func_2) result(total_integrand)

            implicit none

            integer, intent(in) :: idx_1, idx_2, n_1, nc_1
            real(wp), intent(in) :: x, y, z, boundary_val_2, mu_1, m_1, mass_c, trans_coords(2)
            logical, intent(in) :: config_a, mixed_int, channel_func_1, single_func_1, &
                channel_func_2, single_func_2

            real(wp) :: R_1, R_2, gamma_ab, small_r_1, gamma_1, total_integrand, &
                integrand_volume, integrand_1, integrand_2
            logical :: conj

            ! Volume component to be integrated
            integrand_volume = y**2._wp * sin(z)

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
            ! Complex conjugate - not currently used
            conj = .true.
            if (channel_func_1) then
                if (single_func_1) then
                    integrand_1 = get_single_phi(idx_1, boundary_val_2, R_1, small_r_1, gamma_1, &
                        mu_1, m_1, mass_c, config_a, conj)
                else
                    integrand_1 = get_mixed_phi(idx_1, boundary_val_2, R_1, R_2, gamma_ab, &
                        config_a, conj)
                end if
            else
                if (single_func_1) then
                    integrand_1 = get_single_psi(idx_1, n_1, nc_1, boundary_val_2, R_1, &
                        small_r_1, gamma_1, mu_1, m_1, mass_c, config_a, conj)
                else
                    integrand_1 = get_mixed_psi(idx_1, n_1, nc_1, boundary_val_2, R_1, R_2, &
                        gamma_ab, config_a, conj)
                end if
            end if

            conj = .false.
            if (channel_func_2) then
                if (single_func_2) then
                    integrand_2 = get_single_phi(idx_2, boundary_val_2, R_1, small_r_1, gamma_1, &
                        mu_1, m_1, mass_c, config_a, conj)
                else
                    integrand_2 = get_mixed_phi(idx_2, boundary_val_2, R_1, R_2, gamma_ab, &
                        config_a, conj)
                end if
            else
                if (single_func_2) then
                    integrand_2 = get_single_psi(idx_2, n_1, nc_1, boundary_val_2, R_1, &
                        small_r_1, gamma_1, mu_1, m_1, mass_c, config_a, conj)
                else
                    integrand_2 = get_mixed_psi(idx_2, n_1, nc_1, boundary_val_2, R_1, R_2, &
                        gamma_ab, config_a, conj)
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
    function integrate_double_Simpson(simpson_n, x, boundary_val_2, n_1, nc_1, channel_num, &
        basis_num, mu_1, mu_2, m_1, m_2, mass_c, coords, trans_coords, config_a, mixed_int, &
        channel_func_1, single_func_1, channel_func_2, single_func_2) result(total_integral)

            implicit none

            integer, intent(in) :: simpson_n(2), n_1, nc_1, channel_num, basis_num
            real(wp), intent(in) :: x, boundary_val_2, mu_1, mu_2, m_1, m_2, mass_c, coords(2, &
                simpson_n(1)+1, simpson_n(2)+1), trans_coords(2, simpson_n(1)+1, simpson_n(2)+1)
            logical, intent(in) :: config_a, channel_func_1, single_func_1, channel_func_2, &
                single_func_2, mixed_int

            real(wp) :: width(2), y, z, temp_integral, z_integral, total_integral
            integer :: i, j
            logical :: verbose

            ! verbose = .true.
            verbose = .false.

            total_integral = 0._wp
            temp_integral = 0._wp

            ! Get range of y for current x
            width(1) = abs(coords(1, 1, 1) - coords(1, simpson_n(1)+1, 1))
            width(1) = width(1) / real(simpson_n(1), kind=wp)

            ! Integrate each variable in turn, covering full limits
            ! Variables labelled x, y and z for simplicity
            do i = 0, simpson_n(1)

                ! Set value for y at fixed x
                y = coords(1, i+1, 1)

                ! Total z integral given x and y
                z_integral = 0._wp

                ! Get range of z for current x and y
                width(2) = abs(coords(2, i+1, 1) - coords(2, i+1, simpson_n(2)+1))
                width(2) = width(2) / real(simpson_n(2), kind=wp)

                do j = 0, simpson_n(2)

                    if (width(2) == 0._wp) then
                        exit
                    end if

                    ! Set value for z
                    z = coords(2, i+1, j+1)

                    ! Evaluate integral at this point
                    temp_integral = integrand_func(x, y, z, boundary_val_2, mu_1, m_1, mass_c, &
                        trans_coords(:, i+1, j+1), channel_num, basis_num, n_1, nc_1, config_a, &
                        mixed_int, channel_func_1, single_func_1, channel_func_2, single_func_2)

                    if (j == 0 .or. j == simpson_n(2)) then
                        z_integral = z_integral + temp_integral
                    else if (mod(j, 2) == 0) then
                        z_integral = z_integral + 2._wp * temp_integral
                    else
                        z_integral = z_integral + 4._wp * temp_integral
                    end if

                end do

                ! Total z intergral for given x and y
                z_integral = width(2) * z_integral / 3._wp

                ! Use Simpon's rule to add contributions for this subinterval
                temp_integral = z_integral

                if (i == 0 .or. i == simpson_n(1)) then
                    total_integral = total_integral + temp_integral
                else if (mod(i, 2) == 0) then
                    total_integral = total_integral + 2._wp * temp_integral
                else
                    total_integral = total_integral + 4._wp * temp_integral
                end if

            end do
            ! Total integral
            total_integral = width(1) * total_integral / 3._wp

            ! Scale for correct volume element in mixed basis
            if (mixed_int) then
                total_integral = ( (mu_1 * mu_2) / (m_1 * m_2) )**(3._wp/2._wp) * total_integral
            end if

    end function integrate_double_Simpson


    ! Estimate the minimum and maximum values of r_a or r_b, based on known R_a or R_b
    ! or gamma_a or gamma_b, based on known R_a and r_a or R_b and r_b
    function estimate_jacobi_lims(mixed_lims, R_1, small_r_1, mu_1, m_1, mass_c, estimating_r_1, &
        estimating_gamma_1, config_a) result(single_lims)

        implicit none

        real(wp), intent(in) :: R_1, mixed_lims(4), mu_1, m_1, mass_c, small_r_1
        logical, intent(in) :: estimating_r_1, estimating_gamma_1, config_a

        real(wp) :: test_coord(100000), mixed_coords(2), width(2), single_lims(2), test_r_a, &
            temp_gamma
        integer :: i, j, n, count
        logical :: random_range

        ! Use randomly generated coordinates, rather than even distribution across range
        random_range = .true.

        if (random_range) then
            n = size(test_coord)
        else
            n = sqrt(real(size(test_coord), kind=wp))
        end if

        ! Ranges of R_2 and gamma_ab
        width(1) = abs(mixed_lims(2) - mixed_lims(1))
        width(2) = abs(mixed_lims(4) - mixed_lims(3))

        ! Track valid coordinates generated
        count = 0

        if (random_range) then
            do i = 1, n - 4

                ! Generate random R_2 and gamma_ab that change each iteration
                call random_number(mixed_coords)

                mixed_coords(1) = mixed_lims(1) + mixed_coords(1) * width(1)
                mixed_coords(2) = mixed_lims(3) + mixed_coords(2) * width(2)

                if (estimating_r_1) then
                    count = count + 1
                    ! Calculate r_1 from R_1, R_2 and gamma_ab
                    test_coord(count) = calc_small_r_from_mixed(R_1, mixed_coords(1), &
                        mixed_coords(2), mu_1, m_1, mass_c)
                else if (estimating_gamma_1) then

                    ! Attempt to calculate gamma_1 from R_1, r_1 and gamma_ab
                    temp_gamma = calc_gamma_from_integral(mixed_lims, R_1, small_r_1, &
                        mixed_coords(2), mu_1, m_1, mass_c, .true., .true., config_a)

                    ! If R_b is invalid, ignore gamma calculated
                    if (temp_gamma /= -100._wp) then
                        count = count + 1
                        test_coord(count) = temp_gamma
                    end if
                end if
            end do

            do i = n-3, n
                ! Explicitly test ends of limits
                if (i == (n-3) .or. i == (n-2)) then
                    mixed_coords(1) = mixed_lims(1)
                else
                    mixed_coords(1) = mixed_lims(2)
                end if
                if (i == (n-3) .or. i == (n-1)) then
                    mixed_coords(2) = mixed_lims(3)
                else
                    mixed_coords(2) = mixed_lims(4)
                end if

                if (estimating_r_1) then
                    count = count + 1
                    ! Calculate r_1 from R_1, R_2 and gamma_ab
                    test_coord(count) = calc_small_r_from_mixed(R_1, mixed_coords(1), &
                        mixed_coords(2), mu_1, m_1, mass_c)
                else if (estimating_gamma_1) then

                    ! Attempt to calculate gamma from R_1, r_1 and gamma_ab
                    temp_gamma = calc_gamma_from_integral(mixed_lims, R_1, small_r_1, &
                        mixed_coords(2), mu_1, m_1, mass_c, .true., .true., config_a)

                    ! If R_2 is invalid, ignore gamma calculated
                    if (temp_gamma /= -100._wp) then
                        count = count + 1
                        test_coord(count) = temp_gamma
                    end if
                end if
            end do
        else
            do i = 1, n

                ! Evenly sample entire range of each coord
                mixed_coords(1) = mixed_lims(1) + (real(i-1, kind=wp) * &
                    width(1) / real(n, kind=wp))

                do j = 1, n
                    ! Evenly sample entire range of each coord
                    mixed_coords(2) = mixed_lims(3) + (real(j-1, kind=wp) * width(2) / &
                        real(n, kind=wp))

                    if (estimating_r_1) then
                        ! Calculate r_1 from R_1, R_2 and gamma_ab
                        test_coord(count) = calc_small_r_from_mixed(R_1, mixed_coords(1), &
                            mixed_coords(2), mu_1, m_1, mass_c)
                        count = count + 1
                    else if (estimating_gamma_1) then
                        ! Attempt to calculate gamma_1 from R_1, r_1 and gamma_ab
                        temp_gamma = calc_gamma_from_integral(mixed_lims, R_1, small_r_1, &
                            mixed_coords(2), mu_1, m_1, mass_c, .true., .true., config_a)

                        ! If R_2 is invalid, ignore gamma_1 calculated
                        if (temp_gamma /= -100._wp) then
                            count = count + 1
                            test_coord(count) = temp_gamma
                        end if
                    end if
                end do
            end do
        end if

        if (count > 0) then
            single_lims(1) = minval(test_coord(1:count))
            single_lims(2) = maxval(test_coord(1:count))
        else
            single_lims(1) = 0._wp
            single_lims(2) = 0._wp
        end if

    end function estimate_jacobi_lims


    function get_jacobi_lims(mixed_lims, R_1, small_r_1, mu_1, m_1, mass_c, get_r_lims, &
        get_gamma_lims, estimate_lims, config_a) result(lims)

        implicit none

        real(wp), intent(in) :: R_1, small_r_1, mixed_lims(4), mu_1, m_1, mass_c
        logical, intent(in) :: get_r_lims, get_gamma_lims, estimate_lims, config_a

        real(wp) :: lims(2), gamma_tol, temp_gamma

        if (get_r_lims .and. get_gamma_lims) then
            print *, "Only one coordinate flag should be set to true"
            stop
        end if

        gamma_tol = 0.0000001_wp

        if (get_r_lims) then
            if (estimate_lims) then
                lims(1:2) = estimate_jacobi_lims(mixed_lims, R_1, 0._wp, mu_1, m_1, mass_c, &
                    .true., .false., config_a)
            else
                ! If gamma_ab can be greater than pi/2, cos(gamma_ab) < 0 is possible
                if (mixed_lims(4) >= 2._wp*atan(1._wp)) then
                    ! Check negative value large enough to offset, R_2^2, otherwise set R_2 = 0
                    if (abs(2._wp * R_1 * cos(mixed_lims(4)) / mass_c) > &
                    (mixed_lims(2) / m_1)) then
                        lims(1) = calc_small_r_from_mixed(R_1, abs(m_1 * R_1 * &
                            cos(mixed_lims(4)) / mass_c), mixed_lims(4), mu_1, m_1, mass_c)
                    else
                        lims(1) = calc_small_r_from_mixed(R_1, mixed_lims(1), mixed_lims(4), &
                            mu_1, m_1, mass_c)
                    end if
                else
                    lims(1) = calc_small_r_from_mixed(R_1, mixed_lims(1), mixed_lims(4), mu_1, &
                        m_1, mass_c)
                end if

                ! r_1_max depends on current R_1, R_2_max and gamma_ab_min (if less than pi/2)
                lims(2) = calc_small_r_from_mixed(R_1, mixed_lims(2), mixed_lims(3), mu_1, m_1, &
                    mass_c)
            end if

        else if (get_gamma_lims) then
            if (estimate_lims) then
                lims(1:2) = estimate_jacobi_lims(mixed_lims, R_1, small_r_1, mu_1, m_1, mass_c, &
                    .false., .true., config_a)
            else

                ! Note: This option may not work for all R_1 and R_2 limits
                ! gamma_ab must run between 0 and pi
                ! Use with caution

                ! gamma_a is undefined for R_1 = 0 or r_1 = 0
                if (R_1 == 0._wp .or. small_r_1 == 0._wp) then
                    lims(1) = 0._wp
                    lims(2) = 0._wp
                else
                    ! Calculate cos(gamma_1) when R_2 is maximum
                    ! to get the minimum gamma_1
                    lims(1) = mass_c * mu_1 * ((mixed_lims(4) / m_1 )**2._wp - &
                        (R_1 / mass_c)**2._wp - (small_r_1 / mu_1)**2._wp) &
                        / (2._wp * R_1 * small_r_1)

                    ! There may be no valid cos(gamma_a) for maximum R_b
                    ! In this case, gamma_a spans full range
                    if (lims(1) > 1._wp) then
                        lims(1) = 0._wp
                    else if (lims(1) < -1._wp) then
                        ! Allow for rounding errors near pi
                        if (lims(1) > -(1._wp+gamma_tol)) then
                            lims(1) = acos(-1._wp)
                        else
                            lims(1) = 0._wp
                        end if
                    else
                        lims(1) = acos(lims(1))
                    end if

                    ! Get upper limit, starting at lower limit
                    lims(2) = lims(1)

                    ! Test gamma_ab minimum for +- roots
                    temp_gamma = calc_gamma_from_integral(mixed_lims, R_1, small_r_1, &
                        mixed_lims(3), mu_1, m_1, mass_c, .false., .true., config_a)
                    if (temp_gamma /= -100._wp) then
                        lims(2) = temp_gamma
                    end if
                    temp_gamma = calc_gamma_from_integral(mixed_lims, R_1, small_r_1, &
                        mixed_lims(3), mu_1, m_1, mass_c, .false., .false., config_a)
                    if (temp_gamma /= -100._wp .and. temp_gamma > lims(2)) then
                        lims(2) = temp_gamma
                    end if

                    ! Test gamma_ab maximum for +- roots
                    temp_gamma = calc_gamma_from_integral(mixed_lims, R_1, small_r_1, &
                        mixed_lims(4), mu_1, m_1, mass_c, .false., .true., config_a)
                    if (temp_gamma /= -100._wp .and. temp_gamma > lims(2)) then
                            lims(2) = temp_gamma
                    end if
                    temp_gamma = calc_gamma_from_integral(mixed_lims, R_1, small_r_1, &
                            mixed_lims(4), mu_1, m_1, mass_c, .false., .false., config_a)
                    if (temp_gamma /= -100._wp .and. temp_gamma > lims(2)) then
                            lims(2) = temp_gamma
                    end if
                end if
            end if
        end if

    end function get_jacobi_lims


    ! Calculate r_a or r_b from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! For r_a, R_1 = R_a, R_2 = R_b, mu_1 = mu_a and m_1 = m_b
    ! For r_b, R_1 = R_b and R_2 = R_a, mu_1 = mu_b and m_1 = m_a
    function calc_small_r_from_mixed(R_1, R_2, gamma_ab, mu_1, m_1, mass_c) result(r)

        implicit none

        real(wp), intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c

        real(wp) :: r

        r = (R_1 / mass_c)**2._wp  + (R_2 / m_1)**2._wp + 2._wp * &
            (R_1 * R_2 * cos(gamma_ab) / (mass_c * m_1))

        ! Set r to 0 for invalid coordinates, but should not occur
        if (r >= 0) then
            r = mu_1 * sqrt(r)
        else
            print *, "Warning: invalid r: ", r
            r = 0._wp
        end if

    end function calc_small_r_from_mixed


    ! Calculate R_b or R_a from single Jacobi coordinates (R_a, r_a and gamma_a) or
    ! (R_b, r_b and gamma_b)
    ! For R_b:
        ! R_1 = R_a, small_r_1 = r_a, gamma_1 = gamma_a, mu_1 = mu_a, m_1 = m_b, config_a = .true.
    ! For R_a:
        ! R_1 = R_b, small_r_1 = r_b, gamma_1 = gamma_b, mu_1 = mu_b, m_1 = m_a, config_a = .false.
    function calc_R_from_single(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a) result(R_2)

        implicit none

        real(wp), intent(in) :: R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real(wp) :: R_2

        R_2 = (R_1 / mass_c)**2._wp + (small_r_1 / mu_1)**2._wp

        if (config_a) then
            R_2 = R_2 + (2._wp * R_1 * small_r_1 * cos(gamma_1) / (mu_1 * mass_c))
        else
            R_2 = R_2 - (2._wp * R_1 * small_r_1 * cos(gamma_1) / (mu_1 * mass_c))
        end if

        ! Catch numerical errors
        if (R_2 < 0._wp) then
            R_2 = 0._wp
        else
            R_2 = m_1 * sqrt(R_2)
        end if

    end function calc_R_from_single


    ! Calculate gamma_a or gamma_b from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! For gamma_a, R_1 = R_a, R_2 = R_b, mu_1 = mu_a and m = m_b
    ! For gamma_b, R_1 = R_b and R_2 = R_a, mu_1 = mu_b and m = m_a
    function calc_gamma_from_mixed(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a) result(gamma_1)

        implicit none

        real(wp), intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real(wp) :: r, gamma_1, cos_numerator, gamma_tol

        gamma_tol = 0.0000001_wp

        ! Angle undefined for R_1 = 0
        if (R_1 == 0._wp) then
            gamma_1 = 0._wp
            return
        end if

        r = calc_small_r_from_mixed(R_1, R_2, gamma_ab, mu_1, m_1, mass_c)

        ! Angle undefined for r = 0
        if (r == 0._wp) then
            gamma_1 = 0._wp
            return
        end if

        ! Multiply by sgn(gamma_ba) for gamma_a, sgn(gamma_ab) fo gamma_b
        if (config_a) then
            cos_numerator = - mu_1 * ( (R_1 / mass_c) + (R_2 / m_1) * cos(gamma_ab) )
        else
            cos_numerator = mu_1 * ( (R_1 / mass_c) + (R_2 / m_1) * cos(gamma_ab) )
        end if

        gamma_1 = cos_numerator / r

        ! Prevent NaN when taking acos - there is probably a better solution!
        if (gamma_1 > 1._wp) then
            if (gamma_1 <= (1._wp+gamma_tol)) then
                gamma_1 = 1._wp
            else
                print *, "Warning: invalid gamma_1: ", gamma_1
                gamma_1 = 1._wp
                stop
            end if
        end if
        if (gamma_1 < -1._wp) then
            if (gamma_1 >= -(1._wp+gamma_tol)) then
                gamma_1 = -1._wp
            else
                print *, "Warning: invalid gamma_1: ", gamma_1
                gamma_1 = -1._wp
                stop
            end if
        end if

        gamma_1 = acos(gamma_1)

    end function calc_gamma_from_mixed


    ! Calculate gamma_ab from single Jacobi coordinates (R_a, r_a and gamma_a) or
    ! (R_b, r_b and gamma_b)
    function calc_gamma_from_single(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a) &
        result(gamma_ab)

        implicit none

        real(wp), intent(in) :: R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real(wp) :: gamma_ab, gamma_tol, R_2

        ! Angle undefined for R_1 = 0
        if (R_1 == 0._wp) then
            gamma_ab = 0._wp
            return
        end if

        gamma_tol = 0.0000001_wp

        R_2 = calc_R_from_single(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a)

        ! Angle undefined for R_2 = 0
        if (R_2 == 0._wp) then
            gamma_ab = 0._wp
            return
        end if

        if (config_a) then
            gamma_ab = (-1._wp * small_r_1 * cos(gamma_1) / mu_1) - (R_1 / mass_c )
        else
            gamma_ab = (small_r_1 * cos(gamma_1) / mu_1) - (R_1 / mass_c )
        end if

        gamma_ab = (m_1 / R_2) * gamma_ab

        ! Prevent NaN when taking acos - there is probably a better solution!
        if (gamma_ab > 1._wp) then
            if (gamma_ab <= (1._wp+gamma_tol)) then
                gamma_ab = 1._wp
            else
                print *, "Warning: invalid gamma_ab: ", gamma_ab
                gamma_ab = 1._wp
                stop
            end if
        end if
        if (gamma_ab < -1._wp) then
            if (gamma_ab >= -(1._wp+gamma_tol)) then
                gamma_ab = -1._wp
            else
                print *, "Warning: invalid gamma_ab: ", gamma_ab
                gamma_ab = -1._wp
                stop
            end if
        end if

        gamma_ab = acos(gamma_ab)


    end function calc_gamma_from_single


    ! Calculate gamma_a from R_1 = R_a, r_1 = r_a and gamma_ab
    ! or gamma_b from R_1 = R_b, r_1 = r_b and gamma_ab
    ! R_2 is calculated first, and -100._wp is returned if it is invalid
    function calc_gamma_from_integral(lims, R_1, small_r_1, gamma_ab, mu_1, m_1, mass_c, &
        rand_root, pos_root, config_a) result(gamma_1)

        implicit none

        real(wp), intent(in) :: lims(4), R_1, small_r_1, gamma_ab, mu_1, m_1, mass_c
        logical, intent(in) :: rand_root, pos_root, config_a

        logical :: positive_root
        real(wp) :: cos_numerator, cos_denominator, gamma_1, temp_R_2_1, temp_R_2_2, &
            R_2_over_m_1, rand_num, R_2_tol, gamma_tol

        R_2_tol = 0.0000001_wp
        gamma_tol = 0.0000001_wp

        if (rand_root) then
            call random_number(rand_num)
            if (rand_num > 0.5_wp) then
                positive_root = .true.
            else
                positive_root = .true.
            end if
        else
            positive_root = pos_root
        end if

        cos_denominator = small_r_1

        ! gamma_1 is undefined if either R_1 = 0 or r_1 = 0
        if (cos_denominator == 0._wp .or. R_1 == 0._wp) then
            gamma_1 = 0._wp
            return
        end if

        ! Calculate R_2 / m_1
        temp_R_2_1 = (cos_denominator / mu_1)**2._wp + (R_1 / mass_c)**2._wp * &
            (cos(gamma_ab)**2._wp - 1._wp)

        ! Must be positive as will be square rooted
        if (temp_R_2_1 < 0._wp) then
            gamma_1 = -100._wp
            return
        else
            temp_R_2_2 = - (R_1 / mass_c) * cos(gamma_ab)
            if (positive_root) then
                R_2_over_m_1 = temp_R_2_2 + sqrt(temp_R_2_1)
            else
                R_2_over_m_1 = temp_R_2_2 - sqrt(temp_R_2_1)
            end if
        end if

        ! Check R_2 within valid limits
        if ( (R_2_over_m_1 < lims(1) / m_1) .or. (R_2_over_m_1 > lims(2) / m_1) ) then

            ! Try opposite root before aborting
            if (positive_root) then
                R_2_over_m_1 = temp_R_2_2 - sqrt(temp_R_2_1)
            else
                R_2_over_m_1 = temp_R_2_2 + sqrt(temp_R_2_1)
            end if

            if ( (R_2_over_m_1 < lims(1) / m_1) &
            .or. (R_2_over_m_1 > lims(2) / m_1) ) then
                if (R_2_over_m_1 > lims(2) / m_1 &
                .and. R_2_over_m_1 < R_2_tol + lims(2) / m_1) then
                    R_2_over_m_1 = lims(2) / m_1
                else if (R_2_over_m_1 < lims(1) / m_1 &
                .and. R_2_over_m_1 > - R_2_tol + lims(1) / m_1) then
                    R_2_over_m_1 = lims(1) / m_1
                else
                    gamma_1 = -100._wp
                    return
                end if
            end if
        end if

        if (config_a) then
            cos_numerator = - mu_1 * ( (R_1 / mass_c) + R_2_over_m_1 * cos(gamma_ab) )
        else
            cos_numerator = mu_1 * ( (R_1 / mass_c) + R_2_over_m_1 * cos(gamma_ab) )
        end if

        gamma_1 = cos_numerator / cos_denominator

        ! Prevent NaN when taking acos - there is probably a better solution!
        if (gamma_1 > 1._wp) then
            if (gamma_1 <= (1._wp+gamma_tol)) then
                gamma_1 = 1._wp
            else
                print *, "Warning: invalid gamma_1: ", gamma_1
                gamma_1 = 1._wp
                stop
            end if
        end if
        if (gamma_1 < -1._wp) then
            if (gamma_1 >= -(1._wp+gamma_tol)) then
                gamma_1 = -1._wp
            else
                print *, "Warning: invalid gamma_1: ", gamma_1
                gamma_1 = -1._wp
                stop
            end if
        end if

        gamma_1 = acos(gamma_1)

    end function calc_gamma_from_integral


    ! Calculate single Jacobi coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b)
    ! Choice of single coordinates is determined by the config_a flag
    ! from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    function transform_mixed_to_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a) &
        result(single_coords)

        implicit none

        real(wp), intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real(wp) :: single_coords(2)

        ! r_a or r_b
        single_coords(1) = calc_small_r_from_mixed(R_1, R_2, gamma_ab, mu_1, m_1, mass_c)

        ! gamma_a or gamma_b
        single_coords(2) = calc_gamma_from_mixed(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a)

    end function transform_mixed_to_single


    ! Calculate mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! from single Jacobi coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b)
    ! Choice of single coordinates is determined by the config_a flag
    function transform_single_to_mixed(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a) &
        result(mixed_coords)

        implicit none

        real(wp), intent(in) :: R_1, small_r_1, gamma_1, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real(wp) :: mixed_coords(2)

        ! R_b or R_a
        mixed_coords(1) = calc_R_from_single(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, config_a)

        ! gamma_b or gamma_a
        mixed_coords(2) = calc_gamma_from_single(R_1, small_r_1, gamma_1, mu_1, m_1, mass_c, &
            config_a)

    end function transform_single_to_mixed


    ! Transform a grid of mixed Jacobi coordinates (R_a, R_b and gamma_ab) to single Jacobi
    ! coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b) or vice versa
    ! Choice of single coordinates is determined by the config_a flag
    ! Direction of transformation is determined by the mixed_int flag
    ! (mixed_int == .true. means input limits and x are written in mixed coordinates)
    ! Grid defined by limits [R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max]
    ! and number of points in each coordinate [n_R_a, n_R_b, n_gamma_ab]
    ! or equivalent for each set of single Jacobi coordinates
    function transform_grid(x, coords, n, mu_1, m_1, mass_c, config_a, mixed_int) &
        result(trans_coords)

        implicit none

        integer, intent(in) :: n(2)
        real(wp), intent(in) :: x, coords(2, n(1)+1, n(2)+1), mu_1, m_1, mass_c
        logical, intent(in) :: config_a, mixed_int

        integer :: i, j, k
        real(wp) :: y, z, trans_coords(2, n(1)+1, n(2)+1)

        ! Loop over R_b (config_a) or R_a (not config_a) values
        do i = 0, n(1)
            y = coords(1, i+1, 1)
            do j = 0, n(2)
                z = coords(2, i+1, j+1)
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


    function get_grid(x, boundary_val_2, n, mixed_lims, mu_1, m_1, mass_c, config_a, mixed_int) &
        result(coords)

        implicit none

        integer, intent(in) :: n(2)
        real(wp), intent(in) :: x, boundary_val_2, mixed_lims(4), mu_1, m_1, mass_c
        logical, intent(in) :: config_a, mixed_int

        integer :: i, j, k
        real(wp) :: y, z, width(2), coords(2, n(1)+1, n(2)+1), lims(4)

        if (mixed_int) then
            ! R_b or R_a, and gamma_ab limits already defined
            lims = mixed_lims
            width(1) = abs(lims(2) - lims(1)) / real(n(1), kind=wp)
            width(2) = abs(lims(4) - lims(3)) / real(n(2), kind=wp)
        else
            ! Calculate r_a or r_b limits
            ! lims(1:2) = get_jacobi_lims(mixed_lims, x, 0._wp, mu_1, m_1, mass_c, .true., .false., &
            !     .true., config_a)
            ! print *, "Estimate", lims(1:2)

            lims(1:2) = get_jacobi_lims(mixed_lims, x, 0._wp, mu_1, m_1, mass_c, .true., .false., &
                .false., config_a)
            ! print *, "Analytical", lims(1:2)

            width(1) = abs(lims(2) - lims(1)) / real(n(1), kind=wp)
        end if

        ! Loop over R_b or R_a values
        do i = 0, n(1)
            y = lims(1) + real(i, kind=wp) * width(1)

            if (.not. mixed_int) then
                ! Calculate gamma_a or gamma_b limits
                lims(3:4) = get_jacobi_lims(mixed_lims, x, y, mu_1, m_1, mass_c, .false., .true., &
                    .true., config_a)
                ! print *, "Estimate", lims(3:4)

                ! lims(3:4) = get_jacobi_lims(mixed_lims, x, y, mu_1, m_1, mass_c, .false., .true., &
                !     .false., config_a)
                ! print *, "Analytical", lims(3:4)
                width(2) = abs(lims(4) - lims(3)) / real(n(2), kind=wp)
            end if

            do j = 0, n(2)
                z = lims(3) + real(j, kind=wp) * width(2)
                coords(:, i+1, j+1) = (/ y, z /)
            end do
        end do

    end function get_grid


    ! Transform surface amplitudes w_ik from mixed to single Jacobi coordinates
    function transform_amps(old_amps, n_1, nc_1, nt, simpson_n, boundary_val_1, &
        boundary_val_2, mu_1, mu_2, m_1, m_2, mass_c, coords, trans_coords, config_a, mixed_int, &
        sub_comm, sub_comm_size, rank) result(amps)

      implicit none

        integer, intent(in) :: n_1, nc_1, nt, simpson_n(2)
        real(wp), intent(in) :: old_amps(n_1, nt), boundary_val_1, boundary_val_2, mu_1, &
            mu_2, m_1, m_2, mass_c, coords(2, simpson_n(1)+1, simpson_n(2)+1), &
            trans_coords(2, simpson_n(1)+1, simpson_n(2)+1)
        logical, intent(in) :: config_a, mixed_int

        ! MPI variables
        integer, intent(in) :: sub_comm, sub_comm_size, rank

        integer :: i, k, n, m, imin, imax, progress, num_values

        real(wp) :: integral_1(n_1), integral_2(nt), amps(n_1, nt)
        logical :: verbose, channel_func_1, single_func_1, channel_func_2, single_func_2

        ! MPI variables
        integer :: ierr, recv_status(MPI_STATUS_SIZE)

        ! verbose = .true.
        verbose = .false.

        ! Initialise amps array
        do i = 1, n_1
            do k = 1, nt
                amps(i, k) = 0._wp
            end do
        end do

        ! Divide calculation of amps (i runs between 1 and n_1 overall)
        imin = 1 + rank * n_1 / sub_comm_size
        if (rank == sub_comm_size - 1) then
            imax = n_1
        else
            imax = (rank + 1) * n_1 / sub_comm_size
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

            do k = 1, nt
                ! Loop over channel functions (phi or theta)
                ! Take inner product of channel function i in single basis
                ! with channel function n in mixed basis
                channel_func_1 = .true.
                single_func_1 = .true.
                channel_func_2 = .true.
                single_func_2 = .false.
                do n = 1, n_1
                    integral_1(n) = integrate_double_Simpson(simpson_n, boundary_val_1, &
                        boundary_val_2, n_1, nc_1, i, n, mu_1, mu_2, m_1, m_2, mass_c, &
                        coords, trans_coords, config_a, mixed_int, channel_func_1, single_func_1, &
                        channel_func_2, single_func_2)
                end do
                ! Loop over basis functions (psi)
                ! Take inner product of basis function m in mixed basis
                ! with basis function k in single basis
                channel_func_1 = .false.
                single_func_1 = .false.
                channel_func_2 = .false.
                single_func_2 = .true.
                do m = 1, nt
                    integral_2(m) = integrate_double_Simpson(simpson_n, boundary_val_1, &
                        boundary_val_2, n_1, nc_1, m, k, mu_1, mu_2, m_1, m_2, mass_c, &
                        coords, trans_coords, config_a, mixed_int, channel_func_1, single_func_1, &
                        channel_func_2, single_func_2)
                end do

                do n = 1, n_1
                    do m = 1, nt
                        ! print *, integral_1(n), integral_2(m), old_amps(n, m), amps(i, k)
                        amps(i, k) = amps(i, k) + integral_1(n) * integral_2(m) * old_amps(n, m)
                    end do
                end do
            end do
        end do

        ! Send all amps to rank 0
        if (rank == 0) then
            ! Receive from all other ranks
            do i = 1, sub_comm_size - 1
                    imin = 1 + i * n_1 / sub_comm_size
                    if (i == sub_comm_size - 1) then
                            imax = n_1
                    else
                            imax = (i + 1) * n_1 / sub_comm_size
                    end if

                    num_values = (imax - imin + 1) * nt
                    call MPI_Recv(amps(imin:imax, :), num_values, &
                            MPI_DOUBLE_PRECISION, i, 0, sub_comm, recv_status, ierr)
            end do
        else
            ! All non-zero ranks send to rank 0
            num_values = (imax - imin + 1) * nt
            call MPI_Ssend(amps(imin:imax, :), num_values, MPI_DOUBLE_PRECISION, 0, 0, sub_comm, ierr)
        end if

        ! All sends/receives complete
        call MPI_Barrier(sub_comm, ierr)

    end function transform_amps


    ! Calculate surface amplitudes at boundary
    ! Currently calculated through integration, but for appropriate channel functions
    ! should be calculable directly from radial functions at the boundary and coefficients
    function calc_old_amps(n, nc, nt, simpson_n, coords, trans_coords, boundary_val_1, &
        boundary_val_2, mu_1, mu_2, m_1, m_2, mass_c, config_a, mixed_int, sub_comm, &
        sub_comm_size, rank) result(amps)

        implicit none

        integer, intent(in) :: n, nc, nt, simpson_n(2)
        real(wp), intent(in) :: coords(2, simpson_n(1)+1, simpson_n(2)+1), &
            trans_coords(2, simpson_n(1)+1, simpson_n(2)+1), boundary_val_1, boundary_val_2, &
            mu_1, mu_2, m_1, m_2, mass_c
        logical, intent(in) :: config_a, mixed_int

        ! MPI variables
        integer, intent(in) :: sub_comm, sub_comm_size, rank

        real(wp) :: amps(n, nt)
        integer :: i, k, imin, imax, progress, num_values
        logical :: verbose, channel_func_1, single_func_1, channel_func_2, single_func_2

        ! MPI variables
        integer :: ierr, recv_status(MPI_STATUS_SIZE)

        ! verbose = .true.
        verbose = .false.

        ! Take inner product of channel function i in the mixed basis
        ! with basis function k in mixed basis
        channel_func_1 = .true.
        single_func_1 = .false.
        channel_func_2 = .false.
        single_func_2 = .false.

        ! Initialise old amps
        do i = 1, n
            do k = 1, nt
                amps(i, k) = 0._wp
            end do
        end do

        ! Divide calculation of amps (i runs between 1 and n overall)
        imin = 1 + rank * n / sub_comm_size
        if (rank == sub_comm_size - 1) then
            imax = n
        else
            imax = (rank + 1) * n / sub_comm_size
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

            do k = 1, nt
                amps(i, k) = integrate_double_Simpson(simpson_n, boundary_val_1, boundary_val_2, &
                    n, nc, i, k, mu_1, mu_2, m_1, m_2, mass_c, coords, trans_coords, config_a, &
                    mixed_int, channel_func_1, single_func_1, channel_func_2, single_func_2)
                amps(i, k) = amps(i, k) / boundary_val_1
            end do
        end do

        ! Send all amps to rank 0
        if (rank == 0) then
            ! Receive from all other ranks
            do i = 1, sub_comm_size - 1
                    imin = 1 + i * n / sub_comm_size
                    if (i == sub_comm_size - 1) then
                            imax = n
                    else
                            imax = (i + 1) * n / sub_comm_size
                    end if

                    num_values = (imax - imin + 1) * nt
                    call MPI_Recv(amps(imin:imax, :), num_values, &
                            MPI_DOUBLE_PRECISION, i, 0, sub_comm, recv_status, ierr)
            end do
        else
            ! All non-zero ranks send to rank 0
            num_values = (imax - imin + 1) * nt
            call MPI_Ssend(amps(imin:imax, :), num_values, MPI_DOUBLE_PRECISION, 0, 0, sub_comm, ierr)
        end if

        ! All sends/receives complete
        call MPI_Barrier(sub_comm, ierr)

    end function calc_old_amps


    ! Transform surface amplitudes w_ik from mixed to single Jacobi coordinates
    subroutine transform_all_amps(n_1, nc_1, n_2, nc_2, nt, simpson_n, mixed_lims, mu_a, mu_b, &
        m_a, m_b, mass_c, mixed_int, sub_comm)

      implicit none

        integer, intent(in) :: n_1, nc_1, n_2, nc_2, nt, simpson_n(3)
        real(wp), intent(in) :: mixed_lims(6), mu_a, mu_b, m_a, m_b, mass_c
        logical, intent(in) :: mixed_int

        ! MPI variables
        integer, intent(in) :: sub_comm

        integer :: i, j, simpson_n_1(2)
        real(wp) :: old_amps(n_1 + n_2, nt), new_amps(n_1 + n_2, nt), &
            coords_a(2, simpson_n(2)+1, simpson_n(3)+1), &
            trans_coords_a(2, simpson_n(2)+1, simpson_n(3)+1), &
            coords_b(2, simpson_n(1)+1, simpson_n(3)+1), &
            trans_coords_b(2, simpson_n(1)+1, simpson_n(3)+1), boundary_val_1, boundary_val_2, &
            lims_1(4)
        logical :: config_a, verbose
        double precision :: time(9), t_diff

        ! MPI variables
        integer :: sub_comm_size, rank, ierr

        call MPI_Comm_size(sub_comm, sub_comm_size, ierr)
        call MPI_Comm_rank(sub_comm, rank, ierr)

        ! Initialise old amps
        do i = 1, n_1 + n_2
            do j = 1, nt
                old_amps(i, j) = 0._wp
            end do
        end do

        ! Initialise transformed amps
        do i = 1, n_1 + n_2
            do j = 1, nt
                new_amps(i, j) = 0._wp
            end do
        end do

        ! verbose = .true.
        verbose = .false.

        ! At the R_a boundary
        config_a = .true.
        simpson_n_1 = simpson_n(2:3)
        boundary_val_1 = mixed_lims(2)
        boundary_val_2 = mixed_lims(4)
        lims_1 = mixed_lims(3:6)

        if (rank == 0 .and. verbose) then
            print *, "Calculating grid points in chosen basis for first arrangement..."
        end if

        ! Start the timer before grid points calculated
        call MPI_Barrier(sub_comm, ierr)
        time(1) = MPI_Wtime()

        ! Calculate all grid points for Simpson's integration of arrangement A
        coords_a = get_grid(boundary_val_1, boundary_val_2, simpson_n_1, lims_1, mu_a, m_b, &
            mass_c, config_a, mixed_int)

        ! Check the timer after grid points calculated, before grid transformation
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(2) = MPI_Wtime()
            t_diff = time(2) - time(1)
            print *, "Grid calculation time: ", t_diff
        end if

        if (rank == 0 .and. verbose) then
            print *, "Transforming grid points for first arrangement..."
        end if

        ! If mixed_int, calculate single Jacobi coordinates (R_a, r_a, gamma_a)
        ! for all mixed coordinates, else calculate mixed Jacobi coordinates (R_a, R_b, gamma_ab)
        ! for all single coordinates
        trans_coords_a = transform_grid(boundary_val_1, coords_a, simpson_n_1, mu_a, m_b, mass_c, &
            config_a, mixed_int)

        ! Check the timer after grid points transformation, before amplitudes calculated
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(3) = MPI_Wtime()
            t_diff = time(3) - time(2)
            print *, "Grid transformation time: ", t_diff
        end if

        if (rank == 0 .and. verbose) then
            print *, "Calculating surface amplitudes in mixed Jacobi coordinates for first &
                arrangement..."
        end if

        ! Calculate amplitudes in mixed Jacobi coordinates at R_a boundary
        old_amps(:n_1, :) = calc_old_amps(n_1, nc_1, nt, simpson_n_1, coords_a, trans_coords_a, &
            boundary_val_1, boundary_val_2, mu_a, mu_b, m_b, m_a, mass_c, config_a, mixed_int, &
            sub_comm, sub_comm_size, rank)

        ! Check the timer after surface amplitude calculation, before amplitudes transformed
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(4) = MPI_Wtime()
            t_diff = time(4) - time(3)
            print *, "Surface amplitude calculation time: ", t_diff
        end if

        if (rank == 0 .and. verbose) then
            print *, "Transforming surface amplitudes into single Jacobi coordinates for first &
            arrangement..."
        end if

        ! Transform surface amplitudes at R_a boundary
        new_amps(:n_1, :) = transform_amps(old_amps(:n_1, :), n_1, nc_1, nt, simpson_n_1, &
            boundary_val_1, boundary_val_2, mu_a, mu_b, m_b, m_a, mass_c, coords_a, &
            trans_coords_a, config_a, mixed_int, sub_comm, sub_comm_size, rank)

        ! Check the timer after surface amplitude transformation, before grid points calculated
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(5) = MPI_Wtime()
            t_diff = time(5) - time(4)
            print *, "Surface amplitude transformation time: ", t_diff
        end if

        ! At the R_b boundary
        config_a = .false.
        simpson_n_1(1) = simpson_n(1)
        simpson_n_1(2) = simpson_n(3)
        boundary_val_1 = mixed_lims(4)
        boundary_val_2 = mixed_lims(2)
        lims_1(1:2) = mixed_lims(1:2)
        lims_1(3:4) = mixed_lims(5:6)

        if (rank == 0 .and. verbose) then
            print *, "Calculating grid points in chosen basis for second arrangement..."
        end if

        ! Calculate all grid points for Simpson's integration of arrangement B
        coords_b = get_grid(boundary_val_1, boundary_val_2, simpson_n_1, lims_1, mu_b, m_a, &
            mass_c, config_a, mixed_int)

        ! Check the timer after grid points calculated, before grid transformation
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(6) = MPI_Wtime()
            t_diff = time(6) - time(5)
            print *, "Grid calculation time: ", t_diff
        end if

        if (rank == 0 .and. verbose) then
            print *, "Transforming grid points for second arrangement..."
        end if

        ! If mixed_int, calculate single Jacobi coordinates (R_b, r_b, gamma_b)
        ! for all mixed coordinates, else calculate mixed Jacobi coordinates (R_a, R_b, gamma_ab)
        ! for all single coordinates
        trans_coords_b = transform_grid(boundary_val_1, coords_b, simpson_n_1, mu_b, m_a, mass_c, &
            config_a, mixed_int)

        ! Check the timer after grid points transformation, before amplitudes calculated
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(7) = MPI_Wtime()
            t_diff = time(7) - time(6)
            print *, "Grid transformation time: ", t_diff
        end if

        if (rank == 0 .and. verbose) then
            print *, "Calculating surface amplitudes in mixed Jacobi coordinates for second &
            arrangement..."
        end if

        ! Calculate amplitudes in mixed Jacobi coordinates at R_b boundary
        old_amps(n_1+1:, :) = calc_old_amps(n_2, nc_2, nt, simpson_n_1, coords_b, trans_coords_b, &
            boundary_val_1, boundary_val_2, mu_b, mu_a, m_a, m_b, mass_c, config_a, mixed_int, &
            sub_comm, sub_comm_size, rank)

        ! Check the timer after surface amplitude calculation, before amplitudes transformed
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(8) = MPI_Wtime()
            t_diff = time(8) - time(7)
            print *, "Surface amplitude calculation time: ", t_diff
        end if

        if (rank == 0 .and. verbose) then
            print *, "Transforming surface amplitudes into single Jacobi coordinates for first &
                arrangement..."
        end if

        ! Transform surface amplitudes at R_b boundary
        new_amps(n_1+1:, :) = transform_amps(old_amps(n_1+1:, :), n_2, nc_2, nt, simpson_n_1, &
            boundary_val_1, boundary_val_2, mu_b, mu_a, m_a, m_b, mass_c, coords_b, &
            trans_coords_b, config_a, mixed_int, sub_comm, sub_comm_size, rank)

        ! Check the timer after surface amplitude transformation, before end of subroutine
        call MPI_Barrier(sub_comm, ierr)
        if (rank == 0 .and. verbose) then
            time(9) = MPI_Wtime()
            t_diff = time(9) - time(8)
            print *, "Surface amplitude transformation time: ", t_diff
        end if

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