program transform

    use mpi

    implicit none

    integer :: simpson_n(3), i, k, num_channels, num_basis_funcs, nc1
    real :: lims(6), integral, mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b, &
        values(9)
    logical :: save_inputs

    ! MPI variables
    integer :: comm, rank, comm_size, ierr, group, sub_group, sub_ranks(1), sub_comm
    double precision :: time(3), t_diff

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
    simpson_n = (/ 5, 5, 5 /)

    ! Limits of integration (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
    ! At R_a boundary, R_a set to R_a_max and remaining four limits used
    ! At R_b boundary, R_b set to R_b_max and remaining four limits used
    ! Note: 4.*atan(1.) = pi
    lims = (/ 0., 1., 0., 2., 0., 4.*atan(1.) /)

    ! Calculate mass relations (three masses defined within)
    call calc_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

    if (rank == 0 .and. save_inputs) then
            values = (/ mass_a, mass_b, mass_c, lims(1), lims(2), lims(3), lims(4), lims(5), &
                lims(6) /)
            open (unit=41, file='outputs/values', form='unformatted')
            write(41) values
            close (unit=41)
    end if

    ! Check the timer before integration
    call MPI_Barrier(comm, ierr)
    time(2) = MPI_Wtime()

    num_channels = 5
    num_basis_funcs = 5
    nc1 = 3

    call transform_all_amps(num_channels, num_basis_funcs, nc1, simpson_n, lims, mu_a, mu_b, m_a, &
        m_b, mass_c, comm)

    ! Check the timer after transformation, at end of program
    call MPI_Barrier(comm, ierr)
    time(3) = MPI_Wtime()
    if (rank == 0) then
            t_diff = time(3) - time(2)
            print *, "Integration time: ", t_diff

            t_diff = time(3) - time(1)
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


    function get_single_phi(i, x, r, gamma, config_a) result(func)

        implicit none

        integer, intent(in) :: i
        real, intent(in) :: x, r, gamma
        logical, intent(in) :: config_a

        real :: func

        func = 1.

    end function get_single_phi


    function get_mixed_phi(n, x, y, z, config_a) result(func)

        implicit none

        integer, intent(in) :: n
        real, intent(in) :: x, y, z
        logical, intent(in) :: config_a

        real :: func

        func = 1.

    end function get_mixed_phi


    function get_single_psi(m, x, r, gamma, config_a) result(func)

        implicit none

        integer, intent(in) :: m
        real, intent(in) :: x, r, gamma
        logical, intent(in) :: config_a

        real :: func

        func = 1.

    end function get_single_psi


    function get_mixed_psi(k, x, y, z, config_a) result(func)

        implicit none

        integer, intent(in) :: k
        real, intent(in) :: x, y, z
        logical, intent(in) :: config_a

        real :: func

        func = 1.

    end function get_mixed_psi


    ! Calculates total integrand
    ! For config_a x = R_a, y = R_b and z = gamma_ab
    ! For config_a x = R_b, y = R_a and z = gamma_ab
    ! x is not integrated over, but the integrand may be a function of x
    function integrand_func(x, y, z, mu_a, mass_c, m_b, single_coords, i, j, config_a, &
        channel_func) result(total_integrand)

            implicit none

            integer, intent(in) :: i, j
            real, intent(in) :: x, y, z, mu_a, mass_c, m_b, single_coords(2)
            logical, intent(in) :: config_a, channel_func

            real :: total_integrand, integrand_volume, integrand_1, integrand_2, r, gamma

            ! Single Jacobi coordinates
            r = single_coords(1)
            gamma = single_coords(2)

            ! Volume component to be integrated
            integrand_volume = y**2. * sin(z)

            ! Function to be integrated:
            if (channel_func) then
                integrand_1 = get_single_phi(i, x, r, gamma, config_a)
                integrand_2 = get_mixed_phi(j, x, y, z, config_a)
            else
                integrand_1 = get_single_psi(i, x, r, gamma, config_a)
                integrand_2 = get_mixed_psi(j, x, y, z, config_a)
            end if
            total_integrand = integrand_volume * integrand_1 * integrand_2

    end function integrand_func


    ! Use Simpson's rule to integrate f(x,y,z) over two variables: y and z
    ! x is input at a fixed value
    function integrate_double_Simpson(n, lims, x, mu_1, m_1, mass_c, single_coords, config_a, &
        channel_func) result(total_integral)

            implicit none

            integer, intent(in) :: n(2)
            real, intent(in) :: lims(4), x, mu_1, m_1, mass_c, &
                single_coords(2, n(1)+1, n(2)+1)
            logical, intent(in) :: config_a, channel_func

            real :: width(2), y, z, temp_integral, z_integral, total_integral, test_coords(3)
            integer :: i, j, progress
            logical :: verbose

            verbose = .true.

            total_integral = 0.
            temp_integral = 0.

            width(1) = abs(lims(2) - lims(1)) / real(n(1))

            ! print *, "coords", single_coords

            ! Integrate each variable in turn, covering full limits
            ! Variables labelled x, y and z for simplicity
            do i = 0, n(1)

                    ! Set value for R_b at fixed R_a
                    y = lims(1) + real(i) * width(1)

                    ! Total gamma_ab integral given R_a and R_b
                    z_integral = 0.

                    do j = 0, n(2)

                            ! Set value for gamma_ab
                            width(2) = abs(lims(4) - lims(3)) / real(n(2))
                            z = lims(3) + real(j) * width(2)

                            ! test_coords = transform_mixed_to_single(x, y, z, mu_1, m_1, &
                            !     mass_c, .true.)
                            ! print *, "TEST", test_coords
                            ! print *, "INPUT", single_coords(:, i+1, j+1)

                            ! Evaluate integral at this point
                            temp_integral = integrand_func(x, y, z, mu_a, mass_c, m_b, &
                                single_coords(:, i+1, j+1), i+1, j+1, config_a, channel_func)

                            if (j == 0 .or. j == n(2)) then
                                    z_integral = z_integral + temp_integral
                            else if (mod(j, 2) == 0) then
                                    z_integral = z_integral + 2. * temp_integral
                            else
                                    z_integral = z_integral + 4. * temp_integral
                            end if

                    end do

                    ! Total gamma_ab, gamma_a or z intergral for set
                    ! (R_a, R_b), (R_a, r_a) or (x, y)
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

    end function integrate_double_Simpson


    ! Calculate r_a or r_b from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! For r_a, R_1 = R_a, R_2 = R_b, mu_1 = mu_a and m = m_b
    ! For r_b, R_1 = R_b and R_2 = R_a, mu_1 = mu_b and m = m_a
    function calc_r(R_1, R_2, gamma_ab, mu_1, m_1, mass_c) result(r)

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

    end function calc_r


    ! Calculate gamma_a or gamma_bfrom mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    ! For gamma_a, R_1 = R_a, R_2 = R_b, mu_1 = mu_a and m = m_b
    ! For gamma_b, R_1 = R_b and R_2 = R_a, mu_1 = mu_b and m = m_a
    function calc_gamma(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, calc_gamma_a) result(gamma)

        implicit none

        real, intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c
        logical, intent(in) :: calc_gamma_a

        real :: r, gamma, cos_numerator, err

        err = 0.000001

        if (R_1 == 0) then
            gamma = 0
            return
        end if

        r = calc_r(R_1, R_2, gamma_ab, mu_1, m_1, mass_c)

        if (r == 0.) then
            gamma = 0
            return
        end if

        ! Multiply by sgn(gamma_ba) for gamma_a, sgn(gamma_ab) fo gamma_b
        if (calc_gamma_a) then
            cos_numerator = - mu_1 * ( (R_1 / mass_c) + (R_2 / m_1) * cos(gamma_ab) )
        else
            cos_numerator = mu_1 * ( (R_1 / mass_c) + (R_2 / m_1) * cos(gamma_ab) )
        end if

        gamma = cos_numerator / r

        ! Temporary to prevent NaN when taking acos - there is probably a better solution!
        if (gamma > 1. .and. gamma <= (1.+err)) then
            gamma = 1.
        end if
        if (gamma < -1. .and. gamma >= -(1.+err)) then
            gamma = -1.
        end if

        gamma = acos(gamma)

    end function calc_gamma


    ! Calculate single Jacobi coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b)
    ! Choice of single coordinates is determined by the config_a flag
    ! from mixed Jacobi coordinates (R_a, R_b and gamma_ab)
    function transform_mixed_to_single(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a) &
        result(single_coords)

        implicit none

        real, intent(in) :: R_1, R_2, gamma_ab, mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        real :: single_coords(3)

        ! R_a or R_b
        single_coords(1) = R_1

        ! r_a or r_b
        single_coords(2) = calc_r(R_1, R_2, gamma_ab, mu_1, m_1, mass_c)

        ! gamma_a or gamma_b
        single_coords(3) = calc_gamma(R_1, R_2, gamma_ab, mu_1, m_1, mass_c, config_a)

    end function transform_mixed_to_single


    ! Transform a grid of mixed (R_a, R_b and gamma_ab) to
    ! single Jacobi coordinates (R_a, r_a and gamma_a) or (R_b, r_b and gamma_b)
    ! Choice of single coordinates is determined by the config_a flag
    ! Grid defined by limits [R_a_min, R_a_max, R_b_min, R_b_max, gammab_a_min, gamma_ab_max]
    ! and number of points in each coordinate [n_R_a, n_R_b, n_gamma_ab]
    function transform_grid(lims, n, mu_1, m_1, mass_c, config_a) result(single_coords)

        implicit none

        integer, intent(in) :: n(3)
        real, intent(in) :: lims(6), mu_1, m_1, mass_c
        logical, intent(in) :: config_a

        integer :: i, j, k
        real :: R_a(n(1)+1), R_b(n(2)+1), gamma_ab(n(3)+1), width(3), &
            single_coords(3, n(1)+1, n(2)+1, n(3)+1)

        ! Spacing between points for each coordinate
        if (config_a) then
            width(2) = abs(lims(4) - lims(3)) / real(n(2))
            width(3) = abs(lims(6) - lims(5)) / real(n(3))
        else
            width(1) = abs(lims(2) - lims(1)) / real(n(1))
            width(3) = abs(lims(6) - lims(5)) / real(n(3))
        end if

        ! Loop over R_a values
        do i = 0, n(1)

            if (config_a) then
                R_a(i+1) = lims(2)
            else
                R_a(i+1) = lims(1) + real(i) * width(1)
            end if

            ! Loop over R_b values
            do j = 0, n(2)

                if (config_a) then
                    R_b(j+1) = lims(3) + real(j) * width(2)
                else
                    R_b(j+1) = lims(4)
                end if

                ! Loop over gamma_ab values
                do k = 0, n(3)

                    gamma_ab(k+1) = lims(5) + real(k) * width(3)

                    ! Calculate single Jacobi coordinates at these mixed coordinates
                    if (config_a) then
                        single_coords(:, i+1, j+1, k+1) = transform_mixed_to_single(R_a(i+1), &
                            R_b(j+1), gamma_ab(k+1), mu_1, m_1, mass_c, config_a)
                    else
                        single_coords(:, i+1, j+1, k+1) = transform_mixed_to_single(R_b(j+1), &
                            R_a(i+1), gamma_ab(k+1), mu_1, m_1, mass_c, config_a)
                    end if

                end do
            end do
        end do

    end function transform_grid


    ! Transform surface amplitudes w_ik from mixed to single Jacobi coordinates
    function transform_amps(old_amps, num_channels, num_basis_funcs, simpson_n, lims, mu_1, m_1, &
        mass_c, single_coords, config_a, sub_comm) result(amps)

        integer, intent(in) :: num_channels, num_basis_funcs, simpson_n(3)
        real, intent(in) :: old_amps(num_channels, num_basis_funcs), lims(6), mu_1, m_1, mass_c, &
            single_coords(3, simpson_n(1)+1, simpson_n(2)+1, simpson_n(3)+1)
        logical, intent(in) :: config_a

        ! MPI variables
        integer, intent(in) :: sub_comm

        integer :: i, k, n, m, imin, imax, progress, num_values, simpson_n_1(2)
        real :: boundary_val, integral_1(num_channels), integral_2(num_basis_funcs), &
            amps(num_channels, num_basis_funcs), lims_1(4)
        logical :: verbose, channel_func

        ! MPI variables
        integer :: sub_comm_size, ierr, request, recv_status(MPI_STATUS_SIZE)

        call MPI_Comm_size(sub_comm, sub_comm_size, ierr)

        verbose = .true.

        if (config_a) then
            boundary_val = lims(2)
            simpson_n_1(1) = simpson_n(2)
            simpson_n_1(2) = simpson_n(3)
            lims_1(1:2) = lims(3:4)
            lims_1(3:4) = lims(5:6)
        else
            ! R_b boundary
            boundary_val = lims(4)
            simpson_n_1(1) = simpson_n(1)
            simpson_n_1(2) = simpson_n(3)
            lims_1(1:2) = lims(1:2)
            lims_1(3:4) = lims(5:6)
        end if

        ! Initialise amps array
        do i = 1, num_channels
            do k = 1, num_basis_funcs
                amps(i, k) = 0.
            end do
        end do

        ! Divide calculation of amps (i runs between 1 and num_channels overall)
        imin = 1 + rank * num_channels / sub_comm_size
        if (rank == sub_comm_size - 1) then
                imax = num_channels
        else
                imax = (rank + 1) * num_channels / sub_comm_size
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
                ! Loop over channel functions (phi)
                channel_func = .true.
                do n = 1, num_channels
                    if (config_a) then
                        integral_1(n) = integrate_double_Simpson(simpson_n_1, lims_1, &
                            boundary_val, mu_1, m_1, mass_c, single_coords(2:3, 1, :, :), &
                            config_a, channel_func)
                    else
                        integral_1(n) = integrate_double_Simpson(simpson_n_1, lims_1, &
                            boundary_val, mu_1, m_1, mass_c, single_coords(2:3, :, 1, :), &
                            config_a, channel_func)
                    end if
                end do
                ! Loop over basis functions (psi)
                channel_func = .false.
                do m = 1, num_basis_funcs
                    if (config_a) then
                        integral_2(m) = integrate_double_Simpson(simpson_n_1, lims_1, &
                            boundary_val, mu_1, m_1, mass_c, single_coords(:, 1, :, :), &
                            config_a, channel_func)
                    else
                        integral_2(m) = integrate_double_Simpson(simpson_n_1, lims_1, &
                            boundary_val, mu_1, m_1, mass_c, single_coords(:, 1, :, :), &
                            config_a, channel_func)
                    end if
            end do

                do n = 1, num_channels
                    do m = 1, num_basis_funcs
                        amps(i, k) = amps(i, k) + integral_1(n) * integral_2(m) * old_amps(n, m)
                    end do
                end do
                ! print *, "amp(", i, k, ") = ", amps(i, k)
            end do
        end do

        ! Send all amps to rank 0
        if (rank == 0) then
                ! Receive from all other ranks
                do i = 1, sub_comm_size - 1
                        imin = 1 + i * num_channels / sub_comm_size
                        if (i == sub_comm_size - 1) then
                                imax = num_channels
                        else
                                imax = (i + 1) * num_channels / sub_comm_size
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


    ! Placeholder - will read surface amplitudes calculated from inner region
    ! Currently sets all values to 1.
    function get_old_amps(num_channels, num_basis_funcs) result(amps)

        implicit none

        integer, intent(in) :: num_channels, num_basis_funcs

        real :: amps(num_channels, num_basis_funcs)
        integer :: i, j

        do i = 1, num_channels
            do j = 1, num_basis_funcs
                amps(i, j) = 1.
            end do
        end do

    end function get_old_amps


    ! Transform surface amplitudes w_ik from mixed to single Jacobi coordinates
    subroutine transform_all_amps(num_channels, num_basis_funcs, nc1, simpson_n, lims, mu_a, &
        mu_b,m_a, m_b, mass_c, sub_comm)

        integer, intent(in) :: num_channels, num_basis_funcs, simpson_n(3)
        real, intent(in) :: lims(6), mu_a, mu_b, m_a, m_b, mass_c

        ! MPI variables
        integer, intent(in) :: sub_comm

        integer :: i, j, nc1, nc2, simpson_n_a(3), simpson_n_b(3)
        real :: old_amps(num_channels, num_basis_funcs), new_amps(num_channels, num_basis_funcs), &
            single_coords_a(3, 1, simpson_n(2)+1, simpson_n(3)+1), &
            single_coords_b(3, simpson_n(1)+1, 1, simpson_n(3)+1)
        logical :: config_a

        old_amps = get_old_amps(num_channels, num_basis_funcs)
        nc2 = num_channels - nc1

        ! Initialise transformed amps
        do i = 1, num_channels
            do j = 1, num_basis_funcs
                new_amps(i, j) = 0.
            end do
        end do

        ! At R_a boundary
        config_a = .true.
        simpson_n_a = simpson_n
        simpson_n_a(1) = 0

        ! Calculate single Jacobi coordinates (R_a, r_a, gamma_a) for all mixed coordinates
        single_coords_a = transform_grid(lims, simpson_n_a, mu_a, m_b, mass_c, config_a)

        ! Transform surface amplitudes at R_a boundary
        new_amps(:nc1, :) = transform_amps(old_amps, nc1, num_basis_funcs, simpson_n_a, &
            lims, mu_a, m_b, mass_c, single_coords_a, config_a, sub_comm)


        ! At R_b boundary
        config_a = .false.
        simpson_n_b = simpson_n
        simpson_n_b(2) = 0

        ! Calculate single Jacobi coordinates (R_b, r_b, gamma_b) for all mixed coordinates
        single_coords_b = transform_grid(lims, simpson_n_b, mu_b, m_a, mass_c, config_a)

        ! Transform surface amplitudes at R_b boundary
        new_amps(nc1+1:, :) = transform_amps(old_amps, nc2, num_basis_funcs, simpson_n_b, &
            lims, mu_b, m_a, mass_c, single_coords_b, config_a, sub_comm)

        ! Print transformed amplitudes - save in future
        if (rank == 0) then
            print *, new_amps
        end if

    end subroutine transform_all_amps

end program transform