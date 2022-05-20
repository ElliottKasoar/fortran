module precisn
        implicit none
        private
        integer, parameter::   wp = selected_real_kind(12)
        public wp
    end module precisn

program jacobi

    use mpi
    use precisn, only: wp

    implicit none

    integer :: mixed_jacobi_n(3), single_jacobi_n(3), mc_n, test_coord_n(3)
    real(wp) :: mixed_jacobi_lims(6), R_a_integral, R_b_integral, gamma_ab_integral, &
        mixed_integral, mixed_Simpson_integral, single_Simpson_integral, &
        single_MC_integral, mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b, &
        test_Simpson_integral, test_coord_lims(6), values(9)
    logical :: save_inputs

    ! MPI variables
    integer :: comm, rank, comm_size, ierr, group, sub_group, sub_ranks(1), sub_comm
    double precision :: time(7), t_diff

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

    ! Number of points to generate for MC integration
    mc_n = 10000

    ! Number of Simpson cells to split mixed Jacobi ingegral into (R_a_n, R_b_n, gamma_ab_n)
    mixed_jacobi_n = (/ 400, 400, 400 /)

    ! Number of Simpson cells to split single Jacobi integral into (R_a_n, r_a_n, gamma_a_n)
    single_jacobi_n = (/ 500, 500, 500 /)

    ! Number of Simpson cells to split test coordinate integral into (x_n, y_n, x_n)
    test_coord_n = (/ 1000, 1000, 1000 /)

    ! Limits of integration (R_a_min, R_a_max, R_b_min, R_b_max, gamma_ab_min, gamma_ab_max)
    ! Note: 4._wp*atan(1._wp) = pi
    mixed_jacobi_lims = (/ 0._wp, 3._wp, 0._wp, 5._wp, 0._wp, 4._wp*atan(1._wp) /)

    ! Limits of integration (x_min, x_max, y_min, y_max, z_min, z_max)
    test_coord_lims = (/ 0._wp, 1._wp, 0._wp, 1._wp, 0._wp, 1._wp /)

    ! Calculate mass relations (three masses defined within)
    call calculate_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

    if (rank == 0 .and. save_inputs) then
        values = (/ mass_a, mass_b, mass_c, mixed_jacobi_lims(1), mixed_jacobi_lims(2), &
            mixed_jacobi_lims(3), mixed_jacobi_lims(4), mixed_jacobi_lims(5), &
            mixed_jacobi_lims(6) /)
        open (unit=41, file='outputs/values', form='unformatted')
        write(41) values
        close (unit=41)
    end if

    ! Check the timer before separated mixed integration
    call MPI_Barrier(comm, ierr)
    time(2) = MPI_Wtime()

    ! Calculate separable integrals for R_a, R_b and gamma_ab using Simpson's rule
    if (rank == 0) then
        R_a_integral = integrate_single_Simpson(mixed_jacobi_n(1), mixed_jacobi_lims(1), &
            mixed_jacobi_lims(2), .true.)
        R_b_integral = integrate_single_Simpson(mixed_jacobi_n(2), mixed_jacobi_lims(3), &
            mixed_jacobi_lims(4), .true.)
        gamma_ab_integral = integrate_single_Simpson(mixed_jacobi_n(3), &
            mixed_jacobi_lims(5), mixed_jacobi_lims(6), .false.)

        ! Integrals are separable so multiply for total integral
        mixed_integral = R_a_integral * R_b_integral * gamma_ab_integral
        print *, "Separated mixed Jacobi integral using Simpson's rule = ", mixed_integral
    end if

    ! Check the timer after separated mixed integration, before unseparated mixed integration
    call MPI_Barrier(comm, ierr)
    time(3) = MPI_Wtime()
    if (rank == 0) then
        t_diff = time(3) - time(2)
        print *, "Integration time: ", t_diff
    end if

    mixed_Simpson_integral = integrate_triple_Simpson(mixed_jacobi_n, mixed_jacobi_lims, mu_a, &
        mu_b, m_a, m_b, mass_c, .true., .false., .false., comm)
    if (rank == 0) then
        ! Calculate integral for R_a, R_b and gamma_ab using Simpson's rule on each in turn
        print *, "Unseparated mixed Jacobi integral using Simpson's rule = ", &
            mixed_Simpson_integral
    end if

    ! Check the timer after unseparated mixed integration, before single Simpson integration
    call MPI_Barrier(comm, ierr)
    time(4) = MPI_Wtime()
    if (rank == 0) then
        t_diff = time(4) - time(3)
        print *, "Integration time: ", t_diff
    end if

    ! Calculate integral for R_a, r_a and gamma_a using Simpson's rule on each in turn
    single_Simpson_integral = integrate_triple_Simpson(single_jacobi_n, mixed_jacobi_lims, &
        mu_a, mu_b, m_a, m_b, mass_c, .false., .true., .false., comm)
    if (rank == 0) then
        print *, "Single Jacobi integral using Simpson's rule = ", single_Simpson_integral
    end if

    ! Check the timer after single Simpson integration, before single MC integration
    call MPI_Barrier(comm, ierr)
    time(5) = MPI_Wtime()
    if (rank == 0) then
        t_diff = time(5) - time(4)
        print *, "Integration time: ", t_diff
    end if

    ! Calculate integral for R_a, r_a and gamma_a using Monte Carlo integration
    ! single_MC_integral = integrate_MC(mc_n, mixed_jacobi_lims, mu_a, mu_b, m_a, m_b, &
    !     mass_c, comm)
    ! if (rank == 0) then
    !     print *, "Single Jacobi integral using Monte Carlo = ", single_MC_integral
    ! end if

    ! Check the timer after single MC integration, before test integration
    call MPI_Barrier(comm, ierr)
    time(6) = MPI_Wtime()
    if (rank == 0) then
        t_diff = time(6) - time(5)
        print *, "Integration time: ", t_diff
    end if

    ! Calculate integral for x, y and z using Simpson's rule on each in turn
    ! test_Simpson_integral = integrate_triple_Simpson(test_coord_n, test_coord_lims, &
    !     mu_a, mu_b, m_a, m_b, mass_c, .false., .false., .true., comm)
    ! if (rank == 0) then
    !     print *, "Test integral using Simpson's rule = ", test_Simpson_integral
    ! end if

    ! Check the timer after test integration, at end of program
    call MPI_Barrier(comm, ierr)
    time(7) = MPI_Wtime()
    if (rank == 0) then
        t_diff = time(7) - time(6)
        print *, "Integration time: ", t_diff

        t_diff = time(7) - time(1)
        print *, "Program time: ", t_diff
    end if

    call MPI_Finalize(ierr)

contains

    subroutine calculate_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

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

    end subroutine calculate_masses


    ! Uses Simpson's rule to integrate x^2 or sin(x) between a and b
    ! n is the number of subintervals, which must be positive and even
    function integrate_single_Simpson(n, a, b, is_x_squared) result(integral)

        implicit none

        integer, intent(in) :: n
        real(wp), intent(in) :: a, b
        logical, intent(in) :: is_x_squared

        real(wp) :: width, x, integral
        real(wp) :: y(0:n)
        integer :: i

        if (mod(n, 2) /= 0 .or. n < 0) then
            print *, "n must be an even, positive integer"
            stop
        end if

        width = abs(b - a) / real(n, kind=wp)

        integral = 0._wp

        do i = 0, n
            x = a + real(i, kind=wp) * width
            y(i) = mixed_jacobi_integrand_func(x, is_x_squared)

            if (i == 0 .or. i == n) then
                integral = integral + y(i)
            else if (mod(i, 2) == 0) then
                integral = integral + 2._wp * y(i)
            else
                integral = integral + 4._wp * y(i)
            end if
        end do

        integral = width * integral / 3._wp

    end function integrate_single_Simpson


    ! Returns either x^4 or sin(x)cos^2(x) based on flag passed
    ! x = R_a, R_b or gamma_ab
    function mixed_jacobi_integrand_func(x, is_x_squared) result(integrand)

        implicit none

        real(wp), intent(in) :: x
        logical, intent(in) :: is_x_squared

        real(wp) :: integrand

        if (is_x_squared) then
            integrand = x**2._wp * x**2._wp
        else
            integrand = sin(x) * cos(x)**2._wp
        end if

    end function mixed_jacobi_integrand_func


    ! Returns R_a^2 R_b^2 cos^2(gamma_ab) R_a^2 r_a^2 sin(gamma_a)
    ! R_b and gamma_ab are rewritten in terms of R_a, r_a and gamma_a
    ! R_a = x, r_a = y and gamma_a = z
    function single_jacobi_integrand_func(x, y, z, mu_a, mass_c, m_b) result(integrand)

        implicit none

        real(wp), intent(in) :: x, y, z, mu_a, mass_c, m_b

        real(wp) :: integrand, integrand_volume, integrand_R_a, integrand_R_b, &
            integrand_gamma_ab

        ! Volume component to be integrated
        integrand_volume = x**2._wp * y**2._wp * sin(z)

        ! Function to be integrated:
        ! R_a^2
        integrand_R_a = x**2._wp

        ! R_b^2
        integrand_R_b = m_b**2._wp * ( (x / mass_c)**2._wp + (y / mu_a)**2._wp + &
            2._wp * x * y * cos(z) / (mu_a * mass_c) )

        ! cos(gamma_ab)^2
        ! Undefined for R_a or R_b = 0, but integrand = 0 anyway
        ! R_a^2 and R_b^2 must also be positive, so avoid rounding errors
        if (integrand_R_a <= 0._wp .or. integrand_R_b <= 0._wp) then
            integrand_gamma_ab = 0._wp
        else
            integrand_gamma_ab = (m_b / sqrt(integrand_R_b)) * &
            (-((y * cos(z) / mu_a) + x / mass_c ))

            integrand_gamma_ab = integrand_gamma_ab**2._wp
        end if

        integrand = integrand_volume * integrand_R_a * integrand_R_b * integrand_gamma_ab

    end function single_jacobi_integrand_func


    ! Calculate r_a from R_a, R_b and gamma_ab
    function calc_r_a(R_a, R_b, gamma_ab, mu_a, m_b, mass_c) result(coord)

        implicit none

        real(wp), intent(in) :: R_a, R_b, gamma_ab, mu_a, m_b, mass_c

        ! coord = r_a
        real(wp) :: coord

        coord = mu_a * sqrt( (R_a / mass_c)**2._wp  + (R_b / m_b)**2._wp &
            + 2._wp * (R_a * R_b * cos(gamma_ab) / (mass_c * m_b)) )

    end function calc_r_a


    ! Calculate gamma_a from R_a, r_a and gamma_ab
    ! R_b is calculated first, and -100._wp is returned if it is invalid
    function calc_gamma_a(small_r_a, mu_a, R_a, mass_c, gamma_ab, m_b, lims, rand_root, &
        pos_root) result(gamma_a)

        implicit none

        real(wp), intent(in) :: small_r_a, mu_a, R_a, mass_c, gamma_ab, m_b, lims(6)
        logical, intent(in) :: rand_root, pos_root

        logical :: positive_root

        real(wp) :: cos_numerator, gamma_a, temp_R_b_1, temp_R_b_2, R_b_over_m_b, rand_num, &
            gamma_tol, R_b_tol

            gamma_tol = 0.0000001_wp
            R_b_tol = 0.0000001_wp

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

        ! gamma_a is undefined if either R_a = 0 or r_a = 0
        if (small_r_a == 0._wp .or. R_a == 0._wp) then
            gamma_a = 0._wp
            return
        end if

        ! Calculate R_b / m_b
        temp_R_b_1 = (small_r_a / mu_a)**2._wp + (R_a / mass_c)**2._wp * &
            (cos(gamma_ab)**2._wp - 1._wp)

        ! Must be positive as will be square rooted
        if (temp_R_b_1 < 0._wp) then
            gamma_a = -100._wp
            return
        else
            temp_R_b_2 = - (R_a / mass_c) * cos(gamma_ab)
            if (positive_root) then
                R_b_over_m_b = temp_R_b_2 + sqrt(temp_R_b_1)
            else
                R_b_over_m_b = temp_R_b_2 - sqrt(temp_R_b_1)
            end if
        end if

        ! Check R_b within valid limits
        if ( (R_b_over_m_b < lims(3) / m_b) .or. (R_b_over_m_b > lims(4) / m_b) ) then

            ! Try opposite root before aborting
            if (positive_root) then
                R_b_over_m_b = temp_R_b_2 - sqrt(temp_R_b_1)
            else
                R_b_over_m_b = temp_R_b_2 + sqrt(temp_R_b_1)
            end if

            if ( (R_b_over_m_b < lims(3) / m_b) &
            .or. (R_b_over_m_b > lims(4) / m_b) ) then
                if (R_b_over_m_b > lims(4) / m_b &
                .and. R_b_over_m_b < R_b_tol + lims(4) / m_b) then
                    R_b_over_m_b = lims(4) / m_b
                else if (R_b_over_m_b < lims(3) / m_b &
                .and. R_b_over_m_b > - R_b_tol + lims(3) / m_b) then
                    R_b_over_m_b = lims(3) / m_b
                else
                    gamma_a = -100._wp
                    return
                end if
            end if
        end if

        cos_numerator = - mu_a * ( (R_a / mass_c) + R_b_over_m_b * cos(gamma_ab) )

        gamma_a = cos_numerator / small_r_a

        ! Prevent NaN when taking acos - there is probably a better solution!
        if (gamma_a > 1._wp) then
                if (gamma_a <= (1._wp+gamma_tol)) then
                    gamma_a = 1._wp
                else
                    print *, "Warning: invalid gamma_a: ", gamma_a
                    gamma_a = 1._wp
                    stop
                end if
            end if
            if (gamma_a < -1._wp) then
                if (gamma_a >= -(1._wp+gamma_tol)) then
                    gamma_a = -1._wp
                else
                    print *, "Warning: invalid gamma_a: ", gamma_a
                    gamma_a = -1._wp
                    stop
                end if
            end if

        gamma_a = acos(gamma_a)

    end function calc_gamma_a


    ! Use Simpson's rule to integrate three variables
    ! The mixed_jacobi flag uses constant limits, while the single_jacobi and
    ! test_coords flags use limits that depend on other variables
    ! These dependencies are defined in separate functions
    function integrate_triple_Simpson(n, mixed_lims, mu_a, mu_b, m_a, m_b, mass_c, &
        mixed_jacobi, single_jacobi, test_coords, sub_comm) result(total_integral)

        implicit none

        real(wp), intent(in) :: mixed_lims(6), mu_a, mu_b, m_a, m_b, mass_c
        integer, intent(in) :: n(3)
        logical, intent(in) :: mixed_jacobi, single_jacobi, test_coords

        ! MPI variables
        integer, intent(in) :: sub_comm

        real(wp) :: width(3), x, y, z, lims(6), temp_integral, z_integral, y_integral, &
            partial_integral, total_integral, prev_lims(2), r_lims(3, n(1)+1), &
            gamma_lims(4, n(1)+1, n(2)+1)
        integer :: i, j, k, partial_count, count, imin, imax, progress, num_r_values, &
            num_gamma_values
        logical :: verbose, save_lims

        ! MPI variables
        integer :: sub_comm_size, ierr, r_tag, gamma_tag, recv_status_1(MPI_STATUS_SIZE), &
                recv_status_2(MPI_STATUS_SIZE)

        verbose = .true.
        save_lims = .true.

        call MPI_Comm_size(sub_comm, sub_comm_size, ierr)

        if ((mixed_jacobi .and. single_jacobi) .or. (mixed_jacobi .and. test_coords) &
        .or. (single_jacobi .and. test_coords)) then
            print *, "Only one coordinate flag should be set to true"
            stop
        end if

        ! Width's of Simpson's rule subintervals for each coordinate
        ! Still use full mixed R_a limits for single Jacobi and test coords
        lims(1) = mixed_lims(1)
        lims(2) = mixed_lims(2)
        width(1) = abs(lims(2) - lims(1)) / real(n(1), kind=wp)
        if (mixed_jacobi) then
            lims(3) = mixed_lims(3)
            lims(4) = mixed_lims(4)
            lims(5) = mixed_lims(5)
            lims(6) = mixed_lims(6)
            width(2) = abs(lims(4) - lims(3)) / real(n(2), kind=wp)
            width(3) = abs(lims(6) - lims(5)) / real(n(3), kind=wp)
        end if

        partial_integral = 0._wp
        total_integral = 0._wp
        temp_integral = 0._wp

        partial_count = 0
        count = 0

        imin = rank * n(1) / sub_comm_size
        if (rank == sub_comm_size - 1) then
            imax = n(1)
        else
            imax = -1 + (rank + 1) * n(1) / sub_comm_size
        end if

        ! print *, "On rank ", rank, "Loop min = ", imin, "Loop max = ", imax
        ! Integrate each variable in turn, covering full limits
        ! Variables labelled x, y and z for simplicity
        do i = imin, imax

            if ((imax - imin >= 4) .and. verbose) then
                progress = 100 * (i - imin) / (imax - imin)
                if (mod((i - imin), (imax - imin + 1) / 5) == 0 &
                .or. i == imax) then
                    print *, "On rank", rank, "...Progress... ", progress, "%"
                end if
            end if

            ! Set value for R_a (mixed and single) or x (test)
            x = lims(1) + real(i, kind=wp) * width(1)

            ! Total R_b, r_a or y integral for set R_a, R_a or x
            y_integral = 0._wp

            ! For single Jacobi, calculate r_a limits from R_a
            if (single_jacobi) then
                lims(3:4) = get_limits(.true., .false., .false., &
                    x, 0._wp, mixed_lims, mu_a, m_b, mass_c)
                width(2) = abs(lims(4) - lims(3)) / real(n(2), kind=wp)

                if (save_lims) then
                    r_lims(1, i+1) = x
                    r_lims(2:3, i+1) = lims(3:4)
                end if
            end if

            if (test_coords) then
                lims(3:4) = get_test_limits(.true., .false., x, 0._wp)
                width(2) = abs(lims(4) - lims(3)) / real(n(2), kind=wp)
            end if

            do j = 0, n(2)

                if (width(2) == 0._wp) then
                    exit
                end if

                ! Set value for R_b (mixed), r_a (single) or y (test)
                y = lims(3) + real(j, kind=wp) * width(2)

                ! Total gamma_ab, gamma_a or z intergral for set
                ! (R_a, R_b), (R_a, r_a) or (x, y)
                z_integral = 0._wp

                ! For single Jacobi, calculate gamma_a limits from R_a and r_a
                if (single_jacobi) then
                    ! prev_lims = lims(5:6)
                    ! lims(5:6) = get_limits(.false., .true., .true., &
                    !     x, y, mixed_lims, mu_a, m_b, mass_c)
                    ! print *, "gamma_a min, max estimates: ", lims(5:6)
                    lims(5:6) = get_limits(.false., .true., .false., &
                        x, y, mixed_lims, mu_a, m_b, mass_c)
                    ! print *, "gamma_a min, max analytical: ", lims(5:6)

                    ! If unable to find limits, estimate as previous?
                    if (lims(5) == 0._wp .and. lims(6) == 0._wp .and. (i+j+k)>0) then
                        ! lims(5:6) = prev_lims
                        partial_count = partial_count + 1
                    end if

                    if (save_lims) then
                        gamma_lims(1, i+1, j+1) = x
                        gamma_lims(2, i+1, j+1) = y
                        gamma_lims(3:4, i+1, j+1) = lims(5:6)
                    end if

                    width(3) = abs(lims(6) - lims(5)) / real(n(3), kind=wp)
                end if

                if (test_coords) then
                    lims(5:6) = get_test_limits(.false., .true., x, y)
                    width(3) = abs(lims(6) - lims(5)) / real(n(3), kind=wp)
                end if

                do k = 0, n(3)

                    if (width(3) == 0._wp) then
                        exit
                    end if

                    ! Set value for gamma_ab (mixed), gamma_a (single)
                    ! or z (test)
                    z = lims(5) + real(k, kind=wp) * width(3)

                    ! Use Simpson's rule to add contributions
                    ! for this subinterval
                    if (mixed_jacobi) then
                        temp_integral = mixed_jacobi_integrand_func(z, &
                            .false.)
                    else if (single_jacobi) then
                        temp_integral = single_jacobi_integrand_func(x, &
                            y, z, mu_a, mass_c, m_b)
                    else if (test_coords) then
                        temp_integral = test_integrand_func(x, y, z)
                    end if

                    if (k == 0 .or. k == n(3)) then
                        z_integral = z_integral + temp_integral
                    else if (mod(k, 2) == 0) then
                        z_integral = z_integral + 2._wp * temp_integral
                    else
                        z_integral = z_integral + 4._wp * temp_integral
                    end if

                end do

                ! Total gamma_ab, gamma_a or z intergral for set
                ! (R_a, R_b), (R_a, r_a) or (x, y)
                z_integral = width(3) * z_integral / 3._wp

                ! Use Simpon's rule to add contributions for this subinterval
                ! Includes multiplication by total gamma_ab, gamma_a or z integral
                if (mixed_jacobi) then
                    temp_integral = z_integral * &
                        mixed_jacobi_integrand_func(y, .true.)
                else if (single_jacobi) then
                    temp_integral = z_integral
                else if (test_coords) then
                    temp_integral = z_integral
                end if

                if (j == 0 .or. j == n(2)) then
                    y_integral = y_integral + temp_integral
                else if (mod(j, 2) == 0) then
                    y_integral = y_integral + 2._wp * temp_integral
                else
                    y_integral = y_integral + 4._wp * temp_integral
                end if

            end do

            ! Total R_b, r_a or y integral for set R_a, R_a or x
            y_integral = width(2) * y_integral / 3._wp

            ! Use Simpon's rule to add contributions for this subinterval
            ! Includes multiplication by total R_b/r_a/y integral
            if (mixed_jacobi) then
                temp_integral = y_integral * mixed_jacobi_integrand_func(x, .true.)
            else if (single_jacobi) then
                temp_integral = y_integral
            else if (test_coords) then
                temp_integral = y_integral
            end if

            if (i == 0 .or. i == n(1)) then
                partial_integral = partial_integral + temp_integral
            else if (mod(i, 2) == 0) then
                partial_integral = partial_integral + 2._wp * temp_integral
            else
                partial_integral = partial_integral + 4._wp * temp_integral
            end if

        end do

        ! Sum integrals from all processes
        call MPI_Reduce(partial_integral, total_integral, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
            sub_comm, ierr)
        call MPI_Reduce(partial_count, count, 1, MPI_INT, MPI_SUM, 0, sub_comm, ierr)

        ! Send limits to rank 0 to save
        r_tag = 0
        gamma_tag = 1
        if (single_jacobi .and. save_lims) then
            if (rank == 0) then
                ! Receive from all other ranks
                do i = 1, sub_comm_size - 1
                    imin = i * n(1) / sub_comm_size
                    if (i == sub_comm_size - 1) then
                        imax = n(1)
                    else
                        imax = -1 + (i + 1) * n(1) / sub_comm_size
                    end if

                    num_r_values = (imax - imin + 1) * 3
                    call MPI_Recv(r_lims(:, imin+1:imax+1), num_r_values, &
                        MPI_DOUBLE_PRECISION, i, r_tag, sub_comm, recv_status_1, ierr)

                    num_gamma_values = (imax - imin + 1) * 4 * (n(2) + 1)
                    call MPI_Recv(gamma_lims(:, imin+1:imax+1, :), &
                        num_gamma_values, MPI_DOUBLE_PRECISION, i, gamma_tag, &
                        sub_comm, recv_status_2, ierr)
                end do
            else
                ! All non-zero ranks send to rank 0
                num_r_values = (imax - imin + 1) * 3
                call MPI_Ssend(r_lims(:, imin+1:imax+1), num_r_values, MPI_DOUBLE_PRECISION, &
                    0, r_tag, sub_comm, ierr)

                num_gamma_values = (imax - imin + 1) * 4 * (n(2) + 1)
                call MPI_Ssend(gamma_lims(:, imin+1:imax+1, :), num_gamma_values, &
                    MPI_DOUBLE_PRECISION, 0, gamma_tag, sub_comm, ierr)
            end if
        end if

        ! All sends/receives complete
        call MPI_Barrier(sub_comm, ierr)

        if (rank == 0) then
            ! Total (R_a, R_a or x) integral
            total_integral = width(1) * total_integral / 3._wp

            ! For single Jacobi coordinates, multiply by prefactor
            if (single_jacobi) then
                total_integral = ( (mu_a * mu_b) / (m_a * m_b) )**(-3._wp/2._wp) * &
                    total_integral
            end if

            if (verbose) then
                print *, "Number of limits not found: ", count
            end if

            if (single_jacobi .and. save_lims) then
                open (unit=42, file='outputs/r_lims', form='unformatted')
                write(42) r_lims
                close (unit=42)
                open (unit=43, file='outputs/gamma_lims', form='unformatted')
                write(43) gamma_lims
                close (unit=43)
            end if
        end if

    end function integrate_triple_Simpson


    ! Calculate integral with Monte Carlo
    function integrate_MC(n, mixed_lims, mu_a, mu_b, m_a, m_b, mass_c, sub_comm) &
        result(total_integral)

        implicit none

        real(wp), intent(in) :: mixed_lims(6), mu_a, mu_b, m_a, m_b, mass_c
        integer, intent(in) :: n

        ! MPI variables
        integer, intent(in) :: sub_comm

        real(wp) :: width(3), coords(3), single_lims(6), temp_integral, partial_integral, &
            total_integral, r_lims(3, n+1), gamma_lims(4, n+1)
        integer :: i, imin, imax, progress, partial_count, count, num_r_values, &
            num_gamma_values
        logical :: verbose, save_lims

        ! MPI variables
        integer :: sub_comm_size, ierr, r_tag, gamma_tag, recv_status_1(MPI_STATUS_SIZE), &
                recv_status_2(MPI_STATUS_SIZE)

        verbose = .true.
        save_lims = .false.

        call MPI_Comm_size(sub_comm, sub_comm_size, ierr)

        partial_integral = 0._wp
        total_integral = 0._wp

        partial_count = 0
        count = 0

        width = (/ 0._wp, 0._wp, 0._wp /)
        coords = (/ 0._wp, 0._wp, 0._wp /)
        single_lims = (/ 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp /)

        ! R_a_min and R_a_max can use full mixed R_a limits
        single_lims(1) = mixed_lims(1)
        single_lims(2) = mixed_lims(2)

        width(1) = abs(single_lims(2) - single_lims(1))

        imin = rank * n / sub_comm_size
        if (rank == sub_comm_size - 1) then
            imax = n
        else
            imax = -1 + (rank + 1) * n / sub_comm_size
        end if

        ! print *, "On rank ", rank, "Loop min = ", imin, "Loop max = ", imax
        do i = imin, imax

            if ((imax - imin >= 4) .and. verbose) then
                progress = 100 * (i - imin) / (imax - imin)
                if (mod((i - imin), (imax - imin + 1) / 5) == 0 &
                .or. i == imax) then
                    print *, "On rank", rank, "...Progress... ", progress, "%"
                end if
            end if

            temp_integral = 0._wp

            ! Generate random points in range 0 to 1
            call random_number(coords)

            ! Random R_a
            coords(1) = single_lims(1) + coords(1) * width(1)

            ! r_a limits from current R_a
            single_lims(3:4) = get_limits(.true., .false., .false., &
                coords(1), 0._wp, mixed_lims, mu_a, m_b, mass_c)

            if (save_lims) then
                r_lims(1, i+1) = coords(1)
                r_lims(2:3, i+1) = single_lims(3:4)
            end if

            width(2) = abs(single_lims(4) - single_lims(3))

            ! Random r_a in calculated range
            coords(2) = single_lims(3) + coords(2) * width(2)

            ! gamma_a limits from current R_a and r_a
            single_lims(5:6) = get_limits(.false., .true., .true., &
                coords(1), coords(2), mixed_lims, mu_a, m_b, mass_c)

            if (save_lims) then
                gamma_lims(1, i+1) = coords(1)
                gamma_lims(2, i+1) = coords(2)
                gamma_lims(3:4, i+1) = single_lims(5:6)
            end if

            width(3) = abs(single_lims(6) - single_lims(5))

            ! If unable to find limits, estimate as previous?
            if (single_lims(5) == 0._wp .and. single_lims(6) == 0._wp .and. i > 0) then
                ! lims(5:6) = prev_lims
                partial_count = partial_count + 1
            end if

            ! Random gamma_a in calculated range
            coords(3) = single_lims(5) + coords(3) * width(3)

            ! Evaluate integral element
            temp_integral = single_jacobi_integrand_func(coords(1), coords(2), &
                coords(3), mu_a, mass_c, m_b)
            temp_integral = width(1) * width(2) * width(3) * temp_integral
            partial_integral = partial_integral + temp_integral

        end do

        ! Sum integrals from all processes
        call MPI_Reduce(partial_integral, total_integral, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
            sub_comm, ierr)
        call MPI_Reduce(partial_count, count, 1, MPI_INT, MPI_SUM, 0, sub_comm, ierr)

        ! Send limits to rank 0 to save
        r_tag = 0
        gamma_tag = 1
        if (save_lims) then
            if (rank == 0) then
                ! Receive from all other ranks
                do i = 1, sub_comm_size - 1
                    imin = i * n / sub_comm_size
                    if (i == sub_comm_size - 1) then
                        imax = n
                    else
                        imax = -1 + (i + 1) * n / sub_comm_size
                    end if

                    num_r_values = (imax - imin + 1) * 3
                    call MPI_Recv(r_lims(:, imin+1:imax+1), num_r_values, &
                        MPI_DOUBLE_PRECISION, i, r_tag, sub_comm, recv_status_1, ierr)

                    num_gamma_values = (imax - imin + 1) * 4
                    call MPI_Recv(gamma_lims(:, imin+1:imax+1), &
                        num_gamma_values, MPI_DOUBLE_PRECISION, i, gamma_tag, &
                        sub_comm, recv_status_2, ierr)
                end do
            else
                ! All non-zero ranks send to rank 0
                num_r_values = (imax - imin + 1) * 3
                call MPI_Ssend(r_lims(:, imin+1:imax+1), num_r_values, MPI_DOUBLE_PRECISION, &
                    0, r_tag, sub_comm, ierr)

                num_gamma_values = (imax - imin + 1) * 4
                call MPI_Ssend(gamma_lims(:, imin+1:imax+1), num_gamma_values, &
                    MPI_DOUBLE_PRECISION, 0, gamma_tag, sub_comm, ierr)
            end if
        end if

        ! All sends/receives complete
        call MPI_Barrier(sub_comm, ierr)

        if (rank == 0) then

            total_integral = total_integral / real(n, kind=wp)
            total_integral = ( (mu_a * mu_b) / (m_a * m_b) )**(-3._wp/2._wp) * total_integral

            if (verbose) then
                print *, "Number of limits not found: ", count
            end if

            if (save_lims) then
                open (unit=42, file='outputs/r_lims', form='unformatted')
                write(42) r_lims
                close (unit=42)
                open (unit=43, file='outputs/gamma_lims', form='unformatted')
                write(43) gamma_lims
                close (unit=43)
            end if
        end if

    end function integrate_MC


    ! Estimate the minimum and maximum values of r_a, based on known R_a
    ! or gamma_a, based on known R_a and r_a
    function estimate_jacobi_lims(R_a, mixed_lims, mu_a, m_b, mass_c, estimating_r_a, &
        estimating_gamma_a, small_r_a) result(single_lims)

        implicit none

        real(wp), intent(in) :: R_a, mixed_lims(6), mu_a, m_b, mass_c, small_r_a
        logical, intent(in) :: estimating_r_a, estimating_gamma_a

        real(wp) :: test_coord(100000), mixed_coords(3), width(3), single_lims(2), test_r_a, &
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

        ! Generate random R_b and gamma_ab (R_a is fixed)
        call random_number(mixed_coords(2:3))

        ! Ranges of R_a, R_b and gamma_ab
        width(1) = abs(mixed_lims(2) - mixed_lims(1))
        width(2) = abs(mixed_lims(4) - mixed_lims(3))
        width(3) = abs(mixed_lims(6) - mixed_lims(5))

        mixed_coords(1) = R_a

        ! Track valid coordinates generated
        count = 0

        if (random_range) then
            do i = 1, n - 4

                ! Generate random R_b and gamma_ab that change each iteration
                call random_number(mixed_coords(2:3))

                mixed_coords(2) = mixed_lims(3) + mixed_coords(2) * width(2)
                mixed_coords(3) = mixed_lims(5) + mixed_coords(3) * width(3)

                if (estimating_r_a) then
                    count = count + 1
                    test_coord(count) = calc_r_a(mixed_coords(1), &
                    mixed_coords(2), mixed_coords(3), mu_a, m_b, mass_c)
                else if (estimating_gamma_a) then

                    ! Attempt to calculate gamma from R_a, r_a and gamma_ab
                    temp_gamma = calc_gamma_a(small_r_a, mu_a, &
                        mixed_coords(1), mass_c, mixed_coords(3), m_b, &
                        mixed_lims, .true., .true.)

                    ! If R_b is invalid, ignore gamma calculated
                    if (temp_gamma /= -100._wp) then
                        count = count + 1
                        test_coord(count) = temp_gamma
                    end if
                end if
            end do

            do i = n - 3, n
                if (i == (n-3) .or. i == (n-2)) then
                    mixed_coords(2) = mixed_lims(3)
                else
                    mixed_coords(2) = mixed_lims(3) + width(2)
                end if
                if (i == (n-3) .or. i == (n-1)) then
                    mixed_coords(3) = mixed_lims(5)
                else
                    mixed_coords(3) = mixed_lims(5) + width(3)
                end if

                if (estimating_r_a) then
                    count = count + 1
                    test_coord(count) = calc_r_a(mixed_coords(1), &
                        mixed_coords(2), mixed_coords(3), mu_a, m_b, &
                        mass_c)
                else if (estimating_gamma_a) then

                    ! Attempt to calculate gamma from R_a, r_a and gamma_ab
                    temp_gamma = calc_gamma_a(small_r_a, mu_a, &
                        mixed_coords(1), mass_c, mixed_coords(3), m_b, &
                        mixed_lims, .true., .true.)

                    ! If R_b is invalid, ignore gamma calculated
                    if (temp_gamma /= -100._wp) then
                        count = count + 1
                        test_coord(count) = temp_gamma
                    end if
                end if
            end do
        else
            do i = 1, n

                ! Evenly sample entire range of each coord
                mixed_coords(2) = mixed_lims(3) + (real(i-1, kind=wp) * width(2) / real(n, kind=wp))

                do j = 1, n
                    ! Evenly sample entire range of each coord
                    mixed_coords(3) = mixed_lims(5) + (real(j-1, kind=wp) * &
                        width(3) / real(n, kind=wp))

                    if (estimating_r_a) then
                        test_coord(count) = calc_r_a(mixed_coords(1), &
                            mixed_coords(2), mixed_coords(3), mu_a, &
                            m_b, mass_c)
                        count = count + 1
                    else if (estimating_gamma_a) then
                        ! Attempt to calculate gamma from
                        ! R_a, r_a and gamma_ab
                        temp_gamma = calc_gamma_a(small_r_a, mu_a, &
                            mixed_coords(1), mass_c, mixed_coords(3), &
                            m_b, mixed_lims, .true., .true.)

                        ! If R_b is invalid, ignore gamma calculated
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


    ! Get the limits of r_a (get_r_lims) for a known R_a
    ! or the limits of gamma_a (get_gamma_lims) for a known R_a and r_a
    ! The limits can be estimated (estimate_lims) or calculated analytically
    ! Both limits are set to 0._wp if valid coordinates cannot be found
    function get_limits(get_r_lims, get_gamma_lims, estimate_lims, R_a, small_r_a, &
        mixed_lims, mu_a, m_b, mass_c) result(lims)

        implicit none

        real(wp), intent(in) :: R_a, small_r_a, mixed_lims(6), mu_a, m_b, mass_c
        logical, intent(in) :: get_r_lims, get_gamma_lims, estimate_lims

        real(wp) :: lims(2), gamma_tol, temp_gamma

        if (get_r_lims .and. get_gamma_lims) then
            print *, "Only one coordinate flag should be set to true"
            stop
        end if

        gamma_tol = 0.0000001_wp

        if (get_r_lims) then
            if (estimate_lims) then
                lims(1:2) = estimate_jacobi_lims(R_a, mixed_lims, mu_a, m_b, &
                    mass_c, .true., .false., 0._wp)
            else
                ! r_a_min uses R_a_min, R_b_min, gamma_ab_max
                ! (r smaller with larger gamma
                ! although no effect for R_a = R_b = 0)
                lims(1) = calc_r_a(mixed_lims(1), mixed_lims(3), mixed_lims(6), &
                    mu_a, m_b, mass_c)

                ! r_a_max depends on max current (single) R_a
                ! Also uses full (mixed) R_b_max and gamma_ab_min
                lims(2) = calc_r_a(R_a, mixed_lims(4), mixed_lims(5), mu_a, m_b, &
                    mass_c)
            end if

        else if (get_gamma_lims) then
            if (estimate_lims) then
                lims(1:2) = estimate_jacobi_lims(R_a, mixed_lims, mu_a, m_b, &
                    mass_c, .false., .true., small_r_a)
            else

                ! Note: This option may not work for all R_a and R_b limits
                ! gamma_ab must run between 0 and pi
                ! Use with caution

                ! gamma_a is undefined for R_a = 0 or r_a = 0
                if (R_a == 0._wp .or. small_r_a == 0._wp) then
                    lims(1) = 0._wp
                    lims(2) = 0._wp
                else
                    ! Calculate cos(gamma_a) when R_b is maximum
                    ! to get the minimum gamma_a
                    lims(1) = mass_c * mu_a * ((mixed_lims(4) / m_b )**2._wp &
                        - (R_a / mass_c)**2._wp - (small_r_a / mu_a)**2._wp) &
                        / (2._wp * R_a * small_r_a)

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
                    temp_gamma = calc_gamma_a(small_r_a, mu_a, R_a, &
                        mass_c, mixed_lims(5), m_b, mixed_lims, .false., &
                        .true.)
                    if (temp_gamma /= -100._wp) then
                        lims(2) = temp_gamma
                    end if
                    temp_gamma = calc_gamma_a(small_r_a, mu_a, R_a, &
                        mass_c, mixed_lims(5), m_b, mixed_lims, .false., &
                        .false.)
                    if (temp_gamma /= -100._wp .and. temp_gamma > lims(2)) then
                        lims(2) = temp_gamma
                    end if

                    ! Test gamma_ab maximum for +- roots
                    temp_gamma = calc_gamma_a(small_r_a, mu_a, R_a, &
                        mass_c, mixed_lims(6), m_b, mixed_lims, .false., &
                        .true.)
                    if (temp_gamma /= -100._wp .and. temp_gamma > lims(2)) then
                        lims(2) = temp_gamma
                    end if
                    temp_gamma = calc_gamma_a(small_r_a, mu_a, R_a, &
                        mass_c, mixed_lims(6), m_b, mixed_lims, .false., &
                        .false.)
                    if (temp_gamma /= -100._wp .and. temp_gamma > lims(2)) then
                        lims(2) = temp_gamma
                    end if
                end if
            end if
        end if

    end function get_limits


    ! Get limits for the test integrand
    ! y limits are 0 and (1-x)
    ! z limits are 0 and (1-x-y)
    function get_test_limits(get_y_lims, get_z_lims, x, y) result(lims)

        implicit none

        real(wp), intent(in) :: x, y
        logical, intent(in) :: get_y_lims, get_z_lims

        real(wp) :: lims(2)

        if (get_y_lims .and. get_z_lims) then
            print *, "Only one coordinate flag should be set to true"
            stop
        end if

        if (get_y_lims) then
            ! Depends on x
            lims(1:2) = (/ 0._wp, 1._wp - x/)
        else if (get_z_lims) then
            ! Depends on x and y
            lims(1:2) = (/ 0._wp, 1._wp - x - y /)
        end if

    end function get_test_limits

    ! Returns (x+y+z)^-0.5
    ! For 1/denom, requires adjustment for 1/0, and so a much larger n is needed
    ! Also tested with (x+y+z)^-0
    ! Expected results:
    !       (x+y+z)^-0.5: 0.2
    !       (x+y+z)^0.5: 0.142857142857143
    function test_integrand_func(x, y, z) result(integrand)

        implicit none

        real(wp), intent(in) :: x, y, z
        real(wp) :: integrand, denominator

        denominator = sqrt(x + y + z)
        if (denominator  < 0.00000001) then
            denominator =  0.00000001
        end if
        integrand = 1._wp / denominator

    end function test_integrand_func

end program jacobi