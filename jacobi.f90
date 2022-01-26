program jacobi

        implicit none
        integer :: mixed_jacobi_n(3), single_jacobi_n(3), mc_n, test_coord_n(3)
        real :: mixed_jacobi_lims(6), R_a_integral, R_b_integral, theta_ab_integral, &
                mixed_integral, mixed_Simpson_integral, single_Simpson_integral, &
                single_MC_integral, mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b, &
                test_Simpson_integral, test_coord_lims(6)

        ! Number of points to generate for MC integration
        mc_n = 10000

        ! Number of Simpson cells to split mixed Jacobi ingegral into (R_a_n, R_b_n, theta_ab_n)
        mixed_jacobi_n = (/ 50, 50, 50 /)

        ! Number of Simpson cells to split single Jacobi integral into (R_a_n, r_a_n, theta_a_n)
        single_jacobi_n = (/ 50, 50, 50 /)

        ! Number of Simpson cells to split test coordinate integral into (x_n, y_n, x_n)
        test_coord_n = (/ 1000, 1000, 1000 /)

        ! Limits of integration (R_a_min, R_a_max, R_b_min, R_b_max, theta_ab_min, thata_ab_max)
        mixed_jacobi_lims = (/ 0., 5., 0., 7., 0., 4.*atan(1.) /)

        ! Limits of integration (x_min, x_max, y_min, y_max, z_min, z_max)
        test_coord_lims = (/ 0., 1., 0., 1., 0., 1. /)

        ! Calculate separable integrals for R_a, R_b and theta_ab using Simpson's rule
        R_a_integral = integrate_single_Simpson(mixed_jacobi_n(1), mixed_jacobi_lims(1), &
                mixed_jacobi_lims(2), .true.)
        R_b_integral = integrate_single_Simpson(mixed_jacobi_n(2), mixed_jacobi_lims(3), &
                mixed_jacobi_lims(4), .true.)
        theta_ab_integral = integrate_single_Simpson(mixed_jacobi_n(3), mixed_jacobi_lims(5), &
                mixed_jacobi_lims(6), .false.)

        ! Integrals are separable so multiply for total integral
        mixed_integral = R_a_integral * R_b_integral * theta_ab_integral
        print *, "Separated mixed Jacobi integral using Simpson's rule = ", mixed_integral

        ! Calculate integral for R_a, R_b and theta_ab using Simpson's rule on each in turn
        mixed_Simpson_integral = integrate_triple_Simpson(mixed_jacobi_n, mixed_jacobi_lims, &
                mu_a, mu_b, m_a, m_b, mass_c, .true., .false., .false.)
        print *, "Unseparated mixed Jacobi integral using Simpson's rule = ", mixed_Simpson_integral

        ! Calculate mass relations (three masses defined within)
        call calculate_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

        ! Calculate integral for R_a, r_a and theta_a using Simpson's rule on each in turn
        single_Simpson_integral = integrate_triple_Simpson(single_jacobi_n, mixed_jacobi_lims, &
                mu_a, mu_b, m_a, m_b, mass_c, .false., .true., .false.)
        print *, "Single Jacobi integral using Simpson's rule = ", single_Simpson_integral

        ! Calculate integral for R_a, r_a and theta_a using Monte Carlo integration
        single_MC_integral = integrate_MC(mc_n, mixed_jacobi_lims, mu_a, mu_b, m_a, m_b, mass_c)
        print *, "Single Jacobi integral using Monte Carlo = ", single_MC_integral

        ! Calculate integral for x, y and z using Simpson's rule on each in turn
        test_Simpson_integral = integrate_triple_Simpson(test_coord_n, test_coord_lims, &
                mu_a, mu_b, m_a, m_b, mass_c, .false., .false., .true.)
        print *, "Test integral using Simpson's rule = ", test_Simpson_integral

contains

        subroutine calculate_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

                implicit none

                real, intent(inout) :: mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b

                ! Atom masses
                mass_a = 2.
                mass_b = 3.
                mass_c = 4.

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

        end subroutine calculate_masses


        ! Uses Simpson's rule to integrate x^2 or sin(x) between a and b
        ! n is the number of subintervals, which must be positive and even
        function integrate_single_Simpson(n, a, b, is_x_squared) result(integral)

                implicit none

                integer, intent(in) :: n
                real, intent(in) :: a, b
                logical, intent(in) :: is_x_squared

                real :: width, x, integral
                real :: y(0:n)
                integer :: i

                if (mod(n, 2) /= 0 .or. n < 0) then
                        print *, "n must be an even, positive integer"
                        stop
                end if

                width = abs(b - a) / real(n)

                integral = 0.

                do i = 0, n
                        x = a + real(i) * width
                        y(i) = jacobi_integrand_func(x, is_x_squared)

                        if (i == 0 .or. i == n) then
                                integral = integral + y(i)
                        else if (mod(i, 2) == 0) then
                                integral = integral + 2. * y(i)
                        else
                                integral = integral + 4. * y(i)
                        end if
                end do

                integral = width * integral / 3.

                if (is_x_squared) then
                        !print *, "The integral of x^2 from ", a, " to ", b, " is ", integral
                else
                        !print *, "The integral of sin(x) from ", a, " to ", b, " is ", integral
                end if

        end function integrate_single_Simpson


        ! Returns either x^2 or sin(x) based on flag passed
        function jacobi_integrand_func(x, is_x_squared) result(integrand)

                implicit none

                real, intent(in) :: x
                logical, intent(in) :: is_x_squared

                real :: integrand

                if (is_x_squared) then
                        integrand = x ** 2
                else
                        integrand = sin(x)
                end if

        end function jacobi_integrand_func


        function calculate_r_a(mu_a, R_a, M_c, R_b, m_b, gamma_ab) result(coord)

                implicit none

                real, intent(in) :: mu_a, R_a, M_c, R_b, m_b, gamma_ab

                ! coord = r_a
                real :: coord

                coord = mu_a * sqrt( (R_a / M_c)**2  + (R_b / m_b)**2 &
                        + 2. * (R_a * R_b * cos(gamma_ab) / (M_c * m_b)) )

        end function calculate_r_a


        function calculate_theta_a(small_r_a, mu_a, R_a, M_c, R_b, m_b, gamma_ab) result(coord)

                implicit none

                real, intent(in) :: small_r_a, mu_a, R_a, M_c, R_b, m_b, gamma_ab

                ! coord = theta_a
                real :: temp_r_a, cos_numerator, cos_denominator, cos_coord, coord

                temp_r_a = small_r_a
                temp_r_a = calculate_r_a(mu_a, R_a, M_c, R_b, m_b, gamma_ab)

                ! Temporary to prevent division by 0 - there is probably a better solution!
                if (temp_r_a < 0.00000001) then
                        cos_denominator = 0.00000001
                else
                        cos_denominator = temp_r_a
                end if

                cos_numerator = - mu_a * ( (R_a / M_c) + (R_b / m_b) * cos(gamma_ab) )

                cos_coord = cos_numerator / cos_denominator
                ! print *, "cos_theta = ", cos_coord

                ! Temporary to prevent NaN when taking acos - there is probably a better solution!
                if (cos_coord > 1.) then
                        cos_coord = 1.
                end if
                if (cos_coord < -1.) then
                        cos_coord = -1.
                end if

                coord = acos(cos_coord)
                !coord = cos_coord

                ! print *, "r_a, numerator = ", cos_denominator, cos_denominator
                ! print *, "theta = ", coord

        end function calculate_theta_a


        ! Use Simpson's rule to integrate three variables
        ! The mixed_jacobi flag uses constant limits, while the single_jacobi and
        ! test_coords flags use constant limits, defined in other functions
        ! Otherwise the limits are dependent on previous variables
        function integrate_triple_Simpson(n, mixed_lims, mu_a, mu_b, m_a, m_b, mass_c, &
                mixed_jacobi, single_jacobi, test_coords) result(total_integral)

                implicit none

                real, intent(in) :: mixed_lims(6), mu_a, mu_b, m_a, m_b, mass_c
                integer, intent(in) :: n(3)
                logical, intent(in) :: mixed_jacobi, single_jacobi, test_coords

                real :: width(3), x, y, z, lims(6), &
                        temp_integral, z_integral, y_integral, total_integral
                integer :: i, j, k

                if ((mixed_jacobi .and. single_jacobi) .or. (mixed_jacobi .and. test_coords) &
                        .or. (single_jacobi .and. test_coords)) then
                        print *, "Only one coordinate flag should be set to true"
                        stop
                end if

                ! Width's of Simpson's rule subintervals for each coordinate
                ! Still use full mixed R_a limits for single Jacobi
                lims(1) = mixed_lims(1)
                lims(2) = mixed_lims(2)
                width(1) = abs(lims(2) - lims(1)) / real(n(1))
                if (mixed_jacobi) then
                        lims(3) = mixed_lims(3)
                        lims(4) = mixed_lims(4)
                        lims(5) = mixed_lims(5)
                        lims(6) = mixed_lims(6)
                        width(2) = abs(lims(4) - lims(3)) / real(n(2))
                        width(3) = abs(lims(6) - lims(5)) / real(n(3))
                end if

                total_integral = 0.
                temp_integral = 0.

                ! Integrate each variable in turn, covering full limits
                ! Variables labelled x, y and z for simplicity
                do i = 0, n(1)

                        ! Set value for R_a (mixed and single)
                        x = lims(1) + real(i) * width(1)

                        ! Total R_b interegral for set R_a
                        y_integral = 0.

                        ! For single Jacobi, calculate r_a limits from R_a
                        if (single_jacobi) then
                                lims(3:4) = get_limits(.true., .false., .true., &
                                        x, 0., mixed_lims, mu_a, m_b, mass_c)
                        end if

                        if (test_coords) then
                                lims(3:4) = get_test_limits(.true., .false., x, 0.)
                        end if

                        width(2) = abs(lims(4) - lims(3)) / real(n(2))

                        do j = 0, n(2)

                                ! Set value for R_b (mixed) or r_a (single)
                                y = lims(3) + real(j) * width(2)

                                ! Total theta_ab interegral for set R_a, R_b
                                z_integral = 0.

                                ! For single Jacobi, calculate theta_a limits from R_a and r_a
                                if (single_jacobi) then
                                        lims(5:6) = get_limits(.false., .true., .true., &
                                                x, y, mixed_lims, mu_a, m_b, mass_c)

                                end if

                                if (test_coords) then
                                        lims(5:6) = get_test_limits(.false., .true., x, y)
                                end if

                                width(3) = abs(lims(6) - lims(5)) / real(n(3))

                                do k = 0, n(3)

                                        ! Set value for theta_ab (mixed) or theta_a (single)
                                        z = lims(5) + real(k) * width(3)

                                        ! Use Simpson's rule to add contributions for this subinterval
                                        if (single_jacobi .or. mixed_jacobi) then
                                                temp_integral = jacobi_integrand_func(z, .false.)
                                        else if (test_coords) then
                                                temp_integral = test_integrand_func(x, y, z)
                                        end if

                                        if (k == 0 .or. k == n(3)) then
                                                z_integral = z_integral + temp_integral
                                        else if (mod(k, 2) == 0) then
                                                z_integral = z_integral + 2. * temp_integral
                                        else
                                                z_integral = z_integral + 4. * temp_integral
                                        end if

                                end do

                                ! Total theta_ab intergral for set R_a, R_b
                                z_integral = width(3) * z_integral / 3.

                                ! Use Simpon's rule to add contributions for this subinterval
                                ! Inlcudes multiplication by total theta_ab integral
                                if (single_jacobi .or. mixed_jacobi) then
                                        temp_integral = z_integral * jacobi_integrand_func(y, .true.)
                                else if (test_coords) then
                                        temp_integral = z_integral
                                end if

                                if (j == 0 .or. j == n(2)) then
                                        y_integral = y_integral + temp_integral
                                else if (mod(j, 2) == 0) then
                                        y_integral = y_integral + 2. * temp_integral
                                else
                                        y_integral = y_integral + 4. * temp_integral
                                end if

                        end do

                        ! Total R_b integral for set R_a
                        y_integral = width(2) * y_integral / 3.

                        ! Use Simpon's rule to add contributions for this subinterval
                        ! Includes multiplication by total R_b integral
                        if (single_jacobi .or. mixed_jacobi) then
                                temp_integral = y_integral * jacobi_integrand_func(x, .true.)
                        else if (test_coords) then
                                temp_integral = y_integral
                        end if

                        if (i == 0 .or. i == n(1)) then
                                total_integral = total_integral + temp_integral
                        else if (mod(i, 2) == 0) then
                                total_integral = total_integral + 2. * temp_integral
                        else
                                total_integral = total_integral + 4. * temp_integral
                        end if

                end do

                ! Total (R_a) integral
                total_integral = width(1) * total_integral / 3.

                if (single_jacobi) then
                        total_integral = ( (mu_a * mu_b) / (m_a * m_b) )**(-3/2) * total_integral
                end if

        end function integrate_triple_Simpson


        ! Calculate integral with Monte Carlo
        function integrate_MC(n, mixed_lims, mu_a, mu_b, m_a, m_b, mass_c) result(total_integral)

                implicit none

                real, intent(in) :: mixed_lims(6), mu_a, mu_b, m_a, m_b, mass_c
                integer, intent(in) :: n

                real :: width(3), coords(3), single_lims(6), &
                temp_integral, total_integral
                integer :: i

                total_integral = 0.

                width = (/ 0, 0, 0 /)
                coords = (/ 0, 0, 0 /)
                single_lims = (/ 0, 0, 0, 0, 0, 0 /)

                ! R_a_min and R_a_max can use full mixed R_a limits
                single_lims(1) = mixed_lims(1)
                single_lims(2) = mixed_lims(2)

                width(1) = abs(single_lims(2) - single_lims(1))

                do i = 0, n

                        temp_integral = 0.

                        ! Generate random points in range 0 to 1
                        call random_number(coords)
                        ! print *, i, ". Random numbers: ", coords(1), coords(2), coords(3)

                        ! Random R_a
                        coords(1) = single_lims(1) + coords(1) * width(1)

                        single_lims(3:4) = get_limits(.true., .false., .true., &
                                coords(1), 0., mixed_lims, mu_a, m_b, mass_c)

                        width(2) = abs(single_lims(4) - single_lims(3))

                        ! Random r_a
                        coords(2) = single_lims(3) + coords(2) * width(2)

                        single_lims(5:6) = get_limits(.false., .true., .true., &
                                coords(1), coords(2), mixed_lims, mu_a, m_b, mass_c)

                        width(3) = abs(single_lims(6) - single_lims(5))

                        ! Random theta_a
                        coords(3) = single_lims(5) + coords(3) * width(3)

                        temp_integral = evaluate_function(coords(1), coords(2), coords(3))

                        if (total_integral > 1000000000) then
                                print *, "Function evaluated: ", temp_integral
                        end if
                        temp_integral = width(1) * width(2) * width(3) * temp_integral

                        if (total_integral > 1000000000) then
                                print *, "Coords: ", coords
                                print *, "Widths: ", width
                                print *, "Temp int: ", temp_integral
                        end if

                        total_integral = total_integral + temp_integral

                end do

                total_integral = total_integral / real(n)

                total_integral = ( (mu_a * mu_b) / (m_a * m_b) )**(-3/2) * total_integral

        end function integrate_MC


        ! Evaluates hardcoded function based on x, y and z
        function evaluate_function(x, y, z) result(total_func)

                implicit none

                real, intent(in) :: x, y, z

                real :: total_func, F

                total_func = 0.

                ! Function being integrated F(x,y,z)
                F = 1
                ! Multiply F by integration variables in integral 
                total_func = x**2 * y**2 * sin(z) * F
                ! total_func = sqrt(x+y+z)

                ! print *, "R_a, r_a, gamma_a: ", x, y, z
                ! print *, "Total func: ", total_func

        end function evaluate_function


        ! Estimate the minimum and maximum values of r_a, based on known R_a
        function estimate_jacobi_lims(R_a, mixed_lims, mu_a, m_b, mass_c, estimating_r_a, &
                estimating_theta_a, small_r_a) result(single_lims)

                implicit none

                real, intent(in) :: R_a, mixed_lims(6), mu_a, m_b, mass_c, small_r_a
                logical, intent(in) :: estimating_r_a, estimating_theta_a

                real :: test_coord(10000), mixed_coords(3), width(3), single_lims(2)
                integer :: i, j, n
                logical :: random_range

                ! Use randomly generated coordinates, rather than even distribution across range
                random_range = .true.

                if (random_range) then
                        n = size(test_coord)
                else
                        n = sqrt(real(size(test_coord)))
                end if

                ! Generate random R_b and theta_ab (R_a is fixed)
                call random_number(mixed_coords(2:3))

                ! Ranges of R_a, R_b and theta_ab
                width(1) = abs(mixed_lims(2) - mixed_lims(1))
                width(2) = abs(mixed_lims(4) - mixed_lims(3))
                width(3) = abs(mixed_lims(6) - mixed_lims(5))

                mixed_coords(1) = R_a

                if (random_range) then
                        do i = 1, n

                                ! Generate random R_b and theta_ab that change each iteration
                                call random_number(mixed_coords(2:3))

                                mixed_coords(2) = mixed_lims(3) + mixed_coords(2) * width(2)
                                mixed_coords(3) = mixed_lims(5) + mixed_coords(3) * width(3)

                                if (estimating_r_a) then
                                        test_coord(i) = calculate_r_a(mu_a, mixed_coords(1), mass_c, &
                                                mixed_coords(2), m_b, mixed_coords(3))
                                else if (estimating_theta_a) then
                                        test_coord(i) = calculate_theta_a(small_r_a, mu_a, mixed_coords(1), mass_c, &
                                mixed_coords(2), m_b, mixed_coords(3))

                                end if

                        end do
                else
                        do i = 1, n

                                ! Evenly sample entire range of each coord
                                mixed_coords(2) = mixed_lims(3) + (i-1) * width(2) / real(n)

                                do j = 1, n
                                        ! Evenly sample entire range of each coord
                                        mixed_coords(3) = mixed_lims(5) + (j-1) * width(3) / real(n)

                                        if (estimating_r_a) then
                                                test_coord(j + (i-1)*(n-1)) = calculate_r_a(mu_a, mixed_coords(1), mass_c, &
                                                mixed_coords(2), m_b, mixed_coords(3))
                                        else if (estimating_theta_a) then
                                                test_coord(j + (i-1)*(n-1)) = calculate_theta_a(small_r_a, mu_a, &
                                                mixed_coords(1), mass_c, mixed_coords(2), m_b, mixed_coords(3))
                                        end if

                                end do
                        end do
                end if

                single_lims(1) = minval(test_coord)
                single_lims(2) = maxval(test_coord)
                ! print *, single_lims

        end function estimate_jacobi_lims


        function get_limits(get_r_lims, get_theta_lims, estimate_lims, R_a, small_r_a, mixed_lims, mu_a, m_b, mass_c) result(lims)

                implicit none

                real, intent(in) :: R_a, small_r_a, mixed_lims(6), mu_a, m_b, mass_c
                logical, intent(in) :: get_r_lims, get_theta_lims, estimate_lims

                real :: lims(2)

                if (get_r_lims .and. get_theta_lims) then
                        print *, "Only one coordinate flag should be set to true"
                        stop
                end if

                if (get_r_lims) then
                        if (estimate_lims) then
                                lims(1:2) = estimate_jacobi_lims(R_a, mixed_lims, mu_a, m_b, mass_c, &
                                        .true., .false., 0.)
                        else
                                ! r_a_min uses R_a_min, R_b_min, theta_ab_max - r smaller with larger theta)
                                lims(1) = calculate_r_a(mu_a, R_a, mass_c, &
                                        mixed_lims(3), m_b, mixed_lims(6))

                                ! r_a_max depends on current (single) R_a, but uses full (mixed) R_b_max and theta_ab_min
                                lims(2) = calculate_r_a(mu_a, R_a, mass_c, &
                                        mixed_lims(4), m_b, mixed_lims(5))
                        end if

                else if (get_theta_lims) then
                        if (estimate_lims) then
                                lims(1:2) = estimate_jacobi_lims(R_a, mixed_lims, mu_a, m_b, mass_c, &
                                .false., .true., small_r_a)
                        else
                                ! theta_a_min depends on current R_a (and r_a), but uses full R_b_max and theta_ab_min
                                lims(1) = calculate_theta_a(small_r_a, mu_a, R_a, mass_c, &
                                        mixed_lims(4), m_b, mixed_lims(6))

                                ! For R_b = 0, theta_a_max = pi
                                lims(2) = calculate_theta_a(small_r_a, mu_a, R_a, mass_c, &
                                        mixed_lims(3), m_b, mixed_lims(5))
                        end if
                end if

        end function get_limits

        function get_test_limits(get_y_lims, get_z_lims, x, y) result(lims)

                implicit none

                real, intent(in) :: x, y
                logical, intent(in) :: get_y_lims, get_z_lims

                real :: lims(2)

                if (get_y_lims .and. get_z_lims) then
                        print *, "Only one coordinate flag should be set to true"
                        stop
                end if

                if (get_y_lims) then
                        ! Depends on x
                        lims(1:2) = (/ 0., 1. - x/)
                else if (get_z_lims) then
                        ! Depends on x and y
                        lims(1:2) = (/ 0., 1. - x - y /)
                end if

        end function get_test_limits

        ! Returns (x+y+z)^-0.5
        ! For 1/denom, requires adjustment for 1/0, and so a much larger n is needed
        ! Also tested with (x+y+z)^-0.5
        function test_integrand_func(x, y, z) result(integrand)

                implicit none

                real, intent(in) :: x, y, z
                real :: integrand, denominator

                denominator = sqrt(x + y + z)
                if (denominator  < 0.00000001) then
                        denominator =  0.00000001
                end if
                integrand = 1 / denominator

        end function test_integrand_func

end program jacobi