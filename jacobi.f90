program integration

        implicit none
        integer :: mixed_jacobi_n(3), single_jacobi_n(3), mc_n
        real :: mixed_jacobi_lims(6), R_a_integral, R_b_integral, theta_ab_integral, mixed_integral, &
                single_integral, triple_integral, mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b

        ! Number of points to generate for MC integration
        mc_n = 2000

        ! Number of Simpson cells to split mixed Jacobi ingegral into (R_a_n, R_b_n, theta_ab_n)
        mixed_jacobi_n = (/ 100, 100, 100 /)

        ! Number of Simpson cells to split single Jacobi integral into (R_a_n, r_a_n, theta_a_n)
        single_jacobi_n = (/ 100, 100, 100 /)

        ! Limits of integration (R_a_min, R_a_max, R_b_min, R_b_max, theta_ab_min, thata_ab_max)
        mixed_jacobi_lims = (/ 0., 5., 0., 7., 0., 4.*atan(1.) /)
        ! mixed_jacobi_lims = (/ 2., 3., 2., 3., 2., 3. /)

        ! Calculate seperable integrals for R_a, R_b and theta_ab using Simpson's rule
        R_a_integral = integrate_single_Simpson(mixed_jacobi_n(1), mixed_jacobi_lims(1), &
                mixed_jacobi_lims(2), .true.)
        R_b_integral = integrate_single_Simpson(mixed_jacobi_n(2), mixed_jacobi_lims(3), &
                mixed_jacobi_lims(4), .true.)
        theta_ab_integral = integrate_single_Simpson(mixed_jacobi_n(3), mixed_jacobi_lims(5), &
                mixed_jacobi_lims(6), .false.)

        ! Integrals are seperable so multiply for total integral
        mixed_integral = R_a_integral * R_b_integral * theta_ab_integral
        print *, "Mixed Jacobi integral = ", mixed_integral

        ! Calculate integral for R_a, R_b and theta_ab using Simpon's rule on each in turn
        triple_integral = calculate_triple_integral(mixed_jacobi_n, mixed_jacobi_lims)
        print *, "Triple mixed Jacobi integral = ", triple_integral

        ! Calculate mass relations (three masses defined within)
        call calculate_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

        ! Calculate integral for R_a, r_a and theta_a using Simpon's rule on each in turn
        single_integral = calculate_single_S_integral(single_jacobi_n, mixed_jacobi_lims, mu_a, mu_b, m_a, m_b, mass_c)
        print *, "Single Jacobi integral = ", single_integral

        ! Calculate integral for R_a, r_a and theta_a using Monte Carlo integration
        single_integral = calculate_single_mc_integral(mc_n, mixed_jacobi_lims, mu_a, mu_b, m_a, m_b, mass_c)
        print *, "Single Jacobi integral = ", single_integral

        ! call test_limits(mixed_jacobi_lims)

contains

        subroutine calculate_masses(mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b)

                implicit none

                real, intent(inout) :: mass_a, mass_b, mass_c, mass_total, mu_a, mu_b, m_a, m_b

                ! Atom masses
                mass_a = 2.
                mass_b = 3.
                mass_c = 4.

                ! Sum of masses
                mass_total = 0
                mass_total = mass_a + mass_b + mass_c

                ! Internal reduced masses
                m_a = 0
                m_b = 0
                m_a = mass_b * mass_c / (mass_b + mass_c)
                m_b = mass_a * mass_c / (mass_a + mass_c)

                ! Reduced channel masses
                mu_a = 0
                mu_b = 0
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

                width = abs(b-a)/n

                integral = 0

                do i = 0, n
                        x = a + i * width
                        y(i) = integrand_func(x, is_x_squared)

                        if (i == 0 .or. i == n) then
                                integral = integral + y(i)
                        else if (mod(i, 2) == 0) then
                                integral = integral + 2 * y(i)
                        else
                                integral = integral + 4 * y(i)
                        end if
                end do

                integral = width * integral / 3

                if (is_x_squared) then
                        !print *, "The integral of x^2 from ", a, " to ", b, " is ", integral
                else
                        !print *, "The integral of sin(x) from ", a, " to ", b, " is ", integral
                end if

        end function integrate_single_Simpson

        ! Returns either x^2 or sin(x) based on flag passed
        function integrand_func(x, is_x_squared) result(integrand)

                implicit none

                real, intent(in) :: x
                logical, intent(in) :: is_x_squared

                real :: integrand

                if (is_x_squared) then
                        integrand = x ** 2
                else
                        integrand = sin(x)
                end if

        end function integrand_func

        function calculate_r_a(mu_a, R_a, M_c, R_b, m_b, gamma_ab) result(coord)

                implicit none

                real, intent(in) :: mu_a, R_a, M_c, R_b, m_b, gamma_ab

                ! r_a
                real :: coord

                coord = mu_a * sqrt( (R_a / M_c)**2  + (R_b / m_b)**2 &
                        + 2 * (R_a * R_b * cos(gamma_ab) / (M_c * m_b)) )

        end function calculate_r_a

        function calculate_theta_a(small_r_a, mu_a, R_a, M_c, R_b, m_b, gamma_ab) result(coord)

                implicit none

                real, intent(in) :: small_r_a, mu_a, R_a, M_c, R_b, m_b, gamma_ab

                !theta_a
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

        ! Use Simpson's rule to integrate three independent variables
        function calculate_triple_integral(n, lims) result(total_integral)

                implicit none

                real, intent(in) :: lims(6)
                integer, intent(in) :: n(3)

                real :: width(3), x, y, z, corners(9), trap_volume, mid_volume, &
                temp_integral, z_integral, y_integral, total_integral
                integer :: i, j, k

                ! Width's of Simpson's rule subintervals for each coordinate
                width(1) = abs(lims(2) - lims(1)) / n(1)
                width(2) = abs(lims(4) - lims(3)) / n(2)
                width(3) = abs(lims(6) - lims(5)) / n(3)

                total_integral = 0
                temp_integral = 0

                ! Integrate each variable in turn, covering full limits
                ! Variables labelled x, y and z for simplicity
                do i = 0, n(1)

                        ! Set value for R_a
                        x = lims(1) + i * width(1)

                        ! Total R_b interegral for set R_a
                        y_integral = 0

                        do j = 0, n(2)

                                ! Set value for R_b
                                y = lims(3) + j * width(2)

                                ! Total theta_ab interegral for set R_a, R_b
                                z_integral = 0

                                do k = 0, n(3)

                                        ! Set value for theta_ab
                                        z = lims(5) + k * width(3)

                                        ! Use Simpon's rule to add contributions for this subinterval
                                        temp_integral = integrand_func(z, .false.)

                                        if (k == 0 .or. k == n(3)) then
                                                z_integral = z_integral + temp_integral
                                        else if (mod(k, 2) == 0) then
                                                z_integral = z_integral + 2 * temp_integral
                                        else
                                                z_integral = z_integral + 4 * temp_integral
                                        end if

                                end do

                                ! Total theta_ab intergral for set R_a, R_b
                                z_integral = width(3) * z_integral / 3

                                ! Use Simpon's rule to add contributions for this subinterval
                                ! Inlcudes multiplication by total theta_ab integral
                                temp_integral = z_integral * integrand_func(y, .true.)

                                if (j == 0 .or. j == n(2)) then
                                        y_integral = y_integral + temp_integral
                                else if (mod(j, 2) == 0) then
                                        y_integral = y_integral + 2 * temp_integral
                                else
                                        y_integral = y_integral + 4 * temp_integral
                                end if

                        end do

                        ! Total R_b integral for set R_a
                        y_integral = width(2) * y_integral / 3

                        ! Use Simpon's rule to add contributions for this subinterval
                        ! Includes multiplication by total R_b integral
                        temp_integral = y_integral * integrand_func(x, .true.)

                        if (i == 0 .or. i == n(1)) then
                                total_integral = total_integral + temp_integral
                        else if (mod(i, 2) == 0) then
                                total_integral = total_integral + 2 * temp_integral
                        else
                                total_integral = total_integral + 4 * temp_integral
                        end if

                end do

                ! Total (R_a) integral
                total_integral = width(1) * total_integral / 3

        end function calculate_triple_integral

        ! Use Simpson's rule to integrate three dependent variables
        function calculate_single_S_integral(n, mixed_lims, mu_a, mu_b, m_a, m_b, mass_c) result(total_integral)

                implicit none

                real, intent(in) :: mixed_lims(6), mu_a, mu_b, m_a, m_b, mass_c
                integer, intent(in) :: n(3)

                real :: width(3), x, y, z, single_lims(6), &
                temp_integral, z_integral, y_integral, total_integral
                integer :: i, j, k

                total_integral = 0
                temp_integral = 0

                ! R_a_min and R_a_max can use full mixed R_a limits
                single_lims(1) = mixed_lims(1)
                single_lims(2) = mixed_lims(2)

                width(1) = abs(single_lims(2) - single_lims(1)) / n(1)

                do i = 0, n(1)

                        x = single_lims(1) + i * width(1)
                        y_integral = 0

                        ! r_a_min uses R_a_min, R_b_min, theta_ab_max - r smaller with larger theta)
                        single_lims(3) = calculate_r_a(mu_a, x, mass_c, &
                                mixed_lims(3), m_b, mixed_lims(6))

                        ! r_a_max depends on current (single) R_a, but uses full (mixed) R_b_max and theta_ab_min
                        single_lims(4) = calculate_r_a(mu_a, x, mass_c, &
                                mixed_lims(4), m_b, mixed_lims(5))

                        width(2) = abs(single_lims(4) - single_lims(3)) / n(2)
                        ! print *, "r_a: ", single_lims(3), single_lims(4)

                        do j = 0, n(2)

                                y = single_lims(3) + j * width(2)

                                z_integral = 0

                                ! theta_a_min depends on current R_a (and r_a), but uses full R_b_max and theta_ab_min
                                single_lims(5) = calculate_theta_a(y, mu_a, x, mass_c, &
                                mixed_lims(4), m_b, mixed_lims(6))

                                ! For R_b = 0, theta_a_max = pi
                                single_lims(6) = calculate_theta_a(y, mu_a, x, mass_c, &
                                        mixed_lims(3), m_b, mixed_lims(5))

                                width(3) = abs(single_lims(6) - single_lims(5)) / n(3)
                                ! print *, "theta_a: ", single_lims(5), single_lims(6)

                                do k = 0, n(3)

                                        z = single_lims(5) + k * width(3)

                                        temp_integral = integrand_func(z, .false.)

                                        if (k == 0 .or. k == n(3)) then
                                                z_integral = z_integral + temp_integral
                                        else if (mod(k, 2) == 0) then
                                                z_integral = z_integral + 2 * temp_integral
                                        else
                                                z_integral = z_integral + 4 * temp_integral
                                        end if

                                end do

                                z_integral = width(3) * z_integral / 3

                                temp_integral = z_integral * integrand_func(y, .true.)

                                if (j == 0 .or. j == n(2)) then
                                        y_integral = y_integral + temp_integral
                                else if (mod(j, 2) == 0) then
                                        y_integral = y_integral + 2 * temp_integral
                                else
                                        y_integral = y_integral + 4 * temp_integral
                                end if

                        end do

                        y_integral = width(2) * y_integral / 3

                        temp_integral = y_integral * integrand_func(x, .true.)

                        if (i == 0 .or. i == n(1)) then
                                total_integral = total_integral + temp_integral
                        else if (mod(i, 2) == 0) then
                                total_integral = total_integral + 2 * temp_integral
                        else
                                total_integral = total_integral + 4 * temp_integral
                        end if

                end do

                total_integral = width(1) * total_integral / 3

                total_integral = ( (mu_a * mu_b) / (m_a * m_b) )**(-3/2) * total_integral

        end function calculate_single_S_integral

        ! Calculate integral with Monte Carlo
        function calculate_single_mc_integral(n, mixed_lims, mu_a, mu_b, m_a, m_b, mass_c) result(total_integral)

                implicit none

                real, intent(in) :: mixed_lims(6), mu_a, mu_b, m_a, m_b, mass_c
                integer, intent(in) :: n

                real :: width(3), coords(3), single_lims(6), &
                temp_integral, total_integral
                integer :: i

                total_integral = 0

                width = (/ 0, 0, 0 /)
                coords = (/ 0, 0, 0 /)
                single_lims = (/ 0, 0, 0, 0, 0, 0 /)

                ! R_a_min and R_a_max can use full mixed R_a limits
                single_lims(1) = mixed_lims(1)
                single_lims(2) = mixed_lims(2)

                width(1) = abs(single_lims(2) - single_lims(1))
                ! print *, "R_a: ", single_lims(1), single_lims(2)

                do i = 0, n

                        temp_integral = 0

                        ! Generate random points in range 0 to 1
                        call random_number(coords)
                        ! print *, i, ". Random numbers: ", coords(1), coords(2), coords(3)

                        ! Random R_a
                        coords(1) = single_lims(1) + coords(1) * width(1)

                        ! r_a_min uses R_a_min, R_b_min, theta_ab_max - r smaller with larger theta)
                        ! single_lims(3) = calculate_r_a(mu_a, coords(1), mass_c, &
                                ! mixed_lims(3), m_b, mixed_lims(6))

                        ! r_a_max depends on current (single) R_a, but uses full (mixed) R_b_max and theta_ab_min
                        ! single_lims(4) = calculate_r_a(mu_a, coords(1), mass_c, &
                                ! mixed_lims(4), m_b, mixed_lims(5))

                        ! print *, "r_a limits: ", single_lims(3), single_lims(4)

                        ! theta_a_min depends on current R_a (and r_a), but uses full R_b_max and theta_ab_min
                        ! single_lims(5) = calculate_theta_a(coords(2), mu_a, coords(1), mass_c, &
                        !         mixed_lims(4), m_b, mixed_lims(6))

                        ! For R_b = 0, theta_a_max = pi
                        ! single_lims(6) = calculate_theta_a(coords(2), mu_a, coords(1), mass_c, &
                        !         mixed_lims(3), m_b, mixed_lims(5))

                        ! print *, i, ". Old: ", single_lims(3:6)

                        single_lims(3:4) = estimate_jacobi_lims(coords(1), mixed_lims, mass_c, mu_a, m_b, &
                                .true., .false., 0.)
                        width(2) = abs(single_lims(4) - single_lims(3))

                        ! Random r_a
                        coords(2) = single_lims(3) + coords(2) * width(2)

                        single_lims(5:6) = estimate_jacobi_lims(coords(1), mixed_lims, mass_c, mu_a, m_b, &
                                .false., .true., coords(2))
                        width(3) = abs(single_lims(6) - single_lims(5))

                        ! print *, i, ". New: " , single_lims(3:6)

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

                total_integral = total_integral / n

                total_integral = ( (mu_a * mu_b) / (m_a * m_b) )**(-3/2) * total_integral

        end function calculate_single_mc_integral

        ! Evaluates hardcoded function based on x, y and z
        function evaluate_function(x, y, z) result(total_func)

                implicit none

                real, intent(in) :: x, y, z

                real :: total_func, F

                total_func = 0

                ! Function being integrated F(x,y,z)
                F = 1
                ! Multiply F by integration variables in integral 
                total_func = x**2 * y**2 * sin(z) * F
                ! total_func = sqrt(x+y+z)

                ! print *, "R_a, r_a, gamma_a: ", x, y, z
                ! print *, "Total func: ", total_func

        end function evaluate_function

        ! Estimate the minimum and maximum values of r_a, based on known R_a
        function estimate_jacobi_lims(R_a, mixed_lims, mass_c, mu_a, m_b, estimating_r_a, &
                estimating_theta_a, small_r_a) result(single_lims)

                implicit none

                real, intent(in) :: R_a, mixed_lims(6), mass_c, mu_a, m_b, small_r_a
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
                                mixed_coords(2) = mixed_lims(3) + (i-1) * width(2) / n

                                do j = 1, n
                                        ! Evenly sample entire range of each coord
                                        mixed_coords(3) = mixed_lims(5) + (j-1) * width(3) / n

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

end program integration