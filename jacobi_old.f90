program integration

        implicit none
        integer :: mixed_jacobi_n(3), single_jacobi_n(3), mc_n
        real :: mixed_jacobi_limits(6), R_a_integral, R_b_integral, theta_ab_integral, mixed_integral, &
                single_integral

        ! Number of points to generate for MC
        mc_n = 100000

        ! Number of Simpon cells to split mixed Jacobi ingegral into (R_a_n, R_b_n, theta_ab_n)
        mixed_jacobi_n = (/ 100, 100, 1000 /)

        ! Number of Simpon cells to split single Jacobi integral into (R_a_n, r_a_n, theta_a_n)
        single_jacobi_n = (/ 100, 100, 1000 /)

        ! Limits of integration (R_a_min, R_a_max, R_b_min, R_b_max, theta_ab_min, thata_ab_max)
        mixed_jacobi_limits = (/ 0., 3., 0., 5., 0., 2.*atan(1.) /)

        if (mod(mixed_jacobi_n(1),2) == 0 .and. mixed_jacobi_n(1) > 0 &
                .and. mod(mixed_jacobi_n(2),2) == 0 .and. mixed_jacobi_n(2) > 0 &
                .and. mod(mixed_jacobi_n(3),2) == 0 .and. mixed_jacobi_n(3) > 0) then
                R_a_integral = calculate_integral(mixed_jacobi_n(1), mixed_jacobi_limits(1), &
                        mixed_jacobi_limits(2), .true.)
                R_b_integral = calculate_integral(mixed_jacobi_n(2), mixed_jacobi_limits(3), &
                        mixed_jacobi_limits(4), .true.)
                theta_ab_integral = calculate_integral(mixed_jacobi_n(3), mixed_jacobi_limits(5), &
                        mixed_jacobi_limits(6), .false.)

                mixed_integral = R_a_integral * R_b_integral * theta_ab_integral

                print *, "Mixed Jacobi integral = ", mixed_integral
        else
                print *, "n must be an even, positive integer"
        end if

        single_integral = calculate_single_integral(mixed_jacobi_limits, single_jacobi_n, mc_n)
        print *, "Single Jacobi integral = ", single_integral

contains

        function calculate_integral(n, a, b, is_x_squared) result(integral)
                implicit none
                integer, intent(in) :: n
                real, intent(in) :: a, b
                logical, intent(in) :: is_x_squared

                real :: width, x, integral
                real :: y(0:n)
                integer :: i

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
        end function calculate_integral

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

        function calculate_single_integral(mixed_coords_lims, n, mc_n) result(total_integral)
  
                implicit none
                real, intent(in) :: mixed_coords_lims(6)
                integer, intent(in) :: n(3), mc_n

                real :: single_coords_lims(6), single_coords_max_lims(6), frac, & 
                        temp_integral, total_integral, mass_a, mass_b, mass_c, mass_total, mu_a, &
                        mu_b, m_a, m_b, temp_r_a, R_a_integral, small_r_a_integral, theta_a_integral

                ! Atom masses
                mass_a = 2.
                mass_b = 3.
                mass_c = 4.

                ! Sum of masses
                mass_total = mass_a + mass_b + mass_c
                
                ! Internal reduced masses
                m_a = mass_b * mass_c / (mass_b + mass_c)
                m_b = mass_a * mass_c / (mass_a + mass_c)

                ! Reduced channel masses
                mu_a = mass_a * mass_b * mass_c / (mass_total * m_a)
                mu_b = mass_a * mass_b * mass_c / (mass_total * m_b)

                ! Limits of integration for single: (R_a_min, R_a_max, r_a_min, r_a_max, theta_a_min, thata_a_max)
                ! See also: mixed: (R_a_min, R_a_max, R_b_min, R_b_max, theta_ab_min, thata_ab_max)
                single_coords_max_lims = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

                ! R_a
                single_coords_max_lims(1) = mixed_coords_lims(1)
                single_coords_max_lims(2) = mixed_coords_lims(2)
               
                ! r_a_min (uses R_a_min, R_b_min, theta_ab_max - r smaller with larger theta)
                single_coords_max_lims(3) = calculate_r_a(mu_a, mixed_coords_lims(1), mass_c, &
                        mixed_coords_lims(3), m_b, mixed_coords_lims(6))

                ! r_a_max (uses R_a_max, R_b_max, theta_ab_min - r larger with smaller theta)
                single_coords_max_lims(4) = calculate_r_a(mu_a, mixed_coords_lims(2), mass_c, &
                        mixed_coords_lims(4), m_b, mixed_coords_lims(5))

                ! theta_a_min
                ! When cos(theta_a) is max (uses R_a_max, R_b_max, theta_ab_max, r_a)
                single_coords_max_lims(5) = calculate_theta_a(mu_a, mixed_coords_lims(2), mass_c, &
                        mixed_coords_lims(4), m_b, mixed_coords_lims(6))

                ! theta_a_max
                ! When cos(theta) is min (uses R_a_max, R_b_max, theta_ab_min, r_a)
                single_coords_max_lims(6) = calculate_theta_a(mu_a, mixed_coords_lims(2), mass_c, &
                        mixed_coords_lims(4), m_b, mixed_coords_lims(5))

                print *, "Maximum limits: ", single_coords_max_lims
                
                if (mod(n(1),2) == 0 .and. n(1) > 0 &
                        .and. mod(n(2),2) == 0 .and. n(2) > 0 &
                        .and. mod(n(3),2) == 0 .and. n(3) > 0) then

                        R_a_integral = calculate_integral(n(1), single_coords_max_lims(1), &
                                single_coords_max_lims(2), .true.)

                        small_r_a_integral = calculate_integral(n(2), single_coords_max_lims(3), &
                                single_coords_max_lims(4), .true.)
               
                        theta_a_integral = calculate_integral(n(3), single_coords_max_lims(5), &
                                single_coords_max_lims(6), .false.)
                end if

                temp_integral = R_a_integral * small_r_a_integral * theta_ab_integral

                temp_integral = ( (mu_a * mu_b) / (m_a * m_b) )**(-3/2) * temp_integral

                print *, "Single Jacobi integral with maximum limits: ", temp_integral

                frac = mc(mc_n, single_coords_max_lims, mixed_coords_lims, mu_a, mass_c, m_b)

                total_integral = frac * temp_integral

        end function calculate_single_integral

        function mc(n, single_coords_max_lims, mixed_coords_lims, mu_a, mass_c, m_b) result(frac)
       
                implicit none

                real, intent(in) :: single_coords_max_lims(6), mixed_coords_lims(6), mu_a, mass_c, m_b
                integer, intent(in) :: n
                real :: single_coords_lims(6), single_coords(3), coords_in_v, frac
                integer :: i

                ! Number within correct range
                coords_in_v = 0

                ! Limits of integration for single: (R_a_min, R_a_max, r_a_min, r_a_max, theta_a_min, thata_a_max)
                ! See also: mixed: (R_a_min, R_a_max, R_b_min, R_b_max, theta_ab_min, thata_ab_max)
                single_coords_lims = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

                ! Create n points and test if each is within correct range
                do i = 1, n
                        !generate random points in range
                        call random_number(single_coords)

                        ! Coordinates within maximum range (all have 0 minimum, else would shift range)
                        single_coords(1) = single_coords(1) * single_coords_max_lims(2)
                        single_coords(2) = single_coords(2) * single_coords_max_lims(4)
                        single_coords(3) = single_coords(3) * single_coords_max_lims(6)

                        !print *, "Coordinates generated: ", single_coords

                        ! R_a_min and R_a_max can use full mixed R_a limits
                        single_coords_lims(1) = mixed_coords_lims(1)
                        single_coords_lims(2) = mixed_coords_lims(2)

                        ! r_a_min uses R_a_min, R_b_min, theta_ab_max - r smaller with larger theta)
                        single_coords_lims(3) = calculate_r_a(mu_a, mixed_coords_lims(1), mass_c, &
                                mixed_coords_lims(3), m_b, mixed_coords_lims(6))

                        ! r_a_max depends on current (single) R_a, but uses full (mixed) R_b_max and theta_ab_min
                        single_coords_lims(4) = calculate_r_a(mu_a, single_coords(1), mass_c, &
                                mixed_coords_lims(4), m_b, mixed_coords_lims(5))
              
                        ! theta_a_min depends on current R_a and r_a, but uses full R_b_max and theta_ab_min
                        single_coords_lims(5) = calculate_theta_a(mu_a, mixed_coords_lims(1), mass_c, &
                                mixed_coords_lims(4), m_b, mixed_coords_lims(6))

                        ! For R_b = 0, theta_a_max = pi
                        single_coords_lims(6) = calculate_theta_a(mu_a, single_coords(1), mass_c, &
                                mixed_coords_lims(3), m_b, mixed_coords_lims(5))

                        !print *, single_coords(1), mixed_coords_lims(3), mixed_coords_lims(5), single_coords(2)
                        !print *, "Limits: ", single_coords_lims

                        if (single_coords(1) > single_coords_lims(1) & 
                            .and. single_coords(1) < single_coords_lims(2) &
                            .and. single_coords(2) > single_coords_lims(3) & 
                            .and. single_coords(2) < single_coords_lims(4) & 
                            .and. single_coords(3) > single_coords_lims(5) & 
                            .and. single_coords(3) < single_coords_lims(6)) then  
                                coords_in_v = coords_in_v + 1
                        end if

                end do

                frac = coords_in_v / n

                print *, coords_in_v
                print *, "Fraction: ", frac
        
        end function mc

        function calculate_r_a(mu_a, R_a, M_c, R_b, m_b, gamma_ab) result(coord)

                implicit none
                real, intent(in) :: mu_a, R_a, M_c, R_b, m_b, gamma_ab
                
                ! r_a
                real :: coord

                coord = mu_a * sqrt( (R_a / M_c)**2  + (R_b / m_b)**2 &
                        + 2 * (R_a * R_b * cos(gamma_ab) / (M_c * m_b)) )

        end function calculate_r_a

        function calculate_theta_a(mu_a, R_a, M_c, R_b, m_b, gamma_ab) result(coord)

                implicit none
                real, intent(in) :: mu_a, R_a, M_c, R_b, m_b, gamma_ab
               
                !theta_a
                real :: temp_r_a, cos_numerator, cos_denominator, cos_coord, coord

                temp_r_a = calculate_r_a(mu_a, R_a, M_c, R_b, m_b, gamma_ab)

                if (temp_r_a < 0.00000001) then
                        cos_denominator = temp_r_a + 0.00000001
                else
                        cos_denominator = temp_r_a
                end if

                cos_numerator = - mu_a * ( (R_a / M_c) + (R_b / m_b) * cos(gamma_ab) )

                cos_coord = cos_numerator / cos_denominator

                !print *, "r_a, numerator = ", denominator, numerator
                !print *, "cos_theta = ", cos_coord
                
                coord = acos(cos_coord)
                !coord = cos_coord

        end function calculate_theta_a

end program integration
