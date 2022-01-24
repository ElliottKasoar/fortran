program integration

        implicit none
        integer :: R_a_n, R_b_n, theta_ab_n
        real :: R_a_min, R_a_max, R_b_min, R_b_max, theta_ab_min, theta_ab_max, &
                R_a_integral, R_b_integral, theta_ab_integral, total_integral

        ! Number of cells to split into
        R_a_n = 100
        R_b_n = 100
        theta_ab_n = 1000
        
        R_a_min = 0
        R_a_max = 5
        R_b_min = 0
        R_b_max = 7
        theta_ab_min = 0
        theta_ab_max = 4.D0*atan(1.D0)

        if (mod(R_a_n,2) == 0 .and. R_a_n > 0 &
                .and. mod(R_b_n,2) == 0 .and. R_b_n > 0 &
                .and. mod(theta_ab_n,2) == 0 .and. theta_ab_n > 0) then
                R_a_integral = calculate_integral(R_a_n, R_a_min, R_a_max, .true.)
                R_b_integral = calculate_integral(R_b_n, R_b_min, R_b_max, .true.)
                theta_ab_integral = calculate_integral(theta_ab_n, theta_ab_min, theta_ab_max, .false.)

                total_integral = R_a_integral * R_b_integral * theta_ab_integral

                print *, "The integral is ", total_integral
        else
                print *, "n must be an even, positive integer"
        end if

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
                        print *, "The integral of x^2 from ", a, " to ", b, " is ", integral
                else
                        print *, "The integral of sin(x) from ", a, " to ", b, " is ", integral
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

end program integration
