program integration

        implicit none
        integer :: n
        real :: a, b

        n = 6
        a = 0
        b = 5

        if (mod(n,2) == 0 .and. n > 0) then
                call calculate_integral(n, a, b)
        else
                print *, "n must be an even, positive integer"
        end if

contains

        subroutine calculate_integral(n, a, b)
                implicit none
                integer, intent(in) :: n
                real, intent(in) :: a, b

                real :: width, x, total
                real :: y(0:n)
                integer :: i

                width = abs(b-a)/n

                total = 0

                do i = 0, n
                        x = a + i * width
                        y(i) = integrand_func(x)

                        if (i == 0 .or. i == n) then
                                total = total + y(i)
                        else if (mod(i, 2) == 0) then
                                total = total + 2 * y(i)
                        else
                                total = total + 4 * y(i)
                        end if
                end do

                total = width * total / 3

                print *, "The integral of x^2 from ", a, " to ", b, " is ", total

        end subroutine calculate_integral

        function integrand_func(x) result(integrand)
                implicit none
                real, intent(in) :: x

                real :: integrand

                integrand = x ** 2

        end function integrand_func

end program integration
