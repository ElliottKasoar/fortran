program average

        implicit none
        integer :: arr(6)
        integer :: n
        real :: total, average_num

        n = 6
        arr = [1, 6, 8, 10, 12, 40]
        print *, add_ints(n, arr)

        total = 0.
        average_num = 0.
        call add_ints_sub(n, arr, total)
        average_num = total / real(n)
        print *, "Sum of values: ", total
        print *, "Number of values: ", n
        print *, "Average of values (real): ", average_num
        average_num = total / n
        print *, "Average of values: ", average_num

contains

        subroutine add_ints_sub(n, arr, total)
                implicit none
                integer, intent(in) :: n
                integer, intent(in) :: arr(n)
                real, intent(inout) :: total
                integer :: i

                do i = 1, n
                       total = total + arr(i)
                end do
        end subroutine add_ints_sub

        function add_ints(n, arr) result(total)
                implicit none
                integer, intent(in) :: n
                integer, intent(in) :: arr(n)

                integer :: i
                real :: total

                total = 0.

                do i = 1, n
                       total = total + arr(i)
                end do
        end function add_ints

end program average
