program sum

        implicit none
        integer :: arr(6)
        integer :: n, total

        n = 6
        arr = [2, 6, 8, 10, 12, 40]
        print *, add_ints(n, arr)

        total = 0
        call add_ints_sub(n, arr, total)
        print *, total

contains

        subroutine add_ints_sub(n, arr, total)
                implicit none
                integer, intent(in) :: n
                integer, intent(in) :: arr(n)
                integer, intent(inout) :: total
                integer :: i

                do i = 1, n
                       total += arr(i)
                end do
        end subroutine add_ints_sub

        function add_ints(n, arr) result(total)
                implicit none
                integer, intent(in) :: n
                integer, intent(in) :: arr(n)

                integer :: total, i

                total = 0

                do i = 1, n
                       total = total + arr(i)
                end do
        end function add_ints

end program sum
