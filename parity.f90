program parity
        implicit none
        integer :: x(5)
        integer :: i

        x = (/ -1, -10, 2, 0, 5 /)

        do i = 1,5
                if (x(i) == 0) then
                        exit
                end if
                
                if (mod(x(i),2) /= 0) then
                        print *, x(i), " is odd"
                else
                        print *, x(i), " is even"

                end if
        end do
end program parity

