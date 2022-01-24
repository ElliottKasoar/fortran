program format
        implicit none

        real :: pi
        character(len=7) :: x

        pi = 4.0 * atan2(1.0,1.0)
        print "(E12.2)", pi
        print "(E12.4)", pi
        print "(E12.6)", pi
        print "(F12.2)", pi
        print "(F12.4)", pi
        print "(F12.6)", pi
        print "(G12.2)", pi
        print "(G12.4)", pi
        print "(G12.6)", pi
        print *, pi

        x = "(e12.2)"
        !x = "abc"
        !print *, x(:)
        print x, pi
end program format
