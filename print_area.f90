program print_area

        use calc_area

        implicit none
        real :: radius

        print *, 'Enter radius of circle:'
        read(*,*) radius

        print *, "Area: ", area(radius)

end program print_area
