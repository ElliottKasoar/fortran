module calc_area
        implicit none
        private
        real, public, parameter :: pi = 3.14159
        public :: area

contains
        function area(r) result(calculated_area)
              real, intent(in) :: r

              real :: calculated_area

              calculated_area = pi * r**2

      end function area

end module calc_area
