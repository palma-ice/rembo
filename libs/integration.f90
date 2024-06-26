module integration

    use precision, only : wp, TOL_UNDERFLOW

    implicit none

    private
    public :: integrate_trapezoid1D_pt

contains

    function integrate_trapezoid1D_pt(var,zeta) result(var_int)
        ! Integrate a variable from the base to height zeta(nk).
        ! The value of the integral using the trapezoid rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a point of integrated value of var at level zeta(nk).

        implicit none

        real(wp), intent(IN) :: var(:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int

        ! Local variables 
        integer :: k, nk
        real(wp) :: var_mid 
        
        nk = size(var,1)

        ! Initial value is zero
        var_int = 0.0_wp 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
            var_mid = 0.5_wp*(var(k)+var(k-1))
            if (abs(var_mid) .lt. TOL_UNDERFLOW) var_mid = 0.0_wp 
            var_int = var_int + var_mid*(zeta(k) - zeta(k-1))
        end do

        return

    end function integrate_trapezoid1D_pt

end module integration