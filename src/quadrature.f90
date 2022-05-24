!-----------------------------------------------------------------------
!Module: quadrature
!-----------------------------------------------------------------------
!! By Louis Andre
!!
!! Description of module and contents: This module houses the functions that
!! initialize and execute the volume integral approximation method called Boole's Quadrature,
!! as well as an alternative integral approximation method using random rumber
!! generators, aptly names monte_carlo_quad. 
!! This monte_carlo_quad method takes an integral and exapnds it into a sum of fx() functions,
!! multiplies by the volume of the space over which the integration occurs, and 
!! divides by the number of random points withing the volume.
!! 
!! This module receives data and functions from the neutron flux module and executes integral 
!! approximations of those functions over a space specified by the data input.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! 
!! monte_carlo_quad
!!----------------------------------------------------------------------
!! Included functions:
!!
!! booles_quadrature
!! booles_rule
!-----------------------------------------------------------------------
module quadrature
use types

implicit none

private
public :: booles_quadrature, monte_carlo_quad

!-----------------------------------------------------------------------
!Interface: func
!-----------------------------------------------------------------------
!! This defines a new type of procedure in order to allow callbacks
!! in the Monte Carlo quadrature subroutine of an arbitrary function that is given
!! as input and declared as a procedure
!!
!! The arbitrary function receives two rank 1 arrays of arbitrary size.
!! The first array contains an n-dimensional vector representing the
!! point sampled by the Monte Carlo method. The second is a "work array"
!! that contains parameters  necessary to calculate the function to be
!! integrated.
!!----------------------------------------------------------------------
interface
    real(dp) function func(x, data)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: x(:), data(:)
       
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Function: booles_quadrature
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This function takes the incomming array and cuts it into
!! smaller slices to be sent to the booles_rule function, which calculates
!! the expanson approximation to the integral. 
!! This function sums the results of those slices after they return from
!! the booles_rule function, as well as does the necessary checks in order
!! to verify the incomming array is of the correct dimension for Booles 
!! quadrature to function properly. 
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size, i

    fx_size = size(fx)

    

    ! As the diagram below shows, only certain number of grid points
    ! fit the scheme of Boole's quadrature. The following If statement 
    ! ensures that the input array satisfies the condition of being
    ! divisible by 4 when one is subtracted.

    ! |--interval 1---|--interval 2---|--interval 3---|
    ! 1   2   3   4   5   6   7   8   9   10  11  12  13
    ! |---|---|---|---|---|---|---|---|---|---|---|---|
    ! x0  x1  x2  x3  x4
    !                 x0  x1  x2  x3  x4
    !                                 x0  x1  x2  x3  x4


     if (modulo(fx_size-1, 4) /= 0)  then 
         print *, 'fx array size plus one in booles_quadrature has to be divisible by 4'
         stop
     endif

    ! We could implement the full integration here, however to make a cleaner,
    ! easy to read (and debug or maintain) code we will define a smaller
    ! function that returns Boole's five point rule and pass slices (1:5), (5:9),
    ! (9:13), ... of fx to such function to then add all the results. 

    s = 0._dp

     do i = 1, ((fx_size-1)/4)
        s = s + booles_rule(fx( ((4*i) - 3):( (4*i) + 1)), delta_x)
        
     enddo
     
end function booles_quadrature

!-----------------------------------------------------------------------
!! Function: booles_rule
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This function unpacks the slices sent from booles_quadrature
!! and uses the values in the Boole's integral approximation partial summation
!! represented as the variable "s" below. This value is then sent back to 
!! booles_quadrature in order to take their place as a piece if the summation.
!! 
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size
    real(dp) :: fx0, fx1, fx2, fx3, fx4

    fx_size = size(fx)

    ! Additional check to make sure we received 5 values for the partial sum below. 

     if (fx_size /= 5) then
        print*, 'Error in calculation in booles_rule, array fx must have 5 points and is currently ', fx_size,'.' 
        stop  
     endif
    
     fx0 = fx(1)
     fx1 = fx(2)
     fx2 = fx(3)
     fx3 = fx(4)
     fx4 = fx(5)

     
    
     s = ((delta_x*2)/45._dp)*(7*fx0 + 32*fx1 + 12*fx2 + 32*fx3 + 7*fx4)

end function booles_rule

!-----------------------------------------------------------------------
!! Subroutine: monte_carlo_quad
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This is the "meat and potatoes"of the monte_carlo_quad
!! method os approximating integrals through sums and random point generation. 
!! The method is based on prebability of a point landing within the volume 
!! over whicb an integral is evaluated. Summation of the function evaluated
!! at random points within the volume, multiplied by the volume and 
!! divided by the number of sampling points, gives an approximation to 
!! the integral. 
!! ----------------------------------------------------------------------
!! Input:
!!
!! f            procedure   function to be integrated
!! a            real        array containing the lower limits of the integral
!! b            real        array containing the upper limits of the integral
!! data         real        array containing parameters necessary to calculate the function f
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Monte Carlo integral
!! sigma_s      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine monte_carlo_quad(f, a, b, data, n_samples, s, sigma_s)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: a(:), b(:), data(:)
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: s, sigma_s
    real(dp) :: avg_square, square_avg, var_f, volume

    integer :: i, vector_size, b_size
    real(dp), allocatable :: x_vector(:), fx(:), delta_ab(:)
    

    
    vector_size = size(a)
    b_size = size(b)

    ! Defining a Monte Carlo routine that works for an arbitrary number of 
    ! dimensions in the integral. The advantage of Monte Carlo integration,
    ! it's very efficient for high dimensional integrals)

    ! Since a and b give the lower and upper limits they need to have the same size.
    ! Make a check to see if they do have the same size

    if (vector_size /= b_size) then 
          print *, 'a and b arrays in monte_carlo_quad have to be the same size'
         stop       
    endif

    ! Here we allocate memory for the vector containing the sample points and 
    ! for a vector that contains the evaluated function
    allocate(x_vector(1:vector_size))
    allocate(fx(1:n_samples))
    allocate(delta_ab(1:vector_size))

    ! Evaluates the volume of the region using dimentions stored in "a" and "b" arrays. 
    delta_ab = b - a

    volume = 1
    do i=1,vector_size 
        
        volume = volume*(delta_ab(i))
        
    enddo

    s = 0
    square_avg = 0
    avg_square = 0

    do i=1,n_samples
        !generates an array with random numbers in the [0,1) interval
        call random_number(x_vector) 
        ! rescaling to the integration volume [a,b), random points will now represent
        ! a random three dimensional point in the volume of the region. 
        x_vector = a + x_vector*delta_ab  
        ! This random sampling point is now evaluated in the function sent to this subroutine. 
        ! Allowing us to use this subroutine for any function, very useful. 
        fx(i) = f(x_vector,data)
        ! This computes a term of the summation using the function evaluated at the random 
        ! point declared above. This is multiplied by the volume divided by the number
        ! of random points sampled. This is then stored in "s" to be added to the next term 
        ! of the sum when the loop cycles again. 
        s = s + (volume/((float(n_samples))))*(fx(i))

        ! This computes the average of the function squared, to be used in variance calculation. 
        square_avg = square_avg + (1/(float(n_samples)))*(fx(i))*(fx(i))
        ! This computes the square of the average of the function, to be used in the variance calculation.
        avg_square = avg_square + (1/ (float(n_samples)) )*(fx(i))

    enddo
    ! This calculates the variance in the function evaluated at the numerous sample points. 
    var_f = square_avg - (avg_square**2)
    ! This calculates the standard deviation from the mean of the function evaluated at the 
    ! numerous sample points. This is the error calculated for any monte carlo method in the program. 
    ! Also called the uncertainty. 
    sigma_s = volume * (sqrt(var_f)/sqrt( float(n_samples) ))

end subroutine monte_carlo_quad

end module quadrature