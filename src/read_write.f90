!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This module takes input from the user in order to 
!! build the reactor and calculate the measured flux with various 
!! defects and/or detector locations. 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use neutron_flux, only : box_flux_booles, large_x0_flux, box_flux_monte_carlo, total_flux_booles, hollow_box_flux_mc

implicit none

private
public :: read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine takes in information provided by the user
!! in order to measure the flux generated a, user defined, variable
!! distance form a nuclear reactor with a user defined size. Listed 
!! below are the inputs taken from the user, with a descrption of its use.  
!!----------------------------------------------------------------------
!! Output:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(out) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(out) ::  n_grid, m_samples

    !
    print *, 'This program calculates neutron flux of a reactor ....'
    print*, 'You get to define the size of the reactor.'
    print*, 'Flux will be measured a varying distance from reactor'
    print*, 'You define the range of distances the detector can be from the reactor.'

    ! Function are called below to ask for this user input, each input goes
    ! through a do loop in order to verify the inout is of the correct type. 
    depth = read_real('depth D')
    width = read_real('width W')
    height = read_real('height H')
    y_zero = read_real('y coordinate of detector position')
    x_min = read_real('minumum x value of detector position')
    x_max = read_real('maximun x value of detector position')
    x_step = read_real('incriment size for x coordinates of detector position')

    ! The function read_real used above returns a double precision real,
    ! however n_grin and m_samples are integers, that means that we need
    ! another function to get integers. For that we'll define read_integer
    ! as well
    ! These functions are called to ask user for integer values and run them through
    ! a similar do loop to verify they of are the correct type. 
    n_grid = read_integer('number of lattice points N')
    m_samples = read_integer('number Monte Carlo samples M')

    print *,
end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Function containing a do loop in order to verfiy the input
!! given by the user is of the correct type. These will be real numbers.
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        real        A positive non negative number given by the user
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

     do
        read(*,'(a)',iostat=ierror) string
        ! In input is not empty, proceed
        if(string.ne.'') then
            read(string,*,iostat=ierror) x
            ! If input can be made into a number, proceed
            if (ierror == 0) then
                ! if number is positive, we can exit the loop.
                if (x > 0) exit       
                    print *, "'"//trim(string)//"'"// 'cannot be negative or zero, please provide a positive number' 
            else
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif          
        else
            print *, 'that was an empty input, please provide a positive, non-zero, number'
        
        endif
    enddo


end function read_real

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Function containing a do loop in order to verfiy the input
!! given by the user is of the correct type. These will be integers. 
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        integer     A positive non negative number given by the user
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *,
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'
  
     do
        read(*,'(a)',iostat=ierror) string
        if(string.ne.'') then
            read(string,*,iostat=ierror) x
            if (ierror == 0) then
                if (x > 0) exit       
                    print *, "'"//trim(string)//"'"// 'cannot be negative or zero, please provide a positive numebr' 
            else
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif          
        else
            print *, 'that was an empty input, please provide a positive, non-zero, number'
        
        endif
    enddo
end function read_integer

!-----------------------------------------------------------------------
!Subroutine: write_neutron_flux
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine writes the fluxes generated later in
!! the program (for the solid box) into a file named results_basic.dat.
!! This also calles for the calculation of these fluxes for each detector 
!! location. 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(in) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(in) :: n_grid, m_samples

    real(dp) :: x_zero, box_booles, box_mc, box_large_x0, sigma_box
    character(len=*), parameter :: file_name = 'results_basic.dat'
    integer :: unit

    open(newunit=unit,file=file_name)
    write(unit,'(5a28)') 'x_0', 'booles', 'large x_0', 'monte carlo', 'MC uncertainty'
    x_zero = x_min
    do 
        if(x_zero > x_max) exit
        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
        box_large_x0 = large_x0_flux(depth, width, height, x_zero, y_zero)
        call box_flux_monte_carlo(depth, width, height, x_zero, y_zero, m_samples, box_mc, sigma_box)
        write(unit,'(5e28.16)') x_zero, box_booles, box_large_x0, box_mc, sigma_box
        x_zero = x_zero + x_step
    enddo
    close(unit)
    print *, 'The fluxes were written in the '//file_name//' file'
    print *,
end subroutine write_neutron_flux

!-----------------------------------------------------------------------
!! Subroutine: read_advanced_input
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine takes the input from the user for the 
!! second part of the program, calculation the flux of a box reactor with a
!! varialbe size hollow spherical region. Same as the subroutine for the
!! solid box portion of the program, this routine sends the input to get checked
!! in the read_real routine.
!!----------------------------------------------------------------------
!! Output:
!!
!! x_zero       real        x coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!-----------------------------------------------------------------------
 subroutine read_advanced_input(x_zero, r_min, r_max, r_step)
     implicit none
     real(dp), intent(out) :: x_zero, r_min, r_max, r_step

print*, 'This portion of the program calculates flux of a reactor with a hollow spherical chamber inside.'
print*, 'This is the same box reactor as in the first part of the program, only with the hollow sphere inside. '
print*, 'You determing the varying size of the hollow region withing the reactor.'

x_zero = read_real('x coordinate of detector position.')
r_min = read_real('minumum radius of hollow sphere.')
r_max = read_real('maximum redius of hollow sphere.')
r_step = read_real('increment size for the radius of hollow sphere.')

 end subroutine read_advanced_input

!-----------------------------------------------------------------------
!Subroutine: write_advanced_flux
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine writes the fluxes calculated for the hollow
!! box reactor into a file named results_advanced.dat. It also called for the 
!! individual flux's to be calculated with a radius that varies from the user's
!! minumum value to the user's maximun value via user defind steps. 
!!----------------------------------------------------------------------
!! Input:
!! 
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
 subroutine write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)
     implicit none
     real(dp), intent(in) :: depth, width, height, x_zero, y_zero, r_min, r_max, r_step
     integer, intent(in) :: n_grid, m_samples

     real(dp) :: radius, box_booles, hollow_booles, hollow_mc, sigma_hollow
     character(len=*), parameter :: file_name = 'results_advanced.dat'
     integer :: unit

     open(newunit=unit, file=file_name)
     write(unit,'(5a28)') 'radius', 'box booles', 'hollow booles', 'hollow monte carlo', 'MC uncertainty'
     radius = r_min
    
    do 
    
     if (radius > r_max) exit
      box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
      hollow_booles = total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid)
      call hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, m_samples, hollow_mc, sigma_hollow)
      write(unit,'(5e28.16)') radius, box_booles, hollow_booles, hollow_mc, sigma_hollow
      radius = radius + r_step
    enddo
    close(unit)
    print *, 'The fluxes were written in the '//file_name//' file'
    print *,
 end subroutine write_advanced_flux

end module read_write
