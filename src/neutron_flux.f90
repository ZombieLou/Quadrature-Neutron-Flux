!-----------------------------------------------------------------------
!Module: neutron_flux
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Module contains subroutines that calculate flux values
!! for the box reactor with and without a hollow spherical chamber with,
!! centered at the senter of the box, with radius = "radius"
!! Individual descriptions of each subroutine listed above each subroutine below. 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! box_flux_monte_carlo
!! hollow_box_flux_mc
!!----------------------------------------------------------------------
!! Included functions:
!!
!! box_flux_booles
!! sphere_flux_booles
!! total_flux_booles
!! sphere_flux_kernel
!! flux_kernel
!! flux_kernel_vector
!! hollow_box_flux_kernel
!! large_x0_flux
!-----------------------------------------------------------------------
module neutron_flux
use types
use quadrature, only : booles_quadrature, monte_carlo_quad

implicit none

private
public :: box_flux_booles, large_x0_flux, box_flux_monte_carlo , hollow_box_flux_mc, total_flux_booles

contains



!-----------------------------------------------------------------------
!! Function: box_flux_booles
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine is meant to calculate the triple integral for the volume 
!! of the box for which we calculate the flux. Each step of the loop fills an array
!! to be sent to the quadrature module in order to calculate that section of the 
!! triple integral. Then flux is calculated with final array sent to 
!! booles quadrature.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function box_flux_booles(depth, width, height, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_grid
    
    real(dp) :: delta_x, delta_y, delta_z
    real(dp), allocatable :: f_x(:), g_xy(:), h_xyz(:)
    integer :: n_bins, i_x, i_y, i_z
    real(dp) :: x, y, z


    ! bins:              1   2   3   4   5   6   7   8
    ! x interval:      |---|---|---|---|---|---|---|---|
    ! grid points:     1   2   3   4   5   6   7   8   9 
    
    ! interval length: |-------------------------------|
    !                  0                               depth
    ! delta x length:  |---|
    !                  0   delta_x

    !Here n_bins is the amount of intervals/delta_x's we will have.

    !The size of the steps in x, y, and z (called delta_x, delta_y, and delta_z)
    !are definied by the entire possible distance divided by the amount on intervals.

     n_bins = n_grid - 1

     delta_x = depth/(float(n_bins))
     delta_y = width/(float(n_bins))
     delta_z = height/(float(n_bins))


    ! Here we allocate arrays containing the portions of Boole's Rule 
    ! to be sent to the Boole's Quadrature subroutine 
    allocate(  f_x(1:n_grid))
    allocate( g_xy(1:n_grid))
    allocate(h_xyz(1:n_grid))

    ! Before implimenting the do loop for the integration, the starting
    ! x value must be set in order to start at x=0.
 x = -delta_x

    do i_x = 1, n_grid
        !Outside loop for x-direction (depth), terminated if x
        !reaches the max distance over x, which is the user input "depth".
        
         if (x < depth) then
            x = x + delta_x
         else 
            exit
         endif
        ! y set to negative delta_y (step in y) in order for integration
        ! to start at y=0.
        y=-delta_y

        do i_y = 1, n_grid
        !Middle loop for y-direction (width), terminated if y
        !reaches the max distance over y, which is user input "width".

        if (y<width) then
            y = y + delta_y
        else
            exit
        endif
        ! z set to negative delta_z (step in z) in order for integration
        ! to start at z=0.
        z = -delta_z

            do i_z = 1, n_grid
        !Inside loop for z-direction (heigth), terminated if z
        !reaches the max distance over z, which is user input "heigth".
                
                if (z<height) then
                z = z + delta_z
            else 
                exit
            endif
                ! Now you can fill the h_xyz array with an evaluation of the function
                ! to integrate. To make things cleaner will define a function 
                ! below that returns the value we want

                !flux_kernel contains integrand of the triple integral.  
                !array filled to be send to booles quadrature for integral approximation
                h_xyz(i_z) = flux_kernel(x, y, z, x_zero, y_zero)

            enddo
            ! filling array by sending prevous array and delta_z to booles quadrature
            ! which will approximate integral over z.
             g_xy(i_y) = booles_quadrature(h_xyz, delta_z)

        enddo
            ! filling array by sending prevous array and delta_y to booles quadrature
            ! which will approximate integral over y.   
        f_x(i_x) = booles_quadrature(g_xy, delta_y)

    enddo
        ! filling array by sending prevous array and delta_x to booles quadrature
        ! which will approximate integral over x. This is last step is triple
        ! integral, and will provide the flux for the solid box case at 
        ! a particular detector position. 

     flux = booles_quadrature(f_x, delta_x)

end function box_flux_booles

!-----------------------------------------------------------------------
!! Function: flux_kernel
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This is the subroutine that computed the integrand 
!! from box_flux_booles.
!!
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x            real        x coordinate of the small integration volume
!! y            real        y coordinate of the small integration volume
!! y            real        z coordinate of the small integration volume
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function flux_kernel(x, y, z, x0, y0) result(k)
    implicit none
    real(dp), intent(in) :: x, y, z, x0, y0
    
     k = (1/(4._dp*pi*( ((x+x0)**2) + ((y-y0)**2) + ((z)**2) ) ) )

end function flux_kernel

!-----------------------------------------------------------------------
!! Function: large_x0_flux
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine calculates the flux assuming the source is
!! a point source (aka, far away from detector), hense the "large x0"
!!----------------------------------------------------------------------
!! Input:
!!
!! d            real        Depth of the rectangular nuclear reactor
!! w            real        Width of the rectangular nuclear reactor
!! h            real        Height of the rectangular nuclear reactor
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!!----------------------------------------------------------------------
real(dp) function large_x0_flux(d, w, h, x0, y0) result(flux)
    implicit none
    real(dp), intent(in) :: d, w, h, x0, y0
    
    flux = (d*w*h)/( (4._dp*pi)*( ((x0+(d/2._dp))**2) + ((y0-(w/2._dp))**2) + ((h/2._dp)**2)  ) )

end function large_x0_flux

!-----------------------------------------------------------------------
!! Subroutine: box_flux_monte_carlo
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine packs information into arrays a(:), b(:), and data(:).
!! These are sent to monte_carlo_quad in order to calculate flux. This
!! is done in this way in order to keep monte_carlo_quad versatile. 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine box_flux_monte_carlo(depth, width, height, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: flux, sigma_f
    
    real(dp) :: a(1:3), b(1:3), data(1:2)

    ! a is the lower integration limit in the x, y, z coordinates. 
    ! since the origin was placed at the corner of the nuclear reactor the
    ! lower limit is zero in all coordinates
    a = 0._dp

    ! b is the upper integration limit in the x, y, z coordinates.
    b(1) = depth
    b(2) = width
    b(3) = height

    ! This is the 'work array' and contains parameters
    ! (other than the sample point) needed to evaluate the function to
    ! integrate
    data(1) = x_zero
    data(2) = y_zero

    ! monte_carlo_quad is called here, sent flux_kernel_vector as monte_carlo_quad's "f",
    ! this calculates flux integration using flux_kernel's integrand.
    call monte_carlo_quad(flux_kernel_vector, a, b, data, n_samples, flux, sigma_f)
end subroutine box_flux_monte_carlo

!-----------------------------------------------------------------------
!! Function: flux_kernel_vector
!-----------------------------------------------------------------------
!! By: 
!!
!! Description: This function unpacks the data sent from monte_carlo_quad,
!! this data is then sent to flux_kernel in order to obtain the value
!! of the integrand (fx(i) in monte_carlo_quad).
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
! Because of the interface defined in the quadrature module the 
! Monte Carlo subroutine expects a kernel function that receives two
! arrays, the first one contains the sampling point, the second one
! contains the parameters needed to calculate the kernel. 
real(dp) function flux_kernel_vector(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0

    x = x_vector(1)
    y = x_vector(2)
    z = x_vector(3)
    x0 = data(1)
    y0 = data(2)

    ! We're going to use the function we already defined for the 
    ! Boole's integration.
    k = flux_kernel(x, y, z, x0, y0)
end function flux_kernel_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ADVANCED PART STARTS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For Boole's method we need to calculate
! the flux of the solid box and subtract the flux from a solid sphere.
! Let's start defining the function that calculates the flux of a 
! solid sphere

!-----------------------------------------------------------------------
!! Function: sphere_flux_booles
!-----------------------------------------------------------------------
!! By:
!!
!! Description: This function uses the booles quadrature method of 
!! approximating integrals in order to approximate the volume 
!! integral in spherical coordinates which gives the flux produced by the
!! hollow spherical region within the box if it was indeed a reactor of its own.
!! This is because we will subtract this theoetical value from the 
!! solid box in order to get the hollow box flux.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! distance     real        Distance from the center of the reactor to the detector
!! radius       real        Radius of the spherical reactor
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function sphere_flux_booles(distance, radius, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: distance, radius
    integer, intent(in) :: n_grid

    real(dp) :: delta_r, delta_theta
    real(dp), allocatable :: f_r(:), g_rtheta(:)
    integer :: n_bins, i_r, i_theta
    real(dp) :: r, theta

    ! Base the rest of the function on the one for the solid box.
    ! Here we're integrating only two variables (r and theta).
    ! r is integrated from 0 to radius while theta is integrated
    ! from 0 to pi 

    ! "distance" is the distance from the center of the sphere 
    ! to the position of the detector.
    ! It will be given as a input to this function and passed
    ! to the sphere_flux_kernel function defined below.
    ! As before, the number of intervals "n_bins" is defined
    ! here in order to obtai correct sizes of delta_r and 
    ! delta_theta. 
     n_bins = n_grid - 1
     print*, 'n_grid', n_grid
    ! Dividing delta_r and delta_theta into proper chunks,
    ! in order to index through to radius and pi respectively. 
     delta_r = radius/(float(n_bins))
     delta_theta = pi/(float(n_bins))
    ! Allocating arrays for the douple loop for integration over all
    ! space in spherical.
     allocate(g_rtheta(1:n_grid))
     allocate(f_r(1:n_grid))
! As before, setting variable radius of hollow spherical chamber to start 
! at zero requires having r start at -delta_r. 
r = -delta_r

    do i_r = 1, n_grid
    !Inside loop for radial direction, terminated if r
    !reaches the max distance over r, which is user input "r_max".
        if(r < radius) then
            r = r + delta_r
            print*, 'r=', r
        else
            exit
        endif
        ! theta set to negative delta_theta (step in theta) 
        ! in order for integration to start at theta=0.
             theta= -delta_theta

            do i_theta = 1, n_grid
             !Inside loop for theta direction, terminated if theta
             !reaches the max distance over theta, which is pi.   
            if (theta < pi) then
                theta = theta + delta_theta
                print*, 'theta=', theta
            else
                exit
            endif
                 
                ! Now the array of values of the evaluated integrand located in
                ! sphere_flux_kernel subroutine for various theta is filled.
                g_rtheta(i_theta) = sphere_flux_kernel(r, theta, distance)

            enddo
            ! Resulting array g_rtheta above is sent to booles quadrature
            ! along with delta_theta in order to approximate integral over theta.
             f_r(i_r) = booles_quadrature(g_rtheta, delta_theta)

        enddo
        ! Resulting array f_r from above is sent to booles quadrature along
        ! with delta_r in order to approximate integral over r.
        flux = booles_quadrature(f_r, delta_r)




end function sphere_flux_booles

!-----------------------------------------------------------------------
!! Function: sphere_flux_kernel
!-----------------------------------------------------------------------
!! By: 
!! Description: This is a funtion meant to compute the integrand for 
!! the three dimensional integral in spherical done in function sphere_flux_booles.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! r_prime      real        r coordinate of the small integration volume
!! theta        real        theta coordinate of the small integration volume
!! big_r        integer     distance from the center of the sphere to the detector
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function sphere_flux_kernel(r_prime, theta, big_r) result(k)
    implicit none
    real(dp), intent(in) :: r_prime, theta, big_r

    k =  (sin(theta)*r_prime*r_prime) / ( (2._dp)*( (r_prime**2._dp) + (big_r**2._dp) - 2._dp*r_prime*big_r*cos(theta) ) ) 
    
end function sphere_flux_kernel

!-----------------------------------------------------------------------
!! Function: total_flux_booles
!-----------------------------------------------------------------------
!! By:
!!
!! Description: This function calculates the distance form the detector to
!! the center of the box (which is the center of the hollow spherical region).
!! It then sends informtaion to box_flux_booles and sphere_flux_booles in
!! order to retrieve the flux of the solic box and the solid sphere.
!!
!! It then subtracts these two values in order to calculate the flux
!! of the box with a hollow spherical region centered inside. This is done
!! by subtracting the sphere's flux from the box's flux.
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_grid

    real(dp) distance, box_flux, sphere_flux

    ! Now that we have a function to calculate the flux of the solid box and
    ! another one for the solid sphere we just need to use both functions 
    ! and calculate the difference.

    ! distance is the distance between the position of the detector (x_zero, y_zero)
    ! and the center of the sphere (which is also the center of the box)
     distance = sqrt( ((x_zero + (depth/2._dp))**2) + ((y_zero - (width/2._dp) )**2) + ((height/2._dp)**2) )

     box_flux = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid) 
     !Flux calculated for solid box at given detector position
     print*, 'Box Flux Value', box_flux
     sphere_flux = sphere_flux_booles(distance,radius,n_grid)
     !Flux calculated from sphere 
     print*, 'Sphere Flux Value', sphere_flux

     flux = box_flux - sphere_flux
     !Flux of the box with a sphere cut out of middle (box minus sphere)

    

end function total_flux_booles

! The Monte Carlo approach is simpler.
! We just need to define a new kernel function that is zero if the 
! sampling point is inside the sphere and the original kernel if
! the sampling point is outside of the sphere.

! This new kernel will take to arrays, one with the coordinates
! of the sampling point and one with all the other needed parameters
! this time it's more than just the position of the detector.

!-----------------------------------------------------------------------
!! Function: hollow_box_flux_kernel
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Explaination: This subroutine unpacks the data and x_vector arrays,
!! calculates the distance from the random sampling point to the center of
!! of the sphere, and sends some of the unpacked information to the
!! flux kernel subroutine used in the solid box portion of the program. 
!! The flux kernel is sent data only if the distance from the random point in
!! the box to the center of the sphere is greate than the radius of the spherical
!! hollow chamber. This ensures flux values are only calculated for point within the solid
!! portion of the box (inside box but outside spherical chamber)
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the sphere's radius and x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function hollow_box_flux_kernel(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0, radius, x_start, y_start, z_start 
    real(dp) :: distance_to_center
    ! unpacking x_vector which houses random point in box.
     x = x_vector(1)
     y = x_vector(2)
     z = x_vector(3)
    ! unpacking detector location and radius of spherical chamber
     x0 = data(1)
     y0 = data(2)
     radius = data(3)
    ! unpacking and calculating coordinate of center of sphere
    ! and box in cartesian coordinates.
     x_start = data(4)/2._dp
     y_start = data(5)/2._dp
     z_start = data(6)/2._dp





    ! We need to determine whether or not the sampling point is inside the 
    ! sphere. For that ywe can calculate the distance from the  random sampling point
    ! (THIS IS NOT THE POSITION OF THE DETECTOR) and the center of the sphere 
    ! and compare it with the sphere's radius
    ! Origin located at corner of Box

     distance_to_center = sqrt( ((x-x_start)**2) + ((y-y_start)**2) + ((z-z_start)**2) )
    ! If statement in order to ensure flux calculation only occurs within box but outside sphere. 
     if (distance_to_center <= radius) then
    ! no k values if within hollow space
         k = 0._dp
     else
    ! Regular k value if within solid portion of box.
         k = flux_kernel(x, y, z, x0, y0)
         
     end if

end function hollow_box_flux_kernel


!-----------------------------------------------------------------------
!! Subroutine: hollow_box_flux_mc
!-----------------------------------------------------------------------
!! By: 
!!
!! Description: This subroutine packs up information to be sent to
!! monte_carlo_quad, since monte_carlo_quad needs some information
!! passed in arrays (called "work arrays") in order to keep
!! monte_carlo_quad versatile. 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) ::  flux, sigma_f

    real(dp) :: a(1:3), b(1:3), data(1:6), deltaAB(1:3)

    ! Declare "a" array empty for proper scaling in monte_carlo_quad
     a = 0._dp
    ! Fill "b" array with dimention information.
     b(1) = depth
     b(2) = width
     b(3) = height
    ! Declare array deltaAB to have values (b-a)
     deltaAB = b - a

    ! Filling data array with detector position, radius of spherical chamber,
    ! and the depth, width, and height of the box (deltaAB(1), (2), and (3)).

     data(1) = x_zero
     data(2) = y_zero
     data(3) = radius
     data(4) = deltaAB(1)
     data(5) = deltaAB(2)
     data(6) = deltaAB(3)
    ! monte_carlo_quad called here, given function hollow_box_flux_kernel as its "f"
    ! in order to ensure calculation of fx(i) happens only in solid portion of box,
    ! outside of hollow region. 
    
     call monte_carlo_quad(hollow_box_flux_kernel, a, b, data, n_samples, flux, sigma_f )
end subroutine hollow_box_flux_mc

end module neutron_flux