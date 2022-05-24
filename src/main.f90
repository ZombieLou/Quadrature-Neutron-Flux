! Program: nuclear_reactor
! By: Louis Andre
!-----------------------------------------------------------------------------
! This program calculates the flux for a reactor of variable size (box dimensions), 
! via a detector placed a varialbe distance from one of the faces of the box reactor. 

! The second section of the program calculates the flux for the same box except with 
! a spherical hollow chamber in the center of the box (varying the radius of the region). 
! With variable hollow chamber, the flux is calculated a fixed distance away from the 
! reactor.


! You're free to give different values when running the code, the ones 
! bellow are just a suggestion.

! Basic part of the project
! depth = 40.0_dp
! width = 100.0_dp
! height = 60.0_dp
! y_zero = 30._dp
! x_min = 5_dp
! x_max = 200._dp
! x_step = 5._dp
! n_grid = 25
! m_samples = 10000

! Advanced part of the project
! x_zero = 80._dp
! r_min = 1.0_dp
! r_max = 19._dp
! r_step = 1.0_dp
!-----------------------------------------------------------------------------
program nuclear_reactor
use types 

use read_write, only : read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

implicit none

real(dp) :: depth, width, height, y_zero, x_min, x_max, x_step
integer :: n_grid, m_samples
real(dp) :: x_zero, r_min, r_max, r_step


! Basic part of the project, not including the reacor with the spherical hollow region within the volume. 
call read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
call write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)

! Advanced part of the program, allowing for the calculaton of the flux with a varable size
! spherical hollow chamber within the reactor. 
call read_advanced_input(x_zero, r_min, r_max, r_step)
call write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)

end program nuclear_reactor