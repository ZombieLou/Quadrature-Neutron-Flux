# What this Program Does

   The program has two major parts: First, the user defines a box shaped nuclear reactor with a detector placed a variable distance away from one of the faces. The distance of the detector from the reactor is defined by the user (user defines minimum and maximum distance of the detector as well as the distance of each step). The program calculates the flux at the detector location, which starts at the user defined minimum location and moves one step at a time to the maximum location, while calculating flux at each step. 
    The second part of the program calculates the flux from the same box but with a hollow spherical region within the box, centered on the center of the box. The radius of the spherical region starts at a minimum value and steps to a maximum value taking flux measurements at each step. The detector is stationary at a set distance from the reactor while the hollow region radius increases. These minimum, maximum, and step radius values are provided by the user, as well as the detector distance. 

# How to Use

   To use the program, make sure all `.f90` files are in the same directory and type "make" into the terminal, this will compile the files and create an executable called `nuclear-reactor`. In the Terminal type "./nuclear-reactor" to run the program. You will be  prompted to enter values for the height, width, and depth of the box, as well as: the y-coordinate of the detector, the minimum and maximum distance of the detector from the reactor and the step size, the number of intervauls over which the integrations will be approximated (sums over finite number of posiions along limits of integration), and the number of sampling points for the monte carlo method. It will tell you the fluxes (for each detector distance using each integration method) were written into a file called `results_basic.dat`. Open this file if you wish to see the raw data of the calculations. Otherwise you can visualize them by running the jupyter notebook `.ipynb` file.
    The program will then request another set of inputs from the user for the hollow reactor section of the program. You will be prompted to enter the position of the detector (stationary for this part of the program), the minumum and maximum radius of the spherical hollow region, and the step size of this radius. The fluxes will be measured for the box with a hollow region growing from the minimum radius to the maximum radius by the step size of the radius given by the user. It will tell you the fluxes were written into a file called `results_advanced.dat`. Open this file if you wish to see the raw data of the calculations. Otherwise you can visualize them by running the jupyter notebook `.ipynb` file.
    
# Contents
   
   This program contains the files `neutron_flux.f90`, `read_write.f90`, `types.f90`, `main.f90`, `quadrature.f90`, as well as `plots.ipynb`, this `readme.md`, and the `makefile`. 
`neutron_flux.f90` calculates the fluxes for the first and second parts of the program
`quadrature.f90` uses mathematical methods of quadrature to evaluate integration needed for flux values.
`read_write.f90` takes inputs from user and prints data to results files, also initializing calculations and stepping through detector distance and radius.
`types.f90` contains types of inputs to be used (real double precision, and pi ~ 3.14159...)
`main.f90` calls to read/write module to execute the program. 
`plots.ipynb` can be opened with jupyter notebook, contains graphs of flux over distances and radius measures. Hit "Shift + Enter" on keyboard repeatedly to initialize all graphs when notebook is open. 
`makefile` allows user to type "make" into terminal to compile all code
`readme.md` is this file, containing instructions for the use of this program. 
   `
