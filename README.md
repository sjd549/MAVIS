# JULIA
Hpem ELectronic ENgine Analysis: data plotting and analysis software for TECPLOT output files.

Files recognised: 

TECPLOT2D.PDT, kin.pdt, movie1.pdt, movie_icp.pdt


Requires the following software/modules to function:

python-pip, python-numpy, python-matplotlib, ffmpeg, findtools.




HELENA is designed for use in it's own seperate folder. When executed it will search for directories and sub-directories which contain HPEM output files and catagorize these into seprate 'simulations' based on which directory they came from. 
Each simulation will then be processed in turn with the output from requested diagnostics being saved in seperate directories within the simulation directory. Diagnostics which perform comparisons between simulations will be saved in the upper HELENA directory.

E.g. to compare a set of 5 simulations:

Copy the 5 directories containing the output files into your HELENA directory.
Edit the HELENA switchboard, requesting your diagnostics, parameters and image settings.
Run HELENA and allow the diagnostics to be completed.
Obtain your processed images and graphs from within the simulation directories.









