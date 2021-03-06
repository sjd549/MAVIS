## MAVIS
# Mega Analysis & VIsualisation Suite 

Data processing and analysis software for MEGA output files.

Files recognised: 

MEGA Inputs: 	sim128shot#.txt, gshot#, seq.in
MEGA Outputs: 	energy-n.txt, energy-phys.txt, seq.harmonics, seq.moments, FI_markers


MAVIS requires the following software/modules to function:

	python-pip
	python-numpy
	python-scipy
	python-matplotlib
	ffmpeg
	tqdm

Users may install these manually or the required libraries 
can be installed automatically via the -f or --first flags:

./MAVIS -f
./MAVIS --first

MAVIS will require sudo user privilages to perform this action.








MAVIS is designed for use in it's own seperate folder. When executed it will search for directories and sub-directories which contain HPEM output files and catagorize these into seprate 'simulations' based on which directory they came from. 
Each simulation will then be processed in turn with the output from requested diagnostics being saved in seperate directories within the simulation directory. Diagnostics which perform comparisons between simulations will be saved in the root MAVIS directory.

E.g. to compare a set of three simulations:

	1) Copy the three directories containing the output files into your MAVIS directory.

	2) Edit the MAVIS switchboard, requesting your diagnostics, parameters and image settings.

	3) Run MAVIS and wait for diagnostics to complete.

	4) Obtain processed images and data from newly created sub-directories.









