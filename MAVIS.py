#!/usr/bin/env python

#######################################
#		   Point of Contact		      #
#								      #
#	      Dr. Scott J. Doyle	      #
#	     University of Seville	      #
# Dept. Atomic and Molecular Physics  #
#	  Calle Tomas Alva Edison, 7	  #
#	   Seville, Andalusia, Spain      #
#	   Scott.Doyle@Physics.org		  #
#									  #
#######################################
#               'MAVIS'               #
# Mega Analysis & VIsualisation Suite #
#######################################



#====================================================================#
				 #PROGRAM FLAGS AND MODULE IMPORTS#
#====================================================================#

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--first", action="store_true", dest="install", default=False, help="Install prompt for required python modules")
(options, args) = parser.parse_args()

if 'True' in str(options): 
	import os, sys
	import os.path

	print ''
	print 'First time use requires installation of additional python modules'
	print 'Please type your password when prompted to allow installation:'
	print ''
	try:
		os.system('sudo apt-get install python-pip')
		os.system('sudo apt-get install python-matplotlib')
		os.system('sudo apt-get install python-numpy')
		os.system('sudo apt-get install python-scipy')
		os.system('sudo apt-get install ffmpeg')
		os.system('pip install tqdm')
	except:
		print ''
		print 'Error installing required packages'
		print 'Please attempt manual installation'
		print ''
	#endtry
	print ''
	print ''
#endif

#==============#

#Import core modules
import matplotlib.cm as cm
import matplotlib
import numpy as np
import scipy as sp
import math as m
import subprocess
import warnings
import os, sys
import os.path
import glob

#Import additional modules
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.io import FortranFile as ff
from scipy.signal import savgol_filter
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs
from matplotlib import ticker
from scipy import ndimage
from tqdm import tqdm
from pylab import *
 

#====================================================================#
				  		 #LOW LEVEL INPUTS#
#====================================================================#

#Various debug and streamlining options.
DebugMode = False						#Produces debug outputs for relevent diagnostics.

#Warning suppressions
warnings.simplefilter(action='ignore', category=FutureWarning)
#Fix any future warnings, mostly related to Savgol_filter treatment of multidimentional arrays
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images
#Fix "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

#List of recognized data extensions for file readin
FileExtensions = ['.hamonics','moments','.txt','.in','.nam','.dat','.out']

#Default poloidal mesh repository location (local or root)
DirRepository = os.path.abspath(".")+'/Repository'

#Numerical Calculation Methods:
GlobSheathMethod = 'AbsDensity'			#Set Global Sheath Calculation Method.
#Choices: ('AbsDensity','IntDensity')
GlobThrustMethod = 'AxialMomentum'		#Set Global Thrust Calculation Method. 
#Choices:('ThermalVelocity','AxialMomentum')
GlobMeanCalculation = 'MeanFraction'	#Definition of 'mean' EDF value
#Choices: ('MeanEnergy','MeanFraction')

#Data Filtering and Smoothing Methods:
KineticFiltering = True						#Pre-fit kinetic data employing a SavGol filter
PlotKineticFiltering = False				#Plot Filtered Profiles, or employ only in trends.
Glob_SavWindow, Glob_SavPolyOrder = 25, 3	#Window > FeatureSize, Polyorder ~= Smoothness


####################

#Commonly used variable sets.
Ctrl = ['kst','t']
Axes = ['r_psi','gpsi_nrm','q_psi']
Phys = ['vrad','vtheta','vphi','brad','btheta','bphi','erad','etheta','ephi','prs','rho','dns_a','mom_a', 'ppara_a','pperp_a','qpara_a','qperp_a']
Kin = ['R_gc','Z_gc','Phi_gc','p_gc','pphi_gc','mu_gc','E_gc','Lambda_gc','psip','ffff']

#Archived variable sets
#Phys = []

####################

#Common Diagnostic Settings:

#===== AUG#34570 =====#
#radialprofiles = [90,180]
#poloidalprofiles = [0.20,0.40,0.65]
#toroidalprofiles = []
#trendlocation = []

#setting_SEQ = [0,1]
#setting_ntor = [0,-2]
#setting_kstep = [00,40,5]

####################




#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested Variables and Plotting Locations:
variables = ['prs','brad','vrad']		#Requested variables to plot
#['dns_a','mom_a', 'ppara_a','pperp_a']	#Phys

radialprofiles = []						#1D Radial Profiles (fixed theta, phi) :: Poloidal Angle [deg]
poloidalprofiles = [0.20,0.40,0.65]			#1D Poloidal Profiles (fixed rho_pol, phi) :: Norm. Radius [-]
toroidalprofiles = []					#1D Toroidal Profiles (fixed rho_pol, theta) :: Toroidal Angle [deg]
trendlocation = [] 						#Cell location For Trend Analysis [R,theta,phi], ([] = min/max)

#Various Diagnostic Settings:
setting_SEQ = [0,14]					#Simulation SEQ to load		- [Min,Max], [Int], [] to plot SEQ001
setting_ntor = [0,-2]					#ntor range to plot 		- [Min,Max], [Int], [] to plot all
setting_kstep = [00,40,1]				#kstep index range to plot 	- [Min,Max,Step], [Int], [] to plot all


#Requested diagnostics and plotting routines:
savefig_1Denergy = False				#Plot 1D MHD energies (1 Sim) 		(xxx.energy_p)	- Working
savefig_1Denergytrends = False			#Plot 1D MHD energies (multi-Sim) 	(xxx.energy_n)	- Working

savefig_1Dequilibrium = False			#Plot 1D radial/poloidal profiles	(xxx.harmonics) - Working	-ASCII
savefig_2Dequilibrium = True			#Plot 2D poloidal x-sections		(xxx.harmonics)	- Working
savefig_2Dequilmovie = True			#Plot 2D poloidal x-section movies	(xxx.harmonics)	- Working	-ASCII

savefig_2Dcontinuum = False				#Plot 2D harmonic continuum		 	(xxx.harmonics)	- Working
savefig_2Dpolspectrum = False			#Plot 2D poloidal spectra 			(xxx.harmonics)	- Working	-ASCII
SpectralVariable = 'brad'; QuickPROES = False
ContinuumVariable = 'vrad'

savefig_1Dkinetics = False				#Plot 1D kinetic distributions	 	(gc_a_kstepxxx)	- Working	!NEED FUNCS
savefig_2Dkinetics = False				#Plot 2D kinetic distributions	 	(gc_a_kstepxxx)	- Working	!NEED FUNCS


#Requested diagnostic terminal outputs:
print_generaltrends = False				#Verbose Trend Outputs								- In Development


#Write processed data to ASCII files:
write_ASCII = True						#All diagnostic outputs written to ASCII.dat		- In Development
write_ASCIIFormat = 'RSV'				#Choose ASCII file output format ('RSV', 'CSV')		- In Development


#Image plotting options:
image_extension = '.png'				#Extensions ('.png', '.jpg', '.eps')				- Working
image_aspectratio = [10,10]				#[x,y] in cm 										- Working
image_rhocrop = []						#Crop image radius/rho_pol [min,max]				- In Development
image_thetacrop = []					#Crop image poloidal angle [min,max]g				- In Development
image_mpolcrop = [-16,16]				#Crop image poloidal modes [min,max]				- Working (ish)
image_ntorcrop = []						#Crop image toroidal modes [min,max]				- In Development
image_cbarlimit = []					#[min,max] colourbar limits							- Working

image_cmap = 'plasma'					#Toggle global colourmap: 'plasma','IDL_Gamma_II' 	- Working (ish)
image_plotsymmetry = True				#Toggle radial symmetry in images					- In Development
image_contourplot = True				#Toggle contour Lines in images						- Working
image_plotgrid = False					#Plot major/minor gridlines in all figures			- Working
image_plotmesh = False					#Plot material mesh outlines in images				- In Development

image_normalise = False					#Normalise image/profiles to local maximum			- In Development
image_logaxes = [False,False]			#Apply logarithmic axes to image/profiles [x,y]		- In Development

#Overrides the automatic image labelling:
titleoverride = []
legendoverride = []
xaxisoverride = []
xlabeloverride = []
ylabeloverride = []
cbaroverride = []

#=====================================================================#
#=====================================================================#



#============================#
#        ####TODO####        #
#============================#
#
#FINISH THE DENORMALISATION ROUTINES IN ExtractMEGA_Harmonics() SUBROUTINE
#	Need to know how to denormalise each one, also needs a denormalisation toggle in the switchboard
#	Tie this option to any other applicable function, e.g. variablelabelmaker, etc...
#	Needs testing for all diagnostics and any special cases identified
#
#ADD A POLOIDAL HARMONIC 1D PLOTTER (see Liu2010, 2016 papers for examples)
#	This will need to be able to select a poloidal harmonic range and compare between sim folders
#	This will be a critical diagnostic tool for isolating MHD plasma response from Kinetic plasma response
#
#FINISH ExtractMEGA_moments() SUBROUTINE AND ADD MOMENTS DIAGNOSTICS
#	Ideally each existing diagnostic would be able to use either moments OR harmonics
#	If the above is impossible, create copies of existing diagnostics using same flags
#
#NEED TO ADD ENERGY PROFILE SUBPLOT TO THE EQUILMOVIE FIGURES
#	This needs to be an optional extra for both the spectral and equilibrium movie figures
#
#NEED TO ADD VARIABLE LOOP TO THE savefig_polspectrum DIAGNOSTIC
#	Loop over all regular variables by default
#
#NEED TO MAKE savefig_continuum AND savefig_equilibrium BOTH USE THE setting_ntor RANGE WHEN PLOTTING
#	Same applies to all the other switchboard ranges which are currently fudged
#
#ADD A SEQ.nam READIN FUNCTION USING LAMBDA FUNCTIONS
#	Should read each line and assign main variables (removing any trailing ! or !!! comments)
#
#ADD DOTTED LINE TO 1D (2D?) EQUILIBRIUM IMAGES TO SHOW THE VACUUM MASK EXTENT
#	Will require a SEQ.in (.nam) readin function to work, would be quite useful to have anyway!
#
#ADD ONE COPY OF 2D AXES TO ASCII FOLDERS BY DEFAULT AND HOMOGONIZE THE HEADER STRUCTURE
#	Need to save if variables are normalised or not in header, also maybe 'CSV', 'RSV'?
#
#ADD OPTION TO PLOT INPUT 2D FLUX SURFACE FUNCTION AND RADIAL FITS
#  	savefig_gfilefits --> output all figures into the pre-existing 'equil' folder
#
#
#
#
#ALTER setting_kstep TO USE THE ACTUAL KSTEP VALUES AND MAKE A FUNCTION TO TRANSLATE INTO SEQ AND KSTEP INDEX RANGES
#	Require making an icp.nam readin function and extracting all of the write steps and other inputs
#	Require making an additional "kstep_translation" function where the input setting_kstep is 
#	translated into an output set of kstep indices and associated SEQ indices 
#
#FIX KSTEPMOD CALCULATION - CURRENT CALCULATION ONLY WORKS IF KSTEP RANGE IS THE SAME FOR ALL SEQs
#		KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
#
#ERROR WARNINGS FOR SEQ_RANGE AND KSTEP_RANGE IN THE INPUT DECK
#
#SPEED UP READ-IN OF LARGE DATA SETS AS MUCH AS POSSIBLE TO ENABLE FASTER RESPONSE CALCULATIONS
#	Remove "While" in ReadMEGA_Harmonics() function, it currently repeats every previous KStep
#	i.e. the current implimentation has to do KStep! (additive factorial) iterations.
#	Stop reading in the Crdr and Crdz axes on each loop, this is wasteful
#		Note for above, read these axes from the psi.dat file OR from seperatrix file (update icp.com to copy)
#
#ADD SEPERATE OUTPUT FOLDERS FOR EACH VARIABLE FOR 1D EQUIL PROFILES - SAME AS 2D EQUIL PROFILES
#	In general, I guess it's better to have folders of variables rather than folders of ntor
#	There are more variables and fewer ntor, so each folder will have less 'clutter' that way.
#
#ADD DOTTED LINE OR SHADED AREA TO INDICATE SEQ NUMBERS IN THE ENERGY DIAGRAMS (AND JUST GENERALLY)
#	Use KStepMod as the KStepArray indices to apply the line, make an image_SEQline input for switchboard
#
#FIX ISSUE WHERE "outputdata is referenced before assignment" IF FILENAME HAS [] or {} IN IT
#	Issue arises in ExtractMEGA_Energy (and others) where glob.glob() can't handle brackets in file directory 
#
#ADD rhocrop, thetacrop, mpolcrop, ntorcrop OPTIONS TO ALL APPLICABLE DIAGNOSTICS
#	Will require a cropping function which considers any image rotation and corrects for zero crop input '[]'
#
#ADD INFORMATION REGARDING SIMULATION RESOLUTION SEQ, KSTEP, ETC... BELOW MAVIS SPLASH
#
#ADD OPTION TO HOMOGONISE THE COLOURBAR FOR ALL MOVIES (EQUIL / RESPONSE / KINETIC)
#
#ADD GEQDSK READ AND WRITE FUNCTIONS (WITH OPTIONS FOR MEGA / FIESTA FORMATS)
#	Create a new diagnostic (savefig_simsetup) which plots gfile inputs for analysis
#
#Extract lpsi and mpol from data without having to explicitly hard-code it (see Read_Harmonics function)
#
#ADD OPTION TO PLOT DIFFERENT TOROIDAL ANGLES
#


#=====================================================================#
#=====================================================================#

#USEFUL INFORMATION:
#
#lpsi = radial spatial resolution					#[cells]
#ltheta = poloidal angle spatial resolution			#[cells]
#lphi = toroidal angle spatial resolution			#[cells]
#
#mpol = number of poloidal harmonics considered 	#[int]	- low-pass filter limited?
#ntor = number of toroidal harmonics considered 	#[int]	- low-pass filter limited?
#
#PoloidalAngle = 2*np.pi*(float(ltheta)/float(ltheta_total)) 	#[Rad]
#ToroidalAngle = 2*np.pi*(float(lphi)/float(lphi_total))		#[Rad]
#






















#====================================================================#
						#COMMON I/O FUNCTIONS#
#====================================================================#

def DirectoryContents(AbsPath):
#Takes directory path and returns all contents and sub-folders
#Inputs: one string; absolute path to desired folder
#Returns: two lists; one of all sub-folders and one of all non-folder files
#Example: HomeDirContents,HomeDirFolders = DirectoryContents(os.path.abspath("."))

	#Obtain contents of supplied directory and initiate folders list
	DirContents = os.listdir( AbsPath )		#List of all contents (inc non-folders).
	DirFolders = list() 					#List of all folders.

	#Determine any sub-folders within supplied directory
	for i in range(0,len(DirContents)):
		if os.path.isdir(AbsPath+'/'+DirContents[i]) == True:
			#Append slash to returned folder directories to differentiate them
			DirFolders.append('/'+DirContents[i]+'/')
		#endif
	#endfor

	return(DirFolders,DirContents)
#enddef

#=========================#
#=========================#

def CreateNewFolder(Dir,DirString):
#Creates a new folder if one does not already exist.
#Takes destination dir and namestring, returns new directory.

	#Try loop avoids soft-crash if folder already exists
	try:
		NewFolderDir = Dir+DirString+'/'
		os.mkdir(NewFolderDir, 0755);
	except:
		Folder_Already_Exists = 1
	#endtry

	return(NewFolderDir)
#enddef

#=========================#
#=========================#

def FolderNameTrimmer(DirString,Index=1):
#Takes folder names and returns item after requested underscore index.
#Note, index > 1 will return between two underscores, not the entire string.

	#Define delimeter specifying trim location
	delimiter = '_'

	#Try avoids error if folder name doesn't contain delimiter
	try:
		for i in range(0,Index):
			SplitLoc = str(DirString[::-1]).index(delimiter)
			SplitIdx = len(DirString) - SplitLoc
			NameString = DirString[SplitIdx:-1]
			DirString = DirString[:SplitIdx-1]
		#endfor
	except:
		NameString = str(DirString[2:-1])
	#endtry

	return(NameString)
#enddef

#=========================#
#=========================#

def ReadCmap_ASCII(CmapDir,CmapName):
#Takes filename directory to ASCII colourmap (Cmap.txt)
#Returns matplotlib Cmap object for use in plotting:
#	e.g. plt.imshow(Image, aspect='auto', cmap=Cmap_Blues)
#Example: Cmap_Blues = ReadCmap_ASCII(os.getcwd()+'/Blues.txt','Blues')

	#Import custom colourmaps module
	from matplotlib import colors as c

	#Load colourmap from supplied directory
	map = np.loadtxt(CmapDir, delimiter=',')
	Cmap = c.ListedColormap(map.T, name=CmapName)

	return(Cmap)
#enddef

#=========================#
#=========================#

def ReadRawData(Dirlist,NameString,ListIndex):
#Takes directory list and data filename type (e.g. .png, .txt)
#Returns datalist of contents and length of datalist.
#Example: rawdata, datalength = ReadRawData(Dir,'.dat',l)

	#Attempt to open file and extract data
	try:
		DataFileDir = filter(lambda x: NameString in x, Dirlist)
		Rawdata = open(DataFileDir[ListIndex]).readlines()
		nn_data = len(Rawdata)
	except:
		print 'Unable to extract '+str(NameString)
		exit()
	#endtry

	return(Rawdata,nn_data)
#enddef

#=========================#
#=========================#

def WriteFile_ASCII(data,filename,structure='w',Orientation='CSV'):
#Takes a 1D or 2D array and writes to a datafile in ASCII format.
#Imputs: Data, Filename, 'w'rite or 'a'ppend, and orientation (CSV or RSV).
#Example: WriteFile_ASCII(Image, "Filename", 'w', 'CSV')

	#Determine dimensionality of profile.
	#If first dimension is a 1D list ==> 2D array
	if isinstance(data[0], (list, np.ndarray) ) == True:
		#Open new textfile and output 2D image data.
		datafile = open(filename, structure)
		for m in range(0,len(data)):
			for n in range(0,len(data[m])):
				datafile.write(str(data[m][n]))
				datafile.write(' ')
			#endfor
			datafile.write('\n')
		#endfor
		datafile.close()

	#If lowest dimention is scalar: ==> 1D array.
	elif isinstance(data, (list, np.ndarray) ) == True:
		#Open new textfile and output 2D image data.
		datafile = open(filename, structure)
		for n in range(0,len(data)):
			datafile.write(str(data[n]))
			datafile.write(' ')
		#endfor
		datafile.close()

	return()
#enddef

#=========================#
#=========================#

def ReadFile_ASCII(Filename,HeaderIdx=0,Dimension='2D',Orientation='CSV'):
#Reads 1D or 2D data from textfile in ASCII format, returns data and header.
#Input filename, header length, data dimension and orientation (CSV or RSV).
#Example: OutputData,Header = ReadFile_ASCII('/Data.txt', 0, '2D', 'CSV')

	#Define any required lists
	OutputData,Header = list(),list()

	#If data is saved 'Row-wise', use default readin routine.
	if Orientation == 'RSV':
		#Determine dimensionality of profile.
		if Dimension in ['1D','2D']:
			#Read in 2D data from ASCII formatted file.
			datafile = open(Filename)
			RawData = datafile.readlines()

			#Extract header and raw data
			for m in range(0,HeaderIdx): Header.append(RawData[m])
			RawData = RawData[HeaderIdx::]

			#Read each row, split it (space delimited) and save.
			for m in range(HeaderIdx,len(RawData)):
				Row = RawData[m].split()
				for n in range(0,len(Row)):
					try: Row[n] = float(Row[n])
					except: Row[n] = str(Row[n])
				#endfor
				OutputData.append(Row)
			#endfor
		#endif

	#=====#

	#If data is saved 'column-wise', transpose the arrays to correct.
	elif Orientation == 'CSV':
		#Determine dimensionality of profile.
		if Dimension in ['1D','2D']:
			#Read in 2D data from ASCII formatted file.
			datafile = open(Filename)
			RawData = datafile.readlines()

			#Extract header and raw data
			for m in range(0,HeaderIdx): Header.append(RawData[m])
			RawData = RawData[HeaderIdx::]


			#AD-HOC FIX FOR EMPTY MARKER FILES - REMOVE ONCE write_ep() SAVES HEADER
			if len(RawData) in [0,1]:
				return (np.zeros([10,1]).tolist(), np.zeros([10,1]).tolist())
			#endif


			#Enlarge output data array by number of columns
			NumColumns = len(RawData[HeaderIdx+1].split())
			for m in range(0,NumColumns):
				OutputData.append(list())
			#endfor

			#Read each row, split it and save into relevant column of output data.
			for i in range(HeaderIdx,len(RawData)):
				Row = RawData[i].split()
				for j in range(0,len(Row)):
					try: Row[j] = float(Row[j])
					except: Row[j] = str(Row[j])
				#endfor
				for k in range(0,NumColumns):
					OutputData[k].append(Row[k])
				#endfor
			#endfor
		#endif
	#endif

	#=====#

	#Orientation doesn't matter if 0D (scalar data).
	elif Dimension == '0D':
		#Read in 0D data from ASCII formatted file.
		datafile = open(Filename)

		for m in range(0,HeaderIdx): Header.append(RawData[m])
		RawData = datafile.readlines()[HeaderIdx::]
		Row = RawData.split()

		for m in range(0,len(Row)):
			OutputData.append(float(Row[m]))
		#endfor
	#endif

	return(OutputData,Header)
#enddef

#=========================#
#=========================#

def Read_MEGAHarmonics(Filename,Variable,mpol,ntor,lpsi,kstep=np.nan):
#Reads MEGA xxx.harmonics FORTRAN binary output file and extracts data into a 3D or 4D object.
#Data is read for a single variable [variable] over all timesteps [KStep], for a single SEQ
#Inputs: SEQ.harmonics filepath [no Root], Variable name string [Variable] as defined in DataFormat,
#Inputs: Radial mesh resolution [lpsi], poloidal mesh resolution [mpol], 
#Inputs: Total number of toroidal harmonics [ntor] including positive, negative and n=0
#Returns: Data object for requested variable with structure: Data[kstep][mpol][ntor][lpsi][A,B]
#Example: HarmonicsData = Read_MEGAHarmonics('FolderName/data/001.harmonics','bphi',64,5,201]

	#Compute flattened 3D data array length based upon mesh resolution
	n_elem = (mpol+1)*ntor*lpsi*2

	#Define FORTRANFile save data format
	#KStep (kst) is an undefined length 1D integer array		[-]
	#t (SI Time) is an undefined length 1D float array 			[IonGyroFreq*ms]
	#r_psi, gpsi_nrm, q_psi are [lpsi] length 1D float arrays 	[various]
	#All other variables are [n_elem] length 1D float arrays 	[various]
	DataFormat = np.dtype([\
	('kst',np.int32,1),\
	('t',np.float64,1),\
	('r_psi',np.float64,lpsi),\
	('gpsi_nrm',np.float64,lpsi),\
	('q_psi',np.float64,lpsi),\
	('vrad',np.float64,n_elem),\
	('vtheta',np.float64,n_elem),\
	('vphi',np.float64,n_elem),\
	('brad',np.float64,n_elem),\
	('btheta',np.float64,n_elem),\
	('bphi',np.float64,n_elem),\
	('erad',np.float64,n_elem),\
	('etheta',np.float64,n_elem),\
	('ephi',np.float64,n_elem),\
	('prs',np.float64,n_elem),\
	('rho',np.float64,n_elem),\
	('dns_a',np.float64,n_elem),\
	('mom_a',np.float64,n_elem),\
	('ppara_a',np.float64,n_elem),\
	('pperp_a',np.float64,n_elem),\
	('qpara_a',np.float64,n_elem),\
	('qperp_a',np.float64,n_elem),\
	])

	#Initiate data object to store SEQ.harmonic data
	if np.isnan(kstep) == True:
		#Initiate output data object and set appropriate internal structures
		Data = lambda:0
		Data.kst = np.array(([]),int)									#1D KStep Array 	[-]
		Data.time = np.array(([]),float)								#1D Time Array 		[IonGyroFreq*ms]
		Data.data = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)		#1D-3D Data Arrays	[various]
	elif np.isnan(kstep) == False:
		Data = lambda:0
		Data.kst = np.array(([]),int)									#1D KStep Array 	[-]
		Data.time = np.array(([]),float)								#1D Time Array 		[IonGyroFreq*ms]
		Data.vrad    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Radial   Velocity 	Array	[-]
		Data.vtheta  = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Poloidal Velocity 	Array	[-]
		Data.vphi    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Toroidal Velocity 	Array	[-]
		Data.brad    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Radial   B-Field  	Array	[-]
		Data.btheta  = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Poloidal B-Field  	Array	[-]
		Data.bphi    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Toroidal B-Field  	Array	[-]
		Data.erad    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Radial   E-Field  	Array	[-]
		Data.etheta  = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Poloidal E-Field  	Array	[-]
		Data.ephi    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Toroidal E-Field  	Array	[-]
		Data.prs     = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D MHD Pressure	  	Array	[-]
		Data.rho     = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D MHD density		  	Array	[-]
		Data.dns_a   = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D KIN Marker Density	Array	[-]
		Data.mom_a   = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D KIN Marker Momentum	Array	[-]
		Data.ppara_a = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D KIN Para Momentum  	Array	[-]
		Data.pperp_a = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D KIN Perp Momentum  	Array	[-]
		Data.qpara_a = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D KIN Para Safety		Array	[-]
		Data.qperp_a = np.empty(([0,mpol+1,ntor+lpsi,2]),np.float64)	#3D KIN Perp Safefy 	Array	[-]
	#endif

	#Open SEQ.harmonics file and ensure data exists
	try:
		FORTRANFile = ff(Filename,'r')							#Open 001.Harmonics FORTRAN format file
		RawData = FORTRANFile.read_record(DataFormat)			#Read RawData from file in Format
	except:
		print('\n \n  Data file "'+Filename+'" not found or formatted incorrectly - check FORTRAN dtype. \n')
		FORTRANFile.close(); exit()
	#endtry

	#=====#=====# 	#=====#=====# 	#=====#=====#

	#IF kstep == -1 is supplied, ONLY READ TEMPORAL AND SPATIAL AXES into a 1D HarmonicsData object
	#Read 1D rho_pol and q_psi on first iteration, and save 1D kstep and time arrays for plotting
	if kstep == -1:													#elif: Dimension == '1D'
		index = 0													#Initiate search index to zero
		while(True):
			#FORTRANFile.read_record automatically steps through all KSteps in the supplied SEQ.harmonics file.
			#RawData for KStep[i] is of shape RawData[Variable][datapoint] where all variables are flattened to 1D
			try: RawData = FORTRANFile.read_record(DataFormat)		#Read KStep=0, if open then read KStep=1, etc...
			except: FORTRANFile.close(); break						#If no KSteps remain, close FORTRAN file
			#endtry

			#Only extract rho_pol and q_psi on first Kstep
			if index == 0:
				Data.rho_pol = np.sqrt(abs(RawData['gpsi_nrm'][0]))
				Data.q_psi   = RawData['q_psi'][0]
			#Always extract 1D Kstep and Time arrays
			Data.kst     = np.append(Data.kst,  RawData['kst'][0])
			Data.time    = np.append(Data.time, RawData['t'][0])#*1e3/wa
			#Print kstep for debug purposes if requested
			if DebugMode == True: print(str(index)+'-'+str(Data.kst[index]))
			#endif

			#Delete RawData for current KStep and increment KStep counter
	  		del RawData; index += 1
		#endwhile

	#=====#=====#

	#IF Kstep > 0 is supplied, read ALL VARIABLES FOR SINGLE KSTEP into a 3D HarmonicsData object
	#Read data for each Kstep from the supplied SEQ.harmonics output folder
	#Save data for all variables once the requested Kstep is reached
	elif np.isnan(kstep) == False:									#elif: Dimension == '3D'
		index = 0													#Initiate search index to zero
		while(index <= kstep):
			#FORTRANFile.read_record automatically steps through all KSteps in the supplied SEQ.harmonics file.
			#RawData for KStep[i] is of shape RawData[Variable][datapoint] where all variables are flattened to 1D
			try: RawData = FORTRANFile.read_record(DataFormat)		#Read KStep=0, if open then read KStep=1, etc...
			except: FORTRANFile.close(); break						#If no KSteps remain, close FORTRAN file
			#endtry

			#Only extract rho_pol and q_psi on first Kstep
			if index == 0:
				Data.rho_pol = np.sqrt(abs(RawData['gpsi_nrm'][0]))
				Data.q_psi   = RawData['q_psi'][0]
			#Always extract 1D Kstep and Time arrays, until the requested Kstep - Shape = Data.kst[Real]
			if index >= 0:
				Data.kst     = np.append(Data.kst,  RawData['kst'][0])
				Data.time    = np.append(Data.time, RawData['t'][0])		#Normalised to ion gyro freq (1e3/Omega_i)
				#Print kstep for debug purposes if requested
				if DebugMode == True: print(str(index)+'-'+str(Data.kst[index]))
			#If index matches requested kstep, retrieve data for all variables and add to object
			if index == kstep:
				Data.vrad    = np.reshape(RawData['vrad'   ],(mpol+1,ntor,lpsi,2),order='F')
				Data.vtheta  = np.reshape(RawData['vtheta' ],(mpol+1,ntor,lpsi,2),order='F')
				Data.vphi    = np.reshape(RawData['vphi'   ],(mpol+1,ntor,lpsi,2),order='F')
				Data.brad    = np.reshape(RawData['brad'   ],(mpol+1,ntor,lpsi,2),order='F')
				Data.btheta  = np.reshape(RawData['btheta' ],(mpol+1,ntor,lpsi,2),order='F')
				Data.bphi    = np.reshape(RawData['bphi'   ],(mpol+1,ntor,lpsi,2),order='F')
				Data.erad    = np.reshape(RawData['erad'   ],(mpol+1,ntor,lpsi,2),order='F')
				Data.etheta  = np.reshape(RawData['etheta' ],(mpol+1,ntor,lpsi,2),order='F')
				Data.ephi    = np.reshape(RawData['ephi'   ],(mpol+1,ntor,lpsi,2),order='F')
				Data.prs     = np.reshape(RawData['prs'    ],(mpol+1,ntor,lpsi,2),order='F')
				Data.rho     = np.reshape(RawData['rho'    ],(mpol+1,ntor,lpsi,2),order='F')
				Data.dns_a   = np.reshape(RawData['dns_a'  ],(mpol+1,ntor,lpsi,2),order='F')
				Data.mom_a   = np.reshape(RawData['mom_a'  ],(mpol+1,ntor,lpsi,2),order='F')
				Data.ppara_a = np.reshape(RawData['ppara_a'],(mpol+1,ntor,lpsi,2),order='F')
				Data.pperp_a = np.reshape(RawData['pperp_a'],(mpol+1,ntor,lpsi,2),order='F')
				Data.qpara_a = np.reshape(RawData['qpara_a'],(mpol+1,ntor,lpsi,2),order='F')
				Data.qperp_a = np.reshape(RawData['qperp_a'],(mpol+1,ntor,lpsi,2),order='F')
			#endif

			#Delete RawData for current KStep and increment KStep counter
	  		del RawData; index += 1
		#endwhile

		#If final index is below requested kstep no data will be saved, causing crashes.
		#Inform user and exit program (adding dummy data might allow for soft crashes)
		if index < kstep:
			MaxKstep = (setting_SEQ[0]+1)*index*(Data.kst[1]-Data.kst[0])
			print('')
			print('#===========================================#')
			print('Requested Kstep exceeds SEQ.harmonics range')
			print('Max Kstep Index: '+str(index)+'   ::   Max Kstep: '+str(MaxKstep))
			print('Please reduce requested Kstep index and retry')
			print('#===========================================#')
			print('')
			exit()
		#endif

	#=====#=====#

	#IF NO kstep supplied, read SINGLE VARIABLE FOR ALL KSTEP values into a 4D HarmonicsData object
	#Read data for each Kstep from the supplied SEQ.harmonics output folder
	#Save data for requested variable for each Kstep in the process
	elif np.isnan(kstep) == True:									#elif: Dimension == '4D'
		kstep = 0													#Initiate kstep index to zero
		while(True):
			#FORTRANFile.read_record automatically steps through all KSteps in the supplied SEQ.harmonics file.
			#RawData for KStep[i] is of shape RawData[Variable][datapoint] where all variables are flattened to 1D
			#Note: Rawdata Variable field is supplied as a string rather than as an index
			try: RawData = FORTRANFile.read_record(DataFormat)		#Read KStep=0, if open then read KStep=1, etc...
			except: FORTRANFile.close(); break						#If no KSteps remain, close FORTRAN file
			#endtry

			#Only extract rho_pol and q_psi on first Kstep
			if kstep == 0:
				Data.rho_pol = np.sqrt(abs(RawData['gpsi_nrm'][0]))			#[-]
				Data.q_psi   = RawData['q_psi'][0]							#[-]
			#Always extract 1D Kstep and Time arrays - Shape = Data.kst[Real]
			Data.kst  = np.append(Data.kst,  RawData['kst'][0])				#[-]
			Data.time = np.append(Data.time, RawData['t'][0])				#[IonGyroFreq*ms]
			#Print kstep for debug purposes if requested
			if DebugMode == True: print(str(i)+'-'+str(data.kst[i]))

			#Extract requested variable and reshape from 1D array into 4D array
			#Data structures of shape: Data.data[kstep][mpol][ntor][lpsi][A/B]
			Data.data=np.concatenate((Data.data, np.reshape(RawData[Variable],(1,mpol+1,ntor,lpsi,2),order='F')))

			#Delete RawData for current KStep and increment KStep counter
	  		del RawData; kstep += 1
		#endwhile
	#endif

	#=====#=====#
	#=====#=====#

	#Debug outputs: 
	if DebugMode == True:
		#Print filename, KStep Range and Data Shape 
		print( '\n  '+Filename.split('/')[-1]+' Data KStep: '+str(0)+'-'+str(Data.kst[-1])    )
		try: print( '  3D image shape [mpol,ntor,lpsi]: '+str(shape(Data.data[0,:,:,:,0]))    )
		except: print( '  3D image shape [mpol,ntor,lpsi]: '+str(shape(Data.vrad[:,:,:,0])) )

		#Extract 2D and 3D images for the requested input variable (Variable)
		try: 	image3D = Data.data[0,:,:,:,0]		#image3D[mpol,ntor,lpsi] for :: SEQ, kstep[0] - One variable 
		except: image3D = Data.vrad[:,:,:,0]		#image3D[mpol,ntor,lpsi] for :: SEQ, kstep[0] - All variables
		image2D = image3D[:,0,:]					#image2D[mpol,0,lphi] for    :: SEQ, kstep[0], ntor[0]

		#Create figure and define Title, Legend, Axis Labels etc...
		fig,ax = figure(image_aspectratio,1)
		Title = 'Harmonic plot of '+Variable+' at Kstep=0 for simulation \n'+str(Filename.split('/')[-3])
		Xlabel,Ylabel = 'Poloidal Resolution $m_{\\theta}$ [cells]', 'Toroidal Resolution $l_{\phi}$ [cells]'
		Legend = list()

		#Plot example data for debugging purposes
		im = ax.imshow(image2D, aspect='auto', origin='bottom')
		cbar = Colourbar(ax,im,Variable+' [-]',5)
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
		plt.show()
		plt.close('all')
	#endif

	#=====#=====#
	#=====#=====#

	return(Data)
#enddef

#=========================#
#=========================#

def Read_MEGAMoments(DataDir,Variable,ntor,kstep=np.nan,SEQ=0):

#	legenda = {'br','bz','bphi','vr','vz','vphi','er','ez','ephi','epara','rho','prs','a '};

	nr = 512; nz = 512;
#	rmin = 1.02; rmax = 2.35; zmin = -1.26; zmax = 1.15; 		%jesus
	rmin = 1.04; rmax = 2.205; zmin = -1.224; zmax = 1.05;
#	[R,Z] = meshgrid(linspace(rmin,rmax,nr),linspace(zmin,zmax,nz));

#	wci = 1.14749E+08; %rad/s
#	wa = 5.80036E+06/(2*pi*1.71863E+00); % 1/s

#	filepath = '/tokp/scratch/javgomar/MEGA_marconi/n4/2/off/';
#	filepath = '/tokp/scratch/javgomar/MEGA_marconi/n4/pr/512/rot_m/';
#	filepath = '/tokp/scratch/javgomar/MEGA_marconi/plasma_response/aug2020/vtor/2/'
#	filepath = '/tokp/scratch/javgomar/mega_interactive/marsf/0/data/'
#	fid = fopen([filepath '001.moments'],'rb');
#	for i in range(0,1):
#		head = fread(fid,1,'int32'); 								# needed by matlab but put in python
#		data(i).kst = fread(fid,1,'int32');
#		data(i).time= fread(fid,1,'float64');
#		aux = fread(fid,nr*nz*16,'float64');
#		data(i).br =  reshape(aux(0*nr*nz+1:1*nr*nz),[nr,nz]);
#		data(i).bz =  reshape(aux(1*nr*nz+1:2*nr*nz),[nr,nz]);
#		data(i).bp =  reshape(aux(2*nr*nz+1:3*nr*nz),[nr,nz]);
#		data(i).vr =  reshape(aux(3*nr*nz+1:4*nr*nz),[nr,nz]);
#		data(i).vz =  reshape(aux(4*nr*nz+1:5*nr*nz),[nr,nz]);
#		data(i).vp =  reshape(aux(5*nr*nz+1:6*nr*nz),[nr,nz]);
#		data(i).er =  reshape(aux(6*nr*nz+1:7*nr*nz),[nr,nz]);
#		data(i).ez =  reshape(aux(7*nr*nz+1:8*nr*nz),[nr,nz]);
#		data(i).ep =  reshape(aux(8*nr*nz+1:9*nr*nz),[nr,nz]);
#		data(i).epara =  reshape(aux(9*nr*nz+1:10*nr*nz),[nr,nz]);
#		data(i).rho =  reshape(aux(10*nr*nz+1:11*nr*nz),[nr,nz]);
#		data(i).prs =  reshape(aux(11*nr*nz+1:12*nr*nz),[nr,nz]);
#		data(i).ppara_i =  reshape(aux(12*nr*nz+1:13*nr*nz),[nr,nz]);
#		data(i).pperp_i =  reshape(aux(13*nr*nz+1:14*nr*nz),[nr,nz]);
#		data(i).ppara_a =  reshape(aux(14*nr*nz+1:15*nr*nz),[nr,nz]);
#		data(i).pperp_a =  reshape(aux(15*nr*nz+1:16*nr*nz),[nr,nz]);
#	    data(i).vac =  reshape(aux(16*nr*nz+1:17*nr*nz),[nr,nz]);
#	    data(i).jr =  reshape(aux(17*nr*nz+1:18*nr*nz),[nr,nz]);
#	    data(i).jz =  reshape(aux(18*nr*nz+1:19*nr*nz),[nr,nz]);
#	    data(i).jp =  reshape(aux(19*nr*nz+1:20*nr*nz),[nr,nz]);
#		head = fread(fid,1,'int32'); 								# needed by matlab but put in python
#	    data(ii).kst
#	   [num2str(i) ' - ' num2str(data(i).time)]
#enddef

#=========================#
#=========================#

def ExtractMEGA_DataShape(HarmonicsData):
#Extracts MEGA SEQ.harmonics data shapes for use with diagnostics
#Determines if HarmonicsData is 3D or 4D and returns appropriate length array
#Inputs: HarmonicsData3D or 4D of shape [kstep][mpol][ntor][lpsi][A/B] - [kstep] optional
#Returns: DataShapes 1D list with contents: [mpol,ntor,lpsi,ltheta,kstep] - kstep optional
#Example: DataShapes = ExtractMEGA_DataShape(HarmonicsData)
#
#lpsi = radial spatial resolution					#[cells]	- set by ????
#ltheta = poloidal angle spatial resolution			#[cells]	- set by ????
#mpol = number of poloidal harmonics considered 	#[int]		- low-pass filter limited?
#ntor = number of toroidal harmonics considered 	#[int]		- low-pass filter limited?
#

	#==========#	#==========#	#==========#

	#Compare if HarmonicsData contains any of the requested variables by name (attribute)
	#N.B. Python 3.7.x objects seem to initiate with 31 generic attributes by default.
#	DataAttributes = dir(HarmonicsData)
#	Intersection = set(variables).intersection(set(DataAttributes))		

	#Extract data array sizes for a 4D (temporally resolved) HarmonicsData object
	#Object created using: ExtractMEGA_Harmonics(Dir[l],'brad',ntor_tot)
#	if len(Intersection) == 0:
#		kstep_res = HarmonicsData.data.shape[0]		#kstep resolution
#		mpol_res = HarmonicsData.data.shape[1]		#mpol resolution
#		ntor_res = HarmonicsData.data.shape[2]		#ntor resolution
#		lpsi_res = HarmonicsData.data.shape[3]		#lpsi resolution
#		ltheta_res = 256							#ltheta resolution 	#!!!!!!!
	#endif

	#Extract data array sizes for a 3D (single kstep) HarmonicsData object
	#NOTE, this assumes brad is an attribute within HarmonicsData - need a 'generic' attribute tag if possible
	#Object created using: ExtractMEGA_Harmonics(Dir[l],'All',ntor_tot,KstepIndex,SEQ)
#	if len(Intersection) >= 1:
#		mpol_res = HarmonicsData.brad.shape[0]		#mpol resolution
#		ntor_res = HarmonicsData.brad.shape[1]		#ntor resolution
#		lpsi_res = HarmonicsData.brad.shape[2]		#lpsi resolution
#		ltheta_res = 256							#ltheta resolution 	#!!!!!!!
	#endif

	#==========#	#==========#	#==========#

	#Extract data array sizes for a 4D (temporally resolved) HarmonicsData object
	#Object created using: ExtractMEGA_Harmonics(Dir[l],'brad',ntor_tot)
	try:
		kstep_res = HarmonicsData.data.shape[0]		#kstep resolution
		mpol_res = HarmonicsData.data.shape[1]		#mpol resolution
		ntor_res = HarmonicsData.data.shape[2]		#ntor resolution
		lpsi_res = HarmonicsData.data.shape[3]		#lpsi resolution
		ltheta_res = 256							#ltheta resolution 	#!!!!!!!

	#Extract data array sizes for a 3D (single kstep) HarmonicsData object
	#NOTE, this assumes brad is an attribute within HarmonicsData - need a 'generic' attribute tag if possible
	#Object created using: ExtractMEGA_Harmonics(Dir[l],'All',ntor_tot,KstepIndex,SEQ)
	except:
		mpol_res = HarmonicsData.brad.shape[0]		#mpol resolution
		ntor_res = HarmonicsData.brad.shape[1]		#ntor resolution
		lpsi_res = HarmonicsData.brad.shape[2]		#lpsi resolution
		ltheta_res = 256							#ltheta resolution 	#!!!!!!!
	#endtry

	#DataShapes contains HarmonicsData resolutions of form: [mpol,ntor,lpsi,ltheta]	[kstep] <-- optional
	try: DataShapes = [mpol_res,ntor_res,lpsi_res,ltheta_res,kstep_res]
	except: DataShapes = [mpol_res,ntor_res,lpsi_res,ltheta_res]

	#Debug outputs: print explicit data resolutions to terminal
	if DebugMode == True:
		print( 'mpol:',mpol_res )
		print( 'ntor:',ntor_res )
		print( 'lpsi:',lpsi_res )
		print( 'ltheta:',ltheta_res )
		try: print( 'kstep_len',kstep_res )
		except: dummyvar = 1
	#endif

	return(DataShapes)
#enddef

#=========================#
#=========================#

def ExtractMEGA_DataRanges(Dir, DataFile='energy_n'):
#Extracts MEGA temporal axes and toroidal mode number information from specified output file
#Used for extracting axes to plot against or to determine toroidal mode number indices
#Inputs: 
#	Dir - 0D string containing simulation directory (local or root)
#	DataFile - 0D string identifying the MEGA output file to be used for the time axis
#				DataFile options: 'energy_n', 'harmonics', 'moments'
#Outputs:
#	KStepArray - 1D array containing all simulation ksteps over all detected SEQs
#	TimeArray - 1D array containing all simulation times [s] over all detected SEQs
#	ntorArray - 1D array containing: ntor0, ntor_pos, and ntor_tot
#				ntor0 - 0D integer indicating which ntor value represents the n=0 equilibrium data
#				ntor_pos - 0D integer indicating the number of positive modes (Ignoring n=0)
#				ntor_tot - 0D integer indicating the total number of positive and negative modes (Including n=0)
#Example: KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir, DataFile='harmonics')

	#Extract Filename.txt paths for all SEQ for given data filename
	DataSubDir = 'data/'
	Files = sorted(glob.glob( Dir+DataSubDir+'*'+DataFile+'*' ))
	SEQArray = range(0,len(Files))

	#Extract kstep, time and toroidal harmonic data from energy_n.txt
	#energy_n data structure: [variable][timestep]
	Energy_n,Header_n = ExtractMEGA_Energy(Dir, 'energy_n')
	#Determine poloidal and toroidal harmonic ranges
	ntor_tot = ((len(Energy_n)-3)*2)+1			#Total number of positive and negative modes (Including n=0)
	ntor_pos = int(float(ntor_tot-1)/2.0)		#Number of positive modes (Ignoring n=0)
	ntor0 = int(ceil(ntor_tot/2))				#ntor = 0, baseline equilibrium data

	#Determine KStep range, Time range and related intervals
	#Use wchck time intervals
	if DataFile == 'energy_n':
		KStepArray = Energy_n[0]					#KStep Array			[-]
		TimeArray = Energy_n[1]						#Time Array				[ms]

	#Use wharm time intervals
	elif DataFile == 'harmonics':
		#READ ONLY TEMPORAL AND SPATIAL AXES from SEQ.Harmonics, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		HarmonicsData = ExtractMEGA_Harmonics(Dir, Variable='NaN', ntor=ntor_tot, Dimension='1D')
		rho_pol = HarmonicsData.rho_pol				#Normalised radius		[-]

		KStepArray = HarmonicsData.kst				#KStep Array			[-]
		TimeArray = HarmonicsData.time				#Time Array				[ms]

	#Use wsnapshot time intervals
	elif DataFile == 'moments':
#		KStepArray = MomentsData.kst				#KStep Array			[-]
#		TimeArray = MomentsData.time				#Time Array				[ms]
		a = 1										#TO BE COMPLETED

	#Use wep time intervals
	elif DataFile == 'markers':
#		KStepArray = KineticData.kst				#KStep Array			[-]
#		TimeArray = KineticData.time				#Time Array				[ms]
		a = 1										#TO BE COMPLETED
	#endif

	return(SEQArray, KStepArray, TimeArray, [ntor0,ntor_pos,ntor_tot])
#enddef

#=========================#
#=========================#

def ExtractMEGA_SharedDataRange(DirList):
#Determines the maximum shared KStepArray length between all simulation folders
#Returns max shared KStep index and associated Dir index for the relevant simulation folder
#Inputs:
#	DirList - 1D string array containing all simulation folder root directories
#Outputs:
#	MaxSharedKStep - 0D Scalar indicating the maximum shared KStep array length
#	MaxSharedDirIdx - 0D Scalar indicating the Dir[Idx] for the smallest simulation folder
#Example: MaxSharedKStep,MaxSharedDirIdx = ExtractMEGA_SharedDataRange(Dir[l])

	#Initiate any required lists
	TimeArray_Lengths = list()

	#For each detected simulation folder, record the KStep array length
	for l in range(0,len(DirList)):

		#Extract Kstep [-], Time [ms] & toroidal harmonics from energy_n.txt
		SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[l], DataFile='energy_n')
		TimeArray_Lengths.append( len(TimeArray) )
	#endfor

	#Extract the minimum shared TimeArray length and associated Dir Array index
	MaxSharedKStep = min(TimeArray_Lengths)						#Maximum shared KStep array length
	MaxSharedDirIdx = TimeArray_Lengths.index(MaxSharedKStep)	#Dir[Idx] of simulation with minimum length
	
	return(MaxSharedKStep, MaxSharedDirIdx)
#enddef

#=========================#
#=========================#

def ExtractMEGA_PoloidalGrid(Dir,HarmonicsData):
#Extracts poloidal axes from repository .dat files using harmonics data resolution
#Inputs: HarmonicsData object of shape [mpol][ntor][lpsi][A/B], Repository directory
#Returns: Radial (crdr) and axial (crdz) magnetic axes
#Example: crdr,crdz = ExtractMEGA_PoloidalGrid('Repository',HarmonicsData):

	#NOTE	:: NEED TO EXTRACT DATA FROM SEPERATRIX FILE AND COMPARE TO Crdr and Crdz ARRAYS
	#		:: CAN ALSO EXTRACT FROM psi.dat FILE, WHICH SHOULD CONTAIN grrsim4 ARRAYS (If I remember correctly)
	#		:: SEE toMEGA.f OUTPUT FILES FOR FURTHER DETAILS

	#Extract data shape from supplied data object
	DataShape = ExtractMEGA_DataShape(HarmonicsData)

	mpol_res = DataShape[0]
	ntor_res = DataShape[1]
	lpsi_res = DataShape[2]
	ltheta_res = DataShape[3]
#	kstep_res = DataShape[4]

	#Load poloidal mesh grid from repository
	rho_a = np.loadtxt(Dir+'/rho_a.dat')
	Crdr = np.loadtxt(Dir+'/crdmagr.dat').reshape((lpsi_res,ltheta_res),order='F')*rho_a
	Crdz = np.loadtxt(Dir+'/crdmagz.dat').reshape((lpsi_res,ltheta_res),order='F')*rho_a
	Crdr = np.concatenate((Crdr,Crdr[:,0:1]),axis=1)
	Crdz = np.concatenate((Crdz,Crdz[:,0:1]),axis=1)

	return(Crdr,Crdz)
#enddef

#=========================#
#=========================#

def ExtractMEGA_Normalisations(Dir):
#Takes simulation folder directory (absolute path) and returns Sim128 normalisation constants
#Example: Variables,Values,Units = ReadNormConstants(Dir[l])
#
# NOTE: Duplicated variable names in output file --- ADDRESS BY SPLITTING Sim128 FILE INTO SECTIONS
#'D beam inj. vlc.','Crit. vlc. axis','SlowD rate axis'				--- ON TOP AND POST NORM SETS
#'psimax','major_r','rleng','left','right','zleng','raxis','zaxis'	--- ON PRE AND POST NORM SETS
#

	#Normalisation constants are stored within: sim128-aug<Shot>.<num>.txt
	#Location of sim128-aug is typically within folder named 'equil':
	try: sim128File = sorted(glob.glob(Dir+'equil/*sim128-aug*txt'))[0]
	except: sim128File = sorted(glob.glob(Dir+'*sim128-aug*txt'))[0]
	sim128Data = open(sim128File).readlines()

	#Manually define variable names in sim128 file --- ADD FUNCTIONALITY TO ENABLE USER SELECTION
	TopVariables = ['Mag.fld. at axis','Bulk density','Alfven velocity','D gyro frequency','Alfv gyro radius','SlowD time axis']
	PreNormVariables = ['psimax','major_r','rleng','left','right','zleng','bottom_sim','top_sim','raxis','zaxis']
	PostNormVariables = ['psimax','major_r','rleng','left','right','zleng','raxis','zaxis','D beam inj. vlc.','Crit. vlc. axis','SlowD rate axis','maximum ion temperature','maximum elc temperature',]
	InputVariables = TopVariables + PreNormVariables + PostNormVariables
	
	#Initialise variable, values and units output lists
	Variables,Values,Units = list(),list(),list()

	#Identify variable name in sim128 file and strip value and unit
	for i in range(0,len(InputVariables)):
		Variable = InputVariables[i]
	
		#--- ADD FUNCTIONALITY TO RE-NAME VARIABLES FROM PRE-WRITTEN LIST
		Variables.append(Variable)

		Value = filter(lambda x:Variable in x, sim128Data)[0].strip(' \t\n\r,='+Variable)
		try: Value = float(Value)
		except: Value = float(Value[0:11])
		Values.append(Value)

		Unit = filter(lambda x:Variable in x, sim128Data)[0].strip(' \t\n\r,='+Variable)
		try: Unit = '['+Unit.split('[')[1]
		except: Unit = '[-]'
		Units.append(Unit)
	#endfor

	#Print debug output to terminal if requested
	if DebugMode == True:
		for i in range(0,len(Variables)): print Variables[i], Values[i], Units[i]
	#endif

	return(Variables,Values,Units)
#enddef

#=========================#
#=========================#

def ExtractMEGA_Energy(Dir,Filename='energy_n'):
#Reads and concatenates MEGA energy.txt output files
#Takes simulation directory (absolute path) and filename (energy_n, energy_phys)
#Returns output data and header, data of form: [Variable][Timestamp]
#Example: OutputData,Header = ExtractMEGA_Energy('LocalDataFolder/','energy_phys'])

	#Extract Filename.txt paths for all SEQ for given data filename
	DataSubDir = 'data/'
	Files = sorted(glob.glob( Dir+DataSubDir+'*'+Filename+'*' ))

	#For each output file in the current simulation directory:
	for SEQ in range(0,len(Files)):
		#Extract header and output data for first SEQ
		if SEQ == 0:
			Header = ReadFile_ASCII(Files[SEQ],1,'2D','CSV')[1]
			OutputData = ReadFile_ASCII(Files[SEQ],1,'2D','CSV')[0]
		#Extract output data for subSEQuent SEQ's and append to each variable
		elif SEQ > 0:
			TempData = ReadFile_ASCII(Files[SEQ],1,'2D','CSV')[0]
			for j in range(0,len(TempData)):
				OutputData[j] = np.concatenate( (OutputData[j],TempData[j]) )
			#endfor
		#endif

		#Debug outputs: Print datafile name, number of variables, length of variable arrays
		if DebugMode == True:
			print Files[l].split('/')[-1]
			print len(OutputData), len(OutputData[0])
		#endif
	#endfor
	
	return(OutputData,Header)
#endif

#=========================#
#=========================#

def ExtractMEGA_Harmonics(Dir,Variable,ntor,kstep=np.nan,SEQ=0,Dimension='1D'):
#Details on the FORTRAN file format can be found below:
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.FortranFile.read_record.html

	#Extract harmonic output files for all SEQ in requested directory
	DataSubDir = 'data/'
	DataFiles = sorted(glob.glob(Dir+DataSubDir+"/*harm*"))

	#Extract data radial (psi) and poloidal mode number (mpol) resolutions 		
	lpsi,mpol = 201, 64											#!!! HARD-CODED FOR NOW !!!

	#=====#=====#	#=====#=====#

	#IF Kstep == -1, is supplied ONLY READ TEMPORAL AND SPATIAL AXES into a 1D HarmonicsData object
	#Extract Kstep and Time arrays for plotting purposes, ignoring all other data.
	if Dimension == '1D':
		Variable = 'NaN'		#Variable is a dummy input here - Not used.
		#Object Shape: HarmonicsData[kstep][mpol][ntor][lpsi][A/B]
		HarmonicData = list()
		for i in range(0,len(DataFiles)):
			HarmonicData.append( Read_MEGAHarmonics(DataFiles[i],Variable,mpol,ntor,lpsi,kstep=-1) )
		#endfor

		#Concatenate data from all SEQs into one continuous array for each variable within HarmonicData object
		for i in range(0,len(HarmonicData)-1):
			pop = HarmonicData.pop(1)										#'Pops' and removes 1st SEQ data array
			HarmonicData[0].kst  = np.append(HarmonicData[0].kst, pop.kst)	#Appends to zero'th SEQ data array
			HarmonicData[0].time = np.append(HarmonicData[0].time, pop.time)
			del pop															#Refresh pop array and repeat
		#endfor
		HarmonicData = HarmonicData[0]		#Replace data object with fully appended (i.e. flattened) data array
	#endif

	#=====#=====#

	#IF Kstep index > 0 is supplied, read ALL VARIABLES FOR SINGLE KSTEP into a 3D HarmonicsData object
#	elif isinstance(kstep, int) and isinstance(SEQ, int):
	elif Dimension == '3D': 
		#Object Shape: HarmonicsData[mpol][ntor][lpsi][A/B]
		Variable = 'NaN'		#Variable is not used in extraction of 3D data, but is required for denormalisation
		HarmonicData = Read_MEGAHarmonics(DataFiles[SEQ],Variable,mpol,ntor,lpsi,kstep)
	#endif

	#=====#=====#

	#IF NO kstep index supplied, read SINGLE VARIABLE FOR ALL KSTEP values into a 4D HarmonicsData object
#	elif np.isnan(kstep) == True or np.isnan(SEQ) == True:
	elif Dimension == '4D':
		#Object Shape: HarmonicData [kstep][mpol][ntor][lpsi][A/B]
		HarmonicData = list()
		for i in tqdm(range(0,len(DataFiles))):
			HarmonicData.append(Read_MEGAHarmonics(DataFiles[i],Variable,mpol,ntor,lpsi))
		#endfor

		#Concatenate data from all SEQs into one continuous array for each variable within HarmonicData object
		for i in range(0,len(HarmonicData)-1):
			pop = HarmonicData.pop(1)										#'Pops' and removes 1st SEQ data array index
			HarmonicData[0].kst  = np.append(HarmonicData[0].kst, pop.kst)	#Appends to zero'th SEQ data array
			HarmonicData[0].time = np.append(HarmonicData[0].time, pop.time)
			HarmonicData[0].data = np.concatenate((HarmonicData[0].data, pop.data))
			del pop															#Refresh pop array and repeat
		#endfor
		HarmonicData = HarmonicData[0]		#Replace data object with fully appended (i.e. flattened) data array

	#=====#=====#	#=====#=====#

	DenormaliseAtReadin = True
	#De-normalise data if requested			# TO BE COMPLETED - MAY MOVE INTO EXTRACT FUNCTION???
	if DenormaliseAtReadin == True:

		#Remove '/data/' from directory --> Dir now points to simulation root folder
		#		  Reverse       split into 2    Keep preamble  re-reverse   +'/' on end
		NormDir = DataDir[::-1].split('/', 2)   [-1]           [::-1]       +'/'
		#This is disgusting and I apologize to anyone who has to read this...

		#Extract relevant normalisation factors for current simulation folder
		NormVariables,NormValues,NormUnits = ExtractMEGA_Normalisations(NormDir)

		RMax = NormValues[NormVariables.index('right')]							#	Update to right_sim
		RMin = NormValues[NormVariables.index('left')]							#	Update to left_sim
		Raxis = NormValues[NormVariables.index('raxis')]						#
		Rlen = NormValues[NormVariables.index('rleng')]							#

		ZMax = NormValues[NormVariables.index('top_sim')]						#
		ZMin = NormValues[NormVariables.index('bottom_sim')]					#
		Zaxis = NormValues[NormVariables.index('zaxis')]						#
		Zlen = NormValues[NormVariables.index('zleng')]							#

		AlfvenVelocity = NormValues[NormVariables.index('Alfven velocity')] 	#B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = NormValues[NormVariables.index('D gyro frequency')]		#
		IonDensity = NormValues[NormVariables.index('Bulk density')]			#
		B0 = NormValues[NormVariables.index('Mag.fld. at axis')]				#

		Mass_Deuterium = 3.34e-27												#
		eps = 0.5/Raxis															#

		#=====#=====#

#		Phys = ['vrad','vtheta','vphi','brad','btheta','bphi','erad','etheta','ephi','prs','rho']
#		Kin = ['dns_a','mom_a', 'ppara_a','pperp_a','qpara_a','qperp_a']

		#Denormalise temporal and spatial axes
		HarmonicData.kst = HarmonicData.kst								# [-]
		HarmonicData.time = HarmonicData.time * (1e3/IonGyroFreq)		# [ms]

		#Denormalise all variables		# TO BE COMPLETED
		try: 
			HarmonicData.brad = HarmonicData.brad*B0*1000			# [mT]
			HarmonicData.btheta = HarmonicData.btheta*B0*1000		# [mT]
			HarmonicData.bphi = HarmonicData.bphi*B0*1000			# [mT]

			HarmonicData.erad = HarmonicData.erad*1.0				# WHAT IS E-FIELD NORMALISED TO
			HarmonicData.etheta = HarmonicData.etheta*1.0			# WHAT IS E-FIELD NORMALISED TO
			HarmonicData.ephi = HarmonicData.ephi*1.0				# WHAT IS E-FIELD NORMALISED TO

			HarmonicData.vrad = HarmonicData.vrad*1.0				# WHAT IS VELOCITY NORMALISED TO
			HarmonicData.vtheta = HarmonicData.vtheta*1.0			# WHAT IS VELOCITY NORMALISED TO
			HarmonicData.vphi = HarmonicData.vphi*1.0				# WHAT IS VELOCITY NORMALISED TO		

			HarmonicData.rho = HarmonicData.rho*IonDensity			# [m-3]
			HarmonicData.prs = HarmonicData.prs*1.0					# WHAT IS PRESSURE NORMALISED TO
		except: 
			NormFactor = 1.0	#NEED TO DETERMINE THIS FROM DATA VARIABLE
			HarmonicData.data = HarmonicData.data*NormFactor			# [-]
		#endtry
	#endif


	return(HarmonicData)
#enddef

#=========================#
#=========================#

def ExtractMEGA_Moments(DataDir,Variable,ntor,kstep=np.nan,SEQ=0,Dimension='1D'):
#Details on the FORTRAN file format can be found below:
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.FortranFile.read_record.html

	#TO BE WRITTEN
	a = 1

	return(a)
#enddef

#=========================#
#=========================#

def ExtractMEGA_Markers(Dir,KStep,MarkerFileStep=1):
# Reads data from kinetic marker output files (gc_a_kstep000xxxx-00xxx.txt)
# Reads all variables and concatenates output data from all cores into single 2D array
# Inputs: 
#	Dir - Directory String to marker output file folder from root
# 	KStep - KStep value (NOT Index) of output files to be read-in
# 	MarkerFileStep - Optional speedup input, will read every 'n'th output file
# Outputs: 
#	KineticsData - 2D Array of shape [variable][marker(n)]
#	Variables - R, Z, Lambda, E, p, Mu, pphi, fff, fnrml, psip, phi
# Example :: KineticsData,Header_Kin = ExtractMEGA_Markers(Dir[l],KStep=1000,MarkerFileStep=8)

	#Initiate kinetic data array
	KineticsData = list()

	#Define marker folder location and filename format for supplied KStep
	Folder = Dir+'markers/'
	Filename = 'gc_a_kstep'+str('%07.f'%KStep)+'-*.txt'
	#Sort all marker output files numerically by core number
	MarkerFiles = sorted(glob.glob(Folder+Filename))

	#Exit cleanly if no files are found:
	if len(MarkerFiles) == 0:
		print ''
		print '-------------------------------------------'
		print 'No Marker Files Detected, Aborting Analysis'
		print '-------------------------------------------'
		print ''
		exit()
	#endif

	#Cycle through all marker output files (nfiles = ncores) for the current KStep and read data
	for j in tqdm( range(0,len(MarkerFiles),MarkerFileStep) ):

		#Set current marker output file	(j = core number)
		Filename = MarkerFiles[j]

		#Read marker data for current NCore output file
		#MarkerData :: 2D array of shape [variable,marker(n)]
		#Variables :: 0:R, 1:Z, 2:Lambda, 3:E, 4:p, 5:Mu, 6:pphi, 7:fff*fnrml, 8:psip, 9:phi
		#See write_ep() subroutine for full details
		MarkerData,Header = ReadFile_ASCII(Filename, 0, '2D', 'CSV')
#		print np.asarray(MarkerData).shape


		#Concatenate MarkerData from all cores into KineticsData for each timestep
		#KineticsData :: 2D Array of shape [variable,marker(n)]
		if len(KineticsData) == 0: 
			KineticsData = MarkerData					#Initiate KineticsData array on first iteration
		elif len(KineticsData) > 0: 
			KineticsData = [KineticsData[x]+MarkerData[x] for x in range(0,len(KineticsData))]
		#endif

		#Print debug outputs to terminal if requested
		if DebugMode == True: 
			print( Filename.split('/')[-1] )
			print( 'Kin Num Variables: '+str(len(KineticsData)) )
			print( 'Kin Variable Length: '+str(len(KineticsData[0])) )
			print( '' )
		#endif
	#endfor

	return(KineticsData,Header)
#enddef

#=========================#
#=========================#

def Extract_PoloidalImage(HarmonicsData,VariableString,ntor):
#Reduces 3D HarmonicsData [mpol][ntor][lpsi][A/B] into 2D poloidal image [mpol,ltheta]
#Reduces 3D HarmonicsData by extracting only 1 ntor and averaging over mpol
#Also merges the real and imaginary components of HarmonicsData
#Inputs: HarmonicsData object [mpol][ntor][lpsi][A/B], 
#Inputs: VariableString matching HarmonicsData attribute, ntor Index (not mode number)
#Outputs: 2D PoloidalData2D array of shape [lpsi,ltheta]
#Example: PoloidalImage = Extract_PoloidalImage(HarmonicsData,'vrad',1)		

	#Extract data shape from supplied data object
	DataShape = ExtractMEGA_DataShape(HarmonicsData)
	mpol_res = DataShape[0]
	ntor_res = DataShape[1]
	lpsi_res = DataShape[2]
	ltheta_res = DataShape[3]
#	kstep_res = DataShape[4]

	#Extract variable attribute from HarmonicsData object and initiate empty 2D data array
	HarmonicsData3D = getattr(HarmonicsData, VariableString)
	PoloidalData2D = np.ndarray((lpsi_res,ltheta_res),np.float64)

	#For each ltheta cell, compute poloidal angle and sum over all poloidal harmonics
	for ltheta in range (0,ltheta_res,1):

		#Poloidal angle is computed as fraction of ltheta, rotating clockwise from the midplane
		PoloidalAngle = 2*np.pi*(float(ltheta)/float(ltheta_res)) 	#Poloidal angle [Rad]
		ToroidalAngle = 0											#Toroidal angle [Rad]		#Not Implimented

		#For each poloidal mode number, extract harmonics data for the requested ntor
		aux = 0			#Initiate accumulator to zero for each ltheta
		for mpol in range (0,mpol_res,1):
#      			    	HarmonicsData3D[mpol]   [ntor]  [lpsi]  [A/B]
			aux = aux + HarmonicsData3D[mpol,   ntor,   :,      0]* np.cos(PoloidalAngle*mpol+ToroidalAngle*ntor)	#1D radial slice
			aux = aux + HarmonicsData3D[mpol,   ntor,   :,      1]* np.sin(PoloidalAngle*mpol+ToroidalAngle*ntor)	#1D radial slice
			#Append 1D radial slice (lpsi) at current poloidal angle (ltheta) to PoloidalData2D
			PoloidalData2D[:,ltheta] = aux
		#endfor
	#endfor

	#Add a zero degree radial slice by appending the final datapoint (poloidal angle is cyclic)
	#PoloidalData2D has shape: [lpsi,ltheta] - data extracted for a single ntor
	PoloidalData2D = np.concatenate((PoloidalData2D, PoloidalData2D[:,0:1]),axis=1)

	return(PoloidalData2D)
#enddef

#=========================#
#=========================#

def Extract_RadialProfile(HarmonicsData,variable,ntorIdx,theta):
#Extracts radially resolved profiles for a single output variable at a single poloidal angle
#Rounds poloidal angle down - i.e. anticlockwise - as defined from vertical zero.
#Inputs: 
#	HarmonicsData - 4D Object of shape [mpol][ntor][lpsi][A/B], for data at single KStep/Time 
#	VariableString - 0D string matching the desired HarmonicsData variable attribute, e.g. 'brad' 
#	ntorIdx - 0D integer determining the toroidal mode number Index, not absolute mode number
#	Theta - 0D float, determining the poloidal angle of the radial profile [degrees]
#Outputs: 
#	RadialProfile - 1D array of shape [lpsi] (rho_pol), containing the poloidal mode merged variable amplitudes
#Example: Profile = Extract_RadialProfile(HarmonicsData,'brad',ntorIdx=2,theta=64)

	#Select variable and Merge 3D Data into 2D poloidal slice
	#PoloidalImage :: 2D array of Shape [lpsi][ltheta] ~ [R][theta]
	#Image[n][:] = full poloidal profile (clockwise from vertical) for R = Rho_pol[n]
	PoloidalImage = Extract_PoloidalImage(HarmonicsData,variable,ntorIdx)	

	#Rotate Polodialimage into shape: [ltheta][lpsi]
	PoloidalImageTrans = PoloidalImage.transpose()		

	#Define requested poloidal angle
	PoloidalAngle = theta						#Poloidal angle 	 [Degrees]

	#Define one degree of rotation in terms of poloidal cells
	FullRotation = len(PoloidalImageTrans)		#Poloidal angle 	 [Cells]
	OneDegree = float(FullRotation)/float(360)	#Poloidal resolution [Cells / Degree]

	#Define requested poloidal index - rounded down (i.e. rounded anti-clockwise)
	thetaIdx = int(PoloidalAngle*OneDegree)		#Poloidal index 	 [Cell]
	thetaIdx = thetaIdx % FullRotation			#Apply Periodicity

	#Extract 1D radial profile for requested variable
	#Radial Profile :: 1D array of shape [Rho_pol]
	RadialProfile = PoloidalImageTrans[thetaIdx]

	#Print useful debugging outputs to terminal if requested
	if DebugMode == True: 
		print(Dir[l],variable)
		print('Poloidal Resolution: ', FullRotation,' [Cells]')
		print('Requested theta: ',Requested_Angle,' [Degrees]')
		print('Theta Index: ', thetaIdx,' [Cell]')
		print('Radial Profile Length: ', len(Radialprofile),' [Cells]')
	#endif

	return(RadialProfile)
#enddef

#=========================#
#=========================#

def Extract_PoloidalProfile(HarmonicsData,variable,ntorIdx,Radius):
#Extracts poloidally resolved profiles for a single output variable at a single radius
#Rounds radius to closest avaliable rho_pol - automatically limited to 0.0 < rho_pol < 1.0
#Inputs: 
#	HarmonicsData - 4D Object of shape [mpol][ntor][lpsi][A/B], for data at single KStep/Time 
#	VariableString - 0D string matching the desired HarmonicsData variable attribute, e.g. 'brad' 
#	ntorIdx - 0D integer determining the toroidal mode number Index, not absolute mode number
#	Radius - 0D float, determining the radius from Rgeo of the poloidal profile [rho_pol]
#Outputs: 
#	PoloidalProfile - 1D array of shape [ltheta], containing the poloidal mode merged variable amplitudes
#Example: Profile = Extract_PoloidalProfile(HarmonicsData,'brad',ntorIdx=2,Radius=0.2)

	#Extract data resolution and poloidal axes from repository .dat files
	#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
	Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)
	DataShape = ExtractMEGA_DataShape(HarmonicsData)#; print DataShape
	mpol_res = DataShape[0]; ntor_res = DataShape[1]
	lpsi_res = DataShape[2]; ltheta_res = DataShape[3]

	#Determine radial location from user supplied switchboard values
	#Radius is in normalised radius [rho_pol], while RadialLoc is in [m]
	RadialLoc = min(rho_pol, key=lambda x:abs(x-Radius))	#[m]
	RadialIdx = rho_pol.tolist().index(RadialLoc)			#[-]

	#Compute poloidal angle axis from length of radial array
	#i.e. assuming a clockwise angle from vertical zero
	#Crdr :: 2D array of Shape [lpsi][ltheta] ~ [R][theta]
	Crdtheta = list()
	for i in range(0,len(Crdr[RadialIdx])):
#		Crdtheta = np.arctan(Crdz[RadialIdx]/Crdr[RadialIdx])	#tan(theta) = Z/R  
		dtheta = 360.0/len(Crdr[RadialIdx])
		Crdtheta.append( i*dtheta )
	#endfor

	#==========##==========#

	#Select variable and Merge 3D Data into 2D poloidal slice
	#PoloidalImage :: 2D array of Shape [lpsi][ltheta] ~ [R][theta]
	#Image[n][:] = full poloidal profile (clockwise from vertical) for R = Rho_pol[n]
	PoloidalImage = Extract_PoloidalImage(HarmonicsData,variable,ntorIdx)
	
	#Extract single poloidal profile at Radius = RadialLoc
	ThetaArray = Crdtheta[:]
	DataArray = PoloidalImage[RadialIdx][:]

	return(ThetaArray, DataArray)

	#==========##==========#
	#==========##==========#

	ReturnInboardOutboard = False
	if ReturnInboardOutboard == True:
		#One full poloidal rotation corresponds to len(Image[n][:]) cells
		FullRot = len(PoloidalImage[0])					#Total number of Poloidal Cells
		HalfRot = FullRot/2								#Inboard/outboard Poloidal Cells

		#Extract inboard and outboard poloidal profiles for Radius = RadialLoc
		#Both PoloidalImage and Crdr are of length [lpsi] and rotate clockwise from vertical
		#PoloidalImage[R][0] = '12 O'clock' position, PoloidalImage[R][Half] = '6 o'clock' position
		#PoloidalImage :: 2D array of Shape [lpsi][ltheta] ~ [R][theta]
		#Crdtheta :: 1D array of Shape [ltheta]
		ThetaInboard = Crdtheta[0:HalfRot]							#X-Axis	(Theta at Radius R)
		DataInboard = PoloidalImage[RadialIdx][0:HalfRot]			#Y-Axis	(Variable Amplitude)

		ThetaOutboard = Crdtheta[HalfRot:FullRot]					#X-Axis	(Theta at Radius R)
		DataOutboard = PoloidalImage[RadialIdx][HalfRot:FullRot]	#Y-Axis	(Variable Amplitude)

		#Select variable and Merge 3D Data into 2D poloidal slice
		#Inboard/Outboard = 2D array containing:
		#		- poloidal angle axis (clockwise from vertical) for R = Radius
		#		- poloidal variables[i] amplitude profile (clockwise from vertical) for R = Radius
		Inboard,Outboard = Extract_PoloidalProfile(HarmonicsData,variables[i],ntorIdx,Radius)

		#Set the same colour for both inboard and outboard portions of the profile
		ColourCycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
		Colour = ColourCycle[j % len(ColourCycle)]

		#Plot the inboard and outboard portions of the poloidal profile at Radius = RadialLoc
		#The inboard portion of the plot is given a different linestyle to aid viewing
		Inboard = plt.plot(Inboard[0], Inboard[1], color=Colour,ls='-')
		Outboard = plt.plot(Outboard[0], Outboard[1], color=Colour,ls='--',label='_nolegend_')
		plt.show()

		return(Inboard,Outboard)
	#endif
#enddef

#=========================#
#=========================#

def Set_SEQRange(setting_SEQ):
#Convert user switchboard 'setting_SEQ' settings into SEQ loop index ranges
#Attempts to set SEQ index range to user inputs, else defaults to SEQ = 001.
#Inputs:
#	setting_SEQ - 1D array (len 2) containing the switchboard inputs for KStep range and step size
#Outputs:
#	SEQRange - 1D array (len 2) containing the min and max Kstep indices to loop over
#Example: SEQRange = Set_SEQRange(setting_SEQ)
#		  for i in range( SEQRange[0],SEQRange[1]): <Loop>

	#Apply user SEQ range if requested...
	if len(setting_SEQ) == 2: 
		SEQRange = [setting_SEQ[0],setting_SEQ[1]+1]

	#...else default to SEQ = 001
	else:
		SEQRange = [0,1]
		print('--------------------------')
		print('SEQ Range Being Set to 001')
		print('--------------------------')
	#endif

	return(SEQRange)
#enddef

#=========================#
#=========================#

def Set_KStepRange(KStepArray,setting_kstep):
#Convert user switchboard 'setting_kstep' settings into KStep loop index ranges
#Attempts to set Kstep index range and step, else defaults to max range.
#Inputs:
#	KStepArray - 1D array containing the KStep values for the associated data to be looped over
#	setting_kstep - 1D array (len 2 or 3) containing the switchboard inputs for KStep range and step size
#Outputs:
#	KStepRange - 1D array (len 2) containing the min and max Kstep indices to loop over
#	KStepStep - 0D scalar containing the KStep loop interval size.
#Example: KStepRange,KStepStep = Set_KStepRange(KStepArray,setting_kstep)
# 		  for i in range( KStepRange[0],KStepRange[1], KStepStep ): <Loop>

	#Apply user Kstep range and step if possible...
	if len(setting_kstep) == 3: 
		KStepRange = [setting_kstep[0],setting_kstep[1]]
		KStepStep = setting_kstep[2]

	#...or apply user Kstep range with default step size...
	elif len(setting_kstep) == 2:
		KStepRange = setting_kstep
		KStepStep = 1

	#...else default to max range and default step size
	else: 
		KStepRange = [0,len(KStepArray)]
		KStepStep = 1
		print('--------------------------------')
		print('KStep Range Being Set to Maximum')
		print('--------------------------------')
	#endif

	return(KStepRange, KStepStep)
#enddef

#=========================#
#=========================#

def Set_ntorIdx(ntor,ntorArray):
#Takes desired toroidal mode number and returns the index relating to said number
#Inputs:
#	ntor - 0D integer determining the actual toroidal mode number desired
#	ntorArray - 1D array containing: ntor0, ntor_pos, and ntor_tot for the data in question
#				ntor0 - 0D integer indicating which ntor value represents the n=0 equilibrium data
#				ntor_pos - 0D integer indicating the number of positive modes (Ignoring n=0)
#				ntor_tot - 0D integer indicating the total number of positive and negative modes (Including n=0)
#Outputs:
#	ntorIdx - 0D integer which determines the index of the requested ntor input
#Notes:
#	Designed to take `ntorArray` directly from ExtractMEGA_DataRanges() function and relates to associated data
#	e.g. extract ntorArray from 'harmonics' data, plug into this function to determine ntor index for 'harmonics'
#Example: ntorIdx = Set_ntorIdx(ntor,ntorArray)

	#Unpack input ntorArray
	ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
	ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
	ntor0 = ntorArray[0]						#ntor = 0, baseline equilibrium data

	#Create 2D array containing [ntor,Set_ntorIdx] for referencing data to be extracted
	ntor_indices = list()
	for i in range(0,ntor_tot):
		#ntor_indices contains the following [ntor, Set_ntorIdx], for positive, negative and n=0 modes
		if i-ntor_pos >= setting_ntor[0] and i-ntor_pos <= setting_ntor[1]:
			ntor_indices.append( [i-ntor_pos,i] )
		else:
			ntor_indices.append( [i-ntor_pos,np.nan] )
		#endif
	#endfor

	#Print debug output to terminal if requested
	if DebugMode == True: 
		print('')
		print('Toroidal Mode Numbers:',ntor_indices)
		print('')

	#ntor range set by ntor_indices[ntor, Set_ntorIdx], for pos, neg & n=0 modes
#	ntor = 0													#override requested ntor mode number
	ntorIdx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number

	return(ntorIdx)
#enddef

#====================================================================#
#====================================================================#






#====================================================================#
					#COMMON PLOTTING FUNCTIONS#
#====================================================================#

def Matplotlib_GlobalOptions():
#Takes global inputs from switchboard, returns nothing
#Alters global image options, run before any diagnostics
#Attempts to revert matplotlib changes made in 2.0 onwards.
#See: https://matplotlib.org/users/dflt_style_changes.html

#	mpl.style.use('classic')								#Resets to classic 1.x.x format

	#Import external colourmaps - NOT CURRENTLY WORKING, SEE BELOW:
	Dir = os.getcwd()+'/Repository/IDL_gamma_II.txt'
	IDL_Gamma_II = ReadCmap_ASCII(Dir,'IDL_Gamma_II')		#Import PSFT group colourmap
#	Need code to add IDL_Gamma_II to mpl.rcParams default image.cmap list

	#Image options			
	mpl.rcParams['figure.figsize'] = [10.0,10.0]			#Sets default figure size
	mpl.rcParams['figure.dpi'] = 100						#Sets viewing dpi
	mpl.rcParams['savefig.dpi'] = 100						#Sets saved dpi
	mpl.rcParams['image.interpolation'] = 'bilinear'		#Applies bilinear image 'smoothing'
	mpl.rcParams['image.resample'] = True					#Resamples data before colourmapping
	mpl.rcParams['image.cmap'] = image_cmap					#Select global colourmap 
	#Common cmap overrides: 'plasma', 'gnuplot', 'jet',

	#Axis options
	mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'	#View limits coencide with axis ticks
	mpl.rcParams['axes.xmargin'] = 0						#Set default x-axis padding
	mpl.rcParams['axes.ymargin'] = 0						#Set default y-axis padding
	mpl.rcParams['errorbar.capsize'] = 3					#Set error bar end cap width
	mpl.rcParams['font.size'] = 12							#Set global fontsize
	mpl.rcParams['legend.fontsize'] = 'large'				#Set legend fontsize
	mpl.rcParams['figure.titlesize'] = 'medium'				#Set title fontsize

	#Line and Colour options
#	from cycler import cycler								#See below
	mpl.rcParams['axes.prop_cycle']=cycler(color='krbgcym')	#Set Default colour rotation
	mpl.rcParams['lines.linewidth'] = 1.0					#Set Default linewidth

	#Maths and Font options
	mpl.rcParams['mathtext.fontset'] = 'cm'					#Sets 'Latex-like' maths font
	mpl.rcParams['mathtext.rm'] = 'serif'					#Sets default string font

	return()
#enddef
Matplotlib_GlobalOptions()									#Must be run before diagnostics

#=========================#
#=========================#

def figure(subplots=[1,1],gridspec=[],aspectratio=[],shareX=False,shareY=False):
#Create figure and axes with variable aspect ratio, sub-plots and configurations.
#Takes image aspect ratio [x,y], number of subplots [rows, columns] and row/column sharing boolians
#Returns figure and axes seperately.
#Example: fig,ax = figure([2,1],[1,1],image_aspectratio,shareX=False,shareY=False)

	#If integer subplot supplied, convert to list.
	if isinstance(subplots,int) == True: 
		subplots = [subplots,subplots]
	#endif

	#Extract row/column values for easier reading
	XAspect,YAspect = aspectratio[0],aspectratio[1]
	nRows,nCols = subplots[0],subplots[1]

	#Create figure and axis
	if len(aspectratio) == 2:
		fig, ax = plt.subplots(nrows=nRows,ncols=nCols,figsize=(XAspect,YAspect),sharex=shareX,sharey=shareY)
	else:
		fig, ax = plt.subplots(nrows=nRows,ncols=nCols,figsize=(10,10),sharex=shareX,sharey=shareY)
	#endif

	#if gridspec is supplied, set relative panel heights accordingly
	#Panel heights are defined relative to gridspec index [1]
	if len(gridspec) > 0:
		GridSpecArray = gs.GridSpec(subplots[0],subplots[1], height_ratios=gridspec)
		for i in range(0,len(gridspec)): ax[i] = plt.subplot(GridSpecArray[i])
		#endfor
	#endif

	return(fig,ax)
#enddef

#=========================#
#=========================#

def ImageOptions(fig,ax='NaN',Xlabel='',Ylabel='',Title='',Legend=[]):
#Applies plt.options to current figure based on user input.
#Returns nothing, open figure is required, use figure().
#For best results call immediately before saving/displaying figure.
#Example: ImageOptions(fig,plt.gca(),Xlabel,Ylabel,Title,Legend)

	#If no axis is supplied, use current open axis
	if ax == 'NaN': ax = plt.gca()

	#Apply user overrides to plots.
	if len(titleoverride) > 0:
		Title = titleoverride
	if len(legendoverride) > 0:
		Legend = legendoverride
	if len(xlabeloverride) > 0:
		Xlabel = xlabeloverride[0]
	if len(ylabeloverride) > 0:
		Ylabel = ylabeloverride[0]
	#endif

	#Set title and legend if one is supplied.
	if len(Title) > 0:
		if '\n' in Title: ax.set_title(Title, fontsize=14, y=1.03)
		else: ax.set_title(Title, fontsize=18, y=1.03)
	if len(Legend) > 0:
		ax.legend(Legend, fontsize=16, frameon=False)
	#endif

	#Set labels and ticksize.
	ax.set_xlabel(Xlabel, fontsize=24)
	ax.set_ylabel(Ylabel, fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	#Force scientific notation for all axes.
	try: ax.xaxis.get_major_locator().set_params(axis='both',style='sci',scilimits=(-2,3))
	except: Fails_If_Axis_Ticks_Contain_Strings = True
	#endtry

	#Set grid, default is off.
	if image_plotgrid == True: ax.grid(True)
	#endif

	#Plot mesh outline if requested.	### HACKY ###
	if image_plotmesh == True:
		mesh_auto_plot = 1 				# !!AUTO PLOT MESH NOT IMPLIMENTED!! #
	elif image_plotmesh == 'ASDEX' and Crop == True:	
		ManualASDEXMesh(ax)
	#endif

	#Arrange figure such that labels, legends and titles fit within frame.
	try: fig.tight_layout()
	except: Fails_For_Incorrect_Padding = True

	return()
#enddef

#=========================#
#=========================#

def Colourbar(ax='NaN',image='NaN',Label='',Ticks=5,Lim=[]):
#Creates and plots a colourbar with given label and binsize.
#Takes image axis, label string, number of ticks and limits
#Allows pre-defined colourbar limits in form [min,max].
#Returns cbar axis if further changes are required.
#Example: cbar = Colourbar(ax[0],im,'Label',5,Lim=[0,1])

	#Determine supplied image and axis, replacing or breaking if required
	if image == 'NaN': 
		try: image = im
		except:	print('\n Colourbar Image Not Supplied \n')
	if ax == 'NaN': 
		try: ax = plt.gca()
		except: print('\n Colourbar Axes Not Supplied \n')
	#endif

	#Apply label override if supplied
	if len(cbaroverride) > 0: Label = str(cbaroverride)
	
	#Set default font and spacing options and modify if required
	Rotation,Labelpad = 270,30
	LabelFontSize,TickFontsize = 24,18
	if '\n' in Label: Labelpad += 30		#Pad label for multi-line names

	#Create and define colourbar axis
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	cbar = plt.colorbar(image, cax=cax)

	#Set number of ticks, label location and define scientific notation.
	cbar.set_label(Label, rotation=Rotation,labelpad=Labelpad,fontsize=LabelFontSize)
	cbar.formatter.set_scientific(True)
	cbar.formatter.set_powerlimits((-2,3))
	cbar.locator = ticker.MaxNLocator(nBins=Ticks)
	cbar.ax.yaxis.offsetText.set(size=TickFontsize)
	yticks(fontsize=TickFontsize)
	cbar.update_ticks()  

	#Apply colourbar limits if specified.  (lim=[min,max])
	if len(Lim) == 2: im.set_clim(vmin=Lim[0], vmax=Lim[1])

	return(cbar)
#enddef

#=========================#
#=========================#

def InvisibleColourbar(ax='NaN'):
#Creates an invisible colourbar to align subplots without colourbars.
#Takes image axis, returns colourbar axis if further edits are required
#Example: cax = InvisibleColourbar(ax[0])

	#Determine supplied axis, replacing or breaking if required
	if ax == 'NaN': 
		try: ax = plt.gca()
		except: print('\n Colourbar Axes Not Supplied \n')
	#endif

	#Create colourbar axis, ideally should 'find' values of existing cbar! 
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.1)

	#Set new cax to zero size and remove ticks.
	try: cax.set_facecolor('none')				#matplotlib v2.x.x method
	except: cax.set_axis_bgcolor('none')		#matplotlib v1.x.x method
	for axis in ['top','bottom','left','right']:
		cax.spines[axis].set_linewidth(0)
		cax.set_xticks([])
		cax.set_yticks([])
	#endfor

	return(cax)
#enddef

#=========================#
#=========================#

def VariableLabelMaker(variables,Units=[]):
#Makeshift way of creating units for each legend entry.
#Example: VariableLegends = VariableLabelMaker(variables,Units)

	#Convert to single element list if string is supplied
	if type(variables) is not list:
		variables = [variables]
	#endif

	VariableLegends = list()
	for i in range(0,len(variables)):
		#Explicit Control Parameters
		if variables[i] == 'kst':
			Variable = 'kstep'
			VariableUnit = '[iter]'
		elif variables[i] == 't':
			Variable = 'time'
			VariableUnit = '[ms]'
		
		#Explicit Axes
		elif variables[i] == 'r_psi':
			Variable = 'Radial flux surface $r_{\psi}$'					# Check this?   sqrt(r_psi) = rho_pol (I think...)
			VariableUnit = '[-]'
		elif variables[i] == 'gpsi_nrm':
			Variable = '$g_{psi}$ norm'									# UNKNOWN VARIABLE
			VariableUnit = '[-]'
		elif variables[i] == 'q_psi':
			Variable = 'Safety Factor $q_{\psi}$'
			VariableUnit = '[-]'

		#Explicit MHD fieldss
		elif variables[i] == 'brad':
			Variable = 'Radial B-field Perturbation $\delta B_{r}$'
			VariableUnit = '[mT]'
		elif variables[i] == 'btheta':
			Variable = 'Poloidal B-field Perturbation $\delta B_{\\theta}$'
			VariableUnit = '[mT]'
		elif variables[i] == 'bphi':
			Variable = 'Toroidal B-field Perturbation $\delta B_{\phi}$'
			VariableUnit = '[mT]'
		elif variables[i] == 'erad':
			Variable = 'Radial E-field Perturbation $\delta E_{r}$'
			VariableUnit = '[V m$^{-1}$]'
		elif variables[i] == 'etheta':
			Variable = 'Poloidal E-field Perturbation $\delta E_{\\theta}$'
			VariableUnit = '[V m$^{-1}$]'
		elif variables[i] == 'ephi':
			Variable = 'Toroidal E-field Perturbation $\delta E_{\phi}$'
			VariableUnit = '[V m$^{-1}$]'
		#endif

		#Explicit MHD Velocities
		elif variables[i] == 'vrad':
			Variable = 'Radial Velocity Perturbation $\delta v_{r}$'				#Momentum?
			VariableUnit = '[m s$^{-1}$]'
		elif variables[i] == 'vtheta':
			Variable = 'Poloidal Velocity Perturbation $\delta v_{\\theta}$'		#Momentum?
			VariableUnit = '[m s$^{-1}$]'
		elif variables[i] == 'vphi':
			Variable = 'Toroidal Velocity Perturbation $\delta v_{\phi}$'			#Momentum?
			VariableUnit = '[m s$^{-1}$]'

		#Explicit MHD Densities and Pressures
		elif variables[i] == 'rho':
			Variable = 'Plasma Density Perturbation $\delta \n_{e}$'
			VariableUnit = '[m$^{-3}$]'
		elif variables[i] == 'prs':
			Variable = 'Pressure Perturbation $\delta P_{rs}$'
			VariableUnit = '[Pa]'

		#Explicit Fast Particle Momentum
		elif variables[i] == 'mom_a':
			Variable = 'Fast Ion Momentum p$_{FI}$'
			VariableUnit = '[kg m s$^{-1}$]'
		elif variables[i] == 'dns_a':
			Variable = 'Fast Ion density n$_{FI}$'
			VariableUnit = '[-]'
		elif variables[i] == 'ppara_a':
			Variable = 'Fast ion para. momentum ppara$_{FI}$'
			VariableUnit = '[-]'
		elif variables[i] == 'pperp_a':
			Variable = 'Fast Ion perp. momentum pperp$_{FI}$'
			VariableUnit = '[-]'
		elif variables[i] == 'qpara_a':
			Variable = 'Fast Ion para. charge qpara$_{FI}$'
			VariableUnit = '[-]'
		elif variables[i] == 'qperp_a':
			Variable = 'Fast Ion perp. charge qperp$_{FI}$'
			VariableUnit = '[-]'

		#Default if no fitting variable found.
		else:
			Variable = 'Variable'
			VariableUnit = '[-]'
		#endif

		#Append variable and unit if required
		if len(Units) == 0:
			VariableLegends.append(Variable+' '+VariableUnit)
		else:
			VariableLegends.append(Variable+' '+str(Units))
		#endif
	#endfor

	#Reduce list to string if it only contains one element
	if len(VariableLegends) == 1:
		VariableLegends = VariableLegends[0]
	#endif

	#=====#=====#

	#DebugMode: Print VariableLegends
	if DebugMode == 'True':
		print('')
		print(variables)
		print(VariableLegends)
	#endif
	
	#=====#=====#

	return(VariableLegends)
#enddef

#====================================================================#
#====================================================================#






#====================================================================#
				  #COMMON DATA ANALYSIS FUNCTIONS#
#====================================================================#

def VectorDerivative(XArray,YArray,Order=1,Trend='lin'):
#Determines MHD toroidal mode linear growth rates through analysis of 1st and 2nd energy derivatives
#Derivatives are taken of log(E)/dt such that the 1st and 2nd derivatives are linear (flat)
#Solves: Eend = Estart*exp{gamma*dt} over time indices where 2nd derivative is close to zero
#Inputs: 
#	XArray - 1D array of variables over which to provide differentiation (typically time)
#	YArray - 1D array of variables to be differentiated (typically energy, fields etc...)
#	Order = 0D scalar determining the order of derivative (Only defined for Order > 0)
#	Treand = String determining function shape: 'lin', 'exp',
#Outputs:
#	DxDyArray - 1D array of containing n'th order derivative D^nx/Dy^n  :: (x2-x1)/(y2-y1) for n=1
# 	DxArray - 1D array containing nth derivative of XArray, i.e. difference between successive indices
# 	DyArray - 1D array containing nth derivative of YArray, i.e. difference between successive indices
#Warnings: 
#	Function only works for 1D Arrays
#Example: DxDyArray = VectorDerivative(TimeArray,EnergyArray,Order=1 )[0]

	#Compute i'th derivative of supplied arrays
	for i in range(0,Order):

		#Calculate derivative arrays - i.e. compute difference between successive indices
		DxArray = np.diff(XArray).tolist()
		DyArray = np.diff(YArray).tolist()
#		print len(DxArray)					#Check length of Arrays for debugging higher orders

		#Calculate gradient array - i.e. derivative of Y to X
		DxDyArray = [DyArray[j]/DxArray[j] for j in range(0,len(DxArray))]

		#Sum derivatives up to each index in SEQuence to reconstruct XArray and set YArray to DxDy
		#Only required for calculation of higher order derivatives (Order > 1) - ARRAYS NOT RETURNED
		XArray = [sum(DxArray[0:j]).tolist() for j in range(0,len(DxArray))]
		if Trend == 'exp': 		YArray = np.log(DxDyArray)
		elif Trend == 'lin':	YArray = DxDyArray
	#endfor

	return(DxDyArray,DxArray,DyArray)
#enddef

#=========================#
#=========================#

def ComputeTAEThresholds(HarmonicData,Harmonic,eps,Va,ax='NaN'):
#Compute upper and lower Alfven eigenmode threshold frequencies
#UpperThreshold,LowerThreshold = ComputeTAEThresholds(Harmonic,mpol,lpsi,qpsi,rho_pol,eps,AlfvenVelocity,subfig)

	#eps = ????				[-]
	#Va = Alfven Velocity 	[m/s]

	#Extract required data
	data = HarmonicsData.data
	kmax = data.shape[0]					#Maximum kstep of data
	mpol = data.shape[1]					#Number of poloidal modes
	ntor = data.shape[2]					#Number of toroidal modes
	lpsi = data.shape[3]					#Radial Array Dimension [lpsi]

	#Initiate TAE threshold arrays
	UpperThresholds = list()
	LowerThresholds = np.zeros([lpsi,mpol-1])

	#Extract rho_pol and safety factor arrays, and initiate empty threshold arrays
	rho_pol = HarmonicsData.rho_pol
	q_psi = abs(HarmonicsData.q_psi)		#Safety Factor radial profile [lpsi]
	K = np.zeros([lpsi,mpol])				#Initiate empty 'K' array

	#Calculates the Alfven eigenmode thresholds for all simulated poloidal mode numbers
	#PROVIDE SUMMARY OF MATHS AND ASSUMPTIONS 
	#PROVIDE REFERENCE FOR THESE DERIVATION(S)
	for m in range(0,mpol):
		K[:,m] = (m-abs(Harmonic)*q_psi)/(q_psi*R0)
	#endfor
	for m in range(0,mpol-1,1):
		diff  = np.sqrt(((K[:,m]**2-K[:,m+1]**2)**2+4*eps**2*K[:,m]**2*K[:,m+1]**2)) 
		UpperThresholds.append( np.sqrt((K[:,m]**2+K[:,m+1]**2+diff)/(2*(1-eps**2))) )
		LowerThresholds[:,m] = np.sqrt((K[:,m]**2+K[:,m+1]**2-diff)/(2*(1-eps**2)))
	#endfor

	#Plot TAE Thresholds if axis is supplied
	if ax != 'NaN':
		for m in range(0,mpol-1,1):
			Yaxis = UpperThresholds[m]*Va/(2*np.pi)/1000
			ax.plot(rho_pol,Yaxis, 'w--', lw=1.5)
		#endfor
		Yaxis = np.amin(LowerThresholds,axis=1)*Va/(2*np.pi)/1000
		subfig.plot(rho_pol,Yaxis, 'w--', lw=1.5)
	else:
		print('Image Axis Not Supplied - TAE Thresholds Not Plotted')
	#endif

	return(UpperThresholds,LowerThresholds)
#endfor

#=========================#
#=========================#

def ComputeMHDGrowthRates(EnergyArray,TimeArray,Threshold=100):
#Determines MHD toroidal mode linear growth rates through analysis of 1st and 2nd energy derivatives
#Derivatives are taken of log(E)/dt such that the 1st and 2nd derivatives are linear (flat)
#Solves: Eend = Estart*exp{gamma*dt} over time indices where 2nd derivative is close to zero
#Inputs: 
#	EnergyArray - 1D array of temporally resolved MHD energies for a single mode number
#	TimeArray - 1D array of SI times [s] relating to MHD energies provided in EnergyArray
#	Threshold - 0D float determining the 'sensitivity' of the function, higher values are more sensitive.
#				Specifically, Threshold sets the maximum 'distance from zero' when creating LinearRegion
#Outputs:
#	gamma - 0D scalar 'linear' growth rate [s-1] for supplied MHD mode number
# 	Delta1Energy - 1D array containing 1st derivative of provided EnergyArray to TimeArray
# 	Delta2Energy - 1D array containing 2nd derivative of provided EnergyArray to TimeArray
#Example: gamma, dEdt, d2Edt2 = ComputeMHDGrowthRates(Energy_n[ntor],TimeArray)

	#Use log(Energy) for all calculations
	LogEnergyArray = np.log(EnergyArray)

	#Compute 1st derivative of energy:		d log(E) / dt
	Delta1Energy = list()
	Delta1Energy.append(TimeArray[0:-1])											#Add time array
	Delta1Energy.append( VectorDerivative(TimeArray,LogEnergyArray,1,'exp' )[0] )	#Add n'th harmonic array
	
	#Compute 2nd derivative of energy:		d^2 log(E) / dt^2
	Delta2Energy = list()
	Delta2Energy.append(TimeArray[0:-2])											#Add time array
	Delta2Energy.append( VectorDerivative(TimeArray,LogEnergyArray,2,'exp' )[0] )	#Add n'th harmonic array

	#==========##==========#

	#Smoothing the 2nd derivative array to remove any kinetic noise (Savitzk-Golay filter)
	#Smoothing helps determination of LinearRegion threshold values for non-constant growth rates
	if KineticFiltering == True:
		WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
		Delta2Energy_Smooth = (savgol_filter(Delta2Energy, WindowSize, PolyOrder)).tolist()
	#endif

	#Plot kinetic smoothing comparison for debugging purposes if required
	if DebugMode == True:
		plt.plot(Delta2Energy[0],Delta2Energy[1], 'k--', lw=1)				#Unsmoothed
		plt.plot(Delta2Energy[0],Delta2Energy_Smooth[1], 'r-', lw=2)		#Smoothed
		plt.legend(['Unsmoothed','Smoothed'])
		plt.xlabel('Time [s]')
		plt.ylabel('Energy [-]')
		plt.show()
	#endif

	#==========##==========#

	#Determine temporal extent of linear growth region, i.e. where 2nd derivative is close to zero
	#Threshold is the maximum allowed distance from zero, Threshold is set by function input.
	LinearRegion = list()
	for i in range(0,len(Delta2Energy_Smooth[1])):
		if abs(Delta2Energy_Smooth[1][i]) < Threshold: 	LinearRegion.append(1)
		else: 											LinearRegion.append(0)
	#endfor

	#Smooth Threshold to remove most of the kinetic noise (Savitzk-Golay filter)
	WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
	LinearRegion = (savgol_filter(LinearRegion, WindowSize, PolyOrder)).tolist()
	#endif

	#Create 'clean' binary mask for linear growth region, default threshold set to 0.5
	for i in range(0,len(LinearRegion)):
		if LinearRegion[i] > 0.5: LinearRegion[i] = 1.0
		else: LinearRegion[i] = 0.0
	#endfor

	#Compute 'linear' phase growth rate (gamma [s-1]) over full linear phase.
	#Assumes exponential growth where: Eend = Estart*exp{gamma*dt}
	try:
		#Determine linear phase start/end indices and times
		StartIdx = LinearRegion.index(1)
		EndIdx = len(LinearRegion) - LinearRegion[::-1].index(1)
		tstart = TimeArray[StartIdx]				#[s]
		tend = TimeArray[EndIdx]					#[s]
		dt = tend-tstart							#[s]		#Can't be zero

		#Determine linear phase start/end energies
		Estart = EnergyArray[StartIdx]				#[-]		#Can't be zero
		Eend = EnergyArray[EndIdx]					#[-]

		#Compute growth rate: gamma = ln(Eend/Estart)/dt
		gamma = np.log(Eend/Estart)/dt				#[s-1]
		gamma = round(gamma,2)

#	 	THE ABOVE METHOD ISN'T GREAR AS IT ASSUMES A CONSTANT GROWTH RATE OVER THE LINEAR REGION
#		BETTER METHOD IS TO AVERAGE ALL THE `1st DERIVATIVES` USING LinearRegion AS A MASK
#		PROBLEM IS THAT Delta1Energy = (E2-E1) / dt, while growth rate = (E2/E1) / dt  (i.e. diff vs ratio)
#		NEED TO FIX THIS LINE... BUT OTHERWISE THIS SHOULD BE A MORE RELIABLE METHOD.
		Delta1Energy_Masked = ma.masked_array( Delta1Energy, mask=LinearRegion.append(0) )
##		print np.nanmean(Delta1Energy_Masked)
##		print gamma
##		exit()

	#If no linear phase is found, growth rate gamma is set to np.nan
	except:
		tstart = np.nan; tend = np.nan
		Estart = np.nan; Eend = np.nan
		gamma = np.nan
	#endtry

	#Print debug outputs to terminal if required
	if DebugMode == True:
		print('')
		print( round(tstart,3), round(tend,3) )
		print( round(Estart,3), round(Eend,3) )
		print( round(gamma,2), '[s-1]')
		print('')
	#endif

	return(gamma,Delta1Energy,Delta2Energy)
#enddef

#=========================#
#=========================#

def ComputeEquilLCFS(Equilibrium,Threshold=0.0):
#Determines which cells of the equilibrium are within the LCFS and which are outside
#Those inside are unchanged, while values outside are replaced with np.nan()
#Inputs are 2D equilibrium array and LCFS phi threshold (default 0.0)
#Returns 2D EquilibriumLCFS array containing $\Phi(R,Z)$ shaped as: Equilibrium[Row][Column]

	#Initiate required lists
	LCFSEQuil = list()

	#By definition LCFS occours where flux surface equals zero
	for i in range(0,len(Equilibrium)):
		LCFSEQuil.append(list())
		for j in range(0,len(Equilibrium[i])):
			if Equilibrium[i][j] >= Threshold:
				LCFSEQuil[i].append(Equilibrium[i][j])
			else:
				LCFSEQuil[i].append(np.nan)
			#endif
		#endfor
	#endfor

	return(LCFSEQuil)
#enddef

#=========================#
#=========================#

def Normalise(profile,NormFactor=0):
#Takes 1D or 2D array and returns array normalised to maximum value.
#If NormFactor is defined, array will be normalised to this instead.
#Returns normalised image/profile and the max/min normalisation factors.
#Example: NormProfile,Min,Max = Normalise(profile,NormFactor=0)

	#Initiate any required output lists
	NormalisedImage = list()

	#determine dimensionality of profile and select normaliztion method.
	if isinstance(profile[0], (list, np.ndarray) ) == True:

		#Obtain max and min normalization factors for 2D array.
		FlatImage = [item for sublist in profile for item in sublist]
		MaxNormFactor,MinNormFactor = max(FlatImage),min(FlatImage)

		#Fix for division by zero and other infinity related things...
		if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
		if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0
		#endif

		#Normalize 2D array to local maximum.
		if NormFactor == 0: NormFactor = MaxNormFactor
		for i in range(0,len(profile)):
			NormalisedImage.append( [x/NormFactor for x in profile[i]] )
		#endfor
		profile = NormalisedImage
		return(profile,MaxNormFactor,MinNormFactor)

	#Lowest dimention is still list.
	elif isinstance(profile, (list, np.ndarray) ) == True:

		#Obtain max and min normalization factors for 1D profile.
		MaxNormFactor,MinNormFactor = max(profile),min(profile)

		#Fix for division by zero and other infinity related things...
		if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
		if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0

		#Normalize 1D array to local maximum.
		if NormFactor == 0: NormFactor = MaxNormFactor
		for i in range(0,len(profile)):
			profile[i] = profile[i]/NormFactor
		#endfor
	#endif

	return(profile,MinNormFactor,MaxNormFactor)
#enddef

#====================================================================#
#====================================================================#

























#====================================================================#
					#WELCOME TEXT AND INFORMATION#
#====================================================================#

print(']')
print('-------------------------------------------------------')
print(' .___  ___.      ___   ____    ____  __       _______. ')
print(' |   \/   |     /   \  \   \  /   / |  |     /       | ')
print(' |  \  /  |    /  ^  \  \   \/   /  |  |    |   (----` ')
print(' |  |\/|  |   /  /_\  \  \      /   |  |     \   \     ')
print(' |  |  |  |  /  _____  \  \    /    |  | .----)   |    ')
print(' |__|  |__| /__/     \__\  \__/     |__| |_______/     ')
print('                                                 v0.6.7')
print('-------------------------------------------------------')
print('')
print('The following diagnostics were requested:')
print('-----------------------------------------')
if True in [savefig_1Denergy,savefig_1Denergytrends]:
	print('# Energy Convergence Analysis')
if True in [savefig_1Dequilibrium]:
	print('# 1D Equilibrium Analysis')
if True in [savefig_2Dequilibrium,savefig_2Dequilmovie]:
	print('# 2D Equilibrium Analysis')
if True in [savefig_2Dcontinuum]:
	print('# 2D Continuum Analysis')
if True in [savefig_2Dpolspectrum]:
	print('# 2D Spectral Analysis')
if True in [savefig_1Dkinetics,savefig_2Dkinetics]:
	print('# Kinetic Distribution Analysis')
print('-----------------------------------------')
print('')

#====================================================================#
#====================================================================#


#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Create lists for basic processing
Dir = list()			#List of simulation folders 					struc:[folder]
DirFiles = list()		#List of data files in each simulation folder 	struc:[folder][filenames]
NumFolders = 0			#Number of simulation folders

#Mesh sizes and axes
Raxis = list()				#raxis		[m]
Zaxis = list()				#zaxis		[m]
PhiMode = list()			#phimax		[Rad]
RGeo = list()				#major_r	[m]
InboardLCFS = list()		#left 		[m]
OutputLCFS = list()			#right 		[m]

R_mesh = list()
Z_mesh = list()
Phi_mesh = list()

Depth = list()
Radius = list()
Height = list()
dr = list()
dz = list()
dphi = list()

#Lists to store normalisation factors
Variables = list()				#[1D Array] of normalisation factor variable names	- Strings
Values = list()					#[1D Array] of normalisation factor variable values - Floats
Units = list()					#[1D Array] of normalisation factor variable units - Strings

#Lists to store extracted data
HarmonicsData = list()			#[4D Array] of shape Data[kstep][mpol][ntor][lpsi][???] for each variable
MomentsData = list()			#[4D Array] of shape Data[kstep][mpol][ntor][lpsi][???] for each variable
KineticsData = list()			#[2D Array] of shape Data[variables][markers(n)] concatinated for all nodes
EnergyData_phys = list()		#[3D Array] of shape Data[folder][variable][Kstep] for energy_phys.txt
EnergyData_n = list()			#[3D Array] of shape Data[folder][variable][Kstep] for energy_n.txt

#====================================================================#
#====================================================================#


#====================================================================#
						#OBTAIN FILE DIRECTORIES#
#====================================================================#

#Obtain system RAM. (and rename enviroment variable)
mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
mem_gib = mem_bytes/(1024.**3)
ext = image_extension

#Obtain home directory (location of MAVIS) and contents of said directory
Root = os.path.abspath(".")
HomeDirFolders,HomeDirContents = DirectoryContents(Root)

#For each sub-directory in HomeDirFolders:
for l in range(0,len(HomeDirFolders)):

	#Obtain sub-directory names and contents (Simulation folder names and contents)
	SubDirFolders,SubDirContents = DirectoryContents(Root+HomeDirFolders[l])

	#Determine which sub-direcotires contain a '/data/' folder (i.e. MEGA simulation folders)
	if '/data/' in SubDirFolders:
		#Add folder to global simulation list
		Dir.append(Root+HomeDirFolders[l])
		DirFiles.append(list())
		NumFolders += 1
	#Discard folder if it doesn't contain data
	else:
		#Print debug outputs to terminal if requested
		if DebugMode == True:
			print 'Discarding Directory: ', HomeDirFolders[l]
		#endif
	#endif
#endfor
Dir = sorted(Dir)			#Sort MEGA simulation directories into alphanumerical order


for l in range(0,len(Dir)):
	#Extract contents from 'l'th simulation folder and data/ subfolder
	SimDirContents = DirectoryContents(Dir[l])[1]		#Documents in 'Simulation' Folder
	DataDir = '/'+Dir[l]+'data/'						#'data' folder: Root+Dir[l]+'data/'
	DataDirContents = DirectoryContents(DataDir)[1]		#Documents in 'data' Folder

	#Save content files from simulation folder that fit requested data output file extensions
	for j in range(0,len(SimDirContents)):
		Filename = SimDirContents[j]
		if any([x in Filename for x in FileExtensions]):
			Prefix = Dir[l]
			DirFiles[l].append(Prefix+Filename)
			#endif
		#endif
	#endfor

	#Save content files from /data/ subfolder that fit requested data output file extensions
	for j in range(0,len(DataDirContents)):
		Filename = DataDirContents[j]
		if any([x in Filename for x in FileExtensions]):
			Prefix = Dir[l]+'data/'						#Note: Dir ends with a '/'
			DirFiles[l].append(Prefix+Filename)
			#endif
		#endif
	#endfor
#endfor

#If no folders detected end analysis script; else continue to analysis.
if NumFolders > 0:
	print '------------------------------------------'
	print 'Initial Setup Complete, Starting Analysis:'
	print '------------------------------------------'
elif NumFolders == 0:
	print '-------------------------------------------'
	print 'No Ouput Files Detected, Aborting Analysis.'
	print '-------------------------------------------'
	print ''
	exit()
#endif

#=====================================================================#
#=====================================================================#































#====================================================================#
				  #ENERGY & CONVERGENCE DIAGNOSTICS#
#====================================================================#

#====================================================================#
				  	#SPECTRAL ENERGY CONVERGENCE#
#====================================================================#

if savefig_1Denergy == True:

	#For each detected simulation folder
	for l in tqdm(range(0,len(Dir))):

		#Create global 1D diagnostics folder and extract current simulation name
		DirEnergy = CreateNewFolder(Dir[l],'1DEnergy_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Kstep [-], Time [ms] & toroidal harmonics from energy_n.txt
		SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[l], DataFile='energy_n')
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
		ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
		ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
		ntor0 = ntorArray[0]						#ntor = 0, baseline equilibrium data

		#Extract Energy_n outputs and header for plotting
		#energy_n: [ntor][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		Energy_n = Energy_n[2::]					#Remove KStep and Time arrays from array

		#Extract Energy_Phys outputs and header for plotting
		#Energy_Phys: [variable][timestep]
		Energy_Phys,Header_Phys = ExtractMEGA_Energy(Dir[l],'energy_phys')

		#Extract normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
#		print Variables[1],Values[1],Units[1]

		#Compute 1st and 2nd energy derivatives and determine MHD linear growth rates
		#Solves: Eend = Estart*exp{gamma*dt} where 2nd derivative is close to zero
		gamma_Array = list()
		dEnergydt_Array, d2Energydt2_Array = list(),list()
		for i in range(0,len(Energy_n)):
			gamma,dEdt,d2Edt2 = ComputeMHDGrowthRates(Energy_n[i],TimeArray)
			gamma_Array.append(gamma)
			dEnergydt_Array.append(dEdt)
			d2Energydt2_Array.append(d2Edt2)
		#endfor

		#==========##==========#
		#==========##==========#

		#Create figure for energy_n outputs
		fig,ax = figure(subplots=[2,1], aspectratio=image_aspectratio)

		#Energy_n Ax[0] Title, Legend, Axis Labels etc...
		Title = 'Spectrally Resolved Energy Evolution for \n '+DirString
		Xlabel,Ylabel = '', 'Energy $\epsilon_{n}$ (Log$_{10}$) [-]'
		Legend = list()

		#Plot total energy for each harmonic component
		for i in range(0,len(Energy_n)):
			ax[0].plot(TimeArray,np.log10(Energy_n[i]), lw=2)
			Legend.append( '$\gamma'+'_{'+str(i)+'}$ = '+str(gamma_Array[i])+' [s$^{-1}$]' )
		#endfor
		ImageOptions(fig,ax[0],Xlabel,Ylabel,Title,Legend)

		#Energy_n Ax[1] Title, Legend, Axis Labels etc...
		Title = 'Spectrally Resolved Energy Evolution for \n '+DirString
		Xlabel,Ylabel = 'Time [ms]', '$\Delta$ Energy $\\frac{d \epsilon_{n}}{d t}$ (Log$_{10}$) [-]'
		Legend = list()

		#Plot 1st derivative of energy for each harmonic component
		for i in range(0,len(Energy_n)):
			ax[1].plot(dEnergydt_Array[i][0],np.log10(dEnergydt_Array[i][1]), lw=2)
			Legend.append( 'n$_{tor}$ = '+str(i) )
		#endfor
		ImageOptions(fig,ax[1],Xlabel,Ylabel,'',Legend)

		#Save and close open figure(s)
		plt.savefig(DirEnergy+'SpectralEnergy_'+SubString+ext)
#		plt.show()
		plt.close('all')

		#==========##==========#
		#==========##==========#

		#Create figure for energy_phys outputs
		fig,ax = figure(subplots=[3,1], aspectratio=image_aspectratio)

		#Energy_phys[0,1,2] Title, Legend, Axis Labels etc...
		Title = 'Spectrally Integrated Energy Evolution for \n '+DirString
		Xlabel,Ylabel = 'Time [ms]', 'Energy [-]'

		#Plot total thermal, kinetic and magnetic MHD (fluid solver) energy over time
		ax[0].plot(Energy_Phys[1],Energy_Phys[2],'k-',lw=2)		#Kinetic
		ax[0].plot(Energy_Phys[1],Energy_Phys[3],'r-',lw=2)		#Magnetic
		ax[0].plot(Energy_Phys[1],Energy_Phys[4],'b-',lw=2)		#Thermal
		Legend = ['Kinetic','Magnetic','Thermal']				#Header_Phys[2::]
		ImageOptions(fig,ax[0],'','MHD '+Ylabel,Title,Legend)

		#Plot parallel, perpendicular and total fast ion (kinetic solver) energy over time
		ax[1].plot(Energy_Phys[1],Energy_Phys[5],'k-',lw=2)				#Energy parallel to current flow
		ax[1].plot(Energy_Phys[1],Energy_Phys[6],'r-',lw=2)				#Energy perpendicular to current flow
		ax[1].plot(Energy_Phys[1],Energy_Phys[7],'b-',lw=2)				#Total Energy (only for df, not full-f)
		Legend = ['$\hat{J}$ Parallel','$\hat{J}$ Perpendicular','Total']	#Header_Phys[2::]
		ImageOptions(fig,ax[1],'','Fast Ion '+Ylabel,'',Legend)

		#Plot ?Transferred? and ?Total? energy over time
		ax[2].plot(Energy_Phys[1],Energy_Phys[8],'k-',lw=2)		#Transferred
		ax[2].plot(Energy_Phys[1],Energy_Phys[9],'r-',lw=2)		#Total
		Legend = ['Transferred','Total']						#Header_Phys[2::]
		ImageOptions(fig,ax[2],Xlabel,Ylabel,'',Legend)

		#Save and close open figure(s)
		plt.savefig(DirEnergy+'TotalEnergy_'+SubString+ext)
#		plt.show()
		plt.close('all')
	#endfor - Dir loop
#endif


#====================================================================#
				 	 	#SPECTRAL ENERGY TRENDS#
#====================================================================#

if savefig_1Denergytrends == True:

	#Create global 1D diagnostics folder and extract current simulation name
	DirTrends = CreateNewFolder(os.getcwd(),'/1D_Trends/')
	DirEnergy = CreateNewFolder(DirTrends,'/1DEnergy_Trends/')

	#Initiate any required lists
	dEnergydt_Array, d2Energydt2_Array = list(),list()
	gamma_Array = list()

	#Extract maximum shared KStep (index, not value) and associated folder Dir[Idx]
	MaxSharedKStep,MaxSharedDirIdx = ExtractMEGA_SharedDataRange(Dir)

	#Extract Kstep [-], Time [ms] & ntor arrays from energy_n.txt - use simulation with highest shared kstep range
	SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[MaxSharedDirIdx], DataFile='energy_n')
	DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
	DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
	KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
	ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
	ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
	ntor0 = ntorArray[0]						#ntor = 0, baseline equilibrium data

	#For each toroidal mode number
	for i in range(0,ntor_pos+1):

		#Create figure for energy_n trend comparison
		fig,ax = figure(subplots=[2,1], aspectratio=image_aspectratio)

		#Append new dimension to lists used when comparing growth rates
		dEnergydt_Array.append([]), d2Energydt2_Array.append([])
		gamma_Array.append([])

		#Refresh lists used when comparing energy profiles
		Legend0,Legend1 = list(),list()
		TrendAxis = list()

		#For each detected simulation folder
		for l in tqdm(range(0,len(Dir))):

			#Write simulation folder name strings
			DirString = Dir[l].split('/')[-2]
			SubString = DirString.split('_')[-1]
			TrendAxis.append(SubString)

			#Extract Energy_n outputs for all ntor in current folder
			#Energy_n: [ntor][timestep]
			Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
			#EnergyProfile: [timestep]
			EnergyProfile = Energy_n[i+2]					#Energy profile for ntor[i] (i+2, skips Kstep/time)
			EnergyProfile = EnergyProfile[0:MaxSharedKStep]	#Reduce KStep range to minimum shared range

			#Extract normalisation factors for current simulation folder
			Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
#			print Variables[1],Values[1],Units[1]

			#Compute 1st and 2nd energy derivatives and determine MHD linear growth rates
			#Solves: Eend = Estart*exp{gamma*dt} where 2nd derivative is close to zero
			gamma,dEdt,d2Edt2 = ComputeMHDGrowthRates(EnergyProfile,TimeArray)
			#dEnergydt_Array, d2Energydt2_Array 3D arrays of shape: [SimFolder][ntor][kstep]
			#gamma_Array 2D array of shape: [SimFolder][ntor]
			dEnergydt_Array[i].append([dEdt]); d2Energydt2_Array[i].append([d2Edt2])
			gamma_Array[i].append(gamma)

			#============##============#
			#Energy Profile Comparisons#
			#============##============#

			#Energy_n Ax[0] Title, Legend, Axis Labels etc...
			Title = 'MHD Energy Evolution for ntor = '+str(i)+' \n '+DirString
			Xlabel,Ylabel = '', 'Energy $\epsilon_{n}$ (Log$_{10}$) [-]'

			#Plot total energy for each folder (for ntor[i])
			ax[0].plot(TimeArray,np.log10(EnergyProfile), lw=2)
			Legend0.append( SubString+': $\gamma'+'_{'+str(i)+'}$ = '+str(gamma)+' [s$^{-1}$]' )
			#endfor
			ImageOptions(fig,ax[0],Xlabel,Ylabel,Title,Legend0)

			#Energy_n Ax[1] Title, Legend, Axis Labels etc...
			Xlabel,Ylabel = 'Time [ms]', '$\Delta$ Energy $\\frac{d \epsilon_{n}}{d t}$ (Log$_{10}$) [-]'

			#Plot 1st derivative of energy for each folder (for ntor[i])
			ax[1].plot(dEdt[0],np.log10(dEdt[1]), lw=2)
			Legend1.append( SubString )
			#endfor
			ImageOptions(fig,ax[1],Xlabel,Ylabel,'',Legend1)
		#endfor - Dir loop

		#Save and close 1D energy evolution profile figures for current ntor[i]
		plt.savefig(DirEnergy+'SpectralEnergy_n='+str(i)+'_Trends'+ext)
#		plt.show()
		plt.close('all')
	#endfor - ntor loop
	
	#==========##==========#
	#Growth Rate Comparison#
	#==========##==========#

	#Create figure for energy_n growth rate comparison
	fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)

	#Energy_n Ax[0] Title, Legend, Axis Labels etc...
	Title = 'Linear Growth Rate $\gamma$ Comparison for \n '+DirString
	Xlabel,Ylabel = 'Varied Parameter [-]', 'Growth Rate $\gamma_{n}$ [s$^{-1}$]'
	Legend = list()

	#Plot growth rates with respect to simulation folder for each toroidal mode number
	#gamma_Array 2D array of shape: [SimFolder][ntor]
	for i in range(0,len(gamma_Array)):
		ax.plot(TrendAxis, gamma_Array[i], marker='o', markerfacecolor='none', ms=14, lw=2)
		Legend.append('n$_{tor}$ = '+str(i))
	#endfor
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

	#Save and close open figure(s)
	plt.savefig(DirEnergy+'GrowthRate_Trends'+ext)
#	plt.show()
	plt.close('all')

	#endfor
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
#endif

#==========##==========##==========#
#==========##==========##==========#

if any([savefig_1Denergy,savefig_1Denergytrends]) == True:
	print '---------------------------'
	print '1D Energy Analysis Complete'
	print '---------------------------'
#endif

#====================================================================#
#====================================================================#

































#====================================================================#
					   #EQUILIBRIUM DIAGNOSTICS#
#====================================================================#

#====================================================================#
				 	   #1D EQUILIBRIUM PROFILES#
#====================================================================#

if savefig_1Dequilibrium == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - all need looped over... - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		SEQ = setting_SEQ[1]				#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!


		#Create global 1D diagnostics folder and extract current simulation name
		DirEquil1D = CreateNewFolder(Dir[l],'1DEquil_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Kstep [-] & Time [ms] arrays from SEQ.harmonics & toroidal harmonics from energy_n.txt
		SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[l], DataFile='harmonics')
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
		ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
		ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
		ntor0Idx = ntorArray[0]						#ntor = 0 index, contains (var0 + dvar) data

		#Extract toroidal mode number array index (ntorIdx) from requested mode number (ntor)
		ntorIdx = Set_ntorIdx(ntor,ntorArray)

		#Set Kstep ranges as requested - else default to max range
		KStepRange,KStepStep = Set_KStepRange(KStepArray,setting_kstep)
		KStepIdx = KStepRange[1]-1					#Requested KStep index	[-]

		#Set TimeIndex and employ to extract KStep and Time
		IdxOffset = SEQ*KStepMod					#[-]
		KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
		Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

		#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l],'All',ntor_tot,KStepIdx,SEQ,'3D')
		rho_pol = HarmonicsData.rho_pol				#Normalised radius		[-]

		#Extract data resolution and poloidal axes from repository .dat files
		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)
		DataShape = ExtractMEGA_DataShape(HarmonicsData)
		mpol_res = DataShape[0]; ntor_res = DataShape[1]
		lpsi_res = DataShape[2]; ltheta_res = DataShape[3]


		#For each requested variable
		for i in tqdm(range(0,len(variables))):

			#Create Variablelabel with units
			VariableLabel = VariableLabelMaker(variables[i])

			#==========##===========#
			#	 RADIAL PROFILES	#
			#==========##===========#
			if len(radialprofiles) > 0:

				#Create new folder to store radial profiles
				DirEquilRadial = CreateNewFolder(DirEquil1D,'Radial_Profiles/')

				#Create figure and define Title, Legend, Axis Labels etc...
				fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)
				ntorString = ', n='+str(ntor); mpolString=', m='+str(-mpol_res+1)+','+str(mpol_res-1)
				TimeString = ', t='+str(round(Time,3))+' ms'
				Title = VariableLabel+ntorString+mpolString+TimeString+' \n Simulation: '+DirString
				Xlabel,Ylabel = 'Radius $R$ [m]', VariableLabel
				Legend = list()

				#Plot 1D radially resolved profiles for current simulation folder
				#Radial profiles employ fixed poloidal (theta) and toroidal (phi) angles
				for j in range(0,len(radialprofiles)):

					#Define poloidal angle theta and append to legend list
					theta = radialprofiles[j]
					Legend.append('$\\theta$ = '+str(theta)+'$^{\circ}$')

					#Extract radially resolved profile and plot
					#RadialProfile has origin at Rgeo, extending at angle theta clockwise to vertical
					RadialProfile = Extract_RadialProfile(HarmonicsData,variables[i],ntorIdx,theta)
					ax.plot(rho_pol,RadialProfile, lw=2)

					#Save ASCII data to sub-folder
					if write_ASCII == True:
						#Create directory to hold ASCII data
						DirASCII = CreateNewFolder(DirEquilRadial,'1DEquil_Data/')

						#Set ASCII data file name string and header
						SaveString = variables[i]+'_n'+str(ntor)+'_theta='+str(theta)+'_t='+str(round(Time,3))+'.dat'
						Header = [VariableLabel,'   ', '@theta=',theta,'[Deg]', '   R',lpsi_res,  '\n']

						#Write 1D data header, then 2D PoloidalImage
						WriteFile_ASCII(Header, DirASCII+SaveString, 'w', 'RSV')
						WriteFile_ASCII(RadialProfile, DirASCII+SaveString, 'a', write_ASCIIFormat)
					#endif
				#endfor

				#Beautify 1D equilibrium profiles figure
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

				#Save radial equilibrium profiles for current simulation folder
				SaveString = variables[i]+'_Radial_n'+str(ntor)+'_t='+str(round(Time,3))+ext
				plt.savefig(DirEquilRadial+SaveString)
#				plt.show()
				plt.close('all')
			#end - radial profile branch

			#==========##===========#
			#	POLOIDAL PROFILES	#
			#==========##===========#
			if len(poloidalprofiles) > 0:

				#Create new folder to store poloidal profiles
				DirEquilPoloidal = CreateNewFolder(DirEquil1D,'Poloidal_Profiles/')

				#Create figure and define Title, Legend, Axis Labels etc...
				fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)
				ntorString = ', n='+str(ntor); mpolString=', m='+str(-mpol_res+1)+','+str(mpol_res-1)
				TimeString = ', t='+str(round(Time,3))+' ms'
				Title = VariableLabel+ntorString+mpolString+TimeString+' \n Simulation: '+DirString
				Xlabel,Ylabel = 'Poloidal Angle $\\theta$ [Deg]', VariableLabel
				Legend = list()

				#Plot 1D poloidally resolved profiles for current simulation folder
				#Poloidal profiles employ fixed radial (rho_phi) and toroidal (phi) angles
				for j in range(0,len(poloidalprofiles)):

					#Determine radial location from user supplied switchboard values
					#Radius is in normalised radius [rho_pol], while RadialLoc is in [m]
					Radius = poloidalprofiles[j]							#[rho_pol]
					RadialLoc = min(rho_pol, key=lambda x:abs(x-Radius))	#[m]
					RadialIdx = rho_pol.tolist().index(RadialLoc)			#[-]
					
					#Append variable name and SI radial location to legend
					Legend.append('$\\rho_{pol}$ = '+str(round(RadialLoc,2)))

					#Extract poloidally resolved profile and plot 
					#ThetaAxis and ThetaProfile rotate clockwise from vertical at R = Radius
					ThetaAxis,ThetaProfile = Extract_PoloidalProfile(HarmonicsData,variables[i],ntorIdx,Radius)
					ax.plot(ThetaAxis,ThetaProfile, lw=2)

					#Save ASCII data to sub-folder
					if write_ASCII == True:
						#Create directory to hold ASCII data
						DirASCII = CreateNewFolder(DirEquilPoloidal,'1DEquil_Data/')

						#Set ASCII data file name string and header
						SaveString = variables[i]+'_n'+str(ntor)+'_R='+str(Radius)+'_t='+str(round(Time,3))+'.dat'
						Header = [VariableLabel,'   ', '@R=',Radius, 'rho_pol', '   theta', ltheta_res,  '\n']

						#Write 1D data header, then 2D PoloidalImage
						WriteFile_ASCII(Header, DirASCII+SaveString, 'w', 'RSV')
						WriteFile_ASCII(ThetaProfile, DirASCII+SaveString, 'a', write_ASCIIFormat)
					#endif
				#endfor

				#Beautify 1D equilibrium profiles figure
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
				ax.xaxis.set_major_locator(ticker.MultipleLocator(60))
				ax.set_xlim(0,360)

				#Save poloidal equilibrium profiles for current simulation folder
				SaveString = variables[i]+'_Poloidal_n'+str(ntor)+'_t='+str(round(Time,3))+ext
				plt.savefig(DirEquilPoloidal+SaveString)
#				plt.show()
				plt.close('all')
			#endif - poloidal profile branch
		#endfor	- variable loop
	#endfor	- dir loop
#endif


#==========##==========##==========#
#==========##==========##==========#

if any([savefig_1Dequilibrium]) == True:
	print '--------------------------------'
	print '1D Equilibrium Analysis Complete'
	print '--------------------------------'
#endif

#====================================================================#
#====================================================================#


#====================================================================#
				 	      #2D POLOIDAL PLOTS#
#====================================================================#

if savefig_2Dequilibrium == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - all need looped over... - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		SEQ = setting_SEQ[1]				#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
#		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!

		#Create global 2D diagnostics folder and extract current simulation name
		DirEquil2D = CreateNewFolder(Dir[l],'2DEquil_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Kstep [-] & Time [ms] arrays from SEQ.harmonics & toroidal harmonics from energy_n.txt
		SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[l], DataFile='harmonics')
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
		ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
		ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
		ntor0Idx = ntorArray[0]						#ntor = 0 index, contains (var0 + dvar) data

		#Set Kstep ranges as requested - else default to max range
		KStepRange,KStepStep = Set_KStepRange(KStepArray,setting_kstep)
		KStepIdx = KStepRange[1]-1					#Requested KStep index	[-]

		#Set TimeIndex and employ to extract KStep and Time
		IdxOffset = SEQ*KStepMod					#[-]
		KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
		Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

		#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l],'All',ntor_tot,KStepIdx,SEQ,'3D')

		#Extract relevant spatial normalisation factors
		NormVariables,NormValues,NormUnits = ExtractMEGA_Normalisations(Dir[l])
		ZMin = NormValues[NormVariables.index('bottom_sim')]; ZMax = NormValues[NormVariables.index('top_sim')]
		Zgeo = NormValues[NormVariables.index('zaxis')]; Rgeo = NormValues[NormVariables.index('raxis')]

		#Extract data resolution and poloidal axes from repository .dat files
		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)
		Crdz = [Crdz[x]+ZMin for x in range(0,len(Crdz))]		#Offset vertical axis such that Z0 = Zgeo
		DataShape = ExtractMEGA_DataShape(HarmonicsData)
		mpol_res = DataShape[0]; ntor_res = DataShape[1]
		lpsi_res = DataShape[2]; ltheta_res = DataShape[3]


		#For each requested variable
		for i in tqdm(range(0,len(variables))):

			#Create Variablelabel with units
			variable = variables[i]
			VariableLabel = VariableLabelMaker(variable)

			#Create fig of desired size - increasing Xlim with the number of harmonics
			Xaspect, Yaspect = int(10*(float(ntor_tot)/1.75)), 10
			fig,ax = figure(subplots=[1,ntor_tot], aspectratio=[Xaspect,Yaspect])

			for j in range(0,ntor_tot):

				#Set toroidal mode number array index (ntorIdx) and mode number (ntor)
				ntor = -ntor_pos+j
				ntorIdx = j

				#Merge 3D Harmonics Data into 2D poloidal slice for variables[i]
				#PoloidalImage Shape: [lpsi][ltheta] ~~ [R][theta], like an onion (or ogre).
				#i.e. PoloidalImage[n][:] plots a full poloidal profile (clockwise from vertical) for R = Rho_pol[n]
				PoloidalImage = Extract_PoloidalImage(HarmonicsData,variable,ntorIdx)

				#Define Title, Legend, Axis Labels etc...
				SupTitle = VariableLabel+', n='+str(ntor)+', t='+str(round(Time,3))+' ms \n Simulation: '+DirString
				Xlabel,Ylabel = 'Major Radius $R$ [m]', 'Height $Z$ [m]'
				Legend = list()

				#Plot 2D poloidally resolved figure and beautify
				im = ax[ntorIdx].contourf(Crdr, Crdz, PoloidalImage, 100)#; plt.axis('scaled')
				im2 = ax[ntorIdx].contour(Crdr, Crdz, PoloidalImage, 20)#; plt.axis('scaled')

				#Beautify plots - taking account of panel location
				if ntorIdx == 0 and ntor_tot > 1: 						#If first panel with more panels to right
					cbar = Colourbar(ax[ntorIdx],im,'',5)
					ImageOptions(fig,ax[ntorIdx],Xlabel,Ylabel,'$n_{tor}$='+str(ntor),'')
				elif ntorIdx == 0 and ntor_tot == 1:					#If first panel with no panels to right
					cbar = Colourbar(ax[ntorIdx],im,'',5)
	 				ax[ntorIdx].axes.get_yaxis().set_visible(False)
					ImageOptions(fig,ax[ntorIdx],Xlabel,'','$n_{tor}$='+str(ntor),'')
				elif ntorIdx > 0 and ntorIdx < ntor_tot-1: 				#If middle panel with more panels to right
					cbar = Colourbar(ax[ntorIdx],im,'',5)
	 				ax[ntorIdx].axes.get_yaxis().set_visible(False)
					ImageOptions(fig,ax[ntorIdx],Xlabel,'','$n_{tor}$='+str(ntor),'')
				elif ntorIdx == ntor_tot-1 and ntor_tot > 1:			#If last panel with more panels to left
					cbar = Colourbar(ax[ntorIdx],im,VariableLabel,5)
	 				ax[ntorIdx].axes.get_yaxis().set_visible(False)
					ImageOptions(fig,ax[ntorIdx],Xlabel,'','$n_{tor}$='+str(ntor),'')
				#endif

				#OVERLAY 1D PROFILE OUTLINES ONTO THESE SINGLE KSTEP EQUIL IMAGES
				#MAKES THEM USEFUL FOR QUICK DIAGNOSIS
				#for i in range(0,len(radialprofiles)):
					#Stuff
				#endfor
				#for i in range(0,len(poloidalprofiles)):
					#Stuff
				#endfor

			#endfor

			#Save 2D poloidally resolved figure for current simulation
			SaveString = variable+'_t='+str(round(Time,3))+ext
			plt.savefig(DirEquil2D+SaveString)
#			plt.show()
			plt.close('all')

			#==========#

#			TO BE UPDATED TO INCLUDE OUTPUTS FOR THE NTOR RANGE
#			if write_ASCII == True:
				#Create directory to hold ASCII data
#				DirASCII = CreateNewFolder(DirEquil2D,'2DEquil_Data/')

				#Set ASCII data file name string and header
#				SaveString = variable+'_n'+str(ntor)+'_t='+str(round(Time,3))+'.dat'
#				Header = [VariableLabel,'   ', 'R',lpsi_res, 'theta', ltheta_res,  '\n']

				#Write 1D data header, then 2D PoloidalImage
#				WriteFile_ASCII(Header, DirASCII+SaveString, 'w', 'RSV')
#				WriteFile_ASCII(PoloidalImage, DirASCII+SaveString, 'a', write_ASCIIFormat)
			#endif

		#endfor - Variable loop
	#endfor - dir loop
#endif


#====================================================================#
				 	      #2D POLOIDAL MOVIES#
#====================================================================#

if savefig_2Dequilmovie == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - all need looped over... - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!

		#Create global 2D diagnostics folder and extract current simulation name
		DirEquil2D = CreateNewFolder(Dir[l],'2DEquil_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Kstep [-] & Time [ms] arrays from SEQ.harmonics & toroidal harmonics from energy_n.txt
		SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[l], DataFile='harmonics')
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
		ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
		ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
		ntor0Idx = ntorArray[0]						#ntor = 0 index, contains (var0 + dvar) data

		#Extract toroidal mode number array index (ntorIdx) from requested mode number (ntor)
		ntorIdx = Set_ntorIdx(ntor,ntorArray)

		#Set SEQ and Kstep ranges as requested - else default to max range
		KStepRange,KStepStep = Set_KStepRange(KStepArray,setting_kstep)
		SEQRange = Set_SEQRange(setting_SEQ)

		#Extract relevant spatial normalisation factors
		NormVariables,NormValues,NormUnits = ExtractMEGA_Normalisations(Dir[l])
		ZMin = NormValues[NormVariables.index('bottom_sim')]; ZMax = NormValues[NormVariables.index('top_sim')]
		Zgeo = NormValues[NormVariables.index('zaxis')]; Rgeo = NormValues[NormVariables.index('raxis')]


		for j in range(SEQRange[0],SEQRange[1]):
			#Set SEQIndex for current simulation folder
			SEQ = j

			#Extract and plot data for each timestep
			for i in tqdm( range(KStepRange[0],KStepRange[1],KStepStep) ):

				#Set TimeIndex and employ to extract KStep and Time
				KStepIdx = i; 								#[-]
				IdxOffset = SEQ*KStepMod					#[-]
				KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
				Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

				#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
				#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
				#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
				#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
				HarmonicsData = ExtractMEGA_Harmonics(Dir[l],'All',ntor_tot,KStepIdx,SEQ,'3D')

				#Extract data resolution and poloidal axes from repository .dat files
				#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
				Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)
				Crdz = [Crdz[x]+ZMin for x in range(0,len(Crdz))]		#Offset vertical axis such that Z0 = Zgeo
				DataShape = ExtractMEGA_DataShape(HarmonicsData)
				mpol_res = DataShape[0]; ntor_res = DataShape[1]
				lpsi_res = DataShape[2]; ltheta_res = DataShape[3]

				#For each requested variable at the current Kstep
				for j in range(0,len(variables)):

					#Create global 2D diagnostics folder and extract current simulation name
					DirMovie = CreateNewFolder(DirEquil2D,variables[j]+'_n'+str(ntor))

					#Select variable and Merge 3D Data into 2D poloidal slice
					PoloidalImage = Extract_PoloidalImage(HarmonicsData,variables[j],ntorIdx)

					#==========#

					#Create figure and define Title, Legend, Axis Labels etc...
					fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)

					#Extract Variablelabel and define figure texts
					VariableLabel = VariableLabelMaker(variables[j])
					Title = VariableLabel+', n='+str(ntor)+', t='+str(round(Time,3))+' ms \n Simulation: '+DirString
					Xlabel,Ylabel = 'Radius $R$ [m]', 'Height $Z$ [m]'
					Legend = list()

					#Plot 2D poloidally resolved figure and beautify
					im = ax.contourf(Crdr, Crdz, PoloidalImage, 100); plt.axis('scaled')
					im2 = ax.contour(Crdr, Crdz, PoloidalImage, 20); plt.axis('scaled')
					cbar = Colourbar(ax,im,VariableLabel,5)
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

					#Save 2D poloidally resolved figure for current simulation
					SaveString = variables[j]+'_n'+str(ntor)+'_kstep'+str('%07.f'%KStep)+ext
					plt.savefig(DirMovie+SaveString)
#					plt.show()
					plt.close('all')

					#==========#

					if write_ASCII == True:
						#Create directory to hold ASCII data
						DirASCII = CreateNewFolder(DirEquil2D,'2DEquil_Data/')
						DirASCII_Var = CreateNewFolder(DirASCII,variables[j]+'/')

						#Set ASCII data file name string and header
						SaveString = variables[j]+'_n'+str(ntor)+'_t='+str(round(Time,3))+'.dat'
						Header = [VariableLabel,'   ', 'R',lpsi_res, 'theta', ltheta_res,  '\n']

						#Write 1D data header, then 2D PoloidalImage
						WriteFile_ASCII(Header, DirASCII_Var+SaveString, 'w', 'RSV')
						WriteFile_ASCII(PoloidalImage, DirASCII_Var+SaveString, 'a', write_ASCIIFormat)
					#endif
				#endfor - Variable loop

				#!!! AUTO CREATE MOVIES FOR EACH VARIABLE HERE !!!
				#!!! AUTO CREATE MOVIES FOR EACH VARIABLE HERE !!!

			#endfor - Kstep loop
		#endfor - SEQ loop
	#endfor - dir loop
#endif

#==========##==========##==========#
#==========##==========##==========#

if any([savefig_2Dequilibrium,savefig_2Dequilmovie]) == True:
	print '--------------------------------'
	print '2D Equilibrium Analysis Complete'
	print '--------------------------------'
#endif

#====================================================================#
#====================================================================#
































#====================================================================#
			 #TEMPORALLY/SPECTRALLY RESOLVED DIAGNOSTICS#
#====================================================================#

#====================================================================#
					 #POLOIDAL SPECTRUM ANALYSIS#
#====================================================================#

if savefig_2Dpolspectrum == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - all need looped over... - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!
		variable = SpectralVariable			#requested response variable 		!!! Need to impliment vrad etc...

		#Initiate any required lists
		DataAmpPROES_pol,DataAmpPROES_rad = list(),list()
		XaxisPROES = list()

		#Create global 2D diagnostics folder and extract current simulation name
		DirSpectral = CreateNewFolder(Dir[l],'2DSpectral_Plots/')						#Spatio-Temporal Folder
		DirSpectral_ntor = CreateNewFolder(DirSpectral,variable+'_ntor='+str(ntor))		#Spatio-Temporal Images Folder	
		DirString = Dir[l].split('/')[-2]									#Full Simulation Name
		SubString = DirString.split('_')[-1]								#Simulation Nickname

		#Extract Kstep [-] & Time [ms] arrays from SEQ.harmonics & toroidal harmonics from energy_n.txt
		SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[l], DataFile='harmonics')
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
		ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
		ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
		ntor0 = ntorArray[0]						#ntor = 0, baseline equilibrium data

		#Extract Energy_n outputs and header for plotting
		#energy_n: [ntor][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		Energy_TimeArray = Energy_n[1]				#Extract full time array [ms] for plotting
		Energy_n = Energy_n[2::]					#Remove KStep and Time arrays from array

		#Extract toroidal mode number array index (ntorIdx) from requested mode number (ntor)
		ntorIdx = Set_ntorIdx(ntor,ntorArray)

		#Set Kstep ranges as requested - else default to max range
		KStepRange,KStepStep = Set_KStepRange(KStepArray,setting_kstep)
		SEQRange = Set_SEQRange(setting_SEQ)


		#Extract Variablelabel for chosen variable
		VariableLabel = VariableLabelMaker(variable)			#Units='Perturbation [-]'

		for j in tqdm( range(SEQRange[0],SEQRange[1])  ):
			#Set SEQIndex for current simulation folder
			SEQ = j

			#Extract and plot data for each timestep
			for i in range(KStepRange[0],KStepRange[1],KStepStep):

				#Set TimeIndex and employ to extract KStep and Time
				KStepIdx = i								#[-]
				IdxOffset = SEQ*KStepMod					#[-]
				KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
				Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

				#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
				#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
				#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
				#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
				HarmonicsData = ExtractMEGA_Harmonics(Dir[l],'All',ntor_tot,KStepIdx,SEQ,'3D')
				rho_pol = HarmonicsData.rho_pol; q_psi = HarmonicsData.q_psi

				#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
				DataShape = ExtractMEGA_DataShape(HarmonicsData)#; print DataShape
				mpol_res = DataShape[0]; ntor_res = DataShape[1]
				lpsi_res = DataShape[2]; ltheta_res = DataShape[3]

				#Extract radial magnetic field (brad) from SEQ.harmonic object
				#Data is of shape: Data[mpol,ntor,lpsi,A/B]
				Data = getattr(HarmonicsData, variable)

				#Combine spectral components A and B in quadrature to obtain variable amplitude
				#Pos corresponds to resonant poloidal modes 	i.e. +m on RHS of image
				#Neg corresponds to non-resonant poloidal modes i.e. -m on LHS of image
				#One of ntor_pos-ntor or ntor_pos+ntor will equal 0, representing the equilibrium values.
				DataAmpPos = np.sqrt( (Data[:, ntor_pos-ntor,:,0]**2) + (Data[:, ntor_pos-ntor,:,1]**2) )
				DataAmpNeg = np.sqrt( (Data[:, ntor_pos+ntor,:,0]**2) + (Data[:, ntor_pos+ntor,:,1]**2) )
				DataAmpNeg = np.flip( DataAmpNeg,axis=0)		#Flip LHS of image for plotting

				#Concat positive and negative ntor to obtain full poloidal harmonic spectrum
				#DataAmp is of shape: [2*mpol-1][lpsi]
				DataAmp = np.concatenate((DataAmpNeg,DataAmpPos[1:,:]),axis=0)

				#Create Image array and Axes, rotate such that mpol spectrum is on X-axis.
				#Image is of shape: [lpsi][2*mpol+1]	(i.e. Image is orientated [Y,X])
				Image = DataAmp.transpose()
				Xaxis =	[x-int(mpol_res-1) for x in range(0,2*mpol_res-1,1)]	#Poloidal Mode Numbers	[mpolAxis] 
				Yaxis = rho_pol													#Radial Location		[lpsiAxis]

				#If QuickPROES not used, plot a poloidal spectrum for each kstep value
				if QuickPROES == False:
					#Create figure and define Title, Legend, Axis Labels etc...
					AspectRatio = [image_aspectratio[0],image_aspectratio[1]*1.25]
					fig,ax = figure(subplots=[2,1], gridspec=[2,1], aspectratio=AspectRatio)

					#Plot poloidal spectrum figure (R,mpol)
					extent = [Xaxis[0],Xaxis[-1], Yaxis[0],Yaxis[-1]]						#[mpolAxis, lpsiAxis]
					im = ax[0].imshow(Image, extent=extent, aspect='auto', origin='bottom')	#Image orientated [Y,X]
					co = ax[0].contour(Image, extent=extent, origin='lower', levels=10)		#Image orientated [Y,X]
#					im = ax.contourf(Xaxis, Yaxis, Image, 50)
					res = ax[0].plot(-ntor*q_psi, rho_pol, 'w--', lw=2)
					cbar = Colourbar(ax[0],im,VariableLabel,5)
					#####
					Title = 'Poloidal Spectrum: n='+str(ntor)+', m='+str(-mpol_res+1)+','+str(mpol_res-1)+', t='+str(round(Time,3))+' [ms] \n Simulation: '+DirString
					Xlabel,Ylabel = 'Poloidal Harmonic $m_{pol}$ [-]', 'Radial Magnetic Coordinate $\\rho_{pol}$ [-]'
					ImageOptions(fig,ax[0],Xlabel,Ylabel,Title,'')
					ax[0].set_xlim(image_mpolcrop[0],image_mpolcrop[1])

					#Plot total energy for each harmonic component (where "except:" accounts for energy_n -n values)
					try: ax[1].plot(Energy_TimeArray,np.log10(Energy_n[ntorIdx]), lw=2)
					except: ax[1].plot(Energy_TimeArray,np.log10(Energy_n[ntorIdx-ntor_pos]), lw=2)
					ax[1].axvline(TimeArray[KStepIdx+IdxOffset],0,1)
					cbar = InvisibleColourbar(ax[1])
					###
					Xlabel,Ylabel = 'Time [ms]', 'Mode Energy $n_{tor}$ [-]'
					Legend = ['n$_{tor}$ = '+str(ntor)]
					ImageOptions(fig,ax[1],Xlabel,Ylabel,'',Legend)

					#Save poloidal spectrum figure for current SEQ and Kstep
					SaveString = 'PoloidalSpectrum_'+variable+'_n'+str(ntor)+'_kstep'+str('%07.f'%KStep)+ext
					plt.savefig(DirSpectral_ntor+SaveString)
#					plt.show()
					plt.close('all')

					if write_ASCII == True:
						DirASCII = CreateNewFolder(DirSpectral_ntor,"ASCII_Data")			#Spatio-Temporal Data Folder
						
						#Save Yaxis (rho_pol) and safety factor for future plotting
						WriteFile_ASCII(rho_pol, DirASCII+'rho_pol', 'w', write_ASCIIFormat)
						WriteFile_ASCII(q_psi, DirASCII+'q_psi', 'w', write_ASCIIFormat)

						#Write 1D data header, then 2D Radially resolved Spatio-Temporal Image
						SaveString = 'PolSpectrum_'+variable+'_n'+str(ntor)+'_t='+str(round(Time,3))+'.dat'
						Header = [VariableLabel,'   ', 'mpol',extent[0],extent[1], 'rho_pol',extent[2],extent[3], '\n']
						WriteFile_ASCII(Header, DirASCII+SaveString, 'w', 'RSV')
						WriteFile_ASCII(Image, DirASCII+SaveString, 'a', write_ASCIIFormat)
					#endif
				#endif

				#==========##==========#
				#==========##==========#

				#Collapse figure radially to create 'PROES'-like temporal evolution figure
				#Integrate through all radii, maintaining poloidal spectrum (mpol) resolution
				DataAmp1D_pol = list()
				for k in range(0,len(DataAmp)):
					#DataAmp is of shape: [2*mpol-1][lpsi]
					#DataAmp1D is of shape: [2*mpol-1]
					DataAmp1D_pol.append(sum(DataAmp[k][:]))
				#endfor

				#Collapse figure poloidally to create 'PROES'-like temporal evolution figure
				#Integrate through all poloidal modes, maintaining radial (rho_pol) resolution
				DataAmp1D_rad = list()
				DataAmp = DataAmp.transpose()
				for k in range(0,len(DataAmp)):
					#Transposed DataAmp is of shape: [lpsi][2*mpol-1]
					#DataAmp1D is of shape: [lpsi]
					DataAmp1D_rad.append(sum(DataAmp[k][:]))
				#endfor

				#Append 1D spatial arrays into 2D spati-temporal 'PROES-like' image arrays
				#DataAmpPROES_pol: 2D array of shape [kstep][2*mpol+1]
				#DataAmpPROES_rad: 2D array of shape [kstep][lpsi]
				DataAmpPROES_pol.append(DataAmp1D_pol)
				DataAmpPROES_rad.append(DataAmp1D_rad)
				XaxisPROES.append(Time)

			#endfor - Kstep loop
		#endfor - SEQ loop

		#================##=================#
		#	TEMPORALLY RESOLVED PROFILES	#
		#================##=================#

		#Plot spatio-temporally resolved poloidal spectrum figures
		if len(XaxisPROES) > 1:

			#Create 'PROES-like' Image array, rotated such that time is on X-axis.
			#PROESImage is of shape: [lpsi OR 2*mpol+1][kstep]	(i.e. Image is orientated [Y,X])
			PROESImage_pol = np.asarray(DataAmpPROES_pol).transpose()
			PROESImage_rad = np.asarray(DataAmpPROES_rad).transpose()

			#Set image extents for each orientation
			extent_pol = [XaxisPROES[0],XaxisPROES[-1], Xaxis[0],Xaxis[-1]]		#[2*mpol+1][kstep]
			extent_rad = [XaxisPROES[0],XaxisPROES[-1], Yaxis[0],Yaxis[-1]]		#[lpsi][kstep]

			#==========##==========#

			#Radially resolved 'PROES-like' figure: Title, Legend, Axis Labels etc...
			fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)
			Title = 'Poloidally Collapsed: n='+str(ntor)+', m='+str(-mpol_res+1)+','+str(mpol_res-1)+', t='+str(round(Time,3))+' [ms] \n Simulation: '+DirString
			Xlabel,Ylabel = 'Time $t$ [ms]', 'Radial Magnetic Coordinate $\\rho_{pol}$ [-]'
			Legend = list()

			#Plot temporally resolved, poloidally collapsed, figure (R,time)
			im = ax.imshow(PROESImage_rad, extent=extent_rad, aspect='auto', origin='bottom')	#Image orientated [Y,X]
			co = ax.contour(PROESImage_rad, extent=extent_rad, origin='lower', levels=20)		#Image orientated [Y,X]
#			im = plt.contourf(XaxisPROES, Yaxis, PROESImage_rad, 50)
			cbar = Colourbar(ax,im,VariableLabel,5)
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
			ax.set_ylim(0,1)									#ax.set_ylim(image_rhocrop[0],image_rhocrop[1])

			#Save temporal spectrum figure for current simulation directory
			SaveString = 'RadialSpectrum_'+variable+'_n'+str(ntor)+'_t='+str(round(Time,3))+ext
			plt.savefig(DirSpectral+SaveString)
#			plt.show()
			plt.close('all')

			#==========##==========#

			#Poloidially resolved 'PROES-like' figure: Title, Legend, Axis Labels etc...
			fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)
			Title = 'Radially Collapsed: n='+str(ntor)+', m='+str(-mpol_res+1)+','+str(mpol_res-1)+', t='+str(round(Time,3))+' [ms] \n Simulation: '+DirString
			Xlabel,Ylabel = 'Time $t$ [ms]', 'Poloidal Harmonic $m_{pol}$ [-]'
			Legend = list()

			#Plot temporally resolved, radially collapsed, figure (R,time)
			im = ax.imshow(PROESImage_pol, extent=extent_pol, aspect='auto', origin='bottom')	#Image orientated [Y,X]
			co = ax.contour(PROESImage_pol, extent=extent_pol, origin='lower', levels=20)		#Image orientated [Y,X]
#			im = plt.contourf(XaxisPROES, Xaxis, PROESImage_pol, 50)
			cbar = Colourbar(ax,im,VariableLabel,5)
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
			ax.set_ylim(image_mpolcrop[0],image_mpolcrop[1])

			#Save temporal spectrum figure for current simulation directory
			SaveString = 'PoloidalSpectrum_'+variable+'_n'+str(ntor)+'_t='+str(round(Time,3))+ext
			plt.savefig(DirSpectral+SaveString)
#			plt.show()
			plt.close('all')

			#==========##==========#
			#==========##==========#

			if write_ASCII == True:
				DirASCII = CreateNewFolder(DirSpectral,"ASCII_Data")		#Spatio-Temporal Data Folder

				#Write 1D data header, then 2D Radially resolved Spatio-Temporal Image
				SaveString = 'RadialSpectrum_'+variable+'_n'+str(ntor)+'_t='+str(round(Time,3))+'.dat'
				Header = [VariableLabel,'   ', 'time',extent_rad[0],extent_rad[1], 'rho_pol',extent_rad[2],extent_rad[3], '\n']
				WriteFile_ASCII(Header, DirASCII+SaveString, 'w', 'RSV')
				WriteFile_ASCII(PROESImage_rad, DirASCII+SaveString, 'a', write_ASCIIFormat)

				#Write 1D data header, then 2D Poloidally resolved Spatio-Temporal Image
				SaveString = 'PoloidalSpectrum_'+variable+'_n'+str(ntor)+'_t='+str(round(Time,3))+'.dat'
				Header = [VariableLabel,'   ', 'time',extent_pol[0],extent_pol[1], 'mpol',extent_pol[2],extent_pol[3], '\n']
				WriteFile_ASCII(Header, DirASCII+SaveString, 'w', 'RSV')
				WriteFile_ASCII(PROESImage_pol, DirASCII+SaveString, 'a', write_ASCIIFormat)
			#endif
		#endif - PROES plotting branch
	#endfor - Dir loop
#endif - Diag loop


#==========##==========##==========#
#==========##==========##==========#

if any([savefig_2Dpolspectrum]) == True:
	print '--------------------------------------'
	print '2D Poloidal Spectrum Analysis Complete'
	print '--------------------------------------'
#endif

#====================================================================#
#====================================================================#































#====================================================================#
				   #SPECTRAL & HARMONIC DIAGNOSTICS#
#====================================================================#

#====================================================================#
				   #2D CONTINUUM & FOURIER ANALYSIS#
#====================================================================#

if savefig_2Dcontinuum == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - all need looped over... - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		SEQ = setting_SEQ[1]			#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
#		ntor = setting_ntor[1]			#requested ntor mode number			!!! NEEDS A FUNCTION !!!
		variable = ContinuumVariable	#requested continuum variable 		!!! Need to impliment btheta, bphi etc...

		#Create global 2D diagnostics folder and extract current simulation name
		DirContinuum = CreateNewFolder(Dir[l],'2DContinuum_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Kstep [-] & Time [ms] arrays from SEQ.harmonics & toroidal harmonics from energy_n.txt
		SEQArray, KStepArray, TimeArray, ntorArray = ExtractMEGA_DataRanges(Dir[l], DataFile='harmonics')
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/len(SEQArray)	#KStep indices per SEQ 	[-]
		ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
		ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
		ntor0 = ntorArray[0]						#ntor = 0, baseline equilibrium data

		#Extract toroidal mode number array index (ntorIdx) from requested mode number (ntor)
#		ntorIdx = Set_ntorIdx(ntor,ntorArray)

		#Extract Variablelabel for chosen variable
		VariableLabel = VariableLabelMaker(variable,Units=' \n Perturbation [-]')


		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		AlfvenVelocity = Values[Variables.index('Alfven velocity')] #B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = Values[Variables.index('D gyro frequency')]
		IonDensity = Values[Variables.index('Bulk density')]
		B0 = Values[Variables.index('Mag.fld. at axis')]
		R0 = Values[Variables.index('raxis')]
		m_D = 3.34e-27
		eps = 0.5/R0

		#Extract ONE VARIABLES FOR ALL KSTEP from Harmonics, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.data: [4D Array] of shape [kstep][mpol][ntor][lpsi][A/B] for [variable]
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l],variable,ntor_tot,Dimension='4D')
		rho_pol = HarmonicsData.rho_pol; q_psi = HarmonicsData.q_psi;
		Data = HarmonicsData.data

		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		DataShape = ExtractMEGA_DataShape(HarmonicsData)
		mpol, ntor = DataShape[0], DataShape[1]
		lpsi, ltheta = DataShape[2], DataShape[3]
		kmax, dt = DataShape[4], (TimeArray[1]-TimeArray[0])


		#TO DO 	::
			#	:: CHECK THE RHO POL X-AXIS HAS BEEN APPLIED CORRECTLY, MAY BE rho_pol = sqrt(MEGA(Rho_pol))		???
			#	:: ENABLE SELECTION OF VARIABLES TO BE PLOTTED - UPDATE TITLE AND SAVESTRING ACCORDINGLY
			#	:: ENABLE SELECTION OF TOROIDAL MODE NUMBERS TO BE PLOTTED - UPDATE TITLE ACCORDINGLY
			#	:: COMPUTE A 'GROWTH' FUNCTION AND HAVE A TOGGLEABLE SCALING SUCH THAT THE PERTURBATIONS ARE 'FLAT' (PABLO)
			#	:: ENABLE FOURIER TRANSFORM TO BE PERFORMED OVER A USER-DEFINED TIMESCALE
			#	:: UPDATE THE ComputeTAEThresholds() FUNCTION WITH COMMENTS, CITATIONS, AND CHECK MATHS				(PABLO)
			#	:: TRANSLATE AND COMMENT ANY REMAINING ORIGINAL JAVI CODE BELOW
			#	::

		print kmax, mpol, ntor, lpsi

		#Sum Re component of toroidal (n) and poloidal (m) modes for all ksteps
		vcos = list()
		for n in range(0,ntor_tot):
			vcos.append( np.zeros([]) )				
			for m in range(0,mpol):						#Data structure: [kstep][mpol][ntor][lpsi][Re/Im] 
				vcos[n] = vcos[n] + Data[:,m,n,:,0]		#vcos structure: [ntor][kstep][lpsi]
			#endfor
			vcos[n][np.isnan(vcos[n])] = 0				#Remove any NaNs
		#endfor

#		print len(vcos), len(vcos[0]), len(vcos[0][0])
#		print vcos[0]
#		plt.plot(vcos[0])
#		plt.show()
#		exit()

		vcos_fft,vcos_len = list(),list()
		#Extract fourier components from vcos
		for n in range(0,ntor_tot):
			vcos_fft.append( np.fft.fft(vcos[n],axis=0) )	#Take fourier components of variable
		  	vcos_fft[n][0,:] = vcos_fft[n][0,:]*0.0			#Discard imaginary components 					???
			vcos_len.append( int(len(vcos_fft[n])/2)-1 )	#Determine lowpass filter frequency threshold 	???
			vcos_fft[n] = vcos_fft[n][0:vcos_len[n],:]		#Discard upper frequencies (lowpass filter)		???
		#endfor

		#==========##==========#

		#Create fig of desired size - increasing Xlim with the number of harmonics
		Xaspect, Yaspect = int(10*(float(ntor)/1.75)), 12
		fig,ax = figure(subplots=[2,ntor_tot], aspectratio=[Xaspect,Yaspect])

		#For each toroidal harmonic:
		for i in range(0,ntor_tot):

			#Temporal evolution plotted on the top row (row 0)
			if ntor_tot == 1: subfig = ax[0]
			elif ntor_tot > 1: subfig = ax[0,i]
			#endif

			Harmonic = -ntor_pos+i
			#Construct temporal figure axes and meshgrid (not used)
			Xaxis = rho_pol										#[-]
			Yaxis = TimeArray									#[ms]
			extent = [Xaxis[0],Xaxis[-1], Yaxis[0],Yaxis[-1]]
			X,Y = np.meshgrid(Xaxis,Yaxis)						#im = subfig.contourf(X,Y,vcos[i])

			#Plot harmonic temporal evolution
			im = subfig.imshow(vcos[i][::-1], extent=extent, aspect='auto')
			co = subfig.contour(vcos[i], extent=extent, levels=10)
			
			#Add colourbar and beautify plot - taking account of panel location
			if i == 0 and ntor_tot > 1: 					#If first panel with more panels to right
				cbar = Colourbar(subfig,im,'',5)
				ImageOptions(fig,subfig,'','Time [ms]','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
			elif i == 0 and ntor_tot == 1:					#If first panel with no panels to right
				cbar = Colourbar(subfig,im,VariableLabel,5)	#Single Panel colourbar (for reference)
				ImageOptions(fig,subfig,'','Time [ms]','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
			elif i > 0 and i < ntor_tot-1: 					#If middle panel with more panels to right
				cbar = Colourbar(subfig,im,'',5)
				ImageOptions(fig,subfig,'','','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
 				im.axes.get_yaxis().set_visible(False)
			elif i == ntor_tot-1 and ntor_tot > 1:			#If last panel with more panels to left
				cbar = Colourbar(subfig,im,VariableLabel,5)	#Right-most colourbar (for reference)
				ImageOptions(fig,subfig,'','','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
 				im.axes.get_yaxis().set_visible(False)
			#endif

			#==========#

			#Alfven continuum (Fourier) analysis plotted on the bottom row (row 1)
			if ntor_tot == 1: subfig = ax[1]
			elif ntor_tot > 1: subfig = ax[1,i]

			#Construct frequency figure axes and meshgrid (not used)
			Xaxis = rho_pol										#[-]
			Yaxis = np.linspace(0,0.5/dt,vcos_len[i])			#[kHz]
			extent = [Xaxis[0],Xaxis[-1], Yaxis[0],Yaxis[-1]]
			X,Y = np.meshgrid(Xaxis,Yaxis)						#im = subfig.contourf(X,Y,vcos_fft[i])

			#Plot Fourier amplitude spectrum
			im = subfig.imshow(real(vcos_fft[i])[::-1], extent=extent, aspect='auto')
			co = subfig.contour(real(vcos_fft[i]), extent=extent, levels=10)

			#Add colourbar and beautify plot - taking account of panel location
			if i == 0 and ntor_tot > 1: 					#If first panel with more panels to right
				cbar = Colourbar(subfig,im,'',5)
				ImageOptions(fig,subfig,'Normalised Minor Radius $\\rho_{pol}$','Frequency [kHz]','','')
			elif i == 0 and ntor_tot == 1:					#If first panel with no panels to right
				cbar = Colourbar(subfig,im,VariableLabel,5)	#Single Panel colourbar (for reference)
				ImageOptions(fig,subfig,'Normalised Minor Radius $\\rho_{pol}$','Frequency [kHz]','','')
			elif i > 0 and i < ntor_tot-1:   				#If middle panel with more panels to right
				cbar = Colourbar(subfig,im,'',5)
				ImageOptions(fig,subfig,'Normalised Minor Radius $\\rho_{pol}$','','','')
 				im.axes.get_yaxis().set_visible(False)
			elif i == ntor_tot-1 and ntor_tot > 1:			#If last panel with more panels to left
				cbar = Colourbar(subfig,im,VariableLabel,5) #Right-most colourbar (for reference)
				ImageOptions(fig,subfig,'Normalised Minor Radius $\\rho_{pol}$','','','')
 				im.axes.get_yaxis().set_visible(False)
			#endif
			subfig.set_ylim([0,200])

			#Compute and plot TAE thresholds
			UpperThresholds,LowerThresholds = ComputeTAEThresholds(HarmonicsData,Harmonic,eps,AlfvenVelocity,subfig)
		#endfor

		#Minimise spacing between subplots (Spacing may need to be a function of ntor_tot)
		plt.subplots_adjust(wspace=0.2, hspace=0.1)

		#Save continuum figure for variable[j] and simulation folder [l]
		SaveString = variable+'_Continuum_'+SubString+ext
		plt.savefig(DirContinuum+SaveString)
#		plt.show()
		plt.close('all')
	#endfor
#endif

#==========##==========##==========#
#==========##==========##==========#

if any([savefig_2Dcontinuum]) == True:
	print '-----------------------------'
	print '2D Spectral Analysis Complete'
	print '-----------------------------'
#endif

#====================================================================#
#====================================================================#































#====================================================================#
				  #KINETIC & PHASESPACE DIAGNOSTICS#
#====================================================================#

#====================================================================#
			  			 #1D IEDF ANALYSIS#
#====================================================================#

if savefig_1Dkinetics == True:

	#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
	print Dir[l].split('/')[-2]
	nBins = 100					#Kinetics Histogram Bins		- Move to Low-Level Inputs
	KStepMin = 0				#KStepMin						- Automate readin - Use Switchboard?
	KStepMax = 200000			#KStepMax						- Automate readin - Use Switchboard?
	KWep = 10000				#Write_ep save interval (kwep)	- Automate readin - Use Switchboard?
	KMarker = 8					#Marker file readin interval	- Move to Low-Level Inputs

	#Cycle through all simulation folders
	for l in range(0,len(Dir)):

		#Create global kinetics folder and extract current simulation name
		DirKinetics = CreateNewFolder(Dir[l],'2DKinetic_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]


		#KINETICS VARIABLE LOOP GOES HERE

		#Initiate KineticPROES
		KineticPROES = list()

		#Initiate figure and set axes
		fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)
		Xlabel,Ylabel = 'Energy $\epsilon_{i}$ [keV]', 'Ion Energy Distribution Function $f(\epsilon_{i})$ [-]'
		Legend = list()

		#Cycle through all Kstep for given kwep.
		for i in range(KStepMin,KStepMax+1,KWep):

			#Set current KStep
			KStep = i
			#Concatenate variables into KineticsData - Override KineticsData on first iteration
			#KineticsData :: 2D Array of shape [variable,marker(n)]
			#Variables :: R, Z, Lambda, E, p, Mu, pphi, fff*fnrml, psip, phi
			KineticsData,Header_Kin = ExtractMEGA_Markers(Dir[l],KStep,KMarker)


			#Select variable to be plotted (X axis) and histogram into nBins
			XData = KineticsData[3]									#'E_gc'
			HistData1D,XAxis = np.histogram(XData, bins=nBins)

			#Normalise distribution function
			HistSum1D = sum(HistData1D); NormFactor1D = HistSum1D
			HistData1D = [float(HistData1D[x])/float(NormFactor1D) for x in range(0,len(HistData1D))]
			if DebugMode == True: print( "IEDF Integral: ",str(sum(HistData1D)) )

			#Plot figure for current KStep
			ax.plot(XAxis[0:-1], HistData1D, lw=2)

			#Append 1D data to KineticPROES
			KineticPROES.append(HistData1D)
		#endfor

		#Apply image options and save figure - One variable, all KSteps
		ImageOptions(fig,ax,Xlabel,Ylabel,'',Legend)
		plt.show()
		plt.close('all')

		#If more than one KStep was processed, create a temporal IEDF image
		if len(KineticPROES) > 1:
			#Compute mean and modal values
			MeanArray,ModeArray = list(),list()

			#Initiate figure and set axes
			fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)
			Xlabel,Ylabel = 'Energy $\epsilon_{i}$ [keV]','Time $t$ [ms]'
			Legend = list()

			#Plot figure for current KStep
			im = ax.imshow(KineticPROES, aspect='auto', origin='bottom')
			cbar = Colourbar(ax,im,'Ion Energy Distribution Function $f(\epsilon_{i})$ [-]',5)
			ImageOptions(fig,ax,Xlabel,Ylabel,'',Legend)
			plt.show()
		#endif
	#endfor
#endif

#====================================================================#
			  			 #2D IEDF ANALYSIS#
#====================================================================#

if savefig_2Dkinetics == True:

	#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
	print Dir[l].split('/')[-2]
	KMarker = 1					#Marker file readin interval	- Move to Low-Level Inputs
	nBins = 100					#Kinetics Histogram Bins		- Move to Low-Level Inputs
	KStepMin = 000000			#KStepMin						- Automate readin - Use Switchboard?
	KStepMax = 3000000			#KStepMax						- Automate readin - Use Switchboard?
	KWep = 100000				#Write_ep save interval (kwep)	- Automate readin - Use icp.nam readin function?

	Labels = ['Radius $R$ [m]','Height $Z$ [m]','pitch angle $\lambda = \\frac{v_{para}}{v}$ [-]','Energy $\epsilon_{i}$ [keV]','Momentum $p$ \n [kg m${^2}$ s$^{-1}$]','Magnetic Moment $\mu$ [N m T$^{-1}$]','Canonical Momentum $p_{\phi}$ \n [kg m${^2}$ s$^{-1}$]','fff*fnrml [-]','psip [-]','Toroidal Angle $\phi$ [Rads]']
	#N.B. Magnetic moment units: '[A m$^{2}$]', '[J T$^{-1}$]', '[N m T$^{-1}$]'
	#Uppercase Lambda: 'Lambda $\Lambda = \\frac{\mu B}{E}$ [-]'
	#lower lambda = v_para / velocity    ==  	pitch angle

	#Cycle through all simulation folders
	for l in range(0,len(Dir)):

		#Create global kinetics folder and extract current simulation name
		DirKinetics = CreateNewFolder(Dir[l],'2DKinetic_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]


		#KINETICS VARIABLE LOOP GOES HERE

		#Cycle through all Kstep for given kwep.
		for i in range(KStepMin,KStepMax+1,KWep):

			#Set current KStep
			KStep = i

			#Concatenate variables into KineticsData - Override KineticsData on first iteration
			#KineticsData :: 2D Array of shape [variable,marker(n)]
			#Variables :: 0:R, 1:Z, 2:lambda, 3:E, 4:p, 5:Mu, 6:pphi, 7:fff*fnrml, 8:psip, 9:phi
			KineticsData,Header_Kin = ExtractMEGA_Markers(Dir[l],KStep,KMarker)

			#Select variables to be plotted (X,Y axis)
			#Note: "Fast-Ion phase space" refers to plotting Energy vs Canonical Momentum (Y, X)
			XData = KineticsData[0]				# [6] 'pphi_gc'	[Typically Physical Variable]		[3],[6]
			YData = KineticsData[1]				# [3] 'E_gc'	[Typically Spatial Variable]		[0],[1],[9]
			MarkerWeight = KineticsData[7]

			Xlabel = Labels[0]	# Temporary fudge
			Ylabel = Labels[1]	# Temporary fudge

			#	Not sure what this is, but Javi plotted them...
			#	Lambda1 = 2.1*KineticsData[5,:]/(KineticsData[3,:]*1.6e-16)
			# 	HistData2D,XAxis2D,YAxis2D = np.histogram2d(KineticsData[6,:], Lambda1, bins=(XRange,YRange)) 
			# 	HistData2D,XAxis2D,YAxis2D = np.histogram2d(KineticsData[3,:], Lambda1, bins=(XRange,YRange))

			#Extract min/max values and create histogram ranges
			Xmin,Xmax = min(XData), max(XData)
			Ymin,Ymax = min(YData), max(YData)
			XRange = np.linspace(Xmin,Xmax,nBins)
			YRange = np.linspace(Ymin,Ymax,nBins)

			#Select 2D variables to be plotted (X,Y) and histogram over supplied ranges
			HistData2D,XAxis2D,YAxis2D = np.histogram2d(XData,YData, bins=(XRange,YRange))
			HistData2D = np.rot90(HistData2D,1)									#1 is pointy side down, 3 is pointy side up
#			extent = [XAxis2D[0],XAxis2D[-1], YAxis2D[0],YAxis2D[-1]]			#Axes are orientation dependant
			extent = [min(XAxis2D),max(XAxis2D), min(YAxis2D),max(YAxis2D)]		#Axes are orientation indepdenent (ish)

			#Select 1D variable to be plotted (X axis) and histogram into nBins
			XHistData1D,XAxis1D = np.histogram(XData, bins=nBins)

			#Select 1D variable to be plotted (Y axis) and histogram into nBins
			YHistData1D,YAxis1D = np.histogram(YData, bins=nBins)


			#Normalise 2D distribution function
			HistSum2D = sum(sum(HistData2D)); NormFactor2D = HistSum2D
			for x in range(0,len(HistData2D)):
				for y in range(0,len(HistData2D[x])):
					HistData2D[x,y] = float(HistData2D[x,y])/float(NormFactor2D)
				#endfor
			#endfor
			if DebugMode == True: print( "2D IEDF Integral: ",str(sum(HistData2D)) )

			#Normalise 1D (X axis) distribution function
			XHistSum1D = sum(XHistData1D); XNormFactor1D = XHistSum1D
			XHistData1D = [float(XHistData1D[x])/float(XNormFactor1D) for x in range(0,len(XHistData1D))]
			if DebugMode == True: print( "X 1D IEDF Integral: ",str(sum(XHistData1D)) )

			#Normalise 1D (Y axis) distribution function
			YHistSum1D = sum(YHistData1D); YNormFactor1D = YHistSum1D
			YHistData1D = [float(YHistData1D[x])/float(YNormFactor1D) for y in range(0,len(YHistData1D))]
			if DebugMode == True: print( "Y 1D IEDF Integral: ",str(sum(YHistData1D)) )

			#==========#

			fig,ax = figure(subplots=[1,1], aspectratio=image_aspectratio)	

#			Title = VariableLabels[j]+', t='+str(Time)+' \n Simulation: '+DirString
			Title = 'Kinetic Markers, Kstep='+str(KStep).zfill(7)+' \n Simulation: '+DirString
#			Xlabel,Ylabel = 'Energy $\epsilon_{i}$ [keV]','Canonical Momentum $p_{\phi}$ \n [kg m${^2}$ s$^{-1}$]'
			Legend = list()

			#Plot 2D toroidally resolved IEDF
			im1 = ax.imshow(HistData2D, extent=extent, aspect='auto')
			cbar1 = Colourbar(ax,im1,'Marker Count [normalised]',5)
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
			ax.set_xlim(1.2,max(XAxis1D))
			ax.set_ylim(min(YAxis1D),max(YAxis1D))


			if True == False:
				#Initiate figure and set axes
				fig,ax = figure(subplots=[2,1], aspectratio=image_aspectratio, shareX=True)	

	#			Title = VariableLabels[j]+', t='+str(Time)+' \n Simulation: '+DirString
				Title = 'Kinetic Markers, Kstep='+str(KStep).zfill(7)+' \n Simulation: '+DirString
	#			Xlabel,Ylabel = 'Energy $\epsilon_{i}$ [keV]','Canonical Momentum $p_{\phi}$ \n [kg m${^2}$ s$^{-1}$]'
				Legend = list()

				#Plot 2D toroidally resolved IEDF
				im1 = ax[0].imshow(HistData2D, extent=extent, aspect='auto')
				cbar1 = Colourbar(ax[0],im1,'IEDF $f(\epsilon_{i})$ [-]',5)
				ImageOptions(fig,ax[0],'',Ylabel,Title,Legend)
				ax[0].set_xlim(min(XAxis1D),max(XAxis1D))								#USE FIXED FULL RANGE
				ax[0].set_ylim(min(YAxis1D),max(YAxis1D))								#USE FIXED FULL RANGE
				
				#Plot 1D Y-axis integrated IEDF		(i.e. Y-axis collapsed into 1D)						
				im2 = ax[1].plot(XAxis1D[0:-1], XHistData1D, lw=2)
				cbar2 = InvisibleColourbar(ax[1])
				ImageOptions(fig,ax[1],Xlabel,'IEDF $f(\epsilon_{i})$ [-]','',Legend)
				ax[1].set_xlim(min(XAxis1D),max(XAxis1D))								#USE FIXED FULL RANGE
			#endif

			#Save temporal response figure for current simulation directory
			SaveString = 'Kinetics_kstep'+str(KStep).zfill(7)+ext
			plt.savefig(DirKinetics+SaveString)
#			plt.show()
			plt.close('all')

			#==========#

			if write_ASCII == True:
				#Create directory to hold ASCII data
				DirASCII = CreateNewFolder(DirKinetics,'Kinetics_Data/')
#				DirASCII_Var = CreateNewFolder(DirASCII,variables[j]+'/')

				#Set ASCII data file name string and header
#				SaveString = variables[j]+'_n'+str(ntor)+'_t='+str(round(Time,3))+'.dat'
				SaveString = 'kstep='+str(KStep).zfill(7)+'.dat'
				Header = ['VariableLabel','   ', 'X',[Xmin,Xmax], 'Y', [Ymin,Ymax],  '\n']

				#Write 1D data header, then 2D PoloidalImage
				WriteFile_ASCII(Header, DirASCII+SaveString, 'w', 'RSV')
				WriteFile_ASCII(HistData2D.T, DirASCII+SaveString, 'a', write_ASCIIFormat)
			#endif
		#endfor - kstep loop
	#endfor - simulation folder loop
#endif - Diagnostic loop

#==========##==========##==========#
#==========##==========##==========#

if any([savefig_1Dkinetics,savefig_2Dkinetics]) == True:
	print '-------------------------'
	print 'Kinetic Analysis Complete'
	print '-------------------------'
#endif

#====================================================================#
#====================================================================#

#===================================================================#
#===================================================================#
#																	#
#							END OF SCRIPT							#		
#																	#
#===================================================================#
#===================================================================#

exit()


















































#====================================================================#
							# CODE DUMP #
#====================================================================#

# UNUSED OR OUTDATED SNIPPITS OF CODE ARE STORED HERE.




#=========================#
#=========================#

			#DOUBLE 1D HISTOGRAM WITH 2x2 FIGURE PLOT FOR savefig_2Dkinetics DIAGNOSTIC

KineticSubplots = False
if KineticSubplots == True:
	#Initiate figure and set axes
	fig,ax = figure(subplots=[2,2], aspectratio=image_aspectratio)		#shareX=True
	ax[1,1].axis('off')		

#	Title = VariableLabels[j]+', t='+str(Time)+' \n Simulation: '+DirString
	Title = 'Kinetic Markers, Kstep='+str(KStep).zfill(7)+' \n Simulation: '+DirString
#	Xlabel,Ylabel = 'Energy $\epsilon_{i}$ [keV]','Canonical Momentum $p_{\phi}$ \n [kg m${^2}$ s$^{-1}$]'
	Legend = list()

	#Set global figure options
	fig.suptitle(Title, y=1.01)

	#Plot 2D toroidally resolved IEDF
	im1 = ax[0,0].imshow(HistData2D, extent=extent, aspect='auto')
#	ln1 = ax[0,0].plot(XAxis,np.zeros(len(XAxis)), 'r--', lw=2)
	cbar1 = Colourbar(ax[0,0],im1,'IEDF $f(\epsilon_{i})$ [-]',5)
	ImageOptions(fig,ax[0,0],'',Xlabel,'',Legend)
	
	#Plot 1D Y-axis integrated IEDF		(i.e. Y-axis collapsed into 1D)						
	im2 = ax[1,0].plot(XAxis1D[0:-1], XHistData1D, lw=2)
	cbar2 = InvisibleColourbar(ax[1,0])
	ImageOptions(fig,ax[1,0],Ylabel,'IEDF $f(\epsilon_{i})$ [-]','',Legend)
	ax[1,0].set_xlim(min(XAxis1D),max(XAxis1D))									#USE FIXED FULL RANGE
	
	#Plot 1D X-axis integrated IEDF		(i.e. X-axis collapsed into 1D)					
	im3 = ax[0,1].plot(YHistData1D, YAxis1D[0:-1],  lw=2)
	cbar3 = InvisibleColourbar(ax[0,1])
	ImageOptions(fig,ax[0,1],'IEDF $f(\epsilon_{i})$ [-]','','',Legend)
	ax[0,1].set_ylim(min(YAxis1D),max(YAxis1D))									#USE FIXED FULL RANGE
#endif

#Save temporal response figure for current simulation directory
SaveString = 'Kinetics_kstep'+str(KStep).zfill(7)+ext
plt.savefig(DirKinetics+SaveString)
#			plt.show()
plt.close('all')

#=========================#
#=========================#

#Find the resonant surfaces where q_psi = mpol/ntor
def CalcRationalSurfaces(HarmonicsData3D,ntor,Threshold=0.001):

	mpol_res = len(HarmonicsData.brad[:,0,0,0])

	ResonantSurfaces = [[],[]]
	for lpsi in range(0,len(HarmonicsData.q_psi)):
		for mpol in range(0,int(mpol_res/2)):
				
			SafetyFactor = abs(HarmonicsData.q_psi[lpsi])
			try: Resonance = float(ntor)/float(mpol)
			except: Resonance = np.nan

			ResonantDifference = abs(SafetyFactor - Resonance) % 1

			if ResonantDifference < Threshold or abs(ResonantDifference-1) < Threshold:
				ResonantSurfaces[0].append(mpol)
				ResonantSurfaces[1].append(HarmonicsData.rho_pol[lpsi])
				break
			#endif
		#endfor
	#endfor

	return(ResonantSurfaces)
#enddef

#=========================#
#=========================#









