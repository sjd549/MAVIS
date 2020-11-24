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
import os, sys
import os.path
import glob

#Import additional modules
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.io import FortranFile as ff
from scipy.signal import savgol_filter
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
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
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images
#Fix "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

#List of recognized data extensions for file readin
FileExtensions = ['.hamonics','moments','.txt','.in','.nam','.dat','.out']

#Default mesh repository location (local or root)
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

#Archived variable sets
#Phys = []

####################

#Commonly Used Diagnostic Settings:
#### ASDEX ####
#electrodeloc =		[29,44] 					#Reverse [29,62]
#waveformlocs =		[[16,29],[16,44],[16,64],[0,44]]
#DOFWidth =			R;16,Z;21
#TrendLoc =			H[0];R[29,44,64,75]
#ThrustLoc =		75, 						#stdESCT=76, smlESCT=48/54
#SheathROI =		[34,72]
#SourceWidth =		[0.21]						
#Crop =				R[0.65];Z[1.0,4.0] 

####################




#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested Variables and Plotting Locations.
variables = Phys						#Requested variables to plot

radialprofiles = []						#1D Radial Profiles to be plotted (fixed Z,Phi) --
azimuthalprofiles = []					#1D Azimuthal Profiles to be plotted (fixed R,phi) |
toroidalprofiles = []					#1D Toroidal Profiles to be plotted (fixed R,Z) 
trendlocation = [] 						#Cell location For Trend Analysis [R,Z], ([] = min/max)


#Various Diagnostic Settings.
setting_seq = [-1]						#simulation seq to load		- [Int], [-1] to load last seq
setting_ntor = [-1,1]					#ntor range to plot 		- [Min,Max], [Int], [] to plot all
setting_kstep = [399]					#kstep index range to plot 	- [Min,Max], [Int], [] to plot all


#Requested diagnostics and plotting routines.
savefig_2Dequilibrium = False			#Plot 2D equilibrium figures		(xxx.harmonics)	- Working
savefig_2Dtemporal = False				#Plot 2D equilibrium movies			(xxx.harmonics)	- Working
savefig_2Dharmonics = False				#Plot 2D harmonic perturbations 	(xxx.harmonics)	- Development ?Useless?
savefig_2Dfourier = False				#Plot 2D harmonic fourier analysis	(xxx.harmonics)	- Working

savefig_2Dresponse = True				#Plot 2D plasma response 			(xxx.harmonics)	- Development

savefig_1Dtotalenergy = False			#Plot 1D total energy trends 		(xxx.energy_p)	- Working
savefig_1Dspectralenergy = False		#Plot 1D spectral energy trends 	(xxx.energy_n)	- Working


#Steady-State diagnostics terminal output toggles.
print_generaltrends = False				#Verbose Min/Max Trend Outputs						- Development


#Write processed data to ASCII files.
write_ASCII = True						#All diagnostic output written to ASCII.


#Image plotting options.
image_extension = '.png'				#Extensions ('.png', '.jpg', '.eps')
image_aspectratio = [10,10]				#[x,y] in cm 
image_radialcrop = []					#[R1,R2] in cm
image_axialcrop = []					#[Z1,Z2] in cm
image_cbarlimit = []					#[min,max] colourbar limits	

image_plotsymmetry = True				#Toggle radial symmetry
image_contourplot = True				#Toggle contour Lines in images
image_plotgrid = False					#Plot major/minor gridlines on profiles
image_plotmesh = False					#Plot material mesh outlines
image_rotate = True						#Rotate image 90 degrees to the right.

image_normalise = False					#Normalise image/profiles to local max
image_logplot = False					#Plot ln(Data), against linear axis.

#Overrides the automatic image labelling.
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

#Write a 1D midplane profile diagnostic (or just generally a profile diagnostic...)
#Include multi-var capability
#
#Write a poloidal plasma response diagnostic with mpol vs lphi (also name this...)
#
#
#
#Update Read_Harmonics functions to accept 'all' as a variable input
#This will read-in all variables and save sequentally [0,1,2,3,4] like HELENA
#
#Extract lpsi and mpol from data without having to explicitly hard-code it in within the Read_Harmonics functions
#
#Maybe also multiply everything by its normalisation factor when reading in by default?
#Make an option in the switchboard (or deep options) but this would make things much simpler...
#Tie this function to the variable label maker function (use next to eachother)
#
#Enable choice of normalised or non-normalised units in the 1D and 2D plotting routines
#Clean up the usage and definition of the unit normalisations as read-in from the MEGA file
#
#Introduce a method of selecting which time index to use for plotting
#Allow for time index, KStep or SI time (will need a function to determine this)
#
#Determine what purpose the final index has in the Harmonics data output
#HarmonicsData[kstep][mpol][ntor][lpsi][???] <--- Hard-coded to zero for ??? index for now
#
#Create an ImageExtent function to create normalised axes from data size
#


# DIAGNOSTICS TO DO
#	 	- savefig_equil 			PLOT THE REQUESTED ntor RANGE IN SINGLE FIGURE
#		- savefig_temporal 			PLOT MOVIES FOR SINGLE ntor OVER REQUESTED kstep RANGE
#		- savefig_temporal 			HOMOGONISE THE COLOURBAR FOR MOVIE PLOTS (MAKE USER OPTION)
#		- savefig_equil/temporal 	PLOT DIFFERENT TOROIDAL ANGLES?
#		- savefig_equil				CALCULATE FLUX SURFACE FUNCTION AND PLOT
#		- COMPUTE EQUILIBRIUM AND SHAPING PARAMETERS, PLOT ON LINE GRAPH OVER TIME
#		- DECIDE WHAT TO DO WITH HARMONIC DATA - KEEP 2D TEMPORALLY RESOLVED STUFF OR NOT?


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

#Takes absolute directory path name as string
#Returns list of sub-folders and list of all other contents
#Example: HomeDirFolders,HomeDirContents = DirectoryContents(os.path.abspath("."))
def DirectoryContents(AbsPath):
	#Obtain contents of supplied directory and initiate folders list
	DirContents = os.listdir( AbsPath )		#List of all contents (inc non-folders).
	DirFolders = list() 					#List of all folders.

	#Determine sub-folders within supplied directory and correct 'grammar'.
	for i in range(0,len(DirContents)):
		if os.path.isdir(AbsPath+'/'+DirContents[i]) == True:
			DirFolders.append('/'+DirContents[i]+'/')
		#endif
	#endfor

	return(DirFolders,DirContents)
#enddef

#=========================#
#=========================#

#Creates a new folder if one does not already exist.
#Takes destination dir and namestring, returns new directory.
def CreateNewFolder(Dir,DirString):
	try:
		NewFolderDir = Dir+DirString+'/'
		os.mkdir(NewFolderDir, 0755);
	except:
		a = 1
	#endtry
	return(NewFolderDir)
#enddef

#=========================#
#=========================#

#Takes folder names and returns item after requested underscore index.
#Note, index > 1 will return between two underscores, not the entire string.
def FolderNameTrimmer(DirString,Index=1):
	try:
		for i in range(0,Index):
			underscoreloc = str(DirString[::-1]).index('_')
			cutoff = (len(DirString)-underscoreloc)
			NameString = DirString[cutoff:-1]
			DirString = DirString[:cutoff-1]
		#endfor
	except:
		NameString = str(DirString[2:-1])
	#endtry

	return(NameString)
#enddef

#=========================#
#=========================#

#Takes directory list and data filename type (e.g. .png, .txt)
#Returns datalist of contents and length of datalist.
#rawdata, datalength = ExtractRawData(Dir,'.dat',l)
def ExtractRawData(Dirlist,NameString,ListIndex):
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

#Takes a 1D or 2D array and writes to a datafile in ASCII format.
#Imputs: Data, Filename, 'w'rite or 'a'ppend, and orientation (CSV or RSV).
#Example: WriteFile_ASCII(Image, "Filename", 'H', 'w')
def WriteFile_ASCII(data,filename,structure='w',Orientation='CSV'):

	#Determine dimensionality of profile.
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

	#Lowest dimention is scalar: ==> 1D array.
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

#Reads 1D or 2D data from textfile in ASCII format, returns data and header.
#Input filename, header length, data dimension and orientation (CSV or RSV).
#Example: OutputData,Header = ReadFile_ASCII('/Data.txt', 1, '2D', CSV)
def ReadFile_ASCII(Filename,HeaderIdx=0,Dimension='2D',Orientation='CSV'):
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

#Extracts MEGA seq.harmonics data shapes for use with diagnostics
#Determines if HarmonicsData is 3D or 4D and returns appropriate length array
#Inputs: HarmonicsData3D or 4D of shape [kstep][mpol][ntor][lpsi][A/B] - [kstep] optional
#Returns: DataShapes 1D list with contents: [mpol,ntor,lpsi,ltheta,kstep] - kstep optional
#Example: DataShapes = ExtractMEGA_DataShape(HarmonicsData)
def ExtractMEGA_DataShape(HarmonicsData):
	#lpsi = radial spatial resolution					#[cells]	- set by ????
	#ltheta = poloidal angle spatial resolution			#[cells]	- set by ????
	#mpol = number of poloidal harmonics considered 	#[int]		- low-pass filter limited?
	#ntor = number of toroidal harmonics considered 	#[int]		- low-pass filter limited?

	#Compare if HarmonicsData contains any of the requested variables by name (attribute)
	#N.B. Python 3.7.x objects seem to initiate with 31 generic attributes by default.
	DataAttributes = dir(HarmonicsData)
	Intersection = set(variables).intersection(set(DataAttributes))

	#Extract data array sizes for a 4D (temporally resolved) HarmonicsData object
	#Object created using: ExtractMEGA_Harmonics(Dir[l]+'data/','brad',ntor_tot)
	if len(Intersection) == 0:
		kstep_res = HarmonicsData.data.shape[0]		#kstep resolution
		mpol_res = HarmonicsData.data.shape[1]		#mpol resolution
		ntor_res = HarmonicsData.data.shape[2]		#ntor resolution
		lpsi_res = HarmonicsData.data.shape[3]		#lpsi resolution
		ltheta_res = 256							#ltheta resolution 	#!!!!!!!
	#endif

	#Extract data array sizes for a 3D HarmonicsData object
	#Object created using: ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KstepIndex,seq)
	if len(Intersection) >= 1:
		mpol_res = HarmonicsData.vrad.shape[0]		#mpol resolution
		ntor_res = HarmonicsData.vrad.shape[1]		#ntor resolution
		lpsi_res = HarmonicsData.vrad.shape[2]		#lpsi resolution
		ltheta_res = 256							#ltheta resolution 	#!!!!!!!
	#endif

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

#Extracts data resolution from HarmonicsData and poloidal axes from repository .dat files
#Inputs: HarmonicsData object of shape [mpol][ntor][lpsi][A/B], Repository directory
#Returns: GridResolutions of shape [mpol,ntor,lpsi,ltheta] and crdr, crdz poloidal axes
#Example: GridRes,crdr,crdz = ExtractMEGA_PoloidalGrid('Repository',HarmonicsData):
def ExtractMEGA_PoloidalGrid(Dir,HarmonicsData):

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

#Takes simulation folder directory (absolute path) and returns Sim128 normalisation constants
#Example: Variables,Values,Units = ReadNormConstants(Dir[l])
def ExtractMEGA_Normalisations(Dir):

	# NOTE: Duplicated variable names in output file --- ADDRESS BY SPLITTING Sim128 FILE INTO SECTIONS
	#'D beam inj. vlc.','Crit. vlc. axis','SlowD rate axis'				--- ON TOP AND POST NORM SETS
	#'psimax','major_r','rleng','left','right','zleng','raxis','zaxis'	--- ON PRE AND POST NORM SETS

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

#Reads and concatenates MEGA energy.txt output files
#Takes simulation directory (absolute path) and filename (energy_n, energy_phys)
#Returns output data and header, data of form: [Variable][Timestamp]
#Example: OutputData,Header = ExtractMEGA_Energy('LocalDataFolder/','energy_phys'])
def ExtractMEGA_Energy(Dir,Filename='energy_n'):

	#Extract Filename.txt paths for all SEQ for given data filename
	Files = sorted(glob.glob(Dir+'data/*'+Filename+'*'))

	#For each output file in the current simulation directory:
	for SEQ in range(0,len(Files)):
		#Extract header and output data for first SEQ
		if SEQ == 0:
			Header = ReadFile_ASCII(Files[SEQ],1,'2D','CSV')[1]
			OutputData = ReadFile_ASCII(Files[SEQ],1,'2D','CSV')[0]
		#Extract output data for subsequent SEQ's and append to each variable
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

#Reads MEGA xxx.harmonics output file and extracts data into 5D object.
#Data is read for a single variable [variable] over all timesteps [KStep], for a single SEQ
#Inputs: SEQ.harmonics filepath [no Root], Variable name string [Variable] as defined in DataFormat,
#Inputs: Radial mesh resolution [lpsi], poloidal mesh resolution [mpol], 
#Inputs: Total number of toroidal harmonics [ntor] including positive, negative and n=0
#Returns: Data object for requested variable with structure: Data[kstep][mpol][ntor][lpsi][Real,Complex]
#Example: HarmonicsData = Read_MEGAHarmonics('FolderName/data/001.harmonics','bphi',64,5,201]
def Read_MEGAHarmonics(Filename,Variable,mpol,ntor,lpsi,kstep=np.nan):

	#Compute flattened 3D data array length based upon mesh resolution
	n_elem = (mpol+1)*ntor*lpsi*2

	#Define FORTRANFile save data format
	#KStep (kst) is an undefined length 1D integer array		[-]
	#t (SI Time) is an undefined length 1D float array 			[1e8*IonGyroFreq*ms]
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
		Data.vrad    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Radial   Velocity Array	[-]
		Data.vtheta  = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Poloidal Velocity Array	[-]
		Data.vphi    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Toroidal Velocity Array	[-]
		Data.brad    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Radial   B-Field  Array	[-]
		Data.btheta  = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Poloidal B-Field  Array	[-]
		Data.bphi    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Toroidal B-Field  Array	[-]
		Data.erad    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Radial   E-Field  Array	[-]
		Data.etheta  = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Poloidal E-Field  Array	[-]
		Data.ephi    = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Toroidal E-Field  Array	[-]
		Data.prs     = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D ????????????????  Array	[-]
		Data.rho     = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D ????????????????  Array	[-]
		Data.dns_a   = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D ????????????????  Array	[-]
		Data.mom_a   = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D ????????????????  Array	[-]
		Data.ppara_a = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Para Pressure???  Array	[-]
		Data.pperp_a = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Perp Pressure???  Array	[-]
		Data.qpara_a = np.empty(([0,mpol+1,ntor,lpsi,2]),np.float64)	#3D Para Safety Fac?  Array	[-]
		Data.qperp_a = np.empty(([0,mpol+1,ntor+lpsi,2]),np.float64)	#3D Perp Safefy Fac?  Array	[-]
	#endif

	#Open SEQ.harmonics file and ensure data exists
	try:
		FORTRANFile = ff(Filename,'r')							#Open 001.Harmonics FORTRAN format file
		RawData = FORTRANFile.read_record(DataFormat)			#Read RawData from file in Format
	except:
		print('\n \n  Data file "'+Filename+'" not found or formatted incorrectly - check FORTRAN dtype. \n')
		FORTRANFile.close(); exit()
	#endtry

	#=====#=====#

	#Read data from each KStep for the supplied SEQ.harmonics output folder
	if np.isnan(kstep) == True:
		kstep = 0
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
			#endif

			#Extract 1D kstep and time arrays of shape: Data.kst[Real]
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

	#=====#=====#

	#Read data from each KStep for the supplied SEQ.harmonics output folder
	elif np.isnan(kstep) == False:
		index = 0
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
			#Always extract Kstep and Time arrays, until the requested Kstep
			if index >= 0:
				Data.kst     = np.append(Data.kst,  RawData['kst'][0])
				Data.time    = np.append(Data.time, RawData['t'][0])#*1e3/wa
				#Print kstep for debug purposes if requested
				if DebugMode == True: print(str(index)+'-'+str(Data.kst[index]))
			#If index matches requested kstep, retrieve data for all variables
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

#Details on the FORTRAN file format can be found below:
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.FortranFile.read_record.html
def ExtractMEGA_Harmonics(DataDir,Variable,ntor,kstep=np.nan,SEQ=0):

	#Extract harmonic output files for all SEQ in requested directory
	DataFiles = sorted(glob.glob(DataDir+"/*harm*"))

	#Extract data poloidal and toroidal resolutions 		
	lpsi,mpol = 201, 64											#!!! HARD-CODED FOR NOW !!!

	#=====#=====#	#=====#=====#

	#HarmonicsData contains single variable for all Kstep and all SEQ
	if np.isnan(kstep) == True or np.isnan(SEQ) == True:
		#Object Shape: HarmonicData [kstep][mpol][ntor][lpsi][A/B]
		HarmonicData = list()
		for i in tqdm(range(0,len(DataFiles))):
			HarmonicData.append(Read_MEGAHarmonics(DataFiles[i],Variable,mpol,ntor,lpsi))
		#endfor

		#Concatenate data from all SEQs into one continuous array for each variable within HarmonicData object
		for i in range(0,len(HarmonicData)-1):
			aux = HarmonicData.pop(1)										#'Pops' and removes 1st SEQ data array
			HarmonicData[0].kst  = np.append(HarmonicData[0].kst, aux.kst)	#Appends to zero'th SEQ data array
			HarmonicData[0].time = np.append(HarmonicData[0].time, aux.time)
			HarmonicData[0].data = np.concatenate((HarmonicData[0].data, aux.data))
			del aux															#Refresh aux array and repeat
		#endfor
		HarmonicData = HarmonicData[0]		#Replace data object with fully appended (i.e. flattened) data array

	#=====#=====#

	#HarmonicsData contains all variables at a single kstep and single SEQ
	elif isinstance(kstep, int) and isinstance(SEQ, int):							#if Variable == 'All'
		#Object Shape: HarmonicsData[mpol][ntor][lpsi][A/B]
		HarmonicData = Read_MEGAHarmonics(DataFiles[SEQ],Variable,mpol,ntor,lpsi,kstep)
	#endif

	#=====#=====#

	return(HarmonicData)
#enddef

#=========================#
#=========================#

#Reduces 3D HarmonicsData [mpol][ntor][lpsi][A/B] into 2D poloidal image [mpol,ltheta]
#Reduces 3D HarmonicsData by extracting only 1 ntor and averaging over mpol
#Also merges the real and imaginary components of HarmonicsData
#Inputs: HarmonicsData object [mpol,ntor,lpsi,A/B], 
#Inputs: VariableString matching HarmonicsData attribute, ntor Index (not mode number)
#Outputs: PoloidalData2D of shape [lpsi,ltheta]
#Example: PoloidalImage = MergePoloidal(HarmonicsData,'vrad',1)
def MergePoloidal(HarmonicsData,VariableString,ntor):		

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
		ToroidalAngle = 0											#Toroidal angle [Rad]

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

#====================================================================#
#====================================================================#






#====================================================================#
					#COMMON PLOTTING FUNCTIONS#
#====================================================================#

#Takes global inputs from switchboard, returns nothing
#Alters global image options, run before any diagnostics
#Attempts to revert matplotlib changes made in 2.0 onwards.
#See: https://matplotlib.org/users/dflt_style_changes.html
def Matplotlib_GlobalOptions():

#	mpl.style.use('classic')								#Resets to classic 1.x.x format
	
	#Image options			
	mpl.rcParams['figure.figsize'] = [10.0,10.0]			#Sets default figure size
	mpl.rcParams['figure.dpi'] = 100						#Sets viewing dpi
	mpl.rcParams['savefig.dpi'] = 100						#Sets saved dpi
	mpl.rcParams['image.interpolation'] = 'bilinear'		#Applies bilinear image 'smoothing'
	mpl.rcParams['image.resample'] = True					#Resamples data before colourmapping
	mpl.rcParams['image.cmap'] = 'plasma'					#Select global colourmap 
	#'jet','plasma','gnuplot'

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
#	mpl.rcParams['axes.prop_cycle']=cycler(color='krbgcym')	#Set Default colour rotation
	mpl.rcParams['lines.linewidth'] = 1.0					#Set Default linewidth

	#Maths and Font options
	mpl.rcParams['mathtext.fontset'] = 'cm'					#Sets 'Latex-like' maths font
	mpl.rcParams['mathtext.rm'] = 'serif'					#Sets default string font

	return()
#enddef
Matplotlib_GlobalOptions()									#Must be run before diagnostics

#=========================#
#=========================#

#Create figure and axes with variable aspect ratio, sub-plots and configurations.
#Takes image aspect ratio [x,y], number of subplots [rows, columns] and row/column sharing boolians
#Returns figure and axes seperately.
#fig,ax = figure(image_aspectratio,[1,1],shareX=False,shareY=False)
def figure(aspectratio=[],subplots=[1,1],shareX=False,shareY=False):

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
	return(fig,ax)
#enddef

#=========================#
#=========================#

#Applies plt.options to current figure based on user input.
#Returns nothing, open figure is required, use figure().
#For best results call immediately before saving/displaying figure.
#ImageOptions(fig,plt.gca(),Xlabel,Ylabel,Title,Legend)
def ImageOptions(fig,ax='NaN',Xlabel='',Ylabel='',Title='',Legend=[]):
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
		ax.set_title(Title, fontsize=14, y=1.03)
	if len(Legend) > 0:
		ax.legend(Legend, fontsize=16, frameon=False)
	#endif

	#Set labels and ticksize.
	ax.set_xlabel(Xlabel, fontsize=24)
	ax.set_ylabel(Ylabel, fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	#Force scientific notation for all axes.
	try: ax.xaxis.get_major_locator().set_params(style='sci',scilimits=(-2,3),axis='both')
	except: Fails_If_Axes_Contain_Strings = True
	#endtry

	#Set grid, default is off.
	if image_plotgrid == True: ax.grid(True)
	#endif

	#Plot mesh outline if requested.	### HACKY ###
	if image_plotmesh == True:
		mesh_auto_plot = 1 # !!AUTO PLOT MESH NOT IMPLIMENTED!! #
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

#Creates and plots a colourbar with given label and binsize.
#Takes image axis, label string, number of ticks and limits
#Allows pre-defined colourbar limits in form [min,max].
#Returns cbar axis if further changes are required.
#cbar = Colourbar(ax[0],im,'Label',5,Lim=[0,1])
def Colourbar(ax='NaN',image='NaN',Label='',Ticks=5,Lim=[]):

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
	if '\n' in [Label]: Labelpad += 25		#Pad label for multi-line names

	#Create and define colourbar axis
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	cbar = plt.colorbar(image, cax=cax)

	#Set number of ticks, label location and define scientific notation.
	cbar.set_label(Label, rotation=Rotation,labelpad=Labelpad,fontsize=LabelFontSize)
	cbar.formatter.set_scientific(True)
	cbar.formatter.set_powerlimits((-2,3))
	cbar.locator = ticker.MaxNLocator(nbins=Ticks)
	cbar.ax.yaxis.offsetText.set(size=TickFontsize)
	yticks(fontsize=TickFontsize)
	cbar.update_ticks()  

	#Apply colourbar limits if specified.  (lim=[min,max])
	if len(Lim) == 2: im.set_clim(vmin=Lim[0], vmax=Lim[1])

	return(cbar)
#enddef

#=========================#
#=========================#

#Creates an invisible colourbar to align subplots without colourbars.
#Takes image axis, returns colourbar axis if further edits are required
#cax = InvisibleColourbar(ax[0])
def InvisibleColourbar(ax='NaN'):

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

#Makeshift way of creating units for each legend entry.
#Example: VariableLegends = VariableLabelMaker(variables)
def VariableLabelMaker(variables):

	#!!! NEED TO CHECK IF LIST FIRST !!!
	#if list then do the normal procedure,
	#else make into single element list and perform procedure
	#!!! NEED TO CHECK IF LIST FIRST !!!

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
			Variable = '$r_{\psi}$'					# UNKNOWN VARIABLE
			VariableUnit = '[-]'
		elif variables[i] == 'gpsi_nrm':
			Variable = '$g_{psi}$ norm'				# UNKNOWN VARIABLE
			VariableUnit = '[-]'
		elif variables[i] == 'q_psi':
			Variable = 'Safety Factor $q_{\psi}$'
			VariableUnit = '[-]'

		#Explicit Densities
		elif variables[i] == 'rho':
			Variable = 'Plasma Density $n_{e}$'
			VariableUnit = '[m$^{-3}$]'

		#Explicit Velocity, Flux and Momentum
		elif variables[i] == 'vrad':
			Variable = 'Radial Velocity'
			VariableUnit = '[m s$^{-1}$]'
		elif variables[i] == 'vtheta':
			Variable = 'Poloidal Velocity'
			VariableUnit = '[m s$^{-1}$]'
		elif variables[i] == 'vpsi':
			Variable = 'Toroidal Velocity'
			VariableUnit = '[m s$^{-1}$]'
		elif variables[i] == 'mom_a':
			Variable = 'Fast Ion Momentum P$_{kin}$'
			VariableUnit = '[kg m s$^{-1}$]'
		elif variables[i] == 'dns_a':
			Variable = 'dns$_{a}$'					# UNKNOWN VARIABLE - KINETICS, FAST IONS?
			VariableUnit = '[-]'
		elif variables[i] == 'dns_a':
			Variable = 'dns$_{a}$'					# UNKNOWN VARIABLE - KINETICS, FAST IONS?
			VariableUnit = '[-]'
		elif variables[i] == 'ppara_a':
			Variable = 'ppara$_{a}$'				# UNKNOWN VARIABLE - KINETICS, FAST IONS?
			VariableUnit = '[-]'
		elif variables[i] == 'pperp_a':
			Variable = 'pperp$_{a}$'				# UNKNOWN VARIABLE - KINETICS, FAST IONS?
			VariableUnit = '[-]'
		elif variables[i] == 'qpara_a':
			Variable = 'qpara$_{a}$'				# UNKNOWN VARIABLE - KINETICS, FAST IONS?
			VariableUnit = '[-]'
		elif variables[i] == 'qperp_a':
			Variable = 'qperp$_{a}$'				# UNKNOWN VARIABLE - KINETICS, FAST IONS?
			VariableUnit = '[-]'

		#Explicit Electrodynamic Properties
		elif variables[i] == 'brad':
			Variable = 'Radial B-field'
			VariableUnit = '[T]'
		elif variables[i] == 'btheta':
			Variable = 'Poloidal B-field'
			VariableUnit = '[T]'
		elif variables[i] == 'bphi':
			Variable = 'Toroidal B-field'
			VariableUnit = '[T]'
		elif variables[i] == 'erad':
			Variable = 'Radial E-field'
			VariableUnit = '[V m^{-1}$]'
		elif variables[i] == 'etheta':
			Variable = 'Poloidal E-field'
			VariableUnit = '[V m^{-1}$]'
		elif variables[i] == 'ephi':
			Variable = 'Toroidal E-field'
			VariableUnit = '[V m^{-1}$]'
		#endif

		#Default if no fitting variable found.
		else:
			Variable = 'Variable'
			VariableUnit = '[Unit]'
		#endif
		VariableLegends.append(Variable+' '+VariableUnit)
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


#Compute upper and lower Alfven eigenmode threshold frequencies
#UpperThreshold,LowerThreshold = TAEThresholds(Harmonic,mpol,lpsi,qpsi,rho_pol,eps,AlfvenVelocity,subfig)
def TAEThresholds(HarmonicData,Harmonic,eps,Va,ax='NaN'):

	#eps = ????				[-]
	#Va = Alfven Velocity 	[m/s]

	#Extract required data
	data = HarmonicsData.data
	kmax = data.shape[0]
	mpol = data.shape[1]
	ntor = data.shape[2]
	lpsi = data.shape[3]

	#Initiate TAE threshold arrays
	UpperThresholds = list()
	LowerThresholds = np.zeros([lpsi,mpol-1])

	#Extract rho_pol and safety factor arrays, and initiate empty threshold arrays
	rho_pol = HarmonicsData.rho_pol
	q = abs(HarmonicsData.q_psi)
	K = np.zeros([lpsi,mpol])

	#PROVIDE SUMMARY OF MATHS AND ASSUMPTIONS 
	#PROVIDE REFERENCE FOR THESE DERIVATION(S)
	for m in range(0,mpol):
		K[:,m] = (m-abs(Harmonic)*q)/(q*R0)
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

#Determines which cells of the equilibrium are within the LCFS and which are outside
#Those inside are unchanged, while values outside are replaced with np.nan()
#Inputs are 2D equilibrium array and LCFS phi threshold (default 0.0)
#Returns 2D EquilibriumLCFS array containing $\Phi(R,Z)$ shaped as: Equilibrium[Row][Column]
def ComputeEquilLCFS(Equilibrium,Threshold=0.0):

	#Initiate required lists
	LCFSEquil = list()

	#By definition LCFS occours where flux surface equals zero
	for i in range(0,len(Equilibrium)):
		LCFSEquil.append(list())
		for j in range(0,len(Equilibrium[i])):
			if Equilibrium[i][j] >= Threshold:
				LCFSEquil[i].append(Equilibrium[i][j])
			else:
				LCFSEquil[i].append(np.nan)
			#endif
		#endfor
	#endfor

	return(LCFSEquil)
#enddef

#=========================#
#=========================#

#Takes 1D or 2D array and returns array normalised to maximum value.
#If NormFactor is defined, array will be normalised to this instead.
#Returns normalised image/profile and the max/min normalisation factors.
#NormProfile,Min,Max = Normalise(profile,NormFactor=0)
def Normalise(profile,NormFactor=0):
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
print('                                                 v0.0.8')
print('-------------------------------------------------------')
print('')
print('The following diagnostics were requested:')
print('-----------------------------------------')
if True in [savefig_1Dtotalenergy,savefig_1Dspectralenergy]:
	print('# 1D Energy Analysis')
if True in [savefig_2Dequilibrium,savefig_2Dtemporal]:
	print('# 2D Equilibrium Analysis')
if True in [savefig_2Dharmonics,savefig_2Dfourier]:
	print('# 2D Spectral Analysis')
print('-----------------------------------------')
print('')

#=====================================================================#
#=====================================================================#










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
EnergyData_phys = list()		#[3D Array] of shape Data[folder][variable][Kstep] for energy_phys.txt
EnergyData_n = list()			#[3D Array] of shape Data[folder][variable][Kstep] for energy_n.txt

#=====================================================================#
#=====================================================================#










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

		#Extract contents from 'l'th simulation folder and /data/ subfolder
		DataDir = Root+HomeDirFolders[l]+'/data/'
		DataDirContents = DirectoryContents(DataDir)[1]		#Data Folder level
		SimDirContents = SubDirContents						#Simulation Folder Level

		#Save content files from simulation folder that fit requested data output file extensions
		for j in range(0,len(SimDirContents)):
			Filename = SimDirContents[j]
			if any([x in Filename for x in FileExtensions]):
				Prefix = Dir[-1]
				DirFiles[-1].append(Prefix+Filename)
				#endif
			#endif
		#endfor

		#Save content files from /data/ subfolder that fit requested data output file extensions
		for j in range(0,len(DataDirContents)):
			Filename = DataDirContents[j]
			if any([x in Filename for x in FileExtensions]):
				Prefix = Dir[-1]+'data/'						#Note: Dir ends with a '/'
				DirFiles[-1].append(Prefix+Filename)
				#endif
			#endif
		#endfor
	else:
		#Print debug outputs to terminal if requested
		if DebugMode == True:
			print 'Discarding directory: ', HomeDirFolders[l]
		#endif
	#endif
#endfor

#If no folders detected end analysis script; else continue to analysis.
if NumFolders > 0:
	print '----------------------------------------'
	print 'Data Readin Complete, Starting Analysis:'
	print '----------------------------------------'
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
					   #GENERAL ENERGY ANALYSIS#
#====================================================================#

#====================================================================#
				 	    #TOTAL ENERGY DIAGNOSTIC#
#====================================================================#

if savefig_1Dtotalenergy == True:

	#For each detected simulation folder
	for l in tqdm(range(0,len(Dir))):
		#Create global 1D diagnostics folder and extract current simulation name
		DirEnergy = CreateNewFolder(Dir[l],'1DEnergy_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_Phys outputs and header for plotting
		#Energy_Phys: [folder][variable][timestep]
		Energy_Phys,Header_Phys = ExtractMEGA_Energy(Dir[l],'energy_phys')

		#Extract normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
#		print Variables[1],Values[1],Units[1]

		#==========##==========#

		#Create fig of desired size.
		fig,ax = figure(image_aspectratio,[3,1])

		#Define Title, Legend, Axis Labels etc...
		Title = 'Spectrally Integrated Energy Evolution for '+DirString
		Xlabel,Ylabel = 'Time [ms]', 'Energy [-]'

		#Plot total thermal, kinetic and magnetic energy over time
		ax[0].plot(Energy_Phys[1],Energy_Phys[2],'k-',lw=2)		#Kinetic
		ax[0].plot(Energy_Phys[1],Energy_Phys[3],'r-',lw=2)		#Magnetic
		ax[0].plot(Energy_Phys[1],Energy_Phys[4],'b-',lw=2)		#Thermal
		Legend = ['Kinetic','Magnetic','Thermal']		#Header_Phys[2::]
		ImageOptions(fig,ax[0],'',Ylabel,Title,Legend)

		#Plot ?co?, ?cntr? and ?total? energy over time
		ax[1].plot(Energy_Phys[1],Energy_Phys[5],'k-',lw=2)		# ?co?
		ax[1].plot(Energy_Phys[1],Energy_Phys[6],'r-',lw=2)		# ?cntr?
		ax[1].plot(Energy_Phys[1],Energy_Phys[7],'b-',lw=2)		# ?total?
		Legend = ['co','cntr','total']			#Header_Phys[2::]
		ImageOptions(fig,ax[1],'',Ylabel,'',Legend)

		#Plot ?Transferred? and ?Total? energy over time
		ax[2].plot(Energy_Phys[1],Energy_Phys[8],'k-',lw=2)		#Transferred
		ax[2].plot(Energy_Phys[1],Energy_Phys[9],'r-',lw=2)		#Total
		Legend = ['Transferred','Total']		#Header_Phys[2::]
		ImageOptions(fig,ax[2],Xlabel,Ylabel,'',Legend)

		#Save and close open figure(s)
		plt.savefig(DirEnergy+'TotalEnergy_'+SubString+ext)
#		plt.show()
		plt.close('all')
	#endfor
#endif

#====================================================================#
				 	 #SPECTRAL ENERGY DIAGNOSTIC#
#====================================================================#

if savefig_1Dspectralenergy == True:

	#For each detected simulation folder
	for l in tqdm(range(0,len(Dir))):
		#Create global 1D diagnostics folder and extract current simulation name
		DirEnergy = CreateNewFolder(Dir[l],'1DEnergy_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header for plotting
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')

		#Extract normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
#		print Variables[1],Values[1],Units[1]

		#Compute rate of change of energy for each harmonic
		DeltaEnergy_n = list()
		DeltaEnergy_n.append(Energy_n[0][1::])		#Add kstep array
		DeltaEnergy_n.append(Energy_n[1][1::])		#Add time array
		for i in range (2,len(Energy_n)):
			DeltaEnergy_n.append( list() )			#Add i'th harmonic array
			for j in range(1,len(Energy_n[i])):
				DeltaEnergy = (Energy_n[i][j]-Energy_n[i][j-1])
				DeltaTime = (Energy_n[1][j]-Energy_n[1][j-1])
				DeltaEnergy_n[i].append( DeltaEnergy/DeltaTime )
			#endfor
		#endfor

		#==========##==========#

		#Create fig of desired size.
		fig,ax = figure(image_aspectratio,[2,1])

		#Define Title, Legend, Axis Labels etc...
		Title = 'Spectrally Resolved Energy Evolution for '+DirString
		Xlabel,Ylabel = 'Time [ms]', 'Energy (Log$_{10}$) [-]'
		Legend = list()

		#Plot total energy of each harmonic component
		for i in range(2,len(Energy_n)):
			ax[0].plot(Energy_n[1],np.log10(Energy_n[i]), lw=2)
			Legend.append( 'n = '+str(i-2) )
		#endfor
		ImageOptions(fig,ax[0],'',Ylabel,Title,Legend)

		#Plot rate of change of energy of each harmonic component
		for i in range(2,len(Energy_n)):
			ax[1].plot(DeltaEnergy_n[1],np.log10(DeltaEnergy_n[i]), lw=2)
			Legend.append( 'n = '+str(i-2) )
		#endfor
		ImageOptions(fig,ax[1],Xlabel,'Delta '+Ylabel,'',Legend)

		#Save and close open figure(s)
		plt.savefig(DirEnergy+'SpectralEnergy_'+SubString+ext)
#		plt.show()
		plt.close('all')
	#endfor
#endif

#==========##==========##==========#
#==========##==========##==========#

if any([savefig_1Dtotalenergy,savefig_1Dspectralenergy]) == True:
	print '---------------------------'
	print '1D Energy Analysis Complete'
	print '---------------------------'
#endif

#====================================================================#
#====================================================================#































#====================================================================#
					 #GENERAL EQUILIBRIUM ANALYSIS#
#====================================================================#

#====================================================================#
				 	      #2D POLOIDAL PLOTS#
#====================================================================#

if savefig_2Dequilibrium == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
		print Dir[l]
		seq = setting_seq[0]				#requested SEQ file index (001 = 0)
		ntor = 0							#requested ntor mode number

		#Create global 2D diagnostics folder and extract current simulation name
		DirEquil = CreateNewFolder(Dir[l],'2DPoloidal_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		KstepArray = Energy_n[0]				#KStep Array	[-]
		TimeArray = Energy_n[1]					#Time Array		[ms]
		ntor_tot = ((len(Energy_n)-3)*2)+1		#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)	#Number of positive modes (Ignoring n=0)

		#### - CAN BE A FUNCTION
		#Create 2D array containing [ntor,ntor_index] for referencing data to be extracted
		ntor_indices = list()
		for i in range(0,ntor_tot):
			#ntor_indices contains the following [ntor, ntor_index], for positive, negative and n=0 modes
			if i-ntor_pos >= setting_ntor[0] and i-ntor_pos <= setting_ntor[1]:
				ntor_indices.append( [i-ntor_pos,i] )
			else:
				ntor_indices.append( [i-ntor_pos,np.nan] )
			#endif
		#endfor
#		print ntor_indices

		#ntor range set by ntor_indices[ntor, ntor_index], for pos, neg & n=0 modes
#		ntor = 0													#requested ntor mode number
		ntor_idx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number
		#### - CAN BE A FUNCTION

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		IonGyroFreq = Values[Variables.index('D gyro frequency')]

		#Set TimeIndex and employ to extract KStep and Time		--- LINK THIS TO SWITCHBOARD
		KstepIndex = setting_kstep[0]
		Kstep = KstepArray[KstepIndex]			#[-]
		Time = TimeArray[KstepIndex]			#[ms]

		#Extract Harmonics outputs for plotting, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.data: [4D Array] of shape [mpol][ntor][lpsi][A/B] for all variables
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KstepIndex,seq)

		#Extract data resolution and poloidal axes from repository .dat files
		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		DataShape = ExtractMEGA_DataShape(HarmonicsData)
		Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)

		#Create Variablelabels
		VariableLabels = VariableLabelMaker(variables)

		#For each requested variable
		for i in tqdm(range(0,len(variables))):

			#Select variable and Merge 3D Data into 2D poloidal slice
			Image = MergePoloidal(HarmonicsData,variables[i],ntor_idx)	

			#Create figure and define Title, Legend, Axis Labels etc...
			fig,ax = figure(image_aspectratio,1)
			Title = VariableLabels[i]+', n='+str(ntor)+', t='+str(round(Time,3))+' ms \n Simulation: '+DirString
			Xlabel,Ylabel = 'Radius $R$ [m]', 'Height $Z$ [m]'
			Legend = list()

			#Plot 2D harmonics figure and beautify
			im = ax.contourf(Crdr, Crdz, Image, 100); plt.axis('scaled')
			cbar = Colourbar(ax,im,VariableLabels[i],5)
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

			#Save 2D harmonics figure for current simulation
			plt.savefig(DirEquil+variables[i]+'_n'+str(ntor)+'_t='+str(round(Time,3))+ext)
#			plt.show()
			plt.close('all')
		#endfor
	#endfor
#endif

#====================================================================#
				 	      #2D POLOIDAL MOVIES#
#====================================================================#

if savefig_2Dtemporal == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
		print Dir[l]
		seq = setting_seq[0]				#requested SEQ file index (001 = 0)
		ntor = 0							#requested ntor mode number

		#Create global 2D diagnostics folder and extract current simulation name
		DirEquil = CreateNewFolder(Dir[l],'2DPoloidal_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		KstepArray = Energy_n[0]				#KStep Array	[-]
		TimeArray = Energy_n[1]					#Time Array		[ms]
		ntor_tot = ((len(Energy_n)-3)*2)+1		#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)	#Number of positive modes (Ignoring n=0)

		#### - CAN BE A FUNCTION
		#Create 2D array containing [ntor,ntor_index] for referencing data to be extracted
		ntor_indices = list()
		for i in range(0,ntor_tot):
			#ntor_indices contains the following [ntor, ntor_index], for positive, negative and n=0 modes
			if i-ntor_pos >= setting_ntor[0] and i-ntor_pos <= setting_ntor[1]:
				ntor_indices.append( [i-ntor_pos,i] )
			else:
				ntor_indices.append( [i-ntor_pos,np.nan] )
			#endif
		#endfor
#		print ntor_indices

		#ntor range set by ntor_indices[ntor, ntor_index], for pos, neg & n=0 modes
#		ntor = 0													#requested ntor mode number
		ntor_idx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number
		#### - CAN BE A FUNCTION

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		IonGyroFreq = Values[Variables.index('D gyro frequency')]

		#Extract Variablelabels
		VariableLabels = VariableLabelMaker(variables)

		#Extract and plot data for each timestep
		for i in tqdm(range(0,len(KstepArray))):

			#Set TimeIndex and employ to extract KStep and Time
			KstepIndex = i
			KStep = KstepArray[KstepIndex]			#[-]
			Time = TimeArray[KstepIndex]			#[ms]

			#Extract Harmonics outputs for plotting, it contains:
			#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
			#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
			#HarmonicsData.data: [3D Array] of shape [mpol][ntor][lpsi][A/B] for all variables
			HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KstepIndex,seq)

			#Extract data resolution and poloidal axes from repository .dat files
			#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
			DataShape = ExtractMEGA_DataShape(HarmonicsData)
			Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)

			#For each requested variable at the current Kstep
			for j in range(0,len(variables)):

				#Create global 2D diagnostics folder and extract current simulation name
				DirMovie = CreateNewFolder(DirEquil,variables[j]+'_n'+str(ntor)+'_Movie/')

				#Select variable and Merge 3D Data into 2D poloidal slice
				Image = MergePoloidal(HarmonicsData,variables[i],ntor_idx)

				#Create figure and define Title, Legend, Axis Labels etc...
				fig,ax = figure(image_aspectratio,1)
				Title = VariableLabels[j]+', n='+str(ntor)+', t='+str(round(Time,3))+' ms \n Simulation: '+DirString
				Xlabel,Ylabel = 'Radius $R$ [m]', 'Height $Z$ [m]'
				Legend = list()

				#Plot 2D harmonics figure and beautify
				im = ax.contourf(Crdr, Crdz, Image, 100); plt.axis('scaled')
				cbar = Colourbar(ax,im,VariableLabels[j],5)
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

				#Save 2D harmonics figure for current simulation
				plt.savefig(DirMovie+variables[j]+'_n'+str(ntor)+'_kstep'+str(Kstep)+ext)
#				plt.show()
				plt.close('all')
			#endfor
		#endfor
	#endfor
#endif


#==========##==========##==========#
#==========##==========##==========#

if any([savefig_2Dequilibrium,savefig_2Dtemporal]) == True:
	print '--------------------------------'
	print '2D Equilibrium Analysis Complete'
	print '--------------------------------'
#endif

#====================================================================#
#====================================================================#





















#====================================================================#
					#GENERAL PLASMA RESPONSE ANALYSIS#
#====================================================================#

#====================================================================#
				 	      #PLASMA RESPONSE PLOTS#
#====================================================================#

if savefig_2Dresponse == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
		print Dir[l]
		seq = setting_seq[0]				#requested SEQ file index (001 = 0)
		ntor = 1							#requested ntor mode number

		#Create global 2D diagnostics folder and extract current simulation name
		DirHarmonics = CreateNewFolder(Dir[l],'2DHarmonic_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		KstepArray = Energy_n[0]				#KStep Array	[-]
		TimeArray = Energy_n[1]					#Time Array		[ms]
		ntor_tot = ((len(Energy_n)-3)*2)+1		#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)	#Number of positive modes (Ignoring n=0)
		ntor0 = int(ceil(ntor_tot/2))			#ntor = 0, baseline equilibrium data

		#### - CAN BE A FUNCTION
		#Create 2D array containing [ntor,ntor_index] for referencing data to be extracted
		ntor_indices = list()
		for i in range(0,ntor_tot):
			#ntor_indices contains the following [ntor, ntor_index], for positive, negative and n=0 modes
			if i-ntor_pos >= setting_ntor[0] and i-ntor_pos <= setting_ntor[1]:
				ntor_indices.append( [i-ntor_pos,i] )
			else:
				ntor_indices.append( [i-ntor_pos,np.nan] )
			#endif
		#endfor
#		print ntor_indices

		#ntor range set by ntor_indices[ntor, ntor_index], for pos, neg & n=0 modes
#		ntor = 1													#requested ntor mode number
		ntor_idx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number
		#### - CAN BE A FUNCTION

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		AlfvenVelocity = Values[Variables.index('Alfven velocity')] #B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = Values[Variables.index('D gyro frequency')]
		IonDensity = Values[Variables.index('Bulk density')]
		B0 = Values[Variables.index('Mag.fld. at axis')]
		R0 = Values[Variables.index('raxis')]
		m_D = 3.34e-27
		eps = 0.5/R0

		#Extract Variablelabels
		VariableLabels = VariableLabelMaker(variables)

		#Extract and plot data for each timestep
		for i in range(398,399):

			#Set TimeIndex and employ to extract KStep and Time
			KstepIndex = i
			Kstep = KstepArray[KstepIndex]			#[-]
			Time = TimeArray[KstepIndex]			#[ms]

			#Extract Harmonics outputs for plotting, it contains:
			#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
			#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
			#HarmonicsData.data: [3D Array] of shape [mpol][ntor][lpsi][A/B] for all variables
			HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KstepIndex,seq)

			#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
			DataShape = ExtractMEGA_DataShape(HarmonicsData)#; print DataShape
			mpol_res = DataShape[0]; ntor_res = DataShape[1]
			lpsi_res = DataShape[2]; ltheta_res = DataShape[3]

			#Extract radial magnetic field (brad) from SEQ.harmonic object
			#Data is of shape: Data[mpol,ntor,lpsi,A/B]
			Data = getattr(HarmonicsData, 'brad')		#NEED TO IMPLIMENT DIFFERENT VARIABLES, vrad for example

			#Combine spectral components A and B in quadrature to obtain the B-field amplitude
			#Pos corresponds to resonant poloidal modes 	i.e. +m on RHS of image
			#Neg corresponds to non-resonant poloidal modes i.e. -m on LHS of image
			#One of ntor_pos-ntor or ntor_pos+ntor will equal 0, the equil values.
			DataAmpPos = np.sqrt( (Data[:, ntor_pos-ntor,:,0]**2) + (Data[:, ntor_pos-ntor,:,1]**2) )
			DataAmpNeg = np.sqrt( (Data[:, ntor_pos+ntor,:,0]**2) + (Data[:, ntor_pos+ntor,:,1]**2) )
			DataAmpNeg = np.flip( DataAmpNeg,axis=0)		#Flip LHS of image for plotting

			#Concat positive and negative ntor to obtain full poloidal harmonic spectrum
			#DataAmp is of shape: DataAmp[lpsi,mpol]
			DataAmp = np.concatenate((DataAmpNeg,DataAmpPos[1:,:]),axis=0)

			#Create Image array and Axes, rotate such that mpol spectrum is on X-axis.
			#Image is of shape: [mpol,lpsi] 
			Image = DataAmp.transpose()*B0*1e4								#[G]? - must use Brad
			Xaxis =	[x-int(mpol_res-1) for x in range(0,2*mpol_res-1,1)]	#Poloidal Mode Numbers
			Yaxis = HarmonicsData.rho_pol									#Radial Location

			#Create figure and define Title, Legend, Axis Labels etc...
			fig,ax = figure(image_aspectratio,1)
			Title = 'Plasma Response: n='+str(ntor)+', t='+str(round(Time,3))+' [ms] \n Simulation: '+DirString
			Xlabel,Ylabel = 'Poloidal Harmonic $m_{pol}$ [-]', 'Normalised Minor Radius $\\rho_{pol}$ [-]'
			Legend = list()

			#Plot figure
			im = ax.contourf(Xaxis, Yaxis, Image, 50)
			res = ax.plot(-ntor*HarmonicsData.q_psi,HarmonicsData.rho_pol, 'w--', lw=2)
			cbar = Colourbar(ax,im,'Radial Magnetic Perturbation B$_{rad}$ [G]',5)
			#####
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
			ax.set_xlim(-16,16)

			#Save 2D harmonics figure for current simulation
			plt.savefig(DirHarmonics+'PlasmaResponse_n'+str(ntor)+'_kstep='+str(Kstep)+ext)
			plt.show()
			plt.close('all')
		#endfor
	#endfor
#endif



#==========##==========##==========#
#==========##==========##==========#

if any([savefig_2Dresponse]) == True:
	print '------------------------------------'
	print '2D Plasma Response Analysis Complete'
	print '------------------------------------'
#endif

#====================================================================#
#====================================================================#







































#====================================================================#
				  #GENERAL SPECTRAL & HARMONIC ANALYSIS#
#====================================================================#

#====================================================================#
				 	      #2D HARMONIC PLOTS#
#====================================================================#

if savefig_2Dharmonics == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		print Dir[l]

		#Create global 2D diagnostics folder and extract current simulation name
		DirHarmonics = CreateNewFolder(Dir[l],'2DHarmonic_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		KstepArray = Energy_n[0]				#KStep Array	[-]
		TimeArray = Energy_n[1]					#Time Array		[ms]
		ntor_tot = ((len(Energy_n)-3)*2)+1		#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)	#Number of positive modes (Ignoring n=0)

		#### - CAN BE A FUNCTION
		#Create 2D array containing [ntor,ntor_index] for referencing data to be extracted
		ntor_indices = list()
		for i in range(0,ntor_tot):
			if i-ntor_pos >= setting_ntor[0] and i-ntor_pos <= setting_ntor[1]:
				#ntor_indices contains the following [ntor, ntor_index], for positive, negative and n=0 modes
				ntor_indices.append( [i-ntor_pos,i] )
			#endif
		#endfor
#		print ntor_indices

		#ntor range set by ntor_indices[ntor, ntor_index], for pos, neg & n=0 modes
		ntor = 1													#requested ntor mode number
		ntor_idx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number
		#### - CAN BE A FUNCTION

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		IonGyroFreq = Values[Variables.index('D gyro frequency')]

		#Extract Variablelabels
		VariableLabels = VariableLabelMaker(variables)


		#For each requested variable
		for i in range(0,len(variables)):

			#Create new folder for each variable
			DirVariable = CreateNewFolder(DirHarmonics,variables[i])

			#Extract Harmonics outputs for plotting, it contains:
			#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
			#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
			#HarmonicsData.data: [4D Array] of shape [kstep][mpol][ntor][lpsi][A/B] for variables[i]
			HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/',variables[i],ntor_tot)
			HarmonicsData.data = HarmonicsData.data[:,:,:,:,0]		#Extract Real part of HarmonicsData

			#Set TimeIndex and employ to extract KStep and Time ----- TimeIndex should be cleaner
			TimeIndex = setting_kstep[0]								
			KStep = KstepArray.kst[TimeIndex]
			Time = TimeArray[TimeIndex]

			#Extract normalised axes and set image extent
			## !!! MAKE THESE A FUNCTION extent = ImageExtent() !!! ##		
			rho_pol = HarmonicsData.rho_pol							#Normalised poloidal angle
			phi_tor = range(1,shape(HarmonicsData.data)[1]+1)		#Array of toroidal cells
			phi_tor = [x/float(len(phi_tor)) for x in phi_tor]		#Normalised toroidal angle
			extent = [rho_pol[0],rho_pol[-1], phi_tor[0],phi_tor[-1]]
			
			#Plot images of variables[i] for all requested ntor
			#ntor range set by ntor_indices[ntor, ntor_index], for pos, neg & n=0 modes
			for j in range(0,ntor_indices[-1][1]):

				#Extract image from HarmonicsData
				ntor = ntor_indices[j][0]						#ntor mode number
				ntor_idx = ntor_indices[j][1]					#index referring to ntor mode number
				Image = HarmonicsData.data[TimeIndex,:,ntor_idx,:]

				#Create figure and define Title, Legend, Axis Labels etc...
				fig,ax = figure(image_aspectratio,1)
				Title = VariableLabels[j]+', ntor='+str(ntor)+', t='+str(Time)+' \n Simulation: '+DirString
				Xlabel,Ylabel = 'Poloidal Extent $\\rho_{pol}$ [-]', 'Toroidal Extent $\phi_{tor}$ [-]'
#				Xlabel,Ylabel = 'Poloidal Resolution $m_{\\theta}$ [cells]', 'Toroidal Resolution $l_{\phi}$ [cells]'
				Legend = list()

				#Plot 2D harmonics figure and beautify
				im = ax.imshow(Image, extent=extent, aspect='auto', origin='bottom')
				cbar = Colourbar(ax,im,VariableLabels[j]+' [-]',5)
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

				#Save 2D harmonics figure for current simulation
				plt.savefig(DirVariable+variables[j]+'_n'+str(ntor)+ext)
#				plt.show()
				plt.close('all')
			#endfor
		#endfor
	#endfor
#endif

#====================================================================#
				 	    #2D FOURIER ANALYSIS#
#====================================================================#

if savefig_2Dfourier == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		print Dir[l]

		#Create global 2D diagnostics folder and extract current simulation name
		DirHarmonics = CreateNewFolder(Dir[l],'2DHarmonic_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		KstepArray = Energy_n[0]				#KStep Array	[-]
		TimeArray = Energy_n[1]					#Time Array		[ms]
		ntor_tot = ((len(Energy_n)-3)*2)+1		#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)	#Number of positive modes (Ignoring n=0)

		#Create 2D array containing [ntor,ntor_index] for referencing data to be extracted
		ntor_indices = list()
		for i in range(0,ntor_tot):
			if i-ntor_pos >= setting_ntor[0] and i-ntor_pos <= setting_ntor[1]:
				#ntor_indices contains the following [ntor, ntor_index], for positive, negative and n=0 modes
				ntor_indices.append( [i-ntor_pos,i] )
			#endif
		#endfor
#		print ntor_indices

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		AlfvenVelocity = Values[Variables.index('Alfven velocity')] #B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = Values[Variables.index('D gyro frequency')]
		IonDensity = Values[Variables.index('Bulk density')]
		B0 = Values[Variables.index('Mag.fld. at axis')]
		R0 = Values[Variables.index('raxis')]
		m_D = 3.34e-27
		eps = 0.5/R0

		#Extract Harmonics outputs for plotting, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.data: [4D Array] of shape [kstep][mpol][ntor][lpsi][A/B] for variables[i]
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','bphi',ntor_tot)
		Data = HarmonicsData.data

		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		DataShape = ExtractMEGA_DataShape(HarmonicsData)
		mpol, ntor = DataShape[0], DataShape[1]
		lpsi, ltheta = DataShape[2], DataShape[3]
		kmax, dt = DataShape[4], (TimeArray[1]-TimeArray[0])
		ntor2 = int(0.5*(ntor-1))

		#Extract Variablelabel
		VariableLabel = 'BField Perturbation $\delta B_{\phi}$ [T]'		#Needs a seperate function


		#BELOW TO STILL BE TRANSLATED
		#ALSO NEED TO ADD COLOURBAR TO THE FOURIER PLOTS!!!
		print kmax, mpol, ntor, lpsi, ntor2


		#Create Vcos list and add zeros equal to the toroidal resolution (Why tho?)
		vcos = list()
		for n in range(0,ntor2):
			vcos.append( np.zeros([]) )
			#Extract Bpol data for each respective toroidal (n) and poloidal (m) cells (adding vcos zeros?)
			for m in range(0,mpol):
				vcos[n] = vcos[n] + Data[:,m,n,:,0]	#data structure: [kstep][mpol][ntor][lpsi][A/B] 
			#endfor
			vcos[n][np.isnan(vcos[n])] = 0			#Remove any NaNs
		#endfor

		vcos_fft,vcos_len = list(),list()
		#Extract fourier components from vcos
		for n in range(0,ntor2):
			vcos_fft.append( np.fft.fft(vcos[n],axis=0) )	# 												???
		  	vcos_fft[n][0,:] = vcos_fft[n][0,:]*0.0			#Discard imaginary components 					???
			vcos_len.append( int(len(vcos_fft[n])/2)-1 )	#Determine lowpass filter frequency threshold 	???
			vcos_fft[n] = vcos_fft[n][0:vcos_len[n],:]		#Discard upper frequencies (lowpass filter)
		#endfor

		#==========##==========#

		#Create fig of desired size - increasing Xlim with the number of harmonics
		Xlim,Ylim = int(10*(float(ntor_pos)/1.75)), 12
		fig,ax = figure([Xlim,Ylim],[2,ntor2])

		#For each toroidal harmonic:
		for i in range(0,ntor2):

			#Temporal evolution plotted on the top row (row 0)
			if ntor2 == 1: subfig = ax[0]
			elif ntor2 > 1: subfig = ax[0,i]
			Harmonic = -ntor2+i									# Why is ntor reversed?
			#Construct figure axes and meshgrid (not used)
			Xaxis = HarmonicsData.rho_pol						#[-]
			Yaxis = (HarmonicsData.time/IonGyroFreq)*1e3		#[ms]
			extent = [Xaxis[0],Xaxis[-1], Yaxis[0],Yaxis[-1]]
			X,Y = np.meshgrid(Xaxis,Yaxis)						#im = subfig.contourf(X,Y,vcos[i])

			#Plot harmonic temporal evolution
			im = subfig.imshow(vcos[i][::-1], extent=extent, aspect='auto')
			co = subfig.contour(vcos[i], extent=extent, levels=10)
			
			#Add colourbar and beautify plot - taking account of panel location
			if i == 0 and ntor2 > 1: 					#If first panel with more panels to right
				ImageOptions(fig,subfig,'','Time [ms]','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
			elif i == 0 and ntor2 == 1:					#If first panel with no panels to right
				cbar = Colourbar(subfig,im,VariableLabel,5)
				ImageOptions(fig,subfig,'','Time [ms]','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
			elif i > 0 and i < ntor2-1: 				#If middle panel with more panels to right
				ImageOptions(fig,subfig,'','','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
 				im.axes.get_yaxis().set_visible(False)
			elif i == ntor2-1 and ntor2 > 1:			#If last panel with more panels to left
				cbar = Colourbar(subfig,im,VariableLabel,5)
				ImageOptions(fig,subfig,'','','n='+str(Harmonic),'')
 				im.axes.get_xaxis().set_visible(False)
 				im.axes.get_yaxis().set_visible(False)
			#endif

			#==========#

			#Fourier analysis plotted on the bottom row (row 1)
			if ntor2 == 1: subfig = ax[1]
			elif ntor2 > 1: subfig = ax[1,i]
			#Construct figure axes and meshgrid (not used)
			Xaxis = HarmonicsData.rho_pol						#[-]
			Yaxis = np.linspace(0,0.5/dt,vcos_len[i])			#[kHz]
			extent = [Xaxis[0],Xaxis[-1], Yaxis[0],Yaxis[-1]]
			X,Y = np.meshgrid(Xaxis,Yaxis)						#im = subfig.contourf(X,Y,vcos_fft[i])

			#Plot fourier amplitude spectrum
			im = subfig.imshow(real(vcos_fft[i])[::-1], extent=extent, aspect='auto')
			co = subfig.contour(real(vcos_fft[i]), extent=extent, levels=10)

			#Add colourbar and beautify plot - taking account of panel location
			if i == 0 and ntor2 > 1: 					#If first panel with more panels to right
				ImageOptions(fig,subfig,'Poloidal Extent $\\rho_{pol}$','Frequency [kHz]','','')
			elif i == 0 and ntor2 == 1:					#If first panel with no panels to right
				cbar = Colourbar(subfig,im,VariableLabel,5)
				ImageOptions(fig,subfig,'Poloidal Extent $\\rho_{pol}$','Frequency [kHz]','','')
			elif i > 0 and i < ntor2-1:   				#If middle panel with more panels to right
				ImageOptions(fig,subfig,'Poloidal Extent $\\rho_{pol}$','','','')
 				im.axes.get_yaxis().set_visible(False)
			elif i == ntor2-1 and ntor2 > 1:			#If last panel with more panels to left
				cbar = Colourbar(subfig,im,VariableLabel,5)
				ImageOptions(fig,subfig,'Poloidal Extent $\\rho_{pol}$','','','')
 				im.axes.get_yaxis().set_visible(False)
			#endif
			subfig.set_ylim([0,200])

			#Compute and plot TAE thresholds
			UpperThresholds,LowerThresholds = TAEThresholds(HarmonicsData,Harmonic,eps,AlfvenVelocity,subfig)
		#endfor

		#Save 2D harmonics figure for current simulation
		plt.savefig(DirHarmonics+'FourierAnalysis_'+SubString+ext)
		plt.show()
		plt.close('all')
	#endfor
#endif

#==========##==========##==========#
#==========##==========##==========#

if any([savefig_2Dharmonics,savefig_2Dfourier]) == True:
	print '-----------------------------'
	print '2D Spectral Analysis Complete'
	print '-----------------------------'
#endif

#====================================================================#
#====================================================================#



















































#====================================================================#
							# CODE DUMP #
#====================================================================#

# UNUSED SNIPPITS OF CODE WILL BE PUT IN HERE.
# MOST CODE IS BOUND FOR REMOVAL BUT SOME MAY STILL HAVE USE.


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









#====================================================================#
				 	 #READING DATA INTO MEMORY#
#====================================================================#

if True == False:
	print'-----------------------'
	print'Beginning Data Read-in.'
	print'-----------------------'

	#Extraction and organization of data from .txt files.
	for l in tqdm(range(0,numfolders)):

		#Extraction and organization of data from .txt files. - NOT CURRENTLY USED
		if True == False:
	 
			#energy_n: [folder][variable][timestep]
			Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
			#Energy_Phys: [folder][variable][timestep]
			Energy_Phys,Header_Phys = ExtractMEGA_Energy(Dir[l],'Energy_Phys')
		#endif

	#===================##===================#
	#===================##===================#

		#Extraction and organization of data from .harmonics files. - NOT CURRENTLY USED
		if True == False:

			#Load data from TECPLOT_KIN file and unpack into 1D array.
			rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
			rawdata_kin.append(rawdata)
		#endif


	#===================##===================#
	#===================##===================#

		#Extraction and organization of data from .moments files. - NOT CURRENTLY USED
		if True == False:

			#Load data from TECPLOT_KIN file and unpack into 1D array.
			rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
			rawdata_kin.append(rawdata)
		#endif


	#===================##===================#
	#===================##===================#

	#Empty and delete any non-global data lists.
	HomeDir,DirContents = list(),list()
	del HomeDir,DirContents

	#Alert user that readin process has ended and continue with selected diagnostics.
	if any([savefig_1Dtotalenergy,savefig_1Dspectralenergy]) == True:
		print '----------------------------------------'
		print 'Data Readin Complete, Starting Analysis:'
		print '----------------------------------------'
	else:
		print '------------------'
		print 'Analysis Complete.'
		print '------------------'
		exit()
	#endif
#endif

#=====================================================================#
#=====================================================================#















