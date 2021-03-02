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
Kin = ['R_gc','Z_gc','Phi_gc','p_gc','pphi_gc','mu_gc','E_gc','Lambda_gc','psip','ffff']

#Archived variable sets
#Phys = []

####################

#Common Diagnostic Settings:

#===== AUG#34570 =====#
#radialprofiles = [90,120,130]
#poloidalprofiles = []
#toroidalprofiles = []
#trendlocation = []

#setting_seq = [0,0]
#setting_ntor = [0,2]
#setting_kstep = [002,003]

####################




#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested Variables and Plotting Locations:
variables = ['brad']	#Phys			#Requested variables to plot			#Phys+Kin

radialprofiles = [90,120]				#1D Radial Profiles (fixed theta, phi) :: Poloidal Angle [deg]
poloidalprofiles = []					#1D Poloidal Profiles (fixed rho_pol, phi) :: Norm. Radius [-]
toroidalprofiles = []					#1D Toroidal Profiles (fixed rho_pol, theta) :: Toroidal Angle [deg]
trendlocation = [] 						#Cell location For Trend Analysis [R,theta,phi], ([] = min/max)

#Various Diagnostic Settings:
setting_seq = [0,0]						#Simulation seq to load		- [Min,Max], [Int], [-1] to load last
setting_ntor = [0,2]					#ntor range to plot 		- [Min,Max], [Int], [] to plot all
setting_kstep = [30,399,1]				#kstep index range to plot 	- [Min,Max,Step], [Int], [] to plot all


#Requested diagnostics and plotting routines:
savefig_1Dtotalenergy = False			#Plot 1D total energy trends 		(xxx.energy_p)	- Working
savefig_1Dspectralenergy = False			#Plot 1D spectral energy trends 	(xxx.energy_n)	- Working

savefig_1Dequilibrium = True			#Plot 1D equilibrium profiles		(xxx.harmonics) - Working
savefig_2Dequilibrium = False			#Plot 2D equilibrium figures		(xxx.harmonics)	- Working
savefig_2Dequilmovie = False			#Plot 2D equilibrium movies			(xxx.harmonics)	- Working
savefig_2Dresponse = False				#Plot 2D plasma response 			(xxx.harmonics)	- Working

savefig_2Dharmonics = False				#Plot 2D harmonic continnum		 	(xxx.harmonics)	- Working

savefig_1Dkinetics = False				#Plot 2D kinetic distributions	 	(gc_a_kstepxxx)	- Working
savefig_2Dkinetics = False				#Plot 1D kinetic distributions	 	(gc_a_kstepxxx)	- Working


#Requested diagnostic terminal outputs:
print_generaltrends = False				#Verbose Min/Max Trend Outputs						- In Development


#Write processed data to ASCII files:
write_ASCII = True						#All diagnostic output written to ASCII.


#Image plotting options:
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

#Overrides the automatic image labelling:
titleoverride = []
legendoverride = []
xaxisoverride = []
xlabeloverride = []
ylabeloverride = []
cbaroverride = []

#=====================================================================#
#=====================================================================#


#					TO DO FOR TODAY:
#
#Combine all equil plots into a folder [Equil Plots] -->  1DProfiles -- > Radial
#																		  Poloidal
#
#														  2DProfiles -- > Poloidal
#
#Get all SEQ, Ntor loops working as they should - Functional Switchboard
#
#FIX ISSUE WHERE "outputdata is referenced before assignment" IF FILENAME HAS A SPACE IN IT
#
#REMOVE RELIANCE ON ENERGY_N FOR KSTEP AND TIME AXES... WON'T WORK IF KHARM != KCHK
#
#MAKE ALL 'BACKEND' INPUTS INTO FUNCTIONS - ENERGY KSTEP READIN, ETC...
#
#FINISH MODE GROWTH RATE DIAGNOSTIC - V1 JUST NEEDS TO PLOT 'SOMETHING'.
#
#FINISH 1D PLOTTING ROUTINE WITH POLOIDAL EXTRACTION FUNCTION
#
#FINISH PLASMA RESPONSE ROUTINE WITH POLOIDALLY RESOLVED VERSION
#




#============================#
#        ####TODO####        #
#============================#

#REMOVE RELIANCE ON ENERGY_N FOR KSTEP AND TIME AXES... WON'T WORK IF KHARM != KCHK
#
#FIX ISSUE WHERE "outputdata is referenced before assignment" IF FILENAME HAS A SPACE IN IT
#
#FINISH 1D PLOTTING ROUTINE WITH POLOIDAL EXTRACTION FUNCTION
#
#FINISH PLASMA RESPONSE ROUTINE WITH POLOIDALLY RESOLVED VERSION
#
#MAKE ALL 'BACKEND' INPUTS INTO FUNCTIONS - ENERGY KSTEP READIN, ETC...
#
#ADD INFORMATION REGARDING SIMULATION RESOLUTION SEQ, KSTEP, ETC... BELOW MAVIS SPLASH
#
#ADD DIAGNOSTIC COMPARING ENERGY CONVERGENCE BETWEEN MULTIPLE SIMULATIONS  (FINISH GAMMA DIAGNOSTIC)
#
#ADD ABILITY TO SAVE ASCII DATA FOR 2D IMAGES (PARTICULARILY PLASMA RESPONSE)
#
#ADD OPTION TO HOMOGONISE THE COLOURBAR FOR EQUIL/RESPONSE/KINETIC MOVIES
#
##CHECK NORMALISATION FOR ALL EXISTING DIAGNOSTICS AND DETAIL IN READ-IN FUNCTION FOR LATER
#ADD DE-NORMALISATION FUNCTION WITH OPTION TO APPLY AT READ-IN STAGE 
#Maybe also multiply everything by its normalisation factor when reading in by default?
#Make an option in the switchboard (or deep options) but this would make things much simpler...
#Tie this function to the variable label maker function (use next to eachother)
#
#Enable choice of normalised or non-normalised units in the 1D and 2D plotting routines
#Clean up the usage and definition of the unit normalisations as read-in from the MEGA file
#
#Extract lpsi and mpol from data without having to explicitly hard-code it (see Read_Harmonics function)
#
#Create an ImageExtent function to create normalised axes from data size
#
# savefig_equilibrium/equilmovie 	OPTION TO PLOT DIFFERENT TOROIDAL ANGLES?
#
# savefig_equilibrium				OPTION TO PLOT 2D FLUX SURFACE FUNCTION?


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

def ExtractRawData(Dirlist,NameString,ListIndex):
#Takes directory list and data filename type (e.g. .png, .txt)
#Returns datalist of contents and length of datalist.
#rawdata, datalength = ExtractRawData(Dir,'.dat',l)

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
#Example: WriteFile_ASCII(Image, "Filename", 'H', 'w')

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

def ExtractMEGA_DataShape(HarmonicsData):
#Extracts MEGA seq.harmonics data shapes for use with diagnostics
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

	#Extract data array sizes for a 3D (single kstep) HarmonicsData object
	#NOTE, this assumes brad is an attribute within HarmonicsData - need a 'generic' attribute tag if possible
	#Object created using: ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KstepIndex,seq)
	if len(Intersection) >= 1:
		mpol_res = HarmonicsData.brad.shape[0]		#mpol resolution
		ntor_res = HarmonicsData.brad.shape[1]		#ntor resolution
		lpsi_res = HarmonicsData.brad.shape[2]		#lpsi resolution
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

def ExtractMEGA_PoloidalGrid(Dir,HarmonicsData):
#Extracts data resolution from HarmonicsData and poloidal axes from repository .dat files
#Inputs: HarmonicsData object of shape [mpol][ntor][lpsi][A/B], Repository directory
#Returns: GridResolutions of shape [mpol,ntor,lpsi,ltheta] and crdr, crdz poloidal axes
#Example: GridRes,crdr,crdz = ExtractMEGA_PoloidalGrid('Repository',HarmonicsData):

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
				Data.time    = np.append(Data.time, RawData['t'][0])#*1e3/wa
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
			MaxKstep = (setting_seq[0]+1)*index*(Data.kst[1]-Data.kst[0])
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

def ExtractMEGA_Harmonics(DataDir,Variable,ntor,kstep=np.nan,SEQ=0,Dimension='1D'):
#Details on the FORTRAN file format can be found below:
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.FortranFile.read_record.html

	#Extract harmonic output files for all SEQ in requested directory
	DataFiles = sorted(glob.glob(DataDir+"/*harm*"))

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
		Variable = 'NaN'		#Variable is a dummy input here - Not used.
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
			pop = HarmonicData.pop(1)										#'Pops' and removes 1st SEQ data array
			HarmonicData[0].kst  = np.append(HarmonicData[0].kst, pop.kst)	#Appends to zero'th SEQ data array
			HarmonicData[0].time = np.append(HarmonicData[0].time, pop.time)
			HarmonicData[0].data = np.concatenate((HarmonicData[0].data, pop.data))
			del pop															#Refresh pop array and repeat
		#endfor
		HarmonicData = HarmonicData[0]		#Replace data object with fully appended (i.e. flattened) data array

	#=====#=====#	#=====#=====#

	DenormaliseAtReadin = True
	#De-normalise data if requested
	if DenormaliseAtReadin == True:

		#Remove '/data/' from directory --> Dir now points to simulation root folder
		#		  Reverse       split into 2    Keep preamble  re-reverse   +'/' on end
		NormDir = DataDir[::-1].split('/', 2)   [-1]           [::-1]       +'/'
		#This is gross and I apologize to anyone reading this...

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(NormDir)
		AlfvenVelocity = Values[Variables.index('Alfven velocity')] 	#B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = Values[Variables.index('D gyro frequency')]		#
		IonDensity = Values[Variables.index('Bulk density')]			#
		B0 = Values[Variables.index('Mag.fld. at axis')]				#
		R0 = Values[Variables.index('raxis')]							#
		m_D = 3.34e-27													#
		eps = 0.5/R0													#

		#Denormalise temporal and spatial axes
		HarmonicData.kst = HarmonicData.kst								# [-]
		HarmonicData.time = HarmonicData.time * (1e3/IonGyroFreq)		# [ms]

		#Denormalise all variables
		# TO BE COMPLETED - MAY MOVE INTO EXTRACT FUNCTION???
	#endif


	return(HarmonicData)
#enddef

#=========================#
#=========================#

def MergePoloidal(HarmonicsData,VariableString,ntor):
#Reduces 3D HarmonicsData [mpol][ntor][lpsi][A/B] into 2D poloidal image [mpol,ltheta]
#Reduces 3D HarmonicsData by extracting only 1 ntor and averaging over mpol
#Also merges the real and imaginary components of HarmonicsData
#Inputs: HarmonicsData object [mpol][ntor][lpsi][A/B], 
#Inputs: VariableString matching HarmonicsData attribute, ntor Index (not mode number)
#Outputs: 2D PoloidalData2D array of shape [lpsi,ltheta]
#Example: PoloidalImage = MergePoloidal(HarmonicsData,'vrad',1)		

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

#=========================#
#=========================#

def Extract_RadialProfile(HarmonicsData,variable,ntorIdx,theta):
#Extracts radially resolved profiles for single variable at single poloidal angle
#Rounds poloidal angle down - i.e. anticlockwise - as defined from vertical zero.
#Inputs: HarmonicsData object [mpol][ntor][lpsi][A/B], 
#Inputs: VariableString matching HarmonicsData attribute, 
#Inputs: ntor Index (not mode number)
#Inputs: poloidal angle (theta) in degrees
#Outputs: 1D RadialProfile of shape [lpsi] (i.e. rho_pol)
#Example: Profile=Extract_RadialProfile(HarmonicsData,'brad',ntorIdx=2,theta=64)

	#Select variable and Merge 3D Data into 2D poloidal slice
	#PoloidalImage :: 2D array of Shape [lpsi][ltheta] ~ [R][theta]
	#Image[n][:] = full poloidal profile (clockwise from vertical) for R = Rho_pol[n]
	PoloidalImage = MergePoloidal(HarmonicsData,variable,ntorIdx)	

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

def ExtractMEGA_Markers(Dir,KStep,MarkerFileStep=1):
# Reads data from kinetic marker output files (gc_a_kstep000xxxx-00xxx.txt)
# Reads all variables and concatenates output data from all cores into single 2D array
#
# Inputs: Dir - Directory String to marker output file folder from root
# Inputs: KStep - KStep value (NOT Index) of output files to be read-in
# Inputs: MarkerFileStep - Optional speedup input, will read every 'n'th output file
#
# Outputs: KineticsData - 2D Array of shape [variable][marker(n)]
#		  Variables - R, Z, Lambda, E, p, Mu, pphi, fff, fnrml, psip, phi
#
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

	#Cycle through all marker files and read data
	for j in tqdm( range(0,len(MarkerFiles),MarkerFileStep) ):

		#Set current marker output file
		Filename = MarkerFiles[j]

		#Read marker data for current NCore output file
		#MarkerData :: 2D array of shape [variable,marker(n)]
		#Variables :: 0:R, 1:Z, 2:Lambda, 3:E, 4:p, 5:Mu, 6:pphi, 7:fff*fnrml, 8:psip, 9:phi
		#See write_ep() subroutine for full details
		MarkerData,Header = ReadFile_ASCII(Filename, 0, '2D', 'CSV')
#		print np.asarray(MarkerData).shape


		#Concatenate variables into KineticsData - Override KineticsData on first iteration
		#KineticsData :: 2D Array of shape [variable,marker(n)]
		if len(KineticsData) == 0: 
			KineticsData = MarkerData
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
#Example: VariableLegends = VariableLabelMaker(variables,Units)
def VariableLabelMaker(variables,Units=[]):

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
			Variable = '$r_{\psi}$'					# UNKNOWN VARIABLE
			VariableUnit = '[-]'
		elif variables[i] == 'gpsi_nrm':
			Variable = '$g_{psi}$ norm'				# UNKNOWN VARIABLE
			VariableUnit = '[-]'
		elif variables[i] == 'q_psi':
			Variable = 'Safety Factor $q_{\psi}$'
			VariableUnit = '[-]'

		#Explicit MHD fieldss
		elif variables[i] == 'brad':
			Variable = 'Radial B-field $B_{r}$'
			VariableUnit = '[T]'
		elif variables[i] == 'btheta':
			Variable = 'Poloidal B-field $B_{\\theta}$'
			VariableUnit = '[T]'
		elif variables[i] == 'bphi':
			Variable = 'Toroidal B-field $B_{\phi}$'
			VariableUnit = '[T]'
		elif variables[i] == 'erad':
			Variable = 'Radial E-field $E_{r}$'
			VariableUnit = '[V m$^{-1}$]'
		elif variables[i] == 'etheta':
			Variable = 'Poloidal E-field $E_{\\theta}$'
			VariableUnit = '[V m$^{-1}$]'
		elif variables[i] == 'ephi':
			Variable = 'Toroidal E-field $E_{\phi}$'
			VariableUnit = '[V m$^{-1}$]'
		#endif

		#Explicit MHD Velocities
		elif variables[i] == 'vrad':
			Variable = 'Radial Velocity $v_{r}$'				#Momentum?
			VariableUnit = '[m s$^{-1}$]'
		elif variables[i] == 'vtheta':
			Variable = 'Poloidal Velocity $v_{\\theta}$'		#Momentum?
			VariableUnit = '[m s$^{-1}$]'
		elif variables[i] == 'vphi':
			Variable = 'Toroidal Velocity $v_{\phi}$'			#Momentum?
			VariableUnit = '[m s$^{-1}$]'

		#Explicit MHD Densities and Pressures
		elif variables[i] == 'rho':
			Variable = 'Plasma Density $n_{e}$'
			VariableUnit = '[m$^{-3}$]'
		elif variables[i] == 'prs':
			Variable = 'Pressure $P_{rs}$'
			VariableUnit = '[Pa]'

		#Explicit Fast Particle Momentum
		elif variables[i] == 'mom_a':
			Variable = 'Fast Ion Momentum P$_{kin}$'
			VariableUnit = '[kg m s$^{-1}$]'
		elif variables[i] == 'dns_a':
			Variable = 'Fast Ion dns$_{a}$'				# UNKNOWN VARIABLE - KINETICS
			VariableUnit = '[-]'
		elif variables[i] == 'dns_a':
			Variable = 'Fast Ion dns$_{a}$'				# UNKNOWN VARIABLE - KINETICS
			VariableUnit = '[-]'
		elif variables[i] == 'ppara_a':
			Variable = 'Fast Ion ppara$_{a}$'			# UNKNOWN VARIABLE - KINETICS
			VariableUnit = '[-]'
		elif variables[i] == 'pperp_a':
			Variable = 'Fast Ion pperp$_{a}$'			# UNKNOWN VARIABLE - KINETICS
			VariableUnit = '[-]'
		elif variables[i] == 'qpara_a':
			Variable = 'Fast Ion qpara$_{a}$'			# UNKNOWN VARIABLE - KINETICS
			VariableUnit = '[-]'
		elif variables[i] == 'qperp_a':
			Variable = 'Fast Ion qperp$_{a}$'			# UNKNOWN VARIABLE - KINETICS
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

	#Calculates the Alfven eigenmode thresholds for all simulated poloidal mode numbers
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
print('                                                 v0.4.0')
print('-------------------------------------------------------')
print('')
print('The following diagnostics were requested:')
print('-----------------------------------------')
if True in [savefig_1Dtotalenergy,savefig_1Dspectralenergy]:
	print('# Energy Convergence Analysis')
if True in [savefig_1Dequilibrium]:
	print('# 1D Equilibrium Analysis')
if True in [savefig_2Dequilibrium,savefig_2Dequilmovie]:
	print('# 2D Equilibrium Analysis')
if True in [savefig_2Dharmonics]:
	print('# 2D Spectral Analysis')
if True in [savefig_2Dresponse]:
	print('# 2D Plasma Response Analysis')
if True in [savefig_1Dkinetics,savefig_2Dkinetics]:
	print('# Kinetic Distribution Analysis')
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
KineticsData = list()			#[2D Array] of shape Data[variables][markers(n)] concatinated for all nodes
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
				  #ENERGY & CONVERGENCE DIAGNOSTICS#
#====================================================================#

#====================================================================#
				 	    #TOTAL ENERGY ANALYSIS#
#====================================================================#

if savefig_1Dtotalenergy == True:

	#For each detected simulation folder
	for l in tqdm(range(0,len(Dir))):
		#Create global 1D diagnostics folder and extract current simulation name
		DirEnergy = CreateNewFolder(Dir[l],'1DEnergy_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_Phys outputs and header for plotting
		#Energy_Phys: [variable][timestep]
		Energy_Phys,Header_Phys = ExtractMEGA_Energy(Dir[l],'energy_phys')

		#Extract normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
#		print Variables[1],Values[1],Units[1]

		#==========##==========#

		#Create fig of desired size.
		fig,ax = figure(image_aspectratio,[3,1])

		#Define Title, Legend, Axis Labels etc...
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
	#endfor
#endif

#====================================================================#
				 	 #SPECTRAL ENERGY DIAGNOSTIC#
#====================================================================#


#Warning: Output is only defined for Order > 0!
#Example: DxDyArray = VectorDerivative(TimeArray,EnergyArray,Order=1 )[0]
def VectorDerivative(XArray,YArray,Order=1):

	#Compute i'th derivative of supplied arrays
	for i in range(0,Order):

		#Calculate derivative arrays - i.e. compute difference between successive indices
		DxArray = np.diff(XArray).tolist()
		DyArray = np.diff(YArray).tolist()
		print len(DxArray)

		#Calculate gradient array - i.e. derivative of Y to X
		DxDyArray = [DyArray[j]/DxArray[j] for j in range(0,len(DxArray))]

		#Sum derivatives up to each index in sequence to reconstruct XArray and set YArray to DxDy
		#Only required for calculation of higher order derivatives (Order > 1) - ARRAYS NOT RETURNED
		XArray = [sum(DxArray[0:j]).tolist() for j in range(0,len(DxArray))]
		YArray = np.log(DxDyArray)		# HACKY, ONLY TRUE FOR EXPONENTIAL GROWTH
	#endfor

	return(DxDyArray,DxArray,DyArray)
#enddef





if savefig_1Dspectralenergy == True:

	#For each detected simulation folder
	for l in tqdm(range(0,len(Dir))):

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
#		print Dir[l].split('/')[-2]
		seq = setting_seq[1]				#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
		ntor = 2 #setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!
		KstepIdx = setting_kstep[1]			#requested kstep index 				!!! NEEDS A FUNCTION !!!


		#Create global 1D diagnostics folder and extract current simulation name
		DirEnergy = CreateNewFolder(Dir[l],'1DEnergy_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header for plotting
		#energy_n: [variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		KStepArray = Energy_n[0]					#KStep Array			[-]
		TimeArray = Energy_n[1]						#Time Array				[ms]
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/(seq+1)			#KStep indices per seq 	[-]

		#Extract normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
#		print Variables[1],Values[1],Units[1]

		#Compute 1st derivative of energy for each harmonic
		DeltaEnergy_n = list()
		DeltaEnergy_n.append(KStepArray[0:-1])			#Add kstep array
		DeltaEnergy_n.append(TimeArray[0:-1])			#Add time array
		for i in range(2,len(Energy_n)):
			EnergyArray = np.log(Energy_n[i])			#Add i'th harmonic array
			DeltaEnergy_n.append( VectorDerivative(TimeArray,EnergyArray,1 )[0] )
			#endfor
		#endfor

		#Compute 2nd derivative of energy for each harmonic
		Delta2Energy_n = list()
		Delta2Energy_n.append(KStepArray[0:-2])			#Add kstep array
		Delta2Energy_n.append(TimeArray[0:-2])			#Add time array
		for i in range(2,len(Energy_n)):
			EnergyArray = np.log(Energy_n[i])			#Add i'th harmonic array
			Delta2Energy_n.append( VectorDerivative(TimeArray,EnergyArray,2 )[0] )
			#endfor
		#endfor



		#==========##==========#

		#NOTES 
#		- 	VECTOR DERIVATIVE NEEDS COMMENTED AND CLEANED UP (QUITE CONFUSING AT THE MOMENT)
#		- 	THE 1ST AND 2ND DERIVATIVES NEED CLEANED UP, IDEALLY JUST WITH ONE OR TWO LINES
#		-	NEED TO CHECK PHYSICAL 'MEANING' OF THE DERIVATIVE, I THINK FLAT DERIVATIVES MAKE SENSE HERE?
#		-	NEED TO CLEAN UP THE LinearRegion AND SMOOTHING STUFF, POTENTIALLY MAKE A FUNCTION?
#		- 	THE BELOW METHOD IS BAD AS IT ONLY USES THE FIRST AND LAST ELEMENTS OF LINEAR REGION
#			BETTER METHOD IS TO AVERAGE ALL THE 1st DERIVATIVES USING LinearRegion AS A MASK
#			NEED TO CHECK THE MATHS ON THAT THOUGH, JUST TO ENSURE IT ALL MAKES SENSE...	

		#Smooth 2nd derivative trends
		#Smooth kinetic data prior to analysis if requested (Savitzk-Golay filter)
		if KineticFiltering == True:
			WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
			Delta2Energy_n[ntor+2] = (savgol_filter(Delta2Energy_n[ntor+2], WindowSize, PolyOrder)).tolist()
		#endif

#		if DebugMode == True:
#			WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
#			Delta2Energy_n_DEBUG = (savgol_filter(Delta2Energy_n[2], WindowSize, PolyOrder)).tolist()
#
#			plt.plot(Delta2Energy_n[1],Delta2Energy_n[2], 'k-', lw=2)			#Unsmoothed
#			plt.plot(Delta2Energy_n[1],Delta2Energy_n_DEBUG, 'r-', lw=2)		#Smoothed
#			plt.legend(['Unsmoothed','Smoothed'])
#			plt.show()
		#endif

		#Determine temporal extent of linear growth region, i.e. where 2nd derivative is close to zero
		Threshold = 100
		LinearRegion = list()
		for i in range(0,len(Delta2Energy_n[0])):
			if abs(Delta2Energy_n[ntor+2][i]) < Threshold: 	LinearRegion.append(1)
			else: 											LinearRegion.append(0)
		#endfor

		#Smooth Threshold (Savitzk-Golay filter)
		WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
		LinearRegion = (savgol_filter(LinearRegion, WindowSize, PolyOrder)).tolist()
		#endif

		#Cut Linear Region at 0.5
		for i in range(0,len(LinearRegion)):
			if LinearRegion[i] > 0.5: LinearRegion[i] = 1.0
			else: LinearRegion[i] = 0.0
		#endfor

		#Compute linear phase growth rate gamma, where Eend = Estart*exp{gamma*dt}
		#gamma = ln(Eend/Estart)/dt
		StartIdx = LinearRegion.index(1)
		EndIdx = len(LinearRegion) - LinearRegion[::-1].index(1)

		tstart = TimeArray[StartIdx]
		tend = TimeArray[EndIdx]
		dt = tend-tstart

		Estart = Energy_n[ntor+2][StartIdx]					#Can't be zero
		Eend = Energy_n[ntor+2][EndIdx]
		Eratio = Eend/Estart

		gamma = np.log(Eend/Estart)/dt
		print('')
		print( round(tstart,3), round(tend,3) )
		print( round(Estart,3), round(Eend,3) )
		print( round(gamma,2), '[s-1]')
		print('')

#		exit()

		#==========##==========#



		#Create fig of desired size.
		fig,ax = figure(image_aspectratio,[2,1])

		#Define Title, Legend, Axis Labels etc...
		Title = 'Spectrally Resolved Energy Evolution for \n '+DirString
		Xlabel,Ylabel = '', 'Energy $\epsilon_{n}$ (Log$_{10}$) [-]'
		Legend = list()

		#Plot total energy for each harmonic component
		for i in range(2,len(Energy_n)):
			ax[0].plot(TimeArray,np.log10(Energy_n[i]), lw=2)
			Legend.append( 'n = '+str(i-2) )
		#endfor
		ImageOptions(fig,ax[0],Xlabel,Ylabel,Title,Legend)


		#Define Title, Legend, Axis Labels etc...
		Title = 'Spectrally Resolved Energy Evolution for \n '+DirString
		Xlabel,Ylabel = 'Time [ms]', '$\Delta$ Energy $\\frac{d \epsilon_{n}}{d t}$ (Log$_{10}$) [-]'
		Legend = list()

		#Plot 1st derivative of energy for each harmonic component
		for i in range(2,len(Energy_n)):
			ax[1].plot(TimeArray[0:-1],np.log10(DeltaEnergy_n[i]), lw=2)
			Legend.append( 'n = '+str(i-2) )
		#endfor
		ax[1].plot(Delta2Energy_n[1],LinearRegion)				#if DebugMode == True:
		ImageOptions(fig,ax[1],Xlabel,Ylabel,'',Legend)

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


















#Input: Dir is the simulation root directory folder
#Input: DataFile is a string determining which output file is used for the time axes
#		DataFile options: 'energy_n', 'harmonics', 'moments'
#Output:
#Example: KStepArray, TimeArray, ntorArray = ExtractMEGA_Axes(Dir, DataFile='harmonics')
def ExtractMEGA_Axes(Dir, DataFile='energy_n'):

	#Extract kstep, time and toroidal harmonic data from energy_n.txt
	#energy_n data structure: [variable][timestep]
	Energy_n,Header_n = ExtractMEGA_Energy(Dir[l], 'energy_n')
	#Determine poloidal and toroidal harmonic ranges
	ntor_tot = ((len(Energy_n)-3)*2)+1			#Total number of positive and negative modes (Including n=0)
	ntor_pos = int(float(ntor_tot-1)/2.0)		#Number of positive modes (Ignoring n=0)
	ntor0 = int(ceil(ntor_tot/2))				#ntor = 0, baseline equilibrium data

	#READ ONLY TEMPORAL AND SPATIAL AXES from SEQ.Harmonics, it contains:
	#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
	#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
	HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/', Variable='NaN', ntor=ntor_tot, Dimension='1D')
	rho_pol = HarmonicsData.rho_pol				#Normalised radius		[-]

	#Determine KStep range, Time range and related intervals
	#Use wchck time intervals
	if DataFile == 'energy_n':
		KStepArray = Energy_n[0]					#KStep Array			[-]
		TimeArray = Energy_n[1]						#Time Array				[ms]

	#Use wharm time intervals
	elif DataFile == 'harmonics':
		KStepArray = HarmonicsData.kst				#KStep Array			[-]
		TimeArray = HarmonicsData.time				#Time Array				[ms]

	#Use wsnapshot time intervals
	elif DataFile == 'moments':
		a = 1										#TO BE COMPLETED
	#endif

	return(KStepArray, TimeArray, [ntor0,ntor_pos,ntor_tot])
#enddef



#Takes toroidal mode number and returns index relating to said number.
#Example: ntorIdx = ntor_Index(ntor,ntorArray)
def ntor_Index(ntor,ntorArray):

	#Unpack input ntorArray
	ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
	ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
	ntor0 = ntorArray[0]						#ntor = 0, baseline equilibrium data

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

	#Print debug output to terminal if requested
	if DebugMode == True: 
		print('')
		print('Toroidal Mode Numbers:',ntor_indices)
		print('')

	#ntor range set by ntor_indices[ntor, ntor_index], for pos, neg & n=0 modes
#	ntor = 0													#override requested ntor mode number
	ntorIdx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number

	return(ntorIdx)
#enddef




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
		seq = setting_seq[1]				#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!
		KStepIdx = setting_kstep[1]			#requested kstep index 				!!! NEEDS A FUNCTION !!!


		#Create global 1D diagnostics folder and extract current simulation name
		DirEquil1D = CreateNewFolder(Dir[l],'1DEquil_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Kstep [-] & Time [ms] arrays from SEQ.harmonics & toroidal harmonics from energy_n.txt
		KStepArray, TimeArray, ntorArray = ExtractMEGA_Axes(Dir, DataFile='harmonics')
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/(seq+1)			#KStep indices per seq 	[-]
		ntor_tot = ntorArray[2]						#Total number of positive & negative modes (Inc n=0)
		ntor_pos = ntorArray[1]						#Number of positive modes (Ignoring n=0)
		ntor0 = ntorArray[0]						#ntor = 0, baseline equilibrium data

		#Extract toroidal mode number array index (ntorIdx) from requested mode number (ntor)
		ntorIdx = ntor_Index(ntor,ntorArray)

		#Set TimeIndex and employ to extract KStep and Time
		IdxOffset = seq*KStepMod					#[-]
		KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
		Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

		#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KStepIdx,seq,'3D')
		rho_pol = HarmonicsData.rho_pol				#Normalised radius		[-]

		#Extract data resolution and poloidal axes from repository .dat files
		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)
		DataShape = ExtractMEGA_DataShape(HarmonicsData)#; print DataShape
		mpol_res = DataShape[0]; ntor_res = DataShape[1]
		lpsi_res = DataShape[2]; ltheta_res = DataShape[3]


		#For each requested variable
		for i in tqdm(range(0,len(variables))):

			#Create Variablelabel with units
			VariableLabel = VariableLabelMaker(variables[i])

			#==========#	#==========#

			#Create figure and define Title, Legend, Axis Labels etc...
#			fig,ax = figure(image_aspectratio,1)
#			Title = VariableLabel+', n='+str(ntor)+', m='+str(-mpol_res+1)+','+str(mpol_res-1)+', t='+str(round(Time,3))+' ms \n Simulation: '+DirString
#			Xlabel,Ylabel = 'Radius $R$ [m]', VariableLabel
#			Legend = list()

			#Plot 1D poloidally resolved profiles for current simulation folder
			#Poloidal profiles employ fixed radial (rho_phi) and toroidal (phi) angles
#			for j in range(0,len(poloidalprofiles)):

				#Create new folder to store poloidal profiles
#				if j == 0: DirEquilPoloidal = CreateNewFolder(DirEquil1D,'Poloidal/')

				#### CAN BE A FUNCTION (Extract_PoloidalProfile())
				#Poloidally resolved profiles
	#			Full = len(Image[0])					#Outboard Poloidal Cells
	#			Half = Full/2							#Inboard Poloidal Cells
	#			plt.plot(Crdr[50][0:Half],Image[50][0:Half], 'k-')		#Outboard poloidal profile at R = Rho_pol[50]
	#			plt.plot(Crdr[50][Half:Full],Image[50][Half:Full], 'r-')#Inboard poloidal profile at R = Rho_pol[50]
	##			plt.plot(Crdr[50],Image[50])							#Full Poloidal profile at R = Rho_pol[50]
	#			plt.xlabel('Radius $R$ [m]')
	#			plt.ylabel('Height $Z$ [m]')
	#			plt.show()
				#### CAN BE A FUNCTION (Extract_PoloidalProfile())
			#endfor

			#Beautify 1D equilibrium profiles figure
#			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

			#Save poloidal equilibrium profiles for current simulation folder
#			SaveString = variables[i]+'_Poloidal_n'+str(ntor)+'_t='+str(round(Time,3))+ext
#			plt.savefig(DirEquilPoloidal+SaveString)
#			plt.show()
#			plt.close('all')

			#==========#	#==========#

			#Create figure and define Title, Legend, Axis Labels etc...
			fig,ax = figure(image_aspectratio,1)
			ntorString = ', n='+str(ntor); mpolString=', m='+str(-mpol_res+1)+','+str(mpol_res-1)
			TimeString = ', t='+str(round(Time,3))+' ms'
			Title = VariableLabel+ntorString+mpolString+TimeString+' \n Simulation: '+DirString
			Xlabel,Ylabel = 'Radius $R$ [m]', VariableLabel
			Legend = list()

			#Plot 1D radially resolved profiles for current simulation folder
			#Radial profiles employ fixed poloidal (theta) and toroidal (phi) angles
			for j in range(0,len(radialprofiles)):

				#Create new folder to store radial profiles
				if j == 0: DirEquilRadial = CreateNewFolder(DirEquil1D,'Radial/')

				#Define poloidal angle theta and append to legend list
				theta = radialprofiles[j]
				Legend.append('$\\theta$ = '+str(theta)+'$^{\circ}$')

				#Extract radially resolved profile and plot 
				RadialProfile = Extract_RadialProfile(HarmonicsData,variables[i],ntorIdx,theta)
				ax.plot(rho_pol,RadialProfile, lw=2)

				#Save ASCII data to sub-folder
#				Write_data_to_file(RadialProfile)			### TO BE ADDED ###
			#endfor

			#Beautify 1D equilibrium profiles figure
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

			#Save radial equilibrium profiles for current simulation folder
			SaveString = variables[i]+'_Radial_n'+str(ntor)+'_t='+str(round(Time,3))+ext
			plt.savefig(DirEquilRadial+SaveString)
#			plt.show()
			plt.close('all')
		#endfor	- variable loop
	#endfor	- directory loop
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

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		seq = setting_seq[1]				#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!
		KstepIdx = setting_kstep[1]			#requested kstep index 				!!! NEEDS A FUNCTION !!!

		#Create global 2D diagnostics folder and extract current simulation name
		DirEquil = CreateNewFolder(Dir[l],'2DPoloidal_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#### - CAN BE A FUNCTION
		#Extract kstep, time and toroidal harmonic data from energy_n.txt
		#energy_n data structure: [variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		#Determine KStep range, Time range and related intervals
		KStepArray = Energy_n[0]					#KStep Array			[-]
		TimeArray = Energy_n[1]						#Time Array				[ms]
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/(seq+1)			#KStep indices per seq 	[-]
		#Determine poloidal and toroidal harmonic ranges
		ntor_tot = ((len(Energy_n)-3)*2)+1			#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)		#Number of positive modes (Ignoring n=0)
		ntor0 = int(ceil(ntor_tot/2))				#ntor = 0, baseline equilibrium data
		#### - CAN BE A FUNCTION

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
		ntorIdx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number
		#### - CAN BE A FUNCTION

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		IonGyroFreq = Values[Variables.index('D gyro frequency')]

		#Set TimeIndex and employ to extract KStep and Time
		KStepIdx = setting_kstep[1]; seqIdx = seq	#Duplicated from Development Settings Above
		IdxOffset = seqIdx*KStepMod					#[-]
		KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
		Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

		#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KStepIdx,seqIdx,'3D')

		#Extract data resolution and poloidal axes from repository .dat files
		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)
		DataShape = ExtractMEGA_DataShape(HarmonicsData)#; print DataShape
		mpol_res = DataShape[0]; ntor_res = DataShape[1]
		lpsi_res = DataShape[2]; ltheta_res = DataShape[3]

		#For each requested variable
		for i in tqdm(range(0,len(variables))):

			#Create Variablelabel with units
			VariableLabel = VariableLabelMaker(variables[i])

			#Select variable and Merge 3D Data into 2D poloidal slice
			#Image Shape: [lpsi][ltheta] ~~ [R][theta], like an onion (or ogre).
			#i.e. Image[n][:] plots a full poloidal profile (clockwise from vertical) for R = Rho_pol[n]
			Image = MergePoloidal(HarmonicsData,variables[i],ntorIdx)

			#Create figure and define Title, Legend, Axis Labels etc...
			fig,ax = figure(image_aspectratio,1)
			Title = VariableLabel+', n='+str(ntor)+', t='+str(round(Time,3))+' ms \n Simulation: '+DirString
			Xlabel,Ylabel = 'Major Radius $R$ [m]', 'Height $Z$ [m]'
			Legend = list()

			#Plot 2D poloidally resolved figure and beautify
			im = ax.contourf(Crdr, Crdz, Image, 100); plt.axis('scaled')
			im2 = ax.contour(Crdr, Crdz, Image, 20); plt.axis('scaled')
			cbar = Colourbar(ax,im,VariableLabel,5)
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

			#Save 2D poloidally resolved figure for current simulation
			SaveString = variables[i]+'_n'+str(ntor)+'_t='+str(round(Time,3))+ext
			plt.savefig(DirEquil+SaveString)
#			plt.show()
			plt.close('all')
		#endfor
	#endfor
#endif


#====================================================================#
				 	      #2D POLOIDAL MOVIES#
#====================================================================#

if savefig_2Dequilmovie == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		seq = setting_seq[1]				#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!

		#Create global 2D diagnostics folder and extract current simulation name
		DirEquil = CreateNewFolder(Dir[l],'2DPoloidal_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#### - CAN BE A FUNCTION
		#Extract kstep, time and toroidal harmonic data from energy_n.txt
		#energy_n data structure: [variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		#Determine KStep range, Time range and related intervals
		KStepArray = Energy_n[0]					#KStep Array			[-]
		TimeArray = Energy_n[1]						#Time Array				[ms]
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/(seq+1)			#KStep indices per seq 	[-]
		#Determine poloidal and toroidal harmonic ranges
		ntor_tot = ((len(Energy_n)-3)*2)+1			#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)		#Number of positive modes (Ignoring n=0)
		ntor0 = int(ceil(ntor_tot/2))				#ntor = 0, baseline equilibrium data
		#### - CAN BE A FUNCTION

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
		ntorIdx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number
		#### - CAN BE A FUNCTION

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		IonGyroFreq = Values[Variables.index('D gyro frequency')]


		#Apply user Kstep range if requested - else default to max range
		if len(setting_kstep) == 3: 
			KStepRange = [setting_kstep[0],setting_kstep[1]]
			KStepStep = setting_kstep[2]
		elif len(setting_kstep) == 2:
			KStepRange = setting_kstep
			KStepStep = 1
		else: KStepRange = [0,len(KStepArray)]
		#endif

		#Extract and plot data for each timestep
		for i in tqdm( range(KStepRange[0],KStepRange[1],KStepStep) ):

			#Set TimeIndex and employ to extract KStep and Time
			KStepIdx = i; seqIdx = seq					#Add these to seq loop.
			IdxOffset = seqIdx*KStepMod					#[-]
			KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
			Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

			#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
			#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
			#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
			#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
			HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KStepIdx,seqIdx,'3D')

			#Extract data resolution and poloidal axes from repository .dat files
			#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
			DataShape = ExtractMEGA_DataShape(HarmonicsData)
			Crdr,Crdz = ExtractMEGA_PoloidalGrid(DirRepository,HarmonicsData)

			#For each requested variable at the current Kstep
			for j in range(0,len(variables)):

				#Create global 2D diagnostics folder and extract current simulation name
				DirMovie = CreateNewFolder(DirEquil,variables[j]+'_n'+str(ntor)+'_Movie/')

				#Select variable and Merge 3D Data into 2D poloidal slice
				Image = MergePoloidal(HarmonicsData,variables[j],ntorIdx)


				#Create figure and define Title, Legend, Axis Labels etc...
				fig,ax = figure(image_aspectratio,1)

				#Extract Variablelabel and define figure texts
				VariableLabel = VariableLabelMaker(variables[j])
				Title = VariableLabel+', n='+str(ntor)+', t='+str(round(Time,3))+' ms \n Simulation: '+DirString
				Xlabel,Ylabel = 'Radius $R$ [m]', 'Height $Z$ [m]'
				Legend = list()

				#Plot 2D poloidally resolved figure and beautify
				im = ax.contourf(Crdr, Crdz, Image, 100); plt.axis('scaled')
				im2 = ax.contour(Crdr, Crdz, Image, 20); plt.axis('scaled')
				cbar = Colourbar(ax,im,VariableLabel,5)
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

				#Save 2D poloidally resolved figure for current simulation
				SaveString = variables[j]+'_n'+str(ntor)+'_kstep'+str('%07.f'%KStep)+ext
				plt.savefig(DirMovie+SaveString)
#				plt.show()
				plt.close('all')
			#endfor

			#!!! AUTO CREATE MOVIES FOR EACH VARIABLE HERE !!!
			#!!! AUTO CREATE MOVIES FOR EACH VARIABLE HERE !!!

		#endfor
	#endfor
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
					 #PLASMA RESPONSE DIAGNOSTICS#
#====================================================================#

#====================================================================#
					 #POLOIDAL RESPONSE ANALYSIS#
#====================================================================#

if savefig_2Dresponse == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		seq = setting_seq[1]				#requested SEQ file index (001 = 0)	!!! NEEDS A FUNCTION !!!
		ntor = setting_ntor[1]				#requested ntor mode number			!!! NEEDS A FUNCTION !!!
		variable = 'brad'					#requested response variable 		!!! Need to impliment vrad etc...

		#Initiate any required lists
		DataAmpPROES = list()
		XaxisPROES = list()

		#Create global 2D diagnostics folder and extract current simulation name
		DirResponse = CreateNewFolder(Dir[l],'2DResponse_Plots/')			#Response Folder
		DirResponse_ntor = CreateNewFolder(DirResponse,'ntor='+str(ntor))	#Toroidal Mode Folder	
		DirString = Dir[l].split('/')[-2]									#Full Simulation Name
		SubString = DirString.split('_')[-1]								#Simulation Nickname

		#### - CAN BE A FUNCTION
		#Extract kstep, time and toroidal harmonic data from energy_n.txt
		#energy_n data structure: [variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		#Determine KStep range, Time range and related intervals
		KStepArray = Energy_n[0]					#KStep Array			[-]
		TimeArray = Energy_n[1]						#Time Array				[ms]
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
		KStepMod = len(KStepArray)/(seq+1)			#KStep indices per seq 	[-]
		#Determine poloidal and toroidal harmonic ranges
		ntor_tot = ((len(Energy_n)-3)*2)+1			#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)		#Number of positive modes (Ignoring n=0)
		ntor0 = int(ceil(ntor_tot/2))				#ntor = 0, baseline equilibrium data
		#### - CAN BE A FUNCTION

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
		ntorIdx = [item[0] for item in ntor_indices].index(ntor)	#index referring to ntor mode number
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

		#Extract Variablelabel for chosen variable
		VariableLabel = VariableLabelMaker(variable,Units='Perturbation [-]')

		#Apply user Kstep range if requested - else default to max range
		if len(setting_kstep) == 2: KStepRange = setting_kstep
		else: KStepRange = [0,len(KStepArray)]
		#endif

		#Apply user seq range if requested - else use supplied single seq
		if len(setting_seq) == 2: seqRange = [0,setting_seq[1]+1]
		else: seqRange = [setting_seq,setting_seq+1]
		#endif

		for j in range(seqRange[0],seqRange[1]):
			#Set SeqIndex for current simulation folder
			seqIdx = j

			#Extract and plot data for each timestep
			for i in tqdm(range(KStepRange[0],KStepRange[1])):

				#Set TimeIndex and employ to extract KStep and Time
				KStepIdx = i
				IdxOffset = seqIdx*KStepMod					#[-]
				KStep = KStepArray[KStepIdx+IdxOffset]		#[-]
				Time = TimeArray[KStepIdx+IdxOffset]		#[ms]

				#Extract ALL VARIABLES FOR SINGLE KSTEP from Harmonics, it contains:
				#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
				#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
				#HarmonicsData.Variables[i]: [3D Array] of shape [mpol][ntor][lpsi][A/B] for a single kstep
				HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/','All',ntor_tot,KStepIdx,seqIdx,'3D')

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
				#DataAmp is of shape: DataAmp[lpsi,mpol]
				DataAmp = np.concatenate((DataAmpNeg,DataAmpPos[1:,:]),axis=0)

				#Create Image array and Axes, rotate such that mpol spectrum is on X-axis.
				#Image is of shape: [mpol][lpsi] 
				Image = DataAmp.transpose()*B0*1e4								#[G]? - must use Brad
				Xaxis =	[x-int(mpol_res-1) for x in range(0,2*mpol_res-1,1)]	#Poloidal Mode Numbers
				Yaxis = HarmonicsData.rho_pol									#Radial Location

				#Create figure and define Title, Legend, Axis Labels etc...
				fig,ax = figure(image_aspectratio,1)
				Title = 'Plasma Response: n='+str(ntor)+', m='+str(-mpol_res+1)+','+str(mpol_res-1)+', t='+str(round(Time,3))+' [ms] \n Simulation: '+DirString
				Xlabel,Ylabel = 'Poloidal Harmonic $m_{pol}$ [-]', 'Radial Magnetic Coordinate $\\rho_{pol}$ [-]'
				Legend = list()

				#Plot poloidal response spectrum figure (R,mpol)
				im = ax.contourf(Xaxis, Yaxis, Image, 50)
				res = ax.plot(-ntor*HarmonicsData.q_psi,HarmonicsData.rho_pol, 'w--', lw=2)
				cbar = Colourbar(ax,im,VariableLabel,5)
				#####
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)
				ax.set_xlim(-16,16)

				#Save poloidal response figure for current seq and Kstep
				SaveString = 'PlasmaResponse_'+variable+'_n'+str(ntor)+'_kstep'+str('%07.f'%KStep)+ext
				plt.savefig(DirResponse_ntor+SaveString)
#				plt.show()
				plt.close('all')

				#==========#

				#Collapse figure poloidally - integrate through all poloidal modes
				DataAmp1D = list()
				DataAmp = DataAmp.transpose()
				for k in range(0,len(DataAmp)):
					DataAmp1D.append(sum(DataAmp[k][:]))
				#endfor
				
				#Append 1D array to 2D PROES image array
				#DataAmpPROES: 2D array of shape [kstep][lpsi]
				DataAmpPROES.append(DataAmp1D)
				XaxisPROES.append(Time)
			#endfor
		#endfor

		#If more than one KStep was processed, create a temporal plasma response image
		if len(DataAmpPROES) > 1:
			#Create figure and define Title, Legend, Axis Labels etc...
			fig,ax = figure(image_aspectratio,1)
			Title = 'Plasma Response: n='+str(ntor)+', m='+str(-mpol_res+1)+','+str(mpol_res-1)+', t='+str(round(Time,3))+' [ms] \n Simulation: '+DirString
			Xlabel,Ylabel = 'Time $t$ [ms]', 'Radial Magnetic Coordinate $\\rho_{pol}$ [-]'
			Legend = list()

			#Plot temporally resolved, poloidally collapsed, response figure (R,time)
			DataAmpPROES = np.asarray(DataAmpPROES).transpose()			#Transpose to align time on X-axis
			im = plt.contourf(XaxisPROES, Yaxis, DataAmpPROES, 50)
			cbar = Colourbar(ax,im,VariableLabel,5)
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

			#Save temporal response figure for current simulation directory
			SaveString = 'PlasmaResponse_'+variable+'_n'+str(ntor)+'_t='+str(round(Time,3))+ext
			plt.savefig(DirResponse+SaveString)
	#		plt.show()
			plt.close('all')
		#endif
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
				   #SPECTRAL & HARMONIC DIAGNOSTICS#
#====================================================================#

#====================================================================#
				   #2D HARMONIC & FOURIER ANALYSIS#
#====================================================================#

if savefig_2Dharmonics == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		#DEVELOPMENT SETTINGS - settings_inputs to be moved to switchboard
		print Dir[l].split('/')[-2]
		variable = 'bphi'				#requested response variable 		!!! Need to impliment btheta, vrad etc...

		#Create global 2D diagnostics folder and extract current simulation name
		DirHarmonics = CreateNewFolder(Dir[l],'2DHarmonic_Plots/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#### - CAN BE A FUNCTION
		#Extract kstep, time and toroidal harmonic data from energy_n.txt
		#energy_n data structure: [variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		#Determine KStep range, Time range and related intervals
		KStepArray = Energy_n[0]					#KStep Array			[-]
		TimeArray = Energy_n[1]						#Time Array				[ms]
		DeltaKstep = KStepArray[1]-KStepArray[0]	#KStep Interval 		[-]
		DeltaTime = TimeArray[1]-TimeArray[0]		#Time Interval 			[ms]
#		KStepMod = len(KStepArray)/(seq+1)			#KStep indices per seq 	[-]
		#Determine poloidal and toroidal harmonic ranges
		ntor_tot = ((len(Energy_n)-3)*2)+1			#Number of positive and negative modes (Including n=0)
		ntor_pos = int(float(ntor_tot-1)/2.0)		#Number of positive modes (Ignoring n=0)
		ntor0 = int(ceil(ntor_tot/2))				#ntor = 0, baseline equilibrium data
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
		HarmonicsData = ExtractMEGA_Harmonics(Dir[l]+'data/',variable,ntor_tot,Dimension='4D')
		Data = HarmonicsData.data

		#DataShape contains data resolution of form: [mpol,ntor,lpsi,ltheta]
		DataShape = ExtractMEGA_DataShape(HarmonicsData)
		mpol, ntor = DataShape[0], DataShape[1]
		lpsi, ltheta = DataShape[2], DataShape[3]
		kmax, dt = DataShape[4], (TimeArray[1]-TimeArray[0])
		ntor2 = int(0.5*(ntor-1))					#Positive ntor (excluding n=0)

		#Extract Variablelabels
		VariableLabel = VariableLabelMaker(variable)

		#BELOW TO STILL BE TRANSLATED
		#ALSO NEED TO ADD COLOURBAR TO THE FOURIER PLOTS!!!
		print kmax, mpol, ntor, lpsi, ntor2


		#Sum Re component of toroidal (n) and poloidal (m) modes for all ksteps
		vcos = list()
		for n in range(0,ntor2):
			vcos.append( np.zeros([]) )				
			for m in range(0,mpol):						#Data structure: [kstep][mpol][ntor][lpsi][Re/Im] 
				vcos[n] = vcos[n] + Data[:,m,n,:,0]		#vcos structure: [ntor][kstep][lpsi]
			#endfor
			vcos[n][np.isnan(vcos[n])] = 0			#Remove any NaNs
		#endfor

#		print len(vcos), len(vcos[0]), len(vcos[0][0])
#		print vcos[0]
#		plt.plot(vcos[0])
#		plt.show()
#		exit()

		vcos_fft,vcos_len = list(),list()
		#Extract fourier components from vcos
		for n in range(0,ntor2):
			vcos_fft.append( np.fft.fft(vcos[n],axis=0) )	#Take fourier components of variable
		  	vcos_fft[n][0,:] = vcos_fft[n][0,:]*0.0			#Discard imaginary components 					???
			vcos_len.append( int(len(vcos_fft[n])/2)-1 )	#Determine lowpass filter frequency threshold 	???
			vcos_fft[n] = vcos_fft[n][0:vcos_len[n],:]		#Discard upper frequencies (lowpass filter)		???
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
			#endif

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
				ImageOptions(fig,subfig,'Radius $R$','Frequency [kHz]','','')
			elif i == 0 and ntor2 == 1:					#If first panel with no panels to right
				cbar = Colourbar(subfig,im,VariableLabel,5)
				ImageOptions(fig,subfig,'Radius $R$','Frequency [kHz]','','')
			elif i > 0 and i < ntor2-1:   				#If middle panel with more panels to right
				ImageOptions(fig,subfig,'Radius $R$','','','')
 				im.axes.get_yaxis().set_visible(False)
			elif i == ntor2-1 and ntor2 > 1:			#If last panel with more panels to left
				cbar = Colourbar(subfig,im,VariableLabel,5)
				ImageOptions(fig,subfig,'Radius $R$','','','')
 				im.axes.get_yaxis().set_visible(False)
			#endif
			subfig.set_ylim([0,200])

			#Compute and plot TAE thresholds
			UpperThresholds,LowerThresholds = TAEThresholds(HarmonicsData,Harmonic,eps,AlfvenVelocity,subfig)
		#endfor

		#Save 2D harmonics figure for current simulation
		plt.savefig(DirHarmonics+'FourierAnalysis_'+SubString+ext)
#		plt.show()
		plt.close('all')
	#endfor
#endif

#==========##==========##==========#
#==========##==========##==========#

if any([savefig_2Dharmonics]) == True:
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
		fig,ax = figure(image_aspectratio,1)
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
			fig,ax = figure(image_aspectratio,1)
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
	nBins = 100					#Kinetics Histogram Bins		- Move to Low-Level Inputs
	KStepMin = 200000				#KStepMin						- Automate readin - Use Switchboard?
	KStepMax = 400000			#KStepMax						- Automate readin - Use Switchboard?
	KWep = 10000				#Write_ep save interval (kwep)	- Automate readin - Use Switchboard?
	KMarker = 1					#Marker file readin interval	- Move to Low-Level Inputs


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
			#Variables :: R, Z, Lambda, E, p, Mu, pphi, fff*fnrml, psip, phi
			KineticsData,Header_Kin = ExtractMEGA_Markers(Dir[l],KStep,KMarker)

			#Select variables to be plotted (X,Y axis)
			XData = KineticsData[0]				# [3] 'E_gc'	[Typically Spatial Variable]
			YData = KineticsData[6]				# [6] 'pphi_gc'	[Typically Physical Variable]

			#Extract min/max values and create histogram ranges
			Xmin,Xmax = min(XData), max(XData)
			Ymin,Ymax = min(YData), max(YData)
			XRange = np.linspace(Xmin,Xmax,nBins)
			YRange = np.linspace(Ymin,Ymax,nBins)

			#Select 2D variables to be plotted (X,Y) and histogram over supplied ranges
			HistData2D,XAxis2D,YAxis2D = np.histogram2d(XData,YData, bins=(XRange,YRange))
			extent = [min(XAxis2D),max(XAxis2D), min(YAxis2D),max(YAxis2D)]				#USE FIXED FULL RANGE

			#Select 1D variable to be plotted (X axis) and histogram into nBins
			HistData1D,XAxis1D = np.histogram(XData, bins=nBins)

			#Normalise 2D distribution function
			HistSum2D = sum(sum(HistData2D)); NormFactor2D = HistSum2D
			for x in range(0,len(HistData2D)):
				for y in range(0,len(HistData2D[x])):
					HistData2D[x,y] = float(HistData2D[x,y])/float(NormFactor2D)
				#endfor
			#endfor
			if DebugMode == True: print( "2D IEDF Integral: ",str(sum(HistData2D)) )

			#Normalise 1D distribution function
			HistSum1D = sum(HistData1D); NormFactor1D = HistSum1D
			HistData1D = [float(HistData1D[x])/float(NormFactor1D) for x in range(0,len(HistData1D))]
			if DebugMode == True: print( "1D IEDF Integral: ",str(sum(HistData1D)) )

			#Initiate figure and set axes
			fig,ax = figure(image_aspectratio,[2,1],shareX=True)
#			Title = VariableLabels[j]+', ntor='+str(ntor)+', t='+str(Time)+' \n Simulation: '+DirString
			Title = 'Kinetic Markers, Kstep='+str(i)+' \n Simulation: '+DirString
#			Xlabel,Ylabel = 'Energy $\epsilon_{i}$ [keV]','Angular Momentum $p_{\phi}$ \n [kg m${^2}$ s$^{-1}$]'
			Xlabel,Ylabel = 'Radius R','Angular Momentum $p_{\phi}$ \n [kg m${^2}$ s$^{-1}$]'
#			Xlabel,Ylabel = 'Radius R','Energy $\epsilon_{i}$ [keV]'
#			Xlabel,Ylabel = 'Radius R','Height Z'
			Legend = list()

			#Plot 2D toroidally resolved IEDF
			im1 = ax[0].imshow(HistData2D.T, extent=extent, aspect='auto')
#			ln1 = ax[0].plot(XAxis,np.zeros(len(XAxis)), 'r--', lw=2)
			cbar1 = Colourbar(ax[0],im1,'IEDF $f(\epsilon_{i})$ [-]',5)
			ImageOptions(fig,ax[0],'',Ylabel,Title,Legend)
			
			#Plot 1D toroidally integrated IEDF								
			im2 = ax[1].plot(XAxis1D[0:-1], HistData1D, lw=2)
			cbar2 = InvisibleColourbar(ax[1])
			ImageOptions(fig,ax[1],Xlabel,'IEDF $f(\epsilon_{i})$ [-]','',Legend)
			ax[1].set_xlim(min(XAxis1D),max(XAxis1D))									#USE FIXED FULL RANGE

			#Save temporal response figure for current simulation directory
			SaveString = 'Kinetics_'+str(i)+'_'+ext
			plt.savefig(DirKinetics+SaveString)
#			plt.show()
			plt.close('all')
		#endfor
	#endfor
#endif

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
			rawKineticsData.append(rawdata)
		#endif


	#===================##===================#
	#===================##===================#

		#Extraction and organization of data from .moments files. - NOT CURRENTLY USED
		if True == False:

			#Load data from TECPLOT_KIN file and unpack into 1D array.
			rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
			rawKineticsData.append(rawdata)
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











#====================================================================#
				 	 #KINETIC IEDF JAVI-METHOD#
#====================================================================#

savefig_kineticsJAVI = False
if savefig_kineticsJAVI == True:

	#Set KStep range and interval (Should determine automatically)
	KStepMax = 2000
	KStepMod = 1000
	KStep = 1000 		#Development Input - Remove when KStep loop is added

	#Set number of cores and interval (Should determine automatically)
	NCoreTot = 64		#Not actually needed
	NCoreMod = 8

	#KSTEP LOOP GOES HERE

	#Define marker folder location and filename for current KStep
	Folder = Dir[0]+'markers/'
	Filename = 'gc_a_kstep'+str('%07.f'%KStep)+'-*.txt'
	#Sort all marker communication files numerically by core number
	CommFiles = sorted(glob.glob(Folder+Filename))
	print(Filename)

	#Initiate kinetic data array of size num variables saved in write_ep() subroutine
	#KineticsData :: 2D Array of shape [marker(n),variable(n)]
	KineticsData = np.array([]).reshape(0,10)

	#Concatenate marker data from all cores into KineticsData
	for i in tqdm( range(len(CommFiles)/NCoreMod) ):
		Data_Fragment = np.loadtxt(CommFiles[i])
		if Data_Fragment.ndim == 2:
			KineticsData = np.concatenate((KineticsData,Data_Fragment))
		elif len(Data_Fragment) == 1:
			KineticsData = np.concatenate((KineticsData,Data_Fragment.T))
		#endif
	#endfor

	#==========##==========##==========#

	#1D TEST ZONE
	if True == False:
		nBins = 100
		XData1D = KineticsData[:,3]						#'E_gc'

		#Histogram data over bin range
		Hist = np.histogram(XData1D, bins=nBins)[0]
		#Normalise distribution function
		HistSum = sum(Hist)
		NormFactor = HistSum
		Hist = [float(Hist[x])/float(NormFactor) for x in range(0,len(Hist))]
		if DebugMode == True: print( "Integrated Distribution: ",str(sum(Hist)) )

		#Create Energy axis :: Scale dEnergy from min to max and create range of len(nBins)
		Xmin,Xmax = min(XData1D),max(XData1D); 	dBin = (Xmax-Xmin)/nBins
		XAxis = [float(x)*dBin for x in range(1,nBins+1)]

		#Initiate figure and set axes
		fig,ax = figure(image_aspectratio,1)
		Xlabel,Ylabel = 'Energy $\epsilon_{i}$ [keV]', 'Ion Energy Distribution Function $f(\epsilon_{i})$ [-]'
		Legend = list()

		#Plot figure
		ax.plot(XAxis,Hist, lw=2)
		ImageOptions(fig,ax,Xlabel,Ylabel,'',Legend)
		#
		plt.show()
	#endif
	#1D TEST ZONE

	#==========##==========##==========#

	if True == True:
		#
		XData = KineticsData[:,6]						#'pphi_gc'
		Xmin,Xmax = min(XData), max(XData)
		YData = KineticsData[:,3]						#'E_gc'
		Ymin,Ymax = min(YData), max(YData)

		#
		nBins = 100
		XRange = np.linspace(Xmin,Xmax,nBins)
		YRange = np.linspace(Ymin,Ymax,nBins)

		#
		Hist,XAxis,YAxis = np.histogram2d(XData, YData, bins=(XRange, YRange))
		extent = [min(XAxis),max(XAxis), min(YAxis),max(YAxis)]

		#Normalise distribution function
		HistSum = sum(sum(Hist)); NormFactor = HistSum
		for x in range(0,len(Hist)):
			for y in range(0,len(Hist[x])):
				Hist[x,y] = float(Hist[x,y])/float(NormFactor)
			#endfor
		#endfor
		if DebugMode == True: print( "IEDF Integral: ",str(sum(Hist)) )

		#Initiate figure and set axes
		fig,ax = figure(image_aspectratio,1)
	#	Title = VariableLabels[j]+', ntor='+str(ntor)+', t='+str(Time)+' \n Simulation: '+DirString
		Xlabel,Ylabel = 'Toroidal Angular Mom. [kg m${^2}$ s$^{-1}$]', 'Energy $\epsilon_{i}$ [keV]'
		Legend = list()

		#Plot figure
		im = ax.imshow(Hist, extent=extent, aspect='auto')
		cbar = Colourbar(ax,im,'Ion Energy Distribution Function $f(\epsilon_{i})$ [-]',5)	#$\int f(\epsilon_{i}) dt$
		ImageOptions(fig,ax,Xlabel,Ylabel,'',Legend)
		#
		plt.show()
	#endif
#endif

#=====================================================================#
#=====================================================================#



