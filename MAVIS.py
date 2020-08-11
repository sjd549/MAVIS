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
import numpy as np
import scipy as sp
import math as m
import subprocess
import os, sys
import os.path
import glob

#Enforce matplotlib to avoid instancing undisplayed windows
#matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
import matplotlib
matplotlib.use('Agg')

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
DebugMode = False					#Produces debug outputs for relevent diagnostics.

#Warning suppressions
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images
#Fix "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

#List of recognized data extensions for file readin
FileExtensions = ['.hamonics','moments','.txt','.in','.nam','.dat','.out']

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

#Define units for particular variables
PressureUnit = 'Torr'						#'Torr','mTorr','Pa'
BFieldUnit	=  'Gauss'						#'Gauss','Tesla'


####################

#Commonly used variable sets.
Phys = ['P-POT','TE','EF-TOT','EAMB-Z','EAMB-R','RHO','BR','BRS','BZ','BZS','BT','VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','EFLUX-R','EFLUX-Z','JZ-NET','JR-NET','TG-AVE','PRESSURE','POW-RF','POW-RF-E','POW-ICP','EB-ESORC','COLF']

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
Variables = []
MultiVar = []							#Additional variables plotted ontop of [Variables]
radialprofiles = []						#1D Radial Profiles to be plotted (fixed Z,Phi) --
azimuthalprofiles = []					#1D Azimuthal Profiles to be plotted (fixed R,phi) |
toroidalprofiles = []					#1D Toroidal Profiles to be plotted (fixed R,Z) 
TrendLocation = [] 						#Cell location For Trend Analysis [R,Z], ([] = min/max)


#Various Diagnostic Settings.
phasecycles = 1.01						#Number of waveform phase cycles to be plotted. (float)
DoFWidth = 21							#PROES Depth of Field (symmetric about image plane) (cells)
ThrustLoc = 75							#Z-axis cell for thrust calculation  (cells)
SheathROI = [34,72]						#Sheath Region of Interest, (Start,End) [cells]
SourceWidth = [16]						#Source Dimension at ROI, leave empty for auto. [cells]
EDF_Threshold = 0.01					#Upper Recognised EEDF/IEDF energy fraction (Plot all: 0.0)


#Requested diagnostics and plotting routines.
savefig_convergence = False				#Requires movie_icp.pdt
savefig_plot2D = False					#Requires TECPLOT2D.PDT

savefig_monoprofiles = False			#Single-Variables; fixed height/radius
savefig_multiprofiles = False			#Multi-Variables; same folder
savefig_comparelineouts = False			#Multi-Variables; all folders
savefig_trendphaseaveraged = False		#Single-Variables; fixed cell location (or max/min)
savefig_trendphaseresolved = True		#Single-Variables; Phase-resolved data.
savefig_pulseprofiles = False			#Single-Variables; plotted against real-time axis

savefig_phaseresolve1D = False			#1D Phase Resolved Images
savefig_phaseresolve2D = False			#2D Phase Resolved Images
savefig_PROES =	False					#Simulated PROES Diagnostic

savefig_IEDFangular = False				#2D images of angular IEDF; single folders.
savefig_IEDFtrends = False				#1D IEDF trends; all folders.
savefig_EEDF = False					#NO PLOTTING ROUTINE		#IN DEVELOPMENT#

#Write processed data to ASCII files.
write_ASCII = True						#All diagnostic output written to ASCII.


#Steady-State diagnostics terminal output toggles.
print_generaltrends = False				#Verbose Min/Max Trend Outputs.
print_Knudsennumber = False				#Print cell averaged Knudsen Number
print_soundspeed = False				#Print cell averaged sound speed
print_totalpower = False				#Print all requested total powers
print_DCbias = False					#Print DC bias at electrodeloc
print_thrust = False					#Print neutral, ion and total thrust
print_sheath = False					#Print sheath width at electrodeloc


#Image plotting options.
image_extension = '.png'				#Extensions ('.png', '.jpg', '.eps')
image_aspectratio = [10,10]				#[x,y] in cm [Doesn't rotate dynamically]
image_radialcrop = [0.65]				#[R1,R2] in cm
image_axialcrop = [1.0,4.0]				#[Z1,Z2] in cm
image_cbarlimit = []					#[min,max] colourbar limits	

image_plotsymmetry = True				#Toggle radial symmetry
image_numericaxis = False				#### NOT IMPLIMENTED ####
image_contourplot = True				#Toggle contour Lines in images
image_1Doverlay = False					#Overlay location(s) of radialineouts/heightlineouts
image_plotgrid = False					#Plot major/minor gridlines on profiles
image_plotmesh = 'PRCCP'				#Plot material mesh outlines ('Auto','PRCCP','HyperionII','EVgeny')
image_rotate = True						#Rotate image 90 degrees to the right.

image_normalize = False					#Normalize image/profiles to local max
image_logplot = False					#Plot ln(Data), against linear axis.
image_sheath = True						#Plot sheath width onto 2D images.


#Overrides the automatic image labelling.
titleoverride = []
legendoverride = []
xaxisoverride = []
xlabeloverride = []
ylabeloverride = []
cbaroverride = ['NotImplimented']

#=====================================================================#
#=====================================================================#



#============================#
#        ####TODO####        #
#============================#

# Stuff
# Stuff
# Stuff

























#====================================================================#
					#FUNDAMENTAL I/O FUNCTIONS#
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


#Takes simulation folder directory (absolute path) and returns Sim128 normalisation constants
#Example: Variables,Values,Units = ReadNormConstants(Dir)
def ReadNormConstants(Dir):

	# NOTE: Duplicated variable names in output file --- ADDRESS BY SPLITTING Sim128 FILE INTO SECTIONS
	#'D beam inj. vlc.','Crit. vlc. axis','SlowD rate axis'				--- ON TOP AND POST NORM SETS
	#'psimax','major_r','rleng','left','right','zleng','raxis','zaxis'	--- ON PRE AND POST NORM SETS

	#Normalisation constants are stored within: sim128-aug<Shot>.<num>.txt
	sim128File = sorted(glob.glob(Dir+'*sim128-aug*txt'))[0]
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

	return(Variables,Values,Units)
#enddef


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


#Takes a 1D or 2D array and writes to a datafile in ASCII format.
#Three imputs, Data to be written, Filename, 'w'rite or 'a'ppend.
#WriteDataToFile(Image, FolderNameTrimmer(Dirlist[l])+Variablelist[k])
def WriteDataToFile(data,filename,structure='w'):

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


#Reads 1D or 2D data from textfile in ASCII format.
#One input, filename string, returns data array.
def ReadDataFromFile(Filename,Dimension='1D'):
	datafile = open(Filename)
	OutputData = list()

	#Determine dimensionality of profile.
	if Dimension == '2D':
		#Read in 2D data from ASCII formatted file.	
		RawData = datafile.readlines()
		for m in range(0,len(RawData)):
			Row = RawData[m].split()
			for n in range(0,len(Row)):
				#Convert to float if possible.
				try: Row[n] = float(Row[n])
				except: Row[n] = Row[n]
			#endfor
			OutputData.append(Row)
		#endfor

	#Lowest dimention is scalar: ==> 1D array.
	elif Dimension == '1D':
		#Read in 1D data from ASCII formatted file.
		Row = datafile.readline().split()
		for m in range(0,len(Row)):
			OutputData.append(float(Row[m]))
		#endfor
	#endif

	return(OutputData)
#enddef


#Takes array of strings and compares to variable string.
#Returns true if any element of stringarray is in variable.
def IsStringInVariable(variable,stringarray):

	boolian = False
	#Check if each element of string is inside variable.
	for i in range(0,len(stringarray)):
		if stringarray[i] in variable:
			boolian = True
		#endif
	#endfor
	return(boolian)
#enddef

#=====================================================================#
#=====================================================================#



















#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Create lists for basic processing
Dir = list()			#List of simulation folders 					struc:[folder]
DirFiles = list()		#List of data files in each simulation folder 	struc:[folder][filenames]
NumFolders = 0			#Number of simulation folders

Globalvarlist = list()
Globalnumvars = list()

#Create mesh_size lists and SI conversion
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

#Lists for icp.nam variables
VRFM,VRFM2 = list(),list()
FREQM,FREQM2 = list(),list()
FREQC        = list()
FREQGLOB,IRFPOW = list(),list()
MAXFREQ,MINFREQ = list(),list()
PRESOUT = list()
IETRODEM = list()
IMOVIE_FRAMES = list()

#Lists for icp.dat variables
header_icpdat = list()			#[SpeciesName, Charge, MolecularWeight, StickingCoeff,
								# Transport, ReturnFrac, ReturnName]
AtomicSpecies = list()			#All species contained within chemistry set
FluidSpecies  = list() 			#All 'bulk' fluid species (for fluid dynamics analysis)
NeutSpecies	= list()			#All neutral and metastable species
PosSpecies = list()				#All positive ion species
NegSpecies = list()				#All negative ion species

#Lists to store raw data
rawdata_2D = list()
rawdata_kin = list()
rawdata_phasemovie = list()
rawdata_itermovie = list()
rawdata_IEDF = list()
rawdata_mcs = list()

Data = list()					#Data[folder][Variable][Datapoint]
DataIEDF = list()				#Data[folder][Variable][Datapoint]
DataEEDF = list()				#Data[folder][Variable][Datapoint]
IterMovieData = list()			#ITERMovieData[folder][timestep][variable][datapoints]
PhaseMovieData = list()			#PhaseMovieData[folder][timestep][variable][datapoints]

Moviephaselist = list()			#'CYCL = n'
MovieIterlist = list()			#'ITER = n'
EEDF_TDlist = list()			#'???'

header_itermovie = list()
header_phasemovie = list()
header_IEDFlist = list()
header_2Dlist = list()

#=====================================================================#
#=====================================================================#










#====================================================================#
					#WELCOME TEXT AND INFORMATION#
#====================================================================#

print ''
print '-------------------------------------------------------'
print ' .___  ___.      ___   ____    ____  __       _______. '
print ' |   \/   |     /   \  \   \  /   / |  |     /       | '
print ' |  \  /  |    /  ^  \  \   \/   /  |  |    |   (----` '
print ' |  |\/|  |   /  /_\  \  \      /   |  |     \   \     '
print ' |  |  |  |  /  _____  \  \    /    |  | .----)   |    '
print ' |__|  |__| /__/     \__\  \__/     |__| |_______/     '
print '                                                 v0.0.1'
print '-------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if True in [savefig_phaseresolve2D,savefig_PROES]:
	print'# 2D Phase-Resolved Movie Processing'
if True in [savefig_phaseresolve1D]:
	print'# 1D Phase-Resolved Profile Processing'
if True in [savefig_monoprofiles,savefig_multiprofiles,savefig_comparelineouts,savefig_pulseprofiles]:
	print'# 1D Steady-State Profile Processing'
if True in [print_generaltrends,print_Knudsennumber,print_soundspeed, print_totalpower,print_DCbias,print_thrust]:
	print'# 1D Specific Trend Analysis'
if savefig_trendphaseaveraged == True:
	print'# 1D Steady-State Trend Processing'
if savefig_trendphaseresolved == True:
	print'# 1D Phase-Resolved Trend Processing'
if True in [savefig_IEDFangular,savefig_IEDFtrends,savefig_EEDF]:
	print'# Angular Energy Distribution Processing'
print '-----------------------------------------'
print ''

#=====================================================================#
#=====================================================================#













#====================================================================#
					#OBTAINING FILE DIRECTORIES#
#====================================================================#

#Obtain system RAM. (and rename enviroment variable)
mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
mem_gib = mem_bytes/(1024.**3)
ext = image_extension


#Obtain home directory and contents
Root = os.path.abspath(".")
HomeDirFolders,HomeDirContents = DirectoryContents(Root)

#For each sub-directory in HomeDirFolders:
for i in range(0,len(HomeDirFolders)):

	#Obtain sub-directory contents (Simulation folder level)
	SubDirFolders,SubDirContents = DirectoryContents(Root+HomeDirFolders[i])

	#Determine folders containing '/data/' folder - i.e. MEGA simulation folders
	if '/data/' in SubDirFolders:
		#Add folder to global simulation list
		Dir.append(Root+HomeDirFolders[i])
		DirFiles.append(list())
		NumFolders += 1

		#Extract contents from current simulation folder and /data/ subfolder
		DataDir = Root+HomeDirFolders[i]+'/data/'
		DataDirContents = DirectoryContents(DataDir)[1]		#Data Folder level
		SimDirContents = SubDirContents						#Simulation Folder Level

		#Save content files from simulation folder that fit requested file extensions
		for j in range(0,len(SimDirContents)):
			Filename = SimDirContents[j]
			if any([x in Filename for x in FileExtensions]):
				Prefix = Dir[-1]
				DirFiles[-1].append(Prefix+Filename)
				#endif
			#endif
		#endfor

		#Save content files from /data/ subfolder that fit requested file extensions
		for j in range(0,len(DataDirContents)):
			Filename = DataDirContents[j]
			if any([x in Filename for x in FileExtensions]):
				Prefix = Dir[-1]+'data/'
				DirFiles[-1].append(Prefix+Filename)
				#endif
			#endif
		#endfor
	else:
		NOT_SIMULATION_FOLDER = 1
	#endif
#endfor
#If no folders detected, end analysis script.
if NumFolders == 0:
	print '-------------------------------------------'
	print 'No Ouput Files Detected, Aborting Analysis.'
	print '-------------------------------------------'
	print ''
	exit()
#endif

#==========##========================##==========#
#==========##========================##==========#










l=0

print Dir[l]

Variables,Values,Units = ReadNormConstants(Dir[l])
print Variables[1],Values[1],Units[1]
print ''

files = sorted(glob.glob(Dir[l]+'data/*energy_n*'))
print files
print ''

aux = np.loadtxt(files[0],skiprows=1)
energy_phys = np.empty((0,aux.shape[1]))
for afile in files:
  print afile.split('/')[-1]
  aux = np.loadtxt(afile,skiprows=1)
  energy_phys = np.concatenate((energy_phys,aux))
#endfor
print len(energy_phys)
print ''

exit()






























#====================================================================#
				 	 #READING DATA INTO MEMORY#
#====================================================================#

print'-----------------------'
print'Beginning Data Read-in.'
print'-----------------------'

#Extraction and organization of data from .PDT files.
for l in tqdm(range(0,numfolders)):

	#Load data from TECPLOT2D file and unpack into 1D array.
	rawdata, nn_2D = ExtractRawData(Dir,'TECPLOT2D.PDT',l)
	rawdata_2D.append(rawdata)

	#Read through all variables for each file and stop when list ends.
	Variablelist,HeaderEndMarker = ['Radius','Height'],'ZONE'
	for i in range(2,nn_2D):
		if HeaderEndMarker in str(rawdata_2D[l][i]): break
		else: Variablelist.append(str(rawdata_2D[l][i][:-2].strip(' \t\n\r\"')))
		#endif
	#endfor
	numvariables_2D,header_2D = len(Variablelist),len(Variablelist)+2
	header_2Dlist.append(header_2D)

	#Seperate total 1D data array into sets of data for each variable.
	CurrentFolderData = SDFileFormatConvertorHPEM(rawdata_2D[l],header_2D,numvariables_2D)

	#Save all variables for folder[l] to Data.
	#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
	Data.append(CurrentFolderData)


#===================##===================#
#===================##===================#

	#Kinetics data readin - NOT CURRENTLY USED
	if True == False:

		#Load data from TECPLOT_KIN file and unpack into 1D array.
		rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
		rawdata_kin.append(rawdata)
	#endif


#===================##===================#
#===================##===================#

	#IEDF/NEDF file readin.
	if True in [savefig_IEDFangular,savefig_IEDFtrends]:

		#Define arguments and autorun conv_prof.exe if possible.
		#### THIS IS HACKY, WON'T ALWAYS WORK, ARGS LIST NEEDS AUTOMATING ####
		IEDFVarArgs = ['1','1','1','1','1'] 	#Works for 2 species 1 surface.
		ExtraArgs = ['1','1','1','1','1','1','1','1','1','1']#[]	#Hack For Additional Species
		Args = ['pcmc.prof','title','1','1','1'] + IEDFVarArgs + ExtraArgs + ['0','0']
		DirAdditions = ['iprofile_tec2d.pdt','nprofile_tec2d.pdt','iprofile_tec1d.pdt', 'nprofile_tec1d.pdt','iprofile_zones_tec1d.pdt','nprofile_zones_tec1d.pdt']
		try: AutoConvProfData('./conv_prof.exe',Args,DirAdditions)
		except: print 'ConvProf Failure:'+Dirlist[l]

		#Load data from IEDFprofile file and unpack into 1D array.
		rawdata, nn_IEDF = ExtractRawData(Dir,'iprofile_tec2d.pdt',l)
		rawdata_IEDF.append(rawdata)

		#Read through all variables for each file and stop when list ends.
		IEDFVariablelist,HeaderEndMarker = ['Theta [deg]','Energy [eV]'],'ZONE'
		for i in range(2,nn_IEDF):
			#Grab EDFangle(I),EDFbins(J) values from the ZONE line, these outline the datablock size.
			if HeaderEndMarker in str(rawdata_IEDF[l][i]): 
				I = int(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[0]))
				J = int(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[1]))
				EDFangle, EDFbins = I,J
				break
			else: IEDFVariablelist.append(str(rawdata_IEDF[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		numvariables_IEDF,header_IEDF = len(IEDFVariablelist),len(IEDFVariablelist)+2
		header_IEDFlist.append(header_IEDF)

		#Seperate total 1D data array into sets of data for each variable.
		CurrentFolderData = SDFileFormatConvertorHPEM(rawdata_IEDF[l],header_IEDF,numvariables_IEDF,0,I,J)

		#Save all variables for folder[l] to Data.
		#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
		DataIEDF.append(CurrentFolderData)
	#endif


#===================##===================#
#===================##===================#

	#EEDF data readin.
	if savefig_EEDF == True:

		#Load data from MCS.PDT file and unpack into 1D array.
		rawdata, nn_mcs = ExtractRawData(Dir,'boltz_tec.pdt',l)
		rawdata_mcs.append(rawdata)

		#Unpack each row of data points into single array of floats.
		#Removing 'spacing' between the floats and ignoring variables above data.
		Energy,Fe = list(),list()
		for i in range(3,len(rawdata_mcs[l])):
			if 'ZONE' in rawdata_mcs[l][i]:
				EEDF_TDlist.append( rawdata_mcs[l][i].split('"')[-2].strip(' ') )
				DataEEDF.append([Energy,Fe])
				Energy,Fe = list(),list()
			#endif
			try:
				Energy.append( float(rawdata_mcs[l][i].split()[0]) )
				Fe.append( float(rawdata_mcs[l][i].split()[1]) )
			except:
				NaN_Value = 1
			#endtry
		#endfor
		a,b = 0,5
		for i in range(a,b):
			plt.plot(DataEEDF[i][0],DataEEDF[i][1], lw=2)
		plt.legend(EEDF_TDlist[a:b])
		plt.xlabel('Energy [eV]')
		plt.ylabel('F(e) [eV-3/2]')
		plt.show()
	#endif


#===================##===================#
#===================##===================#

	if True in [savefig_convergence,savefig_pulseprofiles]:

		#Load data from movie_icp file and unpack into 1D array.
		rawdata,nn_itermovie = ExtractRawData(Dir,'movie_icp.pdt',l)
		rawdata_itermovie.append(rawdata)

		#Read through all variables for each file and stop when list ends. 
		#movie_icp has geometry at top, therefore len(header) != len(variables).
		#Only the first encountered geometry is used to define variable zone.
		VariableEndMarker,HeaderEndMarker = 'GEOMETRY','ITER'
		variablelist,numvar = list(),0
		for i in range(2,nn_itermovie):
			if HeaderEndMarker in str(rawdata[i]): 
				header_iter = i+1		# +1 to skip to data row.
				break
			if VariableEndMarker in str(rawdata[i]) and numvar == 0:
				numvar = (i-3)	# -3 to not include R,Z and remove overflow.
			if len(rawdata[i]) > 1 and numvar == 0: 
				variablelist.append(str(rawdata_itermovie[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		header_itermovie.append(header_iter)

		#Rough method of obtaining the movie_icp iter locations for data extraction.
		Iterloc = list()
		MovieIterlist.append(list())
		for i in range(0,len(rawdata)):
			if "ITER=" in rawdata[i]:
				Iterloc.append(i+1)

				IterStart=rawdata[i].find('ITER')
				MovieIterlist[l].append(rawdata[i][IterStart:IterStart+9])
			#endif
		#endfor

		#Cycle through all iterations for current datafile, appending per cycle.
		CurrentFolderData,CurrentFolderIterlist = list(),list()
		for i in range(0,len(Iterloc)):
			if i == 0:
				CurrentIterData = SDFileFormatConvertorHPEM(rawdata,Iterloc[i],numvar+2,offset=2)
				CurrentFolderData.append(CurrentIterData[0:numvar])
			else:
				CurrentIterData = SDFileFormatConvertorHPEM(rawdata,Iterloc[i],numvar)
				CurrentFolderData.append(CurrentIterData)
			#endif
		#endfor
		IterMovieData.append(CurrentFolderData)
	#endif


#===================##===================#
#===================##===================#

	Batch=False
	if True in [savefig_phaseresolve2D,savefig_phaseresolve1D,savefig_PROES] and Batch==True:

		#Load data from movie_icp file and unpack into 1D array.
		rawdata,nn_phasemovie = ExtractRawData(Dir,'movie1.pdt',l)
		rawdata_phasemovie.append(rawdata)

		#Read through all variables for each file and stop when list ends. 
		#Movie1 has geometry at top, therefore len(header) != len(variables).
		#Only the first encountered geometry is used to define variable zone.
		VariableEndMarker,HeaderEndMarker = 'GEOMETRY','ZONE'
		variablelist,numvar = list(),0
		for i in range(2,nn_phasemovie):
			if HeaderEndMarker in str(rawdata_phasemovie[l][i]): 
				header_phase = i+2		# +2 to skip to data row.
				break
			if VariableEndMarker in str(rawdata_phasemovie[l][i]) and numvar == 0:
				numvar = (i-3)	# -3 to not include R,Z and remove overflow.
			if len(rawdata_phasemovie[l][i]) > 1 and numvar == 0: 
				variablelist.append(str(rawdata_phasemovie[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		header_phasemovie.append(header_phase)

		#Rough method of obtaining the movie1.pdt cycle locations for data extraction.
		cycleloc = list()
		for i in range(0,len(rawdata_phasemovie[l])):
			if "CYCL=" in rawdata_phasemovie[l][i]:
				cycleloc.append(i+1)
			#endif
		#endfor

		#Cycle through all phases for current datafile, appending per cycle.
		CurrentFolderData,CurrentFolderPhaselist = list(),list()
		for i in range(0,len(cycleloc)):
			if i == 0:
				CurrentPhaseData = SDFileFormatConvertorHPEM(rawdata,cycleloc[i],numvar+2,offset=2)
				CurrentFolderData.append(CurrentPhaseData[0:numvar])
			else:
				CurrentPhaseData = SDFileFormatConvertorHPEM(rawdata,cycleloc[i],numvar)
				CurrentFolderData.append(CurrentPhaseData)
			#endif
			CurrentFolderPhaselist.append('CYCL = '+str(i+1))
		#endfor
		Moviephaselist.append(CurrentFolderPhaselist)
		PhaseMovieData.append(CurrentFolderData)
	#endif
#endfor


#===================##===================#
#===================##===================#
#===================##===================#


#Create global list of all variable names and find shortest list.
for l in range(0,numfolders):
	#Alphabetize the Variablelist and keep global alphabetized list.
	tempvarlist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])[1]
	tempvarlist = sort(tempvarlist)
	numvars = len(tempvarlist)

	Globalvarlist.append(tempvarlist)
	Globalnumvars.append(numvars)
#endfor

#Find the folder with the fewest avaliable variables.
val, idx = min((val, idx) for (idx, val) in enumerate(Globalnumvars))
Comparisonlist = Globalvarlist[idx]


#===================##===================#
#===================##===================#


#Empty and delete any non-global data lists.
tempdata,tempdata2 = list(),list()
data_array,templineout = list(),list()
Energy,Fe,rawdata_mcs = list(),list(),list()
Variablelist,variablelist = list(),list()
HomeDir,DirContents = list(),list()
del RADIUS,RADIUST,HEIGHT,HEIGHTT,DEPTH,SYM
del data_array,tempdata,tempdata2,templineout
del Variablelist,variablelist
del Energy,Fe,rawdata_mcs
del HomeDir,DirContents


#Alert user that readin process has ended and continue with selected diagnostics.
if any([savefig_plot2D, savefig_phaseresolve2D, savefig_convergence, savefig_monoprofiles, savefig_multiprofiles, savefig_comparelineouts, savefig_pulseprofiles, savefig_trendphaseresolved, savefig_phaseresolve1D, savefig_PROES, savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_totalpower, print_DCbias, print_thrust, savefig_IEDFangular, savefig_IEDFtrends, savefig_EEDF]) == True:
	print '----------------------------------------'
	print 'Data Readin Complete, Starting Analysis:'
	print '----------------------------------------'
else:
	print '------------------'
	print 'Analysis Complete.'
	print '------------------'
#endif


#=====================================================================#
#=====================================================================#



























#====================================================================#
				  #COMMONLY USED PLOTTING FUNCTIONS#
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
#	mpl.rcParams['axes.prop_cycle']=cycler(color='bgrcmyk')	#Set default colour names
	mpl.rcParams['lines.linewidth'] = 1.0					#Set Default linewidth

	#Maths and Font options
	mpl.rcParams['mathtext.fontset'] = 'cm'					#Sets 'Latex-like' maths font
	mpl.rcParams['mathtext.rm'] = 'serif'					#Sets default string font

	return()
#enddef
Matplotlib_GlobalOptions()	#MUST BE RUN BEFORE ANY DIAGNOSTICS!!!!

#=========================#
#=========================#











#====================================================================#
				 #GENERAL TREND PLOTTING ANALYSIS#
#====================================================================#




#=====================================================================#
#=====================================================================#







