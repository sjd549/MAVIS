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

radialprofiles = []						#1D Radial Profiles to be plotted (fixed Z,Phi) --
azimuthalprofiles = []					#1D Azimuthal Profiles to be plotted (fixed R,phi) |
toroidalprofiles = []					#1D Toroidal Profiles to be plotted (fixed R,Z) 
TrendLocation = [] 						#Cell location For Trend Analysis [R,Z], ([] = min/max)


#Various Diagnostic Settings.
HarmonicRange = [0,2]					#Harmonic range to be plotted [Min,Max]

#Requested diagnostics and plotting routines.
savefig_energyphys1D = True				#Plot 1D physical trends (xxx.energy_phys)
savefig_energyharm1D = True				#Plot 1D harmonic trends (xxx.energy_n)

savefig_harmonics2D = True				#Plot 2D harmonic trends (xxx.harmonics)

savefig_equilibrium2D = True			#Plot 2D ...


#Steady-State diagnostics terminal output toggles.
print_generaltrends = False				#Verbose Min/Max Trend Outputs.


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
image_1Doverlay = False					#Overlay location(s) of radialineouts/heightlineouts
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
cbaroverride = ['NotImplimented']

#=====================================================================#
#=====================================================================#






#============================#
#        ####TODO####        #
#============================#

#Stuff
#Stuff
#Stuff

# Possible Diagnostics?
#savefig_monoprofiles = False			#Single-Variables; fixed height/radius
#savefig_multiprofiles = False			#Multi-Variables; same folder
#savefig_comparelineouts = False		#Multi-Variables; all folders
#savefig_trendphaseaveraged = False		#Single-Variables; fixed cell location (or max/min)
#savefig_trendphaseresolved = True		#Single-Variables; Phase-resolved data.
#savefig_pulseprofiles = False			#Single-Variables; plotted against real-time axis

#savefig_phaseresolve1D = False			#1D Phase Resolved Images
#savefig_phaseresolve2D = False			#2D Phase Resolved Images
#savefig_PROES = False					#Simulated PROES Diagnostic

#savefig_IEDFangular = False			#2D images of angular IEDF; single folders.
#savefig_IEDFtrends = False				#1D IEDF trends; all folders.
#savefig_EEDF = False					#NO PLOTTING ROUTINE		#IN DEVELOPMENT#


























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

#Takes a 1D or 2D array and writes to a datafile in ASCII format.
#Imputs: Data, Filename, 'w'rite or 'a'ppend, and orientation (CSV or RSV).
#Example: WriteDataToFile(Image, "Filename", 'H', 'w')
def WriteDataToFile(data,filename,structure='w',Orientation='CSV'):

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

#Reads 1D or 2D data from textfile in ASCII format, returns data and header.
#Input filename, header length, data dimension and orientation (CSV or RSV).
#Example: OutputData,Header = ReadDataFromFile('/Data.txt', 1, '2D', CSV)
def ReadDataFromFile(Filename,HeaderIdx=0,Dimension='2D',Orientation='CSV'):
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

#Takes simulation folder directory (absolute path) and returns Sim128 normalisation constants
#Example: Variables,Values,Units = ReadNormConstants(Dir[l])
def ReadMEGANormalisations(Dir):

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

#=========================#

#Reads and concatenates MEGA energy.txt output files
#Takes simulation directory (absolute path) and filename (energy_n, Energy_Phys)
#Returns output data and header, data of form: [Variable][Timestamp]
#Example: OutputData,Header = ExtractMEGAEnergy(Dir[l],Filename)
def ExtractMEGAEnergy(Dir,Filename='energy_n'):

	#Extract Filename.txt paths for all SEQ for given data filename
	Files = sorted(glob.glob(Dir+'data/*'+Filename+'*'))

	#For each output file in the current simulation directory:
	for SEQ in range(0,len(Files)):
		#Extract header and output data for first SEQ
		if SEQ == 0:
			Header = ReadDataFromFile(Files[SEQ],1,'2D','CSV')[1]
			OutputData = ReadDataFromFile(Files[SEQ],1,'2D','CSV')[0]
		#Extract output data for subsequent SEQ's and append to each variable
		elif SEQ > 0:
			TempData = ReadDataFromFile(Files[SEQ],1,'2D','CSV')[0]
			for j in range(0,len(TempData)):
				OutputData[j] = np.concatenate( (OutputData[j],TempData[j]) )
			#endfor
		#endif

		#Print debugger outputs to terminal if requested
		if DebugMode == True:
			print Files[i].split('/')[-1]
			print len(OutputData), len(OutputData[0])
		#endif
	#endfor
	
	return(OutputData,Header)
#endif

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





#====================================================================#
					#COMMON PLOTTING FUNCTIONS#
#====================================================================#


#Create figure of desired size and with variable axes.
#Returns figure and axes seperately.
#fig,ax = figure(image_aspectratio,1,shareX=False)
def figure(aspectratio=[],subplots=1,shareX=False):
	if len(aspectratio) == 2:
		fig, ax = plt.subplots(subplots, figsize=(aspectratio[0],aspectratio[1]),sharex=shareX)
	else:
		fig, ax = plt.subplots(subplots, figsize=(10,10), sharex=shareX)
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

	#Force scientific notation for all axes, accounting for non-scalar x-axes.
	try: ax.xaxis.get_major_locator().set_params(style='sci',scilimits=(-2,3),axis='both')
	except: Axes_Contain_Strings = True
#	try: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='both')	#Old tickformat.
#	except: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='y')	#Old tickformat.
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

	#Crop image dimensions, use provided dimensions...
#	if isinstance(Crop, (list, np.ndarray) ) == True:
#		CropImage(ax,Extent=Crop,Rotate=Rotate)
	#... or default dimensions if not directly provided.
#	elif any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
#		if Crop == True:
#			CropImage(ax,Rotate=Rotate)
		#endif
	#endif

	#Arrange figure such that labels, legends and titles fit within frame.
	fig.tight_layout()

	return()
#enddef

#=========================#
#=========================#

#Creates and plots a colourbar with given label and binsize.
#Takes image axis, label string, number of ticks and limits
#Allows pre-defined colourbar limits in form [min,max].
#Returns cbar axis if further changes are required.
#cbar = Colourbar(ax[0],'Label',5,Lim=[0,1])
def Colourbar(ax='NaN',Label='',Ticks=5,Lim=[]):
	if ax == 'NaN': ax = plt.gca()

	#Set default font and spacing options and modify if required
	Rotation,Labelpad = 270,30
	LabelFontSize,TickFontsize = 24,18
	if '\n' in Label: Labelpad += 25		#Pad label for multi-line names

	#Create and define colourbar axis
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)

	#Set number of ticks, label location and define scientific notation.
	cbar.set_label(Label, rotation=Rotation,labelpad=Labelpad,fontsize=LabelFontSize)
	cbar.formatter.set_powerlimits((-2,3))
	cbar.locator = ticker.MaxNLocator(nbins=Ticks)
	cbar.ax.yaxis.offsetText.set(size=TickFontsize)
	yticks(fontsize=TickFontsize)

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
	if ax == 'NaN': ax = plt.gca()

	#Create colourbar axis, ideally should 'find' values of existing cbar! 
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)

	#Set new cax to zero size and remove ticks.
	try: cax.set_facecolor('none')				#matplotlib v2.x.x method
	except: cax.set_axis_bgcolor('none')		#matplotlib v1.x.x method
	for axis in ['top','bottom','left','right']:
		cax.spines[axis].set_linewidth(0)
	cax.set_xticks([])
	cax.set_yticks([])

	return(cax)
#enddef

#=========================#
#=========================#



#====================================================================#
				  #COMMON DATA ANALYSIS FUNCTIONS#
#====================================================================#


#Takes 1D or 2D array and returns array normalized to maximum value.
#If NormFactor is defined, array will be normalized to this instead.
#Returns normalized image/profile and the max/min normalization factors.
#NormProfile,Min,Max = Normalize(profile,NormFactor=0)
def Normalize(profile,NormFactor=0):
	NormalizedImage = list()

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
			NormalizedImage.append( [x/NormFactor for x in profile[i]] )
		#endfor
		profile = NormalizedImage
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

#=========================#
#=========================#

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
print '                                                 v0.0.3'
print '-------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if True in [savefig_energyphys1D,savefig_energyharm1D]:
	print'# 1D Energy Analysis'
print '-----------------------------------------'
print ''

#=====================================================================#
#=====================================================================#










#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Create lists for basic processing
Dir = list()			#List of simulation folders 					struc:[folder]
DirFiles = list()		#List of data files in each simulation folder 	struc:[folder][filenames]
NumFolders = 0			#Number of simulation folders

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

header_itermovie = list()
header_phasemovie = list()
header_IEDFlist = list()
header_2Dlist = list()

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
for l in range(0,len(HomeDirFolders)):

	#Obtain sub-directory contents (Simulation folder level)
	SubDirFolders,SubDirContents = DirectoryContents(Root+HomeDirFolders[l])

	#Determine folders containing '/data/' folder - i.e. MEGA simulation folders
	if '/data/' in SubDirFolders:
		#Add folder to global simulation list
		Dir.append(Root+HomeDirFolders[l])
		DirFiles.append(list())
		NumFolders += 1

		#Extract contents from current simulation folder and /data/ subfolder
		DataDir = Root+HomeDirFolders[l]+'/data/'
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
				Prefix = Dir[-1]+'data/'						#Note: Dir ends with a '/'
				DirFiles[-1].append(Prefix+Filename)
				#endif
			#endif
		#endfor
	else:
		CURRENT_FOLDER_IS_NOT_SIMULATION_FOLDER = 1
	#endif
#endfor

#If no folders detected, end analysis script.
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
				 	    #1D ENERGY DIAGNOSTICS#
#====================================================================#

if savefig_energyphys1D == True:

	#For each detected simulation folder
	for l in tqdm(range(0,len(Dir))):
		#Create global 1D diagnostics folder and extract current simulation name
		DirEnergy = CreateNewFolder(Dir[l],'1DEnergy_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_Phys outputs and header for plotting
		#Energy_Phys: [folder][variable][timestep]
		Energy_Phys,Header_Phys = ExtractMEGAEnergy(Dir[l],'energy_phys')

		#Extract normalisation factors for current simulation folder
		Variables,Values,Units = ReadMEGANormalisations(Dir[l])
#		print Variables[1],Values[1],Units[1]

		#==========#

		#Create fig of desired size.
		fig,ax = figure(image_aspectratio,3)

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

#==========##==========##==========#
#==========##==========##==========#

if savefig_energyharm1D == True:

	#For each detected simulation folder
	for l in tqdm(range(0,len(Dir))):
		#Create global 1D diagnostics folder and extract current simulation name
		DirEnergy = CreateNewFolder(Dir[l],'1DEnergy_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header for plotting
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGAEnergy(Dir[l],'energy_n')

		#Extract normalisation factors for current simulation folder
		Variables,Values,Units = ReadMEGANormalisations(Dir[l])
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

		#==========#

		#Create fig of desired size.
		fig,ax = figure(image_aspectratio,2)

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

if any([savefig_energyphys1D,savefig_energyharm1D]) == True:
	print '---------------------------'
	print '1D Energy Analysis Complete'
	print '---------------------------'
#endif

#====================================================================#
#====================================================================#






























































#====================================================================#
							# CODE DUMP #
#====================================================================#

# UNUSED SNIPPITS OF CODE WILL BE PUT IN HERE.
# MOST CODE IS BOUND FOR REMOVAL BUT SOME MAY STILL HAVE USE.






#=====================================================================#
#=====================================================================#




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
			Energy_n,Header_n = ExtractMEGAEnergy(Dir[l],'energy_n')
			#Energy_Phys: [folder][variable][timestep]
			Energy_Phys,Header_Phys = ExtractMEGAEnergy(Dir[l],'Energy_Phys')
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
	if any([savefig_energyphys1D,savefig_energyharm1D]) == True:
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















