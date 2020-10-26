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
set_harmonicrange = [0,2]				#Harmonic range to be plotted [Min,Max]

#Requested diagnostics and plotting routines.
savefig_1Dtotalenergy = False			#Plot 1D physical trends (xxx.energy_phys)		-Javi Diagnostic
savefig_1Dspectralenergy = False		#Plot 1D harmonic trends (xxx.energy_n)			-Javi Diagnostic

savefig_2Dequilibria = False			#Plot 2D ... UNDER CONSTRUCTION
savefig_2Dharmonics = True				#Plot 2D harmonic trends 	(xxx.harmonics)
savefig_2Dfourier = False				#Plot 2D fourier analysis	(xxx.harmonics)		-Javi Diagnostic


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

#Enable choice of normalised or non-normalised units in the 1D and 2D plotting routines
#Increase user flexibility for plotting routines e.g: [min,max] harmonics, etc...
#


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
#=========================#

#Takes simulation folder directory (absolute path) and returns Sim128 normalisation constants
#Example: Variables,Values,Units = ReadNormConstants(Dir[l])
def ExtractMEGA_Normalisations(Dir):

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

	#Print debug output to terminal if requested
	if DebugMode == True:
		for i in range(0,len(Variables)): print Variables[i], Values[i], Units[i]
	#endif

	return(Variables,Values,Units)
#enddef

#=========================#
#=========================#

#Reads and concatenates MEGA energy.txt output files
#Takes simulation directory (absolute path) and filename (energy_n, Energy_Phys)
#Returns output data and header, data of form: [Variable][Timestamp]
#Example: OutputData,Header = ExtractMEGA_Energy(Dir[l],Filename)
def ExtractMEGA_Energy(Dir,Filename='energy_n'):

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

#Details on the FORTRAN file format can be found below:
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.FortranFile.read_record.html
def Extract_MEGAHarmonics(Variable,ntor,DataDir):

	def Read_MEGAHarmonics(lpsi,mpol,ntor,Filename,Variable):

		#Define FORTRANFile save data format
		#KStep (kst) is an undefined length 1D integer array		[-]
		#t (SI Time) is an undefined length 1D float array 			[IonGyroFreq*ms]
		#r_psi, gpsi_nrm, q_psi are (lpsi) length 1D float arrays 	[various]
		#All other variables are (n_elem) length 1D float arrays 	[various]
		n_elem = (mpol+1)*(2*ntor+1)*lpsi*2
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

		#Test reading a single 001 file to ensure data exists
		try:
			FORTRANFile = ff(Filename,'r')							#Name FORTRAN format file
			RawData = FORTRANFile.read_record(DataFormat)			#Extract RawData from file in Format
		except:
			print('Error: output data file ',Filename,' not found.')
			FORTRANFile.close(); exit()
		#endtry

		#Initiate output data object and set appropriate internal structures
		Data = lambda:0
		Data.kst = np.array(([]),int)									#1D KStep Array 	[-]
		Data.time = np.array(([]),float)								#1D Time Array 		[IonGyroFreq*ms]
		Data.data = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)	#1D-3D Data Arrays	[various]

		#Read each .harmonics FORTRANFile in sequence for the current simulation folder
		SEQ = 1
		while(True):
			#Attempt to read the i'th sequence (SEQ) .harmonics FORTRANFile using the defined DataFormat
			#Read attempts automatically choose the next FORTRANFile in sequence (SEQ) until none remain
			#RawData for SEQ00i is of shape RawData[Variable][datapoint] where all variables are flattened to 1D
			try: RawData = FORTRANFile.read_record(DataFormat)
			except: FORTRANFile.close(); break

			#For first SEQ=001, extract rho_pol and q_psi
			if SEQ == 1:
				Data.rho_pol = np.sqrt(abs(RawData['gpsi_nrm'][0]))			#[-]
				Data.q_psi   = RawData['q_psi'][0]							#[-]
			#endif
			#Extract 1D kstep and time arrays
			Data.kst  = np.append(Data.kst,  RawData['kst'][0])				#[-]
			Data.time = np.append(Data.time, RawData['t'][0])				#[IonGyroFreq*ms]

			#Extract all other data and reshape from 1D array into 3D array 
			#Data object for input variable [Variable] is of shape Data[kstep][mpol][ntor][lpsi][???]
			Data.data=np.concatenate((Data.data, np.reshape(RawData[Variable],(1,mpol+1,2*ntor+1,lpsi,2),order='F')))

			#Delete RawData for FORTRANFile SEQ00i and increment SEQ counter
      		del RawData
      		SEQ = SEQ+1
    	#endwhile

		#Debug outputs: 
		if DebugMode == True:
			#Print filename, KStep Range and Data Shape 
			print( '\n  '+Filename.split('/')[-1]+' Data KStep: '+str(i)+'-'+str(Data.kst[i-1]) )
			print( '  3D image shape [mpol,ntor,lpsi]: '+str(shape(Data.data[0,:,:,:,0]))      )

			#Extract 2D and 3D images for the requested input variable (Variable)
			image3D = Data.data[0,:,:,:,0]		#image3D[mpol,ntor,lpsi] for :: SEQ00i, Kstep[0]
			image2D = image3D[:,0,:]			#image2D[mpol,0,lphi] for    :: SEQ00i, Kstep[0], ntor[0]

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

		return(Data)
	#enddef
  
	#==========#==========#
	#==========#==========#


	#Extract harmonic output files for all SEQ in requested directory
	DataFiles = sorted(glob.glob(DataDir+"/*harm*"))

	#Extract data poloidal and toroidal resolutions 		
	lpsi,mpol = 201, 64														#!!! HARD-CODED FOR NOW !!!

	#Extract HarmonicData from each SEQ:  HarmonicData [kmax][mpol][ntor][lpsi] 
	HarmonicData = list()
	for i in tqdm(range(0,len(DataFiles))):
		HarmonicData.append(Read_MEGAHarmonics(lpsi,mpol,ntor,DataFiles[i],Variable))
	#endfor

	#Flatten each SEQ entry in HarmonicData, i.e. remove any sub-lists between successive SEQs.
	#LEGACY CODE - NOT SURE IF NEEDED BUT SHOULDN'T CAUSE ANY ISSUES IF DATA IS ALREADY 'FLAT'
	for i in range(0,len(HarmonicData)-1):
		aux = HarmonicData.pop(1)											#Grabs and removes 1st sub-list
		HarmonicData[0].kst  = np.append(HarmonicData[0].kst, aux.kst)		#Appends to zero'th list
		HarmonicData[0].time = np.append(HarmonicData[0].time, aux.time)
		HarmonicData[0].data = np.concatenate((HarmonicData[0].data, aux.data))
		del aux																#Refresh and grab 2nd (now 1st) sub-list
	#endfor
	HarmonicData = HarmonicData[0]			#Replace with fully appended (i.e. flattened) zero'th list

	return(HarmonicData)
#enddef

#=========================#
#=========================#

#Details on the FORTRAN file format can be found below:
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.FortranFile.read_record.html
def read_harmonics_all(ntor,dir1,kstep,seq):

  def load_harmonics_all(lpsi,mpol,ntor,afile,kstep):
    data = lambda:0
    n_elem = (mpol+1)*(2*ntor+1)*lpsi*2
    dtype = np.dtype([\
    ('kst',np.int32,1),('t',np.float64,1),\
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
    data.kst = np.array(([]),int)
    data.time = np.array(([]),float)
    data.vrad    = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.vtheta  = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.vphi    = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.brad    = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.btheta  = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.bphi    = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.erad    = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.etheta  = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.ephi    = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.prs     = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.rho     = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.dns_a   = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.mom_a   = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.ppara_a = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.pperp_a = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.qpara_a = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    data.qperp_a = np.empty(([0,mpol+1,2*ntor+1,lpsi,2]),np.float64)
    f = ff(afile,'r')
    ii = 0
    a = f.read_record(dtype)
    while(ii<=kstep):
      try:
        a = f.read_record(dtype)
      except IOError:
        print("IO error detected, breaking...")
        f.close()
        break
      except TypeError:
        print("Type Error detected, breaking...")
        f.close()
        break
      data.kst     = np.append(data.kst,a['kst'][0])
      data.time    = np.append(data.time,a['t'][0])#*1e3/wa
      print(str(ii)+'-'+str(data.kst[ii]))
      if ii == 0:
        data.rho_pol = np.sqrt(abs(a['gpsi_nrm'][0]))
        data.q_psi   = a['q_psi'][0]
      elif ii==kstep:
        data.vrad    = np.reshape(a['vrad'   ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.vtheta  = np.reshape(a['vtheta' ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.vphi    = np.reshape(a['vphi'   ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.brad    = np.reshape(a['brad'   ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.btheta  = np.reshape(a['btheta' ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.bphi    = np.reshape(a['bphi'   ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.erad    = np.reshape(a['erad'   ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.etheta  = np.reshape(a['etheta' ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.ephi    = np.reshape(a['ephi'   ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.prs     = np.reshape(a['prs'    ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.rho     = np.reshape(a['rho'    ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.dns_a   = np.reshape(a['dns_a'  ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.mom_a   = np.reshape(a['mom_a'  ],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.ppara_a = np.reshape(a['ppara_a'],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.pperp_a = np.reshape(a['pperp_a'],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.qpara_a = np.reshape(a['qpara_a'],(mpol+1,2*ntor+1,lpsi,2),order='F')      
        data.qperp_a = np.reshape(a['qperp_a'],(mpol+1,2*ntor+1,lpsi,2),order='F')      
      del a
      ii = ii+1
    return data
  
  #dir1 = os.getcwd()
  files = sorted(glob.glob(dir1+"/*harm*"))
  lpsi = 201
  mpol = 64
  data = load_harmonics_all(lpsi,mpol,ntor,files[seq],kstep)
  return data
#enddef

#=========================#
#=========================#






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
#	mpl.rcParams['axes.prop_cycle']=cycler(color='bgrcmyk')	#Set default colour names
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
def Colourbar(ax='NaN',image='NaN',Label='',Ticks=5,Lim=[]):

	#Determine supplied image and axis, replacing or breaking if required
	if image == 'NaN': print('Colourbar Image Not Supplied')
	if ax == 'NaN': ax = plt.gca()

	#Set default font and spacing options and modify if required
	Rotation,Labelpad = 270,30
	LabelFontSize,TickFontsize = 24,18
	if '\n' in Label: Labelpad += 25		#Pad label for multi-line names

	#Create and define colourbar axis
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)
	cbar = plt.colorbar(image, cax=cax)

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


#Compute upper and lower Alfven eigenmode threshold frequencies
#UpperThreshold,LowerThreshold = TAEThresholds(Harmonic,mpol,lpsi,qpsi,rho_pol,eps,AlfvenVelocity,subfig)
def TAEThresholds(HarmonicData,Harmonic,eps,AlfvenVelocity,ax='NaN'):

	#Extract required data
	data = HarmonicsData.data
	kmax = data.shape[0]
	mpol = data.shape[1]
	ntor = data.shape[2]
	lpsi = data.shape[3]

	#Initiate TAE threshold arrays
	UpperThresholds = list()
	LowerThresholds = np.zeros([lpsi,mpol-1])

	#Extract safety factor (???) array and initiate threshold arrays
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
			Yaxis = UpperThresholds[m]*AlfvenVelocity/(2*np.pi)/1000
			ax.plot(rho_pol,Yaxis, 'w--', lw=1.5)
		#endfor
		Yaxis = np.amin(LowerThresholds,axis=1)*AlfvenVelocity/(2*np.pi)/1000
		subfig.plot(rho_pol,Yaxis, 'w--', lw=1.5)
	#endif

	return(UpperThresholds,LowerThresholds)
#endfor

#=========================#
#=========================#

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

#=====================================================================#
#=====================================================================#

























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
print('                                                 v0.0.5')
print('-------------------------------------------------------')
print('')
print('The following diagnostics were requested:')
print('-----------------------------------------')
if True in [savefig_1Dtotalenergy,savefig_1Dspectralenergy]:
	print('# 1D Energy Analysis')
if True in [savefig_2Dequilibria]:
	print('# 2D Equilibria Plots')
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
				 	    #1D ENERGY DIAGNOSTICS#
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

#==========##==========##==========#
#==========##==========##==========#

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
				 	    #2D EQUILIBRIA DIAGNOSTICS#
#====================================================================#

if savefig_2Dequilibria == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		print Dir[l]

		#Create global 1D diagnostics folder and extract current simulation name
		DirHarmonics = CreateNewFolder(Dir[l],'2DHarmonic_Profiles/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		NumHarmonics = len(Energy_n)-3		#Ignore kstep, time, n=0

		#Extract Harmonics outputs for plotting, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.data: [4D Array] of shape [kstep][mpol][ntor][lpsi] for variables[i]
		HarmonicsData = Extract_MEGAHarmonics('bphi',NumHarmonics,Dir[l]+'data/')
		## ??? PROBABLY NEED .moments DATA FOR EQUILIBRIA PLOTS ??? ##

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		AlfvenVelocity = Values[Variables.index('Alfven velocity')] #B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = Values[Variables.index('D gyro frequency')]
		IonDensity = Values[Variables.index('Bulk density')]
		B0 = Values[Variables.index('Mag.fld. at axis')]
		R0 = Values[Variables.index('raxis')]
		m_D = 3.34e-27
		eps = 0.5/R0
	#endfor

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
		DirHarmonics = CreateNewFolder(Dir[l],'2DHarmonic_Images/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		NumHarmonics = len(Energy_n)-3		#Ignore kstep, time, n=0

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		AlfvenVelocity = Values[Variables.index('Alfven velocity')] #B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = Values[Variables.index('D gyro frequency')]
		IonDensity = Values[Variables.index('Bulk density')]
		B0 = Values[Variables.index('Mag.fld. at axis')]
		R0 = Values[Variables.index('raxis')]
		m_D = 3.34e-27
		eps = 0.5/R0



		#TESTING - REMOVE ME!!!
		variables = ['bphi']
		set_harmonicrange = [0,4]
		#TESTING - REMOVE ME!!!



		#For each requested variable
		for i in range(0,len(variables)):

			#Create new folder for each variable
			DirVariable = CreateNewFolder(DirHarmonics,variables[i])

			#Extract Harmonics outputs for plotting, it contains:
			#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
			#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
			#HarmonicsData.data: [4D Array] of shape [kstep][mpol][ntor][lpsi] for variables[i]
			HarmonicsData = Extract_MEGAHarmonics(variables[i],NumHarmonics,Dir[l]+'data/')

			#HarmonicsData.data is of shape [kstep][mpol][ntor][lpsi][???]
			HarmonicsData.data = HarmonicsData.data[:,:,:,:,0]		#Final index is unusual, should fix this

			#Extract Variablelabel				
			VariableLabel = variables[i]							#Needs a seperate function to label vars

			#Set TimeIndex and employ to extract KStep and Time
			TimeIndex = 600											#!!! Hard-coded to final entry for now !!!
			KStep = HarmonicsData.kst[TimeIndex]					#Need a method of altering these
			Time = (HarmonicsData.time[TimeIndex]/IonGyroFreq)*1e3	#Need a method of altering these

			#Extract normalised axes and set image extent
			## !!! MAKE THESE A FUNCTION extent = ImageExtent() !!! ##		
			rho_pol = HarmonicsData.rho_pol							#Normalised poloidal angle
			phi_tor = range(1,shape(HarmonicsData.data)[1]+1)		#Array of toroidal cells
			phi_tor = [x/float(len(phi_tor)) for x in phi_tor]		#Normalised toroidal angle
			extent = [rho_pol[0],rho_pol[-1], phi_tor[0],phi_tor[-1]]

			#Plot images of variables[i] for all ntor
			for j in range(set_harmonicrange[0],set_harmonicrange[1]):

				#Extract image from HarmonicsData
				ntor = j											#ntor is symmetrical, needs fixing
				Image = HarmonicsData.data[TimeIndex,:,ntor,:]
				print shape(HarmonicsData.data)

				#Create figure and define Title, Legend, Axis Labels etc...
				fig,ax = figure(image_aspectratio,1)
				Title = VariableLabel+' for ntor='+str(ntor)+' at KStep='+str(KStep)+' \n Simulation: '+DirString
				Xlabel,Ylabel = 'Poloidal Extent $\\rho_{pol}$ [-]', 'Toroidal Extent $\phi_{tor}$ [-]'
#				Xlabel,Ylabel = 'Poloidal Resolution $m_{\\theta}$ [cells]', 'Toroidal Resolution $l_{\phi}$ [cells]'
				Legend = list()

				#Plot example data for debugging purposes
				im = ax.imshow(Image, extent=extent, aspect='auto', origin='bottom')
				cbar = Colourbar(ax,im,VariableLabel+' [-]',5)
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend)

				#Save 2D harmonics figure for current simulation
				plt.savefig(DirVariable+VariableLabel+'_n'+str(ntor)+ext)
#				plt.show()
				plt.close('all')
			#endfor
		#endfor
	#endfor

	if any([savefig_2Dharmonics]) == True:
		print '--------------------------'
		print '2D Harmonic Plots Complete'
		print '--------------------------'
	#endif
#endif












#====================================================================#
				 	    #2D FOURIER ANALYSIS#
#====================================================================#

if savefig_2Dfourier == True:

	#For each detected simulation folder
	for l in range(0,len(Dir)):

		print Dir[l]

		#Create global 2D diagnostics folder and extract current simulation name
		DirHarmonics = CreateNewFolder(Dir[l],'2DHarmonic_Images/')
		DirString = Dir[l].split('/')[-2]
		SubString = DirString.split('_')[-1]

		#Extract Energy_n outputs and header, used to determine harmonic range
		#energy_n: [folder][variable][timestep]
		Energy_n,Header_n = ExtractMEGA_Energy(Dir[l],'energy_n')
		NumHarmonics = len(Energy_n)-3		#Ignore kstep, time, n=0

		#Extract Harmonics outputs for plotting, it contains:
		#HarmonicsData.rho_pol [1D array] :: HarmonicsData.q_psi [1D array]
		#HarmonicsData.kst [1D array]     :: HarmonicsData.time [1D array]    
		#HarmonicsData.data: [4D Array] of shape [kstep][mpol][ntor][lpsi] for variables[i]
		HarmonicsData = Extract_MEGAHarmonics('bphi',NumHarmonics,Dir[l]+'data/')

		#Extract relevant normalisation factors for current simulation folder
		Variables,Values,Units = ExtractMEGA_Normalisations(Dir[l])
		AlfvenVelocity = Values[Variables.index('Alfven velocity')] #B0/np.sqrt(4e-7*np.pi*IonDensity*m_D)
		IonGyroFreq = Values[Variables.index('D gyro frequency')]
		IonDensity = Values[Variables.index('Bulk density')]
		B0 = Values[Variables.index('Mag.fld. at axis')]
		R0 = Values[Variables.index('raxis')]
		m_D = 3.34e-27
		eps = 0.5/R0



		#TO STILL BE TRANSLATED 
		Data = HarmonicsData.data
		dt = (HarmonicsData.time[1]-HarmonicsData.time[0])/IonGyroFreq*1e3
		kmax = Data.shape[0]
		mpol = Data.shape[1]
		ntor = Data.shape[2]
		lpsi = Data.shape[3]
		ltheta = 256
		ntor2 = int(0.5*(ntor-1))

		print kmax, mpol, ntor, lpsi, ntor2
		#TO STILL BE TRANSLATED
		#ALSO NEED TO ADD COLOURBAR TO THE FOURIER PLOTS!!!



		vcos = list()
		#Create Vcos list and add zeros equal to the toroidal resolution?
		for n in range(0,ntor2):
			vcos.append( np.zeros([]) )
			#Extract Bpol data for each respective toroidal (n) and poloidal (m) cells (adding vcos zeros?)
			for m in range(0,mpol):
				vcos[n] = vcos[n] + Data[:,m,n,:,0]	#data structure: [kmax][mpol][ntor][lpsi]  
			#endfor
			vcos[n][np.isnan(vcos[n])] = 0			#Remove any NaNs
		#endfor

		vcos_fft,vcos_len = list(),list()
		#Extract fourier components from vcos
		for n in range(0,ntor2):
			vcos_fft.append( np.fft.fft(vcos[n],axis=0) )	# 											???
		  	vcos_fft[n][0,:] = vcos_fft[n][0,:]*0.0			#Discard imaginary components 				???
			vcos_len.append( int(len(vcos_fft[n])/2)-1 )	#Determine lowpass filter threshold 		???
			vcos_fft[n] = vcos_fft[n][0:vcos_len[n],:]		#Discard upper frequencies (lowpass filter)
		#endfor

		#==========##==========#

		#Create fig of desired size - increasing Xlim with the number of harmonics
		Xlim,Ylim = 10*(NumHarmonics/2), 10
		fig,ax = figure([Xlim,Ylim],[2,ntor2])

		#For each toroidal harmonic:
		for i in range(0,ntor2):

			#Temporal evolution plotted on the top row (row 0)
			if ntor2 == 1: subfig = ax[0]
			elif ntor2 > 1: subfig = ax[0,i]
			Harmonic = -ntor2+i									#Why is it reversed?
			#Construct meshgrid
			Xaxis = HarmonicsData.rho_pol						#[-]
			Yaxis = (HarmonicsData.time/IonGyroFreq)*1e3		#[ms]
			X,Y = np.meshgrid(Xaxis,Yaxis)
			#Plot harmonic temporal evolution
			subfig.contourf(X,Y,vcos[i], cmap='plasma')
			if i == 0: ImageOptions(fig,subfig,'','Time [ms]','n='+str(Harmonic),'')
			if i > 0: ImageOptions(fig,subfig,'','','n='+str(Harmonic),'')

			#==========#

			#Fourier analysis plotted on the bottom row (row 1)
			if ntor2 == 1: subfig = ax[1]
			elif ntor2 > 1: subfig = ax[1,i]
			#Construct meshgrid
			Xaxis = HarmonicsData.rho_pol						#[-]
			Yaxis = np.linspace(0,0.5/dt,vcos_len[i])			#[kHz]
			X,Y = np.meshgrid(Xaxis,Yaxis)
			#Plot fourier amplitude spectrum
			subfig.contourf(X,Y,vcos_fft[i], cmap='plasma')
			if i == 0: ImageOptions(fig,subfig,'Poloidal Angle $\\rho_{pol}$','Frequency [kHz]','','')
			if i > 0: ImageOptions(fig,subfig,'Poloidal Angle $\\rho_{pol}$','','','')
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

if any([savefig_2Dharmonics]) == True:
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















