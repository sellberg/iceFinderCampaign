#!/usr/bin/env python

import numpy as N
import h5py as H
import glob as G
import matplotlib
import matplotlib.pyplot as P
import scipy
import scipy.interpolate as I
import sys, os, re, shutil, subprocess, time
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to view", metavar="xxxx", default="")
parser.add_option("-M", "--min", action="store", type="float", dest="min_value", help="ignore intensities below this q-value in splined angular average (default: 0.06 A-1)", metavar="MIN_VALUE", default="0.06")
parser.add_option("-X", "--max", action="store", type="float", dest="max_value", help="ignore intensities above this q-value in splined angular average (default: 3.48 A-1)", metavar="MAX_VALUE", default="3.48")
parser.add_option("-D", "--delta", action="store", type="float", dest="delta_value", help="spline intensities with this interval in angular average (default: 0.001 A-1)", metavar="DELTA_VALUE", default="0.001")
parser.add_option("-a", "--xaca", action="store_true", dest="xaca", help="copies the xaca files along with the angular averages", default=False)
#parser.add_option("-x", "--xcca", action="store_true", dest="xcca", help="copies the xcca files along with the angular averages", default=False)
parser.add_option("-P", "--nq", action="store", type="int", dest="nQ", help="number of Q bins in correlation file (default: 151)", metavar="N_Q", default="151")
parser.add_option("-Q", "--nphi", action="store", type="int", dest="nPhi", help="number of Phi bins in correlation file (default: 181)", metavar="N_PHI", default="181")
parser.add_option("-o", "--outputdir", action="store", type="string", dest="outputDir", help="output directory (default: output_rxxxx)", metavar="OUTPUT_DIR", default="output")  
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="prints out the frame number as it is processed", default=False)
(options, args) = parser.parse_args()

########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
# SCRATCH
#source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
#ang_avg_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
# RES
#source_dir = "/reg/data/ana12/cxi/cxi74613/res/cleaned_hdf5/"
#ang_avg_dir = "/reg/data/ana12/cxi/cxi74613/res/cleaned_hdf5/"
# FTC
source_dir = "/reg/data/ana12/cxi/cxi74613/ftc/cleaned_hdf5/"
ang_avg_dir = "/reg/data/ana12/cxi/cxi74613/ftc/cleaned_hdf5/"

runtag = "r%s"%(options.runNumber)
write_dir = options.outputDir + '_' + runtag + '/'

h5files = []
foundFiles = []
#Check if write_dir already exists, and which angavg files have already been sorted
if not os.path.exists(write_dir):
	print "There is no directory called %s. Aborting." % (write_dir)
	sys.exit(1)	
else:
	write_anomaly_dir = write_dir 
	originaldir = os.getcwd()
	foundTypes = [options.outputDir + '_' + runtag]
	anomalousTypes = G.glob(write_anomaly_dir+"type[1-9]")
	if (len(anomalousTypes) > 0):
		foundTypes += anomalousTypes
	numTypes = len(foundTypes)
	updateTypes = N.zeros(numTypes, dtype='int64')
	foundTypeNumbers = N.arange(numTypes)
	foundTypeFiles = [[] for i in foundTypeNumbers]
	for i in foundTypeNumbers:
		typeInStringPos = foundTypes[i].find("type")
		if (typeInStringPos == -1):
			foundTypeNumbers[i] = 0
		else:
			foundTypeNumbers[i] = int(foundTypes[i][typeInStringPos+len("type")])
	tcounter = 0 
	for cDir in foundTypes:
		os.chdir(cDir)
		updateTypes[tcounter] = int(input("Spline "+cDir+"/ (1 for yes, 0 for no)? "))
		foundTypeFiles[tcounter] += G.glob("LCLS*angavg.h5")
		foundFiles += G.glob("LCLS*angavg.h5")
		os.chdir(originaldir)
		tcounter += 1
	print "Found %d types (including type0) with %d sorted files." % (numTypes, len(foundFiles))	

if (options.verbose):
	print "Detector at %lf mm away from the interaction region." % (eDD.get_detector_dist_in_meters(runtag)*1000)

#Check if searchDir already exists and what angavg files it includes
searchDir = ang_avg_dir + runtag
if not os.path.exists(searchDir):
	print "There is no directory called %s. Aborting." % (searchDir)
	sys.exit(1)	
else:
	print "Now examining new H5 files in %s/ ..." % (searchDir)   
	searchstring="LCLS+[a-zA-Z0-9\_]+"+runtag+"[a-z0-9\_]+-angavg.h5"
	h5pattern = re.compile(searchstring)
	h5files = [x for x in os.listdir(searchDir) if h5pattern.findall(x)]
	numFiles = len(h5files)
	print "Found %d new H5 files." % (numFiles)

#Check that there are no duplicates in the presorted files
sFound = set(foundFiles)
sH5 = set(h5files)
if (len(sFound) != len(foundFiles)):
	print "Duplicate files exist in the pre-sorted directories. Aborting."
	sys.exit(1)

#Check that all files in searchDir are included in the presorted files
if (sH5.issubset(sFound)):
	print "There are %d pre-sorted files in %s (and its subfolders) that will not get updated." % (len(sFound-sH5), write_dir)
	print "Updating %d pre-sorted files..." % (numFiles)
	fcounter = 0
	for fname in h5files:
		angAvgName = ang_avg_dir + runtag + '/' + fname
		diffractionName = source_dir + runtag + '/' + re.sub("-angavg",'',fname)
		correlationName = source_dir + runtag + '/' + re.sub("-angavg","-xaca",fname)
		if os.path.exists(diffractionName):
			if os.path.exists(correlationName) or not options.xaca:
				#Check if fname is part of type0
				if (fname in set(foundTypeFiles[0])):
					shutil.copyfile(angAvgName, write_dir+fname)
					if options.xaca:
						shutil.copyfile(correlationName, write_dir+re.sub("-angavg","-xaca",fname))
					fcounter += 1
					if (options.verbose and (round(((fcounter*100)%numFiles)/100)==0)):
						print str(fcounter) + " of " + str(numFiles) + " files copied (" + str(fcounter*100/numFiles) + "%)"
				else:
					for i in range(1, numTypes):
						if (fname in set(foundTypeFiles[i])):
							shutil.copyfile(angAvgName, write_dir+'type'+str(foundTypeNumbers[i])+'/'+fname)
							shutil.copyfile(diffractionName, write_dir+'type'+str(foundTypeNumbers[i])+'/'+re.sub("-angavg",'',fname))
							if options.xaca:
								shutil.copyfile(correlationName, write_dir+'type'+str(foundTypeNumbers[i])+'/'+re.sub("-angavg","-xaca",fname))
							fcounter += 1
							if (options.verbose and (round(((fcounter*100)%numFiles)/100)==0)):
								print str(fcounter) + " of " + str(numFiles) + " files copied (" + str(fcounter*100/numFiles) + "%)"
			else:
				print "The correlation file %s does not exist, ignoring %s" % (correlationName, fname)
		else:
			print "The diffraction file %s does not exist, ignoring %s" % (diffractionName, fname)
	print "Updated %d of %d files successfully." % (fcounter, numFiles)
else:
	print "Found %d angavg files in %s/ that are not already sorted. Should sort files first with ./separateAndShowStrongHits.py and ./fineInspection.py, aborting." % (len(sH5-sFound), searchDir)
	sys.exit(1)


#Global parameters
colmax = 1000
colmin = 0
storeFlag = 0

#########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
#########################################################
class img_class (object):
	def __init__(self, inarr, inangavg, inangavgQ , filename, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs, detectorDistance=eDD.get_detector_dist_in_meters(runtag)):
		self.inarr = inarr
		self.filename = filename
		self.inangavg = inangavg
		self.inangavgQ = inangavgQ
		self.wavelength = meanWaveLengthInAngs
		self.detectorDistance = detectorDistance
		global colmax
		global colmin
		global storeFlag
		self.tag = 0
	
	def on_keypress_for_viewing(self,event):
		global colmax
		global colmin
		global storeFlag
		if event.key == 'p':
			pngtag = foundTypes[storeFlag] + '/' + "%s.png" % (self.filename)
			P.savefig(pngtag)
			print "%s saved." % (pngtag)
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = self.inarr.max()
			P.clim(colmin, colmax)
			P.draw()
	
	def on_click(self, event):
		global colmax
		global colmin
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin = self.inarr.min()
				colmax = self.inarr.max()
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			P.clim(colmin, colmax)
			P.draw()
	
	def draw_img_for_viewing(self):
		print "Press 'p' to save PNG."
		global colmax
		global colmin
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		canvas = fig.add_subplot(122)
		canvas.set_title("Angular Average")
		maxAngAvg = (self.inangavg).max()
		numQLabels = len(eDD.iceHInvAngQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in eDD.iceHInvAngQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels
		
		P.plot(self.inangavgQ, self.inangavg)
		P.xlabel("Q (A-1)")
		P.ylabel("I(Q) (ADUs/srad)")
		P.show()


########################################################
# Average and spline types
########################################################

#avgArr = N.zeros((numTypes,1760,1760)) #cxi25410
#avgArr = N.zeros((numTypes,1764,1764)) #old cxi74613
avgArr = N.zeros((numTypes,1762,1762)) #cxi74613
avgRawArr = N.zeros((numTypes,1480,1552))
avgCorrArr = N.zeros((numTypes,options.nQ,options.nPhi))
avgAngAvgQ = N.arange(options.min_value,options.max_value+options.delta_value,options.delta_value)
angAvgLength = int((options.max_value-options.min_value)/options.delta_value)+1
avgAngAvg = N.zeros((numTypes,angAvgLength))
typeOccurences = N.zeros(numTypes)
wavelengths = [[] for i in foundTypeNumbers]

#Loop over each type to spline and sum all hits
for currentlyExamining in range(numTypes):
	dirName = foundTypes[currentlyExamining]
	fcounter = 0
	numFilesInDir = len(foundTypeFiles[currentlyExamining])
	if (updateTypes[currentlyExamining] == 1):
		if (set(foundTypeFiles[currentlyExamining]).issubset(sH5)):
			t1 = time.time()
			print "Now splining H5 files in %s/ ..." % (dirName)
			os.chdir(dirName)
			for fname in foundTypeFiles[currentlyExamining]:
				storeFlag = foundTypeNumbers[currentlyExamining]
				if (options.verbose and (round(((fcounter*100)%numFilesInDir)/100)==0)):
					print str(fcounter) + " of " + str(numFilesInDir) + " files splined (" + str(fcounter*100/numFilesInDir) + "%)"
				diffractionName = source_dir + runtag + "/" + re.sub("-angavg",'',fname)
				correlationName = source_dir + runtag + "/" + re.sub("-angavg","-xaca",fname)
				if os.path.exists(diffractionName):
					if os.path.exists(correlationName) or not options.xaca:
						f = H.File(diffractionName, 'r')
						d = N.array(f['/data/data']).astype(float)
						draw = N.array(f['/data/rawdata']).astype(float)
						currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
						f.close()
						if options.xaca:
							f = H.File(correlationName, 'r')
							dcorr = N.array(f['/data/data']).astype(float)	#currently saved as float (32-bit) although io->writeToHDF5() specifies it as double (64-bit)
							f.close()
							avgCorrArr[currentlyExamining] += dcorr
						angAvgName = fname
						f = H.File(angAvgName, 'r')
						davg = N.array(f['data']['data'])
						f.close()
						f = I.interp1d(davg[0], davg[1])
						davg = f(avgAngAvgQ)
						wavelengths[currentlyExamining].append(currWavelengthInAngs)
						avgArr[currentlyExamining] += d
						avgRawArr[currentlyExamining] += draw
						avgAngAvg[currentlyExamining] += davg
						typeOccurences[currentlyExamining] += 1
						fcounter += 1
					else:
						print "The correlation file %s does not exist, ignoring %s" % (correlationName, fname)
				else:
					print "The diffraction file %s does not exist, ignoring %s" % (diffractionName, fname)
			os.chdir(originaldir)
			t2 = time.time()
			print "Time taken for averaging type" + str(storeFlag) + " = " + str(t2-t1) + " s."
			if (options.verbose):
				print "Mean wavelength = " + str(N.mean(wavelengths[currentlyExamining])) + " A."
				print "Relative change in wavelength = " + str(N.std(wavelengths[currentlyExamining])/N.mean(wavelengths[currentlyExamining]))
				print "max-min wavelength = " + str(N.max(wavelengths[currentlyExamining]) - N.min(wavelengths[currentlyExamining])) + " A."
		else:
			print "Found %d angavg files in %s/ that have not been updated. Should update all files before splining, aborting." % (len(set(foundTypeFiles[currentlyExamining])-sH5), dirName)
			sys.exit(1)


print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."


storeFlag = 0

#Loop over each type to view and save average
for dirName in foundTypes:
	if (typeOccurences[storeFlag] > 0.):
		avgArr[storeFlag] /= typeOccurences[storeFlag]
		avgRawArr[storeFlag] /= typeOccurences[storeFlag]
		avgAngAvg[storeFlag] /= typeOccurences[storeFlag]		
		if options.xaca:
			avgCorrArr[storeFlag] /= typeOccurences[storeFlag]
		if (storeFlag > 0):
			typeTag = runtag+'_type'+str(foundTypeNumbers[storeFlag])
		else:
			typeTag = runtag+'_type0'
		currImg = img_class(avgArr[storeFlag], avgAngAvg[storeFlag], avgAngAvgQ, typeTag, meanWaveLengthInAngs=N.mean(wavelengths[storeFlag]))
		currImg.draw_img_for_viewing()
		f = H.File(dirName +'/'+ typeTag + ".h5", "w")
		entry_1 = f.create_group("/data")
		entry_1.create_dataset("diffraction", data=avgArr[storeFlag])
		entry_1.create_dataset("rawdata", data=avgRawArr[storeFlag])
		if options.xaca:
			entry_1.create_dataset("correlation", data=avgCorrArr[storeFlag])
		entry_1.create_dataset("angavg", data=avgAngAvg[storeFlag])	
		entry_1.create_dataset("angavgQ", data=avgAngAvgQ)	
		f.close()
		print "Successfully updated %s" % (dirName +'/'+ typeTag + ".h5")
	storeFlag += 1


colmax = 1
colmin = -1
storeFlag = 0

#Loop over each type to view correlation
if options.xaca:
	for dirName in foundTypes:
		if (typeOccurences[storeFlag] > 0.):
			if (storeFlag > 0):
				typeTag = runtag+'_type'+str(foundTypeNumbers[storeFlag])
			else:
				typeTag = runtag+'_type0'
			currImg = img_class(avgCorrArr[storeFlag], avgAngAvg[storeFlag], avgAngAvgQ, typeTag+'_xaca', meanWaveLengthInAngs=N.mean(wavelengths[storeFlag]))
			currImg.draw_img_for_viewing()
		storeFlag += 1

