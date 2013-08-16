#!/usr/bin/env python

import numpy as N
import h5py as H
import glob as G
import matplotlib
import matplotlib.pyplot as P
from pylab import *
import scipy
import scipy.interpolate as I
from scipy import *
from scipy import optimize
import sys, os, re, shutil, subprocess, time
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to view", metavar="xxxx", default="")
parser.add_option("-M", "--min", action="store", type="float", dest="min_value", help="ignore intensities below this q-value in splined angular average (default: 0.06 A-1)", metavar="MIN_VALUE", default="0.06")
parser.add_option("-X", "--max", action="store", type="float", dest="max_value", help="ignore intensities above this q-value in splined angular average (default: 3.48 A-1)", metavar="MAX_VALUE", default="3.48")
parser.add_option("-D", "--delta", action="store", type="float", dest="delta_value", help="spline intensities with this interval in angular average (default: 0.001 A-1)", metavar="DELTA_VALUE", default="0.001")
parser.add_option("-p", "--peakfit", action="store_true", dest="peakfit", help="applies Gaussian peak fitting algorithm to the angular averages", default=False)
parser.add_option("-S", "--sonemin", action="store", type="float", dest="S1_min", help="lower limit of range used for S1 peak fitting (default: 1.50 A-1)", metavar="MIN_VALUE", default="1.50")
parser.add_option("-T", "--sonemax", action="store", type="float", dest="S1_max", help="upper limit of range used for S1 peak fitting (default: 2.08 A-1)", metavar="MAX_VALUE", default="2.08")
parser.add_option("-U", "--stwomin", action="store", type="float", dest="S2_min", help="lower limit of range used for S2 peak fitting (default: 2.64 A-1)", metavar="MIN_VALUE", default="2.64")
parser.add_option("-W", "--stwomax", action="store", type="float", dest="S2_max", help="upper limit of range used for S2 peak fitting (default: 3.20 A-1)", metavar="MAX_VALUE", default="3.20")
parser.add_option("-a", "--xaca", action="store_true", dest="xaca", help="saves the xaca files along with the angular averages", default=False)
#parser.add_option("-x", "--xcca", action="store_true", dest="xcca", help="saves the xcca files along with the angular averages", default=False)
parser.add_option("-Q", "--nq", action="store", type="int", dest="nQ", help="number of Q bins in correlation file (default: 151)", metavar="N_Q", default="151")
parser.add_option("-P", "--nphi", action="store", type="int", dest="nPhi", help="number of Phi bins in correlation file (default: 181)", metavar="N_PHI", default="181")
parser.add_option("-e", "--exclude", action="store_true", dest="exclude", help="excludes hits from splining/averaging", default=False)
parser.add_option("-f", "--excludefile", action="store", type="string", dest="excludeFile", help="name of txt file with hits to exlude (default: below100ADUs)", metavar="EXCLUDE_FILENAME", default="below100ADUs")
parser.add_option("-s", "--saveexcluded", action="store_true", dest="saveExcluded", help="flag to save average of excluded hits", default=False)
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
#source_dir = "/reg/d/psdm/cxi/cxi74613/res/cleaned_hdf5/"
#ang_avg_dir = "/reg/d/psdm/cxi/cxi74613/res/cleaned_hdf5/"
# FTC
source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
ang_avg_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"

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
		if options.peakfit:
			updateTypes[tcounter] = int(input("Spline and peak fit "+cDir+"/ (1 for yes, 0 for no)? "))
		else:
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


#Global parameters
colmax = 500
colmin = 0
storeFlag = 0

if options.peakfit:
	#Gaussian peak fitting functions
	fitfunc = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2)) + p[3]*exp(-(x-p[4])**2/(2*p[5]**2))
	errfunc = lambda p, x, y: fitfunc(p, x) - y

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
		#print "Press 'p' to save PNG."
		global colmax
		global colmin
		global p1
		if options.peakfit:
			p0 = [2.2E8, 1.83, 0.25, 1.7E8, 2.98, 0.2]
			index = N.array([((self.inangavgQ > options.S1_min)[i] and (self.inangavgQ < options.S1_max)[i]) or ((self.inangavgQ > options.S2_min)[i] and (self.inangavgQ < options.S2_max)[i]) for i in range(len(self.inangavgQ))])
			[p1, success] = optimize.leastsq(errfunc, p0[:], args=(self.inangavgQ[index],self.inangavg[index]))
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', interpolation='nearest', vmax = colmax, vmin = colmin)
		#cxi74613: invert X-axis to follow CXI-convention
		if not canvas.xaxis_inverted():
			canvas.invert_xaxis()
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		canvas = fig.add_subplot(122)
		if (options.peakfit and success):
			canvas.set_title("AngAvg; S1 = %.3f A-1, S2 = %.3f A-1" % (p1[1], p1[4]), fontsize='medium')
		else:
			canvas.set_title("Angular Average")
		
		maxAngAvg = (self.inangavg).max()
		numQLabels = len(eDD.iceHInvAngQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in eDD.iceHInvAngQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels
		
		if (options.peakfit and success):
			P.plot(self.inangavgQ, self.inangavg, "b-", self.inangavgQ, fitfunc(p1, self.inangavgQ), "g-")
		else:
			P.plot(self.inangavgQ, self.inangavg)
		
		P.xlabel("Q (A-1)")
		P.ylabel("I(Q) (ADU/srad)")
		pngtag = foundTypes[storeFlag] + '/' + "%s.png" % (self.filename)
		P.savefig(pngtag)
		print "%s saved." % (pngtag)
		P.close()
		#P.show()


########################################################
# Average and spline types
########################################################

#avgArr = N.zeros((numTypes,1760,1760)) #cxi25410
#avgArr = N.zeros((numTypes,1764,1764)) #old cxi74613
avgArr = N.zeros((numTypes,1762,1762)) #cxi74613
avgRawArr = N.zeros((numTypes,1480,1552))
avgAngAvgQ = N.arange(options.min_value,options.max_value+options.delta_value,options.delta_value)
angAvgLength = int((options.max_value-options.min_value)/options.delta_value)+1
avgAngAvg = N.zeros((numTypes,angAvgLength))
if options.xaca:
	avgCorrArr = N.zeros((numTypes,options.nQ,options.nPhi))
typeOccurences = N.zeros(numTypes)
wavelengths = [[] for i in foundTypeNumbers]
attenuations = [[] for i in foundTypeNumbers]
avgIntensities = [[] for i in foundTypeNumbers]
maxIntensities = [[] for i in foundTypeNumbers]
if options.peakfit:
	fitint1 = [[] for i in foundTypeNumbers]
	fitpos1 = [[] for i in foundTypeNumbers]
	fitfwhm1 = [[] for i in foundTypeNumbers]
	fitint2 = [[] for i in foundTypeNumbers]
	fitpos2 = [[] for i in foundTypeNumbers]
	fitfwhm2 = [[] for i in foundTypeNumbers]
damaged_events = [] #these are only saved for the non-excluded hits
failed_fits = [] #these are only saved for the non-excluded hits
if (options.exclude and options.saveExcluded):
	#excludedAvgArr = N.zeros((numTypes,1760,1760)) #cxi25410
	#excludedAvgArr = N.zeros((numTypes,1764,1764)) #old cxi74613
	excludedAvgArr = N.zeros((numTypes,1762,1762)) #cxi74613
	excludedAvgRawArr = N.zeros((numTypes,1480,1552))
	excludedAvgAngAvg = N.zeros((numTypes,angAvgLength))
	if options.xaca:
		excludedAvgCorrArr = N.zeros((numTypes,options.nQ,options.nPhi))
	excludedTypeOccurences = N.zeros(numTypes)	
	excludedWavelengths = [[] for i in foundTypeNumbers]
	excludedAttenuations = [[] for i in foundTypeNumbers]
	excludedAvgIntensities = [[] for i in foundTypeNumbers]
	excludedMaxIntensities = [[] for i in foundTypeNumbers]
	if options.peakfit:
		excludedFitint1 = [[] for i in foundTypeNumbers]
		excludedFitpos1 = [[] for i in foundTypeNumbers]
		excludedFitfwhm1 = [[] for i in foundTypeNumbers]
		excludedFitint2 = [[] for i in foundTypeNumbers]
		excludedFitpos2 = [[] for i in foundTypeNumbers]
		excludedFitfwhm2 = [[] for i in foundTypeNumbers]

#Loop over each type to spline and sum all hits
for currentlyExamining in range(numTypes):
	dirName = foundTypes[currentlyExamining]
	fcounter = 0
	dcounter = 0
	if options.peakfit:
		pcounter = 0
	if (options.exclude and options.saveExcluded):
		ecounter = 0
	numFilesInDir = len(foundTypeFiles[currentlyExamining])
	sType = set(foundTypeFiles[currentlyExamining])
	if (updateTypes[currentlyExamining] == 1):
		if (set(foundTypeFiles[currentlyExamining]).issubset(sH5)):
			if (options.exclude and os.path.exists(dirName + '/' + options.excludeFile + ".txt")):
				#Exclude hits
				excludedFiles = []
				print "Reading excluded hits from %s ..." % (dirName + '/' + options.excludeFile + ".txt")
				f = open(dirName + '/' + options.excludeFile + ".txt", 'r')
				for line in f:
					if line[0] != '#':
						line = re.sub("\n", '', line)
						excludedFiles.append(line)
				f.close()
				sExcluded = set(excludedFiles)
				print "Found %d hits to exclude." % (len(sType.intersection(sExcluded)))
				numFilesInDir = len(sType.difference(sExcluded))
			else:
				if options.exclude:
					print "%s does not exist, no hits are excluded." % (dirName + '/' + options.excludeFile + ".txt")
				sExcluded = set([])
			
			t1 = time.time()
			print "Now splining H5 files in %s/ ..." % (dirName)
			os.chdir(dirName)
			storeFlag = currentlyExamining
			for fname in foundTypeFiles[currentlyExamining]:
				diffractionName = source_dir + runtag + '/' + re.sub("-angavg",'',fname)
				correlationName = source_dir + runtag + '/' + re.sub("-angavg","-xaca",fname)
				if os.path.exists(diffractionName):
					if (os.path.exists(correlationName) or not options.xaca):
						#Read cspad-angavg.h5
						angAvgName = fname
						f = H.File(angAvgName, 'r')
						davg = N.array(f['data']['data'])
						f.close()
						if (davg[0].max() < options.max_value):
							print "Error in Q-calibration! Qmax = %s, skipping event." % (davg[0].max())
							if (fname not in sExcluded or not options.exclude):
								damaged_events.append(fname)
								dcounter += 1
								fcounter += 1
							continue
						
						#Read cspad.h5
						f = H.File(diffractionName, 'r')
						d = N.array(f['/data/data']).astype(float)
						draw = N.array(f['/data/rawdata']).astype(float)
						currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
						currAttenuation=f['LCLS']['attenuation'][0]
						f.close()
						
						#Spline angAvg
						f = I.interp1d(davg[0], davg[1])
						davg = f(avgAngAvgQ)
						
						if options.peakfit:
							#Gaussian peak fitting
							#p0 = [2.2E8, 1.83, 0.25, 1.7E8, 2.98, 0.2]
							p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2] #this helped LCLS_2011_Feb28_r0146_173206_9f8d_cspad-angavg.h5 from failing fit
							index = N.array([((avgAngAvgQ > options.S1_min)[i] and (avgAngAvgQ < options.S1_max)[i]) or ((avgAngAvgQ > options.S2_min)[i] and (avgAngAvgQ < options.S2_max)[i]) for i in range(len(avgAngAvgQ))])
							[p1, success] = optimize.leastsq(errfunc, p0[:], args=(avgAngAvgQ[index],davg[index]))
							if (success and p1[1] > options.S1_min and p1[1] < options.S1_max and p1[4] > options.S2_min and p1[4] < options.S2_max):
								if (fname not in sExcluded or not options.exclude):								
									fitint1[currentlyExamining].append(p1[0])
									fitpos1[currentlyExamining].append(p1[1])
									fitfwhm1[currentlyExamining].append(p1[2])
									fitint2[currentlyExamining].append(p1[3])
									fitpos2[currentlyExamining].append(p1[4])
									fitfwhm2[currentlyExamining].append(p1[5])
								elif options.saveExcluded:
									excludedFitint1[currentlyExamining].append(p1[0])
									excludedFitpos1[currentlyExamining].append(p1[1])
									excludedFitfwhm1[currentlyExamining].append(p1[2])
									excludedFitint2[currentlyExamining].append(p1[3])
									excludedFitpos2[currentlyExamining].append(p1[4])
									excludedFitfwhm2[currentlyExamining].append(p1[5])								
							else:
								print "The Gaussian peak fit failed, skipping %s." % (fname)
								if (fname not in sExcluded or not options.exclude):
									failed_fits.append(fname)
									if not os.path.exists("failedFits"):
										os.mkdir("failedFits")
									os.chdir(originaldir)
									currImg = img_class(d, davg, avgAngAvgQ, "failedFits/"+fname, meanWaveLengthInAngs=currWavelengthInAngs)
									currImg.draw_img_for_viewing()
									os.chdir(dirName)
									pcounter += 1
									fcounter += 1
								continue
						
						if options.xaca:
							#Read cspad-xaca.h5
							f = H.File(correlationName, 'r')
							dcorr = N.array(f['/data/data']).astype(float)	#currently saved as float (32-bit) although io->writeToHDF5() specifies it as double (64-bit)
							f.close()
						
						if (fname not in sExcluded or not options.exclude):
							#Sum data to internal arrays
							wavelengths[currentlyExamining].append(currWavelengthInAngs)
							attenuations[currentlyExamining].append(currAttenuation)
							avgIntensities[currentlyExamining].append(draw.mean())
							maxIntensities[currentlyExamining].append(max(davg))
							avgArr[currentlyExamining] += d
							avgRawArr[currentlyExamining] += draw
							avgAngAvg[currentlyExamining] += davg
							if options.xaca:
								avgCorrArr[currentlyExamining] += dcorr
							typeOccurences[currentlyExamining] += 1
							fcounter += 1
						elif options.saveExcluded:
							excludedWavelengths[currentlyExamining].append(currWavelengthInAngs)
							excludedAttenuations[currentlyExamining].append(currAttenuation)
							excludedAvgIntensities[currentlyExamining].append(draw.mean())
							excludedMaxIntensities[currentlyExamining].append(max(davg))
							excludedAvgArr[currentlyExamining] += d
							excludedAvgRawArr[currentlyExamining] += draw
							excludedAvgAngAvg[currentlyExamining] += davg
							if options.xaca:
								excludedAvgCorrArr[currentlyExamining] += dcorr
							excludedTypeOccurences[currentlyExamining] += 1
							ecounter += 1
						
					else:
						print "The correlation file %s does not exist, ignoring %s" % (correlationName, fname)
				else:
					print "The diffraction file %s does not exist, ignoring %s" % (diffractionName, fname)
				if (fname not in sExcluded or not options.exclude):
					if (options.verbose and (round(((fcounter*100)%numFilesInDir)/100)==0)):
						if options.peakfit:
							print str(fcounter) + " of " + str(numFilesInDir) + " files splined and fitted (" + str(fcounter*100/numFilesInDir) + "%)"
						else:
							print str(fcounter) + " of " + str(numFilesInDir) + " files splined (" + str(fcounter*100/numFilesInDir) + "%)"
			
			os.chdir(originaldir)
			t2 = time.time()
			if options.peakfit:
				print "Averaged " + str(fcounter - dcounter - pcounter) + " hits from type" + str(foundTypeNumbers[currentlyExamining]) + " (" + str(dcounter) + " damaged, " + str(pcounter) + " ignored due to failed peak fit)."
			else:
				print "Averaged " + str(fcounter - dcounter) + " hits from type" + str(foundTypeNumbers[currentlyExamining]) + " (" + str(dcounter) + " damaged)."
			if (options.exclude and options.saveExcluded):
				print "Averaged " + str(ecounter) + " excluded hits from type" + str(foundTypeNumbers[currentlyExamining]) + " separately."
			if (t2-t1 < 60):
				print "Time taken for averaging type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(round(t2-t1)) + " s."
			elif (t2-t1 < 3600):
				print "Time taken for averaging type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(int(t2-t1)/60) + " min " + str((round(t2-t1))%60) + " s."
			else:
				print "Time taken for averaging type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(int(t2-t1)/3600) + " h " + str(int((t2-t1)%3600)/60) + " min " + str((round(t2-t1))%60) + " s."
			if (options.verbose and len(wavelengths[currentlyExamining]) > 0):
				print "Mean wavelength = " + str(N.mean(wavelengths[currentlyExamining])) + " A."
				print "Relative change in wavelength = " + str(N.std(wavelengths[currentlyExamining])/N.mean(wavelengths[currentlyExamining]))
				print "max-min wavelength = " + str(N.max(wavelengths[currentlyExamining]) - N.min(wavelengths[currentlyExamining])) + " A."
		else:
			print "Found %d angavg files in %s/ that have not been updated. Should update all files before splining, aborting." % (len(set(foundTypeFiles[currentlyExamining])-sH5), dirName)
			sys.exit(1)


if damaged_events:
	if options.exclude:
		damaged_events_name = write_dir + runtag + "_damaged_events-" + options.excludeFile + ".txt"
	else:
		damaged_events_name = write_dir + runtag + "_damaged_events.txt"
	print "There are %s damaged events that have been ignored." % (len(damaged_events))
	N.array(damaged_events).tofile(damaged_events_name, sep="\n")
	print "Saved damaged events to %s" % (damaged_events_name)


if failed_fits:
	if options.exclude:
		failed_fits_name = write_dir + runtag + "_failed_fits-" + options.excludeFile + ".txt"
	else:
		failed_fits_name = write_dir + runtag + "_failed_fits.txt"
	print "There are %s failed fits that have been ignored." % (len(failed_fits))
	N.array(failed_fits).tofile(failed_fits_name, sep="\n")
	print "Saved failed fits to %s" % (failed_fits_name)


#print "Right-click on colorbar to set maximum scale."
#print "Left-click on colorbar to set minimum scale."
#print "Center-click on colorbar (or press 'r') to reset color scale."
#print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
#print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

storeFlag = 0

#Loop over each type to view and save average
for dirName in foundTypes:
	if (updateTypes[storeFlag] == 1):
		if (storeFlag > 0):
			if options.exclude:
				typeTag = runtag + "_type" + str(foundTypeNumbers[storeFlag]) + "-" + options.excludeFile
				if options.saveExcluded:
					excludedTypeTag = runtag + "_type" + str(foundTypeNumbers[storeFlag]) + "-" + options.excludeFile + "_excludedHits"					
			else:
				typeTag = runtag + "_type" + str(foundTypeNumbers[storeFlag])
		else:
			if options.exclude:
				typeTag = runtag + "_type0" + "-" + options.excludeFile
				if options.saveExcluded:
					excludedTypeTag = runtag + "_type0" + "-" + options.excludeFile + "_excludedHits"					
			else:
				typeTag = runtag + "_type0"
			
		if (typeOccurences[storeFlag] > 0.):
			avgArr[storeFlag] /= typeOccurences[storeFlag]
			avgRawArr[storeFlag] /= typeOccurences[storeFlag]
			avgAngAvg[storeFlag] /= typeOccurences[storeFlag]		
			if options.xaca:
				avgCorrArr[storeFlag] /= typeOccurences[storeFlag]
			currImg = img_class(avgArr[storeFlag], avgAngAvg[storeFlag], avgAngAvgQ, typeTag, meanWaveLengthInAngs=N.mean(wavelengths[storeFlag]))
			currImg.draw_img_for_viewing()
			
			if options.peakfit:
				#Gaussian peak fit statistics
				fitdeltaq = N.array(fitpos2[storeFlag]) - N.array(fitpos1[storeFlag])
				
				fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
				canvas = fig.add_subplot(131)
				canvas.set_title("Histogram")
				P.xlabel("S1 (A-1)")
				P.ylabel("Hist(S1)")
				hist_bins = N.arange(min(fitpos1[storeFlag]), max(fitpos1[storeFlag])+0.004, 0.001) - 0.0015
				pos1_hist, hist_bins = N.histogram(fitpos1[storeFlag], bins=hist_bins)
				pos1_bins = [(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(pos1_hist))]
				P.plot(pos1_bins, pos1_hist)
				P.axvline(p1[1],0,max(pos1_hist),color='r')
				
				canvas = fig.add_subplot(132)
				canvas.set_title("Histogram")
				P.xlabel("deltaQ (A-1)")
				P.ylabel("Hist(deltaQ)")
				hist_bins = N.arange(min(fitdeltaq), max(fitdeltaq)+0.004, 0.001) - 0.0015
				deltaq_hist, hist_bins = N.histogram(fitdeltaq, bins=hist_bins)
				deltaq_bins = [(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(deltaq_hist))]
				P.plot(deltaq_bins, deltaq_hist)
				P.axvline(p1[4]-p1[1],0,max(deltaq_hist),color='r')
				
				canvas = fig.add_subplot(133)
				canvas.set_title("Histogram")
				P.xlabel("S2 (A-1)")
				P.ylabel("Hist(S2)")
				hist_bins = N.arange(min(fitpos2[storeFlag]), max(fitpos2[storeFlag])+0.004, 0.001) - 0.0015
				pos2_hist, hist_bins = N.histogram(fitpos2[storeFlag], bins=hist_bins)
				pos2_bins = [(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(pos2_hist))]
				P.plot(pos2_bins, pos2_hist)
				P.axvline(p1[4],0,max(pos2_hist),color='r')
				
				pngtag = dirName +'/'+ typeTag + "-peak_fit_hist.png"
				P.savefig(pngtag)
				print "%s saved." % (pngtag)
				#P.show()
				P.close()
				
				fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
				canvas = fig.add_subplot(131)
				canvas.set_title("Correlation")
				P.xlabel("deltaQ (A-1)")
				P.ylabel("Average Intensity (ADU/pixel)")
				P.plot(fitdeltaq, avgIntensities[storeFlag], 'r.')
				
				canvas = fig.add_subplot(132)
				canvas.set_title("Correlation")
				P.xlabel("deltaQ (A-1)")
				P.ylabel("Max Intensity (ADU/srad)")
				P.plot(fitdeltaq, maxIntensities[storeFlag], 'r.')
				
				canvas = fig.add_subplot(133)
				canvas.set_title("Correlation")
				P.xlabel("deltaQ (A-1)")
				P.ylabel("Attenuation")
				P.plot(fitdeltaq, attenuations[storeFlag], 'r.')
				
				pngtag = dirName +'/'+ typeTag + "-peak_fit_corr-dq.png"
				P.savefig(pngtag)
				print "%s saved." % (pngtag)
				#P.show()
				P.close()
				
				fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
				canvas = fig.add_subplot(131)
				canvas.set_title("Correlation")
				P.xlabel("Peak1 Intensity (ADU/srad)")
				P.ylabel("Peak1 FWHM (A-1)")
				P.plot(fitint1[storeFlag], fitfwhm1[storeFlag], 'r.')
				
				canvas = fig.add_subplot(132)
				canvas.set_title("Correlation")
				P.xlabel("Peak1 Intensity (ADU/srad)")
				P.ylabel("Max Intensity (ADU/srad)")
				P.plot(fitint1[storeFlag], maxIntensities[storeFlag], 'r.')
				
				canvas = fig.add_subplot(133)
				canvas.set_title("Correlation")
				P.xlabel("Peak2 Intensity (ADU/srad)")
				P.ylabel("Peak2 FWHM (A-1)")
				P.plot(fitint2[storeFlag], fitfwhm2[storeFlag], 'r.')
				
				pngtag = dirName +'/'+ typeTag + "-peak_fit_corr-int.png"
				P.savefig(pngtag)
				print "%s saved." % (pngtag)
				#P.show()
				P.close()
		
		else:
			wavelengths[storeFlag].append([0])
			attenuations[storeFlag].append([0])
			avgIntensities[storeFlag].append([0])
			maxIntensities[storeFlag].append([0])
			if options.peakfit:
				fitint1[storeFlag].append([0])
				fitpos1[storeFlag].append([0])
				fitfwhm1[storeFlag].append([0])
				fitint2[storeFlag].append([0])
				fitpos2[storeFlag].append([0])
				fitfwhm2[storeFlag].append([0])
		
		if (options.exclude and options.saveExcluded and excludedTypeOccurences[storeFlag] > 0.):
			excludedAvgArr[storeFlag] /= excludedTypeOccurences[storeFlag]
			excludedAvgRawArr[storeFlag] /= excludedTypeOccurences[storeFlag]
			excludedAvgAngAvg[storeFlag] /= excludedTypeOccurences[storeFlag]		
			if options.xaca:
				excludedAvgCorrArr[storeFlag] /= excludedTypeOccurences[storeFlag]		
			currImg = img_class(excludedAvgArr[storeFlag], excludedAvgAngAvg[storeFlag], avgAngAvgQ, excludedTypeTag, meanWaveLengthInAngs=N.mean(excludedWavelengths[storeFlag]))
			currImg.draw_img_for_viewing()
		
		#Still save h5 file with array of zeros if no hits
		hdf5tag = dirName +'/'+ typeTag + ".h5"
		f = H.File(hdf5tag, "w")
		entry_1 = f.create_group("/data")
		entry_1.create_dataset("diffraction", data=avgArr[storeFlag])
		entry_1.create_dataset("rawdata", data=avgRawArr[storeFlag])
		if options.xaca:
			entry_1.create_dataset("correlation", data=avgCorrArr[storeFlag])
		entry_1.create_dataset("angavg", data=avgAngAvg[storeFlag])
		entry_1.create_dataset("angavgQ", data=avgAngAvgQ)
		entry_2 = f.create_group("/data/shotStatistics")
		entry_2.create_dataset("wavelength", data=wavelengths[storeFlag])
		entry_2.create_dataset("attenuation", data=attenuations[storeFlag])
		entry_2.create_dataset("intensityAvg", data=avgIntensities[storeFlag])
		entry_2.create_dataset("intensityMax", data=maxIntensities[storeFlag])
		if options.peakfit:
			entry_3 = f.create_group("/data/peakFit")
			entry_3.create_dataset("int1", data=fitint1[storeFlag])
			entry_3.create_dataset("int2", data=fitint2[storeFlag])
			entry_3.create_dataset("pos1", data=fitpos1[storeFlag])
			entry_3.create_dataset("pos2", data=fitpos2[storeFlag])
			entry_3.create_dataset("fwhm1", data=fitfwhm1[storeFlag])
			entry_3.create_dataset("fwhm2", data=fitfwhm2[storeFlag])
		if (options.exclude and options.saveExcluded and excludedTypeOccurences[storeFlag] > 0.):
			entry_4 = f.create_group("/data/excludedHits")
			entry_4.create_dataset("diffraction", data=excludedAvgArr[storeFlag])
			entry_4.create_dataset("rawdata", data=excludedAvgRawArr[storeFlag])
			if options.xaca:
				entry_4.create_dataset("correlation", data=excludedAvgCorrArr[storeFlag])
			entry_4.create_dataset("angavg", data=excludedAvgAngAvg[storeFlag])
			entry_4.create_dataset("angavgQ", data=avgAngAvgQ)
			entry_5 = f.create_group("/data/excludedHits/shotStatistics")
			entry_5.create_dataset("wavelength", data=excludedWavelengths[storeFlag])
			entry_5.create_dataset("attenuation", data=excludedAttenuations[storeFlag])
			entry_5.create_dataset("intensityAvg", data=excludedAvgIntensities[storeFlag])
			entry_5.create_dataset("intensityMax", data=excludedMaxIntensities[storeFlag])
			if options.peakfit:
				entry_6 = f.create_group("/data/excludedHits/peakFit")
				entry_6.create_dataset("int1", data=excludedFitint1[storeFlag])
				entry_6.create_dataset("int2", data=excludedFitint2[storeFlag])
				entry_6.create_dataset("pos1", data=excludedFitpos1[storeFlag])
				entry_6.create_dataset("pos2", data=excludedFitpos2[storeFlag])
				entry_6.create_dataset("fwhm1", data=excludedFitfwhm1[storeFlag])
				entry_6.create_dataset("fwhm2", data=excludedFitfwhm2[storeFlag])
		f.close()
		print "Successfully updated %s" % (hdf5tag)
		
	storeFlag += 1

colmax = 1
colmin = -1
storeFlag = 0

#Loop over each type to view correlation
if options.xaca:
	for dirName in foundTypes:
		if (typeOccurences[storeFlag] > 0.):
			if (storeFlag > 0):
				if options.exclude:
					typeTag = runtag + "_type" + str(foundTypeNumbers[storeFlag]) + "-" + options.excludeFile
					if options.saveExcluded:
						excludedTypeTag = runtag + "_type" + str(foundTypeNumbers[storeFlag]) + "-" + options.excludeFile + "_excludedHits"					
				else:
					typeTag = runtag + "_type" + str(foundTypeNumbers[storeFlag])
			else:
				if options.exclude:
					typeTag = runtag + "_type0" + "-" + options.excludeFile
					if options.saveExcluded:
						excludedTypeTag = runtag + "_type0" + "-" + options.excludeFile + "_excludedHits"					
				else:
					typeTag = runtag + "_type0"
			currImg = img_class(avgCorrArr[storeFlag], avgAngAvg[storeFlag], avgAngAvgQ, typeTag+'_xaca', meanWaveLengthInAngs=N.mean(wavelengths[storeFlag]))
			currImg.draw_img_for_viewing()
			if (options.exclude and options.saveExcluded and excludedTypeOccurences[storeFlag] > 0.):
				currImg = img_class(excludedAvgCorrArr[storeFlag], excludedAvgAngAvg[storeFlag], avgAngAvgQ, excludedTypeTag+'_xaca', meanWaveLengthInAngs=N.mean(excludedWavelengths[storeFlag]))
				currImg.draw_img_for_viewing()				
		storeFlag += 1

