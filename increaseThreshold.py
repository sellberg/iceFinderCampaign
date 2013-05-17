#!/usr/bin/env python

import numpy as N
from numpy import linalg as LA
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
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to view", metavar="rxxxx")
parser.add_option("-t", "--threshold", action="store", type="float", dest="threshold", help="sets threshold for max intensity of angular average below which hits are automatically sorted to sub-type (default:0)", default=0)
parser.add_option("-L", "--lowerBound", action="store", type="int", dest="lowerBound", help="sets lower bound of pixels for max intensity of angular average below which hits are automatically sorted to sub-type (default:200)", default=200)
parser.add_option("-U", "--upperBound", action="store", type="int", dest="upperBound", help="sets upper bound of pixels for max intensity of angular average below which hits are automatically sorted to sub-type (default:1150)", default=1150)
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value (default:2000)", default=2000)
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
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
ang_avg_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
# RES
#source_dir = "/reg/d/psdm/cxi/cxi74613/res/cleaned_hdf5/"
#ang_avg_dir = "/reg/d/psdm/cxi/cxi/cxi74613/res/cleaned_hdf5/"
# FTC
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#ang_avg_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"

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
		updateTypes[tcounter] = int(input("Threshold "+cDir+"/ (1 for yes, 0 for no)? "))
		foundTypeFiles[tcounter] += G.glob("LCLS*angavg.h5")
		foundFiles += G.glob("LCLS*angavg.h5")
		os.chdir(originaldir)
		tcounter += 1
	print "Found %d types (including type0) with %d sorted files." % (numTypes, len(foundFiles))	

if (options.verbose):
	print "Detector at %lf mm away from the interaction region." % (eDD.get_detector_dist_in_meters(runtag)*1000)


#Global parameters
colmax = options.maxIntens
colmin = 0
storeFlag = 0

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, inangavg , filename, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs, detectorDistance=eDD.get_detector_dist_in_meters(runtag)):
		self.origarr = inarr.copy()
		self.inarr = inarr
		self.filename = filename
		self.inangavg = inangavg
		self.wavelength = meanWaveLengthInAngs
		self.detectorDistance = detectorDistance
		self.HIceQ ={}
		global colmax
		global colmin
		global storeFlag
		colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
		colmin = 0
	
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
		for i,j in eDD.iceHInvAngQ.iteritems():
			self.HIceQ[i] = eDD.get_pix_from_invAngsQ_and_detectorDist(runtag,j,self.detectorDistance, wavelengthInAngs=self.wavelength)

		fig = P.figure(num=None, figsize=(13.5, 6), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		(rTemp,cTemp) = self.inarr.shape
		#Approximate, center
		for i,j in self.HIceQ.iteritems():
			circ = P.Circle((rTemp/2, cTemp/2), radius=j)
			circ.set_fill(False)
			circ.set_edgecolor('k')
			canvas.add_patch(circ)
		
		canvas = fig.add_subplot(122)
		canvas.set_title("Angular Average")
		maxAngAvg = (self.inangavg).max()
		
		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels

		P.plot(self.inangavg)

		pngtag = foundTypes[storeFlag] + '/' + "%s.png" % (self.filename)
		P.savefig(pngtag)
		print "%s saved." % (pngtag)
		P.close()
		#P.show()
	
	def draw_img_for_thresholding(self):
		global colmax
		global colmin
		global storeFlag
		for i,j in eDD.iceHInvAngQ.iteritems():
			self.HIceQ[i] = eDD.get_pix_from_invAngsQ_and_detectorDist(runtag,j,self.detectorDistance, wavelengthInAngs=self.wavelength)

		fig = P.figure(num=None, figsize=(13.5, 6), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		(rTemp,cTemp) = self.inarr.shape
		#Approximate, center
		for i,j in self.HIceQ.iteritems():
			circ = P.Circle((rTemp/2, cTemp/2), radius=j)
			circ.set_fill(False)
			circ.set_edgecolor('k')
			canvas.add_patch(circ)
		
		canvas = fig.add_subplot(122)
		canvas.set_title("Angular Average")
		maxAngAvg = (self.inangavg).max()
		
		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels

		P.plot(self.inangavg)
		
		pngtag = foundTypes[storeFlag] + '/' + "below%dADUs/%s.png" % (options.threshold, self.filename)
		P.savefig(pngtag)
		if (options.verbose):
			print "%s saved." % (pngtag)
		
		P.close()


########################################################
# Loop to display and auto-threshold all H5 files
########################################################

avgRawArr = N.zeros((numTypes,1480,1552))
#avgArr = N.zeros((numTypes+1,1760,1760)) #cxi25410
#avgRadAvg = N.zeros((numTypes+1,1233)) #cxi25410
avgArr = N.zeros((numTypes+1,1764,1764)) #cxi74613
avgRadAvg = N.zeros((numTypes+1,1191)) #cxi74613
typeOccurences = N.zeros(numTypes)
threshRawArr = N.zeros((numTypes,1480,1552))
#threshArr = N.zeros((numTypes,1760,1760)) #cxi25410
#threshRadAvg = N.zeros((numTypes,1233)) #cxi25410
threshArr = N.zeros((numTypes,1764,1764)) #cxi74613
threshRadAvg = N.zeros((numTypes,1191)) #cxi74613
threshOccurences = N.zeros(numTypes)
wavelengths = [[] for i in foundTypeNumbers]

#Loop over each type to threshold and sum all hits
for currentlyExamining in range(numTypes):
	dirName = foundTypes[currentlyExamining]
	storeFlag = currentlyExamining
	fcounter = 0
	tcounter = 0
	numFilesInDir = len(foundTypeFiles[currentlyExamining])
	if (updateTypes[currentlyExamining] == 1):
		t1 = time.time()
		print "Now thresholding H5 files in %s/ ..." % (dirName)
		for fname in foundTypeFiles[currentlyExamining]:
			diffractionName = source_dir + runtag + "/" + re.sub("-angavg",'',fname)
			if os.path.exists(diffractionName):
				f = H.File(diffractionName, 'r')
				d = N.array(f['/data/data']).astype(float)
				draw = N.array(f['/data/rawdata']).astype(float)
				currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
				currDetectorDist=(1.E-3)*f['LCLS']['detectorPosition'][0]
				f.close()
			elif os.path.exists(dirName + '/' + re.sub("-angavg",'',fname)):
				diffractionName = dirName + '/' + re.sub("-angavg",'',fname)
				f = H.File(diffractionName, 'r')
				d = N.array(f['/data/data']).astype(float)
				draw = N.array(f['/data/rawdata']).astype(float)
				currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
				currDetectorDist=(1.E-3)*f['LCLS']['detectorPosition'][0]
				f.close()
			else:
				print "%s does not exist, aborting." % (diffractionName)
				sys.exit(1)
			angAvgName = dirName + '/' + fname
			f = H.File(angAvgName, 'r')
			if (len(N.array(f['data']['data'])) == 2):
				davg = N.array(f['data']['data'][1])
			else:
				davg = N.array(f['data']['data'][0])
			f.close()
			baseName = dirName + '/' + re.sub("-angavg.h5",'',fname)
			if ((davg[options.lowerBound:options.upperBound]).max() < options.threshold):
				subdir = foundTypes[currentlyExamining] + '/' + "below%dADUs/" % (options.threshold)
				if (not os.path.exists(subdir)):
					os.mkdir(subdir)
				if (options.verbose):
					print "Moving %s to " % (fname) + subdir
				os.system("mv " + baseName + "* " + subdir)
				threshArr[currentlyExamining] += d
				threshRawArr[currentlyExamining] += draw
				threshRadAvg[currentlyExamining] += davg
				threshOccurences[currentlyExamining] += 1
				currImg = img_class(d, davg, fname, currWavelengthInAngs, currDetectorDist)
				currImg.draw_img_for_thresholding()
				tcounter += 1
			else:
				avgArr[currentlyExamining] += d
				avgRawArr[currentlyExamining] += draw
				avgRadAvg[currentlyExamining] += davg
				typeOccurences[currentlyExamining] += 1
			wavelengths[currentlyExamining].append(currWavelengthInAngs)
			fcounter += 1
			if (options.verbose and (round(((fcounter*100)%numFilesInDir)/100)==0)):
				print str(fcounter) + " of " + str(numFilesInDir) + " files updated (" + str(fcounter*100/numFilesInDir) + "%)"
		t2 = time.time()
		print "Thresholded %d events in type" % (tcounter) + str(foundTypeNumbers[currentlyExamining])
		if (t2-t1 < 60):
			print "Time taken for thresholding type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(round(t2-t1)) + " s."
		elif (t2-t1 < 3600):
			print "Time taken for thresholding type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(int(t2-t1)/60) + " min " + str((round(t2-t1))%60) + " s."
		else:
			print "Time taken for thresholding type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(int(t2-t1)/3600) + " h " + str(int((t2-t1)%3600)/60) + " min " + str((round(t2-t1))%60) + " s."
		if (options.verbose and len(wavelengths[currentlyExamining]) > 0):
			print "Mean wavelength = " + str(N.mean(wavelengths[currentlyExamining])) + " A."
			print "Relative change in wavelength = " + str(N.std(wavelengths[currentlyExamining])/N.mean(wavelengths[currentlyExamining]))
			print "max-min wavelength = " + str(N.max(wavelengths[currentlyExamining]) - N.min(wavelengths[currentlyExamining])) + " A."


#print "Right-click on colorbar to set maximum scale."
#print "Left-click on colorbar to set minimum scale."
#print "Center-click on colorbar (or press 'r') to reset color scale."
#print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
#print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

storeFlag = 0

#Loop over each type to view and save average
for dirName in foundTypes:
	if (typeOccurences[storeFlag] > 0.):
		avgArr[storeFlag] /= typeOccurences[storeFlag]
		avgRawArr[storeFlag] /= typeOccurences[storeFlag]
		avgRadAvg[storeFlag] /= typeOccurences[storeFlag]		
		typeTag = runtag+'_type'+str(foundTypeNumbers[storeFlag])
		currImg = img_class(avgArr[storeFlag], avgRadAvg[storeFlag], typeTag, meanWaveLengthInAngs=N.mean(wavelengths[storeFlag]))
		currImg.draw_img_for_viewing()
		f = H.File(dirName +'/'+ typeTag + ".h5", "w")
		entry_1 = f.create_group("/data")
		entry_1.create_dataset("diffraction", data=avgArr[storeFlag])
		entry_1.create_dataset("rawdata", data=avgRawArr[storeFlag])
		entry_1.create_dataset("angavg", data=avgRadAvg[storeFlag])
		f.close()
		print "Successfuly updated %s" % (dirName +'/'+ typeTag + ".h5")
	if (threshOccurences[storeFlag] > 0.):
		threshArr[storeFlag] /= threshOccurences[storeFlag]
		threshRawArr[storeFlag] /= threshOccurences[storeFlag]
		threshRadAvg[storeFlag] /= threshOccurences[storeFlag]		
		typeTag = "below%dADUs/"%(options.threshold)+runtag+'_type'+str(foundTypeNumbers[storeFlag])+"-below_%dADUs"%(options.threshold)
		currImg = img_class(threshArr[storeFlag], threshRadAvg[storeFlag], typeTag, meanWaveLengthInAngs=N.mean(wavelengths[storeFlag]))
		currImg.draw_img_for_viewing()
		f = H.File(dirName +'/'+ typeTag + ".h5", "w")
		entry_1 = f.create_group("/data")
		entry_1.create_dataset("diffraction", data=avgArr[storeFlag])
		entry_1.create_dataset("rawdata", data=avgRawArr[storeFlag])
		entry_1.create_dataset("angavg", data=avgRadAvg[storeFlag])	
		f.close()
		print "Successfuly updated %s" % (dirName +'/'+ typeTag + ".h5")
	storeFlag += 1

