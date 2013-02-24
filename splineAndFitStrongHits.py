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
parser.add_option("-m", "--min", action="store", type="float", dest="min_value", help="ignore intensities below this q-value in splined angular average (default: 0.06 A-1)", metavar="MIN_VALUE", default="0.06")
parser.add_option("-x", "--max", action="store", type="float", dest="max_value", help="ignore intensities above this q-value in splined angular average (default: 3.48 A-1)", metavar="MAX_VALUE", default="3.48")
parser.add_option("-d", "--delta", action="store", type="float", dest="delta_value", help="spline intensities with this interval in angular average (default: 0.001 A-1)", metavar="DELTA_VALUE", default="0.001")
parser.add_option("-p", "--peakfit", action="store_true", dest="peakfit", help="applies Gaussian peak fitting algorithm to the angular averages", default=False)
parser.add_option("-s", "--sonemin", action="store", type="float", dest="S1_min", help="lower limit of range used for S1 peak fitting (default: 1.50 A-1)", metavar="MIN_VALUE", default="1.50")
parser.add_option("-t", "--sonemax", action="store", type="float", dest="S1_max", help="upper limit of range used for S1 peak fitting (default: 2.08 A-1)", metavar="MAX_VALUE", default="2.08")
parser.add_option("-u", "--stwomin", action="store", type="float", dest="S2_min", help="lower limit of range used for S2 peak fitting (default: 2.64 A-1)", metavar="MIN_VALUE", default="2.64")
parser.add_option("-w", "--stwomax", action="store", type="float", dest="S2_max", help="upper limit of range used for S2 peak fitting (default: 3.20 A-1)", metavar="MAX_VALUE", default="3.20")
parser.add_option("-o", "--outputdir", action="store", type="string", dest="outputDir", help="output directory (default: output_rxxxx)", metavar="OUTPUT_DIR", default="output")  
parser.add_option("-e", "--exclude", action="store_true", dest="exclude", help="excludes hits from splining/averaging", default=False)
parser.add_option("-f", "--excludefile", action="store", type="string", dest="excludeFile", help="name of text file with hits to exlude (default: below100ADUs)", metavar="EXCLUDE_FILENAME", default="below100ADUs")
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
#source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
#ang_avg_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
# RES
#source_dir = "/reg/data/ana12/cxi/cxi25410/res/cleaned_hdf5/"
#ang_avg_dir = "/reg/data/ana12/cxi/cxi25410/res/cleaned_hdf5/"
# FTC
source_dir = "/reg/data/ana12/cxi/cxi25410/ftc/cleaned_hdf5/"
ang_avg_dir = "/reg/data/ana12/cxi/cxi25410/ftc/cleaned_hdf5/"

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
		if (options.peakfit):
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

if (options.peakfit):
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
		print "Press 'p' to save PNG."
		global colmax
		global colmin
		global p1
		if (options.peakfit):
			p0 = [2.2E8, 1.83, 0.25, 1.7E8, 2.98, 0.2]
			index = N.array([((self.inangavgQ > options.S1_min)[i] and (self.inangavgQ < options.S1_max)[i]) or ((self.inangavgQ > options.S2_min)[i] and (self.inangavgQ < options.S2_max)[i]) for i in range(len(self.inangavgQ))])
			[p1, success] = optimize.leastsq(errfunc, p0[:], args=(self.inangavgQ[index],self.inangavg[index]))
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
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

avgArr = N.zeros((numTypes,1760,1760))
avgRawArr = N.zeros((numTypes,1480,1552))
avgAngAvgQ = N.arange(options.min_value,options.max_value+options.delta_value,options.delta_value)
angAvgLength = int((options.max_value-options.min_value)/options.delta_value)+1
avgAngAvg = N.zeros((numTypes,angAvgLength))
typeOccurences = N.zeros(numTypes)
damaged_events = []
wavelengths = [[] for i in foundTypeNumbers]
attenuations = [[] for i in foundTypeNumbers]
avgIntensities = [[] for i in foundTypeNumbers]
maxIntensities = [[] for i in foundTypeNumbers]
if (options.peakfit):
	fitint1 = [[] for i in foundTypeNumbers]
	fitpos1 = [[] for i in foundTypeNumbers]
	fitfwhm1 = [[] for i in foundTypeNumbers]
	fitint2 = [[] for i in foundTypeNumbers]
	fitpos2 = [[] for i in foundTypeNumbers]
	fitfwhm2 = [[] for i in foundTypeNumbers]

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
				diffractionName = source_dir + runtag + "/" + re.sub("-angavg",'',fname)
				if os.path.exists(diffractionName):
					angAvgName = fname
					f = H.File(angAvgName, 'r')
					davg = N.array(f['data']['data'])
					f.close()
					if (davg[0].max() < options.max_value):
						print "Error in Q-calibration! Qmax = %s, skipping event." % (davg[0].max())
						damaged_events.append(fname)
						continue
					
					f = H.File(diffractionName, 'r')
					d = N.array(f['/data/data']).astype(float)
					draw = N.array(f['/data/rawdata']).astype(float)
					currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
					currAttenuation=f['LCLS']['attenuation'][0]
					f.close()
					f = I.interp1d(davg[0], davg[1])
					davg = f(avgAngAvgQ)
					wavelengths[currentlyExamining].append(currWavelengthInAngs)
					attenuations[currentlyExamining].append(currAttenuation)
					avgIntensities[currentlyExamining].append(draw.mean())
					maxIntensities[currentlyExamining].append(max(davg))
					avgArr[currentlyExamining] += d
					avgRawArr[currentlyExamining] += draw
					avgAngAvg[currentlyExamining] += davg
					
					if (options.peakfit):
						#Gaussian peak fitting
						p0 = [2.2E8, 1.83, 0.25, 1.7E8, 2.98, 0.2]
						index = N.array([((avgAngAvgQ > options.S1_min)[i] and (avgAngAvgQ < options.S1_max)[i]) or ((avgAngAvgQ > options.S2_min)[i] and (avgAngAvgQ < options.S2_max)[i]) for i in range(len(avgAngAvgQ))])
						[p1, success] = optimize.leastsq(errfunc, p0[:], args=(avgAngAvgQ[index],davg[index]))
						if (success):
							fitint1[currentlyExamining].append(p1[0])
							fitpos1[currentlyExamining].append(p1[1])
							fitfwhm1[currentlyExamining].append(p1[2])
							fitint2[currentlyExamining].append(p1[3])
							fitpos2[currentlyExamining].append(p1[4])
							fitfwhm2[currentlyExamining].append(p1[5])
						else:
							print "The Gaussian peak fit failed, skipping %s" % (fname)
					
					typeOccurences[currentlyExamining] += 1
					fcounter += 1
				else:
					print "The diffraction file %s does not exist, ignoring %s" % (diffractionName, fname)
				if (options.verbose and (round(((fcounter*100)%numFilesInDir)/100)==0)):
					if (options.peakfit):
						print str(fcounter) + " of " + str(numFilesInDir) + " files splined and fitted (" + str(fcounter*100/numFilesInDir) + "%)"
					else:
						print str(fcounter) + " of " + str(numFilesInDir) + " files splined (" + str(fcounter*100/numFilesInDir) + "%)"
			
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


if damaged_events:
	damaged_events_name = write_dir + runtag + "_damaged_events.txt"
	print "There are %s damaged events that have been ignored." % (len(damaged_events))
	N.array(damaged_events).tofile(damaged_events_name, sep="\n")
	print "Saved damaged events to %s" % (damaged_events_name)


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
		if(storeFlag > 0):
			typeTag = runtag+'_type'+str(foundTypeNumbers[storeFlag])
		else:
			typeTag = runtag+'_type0'
		
		currImg = img_class(avgArr[storeFlag], avgAngAvg[storeFlag], avgAngAvgQ, typeTag, meanWaveLengthInAngs=N.mean(wavelengths[storeFlag]))
		currImg.draw_img_for_viewing()
		
		if (options.peakfit):
			#Gaussian peak fit statistics
			fitdeltaq = N.array(fitpos2[storeFlag]) - N.array(fitpos1[storeFlag])
			
			fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
			canvas = fig.add_subplot(131)
			canvas.set_title("Histogram")
			P.xlabel("S1 (A-1)")
			P.ylabel("Hist(S1)")
			hist_bins = N.arange(min(fitpos1[storeFlag]), max(fitpos1[storeFlag])+0.003, 0.001) - 0.0005
			pos1_hist, hist_bins = N.histogram(fitpos1[storeFlag], bins=hist_bins)
			pos1_bins = [(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(pos1_hist))]
			P.plot(pos1_bins, pos1_hist)
			P.axvline(p1[1],0,max(pos1_hist),color='r')
			
			canvas = fig.add_subplot(132)
			canvas.set_title("Histogram")
			P.xlabel("deltaQ (A-1)")
			P.ylabel("Hist(deltaQ)")
			hist_bins = N.arange(min(fitdeltaq), max(fitdeltaq)+0.003, 0.001) - 0.0005
			deltaq_hist, hist_bins = N.histogram(fitdeltaq, bins=hist_bins)
			deltaq_bins = [(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(deltaq_hist))]
			P.plot(deltaq_bins, deltaq_hist)
			P.axvline(p1[4]-p1[1],0,max(deltaq_hist),color='r')
			
			canvas = fig.add_subplot(133)
			canvas.set_title("Histogram")
			P.xlabel("S2 (A-1)")
			P.ylabel("Hist(S2)")
			hist_bins = N.arange(min(fitpos2[storeFlag]), max(fitpos2[storeFlag])+0.003, 0.001) - 0.0005
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
		
		f = H.File(dirName +'/'+ typeTag + ".h5", "w")
		entry_1 = f.create_group("/data")
		entry_1.create_dataset("diffraction", data=avgArr[storeFlag])
		entry_1.create_dataset("rawdata", data=avgRawArr[storeFlag])
		entry_1.create_dataset("angavg", data=avgAngAvg[storeFlag])
		entry_1.create_dataset("angavgQ", data=avgAngAvgQ)
		entry_1.create_dataset("attenuation", data=attenuations[storeFlag])
		entry_1.create_dataset("intensityAvg", data=avgIntensities[storeFlag])
		entry_1.create_dataset("intensityMax", data=maxIntensities[storeFlag])
		if (options.peakfit):
			entry_2 = f.create_group("/data/peakFit")
			entry_2.create_dataset("int1", data=fitint1[storeFlag])
			entry_2.create_dataset("int2", data=fitint2[storeFlag])
			entry_2.create_dataset("pos1", data=fitpos1[storeFlag])
			entry_2.create_dataset("pos2", data=fitpos2[storeFlag])
			entry_2.create_dataset("fwhm1", data=fitfwhm1[storeFlag])
			entry_2.create_dataset("fwhm2", data=fitfwhm2[storeFlag])
		f.close()
		print "Successfully updated %s" % (dirName +'/'+ typeTag + ".h5")
		
	storeFlag += 1
