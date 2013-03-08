#!/usr/bin/env python
import numpy as N
from numpy import linalg as LA
import h5py as H
import glob as G
import pylab as P
import sys, os, re, shutil, subprocess, time
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to view", metavar="rxxxx")
parser.add_option("-i", "--inspectOnly", action="store_true", dest="inspectOnly", help="inspect output directory", default=False)
parser.add_option("-o", "--outputDir", action="store", type="string", dest="outputDir", help="output directory will be appended by run number (default: output_rxxxx); separate types will be stored in output_rxxxx/anomaly/type[1-3]", default="output")
parser.add_option("-a", "--auto", action="store_true", dest="auto", help="ignores any manual sorting (of hits not automatically classified as type2)",default=False)
parser.add_option("-W", "--waterAveraging", action="store_true", dest="averageWaterTypes", help="average pattern and angavg of water types", default=False)
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value (default:2000)", default=2000)
parser.add_option("-S", "--sortTypes", action="store", type="int", dest="sortTypes", help="default:0. -1(descending total intens), 0(peakyness), 1(ascending total intens).", default=0)
parser.add_option("-T", "--thresholdType2", action="store", type="float", dest="thresholdType2", help="sets threshold for max intensity of angular average below which hits are automatically sorted to type2 (default:0)", default=0)
parser.add_option("-L", "--lowerBoundType2", action="store", type="int", dest="lowerBoundType2", help="sets lower bound of pixels for max intensity of angular average below which hits are automatically sorted to type2 (default:200)", default=200)
parser.add_option("-U", "--upperBoundType2", action="store", type="int", dest="upperBoundType2", help="sets upper bound of pixels for max intensity of angular average below which hits are automatically sorted to type2 (default:1150)", default=1150)
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
write_anomaly_dir = write_dir
if(not os.path.exists(write_anomaly_dir)):
	os.mkdir(write_anomaly_dir)
numTypes = 9
write_anomaly_dir_types = [write_dir]
for i in range(numTypes):
	write_anomaly_dir_types.append(write_anomaly_dir+"type"+str(i+1)+"/") 

#Change into data directory to extract *angavg.h5 files
arr = []
originaldir=os.getcwd()
os.chdir(write_dir)
files = G.glob("LCLS*angavg.h5")
print "reading ang_avgs.."
for i in files:
	f = H.File(i, 'r')
	if (len(N.array(f['data']['data'])) == 2):
		arr.append(N.array(f['data']['data'][1]))
	else:
		arr.append(N.array(f['data']['data'][0]))
	f.close()
os.chdir(originaldir)
masterArr = N.array(arr)
numData = len(masterArr)
if (numData == 0):
	print "no ang_avgs in %s aborting." % write_dir
	sys.exit(1)
angAvgLen = len(masterArr[0])

#Normalize to water ring
normed_arr = N.zeros((numData, angAvgLen))
sorted_arr = N.zeros((numData, angAvgLen))
sortedFileNames = []
unnormed_arr = masterArr.copy()
for i in range(numData):
	temp = masterArr[i]
	max_temp = N.max(temp[530:560])
	min_temp = N.min(temp[50:1153])
	normed_arr[i] = (temp - min_temp) / (max_temp - min_temp) 

#Sorting routines
if(options.sortTypes==-1):
	print "sorting by total intensities in descending order.."
	scoreKeeper = [N.sum(N.abs(i)) for i in unnormed_arr]
	ordering = (N.argsort(scoreKeeper))[-1::-1]
	sorted_arr = normed_arr[ordering]
	sortedFileNames = N.array(files)[ordering]
elif (options.sortTypes==1):
	print "sorting by total intensities in ascending order.."
	scoreKeeper = [N.sum(N.abs(i)) for i in unnormed_arr]
	ordering = N.argsort(scoreKeeper)
	sorted_arr = normed_arr[ordering]
	sortedFileNames = N.array(files)[ordering]
elif (options.sortTypes==0):
	print "sorting by maximum of median filtered ang_avgs.."
	filterLen = 5
	medianFiltered_arr = N.zeros((numData, angAvgLen-filterLen))
	for i in range(numData):
		for j in range(len(normed_arr[i])-filterLen):
			medianFiltered_arr[i][j] = N.median(normed_arr[i][j:j+filterLen])
	scoreKeeper = [N.max(N.abs(i[201:1001]-i[200:1000])) for i in medianFiltered_arr]
	ordering = N.argsort(scoreKeeper)
	sorted_arr = normed_arr[ordering]
	sortedFileNames = N.array(files)[ordering]

#Global parameters
colmax = options.maxIntens
colmin = 0
storeFlag = 0
########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, inangavg , filename, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs, detectorDistance=eDD.get_detector_dist_in_meters(runtag)):
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.inangavg = inangavg
		self.wavelength = meanWaveLengthInAngs
		self.detectorDistance = detectorDistance
		self.HIceQ ={}
		global colmax
		global colmin
		global storeFlag
		self.tag = 0

	def on_keypress_for_tagging(self, event):
		global colmax
		global colmin
		global storeFlag
		if event.key in [str(i) for i in range(numTypes+1)]:
			storeFlag = int(event.key)
			
			if(options.inspectOnly):
				print "Inspection only mode."
			else:
				if(not os.path.exists(write_anomaly_dir_types[storeFlag])):
					os.mkdir(write_anomaly_dir_types[storeFlag])
				pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
				if(self.tag != 0):
					#delete previous assignment
					pngtag = write_anomaly_dir_types[self.tag] + "%s.png" % (self.filename)
					if os.path.isfile(pngtag):
						os.remove(pngtag)
						print "%s removed!" % (pngtag)
					else:
						print "No action taken."
					#Save new assignment if it's store flag not type 0
					if (storeFlag !=0):
							pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
							P.savefig(pngtag)
							print "%s saved." % (pngtag)
							self.tag = storeFlag
				else:
					P.savefig(pngtag)
					print "%s saved." % (pngtag)
					self.tag = storeFlag
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
			P.clim(colmin, colmax)
			P.draw()

	def on_keypress_for_viewing(self, event):
		global colmax
		global colmin
		global storeFlag
		if event.key == 'p':
			pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
			if(options.inspectOnly):
				print "Inspection only mode."
			else:
				P.savefig(pngtag)
				print "%s saved." % (pngtag)
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
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
		if(not options.inspectOnly):
			print "Press 'p' to save PNG."
		global colmax
		global colmin
		P.close('all')
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 2000, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		canvas = fig.add_subplot(122)
		canvas.set_title("angular average")
		maxAngAvg = (self.inangavg).max()
		for i,j in eDD.iceHInvAngQ.iteritems():
			self.HIceQ[i] = eDD.get_pix_from_invAngsQ_and_detectorDist(runtag,j,self.detectorDistance, wavelengthInAngs=self.wavelength)

		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels

		P.plot(self.inangavg)
		P.show()

	def draw_img_for_tagging(self):
		if(not options.inspectOnly and not options.auto):
			print "Press 1-"+ str(numTypes)+ " to save png (overwrites old PNGs); Press 0 to undo (deletes png if wrongly saved)."
		global colmax
		global colmin
		global storeFlag
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_tagging)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 2000, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		canvas = fig.add_subplot(122)
		canvas.set_title("angular average")
		maxAngAvg = (self.inangavg).max()
		for i,j in eDD.iceHInvAngQ.iteritems():
			self.HIceQ[i] = eDD.get_pix_from_invAngsQ_and_detectorDist(runtag,j,self.detectorDistance, wavelengthInAngs=self.wavelength)

		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels
		
		P.plot(self.inangavg)
		#Check for max intensity of angular average below threshold
		if (not options.inspectOnly and (self.inangavg[options.lowerBoundType2:options.upperBoundType2]).max() < options.thresholdType2):
			print "%s automatically characterized as type2." % (self.filename)
			storeFlag = 2
			if(not os.path.exists(write_anomaly_dir_types[storeFlag])):
				os.mkdir(write_anomaly_dir_types[storeFlag])
			self.axes.set_clim(0, 2000)
			pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
			P.savefig(pngtag)
			print "%s saved." % (pngtag)
			self.tag = storeFlag
			P.close()
		elif (options.auto):
			P.close()
		else:
			P.show()
	
	def draw_spectrum(self):
		print "Press 'p' to save PNG."
		global colmax
		global colmin
		localColMax=self.inarr.max()
		localColMin=self.inarr.min()
		aspectratio = 1.5*(self.inarr.shape[1])/(float(self.inarr.shape[0]))
		fig = P.figure(num=None, figsize=(13, 10), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing) 
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_axes([0.05,0.05,0.6,0.9], xlabel="q", ylabel="normalized angular average (will prompt to examine data larger than cutoff)")
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', aspect=aspectratio, vmax = localColMax, vmin = localColMin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim() 
		canvas2 = fig.add_axes([0.7,0.05,0.25,0.9], xlabel="log(sorting score)", ylabel="data")
		canvas2.set_ylim([0,numData])
		canvas2.plot(N.log(N.array(scoreKeeper)[ordering]),range(numData))
		P.show()  	


print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

currImg = img_class(sorted_arr, None, runtag+"_spectrum_sort%s"%(options.sortTypes))
currImg.draw_spectrum()
cutoff = int(input("ice/water cutoff? "))

########################################################
# Loop to display all non-anomalous H5 files. 
########################################################

avgArr = N.zeros((numTypes+1,1760,1760)) #cxi25410
avgRadAvg = N.zeros((numTypes+1,1233)) #cxi25410
#avgArr = N.zeros((numTypes+1,1764,1764)) #cxi74613
#avgRadAvg = N.zeros((numTypes+1,1191)) #cxi74613
typeOccurences = N.zeros(numTypes+1)
waveLengths={}
for i in range(numTypes):
	waveLengths[i] = []

if(options.averageWaterTypes and len(sortedFileNames[:cutoff]) > 0):
	print "averaging water-only types.."
	for fname in sortedFileNames[:cutoff]:
		diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
		f = H.File(diffractionName, 'r')
		d = N.array(f['/data/data'])
		waveLengths[0].append(f['LCLS']['photon_wavelength_A'][0])
		f.close()
		avgArr[0] += d
		angAvgName = write_dir + fname
		f = H.File(angAvgName, 'r')
		if (len(N.array(f['data']['data'])) == 2):
			davg = N.array(f['data']['data'][1])
		else:
			davg = N.array(f['data']['data'][0])
		f.close()
		avgRadAvg[0] += davg
		typeOccurences[0] += 1.

	avgArr[0] /= typeOccurences[0]
	avgRadAvg[0] /= typeOccurences[0]

	storeFlag = 0
	currImg = img_class(avgArr[0], avgRadAvg[0], runtag+"_type0",meanWaveLengthInAngs=N.mean(waveLengths[0]))
	currImg.draw_img_for_viewing()

########################################################
# Loop to display all H5 files with ice anomalies. 
########################################################
print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

anomalies = sortedFileNames[cutoff:]

waveLengths={}
rangeNumTypes = range(1,numTypes+1)
for i in range(numTypes):
	waveLengths[i] = []

#Tag anomalies
fcounter = 0
t1 = time.time()
for fname in anomalies:
	storeFlag=0
	diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
	f = H.File(diffractionName, 'r')
	d = N.array(f['/data/data'])
	currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
	currDetectorDist=(1.E-3)*f['LCLS']['detectorPosition'][0] 
	f.close()
	correlationName = write_dir+re.sub("-angavg","-xaca",fname)
	angAvgName = write_dir + fname
	f = H.File(angAvgName, 'r')
	if (len(N.array(f['data']['data'])) == 2):
		davg = N.array(f['data']['data'][1])
	else:
		davg = N.array(f['data']['data'][0])
	f.close()
	print "%d: wavelength = %lf A, detectorPos = %lf mm"%(cutoff+fcounter,currWavelengthInAngs,1000*currDetectorDist)
	currImg = img_class(d, davg, fname, meanWaveLengthInAngs=currWavelengthInAngs, detectorDistance=currDetectorDist)
	currImg.draw_img_for_tagging()
	
	if((storeFlag in rangeNumTypes) and not options.inspectOnly):
		waveLengths[storeFlag].append(currWavelengthInAngs)
		avgArr[storeFlag] += d
		avgRadAvg[storeFlag] += davg
		typeOccurences[storeFlag] += 1
		if(not os.path.exists(write_anomaly_dir_types[storeFlag])):
			os.mkdir(write_anomaly_dir_types[storeFlag])

		print "mv " + angAvgName + " " + write_anomaly_dir_types[storeFlag]
		os.system("mv " + angAvgName + " " + write_anomaly_dir_types[storeFlag])
		os.system("cp " + diffractionName + " " + write_anomaly_dir_types[storeFlag])
		if(os.path.exists(correlationName)):
			print "mv " + correlationName + " " + write_anomaly_dir_types[storeFlag]
			os.system("mv " + correlationName + " " + write_anomaly_dir_types[storeFlag])
	
	fcounter += 1

t2 = time.time()
if (t2-t1 < 60):
	print "Time taken for averaging type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(round(t2-t1)) + " s."
elif (t2-t1 < 3600):
	print "Time taken for averaging type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(int(t2-t1)/60) + " min " + str((round(t2-t1))%60) + " s."
else:
	print "Time taken for averaging type" + str(foundTypeNumbers[currentlyExamining]) + " = " + str(int(t2-t1)/3600) + " h " + str(int((t2-t1)%3600)/60) + " min " + str((round(t2-t1))%60) + " s."

#View the averages. Tagging disabled.
for i in range(numTypes):
	if (typeOccurences[i] > 0.):
		print "averaging new hits in type%s.."%i
		storeFlag=i
		avgArr[i] /= typeOccurences[i]
		avgRadAvg[i] /= typeOccurences[i]
		typeTag = runtag+'_'+'type'+str(i)
		currImg = img_class(avgArr[i], avgRadAvg[i], typeTag, meanWaveLengthInAngs=N.mean(waveLengths[i]))
		currImg.draw_img_for_viewing()
		(sx,sy) = avgArr[i].shape
		if (not options.inspectOnly):
			f = H.File(write_anomaly_dir_types[i] + typeTag + ".h5", "w")
			entry_1 = f.create_group("/data")
			entry_1.create_dataset("diffraction", data=avgArr[i])
			entry_1.create_dataset("angavg", data=avgRadAvg[i])	
			f.close()
