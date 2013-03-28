#!/usr/bin/env python

# Written by J. Sellberg on 2012-01-09
# Has the same functionality as analyzeSums-droplet_dispenser.py but for the new data format
# Modified by J.S. on 2012-01-25 to include the thresholding analysis
# Modified by J.S. on 2013-03-27 to include the peakfitting analysis and switch thresholding through OptionParser

import numpy as N
from numpy import linalg as LA
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
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value", metavar="MAX_VALUE", default=2000)
parser.add_option("-t", "--threshold", action="store", type="float", dest="threshold", help="sets the intensity threshold for the angular averages: 20 (default), 50, or 100", default=20)
parser.add_option("-a", "--xaca", action="store_true", dest="xaca", help="average the xaca files along with the angular averages", default=False)
parser.add_option("-p", "--peakfit", action="store_true", dest="peakfit", help="applies the peakfitting analysis to the angular averages", default=False)
parser.add_option("-S", "--sonemin", action="store", type="float", dest="S1_min", help="lower limit of range used for S1 peak fitting (default: 1.50 A-1)", metavar="MIN_VALUE", default=1.50)
parser.add_option("-T", "--sonemax", action="store", type="float", dest="S1_max", help="upper limit of range used for S1 peak fitting (default: 2.08 A-1)", metavar="MAX_VALUE", default=2.08)
parser.add_option("-U", "--stwomin", action="store", type="float", dest="S2_min", help="lower limit of range used for S2 peak fitting (default: 2.64 A-1)", metavar="MIN_VALUE", default=2.64)
parser.add_option("-W", "--stwomax", action="store", type="float", dest="S2_max", help="upper limit of range used for S2 peak fitting (default: 3.20 A-1)", metavar="MAX_VALUE", default=3.20)
parser.add_option("-e", "--exclude", action="store_true", dest="exclude", help="average excluded hits", default=False)
parser.add_option("-f", "--excludefile", action="store", type="string", dest="excludeFile", help="name of txt file with hits to exlude (default: below100ADUs)", metavar="EXCLUDE_FILENAME", default="below100ADUs")
parser.add_option("-s", "--saveexcluded", action="store_true", dest="saveExcluded", help="flag to save average of excluded hits", default=False)
(options, args) = parser.parse_args()

# adapted to droplet dispenser (water only)
# with r0102
#runs = [[105],[111],[68],[63,64],[110],[73],[77,79],[88],[97],[102]]
#nhits_water = [[2610],[7860],[3900],[381,3100],[3000],[2280],[1225,1480],[1880],[1820],[1675]]
#dhits_water = [[0],[83],[21],[0,4],[22],[6],[11,0],[0],[0],[]]
#temperatures = [290,260,254,249,244,243,240,235,233,232]
# without r0102
#runs = [[105],[111],[68],[63,64],[110],[73],[77,79],[88],[97]]
#nhits_water = [[2610],[7860],[3900],[381,3100],[3000],[2280],[1225,1480],[1880],[1820]]
# with resorted r0063 (current in iceFinderCampaign)
runs = [[105],[111],[68],[63,64],[110],[73],[77,79],[88],[97]]
nhits_water = [[2610],[7860],[3900],[377,3100],[3000],[2280],[1225,1480],[1880],[1820]]
# hits with damged Q-calibration (THESE ARE NOW INCLUDED WITH THE NEW CORRECTIONS)
dhits_water = [[0],[83],[21],[0,4],[22],[6],[11,0],[0],[0]]
# hits with failed peak fitting (2013-03-27, p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2])
failedFits20_water = [[155],[0],[1],[0,0],[0],[0],[0,0],[0],[0]]
# thresholded hits below 50 ADUs
runs = [[105],[111],[68],[63,64],[110],[73],[77,79],[88],[97]]
t50hits_water = [[0],[0],[0],[0,0],[0],[0],[0,0],[0],[0]]
failedFits50_water = [[155],[0],[1],[0,0],[0],[0],[0,0],[0],[0]]
# thresholded hits below 100 ADUs
runs = [[105],[111],[68],[63,64],[110],[73],[77,79],[88],[97]]
t100hits_water = [[0],[0],[0],[0,0],[0],[0],[0,0],[1453],[1020]]
failedFits100_water = [[155],[0],[1],[0,0],[0],[0],[0,0],[0],[0]]
#colors = ['b','g','r','c','m','y','k']
colors = ['r','g','b','c','m','y','k']
#temperatures = [290,260,254,249,244,243,240,235,233]
temperatures = [283,255,249,245,241,240,237,233,231]
#temperatures = [283,257,250,246,241,241,237,233,231] # FINAL temperatures LicThesis
#distances = [0.750000,10.000000,15.681899,20.681899,30.000000,30.699899,40.699899,60.699899,70.649900]
distances = [0.750000,10.011910,15.798895,20.789092,30.006612,30.738471,40.759723,60.710063,70.662700] # FINAL distances
thresholds = [20, 50, 100]

# SCRATCH
#source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
#sorting_dir = "/reg/d/psdm/cxi/cxi25410/scratch/iceFinderCampaign/"
# RES & FTC
source_dir = "/reg/d/psdm/cxi/cxi25410/ftc/cleaned_hdf5/"
#source_dir = "/reg/d/psdm/cxi/cxi25410/res/cleaned_hdf5/"
sorting_dir = "/reg/d/psdm/cxi/cxi25410/res/iceFinderCampaign/"

original_dir = os.getcwd() + '/'
run_tag = "r0063"

thresholdIndex = N.argmax((N.array(thresholds) == options.threshold))
if not (N.array(thresholds) == options.threshold).any():
	print "Chosen threshold of %d ADUs has not been sorted, aborting." % (options.threshold)
	sys.exit(1)

if (options.exclude and options.threshold == 20):
	print "No hits can be excluded at the original threshold (20 ADUs), aborting."
	sys.exit(1)
if (not options.exclude and options.threshold > 20):
	print "Hits have to be excluded at threshold = %d ADUs (>20 ADUs), enabling the --exclude flag." % options.threshold
	options.exclude = True

# statistics for water hits
runhits_water = [[[] for j in runs] for i in thresholds]
len_runs = [len(i) for i in runs]

for i in N.arange(len(runs)):
	for j in N.arange(len(runs[i])):
		if options.peakfit:
			runhits_water[0][i].append(nhits_water[i][j] - failedFits20_water[i][j])
			runhits_water[1][i].append(nhits_water[i][j] - t50hits_water[i][j] - failedFits50_water[i][j])
			runhits_water[2][i].append(nhits_water[i][j] - t50hits_water[i][j] - t100hits_water[i][j] - failedFits100_water[i][j])
		else:
			runhits_water[0][i].append(nhits_water[i][j])
			runhits_water[1][i].append(nhits_water[i][j] - t50hits_water[i][j])
			runhits_water[2][i].append(nhits_water[i][j] - t50hits_water[i][j] - t100hits_water[i][j])

water_pattern = []
water_correlation = []
water_angavg = []
water_angavgQ = []
water_angavgpos1 = []
water_angavgpos2 = []
water_fitint1 = []
water_fitpos1 = []
water_fitfwhm1 = []
water_fitint2 = []
water_fitpos2 = []
water_fitfwhm2 = []
water_fitdeltaq = []
water_pattern_shape = False
water_correlation_shape = False
water_angavg_shape = False
water_angavgQ_shape = N.array([False])
if (options.exclude and options.saveExcluded):
	excludedWater_pattern = []
	excludedWater_correlation = []
	excludedWater_angavg = []
	excludedWater_angavgQ = []
	excludedWater_angavgpos1 = []
	excludedWater_angavgpos2 = []
	excludedWater_fitint1 = []
	excludedWater_fitpos1 = []
	excludedWater_fitfwhm1 = []
	excludedWater_fitint2 = []
	excludedWater_fitpos2 = []
	excludedWater_fitfwhm2 = []
	excludedWater_fitdeltaq = []

#Global parameters
colmax = options.maxIntens
colmin = 0

if options.peakfit:
	#Gaussian peak fitting functions
	fitfunc = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2)) + p[3]*exp(-(x-p[4])**2/(2*p[5]**2))
	errfunc = lambda p, x, y: fitfunc(p, x) - y

#########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
#########################################################
class img_class (object):
	def __init__(self, inarr, inangavg, inangavg_Q, filename, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs, detectorDistance=eDD.get_detector_dist_in_meters(run_tag)):
		self.origarr = inarr.copy()
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.inangavg = inangavg
		self.inangavg_Q = inangavg_Q
		self.wavelength = meanWaveLengthInAngs
		self.detectorDistance = detectorDistance
		global HIceQ
		global HIceQLabel
		global HIceP
		global colmax
		global colmin
		colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
		colmin = self.inarr.min()
	
	def on_keypress_for_viewing(self,event):
		global HIceQ
		global HIceQLabel
		global HIceP
		global colmax
		global colmin
		if event.key == 'e':
			epstag =  "%s.eps" % (self.filename)
			P.savefig(epstag, format='eps')
			print "%s saved." % (epstag)
		if event.key == 'p':
			pngtag =  "%s.png" % (self.filename)
			P.savefig(pngtag)
			print "%s saved." % (pngtag)
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
			P.clim(colmin, colmax)
			P.draw()
	
	def on_click(self, event):
		global HIceQ
		global HIceQLabel
		global HIceP
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
		
	def draw_img_for_viewing_average(self):
		global HIceQ
		global HIceQLabel
		global HIceP
		global colmax
		global colmin
		print "Press 'p' to save PNG, 'e' to save EPS."
		
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename, fontsize='x-large', fontstretch='condensed')
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 300, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		canvas = fig.add_subplot(122)
		canvas.set_title("Angular Average", fontsize='x-large')
		P.xlabel("Q (Angstroms-1)", fontsize='x-large')
		P.ylabel("Average Intensity (ADUs)", fontsize='x-large')
		maxAngAvg = (self.inangavg).max()
		numQLabels = len(HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in HIceQ.iteritems():
			labelPosition = HIceP[i]*maxAngAvg/numQLabels
			P.axvline(j, 0, colmax, color='k', ls='--')
			P.text(HIceQLabel[i], labelPosition, str(i), rotation="45")
		
		P.plot(self.inangavg_Q, self.inangavg, color='r')
		canvas.set_ylim(bottom=0)
		
		P.show()
	
	def draw_img_for_viewing_pattern(self):
		global colmax
		global colmin
		self.inangavg = self.inangavg*(self.inangavg>0)
		#print "Press 'p' to save PNG, 'e' to save EPS."
		
		fig = P.figure(num=None, figsize=(8, 5), dpi=100, facecolor='w', edgecolor='k')
		fig.subplots_adjust(wspace=0.1)
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title("Water T=%sK"%(self.inangavg_Q[1]) + " (%s Hits)"%(self.inangavg_Q[0]), fontsize=22, fontname='sans-serif', fontweight='roman')
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 100, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		#P.show()
		pngtag =  "%s.png" % (self.filename)
		P.savefig(pngtag)
		print "%s saved." % (pngtag)
		P.close()


# create averages for each distance from premade averages for each run
for i in N.arange(len(runs)):
	temp_water_pattern = []
	temp_water_correlation = []
	temp_water_angavg = []
	temp_water_angavgQ = []
	temp_water_fitint1 = []
	temp_water_fitpos1 = []
	temp_water_fitfwhm1 = []
	temp_water_fitint2 = []
	temp_water_fitpos2 = []
	temp_water_fitfwhm2 = []
	if (options.exclude and options.saveExcluded):
		excludedTemp_water_pattern = []
		excludedTemp_water_correlation = []
		excludedTemp_water_angavg = []
		excludedTemp_water_angavgQ = []
		excludedTemp_water_fitint1 = []
		excludedTemp_water_fitpos1 = []
		excludedTemp_water_fitfwhm1 = []
		excludedTemp_water_fitint2 = []
		excludedTemp_water_fitpos2 = []
		excludedTemp_water_fitfwhm2 = []
	for j in N.arange(len(runs[i])):
		if (runs[i][j] < 100):
			run_tag = "r00%s"%(runs[i][j])
		else:
			run_tag = "r0%s"%(runs[i][j])
		run_dir =  'output_' + run_tag + '/'
		#water
		if (nhits_water[i][j] != 0):
		#if ((nhits_water[i][j]-dhits_water[i][j]) != 0):
			if options.exclude:
				typeTag = run_tag + "_type0-" + options.excludeFile
			else:
				typeTag = run_tag + "_type0"
			if os.path.exists(sorting_dir + run_dir + typeTag + '.h5'):
				print 'found: ' + sorting_dir + run_dir + typeTag + '.h5'
				f = H.File(sorting_dir + run_dir + typeTag + '.h5', 'r')
				if (runhits_water[thresholdIndex][i][j] > 0):
					temp_water_pattern.append(N.array(f['data']['diffraction']))
					water_pattern_shape = N.array(f['data']['diffraction']).shape
					if options.xaca:
						temp_water_correlation.append(N.array(f['data']['correlation']))
						water_correlation_shape = N.array(f['data']['correlation']).shape
					temp_water_angavg.append(N.array(f['data']['angavg']))
					water_angavg_shape = N.array(f['data']['angavg']).shape
					temp_water_angavgQ.append(N.array(f['data']['angavgQ']))
					water_angavgQ_shape = N.array(f['data']['angavgQ'])
					if options.peakfit:
						temp_water_fitint1.append(N.array(f['data']['peakFit']['int1']))
						temp_water_fitpos1.append(N.array(f['data']['peakFit']['pos1']))
						temp_water_fitfwhm1.append(N.array(f['data']['peakFit']['fwhm1']))
						temp_water_fitint2.append(N.array(f['data']['peakFit']['int2']))
						temp_water_fitpos2.append(N.array(f['data']['peakFit']['pos2']))
						temp_water_fitfwhm2.append(N.array(f['data']['peakFit']['fwhm2']))
				elif (water_pattern_shape and (water_correlation_shape or not options.xaca) and water_angavg_shape and water_angavgQ_shape.any()):
					temp_water_pattern.append(N.zeros(water_pattern_shape))
					if options.xaca:
						temp_water_correlation.append(N.zeros(water_correlation_shape))
					temp_water_angavg.append(N.zeros(water_angavg_shape))
					temp_water_angavgQ.append(water_angavgQ_shape)
					print "No water hits for r0%s, paddings zeros." % (runs[i][j])
				else:
					print "No water hits for r0%s and shape is unknown, aborting." % (runs[i][j])
					sys.exit(1)
				#excluded hits
				if (options.exclude and options.saveExcluded):
					if ((runhits_water[0][i][j] - runhits_water[thresholdIndex][i][j]) > 0):
						excludedTemp_water_pattern.append(N.array(f['data']['excludedHits']['diffraction']))
						if options.xaca:
							excludedTemp_water_correlation.append(N.array(f['data']['excludedHits']['correlation']))
						excludedTemp_water_angavg.append(N.array(f['data']['excludedHits']['angavg']))
						excludedTemp_water_angavgQ.append(N.array(f['data']['excludedHits']['angavgQ']))
						if options.peakfit:
							excludedTemp_water_fitint1.append(N.array(f['data']['excludedHits']['peakFit']['int1']))
							excludedTemp_water_fitpos1.append(N.array(f['data']['excludedHits']['peakFit']['pos1']))
							excludedTemp_water_fitfwhm1.append(N.array(f['data']['excludedHits']['peakFit']['fwhm1']))
							excludedTemp_water_fitint2.append(N.array(f['data']['excludedHits']['peakFit']['int2']))
							excludedTemp_water_fitpos2.append(N.array(f['data']['excludedHits']['peakFit']['pos2']))
							excludedTemp_water_fitfwhm2.append(N.array(f['data']['excludedHits']['peakFit']['fwhm2']))
					elif (water_pattern_shape and (water_correlation_shape or not options.xaca) and water_angavg_shape and water_angavgQ_shape.any()):
						excludedTemp_water_pattern.append(N.zeros(water_pattern_shape))
						if options.xaca:
							excludedTemp_water_correlation.append(N.zeros(water_correlation_shape))
						excludedTemp_water_angavg.append(N.zeros(water_angavg_shape))
						excludedTemp_water_angavgQ.append(water_angavgQ_shape)
						print "No excluded water hits for r0%s, paddings zeros." % (runs[i][j])
					else:
						print "No excluded water hits for r0%s and shape is unknown, aborting." % (runs[i][j])
						sys.exit(1)
				
				f.close()
			else:
				print 'file missing: ' + sorting_dir + run_dir + typeTag + '.h5'
	
	#plot temp_angavg
	fig = P.figure(num=None, figsize=(8, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(111)
	canvas.set_title("Water, T = %s K"%(temperatures[i]))
	P.xlabel("Q (Angstroms-1)")
	P.ylabel("Average Intensity (ADUs)")
	for j in N.arange(len(runs[i])):
		P.plot(temp_water_angavgQ[j], temp_water_angavg[j], color=colors[j], label="r0%s"%(runs[i][j]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels)
	
	print "output_runs-T%sK-%dADUs-angavg.png saved."%(temperatures[i],options.threshold)
	P.savefig(original_dir + "output_runs-T%sK-%dADUs-angavg.png"%(temperatures[i],options.threshold))
	#P.savefig(original_dir + "output_runs-T%sK-%dADUs-angavg.eps"%(temperatures[i],options.threshold), format='eps')
	#P.show()
	P.close()
	
	sumwater = float(sum(runhits_water[thresholdIndex][i]))
	if (options.exclude and options.saveExcluded):
		excludesumwater = float(sum(runhits_water[0][i]) - sum(runhits_water[thresholdIndex][i]))
	
	for j in N.arange(len(runs[i])):
		nwater = runhits_water[thresholdIndex][i][j]
		if (options.exclude and options.saveExcluded):
			excludenwater = runhits_water[0][i][j] - runhits_water[thresholdIndex][i][j]
		
		if (sumwater > 0):
			temp_water_angavg[j] *= nwater/sumwater
			temp_water_pattern[j] *= nwater/sumwater
			if options.xaca:
				temp_water_correlation[j] *= nwater/sumwater
		if (options.exclude and options.saveExcluded):
			if (excludesumwater > 0):
				excludedTemp_water_angavg[j] *= excludenwater/excludesumwater
				excludedTemp_water_pattern[j] *= excludenwater/excludesumwater
				if options.xaca:
					excludedTemp_water_correlation[j] *= excludenwater/excludesumwater
	
	water_angavg.append(N.array(temp_water_angavg).sum(axis=0))
	water_angavgQ.append(N.array(temp_water_angavgQ).mean(axis=0))
	water_pattern.append(N.array(temp_water_pattern).sum(axis=0))
	water_correlation.append(N.array(temp_water_correlation).sum(axis=0))
	if options.xaca:
		water_correlation.append(N.array(temp_water_correlation).sum(axis=0))
	if (options.exclude and options.saveExcluded):
		excludedWater_angavg.append(N.array(excludedTemp_water_angavg).sum(axis=0))
		excludedWater_angavgQ.append(N.array(excludedTemp_water_angavgQ).mean(axis=0))
		excludedWater_pattern.append(N.array(excludedTemp_water_pattern).sum(axis=0))
		if options.xaca:
			excludedWater_correlation.append(N.array(excludedTemp_water_correlation).sum(axis=0))
	
	currImg = img_class(water_pattern[i], water_angavg[i], [int(sumwater), temperatures[i]], "output_runs-T%sK-%dADUs-pattern"%(temperatures[i],options.threshold))
	currImg.draw_img_for_viewing_pattern()
	
	
	#Gaussian peak fit statistics
	if options.peakfit:
		temp_water_fitint1 = N.array([shot for run in temp_water_fitint1 for shot in run])
		temp_water_fitpos1 = N.array([shot for run in temp_water_fitpos1 for shot in run])
		temp_water_fitfwhm1 = N.array([shot for run in temp_water_fitfwhm1 for shot in run])
		temp_water_fitint2 = N.array([shot for run in temp_water_fitint2 for shot in run])
		temp_water_fitpos2 = N.array([shot for run in temp_water_fitpos2 for shot in run])
		temp_water_fitfwhm2 = N.array([shot for run in temp_water_fitfwhm2 for shot in run])
		#temp_water_fitpos1 = N.array([temp_water_fitpos1[j] for j in range(len(temp_water_fitpos1)) if (temp_water_fitpos1[j] > 1.7 and temp_water_fitpos1[j] < 2.0] and temp_water_fitpos2[j] > 2.8 and temp_water_fitpos2[j] < 3.2))
		#temp_water_fitpos2 = N.array([temp_water_fitpos2[j] for j in range(len(temp_water_fitpos2)) if (temp_water_fitpos1[j] > 1.7 and temp_water_fitpos1[j] < 2.0] and temp_water_fitpos2[j] > 2.8 and temp_water_fitpos2[j] < 3.2))
		temp_water_fitdeltaq = temp_water_fitpos2 - temp_water_fitpos1
		water_fitint1.append(temp_water_fitint1)
		water_fitpos1.append(temp_water_fitpos1)
		water_fitfwhm1.append(temp_water_fitfwhm1)
		water_fitint2.append(temp_water_fitint2)
		water_fitpos2.append(temp_water_fitpos2)
		water_fitfwhm2.append(temp_water_fitfwhm2)
		water_fitdeltaq.append(temp_water_fitdeltaq)
		if (options.exclude and options.saveExcluded):
			excludedTemp_water_fitint1 = N.array([shot for run in excludedTemp_water_fitint1 for shot in run])
			excludedTemp_water_fitpos1 = N.array([shot for run in excludedTemp_water_fitpos1 for shot in run])
			excludedTemp_water_fitfwhm1 = N.array([shot for run in excludedTemp_water_fitfwhm1 for shot in run])
			excludedTemp_water_fitint2 = N.array([shot for run in excludedTemp_water_fitint2 for shot in run])
			excludedTemp_water_fitpos2 = N.array([shot for run in excludedTemp_water_fitpos2 for shot in run])
			excludedTemp_water_fitfwhm2 = N.array([shot for run in excludedTemp_water_fitfwhm2 for shot in run])
			excludedTemp_water_fitdeltaq = excludedTemp_water_fitpos2 - excludedTemp_water_fitpos1
			if (len(excludedTemp_water_fitint1) == 0):
				excludedWater_fitint1.append([0])
				excludedWater_fitpos1.append([0])
				excludedWater_fitfwhm1.append([0])
				excludedWater_fitint2.append([0])
				excludedWater_fitpos2.append([0])
				excludedWater_fitfwhm2.append([0])
				excludedWater_fitdeltaq.append([0])
			else:
				excludedWater_fitint1.append(excludedTemp_water_fitint1)
				excludedWater_fitpos1.append(excludedTemp_water_fitpos1)
				excludedWater_fitfwhm1.append(excludedTemp_water_fitfwhm1)
				excludedWater_fitint2.append(excludedTemp_water_fitint2)
				excludedWater_fitpos2.append(excludedTemp_water_fitpos2)
				excludedWater_fitfwhm2.append(excludedTemp_water_fitfwhm2)
				excludedWater_fitdeltaq.append(excludedTemp_water_fitdeltaq)
		
		p0 = [2.2E8, 1.83, 0.25, 1.7E8, 2.98, 0.2]
		index = N.array([((water_angavgQ[i] > options.S1_min)[j] and (water_angavgQ[i] < options.S1_max)[j]) or ((water_angavgQ[i] > options.S2_min)[j] and (water_angavgQ[i] < options.S2_max)[j]) for j in range(len(water_angavgQ[i]))])
		[p1, success] = optimize.leastsq(errfunc, p0[:], args=(water_angavgQ[i][index],water_angavg[i][index]))
		
		#FIRST PLOT
		fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
		canvas = fig.add_subplot(131)
		canvas.set_title("Mean = %.3f A-1, Median = %.3f A-1, RMS = %.3f A-1"%(temp_water_fitpos1.mean(), N.median(temp_water_fitpos1), temp_water_fitpos1.std()), fontsize='small')
		P.xlabel("S1 (A-1)")
		P.ylabel("Hist(S1)")
		hist_bins = N.arange(1.70, 2.002, 0.001) - 0.0005
		pos1_hist, hist_bins = N.histogram(temp_water_fitpos1, bins=hist_bins)
		pos1_bins = [(hist_bins[j] + hist_bins[j+1])/2 for j in range(len(pos1_hist))]
		P.plot(pos1_bins, pos1_hist)
		P.xlim([1.6, 2.1])
		if (success):
			P.axvline(p1[1],0,max(pos1_hist),color='r')
			water_angavgpos1.append(p1[1])
	
		canvas = fig.add_subplot(132)
		canvas.set_title("Mean = %.3f A-1, Median = %.3f A-1, RMS = %.3f A-1"%(temp_water_fitdeltaq.mean(), N.median(temp_water_fitdeltaq), temp_water_fitdeltaq.std()), fontsize='small')
		P.xlabel("deltaQ (A-1)")
		P.ylabel("Hist(deltaQ)")
		hist_bins = N.arange(0.9, 1.402, 0.001) - 0.0005
		deltaq_hist, hist_bins = N.histogram(temp_water_fitdeltaq, bins=hist_bins)
		deltaq_bins = [(hist_bins[j] + hist_bins[j+1])/2 for j in range(len(deltaq_hist))]
		P.plot(deltaq_bins, deltaq_hist)
		P.xlim([0.9, 1.4])
		if (success):
			P.axvline(p1[4]-p1[1],0,max(deltaq_hist),color='r')
	
		canvas = fig.add_subplot(133)
		canvas.set_title("Mean = %.3f A-1, Median = %.3f A-1, RMS = %.3f A-1"%(temp_water_fitpos2.mean(), N.median(temp_water_fitpos2), temp_water_fitpos2.std()), fontsize='small')
		P.xlabel("S2 (A-1)")
		P.ylabel("Hist(S2)")
		hist_bins = N.arange(2.80, 3.202, 0.001) - 0.0005
		pos2_hist, hist_bins = N.histogram(temp_water_fitpos2, bins=hist_bins)
		pos2_bins = [(hist_bins[j] + hist_bins[j+1])/2 for j in range(len(pos2_hist))]
		P.plot(pos2_bins, pos2_hist)
		P.xlim([2.7, 3.2])
		if (success):
			P.axvline(p1[4],0,max(pos2_hist),color='r')
			water_angavgpos2.append(p1[4])
	
		print "output_runs-T%sK-%dADUs-peakfit_hist.png saved."%(temperatures[i],options.threshold)
		P.savefig(original_dir + "output_runs-T%sK-%dADUs-peakfit_hist.png"%(temperatures[i],options.threshold))
		#P.savefig(original_dir + "output_runs-T%sK-%dADUs-peakfit_hist.eps"%(temperatures[i],options.threshold), format='eps')
		#P.show()
		P.close()

		#SECOND PLOT
		fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
		canvas = fig.add_subplot(131)
		canvas.set_title("S1 Median = %.3f A-1, deltaQ Median = %.3f A-1"%(N.median(temp_water_fitpos1), N.median(temp_water_fitdeltaq)), fontsize='small')
		P.xlabel("S1 (A-1)")
		P.ylabel("deltaQ (A-1)")
		P.plot(temp_water_fitpos1, temp_water_fitdeltaq, 'r.', N.median(temp_water_fitpos1), N.median(temp_water_fitdeltaq), 'kx')
		
		canvas = fig.add_subplot(132)
		canvas.set_title("S1 Median = %.3f A-1, int Median = %d ADU/srad"%(N.median(temp_water_fitpos1), N.median(temp_water_fitint1)), fontsize='small')
		P.xlabel("S1 (A-1)")
		P.ylabel("S1 peak intensity (ADU/srad)")
		P.plot(temp_water_fitpos1, temp_water_fitint1, 'r.', N.median(temp_water_fitpos1), N.median(temp_water_fitint1), 'kx')
		
		canvas = fig.add_subplot(133)
		canvas.set_title("S1 Median = %.3f A-1, FWHM Median = %.3f A-1"%(N.median(temp_water_fitpos1), N.median(temp_water_fitfwhm1)), fontsize='small')
		P.xlabel("S1 (A-1)")
		P.ylabel("S1 FWHM (A-1)")
		P.plot(temp_water_fitpos1, temp_water_fitfwhm1, 'r.', N.median(temp_water_fitpos1), N.median(temp_water_fitfwhm1), 'kx')
		
		print "output_runs-T%sK-%dADUs-peakfit_corr-s1.png saved."%(temperatures[i],options.threshold)
		P.savefig(original_dir + "output_runs-T%sK-%dADUs-peakfit_corr-s1.png"%(temperatures[i],options.threshold))
		#P.savefig(original_dir + "output_runs-T%sK-peakfit_corr-s1.eps"%(temperatures[i],options.threshold), format='eps')
		#P.show()
		P.close()

		#THIRD PLOT
		fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
		canvas = fig.add_subplot(131)
		canvas.set_title("S2 Median = %.3f A-1, deltaQ Median = %.3f A-1"%(N.median(temp_water_fitpos2), N.median(temp_water_fitdeltaq)), fontsize='small')
		P.xlabel("S2 (A-1)")
		P.ylabel("deltaQ (A-1)")
		P.plot(temp_water_fitpos2, temp_water_fitdeltaq, 'r.', N.median(temp_water_fitpos2), N.median(temp_water_fitdeltaq), 'kx')
		
		canvas = fig.add_subplot(132)
		canvas.set_title("S2 Median = %.3f A-1, int Median = %.0f ADU/srad"%(N.median(temp_water_fitpos2), N.median(temp_water_fitint2)), fontsize='small')
		P.xlabel("S2 (A-1)")
		P.ylabel("S2 peak intensity (ADU/srad)")
		P.plot(temp_water_fitpos2, temp_water_fitint2, 'r.', N.median(temp_water_fitpos2), N.median(temp_water_fitint2), 'kx')
		
		canvas = fig.add_subplot(133)
		canvas.set_title("S2 Median = %.3f A-1, S1 Median = %.3f A-1"%(N.median(temp_water_fitpos2), N.median(temp_water_fitpos1)), fontsize='small')
		P.xlabel("S2 (A-1)")
		P.ylabel("S1 (A-1)")
		P.plot(temp_water_fitpos2, temp_water_fitpos1, 'r.', N.median(temp_water_fitpos2), N.median(temp_water_fitpos1), 'kx')
		
		print "output_runs-T%sK-%dADUs-peakfit_corr-s2.png saved."%(temperatures[i],options.threshold)
		P.savefig(original_dir + "output_runs-T%sK-%dADUs-peakfit_corr-s2.png"%(temperatures[i],options.threshold))
		#P.savefig(original_dir + "output_runs-T%sK-%dADUs-peakfit_corr-s2.eps"%(temperatures[i],options.threshold), format='eps')
		#P.show()
		P.close()


#plot angavg
fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(121)
#canvas.set_title("Water", fontsize='x-large')
#P.xlabel("Q (Angstroms-1)", fontsize='x-large')
#P.ylabel("Average Intensity (ADUs)", fontsize='x-large')
canvas.set_title("Water")
P.xlabel("Q (Angstroms-1)")
P.ylabel("Average Intensity (ADUs/srad)")
for i in N.arange(len(runs)):
	P.plot(water_angavgQ[i], water_angavg[i], color=colors[i%len(colors)], label="T = %s K"%(temperatures[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left', bbox_to_anchor = (1, 1))

print "output_runs-droplet_dispenser-all_T-angavg_Q_%dADUs.png saved." % options.threshold
P.savefig(original_dir + "output_runs-droplet_dispenser-all_T-angavg_Q_%dADUs.png" % options.threshold)
#P.savefig(original_dir + "output_runs-droplet_dispenser-all_T-angavg_Q_%dADUs.eps" % options.threshold, format='eps')
#P.show()
P.close()

#save to file
hdf5tag = "output_runs-droplet_dispenser-all_T+Q_%dADUs.h5" % options.threshold
f = H.File(original_dir + hdf5tag, 'w')
entry_1 = f.create_group("/data")
for i in N.arange(len(runs)):
	entry_2 = f.create_group("/data/%.1fmm"%(distances[i]))
	entry_3 = f.create_group("/data/%.1fmm/water"%(distances[i]))
	entry_3.create_dataset("diffraction", data=water_pattern[i])
	if options.xaca:
		entry_3.create_dataset("correlation", data=water_correlation[i])
	entry_3.create_dataset("angavg", data=water_angavg[i])
	entry_3.create_dataset("angavg_Q", data=water_angavgQ[i])
	if options.peakfit:
		entry_4 = f.create_group("/data/%.1fmm/water/peakFit"%(distances[i]))
		entry_4.create_dataset("int1", data=water_fitint1[i])
		entry_4.create_dataset("int2", data=water_fitint2[i])
		entry_4.create_dataset("pos1", data=water_fitpos1[i])
		entry_4.create_dataset("pos2", data=water_fitpos2[i])
		entry_4.create_dataset("fwhm1", data=water_fitfwhm1[i])
		entry_4.create_dataset("fwhm2", data=water_fitfwhm2[i])
		entry_4.create_dataset("deltaQ", data=water_fitdeltaq[i])
	if (options.exclude and options.saveExcluded):
		entry_5 = f.create_group("/data/%.1fmm/water/excludedHits"%(distances[i]))
		entry_5.create_dataset("diffraction", data=excludedWater_pattern[i])
		if options.xaca:
			entry_5.create_dataset("correlation", data=excludedWater_correlation[i])
		entry_5.create_dataset("angavg", data=excludedWater_angavg[i])
		entry_5.create_dataset("angavg_Q", data=excludedWater_angavgQ[i])
		if options.peakfit:
			entry_6 = f.create_group("/data/%.1fmm/water/excludedHits/peakFit"%(distances[i]))
			entry_6.create_dataset("int1", data=excludedWater_fitint1[i])
			entry_6.create_dataset("int2", data=excludedWater_fitint2[i])
			entry_6.create_dataset("pos1", data=excludedWater_fitpos1[i])
			entry_6.create_dataset("pos2", data=excludedWater_fitpos2[i])
			entry_6.create_dataset("fwhm1", data=excludedWater_fitfwhm1[i])
			entry_6.create_dataset("fwhm2", data=excludedWater_fitfwhm2[i])
			entry_6.create_dataset("deltaQ", data=excludedWater_fitdeltaq[i])
f.close()
print "Successfully updated %s" % hdf5tag
