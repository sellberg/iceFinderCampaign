#!/usr/bin/env python

# Written by J. Sellberg on 2013-02-26
# Compares the output HDF5 files of averageRuns-aerojet.py from various thresholds/sorting algorithms

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
parser.add_option("-a", "--xaca", action="store_true", dest="xaca", help="averages the xaca files along with the angular averages", default=False)
parser.add_option("-p", "--peakfit", action="store_true", dest="peakfit", help="applies the peakfitting analysis to the angular averages", default=False)
parser.add_option("-S", "--sonemin", action="store", type="float", dest="S1_min", help="lower limit of range used for S1 peak fitting (default: 1.50 A-1)", metavar="MIN_VALUE", default=1.50)
parser.add_option("-T", "--sonemax", action="store", type="float", dest="S1_max", help="upper limit of range used for S1 peak fitting (default: 2.08 A-1)", metavar="MAX_VALUE", default=2.08)
parser.add_option("-U", "--stwomin", action="store", type="float", dest="S2_min", help="lower limit of range used for S2 peak fitting (default: 2.64 A-1)", metavar="MIN_VALUE", default=2.64)
parser.add_option("-W", "--stwomax", action="store", type="float", dest="S2_max", help="upper limit of range used for S2 peak fitting (default: 3.20 A-1)", metavar="MAX_VALUE", default=3.20)
(options, args) = parser.parse_args()

files = ["output_runs-aerojet-all_T+Q_20ADUs-all.h5", "output_runs-aerojet-all_T+Q_20ADUs.h5"]
tags = ["all", "without failedFits"]
# resorted 2012-11-18 (r0167, r0171 r0172)
runs = [[114],[121],[118],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
nhits_water = [[12077],[3308],[10143],[15702],[3508,789,4986],[320,203,104,2280,1159,919,207],[177,127,143,191,206],[30,75]]
nhits_ice = [[1],[7],[7],[11],[2,2,8],[51,42,22,266,127,499,140],[786,625,433,749,1705],[2865,630]]
nevents = [[181089],[130074],[208489],[298050],[158756,62320,316507],[114978,194250,16146,309784,110923,211603,56893],[218817,166847,82800,198154,223197],[60809,213783]]
# thresholded hits below 50 ADUs
runs = [[114],[121],[118],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
t50hits_water = [[4028],[536],[1305],[2032],[858,120,788],[78,64,14,545,305,498,102],[120,63,90,95,115],[6,36]]
t50hits_ice = [[0],[4],[3],[8],[0,0,5],[4,1,5,3,5,67,18],[312,285,27,34,830],[2402,167]]
# thresholded hits below 100 ADUs
runs = [[114],[121],[118],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
t100hits_water = [[4493],[512],[1077],[2172],[677,158,573],[116,25,17,299,151,396,102],[66,54,63,101,86],[3,36]]
t100hits_ice = [[0],[1],[3],[0],[0,1,1],[8,6,1,25,25,127,35],[230,189,170,249,417],[319,210]]
# hits with failed peak fitting (2013-02-26, p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2])
runs = [[114],[121],[118],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
failedFits_water = [[832],[234],[1312],[587],[337,19,293],[12,39,3,241,126,90,20],[11,5,2,8,24],[1,6]]
#colors = ['b','g','r','c','m','y','k']
colors = ['r','g','b','c','m','y','k']
#temperatures = [264,235,234,227,224,221,220,219]
temperatures = [270,253,252,232,227,223,221,220]
#temperatures = [293,256,250,232,227,223,221,220] # FINAL temperatures from LicThesis
#distances = [1.085200,11.585200,11.985100,21.585200,30.540800,40.540800,45.540800,50.540800]
distances = [0.668972,11.548956,12.349269,21.407288,30.603057,40.403027,45.595702,50.619378] # FINAL distances

# reference values from supplementary material Huang_resubmitted120910_Supplementary_Notes.pdf
ref_temp = [250, 232, 227, 223, 221, 220]
ref_s1 = N.array([1.972, 1.901, 1.869, 1.851, 1.868, 1.867])
ref_s2 = N.array([3.065, 3.049, 3.060, 3.060, 3.079, 3.050])
ref_dq = ref_s2 - ref_s1

# hexagonal ice peaks
HIceQ = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.222, '112':3.272, '201':3.324}
HIcePos = {'100':10.5, '002':9.5, '101':8.5, '102':6.2, '110':6.2, '103':6.2, '200':7., '112':6., '201':5.}
HIceQLabel = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.324, '112':3.374, '201':3.426}
# for imaging class
HIceP = {'100':9.5, '002':8.5, '101':7.5, '102':6.2, '110':5.7, '103':5.2, '200':7., '112':6., '201':5.}

# SCRATCH
#source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
#sorting_dir = "/reg/d/psdm/cxi/cxi25410/scratch/iceFinderCampaign/"
# RES & FTC
source_dir = "/reg/data/ana12/cxi/cxi25410/ftc/cleaned_hdf5/"
#source_dir = "/reg/data/ana12/cxi/cxi25410/res/cleaned_hdf5/"
sorting_dir = "/reg/data/ana12/cxi/cxi25410/res/iceFinderCampaign/"

original_dir = os.getcwd() + '/'
run_tag = "r0144"

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
		self.inarr = N.array(inarr).mean(axis=0)
		self.inarr = self.inarr*(self.inarr>0)
		self.filename = filename
		self.inangavg = inangavg
		self.inangavgQ = inangavg_Q
		self.wavelength = meanWaveLengthInAngs
		self.detectorDistance = detectorDistance
		global HIceQ
		global HIceQLabel
		global HIceP
		global colmax
		global colmin
		colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
		colmin = 0
	
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
	
	def draw_img_for_viewing_water(self):
		#print "Press 'p' to save PNG."
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
		for i in range(len(self.inangavg)):
			if options.peakfit:
				p0 = [2.2E8, 1.83, 0.25, 1.7E8, 2.98, 0.2]
				index = N.array([((self.inangavgQ[i] > options.S1_min)[j] and (self.inangavgQ[i] < options.S1_max)[j]) or ((self.inangavgQ[i] > options.S2_min)[j] and (self.inangavgQ[i] < options.S2_max)[j]) for j in range(len(self.inangavgQ[i]))])
				[p1, success] = optimize.leastsq(errfunc, p0[:], args=(self.inangavgQ[i][index],self.inangavg[i][index]))
				if success:
					P.plot(self.inangavgQ[i], self.inangavg[i], color="%s"%colors[i], label="S1 = %.3f A-1, S2 = %.3f A-1" % (p1[1], p1[4]))
					P.plot(self.inangavgQ[i], fitfunc(p1, self.inangavgQ[i]), "%s:"%colors[i])
				else:
					P.plot(self.inangavgQ[i], self.inangavg[i], "%s-"%colors[i])
			else:
				P.plot(self.inangavgQ[i], self.inangavg[i], "%s-"%colors[i])
		
		handles, labels = canvas.get_legend_handles_labels()
		leg = canvas.legend(handles, labels, loc='lower right')
		setp(leg.get_texts(), fontsize='medium')
		canvas.set_title("Angular Average")
		P.xlabel("Q (A-1)")
		P.ylabel("I(Q) (ADU/srad)")
		pngtag = original_dir + "peakfit-gdvn_%s.png" % (self.filename)
		P.savefig(pngtag)
		print "%s saved." % (pngtag)
		P.close()
		#P.show()
	
	def draw_img_for_viewing_ice(self):
		#print "Press 'p' to save PNG."
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
		P.ylabel("I(Q) (ADU/srad)")
		pngtag = original_dir + "peakfit-gdvn_%s.png" % (self.filename)
		P.savefig(pngtag)
		print "%s saved." % (pngtag)
		P.close()
		#P.show()


water_pattern = [[] for i in distances]
if options.xaca:
	water_correlation = [[] for i in distances]
water_angavg = [[] for i in distances]
water_angavgQ = [[] for i in distances]
if options.peakfit:
	water_fitint1 = [[] for i in distances]
	water_fitpos1 = [[] for i in distances]
	water_fitfwhm1 = [[] for i in distances]
	water_fitint2 = [[] for i in distances]
	water_fitpos2 = [[] for i in distances]
	water_fitfwhm2 = [[] for i in distances]
	water_fitdeltaq = [[] for i in distances]
ice_pattern = [[] for i in distances]
if options.xaca:
	ice_correlation = [[] for i in distances]
ice_angavg = [[] for i in distances]
ice_angavgQ = [[] for i in distances]

# read and plot averages for each distance from premade files
for i in range(len(files)):
	if os.path.exists(original_dir + files[i]):
		print 'found: ' + original_dir + files[i]
		f = H.File(original_dir + files[i], 'r')
		for j in N.arange(len(distances)):
			water_pattern[j].append(N.array(f['data']['%.1fmm'%distances[j]]['water']['diffraction']))
			if options.xaca:
				water_correlation[j].append(N.array(f['data']['%.1fmm'%distances[j]]['water']['correlation']))
			water_angavg[j].append(N.array(f['data']['%.1fmm'%distances[j]]['water']['angavg']))
			water_angavgQ[j].append(N.array(f['data']['%.1fmm'%distances[j]]['water']['angavg_Q']))
			ice_pattern[j].append(N.array(f['data']['%.1fmm'%distances[j]]['ice']['diffraction']))
			if options.xaca:
				ice_correlation[j].append(N.array(f['data']['%.1fmm'%distances[j]]['ice']['correlation']))
			ice_angavg[j].append(N.array(f['data']['%.1fmm'%distances[j]]['ice']['angavg']))
			ice_angavgQ[j].append(N.array(f['data']['%.1fmm'%distances[j]]['ice']['angavg_Q']))
			
			#Gaussian peak fit statistics
			if options.peakfit:		
				p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2]
				index = N.array([((water_angavgQ[j][i] > options.S1_min)[k] and (water_angavgQ[j][i] < options.S1_max)[k]) or ((water_angavgQ[j][i] > options.S2_min)[k] and (water_angavgQ[j][i] < options.S2_max)[k]) for k in range(len(water_angavgQ[j][i]))])
				[p1, success] = optimize.leastsq(errfunc, p0[:], args=(water_angavgQ[j][i][index],water_angavg[j][i][index]))
				if success:
					water_fitint1[j].append(p1[0])
					water_fitpos1[j].append(p1[1])
					water_fitfwhm1[j].append(p1[2])
					water_fitint2[j].append(p1[3])
					water_fitpos2[j].append(p1[4])
					water_fitfwhm2[j].append(p1[5])
					water_fitdeltaq[j].append(p1[4]-p1[1])
				else:
					"Gaussian peak fit failed for %.1fmm in %s" % (distances[j], files[i])
					sys.exit(1)
		
		f.close()
	else:
		print "Could not find %s, aborting." % (original_dir + files[i])
		sys.exit(1)


for i in N.arange(len(distances)):
	currImg = img_class(water_pattern[i], water_angavg[i], water_angavgQ[i], "%.1fmm"%distances[i])
	currImg.draw_img_for_viewing_water()


#Gaussian peak fit statistics
if options.peakfit:
	#for i in range(len(files)):
	#	txttag = "peakfit-gdvn-%sADUs.txt"%(tags[i])
	for i in range(len(distances)):
		#txttag = "peakfit-gdvn_%.1fmm.txt"%(distances[i])
		#N.array([water_fitint1[i], water_fitpos1[i], water_fitfwhm1[i], water_fitint2[i], water_fitpos2[i], water_fitfwhm2[i]]).tofile(txttag, sep = "\n", format="%lf")
		#print "%s saved."%(txttag)
		csvtag = "peakfit-gdvn_%.1fmm.csv"%(distances[i])
		N.savetxt(csvtag, N.array([water_fitint1[i], water_fitpos1[i], water_fitfwhm1[i], water_fitint2[i], water_fitpos2[i], water_fitfwhm2[i]]), delimiter=",")
		print "%s saved."%(csvtag)
	
	print "include plots of peakfit_2gaussian_finalQ-S1_S2.png, peakfit_2gaussian_finalQ-S1_S2_with_ref-diff.png"
