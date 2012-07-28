#!/usr/bin/env python

# Written by J. Sellberg on 2012-07-24
# Reads the output HDF5 files from sorting and plots the ice data for the supporting material of the Nature manuscript

import numpy as N
from numpy import linalg as LA
import h5py as H
import glob as G
import pylab as P
import sys, os, re, shutil, subprocess, time
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value", default=10000)
parser.add_option("-f", "--file", action="store", type="string", dest="file", help="file name to extract ice data from", default="output_runs-aerojet-all_T+Q_20ADUs.h5")
(options, args) = parser.parse_args()

# with r0165
runs = [[114],[118,121],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
nhits_water = [[12077],[3308],[10143],[15702],[3508,789,4986],[320,203,104,2280,1159,919,207],[188,127,156,208,206],[30,75]] #r0151 actually has 1159 water hits but two had damaged Q-calibration, THESE ARE NOW INCLUDED WITH THE NEW CORRECTIONS
nhits_ice = [[1],[7],[7],[11],[2,2,8],[51,42,22,266,127,499,140],[780,625,425,732,1705],[2865,630]]
nevents = [[181089],[130074],[208489],[298050],[158756,62320,316507],[114978,194250,16146,309784,110923,211603,56893],[218817,166847,82800,198154,223197],[60809,213783]]
#colors = ['b','g','r','c','m','y','k']
colors = ['r','g','b','c','m','y','k']
#temperatures = [293.2,255.7,250.2,232.4,226.8,223.2,221.6,220.3]
temperatures = [293,256,250,232,227,223,222,220]
distances = [0.668972,11.548956,12.349269,21.407288,30.603057,40.403027,45.595702,50.619378]
# hexagonal ice peaks
HIceQ = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.222, '112':3.272, '201':3.324}
#HIcePos = {'100':10.5, '002':9.5, '101':8.5, '102':7.2, '110':9.2, '103':9.5, '200':12.5, '112':11.5, '201':10.5}
HIcePos = {'100':10.5, '002':9.5, '101':8.5, '102':6.2, '110':6.2, '103':6.2, '200':7., '112':6., '201':5.}
HIceQLabel = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.324, '112':3.374, '201':3.426}
# for imaging class
#HIceP = {'100':9.5, '002':8.5, '101':7.5, '102':6.2, '110':5.7, '103':5.2, '200':7., '112':6., '201':5.}
HIceP = {'100':10.5, '002':9.5, '101':8.5, '102':6.2, '110':6.2, '103':6.2, '200':7., '112':6., '201':5.}

source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
sorting_dir = "/reg/d/psdm/cxi/cxi25410/scratch/iceFinderCampaign/"
original_dir = os.getcwd() + '/'
run_tag = "r0166"

colmax = options.maxIntens
colmin = 0

#########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
#########################################################
class img_class (object):
	def __init__(self, inarr, inangavg, inangavg_Q, filename, title, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs, detectorDistance=eDD.get_detector_dist_in_meters(run_tag)):
		self.origarr = inarr.copy()
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.title = title
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
		
		fig = P.figure(num=None, figsize=(6.5, 12), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(211)
		canvas.set_title(self.title, fontsize='x-large')
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 100, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		canvas = fig.add_subplot(212)
		canvas.set_title("Angular Average", fontsize='large')
		P.xlabel(r"Momentum Transfer (${\AA}^{-1}$)", fontsize='large')
		P.ylabel("Radial Intensity (ADUs/srad)", fontsize='large')
		maxAngAvg = (self.inangavg).max()
		numQLabels = len(HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in HIceQ.iteritems():
			labelPosition = HIceP[i]*maxAngAvg/numQLabels
			P.axvline(j, 0, colmax, color='k', ls='--')
			P.text(HIceQ[i], labelPosition, str(i), rotation="45")
		
		P.plot(self.inangavg_Q, self.inangavg, color='r')
		#P.axis([1, 3.5, 0, 300000000])
		P.xlim([1, 3.5])
		canvas.set_ylim(bottom=0)
		
		#self.filename = self.filename + "-ice_only-angavg+pattern_20ADUs"
		self.filename = self.filename + "-ice_only-angavg+pattern_100ADUs"
		P.show()
	
	def draw_img_for_viewing_pattern(self):
		global HIceQ
		global HIceQLabel
		global HIceP
		global colmax
		global colmin
		print "Press 'p' to save PNG, 'e' to save EPS."
		
		fig = P.figure(num=None, figsize=(6, 5.5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title(self.title, fontsize='x-large')
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 100, vmin = 0)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = 150, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		#self.filename = self.filename + "-ice_only-pattern_20ADUs"
		self.filename = self.filename + "-ice_only-pattern_100ADUs"
		P.show()


pattern = []
angavg = []
angavgQ = []
#if os.path.exists(sorting_dir + options.file):
if os.path.exists(options.file):
	print 'found: ' + sorting_dir + options.file
	#f = H.File(sorting_dir + options.file, 'r')
	f = H.File(options.file, 'r')
	# loop over ice data
	for i in range(5,8):
		pattern.append(N.array(f['data']['%.1fmm'%(distances[i])]['ice']['diffraction']))
		angavg.append(N.array(f['data']['%.1fmm'%(distances[i])]['ice']['angavg']))
		angavgQ.append(N.array(f['data']['%.1fmm'%(distances[i])]['ice']['angavg_Q']))
	f.close()
else:
	print 'file missing: ' + sorting_dir + options.file
	sys.exit(1)


fig = P.figure(num=None, figsize=(10.5, 5), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(111)
#canvas.set_title("Angular Average", fontsize='large')
P.xlabel(r"Momentum Transfer (${\AA}^{-1}$)", fontsize='large')
#P.ylabel("Radial Intensity (ADUs/srad)", fontsize='large')
P.ylabel("Radial Intensity (ADUs/pixel)", fontsize='large')
for i in range(3):
	P.plot(angavgQ[i], angavg[i], label="T = %d K"%(temperatures[i+5]), color=colors[i])
#P.axis([1, 3.5, 0, 300000000])
P.xlim([1, 3.5])
handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left')

maxAngAvg = angavg[0].max()
numQLabels = len(HIceQ.keys())+1
for k,j in HIceQ.iteritems():
	labelPosition = HIcePos[k]*maxAngAvg/numQLabels
	P.axvline(j, 0, maxAngAvg, color='k', ls='--')
	P.text(HIceQ[k], labelPosition, str(k), rotation="45")

#P.savefig(original_dir + "output_runs-ice_T-ice_only-angavg_20ADUs.png")
#print "%s saved." % (original_dir + "output_runs-ice_T-ice_only-angavg_20ADUs.png")
#P.savefig(original_dir + "output_runs-ice_T-ice_only-angavg_20ADUs.eps", format='eps')
#print "%s saved." % (original_dir + "output_runs-ice_T-ice_only-angavg_20ADUs.eps")
P.savefig(original_dir + "output_runs-ice_T-ice_only-angavg_100ADUs.png")
print "%s saved." % (original_dir + "output_runs-ice_T-ice_only-angavg_100ADUs.png")
P.savefig(original_dir + "output_runs-ice_T-ice_only-angavg_100ADUs.eps", format='eps')
print "%s saved." % (original_dir + "output_runs-ice_T-ice_only-angavg_100ADUs.eps")
P.show()

for i in range(3):
	currImg = img_class(pattern[i], angavg[i], angavgQ[i], "output_runs-T%dK"%(temperatures[i+5]), "T = %d K"%(temperatures[i+5]))
	currImg.draw_img_for_viewing_pattern()
	currImg.draw_img_for_viewing_average()

