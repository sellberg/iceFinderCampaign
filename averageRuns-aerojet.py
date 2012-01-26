#!/usr/bin/env python

# Written by J. Sellberg on 2012-01-09
# Has the same functionality as analyzeSums-droplet_dispenser.py but for the new data format
# Modified by J.S. on 2012-01-25 to include the thresholding analysis

import numpy as N
from numpy import linalg as LA
import h5py as H
import glob as G
import os 
import re
import pylab as P
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value", default=2000)
(options, args) = parser.parse_args()

# with r0165
runs = [[114],[121],[118],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
nhits_water = [[12077],[3308],[10143],[15702],[3508,789,4986],[320,203,104,2280,1157,919,207],[188,127,156,208,206],[30,75]] #r0151 actually has 1159 water hits but two had damaged Q-calibration
nhits_ice = [[1],[7],[7],[11],[2,2,8],[51,42,22,266,127,499,140],[780,625,425,732,1705],[2865,630]]
# thresholded hits below 50 ADUs
runs = [[114],[121],[118],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
t50hits_water = [[4028],[536],[1305],[2032],[858,120,788],[78,64,14,545,305,498,102],[120,63,90,95,115],[6,36]]
t50hits_ice = [[0],[4],[3],[8],[0,0,5],[4,1,5,3,5,67,18],[312,285,27,34,830],[2402,167]]
# thresholded hits below 100 ADUs
#runs = [[114],[118,121],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
#t100hits_water = [[],[,],[],[,,],[194,,,,,,],[,,,,],[,]]
#t100hits_ice = [[],[,],[],[,,],[12,,,,,,],[,,,,],[,]]
#colors = ['b','g','r','c','m','y','k']
colors = ['r','g','b','c','m','y','k']
#temperatures = [264,235,234,227,224,221,220,219]
temperatures = [270,253,252,232,227,223,221,220]
#distances = [1.085200,11.585200,11.985100,21.585200,30.540800,40.540800,45.540800,50.540800]
distances = [0.668972,11.548956,12.349269,21.407288,30.603057,40.403027,45.595702,50.619378]

# hexagonal ice peaks
HIceQ = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.222, '112':3.272, '201':3.324}
HIcePos = {'100':10.5, '002':9.5, '101':8.5, '102':6.2, '110':6.2, '103':6.2, '200':7., '112':6., '201':5.}
HIceQLabel = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.324, '112':3.374, '201':3.426}
# for imaging class
HIceP = {'100':9.5, '002':8.5, '101':7.5, '102':6.2, '110':5.7, '103':5.2, '200':7., '112':6., '201':5.}

source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
sorting_dir = "/reg/d/psdm/cxi/cxi25410/scratch/iceFinderCampaign/"
original_dir = os.getcwd() + '/'
run_tag = "r0144"

water_pattern = []
water_angavg = []
water_angavgQ = []
ice_pattern = []
ice_angavg = []
ice_angavgQ = []

colmax = options.maxIntens
colmin = 0

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
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = 300, vmin = 0)
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
		global HIceQ
		global HIceQLabel
		global HIceP
		global colmax
		global colmin
		self.inangavg = self.inangavg*(self.inangavg>0)
		#print "Press 'p' to save PNG, 'e' to save EPS."
		
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		fig.subplots_adjust(wspace=0.1)
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title("Water (%s Hits)"%(self.inangavg_Q[0]), fontsize=22, fontname='sans-serif', fontweight='roman')
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 200, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		canvas = fig.add_subplot(122)
		canvas.set_title("Ice (%s Hits)"%(self.inangavg_Q[1]), fontsize=22, fontname='sans-serif', fontweight='roman')
		self.axes = P.imshow(self.inangavg, origin='lower', vmax = colmax, vmin = colmin)
		#self.axes = P.imshow(self.inangavg, origin='lower', vmax = 100, vmin = 0)
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
	temp_water_angavg = []
	temp_water_angavgQ = []
	temp_ice_pattern = []
	temp_ice_angavg = []
	temp_ice_angavgQ = []
	for j in N.arange(len(runs[i])):
		run_tag = "r0%s"%(runs[i][j])
		run_dir =  'output_' + run_tag + '/'
		#water
		if (nhits_water[i][j] != 0):
			if os.path.exists(sorting_dir + run_dir + run_tag + '_type0.h5'):
				print 'found: ' + sorting_dir + run_dir + run_tag + '_type0.h5'
				f = H.File(sorting_dir + run_dir + run_tag + '_type0.h5', 'r')
				temp_water_pattern.append(N.array(f['data']['diffraction']))
				temp_water_angavg.append(N.array(f['data']['angavg']))
				temp_water_angavgQ.append(N.array(f['data']['angavgQ']))
				f.close()
			else:
				print 'file missing: ' + sorting_dir + run_dir + run_tag + '_type0.h5'
		#ice
		if (nhits_ice[i][j] != 0):
			if os.path.exists(sorting_dir + run_dir + 'type1/' + run_tag + '_type1.h5'):
				print 'found: ' + sorting_dir + run_dir + 'type1/' + run_tag + '_type1.h5'
				f = H.File(sorting_dir + run_dir + 'type1/' + run_tag + '_type1.h5', 'r')
				temp_ice_pattern.append(N.array(f['data']['diffraction']))
				temp_ice_angavg.append(N.array(f['data']['angavg']))
				temp_ice_angavgQ.append(N.array(f['data']['angavgQ']))
				f.close()
			else:
				print 'file missing: ' + sorting_dir + run_dir + 'type1/' + run_tag + '_type1.h5'
	
	
	#plot temp_angavg
	fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(121)
	canvas.set_title("Water, T = %s K"%(temperatures[i]))
	P.xlabel("Q (Angstroms-1)")
	P.ylabel("Average Intensity (ADUs)")
	for j in N.arange(len(runs[i])):
		P.plot(temp_water_angavgQ[j], temp_water_angavg[j], color=colors[j], label="r0%s"%(runs[i][j]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels)
	canvas = fig.add_subplot(122)
	canvas.set_title("Ice, T = %s K"%(temperatures[i]))
	P.xlabel("Q (Angstroms-1)")
	P.ylabel("Average Intensity (ADUs)")
	for j in N.arange(len(runs[i])):
		P.plot(temp_ice_angavgQ[j], temp_ice_angavg[j], color=colors[j], label="r0%s"%(runs[i][j]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper left')
	
	maxAngAvg = max([temp_ice_angavg[j].max() for j in N.arange(len(temp_ice_angavg))])
	numQLabels = len(HIceQ.keys())+1
	for k,j in HIceQ.iteritems():
		labelPosition = HIcePos[k]*maxAngAvg/numQLabels
		P.axvline(j, 0, maxAngAvg, color='k', ls='--')
		P.text(HIceQLabel[k], labelPosition, str(k), rotation="45")
	
	print "output_runs-T%sK-angavg.png saved."%(temperatures[i])
	P.savefig(original_dir + "output_runs-T%sK-angavg.png"%(temperatures[i]))
	#P.savefig(original_dir + "output_runs-T%sK-angavg.eps"%(temperatures[i]), format='eps')
	#P.show()
	P.close()
	
	
	sumwater = float(sum(nhits_water[i]) - sum(t50hits_water[i]))
	sumice = float(sum(nhits_ice[i]) - sum(t50hits_ice[i]))
	for j in N.arange(len(runs[i])):
		nwater = nhits_water[i][j] - t50hits_water[i][j]
		nice = nhits_ice[i][j] - t50hits_ice[i][j]
		if (nwater > 0):
			temp_water_angavg[j] *= nwater/sumwater
			temp_water_pattern[j] *= nwater/sumwater
		if (nice > 0):
			temp_ice_angavg[j] *= nice/sumice
			temp_ice_pattern[j] *= nice/sumice
	
	water_angavg.append(N.array(temp_water_angavg).sum(axis=0))
	water_angavgQ.append(N.array(temp_water_angavgQ).sum(axis=0)/len(runs[i]))
	water_pattern.append(N.array(temp_water_pattern).sum(axis=0))
	ice_angavg.append(N.array(temp_ice_angavg).sum(axis=0))
	ice_angavgQ.append(N.array(temp_ice_angavgQ).sum(axis=0)/len(runs[i]))
	ice_pattern.append(N.array(temp_ice_pattern).sum(axis=0))
	
	currImg = img_class(water_pattern[i], ice_pattern[i], [int(sumwater), int(sumice)], "output_runs-T%sK-pattern"%(temperatures[i]))
	currImg.draw_img_for_viewing_pattern()

#plot angavg
fig = P.figure(num=None, figsize=(13.5, 10), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(221)
canvas.set_title("Water", fontsize='x-large')
P.xlabel("Q (Angstroms-1)", fontsize='x-large')
P.ylabel("Average Intensity (ADUs)", fontsize='x-large')
for i in N.arange(len(runs)):
	#if i > 3:
		#P.plot(water_angavgQ[i], water_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))
	P.plot(water_angavgQ[i], water_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left', bbox_to_anchor = (1, 1))

canvas = fig.add_subplot(223)
canvas.set_title("Ice", fontsize='x-large')
P.xlabel("Q (Angstroms-1)", fontsize='x-large')
P.ylabel("Average Intensity (ADUs)", fontsize='x-large')
for i in N.arange(len(runs)):
	#if i > 3:
		#P.plot(ice_angavgQ[i], ice_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))
	P.plot(ice_angavgQ[i], ice_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left', bbox_to_anchor = (1, 1))

maxAngAvg = max([ice_angavg[i].max() for i in N.arange(len(ice_angavg))])
#maxAngAvg = max([ice_angavg[i+4].max() for i in N.arange(3)])
numQLabels = len(HIceQ.keys())+1
for k,j in HIceQ.iteritems():
	labelPosition = HIcePos[k]*maxAngAvg/numQLabels
	P.axvline(j, 0, maxAngAvg, color='k', ls='--')
	P.text(HIceQLabel[k], labelPosition, str(k), rotation="45")

print "output_runs-aerojet-all_T-angavg_Q.png saved."
#P.savefig(original_dir + "output_runs-aerojet-ice_T-angavg_Q.png")
#P.savefig(original_dir + "output_runs-aerojet-ice_T-angavg_Q.eps", format='eps')
P.savefig(original_dir + "output_runs-aerojet-all_T-angavg_Q.png")
#P.savefig(original_dir + "output_runs-aerojet-all_T-angavg_Q.eps", format='eps')
#P.show()
P.close()

#save to file
f = H.File(original_dir + "output_runs-aerojet-all_T+Q_50ADUs.h5", 'w')
entry_1 = f.create_group("/data")
for i in N.arange(len(runs)):
	entry_2 = f.create_group("/data/%.1fmm"%(distances[i]))
	entry_3 = f.create_group("/data/%.1fmm/water"%(distances[i]))
	entry_3.create_dataset("diffraction", data=water_pattern[i])
	entry_3.create_dataset("angavg", data=water_angavg[i])
	entry_3.create_dataset("angavg_Q", data=water_angavgQ[i])
	entry_4 = f.create_group("/data/%.1fmm/ice"%(distances[i]))
	entry_4.create_dataset("diffraction", data=ice_pattern[i])
	entry_4.create_dataset("angavg", data=ice_angavg[i])
	entry_4.create_dataset("angavg_Q", data=ice_angavgQ[i])

f.close()
print "Successfully updated output_runs-aerojet-all_T+Q_50ADUs.h5"
