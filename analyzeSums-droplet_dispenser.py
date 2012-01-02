#!/usr/bin/env python

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
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value", default=10000)
(options, args) = parser.parse_args()

# adapted to droplet dispenser (water only)
runs = [[105],[111],[68],[63,64],[110],[73],[77,79],[88],[97],[102]]
nhits_water = [[2610],[7860],[3900],[381,3100],[3000],[2280],[1225,1480],[1880],[1820],[1675]]
#colors = ['b','g','r','c','m','y','k']
colors = ['r','g','b','c','m','y','k']
temperatures = [290,260,254,249,244,243,240,235,233,232]
# hexagonal ice peaks
HIceQ = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.222, '112':3.272, '201':3.324}
HIcePos = {'100':10.5, '002':9.5, '101':8.5, '102':6.2, '110':6.2, '103':6.2, '200':7., '112':6., '201':5.}
HIceQLabel = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.324, '112':3.374, '201':3.426}
# for imaging class
HIceP = {'100':9.5, '002':8.5, '101':7.5, '102':6.2, '110':5.7, '103':5.2, '200':7., '112':6., '201':5.}
# single shot analysis (not used)
hits = ["LCLS_2011_Feb28_r0144_163006_2feb","LCLS_2011_Feb28_r0144_162957_2322"]

source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
sorting_dir = "/reg/d/psdm/cxi/cxi25410/scratch/iceFinderCampaign/"
original_dir = os.getcwd() + '/'
run_tag = "r0144"

water_angavg = []
water_pattern = []
Q_angavg = []

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
		global colmax
		global colmin
		self.inangavg = self.inangavg*(self.inangavg>0)
		print "Press 'p' to save PNG, 'e' to save EPS."
		
		fig = P.figure(num=None, figsize=(8, 5), dpi=100, facecolor='w', edgecolor='k')
		fig.subplots_adjust(wspace=0.1)
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title("Water T=%sK"%(self.inangavg_Q[1]) + " (%s Hits)"%(self.inangavg_Q[0]), fontsize=22, fontname='sans-serif', fontweight='roman')
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = 500, vmin = 0)
		#self.axes = P.imshow(self.inarr, origin='lower', vmax = 200, vmin = 0)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		P.show()

"""
for file in hits:
	run_tag = "r0144"
	run_dir =  'output_' + run_tag + '/'
	filePath = sorting_dir + run_dir + 'type1/' + file + '_cspad.h5'
	if os.path.exists(filePath):
		print 'found: ' + filePath
		f = H.File(filePath, 'r')
		d = N.array(f['/data/data']).astype(float)
		currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
		currDetectorDist=(1.E-3)*f['LCLS']['detectorPosition'][0]
		f.close()
	else:
		print 'file missing: ' + filePath
	
	angAvgPath = sorting_dir + run_dir + 'type1/' + file + '_cspad-angavg.h5'
	if os.path.exists(angAvgPath):
		print 'found: ' + angAvgPath
		f = H.File(angAvgPath, 'r')
		davg = N.array(f['data']['data'][0])
		f.close()
	else:
		print 'file missing: ' + angAvgPath
	
	angAvgQPath = source_dir + run_tag + '/' + run_tag + "-angavg_Q.h5"
	if os.path.exists(angAvgQPath):
		print 'found: ' + angAvgQPath
		f = H.File(angAvgQPath,'r')
		davgQ = N.array(f['/data/data'][0])
		f.close()
	else:
		print 'file missing: ' + angAvgQPath
	
	currImg = img_class(d, davg, davgQ, file, meanWaveLengthInAngs=currWavelengthInAngs, detectorDistance=currDetectorDist)
	currImg.draw_img_for_viewing_average()
"""

# FAST VERSION (premade sums)
for i in N.arange(len(runs)):
	temp_water_angavg = []
	temp_water_pattern = []
	temp_Q_angavg = []
	for j in N.arange(len(runs[i])):
		if (runs[i][j] < 100):
			run_tag = "r00%s"%(runs[i][j])
		else:
			run_tag = "r0%s"%(runs[i][j])
		run_dir =  'output_' + run_tag + '/'
		#water
		if (nhits_water[i][j] != 0):
			if os.path.exists(sorting_dir + run_dir + run_tag + '_type0.h5'):
				print 'found: ' + sorting_dir + run_dir + run_tag + '_type0.h5'
				f = H.File(sorting_dir + run_dir + run_tag + '_type0.h5', 'r')
				temp_water_angavg.append(N.array(f['data']['angavg']))
				temp_water_pattern.append(N.array(f['data']['diffraction']))
				f.close()
			else:
				print 'file missing: ' + sorting_dir + run_dir + run_tag + '_type0.h5'
		#qcal
		if os.path.exists(source_dir + run_tag + '/' + run_tag + "-angavg_Q.h5"):
			print 'found: ' + source_dir + run_tag + '/' + run_tag + "-angavg_Q.h5"
			f = H.File(source_dir + run_tag + '/' + run_tag + "-angavg_Q.h5",'r')
			temp_Q_angavg.append(N.array(f['/data/data'])[0])
			f.close()
		else:
			print 'file missing: ' + source_dir + run_tag + '/' + run_tag + "-angavg_Q.h5"
	
	
	#plot temp_angavg
	fig = P.figure(num=None, figsize=(8, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(111)
	canvas.set_title("Water, T = %s K"%(temperatures[i]))
	P.xlabel("Q (Angstroms-1)")
	P.ylabel("Average Intensity (ADUs)")
	for j in N.arange(len(runs[i])):
		P.plot(temp_Q_angavg[j], temp_water_angavg[j], color=colors[j], label="r0%s"%(runs[i][j]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels)
	P.savefig(original_dir + "output_runs-T%sK-angavg.png"%(temperatures[i]))
	#P.savefig(original_dir + "output_runs-T%sK-angavg.eps"%(temperatures[i]), format='eps')
	P.show()
	
	
	sumwater = float(sum(nhits_water[i]))
	for j in N.arange(len(runs[i])):
		nwater = nhits_water[i][j]
		if (nwater != 0):
			temp_water_angavg[j] *= nwater/sumwater
			temp_water_pattern[j] *= nwater/sumwater
	
	
	water_angavg.append(N.array(temp_water_angavg).sum(axis=0))
	water_pattern.append(N.array(temp_water_pattern).sum(axis=0))
	Q_angavg.append(N.array(temp_Q_angavg).mean(axis=0))
	
	currImg = img_class(water_pattern[i], water_angavg[i], [int(sumwater), temperatures[i]], "output_runs-T%sK-pattern"%(temperatures[i]))
	currImg.draw_img_for_viewing_pattern()

#plot angavg
fig = P.figure(num=None, figsize=(8, 5), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(111)
canvas.set_title("Water", fontsize='x-large')
P.xlabel("Q (Angstroms-1)", fontsize='x-large')
P.ylabel("Average Intensity (ADUs)", fontsize='x-large')
for i in N.arange(len(runs)):
	#if i > 3:
		#P.plot(Q_angavg[i], water_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))
	P.plot(Q_angavg[i], water_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels)
P.savefig(original_dir + "output_runs-droplet_dispenser-all_T-angavg_Q.png")
#P.savefig(original_dir + "output_runs-droplet_dispenser-all_T-angavg_Q.eps", format='eps')
P.show()

#save to file
f = H.File(original_dir + "output_runs-droplet_dispenser-all_T+Q.h5", 'w')
entry_1 = f.create_group("/data")
for i in N.arange(len(runs)):
	entry_2 = f.create_group("/data/T%sK"%(temperatures[i]))
	entry_3 = f.create_group("/data/T%sK/water"%(temperatures[i]))
	entry_3.create_dataset("diffraction", data=water_pattern[i])
	entry_3.create_dataset("angavg", data=water_angavg[i])
	entry_3.create_dataset("angavg_Q", data=Q_angavg[i])

f.close()
