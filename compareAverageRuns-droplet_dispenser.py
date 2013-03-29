#!/usr/bin/env python

# Written by J. Sellberg on 2013-03-28
# Compares the output HDF5 files of averageRuns-droplet_dispenser.py from various thresholds/sorting algorithms

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

files = ["output_runs-droplet_dispenser-all_T+Q_20ADUs-failedFits.h5", "output_runs-droplet_dispenser-all_T+Q_100ADUs-failedFits.h5"]
tags = ["20ADUs-without_failedFits", "100ADUs-without_failedFits"]
colors = ['r','g','b','c','m','y','k']
#temperatures = [290,260,254,249,244,243,240,235,233]
#temperatures = [283,255,249,245,241,240,237,233,231] # used by averageRuns-droplet_dispenser.py
temperatures = [283,257,250,246,241,241,237,233,231] # FINAL temperatures LicThesis
#distances = [0.750000,10.000000,15.681899,20.681899,30.000000,30.699899,40.699899,60.699899,70.649900]
distances = [0.750000,10.011910,15.798895,20.789092,30.006612,30.738471,40.759723,60.710063,70.662700] # FINAL distances

# reference values from supplementary material Huang_resubmitted120910_Supplementary_Notes.pdf
ref_temp = [257, 250, 246, 241, 240, 237, 233, 231]
ref_s1 = N.array([1.992, 1.950, 1.939, 1.925, 1.915, 1.899, 1.887, 1.872])
ref_s2 = N.array([3.018, 3.016, 3.015, 3.039, 3.025, 3.025, 3.052, 3.051])
ref_dq = ref_s2 - ref_s1

# SCRATCH
#source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
#sorting_dir = "/reg/d/psdm/cxi/cxi25410/scratch/iceFinderCampaign/"
# RES & FTC
source_dir = "/reg/d/psdm/cxi/cxi25410/ftc/cleaned_hdf5/"
#source_dir = "/reg/d/psdm/cxi/cxi25410/res/cleaned_hdf5/"
sorting_dir = "/reg/d/psdm/cxi/cxi25410/res/iceFinderCampaign/"

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
	def __init__(self, inarr, inangavg, inangavg_Q, filename, filetag="", normalize=False, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs, detectorDistance=eDD.get_detector_dist_in_meters(run_tag)):
		self.inarr = N.array(inarr).mean(axis=0)
		self.inarr = self.inarr*(self.inarr>0)
		self.filename = filename
		if (filetag != ""):
			self.filetag = filetag[:] #this is necessary to pass filename by value and not by reference (which would change the input parameters outside of the class)
			for i in range(len(self.filetag)):
				self.filetag[i] += ": "
		else:
			self.filetag = ["" for i in inangavg]
		if normalize:
			self.angavgmax = N.zeros(len(inangavg))
			for i in range(len(inangavg)):
				self.angavgmax[i] = inangavg[i].max()
			self.angavgmaxindex = self.angavgmax.argmax()
			self.inangavg = []
			for i in range(len(inangavg)):
				self.inangavg.append(inangavg[i]*self.angavgmax[self.angavgmaxindex]/self.angavgmax[i])
		else:
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
				p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2]
				index = N.array([((self.inangavgQ[i] > options.S1_min)[j] and (self.inangavgQ[i] < options.S1_max)[j]) or ((self.inangavgQ[i] > options.S2_min)[j] and (self.inangavgQ[i] < options.S2_max)[j]) for j in range(len(self.inangavgQ[i]))])
				[p1, success] = optimize.leastsq(errfunc, p0[:], args=(self.inangavgQ[i][index],self.inangavg[i][index]))
				if success:
					P.plot(self.inangavgQ[i], self.inangavg[i], color="%s"%colors[i%7], label="%sS1 = %.3f A-1, S2 = %.3f A-1" % (self.filetag[i], p1[1], p1[4]))
					P.plot(self.inangavgQ[i], fitfunc(p1, self.inangavgQ[i]), "%s:"%colors[i%7])
				else:
					P.plot(self.inangavgQ[i], self.inangavg[i], "%s-"%colors[i%7])
			else:
				P.plot(self.inangavgQ[i], self.inangavg[i], "%s-"%colors[i%7])
		
		handles, labels = canvas.get_legend_handles_labels()
		if (len(self.inangavg) > 6):
			canvas.legend(handles, labels, loc='upper left', prop={'size':6})
		else:
			canvas.legend(handles, labels, loc='lower right', prop={'size':6})
		canvas.set_title("Angular Average")
		P.xlabel("Q (A-1)")
		P.ylabel("I(Q) (ADU/srad)")
		pngtag = original_dir + "peakfit-dod_%s.png" % (self.filename)
		P.savefig(pngtag)
		print "%s saved." % (pngtag)
		epstag = original_dir + "peakfit-dod_%s.eps" % (self.filename)
		P.savefig(epstag, format='eps')
		print "%s saved." % (epstag)
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
	singleshot_water_fitpos1_mean = [[] for i in files]
	singleshot_water_fitpos1_median = [[] for i in files]
	singleshot_water_fitpos1_std = [[] for i in files]
	singleshot_water_fitpos2_mean = [[] for i in files]
	singleshot_water_fitpos2_median = [[] for i in files]
	singleshot_water_fitpos2_std = [[] for i in files]
	singleshot_water_fitdeltaq_mean = [[] for i in files]
	singleshot_water_fitdeltaq_median = [[] for i in files]
	singleshot_water_fitdeltaq_std = [[] for i in files]

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
			#Gaussian peak fit statistics
			if options.peakfit:
				#Read Gaussian peak fit statistics from single-shot distributions
				singleshot_water_fitpos1_mean[i].append(N.mean(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['pos1']))
				singleshot_water_fitpos1_median[i].append(N.median(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['pos1']))
				singleshot_water_fitpos1_std[i].append(N.std(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['pos1']))
				singleshot_water_fitpos2_mean[i].append(N.mean(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['pos2']))
				singleshot_water_fitpos2_median[i].append(N.median(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['pos2']))
				singleshot_water_fitpos2_std[i].append(N.std(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['pos2']))
				singleshot_water_fitdeltaq_mean[i].append(N.mean(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['deltaQ']))
				singleshot_water_fitdeltaq_median[i].append(N.median(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['deltaQ']))
				singleshot_water_fitdeltaq_std[i].append(N.std(f['data']['%.1fmm'%distances[j]]['water']['peakFit']['deltaQ']))
				
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
	currImg = img_class(water_pattern[i], water_angavg[i], water_angavgQ[i], "%.1fmm"%distances[i], tags)
	currImg.draw_img_for_viewing_water()
	currImg = img_class(water_pattern[i], water_angavg[i], water_angavgQ[i], "%.1fmm-norm"%distances[i], tags, True)
	currImg.draw_img_for_viewing_water()


temp_water_pattern = [[] for i in files]
if options.xaca:
	temp_water_correlation = [[] for i in files]
temp_water_angavg = [[] for i in files]
temp_water_angavgQ = [[] for i in files]

temperature_tags = []
for i in range(len(distances)):
	temperature_tags.append("T%sK"%temperatures[i])
	for j in range(len(files)):
		temp_water_pattern[j].append(water_pattern[i][j])
		if options.xaca:
			temp_water_correlation[j].append(water_correlation[i][j])
		temp_water_angavg[j].append(water_angavg[i][j])
		temp_water_angavgQ[j].append(water_angavgQ[i][j])

for i in range(len(files)):
	currImg = img_class(temp_water_pattern[i], temp_water_angavg[i], temp_water_angavgQ[i], "%s"%tags[i], temperature_tags)
	currImg.draw_img_for_viewing_water()
	currImg = img_class(temp_water_pattern[i], temp_water_angavg[i], temp_water_angavgQ[i], "%s-norm"%tags[i], temperature_tags, True)
	currImg.draw_img_for_viewing_water()

#Gaussian peak fit statistics
if options.peakfit:
	temp_water_fitint1 = [[] for i in files]
	temp_water_fitpos1 = [[] for i in files]
	temp_water_fitfwhm1 = [[] for i in files]
	temp_water_fitint2 = [[] for i in files]
	temp_water_fitpos2 = [[] for i in files]
	temp_water_fitfwhm2 = [[] for i in files]
	temp_water_fitdeltaq = [[] for i in files]
	
	for i in range(len(distances)):
		for j in range(len(files)):
			temp_water_fitint1[j].append(water_fitint1[i][j])
			temp_water_fitpos1[j].append(water_fitpos1[i][j])
			temp_water_fitfwhm1[j].append(water_fitfwhm1[i][j])
			temp_water_fitint2[j].append(water_fitint2[i][j])
			temp_water_fitpos2[j].append(water_fitpos2[i][j])
			temp_water_fitfwhm2[j].append(water_fitfwhm2[i][j])
			temp_water_fitdeltaq[j].append(water_fitdeltaq[i][j])
	
	for i in range(len(files)):
		csvtag = "peakfit-dod_%s.csv"%(tags[i])
		N.savetxt(csvtag, N.array([temp_water_fitint1[i], temp_water_fitpos1[i], temp_water_fitfwhm1[i], temp_water_fitint2[i], temp_water_fitpos2[i], temp_water_fitfwhm2[i], singleshot_water_fitpos1_mean[i], singleshot_water_fitpos1_median[i], singleshot_water_fitpos1_std[i], singleshot_water_fitpos2_mean[i], singleshot_water_fitpos2_median[i], singleshot_water_fitpos2_std[i], singleshot_water_fitdeltaq_mean[i], singleshot_water_fitdeltaq_median[i], singleshot_water_fitdeltaq_std[i]]), delimiter=",")
		print "%s saved."%(csvtag)
	
	#FIRST PLOT
	fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(131)
	canvas.set_title("S1 / S2")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_s1, color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(ref_temp, ref_s2, color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(temperatures, temp_water_fitpos1[i], color="%s"%colors[i], marker='o', label="%s" % (tags[i]))
		P.plot(temperatures, temp_water_fitpos2[i], color="%s"%colors[i], linestyle='--', marker='o')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([1.8, 3.2])
	
	canvas = fig.add_subplot(132)
	canvas.set_title("S1 / S2 mean centered")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_s1 - ref_s1.mean(), color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(ref_temp, ref_s2 - ref_s2.mean(), color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(temperatures, N.array(temp_water_fitpos1[i]) - N.array(temp_water_fitpos1[i])[2:].mean(), color="%s"%colors[i], marker='o', label="%s" % (tags[i]))
		P.plot(temperatures, N.array(temp_water_fitpos2[i]) - N.array(temp_water_fitpos2[i])[2:].mean(), color="%s"%colors[i], linestyle='--', marker='o')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([-0.06, 0.1])
	
	canvas = fig.add_subplot(133)
	canvas.set_title("deltaQ")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_dq, color='k', marker='o', label="Huang_resubmitted120910")
	for i in range(len(files)):
		P.plot(temperatures, temp_water_fitdeltaq[i], color="%s"%colors[i], marker='o', label="%s" % (tags[i]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper right', prop={'size':6})
	canvas.set_ylim([0.95, 1.25])
	
	P.savefig(original_dir + "peakfit-dod_vs_T.png")
	print "peakfit-dod_vs_T.png saved."
	P.savefig(original_dir + "peakfit-dod_vs_T.eps", format='eps')
	print "peakfit-dod_vs_T.eps saved."
	#P.show()
	P.close()


	#SECOND PLOT
	fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(131)
	canvas.set_title("S1 / S2")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_s1, color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(distances[1:], ref_s2, color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(distances, temp_water_fitpos1[i], color="%s"%colors[i], marker='o', label="%s" % (tags[i]))
		P.plot(distances, temp_water_fitpos2[i], color="%s"%colors[i], linestyle='--', marker='o')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([1.8, 3.2])
	
	canvas = fig.add_subplot(132)
	canvas.set_title("S1 / S2 mean centered")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_s1 - ref_s1.mean(), color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(distances[1:], ref_s2 - ref_s2.mean(), color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(distances, N.array(temp_water_fitpos1[i]) - N.array(temp_water_fitpos1[i])[2:].mean(), color="%s"%colors[i], marker='o', label="%s" % (tags[i]))
		P.plot(distances, N.array(temp_water_fitpos2[i]) - N.array(temp_water_fitpos2[i])[2:].mean(), color="%s"%colors[i], linestyle='--', marker='o')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper right', prop={'size':6})
	canvas.set_ylim([-0.06, 0.1])
	
	canvas = fig.add_subplot(133)
	canvas.set_title("deltaQ")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_dq, color='k', marker='o', label="Huang_resubmitted120910")
	for i in range(len(files)):
		P.plot(distances, temp_water_fitdeltaq[i], color="%s"%colors[i], marker='o', label="%s" % (tags[i]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper left', prop={'size':6})
	canvas.set_ylim([0.95, 1.25])
	
	P.savefig(original_dir + "peakfit-dod_vs_dist.png")
	print "peakfit-dod_vs_dist.png saved."
	P.savefig(original_dir + "peakfit-dod_vs_dist.eps", format='eps')
	print "peakfit-dod_vs_dist.eps saved."
	#P.show()
	P.close()


	#THIRD PLOT
	fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(131)
	canvas.set_title("S1 / S2")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_s1, color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(ref_temp, ref_s2, color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(temperatures, singleshot_water_fitpos1_median[i], color="%s"%colors[i], linestyle='-', marker='D', label="%s" % (tags[i]))
		P.plot(temperatures, singleshot_water_fitpos2_median[i], color="%s"%colors[i], linestyle='--', marker='D')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([1.8, 3.2])
	
	canvas = fig.add_subplot(132)
	canvas.set_title("S1 / S2 mean centered")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_s1 - ref_s1.mean(), color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(ref_temp, ref_s2 - ref_s2.mean(), color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(temperatures, N.array(singleshot_water_fitpos1_median[i]) - N.array(singleshot_water_fitpos1_median[i])[2:].mean(), color="%s"%colors[i], linestyle='-', marker='D', label="%s" % (tags[i]))
		P.plot(temperatures, N.array(singleshot_water_fitpos2_median[i]) - N.array(singleshot_water_fitpos2_median[i])[2:].mean(), color="%s"%colors[i], linestyle='--', marker='D')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([-0.06, 0.1])
	
	canvas = fig.add_subplot(133)
	canvas.set_title("deltaQ")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_dq, color='k', marker='o', label="Huang_resubmitted120910")
	for i in range(len(files)):
		P.plot(temperatures, singleshot_water_fitdeltaq_median[i], color="%s"%colors[i], linestyle='-', marker='D', label="%s" % (tags[i]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper right', prop={'size':6})
	canvas.set_ylim([0.95, 1.25])
	
	P.savefig(original_dir + "peakfit-dod_single-shot_median_vs_T.png")
	print "peakfit-dod_single-shot_median_vs_T.png saved."
	P.savefig(original_dir + "peakfit-dod_single-shot_median_vs_T.eps", format='eps')
	print "peakfit-dod_single-shot_median_vs_T.eps saved."
	#P.show()
	P.close()


	#FOURTH PLOT
	fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(131)
	canvas.set_title("S1 / S2")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_s1, color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(distances[1:], ref_s2, color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(distances, singleshot_water_fitpos1_median[i], color="%s"%colors[i], linestyle='-', marker='D', label="%s" % (tags[i]))
		P.plot(distances, singleshot_water_fitpos2_median[i], color="%s"%colors[i], linestyle='--', marker='D')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([1.8, 3.2])
	
	canvas = fig.add_subplot(132)
	canvas.set_title("S1 / S2 mean centered")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_s1 - ref_s1.mean(), color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(distances[1:], ref_s2 - ref_s2.mean(), color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.plot(distances, N.array(singleshot_water_fitpos1_median[i]) - N.array(singleshot_water_fitpos1_median[i])[2:].mean(), color="%s"%colors[i], linestyle='-', marker='D', label="%s" % (tags[i]))
		P.plot(distances, N.array(singleshot_water_fitpos2_median[i]) - N.array(singleshot_water_fitpos2_median[i])[2:].mean(), color="%s"%colors[i], linestyle='--', marker='D')
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper right', prop={'size':6})
	canvas.set_ylim([-0.06, 0.1])
	
	canvas = fig.add_subplot(133)
	canvas.set_title("deltaQ")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_dq, color='k', marker='o', label="Huang_resubmitted120910")
	for i in range(len(files)):
		P.plot(distances, singleshot_water_fitdeltaq_median[i], color="%s"%colors[i], linestyle='-', marker='D', label="%s" % (tags[i]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper left', prop={'size':6})
	canvas.set_ylim([0.95, 1.25])
	
	P.savefig(original_dir + "peakfit-dod_single-shot_median_vs_dist.png")
	print "peakfit-dod_single-shot_median_vs_dist.png saved."
	P.savefig(original_dir + "peakfit-dod_single-shot_median_vs_dist.eps", format='eps')
	print "peakfit-dod_single-shot_median_vs_dist.eps saved."
	#P.show()
	P.close()
	
	
	#FIFTH PLOT
	fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(131)
	canvas.set_title("S1 / S2")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_s1, color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(ref_temp, ref_s2, color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.errorbar(temperatures, singleshot_water_fitpos1_mean[i], yerr=singleshot_water_fitpos1_std[i], fmt="-%ss"%colors[i], label="%s" % (tags[i]))
		P.errorbar(temperatures, singleshot_water_fitpos2_mean[i], yerr=singleshot_water_fitpos2_std[i], fmt="--%ss"%colors[i])
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([1.8, 3.2])
	
	canvas = fig.add_subplot(132)
	canvas.set_title("S1 / S2 mean centered")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_s1 - ref_s1.mean(), color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(ref_temp, ref_s2 - ref_s2.mean(), color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.errorbar(temperatures, N.array(singleshot_water_fitpos1_mean[i]) - N.array(singleshot_water_fitpos1_mean[i])[2:].mean(), yerr=singleshot_water_fitpos1_std[i], fmt="-%ss"%colors[i], label="%s" % (tags[i]))
		P.errorbar(temperatures, N.array(singleshot_water_fitpos2_mean[i]) - N.array(singleshot_water_fitpos2_mean[i])[2:].mean(), yerr=singleshot_water_fitpos2_std[i], fmt="--%ss"%colors[i])
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([-0.06, 0.1])
	
	canvas = fig.add_subplot(133)
	canvas.set_title("deltaQ")
	P.xlabel("T (K)")
	P.ylabel("Q (A-1)")
	P.plot(ref_temp, ref_dq, color='k', marker='o', label="Huang_resubmitted120910")
	for i in range(len(files)):
		P.errorbar(temperatures, singleshot_water_fitdeltaq_mean[i], yerr=singleshot_water_fitdeltaq_std[i], fmt="-%ss"%colors[i], label="%s" % (tags[i]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper right', prop={'size':6})
	canvas.set_ylim([0.95, 1.25])
	
	P.savefig(original_dir + "peakfit-dod_single-shot_mean_vs_T.png")
	print "peakfit-dod_single-shot_mean_vs_T.png saved."
	P.savefig(original_dir + "peakfit-dod_single-shot_mean_vs_T.eps", format='eps')
	print "peakfit-dod_single-shot_mean_vs_T.eps saved."
	#P.show()
	P.close()
	
	
	#SIXTH PLOT
	fig = P.figure(num=None, figsize=(18.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(131)
	canvas.set_title("S1 / S2")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_s1, color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(distances[1:], ref_s2, color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.errorbar(distances, singleshot_water_fitpos1_mean[i], yerr=singleshot_water_fitpos1_std[i], fmt="-%ss"%colors[i], label="%s" % (tags[i]))
		P.errorbar(distances, singleshot_water_fitpos2_mean[i], yerr=singleshot_water_fitpos2_std[i], fmt="--%ss"%colors[i])
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='center right', prop={'size':6})
	canvas.set_ylim([1.8, 3.2])
	
	canvas = fig.add_subplot(132)
	canvas.set_title("S1 / S2 mean centered")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_s1 - ref_s1.mean(), color='k', marker='o', label="Huang_resubmitted120910")
	P.plot(distances[1:], ref_s2 - ref_s2.mean(), color='k', linestyle='--', marker='o')
	for i in range(len(files)):
		P.errorbar(distances, N.array(singleshot_water_fitpos1_mean[i]) - N.array(singleshot_water_fitpos1_mean[i])[2:].mean(), yerr=singleshot_water_fitpos1_std[i], fmt="-%ss"%colors[i], label="%s" % (tags[i]))
		P.errorbar(distances, N.array(singleshot_water_fitpos2_mean[i]) - N.array(singleshot_water_fitpos2_mean[i])[2:].mean(), yerr=singleshot_water_fitpos2_std[i], fmt="--%ss"%colors[i])
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper right', prop={'size':6})
	canvas.set_ylim([-0.06, 0.1])
	
	canvas = fig.add_subplot(133)
	canvas.set_title("deltaQ")
	P.xlabel("distance (mm)")
	P.ylabel("Q (A-1)")
	P.plot(distances[1:], ref_dq, color='k', marker='o', label="Huang_resubmitted120910")
	for i in range(len(files)):
		P.errorbar(distances, singleshot_water_fitdeltaq_mean[i], yerr=singleshot_water_fitdeltaq_std[i], fmt="-%ss"%colors[i], label="%s" % (tags[i]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper left', prop={'size':6})
	canvas.set_ylim([0.95, 1.25])
	
	P.savefig(original_dir + "peakfit-dod_single-shot_mean_vs_dist.png")
	print "peakfit-dod_single-shot_mean_vs_dist.png saved."
	P.savefig(original_dir + "peakfit-dod_single-shot_mean_vs_dist.eps", format='eps')
	print "peakfit-dod_single-shot_mean_vs_dist.eps saved."
	#P.show()
	P.close()
	
	
	#for i in range(len(distances)):
		#txttag = "peakfit-dod_%.1fmm.txt"%(distances[i])
		#N.array([water_fitint1[i], water_fitpos1[i], water_fitfwhm1[i], water_fitint2[i], water_fitpos2[i], water_fitfwhm2[i]]).tofile(txttag, sep = "\n", format="%lf")
		#print "%s saved."%(txttag)
		#csvtag = "peakfit-dod_%.1fmm.csv"%(distances[i])
		#N.savetxt(csvtag, N.array([water_fitint1[i], water_fitpos1[i], water_fitfwhm1[i], water_fitint2[i], water_fitpos2[i], water_fitfwhm2[i]]), delimiter=",")
		#print "%s saved."%(csvtag)
