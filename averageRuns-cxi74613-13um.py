#!/usr/bin/env python

# Written by J. Sellberg on 2013-04-01
# Has the same functionality as averageRuns-aerojet.py but for the new data from Jan 2013
# using the CXI Feb2011-1 nozzle at 300 PSI N2 (gas), 175 PSI He (liquid), driven at 200 kHz (20 Vpp)

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

# regular
runs = [[90, 91, 92, 93, 95, 96, 97, 98], [100, 101, 102, 103, 104], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
#nhits_water = [[166, 279, 178, 1582, 628, 1543, 1234, 1238], [452, 1054, 1041, 1487, 2542], [412, 438, 364, 100, 171, 96, 225, 83, 48, 252, 169, 323, 310, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
#nhits_ice = [[0, 0, 0, 3, 4, 0, 0, 3], [1, 2, 0, 3, 3], [0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
oldnevents = [[20111, 30633, 26921, 197046, 63811, 179230, 133414, 115540], [134509, 120648, 121893, 174011, 310714], [73828, 67206, 92611, 19743, 32708, 18402, 39955, 11670, 6798, 37989, 28012, 48348, 46672, 54121, 111866, 17966, 12967, 29348, 42683, 33433, 38800, 2156, 107000, 4375, 57549], [2928, 58383, 103888, 75020, 94569, 61638, 113293], [14000, 2702, 139053, 105357, 19294, 17119, 65695, 36383]]
# hits with failed peak fitting (2013-02-26, p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2])
#failedFits20_water = [[90, 91, 92, 93, 95, 96, 97, 98], [100, 101, 102, 103, 104], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
# thresholded hits below 50 ADUs
#t50hits_water = [[16, 49, 22, 201, 83, 219, 156, 128], [434, 106, 103, 154, 261], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
#t50hits_ice = [[0, 0, 0, 1, 0, 0, 0, 2], [1, 1, 0, 1, 2], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
# hits with failed peak fitting (2013-03-04, p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2])
#failedFits50_water = [[90, 91, 92, 93, 95, 96, 97, 98], [100, 101, 102, 103, 104], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
# thresholded hits below 100 ADUs
#t100hits_water = [[10, 24, 10, 147, 63, 119, 94, 89], [18, 101, 102, 103, 104], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
#t100hits_ice = [[0, 0, 0, 0, 3, 0, 0, 1], [0, 101, 102, 103, 104], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
# hits with failed peak fitting (2013-02-28, p0 = [1.1E9, 1.83, 0.25, 8.5E8, 2.98, 0.2])
#failedFits100_water = [[90, 91, 92, 93, 95, 96, 97, 98], [100, 101, 102, 103, 104], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
#colors = ['b','g','r','c','m','y','k']
colors = ['r','g','b','c','m','y','k']
temperatures = [230,228,227,225,224] # 12.8 um droplets, 5.5 m/s, gamma = 0.8, 10 mm delay of cooling
distances = [30.3599791314158,35.3687000607684,40.3787409506768,45.3885869580407,50.3988203584402] # FINAL distances
thresholds = [20, 50, 100]

# hexagonal ice peaks
HIceQ = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.222, '112':3.272, '201':3.324}
HIcePos = {'100':10.5, '002':9.5, '101':8.5, '102':6.2, '110':6.2, '103':6.2, '200':7., '112':6., '201':5.}
HIceQLabel = {'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.324, '112':3.374, '201':3.426}
# for imaging class
HIceP = {'100':9.5, '002':8.5, '101':7.5, '102':6.2, '110':5.7, '103':5.2, '200':7., '112':6., '201':5.}

# SCRATCH
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
#sorting_dir = "/reg/d/psdm/cxi/cxi74613/scratch/iceFinderCampaign/"
# RES & FTC
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/res/cleaned_hdf5/"
sorting_dir = "/reg/d/psdm/cxi/cxi74613/res/iceFinderCampaign/"

original_dir = os.getcwd() + '/'
run_tag = "r0144"

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

# statistics for hit rates
nhits_water = [[[] for j in runs] for i in thresholds]
nhits_ice = [[[] for j in runs] for i in thresholds]
nhits_tot = [[[] for j in runs] for i in thresholds]
sumhits_water = [[] for i in thresholds]
sumhits_ice = [[] for i in thresholds]
sumhits_tot = [[] for i in thresholds]
nevents = [[] for i in runs]
ratio20 = [[] for i in runs]
ratio50 = [[] for i in runs]
ratio100 = [[] for i in runs]
hitrate20 = [[] for i in runs]
hitrate50 = [[] for i in runs]
hitrate100 = [[] for i in runs]
ratios = [[] for i in thresholds]
ratio_deviations = [[] for i in thresholds]
hitrates = [[] for i in thresholds]
hitrate_deviations = [[] for i in thresholds]
len_runs = [len(i) for i in runs]

#fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
fig = P.figure(num=None, figsize=(25, 5), dpi=100, facecolor='w', edgecolor='k')
#canvas = fig.add_subplot(121)
canvas = fig.add_subplot(131)
canvas.set_title("Thresholded Ice Ratios", fontsize='x-large')
P.xlabel("Sample-Nozzle Distance (mm)", fontsize='x-large')
P.ylabel("Ice Ratio", fontsize='x-large')

for i in N.arange(len(runs)):
	for j in N.arange(len(runs[i])):
		if (runs[i][j] < 100):
			run_tag = "r00%s"%(runs[i][j])
		else:
			run_tag = "r0%s"%(runs[i][j])
		if options.peakfit:
			nhits_water[0][i].append(eDD.get_type(run_tag, 0) - eDD.get_failedFits_from_type(run_tag, 0))
			nhits_water[1][i].append(eDD.get_type(run_tag, 0) - eDD.get_type_below(run_tag, 0, 50) - eDD.get_failedFits_from_type_above(run_tag, 0, 50))
			nhits_water[2][i].append(eDD.get_type(run_tag, 0) - eDD.get_type_below(run_tag, 0, 50) - eDD.get_type_above_and_below(run_tag, 0, 50, 100) - eDD.get_failedFits_from_type_above(run_tag, 0, 100))
		else:
			nhits_water[0][i].append(eDD.get_type(run_tag, 0))
			nhits_water[1][i].append(eDD.get_type(run_tag, 0) - eDD.get_type_below(run_tag, 0, 50))
			nhits_water[2][i].append(eDD.get_type(run_tag, 0) - eDD.get_type_below(run_tag, 0, 50) - eDD.get_type_above_and_below(run_tag, 0, 50, 100))
		nhits_ice[0][i].append(eDD.get_type(run_tag, 1))
		nhits_ice[1][i].append(eDD.get_type(run_tag, 1) - eDD.get_type_below(run_tag, 1, 50))
		nhits_ice[2][i].append(eDD.get_type(run_tag, 1) - eDD.get_type_below(run_tag, 1, 50) - eDD.get_type_above_and_below(run_tag, 1, 50, 100))
		nhits_tot[0][i].append(float(nhits_water[0][i][j] + nhits_ice[0][i][j]))
		nhits_tot[1][i].append(float(nhits_water[1][i][j] + nhits_ice[1][i][j]))
		nhits_tot[2][i].append(float(nhits_water[2][i][j] + nhits_ice[2][i][j]))
		nevents[i].append(eDD.get_events(run_tag))
		ratio20[i].append(eDD.get_type(run_tag, 1)/nhits_tot[0][i][j])
		hitrate20[i].append(nhits_tot[0][i][j]/nevents[i][j])
		if nhits_tot[1][i][j]:
			ratio50[i].append((eDD.get_type(run_tag, 1) - eDD.get_type_below(run_tag, 1, 50))/nhits_tot[1][i][j])
		else:
			ratio50[i].append(0)
		hitrate50[i].append(nhits_tot[1][i][j]/nevents[i][j])
		if nhits_tot[2][i][j]:
			ratio100[i].append((eDD.get_type(run_tag, 1) - eDD.get_type_below(run_tag, 1, 50) - eDD.get_type_above_and_below(run_tag, 1, 50, 100))/nhits_tot[2][i][j])
		else:
			ratio100[i].append(0)
		hitrate100[i].append(nhits_tot[2][i][j]/nevents[i][j])
	
	ratios[0].append(N.mean(ratio20[i]))
	ratios[1].append(N.mean(ratio50[i]))
	ratios[2].append(N.mean(ratio100[i]))
	ratio_deviations[0].append(N.std(ratio20[i]))
	ratio_deviations[1].append(N.std(ratio50[i]))
	ratio_deviations[2].append(N.std(ratio100[i]))
	hitrates[0].append(sum(nhits_tot[0][i])/sum(nevents[i]))
	hitrates[1].append(sum(nhits_tot[1][i])/sum(nevents[i]))
	hitrates[2].append(sum(nhits_tot[2][i])/sum(nevents[i]))
	hitrate_deviations[0].append(N.std(hitrate20[i]))
	hitrate_deviations[1].append(N.std(hitrate50[i]))
	hitrate_deviations[2].append(N.std(hitrate100[i]))
	sumhits_water[0].append(float(sum(nhits_water[0][i])))
	sumhits_water[1].append(float(sum(nhits_water[1][i])))
	sumhits_water[2].append(float(sum(nhits_water[2][i])))
	sumhits_ice[0].append(float(sum(nhits_ice[0][i])))
	sumhits_ice[1].append(float(sum(nhits_ice[1][i])))
	sumhits_ice[2].append(float(sum(nhits_ice[2][i])))
	sumhits_tot[0].append(sum(nhits_tot[0][i]))
	sumhits_tot[1].append(sum(nhits_tot[1][i]))
	sumhits_tot[2].append(sum(nhits_tot[2][i]))
	if (sumhits_water[0][i] < 0 or sumhits_ice[0][i] < 0 or sumhits_tot[0][i] < 0):
		print "runs:", runs[i]
		print "20 ADUs: sumhits_water[%d] = %d" % (i, sumhits_water[0][i])
		print "20 ADUs: sumhits_ice[%d] = %d" % (i, sumhits_ice[0][i])
		print "20 ADUs: sumhits_tot[%d] = %d" % (i, sumhits_tot[0][i])
	if (sumhits_water[1][i] < 0 or sumhits_ice[1][i] < 0 or sumhits_tot[1][i] < 0):
		print "runs:", runs[i]
		print "50 ADUs: sumhits_water[%d] = %d" % (i, sumhits_water[1][i])
		print "50 ADUs: sumhits_ice[%d] = %d" % (i, sumhits_ice[1][i])
		print "50 ADUs: sumhits_tot[%d] = %d" % (i, sumhits_tot[1][i])
	if (sumhits_water[2][i] < 0 or sumhits_ice[2][i] < 0 or sumhits_tot[2][i] < 0):
		print "runs:", runs[i]
		print "100 ADUs: sumhits_water[%d] = %d" % (i, sumhits_water[2][i])
		print "100 ADUs: sumhits_ice[%d] = %d" % (i, sumhits_ice[2][i])
		print "100 ADUs: sumhits_tot[%d] = %d" % (i, sumhits_tot[2][i])
	
	P.scatter(distances[i], ratios[0][i], color='r', marker='o')
	P.scatter(distances[i], ratios[1][i], color='g', marker='o')
	P.scatter(distances[i], ratios[2][i], color='b', marker='o')

for i in N.arange(3):
	P.plot(distances, ratios[i], color=colors[i], label="%s ADUs"%(thresholds[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left')

#canvas = fig.add_subplot(122)
canvas = fig.add_subplot(132)
canvas.set_title("Standard Deviation of Runs", fontsize='x-large')
P.xlabel("Sample-Nozzle Distance (mm)", fontsize='x-large')
P.ylabel("Standard Deviation", fontsize='x-large')

for i in N.arange(3):
	P.plot(distances, ratio_deviations[i], color=colors[i], label="%s ADUs"%(thresholds[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left')

canvas = fig.add_subplot(133)
canvas.set_title("Number of Hits", fontsize='x-large')
P.xlabel("Sample-Nozzle Distance (mm)", fontsize='x-large')
P.ylabel("Number of Hits", fontsize='x-large')

for i in N.arange(3):
	P.plot(distances, sumhits_water[i], color=colors[i], label="Water %s ADUs"%(thresholds[i]))
	P.plot(distances, sumhits_ice[i], color=colors[i+3], label="Ice %s ADUs"%(thresholds[i]))
	P.plot(distances, sumhits_tot[i], color=colors[i], label="Total %s ADUs"%(thresholds[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, bbox_to_anchor = (1, 1), loc='upper left')

print "output_runs-13um-ice_ratios.png saved."
P.savefig(original_dir + "output_runs-13um-ice_ratios.png")
#P.savefig(original_dir + "output_runs-13um-ice_ratios.png", format='eps')
#P.show()
P.close()

# SECOND PLOT
fig = P.figure(num=None, figsize=(8, 5), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(111)
canvas.set_title("Number of Hits", fontsize='x-large')
P.xlabel("Sample-Nozzle Distance (mm)", fontsize='x-large')
P.ylabel("Number of Hits", fontsize='x-large')

for i in N.arange(3):
	P.fill_between(distances, sumhits_water[i], color=colors[i])
	P.plot(distances, sumhits_ice[i], color=colors[i+3], label="Ice %s ADUs"%(thresholds[i]))
	P.plot(distances, sumhits_tot[i], color=colors[i], label="Total %s ADUs"%(thresholds[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper right')

print "output_runs-13um-ice_hits.png saved."
P.savefig(original_dir + "output_runs-13um-ice_hits.png")
#P.savefig(original_dir + "output_runs-13um-ice_hits.png", format='eps')
#P.show()
P.close()

# THIRD PLOT
fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(121)
canvas.set_title("Hitrates", fontsize='x-large')
P.xlabel("Sample-Nozzle Distance (mm)", fontsize='x-large')
P.ylabel("Hitrate", fontsize='x-large')

for i in N.arange(3):
	P.plot(distances, hitrates[i], color=colors[i], label="%s ADUs"%(thresholds[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper right')

canvas = fig.add_subplot(122)
canvas.set_title("Hitrate Deviations of Runs", fontsize='x-large')
P.xlabel("Sample-Nozzle Distance (mm)", fontsize='x-large')
P.ylabel("Standard Deviation", fontsize='x-large')

for i in N.arange(3):
	P.plot(distances, hitrate_deviations[i], color=colors[i], label="%s ADUs"%(thresholds[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper right')

print "output_runs-13um-ice_hitrates.png saved."
P.savefig(original_dir + "output_runs-13um-ice_hitrates.png")
#P.savefig(original_dir + "output_runs-13um-ice_ratios.png", format='eps')
#P.show()
P.close()

for i in N.arange(3):
	txttag = "output_runs-13um-ice_ratios-%sADUs.txt"%(thresholds[i])
	N.array([distances, len_runs, sumhits_tot[i], sumhits_water[i], sumhits_ice[i], ratios[i], ratio_deviations[i], hitrates[i], hitrate_deviations[i]]).tofile(txttag, sep = "\n", format="%lf")
	print "%s saved."%(txttag)

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
ice_pattern = []
ice_correlation = []
ice_angavg = []
ice_angavgQ = []
ice_pattern_shape = False
ice_correlation_shape = False
ice_angavg_shape = False
ice_angavgQ_shape = N.array([False])
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
	excludedIce_pattern = []
	excludedIce_correlation = []
	excludedIce_angavg = []
	excludedIce_angavgQ = []

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
	temp_water_correlation = []
	temp_water_angavg = []
	temp_water_angavgQ = []
	temp_water_fitint1 = []
	temp_water_fitpos1 = []
	temp_water_fitfwhm1 = []
	temp_water_fitint2 = []
	temp_water_fitpos2 = []
	temp_water_fitfwhm2 = []
	temp_ice_pattern = []
	temp_ice_correlation = []
	temp_ice_angavg = []
	temp_ice_angavgQ = []
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
		excludedTemp_ice_pattern = []
		excludedTemp_ice_correlation = []
		excludedTemp_ice_angavg = []
		excludedTemp_ice_angavgQ = []
	for j in N.arange(len(runs[i])):
		if (runs[i][j] < 100):
			run_tag = "r00%s"%(runs[i][j])
		else:
			run_tag = "r0%s"%(runs[i][j])
		run_dir =  'output_' + run_tag + '/'
		#water
		if (eDD.get_type(run_tag, 0) != 0):
			if options.exclude:
				typeTag = run_tag + "_type0-" + options.excludeFile
			else:
				typeTag = run_tag + "_type0"
			if os.path.exists(sorting_dir + run_dir + typeTag + '.h5'):
				print 'found: ' + sorting_dir + run_dir + typeTag + '.h5'
				f = H.File(sorting_dir + run_dir + typeTag + '.h5', 'r')
				if (nhits_water[thresholdIndex][i][j] > 0):
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
					print "No water hits for r0%s, padding zeros." % (runs[i][j])
				else:
					print "No water hits for r0%s and shape is unknown, aborting." % (runs[i][j])
					sys.exit(1)
				#excluded hits
				if (options.exclude and options.saveExcluded):
					if ((nhits_water[0][i][j] - nhits_water[thresholdIndex][i][j]) > 0):
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
						print "No excluded water hits for r0%s, padding zeros." % (runs[i][j])
					else:
						print "No excluded water hits for r0%s and shape is unknown, aborting." % (runs[i][j])
						sys.exit(1)
				
				f.close()
			else:
				print 'file missing: ' + sorting_dir + run_dir + typeTag + '.h5'
				sys.exit(1)
		elif (water_pattern_shape and (water_correlation_shape or not options.xaca) and water_angavg_shape and water_angavgQ_shape.any()):
			temp_water_pattern.append(N.zeros(water_pattern_shape))
			if options.xaca:
				temp_water_correlation.append(N.zeros(water_correlation_shape))
			temp_water_angavg.append(N.zeros(water_angavg_shape))
			temp_water_angavgQ.append(water_angavgQ_shape)
			print "No water hits for r0%s, padding zeros." % (runs[i][j])
			if (options.exclude and options.saveExcluded):
				excludedTemp_water_pattern.append(N.zeros(water_pattern_shape))
				if options.xaca:
					excludedTemp_water_correlation.append(N.zeros(water_correlation_shape))
				excludedTemp_water_angavg.append(N.zeros(water_angavg_shape))
				excludedTemp_water_angavgQ.append(water_angavgQ_shape)
				print "No excluded water hits for r0%s, padding zeros." % (runs[i][j])
		else:
			temp_water_pattern.append(N.zeros([1764, 1764])) #cxi74613
			if options.xaca:
				temp_water_correlation.append(N.zeros([128, 512])) #cxi74613
			temp_water_angavg.append(N.zeros([3491])) #cxi74613
			temp_water_angavgQ.append(N.arange(0.09,3.581,0.001)) #cxi74613
			print "No excluded water hits for r0%s, padding zeros with standard shape." % (runs[i][j])
			if (options.exclude and options.saveExcluded):
				excludedTemp_water_pattern.append(N.zeros([1764, 1764])) #cxi74613
				if options.xaca:
					excludedTemp_water_correlation.append(N.zeros([128, 512])) #cxi74613
				excludedTemp_water_angavg.append(N.zeros([3491])) #cxi74613
				excludedTemp_water_angavgQ.append(N.arange(0.09,3.581,0.001)) #cxi74613
				print "No excluded water hits for r0%s, padding zeros with standard shape." % (runs[i][j])
		
		#ice
		if (eDD.get_type(run_tag, 1) != 0):
			if options.exclude:
				typeTag = run_tag + "_type1-" + options.excludeFile
			else:
				typeTag = run_tag + "_type1"
			if os.path.exists(sorting_dir + run_dir + 'type1/' + typeTag + '.h5'):
				print 'found: ' + sorting_dir + run_dir + 'type1/' + typeTag + '.h5'
				f = H.File(sorting_dir + run_dir + 'type1/' + typeTag + '.h5', 'r')
				if (nhits_ice[thresholdIndex][i][j] > 0):
					temp_ice_pattern.append(N.array(f['data']['diffraction']))
					ice_pattern_shape = N.array(f['data']['diffraction']).shape
					if options.xaca:
						temp_ice_correlation.append(N.array(f['data']['correlation']))
						ice_correlation_shape = N.array(f['data']['correlation']).shape
					temp_ice_angavg.append(N.array(f['data']['angavg']))
					ice_angavg_shape = N.array(f['data']['angavg']).shape
					temp_ice_angavgQ.append(N.array(f['data']['angavgQ']))
					ice_angavgQ_shape = N.array(f['data']['angavgQ'])
				elif (ice_pattern_shape and (ice_correlation_shape or not options.xaca) and ice_angavg_shape and ice_angavgQ_shape.any()):
					temp_ice_pattern.append(N.zeros(ice_pattern_shape))
					if options.xaca:
						temp_ice_correlation.append(N.zeros(ice_correlation_shape))
					temp_ice_angavg.append(N.zeros(ice_angavg_shape))
					temp_ice_angavgQ.append(ice_angavgQ_shape)
					print "No ice hits for r0%s, padding zeros." % (runs[i][j])
				else:
					print "No ice hits for r0%s and shape is unknown, aborting." % (runs[i][j])
					sys.exit(1)
				#excluded hits
				if (options.exclude and options.saveExcluded): 
					if ((nhits_ice[0][i][j] - nhits_ice[thresholdIndex][i][j]) > 0):
						excludedTemp_ice_pattern.append(N.array(f['data']['excludedHits']['diffraction']))
						if options.xaca:
							excludedTemp_ice_correlation.append(N.array(f['data']['excludedHits']['correlation']))
						excludedTemp_ice_angavg.append(N.array(f['data']['excludedHits']['angavg']))
						excludedTemp_ice_angavgQ.append(N.array(f['data']['excludedHits']['angavgQ']))
					elif (ice_pattern_shape and (ice_correlation_shape or not options.xaca) and ice_angavg_shape and ice_angavgQ_shape.any()):
						excludedTemp_ice_pattern.append(N.zeros(ice_pattern_shape))
						if options.xaca:
							excludedTemp_ice_correlation.append(N.zeros(ice_correlation_shape))
						excludedTemp_ice_angavg.append(N.zeros(ice_angavg_shape))
						excludedTemp_ice_angavgQ.append(ice_angavgQ_shape)
						print "No excluded ice hits for r0%s, padding zeros." % (runs[i][j])
					else:
						print "No excluded ice hits for r0%s and shape is unknown, aborting." % (runs[i][j])
						sys.exit(1)
				f.close()
			else:
				print 'file missing: ' + sorting_dir + run_dir + 'type1/' + typeTag + '.h5'
				sys.exit(1)
		elif (ice_pattern_shape and (ice_correlation_shape or not options.xaca) and ice_angavg_shape and ice_angavgQ_shape.any()):
			temp_ice_pattern.append(N.zeros(ice_pattern_shape))
			if options.xaca:
				temp_ice_correlation.append(N.zeros(ice_correlation_shape))
			temp_ice_angavg.append(N.zeros(ice_angavg_shape))
			temp_ice_angavgQ.append(ice_angavgQ_shape)
			print "No ice hits for r0%s, padding zeros." % (runs[i][j])
			if (options.exclude and options.saveExcluded):
				excludedTemp_ice_pattern.append(N.zeros(ice_pattern_shape))
				if options.xaca:
					excludedTemp_ice_correlation.append(N.zeros(ice_correlation_shape))
				excludedTemp_ice_angavg.append(N.zeros(ice_angavg_shape))
				excludedTemp_ice_angavgQ.append(ice_angavgQ_shape)
				print "No excluded ice hits for r0%s, padding zeros." % (runs[i][j])
		else:
			temp_ice_pattern.append(N.zeros([1764, 1764])) #cxi74613
			if options.xaca:
				temp_ice_correlation.append(N.zeros([128, 512])) #cxi74613
			temp_ice_angavg.append(N.zeros([3491])) #cxi74613
			temp_ice_angavgQ.append(N.arange(0.09,3.581,0.001)) #cxi74613
			print "No excluded ice hits for r0%s, padding zeros with standard shape." % (runs[i][j])
			if (options.exclude and options.saveExcluded):
				excludedTemp_ice_pattern.append(N.zeros([1764, 1764])) #cxi74613
				if options.xaca:
					excludedTemp_ice_correlation.append(N.zeros([128, 512])) #cxi74613
				excludedTemp_ice_angavg.append(N.zeros([3491])) #cxi74613
				excludedTemp_ice_angavgQ.append(N.arange(0.09,3.581,0.001)) #cxi74613
				print "No excluded ice hits for r0%s, padding zeros with standard shape." % (runs[i][j])
	
	
	#plot temp_angavg
	fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
	canvas = fig.add_subplot(121)
	canvas.set_title("Water, T = %s K"%(temperatures[i]))
	P.xlabel("Q (Angstroms-1)")
	P.ylabel("Average Intensity (ADUs)")
	for j in N.arange(len(runs[i])):
		P.plot(temp_water_angavgQ[j], temp_water_angavg[j], color=colors[j % len(colors)], label="r0%s"%(runs[i][j]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels)
	canvas = fig.add_subplot(122)
	canvas.set_title("Ice, T = %s K"%(temperatures[i]))
	P.xlabel("Q (Angstroms-1)")
	P.ylabel("Average Intensity (ADUs)")
	for j in N.arange(len(runs[i])):
		P.plot(temp_ice_angavgQ[j], temp_ice_angavg[j], color=colors[j % len(colors)], label="r0%s"%(runs[i][j]))
	
	handles, labels = canvas.get_legend_handles_labels()
	canvas.legend(handles, labels, loc='upper left')
	
	maxAngAvg = max([temp_ice_angavg[j].max() for j in N.arange(len(temp_ice_angavg))])
	numQLabels = len(HIceQ.keys())+1
	for k,j in HIceQ.iteritems():
		labelPosition = HIcePos[k]*maxAngAvg/numQLabels
		P.axvline(j, 0, maxAngAvg, color='k', ls='--')
		P.text(HIceQLabel[k], labelPosition, str(k), rotation="45")
	
	print "output_runs-T%sK-13um-%dADUs-angavg.png saved."%(temperatures[i],options.threshold)
	P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-angavg.png"%(temperatures[i],options.threshold))
	#P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-angavg.eps"%(temperatures[i],options.threshold), format='eps')
	#P.show()
	P.close()
	
	sumwater = float(sum(nhits_water[thresholdIndex][i]))
	sumice = float(sum(nhits_ice[thresholdIndex][i]))
	if (options.exclude and options.saveExcluded):
		excludesumwater = float(sum(nhits_water[0][i]) - sum(nhits_water[thresholdIndex][i]))
		excludesumice = float(sum(nhits_ice[0][i]) - sum(nhits_ice[thresholdIndex][i]))
	
	for j in N.arange(len(runs[i])):
		nwater = nhits_water[thresholdIndex][i][j]
		nice = nhits_ice[thresholdIndex][i][j]
		if (options.exclude and options.saveExcluded):
			excludenwater = nhits_water[0][i][j] - nhits_water[thresholdIndex][i][j]
			excludenice = nhits_ice[0][i][j] - nhits_ice[thresholdIndex][i][j]
		
		if (sumwater > 0):
			temp_water_angavg[j] *= nwater/sumwater
			temp_water_pattern[j] *= nwater/sumwater
			if options.xaca:
				temp_water_correlation[j] *= nwater/sumwater
		if (sumice > 0):
			temp_ice_angavg[j] *= nice/sumice
			temp_ice_pattern[j] *= nice/sumice
			if options.xaca:
				temp_ice_correlation[j] *= nice/sumice
		if (options.exclude and options.saveExcluded):
			if (excludesumwater > 0):
				excludedTemp_water_angavg[j] *= excludenwater/excludesumwater
				excludedTemp_water_pattern[j] *= excludenwater/excludesumwater
				if options.xaca:
					excludedTemp_water_correlation[j] *= excludenwater/excludesumwater
			if (excludesumice > 0):
				excludedTemp_ice_angavg[j] *= excludenice/excludesumice
				excludedTemp_ice_pattern[j] *= excludenice/excludesumice
				if options.xaca:
					excludedTemp_ice_correlation[j] *= excludenice/excludesumice
	
	water_angavg.append(N.array(temp_water_angavg).sum(axis=0))
	water_angavgQ.append(N.array(temp_water_angavgQ).mean(axis=0))
	water_pattern.append(N.array(temp_water_pattern).sum(axis=0))
	if options.xaca:
		water_correlation.append(N.array(temp_water_correlation).sum(axis=0))
	ice_angavg.append(N.array(temp_ice_angavg).sum(axis=0))
	ice_angavgQ.append(N.array(temp_ice_angavgQ).mean(axis=0))
	ice_pattern.append(N.array(temp_ice_pattern).sum(axis=0))
	if options.xaca:
		ice_correlation.append(N.array(temp_ice_correlation).sum(axis=0))	
	if (options.exclude and options.saveExcluded):
		excludedWater_angavg.append(N.array(excludedTemp_water_angavg).sum(axis=0))
		excludedWater_angavgQ.append(N.array(excludedTemp_water_angavgQ).mean(axis=0))
		excludedWater_pattern.append(N.array(excludedTemp_water_pattern).sum(axis=0))
		if options.xaca:
			excludedWater_correlation.append(N.array(excludedTemp_water_correlation).sum(axis=0))
		excludedIce_angavg.append(N.array(excludedTemp_ice_angavg).sum(axis=0))
		excludedIce_angavgQ.append(N.array(excludedTemp_ice_angavgQ).mean(axis=0))
		excludedIce_pattern.append(N.array(excludedTemp_ice_pattern).sum(axis=0))
		if options.xaca:
			excludedIce_correlation.append(N.array(excludedTemp_ice_correlation).sum(axis=0))	
	
	currImg = img_class(water_pattern[i], ice_pattern[i], [int(sumwater), int(sumice)], "output_runs-T%sK-13um-%dADUs-pattern"%(temperatures[i],options.threshold))
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
	
		print "output_runs-T%sK-13um-%dADUs-peakfit_hist.png saved."%(temperatures[i],options.threshold)
		P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-peakfit_hist.png"%(temperatures[i],options.threshold))
		#P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-peakfit_hist.eps"%(temperatures[i],options.threshold), format='eps')
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
		
		print "output_runs-T%sK-13um-%dADUs-peakfit_corr-s1.png saved."%(temperatures[i],options.threshold)
		P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-peakfit_corr-s1.png"%(temperatures[i],options.threshold))
		#P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-peakfit_corr-s1.eps"%(temperatures[i],options.threshold), format='eps')
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
		
		print "output_runs-T%sK-13um-%dADUs-peakfit_corr-s2.png saved."%(temperatures[i],options.threshold)
		P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-peakfit_corr-s2.png"%(temperatures[i],options.threshold))
		#P.savefig(original_dir + "output_runs-T%sK-13um-%dADUs-peakfit_corr-s2.eps"%(temperatures[i],options.threshold), format='eps')
		#P.show()
		P.close()


#plot angavg
fig = P.figure(num=None, figsize=(13.5, 10), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(221)
#canvas.set_title("Water", fontsize='x-large')
#P.xlabel("Q (Angstroms-1)", fontsize='x-large')
#P.ylabel("Average Intensity (ADUs/srad)", fontsize='x-large')
canvas.set_title("Water")
P.xlabel("Q (Angstroms-1)")
P.ylabel("Average Intensity (ADUs/srad)")
for i in N.arange(len(runs)):
	#if i > 3:
		#P.plot(water_angavgQ[i], water_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))
	P.plot(water_angavgQ[i], water_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left', bbox_to_anchor = (1, 1))

canvas = fig.add_subplot(223)
#P.xlabel("Q (Angstroms-1)", fontsize='x-large')
#P.ylabel("Average Intensity (ADUs/srad)", fontsize='x-large')
canvas.set_title("Ice")
P.xlabel("Q (Angstroms-1)")
P.ylabel("Average Intensity (ADUs/srad)")
for i in N.arange(len(runs)):
	P.plot(ice_angavgQ[i], ice_angavg[i], color=colors[i-4], label="T = %s K"%(temperatures[i]))

handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper left', bbox_to_anchor = (1, 1))

maxAngAvg = max([ice_angavg[i].max() for i in N.arange(len(ice_angavg))])
numQLabels = len(HIceQ.keys())+1
for k,j in HIceQ.iteritems():
	labelPosition = HIcePos[k]*maxAngAvg/numQLabels
	P.axvline(j, 0, maxAngAvg, color='k', ls='--')
	P.text(HIceQLabel[k], labelPosition, str(k), rotation="45")

print "output_runs-gdvn-13um-all_T-angavg_Q_%dADUs.png saved." % options.threshold
P.savefig(original_dir + "output_runs-gdvn-13um-all_T-angavg_Q_%dADUs.png" % options.threshold)
#P.savefig(original_dir + "output_runs-gdvn-13um-all_T-angavg_Q_%dADUs.eps" % options.threshold, format='eps')
#P.show()
P.close()

#save to file
hdf5tag = "output_runs-gdvn-13um-all_T+Q_%dADUs.h5" % options.threshold
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
	entry_4 = f.create_group("/data/%.1fmm/ice"%(distances[i]))
	entry_4.create_dataset("diffraction", data=ice_pattern[i])
	if options.xaca:
		entry_4.create_dataset("correlation", data=ice_correlation[i])
	entry_4.create_dataset("angavg", data=ice_angavg[i])
	entry_4.create_dataset("angavg_Q", data=ice_angavgQ[i])
	if options.peakfit:
		entry_5 = f.create_group("/data/%.1fmm/water/peakFit"%(distances[i]))
		entry_5.create_dataset("int1", data=water_fitint1[i])
		entry_5.create_dataset("int2", data=water_fitint2[i])
		entry_5.create_dataset("pos1", data=water_fitpos1[i])
		entry_5.create_dataset("pos2", data=water_fitpos2[i])
		entry_5.create_dataset("fwhm1", data=water_fitfwhm1[i])
		entry_5.create_dataset("fwhm2", data=water_fitfwhm2[i])
		entry_5.create_dataset("deltaQ", data=water_fitdeltaq[i])
	if (options.exclude and options.saveExcluded):
		entry_6 = f.create_group("/data/%.1fmm/water/excludedHits"%(distances[i]))
		entry_6.create_dataset("diffraction", data=excludedWater_pattern[i])
		if options.xaca:
			entry_6.create_dataset("correlation", data=excludedWater_correlation[i])
		entry_6.create_dataset("angavg", data=excludedWater_angavg[i])
		entry_6.create_dataset("angavg_Q", data=excludedWater_angavgQ[i])
		entry_7 = f.create_group("/data/%.1fmm/ice/excludedHits"%(distances[i]))
		entry_7.create_dataset("diffraction", data=excludedIce_pattern[i])
		if options.xaca:
			entry_7.create_dataset("correlation", data=excludedIce_correlation[i])
		entry_7.create_dataset("angavg", data=excludedIce_angavg[i])
		entry_7.create_dataset("angavg_Q", data=excludedIce_angavgQ[i])
		if options.peakfit:
			entry_8 = f.create_group("/data/%.1fmm/water/excludedHits/peakFit"%(distances[i]))
			entry_8.create_dataset("int1", data=excludedWater_fitint1[i])
			entry_8.create_dataset("int2", data=excludedWater_fitint2[i])
			entry_8.create_dataset("pos1", data=excludedWater_fitpos1[i])
			entry_8.create_dataset("pos2", data=excludedWater_fitpos2[i])
			entry_8.create_dataset("fwhm1", data=excludedWater_fitfwhm1[i])
			entry_8.create_dataset("fwhm2", data=excludedWater_fitfwhm2[i])
			entry_8.create_dataset("deltaQ", data=excludedWater_fitdeltaq[i])
f.close()
print "Successfully updated %s" % hdf5tag
