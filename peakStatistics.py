#!/usr/bin/env python

import numpy as N
import h5py as H
import glob as G
import matplotlib
import matplotlib.pyplot as P
import scipy
import scipy.interpolate as I
import sys, os, re, shutil, subprocess, time
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--threshold", action="store", type="int", dest="threshold", help="sets threshold for max intensity of angular average below which hits are automatically sorted to sub-type (default:0)", default=20)
parser.add_option("-o", "--output", action="store", type="string", dest="output", help="output name (default: output)", metavar="OUTPUTNAME", default="output")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="prints out the frame number as it is processed", default=False)
(options, args) = parser.parse_args()

# INPUT DATA
#runs_aero = [[114],[121],[118],[123],[129,130,133],[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
#runs_dod = [[105],[111],[68],[63,64],[110],[73],[77,79],[88],[97]]
runs = [[144,145,146,147,151,169,170],[167,168,171,172,173],[165,166]]
#distances_aero = [0.668972,11.548956,12.349269,21.407288,30.603057,40.403027,45.595702,50.619378]
#distances_dod = [0.750000,10.011910,15.798895,20.789092,30.006612,30.738471,40.759723,60.710063,70.662700]
distances = [40.403027,45.595702,50.619378]
temperatures = [223,221,220]
colors = ['r','g','b','c','m','y','k']

# OPEN DATA
ice_name = [[] for i in runs]
ice_iavg = [[] for i in runs]
ice_npeaks = [[] for i in runs]

for i in N.arange(len(runs)):
	for j in N.arange(len(runs[i])):

		file_hits = "npeaks/r0%d-cleanedhits.txt" % (runs[i][j])
		file_ice = "hits_from_script/r0%d_hits-type1.txt" % (runs[i][j])
		file_below50 = "%s_r0%d/type1/below50ADUs.txt" % (options.output, runs[i][j])
		file_below100 = "%s_r0%d/type1/below100above50ADUs.txt" % (options.output, runs[i][j])

		fhits = N.loadtxt(file_hits, dtype='str', comments='#', delimiter=', ')
		fice = N.loadtxt(file_ice, dtype='str', comments='#', delimiter=' ')
		
		hits_below50 = []
		if os.path.exists(file_below50):
			fbelow50 = N.loadtxt(file_below50, dtype='str', comments='#', delimiter=' ')
			if (fbelow50.size == 1):
				fbelow50 = [fbelow50.base[0]]
			
			for hit in fbelow50:
				dirInStringPos = hit.find("-angavg")
				if (dirInStringPos != -1):
					hits_below50.append(hit[:dirInStringPos])
				else:
					hits_below50.append(hit)
		
		hits_below100 = []
		if os.path.exists(file_below100):		
			fbelow100 = N.loadtxt(file_below100, dtype='str', comments='#', delimiter=' ')
			if (fbelow100.size == 1):
				fbelow100 = [fbelow100.base[0]]

			for hit in fbelow100:
				dirInStringPos = hit.find("-angavg")
				if (dirInStringPos != -1):
					hits_below100.append(hit[:dirInStringPos])
				else:
					hits_below100.append(hit)
		
		hits = fhits.transpose()[0]
		#if os.path.exists(dirName+re.sub("-angavg",'',fname)):
		hits_name = []
		for hit in hits:
			dirInStringPos = hit.find("/")
			if (dirInStringPos != -1):
				hits_name.append(hit[dirInStringPos+1:])
			else:
				hits_name.append(hit)

		hits_iavg = N.array(fhits.transpose()[1], dtype='float')
		hits_npeaks = N.array(fhits.transpose()[2], dtype='int')		

		sHits = set(hits_name)
		sIce = set(fice)
		sBelow50 = set(hits_below50)
		sBelow100 = set(hits_below100)

		if (sIce.issubset(sHits)):
			if (sBelow50.issubset(sIce) and sBelow100.issubset(sIce)):
				counter = 0
				for hit in hits_name:
					if (options.threshold == 20):
						if (hit in sIce):
							ice_name[i].append(hits_name[counter])
							ice_iavg[i].append(hits_iavg[counter])
							ice_npeaks[i].append(hits_npeaks[counter])
					elif (options.threshold == 50):
						if (hit in sIce and hit not in sBelow50):
							ice_name[i].append(hits_name[counter])
							ice_iavg[i].append(hits_iavg[counter])
							ice_npeaks[i].append(hits_npeaks[counter])
					elif (options.threshold == 100):
						if (hit in sIce and hit not in sBelow50 and hit not in sBelow100):
							ice_name[i].append(hits_name[counter])
							ice_iavg[i].append(hits_iavg[counter])
							ice_npeaks[i].append(hits_npeaks[counter])
					else:
						print "Chosen threshold of %d ADUs has not been sorted, aborting." % (options.threshold)
						sys.exit(1)
					counter += 1
			else:
				print "Found %d ice hits for higher thresholds in r0%d that do not exist in %s, aborting. Rerun: ./makeHitLists" % (len(sBelow50-sIce)+len(sBelow100-sIce), runs[i][j], file_ice)
				sys.exit(1)
		else:
			print "Found %d ice hits in r0%d that do not exist in %s, aborting." % (len(sIce-sHits), runs[i][j], file_hits)
			sys.exit(1)


# STATISTICS
for j in N.arange(len(runs)):
	print "Statistics: d=%.1f mm, %d ADUs" % (distances[j], options.threshold)
	npeaks = N.array(ice_npeaks[j])
	nmean = N.mean(npeaks)
	nstd = N.std(npeaks)
	nmax = N.max(npeaks)
	print "\tmean = %f, std = %f, max = %f" % (nmean, nstd, nmax)

# PLOT DATA
fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(121)
canvas.set_title("Histogram - %d ADUs" % (options.threshold))
P.xlabel("Npeaks")
P.ylabel("Hist(Npeaks)")
for j in N.arange(len(runs)):
	hist_bins = N.arange(max(ice_npeaks[j])+3) - 0.5
	ice_hist, hist_bins = N.histogram(ice_npeaks[j], bins=hist_bins)
	ice_bins = [(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(ice_hist))]
	P.plot(ice_bins, ice_hist, color=colors[j], label="%.1fmm" % (distances[j]))

P.xlim([0, 100])
handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper right')

canvas = fig.add_subplot(122)
canvas.set_title("Probability - %d ADUs" % (options.threshold))
P.xlabel("Npeaks")
P.ylabel("P(Npeaks)")
for j in N.arange(len(runs)):
	hist_bins = N.arange(max(ice_npeaks[j])+3) - 0.5
	ice_hist, hist_bins = N.histogram(ice_npeaks[j], bins=hist_bins, normed=True)
	ice_bins = [(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(ice_hist))]
	P.plot(ice_bins, ice_hist, color=colors[j], label="%.1fmm" % (distances[j]))

P.xlim([0, 100])
handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper right')

pngtag = "peakStatistics-%dADUs.png" % (options.threshold)
P.savefig(pngtag)
print "%s saved."%(pngtag)
#epstag = "peakStatistics-%dADUs.eps" % (options.threshold)
#P.savefig(epstag, format='eps')
#print "%s saved."%(epstag)
P.show()

# SAVE HISTOGRAMS TO TEXTFILE NOT YET IMPLEMENTED SINCE IT HASN'T BEEN NECESSARY
#for i in N.arange(ndod):
#	txttag = "sample-nozzle_distances-droplet_dispenser_%s.txt"%(labels_dod[i])
#	N.array([mean_dod[i],std_dod[i],maxmin_dod[i],error_dod[i]]).tofile(txttag, sep = "\n", format="%lf")
#	print "%s saved."%(txttag)
#	#N.array(["Mean (mm)",mean_dod[i],"Standard Dev. (mm)",std_dod[i],"Max-Min (mm)",maxmin_dod[i],"Error (%)",maxmin_dod[i]]).tofile("sample-nozzle_distances-droplet_dispenser_%s.txt"%(labels_dod[i]), sep = "\n", format="%lf")
