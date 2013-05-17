#!/usr/bin/env python

# Written by J. Sellberg on 2013-04-08
# Wrapper or sortHits for the new data from Jan 2013
# using the HN130110-5 nozzle at 190-200 PSI N2 (gas), 480-490 PSI He (liquid), driven at 900 kHz (20 Vpp)

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
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--remote", action="store_true", dest="remote", help="saves the output text files remotely in $HOME", default=False)
(options, args) = parser.parse_args()

# regular
runs = [[13, 14, 15, 16], [20, 21, 22], [23, 24, 25], [27, 28, 29]]
colors = ['r','g','b','c','m','y','k']
temperatures = [298,243,238,235] # 8.7 um droplets, 30.5 m/s, gamma = 0.8, 10 mm delay of cooling
distances = [0.400029950159991,25.8792070078867,36.0411564485204,46.0943296380132] # FINAL distances

# SCRATCH
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
#sorting_dir = "/reg/d/psdm/cxi/cxi74613/scratch/iceFinderCampaign/"
# RES & FTC
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/res/cleaned_hdf5/"
sorting_dir = "/reg/d/psdm/cxi/cxi74613/res/iceFinderCampaign/"

original_dir = os.getcwd() + '/'

# create hit lists for each type for each run
for i in N.arange(len(runs)):
	for j in N.arange(len(runs[i])):
		if (runs[i][j] < 100):
			run_tag = "r00%s"%(runs[i][j])
		else:
			run_tag = "r0%s"%(runs[i][j])
		run_file =  run_tag + '_strong_hits.txt'
		os.system("./sortHits " + run_file)

