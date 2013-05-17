#!/usr/bin/env python

# Written by J. Sellberg on 2013-04-08
# Wrapper or makeHitLists and makeHitListsRemote for the new data from Jan 2013
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
parser.add_option("-r", "--remote", action="store_true", dest="remote", help="saves the output text files remotely in $HOME", default=False)
(options, args) = parser.parse_args()

# regular
runs = [[90, 91, 92, 93, 95, 96, 97, 98], [100, 101, 102, 103, 104], [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136], [137, 138, 139, 140, 141, 142, 143], [145, 146, 147, 148, 149, 150, 151, 152]]
colors = ['r','g','b','c','m','y','k']
temperatures = [230,228,227,225,224] # 12.8 um droplets, 5.5 m/s, gamma = 0.8, 10 mm delay of cooling
distances = [30.3599791314158,35.3687000607684,40.3787409506768,45.3885869580407,50.3988203584402] # FINAL distances

# SCRATCH
#source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
#sorting_dir = "/reg/d/psdm/cxi/cxi74613/scratch/iceFinderCampaign/"
# RES & FTC
source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
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
		run_dir =  'output_' + run_tag
		if options.remote:
			os.system("./makeHitListsRemote " + run_dir)
		else:
			os.system("./makeHitLists " + run_dir)

