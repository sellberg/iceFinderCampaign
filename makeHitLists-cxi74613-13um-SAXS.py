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

# needs to be changed by Chen Chen
runs = [[153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163]]

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

