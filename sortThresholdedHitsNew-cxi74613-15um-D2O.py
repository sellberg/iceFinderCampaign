#!/usr/bin/env python

# Written by J. Sellberg and Chen Chen on 2014-02-13
# Wrapper or sortHits for the new data from Jan 2013
# using the CXI Feb2011-2 nozzle with D2O at 600 PSI N2 (gas), 320 PSI He (liquid), driven at 100 kHz (30 Vpp)

import numpy as N
import glob as G
import sys, os, re, shutil, subprocess, time
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--remote", action="store_true", dest="remote", help="saves the output text files remotely in $HOME", default=False)
(options, args) = parser.parse_args()

# Needs to be changed by Chen Chen
runs = [[180], [182, 183], [193, 194, 195, 196, 197, 198], [202, 203], [205, 206], [210, 211], [216, 217], [223, 224, 225], [227, 228, 229]]

original_dir = os.getcwd() + '/'

# create hit lists for each type for each run
for i in N.arange(len(runs)):
	for j in N.arange(len(runs[i])):
		if (runs[i][j] < 100):
			run_tag = "r00%s"%(runs[i][j])
		else:
			run_tag = "r0%s"%(runs[i][j])
		run_file =  "output_" + run_tag + '/below300above200ADUs.txt'
		if os.path.isfile(run_file):
			os.system("./sortThresholdedHitsNew " + run_file)
		else:
			print run_file + " does not exist, skipping " + run_tag

