#!/usr/bin/env python

import glob as G
import sys, os, re, shutil, subprocess, time
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="filename", help="file with hits you wish to check type distribution", metavar="FILENAME", default="")
parser.add_option("-t", "--type", action="store", type="int", dest="type", help="type from which you want to find number of hits in file (default: type0)", metavar="TYPENUMBER", default=0)
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to check", metavar="xxxx", default="")
parser.add_option("-o", "--outputdir", action="store", type="string", dest="outputDir", help="output directory (default: output_rxxxx)", metavar="OUTPUT_DIR", default="output")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="prints out additional information", default=False)
(options, args) = parser.parse_args()

runtag = "r%s"%(options.runNumber)
write_dir = options.outputDir + '_' + runtag + '/'

# check if write_dir exists
if not os.path.exists(write_dir):
	print "There is no directory called %s. Aborting." % (write_dir)
	sys.exit(1)
else:
	# read hits from dir
	originaldir = os.getcwd()
	if (options.type == 0):
		write_anomaly_dir = write_dir 
	else:
		write_anomaly_dir = write_dir + "type%d" % (options.type)
		if not os.path.exists(write_anomaly_dir):
			if options.verbose:
				print "There is no directory called %s. Aborting." % (write_anomaly_dir)
			else:
				print 0
			sys.exit(1)
	foundHits = []
	os.chdir(write_anomaly_dir)
	foundHits += G.glob("LCLS*angavg.h5")
	os.chdir(originaldir)
	# read hits from file
	if not (options.filename == ""):
		fileHits = []
		file_path = write_dir + options.filename
		if os.path.exists(file_path):
			f = open(file_path, 'r')
			for line in f:
				if line[0] != '#':
					line = re.sub("\n", '', line)
					fileHits.append(line)
		else:
			print "There is no file called %s. Aborting." % (file_path)
			sys.exit(1)
	else:
		print "There is no file specified. Aborting."
		sys.exit(1)
	# make sets of hits
	sDirHits = set(foundHits)
	sFileHits = set(fileHits)
	sHits = sDirHits.intersection(sFileHits)
	if (options.verbose):
		print "Found %d hits from type%d in %s" % (len(sHits), options.type, options.filename)
		print "Also found %d hits in %s from other types." % (len(sFileHits.difference(sDirHits)), options.filename)
		print "Also found %d hits in type%d not present in %s" % (len(sDirHits.difference(sFileHits)), options.type, options.filename)
	else:
		print len(sHits)

