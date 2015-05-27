#!/opt/local/bin/python2.7
import sys
import os
#import multiprocessing
import pysam
from pybedtools import BedTool

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--bam",   dest="bam_file",	 help="Path of BAM file to parse",			default=False)
parser.add_option("--last",  dest="last_file",	 help="Path of LAST file to parse",			default=False)
parser.add_option("--reg",   dest="region_file", help="Regions of interest BED file",		default=False)
parser.add_option("--qual",  dest="qual_score",  help="Minimum quality for alignments",		default=250)
parser.add_option("--nral",  dest="nr_align",	 help="Number of alignments to display",	default=10)
(options, args) = parser.parse_args()


def check_arguments(options):
	#print("Checking arguments")
	if not os.path.exists(options.bam_file):
		print("Invalid BAM file %s"%(options.bam_file))
		return False

	if not os.path.exists(options.region_file):
		print("Invalid region BED file %s"%(options.region_file))
		return False

	if not os.path.exists(options.last_file):
		print("Invalid LAST file %s"%(options.last_file))
		return False

	return True


def gather_sv_data(options):
	collection = []
	
	# Read regions of interest BED file
	regions = BedTool(options.region_file)

	# Read BAM file
	bamfile = pysam.AlignmentFile(options.bam_file, "rb")

	# Intersect regions
	for reg in regions:
		for read in bamfile.fetch(reg.chrom, reg.start, reg.end):
			#print read
			if read.query_name.endswith("2d"):
				collection.append(read.query_name)
				#print read.reference_id, read.reference_start, read.reference_end
				#print read.query_name, read.query_alignment_start, read.query_alignment_end

	bamfile.close()
	return collection


from itertools import *

def isHeaderLine(line):
	return line.startswith("#")

def gather_alt_mappings(options, collection):
	# Parse LAST file
	with open(options.last_file,'r') as f:
		for line in dropwhile(isHeaderLine, f):
			if isHeaderLine(line):
				continue

			lines_gen = islice(f, 4)

			a = line.strip()
			b = next(lines_gen).strip()
			c = next(lines_gen).strip()
			d = next(lines_gen)

			if (len(a) <= 5):
				print(score)
			score = int(a.split("=")[1])
			ref =  b.split()
			read = c.split()

			if (read[1] in collection and score >= options.qual_score):
				print read[1], score, ref[1], ref[2]




if check_arguments(options):
	collection = gather_sv_data(options)
	gather_alt_mappings(options, collection)

print("DONE")

