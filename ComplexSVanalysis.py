#!/opt/local/bin/python2.7
import sys
import os
#import multiprocessing
from itertools import *

# BAM and BED handling
from pysam import AlignmentFile
from pybedtools import BedTool

# BIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

# Plotting
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

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
	bamfile = AlignmentFile(options.bam_file, "rb")

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




def isHeaderLine(line):
	return line.startswith("#")

def gather_alt_mappings(options, collection):
	# Prepare feature set
	alt_mappings = {}
	for read in collection:
		alt_mappings[read] = SeqRecord(Seq("ATCGTCGTA"), id=read, name=read)

	# Parse LAST file
	with open(options.last_file,'r') as f:
		for line in dropwhile(isHeaderLine, f):
			if isHeaderLine(line):
				continue
			if len(line) <= 5:
				continue

			lines_gen = islice(f, 4)

			score = int(line.strip().split("=")[1])
			ref =  next(lines_gen).strip().split()
			read = next(lines_gen).strip().split()
			empty = next(lines_gen)

			strand = 1
			if read[4] == '-':
				strand = -1

			if (read[1] in collection and score >= options.qual_score):
				feature = SeqFeature(FeatureLocation(int(read[2]), int(read[3]), strand=strand), type="Read", ref=ref[1])
				alt_mappings[read[1]].features.append(feature)

				print "%i %s:%s-%s  -> %s:%s-%s" %(score, read[1], read[2], read[3], ref[1], ref[2], ref[3])

	return alt_mappings


def plot_alt_mappings(options, alt_mappings):
	#GenomeDiagram.FeatureSet()
	color = int(ref[1])

	for read in alt_mappings:
		alt_mappings[read].add_feature(feature, color=color, label=True)



if check_arguments(options):
	collection = gather_sv_data(options)
	alt_mappings = gather_alt_mappings(options, collection)

print("DONE")

