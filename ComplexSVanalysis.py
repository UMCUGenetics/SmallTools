#!/opt/local/bin/python2.7

# GENERAL
import os
from itertools import dropwhile, islice
from math import log

# BAM and BED handling
from pysam import AlignmentFile, AlignedSegment
from pybedtools import BedTool

# PLOTTING
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from optparse import OptionParser

# ------------------------------------------------------------------------------------------------------------------------

class LASTregion:
	"""LAST region class"""

	def __init__(self, array):
		self.loc = array[1]
		self.pos = int(array[2])
		self.length = int(array[3])
		self.strand = array[4]
		self.totlen = int(array[5])
		self.seq = array[6]

	def __repr__(self):
		return "%s:%i-%i" %(self.loc, self.pos, self.pos+self.length)

	def __str__(self):
		return "%s:%i-%i" %(self.loc, self.pos, self.pos+self.length)

	def get_loc(self):
		if self.loc == 'X':
			return 23
		elif self.loc == 'Y':
			return 24
		else:
			return int(self.loc)

class LASTmapping:
	"""LAST mapping class"""
	
	def __init__(self, score, ref, aln):
		self.score = score
		self.ref = LASTregion(ref)
		self.aln = LASTregion(aln)

	def __eq__(self, other):
		return self.score == other.score

	def __lt__(self, other):
		return self.score < other.score

	def __len__(self):
		return self.aln.totlen

	def __repr__(self):
		return "%i %s -> %s" %(self.score, self.aln, self.ref)
	
	def __str__(self):
		return "%i %s -> %s" %(self.score, self.aln, self.ref)


# ------------------------------------------------------------------------------------------------------------------------

parser = OptionParser()
parser.add_option("--bam",   dest="bam_file",	 help="Path of BAM file to parse",			default=False)
parser.add_option("--last",  dest="last_file",	 help="Path of LAST file to parse",			default=False)
parser.add_option("--reg",   dest="region_file", help="Regions of interest BED file",		default=False)
parser.add_option("--qual",  dest="qual_score",  help="Minimum quality for alignments",		default=300)
parser.add_option("--nral",  dest="nr_align",	 help="Number of alignments to display",	default=20)
(options, args) = parser.parse_args()

# ------------------------------------------------------------------------------------------------------------------------

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


def gather_sv_data(options, collection):
	# Read regions of interest BED file
	regions = BedTool(options.region_file)

	# Read BAM file
	bamfile = AlignmentFile(options.bam_file, "rb")

	# Intersect regions
	for reg in regions:
		for read in bamfile.fetch(reg.chrom, reg.start, reg.end):
			#print read
			if read.query_name.endswith("2d"):
				collection[read.query_name] = []
				#print read.reference_id, read.reference_start, read.reference_end
				#print read.query_name, read.query_alignment_start, read.query_alignment_end

	bamfile.close()


def isHeaderLine(line):
	return line.startswith("#")



def gather_alt_mappings(options, collection):
	# Parse LAST file
	with open(options.last_file,'r') as f:
		for line in dropwhile(isHeaderLine, f):
			if isHeaderLine(line):
				continue
			if len(line) <= 5:
				continue

			lines_gen = islice(f, 4)

			score = int(line.strip().split("=")[1])
			ref = next(lines_gen).strip().split()
			aln = next(lines_gen).strip().split()
			empty = next(lines_gen)

			strand = 1
			if aln[4] == '-':
				strand = -1

			if (aln[1] in collection and score >= options.qual_score):
				collection[aln[1]].append(LASTmapping(score, ref, aln))

				#print "%i %s:%s-%s  -> %s:%s-%s" %(score, aln[1], aln[2], aln[3], ref[1], ref[2], ref[3])


def plot_alt_mappings(options, collection):
	#GenomeDiagram.FeatureSet()
	color_list = plt.cm.Set1(np.linspace(0, 1, 24))
	chroms = [str(x) for x in range(1,25)]
	#print len(color_list)

	handles = []
	for j in range(0, 24):
		handles.append(mpatches.Patch(color=color_list[j], label=chroms[j]))
	
	#red_patch = mpatches.Patch(color='red', label='The red data')
	#plt.legend(handles=[red_patch])


	for read in collection:
		#print read, collection[read]

		fig = plt.figure()
		ax = fig.add_subplot(111)

		#alignments = collection[read].sort()
		alignments = collection[read]
		sortaln = sorted(alignments, reverse=True)
		#print(sortaln)

		if len(alignments) == 0:
			continue

		# SET plot limits
		ax.set_xlim(0, sortaln[0].aln.totlen)
		ax.set_ylim(-10, 10)

		for i in range(0, options.nr_align):
			
			score = log(sortaln[i].score)
			if sortaln[i].aln.strand == '-':
				score = score*-1

			ax.broken_barh((sortaln[i].aln.pos, sortaln[i].aln.length), (score-0.1, score+0.1), facecolors=color_list[sortaln[i].ref.get_loc()])

		# Mark-up of plot
		ax.set_yticks(range(-10, 10, 1))
		ax.grid(True)
		ax.set_xlabel('Location within read')
		ax.set_ylabel('Log(score)')

		# Add color legend
		ax.legend(handles=handles, lables=chroms)
		#, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		#plt.show()

		# SAVE file
		fig.savefig(read+'.pdf')

# ------------------------------------------------------------------------------------------------------------------------

if check_arguments(options):
	collection = {}
	gather_sv_data(options, collection)
	gather_alt_mappings(options, collection)
	plot_alt_mappings(options, collection)

print("DONE")

