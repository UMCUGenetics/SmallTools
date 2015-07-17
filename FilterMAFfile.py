#!/opt/local/bin/python2.7

# GENERAL
import os
from itertools import dropwhile, islice
from math import log

# MAF handling
from pysam import AlignmentFile, AlignedSegment


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
		return "%s:%i-%i %s" %(self.loc, self.pos, self.pos+self.length, self.strand)

	def __str__(self):
		return "%s:%i-%i %s" %(self.loc, self.pos, self.pos+self.length, self.strand)

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
parser.add_option("--bam",   dest="bam_file",	   help="Path of BAM file to parse",		default=False)
parser.add_option("--maf",   dest="maf_file",	   help="Path of MAF file to parse",		default=False)
parser.add_option("--qual",  dest="qual_score",  help="Minimum quality for alignments",		default=300)
(options, args) = parser.parse_args()

# ------------------------------------------------------------------------------------------------------------------------

def check_arguments(options):
	#print("Checking arguments")
	if not os.path.exists(options.bam_file):
		print("Invalid BAM file %s"%(options.bam_file))
		return False

	if not os.path.exists(options.maf_file):
		print("Invalid MAF file %s"%(options.maf_file))
		return False

	return True



def isHeaderLine(line):
	return line.startswith("#")

def gather_alt_mappings(options, collection):
	# Parse MAF file
	with open(options.maf_file,'r') as f:
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

			if (score >= options.qual_score):
			  if (aln[1] not in collection):
          collection[aln[1]] = []
				collection[aln[1]].append(LASTmapping(score, ref, aln))

				#print "%i %s:%s-%s  -> %s:%s-%s" %(score, aln[1], aln[2], aln[3], ref[1], ref[2], ref[3])

# ------------------------------------------------------------------------------------------------------------------------

if check_arguments(options):
	collection = {}
	gather_alt_mappings(options, collection)
	write_alt_mappings(options, collection)

print("DONE")

