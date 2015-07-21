#!/opt/local/bin/python2.7

# GENERAL
import os
from itertools import dropwhile, islice
from math import log
from optparse import OptionParser

# ------------------------------------------------------------------------------------------------------------------------

class MAFregion:
	"""MAF region class"""

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
	def to_maf(self):
		return "s %s\t\t%s %s %s %s %s"%(self.loc, self.pos, self.length, self.strand,"00000",self.seq)

class MAFmapping:
	"""MAF mapping class"""
	
	def __init__(self, score, ref, aln):
		self.score = score
		self.ref = MAFregion(ref)
		self.aln = MAFregion(aln)

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
		
	def to_maf(self):
		return "a score=%s\n%s\n%s\n"%(self.score, self.ref.to_maf(), self.aln.to_maf())
		


# ------------------------------------------------------------------------------------------------------------------------

parser = OptionParser()
parser.add_option("--maf",   dest="maf_file",	 help="Path of input MAF file to parse",	default=False)
parser.add_option("--out",   dest="out_file",	 help="Path of output MAF file to write",	default=False)
parser.add_option("--qual",  dest="qual_score",  help="Minimum quality for alignments",		default=300)
parser.add_option("--2d",    dest="twod_only",	 help="Flag to use 2D reads only",		default=True)
(options, args) = parser.parse_args()

# ------------------------------------------------------------------------------------------------------------------------

def check_arguments(options):
	#print("Checking arguments")
	if not os.path.exists(options.maf_file):
		print("Invalid MAF file %s"%(options.maf_file))
		return False

	#print options.qual_score + 0
	if options.qual_score < 0 or not isinstance( options.qual_score, int ):
		print("Invalid Quality Score cut-off %s"%(options.qual_score))
		return False
	
	print "-"*30
	print "Running with the following optios:"
	print options
	print "-"*30
	return True

def isHeaderLine(line):
	return line.startswith("#")

def filter_alt_mappings(options):
	collection = {}
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

			# Store strand information +=1  -=-1
			strand = int(aln[4]+'1')
				
			# Check if alignment passes QC
			qcpass = False
			
			if (score >= options.qual_score):
				# 2D reads always used 
				if "_2D_" in aln[1]:
					qcpass = True
				# 1D reads only when specified
				elif not options.twod_only:
					qcpass = True
			
			# Store alignments that pass
			if qcpass:
				if (aln[1] not in collection):
		          		collection[aln[1]] = []
				collection[aln[1]].append(MAFmapping(score, ref, aln))

			#print "%i %s:%s-%s  -> %s:%s-%s" %(score, aln[1], aln[2], aln[3], ref[1], ref[2], ref[3])
	return collection

def write_alt_mappings(options, collection):
	outf = open(options.out_file, 'w')
	
	inf = open(options.maf_file,'r')
	for line in inf:
		if isHeaderLine(line):
			outf.write(line)
		else:
			break
	inf.close()
	
	for read in collection:
		for aln in colleaction[read]:
			outf.write(aln.to_maf())
	outf.close()
	
# ------------------------------------------------------------------------------------------------------------------------

if check_arguments(options):
	filtered_reads = filter_alt_mappings(options)
	write_alt_mappings(options, filtered_reads)

print("DONE")

