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
	
	def to_bed(self):
		return "%s\t%i\t%i\t%s" %(self.loc, self.pos, self.pos+self.length, self.strand)

	def get_loc(self):
		if self.loc == 'X':
			return 23
		elif self.loc == 'Y':
			return 24
		else:
			return int(self.loc)
	def to_maf(self):
		return "s %s %s %s %s %s %s"%(self.loc, self.pos, self.length, self.strand, self.totlen, self.seq)

class MAFmapping:
	"""MAF mapping class"""
	
	def __init__(self, score, ref, aln):
		self.score = score
		self.ref = MAFregion(ref)
		self.aln = MAFregion(aln)
		self.qual = []

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

	def set_qual(self, qual):
		self.qual = qual
		
	def to_maf(self):
		if self.qual != []:
			return "a score=%s\n%s\n%s\nq %s %s\n\n"%(self.score, self.ref.to_maf(), self.aln.to_maf(), self.aln.loc, self.qual)
		else:
			return "a score=%s\n%s\n%s\n\n"%(self.score, self.ref.to_maf(), self.aln.to_maf(), self.qual)


# ------------------------------------------------------------------------------------------------------------------------

parser = OptionParser()
parser.add_option("--maf",   dest="maf_file",	 help="Path of input MAF file to parse",	default=False)
parser.add_option("--out",   dest="out_file",	 help="Path of output MAF file to write",	default=False)
parser.add_option("--qual",  dest="qual_score",  help="Minimum quality for alignments",		default=300)
parser.add_option("--2d",    dest="twod_only",	 help="Flag to use 2D reads only",			default=False)
parser.add_option("--bq",    dest="hasb_qual",	 help="Flag to indicate precens of base quality scores", default=True)
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

			lines_gen = False
			qual = []

			if options.hasb_qual:
				lines_gen = islice(f, 5)

			else:
				lines_gen = islice(f, 4)

			score = int(line.strip().split("=")[1])
			ref = next(lines_gen).strip().split()
			aln = next(lines_gen).strip().split()
			if options.hasb_qual:
				qual = next(lines_gen).strip().split()
			empty = next(lines_gen)


			# Store strand information +=1  -=-1
			strand = int(aln[4]+'1')
				
			# Check if alignment passes QC
			qcpass = False
			
			if (score >= options.qual_score):
				# 2D reads always used 
				if "_2D_000_2d" in aln[1]:
					qcpass = True
				# 1D reads only when specified
				elif not options.twod_only:
					qcpass = True
			
			# Store alignments that pass
			if qcpass:
				if (aln[1] not in collection):
		          		collection[aln[1]] = []
				newmapping = MAFmapping(score, ref, aln)
				if options.hasb_qual:
					newmapping.set_qual(qual)
				collection[aln[1]].append(newmapping)
				
	return collection


#def find_optimal_path(options, collection):
	# for each sequenced read
#	for read in collection:
		# go through mapping options
#		for mapping in collection[read]:

	# KEEP TRACK OF UNMAPPED REGIONS
	# ITERATE UNTIL X% MAPPED OR NO MORE ADDITIONS POSSIBLE
	#
	# RANGE ==
	# 	mapping.aln.pos - mapping.aln.pos+mapping.aln.length

	

def write_alt_mappings(options, collection):
	outf = open(options.out_file, 'w')
	outf_optimal = open(options.out_file.replace(".maf",".optimal.maf"), 'w')

	# Get header from original file and write to new MAF file
	inf = open(options.maf_file,'r')
	for line in inf:
		if isHeaderLine(line):
			outf.write(line)
		else:
			break
	inf.close()
	
	# Go through reads and write Mapping info
	for read in collection:
		optimal = collection[read][0]

		for mapping in collection[read]:
			outf.write(mapping.to_maf())
			if mapping.score >= optimal.score:
				optimal=mapping

		outf_optimal.write(optimal.ref.to_bed()+"\t"+optimal.aln.loc+"\n")

	outf_optimal.close()
	outf.close()
	
# ------------------------------------------------------------------------------------------------------------------------

if check_arguments(options):
	filtered_reads = filter_alt_mappings(options)
	write_alt_mappings(options, filtered_reads)

print("DONE")

