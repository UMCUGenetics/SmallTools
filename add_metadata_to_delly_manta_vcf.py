#!/usr/bin/env python

"""
add_metadata_to_delly_manta_vcf.py

##INFO=<ID=SOMATIC  is automatically added after i) the delly somatic filtering step ii) manta somatic calling

"""

import sys
import argparse
import re
import vcf

def add_meta2vcf(vcf_file):
	try:
		f = open(vcf_file, 'r')
	except IOError:
		sys.exit("Error: Can't open vcf file: {0}".format(vcf_file))
	else:
		with f:
			vcf_header = []
			vcf_chrom = []
			vcf_variants = []
			countline = 0
			version = ""
			for line in f:
				countline = countline + 1
				line = line.strip('\n')
				if line.startswith('##'):
					##Print original vcf meta-information lines##
					vcf_header.append(line)
					if line.startswith('##cmdline'):
						find_manta = re.findall('(manta_\w+\.\w+\.\w+)', line)
						version = find_manta[0]
					## print header lines and Add meta-information lines with melter info to vcf
				elif line.startswith("#CHROM"):
					vcf_chrom.append(line)
					countline = 0
				else:
					variant = line.split('\t')
					if countline == 1:
						find_delly = re.findall('(EMBL.DELLYv\w+\.\w+.\w+)', variant[7])
						if not version:
							version = find_delly[0]
					if not "DELLY" in variant[7]:
						variant[7] = variant[7]+";SVMETHOD={0}".format(find_manta[0])
					vcf_variants.append("\t".join(variant))
					
					
			print "\n".join(vcf_header)
			print "##INFO=<ID=caller={0}".format(version)
			print "\n".join(vcf_chrom)
			print "\n".join(vcf_variants)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = 'add metadata to a SV VCF file')
	required_named = parser.add_argument_group('required named arguments')
	required_named.add_argument('-v', '--vcf_file', help='path/to/file.vcf', required=True)
	args = parser.parse_args()
	add_meta2vcf(args.vcf_file)

