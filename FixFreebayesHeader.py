#!/usr/local/bin/python
#fix freebayes header
import os
import glob
from optparse import OptionParser
# -------------------------------------------------
parser = OptionParser()
parser.add_option("--vcfdir",	dest="vcfdir",		help="Path to directory containing VCF files",		default=False)
(options, args) = parser.parse_args()
# -------------------------------------------------

SAMPLEHEADER="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s	%s\n"

# -------------------------------------------------

# CHECK AND GENERATE GZ AND TBI
def fix_header(vcffile):
	freader = open(vcffile, 'r')
	fwriter = open(vcffile.replace(".vcf","_fixed.vcf"), 'w')

	samples = []
	header = False
	for line in freader:
		if line.startswith("##"):
			fwriter.write(line)
			if line.startswith("##commandline="):
				##commandline="/home/cog/pprins/run6/bin/freebayes -f /hpc/cog_bioinf/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta -C 3 -t /hpc/cog_bioinf/ENRICH/kinome_design_SS_V2_110811.bed --pooled-discrete --genotype-qualities --min-coverage 5 --no-indels --no-mnps --no-complex /home/cog/pprins/run6/data/freebayes/merged_MBC019R_F3_20130528_rmdup_kinome_design_SS_V2_110811.bam /home/cog/pprins/run6/data/freebayes/merged_MBC019T_F3_20130528_rmdup_kinome_design_SS_V2_110811.bam""
				items = line.strip().split(" ")[-2:-1]
				print items
				samples = [k.split("_")[1] for k in items]
				print samples
		elif not header:
			fwriter.write(SAMPLEHEADER%(samples[0], samples[1]))
			header=True
		else:
			fwriter.write(line)

	freader.close()
	fwriter.close()

# -------------------------------------------------

file_list = glob.glob(os.path.join(options.vcfdir, "*.vcf"))
for vcf_file in file_list:
	fix_header(vcf_file)

os.system("mkdir fixed")
os.system("mv *_fixed.vcf fixed")

# -------------------------------------------------
