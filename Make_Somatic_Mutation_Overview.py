#!/usr/bin/python
import os
import re
import vcf
import glob
import numpy as np

#GENE FORMAT
##chr	start	stop	name
#3	178866311	178952497	PIK3CA

from optparse import OptionParser
# -------------------------------------------------
parser = OptionParser()
parser.add_option("--vcfdir",	dest="vcfdir",		help="Path to directory containing VCF files",	default=False)
parser.add_option("--outdir",	dest="outdir",		help="Path to directory to write output to",	default="./DriverProfile/")
parser.add_option("--genelist",	dest="genelist",	help="File containing Genes to test/plot)",		default=False)

parser.add_option("--bgzip",	dest="bgzip",		help="Path to bgzip binary",				default="bgzip")
parser.add_option("--tabix",	dest="tabix",		help="Path to tabix binary",				default="tabix")

parser.add_option("--t",	dest="nrcpus",		help="Number of CPUs to use per sample",		default=2)

parser.add_option("--dp",	dest="mindepth",	help="Minumum read depth to consider reliable",	default=10)
parser.add_option("--af",	dest="minvaf",		help="Minumum variant allele fraction",			default=0.15)
parser.add_option("--pf",	dest="popfreq",		help="Maximum popultaion frequency",			default=0.05)
(options, args) = parser.parse_args()
# -------------------------------------------------

vocabulary = {"None":-1, "clean":0, "sequence_feature":0, "synonymous_variant":0, "intron_variant":0, "3_prime_UTR_variant":0.5, "5_prime_UTR_variant":0.5, "non_coding_exon_variant":0.5, "TF_binding_site_variant":1.0, "missense_variant":1.5, "splice_region_variant":2, "splice_donor_variant":2, "splice_acceptor_variant":2, "inframe_deletion":2.1, "inframe_insertion":2.1, "disruptive_inframe_deletion":2.5, "disruptive_inframe_insertion":2.5, "5_prime_UTR_premature_start_codon_gain_variant":3, "stop_gained":4, "nonsense_mediated_decay":4, "frameshift_variant":5}
toselect = [k for k,v in vocabulary.items() if v > 1]

# -------------------------------------------------
debug = False
# -------------------------------------------------
def check_arguments():
	if not os.path.exists(options.vcfdir):
		print("Invalid BAM folder %s"%(options.vcfdir))
		return False

	if not os.path.exists(options.outdir):
		print("Creating output folder %s"%(options.outdir))
		try:
			os.makedir(options.outdir)
		except OSError:
			print("Invalid / unable to create, output folder %s"%(options.outdir))
			return False

	print("Running with the following settings:")
	print("------------------------------------")
	print(options)
	print("------------------------------------")
	return True

# -------------------------------------------------

# Extract population frequency from VCF record
# Annoation assumed to be in SNPeff formatting
def find_popfreq(vcf_record):
	popfreq=[0.0]
	#print(vcf_record.INFO)
	freq_fields = ["dbNSFP_ExAC_AF", "dbNSFP_ExAC_Adj_AF", "GoNLv5_Freq"]

	for field in freq_fields:
		if field in vcf_record.INFO:
			#print vcf_record.INFO
			popfreq.append([float(x) for x in vcf_record.INFO[field].split(",")])

	#print(popfreq)
	return(popfreq)

# Determine the most damaging effect of the variant
def find_effects(vcf_record):
	maxeffect="None"
	if debug: print(vcf_record.INFO)

	if "ANN" not in vcf_record.INFO:
		return maxeffect

	# TRAVERSE ALL ANNOTATIONS
	for pred in vcf_record.INFO["ANN"]:
		# SPLIT THE SEPERATE FIELDS WITHIN THE ANNOTATION
		items = pred.split("|")
		allele = items[0]
		effects = items[1].split("&")
		for effect in effects:
			if effect not in vocabulary:
				# A NEW MUTATION EFFECT WAS FOUND
				if debug:
					print "NEW Mutation effect identified:"
					print pred
					print effect

			else:
				# STORE THE MOST DELETERIOUS EFFECT
				if vocabulary[effect] > vocabulary[maxeffect]:
					maxeffect = effect
	if debug:
		print maxeffect
	return(maxeffect)

# ETRACT THE MOST DELETERIOUS MUTATIONS IN A GENE
def select_maximum_effect(effects):
	effectvalues = [vocabulary[eff] for eff in effects]
	if debug: print(effectvalues)
	indices = np.argmax(effectvalues)
	return(effects[indices])

# CHECK AND GENERATE GZ AND TBI
def zip_and_index(vcffile):
	if not os.path.exists(vcffile+".gz"):
		os.system(options.bgzip+" -c "+vcffile+" > "+vcffile+".gz")
	if not os.path.exists(vcffile+".gz"+".tbi"):
		os.system(options.tabix+" "+vcffile+".gz")

# -------------------------------------------------

def main():
	file_list = glob.glob(os.path.join(options.vcfdir, "*.vcf"))
	for vcf_file in file_list:
		zip_and_index(vcf_file)

	genelist=open(options.genelist, 'r').read().split('\n')

	df = {}

	for vcf_file in file_list:
		if (debug): print vcf_file
		vcfread = vcf.Reader(open(vcf_file+".gz",'r'), compressed="gz")
		sample = vcfread.samples[0]
		print sample
		df[sample] = {}

		# FOR EACH GENE OF INTREST
		for gene in genelist:
			if len(gene)<=0:
				continue
			thisgene = dict(zip(["Chr","Start","Stop","SYMBOL"], gene.strip().split('\t')))

			if (debug):
				print(thisgene)

			print(thisgene)
			# FOR EACH TUMOR SAMPLE
			vcf_records = vcfread.fetch(thisgene["Chr"], int(thisgene["Start"])-20, int(thisgene["Stop"])+20)

			effects = []
			# FILTER NON-QC RECORDS
			for vcf_record in vcf_records:
				print vcf_record
				# CHECK TOTAL COVERAGE OF IDENTIFIED ALLELLES
				if isinstance(vcf_record.genotype(sample)['AD'], int):
					#if vcf_record.genotype(sample)['AD'] < options.mindepth:
					# IGNORE SINGLE AD VALUE SAMPLES
					continue
				if sum(vcf_record.genotype(sample)['AD']) < options.mindepth:
					continue

				# CHECK VAF
				#if debug: print sum(vcf_record.genotype(sample)['AD'][1:])*1.0/sum(vcf_record.genotype(sample)['AD'])
				if (sum(vcf_record.genotype(sample)['AD'][1:])*1.0/sum(vcf_record.genotype(sample)['AD'])) < options.minvaf:
					continue

				# CHECK POPULATION FREQUENCY
				# print find_popfreq(vcf_record)
				if max(find_popfreq(vcf_record)) > options.popfreq:
					continue

				effects.append(find_effects(vcf_record))

			if len(effects) <= 0:
				df[sample][thisgene["SYMBOL"]] = "None"
			else:
				df[sample][thisgene["SYMBOL"]] = select_maximum_effect(effects)


		#print(sample, df[sample])

	print "Sample\t"+'\t'.join(df[sample].keys())
	for sp in df:
		print sp+'\t'+'\t'.join(df[sp].values())



# -------------------------------------------------

print("Starting Analysis")

if __name__ == '__main__':
	if check_arguments():
		main()
	else:
		print("Error in provided arguments")

print("DONE")

# -------------------------------------------------

def test():
	for j in range(0, len(gts)):
		effect = "None"
		if gts[j] == "0/0":
			effect = "clean"

		if gts[j] == "0/1" or gts[j] == "1/1":
			effect = posdict[i]["Effects"][0]

		if gts[j] == "0/2" or gts[j] == "2/2":
			effect = posdict[i]["Effects"][1]

		if gts[j] == "1/2":
			print("1/2 detected: ",i,j)
			effect = posdict[i]["Effects"][1]
			# FIXME

		if vocabulary[effect] > vocabulary[varcount[j]]:
			varcount[j] = effect

	measured = [x!="None" for x in varcount]
	affected = [x in toselect for x in varcount ]
	samplenames = list(gtdf.columns.values)
	tumors = ["R" not in x for x in samplenames]
	newdat = dict(zip(samplenames, varcount))
	df[thisgene[3]] = pd.Series(newdat)

	print(df)
	df.to_csv("MutationOverview.txt",sep='\t')
