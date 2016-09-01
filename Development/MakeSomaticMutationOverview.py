#!/usr/bin/python
import os
import re
import vcf

#GENE FORMAT
#3	178866311	178952497	PIK3CA

from optparse import OptionParser
# -------------------------------------------------
parser = OptionParser()
parser.add_option("--vcfdir",	dest="vcfdir",		help="Path to directory containing VCF files",		default=False)
parser.add_option("--outdir",	dest="outdir",		help="Path to directory to write output to",		default="./DriverProfile/")
parser.add_option("--genelist",	dest="genelist",	help="File containing Genes to test/plot)",		default=False)

parser.add_option("--bgzip",	dest="bgzip",		help="Path to bgzip binary",				default="bgzip")
parser.add_option("--tabix",	dest="tabix",		help="Path to tabix binary",				default="tabix")

parser.add_option("--t",	dest="nrcpus",		help="Number of CPUs to use per sample",		default=2)

parser.add_option("--dp",	dest="mindepth",	help="Minumum read depth to consider reliable",		default=10)
parser.add_option("--af",	dest="minvaf",		help="Minumum variant allele fraction",			default=0.15)
parser.add_option("--pf",	dest="popfreq",		help="Maximum popultaion frequency",			default=0.05)
(options, args) = parser.parse_args()
# -------------------------------------------------

vocabulary = {"None":-1, "clean":0, "sequence_feature":0, "synonymous_variant":0, "intron_variant":0, "3_prime_UTR_variant":0.5, "5_prime_UTR_variant":0.5, "non_coding_exon_variant":0.5, "missense_variant":1.5, "splice_region_variant":2, "splice_donor_variant":2, "splice_acceptor_variant":2, "inframe_deletion":2.1, "stop_gained":4, "nonsense_mediated_decay":4, "frameshift_variant":5}
toselect = ["missense_variant", "splice_region_variant", "inframe_deletion", "stop_gained", "nonsense_mediated_decay", "frameshift_variant"]

# -------------------------------------------------




def find_popfreq(vcf_record):
	popfreq=0.0

	print(vcf_record.INFO["dbNSFP_ExAC_AF"])
	ExAC_search = vcf_record.INFO["dbNSFP_ExAC_AF"]
	if ExAC_search:
		popfreq = [float(x) for x in ExAC_search.group(1).split(",")]

	ExAC_search = vcf_record.INFO["dbNSFP_ExAC_Adj_AF"]
	if ExAC_search:
		popfreq = [float(x) for x in ExAC_search.group(1).split(",")]

	GoNL_search = vcf_record.INFO["GoNLv5_Freq"]
	if GoNL_search:
		popfreq = [float(x) for x in GoNL_search.group(1).split(",")]

	return(popfreq)

def find_effects(vcf_record):

	print(vcf_record.INFO["ANN"])
	if vcf_record.INFO["ANN"]:
		# STORE ALL ANNOTATIONS
		ann = vcf_record.INFO["ANN"].split(",")
		for pred in ann:
			# SPLIT THE SEPERATE FIELDS WITHIN THE ANNOTATION
			items = pred.split("|")
			allele = items[0]
			effects = items[1].split("&")
			for effect in effects:
				if effect not in vocabulary:
					# A NEW MUTATION EFFECT WAS FOUND
					print effect
					print ann
				else:
					# STORE THE MOST DELETERIOUS EFFECT
					if vocabulary[effect] > vocabulary[posdict[varcounter]["Effects"][alt.index(allele)]]:
						posdict[varcounter]["Effects"][alt.index(allele)] = effect


def extract_gene(vcffile, thisgene):
	# SUBSET VCF FOR GENE REGION
	tb = tabix.open(vcffile)
	return(tb.query(thisgene["Chr"], thisgene["Start"]-20, thisgene["Stop"]+20))


def zip_and_index(vcffile):
	# CHECK AND GENERATE GZ AND TBI
	if not os.path.exists(vcffile+".gz"):
		os.system(options.bgzip+" -c "+vcffile+" > "+vcffile+".gz")
	if not os.path.exists(vcffile+".gz"+".tbi"):
		os.system(options.tabix+" "+vcffile+".gz")

# -------------------------------------------------

def main():
	file_list = glob.glob(os.path.join(options.vcfdir, "*.bam"))
	for vcf_file in file_list:
		zip_and_index(vcf_file)

	genelist=open('genelist.txt','r')

	# FOR EACH GENE OF INTREST
	for gene in genelist:
		thisgene = dict(zip(["Chr","Start","Stop","SYMBOL"], gene.strip().split('\t')))
		print(thisgene)

		# FOR EACH TUMOR
		for vcf_file in file_list:
			vcf_records = extract_gene(vcf_file, thisgene)

			# FILTER NON-QC RECORDS
			for vcf_record in vcf_records:
				# FIXME
				# CHECK TOTAL COVERAGE
				if sum(vcf_record.SAMPLE.AD) < options.mindepth:
					continue

				# CHECK VAF
				if (sum(vcf_record.SAMPLE.AD[1:])/sum(vcf_record.SAMPLE.AD)) < options.minvaf:
					continue

				# CHECK POPULATION FREQUENCY
				if max(find_popfreq(vcf_record)) > options.popfreq:
					continue

				effect = find_effects(vcf_record)
	genelist.close()






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
	newdat = dict(zip(samplenames,varcount))
	df[thisgene[3]] = pd.Series(newdat)








print(df)
df.to_csv("MutationOverview.txt",sep='\t')
