#!/usr/bin/python

import sys, os
import vcf
import tabix
from time import time 
from time import sleep
import multiprocessing as mp

from optparse import OptionParser

"""
CADD FORMAT
#Chrom	Pos Ref Alt RawScore	PHRED
1 10001 T A 0.176634	5.959
1 10001 T C 0.114925	5.063
"""
parser = OptionParser()
parser.add_option("--vcf",	 dest="vcf_file",	 help="Path to VCF to annotate",	default=False)
parser.add_option("--snv",	 dest="cadd_snvs",	 help="Path to prescored SNVs",		default=False)
parser.add_option("--indel", dest="cadd_indels", help="Path to prescored InDels",	default=False)
parser.add_option("--t",     dest="nr_cpus", 	 help="Number of CPUs to use",	    default=8)
parser.add_option("--out",	 dest="out_file",	 help="Path to output VCF file",	default="out.vcf")
(options, args) = parser.parse_args()

def check_files(options):
	print("Checking files")
	if not os.path.exists(options.vcf_file):
		print("Invalid VCF file %s"%(options.vcf_file))
		return False

	if not os.path.exists(options.cadd_snvs):
		print("Invalid CADD SNV file %s"%(options.cadd_snvs))
		return False
	if not os.path.exists(options.cadd_snvs+".tbi"):
		print("No Index for CADD SNV file %s"%(options.cadd_snvs+".tbi"))
		return False

	if not os.path.exists(options.cadd_indels):
		print("Invalid CADD InDel file %s"%(options.cadd_indels))
		return False
	if not os.path.exists(options.cadd_indels+".tbi"):
		print("No Index for CADD SNV file %s"%(options.cadd_indels+".tbi"))
		return False

	#print("Ready to go")
	return True

VCF_READER = vcf.Reader(open(options.vcf_file, 'r'))
VCF_WRITER = vcf.Writer(open(options.out_file, 'w'), VCF_READER)
VCF_WRITER.close()

def extract_CADD_score(arguments, q):
	vcf_record, caddfile = arguments
	
	tb = tabix.open(caddfile)
	records = tb.query((vcf_record.CHROM).replace("chr",""), vcf_record.POS-1, vcf_record.POS)

	vcf_record.INFO["RAW_CADD"]   = 0
	vcf_record.INFO["PHRED_CADD"] = 0


	for rec in records:
		if rec[3] == vcf_record.ALT[0]:
			vcf_record.INFO["RAW_CADD"]   = rec[4]
			vcf_record.INFO["PHRED_CADD"] = rec[5]
			break
	
	annotated = VCF_WRITER._map(str, [vcf_record.CHROM, vcf_record.POS, vcf_record.ID, vcf_record.REF]) + [VCF_WRITER._format_alt(vcf_record.ALT), str(vcf_record.QUAL) or '.', VCF_WRITER._format_filter(vcf_record.FILTER), VCF_WRITER._format_info(vcf_record.INFO)]

	q.put(annotated)
	return(annotated)


def listener(q):
	'''listens for messages on the q, writes to file. '''
	sys.stdout.write('Starting listener\n')
	
	f = open(options.out_file, 'wb') 
	f.write('#' + '\t'.join(VCF_WRITER.template._column_headers + VCF_WRITER.template.samples) + '\n')
	f.flush()
	
	while 1:
		m = q.get()
		#sys.stdout.write("----\n")
		#sys.stdout.write(str(m)+'\n')
		if m == 'kill':
			if not q.empty():
				sys.stdout.write('ERROR\n')
				break
			sys.stdout.write('DONE\n')
			break
		
		f.write('\t'.join(m)+'\n')
		f.flush()
	f.close()


def main():
	if not check_files(options):
		print("Error")
		exit(0)

	currtime = time()

	#must use Manager queue here, or will not work
	manager = mp.Manager()
	q = manager.Queue()		
	pool = mp.Pool(int(options.nr_cpus))

	#put listener to work first
	watcher = pool.apply_async(listener, (q,))

	print("Filling Queue")
	#fire off workers
	jobs = []
	for vcf_record in VCF_READER:
		arguments = []
		if vcf_record.is_indel:
			arguments = [vcf_record, options.cadd_indels]
		else:
			arguments = [vcf_record, options.cadd_snvs]

		job = pool.apply_async(extract_CADD_score, (arguments, q))
		jobs.append(job)
		


	print("Collecting results")
	# collect results from the workers through the pool result queue
	for job in jobs:
		job.get()
	
	#now we are done, kill the listener
	q.put('kill')
	
	pool.close()
	pool.join()

	print 'time elapsed:', time() - currtime

if __name__ == "__main__":
	main() 


