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
parser.add_option("--vcf",   dest="vcf_file",	 help="Path to VCF to annotate",    default=False)
parser.add_option("--snv",   dest="cadd_snvs",	 help="Path to prescored SNVs",	    default=False)
parser.add_option("--indel", dest="cadd_indels", help="Path to prescored InDels",   default=False)
parser.add_option("--t",     dest="nr_cpus", 	 help="Number of CPUs to use",	    default=8)
parser.add_option("--out",   dest="out_file",	 help="Path to output VCF file",    default="out.vcf")
(options, args) = parser.parse_args()


def check_arguments(options):
	#print("Checking arguments")
	if not os.path.exists(options.vcf_file):
		print("Invalid VCF file %s"%(options.vcf_file))
		return False
	
	# ---- SNV file ---
	if not os.path.exists(options.cadd_snvs):
		print("Invalid CADD SNV file %s"%(options.cadd_snvs))
		return False
	if not os.path.exists(options.cadd_snvs+".tbi"):
		print("No Index for CADD SNV file %s"%(options.cadd_snvs+".tbi"))
		return False

	# ---- InDel file ---
	if not os.path.exists(options.cadd_indels):
		print("Invalid CADD InDel file %s"%(options.cadd_indels))
		return False
	if not os.path.exists(options.cadd_indels+".tbi"):
		print("No Index for CADD InDel file %s"%(options.cadd_indels+".tbi"))
		return False

	# ---- Other settings ----
	try int(options.nr_cpus):
		pass
	except Exception, e:
		print("Invalid nr of cpus defined %s"%(options.nr_cpus))
		return False

	return True

VCF_READER = vcf.Reader(open(options.vcf_file, 'r'))
VCF_WRITER = vcf.Writer(open(options.out_file, 'w'), VCF_READER)
VCF_WRITER.close()

VALID_CHROMOSOMES = {"1":True,"2":True,"3":True,"4":True,"5":True,"6":True,"7":True,"8":True,"9":True,"10":True,"11":True,"12":True,"13":True,"14":True,"15":True,"16":True,"17":True,"18":True,"19":True,"20":True,"21":True,"22":True,"X":True,,"Y":True}

def extract_CADD_score(arguments, q):
	vcf_record, caddfile = arguments
	
	tb = tabix.open(caddfile)

	chromosome = (vcf_record.CHROM).replace("chr","")
	vcf_record.INFO["RAWCADD"]   = 0
	vcf_record.INFO["PHREDCADD"] = 0

	# Specific for CADD files
	# FIXME: get info about chr or not from provided VCF file
	records = tb.query(chromosome, vcf_record.POS-1, vcf_record.POS)

	# Look for matching mutation
	# Works for SNVs, InDels optimisation is ongoing
	for rec in records:
		if rec[3] == vcf_record.ALT[0]:
			# FIXME: Make requested fields optional through arguments
			vcf_record.INFO["RAWCADD"]   = rec[4]
			vcf_record.INFO["PHREDCADD"] = rec[5]
			break
	
	# workaround since multiprocess can't handle VCF record class objects
	# FIXME: use VCF class records rather than this ugly string
	annotated = VCF_WRITER._map(str, [vcf_record.CHROM, vcf_record.POS, vcf_record.ID, vcf_record.REF]) + [VCF_WRITER._format_alt(vcf_record.ALT), str(vcf_record.QUAL) or '.', VCF_WRITER._format_filter(vcf_record.FILTER), VCF_WRITER._format_info(vcf_record.INFO)]

	# Return results to Queue
	q.put(annotated)
	return(annotated)


def listener(q):
	'''listens for messages on the q, writes to file. '''
	#sys.stdout.write('Starting listener\n')
	
	f = open(options.out_file, 'wb')
	#FIXME: get the rest of the header
	f.write("##INFO=<ID=PHREDCADD,Number=1,Type=Float,Description=\"PHRED scaled CADD score\">")
	f.write("##INFO=<ID=RAWCADD,Number=1,Type=Float,Description=\"Raw CADD score\">") 
	f.write('#' + '\t'.join(VCF_WRITER.template._column_headers + VCF_WRITER.template.samples) + '\n')
	f.flush()
	
	while 1:
		m = q.get()
		if m == 'kill':
			if not q.empty():
				# received kill signal without finishing all the processes
				sys.stdout.write('ERROR\n')
				break
			# received kill signal, finished all the processes, done
			sys.stdout.write('DONE\n')
			break
		
		# A vcf record was found, write to file
		f.write('\t'.join(m)+'\n')
		f.flush()
	f.close()


def main():
	if not check_arguments(options):
		print("Error in provided arguments")
		exit(0)

	currtime = time()

	#Init Manager queue
	manager = mp.Manager()
	q = manager.Queue()
	# Init worker pool
	pool = mp.Pool(int(options.nr_cpus))

	#Init Listener
	watcher = pool.apply_async(listener, (q,))

	#print("Filling Queue")
	#fire off workers
	jobs = []
	for vcf_record in VCF_READER:
		chromosome = (vcf_record.CHROM).replace("chr","")
		if chromosome not in VALID_CHROMOSOMES:
			continue
			
		arguments = []
		if vcf_record.is_indel:
			arguments = [vcf_record, options.cadd_indels]
		else:
			arguments = [vcf_record, options.cadd_snvs]

		job = pool.apply_async(extract_CADD_score, (arguments, q))
		jobs.append(job)
		


	#print("Collecting results")
	# collect results from the workers through the pool result queue
	for job in jobs:
		job.get()
	
	# now we are done, kill the listener
	q.put('kill')
	
	pool.close()
	pool.join()

	print 'time elapsed:', time() - currtime

if __name__ == "__main__":
	main() 


