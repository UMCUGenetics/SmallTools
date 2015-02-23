#!/usr/bin/python

import sys, os
import vcf
import multiprocessing as mp

from optparse import OptionParser


parser = OptionParser()
parser.add_option("--vcf",	 dest="vcf_file",	 help="Path to VCF to annotate",	default=False)
(options, args) = parser.parse_args()

def check_files(options):
	print("Checking files")
	if not os.path.exists(options.vcf_file):
		print("Invalid VCF file %s"%(options.vcf_file))
		return False
	return True

VCF_READER = vcf.Reader(open(options.vcf_file, 'r'))


def extract_CADD_score(arguments, q):
	vcf_record, caddfile = arguments
	
	q.put(vcf_record)
	return(vcf_record)


def listener(q):
	'''listens for messages on the q, writes to file. '''
	sys.stdout.write('Starting listener\n')
	
	while 1:
		m = q.get()
		if m == 'kill':
			if not q.empty():
				sys.stdout.write('ERROR\n')
				break
			sys.stdout.write('DONE\n')
			break



def main():
	if not check_files(options):
		print("Error")
		exit(0)

	#must use Manager queue here, or will not work
	manager = mp.Manager()
	q = manager.Queue()		
	pool = mp.Pool(2)

	#put listener to work first
	watcher = pool.apply_async(listener, (q,))

	print("Filling Queue")
	jobs = []
	for vcf_record in VCF_READER:
		arguments = []
		if vcf_record.is_indel:
			arguments = [vcf_record, options.vcf_file]
		else:
			arguments = [vcf_record, options.vcf_file]

		job = pool.apply_async(extract_CADD_score, (arguments, q))
		jobs.append(job)
		


	print("Collecting results")
	for job in jobs:
		job.get()
	
	#now we are done, kill the listener
	q.put('kill')
	
	pool.close()
	pool.join()


if __name__ == "__main__":
	main() 


