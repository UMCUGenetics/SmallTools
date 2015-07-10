#!/usr/bin/python2.7

import pysam
import glob
import os

import subprocess
from multiprocessing import Pool

from optparse import OptionParser
# -------------------------------------------------
parser = OptionParser()
parser.add_option("--sambamba",	 dest="sambamba", help="Path to sambamba/samtools executable", default="sambamba")
parser.add_option("--bamdir",	 dest="bamdir",   help="Path to directory containing BAM files", default=False)
parser.add_option("--outdir",	 dest="outdir",   help="Path to directory to write output to", default="./telomeres/")
parser.add_option("--repsize",   dest="repsize",  help="Number of required matching 6mers (TTAGGG)", default=10)

(options, args) = parser.parse_args()
# -------------------------------------------------

def check_arguments():
    if not os.path.exists(options.bamdir):
        print("Invalid BAM folder %s"%(options.bamdir))
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

def count_telomeric_reads(bamfile):

    # generate Telomere reads file name   
    telofile = bamfile.replace(options.bamdir,options.outdir).replace(".bam","_TelomericReads.sam")
    
    # check if the file was already generated
    if not os.path.exists(telofile):
        # extract telomeric reads and write to file
        cmd = options.sambamba + " view " + bamfile + " | LC_ALL=C grep -E \"" + "TTAGGG"*options.repsize +"|"+ "CCCTAA"*options.repsize + "\"" + " > " + telofile
        print("Generating sam file: "+telofile)
        p = subprocess.Popen(cmd)
        p.wait()

    # count total number of reads
    total_rc = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfile) ])
            
    # count number of telomeric reads by line count
    telomere_rc = sum(1 for line in open(telofile,'r'))
    
    # return sample ID and count stats
    return('\t'.join([bamfile.split("/")[-1].split("_")[0],str(total_rc), str(telomere_rc), str((telomere_rc/(total_rc*1.0))*100000.0)])+'\n')
    
# -------------------------------------------------

print("Starting Analysis")


if __name__ == '__main__':
    # check specified options
    if not check_arguments():
        break
    
    bamfiles = glob.glob(os.path.join(options.bamdir, "*.bam"))
    pool = Pool(processes=len(bamfiles)) 
    
    # check index for all bamfiles
    for bamfile in bamfiles:
        
        baifile = bamfile+".bai"
        # check if index file exists
        if not os.path.exists(baifile):
            print("No index file found for %s, indexing now"%(bamfile))
            subprocess.call(options.sambamba, " index " + bamfile)
    
    # generate Telomere reads file
    counts = pool.map(count_telomeric_reads, bamfiles)

    # open output file
    output = open(os.path.join(options.outdir, "TelomereCounts.txt"),'w')

    # write header
    output.write('\t'.join(["#Sample","TotalReads","TelomericReads","NormalisedFraction"])+'\n')
    # write results
    for count in counts: 
        output.write(count)
    
    output.close()

print("DONE")

# -------------------------------------------------
