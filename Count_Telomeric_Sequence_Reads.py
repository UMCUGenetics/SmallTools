#!/usr/bin/python2.7

import pysam
import os

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--sambamba",	 dest="sambamba", help="Path to sambamba/samtools executable",   default="sambamba")
parser.add_option("--bamdir",	 dest="bamdir",   help="Path to directory containing bam files", default=False)
parser.add_option("--outdir",	 dest="outdir",   help="Path to directory to write output to",   default="./telomeres/")

(options, args) = parser.parse_args()

#TODO implement argparse to get these variables
sambamba = "/path/to/sambamba/sambamba_vx.x"
workdir = "/path/to/working/directory/"
searchdir = "/path/to/search/bams/"

# -------------------------------------------------

def count_telomeric_reads(bamfile, telofile):
    if not os.path.exists(telofile):
        # extract telomeric reads and write to file
        os.system(sambamba + " view " + bamfile + " | LC_ALL=C grep -E \"" + "TTAGGG"*15 +"|"+ "CCCTAA"*15 + "\"" + " > " + telofile + " & ")

    # count total number of reads
    total_rc = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfile) ])
            
    # count number of telomeric reads
    telomere_rc = os.system("wc -l "+telofile)
            
    # print counts
    print(bamfile, total_rc, telomere_rc)
    
    
def main():    
    for root, dirs, files in os.walk(searchdir):
        for filename in files:
            # select bam files
            if filename.endswith(".bam"):
                #print(filename)
                bamfile = os.path.join(root, filename)
                baifile = bamfile+".bai"
                if not os.path.exists(baifile):
                    print("No index file found, indexing")
                    os.system(sambamba + " index " + bamfile)
                    print("Indexing performed")
                    
                telofile = workdir+filename.replace(".bam","_TelomericReads.sam")
            
                count_telomeric_reads(bamfile, telofile)

# -------------------------------------------------
# Execute program

print("Starting analysis")
main()
print("DONE")


