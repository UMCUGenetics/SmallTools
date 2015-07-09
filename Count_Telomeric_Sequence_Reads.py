#!/usr/bin/python2.7

import pysam
import glob
import os

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

def count_telomeric_reads(bamfile, telofile):
    # check if the file was already generated
    if not os.path.exists(telofile):
        # extract telomeric reads and write to file
        os.system(options.sambamba + " view " + bamfile + " | LC_ALL=C grep -E \"" + "TTAGGG"*options.repsize +"|"+ "CCCTAA"*options.repsize + "\"" + " > " + telofile + " & ")

    # count total number of reads
    total_rc = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfile) ])
            
    # count number of telomeric reads
    telomere_rc = os.system("wc -l "+telofile)
            
    # print counts
    return([total_rc, telomere_rc])
    
# -------------------------------------------------

def main():
    # check specified options
    if not check_arguments():
        return 1
    
    # open output file
    output = open(os.path.join(options.outdir, "TelomereCounts.txt"),'w')
    # write header
    output.write('\t'.join(["#Sample","TotalReads","TelomericReads","NormalisedFraction"])+'\n')
        
    # for all bamfiles
    for filename in glob.glob(os.path.join(options.bamdir, "*.bam")):
        print(filename)
        bamfile = os.path.join(options.bamdir, filename)
        baifile = bamfile+".bai"
        
        # check if index file exists
        if not os.path.exists(baifile):
            print("No index file found for %s, indexing now"%(bamfile))
            os.system(options.sambamba + " index " + bamfile)
            print("Indexing performed")
        
        # generate Telomere reads file name   
        telofile = os.path.join(options.outdir, filename.replace(".bam","_TelomericReads.sam"))
        
        # generate Telomere reads file
        print(bamfile,telofile)
        counts = count_telomeric_reads(bamfile, telofile)
        output.write('\t'.join([filename.split("_")[0],counts[0], ((counts[1],counts[1]/counts[0])*1000)])+'\n')
    
    output.close()

# -------------------------------------------------
# Execute program
# -------------------------------------------------
print("Starting analysis")
main()
print("DONE")
# -------------------------------------------------
