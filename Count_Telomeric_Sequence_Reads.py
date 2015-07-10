#!/usr/bin/python2.7

import pysam
import glob
import os

from time import time 
from time import sleep
import subprocess
import multiprocessing as mp

from optparse import OptionParser
# -------------------------------------------------
parser = OptionParser()
parser.add_option("--sambamba",	 dest="sambamba", help="Path to sambamba/samtools executable", default="sambamba")
parser.add_option("--bamdir",	 dest="bamdir",   help="Path to directory containing BAM files", default=False)
parser.add_option("--outdir",	 dest="outdir",   help="Path to directory to write output to", default="./telomeres/")
parser.add_option("--repsize",   dest="repsize",  help="Number of required matching 6mers (TTAGGG)", default=10)
parser.add_option("--s",         dest="nr_samples",help="Number of Samples to analyse simulatiously", default=8)
parser.add_option("--t",         dest="nr_cpus",   help="Number of CPUs to use per sample", default=4)

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

def count_telomeric_reads(bamfile, q):
    # generate Telomere reads file name  
    #print("---- Processing BAM file: "+bamfile)
    telofile = bamfile.replace(options.bamdir,options.outdir).replace(".bam","_TelomericReads.sam")
    
    # check if the file was already generated
    if not os.path.exists(telofile):
        # extract telomeric reads and write to file
        cmd = options.sambamba + " view " + bamfile + " -t "+ options.nr_cpus +" | LC_ALL=C grep -E \"" + "TTAGGG"*options.repsize +"|"+ "CCCTAA"*options.repsize + "\"" + " > " + telofile
        #print("++++ Generating SAM file: "+telofile)
        os.system(cmd)

    # count total number of reads
    total_rc = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfile) ])
    
    sleep(1)
    
    telomere_rc = 0
    if os.path.exists(telofile):
        # count number of telomeric reads by line count
        telomere_rc = sum(1 for line in open(telofile,'r'))
        print("DONE Processing BAM file: "+bamfile)
    else:
        print("Something went wrong with BAM file: "+bamfile)
        
    # return results
    result = '\t'.join([bamfile.split("/")[-1].split("_")[0],str(total_rc), str(telomere_rc), str((telomere_rc/(total_rc*1.0))*100000.0)])+'\n'
    q.put(result)
    return(result)

    
# -------------------------------------------------
def listener(q):
    '''listens for messages on the q, writes to file. '''
    #sys.stdout.write('Starting listener\n')
    
    f = open(os.path.join(options.outdir, "TelomereCounts.txt"), 'wb')
    f.write('\t'.join(["#Sample","TotalReads","TelomericReads","NormalisedFraction"])+'\n')
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
        
        f.write(m)
        f.flush()
    f.close()
# -------------------------------------------------

def main():
    currtime = time()

    #Init Manager queue
    manager = mp.Manager()
    q = manager.Queue()

    # Init worker pool
    pool = mp.Pool(int(options.nr_samples))

    #Init Listener
    watcher = pool.apply_async(listener, (q,))

    bamfiles = glob.glob(os.path.join(options.bamdir, "*.bam"))
    jobs = []

    #fire off workers
    for bamfile in bamfiles:
        
        baifile = bamfile+".bai"
        # check if index file exists
        if not os.path.exists(baifile):
            print("No index file found for %s, indexing now"%(bamfile))
            subprocess.call(options.sambamba, " index " + bamfile)

        job = pool.apply_async(count_telomeric_reads, (bamfile, q))
        jobs.append(job)        

            
    for job in jobs:
        job.get()

    # now we are done, kill the listener
    q.put('kill')
    
    pool.close()
    pool.join()

    print 'time elapsed:', time() - currtime

# -------------------------------------------------

print("Starting Analysis")

if __name__ == '__main__':
    # check specified options
    if check_arguments():
        main()
    else:
        print("Error in provided arguments")

print("DONE")

# -------------------------------------------------
