#!/usr/bin/python2.7

import pysam
import os
sambamba = "/path/to/sambamba/sambamba_vx.x"
workdir = "/path/to/working/directory/"
searchdir = "/path/to/search/bams/"

for root, dirs, files in os.walk(searchdir):
    for filename in files:
        # select bam files
        if filename.endswith(".bam"):
            print(filename)
            telofile = workdir+filename.replace(".bam","_TelomericReads.sam")
            bamfile = os.path.join(root, filename)

            if not os.path.exists(telofile):
                # extract telomeric reads and write to file
                os.system(sambamba + " view " + bamfile + " | LC_ALL=C grep -E \"" + "TTAGGG"*15 +"|"+ "CCCTAA"*15 + "\"" + " > " + telofile + " & ")
            
            else:
                # count total number of reads
                total_rc = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamfile) ])
                # count number of telomeric reads
                telomere_rc = os.system("wc -l "+telofile)
                # print counts
                print(filename, total_rc, telomere_rc)
