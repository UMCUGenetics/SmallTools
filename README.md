Small tools for genetics research 
========

See the [Small tools Manifesto](https://github.com/pjotrp/bioinformatics)  

---
### Estimate telomere lengths in BAM files
1. Count_Telomeric_Sequence_Reads.py:  
   Extracts reads containing telomeric sequences from BAM files using [SAMBAMBA](http://lomereiter.github.io/sambamba/).
   Number of telomeric 6-mers that need to matched can be set.
   Uses multiprocessing to distribute view/count commands over different processes.
   Stores matching telomeric reads in a SAM file for further processing.  
   Uses total read count of the SAM and BAM files to generate a normalised telomere fraction.  
  
---
### Annotate a VCF with CADD scores
1. Annotate_CADD_Scores_In_VCF.py:  
   Annoates variants in a VCF file with [CADD scores](http://cadd.gs.washington.edu/score) using the [precomputed](http://cadd.gs.washington.edu/download) files provided by  the University of Washington.  
   Uses multiprocessing to distribute tabix lookup commands over different processes.  
  
---  
Hope they help.  
If you have suggestions or improvements let me know.  
