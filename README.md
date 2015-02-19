Small tools for genetics research 
========

See the [Small tools Manifesto](https://github.com/pjotrp/bioinformatics)  

---
### Estimate telomere lengths in BAM files
1. count_temolmeres.py:  
   Extracts reads containing telomeric sequences from BAM files using [SAMBAMBA](http://lomereiter.github.io/sambamba/).  
   Stores these in a sam file for further processing.  
   Currently performs a simple line count and provides the total read count of the BAM to generate a normalised telomere count.  
  
---
### Annotate a VCF with CADD scores
1. Annotate_CADD_Scores_In_VCF.py:  
   Annoates variants in a VCF file with [CADD scores](http://cadd.gs.washington.edu/score) using the [precomputed](http://cadd.gs.washington.edu/download) files provided by  the University of Washington.  
   It uses multiprocessing to distribute tabix lookup commands over different processes.  
  
---  
Hope they help.  
If you have suggestions or improvements let me know.  
