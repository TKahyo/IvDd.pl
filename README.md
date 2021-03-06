#IvDd.pl
  
##Description  
#####IvDd.pl is a Perl script to estimate insertional variations of LTR-retrotransposons from the DGV dataset. 
  
##Requirement
### Tools  
 [1] BEDtools (intersectBed): [Quinlan Lab at the University of Utah](http://bedtools.readthedocs.io/en/latest/)	v2.25.0 or later. -F option of intersectBed is required.  

### Data sets  
 [1] DGV Variants (dgv.txt): Downloadable from [the Database of Genomic Variants (DGV)](http://dgv.tcag.ca/dgv/app/about?ref=GRCh37/hg19). Containing 20 columns.  

 [2] Internal proviral sequence data (RM\_ORF\_list.bed): Annotation data of the internal proviral sequences. Downloadable via RepeatMasker track from the Table Browser in  [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html). In the case of HML-2, HERVK-int should be included. Containing 6 columns.  

#  
    chr1    12840257    12845090    HERVK-int    27617    -  
    chr1    13459295    13460029    HERVK-int    4505     +  
        ...  

 [3] LTR list data (RM\_LTR\_list.bed): Annotation data of LTRs. Downloadable via RepeatMasker track from the Table Browser in  [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html). In the case of HML-2, LTR5_Hs should be included. Containing 6 columns.  

#  
    chr1    10485451    10485712    LTR5_Hs    7438    +    
    chr1    10486981    10487683    LTR5_Hs    7438    +  
        ...  
  
##Usage
    perl IvDd.pl -d dgv.txt -t RM_LTR_list.bed -i RM_int_list.bed [other_options]  >result

#####Options
-h|--help　　print help

-b\*　　/path/intersectBed

-d\*　　DGV dataset (.txt)

-del　　Maximum of deletion length parameter (default 1500)

-f　　intersectBed -f option (default 0.9)

-t\*|-i\*　　LTR\_list.bed (-t) and int\_list.bed (-i) from RepeatMasker data

*required  

Results are printed on the standard output.  
Table captions of the result: [6 captions of RM\_LTR\_list.bed] [Chr, Start, End, ID, CNV, Type and Study of dgv.txt] [6 captions of RM\_int\_list.bed (if any)]  
