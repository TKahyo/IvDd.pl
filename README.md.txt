#IvDd.pl
  
##Description  
#####IvDd.pl is a Perl script to estimate insertional variations of LTR-retrotransposons from the DGV dataset. 
  
##Requirement
#####1. Tools  
* BEDtools (intersectBed): [Quinlan Lab at the University of Utah](http://bedtools.readthedocs.io/en/latest/)	v2.25.0 or later. -F option of intersectBed is required. 
#####2. Data sets  
* DGV Variants (dgv.txt): Downloadable from [the Database of Genomic Variants (DGV)](http://dgv.tcag.ca/dgv/app/about?ref=GRCh37/hg19). Containing 20 columns.  

* ORF data (RM\_ORF\_list.bed): Annotation data of open reading frames (ORFs). Downloadable via RepeatMasker track from the Table Browser in  [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html). If LTR5_Hs elements are analyzed, HERVK-int should be selected. Containing 6 columns.  
#  
    chr1    12840257    12845090    HERVK-int    27617    -  
    chr1    13459295    13460029    HERVK-int    4505     +  
        ...  

* LTR list data (RM\_LTR\_list.bed): Annotation data of LTRs. Downloadable via RepeatMasker track from the Table Browser in  [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html). Containing 6 columns.  
#  
    chr1    10485451    10485712    LTR5_Hs    7438    +    
    chr1    10486981    10487683    LTR5_Hs    7438    +
        ...  
  
##Usage
    perl IvDd.pl -d dgv.txt -t RM_LTR_list.bed -i RM_ORF_list.bed [other_options]  >result

#####Options
-h|--help	print help.  

*-b    Full path of intersectBed

*-d    DGV dataset (.txt)

-del    Maximum of deletion length parameter (default 1500)

-f    intersectBed -f option (default 0.9)

*-t|-i    LTR\_list.bed (-t) and ORF\_list.bed (-i) from RepeatMasker data

*required options  

Results are printed on the standard output.  
Table captions of the result: [6 captions of RM\_LTR\_list.bed] [Chr, Start, End, ID, CNV, Type and Study of dgv.txt] [6 captions of RM\_ORF\_list.bed]  


