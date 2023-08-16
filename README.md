# CADD_summary

IMPORTANT:  
-CADD scores for homozygous alleles are not counted twice, but just once.  
-Multiple entries for one positions are discarded, only one is taken into account  

#### Usage: python script.py --vep file.vep --tfile tfile --output cadd_sum_output

python sum_v4.py -h  
usage: sum_v4.py [-h] --vep VEP --tfile TFILE [--CADD CADD] [--threshold THRESHOLD] [--gene GENE] [--start_bp START_BP] [--end_bp END_BP] [--denovo DENOVO] --output OUTPUT  

CADD file summary

optional arguments:  
  -h, --help            show this help message and exit  
  --vep VEP             vep file with CADD scores  
  --tfile TFILE         Prefix of tped and tfam, both files will be read  
  --CADD CADD           Either CADD_PHRED or CADD_RAW score will be estimated, default CADD_PHRED  
  --threshold THRESHOLD CADD threshold  
  --gene GENE           File with gene list to estimate CADD from  
  --start_bp START_BP   Int with start position in bp  
  --end_bp END_BP       Int with end position in bp  
  --denovo DENOVO       Estimate de novo mutations, WORK IN PROGRESS DO NOT USE  
  --output OUTPUT       output name file  
