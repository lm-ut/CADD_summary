# CADD_summary

IMPORTANT:  
-CADD scores for homozygous alleles are not counted twice, just once.  
-Multiple entries for one position are discarded, only one is taken into account  

#### Basic usage: python script.py --vep file.vep --tfile tfile --output cadd_sum_output  

python script.py -h  
usage: script.py [-h] --vep VEP --tfile TFILE [--CADD CADD] [--threshold THRESHOLD] [--gene GENE] [--start_bp START_BP] [--end_bp END_BP] [--denovo DENOVO] --output OUTPUT  
<br>
[optional] arguments:  
  -h, --help            show this help message and exit  
  <strong>--vep</strong> VEP:   vep file with CADD scores  
  <strong>--tfile</strong> TFILE:         Prefix of tped and tfam, both files will be read  
  [<strong>--CADD</strong>] CADD_PHRED/CADD_RAW:           Either CADD_PHRED or CADD_RAW score will be estimated, default CADD_PHRED  
  [<strong>--threshold</strong>] THRESHOLD: CADD threshold, Int  
  [<strong>--gene</strong>] GENE_LIST:           File with gene list to estimate CADD from  
  [<strong>--start_bp</strong>] START_BP:   Int with start position in bp  
  [<strong>--end_bp</strong>] END_BP:       Int with end position in bp  
  [<strong>--denovo</strong>] DENOVO:       Estimate de novo mutations, WORK IN PROGRESS DO NOT USE  
  --output OUTPUT       output name file  
