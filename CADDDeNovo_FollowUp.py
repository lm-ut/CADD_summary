#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import re
import sys

## Usage
#python script.py chr22_Affected_IndexID.csv MotherID FatherID

# Init. File
FILE = sys.argv[1]

# Extract IndexID and CHR
ID = ''.join(re.findall(r"G\w*",FILE))
CHR = ''.join(re.findall(r"chr\d+",FILE))

# Init. MotherID and FatherID
ID1 = sys.argv[2]
ID2 = sys.argv[3]

# Read file (pandas)
cadd_file = pd.read_csv(FILE, sep = '\t', index_col=0)

# Extract only entries where allele is unique to Index
DENOVO = cadd_file[((cadd_file['RES_' + ID + '_1']).notna() |
                         (cadd_file['RES_' + ID + '_2']).notna()) &
                        (cadd_file['RES_' + ID1 + '_1']).isna() &
                        (cadd_file['RES_' + ID1 + '_2']).isna() &
                        (cadd_file['RES_' + ID2 + '_1']).isna() &
                        (cadd_file['RES_' + ID2 + '_2']).isna()]


# Extract only col of interest
subset_cols = ['RES_' + ID + '_1','RES_' + ID + '_2','Gene','CADD_PHRED','Allele', 
'POS']
subset_denovo = DENOVO[subset_cols]

# If there are de novo, create new file,
# Otherwise state that there are no de novo in chr X

if len(subset_denovo)>1:
    subset_denovo.to_csv('DENOVO_'+FILE, sep = '\t')
else:
    print('NO DENOVO FOUND IN: '+CHR)
