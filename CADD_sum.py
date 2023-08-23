#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from itertools import compress
import io
import re
import argparse
import os
import sys


parser = argparse.ArgumentParser(description='CADD file summary')
parser.add_argument('--vep',help='vep file with CADD scores',required=True,default='1')
parser.add_argument('--tfile',help='Prefix of tped and tfam, both files will be read',required=True,default='1')
parser.add_argument('--CADD',help='Either CADD_PHRED or CADD_RAW score will be estimated, default CADD_PHRED',required=False,default='CADD_PHRED')
parser.add_argument('--threshold',help='CADD threshold required',required=False,default='10')
parser.add_argument('--gene',help='File with gene list to estimate CADD from',required=False,default='')
parser.add_argument('--start_bp',help='Int with start position in bp',required=False,default='')
parser.add_argument('--end_bp', help='Int with end position in bp',required=False,default='')
parser.add_argument('--denovo', help='Estimate de novo mutations, WORK IN PROGRESS DO NOT USE', required=False,default='NO')
parser.add_argument('--output', help='output name file', required=True,default='output')

#parser.add_argument('--plot',help='Want to plot? YES or NO, default NO',required=False,default='NO')

args = vars(parser.parse_args())



# Default
CADD_FILE = args['vep']
TPED = args['tfile']+'.tped'
TFAM = args['tfile']+'.tfam'
OUTPUT = args['output']

# Not Default
TARGET = args['CADD']
THRESHOLD = args['threshold']
GENE_LIST = args['gene']
START_BP = args['start_bp']
END_BP = args['end_bp']

# Work in Progress options
DENOVO = args['denovo']



def log_file(OUTPUT):
    with open(OUTPUT+'.log','w') as logfile:
        logfile.write('CADD file used: '+str(CADD_FILE)+'\n')
        logfile.write('plink tfile used: '+str(TPED)+', '+str(TFAM)+'\n')
        logfile.write('output name used: '+str(OUTPUT)+'\n')
        logfile.write('threshold used: '+str(THRESHOLD)+'\n')
        if START_BP != '' and END_BP != '':
            logfile.write('Bp range used: '+str(START_BP)+'-'+str(END_BP)+'\n')
        else:
            logfile.write('Bp range used: none'+'\n')
        if GENE_LIST != '':
            logfile.write('Gene file used: '+ str(GENE_LIST)+'\n')
        else:
            logfile.write('Gene file used: none'+'\n')       

## Assumes that VEP file has info for only one chromosome
def get_chr(cadd_file):
    Chromosome = cadd_file['CHR'].drop_duplicates().item()
    return Chromosome


def return_IDs(tfam):
    tfam = pd.read_csv(TFAM,sep=' ', names=['POP','ID','A','B','C','D'])
    IDs = tfam['ID']
    return IDs


### Retrieve only the non # rows from the CADD output

def clean_cadd_input(cadd_file):
    all_rows_short_output = []

    with open(cadd_file, 'r') as f:
        VEP_file = [l for l in f if not l.startswith('#')]
        for row in VEP_file:
            row_list = re.split(r';|\t', row)
            result = [TARGET in word for word in row_list]
            subset = list(compress(row_list,result))
            short_output = row_list[0:4] + subset
            all_rows_short_output.append(short_output)
        if len(all_rows_short_output) > 0:
            return all_rows_short_output
        else:
            print("Something is wrong at fun clean_cadd_input(). Your file has length: "+ str(len(cadd_file)))
            print("Stopping program")
            sys.exit()  


### extract columns of interest ##

def extract_columns_cadd(clean_cadd_file):

    input_for_summary = pd.DataFrame(clean_cadd_file, columns = ['Format','CHR:POS','Allele','Gene',TARGET])
    regex = TARGET+'='
    input_for_summary[TARGET] = input_for_summary[TARGET].str.replace(regex, '', regex=True)
    input_for_summary[TARGET] = input_for_summary[TARGET].str.replace(r'\n', '', regex=True).astype('float')

    ## Separating CHR and POS
    input_for_summary[['CHR', 'POS']] = input_for_summary["CHR:POS"].apply(lambda x: pd.Series(str(x).split(":")))
    
    ## Cleaning POS
    #input_for_summary[['POS', 'POS2']] = input_for_summary["POS"].apply(lambda x: pd.Series(str(x).split("-")))
    #general_CADD = input_for_summary.drop(columns=['POS2'])

    general_CADD = input_for_summary
    if len(general_CADD) > 0:
        return general_CADD
    else:
        print("Something is wrong at fun extract_columns_cadd(). Your file has length: "+ str(len(cadd_file)))
        print("Stopping program")
        sys.exit()



## Subset based on GENE_LIST parameter ##

def extract_gene_list(GENE_LIST,general_CADD):

    if GENE_LIST != "":
        print('Checking gene_list')
        gene_list = pd.read_csv(GENE_LIST, names = ['Gene'])
        joined_on_gene_list = pd.merge(general_CADD, gene_list, on=['Gene'], how='right')
        CADD_filtered_genes = joined_on_gene_list
        if len(CADD_filtered_genes) < 1:
            print("No genes found")
            sys.exit()
        else:
            return CADD_filtered_genes
    else:
        print("No gene list given, selecting all available entires")
        CADD_filtered_genes = general_CADD
        return CADD_filtered_genes



## Select entries within the bp range ##

def bp_range(cadd_file, start_bp, end_bp):
    
    if (start_bp != '') and (end_bp != ''):
        print('Checking positions from '+str(start_bp)+" to "+str(end_bp))
        cadd_file["POS"] = pd.to_numeric(cadd_file["POS"])
        df = cadd_file[(cadd_file['POS'] >= int(start_bp)) & (cadd_file['POS'] <= int(end_bp))]
        if len(df) > 0:
            return df
        else:
            print("Something is wrong at fun bp_range(). Your file has length: "+ str(len(cadd_file)))
            print("Stopping program")
            sys.exit()
    else:
        print("No bp range given, selecting all available SNPs")
        return cadd_file


def sort_and_clean_cadd(cadd_file, cadd_type, threshold):  

## Sorted file, to remove duplicate entries, 
    sorted_cadd = cadd_file.drop_duplicates(subset=['CHR:POS','Allele'])

## Absent alleles (-) and indel needs to be removed
    cadd_usable_snps = sorted_cadd[sorted_cadd['Allele'].str.match('^[A]$|^[C]$|^[G]$|^[T]$')== True]
    
## Apply threshold
    cadd_threshold = cadd_usable_snps[cadd_usable_snps[cadd_type] > float(threshold)]
    if len(cadd_threshold) < 1:
        print("No entries pass the threshold")
        sys.exit()
    else:
        return cadd_threshold



## Assign label 'Affected' if CADD variant is present ##

def assign_affected(sharedpos_caddfile):
    IDs = return_IDs(TFAM)

    #df.loc[df['column name'] condition, 'new column name'] = 'value if condition is met'
    
    for i in range(0,len(IDs)):
        sharedpos_caddfile.loc[sharedpos_caddfile['Allele'] == sharedpos_caddfile[IDs[i]+'_1'], 'RES_'+IDs[i]+'_1'] = 'Affected'
        sharedpos_caddfile.loc[sharedpos_caddfile['Allele'] == sharedpos_caddfile[IDs[i]+'_2'], 'RES_'+IDs[i]+'_2'] = 'Affected'

    subset_shared_pos_ordered = sharedpos_caddfile.T.groupby(level=0).first().T
    
    return subset_shared_pos_ordered



## Hard-coded for family based analyses

def de_novo(cadd_file):
    print('Looking for denovo mutations')
    IDs = return_IDs(TFAM)
    denovo0 = cadd_file[((cadd_file['RES_'+IDs[0]+'_1']).notna() | 
                                          (cadd_file['RES_'+IDs[0]+'_2']).notna()) & 
                                          (cadd_file['RES_'+IDs[1]+'_1']).isna() &
                                          (cadd_file['RES_'+IDs[1]+'_2']).isna() &
                                          (cadd_file['RES_'+IDs[2]+'_1']).isna() &
                                          (cadd_file['RES_'+IDs[2]+'_2']).isna()]

    denovo1 = cadd_file[((cadd_file['RES_'+IDs[1]+'_1']).notna() | 
                                          (cadd_file['RES_'+IDs[1]+'_2']).notna()) &
                                          (cadd_file['RES_'+IDs[0]+'_1']).isna() &
                                          (cadd_file['RES_'+IDs[0]+'_2']).isna() &
                                          (cadd_file['RES_'+IDs[2]+'_1']).isna() &
                                          (cadd_file['RES_'+IDs[2]+'_2']).isna()]

    denovo2 = cadd_file[((cadd_file['RES_'+IDs[2]+'_1']).notna() | 
                                          (cadd_file['RES_'+IDs[2]+'_2']).notna()) & 
                                          (cadd_file['RES_'+IDs[1]+'_1']).isna() &
                                          (cadd_file['RES_'+IDs[1]+'_2']).isna() &
                                          (cadd_file['RES_'+IDs[0]+'_1']).isna() &
                                          (cadd_file['RES_'+IDs[0]+'_2']).isna()]

    denovo_info = pd.DataFrame(columns = ['ID','DeNovo_N'])

    denovo_info.at[0,'ID'] = IDs[0]
    denovo_info.at[1,'ID'] = IDs[1]
    denovo_info.at[2,'ID'] = IDs[2]

    denovo_info.at[0,'DeNovo_N'] = len(denovo0)
    denovo_info.at[1,'DeNovo_N'] = len(denovo1)
    denovo_info.at[2,'DeNovo_N'] = len(denovo2)
    return denovo_info



## Summary ##    

def get_summary(cadd_file):
    if len(cadd_file) >0:
        CHROMOSOME = get_chr(cadd_file)
        CADD_SUM = cadd_file[TARGET].sum()
        CADD_AVG = cadd_file[TARGET].mean()
        CADD_MIN = cadd_file[TARGET].min()
        CADD_MAX = cadd_file[TARGET].max()
        CADD_N = len(cadd_file[TARGET])


        short_summary = pd.DataFrame(columns = ['CHR','CADD_SUM','CADD_AVG','CADD_MIN','CADD_MAX','N'])
        short_summary.loc[len(short_summary)] = {'CHR': CHROMOSOME, 'CADD_SUM': CADD_SUM, 'CADD_AVG': CADD_AVG, 'CADD_MIN': CADD_MIN, 'CADD_MAX': CADD_MAX, 'N': CADD_N}
        summary_final = short_summary.rename(index = {0: TARGET})
        return summary_final

    else:
        print("File empty")
        sys.exit()



## FIND positions in TPED that overlap with the CADD score file ##

def overlap_tped_pos(tped_file, cadd_file):

    IDs = return_IDs(TFAM)

## Duplicating the samples so that the IDs are haplotype-based, IDs are now IDs_1 and IDs_2 

    IDs_1 = pd.Series([ID for ID in IDs + "_1"], name="ID")
    IDs_2 = pd.Series([ID for ID in IDs + "_2"], name="ID")

    IDs_duplicated = pd.merge(IDs_1, IDs_2, right_index=True, left_index=True)
    One_ColIDs = IDs_duplicated.apply(pd.Series).stack().reset_index(drop=True)
    tmp = One_ColIDs.tolist()

## Setting col names for tped
    colnames = ['CHR','rs','cM','POS']
    colnames.extend(tmp)

## Reading tped and finding overlapping positions
    
    tped = pd.read_csv(tped_file, sep=' ', names=colnames)
    #print(tped)
    cleaned_tped = tped.drop(['CHR','rs','cM'],axis=1)
    cleaned_tped['POS'] = cleaned_tped['POS'].astype(int)

    cadd_file['POS']=cadd_file['POS'].copy().astype(int)
    cleaned_general_CADD = cadd_file[["POS","CHR","Allele","Gene","CADD_PHRED"]]
    
    
    shared_pos = cleaned_tped.merge(cleaned_general_CADD, on=["POS"])
    if len(shared_pos) < 1:
        print("No positions found in tped")
        sys.exit()
    else:
        return shared_pos



## Summary ind-based ##

def ind_based_summary(cadd_file):

    summary_fam = pd.DataFrame(columns = ['CHR','CADD_SUM','CADD_AVG','CADD_MIN','CADD_MAX','N','ID'])
    IDs = return_IDs(TFAM)

    for i in range(0,len(IDs)):
          
        RES_tmp = cadd_file[(cadd_file['RES_'+IDs[i]+'_1'] == 'Affected') | (cadd_file['RES_'+IDs[i]+'_2'] == 'Affected')].copy()
        
        if len(RES_tmp) > 1:
            RES_tmp.dropna(axis=1, how='all', inplace=True)
            RES_tmp_sum = get_summary(RES_tmp)
            RES_tmp_sum['ID'] = IDs[i] 
            
            summary_fam = pd.concat([summary_fam,RES_tmp_sum])
            #summary_fam = summary_fam.append(RES_tmp_sum)
            #print(summary_fam)

        else:
            print(IDs[i]," has no overlapping alleles")
            summary_fam.loc[len(summary_fam)] = {'CHR': 'NaN', 'CADD_SUM': 'NaN', 'CADD_AVG': 'NaN', 'CADD_MIN': 'NaN', 'CADD_MAX': 'NaN','N': 'NaN', 'ID': IDs[i]}
       
    #summary_fam = pd.concat([summary_fam,RES_tmp_sum])
    #print(summary_fam)

    return summary_fam



def main():
    print("You're now using a script that creates summary statistics based on CADD values")
    print("")
    print("I'll now clean the data and check if the cleaning steps went smoothly")
    print("")
    
    IDs = return_IDs(TFAM)
    
    # Cleaning input file
    ## Extracting non-# lines and columns of interes
    all_rows_short_output = clean_cadd_input(CADD_FILE)
    general_CADD_file = extract_columns_cadd(all_rows_short_output)

    CHROMOSOME = get_chr(general_CADD_file)

    ## Extracting genes from list and bp range
    gene_list_CADD = extract_gene_list(GENE_LIST,general_CADD_file)
    
    rangebp_cadd = bp_range(gene_list_CADD, START_BP, END_BP)
    
    ## Applying threshold
    cadd_threshold = sort_and_clean_cadd(rangebp_cadd, TARGET, THRESHOLD)
       
    # Summaries
    ## General CADD summary ##
    summary_general = get_summary(cadd_threshold)
    summary_general['ID'] = 'CADD_PHRED_General'
    #summary_final.to_csv(str(cadd_threshold['CHR'][0])+'_GeneralCADD.csv', sep='\t')
    
    print("General CADD info obtained, now working on the following samples: ",IDs.tolist())
    ## Indiv-based CADD assignations ##
    shared_pos = overlap_tped_pos(TPED,cadd_threshold)
    subset_shared_pos_ordered = assign_affected(shared_pos)

    ## Get summary per each sample ##
    summary_fam = ind_based_summary(subset_shared_pos_ordered)
    
    if DENOVO == 'YES':
        denovo_info = de_novo(subset_shared_pos_ordered)

        ## Create file with summary information
        Final_Summary = pd.merge(summary_fam, denovo_info, on ='ID')
        #Final_Summary = Final_Summary.append(summary_general)
        Final_Summary = pd.concat([Final_Summary,summary_general])
    else:
        Final_Summary = summary_fam
    
    # Writing outputs 
    ## Write out info on affected allele per each family memeber
    print("I'm writing your outputs, please check for any new .csv file in your working directory")
    subset_shared_pos_ordered.to_csv(str(OUTPUT)+'_Affected.csv', sep='\t',na_rep='NULL')    
    Final_Summary.to_csv(str(OUTPUT)+'_SummaryFam.csv', sep='\t', index = False)
    log_file(OUTPUT) 
    print("DONE!")


if __name__ == "__main__":
    main()

