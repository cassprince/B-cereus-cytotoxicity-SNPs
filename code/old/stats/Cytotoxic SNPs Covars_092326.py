# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:19:16 2023

@author: cassp
"""

#Import necessary packages.
import pandas as pd
import numpy as np
from Bio import AlignIO
from statsmodels.sandbox.stats.multicomp import multipletests
import statsmodels.api as sm
import os
import re

#Define the function that applies the logistic regression model to each site in the alignment.
def modelTox(nTP, cutoff1, cutoff2):
    #Create empty dataframe for output.
    dfOutput = pd.DataFrame(0, index = range(0), columns = ["gene", "position", "nuc","name", "total_w", "total_wo", "avg_tox_w", "avg_tox_wo", "sd_tox_w", "sd_tox_wo", "lr_p-val_tox", "lr_p-val_phylo"])
    position = 0
    while position < len(alignment[0]):
        nuc_hyph = alignment[: , position]
    
        #Remove hyphens so they don't have an impact as "negatives".
        cytotoxicityShort = cytotoxicity[(nuc_hyph!="-")]
        testLR = data[(nuc_hyph!="-")]    
        nuc = nuc_hyph[(nuc_hyph!="-")]
    
        #If there is sufficient variation between sequences to correlate cytotox and nucleotide (mostly because f scores were 0 in most cases), pursue the site as possibly containing a SNP. The cutoffs for # of sequences with the given nucleotide in the site are between 20 and 80%. This prevents the perfectly conserved nucleotides from being investigated, as they wouldn't be informative. Also if a nucleotide is very rarely found, it won't be investigated. We want SNPs that aren't extremely rare. Those would likely be unhelpful in application.
        if list(nuc).count(nTP) <= cutoff1 * len(nuc) and list(nuc).count(nTP) >= cutoff2 * len(nuc): 
            y = np.where((nuc == nTP), 1, 0) #Where a given nucleotide exists in the position, the array is assigned a 1. Where there is not, the array is assigned a 0.
            withNTP = cytotoxicityShort[nuc == nTP]
            withoutNTP = cytotoxicityShort[nuc != nTP]
    
            testY = pd.DataFrame(y, dtype = "int")
            testLR = testLR[['cytotoxicity', 'panC_group']]
            testLR = testLR.reset_index(drop=True)
            
            log_reg = sm.Logit(testY, testLR).fit() #Fit the logistic regression model.
    
            #Record general data.
            df_new_row = pd.DataFrame({"gene":[gene], 
                                       "position":[position + 1], 
                                       "nuc":[nTP], 
                                       "name":[gene + '_SNP_' + nTP + '_' + str(position+1)], 
                                       "total_w":[len(withNTP)], "total_wo":[len(withoutNTP)], 
                                       "avg_tox_w":[np.mean(withNTP)], 
                                       "avg_tox_wo":[np.mean(withoutNTP)], 
                                       "sd_tox_w":[np.std(withNTP)], 
                                       "sd_tox_wo":[np.std(withoutNTP)], 
                                       "median_tox_w":[np.median(withNTP)],
                                       "median_tox_wo":[np.median(withoutNTP)],
                                       "lr_p-val_tox":[log_reg.pvalues['cytotoxicity']], 
                                       "lr_p-val_phylo":[log_reg.pvalues['panC_group']]})
            dfOutput = pd.concat([dfOutput, df_new_row], ignore_index = True)
            
            #For the SNP presence or absence matrix, always identify "having the SNP" as having the more cytotoxic SNP.
            if np.mean(withNTP) > np.mean(withoutNTP):
                data.loc[nuc_hyph == nTP, (gene + '_SNP_' + nTP + '_' + str(position+1))] = 1
                data.loc[nuc_hyph != nTP, (gene + '_SNP_' + nTP + '_' + str(position+1))] = 0 
                data.loc[nuc_hyph == "-", (gene + '_SNP_' + nTP + '_' + str(position+1))] = np.nan
            if np.mean(withNTP) < np.mean(withoutNTP):
                data.loc[nuc_hyph == nTP, (gene + '_SNP_' + nTP + '_' + str(position+1))] = 0 
                data.loc[nuc_hyph != nTP, (gene + '_SNP_' + nTP + '_' + str(position+1))] = 1 
                data.loc[nuc_hyph == "-", (gene + '_SNP_' + nTP + '_' + str(position+1))] = np.nan
    
        position += 1
    dfOutput = dfOutput[dfOutput['nuc'] != 0]
    return(dfOutput)

##########------------------------------


#Set working directiory and load in alignment and data.
os.chdir(r'C:\Users\cassp\OneDrive\Documents\GitHub\B-cereus-cytotoxicity-SNPs\blast_results\filtered_megablast_qc50_090126\alignments')

mSheet = pd.read_excel(r"C:\Users\cassp\OneDrive - Cornell University\Biomarkers paper\Mastersheet_082026.xlsx") #import file with cytotoxicity data (same order as fasta)

mSheet['panC_group'] = mSheet['Adjusted_panC_Group(predicted_species)'].replace({'[*]': ''}, regex = True)

mSheet['panC_group'] = mSheet['panC_group'].replace({'Group_clarus' : 0,
                                                     'Group_I(pseudomycoides)' : 1, 
                                                     'Group_II(mosaicus/luti)': 2, 
                                                     'Group_III(mosaicus)': 3, 
                                                     'Group_IV(cereus_sensu_stricto)': 4, 
                                                     'Group_V(toyonensis)': 5, 
                                                     'Group_VI(mycoides/paramycoides)': 6, 
                                                     'Group_VII(cytotoxicus)' : 7, 
                                                     'Group_VIII(mycoides)' : 8})

mSheet['cytotoxicity'] = mSheet['Cytotoxicity ( >0.7 is cytotoxic)']

metadata = mSheet.loc[:, ["Isolate", "cytotoxicity", "panC_group"]]
nucleotides = ["A", "T", "G", "C"]
proteins = ["A", "R", "D", "N", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

file = "nheA_up_align.fasta"
gene = file.replace("_align.fasta", "")
print(gene)

alignmentSeq = AlignIO.read(file, "fasta") #Import alignment file of sequences in fasta format.
alignment = np.array([list(rec) for rec in alignmentSeq])

IDs = []
for record in alignmentSeq:
    #Change names to match mastersheet. If you downloaded the sequences of hits from BTyper3, this should work. If you procured sequences by other means, you may need to change this part so you can get matching names with the mastersheet.
    
    recordID = record.id
    PS_ID = recordID.split("_", 1)[0]
    row = mSheet[mSheet["Isolate"].str.contains(PS_ID)]
    if len(row) != 0:
        IDs.append(PS_ID)

dfIDs = pd.DataFrame(IDs, columns = ["Isolate"]) 

data = dfIDs.merge(metadata)
#data = data.sort_values('Isolate')
cytotoxicity = data['cytotoxicity'].values
cytotoxicity = cytotoxicity.reshape(-1,1)

df_full = pd.DataFrame(0, index = range(0), columns = ["gene", "position", "nuc","name", "total_w", "total_wo", "avg_tox_w", "avg_tox_wo", "sd_tox_w", "sd_tox_wo", "lr_p-val_tox", "lr_p-val_phylo"])

cutoff1 = 0.95
cutoff2 = 0.05

# Run the model on each nucleotide or amino acid. Append data for each nucleotide to same dataframe (df_full). 
if "prot" in file:
    for i in proteins:
        df = modelTox(i, cutoff1, cutoff2)
        df_full = pd.concat([df,df_full], ignore_index = True)
        
if "up" in file:
    for i in nucleotides:
        df = modelTox(i, cutoff1, cutoff2)
        df_full = pd.concat([df,df_full], ignore_index = True)

#Bonferroni correction of p-values.
p_adjusted_cyt = multipletests(df_full["lr_p-val_tox"], alpha=0.05, method='bonferroni')
p_adjusted_phylo = multipletests(df_full["lr_p-val_phylo"], alpha=0.05, method='bonferroni')
corrected_cyt = pd.DataFrame(p_adjusted_cyt[1], columns = ["lr_p-val_tox_bonf"])
corrected_phylo = pd.DataFrame(p_adjusted_phylo[1], columns = ["lr_p-val_phylo_bonf"])

df_final = pd.concat([df_full, corrected_cyt, corrected_phylo], axis = 1)


df_final_filt = df_final[df_final['lr_p-val_tox_bonf'] < 0.05]
data_filt = data[data.columns.intersection(list(df_final_filt['name']) + ["cytotoxicity", "panC_group", "Isolate"])]

#Save the output files as .csv files. dfOutput has the general info about the SNPs. data is the SNP presence/absence matrix.
os.chdir(r'C:\Users\cassp\OneDrive\Documents\GitHub\B-cereus-cytotoxicity-SNPs/data')
df_final_filt.to_csv(f"{gene}_logreg_092326.csv", index = False)
data_filt.to_csv(f"{gene}_snps_092326.csv", index = False)

df_final.to_csv(f"{gene}_unfilt_logreg_092326.csv", index = False)
data.to_csv(f"{gene}_unfilt_snps_092326.csv", index = False)
