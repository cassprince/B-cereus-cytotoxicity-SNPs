# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:19:16 2023

@author: cassp
"""

#Import necessary packages.
import pandas as pd
import numpy as np
import re
from Bio import AlignIO
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, precision_score
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import chisquare
from scipy import stats
import statsmodels.api as sm
import os

#Set working directiory and load in alignment and data.
os.chdir(r'C:\Users\cassp\OneDrive\Documents\GitHub\B-cereus-cytotoxicity-SNPs\blast_results\filtered_megablast_qc50_090126\alignments')

file = 'nheA_up_align.fasta'
mSheet = pd.read_excel(r"C:\Users\cassp\OneDrive - Cornell University\Biomarkers paper\Mastersheet_082026.xlsx") #import file with cytotoxicity data (same order as fasta)

mSheet['panC_group'] = mSheet['Adjusted_panC_Group(predicted_species)'].replace({'Group_clarus' : 0, 'Group_I(pseudomycoides)' : 1, 'Group_II(mosaicus/luti)': 2, 'Group_II(mosaicus/luti)*' : 2, 'Group_III(mosaicus)': 3, 'Group_IV(cereus_sensu_stricto)': 4, 'Group_V(toyonensis)': 5, 'Group_V(toyonensis)*': 5, 'Group_VI(mycoides/paramycoides)': 6, 'Group_VI(mycoides/paramycoides)*' : 6, 'Group_VII(cytotoxicus)' : 7, 'Group_VIII(mycoides)' : 8})

mSheet['cytotoxicity'] = mSheet['Cytotoxicity ( >0.7 is cytotoxic)']

metadata = mSheet.loc[:, ["Isolate", "cytotoxicity", "panC_group"]]

#Create empty dataframe for output.
dfOutput = pd.DataFrame(0, index = range(100000), columns = ["gene", "position", "nuc","name", "total_w", "total_wo", "avg_tox_w", "avg_tox_wo", "sd_tox_w", "sd_tox_wo", "lr_p-val_tox", "lr_p-val_phylo"])

fileName = str(file)
gene = fileName.replace("_all_align.fasta", "")
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



#Define the function that applies the logistic regression model to each site in the alignment.
def modelTox(nTP, cutoff1, cutoff2, p_cutoff, a_cutoff):
    position = 0
    while position < len(alignment[0]):
        print(position)
        nuc_hyph = alignment[: , position]
    
        #Remove hyphens so they don't have an impact as "negatives".
        cytotoxicityShort = cytotoxicity[(nuc_hyph!="-")]
        dataWith = data[(nuc_hyph!="-")]
        testLR = data[(nuc_hyph!="-")]    
        
        nuc = nuc_hyph[(nuc_hyph!="-")]

        print(list(nuc).count(nTP) <= cutoff1 * len(nuc))
        print(list(nuc).count(nTP) >= cutoff2 * len(nuc))

        #If there is sufficient variation between sequences to correlate cytotox and nucleotide (mostly because f scores were 0 in most cases), pursue the site as possibly containing a SNP. The cutoffs for # of sequences with the given nucleotide in the site are between 20 and 80%. This prevents the perfectly conserved nucleotides from being investigated, as they wouldn't be informative. Also if a nucleotide is very rarely found, it won't be investigated. We want SNPs that aren't extremely rare. Those would likely be unhelpful in application.
        if list(nuc).count(nTP) <= cutoff1 * len(nuc) and list(nuc).count(nTP) >= cutoff2 * len(nuc): 
            y = np.where((nuc == nTP), 1, 0) #Where a given nucleotide exists in the position, the array is assigned a 1. Where there is not, the array is assigned a 0.
            withNTP = cytotoxicityShort[nuc == nTP]
            withoutNTP = cytotoxicityShort[nuc != nTP]
            
            
            
    
            testY = pd.DataFrame(y, dtype = "int")
            testLR = testLR[['cytotoxicity', 'panC_group']]
            testLR = testLR.reset_index(drop=True)
            
            
            log_reg = sm.Logit(testY, testLR).fit() #Fit the logistic regression model.
            #pred = list(map(round, log_reg.predict(testLR)))
            #print(log_reg.summary())
            #print(log_reg.pvalues)
            
            #Confusion matrix of actual gene presence vs predicted gene presence based on the LogReg sigmoid line. ***True negatives are top-left, true-positives are bottom-right.***
            #confMatrix = confusion_matrix(testY, pred).flatten()
            #print(confMatrix)
            #acc = accuracy_score(y, pred)
            #prec = precision_score(y, pred)
            #print('Accuracy score:', acc)
            #print('Precision score:', prec)
    
            #In some cases, all of the isolates were being classified as positives or negatives. This is very uninformative, so I filtered for SNPs that didn't have this problem. I also added precision and accuracy cutoffs. Currently they're set at 70%.

            #Record general data.
            dfOutput.loc[position, "gene"] = gene
            dfOutput.loc[position, "position"] = position + 1
            dfOutput.loc[position, "nuc"] = nTP
            dfOutput.loc[position, "name"] = gene + '_SNP_' + str(position+1)
            dfOutput.loc[position, "avg_tox_w"] = np.mean(withNTP)
            dfOutput.loc[position, "avg_tox_wo"] = np.mean(withoutNTP)
            dfOutput.loc[position, "sd_tox_w"] = np.std(withNTP)
            dfOutput.loc[position, "sd_tox_wo"] = np.std(withoutNTP)
            dfOutput.loc[position, "total_w"] = len(withNTP)
            dfOutput.loc[position, "total_wo"] = len(withoutNTP)
            #Record logistic regression specific data.
            #dfOutput.iloc[position, 4:8] = confMatrix
            #dfOutput.loc[position, "Accuracy"] = acc
            #dfOutput.loc[position, "Precision"] = prec
            dfOutput.loc[position, "lr_p-val_tox"] = log_reg.pvalues[0]
            dfOutput.loc[position, "lr_p-val_phylo"] = log_reg.pvalues[1]
            
            #For the SNP presence or absence matrix, always identify "having the SNP" as having the more cytotoxic SNP.
            if np.mean(withNTP) > np.mean(withoutNTP):
                data.loc[nuc_hyph == nTP, (gene + '_SNP_' + str(position+1))] = 1
                data.loc[nuc_hyph != nTP, (gene + '_SNP_' + str(position+1))] = 0 
                data.loc[nuc_hyph == "-", (gene + '_SNP_' + str(position+1))] = np.nan
            if np.mean(withNTP) < np.mean(withoutNTP):
                data.loc[nuc_hyph == nTP, (gene + '_SNP_' + str(position+1))] = 0 
                data.loc[nuc_hyph != nTP, (gene + '_SNP_' + str(position+1))] = 1 
                data.loc[nuc_hyph == "-", (gene + '_SNP_' + str(position+1))] = np.nan
    
        position += 1

cutoff1 = 1
cutoff2 = 0
p_cutoff = 0
a_cutoff = 0

#Run the model on each nucleotide. depending on your alignment, you may need to capitalize the letters. You could also run this program on amino acid alignments. You just need to change the letter to the single-letter amino acid code. Ex. "F" for phenylalanine.
modelTox("a", cutoff1, cutoff2, p_cutoff, a_cutoff)
modelTox("t", cutoff1, cutoff2, p_cutoff, a_cutoff)
modelTox("g", cutoff1, cutoff2, p_cutoff, a_cutoff)
modelTox("c", cutoff1, cutoff2, p_cutoff, a_cutoff)





#Bonferroni correction of p-values.
#p_adjusted_cyt = multipletests(dfOutput["lr_p-val_tox"], alpha=0.05, method='bonferroni')
#p_adjusted_phylo = multipletests(dfOutput["LogReg p-val Phylo"], alpha=0.05, method='bonferroni')
#corrected_cyt = pd.DataFrame(p_adjusted_cyt[1], columns = ["lr_p-val_tox_bonf"])
#corrected_phylo = pd.DataFrame(p_adjusted_phylo[1], columns = ["LogReg Phylo Bonferroni Corrected p-val"])
#reject_cyt = pd.DataFrame(p_adjusted_cyt[0], columns = ["LogReg Cyt Reject null hypothesis?"])
#reject_phylo = pd.DataFrame(p_adjusted_phylo[0], columns = ["LogReg Phylo Reject null hypothesis?"])
#dfOutput = pd.concat([dfOutput, corrected_cyt, corrected_phylo, reject_cyt, reject_phylo], axis = 1)
#dfOutput = dfOutput[dfOutput['Nucleotide'] != 0]

#df_filt = dfOutput[dfOutput["LogReg Cyt Reject null hypothesis?"] == True]


#Save the output files as .csv files. dfOutput has the general info about the SNPs. data is the SNP presence/absence matrix.
#os.chdir(r'C:\Users\cassp\OneDrive\Documents\Kovac Lab\Biomarkers paper\SNP hits\4_21_25')
#df_filt.to_csv(f"{gene}_logreg_4_21_25.csv")
#data.to_csv(f"{gene}_snps_4_21_25.csv")
