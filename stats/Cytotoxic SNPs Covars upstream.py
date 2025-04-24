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
os.chdir(r'C:\Users\cassp\OneDrive\Documents\Kovac Lab\Biomarkers paper\upstream_4_21_25')

file = 'cytkUp_filt.fasta'
mSheet = pd.read_excel(r"C:\Users\cassp\OneDrive\Documents\Kovac Lab\Biomarkers paper\Mastersheet_no_clones.xlsx") #import file with cytotoxicity data (same order as fasta)

mSheet['Adjusted panC Group (predicted species)'] = mSheet['Adjusted panC Group (predicted species)'].replace({'Group_clarus' : 0, 'Group_I(pseudomycoides)' : 1, 'Group_II(mosaicus/luti)': 2, 'Group_III(mosaicus)': 3, 'Group_IV(cereus_sensu_stricto)': 4, 'Group_V(toyonensis)': 5, 'Group_VI(mycoides/paramycoides)': 6, 'Group_VII(cytotoxicus)' : 7, 'Group_VIII(mycoides)' : 8})

metadata = mSheet.loc[:, ["Isolate","Average Cell Viability ( >0.7 is cytotoxic)", "Adjusted panC Group (predicted species)"]]

#Create empty dataframe for output.
dfOutput = pd.DataFrame(0, index = range(100000), columns = ["Gene", "Position", "Nucleotide","name", "Clades Represented", "True Negatives", "False Positives", "False Negatives", "True Positives", "# with", "# without", "Avg Cytotoxicity With", "Avg Cytotoxicity Without", "Stdev Cytotoxicity With", "Stdev Cytotoxicity Without", "Accuracy", "Precision", "LogReg p-val Cyt", "LogReg p-val Phylo"])

fileName = str(file)
fileName = fileName.replace(".fasta", "")
fileName = fileName.replace("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\upstream_4_21_25\\", "")
gene = fileName.replace("Up", "_promoter")
print(gene)

alignmentSeq = AlignIO.read(file, "fasta") #Import alignment file of sequences in fasta format.
alignment = np.array([list(rec) for rec in alignmentSeq])

IDs = []
for record in alignmentSeq:
    #Change names to match mastersheet. If you downloaded the sequences of hits from BTyper3, this should work. If you procured sequences by other means, you may need to change this part so you can get matching names with the mastersheet.
    
    recordID = record.id
    print(recordID)
    PS_ID = recordID.split("___", 1)[1]
    PS_ID = PS_ID.split("_", 1)[0]
    print(PS_ID)
    row = mSheet[mSheet["Isolate"].str.contains(PS_ID)]
    if len(row) != 0:
        IDs.append(PS_ID)
        
    


dfIDs = pd.DataFrame(IDs, columns = ["Isolate"]) 

data = dfIDs.merge(metadata)
#data = data.sort_values('Isolate')
cytotoxicity = data['Average Cell Viability ( >0.7 is cytotoxic)'].values
cytotoxicity = cytotoxicity.reshape(-1,1)

#Define the function that applies the logistic regression model to each site in the alignment.
def modelTox(nTP, cutoff1, cutoff2, p_cutoff, a_cutoff):
    position = 0
    while position < len(alignment[0]):
        #print(position)
        nuc_hyph = alignment[: , position]
    
        #Remove hyphens so they don't have an impact as "negatives".
        cytotoxicityShort = cytotoxicity[(nuc_hyph!="-")]
        dataWith = data[(nuc_hyph!="-")]
        testLR = data[(nuc_hyph!="-")]    
        
        nuc = nuc_hyph[(nuc_hyph!="-")]
        #dataWith = dataWith[nuc == nTP]
        clade = []
        
        #If there is sufficient variation between sequences to correlate cytotox and nucleotide (mostly because f scores were 0 in most cases), pursue the site as possibly containing a SNP. The cutoffs for # of sequences with the given nucleotide in the site are between 20 and 80%. This prevents the perfectly conserved nucleotides from being investigated, as they wouldn't be informative. Also if a nucleotide is very rarely found, it won't be investigated. We want SNPs that aren't extremely rare. Those would likely be unhelpful in application.
        if list(nuc).count(nTP) <= cutoff1 * len(nuc) and list(nuc).count(nTP) >= cutoff2 * len(nuc): 
            y = np.where((nuc == nTP), 1, 0) #Where a given nucleotide exists in the position, the array is assigned a 1. Where there is not, the array is assigned a 0.
            withNTP = cytotoxicityShort[nuc == nTP]
            withoutNTP = cytotoxicityShort[nuc != nTP]
                      
            #Keep track of the clades that have the SNP.
            for isolate, row in dataWith.iterrows():
                if row["Adjusted panC Group (predicted species)"] not in clade:
                    clade.append(row["Adjusted panC Group (predicted species)"])
    
            testY = pd.DataFrame(y, dtype = "int")
            testLR = testLR[['Average Cell Viability ( >0.7 is cytotoxic)', 'Adjusted panC Group (predicted species)']]
            testLR = testLR.reset_index(drop=True)
            #print(position)
            print(position)
            log_reg = sm.Logit(testY, testLR).fit() #Fit the logistic regression model.
            pred = list(map(round, log_reg.predict(testLR)))
            #print(log_reg.summary())
            #print(log_reg.pvalues)
            
            #Confusion matrix of actual gene presence vs predicted gene presence based on the LogReg sigmoid line. ***True negatives are top-left, true-positives are bottom-right.***
            confMatrix = confusion_matrix(testY, pred).flatten()
            acc = accuracy_score(y, pred)
            prec = precision_score(y, pred)
            #print('Accuracy score:', acc)
            #print('Precision score:', prec)
    
            #In some cases, all of the isolates were being classified as positives or negatives. This is very uninformative, so I filtered for SNPs that didn't have this problem. I also added precision and accuracy cutoffs. Currently they're set at 70%.
            if (0 not in confMatrix) and (prec > p_cutoff) and (acc > a_cutoff): 
                #Record general data.
                dfOutput.loc[position, "Gene"] = gene
                dfOutput.loc[position, "Position"] = position + 1
                dfOutput.loc[position, "Nucleotide"] = nTP
                dfOutput.loc[position, "name"] = gene + '_SNP_' + str(position+1)
                dfOutput.loc[position, "Clades Represented"] = str(clade)
                dfOutput.loc[position, "Avg Cytotoxicity With"] = np.mean(withNTP)
                dfOutput.loc[position, "Avg Cytotoxicity Without"] = np.mean(withoutNTP)
                dfOutput.loc[position, "Stdev Cytotoxicity With"] = np.std(withNTP)
                dfOutput.loc[position, "Stdev Cytotoxicity Without"] = np.std(withoutNTP)
                dfOutput.loc[position, "# with"] = len(withNTP)
                dfOutput.loc[position, "# without"] = len(withoutNTP)
                #Record logistic regression specific data.
                dfOutput.iloc[position, 5:9] = confMatrix
                dfOutput.loc[position, "Accuracy"] = acc
                dfOutput.loc[position, "Precision"] = prec
                dfOutput.loc[position, "LogReg p-val Cyt"] = log_reg.pvalues[0]
                dfOutput.loc[position, "LogReg p-val Phylo"] = log_reg.pvalues[1]
                
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

cutoff1 = 0.8
cutoff2 = 0.2
p_cutoff = 0.7
a_cutoff = 0.7

#Run the model on each nucleotide. depending on your alignment, you may need to capitalize the letters. You could also run this program on amino acid alignments. You just need to change the letter to the single-letter amino acid code. Ex. "F" for phenylalanine.
modelTox("a", cutoff1, cutoff2, p_cutoff, a_cutoff)
modelTox("t", cutoff1, cutoff2, p_cutoff, a_cutoff)
modelTox("g", cutoff1, cutoff2, p_cutoff, a_cutoff)
modelTox("c", cutoff1, cutoff2, p_cutoff, a_cutoff)

#Bonferroni correction of p-values.
p_adjusted_cyt = multipletests(dfOutput["LogReg p-val Cyt"], alpha=0.05, method='bonferroni')
p_adjusted_phylo = multipletests(dfOutput["LogReg p-val Phylo"], alpha=0.05, method='bonferroni')
corrected_cyt = pd.DataFrame(p_adjusted_cyt[1], columns = ["LogReg Cyt Bonferroni Corrected p-val"])
corrected_phylo = pd.DataFrame(p_adjusted_phylo[1], columns = ["LogReg Phylo Bonferroni Corrected p-val"])
reject_cyt = pd.DataFrame(p_adjusted_cyt[0], columns = ["LogReg Cyt Reject null hypothesis?"])
reject_phylo = pd.DataFrame(p_adjusted_phylo[0], columns = ["LogReg Phylo Reject null hypothesis?"])
dfOutput = pd.concat([dfOutput, corrected_cyt, corrected_phylo, reject_cyt, reject_phylo], axis = 1)
dfOutput = dfOutput[dfOutput['Nucleotide'] != 0]

df_filt = dfOutput[dfOutput["LogReg Cyt Reject null hypothesis?"] == True]
data_filt = pd.concat([data.iloc[:,0:3], data[dfOutput["name"]]], axis = 1)

#Save the output files as .csv files. dfOutput has the general info about the SNPs. data is the SNP presence/absence matrix.
os.chdir(r'C:\Users\cassp\OneDrive\Documents\Kovac Lab\Biomarkers paper\upstream_4_21_25')
df_filt.to_csv(f"{gene}_logreg_4_21_25.csv")
data.to_csv(f"{gene}_snps_4_21_25.csv")
