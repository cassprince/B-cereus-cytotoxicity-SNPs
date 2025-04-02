# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 14:06:02 2023

@author: cassp
"""

###The logreg results suck, stick with t-test.
###

import pandas as pd
import numpy as np
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, precision_score
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
import statsmodels.api as sm

mSheet = pd.read_excel(r"C:\Users\cassp\OneDrive\Documents\Kovac Lab\Biomarkers paper\Mastersheet_no_clones.xlsx")
cytotoxicity = mSheet.loc[:, "Average Cell Viability ( >0.7 is cytotoxic)"]
testY = pd.DataFrame(cytotoxicity)

yvals_prep = mSheet.loc[:, "vir|nhe":"sphingomyelinase_Sph(gene)"]

dfOutput = pd.DataFrame(0, index = range(14), columns = ["Gene", "True Negatives", "False Positives", "False Negatives", "True Positives", "# with", "# without", "Avg Cytotoxicity With", "Avg Cytotoxicity Without", "Stdev Cytotoxicity With", "Stdev Cytotoxicity Without", "Accuracy", "Precision", "LogReg p-val"])

yvals = pd.DataFrame({'vir_nheA' : []})

yvals["vir_nheA"] = np.where(yvals_prep["vir|nhe"].str.contains("nheA"), 1, 0)
yvals["vir_nheB"] = np.where(yvals_prep["vir|nhe"].str.contains("nheB"), 1, 0)
yvals["vir_nheC"] = np.where(yvals_prep["vir|nhe"].str.contains("nheC"), 1, 0)
yvals["vir_nheABC"] = np.where(yvals_prep["vir|nhe"].str.contains("3/3"), 1, 0)

yvals["vir_hblA"] = np.where(yvals_prep["vir|hbl"].str.contains("hblA"), 1, 0)
yvals["vir_hblB"] = np.where(yvals_prep["vir|hbl"].str.contains("hblB"), 1, 0)
yvals["vir_hblC"] = np.where(yvals_prep["vir|hbl"].str.contains("hblC"), 1, 0)
yvals["vir_hblD"] = np.where(yvals_prep["vir|hbl"].str.contains("hblD"), 1, 0)
yvals["vir_hblACD"] = np.where((yvals_prep["vir|hbl"].str.contains("hblA") & yvals_prep["vir|hbl"].str.contains("hblC") & yvals_prep["vir|hbl"].str.contains("hblD")), 1, 0)
yvals["vir_hblABCD"] = np.where((yvals_prep["vir|hbl"].str.contains("hblA") & yvals_prep["vir|hbl"].str.contains("hblC") & yvals_prep["vir|hbl"].str.contains("hblD") & yvals_prep["vir|hbl"].str.contains("hblB")), 1, 0)

yvals["vir_cytK1"] = np.where(yvals_prep["vir|cytK"].str.contains("cytK-1"), 1, 0)
yvals["vir_cytK2"] = np.where(yvals_prep["vir|cytK"].str.contains("cytK-2"), 1, 0)
yvals["vir_sph"] = np.where(yvals_prep["sphingomyelinase_Sph(gene)"].str.contains("1/1"), 1, 0)


position = 0
for gene in yvals:
    col = yvals[gene]
    withGene = cytotoxicity[col == 1]
    withoutGene = cytotoxicity[col == 0]
    variance = withGene.var()/withoutGene.var()
    
    if 4 > variance > 0.25:
        test = stats.ttest_ind(withGene, withoutGene, equal_var = True)
    else:
        test = stats.ttest_ind(withGene, withoutGene, equal_var = False)
    
    log_reg = sm.Logit(col, testY).fit() #Fit the logistic regression model.
    print(log_reg.summary())
    pred = list(map(round, log_reg.predict(col)))
    confMatrix = confusion_matrix(col, pred).flatten()
    
    #print(confMatrix)
    
    #print(f"Model score: {model.score(cytotoxicity, col)}") #accuracy of fit
    #print(confusion_matrix(col, model.predict(cytotoxicity)))
    #print('report:', (classification_report(col, model.predict(cytotoxicity))), sep='\n')
    dfOutput.loc[position, "Gene"] = gene
    dfOutput.loc[position, "Avg Cytotoxicity With"] = withGene.mean()
    dfOutput.loc[position, "Avg Cytotoxicity Without"] = withoutGene.mean()
    dfOutput.loc[position, "Stdev Cytotoxicity With"] = withGene.std()
    dfOutput.loc[position, "Stdev Cytotoxicity Without"] = withoutGene.std()
    #Record logistic regression specific data.
    dfOutput.iloc[position, 1:5] = confMatrix
    dfOutput.loc[position, "Accuracy"] = accuracy_score(col, pred)
    dfOutput.loc[position, "Precision"] = precision_score(col, pred)
    dfOutput.loc[position, "LogReg p-val Cyt"] = log_reg.pvalues[0]
    #Record two-sample t-test specific data.
    dfOutput.loc[position, "# with"] = len(withGene)
    dfOutput.loc[position, "# without"] = len(withoutGene)
    
    position += 1


p_adjusted_LR = multipletests(dfOutput["LogReg p-val"], alpha=0.05, method='bonferroni')
corrected_LR = pd.DataFrame(p_adjusted_LR[1], columns = ["LogReg Bonferroni Corrected p-val"])
reject_LR = pd.DataFrame(p_adjusted_LR[0], columns = ["LogReg Reject null hypothesis?"])
dfOutput = pd.concat([dfOutput, corrected_LR, reject_LR], axis = 1)

dfOutput.to_csv(r"C:\Users\cassp\OneDrive\Documents\Kovac Lab\Biomarkers paper\gene_presence_stats_statsmodels.csv")
