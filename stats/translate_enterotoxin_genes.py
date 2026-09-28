# -*- coding: utf-8 -*-
"""
Created on Thu May 25 18:41:08 2023

@author: cassp
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import pandas as pd
import matplotlib.pyplot as plt
import os

def find_orf(seq, table):
    lengths = []
    for frame in range(3):
        trans = seq[frame:].translate(table)
        stopCount = trans.count("*")
        indeces = [m.start() for m in re.finditer("\*", str(trans))]
        lengths.append(len(indeces))
        """
        print(" ")
        print("frame: +", frame)
        print(trans)
        print("Number of stop codons", stopCount)
        """
    return(lengths.index(min(lengths)))


os.chdir(r"C:\Users\cassp\OneDrive\Documents\GitHub\B-cereus-cytotoxicity-SNPs\blast_results\filtered_megablast_qc50_090126\fasta_files\gene_seqs")
table = 11

for x in os.listdir(r"C:\Users\cassp\OneDrive\Documents\GitHub\B-cereus-cytotoxicity-SNPs\blast_results\filtered_megablast_qc50_090126\fasta_files\gene_seqs"):
    
    name = re.sub(".fasta", "_prot.fasta", x)
    print(name)
    file = SeqIO.parse(x, "fasta")
    
    seq_list = []
    seq_des_list = []
    
    count = 0
    for seq_record in file:
        seqStr = str(seq_record.seq)
        if "-" in seqStr:
            seqStr = re.sub("-", "", seqStr)
            seq = Seq(seqStr)
        else:
            seq = seq_record.seq
            
        ORF = find_orf(seq, table)
        transORF = seq[ORF:].translate(table)
        rec = SeqRecord(transORF, id = seq_record.id)
        seq_list.append(rec)
        
        #print(seq_record.id)
        #print(transORF)
    
        count += 1

    SeqIO.write(seq_list, "C:\\Users\\cassp\\OneDrive\\Documents\\GitHub\\B-cereus-cytotoxicity-SNPs\\blast_results\\filtered_megablast_qc50_090126\\fasta_files\\translated_seqs\\" + name, "fasta")
