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
from Bio import pairwise2



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


file = SeqIO.parse(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Sequences\prfB_n_alignment.fas", "fasta")
table = 11

seq_list = []
seq_des_list = []
seq_lengths = pd.DataFrame(0, index = range(6000), columns = ["id", "dna_len", "aa_len"])

ref = "ATGATGATGATGATGATGATAGTGATATATATAGAGGAGA"
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
    
    seq_lengths.iloc[count, 0] = seq_record.id
    seq_lengths.iloc[count, 1] = len(seq_record)
    seq_lengths.iloc[count, 2] = len(seq_record)/3
    print(seq_record.id)
    print(transORF)
    alignments = pairwise2.align.globalxx(ref, seqStr, one_alignment_only = True)
    print([a[2] for a in alignments][0])
    print(type([a[4] for a in alignments][0]))
    
    score = pairwise2.align.globalxx(ref, seqStr, score_only = True)
    perc_id = score/len(ref)

    count += 1


seq_lengths = seq_lengths.loc[(seq_lengths != 0).any(axis = 1)] #Drop any rows with all 0s

plt.hist(seq_lengths["aa_len"], bins = 50)
plt.xlabel("number of hits")
plt.ylabel("length of hit (aa)")

#SeqIO.write(seq_list, r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\Bella\lautus_hit_seqs_AA.fasta", "fasta")
