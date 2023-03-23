#Modules
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
from collections import Counter

#Functions
def read_fasta(name):
    fasta_seqs = SeqIO.parse(open(name + '.fa'),'fasta')
    data = []
    for fasta in fasta_seqs:
        data.append([fasta.id, str(fasta.seq).strip().upper()])
            
    return data

def CntSubstr(pattern, string):
    a = [m.start() for m in re.finditer(
        '(?={0})'.format(re.escape(pattern)), string)]
    return a

def sim_score(seq1,seq2):
    return sum (seq1[i] == seq2[i] for i in range(len(seq1)))

#Data
in_fasta_1 = 'sd108'
seqs_df = pd.DataFrame(read_fasta(in_fasta_1), columns = ['name', 'sequence'])

#'NGG' on target gRNA extraction from scaffolds
grna_seq = []
count = 0
pam_pattern = ['GG'] #variable_1
for l in range(len(pam_pattern)):
    curr_pattern = pam_pattern[l]
    for i in range(np.shape(seqs_df)[0]):

        #forward strand
        curr_scaffold = seqs_df['sequence'][i]
        all_index = CntSubstr(curr_pattern, curr_scaffold) 
        for j in range(len(all_index)):
            if all_index[j] > 20:
                grna_seq.append(curr_scaffold[all_index[j]-21:all_index[j]-1])
                count = count + 1

        #reverse strand
        scaffold_seq = Seq(seqs_df['sequence'][i])
        curr_scaffold = str(scaffold_seq.reverse_complement())
        all_index = CntSubstr(curr_pattern, curr_scaffold) 
        for j in range(len(all_index)):
            if all_index[j] > 20:
                grna_seq.append(curr_scaffold[all_index[j]-21:all_index[j]-1])
                count = count + 1
                

all_d = Counter(grna_seq)
all_grna_df = pd.DataFrame.from_dict(all_d, orient='index').reset_index()
all_grna_df.columns = ['sequence', 'frequency']
unique_gg_grna_df = all_grna_df.loc[all_grna_df['frequency'] == 1].reset_index(drop=True)

no_of_max_mismatch = 6 #variable_2
max_match = 20 - no_of_max_mismatch  

grna_of_interest = []
for i in range(np.shape(unique_gg_grna_df)[0]):
    grna_flag = 1
    for j in range(np.shape(all_grna_df)[0]):
        if sim_score(unique_gg_grna_df['sequence'][i],all_grna_df['sequence'][j]) < max_match + 1:
            grna_flag = grna_flag*1
        elif sim_score(unique_gg_grna_df['sequence'][i],all_grna_df['sequence'][j]) == 20:
            grna_flag = grna_flag*1
        else:
            grna_flag = 0
            break
            
    if grna_flag == 1:
        grna_of_interest.append(unique_gg_grna_df['sequence'][i])
        
pd.DataFrame(grna_of_interest).to_csv('grna_with_' + str(no_of_max_mismatch) + 'bp_mismatch.csv',index=False)