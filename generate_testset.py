#%%
from Bio import SeqIO
import random

from Bio.Seq import Seq

#%%
def number_of_sequences(fasta_file):
    """
    Counts the number of sequences in a fasta file.
    """
    return len(list(SeqIO.parse(fasta_file, "fasta")))
#%%
#%%
def subsample_fasta(infile,outfile,target_count):
    """
    Subsamples a fasta file.
    """
    in_count = number_of_sequences(infile)
    reduction_factor = in_count//target_count
    with open(outfile, "w") as outfile:
        for n,record in enumerate(SeqIO.parse(infile, "fasta")):
            if n % reduction_factor == 0:
                SeqIO.write(record, outfile, "fasta-2line")
#%%
subsample_fasta("input/clades.fasta","intermediary/clades_subsampled.fasta",15)
#%%
def subsequence_fasta(infile,outfile,start,end):
    """
    Extracts a subsequence from a fasta file.
    """
    with open(outfile, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            record.id = record.id + " subsequence " + str(start) + "_" + str(end)
            record.name = ""
            record.description = ""
            SeqIO.write(record[start:end], outfile, "fasta-2line")
# %%
subsequence_fasta("intermediary/clades_subsampled.fasta","intermediary/S.fasta",21500,26000)
subsequence_fasta("intermediary/clades_subsampled.fasta","intermediary/ORF1a.fasta",200,14000)
#%%
def sprinkle_N_stretches(infile,outfile,p_to_N, p_from_N):
    """
    Turn stretches into Ns using Markov process
    """
    with open(outfile, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            state = 0
            sequence_in = str(record.seq)
            sequence_out = ""
            for base in sequence_in:
                if state == 0:
                    if random.random() < p_to_N:
                        state = 1
                else:
                    if random.random() < p_from_N:
                        state = 0
                if state == 0:
                    sequence_out += base
                else:
                    sequence_out += "N"
            record.seq = Seq(sequence_out)
            record.id = record.id + " N stretches"
            record.name = ""
            record.description = ""
            SeqIO.write(record, outfile, "fasta-2line")

# %%
sprinkle_N_stretches("intermediary/S.fasta","intermediary/S_N.fasta",0.002,0.02)
# %%
def sprinkle_early_indels(infile,outfile,p_del,last_indel_pos):
    with open(outfile, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            sequence_in = str(record.seq)
            sequence_out = ""
            for pos,base in enumerate(sequence_in):
                if pos < last_indel_pos:
                    if random.random() > p_del:
                        sequence_out += base
                else:
                    sequence_out += base
            record.seq = Seq(sequence_out)
            record.id = record.id + " early indels"
            record.name = ""
            record.description = ""
            SeqIO.write(record, outfile, "fasta-2line")
#%%
sprinkle_early_indels("intermediary/ORF1a.fasta","intermediary/ORF1_early_indels.fasta",0.002,500)
# %%
import glob
with open("output/all_sequences.fasta", "w") as outfile:
    for file in glob.glob("intermediary/*.fasta"):
        with open(file,"r") as infile:
            outfile.write(infile.read())
# %%
