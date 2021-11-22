#%%
from Bio import SeqIO

#%%
def number_of_sequences(fasta_file):
    """
    Counts the number of sequences in a fasta file.
    """
    return len(list(SeqIO.parse(fasta_file, "fasta")))
#%%
print(number_of_sequences("clades.fasta"))
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
subsample_fasta("clades.fasta","clades_subsampled.fasta",50)
#%%
# %%
