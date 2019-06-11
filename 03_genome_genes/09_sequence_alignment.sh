# translate the sequences in the fasta file
t_coffee -other_pg seq_reformat -in 1.fasta -action +translate -output fasta_seq > 1_protein.fasta

# align the sequences
t_coffee 1_protein.fasta -mode mcoffee

# back translate to nucleotides
t_coffee -other_pg seq_reformat -in 1.fasta -in2 1_protein.aln -action +thread_dna_on_prot_aln \
-output fasta_aln > 1_aligned.fasta

# remove gaps from alignment
trimal -in 1_aligned.fasta -out 1_aligned_trimmed.fasta -nogaps -fasta

# 
