# working directory = /lustre/scratch/jmanthey/05_certhia_cds/cds_fasta_files

# translate the sequences in the fasta file
t_coffee -other_pg seq_reformat -in 1.fasta -action +translate -output fasta_seq > 1_protein.fasta

# align the sequences
t_coffee 1_protein.fasta -mode mcoffee

# replace any stop codons (including the final stop codon) with a gap
sed -i 's/X/-/g' 1_protein.aln

# back translate to nucleotides
t_coffee -other_pg seq_reformat -in 1.fasta -in2 1_protein.aln -action +thread_dna_on_prot_aln \
-output fasta_aln > 1_aligned.fasta

# remove gaps from alignment
trimal -in 1_aligned.fasta -out 1_aligned_trimmed.fasta -nogaps -fasta

# create a new directory for this alignment and move to it
mkdir 1
cd 1/

# copy a control file to this directory
cp ../../codeml.ctl .

# copy the alignment to this directory
cp ../1_aligned_trimmed.fasta .

# change the control file input and output names 
sed -i 's/blank.fasta/1_aligned_trimmed.fasta/g' codeml.ctl
sed -i 's/blank_output.txt/1_output.txt/g' codeml.ctl

# run 1st iteration of codeml
codeml codeml.ctl




# copy certhia control file to this directory
cp ../../codeml_2.ctl .

# change the control file input and output names 
sed -i 's/blank.fasta/1_aligned_trimmed.fasta/g' codeml_2.ctl
sed -i 's/blank_output.txt/1_output_certhia.txt/g' codeml_2.ctl

# run certhia iteration of codeml
codeml codeml_2.ctl

