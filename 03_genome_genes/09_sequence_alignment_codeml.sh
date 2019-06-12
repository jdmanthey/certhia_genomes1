#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N c_align_paml
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -t 1-1

# working directory = /lustre/scratch/jmanthey/05_certhia_cds/cds_fasta_files

# translate the sequences in the fasta file
t_coffee -other_pg seq_reformat -in ${SGE_TASK_ID}.fasta -action +translate -output fasta_seq > ${SGE_TASK_ID}_protein.fasta

# align the sequences
t_coffee ${SGE_TASK_ID}_protein.fasta -mode mcoffee

# back translate to nucleotides
t_coffee -other_pg seq_reformat -in ${SGE_TASK_ID}.fasta -in2 ${SGE_TASK_ID}_protein.aln -action +thread_dna_on_prot_aln \
-output fasta_aln > ${SGE_TASK_ID}_aligned.fasta

# remove gaps from alignment
trimal -in ${SGE_TASK_ID}_aligned.fasta -out ${SGE_TASK_ID}_aligned_trimmed.fasta -nogaps -fasta

# create a new directory for this alignment and move to it
mkdir ${SGE_TASK_ID}_work
cd ${SGE_TASK_ID}_work/

# copy a control file to this directory
cp ../../codeml.ctl .

# copy the alignment to this directory
cp ../${SGE_TASK_ID}_aligned_trimmed.fasta .

# change the control file input and output names 
sed -i "s/blank.fasta/${SGE_TASK_ID}_aligned_trimmed.fasta/g" codeml.ctl
sed -i "s/blank_output.txt/${SGE_TASK_ID}_output.txt/g" codeml.ctl

# run 1st iteration of codeml
codeml codeml.ctl

# copy certhia control file to this directory
cp ../../codeml_2.ctl .

# change the control file input and output names for certhia
sed -i "s/blank.fasta/${SGE_TASK_ID}_aligned_trimmed.fasta/g" codeml_2.ctl
sed -i "s/blank_output.txt/${SGE_TASK_ID}_output_certhia.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree.tre/unrooted_tree_certhia.tre/g' codeml_2.ctl

# run certhia iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for ficedula
sed -i "s/${SGE_TASK_ID}_output_certhia.txt/${SGE_TASK_ID}_output_ficedula.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_certhia.tre/unrooted_tree_ficedula.tre/g' codeml_2.ctl

# run ficedula iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for parus
sed -i "s/${SGE_TASK_ID}_output_ficedula.txt/${SGE_TASK_ID}_output_parus.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_ficedula.tre/unrooted_tree_parus.tre/g' codeml_2.ctl

# run parus iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for taeniopygia
sed -i "s/${SGE_TASK_ID}_output_parus.txt/${SGE_TASK_ID}_output_taeniopygia.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_parus.tre/unrooted_tree_taeniopygia.tre/g' codeml_2.ctl

# run parus iteration of codeml
codeml codeml_2.ctl

# add all likelihoods and omegas to a single file
grep '^lnL' ${SGE_TASK_ID}_output.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^omega' ${SGE_TASK_ID}_output.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_certhia.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_certhia.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_ficedula.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_ficedula.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_parus.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_parus.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_taeniopygia.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_taeniopygia.txt >> /lustre/scratch/jmanthey/05_certhia_cds/output/${SGE_TASK_ID}_total_output.txt

# move the alignment fasta to the output directory 
mv ${SGE_TASK_ID}_aligned_trimmed.fasta ../../output
