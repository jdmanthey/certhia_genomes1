#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=fsc_inv
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-100

# move to main working directory
cd /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/02_no_inv

# set the replicate number
rep=${SLURM_ARRAY_TASK_ID}

for prefix in isolation_w_migration isolation secondary_contact speciation_w_geneflow; do

	# move to model observed directory (as opposed to those for bootstraps)
	cd ${prefix}/observed

	# make a directory for this replicate
	mkdir $rep

	# move to that directory
	cd $rep
	
	# copy the table of parameters to estimate here
	cp /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/parameters/${prefix}.tpl .

	# copy the parameter estimate starting values here
	cp /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/parameters/${prefix}.est .

	# copy the observed sfs here
	cp ../../../fsc_no_inv_jointMAFpop1_0.obs ${prefix}_jointMAFpop1_0.obs

	# run fastsimcoal2
	/home/jmanthey/fsc26_linux64/fsc26 -t ${prefix}.tpl -n 200000 -e ${prefix}.est -M -L 30 -q -m
	
	# move back to working directory
	cd ../../..

done

