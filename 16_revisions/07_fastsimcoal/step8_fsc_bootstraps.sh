#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=fsc_boot1
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-100

# move to main working directory
cd /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/02_no_inv

# set the replicate number
rep=${SLURM_ARRAY_TASK_ID}

# move to bootstrap directory
cd bs${rep}

# loop for 20 iterations of fsc
for i in {1..20}; do

	# make a directory for this replicate
	mkdir $i

	# move to that directory
	cd $i
	
	# copy the table of parameters to estimate here
	cp /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/parameters/secondary_contact.tpl .

	# copy the parameter estimate starting values here
	cp /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/parameters/secondary_contact.est .
		
	# copy the bootstrap sfs here
	cp ../secondary_contact_jointMAFpop1_0.obs .
	
	# run fastsimcoal2
	/home/jmanthey/fsc26_linux64/fsc26 -t secondary_contact.tpl -n 200000 -e secondary_contact.est -M -L 30 -q -m
	
	# move back to working directory
	cd ..

done
