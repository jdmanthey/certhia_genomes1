#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=fsc_inv
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

cd /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/02_no_inv

for prefix in isolation_w_migration isolation secondary_contact speciation_w_geneflow; do

	# move to model observed directory (as opposed to those for bootstraps)
	cd ${prefix}/best_ll

	# copy obs file here and rename
	cp ../../fsc_no_inv_jointMAFpop1_0.obs ${prefix}_maxL_jointMAFpop1_0.obs
	
	# run 100 replicates of fsc with these parameters
	for i in ${1..100}; do
		# run fsc
		/home/jmanthey/fsc26_linux64/fsc26 -i ${prefix}_maxL.par -n 100000 -m -q
		# get likelihood value to a file
		tail -n1 ${prefix}_maxL/${prefix}_maxL.lhoods | xargs >> ${prefix}.lhoods
		# remove iteration's output directory
		rm -r ${prefix}_maxL
	done
		
	# change directories
	cd /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/02_no_inv
done

for prefix in isolation_w_migration isolation secondary_contact speciation_w_geneflow; do

	# copy lhoods files here
	cp ${prefix}/best_ll/${prefix}.lhoods .
	
done


