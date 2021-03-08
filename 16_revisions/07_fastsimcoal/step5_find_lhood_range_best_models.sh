cd /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/02_no_inv

cd isolation_w_migration
cp -r observed/6/isolation_w_migration best_ll
cd ../isolation
cp -r observed/17/isolation best_ll
cd ../secondary_contact
cp -r observed/15/secondary_contact best_ll
cd ../speciation_w_geneflow
cp -r observed/9/speciation_w_geneflow best_ll
cd ..


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



# copy to working directory
for prefix in diff_geneflow isolation_w_migration isolation secondary_contact speciation_w_geneflow; do
	cp ${prefix}/best_ll/${prefix}.lhoods .
done





