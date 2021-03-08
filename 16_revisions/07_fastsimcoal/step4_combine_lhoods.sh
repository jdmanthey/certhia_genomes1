
for prefix in isolation_w_migration isolation secondary_contact speciation_w_geneflow; do

	# make the likelihoods file
	grep 'NPOP1' $prefix/observed/1/$prefix/$prefix.bestlhoods > $prefix.likelihoods

	for i in {1..100}; do
	grep -v 'NPOP1' $prefix/observed/$i/$prefix/$prefix.bestlhoods >> $prefix.likelihoods
	done
done

