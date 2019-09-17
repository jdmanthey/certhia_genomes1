# run with up to three migration edges:
src/treemix -i input/certhia_contact_10000.treemix.gz -root Outgroup -o output/certhia_contact_10000.treemix
src/treemix -i input/certhia_contact_10000.treemix.gz -m 1 -g output/certhia_contact_10000.treemix.vertices.gz output/certhia_contact_10000.treemix.edges.gz -o output/certhia_contact_10000_m1.treemix
src/treemix -i input/certhia_contact_10000.treemix.gz -m 1 -g output/certhia_contact_10000_m1.treemix.vertices.gz output/certhia_contact_10000_m1.treemix.edges.gz -o output/certhia_contact_10000_m2.treemix
src/treemix -i input/certhia_contact_10000.treemix.gz -m 1 -g output/certhia_contact_10000_m2.treemix.vertices.gz output/certhia_contact_10000_m2.treemix.edges.gz -o output/certhia_contact_10000_m3.treemix


#bootstraps over 200s snps with migration edges
for i in {1..100}; do
    src/treemix -i input/certhia_contact_10000.treemix.gz -m 1 -g output/certhia_contact_10000_m2.treemix.vertices.gz output/certhia_contact_10000_m2.treemix.edges.gz -bootstrap -k 200 -o output/$i.treemix
done;

# unzip the tree files
for i in { ls *treeout.gz }; do
    gzip -d $i
done;

# in R:
#summarize bootstraps
x <- list.files(pattern="*treeout")
for(a in 1:length(x)) {
	if (a==1) {
		output <- scan(x[a], what="character")[1]
	} else {
		output <- c(output, scan(x[a], what="character")[1])
	}
}
write(output, file="certhia_10000_bootstraps.trees", ncolumns=1)

# in bash 
# summarize bootstraps
sumtrees.py --output=certhia_treemix_summed.tre --min-clade-freq=0.05 certhia_10000_bootstraps.trees
