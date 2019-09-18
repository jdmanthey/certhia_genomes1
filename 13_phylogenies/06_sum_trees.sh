# maximum clade credibility tree (simple summarization) from all gene trees using dendropy:
# gives info about which nodes have support from what proportion of gene trees
sumtrees.py --output=certhia_summed.tre --min-clade-freq=0.05 certhia_total.trees

# coalescent tree of all gene trees using ASTRAL III
# automatically calculates local branch support using quartets, described here: https://doi.org/10.1093/molbev/msw079
java -jar ~/Astral/astral.5.6.3.jar -i certhia_total.trees -o certhia_astral.tre 2> certhia_astral.log
