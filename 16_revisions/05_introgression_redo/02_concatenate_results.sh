# combine the output into a single file
# first add a header
grep 'start' Ca_0001_Tg_2:100000001-100100000__stats.txt > ../window_introgression.txt
# add the  stats to the file
for i in $( ls *txt ); do grep 'outgroup' $i >> ../window_introgression.txt; done
