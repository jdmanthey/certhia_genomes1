# adds the smc meta data header for each file from the original full data files


for i in $( ls -d bootstrap_1/* );
do j=$(basename $i);
head -n1 $j | cat - bootstrap_1/$j > temp && mv temp bootstrap_1/$j;
head -n1 $j | cat - bootstrap_10/$j > temp && mv temp bootstrap_10/$j;
head -n1 $j | cat - bootstrap_2/$j > temp && mv temp bootstrap_2/$j;
head -n1 $j | cat - bootstrap_3/$j > temp && mv temp bootstrap_3/$j;
head -n1 $j | cat - bootstrap_4/$j > temp && mv temp bootstrap_4/$j;
head -n1 $j | cat - bootstrap_5/$j > temp && mv temp bootstrap_5/$j;
head -n1 $j | cat - bootstrap_6/$j > temp && mv temp bootstrap_6/$j;
head -n1 $j | cat - bootstrap_7/$j > temp && mv temp bootstrap_7/$j;
head -n1 $j | cat - bootstrap_8/$j > temp && mv temp bootstrap_8/$j;
head -n1 $j | cat - bootstrap_9/$j > temp && mv temp bootstrap_9/$j;
done

