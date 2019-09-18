source activate smcpp_mod


cd /lustre/scratch/jmanthey/01_certhia_genomics/03_vcf/smc_vcf/

cd split_1

smc++ plot split_1.pdf model.final.json -c -g 2

cd ../split_2

smc++ plot split_2.pdf model.final.json -c -g 2

cd ../split_3

smc++ plot split_3.pdf model.final.json -c -g 2

cd ../split_4

smc++ plot split_4.pdf model.final.json -c -g 2

cd ../split_5

smc++ plot split_5.pdf model.final.json -c -g 2

cd ../split_6

smc++ plot split_6.pdf model.final.json -c -g 2

cd ../split_7

smc++ plot split_7.pdf model.final.json -c -g 2

cd ../split_8

smc++ plot split_8.pdf model.final.json -c -g 2

cd ../split_9

smc++ plot split_9.pdf model.final.json -c -g 2

cd ../split_10

smc++ plot split_10.pdf model.final.json -c -g 2
