## Process to look for selective sweeps between the two certhia lineages

1. Run the commands in the 01* scripts to run discoal simulations. I chose to separate simulations because of memory errors sometimes causing problems. The reasoning behind parameter priors and estimates of demography, etc. are in the discoal_notes.txt file. 
2. Concatenate all the output with the commands in 02* and then use R to get 2000 of the simulations for each scenario (03 script). Make sure to change the header so it is clear there were 2000 simulations (03b script). Zip the simulations (04*) for use in diploS/HIC.
3. Use the demographic analysis reference genome masks to create a mask for analysis here (05*). 
4. Calculate feature vectors for the simulation data (06* scripts) and organize all of the data for training (07*). 
5. Training of the prediction algorithm (08*).
6. Prepare SNP genotype vcfs for use in diploS/HIC (09* scripts).
7. Use R script (10*) to create cluster submission jobs for turning the vcf files into feature vectors for two different sized windows (20 kbp and 50 kbp). 
8. Predict patterns in empirical data from the simulation trained models (11*).
9. Concatenate the output predictions (12*).
10. Plot and explore output in R (13*).
