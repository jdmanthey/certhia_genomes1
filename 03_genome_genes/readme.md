The scripts and steps in this directory were the steps performed for annotating the certhia genome
and assessing the completeness of the genome (based on vertebrate orthologous single-copy genes).

step1 Install and verify maker can use mpi \n
step2 Run maker round 1 with an initial prediction based on proteins from three songbirds \n
step3 Run SNAP trained by the maker round 1 output \n
step4 Run augustus trained by the maker round 1 output \n
step5 Run maker again with the initial maker output and the ab initio models from SNAP and augustus \n

step20 Run BUSCO with the vertebrate orthologous single-copy genes to assess the completeness of the 
genome assembly. 
