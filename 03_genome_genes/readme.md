The scripts and steps in this directory were the steps performed for annotating the certhia genome
and assessing the completeness of the genome (based on vertebrate orthologous single-copy genes).

1. Install and verify maker can use mpi
2. Run maker round 1 with an initial prediction based on proteins from three songbirds.
3. Run SNAP trained by the maker round 1 output
4. Run augustus trained by the maker round 1 output
5. Run maker again with the initial maker output and the ab initio models from SNAP and augustus

20. Run BUSCO with the vertebrate orthologous single-copy genes to assess the completeness of the 
genome assembly. 
