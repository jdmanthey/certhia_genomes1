
cd /lustre/scratch/jmanthey/07_certhia_sim

source activate selection2

python ~/diploSHIC/diploSHIC.py train training_albescens/ training_albescens/ albescensModel

python ~/diploSHIC/diploSHIC.py train training_americana/ training_americana/ americanaModel
