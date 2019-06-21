
cd /lustre/scratch/jmanthey/07b_certhia_diplo

source activate selection2

python ~/diploSHIC/diploSHIC.py train training_albescens/ training_albescens/ albescensModel
diploSHIC loss: 0.762922
diploSHIC accuracy: 0.656000

python ~/diploSHIC/diploSHIC.py train training_americana/ training_americana/ americanaModel
diploSHIC loss: 0.682005
diploSHIC accuracy: 0.690000
