
cd /lustre/scratch/jmanthey/07_certhia_sim

source activate selection2

# run several training sessions (e.g., ~5) and choose the best training model for each species

python ~/diploSHIC/diploSHIC.py train training_albescens/ training_albescens/ albescensModel5
diploSHIC loss: 0.739600
diploSHIC accuracy: 0.688000

python ~/diploSHIC/diploSHIC.py train training_americana/ training_americana/ americanaModel4
diploSHIC loss: 0.647544
diploSHIC accuracy: 0.712000

mv albescensModel5.json albescensModel.json  
mv albescensModel5.weights.hdf5 albescensModel.weights.hdf5

mv americanaModel4.json americanaModel.json
mv americanaModel4.weights.hdf5 americanaModel.weights.hdf5
