#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=stat_windows
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-250

# Set the number of runs that each SLURM task should do
PER_TASK=40

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

module load intel R

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
	echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run

	input_array=$( head -n${run} window_list.txt | tail -n1 )

	Rscript calculate_windows.r /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/01_basic_filter/windows/${input_array} window_popmap.txt

done
