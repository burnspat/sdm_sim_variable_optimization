#!/bin/bash

# ----- INPUTS -----
# - SLURM -
#SBATCH --job-name=sim_model
#SBATCH --output=~/logs/sim_mod_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --mem=25G
#SBATCH --array=1-1312


# the parameters to use for each job
param_file='~/inputs/simulation_method_combos_dredgelast.csv'

# the preprocessed presence-absence + predictors table
dt='~/inputs/borneo_simspec_wpreds_cleaned_nospectra.csv'

# the folder with rasters to use for prediction. If 'NULL', will not do raster predictions
#'~/inputs/multiscale_pred_rasters/'
rast_path='NULL'

# where to save results
save_dir='~/results/'



# ----- PROCESSING -----
SECONDS=0
date_time=`date +%Y%m%d_%H%M%S`
echo "The starting date_time: " $date_time
echo
#echo "SLURM_JOBID: "$SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID: "$SLURM_ARRAY_JOB_ID
echo "SLURM ARRAY TASK ID: "$SLURM_ARRAY_TASK_ID
echo

# adjust the slurm array index to account for the header at the top of the parameter file
row_num=$(( $SLURM_ARRAY_TASK_ID + 1 ))

# grab arguments from the list of all modeling combinations for each simulation iteration
in_params=$(sed "${row_num}q;d" $param_file)
echo "The input parameters: "$in_params
echo

# runid	pa_i	seed  method	vif	mrmr	boruta	dredge	parsimony
runid=$(echo $in_params | cut -d "," -f1)
pa_i=$(echo $in_params | cut -d "," -f2)
seed=$(echo $in_params | cut -d "," -f3)
method=$(echo $in_params | cut -d "," -f4)
vif=$(echo $in_params | cut -d "," -f5)
mrmr=$(echo $in_params | cut -d "," -f6)
boruta=$(echo $in_params | cut -d "," -f7)
dredge=$(echo $in_params | cut -d "," -f8)
parsimony=$(echo $in_params | cut -d "," -f9)


# Load the R environment
module load anaconda3
conda activate /projects/above_gedi/users/pburns/envs/biodiv_mod2

# run the modeling script
echo "Starting R script..."
echo
srun Rscript /projects/above_gedi/users/zk/var_opt/code/simulation_model_fit_eval_predict.R $save_dir $dt $runid $pa_i $seed $method $vif $mrmr $boruta $dredge $parsimony $rast_path

# - ENDING -
date_time=`date +%Y%m%d_%H%M%S`
echo "The ending date_time: " $date_time
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
