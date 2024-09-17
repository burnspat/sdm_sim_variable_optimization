# sdm_sim_variable_optimization
code repo for Ecological Informatics publication "Simulating multi-scale optimization and variable selection in species distribution modeling", by Cushman et al. 2024

We first simulated species occurrence across a portion of Borneo using four known predictor variables and parameters found within the crosswalk folder. For each simulated occurrence point we then extracted a suite of multi-scale predictor variables (simulation_pa_predictors_extract_gee.js) and used various variable selection methods and model algorithms to model the simulated distribution (simulation_model_sbatch.sh). 

We used Northern Arizona University's High Performance Computing System and Slurm Workload Manager to distribute jobs, where each job corresponds to a different simulated species modeling scenario (i.e. different variable selection method or simulation bootstrap). This set of scripts will not run "out of the box" unless Slurm is available on your computer/HPC and you have installed the necessary R packages. You will also need to manually change input and output file paths. The required R packages are listed in simulation_funcs.R. Note that we used a custom conda environment to manage packages. You can clone that environment using this yml filebiodiv_mod.yml, but be aware that it contains extra packages that aren't necessary for this analysis.

The R script simulation_funcs.R contains various custom modeling functions and is sourced by the R script simulation_model_fit_eval_predict.R. In theory you don't need to edit these scripts.

To run the scripts, first edit the arguments in INPUTS section of simulation_model_sbatch.sh. Then run the shell script from the linux command line, like so
$ sh simulation_model_sbatch.sh
