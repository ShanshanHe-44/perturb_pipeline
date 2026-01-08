#!/usr/bin/env bash
#SBATCH --job-name=perturb_pipeline_master
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=perturb_master_%j.out
#SBATCH --error=perturb_master_%j.err
#SBATCH --partition=engreitz,owners,normal

source /home/users/shshanhe/oak_shshanhe/tools/miniforge3/bin/activate
conda activate /home/groups/engreitz/Users/tony/anaconda3/envs/kb

# call current submit.sh
echo "processing analysis ${analysis_name}"
CONFIG="$CONFIG" ./submit.sh