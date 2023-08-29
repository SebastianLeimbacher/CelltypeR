#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --account=rrg-tdurcan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4  # Number of CPU cores per task
#SBATCH --mem=16G  # Memory per node

# Specify the log file path
LOG_FILE=/lustre03/project/6070393/rhalena/Test_celltypeR/logfile_trainRFM.log





module load r/4.2.1
# Set environment variables for R
R_VERSION=4.2.1
R_LIB_PATH=/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2

# Run your R script using Rscript and add the log file to save the print statments
Rscript TrainRFM.R > $LOG_FILE 2>&1

