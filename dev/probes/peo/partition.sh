#!/bin/bash
# WP24 partition probe: stage a partition run, then run.R on it, on the
# 64-bit Armadillo library (the partition build is dense).
# Submit from the login node:
#   ssh orion 'module load slurm; R=/mnt/project/FjellheimLab/martpali/
#     rcomplex-testrun; sbatch -o $R/complex_<exp>/logs/<run>-%j.out
#     $R/peo/scripts/partition.sh <exp> <run>'
# run: organ_partition (LeafvRoot) or kmeans_partition (Shin2021, C3vC4).
#SBATCH --partition=orion
#SBATCH --account=nmbu
#SBATCH --job-name=peo-partition
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
set -euo pipefail
ROOT=/mnt/project/FjellheimLab/martpali/rcomplex-testrun
P=$ROOT/peo
L=$ROOT/complex_$1/logs
set +u
module load R/4.4.2
set -u
export R_LIBS="$P/lib-arma64:$P/lib-0.4.2:$P/lib-deps"
R_LIBS="$R_LIBS:$ROOT/validation-2026-09-17/lib-deps"
R_LIBS="$R_LIBS:/mnt/users/martpali/R/library-4.4"
export PEO_ROOT=$ROOT
# run copies, so a later push cannot change the scripts under the job;
# run_wp24.R is run.R with partition support
cp "$P/scripts/stage_partition.R" "$L/stage_partition-$SLURM_JOB_ID.R"
cp "$P/scripts/run_wp24.R" "$L/run-$SLURM_JOB_ID.R"
Rscript "$L/stage_partition-$SLURM_JOB_ID.R" "$1" "$2"
Rscript "$L/run-$SLURM_JOB_ID.R" "$1" "$2"
