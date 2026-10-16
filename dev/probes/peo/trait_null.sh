#!/bin/bash
# Trait pass on the null edges of a finished run (trait_null.R).
# Submit from the login node:
#   ssh orion 'module load slurm; R=/mnt/project/FjellheimLab/martpali/
#     rcomplex-testrun; sbatch -o $R/complex_<exp>/logs/trait_null-%j.out
#     $R/peo/scripts/trait_null.sh <exp> <run>'
#SBATCH --partition=orion
#SBATCH --account=nmbu
#SBATCH --job-name=peo-trait-null
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
set -euo pipefail
ROOT=/mnt/project/FjellheimLab/martpali/rcomplex-testrun
P=$ROOT/peo
set +u
module load R/4.4.2
set -u
export R_LIBS="$P/lib-0.4.2:$P/lib-deps:$ROOT/validation-2026-09-17/lib-deps"
R_LIBS="$R_LIBS:/mnt/users/martpali/R/library-4.4"
export PEO_ROOT=$ROOT
# run a copy, so a later push cannot change the script under the job
cp "$P/scripts/trait_null.R" \
  "$ROOT/complex_$1/logs/trait_null-$SLURM_JOB_ID.R"
Rscript "$ROOT/complex_$1/logs/trait_null-$SLURM_JOB_ID.R" "$1" "$2"
