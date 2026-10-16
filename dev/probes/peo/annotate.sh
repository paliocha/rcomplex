#!/bin/bash
# Annotate trait_specific, differentiated and lineage_specific hogs.
# Submit from the login node:
#   ssh orion 'module load slurm; R=/mnt/project/FjellheimLab/martpali/
#     rcomplex-testrun; sbatch -o $R/complex_<exp>/logs/annotate-%j.out
#     $R/peo/scripts/annotate.sh <exp> <run>'
#SBATCH --partition=orion
#SBATCH --account=nmbu
#SBATCH --job-name=peo-annotate
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
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
cp "$P/scripts/annotate_hogs.R" \
  "$ROOT/complex_$1/logs/annotate_hogs-$SLURM_JOB_ID.R"
Rscript "$ROOT/complex_$1/logs/annotate_hogs-$SLURM_JOB_ID.R" "$1" "$2"
