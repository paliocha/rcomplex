#!/bin/bash
#SBATCH --partition=orion
#SBATCH --account=nmbu
#SBATCH --job-name=FastOMA-Shin
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=FAIL

# FastOMA v0.5.1 on the six Ensembl Plants release-38 proteomes of
# Shin et al. 2021. Copy of peo-fastoma/run_fastoma.sh with every
# writable path (launch dir, work, output, logs, .nextflow, NXF_HOME)
# under complex_Shin2021/fastoma. The image cache and the OMAmer db are
# read from the PEO run. Head job: Nextflow submits its own jobs.
# Submit with jobs.sh submit fastoma (it sets the log path), or
#   sbatch -o <log> shin_fastoma.sh [-resume]

set -eo pipefail
ROOT=/mnt/project/FjellheimLab/martpali/rcomplex-testrun
PROJ=$ROOT/complex_Shin2021/fastoma
PEO=/mnt/project/FjellheimLab/martpali/AnnualPerennial/peo-fastoma
mkdir -p "$PROJ/logs"
cd "$PROJ"

eval "$($HOME/.local/bin/micromamba shell hook -s bash)"
micromamba activate Nextflow
set -u
export NXF_HOME="$PROJ/.nextflow_home"
echo "Job ${SLURM_JOB_ID} on $(hostname), start $(date)"
echo "Nextflow: $(nextflow -version 2>&1 | grep version)"

export NXF_SINGULARITY_CACHEDIR=$PEO/singularity_cache
export NXF_WORK="$PROJ/work"
ulimit -n 131072 2>/dev/null || ulimit -n 65536 2>/dev/null || true

OMAMER_DB=/mnt/project/FjellheimLab/martpali/FastOMA/omamer_db/LUCA.h5
[[ -f "$OMAMER_DB" ]] || { echo "ERROR: $OMAMER_DB missing"; exit 1; }
echo "OMAmer DB: $OMAMER_DB"
echo "Proteomes: $(ls in_folder/proteome | wc -l)"

NF_CMD="nextflow run dessimozlab/FastOMA -r v0.5.1 \
    -profile slurm,singularity,large \
    -c shin_nextflow.config \
    --input in_folder \
    --omamer_db $OMAMER_DB \
    --output_folder output $*"
echo "Command: $NF_CMD"
eval $NF_CMD
EXIT_CODE=$?
echo "Nextflow exit: $EXIT_CODE, end $(date)"
exit $EXIT_CODE
