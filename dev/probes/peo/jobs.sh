#!/bin/bash
# Dispatcher for the PEO benchmark of rcomplex 0.4.2 on Orion.
#
# Local (Mac):   bash dev/probes/peo/jobs.sh push
#   builds the tarball into prepare_data/cache/peo/ and copies it and
#   these scripts to $ROOT/peo/.
# Login node:    ssh orion 'module load slurm;
#                  bash $ROOT/peo/scripts/jobs.sh submit <step> [args]'
#   submits one job and prints its id. Steps and requests:
#     install                 4 cores,  16 GB,  1 h   -> $ROOT/peo/
#     hogs                    4 cores,  64 GB,  4 h   -> $ROOT/peo/
#     stage <exp>             4 cores,  64 GB,  2 h
#     run <exp> <run>        32 cores, 300 GB, 24 h
#     score <exp>             4 cores,  64 GB,  2 h
#   exp: LeafvRoot, C3vC4 or nodulation. run: leaf, root or organ_block.
#   Logs: $ROOT/complex_<exp>/logs/<step>-%j.out.
# Inside a job the same script runs the step; it refuses to run a step
# outside SLURM, so nothing heavy lands on the login node.
set -euo pipefail
ROOT=/mnt/project/FjellheimLab/martpali/rcomplex-testrun
P=$ROOT/peo
SELF=$P/scripts/jobs.sh
TAR=rcomplex_0.4.2.tar.gz
O=/mnt/project/FjellheimLab/martpali/AnnualPerennial/peo-fastoma/output
LIB=$P/lib-0.4.2
DEPS=$P/lib-deps
OLD_DEPS=$ROOT/validation-2026-09-17/lib-deps
HOME_LIB=/mnt/users/martpali/R/library-4.4

push() {
  local repo=/Users/martin/Documents/R/rcomplex/.claude/worktrees/peo
  local cache=/Users/martin/Documents/R/rcomplex/prepare_data/cache/peo
  Rscript -e "pkgbuild::build('$repo', dest_path = '$cache')"
  ssh -o BatchMode=yes orion "mkdir -p $P/scripts $P/pkg"
  scp -o BatchMode=yes "$repo"/dev/probes/peo/* orion:$P/scripts/
  scp -o BatchMode=yes "$cache/$TAR" orion:$P/pkg/
}

submit() {
  local step=$1 res log
  shift
  case $step in
    install) res="-c 4 --mem=16G -t 01:00:00"; log=$P/install-%j.out ;;
    hogs) res="-c 4 --mem=64G -t 04:00:00"; log=$P/hogs-%j.out ;;
    stage) res="-c 4 --mem=64G -t 02:00:00" ;;
    run) res="-c 32 --mem=300G -t 24:00:00" ;;
    score) res="-c 4 --mem=64G -t 02:00:00" ;;
    *) echo "step: unknown '$step'" >&2; exit 2 ;;
  esac
  if [ -z "${log:-}" ]; then
    log=$ROOT/complex_${1:?exp}/logs/$step${2:+-$2}-%j.out
  fi
  mkdir -p "$(dirname "$log")"
  # shellcheck disable=SC2086
  sbatch --parsable -p orion -A nmbu -J "peo-$step${1:+-$1}${2:+-$2}" \
    $res --mail-type=FAIL -o "$log" "$SELF" "$step" "$@"
}

in_job() {
  if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo "step $1 runs only inside a SLURM job; use: submit $1" >&2
    exit 1
  fi
  set +u
  module load R/4.4.2
  set -u
  export R_LIBS="$LIB:$DEPS:$OLD_DEPS:$HOME_LIB"
  export PEO_ROOT=$ROOT
}

case "${1:-}" in
  push) push ;;
  submit) shift; submit "$@" ;;
  install)
    in_job install
    mkdir -p "$LIB" "$DEPS"
    # ~/.R/Makevars points at a conda compiler; the container has g++.
    : > "$P/Makevars.empty"
    export R_MAKEVARS_USER=$P/Makevars.empty
    cat > "$P/deps.R" <<'REOF'
args <- commandArgs(TRUE)
tmp <- file.path(dirname(args[1]), "desc")
untar(args[1], files = "rcomplex/DESCRIPTION", exdir = tmp)
d <- read.dcf(file.path(tmp, "rcomplex", "DESCRIPTION"))
imp <- gsub("\\s+", " ", paste(d[1, "Imports"], d[1, "LinkingTo"],
  sep = ","))
imp <- trimws(sub("[(].*$", "", strsplit(imp, ",")[[1]]))
imp <- c(imp[nzchar(imp)], "ape", "R.utils")
need <- setdiff(imp, rownames(installed.packages()))
cat("imports:", imp, "\nmissing:", if (length(need)) need else "none",
  "\n")
if (length(need)) {
  install.packages(need, lib = args[2],
    repos = "https://cloud.r-project.org")
}
REOF
    Rscript "$P/deps.R" "$P/pkg/$TAR" "$DEPS"
    R CMD INSTALL --library="$LIB" --no-multiarch "$P/pkg/$TAR"
    # Files R writes in a setgid dir must keep the project group.
    mkdir -p "$P/tmp"
    Rscript -e 'writeLines("x", file.path(commandArgs(TRUE), "gid.txt"))' \
      "$P/tmp"
    echo "group of a file R wrote: $(stat -c %G "$P/tmp/gid.txt")"
    Rscript -e 'cat("rcomplex", as.character(packageVersion("rcomplex")),
      find.package("rcomplex"), "\n")' ;;
  hogs)
    in_job hogs
    python3 "$P/scripts/hogs_at_level.py" "$O/FastOMA_HOGs.orthoxml" \
      "$O/species_tree_checked.nwk" "$P/scripts/sets.tsv" "$ROOT" ;;
  stage)
    in_job stage
    # fread unpacks .gz into TMPDIR: keep SLURM's node-local one. Never
    # a dir under $ROOT: R chmods what it unpacks there, which clears
    # setgid and charges the files to the personal group (quota).
    echo "TMPDIR=${TMPDIR:-unset}"
    Rscript "$P/scripts/stage.R" "$2" ;;
  run)
    in_job run
    Rscript "$P/scripts/run.R" "$2" "$3" ;;
  score)
    in_job score
    E=$ROOT/complex_$2
    mkdir -p "$E/score"
    set +u
    module load slurm
    set -u
    for f in "$E"/runs/*/timing.tsv; do
      id=$(tail -1 "$f" | cut -f4)
      sacct -j "$id" -P -n --format=JobID,State,Elapsed,MaxRSS || true
    done > "$E/score/sacct.tsv"
    Rscript "$P/scripts/score.R" "$2" ;;
  *) sed -n '2,19p' "$0"; exit 2 ;;
esac
