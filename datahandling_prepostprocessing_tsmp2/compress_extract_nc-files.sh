#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --account=slts
#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --job-name=nc-archive
#SBATCH --time=01:25:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage:"
  echo "  sbatch compress_extract_nc-files.sh compress folder*"
  echo "  sbatch compress_extract_nc-files.sh extract  folder*.tar"
  exit 1
fi

MODE="$1"
shift
INPUTS=("$@")

MAX_PARALLEL=${SLURM_NTASKS:-1}

# ---- timer start ----
START_TIME=$(date +%s)
echo "Job started at: $(date)"
echo "Mode: $MODE"
echo "Parallel jobs: $MAX_PARALLEL"
echo "Inputs:"
printf "  %s\n" "${INPUTS[@]}"

case "$MODE" in
  compress)
    echo "Starting gzip step..."

    find "${INPUTS[@]}" -type f -name "*.nc" -print0 \
      | xargs -0 -n 1 -P "$MAX_PARALLEL" gzip

    echo "Gzip finished, starting tar step..."

    for dir in "${INPUTS[@]}"; do
      absdir=$(realpath "$dir")
      parentdir=$(dirname "$absdir")
      basename_dir=$(basename "$absdir")

#      tar -cf "${dir}.tar" "$dir"
      tar -C "$parentdir" -cf "$parentdir/${basename_dir}.tar" "$basename_dir"
    done

    ;;

  extract)
    echo "Starting untar step..."

    for tarfile in "${INPUTS[@]}"; do
      abs_tar=$(realpath "$tarfile")
      tar_dir=$(dirname "$abs_tar")
      tar_base=$(basename "$abs_tar")

#      tar -xf "$tarfile"
      tar -C "$tar_dir" -xf "$abs_tar"
    done

    echo "Untar finished, starting gunzip step..."

    find . -type f -name "*.nc.gz" -print0 \
      | xargs -0 -n 1 -P "$MAX_PARALLEL" gunzip

    ;;

  *)
    echo "ERROR: Unknown mode '$MODE' (use compress or extract)"
    exit 1
    ;;
esac

# ---- timer end ----
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

printf "Job finished at: %s\n" "$(date)"
printf "Total runtime: %02d:%02d:%02d (hh:mm:ss)\n" \
  $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
