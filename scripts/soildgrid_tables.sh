#!/bin/bash
#SBATCH --job-name=soilgrids
#SBATCH --output=soilgrids_download_%j.out
#SBATCH --error=soilgrids_download_%j.err
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G


set -euo pipefail

BASE_DIR="$HOME/scratch/Podostemaceae/niche/soilgrids"
TASK_TABLE="${BASE_DIR}/soilgrids_download_tasks.tsv"

mkdir -p "$BASE_DIR"

# Neotropical extent
# CRS: WGS84 / EPSG:4326
XMIN=-120
XMAX=-30
YMIN=-35
YMAX=35
TILE_SIZE=5

VARS=(
  bdod
  cec
  cfvo
  clay
  nitrogen
  phh2o
  silt
  soc
)

rm -f "$TASK_TABLE"

for VAR in "${VARS[@]}"; do
  for ((X=XMIN; X<XMAX; X+=TILE_SIZE)); do
    for ((Y=YMIN; Y<YMAX; Y+=TILE_SIZE)); do

      X2=$((X + TILE_SIZE))
      Y2=$((Y + TILE_SIZE))

      XLAB=$([[ $X -lt 0 ]] && echo "W${X#-}" || echo "E${X}")
      YLAB=$([[ $Y -lt 0 ]] && echo "S${Y#-}" || echo "N${Y}")

      TILE_ID="${XLAB}_${YLAB}"

      echo -e "${VAR}\t${X}\t${Y}\t${X2}\t${Y2}\t${TILE_ID}" >> "$TASK_TABLE"

    done
  done
done

echo "Task table created:"
echo "$TASK_TABLE"
echo "Number of tasks:"
wc -l "$TASK_TABLE"
