#!/bin/bash
#SBATCH --job-name=soilgrids

set -euo pipefail

BASE_DIR="$HOME/scratch/Podostemaceae/niche/soilgrids"
TASK_TABLE="${BASE_DIR}/soilgrids_download_tasks.tsv"
#TASK_TABLE="${BASE_DIR}/soilgrids_download_tasks_soil_fix.tsv"
OUTDIR="${BASE_DIR}/raw_tiles"

OUTDIR="${BASE_DIR}/raw_tiles"
LOGDIR="${BASE_DIR}/bad_tiles_logs"

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_TABLE")

VAR=$(echo "$LINE" | cut -f1)
TXMIN=$(echo "$LINE" | cut -f2)
TYMIN=$(echo "$LINE" | cut -f3)
TXMAX=$(echo "$LINE" | cut -f4)
TYMAX=$(echo "$LINE" | cut -f5)
TILE_ID=$(echo "$LINE" | cut -f6)

OUTFILE="${OUTDIR}/${VAR}_0_5cm_Q05_${TILE_ID}.tif"

URL="https://maps.isric.org/mapserv?map=/map/${VAR}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=${VAR}_0-5cm_Q0.5&FORMAT=GEOTIFF_INT16&SUBSET=X(${TXMIN},${TXMAX})&SUBSET=Y(${TYMIN},${TYMAX})&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"

echo "Downloading $VAR | $TILE_ID | EPSG:4326 extent: $TXMIN $TYMIN $TXMAX $TYMAX"

# If file already exists and is valid, skip
if [[ -s "$OUTFILE" ]]; then
  if gdalinfo "$OUTFILE" > /dev/null 2>&1; then
    echo "Valid file already exists: $OUTFILE"
    exit 0
  else
    echo "Existing file is invalid. Removing: $OUTFILE"
    rm -f "$OUTFILE"
  fi
fi

MAX_TRIES=2

for TRY in $(seq 1 "$MAX_TRIES"); do

  echo "Attempt $TRY of $MAX_TRIES"

  wget -q -O "$OUTFILE" "$URL"

  if gdalinfo "$OUTFILE" > /dev/null 2>&1; then
    echo "Completed valid TIFF: $OUTFILE"
    exit 0
  fi

  echo "Attempt $TRY failed: $OUTFILE"
  rm -f "$OUTFILE"

done

echo -e "${VAR}\t${TILE_ID}\t${TXMIN}\t${TYMIN}\t${TXMAX}\t${TYMAX}" >> "${LOGDIR}/${VAR}_bad_tiles.tsv"

echo "Persistent bad tile after ${MAX_TRIES} attempts: $OUTFILE"
exit 1
