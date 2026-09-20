#!/bin/bash
#SBATCH --job-name=soilgrids_1km

set -euo pipefail

source ~/.bashrc
conda activate lsdtopo

BASE_DIR=~/scratch/Podostemaceae/niche/
TILES_DIR=${BASE_DIR}/soilgrids/raw_tiles
OUT_DIR=${BASE_DIR}/soilgrids/variables_1km

mkdir -p "$OUT_DIR"

VARIABLES=(
  bdod
  cec
  cfvo
  clay
  nitrogen
  phh2o
  silt
  soc
)

VAR=${VARIABLES[$SLURM_ARRAY_TASK_ID]}

TMP_DIR="${OUT_DIR}/tmp_${VAR}_1km"
LIST_FILE="${OUT_DIR}/${VAR}_list.txt"
VRT_FILE="${OUT_DIR}/Neotropics_${VAR}_1km.vrt"
OUT_TIF="${OUT_DIR}/Neotropics_${VAR}_1km.tif"

mkdir -p "$TMP_DIR"

echo "======================================"
echo "Processing SoilGrids variable: ${VAR}"
echo "======================================"

find "$TILES_DIR" -name "${VAR}_0_5cm_Q05_*.tif" | sort > "$LIST_FILE"

echo "Number of ${VAR} tiles found:"
wc -l "$LIST_FILE"

if [ ! -s "$LIST_FILE" ]; then
  echo "No files found for ${VAR}. Exiting."
  exit 1
fi

rm -f "$VRT_FILE" "$OUT_TIF"
rm -f "$TMP_DIR"/*.tif

echo "Resampling SoilGrids tiles to 1 km EPSG:4326..."

i=0
while read -r TILE; do
  i=$((i + 1))

  BASENAME=$(basename "$TILE" .tif)
  TMP_TIF="${TMP_DIR}/${BASENAME}_1km_epsg4326.tif"

  echo "[$i] Processing $BASENAME"

  gdalwarp \
    -multi \
    -wo NUM_THREADS=ALL_CPUS \
    -t_srs EPSG:4326 \
    -tr 0.008333333333 0.008333333333 \
    -tap \
    -r average \
    -dstnodata -9999 \
    -ot Float32 \
    -of GTiff \
    -co TILED=YES \
    -co COMPRESS=DEFLATE \
    -co BIGTIFF=YES \
    -co PREDICTOR=3 \
    "$TILE" \
    "$TMP_TIF"

#    -srcnodata -32768 \

done < "$LIST_FILE"

echo "Building VRT from 1 km SoilGrids tiles..."

gdalbuildvrt \
  -srcnodata -9999 \
  -vrtnodata -9999 \
  "$VRT_FILE" \
  "$TMP_DIR"/*.tif

echo "Creating final Neotropics raster for ${VAR}..."

gdalwarp \
  -multi \
  -wo NUM_THREADS=ALL_CPUS \
  -t_srs EPSG:4326 \
  -te -120 -35 -30 35 \
  -te_srs EPSG:4326 \
  -tr 0.008333333333 0.008333333333 \
  -tap \
  -r average \
  -srcnodata -9999 \
  -dstnodata -9999 \
  -ot Float32 \
  -of GTiff \
  -co TILED=YES \
  -co COMPRESS=DEFLATE \
  -co BIGTIFF=YES \
  -co PREDICTOR=3 \
  "$VRT_FILE" \
  "$OUT_TIF"

echo "Checking final output:"
gdalinfo -stats "$OUT_TIF" | grep -E "Size is|EPSG|Minimum|Maximum|Mean|VALID_PERCENT|NoData|Pixel Size"

ls -lh "$OUT_TIF"

# Uncomment only if you want to delete intermediate files after each variable finishes.
 rm -rf "$TMP_DIR"
 rm -f "$VRT_FILE"

echo "Done with ${VAR}."

