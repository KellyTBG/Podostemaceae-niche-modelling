#!/bin/bash
#SBATCH --job-name=lsdtopo_vars_1km
#SBATCH --output=lsdtopo_vars_1km_%A_%a.out
#SBATCH --error=lsdtopo_vars_1km_%A_%a.err
#SBATCH --partition=himem
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --mail-user=xx@xx
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source ~/.bashrc
conda activate lsdtopo

BASE_DIR=~/scratch/Podostemaceae/niche/fabdem
TILES_DIR=${BASE_DIR}/tiles_surface
OUT_DIR=${BASE_DIR}/variables_1km

mkdir -p "$OUT_DIR"

VARIABLES=(
  metrics_ASPECT
  metrics_CURV
  metrics_hs
  metrics_PLFMCURV
  metrics_PROFCURV
  hydro_d8_area
  hydro_dinf_area
  hydro_FMD_area
  hydro_MD_area
  hydro_QMD_area
)

VAR=${VARIABLES[$SLURM_ARRAY_TASK_ID]}

TMP_DIR="${OUT_DIR}/tmp_${VAR}_1km"
LIST_FILE="${OUT_DIR}/${VAR}_list.txt"
VRT_FILE="${OUT_DIR}/Neotropics_${VAR}_1km.vrt"
OUT_TIF="${OUT_DIR}/Neotropics_${VAR}_1km.tif"

mkdir -p "$TMP_DIR"

echo "======================================"
echo "Processing variable: ${VAR}"
echo "======================================"

find "$TILES_DIR" -name "*_${VAR}.bil" | sort > "$LIST_FILE"

echo "Number of ${VAR} tiles found:"
wc -l "$LIST_FILE"

if [ ! -s "$LIST_FILE" ]; then
  echo "No files found for ${VAR}. Exiting."
  exit 1
fi

rm -f "$VRT_FILE" "$OUT_TIF"
rm -f "$TMP_DIR"/*.tif

echo "Reprojecting tiles to EPSG:4326 and aggregating to 1 km..."

i=0
while read -r TILE; do
  i=$((i + 1))

  BASENAME=$(basename "$TILE" .bil)
  TMP_TIF="${TMP_DIR}/${BASENAME}_1km_epsg4326.tif"

  echo "[$i] Processing $BASENAME"

  # This step converts each tile from its original UTM zone to EPSG:4326
  # and resamples it to 30 arc-seconds, equivalent to WorldClim 30s (~1 km).

  gdalwarp \
    -multi \
    -wo NUM_THREADS=ALL_CPUS \
    -t_srs EPSG:4326 \
    -tr 0.008333333333 0.008333333333 \
    -tap \
    -r average \
    -srcnodata -9999 \
    -dstnodata -9999 \
    -of GTiff \
    -co TILED=YES \
    -co COMPRESS=DEFLATE \
    -co BIGTIFF=YES \
    -co PREDICTOR=3 \
    "$TILE" \
    "$TMP_TIF"

done < "$LIST_FILE"

echo "Building VRT from 1 km EPSG:4326 tiles..."

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

# Uncomment these two lines only if you want to delete intermediate files after each variable finishes.
# rm -rf "$TMP_DIR"
# rm -f "$VRT_FILE"

echo "Done with ${VAR}."
