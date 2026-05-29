#!/bin/bash
#SBATCH --job-name=lsd_hydro
#SBATCH --output=lsd_hydro_%A_%a.out
#SBATCH --error=lsd_hydro_%A_%a.err
#SBATCH --partition=himem
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --mail-user=xx@xx
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source ~/.bashrc
conda activate lsdtopo

export PATH=/path/to/LSDTopoTools2/src/build:$PATH
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}
ulimit -s unlimited

BASE_DIR=~/scratch/Podostemaceae/niche/fabdem
TILE_TABLE=${BASE_DIR}/tiles_neotropics.tsv
OUT_BASE=${BASE_DIR}/tiles_surface

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TILE_TABLE")

TILE_NAME=$(echo "$LINE" | cut -f1)

OUT_DIR=${OUT_BASE}/${TILE_NAME}
cd "$OUT_DIR"

# Skip completed hydro tiles
if [ -f ${TILE_NAME}_hydro_CN.csv ] || [ -f ${TILE_NAME}_hydro_CN.geojson ]; then
    echo "${TILE_NAME} hydro already processed. Skipping."
    exit 0
fi

cat > ${TILE_NAME}_hydro.driver <<EOF
read fname: ${TILE_NAME}_utm
write fname: ${TILE_NAME}_hydro
channel heads fname: NULL

remove_seas: true
minimum_elevation: 0

min_slope_for_fill: 0.0001

print_d8_drainage_area_raster: true
print_dinf_drainage_area_raster: true
print_QuinnMD_drainage_area_raster: true
print_FreemanMD_drainage_area_raster: true
print_MD_drainage_area_raster: true

threshold_contributing_pixels: 1000
print_stream_order_raster: true
print_channels_to_csv: true
print_junctions_to_csv: true
print_junctions_angles_to_csv: true

convert_csv_to_geojson: true
EOF

lsdtt-basic-metrics ${TILE_NAME}_hydro.driver
