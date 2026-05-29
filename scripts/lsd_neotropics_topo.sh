#!/bin/bash
#SBATCH --job-name=lsd_neotropics_topo
#SBATCH --output=lsd_neotropics_topo_%A_%a.out
#SBATCH --error=lsd_neotropics_topo_%A_%a.err
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
VRT=${BASE_DIR}/Neotropics_FABDEM.vrt
TILE_TABLE=${BASE_DIR}/tiles_neotropics.tsv
OUT_BASE=${BASE_DIR}/tiles_surface

mkdir -p "$OUT_BASE"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TILE_TABLE")

TILE_NAME=$(echo "$LINE" | cut -f1)
WEST=$(echo "$LINE" | cut -f2)
SOUTH=$(echo "$LINE" | cut -f3)
EAST=$(echo "$LINE" | cut -f4)
NORTH=$(echo "$LINE" | cut -f5)
UTM_EPSG=$(echo "$LINE" | cut -f6)

OUT_DIR=${OUT_BASE}/${TILE_NAME}
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

if [ -f ${TILE_NAME}_metrics_SLOPE.bil ]; then
    echo "${TILE_NAME} already processed. Skipping."
    exit 0
fi

gdalwarp \
-te_srs EPSG:4326 \
-te ${WEST} ${SOUTH} ${EAST} ${NORTH} \
-t_srs EPSG:${UTM_EPSG} \
-tr 30 30 \
-r bilinear \
-dstnodata -9999 \
-of ENVI \
-co INTERLEAVE=BIL \
"$VRT" \
${TILE_NAME}_utm.bil

if [ -f ${TILE_NAME}_utm.bil.hdr ]; then
    mv ${TILE_NAME}_utm.bil.hdr ${TILE_NAME}_utm.hdr
fi

gdalinfo ${TILE_NAME}_utm.bil | grep "Size is"

cat > ${TILE_NAME}_surface_metrics.driver <<EOF
read fname: ${TILE_NAME}_utm
write fname: ${TILE_NAME}_metrics

remove_seas: true
minimum_elevation: 0

write_hillshade: true
print_slope: true
print_aspect: true
print_curvature: true
print_planform_curvature: true
print_profile_curvature: true
EOF

lsdtt-basic-metrics ${TILE_NAME}_surface_metrics.driver