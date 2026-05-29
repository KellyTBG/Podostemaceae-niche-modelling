#!/bin/bash
#SBATCH --job-name=fabdem_tile_table
#SBATCH --output=fabdem_tile_table_%j.out
#SBATCH --error=fabdem_tile_table_%j.err
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G

# Create a table of 5-degree tiles covering the Neotropics.
# The table includes tile name, bounding box coordinates, and UTM EPSG code.


OUT="tiles_neotropics.tsv"
TILE_SIZE=5

rm -f "$OUT"

for SOUTH in $(seq -35 $TILE_SIZE 30); do
    NORTH=$((SOUTH + TILE_SIZE))

    for WEST in $(seq -120 $TILE_SIZE -35); do
        EAST=$((WEST + TILE_SIZE))

        CENTER_LON=$(( (WEST + EAST) / 2 ))
        CENTER_LAT=$(( (SOUTH + NORTH) / 2 ))

        ZONE=$(( (CENTER_LON + 180) / 6 + 1 ))

        if [ "$CENTER_LAT" -ge 0 ]; then
            EPSG=$((32600 + ZONE))
            LAT_TAG=$(printf "N%02d" ${CENTER_LAT#-})
        else
            EPSG=$((32700 + ZONE))
            LAT_TAG=$(printf "S%02d" ${CENTER_LAT#-})
        fi

        LON_TAG=$(printf "W%03d" ${WEST#-})
        TILE_NAME="tile_${LON_TAG}_${LAT_TAG}"

        echo -e "${TILE_NAME}\t${WEST}\t${SOUTH}\t${EAST}\t${NORTH}\t${EPSG}" >> "$OUT"
    done
done
