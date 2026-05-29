# **Environmental Variable Acquisition: Topography & Hidrology**
# DEM Processing Workflow for Neotropical Analyses

Developed by Kelly T. Bocanegra-González at the Bedoya Lab (NYBG), May 2026.

## LSDTopoTools and Lsdtopytools

To perform large-scale topographic analyses across the Neotropics, we installed and configured both `LSDTopoTools2` and `lsdtopytools`.

`LSDTopoTools2` is a high-performance geomorphometric analysis software, designed for extracting terrain metrics from Digital Elevation Models (DEMs). The software includes tools for calculating slope, curvature, roughness, drainage networks, hillshade, and other landscape attributes commonly used in geomorphology and ecological modelling.

In addition, `lsdtopytools` was installed as a Python interface for interacting with LSDTopoTools outputs and workflows from Python environments.

## 1. Install LSDTopoTools
https://lsdtopotools.github.io/

Clone: 

```cd ~/scratch/Podostemaceae/niche/fabdem git clone https://github.com/LSDtopotools/LSDTopoTools2.git```

### 1.2 Compile LSDTopoTools

```cd LSDTopoTools2```

```mkdir build```

```cd build```

```cmake ..```

```make -j 8```

### 1.3 Add LSDTopoTools to PATH
The LSDTopoTools executable directory was added to the system `PATH` to allow commandssuch as `lsdtt-basic-metrics` to be executed from any directory without specifying the full executable path.

```export PATH=~/scratch/Podostemaceae/niche/fabdem/LSDTopoTools2/src/build:$PATH```

### 1.4 Check GDAL installation

```gdalinfo --version```

### 1.5 Check available LSDTopoTools commands

```ls ~/scratch/Podostemaceae/niche/fabdem/LSDTopoTools2/src/build```

### 1.6 Confirm `lsdtt-basic-metrics` works

```lsdtt-basic-metrics -h```

This confirmed that LSDTopoTools was compiled correctly and that the executable can becalled from the terminal.

## 2. Install Lsdtopytools

A Python interface for LSDTopoTools (`lsdtopytools`) was also installed to facilitate terrain analysis workflows directly from Python environments.

Installation:

```conda activate lsdtopo```

```pip install lsdtopytools```

```conda activate lsdtopo```

Test installation:

```python -c "import lsdtopytools; print('lsdtopytools installed correctly')"```

`lsdtopytools` provides Python wrappers and utilities for interacting with LSDTopoTools workflows, enabling DEM preprocessing, terrain metric extraction, and integration with Python-based ecological and spatial analyses.

## 3. FABDEM Download
https://research-information.bris.ac.uk/en/datasets/fabdem-v1-2/

FABDEM (Forest And Buildings removed Copernicus DEM) is a global Digital Elevation Model (DEM) derived from the Copernicus GLO-30 dataset, in which the height bias produced by forests and buildings has been removed. The dataset provides a more accurate representation of bare-earth topography and is particularly useful for hydrological, geomorphological, ecological, and niche modelling analyses.
For this project, FABDEM V1-2 tiles covering the Neotropical region were downloaded and prepared for large-scale terrain analysis using LSDTopoTools.

- Extent used (Neotropics area): ```xmin, ymin, xmax, ymax = -120, -35, -30, 35``` 

- Downloaded all FABDEM V1-2 tiles intersecting the Neotropical region.

- Final dataset size: ```text ~46 GB```

script used: *fabdem.sh*

### 3.1 Extract FABDEM ZIP files

```cd ~/scratch/Podostemaceae/niche/fabdem```

```mkdir fabdem_tifs```

```for f in *.zip; do unzip -o "$f" -d fabdem_tifs done```

```find -name "*.tif" | wc -l```

```find . -name "*.tif" > tif_list.txt```

### 3.4 Create a virtual mosaic file (.vrt)

```gdalbuildvrt Neotropics_FABDEM_america.vrt -input_file_list tif_list_america.txt``` #This command creates a virtual mosaic (.vrt) that combines all the TIFF files listed in tif_list_america.txt without generating a new heavy raster file yet.

- output: Neotropics_FABDEM.vrt

```gdalinfo Neotropics_FABDEM.vrt``` #help to verify integity

```gdalinfo Neotropics_FABDEM.vrt | grep -A6 "Corner Coordinates"```#help to verify integrity

## 4. Generation of Neotropical tiles
## LSDTopoTools tile processing

To enable scalable continental analyses, the Neotropical region was subdivided into regular 5° × 5° tiles stored in a tab-delimited table *tiles_neotropics.tsv* generated with the script *tiles_coordinates.sh*.

For each tile, the Neotropical FABDEM virtual raster (Neotropics_FABDEM.vrt) was clipped and reprojected to the corresponding UTM zone in metric coordinates using GDAL (gdalwarp). Tiles were exported in ENVI/BIL format, which is required by LSDTopoTools.

Topographic surface metrics were computed using lsdtt-basic-metrics, including slope, aspect, general curvature, planform curvature, profile curvature, and hillshade rasters. Hydrological analyses were processed independently using a second workflow that calculated drainage area and channel network variables.

Processing was parallelized using SLURM array jobs:

```N=$(wc -l < tiles_neotropics.tsv)``` 

```sbatch --array=1-$(wc -l < tiles_neotropics.tsv)%10 lsd_neotropics_topo.sh```

Script: *lsd_neotropics_topo.sh* 

```N=$(wc -l < tiles_neotropics.tsv)``` 

```sbatch --array=1-${N}%10 lsd_neotropics_hydro.sh```

Script: *lsd_neotropics_hydro.sh*


## 5. Generation of continental topographic and hydrological mosaics

All raster variables generated independently for each 5° × 5° tile were reprojected to a common geographic coordinate system (WGS84; EPSG:4326) and merged into continental mosaics using GDAL (gdalwarp). Final rasters were generated at ~30 m spatial resolution (1 arc-second) and exported as compressed Cloud Optimized GeoTIFFs (COGs). The resulting mosaics included topographic variables (slope, aspect, hillshade, and curvature metrics) and hydrological flow accumulation variables (d8_area, dinf_area, MD_area, FMD_area, and QMD_area).

Note: ~30 m spatial resolution was used because the aim was to model the data alongside other raster layers that do not have a finer resolution.

Script: *mosaic_lsdtopo_1km.sh*

## 6. Soil variables extraction

Soil variables were obtained from the ISRIC SoilGrids global database at 0–5 cm depth (https://soilgrids.org/). Selected variables included bulk density (bdod), cation exchange capacity (cec), coarse fragments (cfvo), clay content (clay), nitrogen, pH in water (phh2o), silt content (silt), and soil organic carbon (soc). These variables were selected after exploratory analyses, including Pearson correlation and PCA, to retain ecologically informative and minimally redundant predictors. *(Notebook: 03_soil.ipy)*.

To enable scalable downloads across the Neotropics, the study area was subdivided into regular 5° × 5° geographic tiles in WGS84 geographic coordinates (EPSG:4326). Download tasks were automatically generated using the script `soilgrids_tables.sh`, which created a table containing all variable–tile combinations required for the analyses.

Raster downloads were then performed using SLURM array jobs through the script `soilgrids_download_5X5.sh`. Each array task downloaded a single SoilGrids variable for a single geographic tile using the SoilGrids Web Coverage Service (WCS), allowing efficient and reproducible large-scale retrieval of environmental rasters across the Neotropics. 

## 7. Generation of continental soil variable mosaics

Finally, use *mosaic_soildgrids_1km.sh* to creat the mosaic for variable.


## 8. CONTINUE WITH Notebook: 03_soil.ipy
