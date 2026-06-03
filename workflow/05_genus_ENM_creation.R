#Genus-Level Podostemaceae ENM Evaluation, Selection, and Plotting
#Luke Sparreo/Bedoya Lab, June 2, 2026

#Connect to Drive to import selected occurrences, variables, and raster data
library(googledrive)
drive_auth(scopes = "https://www.googleapis.com/auth/drive")

#Load occurrence and selected environmental variable matrix
file_info <- drive_ls(
  as_id("1PVws26_O5PI7YMp1zr5ecw9BAIrOxyqW"),
  pattern = "Podostemaceae_environmental_matrix_bio_topo_soil_selected.csv"
)
drive_download(
  as_id(file_info$id),
  path = "Podostemaceae_matrix.csv",
  overwrite = TRUE
)
soil_df <- read.csv("Podostemaceae_matrix.csv")
View(soil_df)

#Examine dataframe by genus
library(dplyr)
library(stringr)
soil_df <- soil_df %>%
  mutate(Genus = str_extract(Species, "^\\S+"))

#Count by genus
genus_counts <- soil_df %>%
  count(Genus, sort = TRUE)

print(genus_counts)

#Load raster data for selected variables
library(terra)
selected_vars <- c(
  "wc2.1_30s_bio_1", "wc2.1_30s_bio_4",
  "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14",
  "Neotropics_SLOPE_1km", "Neotropics_metrics_CURV_1km", "Neotropics_metrics_PROFCURV_1km",
  "Neotropics_hydro_MD_area_1km", "Neotropics_hydro_dinf_area_1km",
  "Neotropics_bdod_1km", "Neotropics_cec_1km", "Neotropics_cfvo_1km",
  "Neotropics_clay_1km", "Neotropics_phh2o_1km",
  "Neotropics_silt_1km"
)

#List files from Drive
climate_files <- drive_ls(as_id("1bSqUuhW4aRJ8ENqOttXHPRZRIW3Q1yYi"), pattern = ".tif")
fabdem_files  <- drive_ls(as_id("1nubDTxE7iExFszHD3HbxrDwCm7NyEU6B"), pattern = ".tif")
soil_files    <- drive_ls(as_id("1eTZCy-5kL7EgWRuF_ouQlml1teW6mEsL"), pattern = ".tif")

#Keep only .tif files
climate_files <- climate_files[grepl("\\.tif$", climate_files$name), ]
fabdem_files  <- fabdem_files[grepl("\\.tif$",  fabdem_files$name), ]
soil_files    <- soil_files[grepl("\\.tif$",    soil_files$name), ]

#Filter to only selected variables
climate_files <- climate_files[gsub("\\.tif$", "", climate_files$name) %in% selected_vars, ]
fabdem_files  <- fabdem_files[gsub("\\.tif$",  "", fabdem_files$name)  %in% selected_vars, ]
soil_files    <- soil_files[gsub("\\.tif$",    "", soil_files$name)    %in% selected_vars, ]

climate_paths <- sapply(1:nrow(climate_files), function(i) {
  local_path <- file.path("/Users/lukesparreo/rasters/climate", climate_files$name[i])
  if (!file.exists(local_path)) drive_download(as_id(climate_files$id[i]), path = local_path, overwrite = FALSE)
  local_path
})
fabdem_paths <- sapply(1:nrow(fabdem_files), function(i) {
  local_path <- file.path("/Users/lukesparreo/rasters/fabdem", fabdem_files$name[i])
  if (!file.exists(local_path)) drive_download(as_id(fabdem_files$id[i]), path = local_path, overwrite = FALSE)
  local_path
})
soil_paths <- sapply(1:nrow(soil_files), function(i) {
  local_path <- file.path("/Users/lukesparreo/rasters/soilgrids", soil_files$name[i])
  if (!file.exists(local_path)) drive_download(as_id(soil_files$id[i]), path = local_path, overwrite = FALSE)
  local_path
})

climate_stack <- rast(climate_paths)
fabdem_stack  <- rast(fabdem_paths)
soil_stack    <- rast(soil_paths)

#Raster extents very slightly different due to rounding difference, crop to smallest
ref_ext <- ext(climate_stack)
fabdem_stack <- crop(fabdem_stack, ref_ext)
soil_stack   <- crop(soil_stack,   ref_ext)

#Stack selected rasters
env_stack_selected <- c(
  climate_stack,
  fabdem_stack,
  soil_stack
)

print(env_stack_selected)

#Extract occurrence coordinates
occs <- soil_df[, c("Longitude", "Latitude")]

#Create buffered convex hulls with function using terra for each genus
library(sf)
library(dplyr)
make_hull <- function(pts) {
  occ_vect    <- vect(pts, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  hull        <- convHull(occ_vect)
  hull_ea     <- project(hull, "EPSG:6933")          # equal-area projection
  hull_buf_ea <- buffer(hull_ea, width = 111320)     # 1 degree ~ 111.32 km
  hull_buf    <- project(hull_buf_ea, "EPSG:4326")   # back to WGS84
  return(hull_buf)
}

#Begin running models
library(ENMeval)

#1: Marathrum
#Filter dataset
occs_marathrum <- soil_df %>%
  filter(Genus == "Marathrum") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))
dim(occs_marathrum)

#Create 1-degree buffer around coordinates
hull_marathrum <- make_hull(occs_marathrum)
env_marathrum <- crop(env_stack_selected, hull_marathrum)
env_marathrum <- mask(env_marathrum, hull_marathrum)

#Run ENM
enmeval_marathrum <- ENMevaluate(
  occs = occs_marathrum,
  envs = env_marathrum,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(
    fc = c("L", "LQ", "LQH", "LQHP"),
    rm = c(1, 2, 4)
  ),
  doClamp = TRUE,
  parallel = FALSE
)

#Examine models
eval_marathrum <- eval.results(enmeval_marathrum)
eval_marathrum
best_marathrum <- eval_marathrum[which.min(eval_marathrum$AICc), ]

#Plot best Marathrum model
eval_pred <- eval.predictions(enmeval_marathrum)

best_row <- which.min(eval_marathrum$AICc)

pred_marathrum <- eval_pred[[best_row]]

#Continuous
par(mfrow = c(1,2))
plot(pred_marathrum,
     main = paste("Marathrum Continuous | FC=", eval_marathrum$fc[best_row],
                  "RM=", eval_marathrum$rm[best_row],
                  "AUC=", round(eval_marathrum$auc.val.avg[best_row], 3)))

points(occs_marathrum, pch = 20, cex = 0.3, col = "red")

#Binary
threshold_10p <- eval_marathrum$or.10p.avg[best_row]

binary_marathrum <- terra::ifel(pred_marathrum >= threshold_10p, 1, 0)

plot(binary_marathrum,
     main = paste("Marathrum Binary (10%) | FC=",
                  eval_marathrum$fc[best_row],
                  "RM=", eval_marathrum$rm[best_row]))

points(occs_marathrum, pch = 20, cex = 0.3, col = "red")

#CONTINUE HERE WITH OTHER GENERA
