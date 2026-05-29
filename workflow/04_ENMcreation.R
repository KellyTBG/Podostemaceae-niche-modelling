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

#Load raster data
library(terra)
data_folder <- drive_ls(as_id("1m00pAcWo4Q3t2ei4YQCEeOvQUkSzQ28o"),
                        pattern = "02_data_processed")
climate_files <- drive_ls(as_id("1bSqUuhW4aRJ8ENqOttXHPRZRIW3Q1yYi"), pattern = ".tif")
fabdem_files  <- drive_ls(as_id("1nubDTxE7iExFszHD3HbxrDwCm7NyEU6B"), pattern = ".tif")
soil_files    <- drive_ls(as_id("1eTZCy-5kL7EgWRuF_ouQlml1teW6mEsL"), pattern = ".tif")

# Keep only .tif files (exclude .aux.xml)
climate_files <- climate_files[grepl("\\.tif$", climate_files$name), ]
fabdem_files  <- fabdem_files[grepl("\\.tif$", fabdem_files$name), ]
soil_files    <- soil_files[grepl("\\.tif$", soil_files$name), ]

# Create local folders
dir.create("rasters/climate",   recursive = TRUE, showWarnings = FALSE)
dir.create("rasters/fabdem",    recursive = TRUE, showWarnings = FALSE)
dir.create("rasters/soilgrids", recursive = TRUE, showWarnings = FALSE)

# Function to download a tif file
download_raster <- function(file_row, subfolder) {
  local_path <- file.path("rasters", subfolder, file_row$name)
  drive_download(as_id(file_row$id), path = local_path, overwrite = TRUE)
  local_path
}

# Download all tifs
climate_paths <- sapply(1:nrow(climate_files), function(i) 
  download_raster(climate_files[i, ], "climate"))
fabdem_paths  <- sapply(1:nrow(fabdem_files),  function(i) 
  download_raster(fabdem_files[i, ],  "fabdem"))
soil_paths    <- sapply(1:nrow(soil_files),    function(i) 
  download_raster(soil_files[i, ],    "soilgrids"))

# Stack each group
climate_stack <- rast(climate_paths)
fabdem_stack  <- rast(fabdem_paths)
soil_stack    <- rast(soil_paths)

print(climate_stack)
print(fabdem_stack)
print(soil_stack)

# Raster extents very slightly different due to rounding difference, crop to smallest
ref_ext <- ext(climate_stack)
fabdem_stack <- crop(fabdem_stack, ref_ext)
soil_stack   <- crop(soil_stack,   ref_ext)

# Stack all layers together
library(ENMeval)
env_stack <- c(climate_stack, fabdem_stack, soil_stack)

print(env_stack)
colnames(soil_df)

#Run models for testing
#Extract occurrence coordinates and selected variables
occs <- soil_df[, c("Longitude", "Latitude")]

selected_vars <- c(
  "wc2.1_30s_bio_1", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5",
  "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14",
  "Neotropics_SLOPE_1km", "Neotropics_metrics_CURV_1km",
  "Neotropics_metrics_PLFMCURV_1km", "Neotropics_metrics_PROFCURV_1km",
  "Neotropics_metrics_hs_1km",
  "Neotropics_hydro_MD_area_1km", "Neotropics_hydro_dinf_area_1km",
  "Neotropics_bdod_1km", "Neotropics_cec_1km", "Neotropics_cfvo_1km",
  "Neotropics_clay_1km", "Neotropics_nitrogen_1km", "Neotropics_phh2o_1km",
  "Neotropics_silt_1km", "Neotropics_soc_1km"
)

# Subset stack to selected variables
env_stack_selected <- env_stack[[selected_vars]]

print(env_stack_selected)

#Run ENMeval with block cross-validation and a range of regularization multipliers
terraOptions(memfrac = 0.9)  #Allow terra to use more available RAM
enmeval_results <- ENMevaluate(
  occs = occs,
  envs = env_stack_selected,
  algorithm = "maxnet",        
  partitions = "block",        #Block cross-validation
  tune.args = list(
    fc = c("L", "LQ", "LQH", "LQHP"),  #Testing feature classes: linear, quadratic, hinge
    rm = c(1, 2, 4)      #Testing regularization multipliers
  )
)

# Check results
eval.results(enmeval_results)

eval_tbl <- eval.results(enmeval_results)
best_aic <- eval_tbl[which.min(eval_tbl$AICc), ]
best_auc <- eval_tbl[which.max(eval_tbl$auc.val.avg), ]
best_or  <- eval_tbl[which.min(eval_tbl$or.10p.avg), ]

best_aic
best_auc
best_or

#Prepare data for running and plotting single best model (Testing: FC = L, RM = 2 and FC = LQ, RM = 4)
#Occurrence coords
occs <- soil_df[, c("Longitude", "Latitude")]

#Background sampling directly from raster
bg <- spatSample(
  env_stack_selected[[1]],
  size = 10000,
  method = "random",
  xy = TRUE,
  na.rm = TRUE
)

#Extract environmental values for points
occs_env <- extract(env_stack_selected, occs)[, -1]

bg_env <- extract(
  env_stack_selected,
  bg[, c("x", "y")]
)[, -1]

#Remove NAs (this is done automatically by ENMEval, but not in Maxnet)
occs_keep <- complete.cases(occs_env)
bg_keep   <- complete.cases(bg_env)

occs_env <- occs_env[occs_keep, ]
bg_env   <- bg_env[bg_keep, ]

#Create presence/absence vector
p <- c(
  rep(1, nrow(occs_env)),
  rep(0, nrow(bg_env))
)

dat <- rbind(occs_env, bg_env)

#Maxnet model (Testing: FC = L, RM = 2 and FC = LQ, RM = 4)
library(maxnet)
maxnet_model <- maxnet(
  p = p,
  data = dat,
  f = maxnet.formula(p, dat, classes = "lq"),
  regmult = 4
)

#Predict across raster
prediction <- predict(
  env_stack_selected,
  maxnet_model,
  type = "cloglog",
  na.rm = TRUE
)

#Plot continuous map
par(mfrow = c(1,1)) 
plot(prediction, main= "Continuous Podostemaceae ENM (FC= LQ, RM= 4, AUC= 0.866, OR10P= 0.162)")
points(occs, pch = 20, cex = 0.1, col = "red")


#Binary thresholding, minimum 10% training presence
binary_map <- prediction >= thresh_10p
plot(binary_map,
     col = c("white", "darkgreen"),
     legend = FALSE, main= "Binary Podostemaceae ENM (FC= LQ, RM= 4, AUC= 0.866, OR10P= 0.162)")

points(occs, pch = 20, cex = 0.1, col = "red")


#Binary thresholding, minimum 20% training presence (stricter suitability)
thresh_20p <- quantile(
  extract(prediction, occs)[,2],
  0.20,
  na.rm = TRUE
)
binary_map <- prediction >= thresh_20p
plot(binary_map,
     col = c("white", "darkgreen"),
     legend = FALSE, main= "Binary Podostemaceae ENM (FC= LQ, RM= 4, AUC= 0.866, OR10P= 0.162)")

points(occs, pch = 20, cex = 0.1, col = "red")

#Examine variable response curves for selected model
par(mfrow = c(5,5))

for (v in names(dat)) {
  
  x <- seq(min(dat[[v]], na.rm=TRUE),
           max(dat[[v]], na.rm=TRUE),
           length.out=100)
  
  newdata <- data.frame(matrix(
    apply(dat, 2, mean, na.rm=TRUE),
    nrow=100,
    ncol=ncol(dat),
    byrow=TRUE
  ))
  
  colnames(newdata) <- colnames(dat)
  newdata[[v]] <- x
  
  y <- predict(maxnet_model, newdata, type="cloglog")
  
  plot(x, y, type="l", main=v, ylab="suitability")
}
