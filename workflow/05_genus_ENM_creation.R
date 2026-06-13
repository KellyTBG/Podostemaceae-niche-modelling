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

#1: Marathrum, 727 occurrences
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
marathrum_stats <- as.data.frame(eval_marathrum)
write.csv(marathrum_stats, "marathrum_enmeval_results.csv", row.names = FALSE)

#Plot best Marathrum model
eval_pred <- eval.predictions(enmeval_marathrum)

best_row <- which.min(eval_marathrum$AICc)

pred_marathrum <- eval_pred[[best_row]]

#Continuous
par(mfrow = c(1,2))
plot(pred_marathrum,
     main = paste("Marathrum Continuous | FC=", eval_marathrum$fc[best_row],
                  "RM=", eval_marathrum$rm[best_row]))

points(occs_marathrum, pch = 20, cex = 0.3, col = "red")

#Binary
threshold_10p <- eval_marathrum$or.10p.avg[best_row]

binary_marathrum <- terra::ifel(pred_marathrum >= threshold_10p, 1, 0)

plot(binary_marathrum,
     main = paste("Marathrum Binary (10%) | FC=",
                  eval_marathrum$fc[best_row],
                  "RM=", eval_marathrum$rm[best_row]))

points(occs_marathrum, pch = 20, cex = 0.3, col = "red")

#2. Podostemum, 475 occurrences
occs_podostemum <- soil_df %>%
  filter(Genus == "Podostemum") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_podostemum <- make_hull(occs_podostemum)
env_podostemum <- crop(env_stack_selected, hull_podostemum)
env_podostemum <- mask(env_podostemum, hull_podostemum)
#Largest extent, need to aggregate for RAM
env_podostemum <- aggregate(env_podostemum, fact = 2, fun = "mean")

enmeval_podostemum <- ENMevaluate(
  occs = occs_podostemum,
  envs = env_podostemum,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_podostemum <- eval.results(enmeval_podostemum)
print(eval_podostemum)
podostemum_stats <- as.data.frame(eval_podostemum)
write.csv(podostemum_stats, "podostemum_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_podostemum$AICc)
eval_pred <- eval.predictions(enmeval_podostemum)

pred_podostemum <- eval_pred[[best_row]]

threshold_10p <- eval_podostemum$or.10p.avg[best_row]
binary_podostemum <- terra::ifel(pred_podostemum >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_podostemum,
     main = paste0("Podostemum Continuous | FC=", eval_podostemum$fc[best_row],
                   " RM=", eval_podostemum$rm[best_row]))
points(occs_podostemum, pch=20, cex=0.3, col="red")

plot(binary_podostemum,
     main = "Podostemum Binary (10%)")
points(occs_podostemum, pch=20, cex=0.3, col="red")

#3. Apinagia, 401 occurrences
occs_apinagia <- soil_df %>%
  filter(Genus == "Apinagia") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_apinagia <- make_hull(occs_apinagia)
env_apinagia <- crop(env_stack_selected, hull_apinagia)
env_apinagia <- mask(env_apinagia, hull_apinagia)

enmeval_apinagia <- ENMevaluate(
  occs = occs_apinagia,
  envs = env_apinagia,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_apinagia <- eval.results(enmeval_apinagia)
print(eval_apinagia)
apinagia_stats <- as.data.frame(eval_apinagia)
write.csv(apinagia_stats, "apinagia_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_apinagia$AICc)
eval_pred <- eval.predictions(enmeval_apinagia)

pred_apinagia <- eval_pred[[best_row]]

threshold_10p <- eval_apinagia$or.10p.avg[best_row]
binary_apinagia <- terra::ifel(pred_apinagia >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_apinagia,
     main = paste0("Apinagia Continuous | FC=", eval_apinagia$fc[best_row],
                   " RM=", eval_apinagia$rm[best_row]))
points(occs_apinagia, pch=20, cex=0.3, col="red")

plot(binary_apinagia,
     main = "Apinagia Binary (10%)")
points(occs_apinagia, pch=20, cex=0.3, col="red")

#4. Tristicha, 381 occurrences
occs_tristicha <- soil_df %>%
  filter(Genus == "Tristicha") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_tristicha <- make_hull(occs_tristicha)
env_tristicha <- crop(env_stack_selected, hull_tristicha)
env_tristicha <- mask(env_tristicha, hull_tristicha)
#Second largest extent, need to aggregate for RAM
env_tristicha <- aggregate(env_tristicha, fact = 2, fun = "mean")

enmeval_tristicha <- ENMevaluate(
  occs = occs_tristicha,
  envs = env_tristicha,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE,
)

eval_tristicha <- eval.results(enmeval_tristicha)
print(eval_tristicha)
tristicha_stats <- as.data.frame(eval_tristicha)
write.csv(tristicha_stats, "tristicha_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_tristicha$AICc)
eval_pred <- eval.predictions(enmeval_tristicha)

pred_tristicha <- eval_pred[[best_row]]

threshold_10p <- eval_tristicha$or.10p.avg[best_row]
binary_tristicha <- terra::ifel(pred_tristicha >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_tristicha,
     main = paste0("Tristicha Continuous | FC=", eval_tristicha$fc[best_row],
                   " RM=", eval_tristicha$rm[best_row]))
points(occs_tristicha, pch=20, cex=0.3, col="red")

plot(binary_tristicha,
     main = "Tristicha Binary (10%)")
points(occs_tristicha, pch=20, cex=0.3, col="red")

#5. Mourera, 295 occurrences
occs_mourera <- soil_df %>%
  filter(Genus == "Mourera") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_mourera <- make_hull(occs_mourera)
env_mourera <- crop(env_stack_selected, hull_mourera)
env_mourera <- mask(env_mourera, hull_mourera)

enmeval_mourera <- ENMevaluate(
  occs = occs_mourera,
  envs = env_mourera,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_mourera <- eval.results(enmeval_mourera)
print(eval_mourera)
mourera_stats <- as.data.frame(eval_mourera)
write.csv(mourera_stats, "mourera_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_mourera$AICc)
eval_pred <- eval.predictions(enmeval_mourera)

pred_mourera <- eval_pred[[best_row]]

threshold_10p <- eval_mourera$or.10p.avg[best_row]
binary_mourera <- terra::ifel(pred_mourera >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_mourera,
     main = paste0("Mourera Continuous | FC=", eval_mourera$fc[best_row],
                   " RM=", eval_mourera$rm[best_row]))
points(occs_mourera, pch=20, cex=0.3, col="red")

plot(binary_mourera,
     main = "Mourera Binary (10%)")
points(occs_mourera, pch=20, cex=0.3, col="red")

#6. Rhyncholacis, 165 occurrences
occs_rhyncholacis <- soil_df %>%
  filter(Genus == "Rhyncholacis") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_rhyncholacis <- make_hull(occs_rhyncholacis)
env_rhyncholacis <- crop(env_stack_selected, hull_rhyncholacis)
env_rhyncholacis <- mask(env_rhyncholacis, hull_rhyncholacis)

enmeval_rhyncholacis <- ENMevaluate(
  occs = occs_rhyncholacis,
  envs = env_rhyncholacis,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_rhyncholacis <- eval.results(enmeval_rhyncholacis)
print(eval_rhyncholacis)
rhyncholacis_stats <- as.data.frame(eval_rhyncholacis)
write.csv(rhyncholacis_stats, "rhyncholacis_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_rhyncholacis$AICc)
eval_pred <- eval.predictions(enmeval_rhyncholacis)

pred_rhyncholacis <- eval_pred[[best_row]]

threshold_10p <- eval_rhyncholacis$or.10p.avg[best_row]
binary_rhyncholacis <- terra::ifel(pred_rhyncholacis >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_rhyncholacis,
     main = paste0("Rhyncholacis Continuous | FC=", eval_rhyncholacis$fc[best_row],
                   " RM=", eval_rhyncholacis$rm[best_row]))
points(occs_rhyncholacis, pch=20, cex=0.3, col="red")

plot(binary_rhyncholacis,
     main = "Rhyncholacis Binary (10%)")
points(occs_rhyncholacis, pch=20, cex=0.3, col="red")


#7. Castelnavia, 124 occurrences
occs_castelnavia <- soil_df %>%
  filter(Genus == "Castelnavia") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_castelnavia <- make_hull(occs_castelnavia)
env_castelnavia <- crop(env_stack_selected, hull_castelnavia)
env_castelnavia <- mask(env_castelnavia, hull_castelnavia)

enmeval_castelnavia <- ENMevaluate(
  occs = occs_castelnavia,
  envs = env_castelnavia,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_castelnavia <- eval.results(enmeval_castelnavia)
print(eval_castelnavia)
castelnavia_stats <- as.data.frame(eval_castelnavia)
write.csv(castelnavia_stats, "castelnavia_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_castelnavia$AICc)
eval_pred <- eval.predictions(enmeval_castelnavia)

pred_castelnavia <- eval_pred[[best_row]]

threshold_10p <- eval_castelnavia$or.10p.avg[best_row]
binary_castelnavia <- terra::ifel(pred_castelnavia >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_castelnavia,
     main = paste0("Castelnavia Continuous | FC=", eval_castelnavia$fc[best_row],
                   " RM=", eval_castelnavia$rm[best_row]))
points(occs_castelnavia, pch=20, cex=0.3, col="red")

plot(binary_castelnavia,
     main = "Castelnavia Binary (10%)")
points(occs_castelnavia, pch=20, cex=0.3, col="red")


#8. Noveloa, 73 occurrences
occs_noveloa <- soil_df %>%
  filter(Genus == "Noveloa") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_noveloa <- make_hull(occs_noveloa)
env_noveloa <- crop(env_stack_selected, hull_noveloa)
env_noveloa <- mask(env_noveloa, hull_noveloa)

enmeval_noveloa <- ENMevaluate(
  occs = occs_noveloa,
  envs = env_noveloa,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_noveloa <- eval.results(enmeval_noveloa)
print(eval_noveloa)
noveloa_stats <- as.data.frame(eval_noveloa)
write.csv(noveloa_stats, "noveloa_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_noveloa$AICc)
eval_pred <- eval.predictions(enmeval_noveloa)

pred_noveloa <- eval_pred[[best_row]]

threshold_10p <- eval_noveloa$or.10p.avg[best_row]
binary_noveloa <- terra::ifel(pred_noveloa >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_noveloa,
     main = paste0("Noveloa Continuous | FC=", eval_noveloa$fc[best_row],
                   " RM=", eval_noveloa$rm[best_row]))
points(occs_noveloa, pch=20, cex=0.3, col="red")

plot(binary_noveloa,
     main = "Noveloa Binary (10%)")
points(occs_noveloa, pch=20, cex=0.3, col="red")


#9. Lophogyne, 71 occurrences
occs_lophogyne <- soil_df %>%
  filter(Genus == "Lophogyne") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_lophogyne <- make_hull(occs_lophogyne)
env_lophogyne <- crop(env_stack_selected, hull_lophogyne)
env_lophogyne <- mask(env_lophogyne, hull_lophogyne)

enmeval_lophogyne <- ENMevaluate(
  occs = occs_lophogyne,
  envs = env_lophogyne,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_lophogyne <- eval.results(enmeval_lophogyne)
print(eval_lophogyne)
lophogyne_stats <- as.data.frame(eval_lophogyne)
write.csv(lophogyne_stats, "lophogyne_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_lophogyne$AICc)
eval_pred <- eval.predictions(enmeval_lophogyne)

pred_lophogyne <- eval_pred[[best_row]]

threshold_10p <- eval_lophogyne$or.10p.avg[best_row]
binary_lophogyne <- terra::ifel(pred_lophogyne >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_lophogyne,
     main = paste0("Lophogyne Continuous | FC=", eval_lophogyne$fc[best_row],
                   " RM=", eval_lophogyne$rm[best_row]))
points(occs_lophogyne, pch=20, cex=0.3, col="red")

plot(binary_lophogyne,
     main = "Lophogyne Binary (10%)")
points(occs_lophogyne, pch=20, cex=0.3, col="red")


#10. Weddellina, 68 occurrences
occs_weddellina <- soil_df %>%
  filter(Genus == "Weddellina") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_weddellina <- make_hull(occs_weddellina)
env_weddellina <- crop(env_stack_selected, hull_weddellina)
env_weddellina <- mask(env_weddellina, hull_weddellina)

enmeval_weddellina <- ENMevaluate(
  occs = occs_weddellina,
  envs = env_weddellina,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_weddellina <- eval.results(enmeval_weddellina)
print(eval_weddellina)
weddellina_stats <- as.data.frame(eval_weddellina)
write.csv(weddellina_stats, "weddellina_enmeval_results.csv", row.names = FALSE)

best_row <- which.min(eval_weddellina$AICc)
eval_pred <- eval.predictions(enmeval_weddellina)

pred_weddellina <- eval_pred[[best_row]]

threshold_10p <- eval_weddellina$or.10p.avg[best_row]
binary_weddellina <- terra::ifel(pred_weddellina >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_weddellina,
     main = paste0("Weddellina Continuous | FC=", eval_weddellina$fc[best_row],
                   " RM=", eval_weddellina$rm[best_row]))
points(occs_weddellina, pch=20, cex=0.3, col="red")

plot(binary_weddellina,
     main = "Weddellina Binary (10%)")
points(occs_weddellina, pch=20, cex=0.3, col="red")


#11. Oserya, 61 occurrences
occs_oserya <- soil_df %>%
  filter(Genus == "Oserya") %>%
  select(Longitude, Latitude) %>%
  filter(complete.cases(.))

hull_oserya <- make_hull(occs_oserya)
env_oserya <- crop(env_stack_selected, hull_oserya)
env_oserya <- mask(env_oserya, hull_oserya)

enmeval_oserya <- ENMevaluate(
  occs = occs_oserya,
  envs = env_oserya,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
  doClamp = TRUE,
  parallel = FALSE
)

eval_oserya <- eval.results(enmeval_oserya)
print(eval_oserya)

best_row <- which.min(eval_oserya$AICc)
eval_pred <- eval.predictions(enmeval_oserya)
oserya_stats <- as.data.frame(eval_oserya)
write.csv(oserya_stats, "oserya_enmeval_results.csv", row.names = FALSE)

pred_oserya <- eval_pred[[best_row]]

threshold_10p <- eval_oserya$or.10p.avg[best_row]
binary_oserya <- terra::ifel(pred_oserya >= threshold_10p, 1, 0)

par(mfrow=c(1,2))
plot(pred_oserya,
     main = paste0("Oserya Continuous | FC=", eval_oserya$fc[best_row],
                   " RM=", eval_oserya$rm[best_row]))
points(occs_oserya, pch=20, cex=0.3, col="red")

plot(binary_oserya,
     main = "Oserya Binary (10%)")
points(occs_oserya, pch=20, cex=0.3, col="red")

#Save best-model rasters to working directory
genera <- list(
  list(name = "Marathrum",    pred = pred_marathrum,    bin = binary_marathrum),
  list(name = "Podostemum",   pred = pred_podostemum,   bin = binary_podostemum),
  list(name = "Apinagia",     pred = pred_apinagia,     bin = binary_apinagia),
  list(name = "Tristicha",    pred = pred_tristicha,    bin = binary_tristicha),
  list(name = "Mourera",      pred = pred_mourera,      bin = binary_mourera),
  list(name = "Rhyncholacis", pred = pred_rhyncholacis, bin = binary_rhyncholacis),
  list(name = "Castelnavia",  pred = pred_castelnavia,  bin = binary_castelnavia),
  list(name = "Noveloa",      pred = pred_noveloa,      bin = binary_noveloa),
  list(name = "Lophogyne",    pred = pred_lophogyne,    bin = binary_lophogyne),
  list(name = "Weddellina",   pred = pred_weddellina,   bin = binary_weddellina),
  list(name = "Oserya",       pred = pred_oserya,       bin = binary_oserya)
)

for (g in genera) {
  writeRaster(
    g$pred,
    filename  = paste0(g$name, "_continuous_best.tif"),
    overwrite = TRUE,
    datatype  = "FLT4S" 
  )
  writeRaster(
    g$bin,
    filename  = paste0(g$name, "_binary10p_best.tif"),
    overwrite = TRUE,
    datatype  = "INT1U" 
  )
  message("Saved: ", g$name)
}

#Plotting hull outlines on Neotropics map for comparability
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

# Hull outlines on Neotropics map
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

land <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_crop(land, st_bbox(c(xmin=-120, xmax=-30, ymin=-40, ymax=35), crs=st_crs(4326)))

rivers <- ne_download(scale = 10, type = "rivers_lake_centerlines", 
                      category = "physical", returnclass = "sf")
rivers <- st_crop(rivers, st_bbox(c(xmin=-120, xmax=-30, ymin=-40, ymax=35), crs=st_crs(4326)))

genera <- c("Marathrum","Podostemum","Apinagia","Tristicha","Mourera",
            "Rhyncholacis","Castelnavia","Noveloa","Lophogyne","Weddellina","Oserya")

#Hull outlines overlayed on Neotropical river map
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

land <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_crop(land, st_bbox(c(xmin=-120, xmax=-30, ymin=-40, ymax=35), crs=st_crs(4326)))

rivers_primary <- ne_download(scale = 10, type = "rivers_lake_centerlines", 
                              category = "physical", returnclass = "sf")

rivers_secondary <- ne_download(scale = 10, type = "rivers_lake_centerlines_scale_rank", 
                                category = "physical", returnclass = "sf")

rivers_primary   <- st_crop(rivers_primary,   st_bbox(c(xmin=-120, xmax=-30, ymin=-40, ymax=35), crs=st_crs(4326)))
rivers_secondary <- st_crop(rivers_secondary, st_bbox(c(xmin=-120, xmax=-30, ymin=-40, ymax=35), crs=st_crs(4326)))

genera <- c("Marathrum","Podostemum","Apinagia","Tristicha","Mourera",
            "Rhyncholacis","Castelnavia","Noveloa","Lophogyne","Weddellina","Oserya")

genus_colors <- c(
  "Marathrum"    = "#E41A1C",
  "Podostemum"   = "#377EB8",
  "Apinagia"     = "#4DAF4A",
  "Tristicha"    = "#FF7F00",
  "Mourera"      = "#984EA3",
  "Rhyncholacis" = "#A65628",
  "Castelnavia"  = "#F781BF",
  "Noveloa"      = "#999999",
  "Lophogyne"    = "#FFFF33",
  "Weddellina"   = "#00CED1",
  "Oserya"       = "#2E8B57"
)

hull_list <- list()
for (g in genera) {
  pts <- soil_df %>%
    filter(Genus == g) %>%
    select(Longitude, Latitude) %>%
    filter(complete.cases(.))
  hull_list[[g]] <- make_hull(pts)
}

hull_sf_list <- lapply(hull_list, function(h) st_as_sf(h))

library(ggplot2)
hull_combined <- do.call(rbind, lapply(genera, function(g) {
  sf_obj <- hull_sf_list[[g]]
  sf_obj$Genus <- g
  sf_obj
}))

ggplot() +
  geom_sf(data = land, fill = "#f0ede4", color = "#aaaaaa", linewidth = 0.3) +
  geom_sf(data = rivers_secondary, color = "#7ab8d4", linewidth = 0.15) +
  geom_sf(data = rivers_primary,   color = "#7ab8d4", linewidth = 0.4) +
  geom_sf(data = hull_combined, 
          aes(color = Genus, fill = Genus),
          linewidth = 0.8, alpha = 0.08) +
  scale_color_manual(values = genus_colors) +
  scale_fill_manual(values = genus_colors) +
  coord_sf(xlim = c(-120, -30), ylim = c(-40, 30), expand = FALSE) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "#d6e8f5"),
    legend.position  = "right",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(size = 11)
  )

#Genus-level statistical test: does best validation AUC correlate with occurrence count?
library(dplyr)

#Occurrence counts
genus_n <- data.frame(
  Genus = c("Marathrum","Podostemum","Apinagia","Tristicha","Mourera",
            "Rhyncholacis","Castelnavia","Noveloa","Lophogyne","Weddellina","Oserya"),
  n_occs = c(727, 475, 401, 381, 295, 165, 124, 73, 71, 68, 61)
)

#Load all CSVs and extract best model AUC (lowest AICc)
stats_dir <- "/Users/lukesparreo/Desktop/PodostemaceaeENMs/ModelStats/"

genera_files <- c(
  "Marathrum"    = paste0(stats_dir, "marathrum_enmeval_results.csv"),
  "Podostemum"   = paste0(stats_dir, "podostemum_enmeval_results.csv"),
  "Apinagia"     = paste0(stats_dir, "apinagia_enmeval_results.csv"),
  "Tristicha"    = paste0(stats_dir, "tristicha_enmeval_results.csv"),
  "Mourera"      = paste0(stats_dir, "mourera_enmeval_results.csv"),
  "Rhyncholacis" = paste0(stats_dir, "rhyncholacis_enmeval_results.csv"),
  "Castelnavia"  = paste0(stats_dir, "castelnavia_enmeval_results.csv"),
  "Noveloa"      = paste0(stats_dir, "noveloa_enmeval_results.csv"),
  "Lophogyne"    = paste0(stats_dir, "lophogyne_enmeval_results.csv"),
  "Weddellina"   = paste0(stats_dir, "weddellina_enmeval_results.csv"),
  "Oserya"       = paste0(stats_dir, "oserya_enmeval_results.csv")
)

best_auc <- lapply(names(genera_files), function(g) {
  df <- read.csv(genera_files[g])
  best_row <- df[which.max(df$auc.val.avg), ]
  data.frame(
    Genus       = g,
    best_auc    = best_row$auc.val.avg,
    best_fc     = best_row$fc,
    best_rm     = best_row$rm
  )
}) %>% bind_rows()

#Merge with occurrence counts
auc_df <- left_join(best_auc, genus_n, by = "Genus")
print(auc_df)

#Spearman correlation
spearman_test <- cor.test(auc_df$n_occs, auc_df$best_auc, method = "spearman")
print(spearman_test)

#Plot
library(ggplot2)
library(ggrepel)

ggplot(auc_df, aes(x = n_occs, y = best_auc, label = Genus)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey60", linetype = "dashed", linewidth = 0.7) +
  geom_point(size = 3, color = "#2c7bb6") +
  geom_text_repel(size = 3, max.overlaps = 20, fontface = "italic") +
  annotate("text", 
           x = max(auc_df$n_occs) * 0.7, 
           y = min(auc_df$best_auc) + 0.01,
           label = paste0("Spearman ρ = ", round(spearman_test$estimate, 3),
                          "\np = ", round(spearman_test$p.value, 3)),
           size = 3.5, hjust = 0) +
  labs(x = "Number of Occurrences",
       y = "Best Model Validation AUC",
       title = "Effect of Occurrence Count on Model Performance")
