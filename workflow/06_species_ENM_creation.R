#Species-Level Podostemaceae ENM Evaluation, Selection, and Plotting
#Luke Sparreo/Bedoya Lab, June 13-30, 2026

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

#Curate dataset per taxonomic revisions and occurance errors
soil_df$Species[soil_df$Species == "Marathrum capillaceum"] <- "Lophogyne capillacea"
soil_df$Species[soil_df$Species == "Lonchostephus elegans"] <- "Mourera elegans"
soil_df$Species[soil_df$Species == "Lophogyne aeruginosa"]  <- "Lophogyne penicillata"
soil_df$Species[soil_df$Species == "Apinagia riedelii"]     <- "Apinagia fucoides"
soil_df$Species[soil_df$Species == "Apinagia batrachifolia"] <- "	Apinagia batrachiifolia"
soil_df$Species[soil_df$Species == "Apinagia yguazuensis"] <- "Apinagia fucoides"
soil_df$Species[soil_df$Species == "Apinagia secundiflora"]  <- "Apinagia richardiana"
soil_df$Species[soil_df$Species == "Apinagia corymbosa"]     <- "Apinagia richardiana"
soil_df$Species[soil_df$Species == "Apinagia exilis"]  <- "Apinagia richardiana"
soil_df$Species[soil_df$Species == "Podostemum flagelliforme"]     <- "Devillea flagelliformis"

soil_df <- soil_df[!soil_df$GBIF_code %in% c(1260131418, 1260131425, 1844433356, 4061593219, 1322407602, 3068331619), ]
View(soil_df)

#Examine dataframe by genus
library(dplyr)
library(stringr)
soil_df <- soil_df %>%
  mutate(Genus = str_extract(Species, "^\\S+"))

genus_counts <- soil_df %>%
  count(Genus, sort = TRUE)
print(genus_counts)

species_counts <- soil_df %>%
  count(Species, sort = TRUE)
print(species_counts)

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

#Species-level ENMs, can be run as a loop because of smaller distributions not requiring RAM reset
make_hull <- function(pts) {
  occ_vect    <- vect(pts, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  hull        <- convHull(occ_vect)
  hull_ea     <- project(hull, "EPSG:6933")
  
  # Training extent: 1-degree buffer
  hull_buf_ea  <- buffer(hull_ea, width = 111320)
  hull_buf     <- project(hull_buf_ea, "EPSG:4326")
  
  # Projection extent: 5-degree buffer
  hull_proj_ea <- buffer(hull_ea, width = 5 * 111320)
  hull_proj    <- project(hull_proj_ea, "EPSG:4326")
  
  return(list(train = hull_buf, proj = hull_proj))
}

species_counts <- soil_df %>%
  count(Species) %>%
  filter(n > 60) %>%
  filter(Species != "", !grepl("^Marathrum$", Species)) %>%
  rename(n_occs = n)
print(species_counts)

contrib_list <- list()

for (sp in species_counts$Species) {
  
  message("\n--- Running ENM for: ", sp, " ---")
  sp_clean <- gsub(" ", "_", sp)
  
  occs_sp <- soil_df %>%
    filter(Species == sp) %>%
    select(Longitude, Latitude) %>%
    filter(complete.cases(.))
  
  if (nrow(unique(occs_sp)) < 3) {
    message("  Skipping — fewer than 3 unique points.")
    next
  }
  
  hulls_sp <- make_hull(occs_sp)
  
  env_train_sp <- crop(env_stack_selected, hulls_sp$train)
  env_train_sp <- mask(env_train_sp, hulls_sp$train)
  env_proj_sp  <- crop(env_stack_selected, hulls_sp$proj)
  
  if (ncell(env_train_sp) > 5e6) {
    message("  Aggregating training raster...")
    env_train_sp <- aggregate(env_train_sp, fact = 2, fun = "mean")
  }
  
  enmeval_sp <- tryCatch(
    ENMevaluate(
      occs       = occs_sp,
      envs       = env_train_sp,
      algorithm  = "maxnet",
      partitions = "block",
      tune.args  = list(fc = c("L","LQ","LQH","LQHP"), rm = c(1,2,4)),
      doClamp    = TRUE,
      parallel   = FALSE
    ),
    error = function(e) { message("  Failed: ", e$message); return(NULL) }
  )
  
  if (is.null(enmeval_sp)) next
  
  eval_sp   <- eval.results(enmeval_sp)
  models_sp <- eval.models(enmeval_sp)
  
  eval_sp$tune.args <- as.character(eval_sp$tune.args)
  eval_sp$AICc[is.na(eval_sp$AICc)] <- 999999
  
  write.csv(
    as.data.frame(eval_sp),
    file.path("Jun29_genus_and_species_ENMs", paste0(sp_clean, "_enmeval_results.csv")),
    row.names = FALSE
  )
  
  best_row   <- which.min(eval_sp$AICc)
  best_tune  <- eval_sp$tune.args[best_row]
  best_model <- models_sp[[best_tune]]
  
  # Extract variable contributions
  coefs     <- best_model$betas
  var_names <- gsub("\\^2|\\(hinge\\)|\\(threshold\\)|:.*", "", names(coefs))
  var_names <- gsub("I\\(|\\)", "", var_names)
  contrib   <- tapply(abs(coefs), var_names, sum)
  contrib   <- contrib / sum(contrib) * 100
  contrib_list[[sp]] <- data.frame(
    Species          = sp,
    Variable         = names(contrib),
    Contribution_pct = round(as.numeric(contrib), 2)
  )
  
  pred_sp <- terra::predict(
    env_proj_sp,
    best_model,
    fun = function(model, ...) predict(model, ..., type = "cloglog", clamp = TRUE),
    na.rm = TRUE
  )
  
  threshold_10p <- eval_sp$or.10p.avg[best_row]
  binary_sp     <- terra::ifel(pred_sp >= threshold_10p, 1, 0)
  
  png(
    file.path("Jun29_genus_and_species_ENMs", paste0(sp_clean, "_ENM_plots.png")),
    width = 1600, height = 700, res = 150
  )
  par(mfrow = c(1, 2))
  plot(pred_sp,
       main = paste0(sp, " Continuous\nFC=", eval_sp$fc[best_row],
                     " RM=", eval_sp$rm[best_row]))
  points(occs_sp, pch = 20, cex = 0.4, col = "red")
  plot(binary_sp,
       main = paste0(sp, " Binary (10%)\nFC=", eval_sp$fc[best_row],
                     " RM=", eval_sp$rm[best_row]))
  points(occs_sp, pch = 20, cex = 0.4, col = "red")
  dev.off()
  
  writeRaster(pred_sp,
              file.path("Jun29_genus_and_species_ENMs", paste0(sp_clean, "_continuous_best.tif")),
              overwrite = TRUE, datatype = "FLT4S")
  writeRaster(binary_sp,
              file.path("Jun29_genus_and_species_ENMs", paste0(sp_clean, "_binary10p_best.tif")),
              overwrite = TRUE, datatype = "INT1U")
  
  message("  Saved: ", sp_clean)
}

# Combine and save contributions
contrib_all <- do.call(rbind, contrib_list)
rownames(contrib_all) <- NULL
write.csv(contrib_all, "Jun29_genus_and_species_ENMs/variable_contributions.csv", row.names = FALSE)

library(tidyr)
contrib_wide <- contrib_all %>%
  pivot_wider(names_from = Variable, values_from = Contribution_pct, values_fill = 0)
write.csv(contrib_wide, "Jun29_genus_and_species_ENMs/variable_contributions_wide.csv", row.names = FALSE)

#Species-level statistical test: does best validation AUC correlate with occurrence count?
library(dplyr)

#Occurrence counts
species_n <- soil_df %>%
  count(Species, name = "n_occs")

#Load all species ENMeval CSVs and extract best model AUC
stats_dir <- "Jun29_genus_and_species_ENMs/"

species_files <- list.files(
  stats_dir,
  pattern = "_enmeval_results\\.csv$",
  full.names = TRUE
)

best_auc <- lapply(species_files, function(f) {
  
  df <- read.csv(f)
  
  species_name <- basename(f)
  species_name <- gsub("_enmeval_results\\.csv", "", species_name)
  species_name <- gsub("_", " ", species_name)
  
  best_row <- df[which.max(df$auc.val.avg), ]
  
  data.frame(
    Species    = species_name,
    best_auc   = best_row$auc.val.avg,
    best_fc    = best_row$fc,
    best_rm    = best_row$rm
  )
  
}) %>% bind_rows()

#Merge with occurrence counts
auc_df <- left_join(best_auc, species_n, by = "Species")
print(auc_df)

#Spearman correlation
spearman_test <- cor.test(
  auc_df$n_occs,
  auc_df$best_auc,
  method = "spearman"
)

print(spearman_test)

#Plot
library(ggplot2)
library(ggrepel)

ggplot(auc_df, aes(x = n_occs, y = best_auc, label = Species)) +
  geom_smooth(method = "lm",
              se = TRUE,
              color = "grey60",
              linetype = "dashed",
              linewidth = 0.7) +
  geom_point(size = 3, color = "#2c7bb6") +
  geom_text_repel(size = 2.5,
                  max.overlaps = Inf,
                  fontface = "italic") +
  annotate("text",
           x = max(auc_df$n_occs) * 0.7,
           y = min(auc_df$best_auc) + 0.01,
           label = paste0(
             "Spearman ρ = ",
             round(spearman_test$estimate, 3),
             "\np = ",
             signif(spearman_test$p.value, 3)
           ),
           size = 3.5,
           hjust = 0) +
  labs(
    x = "Number of Occurrences",
    y = "Best Model Validation AUC",
    title = "Effect of Occurrence Count on Species-Level Model Performance"
  )
