# Load libraries
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mgcv)
library(networkD3)
library(pROC)
library(blockCV)
library(geodata)
library(brms)
library(terra)
library(stringr)
library(progress)


# -----------------------------
# Load baseline bioclim data
# -----------------------------
data_folder <- "data/wc2.1_30s_bio/"
files <- list.files(data_folder, pattern = "^wc2\\.1_30s_bio_.*\\.tif$", full.names = TRUE)
bio_stack <- rast(files)

# -----------------------------
# Project to match climate raster CRS
iceland_shape <- vect("data/Iceland_shapefile/is_1km.shp")
iceland_shape <- project(iceland_shape, crs(bio_stack))

# Crop and mask to Iceland
bio_iceland <- crop(bio_stack, iceland_shape)
bio_regional_crop <- mask(bio_iceland, iceland_shape)

# -----------------------------
# Load DEM and calculate slope & TPI
# -----------------------------
dem_regional <- rast('data/wc2.1_30s_elev/wc2.1_30s_elev.tif')
dem_regional_crop <- crop(dem_regional, iceland_shape)
dem_regional_crop <- mask(dem_regional_crop, iceland_shape)
names(dem_regional_crop) <- "elev"

# Slope
slope_regional <- terrain(dem_regional_crop, v = "slope", unit = "degrees")
names(slope_regional) <- "slope"

# Topographic Position Index (5x5 moving window)
w <- matrix(1, 5, 5)
local_mean <- focal(dem_regional_crop, w = w, fun = mean, na.rm = TRUE)
tpi <- dem_regional_crop - local_mean
names(tpi) <- "tpi"

# -----------------------------
# Load soil data
# -----------------------------
soil_stack <- rast("data/soils/iceland_soilgrid_0-5cm.tif")
soil_stack <- project(soil_stack, bio_stack)
soil_stack <- crop(soil_stack, iceland_shape)
soil_stack <- mask(soil_stack, iceland_shape)

# -----------------------------
# Combine baseline layers into env_stack + include elevation
# -----------------------------
env_stack <- c(bio_regional_crop, slope_regional, tpi, soil_stack, dem_regional_crop)

# -----------------------------
# read in points and sampling raster stack
# -----------------------------

sdm_points <- readRDS("models/occurences.RDS")

# 1. Project points to match the environmental stack
sdm_points_proj <- project(sdm_points, crs(env_stack))

# 2. Get the cell numbers for each point
# This returns a matrix where the second column contains the cell index
pt_cells <- cells(env_stack[[1]], sdm_points_proj)

# 3. Identify and keep only one point per unique cell
# We use the 'cell' column from the pt_cells matrix
unique_indices <- which(!duplicated(pt_cells[, "cell"]))
sdm_thinned_v <- sdm_points_proj[unique_indices, ]

# 4. Remove points that fall on NA values (important for coastal Iceland)
# Extract values from the stack at these points
extracted_vals <- terra::extract(env_stack[[1]], sdm_thinned_v)
sdm_thinned_v <- sdm_thinned_v[!is.na(extracted_vals[, 2]), ]

pa_climate <- terra::extract(env_stack, sdm_thinned_v, xy = TRUE)
pa_climate$presence <- sdm_thinned_v$presence
pa_climate <- na.omit(pa_climate)
pa_climate$presence <- factor(pa_climate$presence, levels = c(0,1))
pa_climate <- pa_climate %>% dplyr::select(-c(ID,x,y))

# -----------------------------
# Prepare predictor matrix
# -----------------------------
predictors <- pa_climate %>% dplyr::select(-presence)

# -----------------------------
# Run PCA
# -----------------------------
# Scale = TRUE standardizes variables (important since units differ)
pca <- prcomp(predictors, scale. = TRUE)

# Inspect variance explained
summary(pca)

# Choose top PCs explaining ~80% of variance
# Example: let's take the first 5 PCs
n_pcs <- 8
pc_scores <- pca$x[, 1:n_pcs]

# Combine with response
gam_data <- data.frame(presence = pa_climate$presence, pc_scores)

# Rename PC columns for convenience
colnames(gam_data)[-1] <- paste0("PC", 1:n_pcs)

# -----------------------------
# Fit GAM with PCs
# -----------------------------
#gam_pca <- mgcv::gam(
#  presence ~ s(PC1, k=20) + s(PC2, k=20) + s(PC3, k=20) + s(PC4, k=20) + s(PC5, k=20) + s(PC6, k=20) + s(PC7, k=20) + s(PC8, k=20),
#  family = binomial(link = "logit"),
# data = gam_data,
#  select = T
#)

#gam.check(gam_pca)

# Optional: plot smooths
#plot(gam_pca, pages=1, shade=TRUE, rug=TRUE)

# -----------------------------
# Fit brms with PCs
# -----------------------------
library(parallel)

# 1. Detect cores and reserve one
n_cores <- parallel::detectCores() - 1

# 2. Fit the model
birch_brms_pca <- brm(
  formula = presence ~ s(PC1, k=20) + s(PC2, k=20) + s(PC3, k=20) + 
    s(PC4, k=20) + s(PC5, k=20) + s(PC6, k=20) + 
    s(PC7, k=20) + s(PC8, k=20),
  data = gam_data,
  family = bernoulli(link = "logit"),
  chains = 4, 
  cores = n_cores,      # Uses all but one core
  iter = 2000, 
  warmup = 1000,
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr"  # Recommended for speed
)

plot(birch_brms_pca)
pp_check(birch_brms_pca, type = "stat", stat = "mean")

# ---------------------------------------------------------
# BASELINE PROJECTION: SAVING 50 POSTERIOR DRAWS
# ---------------------------------------------------------

# 1. Prepare the raster data as a data frame for PCA projection
# Ensure env_df is created within this block to avoid "object not found" errors
# xy = TRUE keeps coordinates, na.rm = TRUE ensures we only predict for Iceland landmass

message("Preparing raster data for projection...")
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)

# 2. Project the raster data into the same PCA space as the model
# Using the pre-trained 'pca' object to transform our Iceland grid
pca_transformed <- predict(pca, newdata = env_df[, -c(1,2)]) # exclude x,y columns

# Keep only the PCs used in the model
pca_df <- as.data.frame(pca_transformed[, 1:n_pcs])
colnames(pca_df) <- paste0("PC", 1:n_pcs)

# 3. Draw 50 samples from the posterior distribution
# posterior_epred gives us the expected probability (0-1) for each cell
message("Generating 50 posterior draws from brms model...")
epred_mat <- posterior_epred(
  birch_brms_pca, 
  newdata = pca_df, 
  ndraws = 50
)

# 4. Reconstruct the 50-layer Prediction Stack
# Optimization: Populate a single multi-layer stack instead of a list of rasters
message("Reconstructing spatial layers...")
raster_template <- env_stack[[1]]
prediction_stack <- rast(raster_template, nlyr = 50)

# Get the cell indices for all our valid data points (non-NA)
# This maps the rows of epred_mat back to their correct physical location in Iceland
cell_indices <- cellFromXY(raster_template, as.matrix(env_df[, c("x", "y")]))

# Assign each draw to a layer
for(i in 1:50) {
  # epred_mat is [draws, observations]
  prediction_stack[[i]][cell_indices] <- epred_mat[i, ]
}

names(prediction_stack) <- paste0("draw_", 1:50)

# 5. Save the output
output_dir <- "models/future_projections_brms_all_draws"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- file.path(output_dir, "birch_posterior_draws_baseline.tif")

writeRaster(
  prediction_stack, 
  filename = output_file,
  overwrite = TRUE
)

message("Finished! Saved 50 layers to: ", output_file)




# ==============================================================================
# SPATIAL CROSS-VALIDATION & PERFORMANCE BENCHMARKING
# ==============================================================================
# This script executes a 5-Fold Spatial Block Cross-Validation (CV) to 
# rigorously test the model's predictive power and transferability.
#
# KEY CONCEPTS:
# 1. SPATIAL DECOUPLING (blockCV):
#    Traditional random CV often overestimates model performance due to spatial 
#    autocorrelation. By using 50km spatial blocks, we ensure that training 
#    and testing data are geographically separated, forcing the model to 
#    generalize across different regions of Iceland.
#
# 2. DIMENSIONAL INTEGRITY:
#    The validation loop applies the global PCA rotation to each spatial fold. 
#    This ensures that the environmental "latent space" remains consistent 
#    regardless of which geographic region is being held out.
#
# 3. DUAL-METRIC EVALUATION:
#    - AUC (Discrimination): Measures how well the model distinguishes 
#      between presence and absence (Relative ranking).
#    - RMSE (Calibration): Measures the magnitude of error between predicted 
#      probabilities and actual outcomes (Absolute accuracy).
#
# 4. CALIBRATION DIAGNOSTICS:
#    Aggregated out-of-sample predictions are binned into deciles to produce 
#    a Calibration Plot. This visualizes whether a "70% probability" actually 
#    corresponds to a 70% observed frequency of birch in the field.
# ==============================================================================

# --- 1. Prepare Data for Spatial Blocking (unchanged) ---
pa_data <- sdm_thinned_v 
pa_sf <- st_as_sf(pa_data, coords = c("x", "y"), crs = crs(env_stack, proj=TRUE))
pa_sf <- st_transform(pa_sf, crs(env_stack))

# --- 2. Spatial Blocking (unchanged) ---
set.seed(42)
sb <- blockCV::spatialBlock(
  speciesData = pa_sf,
  species = "presence",
  rasterLayer = env_stack[[1]],
  theRange = 50000,
  k = 5,
  selection = "random",
  iteration = 100,
  showBlocks = FALSE
)
pa_sf$foldID <- sb$foldID

# --- 3. Initialize CV Storage and Define Predictor Order ---
k_folds <- 5
auc_results <- numeric(k_folds)
rmse_results <- numeric(k_folds)

# ***NEW: Initialize dataframe to store all out-of-sample predictions for the plot***
cv_predictions_all <- data.frame(
  actual_outcome = numeric(),
  predicted_prob = numeric()
)

# Get the exact names of the ORIGINAL predictors used to train the PCA.
pca_predictor_names <- names(pca$center) 
n_pcs <- 8 # Define this once

# --- 4. The 5-Fold Spatial Cross-Validation Loop (FIXED DATA PREPARATION) ---
cat("\n============================================\n")
cat("Starting 5-Fold Spatial Cross-Validation (AUC, RMSE, & Plotting)...\n")
cat("============================================\n")

for (i in 1:k_folds) {
  cat(paste0("  Processing Fold ", i, "...\n"))
  
  # a. Split points based on spatial blocks (training and testing)
  train_sf <- pa_sf[pa_sf$foldID != i, ]
  test_sf  <- pa_sf[pa_sf$foldID == i, ]
  
  # b. Extract environmental data
  train_vect <- terra::vect(train_sf)
  test_vect  <- terra::vect(test_sf)
  train_env <- as.data.frame(terra::extract(env_stack, train_vect, bind = TRUE))
  test_env  <- as.data.frame(terra::extract(env_stack, test_vect, bind = TRUE))
  
  
  # === FIX START: Data Cleaning and PCA Preparation ===
  # Remove ID/foldID columns and NA rows in two clear steps.
  
  # Data Cleaning for Training Set
  train_env_clean <- train_env %>% 
    dplyr::select(-any_of(c("ID", "foldID"))) %>% 
    na.omit() 
  
  # Data Cleaning for Testing Set
  test_env_clean <- test_env %>% 
    dplyr::select(-any_of(c("ID", "foldID"))) %>% 
    na.omit()
  
  # Select the original predictor columns for the PCA transformation
  train_predictors <- train_env_clean %>% dplyr::select(all_of(pca_predictor_names))
  test_predictors  <- test_env_clean %>% dplyr::select(all_of(pca_predictor_names))
  
  # Apply PCA (using the *global* PCA object 'pca' trained on all data)
  train_pcs <- predict(pca, newdata = train_predictors)[, 1:n_pcs]
  test_pcs  <- predict(pca, newdata = test_predictors)[, 1:n_pcs]
  
  # Prepare GAM data frames
  train_gam_data <- data.frame(presence = factor(train_env_clean$presence, levels = c(0,1)),
                               train_pcs)
  test_gam_data  <- data.frame(presence = factor(test_env_clean$presence, levels = c(0,1)),
                               test_pcs)
  # Rename PC columns for GAM formula consistency
  colnames(train_gam_data)[-1] <- paste0("PC", 1:n_pcs)
  colnames(test_gam_data)[-1] <- paste0("PC", 1:n_pcs)
  
  # === FIX END ===
  
  # Fit GAM using only the training fold data
  gam_cv <- mgcv::gam(
    presence ~ s(PC1, k=20) + s(PC2, k=20) + s(PC3, k=20) + s(PC4, k=20) + s(PC5, k=20) + s(PC6, k=20) + s(PC7, k=20) + s(PC8, k=20),
    family = binomial(link = "logit"),
    data = train_gam_data,
    select = TRUE
  )
  
  # Predict on the hold-out test fold
  predictions <- predict(gam_cv, newdata = test_gam_data, type = "response")
  actual_outcome <- as.numeric(as.character(test_gam_data$presence))
  
  # --- e. Evaluate Performance (AUC and RMSE) ---
  
  # AUC
  roc_curve <- pROC::roc(response = test_gam_data$presence, predictor = predictions, quiet = TRUE)
  auc_results[i] <- pROC::auc(roc_curve)
  
  # RMSE
  rmse_results[i] <- sqrt(mean((actual_outcome - predictions)^2, na.rm = TRUE))
  
  # ***NEW: Store predictions for calibration plot***
  fold_results <- data.frame(
    actual_outcome = actual_outcome,
    predicted_prob = predictions
  )
  cv_predictions_all <- rbind(cv_predictions_all, fold_results)
}

# --- 5. Aggregate and Report Results (unchanged) ---
avg_auc <- mean(auc_results, na.rm = TRUE)
sd_auc <- sd(auc_results, na.rm = TRUE)
avg_rmse <- mean(rmse_results, na.rm = TRUE)
sd_rmse <- sd(rmse_results, na.rm = TRUE)

cat("\n============================================\n")
cat("    FINAL SPATIAL CROSS-VALIDATION RESULTS\n")
cat("============================================\n")
cat("Method: 5-Fold Spatial Block CV (50km blocks)\n")
cat("Model: GAM with 8 Principal Components\n")
cat(paste0("Average AUC (Discrimination): ", round(avg_auc, 3), " (SD: ", round(sd_auc, 3), ")\n"))
cat(paste0("Average RMSE (Error Magnitude): ", round(avg_rmse, 3), " (SD: ", round(sd_rmse, 3), ")\n"))
cat("============================================\n")


# --- 6. Calibration Plot Generation (unchanged) ---

# Aggregate data into 10 bins (deciles) for plotting
calibration_data <- cv_predictions_all %>%
  mutate(pred_bin = cut(predicted_prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(pred_bin) %>%
  summarise(
    mean_predicted = mean(predicted_prob, na.rm = TRUE),
    mean_observed = mean(actual_outcome, na.rm = TRUE),
    n_points = n()
  ) %>%
  ungroup()


# Create the plot
calib_plot <- ggplot(calibration_data, aes(x = mean_predicted, y = mean_observed)) +
  # Add the ideal 45-degree line for perfect calibration
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  # Add the observed calibration points (weighted by the number of points in the bin)
  geom_point(aes(size = n_points), color = "darkblue") +
  # Add a smoothed line (often a LOESS line) to show the trend
  geom_smooth(method = "loess", se = FALSE, color = "red", linetype = "solid") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Spatial Cross-Validation Calibration Plot (Aggregated)",
    x = "Mean Predicted Probability (in bin)",
    y = "Mean Observed Frequency (in bin)",
    size = "N points"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(calib_plot)

# --- 7. Final Steps for Reporting (unchanged) ---
# Also report the overall goodness-of-fit for context
cat(paste0("\nGlobal Deviance Explained by final model: ", round(summary(gam_cv)$dev.expl * 100, 1), "%\n"))
cat("\nModel Validation Complete. Next step: Future Climate Projection.\n")


# ==============================================================================
# FUTURE PROJECTION: SAVING ALL POSTERIOR DRAWS
# ==============================================================================
# This script automates the transformation of raw CMIP6 climate projections 
# into 50-layer Bayesian posterior probability stacks.
#
# KEY OPERATIONS:
# 1. LATENT SPACE PROJECTION:
#    Future climate tiles are resampled to the 1km Iceland grid and projected 
#    into the same PCA-transformed environmental space used during model training.
#
# 2. BAYESIAN UNCERTAINTY PRESERVATION:
#    Instead of predicting a single "mean," we extract 50 discrete posterior 
#    draws (epred) for every pixel. This preserves the full distribution of 
#    the brms model's internal uncertainty across the 160 scenarios.
#
# 3. SPATIOTEMPORAL HARMONIZATION:
#    Static physiographic variables (Slope, TPI, Soil) are combined with 
#    dynamic climate layers to ensure "site-specific" refugia characteristics 
#    are captured (e.g., topographic sheltering in Icelandic fjords).
#
# 4. HIGH-THROUGHPUT COMPRESSION:
#    Outputs are saved as multi-layer Cloud Optimized GeoTIFFs (COG-style) 
#    using DEFLATE compression to manage the massive data volume generated 
#    by the 50-draw Bayesian approach.
# ==============================================================================



n_draws <- 50
# Pre-sampling ensures the same 50 "versions" of the model are used for all 160 scenarios
draw_ids <- sample(seq_len(posterior::ndraws(birch_brms_pca)), n_draws)

static_layers_to_add <- c(slope_regional, tpi, soil_stack, dem_regional_crop)

data_path   <- "data/climate_change_projections"
output_path <- "models/future_projections_brms_all_draws"
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

future_files <- list.files(data_path, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

for (f in future_files) {
  cat("--- Processing:", basename(f), "---\n")
  
  try({
    # A. Load future climate layers
    fut_tile <- rast(f)
    
    # B. Match geometry (Resampling to ensure perfect alignment)
    fut_tile_resampled <- resample(fut_tile, slope_regional, method = "bilinear")
    
    # C. Rename climate layers to match training
    names(fut_tile_resampled) <- paste0("wc2.1_30s_bio_", 1:19)
    
    # D. Add static layers
    fut_env_stack <- c(fut_tile_resampled, static_layers_to_add)
    
    # E. Extract values
    fut_values <- terra::values(fut_env_stack, mat = TRUE)
    ok <- complete.cases(fut_values)
    
    if (sum(ok) > 0) {
      
      # ---- PCA TRANSFORMATION ----
      fut_pcs <- predict(
        pca,
        newdata = as.data.frame(fut_values[ok, ])
      )[, 1:n_pcs]
      
      fut_pcs <- as.data.frame(fut_pcs)
      colnames(fut_pcs) <- paste0("PC", 1:n_pcs)
      
      # ---- POSTERIOR PREDICTION (Getting the raw draws) ----
      # epred will be a matrix of [draws (50) x pixels]
      epred <- posterior_epred(
        birch_brms_pca,
        newdata = fut_pcs,
        draw_ids = draw_ids,
        re_formula = NA
      )
      
      # ---- REBUILD MULTI-LAYER RASTER ----
      # We create a stack with one layer for every draw
      template <- rast(slope_regional)
      
      # Create an empty multi-layer SpatRaster
      # We use 'replicate' to create a stack of 50 empty layers
      draw_stack <- rast(lapply(1:n_draws, function(x) template))
      
      # Fill each layer with its corresponding draw
      # We transpose epred so it is [pixels x draws] to match the assignment
      values(draw_stack)[ok] <- t(epred)
      
      # Name layers by the draw ID
      names(draw_stack) <- paste0("draw_", draw_ids)
      
      # ---- SAVE ----
      base_name <- tools::file_path_sans_ext(basename(f))
      
      # Compression is essential here as we are saving 50 layers per scenario
      writeRaster(
        draw_stack,
        filename = file.path(output_path, paste0(base_name, "_50draws.tif")),
        overwrite = TRUE,
        gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
      )
      
      cat("✓ Saved 50-layer stack successfully\n")
      
    } else {
      cat("! Warning: No valid pixels after masking\n")
    }
    
    # Cleanup memory
    rm(fut_tile, fut_tile_resampled, fut_env_stack, fut_values, fut_pcs, epred, draw_stack)
    gc()
    tmpFiles(remove = TRUE)
    
  }, silent = FALSE)
}


# ==============================================================================
# MACROREFUGIA IDENTIFICATION: BAYESIAN MULTI-CRITERIA REALIZATIONS
# ==============================================================================
# This script processes 2-decade discrete time windows to identify 
# macrorefugia for Betula pubescens in Iceland.
#
# KEY OPERATIONS:
# 1. BASELINE INTEGRATION:
#    The baseline stack (birch_posterior_draws_baseline.tif) is used as the 
#    t=0 point for every trajectory to calculate survival and persistence.
#
# 2. DRAW-LEVEL INTEGRITY:
#    Processes Draw N of the baseline alongside Draw N of all future windows.
#
# 3. MCDA LOGIC:
#    Combines Persistence (proportion of time suitable), Mean Probability, 
#    and Survival (how long suitability lasts before the first loss).
# ==============================================================================

# -----------------------------
# CONTROL PANEL
# -----------------------------
tau_thresh      <- 0.5   # Binary suitability threshold
weights         <- c(persist = 0.6, meanp = 0.3, survival = 0.1) 
power_gamma     <- 1.4   # Gamma transformation to highlight core refugia
n_draws         <- 50    # Number of Bayesian posterior draws
viability_limit <- 0.6   

# Paths (Adjusted to your specific Windows directory)
input_dir  <- "models/future_projections_brms_all_draws"
output_dir <- "models/refugia_realizations_stacks"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. Separate Baseline from Future Files
all_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
baseline_path <- all_files[str_detect(all_files, "baseline")]
future_files  <- all_files[!str_detect(all_files, "baseline")]

# 2. Helper to extract GCM_SSP Scenario ID
# Regex looks for the GCM and SSP pattern (e.g., ACCESS-CM2_ssp126)
get_scenario_id <- function(x) {
  str_extract(basename(x), "(?<=bioc_).*?(?=_20\\d{2})")
}

scenarios <- unique(na.omit(get_scenario_id(future_files)))

# 3. Load Baseline Raster once (50 layers)
baseline_rast <- rast(baseline_path)

# -----------------------------
# MAIN PROCESSING LOOP
# -----------------------------
pb <- progress_bar$new(
  format = "  Refugia Realizations [:bar] :percent | ETA: :eta | :current/:total (:message)",
  total = length(scenarios), 
  clear = FALSE, 
  width = 100
)

for (scen in scenarios) {
  pb$tick(tokens = list(message = scen))
  
  # A. Filter and sort future time-steps for this scenario
  scen_files <- future_files[get_scenario_id(future_files) == scen]
  # Sort by year to ensure correct chronological order (2021, 2041, 2061, 2081)
  scen_files <- scen_files[order(str_extract(basename(scen_files), "20\\d{2}"))]
  
  # Load future stacks (each is 50 layers)
  future_list <- lapply(scen_files, rast)
  n_periods <- length(future_list) + 1 # +1 for the baseline
  
  refugia_stack_list <- vector("list", n_draws)
  
  # B. Process each of the 50 Draws
  for (d in 1:n_draws) {
    
    # Construct the chronological trajectory for this draw
    # [Baseline_Draw_D, 2021_Draw_D, 2041_Draw_D...]
    draw_trajectory <- rast(c(
      list(baseline_rast[[d]]), 
      lapply(future_list, function(x) x[[d]])
    ))
    
    # 1. Binary Suitability
    bin_stack <- draw_trajectory >= tau_thresh
    
    # 2. Persistence Proportion (0..1)
    persist_prop <- app(bin_stack, fun = sum, na.rm = TRUE) / n_periods
    
    # 3. Mean Probability (0..1)
    mean_prob <- app(draw_trajectory, fun = mean, na.rm = TRUE)
    
    # 4. Survival: Consecutive windows of suitability from baseline
    calc_surv_func <- function(v) {
      if(all(is.na(v))) return(NA)
      if(v[1] == 0) return(0) # Not suitable at start = 0 survival
      fut_v <- v[-1]
      loss_idx <- which(fut_v == 0)[1]
      if(is.na(loss_idx)) return(length(fut_v)) 
      return(loss_idx - 1)
    }
    
    surv_count <- app(bin_stack, fun = calc_surv_func)
    norm_surv  <- clamp(surv_count / (n_periods - 1), 0, 1)
    
    # 5. MCDA Composite
    score_raw <- (persist_prop * weights["persist"]) +
      (mean_prob    * weights["meanp"]) +
      (norm_surv    * weights["survival"])
    
    # 6. Global Normalization & Power Transform
    # Using global range across the draw to maintain relative variance
    m_val <- global(score_raw, "min", na.rm = TRUE)[1,1]
    M_val <- global(score_raw, "max", na.rm = TRUE)[1,1]
    
    if(!is.na(M_val) && M_val > m_val) {
      score_norm <- (score_raw - m_val) / (M_val - m_val)
    } else {
      score_norm <- score_raw
    }
    
    refugia_stack_list[[d]] <- score_norm^power_gamma
  }
  
  # C. Save 50-layer realization stack for the GCM/SSP
  final_stack <- rast(refugia_stack_list)
  names(final_stack) <- paste0("draw_", 1:n_draws)
  
  writeRaster(
    final_stack,
    filename = file.path(output_dir, paste0(scen, "_refugia_realizations.tif")),
    overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
  )
  
  # Clean memory for next scenario
  rm(future_list, refugia_stack_list, final_stack, draw_trajectory, bin_stack)
  gc()
}


# ==============================================================================
# GRAND UNIFICATION: CONSOLIDATING 2,000 REALIZATIONS
# ==============================================================================
# This script aggregates 40 scenario-stacks (50 draws each) into a single 
# Master Refugia Map for use in secondary MCDA.
#
# METRICS GENERATED:
# 1. ENSEMBLE MEAN: The "Expected" macrorefugia score.
# 2. RELIABILITY: 1 - (Standard Deviation / Mean), penalizing disagreement.
# 3. CONSENSUS: The percentage of models agreeing on high-quality refugia.
# ==============================================================================

library(terra)
library(progress)

# -----------------------------
# CONTROL PANEL
# -----------------------------
input_dir   <- "models/refugia_realizations_stacks"
output_path <- "models/grand_unifying_refugia_map.tif"
agreement_threshold <- 0.6  # Score above which a cell is "Strict Refugium"

# -----------------------------
# AGGREGATION ENGINE
# -----------------------------
realization_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# We use a 'Mean of Means' approach to avoid loading 2,000 layers into RAM at once
# We will track running sums to calculate the global Mean and SD
sum_rast   <- NULL
sum_sq_rast <- NULL
count      <- 0
agreement_count <- NULL

pb <- progress_bar$new(
  format = "  Unifying Scenarios [:bar] :percent | Total Layers: :current/:total",
  total = length(all_files), width = 100
)

for (f in realization_files) {
  # Load the 50-draw stack for this scenario
  r <- rast(f)
  
  # 1. Update Running Sum (for Mean)
  scen_sum <- app(r, fun = sum, na.rm = TRUE)
  if(is.null(sum_rast)) {
    sum_rast <- scen_sum
  } else {
    sum_rast <- sum_rast + scen_sum
  }
  
  # 2. Update Running Sum of Squares (for Standard Deviation)
  scen_sum_sq <- app(r^2, fun = sum, na.rm = TRUE)
  if(is.null(sum_sq_rast)) {
    sum_sq_rast <- scen_sum_sq
  } else {
    sum_sq_rast <- sum_sq_rast + scen_sum_sq
  }
  
  # 3. Update Agreement Count (How many draws > threshold)
  scen_agree <- app(r >= agreement_threshold, fun = sum, na.rm = TRUE)
  if(is.null(agreement_count)) {
    agreement_count <- scen_agree
  } else {
    agreement_count <- agreement_count + scen_agree
  }
  
  count <- count + nlyr(r)
  pb$tick()
}

# -----------------------------
# FINAL CALCULATIONS
# -----------------------------
# 1. Grand Mean
grand_mean <- sum_rast / count

# 2. Grand Standard Deviation (Formula: sqrt(E[X^2] - (E[X])^2))
grand_sd <- sqrt((sum_sq_rast / count) - (grand_mean^2))

# 3. Reliability Index (1 - CV)
# We clamp CV at 1.0 to avoid negative reliability
cv <- grand_sd / (grand_mean + 0.0001) # small epsilon to avoid div by zero
reliability <- clamp(1 - cv, 0, 1)

# 4. Consensus Percentage
consensus_pct <- agreement_count / count

# -----------------------------
# THE GRAND UNIFYING SCORE
# -----------------------------
# We weight the mean by its reliability. 
# A high mean score is only "grand" if the models agree on it.
grand_unifying_score <- grand_mean * reliability

# Export as a multi-band TIFF
master_stack <- c(grand_unifying_score, grand_mean, reliability, consensus_pct)
names(master_stack) <- c("Grand_Unifying_Score", "Ensemble_Mean", "Reliability", "Consensus_Pct")

writeRaster(
  master_stack, 
  filename = output_path, 
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
)

cat("\n--- Grand Unifying Map Exported to models/ ---")
