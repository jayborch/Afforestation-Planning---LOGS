# Load libraries
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(randomForest)
library(mgcv)

# -----------------------------
# read in points
# -----------------------------
sdm_points <- readRDS("models/occurences.RDS")

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
# sampling raster stack
# -----------------------------

pa_climate <- terra::extract(env_stack, sdm_points, xy = TRUE)
pa_climate$presence <- sdm_points$presence
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
gam_pca <- mgcv::gam(
  presence ~ s(PC1, k=12) + s(PC2, k=12) + s(PC3, k=12) + s(PC4, k=12) + s(PC5, k=12) + s(PC6, k=12) + s(PC7, k=12) + s(PC8, k=12),
  family = binomial(link = "logit"),
  data = gam_data,
  select = TRUE
)

gam.check(gam_pca)

# Optional: plot smooths
plot(gam_pca, pages=1, shade=TRUE, rug=TRUE)


# -----------------------------
# 1. Extract raster values
# -----------------------------
# env_stack = your current environmental raster stack
raster_values <- terra::values(env_stack, mat = TRUE)

# Keep only complete cases
complete_rows <- complete.cases(raster_values)
raster_values_complete <- raster_values[complete_rows, ]

# -----------------------------
# 2. Apply PCA rotation
# -----------------------------
# pca = your PCA object used for training
# Use same PCs as in GAM (PC1–PC5)
raster_pcs <- predict(pca, newdata = raster_values_complete)[, 1:8]
colnames(raster_pcs) <- paste0("PC", 1:8)

# Convert to data.frame for GAM
raster_pcs_df <- data.frame(raster_pcs)

# -----------------------------
# 3. Predict with GAM
# -----------------------------
# gam_pca_updated = your trained GAM
predicted_prob <- predict(gam_pca, newdata = raster_pcs_df, type = "response")

# -----------------------------
# 4. Put predictions back into raster
# -----------------------------
pred_raster <- rast(env_stack[[1]])  # use first layer as template
values(pred_raster)[complete_rows] <- predicted_prob

# Optional: mask with original raster to keep NA areas
pred_raster <- mask(pred_raster, env_stack[[1]])

# -----------------------------
# Load projected climate data
# -----------------------------
files <- list.files("data/climate_change_projections", pattern = "\\.tif$", full.names = TRUE)
projected_climate_stack <- rast(files)
projected_crop <- crop(projected_climate_stack, iceland_shape)
projected_crop <- mask(projected_crop, iceland_shape)

# Number of bioclim variables per period
n_vars <- 19
period_names <- c("ACCESS_CM2_ssp370_2021_2040",
                  "ACCESS_CM2_ssp370_2041_2060",
                  "ACCESS_CM2_ssp370_2061_2080",
                  "ACCESS_CM2_ssp370_2081_2100")

# -----------------------------
# Prepare period-specific stacks
# -----------------------------
period_stacks <- list()
for (i in 0:3) {
  start <- i * n_vars + 1
  end <- start + n_vars - 1
  period_stack <- projected_crop[[start:end]]
  
  # Rename layers to match baseline model
  names(period_stack) <- paste0("wc2.1_30s_bio_", 1:19)
  
  # Add STATIC layers (slope, tpi, soil, DEM)
  period_env <- c(period_stack, slope_regional, tpi, soil_stack, dem_regional_crop)
  
  period_stacks[[period_names[i+1]]] <- period_env
}

# -----------------------------
# Predict birch probability using GAM + PCA
# -----------------------------
predictions_list <- list()
for (nm in names(period_stacks)) {
  cat("Predicting:", nm, "\n")
  
  # Extract raster values
  raster_values <- terra::values(period_stacks[[nm]], mat=TRUE)
  complete_rows <- complete.cases(raster_values)
  raster_values_complete <- raster_values[complete_rows, ]
  
  # Apply PCA rotation (same as training)
  raster_pcs <- predict(pca, newdata=raster_values_complete)[, 1:n_pcs]
  colnames(raster_pcs) <- paste0("PC", 1:n_pcs)
  raster_pcs_df <- data.frame(raster_pcs)
  
  # Predict probability
  predicted_prob <- predict(gam_pca, newdata=raster_pcs_df, type="response")
  
  # Fill raster
  pred_raster <- rast(period_stacks[[nm]][[1]])  # template
  values(pred_raster)[complete_rows] <- predicted_prob
  pred_raster <- mask(pred_raster, period_stacks[[nm]][[1]])
  
  predictions_list[[nm]] <- pred_raster
}

# -----------------------------
# Stack predictions for plotting
# -----------------------------
stack_for_plot <- rast(predictions_list)
names(stack_for_plot) <- names(predictions_list)

# -----------------------------
# Convert to long dataframe for ggplot
# -----------------------------
df_list <- lapply(names(stack_for_plot), function(tp) {
  r <- stack_for_plot[[tp]]
  df <- as.data.frame(r, xy=TRUE)
  names(df)[3] <- "probability"
  df$time <- tp
  df
})
df_long <- bind_rows(df_list)

# -----------------------------
# Plot with facets
# -----------------------------
ggplot(df_long, aes(x=x, y=y, fill=probability)) +
  geom_raster() +
  scale_fill_viridis_c(option="viridis", na.value="white", limits=c(0,1)) +
  coord_equal() +
  facet_wrap(~time) +
  labs(fill="Birch Probability", x="Longitude", y="Latitude") +
  theme_minimal()





