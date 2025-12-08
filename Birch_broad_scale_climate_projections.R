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
# 1. Predict baseline birch probability
# -----------------------------
baseline_values <- terra::values(env_stack, mat = TRUE)
complete_rows <- complete.cases(baseline_values)
baseline_values_complete <- baseline_values[complete_rows, ]

# Apply PCA rotation (same as training)
baseline_pcs <- predict(pca, newdata = baseline_values_complete)[, 1:n_pcs]
colnames(baseline_pcs) <- paste0("PC", 1:n_pcs)
baseline_pcs_df <- data.frame(baseline_pcs)

# Predict probability
baseline_pred <- predict(gam_pca, newdata = baseline_pcs_df, type = "response")

# Fill raster
baseline_rast <- rast(env_stack[[1]])  # template
values(baseline_rast)[complete_rows] <- baseline_pred
baseline_rast <- mask(baseline_rast, env_stack[[1]])

# -----------------------------
# 2. Predict future periods
# -----------------------------
predictions_list <- list()
for (nm in names(period_stacks)) {
  cat("Predicting:", nm, "\n")
  
  # Extract raster values
  raster_values <- terra::values(period_stacks[[nm]], mat = TRUE)
  complete_rows <- complete.cases(raster_values)
  raster_values_complete <- raster_values[complete_rows, ]
  
  # Apply PCA rotation (same as training)
  raster_pcs <- predict(pca, newdata = raster_values_complete)[, 1:n_pcs]
  colnames(raster_pcs) <- paste0("PC", 1:n_pcs)
  raster_pcs_df <- data.frame(raster_pcs)
  
  # Predict probability
  predicted_prob <- predict(gam_pca, newdata=raster_pcs_df, type="response")
  
  # Fill raster
  pred_raster <- rast(period_stacks[[nm]][[1]])  # template
  values(pred_raster)[complete_rows] <- predicted_prob
  pred_raster <- mask(pred_raster, period_stacks[[nm]][[1]])
  
  # Store in list
  predictions_list[[nm]] <- pred_raster
}

# -----------------------------
# 3. Combine baseline + future periods into one stack
# -----------------------------
stack_with_baseline <- c(baseline_rast, rast(predictions_list))
names(stack_with_baseline) <- c("Baseline", names(predictions_list))

# Write combined stack to disk
writeRaster(stack_with_baseline, "data/climate_change_projections/predictions/ssp370_ACCESS_CM2_stack.tif", overwrite = T)


# -----------------------------
# 2. Convert to long dataframe for ggplot
# -----------------------------
df_list <- lapply(names(stack_with_baseline), function(tp) {
  r <- stack_with_baseline[[tp]]
  df <- as.data.frame(r, xy=TRUE)
  names(df)[3] <- "probability"
  df$time <- tp
  df
})
df_long <- bind_rows(df_list)

# -----------------------------
# 3. Plot faceted probability map including baseline
# -----------------------------
ggplot(df_long, aes(x=x, y=y, fill=probability)) +
  geom_raster() +
  scale_fill_viridis_c(option="viridis", na.value="white", limits=c(0,1)) +
  coord_equal() +
  facet_wrap(~time) +
  labs(fill="Birch Probability", x="Longitude", y="Latitude") +
  theme_minimal()

# -----------------------------
# 4. Prepare categories for Sankey (baseline included)
# -----------------------------
cat_fun <- function(x) {
  cut(x, breaks=c(-Inf, 0.5, 0.75, Inf),
      labels=c("Unsuitable", "Suitable", "Highly suitable"))
}

cat_stack <- stack_with_baseline
for(i in 1:nlyr(cat_stack)) {
  values(cat_stack[[i]]) <- cat_fun(values(cat_stack[[i]]))
}

df <- as.data.frame(cat_stack, xy=TRUE)

# Optional: sample pixels
set.seed(42)
df <- df[sample(1:nrow(df), 5000), ]

# -----------------------------
# 5. Prepare transitions
# -----------------------------
cat_cols <- names(cat_stack)
df_transitions <- df %>% select(all_of(cat_cols))

edges_list <- list()
for(i in 1:(length(cat_cols)-1)) {
  
  tmp_df <- df_transitions[, c(i, i+1)]
  colnames(tmp_df) <- c("from", "to")
  
  # Convert factors to character
  tmp_df <- tmp_df %>% mutate(from = as.character(from),
                              to   = as.character(to))
  
  # Prepend period info
  tmp_df <- tmp_df %>% mutate(from = paste0(cat_cols[i], "_", from),
                              to   = paste0(cat_cols[i+1], "_", to))
  
  # Count transitions
  tmp_count <- as.data.frame(
    aggregate(rep(1, nrow(tmp_df)),
              by = list(from = tmp_df$from, to = tmp_df$to),
              FUN = sum)
  )
  colnames(tmp_count)[3] <- "value"
  
  edges_list[[i]] <- tmp_count
}

edges_df <- do.call(rbind, edges_list)

# -----------------------------
# 6. Create nodes and links for networkD3
# -----------------------------
nodes <- data.frame(name = unique(c(edges_df$from, edges_df$to)))
edges_df$source <- match(edges_df$from, nodes$name) - 1
edges_df$target <- match(edges_df$to, nodes$name) - 1
links <- edges_df %>% select(source, target, value)

# -----------------------------
# 7. Plot interactive Sankey including baseline
# -----------------------------
sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 30
)

# -----------------------------
# Parameters (tweak as needed)
# -----------------------------
suitable_thresh <- 0.5   # threshold for "suitable" used in persistence/survival
weights <- c(persist = 0.6, meanp = 0.3, survival = 0.1)  # must sum to 1
stopifnot(abs(sum(weights) - 1) < 1e-6)

# Optional transform: "none", "power" (gamma>1 accentuates highs), or "logit_like"
transform_type <- "power"
power_gamma <- 1.4        # >1 emphasizes high scores, <1 compresses
# -----------------------------

n_periods <- nlyr(stack_with_baseline)
if(n_periods < 2) stop("stack_with_baseline must include baseline + >=1 future period")

# 1) Binary suitability matrix (0/1) across periods
bin_stack <- app(stack_with_baseline, fun = function(x) as.integer(x >= suitable_thresh))

# 2) Persistence proportion (0..1)
persist_count <- app(bin_stack, fun = sum, na.rm = TRUE)      # 0..n_periods
persist_prop  <- persist_count / n_periods                    # 0..1

# 3) Mean probability across periods (0..1)
mean_prob <- app(stack_with_baseline, fun = mean, na.rm = TRUE)

# 4) Survival: number of future periods the cell remains suitable after baseline
# Compute vectorized for speed
bin_mat <- as.matrix(values(bin_stack))   # rows=cells, cols=periods
nr <- nrow(bin_mat); nc <- ncol(bin_mat)
surv_vec <- rep(NA_real_, nr)

for(i in seq_len(nr)) {
  v <- bin_mat[i, ]
  if(all(is.na(v))) { surv_vec[i] <- NA; next }
  if(v[1] == 0) { surv_vec[i] <- 0; next }        # baseline not suitable => survival 0
  future_v <- v[-1]
  if(all(is.na(future_v))) { surv_vec[i] <- NA; next }
  # If never lost -> survive all future periods
  if(all(future_v == 1, na.rm = TRUE)) {
    surv_vec[i] <- (n_periods - 1)
  } else {
    zero_idx <- which(future_v == 0)
    if(length(zero_idx) == 0) surv_vec[i] <- (n_periods - 1)
    else surv_vec[i] <- (zero_idx[1] - 1)   # number of future periods survived before loss
  }
}

# Put survival back into raster and normalize 0..1
surv_rast <- rast(stack_with_baseline[[1]])
vals_surv <- rep(NA_real_, ncell(surv_rast))
valid_idx <- which(!is.na(bin_mat[,1]))  # cells with data
vals_surv[valid_idx] <- surv_vec[valid_idx]
values(surv_rast) <- vals_surv

max_surv <- (n_periods - 1)
norm_surv <- surv_rast / max_surv
norm_surv <- clamp(norm_surv, lower=0, upper=1)

# 5) Combine into weighted composite score (continuous)
score_raw <- (persist_prop * weights["persist"]) +
  (mean_prob    * weights["meanp"]) +
  (norm_surv    * weights["survival"])

# 6) Normalize to 0..1 robustly (handle NA)
minv <- global(score_raw, "min", na.rm = TRUE)[1,1]
maxv <- global(score_raw, "max", na.rm = TRUE)[1,1]
if (is.na(minv) | is.na(maxv) | maxv == minv) {
  score_norm <- score_raw  # fallback (unlikely)
} else {
  score_norm <- (score_raw - minv) / (maxv - minv)
}

# 7) Optional monotonic transform to emphasize highs (useful in MCDA)
if(transform_type == "power") {
  score_trans <- app(score_norm, fun = function(x) x^power_gamma)
} else if(transform_type == "logit_like") {
  # tiny logit-like: log(x/(1-x)) scaled back to 0..1 (avoid extremes)
  eps <- 1e-6
  score_trans <- app(score_norm, fun = function(x) {
    x <- pmin(pmax(x, eps), 1-eps)
    z <- log(x/(1-x))
    # rescale z to 0..1 using theoretical min/max for x in (eps,1-eps)
    z_min <- log(eps/(1-eps)); z_max <- log((1-eps)/eps)
    (z - z_min) / (z_max - z_min)
  })
} else {
  score_trans <- score_norm
}

# 9) Quick diagnostics: histogram + map
par(mfrow = c(1,2))
hist(values(score_trans), main = "Persistence score (continuous)", xlab = "score", breaks = 50)
plot(score_trans, main = "Persistence score (continuous)")

terra::writeRaster(score_trans, "data/persistence_score.tif")

