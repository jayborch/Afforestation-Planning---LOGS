library(terra)
library(brms)
library(dplyr)

dem <- rast("data/IslandsDEMv1.0_10x10m_isn2016_zmasl.tif")

# resample elevaion
dem_100m <- aggregate(dem, 
                      fact = 10, 
                      fun = "mean", 
                      na.rm = TRUE)
names(dem_100m) <- "elev"


# create slope
slope <- terrain(dem_100m, v = "slope", unit = "degrees")
names(slope) <- "slope"

# 3. TPI (Topographic Position Index)
# Using a 3x3 window (approx 300m context) to define 'local' position
tpi <- focal(dem_100m, w = 3, fun = function(x) x[5] - mean(x[-5]))
names(tpi) <- "tpi"

# read in topographic wetness index
twi <- rast("data/saga_output/TWI.tif")
twi <- project(twi, dem_100m)
names(twi) <- "twi"

# read in solar
solar <- rast("data/saga_output/Mean_GrowingSeason.tif")
solar <- project(solar, dem_100m)
names(solar) <- "solar"

# read in VARI from GEE
vari <- rast("data/GEE/Iceland_VARI_100m_2024_GrowingSeason.tif")
vari <- project(vari, dem_100m)


broad_scale_SDM <- rast("data/broad_scale_rf_prob.tif")

# Broadscale model (already loaded)
broad_prob <- broad_scale_SDM
broad_proj <- terra::project(broad_prob, "EPSG:8088")

# Fill small NA gaps using a local mean
broad_filled <- terra::focal(
  broad_proj,
  w = 3,
  fun = function(x, ...) {
    if (is.na(x[5])) mean(x, na.rm = TRUE) else x[5]
  },
  na.policy = "omit",
  fillvalue = NA
)

# Define a 5 km smoothing radius
radius_m <- 5000
w <- terra::focalMat(broad_filled, d = radius_m, type = "circle")

# Apply smoothing (mean filter)
broad_smooth <- terra::focal(broad_filled, w = w, fun = mean, na.rm = TRUE)
broad_smooth_aligned <- terra::project(broad_smooth, crs(broad_prob), method = "bilinear")
broad_smooth_aligned <- terra::resample(broad_smooth_aligned, broad_prob, method = "bilinear")
broad_edges_added <- terra::cover(broad_prob, broad_smooth_aligned)
broad_scale_SDM <- project(broad_edges_added, dem_100m)
names(broad_scale_SDM) <- "sdm"

# create slope
env_stack <- c(dem_100m, slope, tpi, twi, solar, vari,  broad_scale_SDM)

# read in points
sdm_points <- readRDS("models/occurences.RDS")

# 1. Extract raster values to points
# 'sdm_points' should be a SpatVector or have lon/lat columns
# If it's a data frame with coords, convert it: 
# pts_vec <- vect(sdm_points, geom=c("lon", "lat"), crs=crs(dem_100m))

# 1. Extract values
# This returns a data frame where the first column is 'ID'
extracted_values <- terra::extract(env_stack, sdm_points, xy = T)

# 2. Convert sdm_points to a standard data frame if it's a SpatVector
points_df <- as.data.frame(sdm_points)

# 3. Combine and Filter
analysis_df <- points_df %>%
  # Bind the extracted values (excluding the 'ID' column from extract)
  cbind(extracted_values[, -1]) %>% 
  # Now filter using across() or standard complete.cases on the df
  dplyr::filter(complete.cases(elev, slope, tpi, twi, solar, VARI, sdm))

# Double check the column names to ensure broad_scale_SDM 
# didn't change name during the stack/extract process
colnames(analysis_df)



# Transform probability to logit scale for the offset
analysis_df$logit_rf <- car::logit(analysis_df$sdm, adjust = 0.01)

# 1. Ensure the formula is correct
birch_forestry_formula <- bf(
  presence ~ 
    logit_rf + 
    solar + 
    twi + 
    s(slope, k = 5) + 
    tpi + 
    VARI + I(VARI^2)
)

# 2. Match the priors to your specific get_prior output
balanced_priors <- c(
  # 1. The Climate Anchor (Macro-scale)
  prior(normal(3, 0.1), class = "b", coef = "logit_rf"),
  
  # 2. VARI Refinement (Spectral n-shape)
  # We use 'IVARIE2' because your get_prior showed that specific name.
  prior(normal(1, 0.5), class = "b", coef = "VARI"),
  prior(normal(-1, 0.5), class = "b", coef = "IVARIE2"), 
  
  # 3. THE SMOOTH MUZZLE (Global)
  # This controls s(solar), s(twi), and s(slope). 
  # Since brms doesn't see b_solar or b_twi, we control them here.
  prior(student_t(3, 0, 0.3), class = "sds"), 
  
  # 4. Slope's Linear Piece
  # Your get_prior specifically showed 'sslope_1' exists.
  prior(normal(0, 0.5), class = "b", coef = "sslope_1"),
  
  # 5. Fixed Effect for TPI
  prior(normal(0, 1), class = "b", coef = "tpi"),
  prior(normal(0, 1), class = "b", coef = "twi"),
  prior(normal(0, 1), class = "b", coef = "solar"),
  
  # 6. Intercept (Skeptical Baseline)
  prior(normal(-7, 1), class = "Intercept")
)
# 3. Running the Model
fit_birch_final <- brm(
  formula = birch_forestry_formula,
  data = analysis_df,
  family = bernoulli(link = "logit"),
  prior = balanced_priors,
  chains = 4, 
  iter = 2000, 
  warmup = 1000, 
  cores = 10,
  seed = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)


pp_check(fit_birch_final, ndraws = 50)

library(terra)
library(dplyr)

# 1. Prepare the Data Frame
# Extracting xy is key for the re-mapping later
pred_df_full <- as.data.frame(env_stack, xy = TRUE) 
pred_df <- na.omit(pred_df_full)

# 2. Transform the RF Prior
pred_df$logit_rf <- car::logit(pred_df$sdm, adjust = 0.001)

# 3. Generate Posterior Predictions
# ndraws = 50 to capture the uncertainty of your climate draws
epred_mat <- posterior_epred(
  fit_birch_final, 
  newdata = pred_df, 
  ndraws = 10,      
  re_formula = NA
)

# 4. Extract Summary Stats
pred_df$suitability  <- colMeans(epred_mat)
pred_df$uncertainty  <- apply(epred_mat, 2, sd)
pred_df$prob_success <- colMeans(epred_mat > 0.6)

# 5. RE-MAP TO RASTER (The "Anti-Clunk" Method)
# We merge the predictions back to the full coordinate list to fill the NAs with 'NA'
final_df <- left_join(pred_df_full[, c("x", "y")], 
                      pred_df[, c("x", "y", "suitability", "uncertainty", "prob_success")], 
                      by = c("x", "y"))

# 6. Create the SpatRaster
# Use the first layer of your stack as a template for CRS and Resolution
result_raster <- rast(final_df, type = "xyz", crs = crs(env_stack))

# Optional: Rename layers for clarity
names(result_raster) <- c("Suitability_Mean", "Uncertainty_SD", "Prob_Success_60")

writeRaster(result_raster, "models/multiscale_SDM_bayes.tif", overwrite = T)