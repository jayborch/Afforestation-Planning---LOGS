# Load libraries
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(randomForest)
library(mgcv)
library(networkD3)
library(pROC)
library(blockCV)
library(geodata)
library(SIBER)

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



# BAYESIAN FUNDAMENTAL NICHE MODEL (FNM) - CLIMATE ONLY
# =============================================================
library(terra)
library(dplyr)

# 1. VARIABLE SELECTION
# -------------------------------------------------------------
# Restricted to Bio10 (Summer Heat) and Bio12 (Annual Rain).
# This represents "Atmospheric Permission"—the physiological ceiling 
# independent of transient edaphic factors which can be remediated.
target_vars <- c("wc2.1_30s_bio_10", "wc2.1_30s_bio_12") 

# Extract presence-only data to characterize the empirical niche
raw_presence <- pa_climate %>% 
  filter(presence == 1) %>% 
  dplyr::select(all_of(target_vars)) %>%
  na.omit() %>%
  as.matrix()

# 2. SCALING & STANDARDIZATION
# -------------------------------------------------------------
mu_obs <- colMeans(raw_presence)
sd_obs <- apply(raw_presence, 2, sd)
data_scaled <- scale(raw_presence, center = mu_obs, scale = sd_obs)

# 3. BAYESIAN CONJUGATE PRIOR CALCULATION
# -------------------------------------------------------------
# To ensure the FNM captures the full physiological range of the species, 
# the Bayesian priors are calibrated to encompass known primary refugia—
# such as Hallormsstaður (continental/dry) and Þórsmörk (maritime/wet)—
# where long-term exclusion of grazing has allowed the species to 
# express its fundamental niche.

# Physiological Prior Settings (Real Units):
# Bio10: 8.5°C (Thermal optimum for carbon fixation in sub-arctic birch)
# Bio12: 1150mm (Midpoint to bridge the hygric gap between South and North refugia)
  mu_prior_real <- c(10, 1150) 
mu_prior_scaled <- (mu_prior_real - mu_obs) / sd_obs

n <- nrow(data_scaled)
p <- ncol(data_scaled) 

# kappa_0 (Prior Strength): Set to 0.1 (Weakly Informative).
# This allows the model to stretch the niche ellipse to encompass high-rainfall 
# areas like Þórsmörk by trusting the spatial density of presence points 
# while still being 'informed' by the theoretical optimum.
kappa_0 <- 0.1           
nu_0    <- p + 2       
Psi_0   <- diag(p)     

# Analytical Posterior Update
x_bar <- colMeans(data_scaled)
S <- cov(data_scaled) * (n - 1)

# The Posterior Mean (The 'Physiological Centroid' of the species)
mu_post <- (kappa_0 * mu_prior_scaled + n * x_bar) / (kappa_0 + n)

# The Posterior Covariance (Defining the breadth and orientation of the niche)
Psi_n <- Psi_0 + S + (kappa_0 * n / (kappa_0 + n)) * (x_bar - mu_prior_scaled) %*% t(x_bar - mu_prior_scaled)
sigma_post <- Psi_n / (nu_0 + n)

# 4. GEOGRAPHIC PROJECTION
# -------------------------------------------------------------
# Subset and scale the environmental raster stack
env_sub <- env_stack[[target_vars]]
env_scaled <- (env_sub - mu_obs) / sd_obs

# Apply the Mahalanobis + Chi-Square logic to every pixel
fnm_suitability <- app(env_scaled, fun = function(x) {
  if(any(is.na(x))) return(NA)
  
  # Calculate Mahalanobis distance D^2 from the Bayesian center
  d2 <- mahalanobis(matrix(x, nrow = 1), center = mu_post, cov = sigma_post)
  
  # Convert D^2 to 0-1 probability using Chi-Square distribution (df = 2)
  # A D^2 of 0 (the centroid) maps to a probability of 1.0.
  p_val <- pchisq(d2, df = p, lower.tail = FALSE)
  return(p_val)
})

names(fnm_suitability) <- "Fundamental_Suitability_Climate"



# 5. VISUALIZATION
# -------------------------------------------------------------
# =============================================================
# INTERACTIVE FNM VISUALIZATION (95% BAYESIAN BOUNDARY)
# =============================================================
# =============================================================
# RASTER-BASED FNM VISUALIZATION (No Polygons)
# =============================================================



library(leaflet)
library(terra)

# 1. Prepare the Binary Raster Boundary (The 0.05 Threshold)
# -------------------------------------------------------------
# Create a mask where 1 = inside 95% niche, NA = outside
# We use NA for 'outside' so the map doesn't cover up the satellite imagery
fnm_boundary_raster <- fnm_suitability >= 0.05
fnm_boundary_raster[fnm_boundary_raster == 0] <- NA 

# Project to Web Mercator for Leaflet
boundary_leaflet <- project(fnm_boundary_raster, "EPSG:3857", method = "near")

# 2. Prepare the Suitability Gradient
# -------------------------------------------------------------
fnm_leaflet <- project(fnm_suitability, "EPSG:3857")

# Palette for the Gradient (Subtle)
pal_grad <- colorNumeric(
  palette = "YlGn", 
  domain = c(0, 1), 
  na.color = "transparent"
)

# Palette for the Boundary (High Contrast Red)
pal_bound <- colorFactor(
  palette = c("#ff0000"), # Pure Red
  domain = c(1),
  na.color = "transparent"
)

# 3. Build the Map
# -------------------------------------------------------------
leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  
  # Layer 1: The Fuzzy Suitability Gradient
  addRasterImage(
    fnm_leaflet, 
    colors = pal_grad, 
    opacity = 0.5, 
    group = "Suitability Gradient"
  ) %>%
  
  # Layer 2: The "Hard" 95% Boundary Raster
  # This shows specifically the pixels that meet the statistical cutoff
  addRasterImage(
    boundary_leaflet, 
    colors = pal_bound, 
    opacity = 0.8, 
    group = "95% Boundary (Pixels)"
  ) %>%
  
  addLayersControl(
    overlayGroups = c("Suitability Gradient", "95% Boundary (Pixels)"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(pal = pal_grad, values = c(0, 1), title = "FNM Prob")



# =============================================================
# BAYESIAN POSTERIOR DISPLACEMENT: THE 'TUG-OF-WAR'
# =============================================================

# 1. DEFINE THE ELLIPSE GENERATOR FUNCTION
# -------------------------------------------------------------
get_ellipse <- function(mu, sigma, level = 0.95) {
  mu_vec <- matrix(as.vector(mu), ncol = 1)
  theta <- seq(0, 2 * pi, length.out = 100)
  circle <- rbind(cos(theta), sin(theta)) 
  radius <- sqrt(qchisq(level, df = 2))
  L <- t(chol(sigma)) 
  ellipse_pts <- mu_vec %*% matrix(1, ncol = 100) + radius * (L %*% circle)
  df <- as.data.frame(t(ellipse_pts))
  colnames(df) <- c("Bio10_scaled", "Bio12_scaled")
  return(df)
}

# 2. GENERATE POSTERIOR SAMPLES (The 'Eggs')
# -------------------------------------------------------------
library(mvnfast)
set.seed(123) # For reproducibility
n_to_draw <- 1000

# We sample the uncertainty of the mean (mu) based on the observed sample size (n)
sampled_eggs <- lapply(1:n_to_draw, function(i) {
  list(
    mu = rmvn(1, mu_post, sigma_post / sqrt(n)), 
    sigma = sigma_post
  )
})

# 3. CONSTRUCT PLOTTING DATAFRAMES
# -------------------------------------------------------------
all_ellipses <- data.frame()
centers <- data.frame() 

for(i in 1:length(sampled_eggs)) {
  curr_mu <- as.vector(sampled_eggs[[i]]$mu)
  curr_sigma <- sampled_eggs[[i]]$sigma
  
  # Store center point (scaled)
  centers <- rbind(centers, data.frame(
    Bio10_scaled = curr_mu[1], 
    Bio12_scaled = curr_mu[2], 
    draw_id = i
  ))
  
  # Generate 95% boundary
  temp_df <- get_ellipse(curr_mu, curr_sigma, level = 0.95)
  temp_df$draw_id <- i
  all_ellipses <- rbind(all_ellipses, temp_df)
}

# 4. BACK-TRANSFORM TO REAL CLIMATE UNITS
# -------------------------------------------------------------
all_ellipses_real <- all_ellipses %>%
  mutate(
    Bio10 = Bio10_scaled * sd_obs[1] + mu_obs[1],
    Bio12 = Bio12_scaled * sd_obs[2] + mu_obs[2]
  )

centers_real <- centers %>%
  mutate(
    Bio10 = Bio10_scaled * sd_obs[1] + mu_obs[1],
    Bio12 = Bio12_scaled * sd_obs[2] + mu_obs[2]
  )

# 5. VISUALIZE THE DISPLACEMENT
# -------------------------------------------------------------
prior_target <- c(x = 10, y = 1150)

ggplot() +
  # BACKGROUND: The actual presence points (Realized Niche)
  geom_point(data = as.data.frame(raw_presence), 
             aes(x = wc2.1_30s_bio_10, y = wc2.1_30s_bio_12), 
             alpha = 0.05, color = "grey60", size = 0.5) +
  
  # THE TUG-OF-WAR: Lines from each draw's center to the Prior
  geom_segment(data = centers_real, 
               aes(x = Bio10, y = Bio12, xend = prior_target["x"], yend = prior_target["y"]),
               color = "firebrick1", linewidth = 0.4, alpha = 0.3) +
  
  # THE UNCERTAINTY CLOUD: Individual Bayesian Draws
  geom_path(data = all_ellipses_real, 
            aes(x = Bio10, y = Bio12, group = draw_id), 
            color = "#2E8B57", alpha = 0.1) +
  
  # THE CONSENSUS: The Master Posterior Mean Ellipse
  stat_ellipse(data = as.data.frame(raw_presence), 
               aes(x = wc2.1_30s_bio_10, y = wc2.1_30s_bio_12), 
               color = "black", linewidth = 1.2) +
  
  # THE TARGET: Theoretical Physiological Prior
  annotate("point", x = prior_target["x"], y = prior_target["y"], 
           color = "red", shape = 17, size = 6) +
  
  labs(title = "Niche Displacement & Historical Validation",
       subtitle = "Green cloud represents model uncertainty; Red lines represent the 'pull' of physiological theory.",
       x = "Summer Mean Temperature (°C)", 
       y = "Annual Precipitation (mm)") +
  theme_minimal()



library(plotly)
library(mvnfast)

# Create a dense grid for the mesh
temp_seq <- seq(7, 11, length.out = 80)
rain_seq <- seq(400, 1600, length.out = 80)
grid <- expand.grid(Bio10 = temp_seq, Bio12 = rain_seq)

# Scale and calculate the Bayesian Density
grid_scaled <- scale(grid, center = mu_obs, scale = sd_obs)
grid$suitability <- dmvn(grid_scaled, mu = mu_post, sigma = sigma_post)

# Reshape into a matrix for the mesh
z_matrix <- matrix(grid$suitability, nrow = 80, ncol = 80)

# Build the Mesh Plot
p <- plot_ly(x = ~temp_seq, y = ~rain_seq, z = ~z_matrix) %>%
  add_trace(
    type = "surface",
    contours = list(
      z = list(show = TRUE, usecolormap = TRUE, highlightcolor = "#ff0000", project = list(z = TRUE))
    ),
    opacity = 0.8,
    colorscale = "Viridis",
    line = list(width = 0.5, color = "white"), # This creates the 'mesh' look
    showscale = FALSE
  ) %>%
  layout(
    title = "3D Fundamental Niche Mesh",
    scene = list(
      xaxis = list(title = "Summer Temp (°C)"),
      yaxis = list(title = "Annual Rain (mm)"),
      zaxis = list(title = "Suitability Index"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2))
    )
  )

p