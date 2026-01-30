# ==============================================================================
# MODEL: WALD-MECHANISTIC BIRCH DISPERSAL (LDD) FOR ICELAND
# Version: 3.0 (Wald/Inverse Gaussian Mechanistic implementation)
# ==============================================================================

library(terra)
library(sf)
library(dplyr)

# 1. ENVIRONMENT & DATA INITIALIZATION
cat("--- STEP 1: INITIALIZING SPATIAL DATA ---\n")
wind_long <- readRDS("data/wind/vedurstofa_wind_rose.rds")
SDM_raw   <- rast("data/broad_scale_rf_prob.tif")
SDM       <- project(SDM_raw, "EPSG:8088", res = 1000)
sdm_vals  <- values(SDM, mat = FALSE)
all_coords <- xyFromCell(SDM, 1:ncell(SDM))

# Spatially enable weather stations
wind_stations <- wind_long %>%
  distinct(id, lon, lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(8088)
station_coords <- st_coordinates(wind_stations)

# ------------------------------------------------------------------------------
# 2. CALIBRATION PARAMETERS (WALD SPECIFIC)
# ------------------------------------------------------------------------------
source_threshold <- 0.2
v_threshold      <- 10    # m/s release threshold
max_dist_m       <- 15000 # Wald allows for longer stable tails (15km)

# Mechanistic Constants for Betula pubescens
H  <- 6           # Average release height in meters (tree height + gust lift)
vt <- 0.58        # Terminal velocity (m/s) for birch samaras
lambda_param <- 2 # Turbulence/Spread parameter (lower = more dispersion)

# ------------------------------------------------------------------------------
# 3. SOURCE IDENTIFICATION
# ------------------------------------------------------------------------------
source_idx    <- which(!is.na(sdm_vals) & sdm_vals >= source_threshold)
source_coords <- all_coords[source_idx, ]
Arrival_vals <- rep(0, ncell(SDM))
Impact_vals  <- rep(0, length(source_idx))

# ------------------------------------------------------------------------------
# 4. DISPERSAL SIMULATION LOOP
# ------------------------------------------------------------------------------
cat("--- STEP 4: EXECUTING WALD SIMULATION ---\n")

for (i in seq_along(source_idx)) {
  src_x <- source_coords[i, 1]; src_y <- source_coords[i, 2]
  S_i   <- sdm_vals[source_idx[i]]
  
  # A. LOCAL WIND FIELD (IDW)
  dists_to_st    <- sqrt((src_x - station_coords[,1])^2 + (src_y - station_coords[,2])^2)
  nearest_st_idx <- order(dists_to_st)[1:3]
  weights        <- 1 / (dists_to_st[nearest_st_idx] + 100)^2 
  w_norm         <- weights / sum(weights)
  target_ids     <- wind_stations$id[nearest_st_idx]
  
  st_data_blended <- wind_long %>%
    filter(id %in% target_ids) %>%
    group_by(direction) %>%
    summarise(
      frequency = sum(frequency * w_norm[match(id, target_ids)]),
      weibull_A = sum(weibull_A * w_norm[match(id, target_ids)]),
      .groups   = 'drop'
    ) %>%
    filter(weibull_A >= v_threshold)
  
  if (nrow(st_data_blended) == 0) next
  
  # B. TARGET CELL FILTERING
  nearby_idx <- which(all_coords[,1] > (src_x - max_dist_m) & all_coords[,1] < (src_x + max_dist_m) &
                        all_coords[,2] > (src_y - max_dist_m) & all_coords[,2] < (src_y + max_dist_m))
  
  dist_m <- sqrt((all_coords[nearby_idx, 1] - src_x)^2 + (all_coords[nearby_idx, 2] - src_y)^2)
  ldd_ring <- which(dist_m <= max_dist_m & dist_m > 0)
  
  if (length(ldd_ring) == 0) next
  
  final_target_idx <- nearby_idx[ldd_ring]
  d <- dist_m[ldd_ring] # Distance in meters for Wald math
  theta <- (atan2(all_coords[final_target_idx, 1] - src_x, 
                  all_coords[final_target_idx, 2] - src_y) * 180 / pi) %% 360
  
  C_local_sum <- rep(0, length(final_target_idx))
  
  # C. PLUME ACCUMULATION (WALD KERNEL)
  for (j in 1:nrow(st_data_blended)) {
    u <- st_data_blended$weibull_A[j] # Mean wind speed for this sector
    dispersal_to <- (st_data_blended$direction[j] + 180) %% 360
    
    # 1. Wald Mechanistic Math
    # mu_prime is the expected mean distance: (Height * Wind) / TerminalVelocity
    mu_prime <- (H * u) / vt
    
    # Probability density of landing at distance 'd'
    # f(d) = sqrt(lambda / 2*pi*d^3) * exp(-lambda*(d-mu)^2 / 2*mu^2*d)
    # Note: We use lambda_param as a scaling factor for turbulence
    term1 <- sqrt(lambda_param / (2 * pi * d^3))
    term2 <- exp(-(lambda_param * (d - mu_prime)^2) / (2 * mu_prime^2 * d))
    K_wald <- term1 * term2
    
    # 2. Directional & Energy Weighting
    storm_energy <- (st_data_blended$frequency[j] / 100) * ((u - v_threshold)^3)
    delta_theta  <- pmin(abs(theta - dispersal_to), 360 - abs(theta - dispersal_to))
    W_dir        <- exp(-(delta_theta^2) / (2 * 30^2)) 
    
    T_j <- sdm_vals[final_target_idx]; T_j[is.na(T_j)] <- 0
    
    C_local_sum <- C_local_sum + (S_i * storm_energy * W_dir * K_wald * T_j)
  }
  
  Arrival_vals[final_target_idx] <- Arrival_vals[final_target_idx] + C_local_sum
  Impact_vals[i] <- sum(C_local_sum) 
  
  if (i %% 500 == 0) cat("Progress:", round(100*i/length(source_idx), 1), "%\n")
}

# ------------------------------------------------------------------------------
# 5. EXPORT FINAL RASTER PRODUCTS
# ------------------------------------------------------------------------------
cat("\n--- STEP 5: GENERATING RASTER OUTPUTS ---\n")

# A. ARRIVAL ZONE MAP (Sink potential)
Arrival_final <- rast(SDM)
values(Arrival_final) <- sqrt(Arrival_vals) # Sqrt compression for better mapping
Arrival_final <- mask(Arrival_final, SDM)
if(global(Arrival_final, "max", na.rm=T)[1,1] > 0) Arrival_final <- Arrival_final / global(Arrival_final, "max", na.rm=T)[1,1]

writeRaster(Arrival_final, "data/wind/birch_storm_arrival_zones.tif", overwrite = TRUE)

# B. IMPACT PRIORITY MAP (Source value)
Impact_final <- rast(SDM); values(Impact_final) <- 0
values(Impact_final)[source_idx] <- Impact_vals 
Impact_final <- mask(Impact_final, SDM)
if(global(Impact_final, "max", na.rm=T)[1,1] > 0) Impact_final <- Impact_final / global(Impact_final, "max", na.rm=T)[1,1]

writeRaster(Impact_final, "data/wind/birch_storm_planting_priority.tif", overwrite = TRUE)

cat("Success! LDD Simulation Complete.\n")


library(terra)
library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)

# 1. SELECT A SINGLE HIGH-VALUE SOURCE
set.seed(42)
# Pick a source with high habitat suitability for a "strong" plume
single_src_idx <- source_idx[which.max(sdm_vals[source_idx])]
src_x <- all_coords[single_src_idx, 1]
src_y <- all_coords[single_src_idx, 2]

# 2. ZOOM SETTINGS (15km buffer)
zoom_dist <- 15000 
n_seeds <- 500

# Crop the background for the zoom
crop_ext <- ext(src_x - zoom_dist, src_x + zoom_dist, 
                src_y - zoom_dist, src_y + zoom_dist)
arrival_zoomed <- crop(Arrival_final, crop_ext)
arrival_df <- as.data.frame(arrival_zoomed, xy = TRUE, na.rm = TRUE)
colnames(arrival_df)[3] <- "prob"

# 3. GENERATE PLUME PARTICLES
# Using a specific wind direction (e.g., 225 degrees / SW wind)
dir <- 225 
dist <- rgamma(n_seeds, shape = 2, scale = 3000) # Wald-like spread
ang  <- rnorm(n_seeds, mean = (dir + 180) %% 360, sd = 12) # Tight cone

dest_x <- src_x + dist * sin(ang * pi / 180)
dest_y <- src_y + dist * cos(ang * pi / 180)

# TJ Check (Suitability)
land_coords <- cbind(dest_x, dest_y)
suit_vals <- terra::extract(SDM, land_coords)
suit_vals <- if(ncol(suit_vals) > 1) suit_vals[, 2] else suit_vals[, 1]
suit_vals[is.na(suit_vals)] <- 0
survives <- runif(n_seeds) < suit_vals

# Format flight data
flight_df <- data.frame(
  id = 1:n_seeds,
  start_x = src_x, start_y = src_y,
  end_x = dest_x, end_y = dest_y,
  survives = survives
) %>%
  pivot_longer(cols = c(start_x, end_x), names_to = "stage", values_to = "x") %>%
  mutate(
    y = ifelse(stage == "start_x", start_y, end_y),
    time = ifelse(stage == "start_x", 1, 2),
    alpha_val = case_when(time == 1 ~ 0.6, survives ~ 1, TRUE ~ 0),
    size_val = case_when(time == 1 ~ 0.8, survives ~ 2, TRUE ~ 0)
  )

# 4. RENDER ZOOMED ANIMATION
anim <- ggplot() +
  geom_raster(data = arrival_df, aes(x = x, y = y, fill = prob), alpha = 0.3) +
  scale_fill_viridis_c(option = "mako", guide = "none") +
  
  # The Source Tree
  annotate("point", x = src_x, y = src_y, color = "red", shape = 17, size = 5) +
  
  # The Particles
  geom_point(data = flight_df, 
             aes(x = x, y = y, group = id, size = size_val, alpha = alpha_val), 
             color = "#A3FF00") +
  
  scale_size_identity() +
  scale_alpha_identity() +
  transition_reveal(time) +
  coord_sf(datum = st_crs(8088), 
           xlim = c(src_x - zoom_dist, src_x + zoom_dist),
           ylim = c(src_y - zoom_dist, src_y + zoom_dist)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "black"),
        plot.title = element_text(color = "white", face = "bold", size = 16)) +
  labs(title = "Mechanistic Seed Plume",
       subtitle = "Red triangle = Parent Tree | Green = Established Seedlings")

animate(anim, nframes = 100, fps = 20, width = 800, height = 800)