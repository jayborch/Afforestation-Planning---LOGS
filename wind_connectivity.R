library(terra)
library(sf)
library(dplyr)

# ----- 1. DATA LOADING -----
cat("--- DATA LOADING ---\n")
wind_long <- readRDS("data/wind/vedurstofa_wind_rose.rds")
SDM_raw <- rast("data/broad_scale_rf_prob.tif")

# Project SDM to 1km grid (EPSG:8088)
SDM <- project(SDM_raw, "EPSG:8088", res = 1000)
sdm_vals <- values(SDM, mat = FALSE)

# Pre-calculate coordinates for the whole grid
all_coords <- xyFromCell(SDM, 1:ncell(SDM))

# Get unique stations and coordinates
wind_stations <- wind_long %>%
  distinct(id, lon, lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(8088)

station_coords <- st_coordinates(wind_stations)

# ----- 2. PARAMETERS -----
source_threshold <- 0.2
max_dist_km      <- 10
max_dist_m       <- max_dist_km * 1000
beta             <- 2.0  # Power-law distance decay
sigma_angle      <- 15   # Sharp 15-degree cone for directional plumes
k_neighbors      <- 3    # Number of stations to blend for IDW smoothing

# ----- 3. PREP VECTORS -----
source_idx <- which(!is.na(sdm_vals) & sdm_vals >= source_threshold)
source_coords <- all_coords[source_idx, ]

Arrival_vals <- rep(0, ncell(SDM))
Impact_vals <- rep(0, length(source_idx))

# ----- 4. DUAL-STREAM SIMULATION (WITH IDW BLENDING) -----
cat("Simulating dispersal for", length(source_idx), "sources using IDW blending...\n")

for (i in seq_along(source_idx)) {
  src_x <- source_coords[i, 1]
  src_y <- source_coords[i, 2]
  S_i   <- sdm_vals[source_idx[i]]
  
  # A. IDW STATION BLENDING
  # Calculate distances to all stations
  dists_to_st <- sqrt((src_x - station_coords[,1])^2 + (src_y - station_coords[,2])^2)
  nearest_st_idx <- order(dists_to_st)[1:k_neighbors]
  
  # Calculate weights (1/d^2)
  # Adding 100m to distance to avoid infinity at station location
  w <- 1 / (dists_to_st[nearest_st_idx] + 100)^2
  w_norm <- w / sum(w)
  
  # Blend wind data from the top K stations
  target_ids <- wind_stations$id[nearest_st_idx]
  
  st_data_blended <- wind_long %>%
    filter(id %in% target_ids) %>%
    group_by(direction) %>%
    summarise(
      # Weighted average of wind frequency and energy
      frequency = sum(frequency * w_norm[match(id, target_ids)]),
      weibull_A = sum(weibull_A * w_norm[match(id, target_ids)]),
      .groups = 'drop'
    )
  
  # B. Spatial Filter (10km bounding box)
  nearby_idx <- which(all_coords[,1] > (src_x - max_dist_m) & all_coords[,1] < (src_x + max_dist_m) &
                        all_coords[,2] > (src_y - max_dist_m) & all_coords[,2] < (src_y + max_dist_m))
  
  if (length(nearby_idx) == 0) next
  
  # C. Calculate Distances
  dist_sq <- (all_coords[nearby_idx, 1] - src_x)^2 + (all_coords[nearby_idx, 2] - src_y)^2
  ldd_ring <- which(dist_sq <= (max_dist_m^2))
  
  if (length(ldd_ring) == 0) next
  
  final_target_idx <- nearby_idx[ldd_ring]
  d_km <- sqrt(dist_sq[ldd_ring]) / 1000
  d_km[d_km < 0.5] <- 0.5  # Cap distance to prevent infinity
  
  # D. Calculate Angles
  theta <- (atan2(all_coords[final_target_idx, 1] - src_x, 
                  all_coords[final_target_idx, 2] - src_y) * 180 / pi) %% 360
  
  C_local_sum <- rep(0, length(final_target_idx))
  
  # E. Iterate through BLENDED wind sectors
  for (j in 1:nrow(st_data_blended)) {
    wind_from <- st_data_blended$direction[j]
    dispersal_to <- (wind_from + 180) %% 360
    
    freq <- st_data_blended$frequency[j] / 100 
    A    <- st_data_blended$weibull_A[j]
    
    if (freq <= 0 || A <= 0) next
    
    storm_energy <- freq * (A^3)
    delta_theta <- pmin(abs(theta - dispersal_to), 360 - abs(theta - dispersal_to))
    W_dir <- exp(-(delta_theta^2) / (2 * sigma_angle^2))
    
    K_ldd <- d_km^(-beta)
    T_j   <- sdm_vals[final_target_idx]; T_j[is.na(T_j)] <- 0
    
    C_local_sum <- C_local_sum + (S_i * storm_energy * W_dir * K_ldd * T_j)
  }
  
  Arrival_vals[final_target_idx] <- Arrival_vals[final_target_idx] + C_local_sum
  Impact_vals[i] <- sum(C_local_sum)
  
  if (i %% 500 == 0) cat("Progress:", round(100*i/length(source_idx), 1), "% | Current Max Arrival:", max(Arrival_vals), "\n")
}

# ----- 5. OUTPUT GENERATION -----
cat("\n--- Preparing Final Outputs ---\n")

# A. EXPORT ARRIVAL MAP
Arrival_rast <- rast(SDM)
values(Arrival_rast) <- Arrival_vals
Arrival_masked <- mask(Arrival_rast, SDM)
Arrival_final <- sqrt(Arrival_masked)
max_arr <- global(Arrival_final, "max", na.rm = TRUE)[1,1]
if(max_arr > 0) Arrival_final <- Arrival_final / max_arr

writeRaster(Arrival_final, "data/wind/birch_downstream_arrival_zones.tif", 
            overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)

# B. EXPORT IMPACT MAP
Impact_rast <- rast(SDM)
values(Impact_rast) <- 0
values(Impact_rast)[source_idx] <- Impact_vals 
Impact_masked <- mask(Impact_rast, SDM)
max_imp <- global(Impact_masked, "max", na.rm = TRUE)[1,1]
if(max_imp > 0) Impact_final <- Impact_masked / max_imp

writeRaster(Impact_final, "data/wind/birch_upstream_planting_priority.tif", 
            overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)

cat("Success! Smooth maps generated using IDW wind blending.\n")