
### Site accessibility model fdor afforestation

# This script generates a landscape resistance surface for accessibility analysis. 
# It combines landcover-based resistance (simplified from Vistlendi classes) 
# and terrain-based resistance (slope-based), normalizes both layers, 
# multiplies them to create a combined resistance raster, and incorporates roads. 
# Finally, it calculates accumulated cost distances from the road network 
# using WhiteboxTools for site accessibility modeling.

### Roads network data
# https://gatt.lmi.is/geonetwork/srv/eng/catalog.search#/metadata/41D3314D-55AF-4FE3-9A9F-B9CB18AC5CAA



library(terra)
library(foreign)
library(dplyr)
library(whitebox)
library(sf)
library(leaflet)

wbt_init()

# --- 1. Define simplified resistance classes using Vistlendi --- T
cost_weights <- data.frame(
  Vistlendi = c(
    "Melar og sandlendi",
    "Skriður og klettar",
    "Eyrar",
    "Moslendi",
    "Hraunlendi",
    "Moldir",
    "Strandlendi",
    "Votlendi",
    "Graslendi",
    "Mólendi",
    "Jarðhitasvæði",
    "Aðrar landgerðir",
    "Skóglendi",
    "Ferskvatn",
    "Jöklar",
    "Fjörur"
  ),
  cost = c(
    12,   # Melar og sandlendi (sands & gravels)
    90,   # Skriður og klettar (scree, boulders, cliffs)
    25,   # Eyrar (river plains, unstable sediment)
    35,   # Moslendi (moss heaths, soft ground)
    70,   # Hraunlendi (lava fields, rough terrain)
    20,   # Moldir (loamy soils – generally OK)
    40,   # Strandlendi (coastal habitats, soft sediments)
    95,   # Votlendi (wetlands, mires – high bogging risk)
    15,   # Graslendi (grasslands – good access)
    30,   # Mólendi (heathlands – moderate)
    60,   # Jarðhitasvæði (geothermal, unstable ground)
    5,    # Aðrar landgerðir (urban, agriculture, roads, tracks)
    45,   # Skóglendi (birch forest – obstacles)
    NA,  # Ferskvatn (lakes/rivers – impassable)
    NA,  # Jöklar (glaciers – impassable)
    50    # Fjörur (tidal flats – soft, risky)
  )
)

# --- 2. Load habitat raster ---
habitats <- rast("data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif")

# Load raster attribute table (.dbf)
rat <- as.data.frame(
  read.dbf("data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif.vat.dbf")
)

# --- 3. Join RAT with simplified resistance values ---
rat_joined <- rat %>%
  left_join(cost_weights, by = "Vistlendi")

# --- 4. Build lookup matrix: Value → cost ---
lookup_mat <- as.matrix(rat_joined[, c("Value", "cost")])

# --- 5. Classify raster to create resistance surface ---
resistance_landclass <- classify(habitats, rcl = lookup_mat)

### Now create resistance from terrain

# design some resistances for slopes for afforestation planning
dem <- rast("data/IslandsDEMv1.0_10x10m_isn2016_zmasl.tif")
slope <- terrain(dem, v = "slope", unit = "degrees")


# --- slope resistance function ---
resistance_slope <- function(slope_deg, a = 20, S0 = 18, b = 3) {
  # a  = max additional resistance (1 + a = asymptotic maximum)
  # S0 = slope at which resistance increases fastest
  # b  = steepness of logistic curve
  1 + a / (1 + exp(-(slope_deg - S0) / b))
}
# --- Create slope sequence ---
slope_vals <- seq(0, 45, by = 0.5)
res_vals <- resistance_slope(slope_vals)

# --- Plot ---
plot(
  slope_vals, res_vals, type = "l", lwd = 2, col = "darkred",
  xlab = "Slope (degrees)", ylab = "Resistance multiplier",
  main = "6-Wheeler Slope Resistance Function"
)

terrain_resistance_raster <- app(slope, resistance_slope)

### combine models together
resistance_landclass_proj <- project(resistance_landclass, terrain_resistance_raster)

# normalize both rasters
mm <- minmax(terrain_resistance_raster)
terrain_norm <- (terrain_resistance_raster - mm[1,1]) / (mm[2,1] - mm[1,1])

mm <- minmax(resistance_landclass_proj)
landcover_norm <- (resistance_landclass_proj - mm[1,1]) / (mm[2,1] - mm[1,1])

# combine by multiplying normalized rasters
combined_resistance <- terrain_norm * landcover_norm

### rasterise road network

roads <- vect("data/roads/is_50v_samgongur_epsg_8088.gpkg/is_50v_samgongur_epsg_8088.gpkg", layer = "samgongur_linur")
roads_proj <- project(roads, crs(combined_resistance))
road_raster <- rasterize(roads_proj, combined_resistance, field = 1)


# --- 3. Continuous Tiling Logic ---
# This creates a seamless grid. Every row in tile_mat is now treated as a tile.
tile_mat <- getTileExtents(combined_resistance, y = c(5000, 5000), buffer = 10000)

# We create the tile_df using every single row provided by the tiling function
tile_df <- data.frame(
  b_xmin = tile_mat[, 1], b_xmax = tile_mat[, 2], # Buffered bounds
  b_ymin = tile_mat[, 3], b_ymax = tile_mat[, 4],
  xmin   = tile_mat[, 1] + 10000, xmax   = tile_mat[, 2] - 10000, # Reconstruct original bounds
  ymin   = tile_mat[, 3] + 10000, ymax   = tile_mat[, 4] - 10000
)

# --- 4. Filtering and Processing ---
dir.create("data/site_accessibility/temp", recursive = TRUE, showWarnings = FALSE)
dir.create("data/site_accessibility/final", recursive = TRUE, showWarnings = FALSE)

# Clean up any old partial results to ensure a fresh, clean mosaic
# unlink("data/site_accessibility/final/*") 

cat("Aligning master rasters (pixel-perfect sync)...\n")
road_raster <- resample(road_raster, combined_resistance, method = "near")

cat("Starting robust tiled processing...\n")

for (i in 1:nrow(tile_df)) {
  f_tile  <- paste0("data/site_accessibility/final/tile_", i, ".tif")
  
  # RESUME: skip if already done
  if (file.exists(f_tile)) next
  
  # 1. Define buffered extent
  e_buff_raw <- ext(as.numeric(tile_df[i, 1:4]))
  
  # 2. Crop Resistance with 'extend=TRUE' for coastlines
  res_sub <- try({
    tmp <- crop(combined_resistance, e_buff_raw, snap = "out", extend = TRUE)
    if(all(is.na(values(tmp)))) return(NULL)
    tmp
  }, silent = TRUE)
  
  if (inherits(res_sub, "try-error") || is.null(res_sub)) next
  
  # 3. PIXEL-SYNC: Force road_sub to match res_sub's exact grid
  road_sub <- crop(road_raster, res_sub, snap = "near")
  
  if (!all(dim(res_sub) == dim(road_sub))) {
    road_sub <- resample(road_sub, res_sub, method = "near")
  }
  ext(road_sub) <- ext(res_sub)
  
  # Skip if no roads are found (essential for cost distance)
  if (all(is.na(values(road_sub)))) next
  
  cat("Processing tile", i, "of", nrow(tile_df), "\n")
  
  # --- File Paths ---
  t_res   <- paste0("data/site_accessibility/temp/res_", i, ".tif")
  t_road  <- paste0("data/site_accessibility/temp/road_", i, ".tif")
  t_accum <- paste0("data/site_accessibility/temp/accum_", i, ".tif")
  t_back  <- paste0("data/site_accessibility/temp/back_", i, ".tif")
  
  writeRaster(res_sub, t_res, overwrite = TRUE, datatype = "FLT4S")
  writeRaster(road_sub, t_road, overwrite = TRUE, datatype = "FLT4S")
  
  # 4. Whitebox Calculation
  wbt_cost_distance(
    cost = t_res, 
    source = t_road, 
    out_accum = t_accum,
    out_backlink = t_back
  )
  
  # 5. TRIM AND SAVE: Remove the 5km buffer
  if (file.exists(t_accum)) {
    accum_r <- rast(t_accum)
    e_orig_raw <- ext(as.numeric(tile_df[i, 5:8]))
    
    # Use extend=TRUE to handle coastal tiles hanging into the void
    final_r <- try(crop(accum_r, e_orig_raw, snap = "near", extend = TRUE), silent = TRUE)
    
    if (!inherits(final_r, "try-error") && !is.null(final_r)) {
      if (!all(is.na(values(final_r)))) {
        writeRaster(final_r, f_tile, overwrite = TRUE, datatype = "FLT4S")
      }
    }
  }
  
  # Cleanup temp files for this tile
  temp_to_clean <- c(t_res, t_road, t_accum, t_back)
  file.remove(temp_to_clean[file.exists(temp_to_clean)])
}

# --- 5. Reassemble the Tiles ---
tile_files <- list.files("data/site_accessibility/final", pattern = "\\.tif$", full.names = TRUE)

if (length(tile_files) > 0) {
  vrt_path <- "data/site_accessibility/iceland_accessibility_10m.vrt"
  vrt(tile_files, vrt_path, overwrite = TRUE)
  
  # Final Plotting
  final_map <- rast(vrt_path)
  # Set 0 to NA for better visualization (ocean/roads) if needed
  # final_map[final_map == 0] <- NA 
  
  plot(final_map, main = "Complete Accessibility Model", col = terrain.colors(50))
}

library(terra)

tile_files <- list.files(
  "data/site_accessibility/final",
  pattern = "\\.tif$",
  full.names = TRUE
)

if (length(tile_files) > 0) {
  
  vrt_path <- "data/site_accessibility/iceland_accessibility_10m.vrt"
  tif_path <- "data/site_accessibility/iceland_accessibility_10m.tif"
  
  # 1. Build virtual raster (no alignment issues)
  vrt(tile_files, vrt_path, overwrite = TRUE)
  
  # 2. Load VRT
  final_map <- rast(vrt_path)
  
  # Optional: set 0 to NA
  # final_map[final_map == 0] <- NA
  
  # 3. Write out a real GeoTIFF
  writeRaster(
    final_map,
    tif_path,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
  )
  
  # Sanity check
  plot(final_map, main = "Complete Accessibility Model", col = terrain.colors(50))
}
