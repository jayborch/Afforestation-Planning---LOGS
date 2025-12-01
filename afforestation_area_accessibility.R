
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

### test sample
reyk <- vect("data/GEE/roi_rykb_subset.shp")
reyk_project <- project(reyk, crs(combined_resistance))

combined_resistance <- crop(combined_resistance, reyk_project)
road_raster <- crop(road_raster, reyk_project)

writeRaster(combined_resistance, "data/site_accessibility/resistance.tif", overwrite = TRUE)
writeRaster(road_raster, "data/site_accessibility/roads.tif", overwrite = TRUE)

# Run accumulated cost
wbt_cost_distance(
  cost = "data/site_accessibility/resistance.tif",
  source = "data/site_accessibility/roads.tif",
  out_accum = "data/site_accessibility/accumulated_cost_to_roads.tif",
  out_backlink = "data/site_accessibility/accumulated_cost_to_roads_backlink.tif", 
  verbose_mode = TRUE
)

cost_accum = rast("data/site_accessibility/accumulated_cost_to_roads.tif")

