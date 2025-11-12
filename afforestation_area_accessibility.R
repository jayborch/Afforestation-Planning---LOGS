

library(terra)
library(foreign)
library(dplyr)
library(whitebox)

wbt_init() 

# first need to create a map of resistance across iceland using land classes

cost_weights <- data.frame(
  Value = 1:65,  
  Vistgerð = c(
    "L1.1 Eyðimelavist",
    "L1.2 Grasmelavist",
    "L1.3 Mosamelavist",
    "L1.4 Víðimelavist",
    "L1.5 Sanda- og vikravist",
    "L3.3 Ljónslappaskriðuvist",
    "L3.1 Urðarskriðuvist",
    "L3.2 Grasvíðiskriðuvist",
    "L4.2 Auravist",
    "L4.1 Eyravist",
    "L5.1 Hélumosavist",
    "L5.2 Melagambravist",
    "L5.3 Hraungambravist",
    "L6.2 Fléttuhraunavist",
    "L6.1 Eyðihraunavist",
    "L6.4 Lynghraunavist",
    "L6.3 Mosahraunavist",
    "L2.1 Moldavist",
    "L7.7 Sjávarkletta- og eyjavist",
    "L7.4 Grashólavist",
    "L7.6 Gulstararfitjavist",
    "L7.3 Strandmelhólavist",
    "L7.5 Sjávarfitjungsvist",
    "L7.2 Malarstrandarvist",
    "L7.1 Sandstrandarvist",
    "L8.1 Dýjavist",
    "L8.14 Gulstararflóavist",
    "L8.4 Hrossanálarvist",
    "L8.11 Brokflóavist",
    "L8.10 Hengistararflóavist",
    "L8.12 Starungsflóavist",
    "L8.9 Starungsmýravist",
    "L8.2 Rekjuvist",
    "L8.6 Runnamýravist á láglendi",
    "L8.5 Runnamýravist á hálendi",
    "L8.8 Rústamýravist",
    "L8.3 Sandmýravist",
    "L8.13 Tjarnastararflóavist",
    "L9.7 Blómgresisvist",
    "L9.3 Bugðupuntsvist",
    "L9.5 Grasengjavist",
    "L9.2 Finnungsvist",
    "L9.6 Língresis- og vingulsvist",
    "L9.4 Snarrótarvist",
    "L9.1 Stinnastararvist",
    "L10.6 Fjalldrapamóavist",
    "L10.5 Fléttumóavist",
    "L10.4 Grasmóavist",
    "L10.7 Lyngmóavist á hálendi",
    "L10.8 Lyngmóavist á láglendi",
    "L10.1 Mosamóavist",
    "L10.2 Flagmóavist",
    "L11.1 Ræktarlönd",               
    "L11.2 Túngrónir reitir",          
    "L12.1 Vegsvæði",                 
    "L12.2 Þéttbýli / byggð",         
    "L12.3 Flugvellir / iðnaðarsvæði",
    "L12.4 Námur",                    
    "L12.5 Vökvafylltar lón",         
    "L12.6 Skógræktarsvæði",          
    "L14.1 Þéttbýli og annað manngert land",
    "L14.4 Alaskalúpína",
    "V1 Vötn",
    "L11 Birkiskógur",
    "L14.3 Skógrækt"
    
  ),
  cost = c(
    # L1 – Sandy/gravel
    5,5,5,5,5,
    # L3 – Scree/landslides
    8,8,8,
    # L4 – Fluvial
    6,6,
    # L5 – Moss heath
    4,4,4,
    # L6 – Lava
    7,7,7,7,
    # L2 – Loam
    3,
    # L7 – Coastal
    6,4,5,5,5,5,5,
    # L8 – Wetlands/mire
    9,9,9,9,9,9,9,9,9,9,9,9,9,9,
    # L9 – Grassland
    2,2,2,2,2,2,2,
    # L10 – Heathlands
    4,4,4,4,4,4,4,
    # L11 – Agriculture
    2,2,
    # L12 – Anthropogenic
    1,15,15,15,15,4,
    # L14 / V1 / Forests (corrected to 4 values)
    15,5,20,4
  )
)
# 1. Load raster
habitats <- rast("data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif")
rat <- as.data.frame(read.dbf("data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif.vat.dbf"))

rat_joined <- rat %>%
  left_join(cost_weights, by = "Vistgerð")

lookup_mat <- as.matrix(rat_joined[, c("Value.x", "cost")])
resistance_landclass <- classify(habitats, rcl = lookup_mat)

### Now create resistance from terrain

# design some resistances for slopes for afforestation planning
dem <- rast("data/IslandsDEMv1.0_10x10m_isn2016_zmasl.tif")
slope <- terrain(dem, v = "slope", unit = "degrees")


# --- slope resistance function ---
resistance_slope <- function(slope_deg, a = 4, S0 = 20, b = 4) {
  # a  = max additional resistance (1 + a = asymptotic maximum)
  # S0 = slope where resistance increases fastest (inflection)
  # b  = steepness of the curve
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
combined_resistance <- terrain_resistance_raster + resistance_landclass_proj

### rasterise road network

roads <- vect("data/roads/is_50v_samgongur_epsg_8088.gpkg/is_50v_samgongur_epsg_8088.gpkg", layer = "samgongur_linur")
roads_proj <- project(roads, crs(combined_resistance))
road_raster <- rasterize(roads_proj, resistance_raster, field = 1)

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

# Load rasters
accum <- rast("data/site_accessibility/accumulated_cost_to_roads.tif")
backlink <- rast("data/site_accessibility/accumulated_cost_to_roads_backlink.tif")

# Sample valid starting points
valid_cells <- which(!is.na(values(accum)))
set.seed(42)
sample_cells <- sample(valid_cells, size = 500)
xy_start <- xyFromCell(accum, sample_cells)

# Load rasters
accum <- rast("data/site_accessibility/accumulated_cost_to_roads.tif")
backlink <- rast("data/site_accessibility/accumulated_cost_to_roads_backlink.tif")

# Identify only valid cells in accumulated cost raster
valid_cells <- which(!is.na(values(accum)))

# Sample particles only from valid cells
n_particles <- 2000
set.seed(42)
sample_cells <- sample(valid_cells, size = n_particles)

# Convert sampled cells to XY coordinates
xy_start <- xyFromCell(accum, sample_cells)

# Particle tracing parameters
max_steps <- 1000  # safety limit to avoid infinite loops
paths <- vector("list", n_particles)

# Precompute backlink values as a vector for faster access
backlink_vals <- values(backlink)

# Particle tracing loop
for (i in seq_len(n_particles)) {
  particle_path <- xy_start[i, , drop = FALSE]
  current_cell <- sample_cells[i]
  steps <- 0
  
  while (!is.na(current_cell) && backlink_vals[current_cell] != 0 && steps < max_steps) {
    # D8 neighbors
    nbrs <- adjacent(backlink, current_cell, directions = 8, pairs = TRUE)
    if (nrow(nbrs) == 0) break
    
    # Follow backlink: find neighbor with matching backlink value
    next_idx <- which(backlink_vals[nbrs[,2]] == backlink_vals[current_cell])
    if (length(next_idx) == 0) break
    
    next_cell <- nbrs[next_idx[1], 2]  # take the first matching neighbor
    xy_next <- xyFromCell(backlink, next_cell)
    particle_path <- rbind(particle_path, xy_next)
    
    current_cell <- next_cell
    steps <- steps + 1
  }
  
  paths[[i]] <- particle_path
}

# Plot the raster as background
plot(accum, main = "Particle Paths to Roads", col = terrain.colors(50))

# Overlay paths
for(i in seq_along(paths)) {
  lines(paths[[i]][,1], paths[[i]][,2], col = rgb(1,0,0,0.5))
}

