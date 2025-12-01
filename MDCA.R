
library(terra)
library(foreign)
library(dplyr)

dem <- rast("data/IslandsDEMv1.0_10x10m_isn2016_zmasl.tif")

### species distribution model-----------------------------------------

SDM <- rast("data/multiscale_betula_reykjanes_10m_embeddings.tif")

### read in site accessibility map and take 95 percentile and normalise-----------------------------------------

# Load raster
road_dist <- rast("data/site_accessibility/accumulated_cost_to_roads.tif")
road_dist <- road_dist[[1]]  # ensure single layer

# 1. Extract raster values into an array/vector
vals <- values(road_dist, mat = FALSE)
vals <- vals[!is.na(vals)]  # remove NAs

# 2. Calculate 95th percentile (or any desired percentile)
cap_value <- quantile(vals, 0.95)

# 3. Cap raster values at this percentile
road_dist_capped <- clamp(road_dist, lower = -Inf, upper = cap_value)

# 4. Add constant and log-transform (handles zeros safely)
road_dist_log <- log1p(road_dist_capped + 1)

# 5. Min-max normalize
min_val <- minmax(road_dist_log)[1]
max_val <- minmax(road_dist_log)[2]
road_dist_norm <- (road_dist_log - min_val) / (max_val - min_val)

road_dist_norm <- 1 - road_dist_norm

###  Population-weighted Proximity Preference-----------------------------------------

pop_weighted_proxity_town <- rast("data/pop_weighted_proxity_town.tif")

min_val <- minmax(pop_weighted_proxity_town)[1]
max_val <- minmax(pop_weighted_proxity_town)[2]
pop_weighted_proxity_town_norm <- (pop_weighted_proxity_town - min_val) / (max_val - min_val)



###  Archaeological remains buffer -----------------------------------------

remains <- rast("data/minjastofnun/minjastofnun_raster.tif")


###  bird areas -----------------------------------------
 
IBAs <- rast("data/bird areas/important_bird_areas.tif") ### soft constraint (penalise)
IBAs <- 1 - IBAs

protected_BAs <- rast("data/bird areas/protected_bird_areas.tif") ### hard constraint (exlcude)+

###  water and infrastructure for masking  -----------------------------------------

# --- 1. Define simplified resistance classes using Vistlendi --- T
mask <- data.frame(
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
    0,   # Melar og sandlendi (sands & gravels)
    0,   # Skriður og klettar (scree, boulders, cliffs)
    0,   # Eyrar (river plains, unstable sediment)
    0,   # Moslendi (moss heaths, soft ground)
    0,   # Hraunlendi (lava fields, rough terrain)
    0,   # Moldir (loamy soils – generally OK)
    0,   # Strandlendi (coastal habitats, soft sediments)
    0,   # Votlendi (wetlands, mires – high bogging risk)
    0,   # Graslendi (grasslands – good access)
    0,   # Mólendi (heathlands – moderate)
    0,   # Jarðhitasvæði (geothermal, unstable ground)
    1,    # Aðrar landgerðir (urban, agriculture, roads, tracks)
    0,   # Skóglendi (birch forest – obstacles)
    1,  # Ferskvatn (lakes/rivers – impassable)
    1,  # Jöklar (glaciers – impassable)
    0    # Fjörur (tidal flats – soft, risky)
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
  left_join(mask, by = "Vistlendi")

# --- 4. Build lookup matrix: Value → cost ---
lookup_mat <- as.matrix(rat_joined[, c("Value", "cost")])

# --- 5. Classify raster to create resistance surface ---
resistance_landclass <- classify(habitats, rcl = lookup_mat)
resistance_landclass <- project(resistance_landclass, dem)

### stack, crop to AOI -----------------------------------------

r_stack <- c(pop_weighted_proxity_town_norm, 
             remains, 
             IBAs,
             protected_BAs,
             resistance_landclass)

# Load regions
regions <- vect("data/boundaries/is_50v_mork_epsg_8088.gpkg/is_50v_mork_epsg_8088.gpkg", 
                layer = "mork_logregluumdaemi")

# Reproject to DEM CRS
regions_proj <- project(regions, crs(dem))

# Subset the region you want (2) reyjkanes
rykb_subset <- regions_proj[regions_proj$objectid == 2, ]


# crop stack
r_stack_cropped <- crop(r_stack, rykb_subset)
r_stack_cropped <- c(SDM, road_dist_norm, r_stack_cropped)
r_stack_cropped <- mask(r_stack_cropped, rykb_subset)

names(r_stack_cropped) <- c("Habitat preference", 
                            "Site acscessibility", 
                            "Pref. dist. from. cities", 
                            "Archaeological remains",
                            "Important bird areas",
                            "Protected bird areas",
                            "fresh water, glaciers, and infrastructure")

### factors to be weighted

factors_to_mask <- c(
  r_stack_cropped[[ which(names(r_stack_cropped) == "Archaeological remains") ]],
  r_stack_cropped[[ which(names(r_stack_cropped) == "Protected bird areas") ]],
  r_stack_cropped[[ which(names(r_stack_cropped) == "fresh water, glaciers, and infrastructure") ]]
)

composite_mask <- sum(factors_to_mask)
composite_mask <-ifel(composite_mask > 0, 1, 0)

### factors to be weighted

stack_to_weight <- c(
  r_stack_cropped[[ which(names(r_stack_cropped) == "Habitat preference") ]],
  r_stack_cropped[[ which(names(r_stack_cropped) == "Site acscessibility") ]],
  r_stack_cropped[[ which(names(r_stack_cropped) == "Pref. dist. from. cities") ]],
  r_stack_cropped[[ which(names(r_stack_cropped) == "Important bird areas") ]]
)

### weights

w <- c(0.6, # SDM
       0.2, # accessibility
       0.1, # distance from cities
       0.1) # IBA

### weights
weighted_stack <- stack_to_weight

weighted_layers <- lapply(1:nlyr(stack_to_weight), function(i) {
  stack_to_weight[[i]] * w[i]
})

weighted_stack <- rast(weighted_layers)

suitability <- sum(weighted_stack)
plot(suitability)

masked_suitability <- ifel(SDM < 0.25, 0, suitability)
masked_suitability <- mask(masked_suitability, composite_mask, maskvalues = 1)


writeRaster(masked_suitability, "data/test5.tif")


