
library(terra)


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


###  Population-weighted Proximity Preference-----------------------------------------

pop_weighted_proxity_town <- rast("data/pop_weighted_proxity_town.tif")

min_val <- minmax(pop_weighted_proxity_town)[1]
max_val <- minmax(pop_weighted_proxity_town)[2]
pop_weighted_proxity_town_norm <- (pop_weighted_proxity_town - min_val) / (max_val - min_val)



###  Archaeological remains buffer -----------------------------------------

remains <- rast("data/minjastofnun/minjastofnun_raster.tif")


###  bird areas -----------------------------------------
 
IBAs <- rast("data/bird areas/important_bird_areas.tif") ### soft constraint (penalise)
protected_BAs <- rast("data/bird areas/protected_bird_areas.tif") ### hard constraint (exlcude)+



### stack, crop to AOI -----------------------------------------

r_stack <- c(pop_weighted_proxity_town_norm, 
             remains, 
             IBAs,
             protected_BAs)

# Load regions
regions <- vect("data/boundaries/is_50v_mork_epsg_8088.gpkg/is_50v_mork_epsg_8088.gpkg", 
                layer = "mork_logregluumdaemi")

dem <- rast("data/IslandsDEMv1.0_10x10m_isn2016_zmasl.tif")

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
                            "Protected bird areas")






