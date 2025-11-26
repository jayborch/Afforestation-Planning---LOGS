
# This script models a population-weighted accessibility or "preference" surface around towns. 
# It assumes that influence or preference of a town decays exponentially with distance. 
# Steps:
# 1. Plots example exponential decay curves for different town populations.
# 2. Loads a DEM and a shapefile of towns with population data.
# 3. For each town, computes a raster where preference decreases exponentially with distance from the town, weighted by population.
# 4. Sums the contributions from all towns to produce a cumulative preference raster.
# 5. Resamples the result to match the original DEM resolution and masks town locations to zero.
# 6. Saves the final raster as a GeoTIFF for further spatial analysis.

library(terra)

# Parameters

lambda <- 0.001
distance <- seq(0, 10000, length.out = 500)  # distance in meters
pop_vals <- c(1000, 2000, 4000)  # example populations
colors <- c("darkgreen", "dodgerblue", "firebrick")

# Plot setup

plot(distance, rep(0, length(distance)), type = "n",
     ylim = c(0, max(pop_vals)),
     xlab = "Distance from town (m)",
     ylab = "Preference",
     main = expression(paste("Preference = Population × ", e^{-lambda %.% distance})))

# Add curves for each population

for (i in seq_along(pop_vals)) {
  pref <- pop_vals[i] * exp(-lambda * distance)
  lines(distance, pref, col = colors[i], lwd = 2)
}

# Add legend

legend("topright", legend = paste("Population =", pop_vals),
       col = colors, lwd = 2, bty = "n")



# --- 1. Load base raster and towns ---
dem <- rast("data/IslandsDEMv1.0_10x10m_isn2016_zmasl.tif")
dem_lowres <- aggregate(dem, fact = 10, fun = mean)  # Lower resolution = less memory

towns <- vect("data/towns/roads.shp")
towns <- project(towns, crs(dem_lowres))

# Extract population column
pop <- towns$ibuarfj_24
pop[is.na(pop)] <- 0

# --- 2. Create empty raster to accumulate preference ---
pref_total <- rast(dem_lowres)
values(pref_total) <- 0  # initialize with zeros

lambda <- 0.001

# --- 3. Loop through each town, add weighted preference to total ---
for (i in 1:nrow(towns)) {
  cat("Processing town", i, "with population", pop[i], "\n")
  
  # Skip towns with 0 population
  if (pop[i] == 0) next
  
  # Rasterize current town
  town_r <- rasterize(towns[i, ], dem_lowres, field = 1)
  
  # Compute distance
  dist_r <- distance(town_r)
  
  # Weighted preference
  pref <- pop[i] * exp(-lambda * dist_r)
  
  # Add to total
  pref_total <- pref_total + pref
  
  # Clean up to save memory
  rm(pref, dist_r, town_r)
  gc()
}

# resample to DEM and then mask towns as zero
pref_total_resampled <- terra::resample(pref_total, dem, method = "bilinear")
town_mask <- rasterize(towns, dem, field=1, background=NA)
pref_total_resampled_masked <- ifel(!is.na(town_mask), 0, pref_total_resampled)


writeRaster(
  pref_total_resampled_masked,
  filename = "data/pop_weighted_proxity_town.tif",
  overwrite = TRUE
)

