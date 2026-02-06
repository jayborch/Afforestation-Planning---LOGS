# ===========================================
# WIND EROSION PREPARATION SCRIPT (ICELAND)
# 1 km resolution, IDW-mapped wind roses
# ===========================================

library(terra)
library(sf)
library(dplyr)
library(purrr)
library(furrr)

plan(multisession)  # parallel processing
SDM_raw   <- rast("data/broad_scale_rf_prob.tif")
SDM       <- project(SDM_raw, "EPSG:8088", res = 1000)

# ==============================================================================
# PHASE 2: SOIL PARAMETER DERIVATION (REVISED FOR ICELANDIC CONDITIONS)
# Purpose: Calculate erodibility while preventing "Zero-out" errors.
# ==============================================================================

# Load the soil property stack
soils <- rast("data/soils/iceland_soilgrid_0-5cm.tif")

# ------------------------------------------------------------------------------
# 1. UNIT CONVERSION (SoilGrids Integer -> Percentage)
# ------------------------------------------------------------------------------
Sa  <- soils[["sand"]] / 10
Si  <- soils[["silt"]] / 10
Cl  <- soils[["clay"]] / 10
soc_pct <- soils[["soc"]] / 100
OM      <- soc_pct * 1.724 

# ------------------------------------------------------------------------------
# 2. ERODIBLE FRACTION (EF) - Adjusted for Icelandic Andosols
# ------------------------------------------------------------------------------
# We reduce the OM penalty from 2.59 to 1.0 to reflect the high erodibility 
# of volcanic ash even when organic matter is present.
EF_raw <- (29.09 + 0.31*Sa + 0.17*Si + 0.33*(Sa/(Cl + 0.1)) - 1.0*OM) / 100
EF_rast <- clamp(EF_raw, lower = 0.01, upper = 1)

# ------------------------------------------------------------------------------
# 3. SOIL CRUST FACTOR (SCF) - Linearized for Volcanic Soils
# ------------------------------------------------------------------------------
# Original RWEQ: 1 / (1 + 0.0066 * (Cl^2) + 0.021 * (OM^2))
# Problem: Squared terms drive SCF to ~0 in Iceland's carbon-rich soils.
# Fix: Remove squares to allow for a more realistic, non-zero crust factor.
SCF_raw  <- 1 / (1 + 0.0066 * Cl + 0.021 * OM)

# Clamp at 0.05 to ensure even the most stable soil has some erosion potential.
SCF_rast <- clamp(SCF_raw, lower = 0.05, upper = 1)

# ------------------------------------------------------------------------------
# 4. SPATIAL ALIGNMENT & MASTER SNAP (Scaling to 1km SDM)
# ------------------------------------------------------------------------------
cat("Aligning and snapping soil factors to 1km SDM grid...\n")

# A. Project and Resample EF
# We use 'bilinear' for smooth transitions and snap it to the SDM template
EF_1km <- project(EF_rast, SDM)
SCF_1km <- project(SCF_rast,SDM)

# Clean up names
names(EF_1km)  <- "Erodible_Fraction"
names(SCF_1km) <- "Soil_Crust_Factor"

# ------------------------------------------------------------------------------
# 3. LOAD WIND DATA & TRANSFORM
# ------------------------------------------------------------------------------
cat("Loading and projecting wind stations...\n")
wind_long <- readRDS("data/wind/vedurstofa_wind_rose.rds")

# We use the SDM as the master CRS/Grid template
wind_stations <- wind_long %>%
  distinct(id, lon, lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs(SDM)) 

station_coords <- st_coordinates(wind_stations)

# ------------------------------------------------------------------------------
# 4. PREP COORDINATES (Updated to SDM 1km Grid)
# ------------------------------------------------------------------------------
# Extract XY for every cell in the 1km SDM template
all_coords <- xyFromCell(SDM, 1:ncell(SDM))

# Identify valid cells (Where SDM isn't NA) to save computation time
valid_cells <- which(!is.na(values(SDM)))

# ------------------------------------------------------------------------------
# 5. WIND FACTOR FUNCTION (RWEQ) - No change needed to the physics
# ------------------------------------------------------------------------------
u_t <- 7  # Threshold wind speed (m/s)

wind_factor_vec <- function(A, k, freq, u_t){
  ifelse(A*4 <= u_t, 0, {
    u <- seq(u_t, A*4, by = 0.1)
    pdf <- (k / A) * (u / A)^(k - 1) * exp(-(u / A)^k)
    integral <- sum((u - u_t)^2 * pdf * 0.1)
    freq * integral
  })
}

# ------------------------------------------------------------------------------
# 6. IDW BLENDING FUNCTION - No change needed
# ------------------------------------------------------------------------------
blend_wind_to_cell <- function(cell_x, cell_y, wind_long, station_coords, n_nearest = 3, fudge = 100){
  dists <- sqrt((station_coords[,1] - cell_x)^2 + (station_coords[,2] - cell_y)^2)
  nearest_idx <- order(dists)[1:n_nearest]
  weights <- 1 / (dists[nearest_idx] + fudge)^2
  w_norm <- weights / sum(weights)
  
  nearest_ids <- wind_stations$id[nearest_idx]
  
  blended <- wind_long %>%
    filter(id %in% nearest_ids) %>%
    group_by(direction) %>%
    summarise(
      frequency = sum(frequency * w_norm[match(id, nearest_ids)]),
      weibull_A = sum(weibull_A * w_norm[match(id, nearest_ids)]),
      weibull_k = sum(weibull_k * w_norm[match(id, nearest_ids)]),
      .groups = "drop"
    )
  return(blended)
}

# ------------------------------------------------------------------------------
# 7. COMPUTE V FOR EACH VALID CELL
# ------------------------------------------------------------------------------
cat("Computing Wind Energy Integral (V) for valid cells...\n")

# parallel processing with furrr
V_values <- future_map_dbl(valid_cells, function(cell_idx){
  xy <- all_coords[cell_idx, ]
  blended <- blend_wind_to_cell(xy[1], xy[2], wind_long, station_coords)
  
  V_cell <- sum(pmap_dbl(list(blended$weibull_A, blended$weibull_k, blended$frequency),
                         ~ wind_factor_vec(..1, ..2, ..3, u_t)))
  return(V_cell)
}, .progress = TRUE)

# ------------------------------------------------------------------------------
# 8. CREATE FINAL V RASTER (Mapped to SDM)
# ------------------------------------------------------------------------------
V_rast <- rast(SDM) # Create empty raster from template
values(V_rast) <- NA
values(V_rast)[valid_cells] <- V_values
names(V_rast) <- "Wind_Energy_V"


# ==============================================================================
# PHASE 4: REFINED HABITAT PHYSICS (Vistgerð Integration)
# Purpose: Use detailed habitat types to set realistic Roughness and Cover.
# ==============================================================================

library(foreign) # Required for read.dbf

# 1. LOAD HABITAT DATA & ATTRIBUTES
habitats <- rast("data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif")
rat <- as.data.frame(read.dbf("data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif.vat.dbf"))

# Clean names for joining
rat_clean <- rat %>%
  dplyr::rename(Vistgerd = Vistgerð)

# ------------------------------------------------------------------------------
# 2. ASSIGN PHYSICS VIA KEYWORD DETECTION
# ------------------------------------------------------------------------------
# We map each Vistgerð to its physical 'Wind-Brake' properties.
rat_physics <- rat_clean %>%
  mutate(
    # Aerodynamic Roughness (z0 in meters)
    # Higher z0 = more wind drag/protection.
    z0 = case_when(
      grepl("Birkiskógur|Skógrækt", Vistgerd) ~ 0.500, # Tall canopy
      grepl("Mosa|Lyng|Gras|Móavist", Vistgerd) ~ 0.050, # Thick low vegetation
      grepl("Hraun", Vistgerd)                  ~ 0.080, # Jagged rock surface
      grepl("Melavist|Skriðuvist", Vistgerd)    ~ 0.010, # Sparse/patchy
      grepl("Sanda|Eyði|Aur|Malar|Mold", Vistgerd) ~ 0.001, # Smooth desert
      TRUE ~ 0.005 
    ),
    # Fractional Vegetation Cover (0 to 1)
    # This determines the COG (Cover on Ground) factor.
    fvc_base = case_when(
      grepl("Mosa|Gras|Lyng|Birkiskógur|Skógrækt|Móavist", Vistgerd) ~ 0.95,
      grepl("Melavist|Fléttuhraunavist", Vistgerd) ~ 0.20,
      grepl("Sanda|Eyði|Moldavist|Eyravist", Vistgerd) ~ 0.05,
      TRUE ~ 0.10
    )
  )

# ------------------------------------------------------------------------------
# 3. SPATIAL ALIGNMENT (1km Master Snap)
# ------------------------------------------------------------------------------
cat("Mapping Vistgerð physics to 1km grid...\n")

# A. Create z0 Raster
z0_mat  <- as.matrix(rat_physics[, c("Value", "z0")])
z0_rast <- classify(habitats, rcl = z0_mat)
z0_1km  <- resample(project(z0_rast, SDM), SDM, method = "bilinear")

# B. Create FVC Raster 
# We use aggregate to get the mean cover of the 5m cells within each 1km cell
fvc_mat  <- as.matrix(rat_physics[, c("Value", "fvc_base")])
fvc_rast <- classify(habitats, rcl = fvc_mat)
fvc_agg  <- aggregate(fvc_rast, fact = 200, fun = mean)
fvc_1km  <- resample(project(fvc_agg, SDM), SDM, method = "bilinear")

# ------------------------------------------------------------------------------
# 4. FINAL RWEQ ASSEMBLY
# ------------------------------------------------------------------------------
# K' (Roughness Factor)
K_prime <- 1 / log(10 / (z0_1km + 0.0001))

# COG (Vegetation Factor) - We use the RWEQ exponential decay
COG_val <- exp(-0.0438 * (fvc_1km * 100))

# Current Erosion Risk (Qmax)
# WF, EF, and SCF were calculated in previous sections
Qmax_current <- 10.97 * (V_rast * EF_1km * SCF_1km * K_prime * COG_val)

writeRaster(Qmax_current, "models/Qmax_erosion.tif", overwrite = T)

