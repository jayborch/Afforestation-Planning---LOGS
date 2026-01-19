# ==============================================================================
# CALIBRATED BIRCH RESTORATION MODEL: ICELAND (Tens of Meters Scale)
# ==============================================================================

library(terra)
library(dplyr)
library(ggplot2)
library(foreign)

# ------------------------------------------------------------------------------
# Terra performance options (important for 10 m national rasters)
# ------------------------------------------------------------------------------
terraOptions(
  memfrac = 0.8,
  progress = 1
)

# ------------------------------------------------------------------------------
# (1) VEGETATION QUALITY: Habitat Suitability Weights
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FULL BIRCH HABITAT WEIGHTS: ICELAND
# ------------------------------------------------------------------------------
birch_weights <- data.frame(
  Vistgerd = c(
    "L1.1 Eyðimelavist", "L1.2 Grasmelavist", "L1.3 Mosamelavist", "L1.4 Víðimelavist",
    "L1.5 Sanda- og vikravist", "L1.6 Landmelhólavist",
    "L3.1 Urðarskriðuvist", "L3.2 Grasvíðiskriðuvist", "L3.3 Ljónslappaskriðuvist",
    "L4.1 Eyravist", "L4.2 Auravist",
    "L5.1 Hélumosavist", "L5.2 Melagambravist", "L5.3 Hraungambravist",
    "L6.1 Eyðihraunavist", "L6.2 Fléttuhraunavist", "L6.3 Mosahraunavist", "L6.4 Lynghraunavist",
    "L7.1 Sandstrandarvist", "L7.2 Malarstrandarvist", "L7.3 Strandmelhólavist", "L7.4 Grashólavist",
    "L7.5 Sjávarfitjungsvist", "L7.6 Gulstararfitjavist", "L7.7 Sjávarkletta- og eyjavist",
    "L8.1 Dýjavist", "L8.2 Rekjuvist", "L8.3 Sandmýravist", "L8.4 Hrossanálarvist",
    "L8.5 Runnamýravist á hálendi", "L8.6 Runnamýravist á láglendi", "L8.7 Rimamýravist",
    "L8.8 Rústamýravist", "L8.9 Starungsmýravist", "L8.10 Hengistararflóavist", "L8.11 Brokflóavist",
    "L8.12 Starungsflóavist", "L8.13 Tjarnastararflóavist", "L8.14 Gulstararflóavist",
    "L9.1 Stinnastararvist", "L9.2 Finnungsvist", "L9.3 Bugðupuntsvist", "L9.4 Snarrótarvist",
    "L9.5 Grasengjavist", "L9.6 Língresisvist", "L9.7 Blómgresisvist",
    "L10.1 Mosamóavist", "L10.2 Flagmóavist", "L10.3 Starmóavist", "L10.4 Grasmóavist",
    "L10.5 Fléttumóavist", "L10.6 Fjalldrapamóavist", "L10.7 Lyngmóavist á hálendi", "L10.8 Lyngmóavist á láglendi",
    "L10.9 Víðimóavist", "L10.10 Víðikjarrvist", "L11 Birkiskógur",
    "L12.1 Mýrahveravist", "L12.2 Móahveravist", "L12.3 Fjallahveravist", "L12.4 Hveraleirsvist",
    "L13.1 Jöklar og urðarjöklar",
    "L14.1 Þéttbýli og annað manngert land", "L14.2 Tún og akurlendi", "L14.3 Skógrækt",
    "L14.4 Alaskalúpína", "L14.5 Uppgræðslur", "L14.6 Skógarkerfill og fleiri áþekkar tegundir",
    "F Fjöruvistir", "FX1.1 Sjávarlón", "V1 Vötn", "V2 Ár"
  ),
  Weight = c(
    0.0, 0.2, 0.2, 0.1, 0.0, 0.1,
    0.3, 0.4, 0.2,
    0.1, 0.3,
    0.3, 0.3, 0.2,
    0.0, 0.3, 0.3, 0.4,
    0.0, 0.2, 0.1, 0.5,
    0.6, 0.7, 0.3,
    0.3, 0.4, 0.3, 0.4,
    0.6, 0.7, 0.7,
    0.6, 0.7, 0.6, 0.7,
    0.7, 0.7, 0.8,
    0.5, 0.6, 0.7, 0.7,
    0.85, 0.85, 0.7,
    0.2, 0.2, 0.4, 0.75,
    0.3, 0.5, 0.5, 0.75,
    0.5, 0.9, 1.0,
    0.3, 0.3, 0.2, 0.1,
    0.0,
    0.0, 0.1, 0.8,
    0.2, 0.2, 0.2,
    0.0, 0.0, 0.0, 0.0
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# (2) LOAD DATA
# ------------------------------------------------------------------------------
habitats <- rast(
  "data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif"
)

# Aggregate to 10 m resolution using majority (modal) rule
habitats <- terra::aggregate(
  habitats,       # original raster
  fact = 2,          # factor: 2x2 cells -> 10 m (5m*2)
  fun = "modal",       # majority class
  na.rm = TRUE       # ignore NA cells
)


rat <- as.data.frame(
  read.dbf(
    "data/natt_vg25r_3_utg_epsg_3057/NI_VG25r_3.utg/ni_vg25r_3utg.tif.vat.dbf"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# (3) JOIN HABITAT WEIGHTS
# ------------------------------------------------------------------------------
rat_joined <- rat %>%
  dplyr::rename(Vistgerd = Vistgerð) %>%
  left_join(birch_weights, by = "Vistgerd") %>%
  mutate(
    Weight = ifelse(is.na(Weight), 0, Weight),
    Value  = as.integer(Value)
  )

# Lookup table for classify()
lut <- rat_joined %>%
  dplyr::select(Value, Weight) %>%
  distinct() %>%
  arrange(Value)

# ------------------------------------------------------------------------------
# (4) HABITAT SUITABILITY SURFACE (VERSION-SAFE)
# ------------------------------------------------------------------------------
suitability_layer <- classify(
  habitats,
  as.matrix(lut)
)

# ------------------------------------------------------------------------------
# (5) EXTRACT EROSION SOURCES
# ------------------------------------------------------------------------------
sand_ids <- rat_joined$Value[
  rat_joined$Vistgerd %in% c(
    "L1.5 Sanda- og vikravist",
    "L7.1 Sandstrandarvist"
  )
]

gravel_ids <- rat_joined$Value[
  rat_joined$Vistgerd %in% c(
    "L1.1 Eyðimelavist",
    "L4.2 Auravist"
  )
]

sand_mask   <- habitats %in% sand_ids
gravel_mask <- habitats %in% gravel_ids

# ------------------------------------------------------------------------------
# (6) BUFFER EROSION SOURCES (LIMIT COMPUTATION)
# ------------------------------------------------------------------------------
# Convert logical mask to NA for distance calculation
# TRUE = erosion source, FALSE = NA (distance will be computed from here)
sand_mask_na   <- ifel(sand_mask, 1, NA)
gravel_mask_na <- ifel(gravel_mask, 1, NA)

# Compute Euclidean distance from TRUE pixels (erosion sources)
dist_from_sand   <- distance(sand_mask_na)
dist_from_gravel <- distance(gravel_mask_na)

# ------------------------------------------------------------------------------
# (8) APPLY CALIBRATED DISTANCE DECAY
# ------------------------------------------------------------------------------
# STRATEGY 1: Mobile Sand (front-line stabilization)
pref_sand <- exp(-0.2 * dist_from_sand)

# STRATEGY 2: Gravel Flats (strategic buffer)
pref_gravel <- exp(-0.05 * dist_from_gravel)

# ------------------------------------------------------------------------------
# (9) COMBINE STRATEGIES (DOMINANT PROCESS)
# ------------------------------------------------------------------------------
# Just take the maximum of the two priority rasters
total_proximity_priority <- app(
  c(pref_sand, pref_gravel),
  fun = max,
  na.rm = TRUE
)

# ------------------------------------------------------------------------------
# (10) FINAL RESTORATION PRIORITY MAP
# ------------------------------------------------------------------------------
restoration_priority <- suitability_layer * total_proximity_priority

writeRaster(
  restoration_priority,
  "data/erosion_front/Iceland_Birch_Restoration_Priority.tif",
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# (11) VISUALIZATION: OPERATIONAL SCALE CURVES
# ------------------------------------------------------------------------------
dist_range <- seq(0, 150, 0.5)

viz_df <- data.frame(
  Distance = rep(dist_range, 2),
  Preference = c(
    exp(-0.2 * dist_range),
    exp(-0.05 * dist_range)
  ),
  Strategy = rep(
    c(
      "Mobile Sand (Tactical Edge)",
      "Gravel Flats (Strategic Buffer)"
    ),
    each = length(dist_range)
  )
)

ggplot(viz_df, aes(x = Distance, y = Preference, color = Strategy)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Calibrated Restoration Preference (Operational Scale)",
    subtitle = "Horizontal dashed line represents 50% priority threshold",
    x = "Distance from Erosion Source (meters)",
    y = "Priority Multiplier"
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"))
