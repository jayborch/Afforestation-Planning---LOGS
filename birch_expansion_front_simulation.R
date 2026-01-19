# =========================================================================================
# PROJECT: Multi-Scale National Afforestation Framework for Iceland
# SCRIPT: Spatio-Temporal Cellular Automaton (CA) for Downscaled Birch Migration
# 
# METHODOLOGY OVERVIEW:
# This script executes a "Mechanistic-Stochastic" simulation of Betula pubescens (Birch) 
# expansion across Iceland from 2020 to 2100. It bridges the gap between theoretical 
# climate suitability and biological dispersal reality.
#
# THE THREE PILLARS OF THIS MODEL:
# 1. BAYESIAN SDM: Uses posterior draws from a BRMS model to account for parameter 
#    uncertainty. Suitability isn't a fixed value, but a range of ecological possibilities.
# 2. CLIMATE ENSEMBLE: Processes 2,000 unique futures (10 GCMs x 4 SSPs x 50 Draws). 
#    This captures the full spectrum of the CMIP6 climate forecasts.
# 3. BIOLOGICAL FILTER: A Cellular Automaton that prevents "teleportation." Expansion 
#    is constrained by maturity lags, dispersal limits, and age-dependent mortality.
#
# KEY ECOLOGICAL LOGIC:
# - DISPERSAL: 1km per 5-year step (adjacent cells only).
# - MATURITY: 15-year lag required before a colonized pixel can act as a seed source.
# - ESTABLISHMENT: Requires high suitability (>0.6) and is influenced by "Seed Pressure"
#   (number of neighbors) and "Frontier Fatigue" (penalizing isolated pioneers).
# - PERSISTENCE: Asymmetric mortality. Saplings (<15yr) die at <0.5 suitability; 
#   Adults (>=15yr) endure sub-optimal conditions down to <0.2 suitability.
#
# SCALE: 1km resolution (Strategic Wave) to be used as a constraint for 10m (Tactical Site) planning.
# =========================================================================================
# =========================================================================
# LIBRARIES
# =========================================================================
# =========================================================================
# LIBRARIES
# =========================================================================
library(terra)
library(dplyr)
library(tidyr)
library(RANN)
library(leaflet)

# =========================================================================
# READ AND PROCESS WIND DATA
# =========================================================================
wind_df <- readRDS("data/wind/vedurstofa_wind_rose.rds")

# Precompute wind weights per direction (normalized)
p <- 2.5
wind_df <- wind_df %>%
  mutate(wind_weight = frequency * (weibull_A^p) * gamma(1 + p / weibull_k)) %>%
  group_by(id) %>%
  mutate(wind_weight = wind_weight / sum(wind_weight)) %>%
  ungroup()

directions <- seq(0, 330, by = 30)

wind_lookup <- wind_df %>%
  select(id, direction, wind_weight) %>%
  pivot_wider(names_from = direction, values_from = wind_weight, values_fill = 0)

wind_matrix_raw <- as.matrix(wind_lookup[,-1])
rownames(wind_matrix_raw) <- wind_lookup$id
colnames(wind_matrix_raw) <- as.character(directions)

# Convert meteorological (FROM) → dispersal (TO)
downwind_dirs <- (as.numeric(colnames(wind_matrix_raw)) + 180) %% 360
downwind_dirs[downwind_dirs == 360] <- 0
wind_matrix <- wind_matrix_raw
colnames(wind_matrix) <- as.character(downwind_dirs)

# =========================================================================
# PREVAILING WIND (DETERMINISTIC)
# =========================================================================
dominant_wind <- apply(wind_matrix, 1, function(x) {
  as.numeric(names(which.max(x)))
})

WIND_GATE_WIDTH <- 60  # degrees (±30°)

# =========================================================================
# SETTINGS & PATHS
# =========================================================================
DRAWS_TO_RUN <- 1
path <- "models/future_projections_brms_all_draws"
out_dir <- "models/sim_results"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

gcms <- c("ACCESS-CM2","BCC-CSM2-MR","CMCC-ESM2","EC-Earth3-Veg",
          "INM-CM5-0","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR",
          "MRI-ESM2-0","UKESM1-0-LL")
ssps <- c("ssp126","ssp245","ssp370","ssp585")

all_files <- list.files(path, pattern = "\\.tif$")
b0_all <- rast(file.path(path, "birch_posterior_draws_baseline.tif"))
n_cells <- ncell(b0_all)

# =========================================================================
# PRECOMPUTE NEIGHBORS & WIND MAPPING
# =========================================================================
cat("Precomputing spatial indices...\n")
all_adj <- adjacent(b0_all[[1]], cells = 1:n_cells, directions = 8, pairs = TRUE)
cell_coords <- xyFromCell(b0_all[[1]], 1:n_cells)

# Neighbor sectors
from_coords <- cell_coords[all_adj[,1], ]
to_coords   <- cell_coords[all_adj[,2], ]

dx <- to_coords[,1] - from_coords[,1]
dy <- to_coords[,2] - from_coords[,2]

sector <- (round(atan2(dy, dx) * 180/pi / 30) * 30) %% 360
sector[sector == 360] <- 0

neighbors_df <- data.frame(
  source = all_adj[,1],
  target = all_adj[,2],
  sector = sector
)

neighbors_list <- split(neighbors_df, neighbors_df$source)

# Nearest wind station per cell
wind_coords <- wind_df %>% distinct(id, lon, lat)
nn <- nn2(wind_coords[, c("lon","lat")], cell_coords, k = 1)
cell_wind_ids <- wind_coords$id[nn$nn.idx[,1]]

# =========================================================================
# MASTER ENSEMBLE LOOP
# =========================================================================
for(g in gcms){
  for(s in ssps){
    
    target_files <- all_files[grep(g, all_files)]
    target_files <- target_files[grep(s, target_files)]
    if(length(target_files) < 4) next
    
    t1_all <- rast(file.path(path, target_files[grep("2021-2040", target_files)]))
    t2_all <- rast(file.path(path, target_files[grep("2041-2060", target_files)]))
    t3_all <- rast(file.path(path, target_files[grep("2061-2080", target_files)]))
    t4_all <- rast(file.path(path, target_files[grep("2081-2100", target_files)]))
    
    for(d in 1:DRAWS_TO_RUN){
      
      cat(sprintf("\n>>> GCM: %s | SSP: %s | Draw: %d <<<\n", g, s, d))
      
      b0 <- b0_all[[d]]
      keyframes <- list(
        list(b0, t1_all[[d]]),
        list(t1_all[[d]], t2_all[[d]]),
        list(t2_all[[d]], t3_all[[d]]),
        list(t3_all[[d]], t4_all[[d]])
      )
      
      occupied <- which(values(b0, mat = FALSE) > 0.9)
      arrival_map <- b0 * NA
      arrival_map[occupied] <- 2020
      
      for(win in 1:4){
        for(step in 1:4){
          
          current_year <- 2020 + (win-1)*20 + step*5
          weight <- step * 0.25
          current_suit <- keyframes[[win]][[1]] * (1-weight) +
            keyframes[[win]][[2]] * weight
          suit_vals <- values(current_suit, mat = FALSE)
          
          # -------------------------
          # MORTALITY
          # -------------------------
          if(length(occupied) > 0){
            age <- current_year - values(arrival_map)[occupied]
            occ_suit <- suit_vals[occupied]
            
            dead <- (age < 15 & occ_suit < 0.5) |
              (age >= 15 & occ_suit < 0.2)
            
            if(any(dead)){
              arrival_map[occupied[dead]] <- NA
              occupied <- occupied[!dead]
            }
          }
          
          # -------------------------
          # WIND-GATED DISPERSAL
          # -------------------------
          arrival_vals <- values(arrival_map, mat = FALSE)
          mature <- occupied[(current_year - arrival_vals[occupied]) >= 15]
          
          if(length(mature) > 0){
            new_cells <- integer(0)
            
            for(src in mature){
              
              wid <- cell_wind_ids[src]
              dom_dir <- dominant_wind[wid]
              
              nbs <- neighbors_list[[as.character(src)]]
              if(is.null(nbs)) next
              
              targets <- nbs$target
              valid <- !targets %in% occupied & !is.na(suit_vals[targets])
              if(!any(valid)) next
              
              v_targets <- targets[valid]
              v_sectors <- nbs$sector[valid]
              
              # wind gating
              ang_diff <- abs(v_sectors - dom_dir)
              ang_diff <- pmin(ang_diff, 360 - ang_diff)
              
              downwind <- ang_diff <= (WIND_GATE_WIDTH / 2)
              if(!any(downwind)) next
              
              dw_targets <- v_targets[downwind]
              
              # seed production
              n_seeds <- rpois(1, lambda = 2)
              if(n_seeds > 0){
                landed <- sample(dw_targets, size = n_seeds, replace = TRUE)
                survived <- landed[runif(length(landed)) < suit_vals[landed]^2]
                if(length(survived) > 0){
                  new_cells <- c(new_cells, survived)
                }
              }
              
              # -------------------------
              # LONG-DISTANCE JUMP (ALIGNED)
              # -------------------------
              jump_prob <- max(wind_matrix[wid,], na.rm = TRUE) * 0.1
              if(runif(1) < jump_prob){
                
                jump_dir <- dom_dir + runif(1, -WIND_GATE_WIDTH/2, WIND_GATE_WIDTH/2)
                jump_dir <- jump_dir %% 360
                
                jump_dist <- rpois(1, lambda = 3) + 2
                dx_j <- round(jump_dist * sin(jump_dir * pi/180))
                dy_j <- round(jump_dist * cos(jump_dir * pi/180))
                
                tx <- colFromCell(b0, src) + dx_j
                ty <- rowFromCell(b0, src) + dy_j
                
                if(tx >= 1 && tx <= ncol(b0) && ty >= 1 && ty <= nrow(b0)){
                  tid <- cellFromRowCol(b0, ty, tx)
                  if(!tid %in% occupied & !is.na(suit_vals[tid])){
                    if(runif(1) < suit_vals[tid]){
                      new_cells <- c(new_cells, tid)
                    }
                  }
                }
              }
            }
            
            if(length(new_cells) > 0){
              new_cells <- unique(new_cells)
              arrival_map[new_cells] <- current_year
              occupied <- c(occupied, new_cells)
            }
          }
        }
      }
      
      fname <- sprintf("%s/arrival_%s_%s_draw%02d.tif", out_dir, g, s, d)
      writeRaster(arrival_map, fname, overwrite = TRUE,
                  gdal = c("COMPRESS=DEFLATE"))
      rm(arrival_map, b0); gc()
    }
    
    rm(t1_all, t2_all, t3_all, t4_all); gc()
  }
}




# ---------------------------------------------------------
# Visualization: Non-Averaged Consensus (Quantiles)
# ---------------------------------------------------------

# 1. Define paths and load the 2,000 files
out_dir <- "models/sim_results"
sim_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)

# Create a 'virtual' stack (memory efficient)
results_stack <- rast(sim_files)

# 2. CALCULATION 1: Probability of Presence (%)
# This remains a binary count: what % of our 2,000 "worlds" does birch survive in?
cat("Calculating Presence Probability...\n")
presence_count <- sum(!is.na(results_stack)) 
prob_presence <- (presence_count / nlyr(results_stack)) * 100

# 3. CALCULATION 2: The "Un-flattened" Arrival Years
# Instead of mean(), we use app() with quantiles to see the actual distribution.
cat("Calculating Arrival Quantiles (This takes longer than mean)...\n")

# Median (50th percentile): The "Most Likely" arrival year
arrival_median <- app(results_stack, fun = function(x) quantile(x, probs = 0.5, na.rm = TRUE))

# Safe Bet (90th percentile): The "High Certainty" arrival year 
# (By this year, 90% of simulations agree birch has arrived)
arrival_safe <- app(results_stack, fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE))


# 4. LEAFLET SETUP
# Palettes
pal_prob <- colorNumeric(palette = "viridis", domain = c(0, 100), na.color = "transparent")
pal_year <- colorNumeric(palette = "inferno", domain = c(2020, 2100), na.color = "transparent")

# 5. RENDER INTERACTIVE MAP
leaflet() %>%
  addProviderTiles(providers$OpenTopoMap, group = "Topography") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  
  # Layer 1: Probability of Presence (%)
  addRasterImage(prob_presence, 
                 colors = pal_prob, 
                 opacity = 0.6, 
                 group = "Probability of Presence") %>%
  
  # Layer 2: Median Arrival (Most Likely Future)
  addRasterImage(arrival_median, 
                 colors = pal_year, 
                 opacity = 0.6, 
                 group = "Median Arrival (50%)") %>%
  
  # Layer 3: Safe Arrival (Conservative Future)
  addRasterImage(arrival_safe, 
                 colors = pal_year, 
                 opacity = 0.6, 
                 group = "Safe Arrival (90%)") %>%
  
  # Legends
  addLegend(pal = pal_prob, values = c(0, 100),
            title = "Presence Prob (%)", 
            position = "bottomleft", group = "Probability of Presence") %>%
  
  addLegend(pal = pal_year, values = c(2020, 2100),
            title = "Arrival Year", 
            position = "bottomright") %>%
  
  # Toggle Controls
  addLayersControl(
    baseGroups = c("Topography", "Satellite"),
    overlayGroups = c("Probability of Presence", "Median Arrival (50%)", "Safe Arrival (90%)"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  hideGroup("Safe Arrival (90%)") # Hide by default to keep map clean

