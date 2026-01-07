library(geodata)
library(terra)

# 1. Define the Epic Parameters
models <- c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2", "EC-Earth3-Veg", 
            "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", 
            "MRI-ESM2-0", "UKESM1-0-LL")

# WorldClim CMIP6 standard SSPs: 126 (Green), 245 (Middle), 370 (High), 585 (Fossil-Fuel)
ssps <- c("126", "245", "370", "585")

# All available time windows
periods <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")

# Your specific path
my_path <- "data/climate_change_projections"

# 2. Execute the Loop
# Total downloads: 10 models * 4 SSPs * 4 periods = 160 tiles (~12 GB total)
for (m in models) {
  for (s in ssps) {
    for (p in periods) {
      
      message(sprintf("Processing: %s | SSP %s | Period %s", m, s, p))
      
      tryCatch({
        # Download 30s tile for Iceland center (Lon -18, Lat 65)
        # Note: If file exists, it will skip automatically.
        cmip6_tile(
          lon = -18, 
          lat = 65,
          model = m,
          ssp = s,
          time = p,
          var = "bioc", # Just Bioclimatic variables
          res = 0.5,
          path = my_path
        )
      }, error = function(e) {
        message(sprintf("--- Error downloading %s %s %s: %s", m, s, p, e$message))
      })
      
    }
  }
}