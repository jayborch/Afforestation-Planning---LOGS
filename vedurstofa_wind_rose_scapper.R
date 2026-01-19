library(httr)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)

base_url <- "https://vindatlas.vedur.is/data/rav_3km_1995-2008_"

directions <- seq(0, 330, by = 30)

parse_lib <- function(text, id) {
  
  lines <- str_split(text, "\n")[[1]]
  lines <- str_trim(lines)
  
  # ---- Coordinates ----
  coord_line <- lines[str_detect(lines, "<coordinates>")]
  coords <- str_match(coord_line,
                      "<coordinates>([^,]+),([^,]+),([^<]+)</coordinates>")
  lon <- as.numeric(coords[2])
  lat <- as.numeric(coords[3])
  elev <- as.numeric(coords[4])
  
  # ---- Find the frequency line ----
  # Frequency is the FIRST line with 12 numeric entries
  numeric_lines <- lines[str_detect(lines, "^[-0-9\\. ]+$")]
  split_nums <- str_split(numeric_lines, "\\s+")
  
  lens <- map_int(split_nums, length)
  idx <- which(lens == 12)[1]
  
  freq <- as.numeric(split_nums[[idx]])
  A    <- as.numeric(split_nums[[idx + 1]])
  k    <- as.numeric(split_nums[[idx + 2]])
  
  tibble(
    id = id,
    lon = lon,
    lat = lat,
    elevation = elev,
    direction = directions,
    frequency = freq,
    weibull_A = A,
    weibull_k = k
  )
}

# ---- GRID LOOP ----
results <- list()

for (row in 1:200) {
  for (col in 1:200) {
    
    id <- sprintf("%03d%03d", row, col)
    url <- paste0(base_url, id, ".lib")
    
    resp <- GET(url)
    
    if (status_code(resp) != 200) next
    
    text <- content(resp, "text", encoding = "UTF-8")
    
    parsed <- tryCatch(
      parse_lib(text, id),
      error = function(e) NULL
    )
    
    if (!is.null(parsed)) {
      results[[id]] <- parsed
      cat("Parsed:", id, "\n")
    }
  }
}

wind_df <- bind_rows(results)

saveRDS(wind_df, "data/wind/vedurstofa_wind_rose.rds")