library(terra)
library(ggplot2)
library(dplyr)

options(scipen = 999)

access <- rast("data/site_accessibility/iceland_accessibility_10m.tif")

cult_forest <- vect("data/forests/raektad_skoglendi_Island_07052025.shp")
cult_forest_resamp <- project(cult_forest, access)

forest_access_values <- extract(access, cult_forest_resamp, fun = mean, na.rm = TRUE)


# 2. Rename columns for easier plotting (assuming 'access' is the first layer)
colnames(forest_access_values) <- c("ID", "accessibility")

# 3. Clean NAs (pixels outside the raster but inside polygon buffers)
forest_access_values <- na.omit(forest_access_values)

forest_access_values <- forest_access_values %>%
  filter(accessibility < 250) 


# 4. Create the Density Plot
ggplot(forest_access_values, aes(x = accessibility)) +
  geom_density(fill = "#2ca25f", color = "#006d2c") +
  xlim(0,5) +
  theme_minimal() +
  labs(
    title = "Accessibility Profile of Cultivated Forests",
    subtitle = "Distribution of 10m accessibility scores across all polygons",
    x = "Accessibility (Total Accumulation)",
    y = "Density"
  ) 