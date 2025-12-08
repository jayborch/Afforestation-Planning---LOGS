library(tools)

# Set your root directory
root <- "C:/Users/jay.philip/Desktop/Afforestation-Planning---LOGS/data/GEE/east_fjords_tiles"   # <-- change this

# List all zip files
zips <- list.files(root, pattern = "\\.zip$", full.names = TRUE)

# Unzip each zip into its own folder
for (z in zips) {
  unzip(z, exdir = file.path(root, tools::file_path_sans_ext(basename(z))))
}

# Now collect all images (recursive = TRUE searches inside all folders)
image_patterns <- c("*.tif", "*.jpg", "*.jpeg", "*.png")

image_files <- unlist(
  lapply(image_patterns, function(p) {
    list.files(root, pattern = glob2rx(p), full.names = TRUE, recursive = TRUE)
  })
)

# Remove images that are already in the root
image_files <- image_files[dirname(image_files) != root]

# Move images to root
for (img in image_files) {
  file.rename(img, file.path(root, basename(img)))
}