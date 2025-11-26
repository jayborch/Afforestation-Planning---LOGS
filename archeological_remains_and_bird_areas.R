

### Archaeological remains buffer

# Monuments, immunity area (https://gatt.natt.is/geonetwork/srv/ice/catalog.search#/metadata/177ede70-2f02-4af2-953b-61f05c7a2a7f)
# This dataset contains information on the privacy of surveyed archaeological remains that has been submitted to the Icelandic National Heritage Institute.
# The dataset shows a margin around the locations of archaeological remains that is the estimated privacy, 15 meters for protected archaeological remains and 100 meters for protected archaeological remains. 
# As the privacy is only estimated, it does not meet the requirements for planning or construction purposes, but is for reference only.
# Please note that the dataset is not an exhaustive overview of archaeological remains, buildings and structures in Iceland, which are protected under the Cultural Heritage Act No. 80/2012.

library(terra)
library(sf)

# Path to GeoPackage
gpkg_path <- "data/minjastofnun/fornleifaskraning_fridhelgi.gpkg"

# List layers (using terra + sf approach)
st_layers(gpkg_path)

# Suppose the layer you want is "fridhelgi_polygons"
vec <- vect(gpkg_path, layer = "fornleifaskraning_buffer")

# rasterize to dem
dem <- rast("data/IslandsDEMv1.0_10x10m_isn2016_zmasl.tif")
vec_proj <- project(vec, crs(dem))
r_vec <- rasterize(vec_proj, dem, field = 1, background = 0)

# Write raster
writeRaster(
  r_vec, 
  filename = "data/minjastofnun/minjastofnun_raster.tif", 
  overwrite = TRUE
)

### Protected Bird Areas

# https://gatt.lmi.is/geonetwork/srv/eng/catalog.search#/metadata/c65e5e41-8b99-4371-9c5d-3bbb9c084bb7

# A total of 121 IBAs are defined. 70 IBAs are seabird colonies (SF), 25 IBAs are primarily intertidal zones and adjacent shallow marine waters (FG), 
# 31 IBAs are inland, predominantly fertile wetlands and surface inland waters. A few IBAs fall under two or three categories. These IBAs are located in 65 of the 74 municipalities of Iceland.
# The number of IBAs designated for each species varies; by far, most breeding sites of Fulmarus glacialis (38). No areas were designated for 40 species. For more than half of those IBA criteria have not been defined; o
# thers do not meet the criteria due to their small populations in Iceland. The feature attribute codes for 'stadaFitju' are:
#   0 = not previously designated as IBAs;
#   1 = previously designated as IBAs, outline unchanged;
#   2 = previously designated as IBAs, outline changed;
#   3 = previously designated as IBAs, outline drawn by IINH;
#   4 = protected or on the Register of areas of conservation interest (could as well be previously designated as IBAs), outline unchanged;
#   5 = protected or on the Register of areas of conservation interest (could as well be previously designated as IBAs), outline changed;
#   6 = .area is part of an area that is protected or on the Register of areas of conservation interest (could as well be previously designated as IBAs) or part of this area is protected or on the Register of areas of conservation interest (could as well be previously designated as IBAs), outline drawn by IINH.]

### there are 2 sets of data. One is hard constraint protected areas, while the second is more of a soft constraint

protected_areas <- vect("data/bird areas/natt_f25v_mikilvaegfuglasvaedi_epsg_3057/NI_F25v_mikilvaegFuglasvaedi_1.utg/NI_F25v_mikilvaegFuglasvaedi_SHP/ni_f25v_mikilvaegFuglasvaedi_fl.shp")
protected_areas_proj <- project(protected_areas, crs(dem))
PAs <- protected_areas_proj[protected_areas_proj$stadaFitju %in% c(4,5,6), ]
r_PAs<- rasterize(PAs, dem, field = 1, background = 0)

# Write raster
writeRaster(
  r_PAs, 
  filename = "data/bird areas/protected_bird_areas.tif", 
  overwrite = TRUE
)

IBA <- protected_areas_proj [protected_areas_proj $stadaFitju %in% c(1,2,3), ]
r_IBA <- rasterize(IBA, dem, field = 1, background = 0)

# Write raster
writeRaster(
  r_IBA, 
  filename = "data/bird areas/important_bird_areas.tif", 
  overwrite = TRUE
)




