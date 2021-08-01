################################################################################
################################################################################
# PREDICTION MAP
################################################################################
################################################################################

# load packages
library(raster)
library(sf)
library(sp)
library(keras)
library(meteo) # for tiling

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/06_prepare_cnn_input.R")

# set paths
basedir <- "H:/Daten/Studium/2_Master/4_Semester"
path_rasters <- paste0(basedir, "/4_Daten/rasters_2cm") # input
path_model   <- paste0(basedir, "/4_Daten/models_2cm/tls_rgb_geo/bla.hd5") # input
path_area    <- paste0(basedir, "/4_Daten/sites/convex/area_polygons.shp") # input
path_values  <- paste0(basedir, "/4_Daten/clips_2cm/raster_values_labels_unscaled.csv") # input
path_labels  <- paste0(basedir, "/4_Daten/model_input_2cm_standardized/tls/label_lookup.csv") # input
path_plots   <- paste0(basedir, "/4_Semester/5_Analyse/Plots") # output
path_tiles   <- paste0(path_rasters, "/tiles") # output

# set parameters
n_pixels <- 25
n_bands <- 13
area <- 6 # because it has all classes (1,6,7,8)
crs_raster_las <- "+proj=utm +zone=32 +ellps=WGS84 +units=m +vunits=m +no_defs"
use_these <- c(
  "ortho", "anisotropy_max", "curvature_max", "linearity_max",
  "linearity_sd", "nDSM", "planarity_mean", "planarity_sd",
  "point_density", "reflectance_mean", "reflectance_sd"
)



# use custom color palette, bright
own_colors_named <- list(
  red = "#ff6f69",
  blue = "#00aedb",
  yellow = "#ffcc5c",
  green = "#88d8b0",
  turquoise = "#4abdac",
  pink = "#ff8b94",
  bright_green = "#cbe885"
)
own_colors <- c(
  own_colors_named$blue, own_colors_named$red, own_colors_named$yellow,
  own_colors_named$bright_green, own_colors_named$green
)

################################################################################
# TILE RASTER
################################################################################

# get all rasters of the area
raster_paths <- list.files(path_rasters, pattern = paste0("area_", area), recursive = T, full.names = T)
raster_paths <- raster_paths[grepl(pattern = "[.]tif", raster_paths)]
raster_paths <- raster_paths[!grepl(pattern = "nDSM_filtering", raster_paths)]

# load areas
area_polys <- st_read(path_area)
st_crs(area_polys) <- CRS("+init=EPSG:25832")
area_polys <- st_transform(area_polys, crs_raster_las)
area_poly <- area_polys[area_polys$ID == area, ]

# get lookup table for normalization
lookup_table <- rescale_values(path_values)

# get raster files which are within the selection
types <- c()
for (i in 1:length(raster_paths)) {
  types[i] <- strsplit(basename(raster_paths[i]), "_area_")[[1]][1]
}
raster_paths <- raster_paths[types %in% use_these]
types <- types[types %in% use_these]

# empty object for rasters
raster_list <- list()

# for every raster
for (i in 1:length(raster_paths)) {
  # load raster & type
  type <- types[i]
  path <- raster_paths[i]
  raster_list[[i]] <- stack(path)

  # clip to area
  raster_list[[i]] <- mask(raster_list[[i]], area_poly) # want too big but rectangle plot? use crop()

  bands <- if (type != "ortho") c(type) else c("R", "G", "B")
  for (j in 1:length(bands)) {
    # standardize using mean / sd from all data
    mean_val <- lookup_table[[paste0(bands[j], "_mean")]]
    sd_val <- lookup_table[[paste0(bands[j], "_sd")]]
    raster_list[[i]][[j]] <- scale(raster_list[[i]][[j]], center = mean_val, scale = sd_val)
  }
}

# stack rasters
raster_stack <- stack(raster_list)
rm(raster_list)

# make tiles
agg <- aggregate(raster_stack[[1]], n_pixels)
polys <- as(agg, "SpatialPolygons")
tiles <- lapply(seq_along(polys), function(i) crop(raster_stack, polys[i]))

# save tiles
saveRDS(tiles, file = paste0(path_tiles, ".rds"))

################################################################################
# MAKE PREDICTIONS
################################################################################

# load tiles as arrays
tiles <- readRDS(paste0(path_tiles, ".rds"))
img_array <- array(dim = c(length(tiles), n_pixels, n_pixels, n_bands))
for (i in 1:length(tiles)) {
  img_array[i, , , ] <- values(tiles[[i]])
}

# replace missing values with 0
img_array[is.na(img_array)] <- 0

# load model
model <- load_model_hdf5(path_model)

# predict
preds <- model %>% predict(img_array)
preds <- apply(preds, 1, which.max)

# replace numeric label with text label
label_lookup <- read.csv(path_labels)
for (i in 1:nrow(label_lookup)) {
  preds[preds == label_lookup$new[i]] <- label_lookup$old[[i]]
}

# save as shapefile
polys <- lapply(tiles, function(x) as(extent(x), "SpatialPolygons"))
polys <- do.call(bind, polys)
polys <- SpatialPolygonsDataFrame(Sr = polys, data = data.frame(
  id = 1:length(polys),
  prediction = preds
))
crs(polys) <- CRS(crs_raster_las)
shapefile(polys, paste0(path_tiles, ".shp"), overwrite = TRUE)

################################################################################
# MAP
################################################################################

# load packages
library(ggmap)
library(rgdal)
library(ggplot2)
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)

# load & prepare data
shapefile <- readOGR(paste0(path_tiles, ".shp"))
shapefile[is.na(shapefile$prediction), ]
df <- tidy(shapefile, region = "id")
shapefile$id <- rownames(shapefile@data)
df <- left_join(df, shapefile@data, by = "id")
df <- df[!is.na(df$prediction), ]

# make map
cairo_pdf(file = paste0(path_plots, "/prediction_map.pdf"), family = "Calibri", width = 8.27, height = 5.83)
ggplot() +
  geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = prediction), color = NA) +
  scale_fill_manual(
    values = own_colors, "Predicted\nGround Cover Class",
    labels = c("Blueberry", "Deadwood", "Forest Floor", "Moss", "Spruce")
  ) +
  theme_light() +
  coord_equal() +
  ylim(c(5379630, 5379665)) +
  xlim(c(447025, 447075)) +
  ylab("Northing\n") +
  xlab("\nEasting") +
  theme(text = element_text(family = "Calibri", size = 16))
dev.off()

################################################################################