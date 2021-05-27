################################################################################
################################################################################
# COMBINING RASTER & VEGATATION
################################################################################
################################################################################

# load packages
library(lidR)
library(sf)

# set paths
path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"  # input
path_vegetation <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/vegetation/Export_ODK_clean_checked.kml"  # input
path_clips <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/clips"  # output
path_nDSM <- paste0(path_rasters, "/nDSM_unscaled")  # input

# set parameter
tile_size <- 0.64

################################################################################
# HELPER FUNCTIONS
################################################################################

check_create_dir <- function(path) {
  # checks if directory exists
  # if not, creates it
  if (!dir.exists(path)) {
    print("... creating new folder")
    dir.create(path)
  } else {
    print("... using existing folder")
  }
}

################################################################################
# FILTER MIDSTORY PLOTS
################################################################################

filter_midstory <- function(nDSM_unscaled_dir, plot_path, tile_size) {
  # check if nDSM values are above 2m
  print("... removing midstory points")
  # empty list for storing "bad" points
  midstory_plots <- c()
  # read as kml + transform CRS
  plots <- st_transform(st_read(plot_path), 25832)
  # get all rasters within nDSM_dir
  nDSM_list <- list.files(nDSM_unscaled_dir, pattern=".tif", recursive=TRUE)
  # calculate edge length (divisible by two)
  edge <- ((tile_size*100)%/%2)/100
  # loop through nDSMs
  for (nDSM_path in nDSM_list) {
    # load raster
    nDSM <- raster(paste0(nDSM_unscaled_dir,"/",nDSM_path))
    crs(nDSM) <- CRS("+init=EPSG:25832")
    # get points within raster extent
    subset <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(nDSM), "SpatialPolygons")), st_crs(plots)))
    # loop through plots
    for (idx in 1:nrow(subset)) {
      plot <- subset[idx,]
      # clip raster with point + edge
      center_x <- round(st_coordinates(plot)[,1], 2) # round on cm
      center_y <- round(st_coordinates(plot)[,2], 2) # round on cm
      rectangle <- extent(c(xmin=center_x-edge, xmax=center_x+edge, ymin=center_y-edge, ymax=center_y+edge))
      clip <- crop(nDSM, rectangle)
      # check height
      if ((quantile(as.vector(clip), 0.95, na.rm=T) >= 2) &  !all(is.na(as.vector(clip)))) {
        midstory_plots <- append(midstory_plots, plot$Description)
      }
    }
  }
  # delete midstory points
  filtered_plots <- plots[!plots$Description %in% midstory_plots,]
  # save filtered points
  st_write(filtered_plots, paste0(substr(plot_path, 1, nchar(plot_path)-4), "_no_midstory.kml"), delete_layer = T)
  return(paste0(substr(plot_path, 1, nchar(plot_path)-4), "_no_midstory.kml"))
}

################################################################################
# FILTER OVERLAPPING PLOTS
################################################################################

filter_overlaps <- function(plot_path, tile_size) {
  # check if nDSM values are above 2m
  print("... removing overlapping points")
  # read as kml + transform CRS
  plots <- st_transform(st_read(plot_path), 25832)
  # calculate edge length (divisible by two)
  edge <- ((tile_size*100)%/%2)/100
  # convert points to polygons with according tilesize
  center_x <- round(st_coordinates(plots)[,1], 2) # round on cm
  center_y <- round(st_coordinates(plots)[,2], 2) # round on cm
  polygon_list <- c()
  for (i in 1:nrow(plots)) {
    polygon <- extent(c(xmin=center_x[i]-edge, xmax=center_x[i]+edge, ymin=center_y[i]-edge, ymax=center_y[i]+edge))
    polygon_list <- rbind(polygon_list, st_as_sf(as(polygon, "SpatialPolygons")))
    # TODO: description must be passed over!!
  }
  polygon_list$Description <- plots$Description
  polygon_list <- st_set_crs(polygon_list, 25832)
  # iteratively remove overlapping polygons
  repeat {
    # get list of overlapping points
    overlapping_polygons <- st_intersects(polygon_list, polygon_list)
    # get the indices of overlapping polygons
    overlapping_indices <- c()
    for (i in 1:length(overlapping_polygons)) {
      if (length(overlapping_polygons[[i]]) > 1) {
        overlapping_indices <- c(overlapping_indices, overlapping_polygons[[i]])
      }
    }
    # stop, if there are no duplicates anymore
    if(length(overlapping_indices) == 0) {
      break
    }
    # otherwise, randomly remove one of the overlapping polygons
    print(paste0("number of overlapping poygons: ", length(unique(as.factor(overlapping_indices)))))
    polygon_idx <- overlapping_indices[floor(runif(1, min = 1, max=length(overlapping_indices)+1))]
    polygon_list <- polygon_list[-polygon_idx,]
  }
  # keep plots which have same geometry as filtered_points
  new_plots <- plots[plots$Description %in% polygon_list$Description,]
  # save filtered points
  st_write(new_plots, paste0(substr(plot_path, 1, nchar(plot_path)-4), "_no_overlap.kml"), delete_layer = T)
  return(paste0(substr(plot_path, 1, nchar(plot_path)-4), "_no_overlap.kml"))
}

################################################################################
# RASTER TILES
################################################################################

raster_clip_all <- function(raster_dir, plot_path, tile_size, output_dir) {
  # clips all rasters to small areas around the vegetation plots
  # creates folder structures similar to the input rasters
  check_create_dir(output_dir)
  print("... clipping all rasters")
  # read as kml + transform CRS
  plots <- st_transform(st_read(plot_path), 25832)
  # get all rasters within raster_dir (without unscaled nDSM)
  raster_list <- list.files(raster_dir, pattern=".tif", recursive=TRUE)
  raster_list <- raster_list[!grepl("nDSM_unscaled", raster_list)]
  # calculate edge length (divisible by two)
  edge <- ((tile_size*100)%/%2)/100
  # loop through rasters
  for (raster_path in raster_list) {
    # load raster
    raster <- raster(paste0(raster_dir,"/",raster_path))
    crs(raster) <- CRS("+init=EPSG:25832")
    # get points within raster extent
    subset <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(raster), "SpatialPolygons")), st_crs(plots)))
    # get raster dir & type & name
    subfolder <- strsplit(raster_path, "/")[[1]][1]
    name <- strsplit(strsplit(raster_path, "/")[[1]][2], "[.]")[[1]][1]
    # loop through plots
    for (idx in 1:nrow(subset)) {
      plot <- subset[idx,]
      # clip raster with point + edge
      center_x <- round(st_coordinates(plot)[,1], 2) # round on cm
      center_y <- round(st_coordinates(plot)[,2], 2) # round on cm
      rectangle <- extent(c(xmin=center_x-edge, xmax=center_x+edge, ymin=center_y-edge, ymax=center_y+edge))
      clip <- crop(raster, rectangle)
      # get plot plot_ID & veg_ID
      plot_id <- strsplit(strsplit(plot$Description, ",")[[1]][1], " ")[[1]][2]
      veg_id <- strsplit(strsplit(plot$Description, ",")[[1]][2], " ")[[1]][3]
      # check if raster is empty, only continue if not
      if (!all(is.na(as.vector(clip)))) {
        # save clip
        check_create_dir(paste0(output_dir, "/", subfolder))
        writeRaster(clip, paste0(output_dir, "/", subfolder, "/", name,
                                 "_", plot_id, "_", veg_id, ".tif"), overwrite=TRUE)
      }
    }
  }
}

################################################################################
# EXECUTION
################################################################################

path_vegetation <- filter_midstory(path_nDSM, path_vegetation, tile_size)
path_vegetation <- filter_overlaps(path_vegetation, tile_size)
raster_clip_all(path_rasters, path_vegetation, tile_size, path_clips)

################################################################################