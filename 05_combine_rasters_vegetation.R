################################################################################
################################################################################
# COMBINING RASTER & VEGATATION
################################################################################
################################################################################

# load packages
library(lidR)
library(sf)

# set paths
path_rasters  <- "D:/Masterarbeit_Zoe/4_Daten/rasters"  # input
path_vegetation <- "D:/Masterarbeit_Zoe/4_Daten/vegetation/Export_ODK_clean_checked.kml"  # input
path_clips <- "D:/Masterarbeit_Zoe/4_Daten/clips"  # output
path_nDSM <- paste0(path_rasters, "/nDSM_filtering")  # input
path_points <- "D:/Masterarbeit_Zoe/4_Daten/points/actual_data/04_understory_stems"  # input

# rasters to be clipped, results from collinearity check
clip_these <- c("ortho", "anisotropy_max", "curvature_max", "linearity_max",
                "linearity_sd", "planarity_mean", "planarity_sd", "nDSM",
                "point_density", "reflectance_mean", "reflectance_sd")

# set parameter
tile_size <- 0.5
crs_raster_las <- as.character(crs(readTLSLAScatalog(path_points)))

################################################################################
# HELPER FUNCTIONS
################################################################################

check_create_dir <- function(path) {
  # checks if directory exists
  # if not, creates it
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

################################################################################

rescale_values <- function(folder) {
  # gets minimum & maximum values for rescaling
  # get all raster
  raster_list <- list.files(folder, pattern=".tif", full.names=TRUE, recursive=TRUE)
  raster_list <- raster_list[!grepl("temp", raster_list)]
  # set up empty lookup list
  lookup_list <- list()
  # get all unique raster types
  types <- c()
  for (i in 1:length(raster_list)) {types[i] <- strsplit(basename(raster_list[i]), "_area_")[[1]][1]}
  types <- unique(types)
  # loop through all unique raster types
  for (type_idx in 1:length(types)) {
    # get all raster paths with that type
    type <- types[type_idx]
    type_paths <- raster_list[grepl(type, raster_list)]
    # set up min / max value
    val_min = c()
    val_max = c()
    # loop through all rasters of that type
    for (idx in 1:length(type_paths)) {
      # load raster & get extreme values
      type_raster <- stack(type_paths[idx])
      val_min_temp = min(minValue(type_raster))
      val_max_temp = max(maxValue(type_raster))
      # depending on value & index, set new min / max value
      if (idx == 1) {val_min <- val_min_temp}
      if (val_min_temp < val_min) {val_min <- val_min_temp}
      if (idx == 1) {val_max <- val_max_temp}
      if (val_max_temp > val_max) {val_max <- val_max_temp}
    }
    # save final min and max value un lookup list
    lookup_list[paste0(type, "_min")] <- val_min
    lookup_list[paste0(type, "_max")] <- val_max
  }
  # return full lookup list
  return(lookup_list)
}

################################################################################
# FILTER MIDSTORY PLOTS
################################################################################

# this section would be a lot shorter, if I would habe measured the vegetation
# height during field work & would have excluded vegetation > 2m height

filter_midstory_nDSM <- function(plots, nDSM_unscaled_dir, tile_size) {
  # check if maximum nDSM values are above 1,9m
  # otherwise: delete
  print("... removing midstory points")
  # empty list for storing "bad" points
  remove_plots <- c()
  # get all rasters within nDSM_dir
  nDSM_list <- list.files(nDSM_unscaled_dir, pattern=".tif", recursive=TRUE)
  # calculate edge length (divisible by two)
  edge <- ((tile_size*100)%/%2)/100
  # loop through nDSMs
  for (nDSM_path in nDSM_list) {
    # load raster
    nDSM <- raster(paste0(nDSM_unscaled_dir, "/", nDSM_path))
    crs(nDSM) <- crs_raster_las
    # get points within raster extent
    subset <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(nDSM), "SpatialPolygons")), st_crs(plots)))
    # loop through plots
    for (idx in 1:nrow(subset)) {
      plot <- subset[idx,]
      # clip raster with point + edge
      center_x <- round(st_coordinates(plot)[,1], 2) # round on cm
      center_y <- round(st_coordinates(plot)[,2], 2) # round on cm
      rectangle <- extent(c(xmin=center_x-edge, xmax=center_x+edge,
                            ymin=center_y-edge, ymax=center_y+edge))
      clip <- crop(nDSM, rectangle)
      # check height & cover
      if ((max(as.vector(clip), na.rm=T) >= 1.9) & !all(is.na(as.vector(clip)))) {  # keeps empty clips
        # TODO: guter Punkt um festzulegen, dass gewisser Pixel-Bedeckungsgrad vorhanden sein muss?
        remove_plots <- rbind(remove_plots, plot)
      }
    }
  }
  # delete midstory points
  keep_plots <- plots[!(plots$Description %in% remove_plots$Description),]
  # save filtered points
  return(list(keep = keep_plots, remove = remove_plots))
}

################################################################################

# filter_midstory_points <- function(plots, point_cloud_dir, tile_size) {
#   # check if there are empty vertical bins
#   # (assumption: flying branches from overstory lead to empty vertical bins)
#   # otherwise: delete
#   print("... removing midstory points")
#   # empty list for storing "bad" points
#   remove_plots <- c()
#   # load point clouds as LAScatalog
#   las <- readTLSLAScatalog(point_cloud_dir)  # from all areas combined
#   # calculate edge length (divisible by two)
#   edge <- ((tile_size*100)%/%2)/100
#   # loop through plots
#   for (idx in 1:nrow(plots)) {
#     plot <- plots[idx,]
#     # clip point cloud with point + edge
#     center_x <- round(st_coordinates(plot)[,1], 2) # round on cm
#     center_y <- round(st_coordinates(plot)[,2], 2) # round on cm
#     rectangle <- extent(c(xmin=center_x-edge, xmax=center_x+edge,
#                           ymin=center_y-edge, ymax=center_y+edge))
#     clip <- clip_roi(las, rectangle)
#     # bin the points, 10cm vertical bins until 2m height
#     bin_boolean <- c()
#     for (i in 1:20) {
#       upper_bound <- (i * 10)/100
#       lower_bound <- (upper_bound*100 - 10)/100
#       bin <- 0 < length(filter_poi(clip, Z >= lower_bound & Z < upper_bound)@data$X)  # TRUE: points, FALSE: no points
#       bin_boolean <- c(bin_boolean, bin)
#     }
#     # check if there are empty bins between
#     if (mean(bin_boolean) == 1) {
#       remove_plots <- rbind(remove_plots, plot)
#     }
#   }
#   # delete midstory points
#   keep_plots <- plots[!(plots$Description %in% remove_plots$Description),]
#   # save filtered points
#   return(list(keep = keep_plots, remove = remove_plots))
# }

################################################################################

filter_midstory_visually <- function(plots, point_cloud_dir, tile_size) {
  # check visually & manually if plot should be kept
  print("... removing midstory points")
  # empty list for storing "bad" / "good" points
  remove_plots <- c()
  keep_plots <- c()
  # load point clouds as LAScatalog
  las <- readTLSLAScatalog(point_cloud_dir)  # from all areas combined
  # calculate edge length (divisible by two)
  edge <- ((tile_size*100)%/%2)/100
  # loop through plots
  for (idx in 1:nrow(plots)) {
    plot <- plots[idx,]
    # clip point cloud with point + edge
    center_x <- round(st_coordinates(plot)[,1], 2) # round on cm
    center_y <- round(st_coordinates(plot)[,2], 2) # round on cm
    rectangle <- extent(c(xmin=center_x-edge, xmax=center_x+edge,
                          ymin=center_y-edge, ymax=center_y+edge))
    clip <- clip_roi(las, rectangle)
    # categorize points into 2m or higher
    below_above <- clip@data$Z < 2
    clip <- add_attribute(clip, below_above, "below_above")
    plot(clip, color="below_above")
    print(paste0("Label: ", plot$Name))
    user_input <- readline(prompt="Keep this point? Enter y or n :  ")
    # check if the user wants the point removed
    if (user_input == "n") {
      remove_plots <- rbind(remove_plots, plot)
    } else if (user_input == "y") {
      keep_plots <- rbind(keep_plots, plot)
    } else if (user_input == "exit") {
      print("... mission aborted")
      return(0)
    }
  }
  # save filtered points
  return(list(keep = keep_plots, remove = remove_plots))
}

################################################################################

filter_midstory_all <- function(plot_path, nDSM_unscaled_dir, point_cloud_dir, tile_size) {
  # execute all filtering functions together
  # load plots
  plots <- st_transform(st_read(plot_path), crs_raster_las)
  # filter with nDSM
  out_ndsm <- filter_midstory_nDSM(plots, nDSM_unscaled_dir, tile_size)
  plots <- out_ndsm[[2]]
  # filter visually & manually
  out_visual <- filter_midstory_visually(plots, point_cloud_dir, tile_size)
  # merge all plots to be kept
  keep_plots <- rbind(out_ndsm[[1]], out_visual[[1]])
  # save & return path to filtered kml
  st_write(keep_plots, paste0(substr(plot_path, 1, nchar(plot_path)-4), "_filtered.kml"), delete_layer = T)
  return(paste0(substr(plot_path, 1, nchar(plot_path)-4), "_filtered.kml"))
}

################################################################################
# FILTER OVERLAPPING PLOTS
################################################################################

filter_overlaps <- function(plot_path, tile_size) {
  # check if vegetation tiles are not overlapping
  # otherwise: delete
  print("... removing overlapping points")
  # read as kml + transform CRS
  plots <- st_transform(st_read(plot_path), crs_raster_las)
  # calculate edge length (divisible by two)
  edge <- ((tile_size*100)%/%2)/100
  # convert points to polygons with according tile size
  center_x <- round(st_coordinates(plots)[,1], 2) # round on cm
  center_y <- round(st_coordinates(plots)[,2], 2) # round on cm
  polygon_list <- c()
  for (i in 1:nrow(plots)) {
    polygon <- extent(c(xmin=center_x[i]-edge, xmax=center_x[i]+edge, ymin=center_y[i]-edge, ymax=center_y[i]+edge))
    polygon_list <- rbind(polygon_list, st_as_sf(as(polygon, "SpatialPolygons")))
  }
  polygon_list$Description <- plots$Description
  polygon_list <- st_set_crs(polygon_list, crs_raster_las)
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
    polygon_idx <- overlapping_indices[floor(runif(1, min = 1, max=length(overlapping_indices)+1))]
    polygon_list <- polygon_list[-polygon_idx,]
  }
  # keep plots which have same geometry as filtered_points
  new_plots <- plots[plots$Description %in% polygon_list$Description,]
  # give feedback on how much was removed:
  print(paste0("number of removed plots: ", nrow(plots)-nrow(new_plots)))
  # save filtered points
  st_write(new_plots, paste0(substr(plot_path, 1, nchar(plot_path)-4), "_no_overlap.kml"), delete_layer = T)
  return(paste0(substr(plot_path, 1, nchar(plot_path)-4), "_no_overlap.kml"))
}

################################################################################
# CREATE & RESCALE RASTER TILES
################################################################################

raster_clip_all <- function(raster_dir, plot_path, output_dir, selection, tile_size, rescale=TRUE) {
  # clips all rasters to small areas around the vegetation plots
  # creates folder structures similar to the input rasters
  check_create_dir(output_dir)
  print("... clipping all rasters")
  # read as kml + transform CRS
  plots <- st_transform(st_read(plot_path), crs_raster_las)
  # get all rasters within raster_dir (without unscaled nDSM)
  raster_list <- list.files(raster_dir, pattern=".tif", recursive=TRUE)
  raster_list <- raster_list[!grepl("nDSM_filtering", raster_list)]
  raster_list <- raster_list[!grepl("temp", raster_list)]
  # calculate edge length (divisible by two)
  edge <- ((tile_size*100)%/%2)/100
  # create rescaling lookup table
  if (rescale) {rescale_lookup <- rescale_values(raster_dir)}
  # loop through rasters
  for (raster_path in raster_list) {
    # get raster dir & type & name
    subfolder <- strsplit(raster_path, "/")[[1]][1]
    name <- strsplit(strsplit(raster_path, "/")[[1]][2], "[.]")[[1]][1]
    type <- strsplit(name, "_area_")[[1]][1]
    # check if raster is of interest
    if (type %in% selection) {
      check_create_dir(paste0(output_dir, "/", subfolder))
      # load raster
      raster <- stack(paste0(raster_dir, "/", raster_path))
      crs(raster) <- crs_raster_las
      # get points within raster extent
      subset <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(raster), "SpatialPolygons")), st_crs(plots)))
      # rescale raster
      if (rescale) {
        rescale_min <- rescale_lookup[[paste0(type, "_min")]]
        rescale_max <- rescale_lookup[[paste0(type, "_max")]]
        raster <- (raster-rescale_min)/(rescale_max-rescale_min)
      }
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
          writeRaster(clip, paste0(output_dir, "/", subfolder, "/", name,
                                   "_", plot_id, "_", veg_id, ".tif"), overwrite=TRUE)
        }
      }
    }
  }
}

################################################################################
# EXECUTION
################################################################################

# TODO: eventuell markieren, dass überlappende immer nur ins gleiche Datenset dürfen,
#       oder nur raus werfen ab gewissen overlap Prozent
#       -> 8 plots would be removed due to overlap -> okay

new_path_vegetation <- filter_midstory_all(path_vegetation, path_nDSM, path_points, tile_size)
newer_path_vegetation <- filter_overlaps(new_path_vegetation, tile_size)
raster_clip_all(path_rasters, newer_path_vegetation, path_clips, clip_these, tile_size)

################################################################################