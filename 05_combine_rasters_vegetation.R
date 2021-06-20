################################################################################
################################################################################
# COMBINING RASTER & VEGATATION
################################################################################
################################################################################

# load packages
library(lidR)
library(sf)
library(rgl)
library(jpeg)

# set paths
path_vegetation_kml <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/vegetation/Export_ODK_clean_checked.kml"  # input
path_vegetation_csv <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/vegetation/Export_ODK_clean_checked.csv"  # input
path_vegetation_img <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/vegetation/images"  # input
path_clips          <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/clips_1cm"  # output
path_points_03      <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/points/03_understory"  # input
path_points_04      <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/points/04_understory_stems"  # input
path_rasters        <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/rasters_1cm"  # input

# rasters to be clipped, results from collinearity check
clip_these <- c("ortho", "anisotropy_max", "curvature_max", "linearity_max",
                "linearity_sd", "nDSM", "planarity_mean", "planarity_sd",
                "point_density", "reflectance_mean", "reflectance_sd")

# set parameter
tile_size <- 0.5
crs_points_raster <- as.character(crs(readTLSLAScatalog(path_points_03)))

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

# filter_midstory_nDSM <- function(plots, nDSM_unscaled_dir, tile_size) {
#   # check if maximum nDSM values are above 1,9m
#   # otherwise: delete
#   print("... removing midstory points")
#   # empty list for storing "bad" points
#   remove_plots <- c()
#   # get all rasters within nDSM_dir
#   nDSM_list <- list.files(nDSM_unscaled_dir, pattern=".tif", recursive=TRUE)
#   # calculate edge length (divisible by two)
#   edge <- ((tile_size*100)%/%2)/100
#   # loop through nDSMs
#   for (nDSM_path in nDSM_list) {
#     # load raster
#     nDSM <- raster(paste0(nDSM_unscaled_dir, "/", nDSM_path))
#     crs(nDSM) <- crs_raster_las
#     # get points within raster extent
#     subset <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(nDSM), "SpatialPolygons")), st_crs(plots)))
#     # loop through plots
#     for (idx in 1:nrow(subset)) {
#       plot <- subset[idx,]
#       # clip raster with point + edge
#       center_x <- round(st_coordinates(plot)[,1], 2) # round on cm
#       center_y <- round(st_coordinates(plot)[,2], 2) # round on cm
#       rectangle <- extent(c(xmin=center_x-edge, xmax=center_x+edge,
#                             ymin=center_y-edge, ymax=center_y+edge))
#       clip <- crop(nDSM, rectangle)
#       # check height & cover
#       if ((max(as.vector(clip), na.rm=T) >= 1.9) & !all(is.na(as.vector(clip)))) {  # keeps empty clips
#         # TODO: guter Punkt um festzulegen, dass gewisser Pixel-Bedeckungsgrad vorhanden sein muss?
#         remove_plots <- rbind(remove_plots, plot)
#       }
#     }
#   }
#   # delete midstory points
#   keep_plots <- plots[!(plots$Description %in% remove_plots$Description),]
#   # save filtered points
#   return(list(keep = keep_plots, remove = remove_plots))
# }

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

filter_midstory_visually <- function(plots_kml, plots_csv, plots_img,
                                     point_cloud_dir_03, point_cloud_dir_04,
                                     tile_size, crs_points) {
  # loading data
  plots <- st_transform(st_read(plots_kml), crs_points)
  # check visually & manually if plot should be kept
  # using both point clouds & images & labels
  print("... removing midstory points")
  print(paste0("... checking ", nrow(plots), " plots manually"))
  # empty list for storing "bad" / "good" points
  remove_plots <- c()
  keep_plots <- c()
  # load data
  las_03 <- readTLSLAScatalog(point_cloud_dir_03)
  las_04 <- readTLSLAScatalog(point_cloud_dir_04)
  df <- read.csv(plots_csv)
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
    clip_03 <- clip_roi(las_03, rectangle)
    clip_04 <- clip_roi(las_04, rectangle)
    # categorize points into 2m or higher
    below_above <- clip_04@data$Z < 2
    clip_04 <- add_attribute(clip_04, below_above, "below_above")
    # show 3D plot % image & label
    if(!is.empty(clip_03) & !is.empty(clip_04)) {
      plot(clip_03, axis=T)
      plot(clip_04, color="below_above", axis=T)
      img_path <- paste0(plots_img, "/", df$filename[paste0("plot_ID: ", df$plot_ID, ", veg_ID: ", df$veg_ID) == plot$Description])
      jj <- readJPEG(img_path, native=TRUE)
      (mar=c(0,0,0,0))
      plot(0:1, 0:1, type="n", ann = FALSE, axes = FALSE)
      rasterImage(jj,0,0,1,1)
      print(paste0("ID: ", plot$Description))
      print(paste0("Label: ", plot$Name))
      print(plots_img)
      # ask user whether to keep or not
      user_input <- readline(prompt="Keep this point? Enter y or n :  ")
      # check if the user wants the point removed
      if (user_input == "n") {
        remove_plots <- rbind(remove_plots, plot)
      } else if (user_input == "y") {
        keep_plots <- rbind(keep_plots, plot)
      } else if (user_input == "exit") {
        print("... mission aborted")
        return(list(keep = keep_plots, remove = remove_plots))
      }
    } else {
      print("... empty")
    }
  }
  # save filtered points
  st_write(keep_plots, paste0(substr(plots_kml, 1, nchar(plots_kml)-4), "_filtered.kml"), delete_layer = T)
  return(paste0(substr(plots_kml, 1, nchar(plots_kml)-4), "_filtered.kml"))
}

################################################################################
# FILTER OVERLAPPING PLOTS
################################################################################

filter_overlaps <- function(plot_path, tile_size, crs_points) {
  # check if vegetation tiles are not overlapping
  # otherwise: delete
  print("... removing overlapping points")
  # read as kml + transform CRS
  plots <- st_transform(st_read(plot_path), crs_points)
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
  polygon_list <- st_set_crs(polygon_list, crs_points)
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

raster_clip_all <- function(raster_dir, plot_path, output_dir, selection, tile_size, crs_points, rescale=TRUE) {
  # clips all rasters to small areas around the vegetation plots
  # creates folder structures similar to the input rasters
  check_create_dir(output_dir)
  print("... clipping all rasters")
  # read as kml + transform CRS
  plots <- st_transform(st_read(plot_path), crs_points)
  # get all rasters within raster_dir (without unscaled nDSM)
  raster_list <- list.files(raster_dir, pattern=".tif", recursive=TRUE)
  raster_list <- raster_list[!grepl("nDSM_filtering", raster_list)]
  raster_list <- raster_list[!grepl("DTM", raster_list)]
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
      crs(raster) <- crs_points
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

# filter everything visually (don't execute this again!)
new_path_vegetation <- filter_midstory_visually(path_vegetation_kml, path_vegetation_csv, path_vegetation_img,
                                                path_points_03, path_points_04, tile_size, crs_points_raster)
new_path_vegetation <- paste0(substr(path_vegetation_kml, 1, nchar(path_vegetation_kml)-4), "_filtered.kml")

# filter overlapping tiles
newer_path_vegetation <- filter_overlaps(new_path_vegetation, tile_size, crs_points_raster)
newer_path_vegetation <- paste0(substr(new_path_vegetation, 1, nchar(new_path_vegetation)-4), "_no_overlap.kml")

# clipping rasters to point tiles
raster_clip_all(path_rasters, newer_path_vegetation, path_clips, clip_these, tile_size, crs_points_raster)

################################################################################
