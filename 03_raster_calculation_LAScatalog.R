################################################################################
################################################################################
# RASTER CALCULATION - LAS CATALOGS
################################################################################
################################################################################

# load packages
library(lidR)  # for point clouds, also loads sp & raster

# set paths
path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"  # output
path_points <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/areaXY/testing.las"  # input
#path_points <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/areaXY/OT01cm.laz"  # input
points_name <- substr(basename(path_points), 1, nchar(basename(path_points))-4)  # for the naming pattern

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/03_raster_calculation_functions.R")

# load data
ctg <- readTLSLAScatalog(path_points)
plot(ctg, chunk=TRUE)

# separate huge LAS file into LAScatalog
check_create_dir(paste0(dirname(path_points), "/01_tiled"))
opt_output_files(ctg) <- paste0(dirname(path_points), "/01_tiled/", points_name, "_{ID}")
opt_chunk_buffer(ctg) <- 0  # otherwise chunks are saved with buffer, bad practice
opt_chunk_size(ctg) <- 25
x_corner <- floor(bbox(ctg)[1,1])  # x left
y_corner <- floor(bbox(ctg)[2,1])  # y bottom
opt_chunk_alignment(ctg) <- c(x_corner,y_corner)
# plot(ctg, chunk = TRUE)
ctg_retiled <- catalog_retile(ctg)

################################################################################
# FILTER UNDERSTORY POINTS
################################################################################

normalize_ctg.LAScluster <- function(las) {
  # load the data
  las <- readLAS(las)
  if (is.empty(las)) return(NULL)
  # classify & normalize
  las <- classify_ground(las, csf(class_threshold = 0.3, cloth_resolution = 0.3, sloop_smooth = TRUE))  # has to be small due to smaller areas & heavy slope 
  dtm <- grid_terrain(las, tin(), res = 0.01)  # if bigger, there are many artifacts
  las <- normalize_height(las, dtm, na.rm = T)
  las <- filter_poi(las, Z >= 0)
  # TODO: delete, dummy RGB data
  las <- add_lasrgb(las, R=as.integer(floor(runif(las@data$X)*255)),
                    G=as.integer(floor(runif(las@data$X)*255)),
                    B=as.integer(floor(runif(las@data$X)*255)))
  # delete buffer & return points
  las <- filter_poi(las, buffer == 0)
  return(las)
}

normalize_ctg.LAScatalog <- function(las) {
  # undo previous selections
  opt_select(las) <-  "*"
  # set paramters
  options <- list(
    need_output_file = TRUE,  # output path necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = TRUE)  # combine outputs
  # execute & return
  output  <- catalog_apply(las, normalize_ctg.LAScluster, .options = options)
  return(output)
}

################################################################################

remove_understory_ctg.LAScluster <- function(las, height) {
  # load the data
  las <- readLAS(las)
  if (is.empty(las)) return(NULL)
  # remove everything below certain height
  las <- filter_poi(las, Z <= height)
  # TODO: try to remove stems & floating stuff?
  # delete buffer & return points
  las <- filter_poi(las, buffer == 0)
  return(las)
}

remove_understory_ctg.LAScatalog <- function(las, height) {
  # undo previous selections
  opt_select(las) <-  "*"
  # set paramters
  options <- list(
    need_output_file = TRUE,  # output path necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = TRUE)  # combine outputs
  # execute & return
  output  <- catalog_apply(las, remove_understory_ctg.LAScluster, height = height, .options = options)
  return(output)
}

################################################################################

lax_for_las <- function(point_dir) {
  # creates for all las files of a dir lax files
  # (lax files provide performance improvements)
  las_files <- list.files(point_dir, pattern=".las", recursive=TRUE)
  for (las_file in las_files) {
    rlas::writelax(paste0(point_dir, "/", las_file))
  }
}

################################################################################

# normalize the point clouds
opt_chunk_buffer(ctg_retiled) <- 0.5
check_create_dir(paste0(dirname(path_points), "/02_normalized"))
opt_output_files(ctg_retiled) <- paste0(dirname(path_points), "/02_normalized/", points_name, "_{ID}")
ctg_normalized <- normalize_ctg.LAScatalog(ctg_retiled)

# remove everything above 2 m height
opt_chunk_buffer(ctg_normalized) <- 0.5
check_create_dir(paste0(dirname(path_points), "/03_understory"))
opt_output_files(ctg_normalized) <- paste0(dirname(path_points), "/03_understory/", points_name, "_{ID}")
ctg_understory <- remove_understory_ctg.LAScatalog(ctg_normalized, height=2)

# create lax files for everything (for faster processing)
lax_for_las(dirname(path_points))

################################################################################
# CALCULATE RASTERS
################################################################################

opt_chunk_buffer(ctg_understory) <- 0.5
opt_chunk_size(ctg_understory) <- 25
x_corner <- floor(bbox(ctg_understory)[1,1])  # x left
y_corner <- floor(bbox(ctg_understory)[2,1])  # y bottom
opt_chunk_alignment(ctg_understory) <- c(x_corner,y_corner)
opt_output_files(ctg_understory) <- ""

raster_create_all_ctg(ctg_understory, 0.01, path_rasters, points_name)

################################################################################

opt_chunk_buffer(ctg_normalized) <- 0.5
opt_chunk_size(ctg_normalized) <- 25
x_corner <- floor(bbox(ctg_normalized)[1,1])  # x left
y_corner <- floor(bbox(ctg_normalized)[2,1])  # y bottom
opt_chunk_alignment(ctg_normalized) <- c(x_corner,y_corner)
opt_output_files(ctg_normalized) <- ""

raster_nDSM_ctg.LAScatalog(ctg_normalized, 0.01, paste0(path_rasters, "/nDSM_unscaled"), points_name, rescale=FALSE)

################################################################################
