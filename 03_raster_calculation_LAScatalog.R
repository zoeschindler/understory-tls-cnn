################################################################################
################################################################################
# RASTER CALCULATION - LAS CATALOGS
################################################################################
################################################################################

# load packages
library(lidR)  # for point clouds, also loads sp & raster

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/03_raster_calculation_functions.R")

# set paths
path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"  # output
check_create_dir(path_rasters)
path_points <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/actual_data/OT01cm.laz"  # input
points_name <- substr(basename(path_points), 1, nchar(basename(path_points))-4)  # for the naming pattern

# set chunk parameters
chunk_size <- 15  # my RAM hates everything above, so I hate everything above
buffer_size <- 1  # to avoid edge effects & not having enough points for interpolation
raster_resolution <- 0.01

# load data
ctg <- readTLSLAScatalog(path_points)
lidR:::catalog_laxindex(ctg)  # lax file
plot(ctg, chunk=TRUE)

################################################################################
# TILING HUGE CLOUD
################################################################################

# set options
check_create_dir(paste0(dirname(path_points), "/01_tiled"))
opt_output_files(ctg) <- paste0(dirname(path_points), "/01_tiled/", points_name, "_tiled_{ID}")
opt_chunk_buffer(ctg) <- 0  # otherwise chunks are saved with buffer
opt_chunk_size(ctg) <- chunk_size
plot(ctg, chunk=TRUE)

# execute
ctg_retiled <- catalog_retile(ctg)
if (is.list(ctg_retiled)) {
  ctg_retiled <- readTLSLAScatalog(dirname(ctg_retiled[[1]]))
  # if a list is returned, open the resulting list
}
lidR:::catalog_laxindex(ctg_retiled)  # lax files
warnings()

################################################################################
# NORMALIZE POINT CLOUDS
################################################################################

# read from folder
ctg_retiled <- readTLSLAScatalog(paste0(dirname(path_points), "/01_tiled"))

# set options
opt_chunk_buffer(ctg_retiled) <- buffer_size
opt_chunk_size(ctg_retiled) <- chunk_size
check_create_dir(paste0(dirname(path_points), "/02_normalized"))
opt_output_files(ctg_retiled) <- paste0(dirname(path_points), "/02_normalized/",
                                        points_name, "_normalized_{ID}")

# execute
ctg_normalized <- normalize_ctg.LAScatalog(ctg_retiled)
warnings()
if (is.list(ctg_normalized)) {
  # if a list is returned, open the resulting list
  ctg_normalized <- readTLSLAScatalog(dirname(ctg_normalized[[1]]))
}
lidR:::catalog_laxindex(ctg_normalized)  # lax files

################################################################################
# FILTER UNDERSTORY POINTS
################################################################################

# basically unnecessary, because I could filter catalog by height
# but nice to have these point clouds for later / other programs

# read from folder
ctg_normalized <- readTLSLAScatalog(paste0(dirname(path_points), "/02_normalized"))

# set options
opt_chunk_buffer(ctg_normalized) <- 0
opt_chunk_size(ctg_normalized) <- 0
check_create_dir(paste0(dirname(path_points), "/03_understory"))
opt_output_files(ctg_normalized) <- paste0(dirname(path_points), "/03_understory/",
                                           points_name, "_understory_{ID}")

# execute
ctg_understory <- remove_understory_ctg.LAScatalog(ctg_normalized, height=2)
warnings()
if (is.list(ctg_understory)) {
  # if a list is returned, open the resulting list
  ctg_understory <- readTLSLAScatalog(dirname(ctg_understory[[1]]))
}
lidR:::catalog_laxindex(ctg_understory)

################################################################################
# CALCULATE RASTERS
################################################################################

# read from folder
ctg_understory <- readTLSLAScatalog(paste0(dirname(path_points), "/03_understory"))

# set options
opt_chunk_buffer(ctg_understory) <- buffer_size
opt_chunk_size(ctg_understory) <- chunk_size
#opt_output_files(ctg_understory) <- ""
check_create_dir(paste0(path_rasters, "/temp"))
opt_output_files(ctg_understory) <- paste0(path_rasters, "/temp/temp_{ID}")
# muss noch geändert werden, wenn das für alle Erstellungen genutzt werden soll, dann in Funktionen selber!

# execute - input for CNN
raster_create_all_ctg(ctg_understory, raster_resolution, path_rasters, points_name, rescale=FALSE)
warnings()
# TODO: kills itself because of RAM in geometry feature calculation

###
# TODO: testing, delete later
# testing single functions

raster_reflectance_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/reflectance"),
                                  points_name, rescale=FALSE)
warnings()  # TODO: does not work

raster_geometry_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/reflectance"),
                               points_name, rescale=FALSE)
warnings()  # TODO: does not work
###

# not rescaled, so I can normalize later, per area or per everything
# will be normalized & used as an input for the CNN

################################################################################

# execute - for filtering vegetation plots
raster_nDSM_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/nDSM_unscaled"),
                           points_name, rescale=FALSE)
warnings()

# will be used for filterung vegetation plots with heights above 2m
# cut point cloud is used because otherwise everything would be excluded due to overstory

################################################################################
# CLIP RASTERS TO AOI
################################################################################

# TODO
# need shapes

################################################################################
# NORMALIZE RASTERS
################################################################################

# TODO
# need clipped rasters

# based on rasters of all areas combined?

# based on each area at once?

################################################################################