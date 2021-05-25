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
#path_points <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/areaXY/testing.las"  # input
path_points <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/areaXY/OT01cm.laz"  # input
points_name <- substr(basename(path_points), 1, nchar(basename(path_points))-4)  # for the naming pattern

# set chunk parameters
chunk_size <- 25
buffer_size <- 0.5
raster_resolution <- 0.01

# load data
ctg <- readTLSLAScatalog(path_points)
lidR:::catalog_laxindex(ctg)  # lax file for big file
plot(ctg, chunk=TRUE)

# separate huge LAS file into LAScatalog
check_create_dir(paste0(dirname(path_points), "/01_tiled"))
opt_output_files(ctg) <- paste0(dirname(path_points), "/01_tiled/", points_name, "_{ID}")
opt_chunk_buffer(ctg) <- 0  # otherwise chunks are saved with buffer
opt_chunk_size(ctg) <- chunk_size
x_corner <- floor(bbox(ctg)[1,1])  # x left
y_corner <- floor(bbox(ctg)[2,1])  # y bottom
opt_chunk_alignment(ctg) <- c(x_corner,y_corner)
# plot(ctg, chunk = TRUE)
ctg_retiled <- catalog_retile(ctg)
lidR:::catalog_laxindex(ctg_retiled)  # lax files for small files
warnings()

################################################################################
# FILTER UNDERSTORY POINTS
################################################################################

# normalize the point clouds
opt_chunk_buffer(ctg_retiled) <- buffer_size
check_create_dir(paste0(dirname(path_points), "/02_normalized"))
opt_output_files(ctg_retiled) <- paste0(dirname(path_points), "/02_normalized/", points_name, "_{ID}")
ctg_normalized <- normalize_ctg.LAScatalog(ctg_retiled)
lidR:::catalog_laxindex(ctg_normalized)
warnings()

# remove everything above 2 m height
opt_chunk_buffer(ctg_normalized) <- buffer_size
check_create_dir(paste0(dirname(path_points), "/03_understory"))
opt_output_files(ctg_normalized) <- paste0(dirname(path_points), "/03_understory/", points_name, "_{ID}")
ctg_understory <- remove_understory_ctg.LAScatalog(ctg_normalized, height=2)
lidR:::catalog_laxindex(ctg_understory)
warnings()

################################################################################
# CALCULATE RASTERS
################################################################################

opt_chunk_buffer(ctg_understory) <- buffer_size
opt_chunk_size(ctg_understory) <- chunk_size
x_corner <- floor(bbox(ctg_understory)[1,1])  # x left
y_corner <- floor(bbox(ctg_understory)[2,1])  # y bottom
opt_chunk_alignment(ctg_understory) <- c(x_corner,y_corner)
opt_output_files(ctg_understory) <- ""

raster_create_all_ctg(ctg_understory, raster_resolution, path_rasters, points_name)
warnings()

################################################################################

opt_chunk_buffer(ctg_normalized) <- buffer_size
opt_chunk_size(ctg_normalized) <- chunk_size
x_corner <- floor(bbox(ctg_normalized)[1,1])  # x left
y_corner <- floor(bbox(ctg_normalized)[2,1])  # y bottom
opt_chunk_alignment(ctg_normalized) <- c(x_corner,y_corner)
opt_output_files(ctg_normalized) <- ""

raster_nDSM_ctg.LAScatalog(ctg_normalized, raster_resolution, paste0(path_rasters, "/nDSM_unscaled"),
                           points_name, rescale=FALSE)
warnings()

################################################################################
