################################################################################
################################################################################
# RASTER CALCULATION - LAS FILES
################################################################################
################################################################################

# load packages
library(lidR)  # for point clouds, also loads sp & raster

# set paths
path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"  # output
path_points   <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/areaXY/testing.las"  # input

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/03_raster_calculation_functions.R")

# load data
cloud_raw <- readTLSLAS(path_points)  # only for small point clouds

################################################################################
################################################################################
# RASTER CALCULATION - LAS FILES
################################################################################
################################################################################

# load packages
library(lidR)  # for point clouds, also loads sp & raster

# set paths
path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"  # output
path_points   <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/areaXY/testing.las"  # input

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/03_raster_calculation_functions.R")

# load data
cloud_raw <- readTLSLAS(path_points)  # only for small point clouds

################################################################################
# FILTER UNDERSTORY POINTS
################################################################################

# normalize height
cloud_raw <- classify_ground(cloud_raw, csf(class_threshold = 0.3, cloth_resolution = 0.3, sloop_smooth = TRUE))
dtm <- grid_terrain(filter_ground(cloud_raw), tin(), res = 0.01)
cloud_norm <- normalize_height(cloud_raw, dtm, na.rm = T)

# delete below ground points
cloud_norm <- filter_poi(cloud_norm, Z >= 0)

# delete points above 2m height
cloud_under <- filter_poi(cloud_norm, Z <= 2)

# save changed data
writeLAS(cloud_norm, paste0(substr(path_points, 1, nchar(path_points)-4), "_normalized.las"))
writeLAS(cloud_under, paste0(substr(path_points, 1, nchar(path_points)-4), "_understory.las"))

################################################################################
# CALCULATE RASTERS
################################################################################

raster_create_all(cloud_under, 0.01, path_rasters, "test")
raster_nDSM(cloud_norm, resolution, paste0(raster_dir, "/nDSM_unscaled"), "test", rescale=FALSE)

################################################################################
