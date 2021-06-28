################################################################################
################################################################################
# RASTER CALCULATION - LAS CATALOGS
################################################################################
################################################################################

# load packages
library(lidR)
library(rlas)
library(sf)
library(future)
library(raster)

# load functions
source("C:/Users/Zoe/Documents/understory_classification/5_Analyse/03_raster_calculation_functions.R")

# set chunk parameters
chunk_size <- 15  # my RAM hates everything above, so I hate everything above
buffer_size <- 1  # to avoid edge effects & not having enough points for interpolation
raster_resolution <- 0.01

# set paths
path_points <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/points/day4.laz"  # input
path_area   <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/sites/convex/area_polygons.shp"  # input
path_rasters <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/rasters_", raster_resolution*100, "cm")  # output
check_create_dir(path_rasters)

# load data
ctg <- readTLSLAScatalog(path_points)
#lidR:::catalog_laxindex(ctg)  # lax file, does not work on workstation, idk why
#plot(ctg, chunk=TRUE)

################################################################################
# CLIPPING CLOUDS TO AOI
################################################################################

# use multiple cores
plan(multisession, workers=11L)

# set options
check_create_dir(paste0(dirname(path_points), "/01_tiled"))
opt_chunk_buffer(ctg) <- 0  # otherwise chunks are saved with buffer
opt_chunk_size(ctg) <- chunk_size
plot(ctg, chunk=TRUE)

# load site polygons
area_polys <- st_read(path_area)
st_crs(area_polys) <- CRS("+init=EPSG:25832")
area_polys <- st_transform(area_polys, as.character(crs(ctg)))

# loop through site polygons
for (idx in 1:nrow(area_polys)) {
  # load polygon
  area <- area_polys[idx,]
  area_id <- idx
  # check if area & catalog are even overlapping
  if (!is.null(intersect(extent(ctg), extent(area)))) {
    # print, which area is processed
    print(paste0("... clipping area ", area_id))
    # set output path
    opt_output_files(ctg) <- paste0(dirname(path_points), "/01_tiled/area_", area_id, "_tiled_{ID}")
    # clip points to area
    ctg_retiled <- area_retile_ctg.LAScatalog(ctg, area)
    warnings()
    if (is.list(ctg_retiled)) {
      # if a list is returned, open the resulting list
      ctg_retiled <- readTLSLAScatalog(dirname(ctg_retiled[[1]]))
    }
    lidR:::catalog_laxindex(ctg_retiled)  # lax files
  }
}

# use single core
plan(sequential)

################################################################################
# TILING HUGE CLOUDS
################################################################################

# # set options
# points_name <- substr(basename(path_points), 1, nchar(basename(path_points))-4)  # for the naming pattern
# check_create_dir(paste0(dirname(path_points), "/01_tiled"))
# opt_output_files(ctg) <- paste0(dirname(path_points), "/01_tiled/", points_name, "_tiled_{ID}")
# opt_chunk_buffer(ctg) <- 0  # otherwise chunks are saved with buffer
# opt_chunk_size(ctg) <- chunk_size
# plot(ctg, chunk=TRUE)
# 
# # execute
# ctg_retiled <- catalog_retile(ctg)
# if (is.list(ctg_retiled)) {
#   ctg_retiled <- readTLSLAScatalog(dirname(ctg_retiled[[1]]))
#   # if a list is returned, open the resulting list
# }
# lidR:::catalog_laxindex(ctg_retiled)  # lax files
# warnings()

################################################################################
# CALCULATE DTM
################################################################################

# use multiple cores
plan(multisession, workers=11L)

# get all area IDs
area_IDs <- list.files(paste0(dirname(path_points), "/01_tiled"), pattern=".las")
area_IDs <- as.numeric(unique(lapply(area_IDs, function(x) strsplit(x, split="_")[[1]][2])))

for (area_ID in area_IDs) {
  # print, which area is processed
  print(paste0("... calculating DTM of area ", area_ID))
  
  # read from folder
  file_list <- list.files(paste0(dirname(path_points), "/01_tiled"), pattern=paste0("area_", area_ID), full.names=TRUE)
  file_list <- file_list[grepl("[.]las", file_list)]
  ctg_retiled <- readTLSLAScatalog(file_list)
  
  # set options
  opt_chunk_buffer(ctg_retiled) <- buffer_size * 2  # otherwise there are holes!
  opt_chunk_size(ctg_retiled) <- chunk_size

  # execute
  ctg_normalized <- dtm_ctg.LAScatalog(ctg_retiled, paste0(path_rasters, "/DTM"),
                                       paste0("area_", area_ID))
  warnings()
}

# use single core
plan(sequential)

################################################################################
# NORMALIZE POINT CLOUDS
################################################################################

# use multiple cores
plan(multisession, workers=5L)

# get all area IDs
area_IDs <- list.files(paste0(dirname(path_points), "/01_tiled"), pattern=".las")
area_IDs <- as.numeric(unique(lapply(area_IDs, function(x) strsplit(x, split="_")[[1]][2])))

for (area_ID in area_IDs) {
  # print, which area is processed
  print(paste0("... normalizing area ", area_ID))
  
  # read from folder
  file_list <- list.files(paste0(dirname(path_points), "/01_tiled"), pattern=paste0("area_", area_ID), full.names=TRUE)
  file_list <- file_list[grepl("[.]las", file_list)]
  ctg_retiled <- readTLSLAScatalog(file_list)
  
  # set options
  opt_chunk_buffer(ctg_retiled) <- buffer_size * 2  # otherwise there are holes in the dtm!
  opt_chunk_size(ctg_retiled) <- chunk_size 
  check_create_dir(paste0(dirname(path_points), "/02_normalized"))
  opt_output_files(ctg_retiled) <- paste0(dirname(path_points), "/02_normalized/area_",
                                          area_ID, "_norm_{ID}")
  
  # execute
  ctg_normalized <- normalize_ctg.LAScatalog(ctg_retiled)
  warnings()
  if (is.list(ctg_normalized)) {
    # if a list is returned, open the resulting list
    ctg_normalized <- readTLSLAScatalog(dirname(ctg_normalized[[1]]))
  }
  lidR:::catalog_laxindex(ctg_normalized)  # lax files
}

# use single core
plan(sequential)

################################################################################
# FILTER UNDERSTORY POINTS
################################################################################

# use multiple cores
plan(multisession, workers=7L)

# get all area IDs
area_IDs <- list.files(paste0(dirname(path_points), "/02_normalized"), pattern=".las")
area_IDs <- as.numeric(unique(lapply(area_IDs, function(x) strsplit(x, split="_")[[1]][2])))

################################################################################

# for CNN input creation
# height: 2m, remove stems: yes

for (area_ID in area_IDs) {
  # print, which area is processed
  print(paste0("... filtering area ", area_ID))
  
  # read from folder
  file_list <- list.files(paste0(dirname(path_points), "/02_normalized"), pattern=paste0("area_", area_ID), full.names=TRUE)
  file_list <- file_list[grepl("[.]las", file_list)]
  ctg_normalized <- readTLSLAScatalog(file_list)
  
  # set options
  opt_chunk_buffer(ctg_normalized) <- buffer_size
  opt_chunk_size(ctg_normalized) <- chunk_size
  check_create_dir(paste0(dirname(path_points), "/03_understory"))
  opt_output_files(ctg_normalized) <- paste0(dirname(path_points), "/03_understory/area_", area_ID, "_understory_{ID}")
  
  # execute
  ctg_understory <- filter_understory_ctg.LAScatalog(ctg_normalized, height=2, remove_stems=TRUE)
  warnings()
  if (is.list(ctg_understory)) {
    # if a list is returned, open the resulting list
    ctg_understory <- readTLSLAScatalog(dirname(ctg_understory[[1]]))
  }
  lidR:::catalog_laxindex(ctg_understory)
}

################################################################################

# for vegetation plot filtering
# height: 3m, remove stems: no

for (area_ID in area_IDs) {
  # print, which area is processed
  print(paste0("... filtering area ", area_ID))
  
  # read from folder
  file_list <- list.files(paste0(dirname(path_points), "/02_normalized"), pattern=paste0("area_", area_ID), full.names=TRUE)
  file_list <- file_list[grepl("[.]las", file_list)]
  ctg_normalized <- readTLSLAScatalog(file_list)
  
  # set options
  opt_chunk_buffer(ctg_normalized) <- buffer_size
  opt_chunk_size(ctg_normalized) <- chunk_size
  check_create_dir(paste0(dirname(path_points), "/04_understory_stems"))
  opt_output_files(ctg_normalized) <- paste0(dirname(path_points), "/04_understory_stems/area_", area_ID, "_understory_stems_{ID}")
  
  # execute
  ctg_understory <- filter_understory_ctg.LAScatalog(ctg_normalized, height=2.5, remove_stems=FALSE)
  warnings()
  if (is.list(ctg_understory)) {
    # if a list is returned, open the resulting list
    ctg_understory <- readTLSLAScatalog(dirname(ctg_understory[[1]]))
  }
  lidR:::catalog_laxindex(ctg_understory)
}

################################################################################

# use single core
plan(sequential)

################################################################################
# CALCULATE RASTERS
################################################################################

# use multiple cores
plan(multisession, workers=7L)
# computer unhappy with multisession & geometry raster calculation though!

################################################################################

# get all area IDs
area_IDs <- list.files(paste0(dirname(path_points), "/03_understory"), pattern=".las")
area_IDs <- as.numeric(unique(lapply(area_IDs, function(x) strsplit(x, split="_")[[1]][2])))

for (area_ID in area_IDs) {
  # print, which area is processed
  print(paste0("... calculating rasters of area ", area_ID))
  
  # read from folder
  file_list <- list.files(paste0(dirname(path_points), "/03_understory"), pattern=paste0("area_", area_ID), full.names=TRUE)
  file_list <- file_list[grepl("[.]las", file_list)]
  ctg_understory <- readTLSLAScatalog(file_list)
  
  # set options
  opt_chunk_buffer(ctg_understory) <- buffer_size
  opt_chunk_size(ctg_understory) <- chunk_size
  
  # execute - input for CNN & filtering vegetation plots
  raster_create_all_ctg(ctg_understory, raster_resolution, path_rasters, paste0("area_", area_ID), rescale=FALSE)
  warnings()
  
  # use many cores for this
  # raster_nDSM_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/nDSM"), paste0("area_", area_ID), rescale=FALSE, saving=TRUE)
  # raster_ortho_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/ortho"), paste0("area_", area_ID), rescale=FALSE, saving=TRUE)
  # raster_point_density_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/point_density"), paste0("area_", area_ID), rescale=FALSE, saving=TRUE)
  # raster_reflectance_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/reflectance"), paste0("area_", area_ID), rescale=FALSE, saving=TRUE)
  
  # use single core for this
  # raster_geometry_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/geometry"), paste0("area_", area_ID), rescale=FALSE, saving=TRUE)
}

################################################################################

# get all area IDs
area_IDs <- list.files(paste0(dirname(path_points), "/04_understory_stems"), pattern=".las")
area_IDs <- as.numeric(unique(lapply(area_IDs, function(x) strsplit(x, split="_")[[1]][2])))

for (area_ID in area_IDs) {
  # print, which area is processed
  print(paste0("... calculating rasters of area ", area_ID))
  
  # read from folder
  file_list <- list.files(paste0(dirname(path_points), "/04_understory_stems"), pattern=paste0("area_", area_ID), full.names=TRUE)
  file_list <- file_list[grepl("[.]las", file_list)]
  ctg_understory <- readTLSLAScatalog(file_list)
  
  # set options
  opt_chunk_buffer(ctg_understory) <- buffer_size
  opt_chunk_size(ctg_understory) <- chunk_size
  
  # execute - input for CNN & filtering vegetation plots
  raster_nDSM_ctg.LAScatalog(ctg_understory, raster_resolution, paste0(path_rasters, "/nDSM_filtering"),
                             paste0("area_", area_ID), rescale=FALSE)
  warnings()
}

################################################################################

# use single core
plan(sequential)

################################################################################
