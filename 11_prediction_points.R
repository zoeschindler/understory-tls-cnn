################################################################################
################################################################################
# PREDICTION POINTS
################################################################################
################################################################################

# load packages
library(lidR)
library(sf)
library(future)

# set paths
basedir <- "C:/Users/Zoe/Documents/understory_classification"
path_shp <- paste0(basedir, "/4_Daten/prediction_map/tiles.shp")
path_ctg <- paste0(basedir, "/4_Daten/points/03_understory")
path_out <- paste0(basedir, "/4_Daten/points/area_6_preds")

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
# ADD ATTRIBUTES
################################################################################

merge_shape_attributes.LAScluster <- function(chunk, shp, attributes = c("prediction", "chance")) {
  # returns LAS with additional attributes from shapes
  # load the data
  las <- readLAS(chunk)
  if (is.empty(las)) {
    return(NULL)
  }
  # add all attributes
  for (attr in attributes) {
    las <- merge_spatial(las, shp, attr)
    las <- add_lasattribute(las, name = attr, desc = attr)
  }
  return(las)
}

merge_shape_attributes.LAScatalog <- function(las, shp, output_dir, attributes = c("prediction", "chance")) {
  # returns LAScatalog with additional attributes from shapes
  # undo previous selections
  opt_select(las) <- "* -t"
  # set paramters
  options <- list(
    need_output_file = TRUE, # output path necessary
    need_buffer = FALSE, # no buffer necessary
    automerge = TRUE # combine outputs
  )
  # create output dir
  check_create_dir(output_dir)
  # execute & return
  output <- catalog_apply(las, merge_shape_attributes.LAScluster,
    shp = shp, attributes = attributes, .options = options
  )
  return(output)
}

################################################################################

plan(multisession, workers = 5L)

# load catalog
ctg_list <- list.files(path_ctg, pattern = "[.]las", full.names = T)
ctg_list <- ctg_list[grepl("area_6", ctg_list)]
ctg <- readTLSLAScatalog(ctg_list)

# load shapes
shp_obj <- st_read(path_shp)
st_crs(shp_obj) <- "+proj=utm +zone=32 +ellps=WGS84 +units=m +vunits=m +no_defs"
shp_obj$prediction <- as.numeric(as.factor(shp_obj$prediction))

# set catalog options
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- paste0(path_out, "/area_6_preds_{ID}")

# add shape variables to catalog tiles
ctg_out <- merge_shape_attributes.LAScatalog(ctg, shp_obj, path_out)
warnings()

# create lax files
lidR:::catalog_laxindex(ctg_out)

################################################################################