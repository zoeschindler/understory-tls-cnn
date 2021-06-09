################################################################################
################################################################################
# RASTER CALCULATION - FUNCTIONS
################################################################################
################################################################################

# line 32: path to cpp-file
# line 53: amount of cores

# load packages
library(lidR)  # for point clouds, also loads sp & raster
library(TreeLS)

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

# fast eigenvalue calculation
# source: https://gis.stackexchange.com/questions/395916/get-eigenvalues-of-large-point-cloud-using-lidr
Rcpp::sourceCpp("D:/Masterarbeit_Zoe/5_Analyse/eigen_decomposition.cpp")

################################################################################

metric_ortho <- function(r, g, b, z) {
  # necessary for raster_ortho
  # returns RGB values of the highest points
  index_max_Z <- which.max(z)
  bands = list(
    red = r[index_max_Z],
    green = g[index_max_Z],
    blue = b[index_max_Z]
  )
  return(bands)
}

################################################################################

add_geometry <- function(las) {
  # necessary for raster_geometry
  # returns geometric features based on eigenvalues
  eigen <- eigen_decomposition(las, 20, 16)  # 20 neighbours, 16 cores
  las <- add_lasattribute(las, eigen[,3]/(eigen[,1] + eigen[,2] + eigen[,3]), "curvature", "curvature")
  las <- add_lasattribute(las, (eigen[,1] - eigen[,2])/eigen[,1], "linearity", "linearity")
  las <- add_lasattribute(las, (eigen[,2] - eigen[,3])/eigen[,1], "planarity", "planarity")
  las <- add_lasattribute(las, eigen[,3]/eigen[,1], "sphericity", "sphericity")
  las <- add_lasattribute(las, (eigen[,1] - eigen[,3])/eigen[,1], "anisotropy", "anisotropy")
  return(las)
}

################################################################################

metric_geometry <- function(curvature, linearity, planarity, sphericity, anisotropy) {
  # necessary for raster_geometry
  # returns statistics of geometric features
  geometries = list(
    curvature_mean=mean(curvature),
    curvature_sd=sd(curvature),
    curvature_max=max(curvature),
    linearity_mean=mean(linearity),
    linearity_sd=sd(linearity),
    linearity_max=max(linearity),
    planarity_mean=mean(planarity),
    planarity_sd=sd(planarity),
    planarity_max=max(planarity),
    sphericity_mean=mean(sphericity),
    sphericity_sd=sd(sphericity),
    sphericity_max=max(sphericity),
    anisotropy_mean=mean(anisotropy),
    anisotropy_sd=sd(anisotropy),
    anisotropy_max=max(anisotropy)
  )
  return(geometries)
}

################################################################################

metric_reflectance <- function(r) {
  # necessary for raster_intensity
  # returns statistics of reflectances
  reflectances = list(
    mean = ifelse(length(r)==0, NA, mean(r)),
    sd = ifelse(length(r)==0, NA, sd(r)),
    max = ifelse(length(r)==0, NA, max(r))
  )
  return(reflectances)
}

################################################################################

rescale_raster <- function(raster) {
  # scales raster values between 0 and 1
  # uses minimum & maximum values of all bands to keep proportions
  print("... rescale raster to values between 0 and 1")
  val_min = min(cellStats(raster, "min"), na.rm=TRUE)
  val_max = max(cellStats(raster, "max"), na.rm=TRUE)
  return((raster-val_min)/(val_max-val_min))
}

################################################################################
# RETILE & FILTER UNDERSTORY POINTS
################################################################################

area_retile_ctg.LAScluster <- function(chunk, area) {
  # return las if extents of las & area overlap
  if (!is.null(intersect(extent(chunk), extent(area)))) {
    las <- readLAS(chunk)
    if (is.empty(las)) return(NULL)
    return(las)
  } else {
    return(NULL)
  }
}

area_retile_ctg.LAScatalog <- function(las, area) {
  # retile the catalog per area
  # undo previous selections
  opt_select(las) <-  "*"
  # set paramters
  options <- list(
    need_output_file = TRUE,  # output path necessary
    need_buffer = FALSE,  # buffer not necessary
    automerge = TRUE)  # combine outputs
  # execute & return
  output  <- catalog_apply(las, area_retile_ctg.LAScluster, area = area, .options = options)
  return(output)
}

################################################################################

normalize_ctg.LAScluster <- function(las) {
  # returns normalized point cloud (LAS file)
  # load the data
  las <- readLAS(las)
  if (is.empty(las)) return(NULL)
  # classify & normalize
  las <- classify_ground(las, csf(class_threshold = 0.3, cloth_resolution = 0.1, sloop_smooth = TRUE))  # has to be small due to smaller areas & heavy slope 
  dtm <- grid_terrain(las, tin(), res = 0.01)  # if bigger, there are many artifacts
  las <- normalize_height(las, dtm, na.rm = T)
  las <- filter_poi(las, Z >= -0.3)
  # delete buffer & return points
  las <- filter_poi(las, buffer == 0)
  gc()  # make RAM space
  return(las)
}

normalize_ctg.LAScatalog <- function(las) {
  # returns normalized point cloud (LAS catalog)
  # undo previous selections
  opt_select(las) <-  "*"
  opt_stop_early(las) <- FALSE  # otherwise it stops when chunks at the border are too small
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

filter_understory_ctg.LAScluster <- function(las, height, remove_stems) {
  # returns point cloud without understory (LAS file)
  # load the data
  las <- readLAS(las)
  if (is.empty(las)) return(NULL)
  # get & remove stem trees
  if (remove_stems) {
    map <- treeMap(las, map.hough())
    las <- treePoints(las, map, trp.crop())
    las <- stemPoints(las, stm.hough(pixel_size = 0.01))
    las <- filter_poi(las, Stem == FALSE)
  }
  # remove everything below certain height
  las <- filter_poi(las, Z <= height)
  if (is.empty(las)) return(NULL)
  # voxel metrics: is there a point in the voxel or not?
  voxels <- voxel_metrics(las, length(X), res=0.1, all_voxels=TRUE)
  voxels$V1 <- ifelse(voxels$V1 > 0, 1, 0)
  voxels$V1[is.na(voxels$V1)] <- 0
  # convert coordinates to [cm]
  # because otherwise R adds decimal places and makes this crash
  voxels$X <- as.integer(round(voxels$X*100))
  voxels$Y <- as.integer(round(voxels$Y*100))
  voxels$Z <- as.integer(round(voxels$Z*100))
  # loop from lowest to highest z value, start at 0.5 m height
  z_loop_vals <- sort(unique(voxels$Z))[9:length(unique(voxels$Z))]
  for (z_val in z_loop_vals) {
    # loop through every non-empty voxel with this z value
    z_loop_vox <- voxels[voxels$Z == z_val & voxels$V1 == 1,]
    if (nrow(z_loop_vox != 0)) {
      for (idx in 1:nrow(z_loop_vox)) {
        vox <- z_loop_vox[idx,]
        # set attribute to empty, if z-1, x+-1, y+-1 is all empty
        neighbours <- voxels$V1[voxels$Z == vox$Z-10 &
                                  voxels$X >= vox$X-10 & voxels$X <= vox$X+10 &
                                  voxels$Y >= vox$Y-10 & voxels$Y <= vox$Y+10]
        if (mean(neighbours) == 0) {
          voxels$V1[voxels$X == vox$X & voxels$Y == vox$Y & voxels$Z == vox$Z] <- 0
        }
      }
    }
  }
  # convert coordinates back to [m]
  voxels$X <- voxels$X/100
  voxels$Y <- voxels$Y/100
  voxels$Z <- voxels$Z/100
  # add voxel attributes to the points
  las <- add_lasattribute(las, 1, "V1", "keep voxels with 1")  # create empty attribute
  # save points which should remain the same
  unchanged_las <- filter_poi(las, Z <= min(z_loop_vals-5)/100)
  # for each (filtered) vertical voxel layer, create a raster
  for (z_val in z_loop_vals) {
    # create raster
    z_subset <- voxels[as.integer(voxels$Z*100)==z_val,]
    z_subset <- as.data.frame(z_subset)[,c(1,2,4)]
    new_raster <- rasterFromXYZ(z_subset)
    crs(new_raster) <- crs(las)
    # add raster values to point cloud
    las_z <- filter_poi(las, Z > ((z_val-5)/100) & Z <= ((z_val+5)/100))
    las_z <- merge_spatial(las_z, new_raster, "V1")
    # remove points with V1 == 0
    las_z <- filter_poi(las_z, V1 == 1)
    unchanged_las <- rbind(unchanged_las, las_z)
  } 
  # delete all points with attribute empty
  las <- unchanged_las
  rm(unchanged_las); gc()
  # delete buffer & return points
  las <- filter_poi(las, buffer == 0)
  return(las)
}

filter_understory_ctg.LAScatalog <- function(las, height=2, remove_stems=TRUE) {
  # returns point cloud without understory (LAS catalog)
  # undo previous selections
  opt_select(las) <-  "*"
  # set paramters
  options <- list(
    need_output_file = TRUE,  # output path necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = TRUE)  # combine outputs
  # execute & return
  output  <- catalog_apply(las, filter_understory_ctg.LAScluster, height = height,
                           remove_stems = remove_stems, .options = options)
  return(output)
}

################################################################################
# SINGLE RASTER FUNCTIONS - LAS FILES
################################################################################

raster_nDSM <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns nDSM raster
  check_create_dir(output_dir)
  print("... creating nDSM")
  # create chm / nDSM
  chm <- grid_metrics(point_cloud, ~max(Z), res = resolution)
  # rescale raster
  if (rescale) {
    chm <- rescale_raster(chm)
  }
  if (saving) {
    # save chm / nDSM
    writeRaster(chm, paste0(output_dir, "/nDSM_", output_name, "_",
                            resolution*100, "cm.tif"), overwrite = TRUE)
  } else {
    # return the raster
    return(chm)
  }
}

################################################################################

raster_ortho <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns orthomosaic raster
  check_create_dir(output_dir)
  print("... creating orthomosaic")
  # create orthomosaic
  ortho <- grid_metrics(point_cloud, ~metric_ortho(R,G,B,Z), res = resolution)
  # rescale raster
  if (rescale) {
    ortho <- rescale_raster(ortho)
  }
  if (saving) {
    # save orthomosaic
    writeRaster(ortho, paste0(output_dir, "/ortho_", output_name, "_",
                              resolution*100, "cm.tif"), overwrite = TRUE)
  } else {
    # return the raster
    return(ortho)
  }
}

################################################################################

raster_geometry <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns geometry rasters
  check_create_dir(output_dir)
  print("... creating geometry rasters")
  # get attributes
  point_cloud <- add_geometry(point_cloud)
  print("... eigenvalues were calculated")
  # create raster
  geom_raster <- grid_metrics(point_cloud, ~metric_geometry(curvature, linearity, planarity, sphericity, anisotropy), res = resolution)
  # create empty raster list
  raster_bands <- c()
  # loop through bands
  for (i in 1:dim(geom_raster)[3]) {
    # extract single bands
    raster_band <- geom_raster[[i]]
    type <- names(raster_band)
    # rescale raster
    if (rescale) {
      raster_band <- rescale_raster(raster_band)
    }
    if (saving) {
      # save rasters
      writeRaster(raster_band, paste0(output_dir, "/", type, "_", output_name,
                                      "_", resolution*100, "cm.tif"), overwrite = TRUE)
    } else {
      # save raster to list
      raster_bands <- c(raster_bands, raster_band)
    }
  }
  if (!saving) {
    # return raster stack
    return(brick(stack(raster_bands)))
  }
}

################################################################################

raster_reflectance <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns intensity rasters
  check_create_dir(output_dir)
  print("... creating reflectance rasters")
  # create intensity raster
  reflect_raster <- grid_metrics(point_cloud, ~metric_reflectance(Reflectance), res = resolution)
  # create empty raster list
  raster_bands <- c()
  # loop through bands
  for (i in 1:dim(reflect_raster)[3]) {
    # extract single bands
    raster_band <- reflect_raster[[i]]
    type <- names(raster_band)
    # rescale raster
    if (rescale) {
      raster_band <- rescale_raster(raster_band)
    }
    if (saving) {
      # save raster
      writeRaster(raster_band, paste0(output_dir, "/reflectance_", type, "_",
                                      output_name, "_", resolution*100, "cm.tif"), overwrite = TRUE)
    } else {
      # save raster to list
      raster_bands <- c(raster_bands, raster_band)
    }
  }
  if (!saving) {
    # return raster stack
    return(brick(stack(raster_bands)))
  }
}

################################################################################

raster_point_density <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE){
  # saves or returns point density raster
  check_create_dir(output_dir)
  print("... creating point density raster")
  # create point density raster
  point_density <- grid_density(point_cloud, res = resolution)
  # rescale raster
  if (rescale) {
    point_density <- rescale_raster(point_density)
  }
  if (saving) {
    # save point density raster
    writeRaster(point_density, paste0(output_dir, "/point_density_", output_name, "_",
                                      resolution*100, "cm.tif"), overwrite = TRUE)
  } else {
    # return the raster
    return(point_density)
  }
}

################################################################################
# SINGLE RASTER FUNCTIONS - LAS CATALOGS
################################################################################

raster_nDSM_ctg.LAScluster <- function(chunk, resolution, output_dir, output_name, rescale, saving) {
  # saves or returns nDSM raster (LAS file)
  # load the data
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  # create nDSM raster
  nDSM <- raster_nDSM(las, resolution, output_dir, output_name, rescale, saving)
  # delete buffer & return points
  nDSM <- crop(nDSM, extent(chunk))
  return(nDSM)
}

raster_nDSM_ctg.LAScatalog <- function(las, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns nDSM raster (LAS catalog)
  # undo previous selections
  opt_select(las) <-  "xyz"
  # delete output path
  opt_output_files(las) <- ""
  # set parameters
  options <- list(
    need_output_file = FALSE,  # output path not necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = TRUE,  # combine outputs
    raster_alignment = resolution)  # align chunks & rasters
  # calculate & merge raster
  output  <- catalog_apply(las, raster_nDSM_ctg.LAScluster, resolution = resolution,
                           output_dir = output_dir, output_name = output_name,
                           rescale = FALSE, saving = FALSE, .options = options)
  # rescale raster
  if (rescale) {
    output <- rescale_raster(output)
  }
  if (saving) {
    # save merged output
    writeRaster(output, paste0(output_dir, "/nDSM_", output_name, "_",
                               resolution*100, "cm.tif"), overwrite = TRUE)
  } else {
    # return the raster
    return(output)
  }
}

################################################################################

raster_ortho_ctg.LAScluster <- function(chunk, resolution, output_dir, output_name, rescale, saving) {
  # saves or returns orthomosaic raster (LAS file)
  # load the data
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  # create ortho raster
  ortho <- raster_ortho(las, resolution, output_dir, output_name, rescale, saving)
  # delete buffer & return points
  ortho <- crop(ortho, extent(chunk))
  return(ortho)
}

raster_ortho_ctg.LAScatalog <- function(las, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns orthomosaic raster (LAS catalog)
  # undo previous selections
  opt_select(las) <-  "*"
  # delete output path
  opt_output_files(las) <- ""
  # set parameters
  options <- list(
    need_output_file = FALSE,  # output path not necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = TRUE,  # combine outputs
    raster_alignment = resolution)  # align chunks & rasters
  # calculate & merge raster
  output  <- catalog_apply(las, raster_ortho_ctg.LAScluster, resolution = resolution,
                           output_dir = output_dir, output_name = output_name,
                           rescale = FALSE, saving = FALSE, .options = options)
  # rescale all raster bands at once to keep relations
  if (rescale) {
    output <- rescale_raster(output)
  }
  if (saving) {
    # save merged output
    writeRaster(output, paste0(output_dir, "/ortho_", output_name, "_",
                               resolution*100, "cm.tif"), overwrite = TRUE)
  } else {
    # return the raster
    return(output)
  }
}

################################################################################

raster_geometry_ctg.LAScluster <- function(chunk, resolution, output_dir, output_name, rescale, saving) {
  # saves or returns geometry raster stack (LAS file)
  # load the data
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  # create geometry rasters
  geometries <- raster_geometry(las, resolution, output_dir, output_name, rescale, saving)
  # delete buffer
  geometries <- crop(geometries, extent(chunk))
  # clear memory
  gc()
  return(geometries)
}

raster_geometry_ctg.LAScatalog <- function(las, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns geometry rasters (LAS catalog)
  # undo previous selections
  opt_select(las) <-  "xyz"
  # allow overwriting
  las@output_options$drivers$Raster$param$overwrite <- TRUE
  # set temporary output path
  temp_dir <- paste0(output_dir, "/temp")
  check_create_dir(temp_dir)
  opt_output_files(las) <- paste0(temp_dir, "/temp_geometry_", output_name, "_{ID}")
  # set parameters
  options <- list(
    need_output_file = TRUE,  # output path necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = FALSE,  # don't combine outputs
    raster_alignment = resolution)  # align chunks & rasters
  # calculate & merge raster
  output  <- catalog_apply(las, raster_geometry_ctg.LAScluster, resolution = resolution,
                           output_dir = output_dir, output_name = output_name,
                           rescale = FALSE, saving = FALSE, .options = options)
  print("... tiles saved")
  # load & merge tiles
  raster_list <- lapply(output, stack)
  raster_list$filename <- paste0(temp_dir, "/geometry_all_", output_name, "_", resolution*100, "cm.tif")
  merged_raster <- do.call(merge, raster_list)
  # get layer names
  names <- names(metric_geometry(0,0,0,0,0))
  # loop through bands
  for (idx in 1:dim(merged_raster)[3]) {
    raster_band <- merged_raster[[idx]]
    type <- names[idx]
    # rescale band
    if (rescale) {
      raster_band <- rescale_raster(raster_band)
    }
    if (saving) {
      # save each band separately
      writeRaster(raster_band, paste0(output_dir, "/", type, "_", output_name, "_",
                                      resolution*100, "cm.tif"), overwrite = TRUE)
    } else {
      # replace with rescaled band
      merged_raster[[idx]] <- raster_band
    }
  }
  # delete temporary folder
  unlink(temp_dir, recursive = TRUE)
  # return bands stacked
  if (!saving) {
    return(merged_raster)
  }
}

################################################################################

raster_reflectance_ctg.LAScluster <- function(chunk, resolution, output_dir, output_name, rescale, saving) {
  # saves or returns reflectance raster stack (LAS file)
  # load the data
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  # create reflectance raster
  reflectance <- raster_reflectance(las, resolution, output_dir, output_name, rescale, saving)
  # delete buffer
  reflectance <- crop(reflectance, extent(chunk))
  # clear memory
  gc()
  return(reflectance)
}

raster_reflectance_ctg.LAScatalog <- function(las, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns geometry rasters (LAS catalog)
  # undo previous selections
  opt_select(las) <-  "*"
  # allow overwriting
  las@output_options$drivers$Raster$param$overwrite <- TRUE
  # set temporary output path
  temp_dir <- paste0(output_dir, "/temp")
  check_create_dir(temp_dir)
  opt_output_files(las) <- paste0(temp_dir, "/temp_reflectance_", output_name, "_{ID}")
  # set parameters
  options <- list(
    need_output_file = TRUE,  # output path necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = FALSE,  # don't combine outputs
    raster_alignment = resolution)  # align chunks & rasters
  # calculate & merge raster
  output  <- catalog_apply(las, raster_reflectance_ctg.LAScluster, resolution = resolution,
                           output_dir = output_dir, output_name = output_name,
                           rescale = FALSE, saving = FALSE, .options = options)
  print("... tiles saved")
  # load & merge tiles
  raster_list <- lapply(output, stack)
  raster_list$filename <- paste0(temp_dir, "/reflectance_all_", output_name, "_", resolution*100, "cm.tif")
  merged_raster <- do.call(merge, raster_list)
  # get layer names
  names <- names(metric_reflectance(0))
  # loop through bands
  for (idx in 1:dim(merged_raster)[3]) {
    raster_band <- merged_raster[[idx]]
    type <- names[idx]
    # rescale band
    if (rescale) {
      raster_band <- rescale_raster(raster_band)
    }
    if (saving) {
      # save each band separately
      writeRaster(raster_band, paste0(output_dir, "/reflectance_", type, "_",
                                      output_name, "_", resolution*100, "cm.tif"), overwrite = TRUE)
    } else {
      # replace with rescaled band
      merged_raster[[idx]] <- raster_band
    }
  }
  # delete temporary folder
  unlink(temp_dir, recursive = TRUE)
  # return bands stacked
  if (!saving) {
    return(merged_raster)
  }
}

################################################################################

raster_point_density_ctg.LAScluster <- function(chunk, resolution, output_dir, output_name, rescale, saving) {
  # saves or returns point density raster (LAS file)
  # load the data
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  # create point_density raster
  point_density <- raster_point_density(las, resolution, output_dir, output_name, rescale, saving)
  # delete buffer & return points
  point_density <- crop(point_density, extent(chunk))
  return(point_density)
}

raster_point_density_ctg.LAScatalog <- function(las, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # saves or returns point density raster (LAS catalog)
  # undo previous selections
  opt_select(las) <-  "xyz"
  # delete output path
  opt_output_files(las) <- ""
  # set parameters
  options <- list(
    need_output_file = FALSE,  # output path not necessary
    need_buffer = TRUE,  # buffer necessary
    automerge = TRUE,  # combine outputs
    raster_alignment = resolution)  # align chunks & rasters
  # calculate & merge raster
  output  <- catalog_apply(las, raster_point_density_ctg.LAScluster, resolution = resolution,
                           output_dir = output_dir, output_name = output_name,
                           rescale = FALSE, saving = FALSE, .options = options)
  # rescale raster
  if (rescale) {
    output <- rescale_raster(output)
  }
  if (saving) {
    # save merged output
    writeRaster(output, paste0(output_dir, "/point_density_", output_name, "_",
                               resolution*100, "cm.tif"), overwrite = TRUE)
  } else {
    # return the raster
    return(output)
  }
}

################################################################################
# CREATE ALL RASTERS - LAS FILES
################################################################################

raster_create_all <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # creates all rasters
  # set rescale = TRUE if it should be normalized between 0 and 1
  # set saving = TRUE if raster should be saved, set saving = FALSE to return raster as object
  raster_nDSM(point_cloud, resolution, paste0(output_dir, "/nDSM"), output_name, rescale, saving)
  raster_ortho(point_cloud, resolution, paste0(output_dir, "/ortho"), output_name, rescale, saving)
  raster_geometry(point_cloud, resolution, paste0(output_dir, "/geometry"), output_name, rescale, saving)
  raster_reflectance(point_cloud, resolution, paste0(output_dir, "/reflectance"), output_name, rescale, saving)
  raster_point_density(point_cloud, resolution, paste0(output_dir, "/point_density"), output_name, rescale, saving)
  print("done!")
}

################################################################################
# CREATE ALL RASTERS - LAS CATALOGS
################################################################################

raster_create_all_ctg <- function(ctg, resolution, output_dir, output_name, rescale=TRUE, saving=TRUE) {
  # creates all rasters
  # set rescale = TRUE if it should be normalized between 0 and 1
  # set saving = TRUE if raster should be saved, set saving = FALSE to return raster as object
  raster_nDSM_ctg.LAScatalog(ctg, resolution, paste0(output_dir, "/nDSM"), output_name, rescale, saving)
  gc()
  raster_ortho_ctg.LAScatalog(ctg, resolution, paste0(output_dir, "/ortho"), output_name, rescale, saving)
  gc()
  raster_point_density_ctg.LAScatalog(ctg, resolution, paste0(output_dir, "/point_density"), output_name, rescale, saving)
  gc()
  raster_reflectance_ctg.LAScatalog(ctg, resolution, paste0(output_dir, "/reflectance"), output_name, rescale, saving)
  gc()
  raster_geometry_ctg.LAScatalog(ctg, resolution, paste0(output_dir, "/geometry"), output_name, rescale, saving)
  gc()
  print("done!")
}

################################################################################
