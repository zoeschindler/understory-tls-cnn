################################################################################
################################################################################
# RASTER CALCULATION
################################################################################
################################################################################

# load packages
library(lidR)  # for point clouds, also loads sp & raster

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
Rcpp::sourceCpp("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/eigen_decomposition.cpp")

################################################################################

metric_ortho <- function(r, g, b, z) {
  # necessary for raster_ortho
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
  eigen <- eigen_decomposition(las, 20, 6)  # 20 neighbours, 6 cores
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

get_geometry_dict <- function(){
  # necessary for raster_geometry
  # needed to name the files properly
  return(list(c("curvature", "mean"),  c("curvature", "sd"),  c("curvature", "max"),
              c("linearity", "mean"),  c("linearity", "sd"),  c("linearity", "max"),
              c("planarity", "mean"),  c("planarity", "sd"),  c("planarity", "max"),
              c("sphericity", "mean"), c("sphericity", "sd"), c("sphericity", "max"),
              c("anisotropy", "mean"), c("anisotropy", "sd"), c("anisotropy", "max")))
}

################################################################################

metric_reflectance <- function(r) {
  # necessary for raster_intensity
  reflectances = list(
    mean = mean(r),
    median = median(r),
    sd = sd(r),
    max = max(r)
  )
  return(reflectances)
}

################################################################################

get_reflectance_dict <- function(){
  # necessary for raster_geometry
  # needed to name the files properly
  return(c("mean", "median", "sd", "max"))
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
# SINGLE RASTER FUNCTIONS
################################################################################

raster_nDSM <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE) {
  # takes in point cloud (normalized or unnormalized)
  # saves nDSM in output folder
  check_create_dir(output_dir)
  print("... creating nDSM")
  # classify ground
  point_cloud <- classify_ground(point_cloud, csf(class_threshold = 0.1, cloth_resolution = 1)) 
  # create dtm
  dtm <- grid_terrain(filter_ground(point_cloud), tin(), res = resolution)
  # normalize point cloud
  point_cloud_norm <- normalize_height(point_cloud, dtm, na.rm = T)
  # create chm / nDSM
  chm <- grid_metrics(point_cloud_norm, ~max(Z), res = resolution) # no smoothing at all
  # chm <- grid_canopy(point_cloud_norm, res = resolution, p2r())
  # rescale raster
  if (rescale) {
    chm <- rescale_raster(chm)
  }
  # save chm / nDSM
  writeRaster(chm, paste0(output_dir, "/nDSM_", output_name, "_",
                          resolution*100, "cm.tif"), overwrite = TRUE)
}

################################################################################

raster_ortho <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE) {
  # takes in point cloud (normalized or unnormalized)
  # saves pseudo-orthophoto in output folder
  check_create_dir(output_dir)
  print("... creating orthomosaic")
  # create pseudo-orthophoto
  ortho <- grid_metrics(point_cloud, ~metric_ortho(R,G,B,Z), res = resolution)
  # rescale raster
  if (rescale) {
    ortho <- rescale_raster(ortho)
  }
  # save pseudo-orthophoto
  writeRaster(ortho, paste0(output_dir, "/ortho_", output_name, "_",
                            resolution*100, "cm.tif"), overwrite = TRUE)
}

################################################################################

raster_geometry <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE) {
  # takes in point cloud (normalized or unnormalized)
  # saves several geometry rasters in output folder
  check_create_dir(output_dir)
  print("... creating geometry rasters")
  # get attributes
  point_cloud <- add_geometry(point_cloud)
  dict <- get_geometry_dict()
  print("... eigenvalues were calculated")
  # create raster
  geom_raster <- grid_metrics(point_cloud, ~metric_geometry(curvature, linearity, planarity, sphericity, anisotropy), res = resolution)
  # loop through bands
  for (i in 1:dim(geom_raster)[3]) {
    # extract single bands
    name <- dict[[i]][1]
    type <- dict[[i]][2]
    raster_band <- geom_raster[[i]]
    # rescale raster
    if (rescale) {
      raster_band <- rescale_raster(raster_band)
    }
    # save rasters
    writeRaster(raster_band, paste0(output_dir, "/", name, "_", type, "_",
                                    output_name, "_", resolution*100, "cm.tif"), overwrite = TRUE)
  }
}

################################################################################

raster_reflectance <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE){
  # takes in point cloud (normalized or unnormalized)
  # saves intensity raster in output folder
  check_create_dir(output_dir)
  print("... creating reflectance rasters")
  # create intensity raster
  dict <- get_reflectance_dict()
  reflect_raster <- grid_metrics(point_cloud, ~metric_reflectance(Reflectance), res = resolution)
  # loop through bands
  for (i in 1:dim(reflect_raster)[3]) {
    # extract single bands
    type <- dict[i]
    raster_band <- reflect_raster[[i]]
    # rescale raster
    if (rescale) {
      raster_band <- rescale_raster(raster_band)
    }
    # save raster
    writeRaster(raster_band, paste0(output_dir, "/reflectance_", type, "_",
                                    output_name, "_", resolution*100, "cm.tif"), overwrite = TRUE)
  }
}

################################################################################

raster_point_density <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE){
  # takes in point cloud (normalized or unnormalized)
  # saves point density raster in output folder
  check_create_dir(output_dir)
  print("... creating point density raster")
  # create point density raster
  point_density <- grid_density(point_cloud, res = resolution)
  # rescale raster
  if (rescale) {
    point_density <- rescale_raster(point_density)
  }
  # save point density raster
  writeRaster(point_density, paste0(output_dir, "/point_density_", output_name, "_",
                                    resolution*100, "cm.tif"), overwrite = TRUE)
}

################################################################################
# CREATE ALL RASTERS
################################################################################

raster_create_all <- function(point_cloud, resolution, raster_dir, output_name) {
  raster_nDSM(point_cloud, resolution, paste0(raster_dir, "/nDSM"), output_name)
  raster_ortho(point_cloud, resolution, paste0(raster_dir, "/ortho"), output_name)
  raster_geometry(point_cloud, resolution, paste0(raster_dir, "/geometry"), output_name)
  raster_reflectance(point_cloud, resolution, paste0(raster_dir, "/reflectance"), output_name)
  raster_point_density(point_cloud, resolution, paste0(raster_dir, "/point_density"), output_name)
  print("done!")
}

################################################################################