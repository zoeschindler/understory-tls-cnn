################################################################################
################################################################################
# PREPROCESSING
################################################################################
################################################################################

# load packages
library(dplyr)
library(lidR)
library(sf)

# set path
points_path  <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/points"
kml_path    <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/scan_positions/2021-05-03 wino np earth.kml"
setwd(points_path)

# *.asc to *.las in CC
# *.kmz to *.kml in Google Earth

################################################################################
# EXTRACT GROUND - 2 M
################################################################################

prepare_scan <- function(path) {
  scan <- readTLSLAS(path)
  scan <- scan %>%
    # classify ground
    classify_ground(., csf(class_threshold = 2.0, cloth_resolution = 0.3)) %>%
    filter_ground(.) %>%
    # rename scalar fields
    add_lasattribute(., .@data$`Scalar field`, "Reflectance", "Reflectance") %>%
    remove_lasattribute(., "Scalar field") %>%
    add_lasattribute(., .@data$`Scalar field #2`, "Deviation", "Deviation") %>%
    remove_lasattribute(., "Scalar field #2") %>%
    add_lasattribute(., .@data$`Scalar field #3`, "Target Index", "Target Index") %>%
    remove_lasattribute(., "Scalar field #3") %>%
    add_lasattribute(., .@data$`Scalar field #4`, "Target Count", "Target Count") %>%
    remove_lasattribute(., "Scalar field #4") %>%
    add_lasattribute(., .@data$`Scalar field #5`, "Source Indicator", "Source Indicator") %>%
    remove_lasattribute(., "Scalar field #5") %>%
    # save
    writeLAS(., paste0(substr(path, 1, nchar(path)-4), "_ground.las"))
  gc()
  return(paste0(substr(path, 1, nchar(path)-4), "_ground.las"))
}

scan2_path  <- prepare_scan("ScanPos002 - SINGLESCANS - 210503_102933 - Cloud.las")
scan3_path  <- prepare_scan("ScanPos003 - SINGLESCANS - 210503_103110 - Cloud.las")
scan4_path  <- prepare_scan("ScanPos004 - SINGLESCANS - 210503_103257 - Cloud.las")
scan5_path  <- prepare_scan("ScanPos005 - SINGLESCANS - 210503_103526 - Cloud.las")
scan6_path  <- prepare_scan("ScanPos006 - SINGLESCANS - 210503_103753 - Cloud.las")

################################################################################
# FILTERING NOISE
################################################################################

filter_noise <- function(path, res) {
  scan <- readTLSLAS(path)
  scan <- scan %>%
    # classify noise
    classify_noise(., ivf(res = res, n = 5)) %>%
    # delete noise
    filter_poi(., Classification != LASNOISE) %>%
    # save
    writeLAS(., paste0(substr(path, 1, nchar(path)-4), "_noise15cm.las"))
  gc()
  return(paste0(substr(path, 1, nchar(path)-4), "_noise15cm.las"))
}

resolution <- 0.1
scan2_path_noise  <- filter_noise(scan2_path, resolution)
scan3_path_noise  <- filter_noise(scan3_path, resolution)
scan4_path_noise  <- filter_noise(scan4_path, resolution)
scan5_path_noise  <- filter_noise(scan5_path, resolution)
scan6_path_noise  <- filter_noise(scan6_path, resolution)

################################################################################
# CLIP EXAMPLES
################################################################################

# clip_scan_circle <- function(path_scan, path_vegetation, vegetation_id, radius) {
#   if (!dir.exists("clips")) {
#     print("... creating new folder")
#     dir.create("clips")
#   } else {
#     print("... using existing folder")
#   }
#   # read as kml + transform CRS
#   vegetation <- sf::st_read(path_vegetation)
#   vegetation <- sf::st_transform(vegetation, 25832)
#   # get plot coordinates
#   vegetation_plot <- as.data.frame(sf::st_coordinates(vegetation[vegetation_id,]))
#   # clip points within radius
#   scan <- readTLSLAS(path_scan) %>%
#     clip_circle(., xcenter=vegetation_plot$X, ycenter=vegetation_plot$Y, radius = radius) %>%
#     writeLAS(., paste0("clips/",substr(path_scan, 1, 10), "_", vegetation_id, "_", vegetation$Name[vegetation_id], ".las"))
#   gc()
#   return(paste0("clips/", substr(path_scan, 1, 10), "_", vegetation_id, "_", vegetation$Name[vegetation_id], ".las"))
# }
# 
# path_vegetation <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/vegetation/Export_ODK_clean_2D.kml"
# path_scan <- c("ScanPos002 - SINGLESCANS - 210503_102933 - Cloud_ground_noise15cm.las",
#                "ScanPos003 - SINGLESCANS - 210503_103110 - Cloud_ground_noise15cm.las",
#                "ScanPos003 - SINGLESCANS - 210503_103110 - Cloud_ground_noise15cm.las",
#                "ScanPos004 - SINGLESCANS - 210503_103257 - Cloud_ground_noise15cm.las",
#                "ScanPos004 - SINGLESCANS - 210503_103257 - Cloud_ground_noise15cm.las",
#                "ScanPos005 - SINGLESCANS - 210503_103526 - Cloud_ground_noise15cm.las",
#                "ScanPos005 - SINGLESCANS - 210503_103526 - Cloud_ground_noise15cm.las",
#                "ScanPos006 - SINGLESCANS - 210503_103753 - Cloud_ground_noise15cm.las",
#                "ScanPos028 - SINGLESCANS - 210503_112304 - Cloud_ground_noise15cm.las",
#                "ScanPos028 - SINGLESCANS - 210503_112304 - Cloud_ground_noise15cm.las",
#                "ScanPos028 - SINGLESCANS - 210503_112304 - Cloud_ground_noise15cm.las",
#                "ScanPos029 - SINGLESCANS - 210503_112505 - Cloud_ground_noise15cm.las",
#                "ScanPos029 - SINGLESCANS - 210503_112505 - Cloud_ground_noise15cm.las",
#                "ScanPos029 - SINGLESCANS - 210503_112505 - Cloud_ground_noise15cm.las")
# vegetation_id <- c(613, #Scan002
#                    594, 593, #Scan003
#                    616, 617, #Scan004
#                    598, 597, #Scan005
#                    599, #Scan006
#                    670, 690, 689, #Scan028
#                    688, 687, 672) #Scan029
# radius <- 1
# 
# for (i in 1:length(path_scan)) {
#   clip_scan_circle(path_scan[i], path_vegetation, vegetation_id[i], radius)
# }

################################################################################
# FILTER REFLECTANCE & DEVIATION
################################################################################

filter_ref_dev <- function(path, reflectance, deviation) {
  scan <- readTLSLAS(path)
  scan <- scan %>%
    # filter with reflectance and deviation
    filter_poi(., Reflectance >= reflectance &  Reflectance < 0 & Deviation <= deviation & Deviation > 0) %>%
    # save
    writeLAS(., paste0(substr(path, 1, nchar(path)-4), "_ref", reflectance, "_dev", deviation, ".las"))
  gc()
  return(paste0(substr(path, 1, nchar(path)-4), "_ref", reflectance, "_dev", deviation, ".las"))
}

reflectance = -20
deviation = 20
scan2_path_noise_filter  <- filter_ref_dev(scan2_path_noise, reflectance, deviation)
scan3_path_noise_filter  <- filter_ref_dev(scan3_path_noise, reflectance, deviation)
scan4_path_noise_filter  <- filter_ref_dev(scan4_path_noise, reflectance, deviation)
scan5_path_noise_filter  <- filter_ref_dev(scan5_path_noise, reflectance, deviation)
scan6_path_noise_filter  <- filter_ref_dev(scan6_path_noise, reflectance, deviation)

################################################################################
# MAXIMUM DISTANCE
################################################################################

clip_scan_circle <- function(path_scan, path_kml, radius) {
  # read las
  scan_las <- readTLSLAS(path_scan)
  # read as kml + transform CRS
  scan_positions <- sf::st_read(path_kml)
  scan_positions_trans <- sf::st_transform(scan_positions, 25832)
  # get scanner coordinates
  scan_name <- substr(path_scan, 1, 10)
  scan_pos <- as.data.frame(sf::st_coordinates(scan_positions_trans[scan_positions_trans$Name == scan_name,]))
  # clip points within radius
  scan <- readTLSLAS(path_scan) %>%
    clip_circle(., xcenter=scan_pos$X, ycenter=scan_pos$Y, radius = radius) %>%
    # save
    writeLAS(., paste0(substr(path_scan, 1, nchar(path_scan)-4), "_circle", radius, ".las"))
  gc()
  return(paste0(substr(path_scan, 1, nchar(path_scan)-4), "_circle", radius, ".las"))
}

radius <- 30
scan2_path_noise_circle  <- clip_scan_circle(scan2_path_noise_filter, kml_path, radius)
scan3_path_noise_circle  <- clip_scan_circle(scan3_path_noise_filter, kml_path, radius)
scan4_path_noise_circle  <- clip_scan_circle(scan4_path_noise_filter, kml_path, radius)
scan5_path_noise_circle  <- clip_scan_circle(scan5_path_noise_filter, kml_path, radius)
scan6_path_noise_circle  <- clip_scan_circle(scan6_path_noise_filter, kml_path, radius)


################################################################################
# MERGE SCANS
################################################################################

random_per_voxel <- function(X, Y, Z, Reflectance, Deviation, T_Index, T_Count) {
  id = floor(runif(n=1, min=1, max=length(Reflectance)+1))
  values = list("original_X"= X[id],
                "original_Y"= Y[id],
                "original_Z"= Z[id],
                "Reflectance"=Reflectance[id],
                "Deviation"=Deviation[id],
                "Target Index"=T_Index[id],
                "Target Count"=T_Count[id])
  return(values)
}

mean_per_voxel <- function(Reflectance, Deviation) {
  values = list("Reflectance"=mean(Reflectance),
                "Deviation"=mean(Deviation))
  return(values)
}

min_dev_per_voxel <- function(X, Y, Z, Reflectance, Deviation, T_Index, T_Count) {
  id = which.min(Deviation)
  values = list("original_X"= X[id],
                "original_Y"= Y[id],
                "original_Z"= Z[id],
                "Reflectance"=Reflectance[id],
                "Deviation"=Deviation[id],
                "Target Index"=T_Index[id],
                "Target Count"=T_Count[id])
  return(values)
}

merge_and_thin <- function(scan_paths, output_name, resolution=0.01, method=c("random","mean","min_deviation")) {
  # merge point clouds
  combo <- readTLSLAS(scan_paths[1])
  for (i in 2:length(scan_paths)) {
    combo <- rbind(combo, readTLSLAS(scan_paths[i]))
  }
  # select method
  if (method=="random") {
    # method 1: random pixel per voxel, keep coordinates 
    random_data <- voxel_metrics(combo, ~random_per_voxel(X, Y, Z, Reflectance, Deviation, `Target Index`, `Target Count`), resolution)
    random_las <- LAS(data=data.frame(X=random_data$original_X, Y=random_data$original_Y, Z=random_data$original_Z))
    random_las <- add_lasattribute(random_las, random_data$Reflectance, "Reflectance", "Reflectance")
    random_las <- add_lasattribute(random_las, random_data$Deviation, "Deviation", "Deviation")
    random_las <- add_lasattribute(random_las, random_data$`Target Index`, "Target Index", "Target Index")
    random_las <- add_lasattribute(random_las, random_data$`Target Count`, "Target Count", "Target Count")
    writeLAS(random_las, paste0("thinning_random_", output_name, ".las"))
    return(paste0("thinning_random_", output_name, ".las"))
  }
  # select method
  if (method=="mean") {
    # method 2: mean per voxel, assign voxel coordinates
    mean_data <- voxel_metrics(combo, ~mean_per_voxel(Reflectance, Deviation), resolution)
    mean_las <- LAS(data=data.frame(X=mean_data$X, Y=mean_data$Y, Z=mean_data$Z))
    mean_las <- add_lasattribute(mean_las, mean_data$Reflectance, "Reflectance", "Reflectance")
    mean_las <- add_lasattribute(mean_las, mean_data$Deviation, "Deviation", "Deviation")
    writeLAS(mean_las, paste0("thinning_mean_", output_name, ".las"))
    return(paste0("thinning_mean_", output_name, ".las"))
  }
  # select method
  if (method=="min_deviation") {
    # method 3: lowest deviation pixel per voxel, keep coordinates
    min_dev_data <- voxel_metrics(combo, ~min_dev_per_voxel(X, Y, Z, Reflectance, Deviation, `Target Index`, `Target Count`), resolution)
    min_dev_las <- LAS(data=data.frame(X=min_dev_data$original_X, Y=min_dev_data$original_Y, Z=min_dev_data$original_Z))
    min_dev_las <- add_lasattribute(min_dev_las, min_dev_data$Reflectance, "Reflectance", "Reflectance")
    min_dev_las <- add_lasattribute(min_dev_las, min_dev_data$Deviation, "Deviation", "Deviation")
    min_dev_las <- add_lasattribute(min_dev_las, min_dev_data$`Target Index`, "Target Index", "Target Index")
    min_dev_las <- add_lasattribute(min_dev_las, min_dev_data$`Target Count`, "Target Count", "Target Count")
    writeLAS(min_dev_las, paste0("thinning_min_deviation_", output_name, ".las"))
    return(paste0("thinning_min_deviation_", output_name, ".las"))
  }
  # if invalid method
  stop("method invalid")
}

circle <- 10
resolution <- 0.005
scan_paths <- c(scan2_path_noise_circle,
                scan3_path_noise_circle,
                scan4_path_noise_circle,
                scan5_path_noise_circle,
                scan6_path_noise_circle)

merge_and_thin(scan_paths, paste0("cirle", circle, "_res", resolution), resolution, method="random")
# merge_and_thin(scan_paths, paste0("cirle", circle, "_res", resolution), resolution, method="mean")
# merge_and_thin(scan_paths, paste0("cirle", circle, "_res", resolution), resolution, method="min_deviation")

################################################################################