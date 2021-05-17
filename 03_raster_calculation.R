################################################################################
################################################################################
# RASTER CALCULATION
################################################################################
################################################################################

# load packages
library(lidR)  # for point clouds, also loads sp & raster
library(Hmisc)  # for cluster analysis

# set paths
#path_points   <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/thinning_mean_cirle10.las"
path_points   <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/_Testwolke/OT8cm_small.las"
path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"
path_vegetation <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/vegetation/Export_ODK_clean_2D.kml"

# load data
cloud_raw <- readTLSLAS(path_points) # TODO passt eine Wolke komplett in RAM?

################################################################################
# POINT CLOUD QUALITY CHECKS
################################################################################

# what are common quality checks for TLS data?

# save changed data

################################################################################
# FILTER UNDERSTORY POINTS
################################################################################

# set parameters
resolution_dtm = 0.1

# normalize height
cloud_raw <- classify_ground(cloud_raw, csf(class_threshold = 0.1, cloth_resolution = 1))
dtm <- grid_terrain(filter_ground(cloud_raw), tin(), res = resolution_dtm)
cloud_norm <- normalize_height(cloud_raw, dtm, na.rm = T)

# delete below ground points
cloud_norm <- filter_poi(cloud_norm, Z > 0)

###

# delete trees

# https://github.com/Jean-Romain/lidR/wiki/Segment-individual-trees-and-compute-metrics#segment-the-trees
  # sieht mies aus. lieber nicht.

# https://jean-romain.github.io/lidRbook/itd-its.html#its-cloud
  # dauert mir viel zu lange, da hab ich keine Geduld für

# https://gis.stackexchange.com/questions/376772/removing-tree-stems-from-a-point-cloud-in-r
# https://github.com/tiagodc/TreeLS
  # kann auch nur Stamm + Reise ausschneiden um Stamm

# https://github.com/apburt/treeseg/
  # nicht in R & benötigt Ubuntu

###

# delete points above 2m height
cloud_under <- filter_poi(cloud_norm, Z <= 2)
#plot(cloud_under, axis = T)

# save changed data
writeLAS(cloud_norm, paste0(substr(path_points, 1, nchar(path_points)-4), "_normalized.las"))
writeLAS(cloud_under, paste0(substr(path_points, 1, nchar(path_points)-4), "_understory.las"))

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

# raster_geometry <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE) {
#   # takes in point cloud (normalized or unnormalized)
#   # saves several geometry rasters in output folder
#   # https://gis.stackexchange.com/questions/396186/calculating-grid-metrics-in-a-loop-using-column-name-within-variable/396191#396191
#   check_create_dir(output_dir)
#   print("... creating geometry rasters")
#   # get attributes
#   point_cloud <- add_geometry(point_cloud)
#   # loop through all geometries
#   for (name in c("curvature", "linearity", "planarity", "sphericity", "anisotropy")) {
#     # loop through all statistics to use
#     for (type in c("mean", "sd", "max")) {
#       # print progress
#       print(paste0("... creating ", name, " (", type,") raster"))
#       # make rasters
#       expr = paste0("~", type, "(", name, ")")
#       expr = eval(parse(text = expr))
#       geom_raster <- grid_metrics(point_cloud, expr, res = resolution)
#       # rescale rasters
#       if (rescale) {
#         geom_raster <- rescale_raster(geom_raster)
#       }
#       # save rasters
#       writeRaster(geom_raster, paste0(output_dir, "/", name, "_", type, "_",
#                                       output_name, "_", resolution*100, "cm.tif"), overwrite = TRUE)
#     }
#   }
# }

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

# raster_reflectance <- function(point_cloud, resolution, output_dir, output_name, rescale=TRUE){
#   # takes in point cloud (normalized or unnormalized)
#   # saves intensity raster in output folder
#   check_create_dir(output_dir)
#   print("... creating reflectance rasters")
#   # create intensity raster
#   reflect <- grid_metrics(point_cloud, ~metric_reflectance(Reflectance), res = resolution)
#   # rescale raster
#   if (rescale) {
#     reflect <- rescale_raster(reflect)
#   }
#   # save intensity raster
#   writeRaster(reflect, paste0(output_dir, "/reflectance_", output_name, "_",
#                               resolution*100, "cm.tif"), overwrite = TRUE)
# }

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

# # testing raster_nDSM
# raster_nDSM(cloud_norm, 0.01, paste0(path_rasters, "/nDSM"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/nDSM/nDSM_test_1cm.tif"))
# 
# # testing raster_ortho
# raster_ortho(cloud_norm, 0.01, paste0(path_rasters, "/ortho"), "test")
# plotRGB(stack("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/ortho/ortho_test_1cm.tif")*255)
# 
# # testing raster_geometry
# raster_geometry(cloud_norm, 0.01, paste0(path_rasters, "/geometry"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/geometry/planarity_mean_test_1cm.tif"))
# 
# # testing raster_reflectance
# raster_reflectance(cloud_norm, 0.01, paste0(path_rasters, "/reflectance"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/reflectance/reflectance_test_1cm.tif"))
# 
# # testing raster_point_density
# raster_point_density(cloud_norm, 0.01, paste0(path_rasters, "/point_density"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/point_density/point_density_test_1cm.tif"))

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

# alle Bänder in einem Bild speichern?
# ...dann wäre aber viel redundant
# ...aber dann könnte ich den in Keras implementierten image augmenter nutzen
#   (kann ich da einstellen, welche Klasse wie oft erweitert wird??)
# beim Bilder Laden beim Modell die Bilder stacken?
# stack(stack1, stack2)

################################################################################

# calculating all necessary rasters
raster_create_all(cloud_norm, 0.01, path_rasters, "Breisach")
raster_nDSM(cloud_raw, resolution, paste0(raster_dir, "/nDSM_unscaled"), output_name, rescale=FALSE)  # has to be uncut cloud

################################################################################
# COLINEARITY GEOMETRY CHECKS
################################################################################

# path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters_copy"  # multiple areas
# path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"  # single area

# for this, all rasters must be computed and normalized
# this part must be done manually!

# get all relevant raster paths
raster_list <- list.files(path_rasters, pattern=".tif", recursive=TRUE)
raster_list <- raster_list[!grepl("nDSM", raster_list)]
raster_list <- raster_list[!grepl("ortho", raster_list)]

# load all those rasters
rasters <- list()
for (i in raster_list) {rasters[i] = raster(paste0(path_rasters, "/", i))}

# get unique area names
area_names <- c()
for (name in raster_list) {
  name <- strsplit(name, "[.]")[[1]][1]
  name <- strsplit(name, "_")
  name <- name[[1]][length(name[[1]])-1]
  area_names <- c(area_names, name)
}
unique_areas <- unique(area_names)

# stack rasters (of same area respectively)
raster_stacks <- c()
for (area in unique_areas) {
  single_stack <- stack(rasters[grepl(paste0("_", area, "_"), raster_list)])
  raster_stacks <- c(raster_stacks, single_stack)
  rm(single_stack)
}

# change the band names so they match & are short
band_names <- names(raster_stacks[[1]])  # assumes all have same rasters available
new_band_names <- c()
for (name in band_names) {
  name <- strsplit(name, "[.]")[[1]][2]
  name <- strsplit(name, "_")
  name <- name[[1]][1:2]
  name <- paste(name[1], name[2], sep="_")
  new_band_names <- c(new_band_names, name)
}
for (i in 1:length(raster_stacks)) {
  names(raster_stacks[[i]]) <- new_band_names
}

# prepare empty data frame for samples
raster_df <- matrix(ncol = length(new_band_names), nrow=0)
colnames(raster_df) <- new_band_names
raster_df <- data.frame(raster_df)

# draw samples from each plot & put them all into one data frame
for (i in 1:length(raster_stacks)) {
  raster_samples <- sampleRandom(raster_stacks[[i]], size=5000, cells=FALSE, sp=TRUE)
  raster_samples_df <- as.data.frame(raster_samples@data, xy=FALSE)
  raster_df <- rbind(raster_df, raster_samples_df)
  rm(raster_samples)
  rm(raster_samples_df)
}

################################################################################

## PCA 1

# PCA biplot:
# - gleiche Richtung = positiv korreliert
# - entgegengesetzt = negativ korreliert
# - 90 Grad = nicht korreliert

pca_all <- prcomp(raster_df)
summary(pca_all)  # with 5 PCs we are above the 90% explained variance
names(pca_all$sdev) <- as.character(1:20)
screeplot(pca_all, las=1, main="", cex.lab=1.5, xlab="Hauptkomponenten")  # to see which PCs explain how much variance
round(pca_all$rotation, 2)  # here we can see how important the variables is for each PC
biplot(pca_all, cex=c(0.5,1))#, xlim=c(-0.05, 0.05), ylim=c(-0.05, 0.05))

################################################################################

## CLUSTER

# threshold: spearmans rho = 0.7
# beim rausnehmen drauf achten, was mehr zu den ersten PCA beiträgt?

# look at the situation
clust_all_1 <- varclus(as.matrix(raster_df), similarity = c("spearman"))
plot(clust_all_1)
abline(h=0.7^2, lty=2, col="deeppink3")

# handle reflectance
clust_all_2 <- varclus(as.matrix(raster_df[-c(17,19)]), similarity = c("spearman"))
plot(clust_all_2)
abline(h=0.7^2, lty=2, col="deeppink3")

# handle planarity & linearity
clust_all_3 <- varclus(as.matrix(raster_df[-c(17,19, 11,8,9)]), similarity = c("spearman"))
plot(clust_all_3)
abline(h=0.7^2, lty=2, col="deeppink3")

# handle curvature & shericity & ansiotrophy
clust_all_4 <- varclus(as.matrix(raster_df[-c(17,19, 11,8,9, 4,16,3,14,15,5,2)]), similarity = c("spearman"))
plot(clust_all_4)
abline(h=0.7^2, lty=2, col="deeppink3")

# checking if everything under threshold
clust_all_4[["sim"]]
clust_all_4[["sim"]]>(0.7^2)
cor(raster_df[-c(17,19, 11,8,9, 4,16,3,14,15,5,2)], method="spearman")

################################################################################

## PCA 2

pca_remains <- prcomp(raster_df[-c(17,19, 11,8,9, 4,16,3,14,15,5,2)])
summary(pca_remains)  # with 5 PCs we are above the 90% explained variance
names(pca_remains$sdev) <- as.character(1:20)
screeplot(pca_remains, las=1, main="", cex.lab=1.5, xlab="Hauptkomponenten")  # to see which PCs explain how much variance
round(pca_remains$rotation, 2)  # here we can see how important the variables is for each PC
biplot(pca_remains, cex=c(0.5,1))#, xlim=c(-0.05, 0.05), ylim=c(-0.05, 0.05))

################################################################################

# save the final rasters in combinations the CNN should be fed with?

################################################################################