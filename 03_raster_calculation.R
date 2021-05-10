################################################################################
################################################################################
# RASTER CALCULATION
################################################################################
################################################################################

# load packages
library(lidR)
library(TreeLS)

# set paths
#path_points   <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/points/thinning_mean_cirle10.las"
path_points   <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/_Testwolke/OT8cm_small.las"
path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"
path_vegetation <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/vegetation/Export_ODK_clean_2D.kml"



# load data
cloud_raw <- readTLSLAS(path_points) # TODO: passt eine Wolke komplett in RAM?
cloud_norm <- cloud_raw # TODO kommt weg später maaal
rm(cloud_raw) # TODO kommt weg später maaal

# https://cran.r-project.org/web/packages/lidR/lidR.pdf

# TODO: RGB hat mehrere Bänder pro Raster --> doof?

################################################################################
# delete invalid / unrealistic points
################################################################################

# what are common quality checks for TLS data?

# save changed data

################################################################################
# filter understory points
################################################################################

# set parameters
resolution_dtm = 0.1

# normalize height
cloud_raw <- classify_ground(cloud_raw, csf(class_threshold = 0.1, cloth_resolution = 1)) # is later used for creating nDSM

dtm <- grid_terrain(filter_ground(cloud_raw), tin(), res = resolution_dtm)
cloud_norm <- normalize_height(cloud_raw, dtm, na.rm = T)
#plot(cloud_norm, axis = T)

# Punkte mit negativem z-Wert raus kicken
cloud_norm <- filter_poi(cloud_norm, Z <= 2)

# delete trees from data when bigger than 2 meters
?find_trees
# https://github.com/Jean-Romain/lidR/wiki/Segment-individual-trees-and-compute-metrics#segment-the-trees
# https://rdrr.io/cran/lidR/man/add_attribute.html
# https://gis.stackexchange.com/questions/376772/removing-tree-stems-from-a-point-cloud-in-r
# https://jean-romain.github.io/lidRbook/itd-its.html


# delete points above 2m height
cloud_under <- filter_poi(cloud_norm, Z <= 2)
#plot(cloud_under, axis = T)

# (mit Nummer des Returns overstory rauswerfen?)

# save changed data
writeLAS(cloud_norm, paste0(substr(path_points, 1, nchar(path_points)-4),
                            "_normalized.las"))
writeLAS(cloud_under, paste0(substr(path_points, 1, nchar(path_points)-4),
                             "_understory.las"))

################################################################################
# helper functions
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

# even faster eigenvalue calculation
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
  eigen <- eigen_decomposition(las, 20, 6) # 20 neighbours, 6 cores
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
  # TODO: welche Statistik will ich haben? mean?
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
# create single rasters
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
  chm <- grid_canopy(point_cloud_norm, res = resolution, p2r())
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
  # https://gis.stackexchange.com/questions/396186/calculating-grid-metrics-in-a-loop-using-column-name-within-variable/396191#396191
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
    # save rasters
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
  writeRaster(point_density, paste0(output_dir, "/density_", output_name, "_",
                                    resolution*100, "cm.tif"), overwrite = TRUE)
}

################################################################################

# # testing raster_nDSM
# raster_nDSM(cloud_norm, 0.01, paste0(path_rasters, "/nDSM"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/nDSM/nDSM_test_1cm.tif"))
# 
# # testing raster_ortho
# raster_ortho(cloud_norm, 0.01, paste0(path_rasters, "/ortho"), "test")
# plotRGB(stack("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/ortho/ortho_test_10cm.tif")*255)
# 
# # testing raster_geometry
# raster_geometry(cloud_norm, 0.01, paste0(path_rasters, "/geometry"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/geometry/planarity_mean_test_10cm.tif"))
# 
# # testing raster_reflectance
# raster_reflectance(cloud_norm, 0.01, paste0(path_rasters, "/reflectance"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/reflectance/reflectance_test_10cm.tif"))
# 
# # testing raster_point_density
# raster_point_density(cloud_norm, 0.01, paste0(path_rasters, "/density"), "test")
# plot(raster("H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/density/density_test_10cm.tif"))

################################################################################
# create all rasters
################################################################################

raster_create_all <- function(point_cloud, resolution, raster_dir, output_name) {
  raster_nDSM(point_cloud, resolution, paste0(raster_dir, "/nDSM"), output_name)
  raster_ortho(point_cloud, resolution, paste0(raster_dir, "/ortho"), output_name)
  raster_geometry(point_cloud, resolution, paste0(raster_dir, "/geometry"), output_name)
  raster_reflectance(point_cloud, resolution, paste0(raster_dir, "/reflectance"), output_name)
  raster_point_density(point_cloud, resolution, paste0(raster_dir, "/density"), output_name)
  print("done!")
}

# alle Bänder in einem Bild speichern?
# ...dann wäre aber viel redundant
# ...aber dann könnte ich den in Keras implementierten image augmenter nutzen
#   (kann ich da einstellen, welche Klasse wie oft erweitert wird??)
# beim Bilder Laden beim Modell die Bilder stacken?
# stack(stack1, stack2)

# # testing raster_create_all
# raster_create_all(cloud_norm, 0.01, path_rasters, "Breisach")

################################################################################
# create raster tiles
################################################################################

path_raster <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/dummy/raster.tif"
path_points <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/dummy/points.kml"

raster_clip_all <- function(raster_dir, plot_path, tile_size, output_dir) {
  # clips all rasters to small areas around the vegetation plots
  check_create_dir(output_dir)
  print("... clipping all rasters")
  # read as kml + transform CRS
  plots <- sf::st_read(plot_path)
  plots <- sf::st_transform(plots, 25832)
  # TODO: add ID?
  
  # get all rasters within raster_dir
  raster_list <- list.files(raster_dir, pattern=".tif", recursive=TRUE)
  # calculate edge length
  edge <- ((tile_size*100)%/%2)/100
  # loop through rasters
  for (raster_path in raster_list) {
    # load raster
    raster <- raster(paste0(raster_dir,"/",raster_path))
    crs(raster) <- CRS("+init=EPSG:25832")
    # get points within raster extent
    # https://gis.stackexchange.com/questions/230900/crop-simple-features-object-in-r
    subset <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(raster), "SpatialPolygons")), st_crs(plots)))
    # loop through plots
    for (idx in 1:nrow(subset)) {
     plot <- subset[idx,]
      # set clipping extent beginning from plot center
      center_x <- round(sf::st_coordinates(plot)[,1], 2) # round on cm
      center_y <- round(sf::st_coordinates(plot)[,2], 2) # round on cm
      rectangle <- extent(c(xmin=center_x-edge, xmax=center_x+edge,
                            ymin=center_y-edge, ymax=center_y+edge))
      # clip raster
      clip <- crop(raster, rectangle)
      # save clip
      writeRaster(clip, paste0(), overwrite=TRUE) # TODO: Bennenung???
    }
  }
}

# TODO: exclude all tiles where the height of the 99% percentile is above 2m height

library(lidR)
plot_path <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/vegetation/Export_ODK_clean_2D.kml"
raster_dir <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters"
resolution <- 0.01
tile_size <- 0.5
output_dir <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters/clips"

# Datei anlegen, die die Koordinaten der unteren linken Ecke aller Kacheln speichert?

################################################################################
# check for correlations
################################################################################

# hierfür sind fertige, normalisierte Raster notwendig, die in dataframe umgewandlet werden

raster_list <- list.files(path_rasters, pattern=".tif", recursive=TRUE)
rasters <- list()
for (i in raster_list) {rasters[i] = raster(paste0(path_rasters, "/", i))}
# Raster-Werte: pro Raster-Typ eine Spalte
# for (i in raster_list) {rasters[i] = as.vector(raster(paste0(path_rasters, "/", i)))}
# Raster: Raster mit rbind drunter hauen, wenn mehrere da
# 

# PROBLEM:mehrere Bänder pro Bild manchmal! (rgb, reflectance)

# PCA biplot:
# - gleiche Richtung = positiv korreliert
# - entgegengesetzt = negativ korreliert
# - 90 Grad = nicht korreliert

# Code von Env Stat Kurs
# TODO: anpassen an neue Daten

    # PCA
    # making pretty graphs: https://www.datacamp.com/community/tutorials/pca-analysis-r
    
    PCA_bird <- prcomp(bird[,-c(1:3)])
    str(PCA_bird)
    summary(PCA_bird) # with 5 PCs we are above the 90% explained variance
    
    names(PCA_bird$sdev) <- as.character(1:12)
    par(mar=c(5,5,1,1))
    screeplot(PCA_bird, las=1, main="", cex.lab=1.5, xlab="Hauptkomponenten")  # to see which PCs explain how much variance
    
    round(PCA_bird$rotation, 2) # here we can see how important the variables is for each PC
    round(PCA_bird$x,2)  # new coordinates of the data points
    biplot(PCA_bird, xlim=c(-0.045, 0.045), ylim=c(-0.045, 0.045))  # idk what is going on here
    
    
    # CLUSTER
    
    clust_1 <- varclus(as.matrix(bird[,-c(1:3)]), similarity = c("spearman"))  # make a cluster
    plot(clust_1)
    abline(h=0.7^2, lty=2, col="tomato")
    # I have 3 variables, which are okay, I have to exclude a lot
    
    clust_2 <- varclus(as.matrix(bird[,-c(1:3, 5:7, 11, 14:15)]), similarity = c("spearman"))
    plot(clust_2)
    abline(h=(0.7^2), lty=3, col="tomato")
    # decision between predictors
    
    clust_3 <- varclus(as.matrix(bird[,-c(1:3, 5:7, 9, 10:11, 14)]), similarity = c("spearman"))
    plot(clust_3)
    abline(h=(0.7^2), lty=3, col="tomato")
    # kick out SAVANNA, because SHRUBS are more important for my bird
    # kick out PRE_YEAR because we need to have PRE_SUMMER in there
    # kick out TMIN_JUL, TMIN_JAN, T_SUMMER, PRE_WINTER for T_WINTER because all climate zones my bird likes have mild winters
    # kick out SHRUBS for TDIFF because the climate variables mostlikely kinda explain SHRUBS also
    
    # checking if everything under threshold
    clust_3[["sim"]]
    clust_3[["sim"]]>(0.7^2)
    cor(bird[,-c(1:3, 5:7, 9, 10:11, 14)], method="spearman")

    
