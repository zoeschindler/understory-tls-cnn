################################################################################
################################################################################
# PREPARE VEGETATION
################################################################################
################################################################################

# load packages
library(jsonlite)
library(dplyr)
library(sf)
library(sp)
library(jpeg)
library(raster)

# set path
area_dir <- "D:/Masterarbeit_Zoe/4_Daten/sites/convex"  # output
home_dir <- "D:/Masterarbeit_Zoe/4_Daten/vegetation"  # input & output
setwd(home_dir)

# clean the data?
clean = TRUE
if (clean) {
  file_name <- "Export_ODK_clean"
} else {
  file_name <- "Export_ODK_raw"
}

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
# READ JSON
################################################################################

json_to_df <- function (path) {
  full_data <- fromJSON(path)
  final <- data.frame()
  for(d in 1:nrow(full_data)){
    data <- full_data[d,]
    # extract the image subtable
    img <- data$collect_ground_vegetation[[1]]$ground_img
    # exclude the image subtable
    data$collect_ground_vegetation[[1]] <- data$collect_ground_vegetation[[1]][names(data$collect_ground_vegetation[[1]]) != 'ground_img']
    # get the start and end of the questionaire (metadata)
    meta <- data[names(data) != 'collect_ground_vegetation']
    # get the groundsamples
    ground_veg <- data$collect_ground_vegetation[[1]]
    # put it all together
    final <- rbind(final, cbind(cbind(ground_veg, img), meta))
  }
  # replace "other" with description
  final$veg_type[final$veg_type == 'other'] <- final$veg_other_text[final$veg_type == 'other']
  # delete unnecessary column
  final$veg_other_text <- NULL
  # set plot_ID
  final$plot_ID <- final %>%
    group_indices(`plot_loc:Longitude`, `plot_loc:Latitude`)
  # return data
  return(final)
}

add_IDs <- function(df) {
  # set plot_ID
  df$plot_ID <- df %>%
    group_indices(`plot_loc:Longitude`, `plot_loc:Latitude`)
  # set vegetation_ID
  df <- df %>%
    arrange(plot_ID) %>%
    mutate(veg_ID = 1:nrow(df))
  return(df)
}

# execute function
part1 <- json_to_df('Ground_Vegetation_Survey_0_1_results.json')
part2 <- json_to_df('Ground_Vegetation_Survey_0_2_results.json')
final <- add_IDs(rbind(part1, part2))
rm(part1, part2)

################################################################################
# CLEANING
################################################################################

if (clean) {
  # clean up data
  summary(as.factor(final$veg_type))
  # deleting points
  final <- final[substr(final$veg_type,1,8) != "Endpunkt",] # end point
  final <- final[final$veg_type != "Schreibtisch",] # test point
  final <- final[final$veg_type != "Stechpalme",] # only once
  final <- final[final$veg_ID != 940,]  # way to far away
  # renaming points
  final$veg_type[final$veg_type=="fir"] <- "spruce" # too dumb for fir / spruce
  final$veg_type[final$veg_type=="sprouce"] <- "spruce" # typo
  final$veg_type[final$veg_type=="Totholz"] <- "dead_wood" # German to English
  final$veg_type[final$veg_type=="Gras"] <- "grass"
  final$veg_type[final$veg_type=="gras"] <- "grass"
  final$veg_type[final$veg_type=="Boden"] <- "forest_floor"
  final$veg_type[final$veg_type=="Boden?"] <- "forest_floor"
  # check remaining types
  summary(as.factor(final$veg_type))
}

################################################################################
# CREATE SITE SHAPES
################################################################################

if (clean) {  # otherwise, really wrong point are included
  check_create_dir(area_dir)
  final_spatial <- st_as_sf(final, coords = c('veg_loc:Longitude', 'veg_loc:Latitude'), crs = 4326)
  final_spatial <- st_transform(final_spatial, 25832)
  poly_list <- list()
  # loop through all unique plot IDs
  plot_IDs <- unique(final_spatial$plot_ID)
  for (plot_ID in plot_IDs) {
    # select all points belonging to the plot
    points <- final_spatial[final_spatial$plot_ID == plot_ID,]
    points_coords <- st_coordinates(points)
    # calculate convex hull
    plot_poly <- chull(st_coordinates(points))
    plot_poly <- points_coords[c(plot_poly, plot_poly[1]),]  # closed polygon
    # convert to spatial polygon
    poly_list[[plot_ID]] <- Polygons(list(Polygon(plot_poly)), ID=plot_ID)
  }
  # combine single polygons
  plot_poly_all <- SpatialPolygons(poly_list)
  crs(plot_poly_all) <- CRS("+init=EPSG:25832")
  # add buffer
  plot_poly_all <- buffer(plot_poly_all, width=2, dissolve=FALSE)
  # save to shapefile
  shapefile(plot_poly_all, paste0(area_dir, "/area_polygons.shp"))
}

################################################################################
# DOWNLOAD IMAGES
################################################################################

# create image folder
if (!dir.exists(paste0(home_dir, "/images"))) {
  print("... creating new folder")
  dir.create(paste0(home_dir, "/images"))
} else {
  print("... using existing folder")
}

# extract & save images
for (i in 1:nrow(final)) {
  if (!is.na(final$url[i])) {
    if (file.exists(paste0("images/", final$filename[i]))) {
      # print("already exists!")
    } else {
      download.file(final$url[i], paste0("images/", final$filename[i]), mode = 'wb')
    }
  }
}

################################################################################
# EXPORT TO KML & CSV
################################################################################

# prepare
final_sp <- final
final_sp$name <- final_sp$veg_type # for the kmls
final_sp$description <- paste0("plot_ID: ", final_sp$plot_ID, ", veg_ID: ", final_sp$veg_ID) # for the kmls

# csv
final_sp_trans_csv <- st_as_sf(final_sp, coords = c('veg_loc:Longitude', 'veg_loc:Latitude', 'veg_loc:Altitude'), crs = 4326)
final_sp_trans_csv <- st_transform(final_sp_trans_csv, 25832)
final_sp_trans_csv <- cbind(as.data.frame(st_coordinates(final_sp_trans_csv)), as.data.frame(final_sp_trans_csv))
final_sp_trans_csv <- final_sp_trans_csv[1:(ncol(final_sp_trans_csv)-3)] # remove geometry & name & description
st_write(final_sp_trans_csv, paste0(file_name, ".csv"), delete_layer = T)

# 2D kml
final_sp_trans_xy <- st_as_sf(final_sp, coords = c('veg_loc:Longitude', 'veg_loc:Latitude'), crs = 4326)
#final_sp_trans_xy <- st_transform(final_sp_trans_xyz, 25832) # wird eh als WGS84 gespeichert
st_write(final_sp_trans_xy, paste0(file_name, ".kml"), delete_layer = T)

################################################################################
# COMPARING IMAGES & LABELS
################################################################################

# forest_floor / moss / grass: what covers more area (>50%)
# forest_floor: bare earth, litter, dead herbs & grasses
# dead_wood: lying deadwood & shrubby deadwood (knee-height), needs to cover at least half of the ground
# spruce & blueberry: when with leaves would be covering more than half of the ground
# higher cover class > lower cover class, e.g. blueberry on moss --> blueberry

################################################################################

# uncomment if this should be done again manually

# # load data
# df <- read.csv("Export_ODK_clean.csv")
# unique(df$veg_type)
# 
# # loop through all points, show image, ask if label is right (I did this 3x)
# plots <- c()
# label <- c()
# for (i in 1:nrow(df)) {
#   print("-----")
#   print(paste0("Plot: ", df$veg_ID[i]))
#   print(paste0("Label: ", df$veg_type[i]))
#   print(paste0("Remarks: ", df$remarks[i]))
#   print(paste0("Filename: ", df$filename[i]))
#   jj <- readJPEG(paste0(home_dir, "/images/", df$filename[i]), native=TRUE)
#   par(mar=c(0,0,0,0))
#   plot(0:1, 0:1, type="n", ann = FALSE, axes = FALSE)
#   rasterImage(jj,0,0,1,1)
#   input <- readline(prompt="Label correct? ")
#   if (input == "exit") {  # to exit, write "exit"
#     break
#   }
#   plots <- c(plots, df$veg_ID[i])
#   label <- c(label, input)
# }
# 
# # save in data frame
# results <- data.frame("plot_ID"=plots, "label"=label)
# 
# # look at & save wrongly labeled data
# wrong_ID <- results$plot_ID[results$label!="y"]
# wrong_df <- df[df$veg_ID %in% wrong_ID,]
# write.csv(wrong_df, "Export_ODK_clean_wrong.csv", row.names=FALSE)
# 
# # combine marked plots from several runs
# df1 <- read.csv("Export_ODK_clean_wrong1.csv")[1:19]
# df2 <- read.csv("Export_ODK_clean_wrong2.csv")[1:19]
# df3 <- read.csv("Export_ODK_clean_wrong3.csv")[1:19]
# df_all <- rbind(df1, df2, df3)
# df_all <- df_all[!duplicated(df_all),]
# write.csv(df_all, "Export_ODK_clean_wrong_all.csv", row.names=FALSE)
# 
# # relabel into mixed / too_small / check_cloud / correct / below spruce
# plots <- c()
# label <- c()
# for (i in 1:nrow(df_all)) {
#   print("-----")
#   print(paste0("Plot: ", df_all$veg_ID[i]))
#   print(paste0("Label: ", df_all$veg_type[i]))
#   print(paste0("Remarks: ", df_all$remarks[i]))
#   print(paste0("Filename: ", df_all$filename[i]))
#   jj <- readJPEG(paste0(home_dir, "/images/", df_all$filename[i]), native=TRUE)
#   par(mar=c(0,0,0,0))
#   plot(0:1, 0:1, type="n", ann = FALSE, axes = FALSE)
#   rasterImage(jj,0,0,1,1)
#   input <- readline(prompt="Label correct? ")
#   if (input == "exit") {  # to exit, write "exit"
#     break
#   }
#   plots <- c(plots, df_all$veg_ID[i])
#   label <- c(label, input)
# }
# 
# # save in data frame & csv
# results <- data.frame("plot_ID"=plots, "label"=label)
# write.csv(results, "Export_ODK_clean_wrong_labels.csv", row.names=FALSE)

################################################################################

if (clean) {
  # read in csv with IDs and new labels of wrongly labeled points
  results <- read.csv("Export_ODK_clean_wrong_labels.csv")
  
  # delete "correct"
  results <- results[results$label!="correct",]
  delete_these <- results$plot_ID[results$label == "too_small" | results$label == "mixed" | results$label == "check_cloud"]
  relabel_these <- results[results$label != "too_small" & results$label != "mixed" & results$label != "check_cloud",]
  
  # change & save csv
  new_csv <- read.csv("Export_ODK_clean.csv")
  new_csv <- new_csv[!(new_csv$veg_ID %in% delete_these),]
  for (i in 1:nrow(relabel_these)) {
    plot_id <- relabel_these$plot_ID[i]
    new_label <- relabel_these$label[i]
    new_csv$veg_type[new_csv$veg_ID == plot_id] <- new_label
  }
  write.csv(new_csv, "Export_ODK_clean_checked.csv")
  
  # change & save kml
  # more complicated because the ID is part of the description column
  new_kml <- st_read("Export_ODK_clean.kml")
  delete_indices <- c()
  for (i in 1:nrow(new_kml)) {
    veg_id <- as.numeric(strsplit(new_kml$Description, " ")[[i]][4])
    if (veg_id %in% delete_these) {
      delete_indices <- c(delete_indices, i)
    }
    if (veg_id %in% relabel_these$plot_ID) {
      new_label <- relabel_these$label[veg_id == relabel_these$plot_ID]
      new_kml$Name[i] <- new_label
    }
  }
  new_kml <- new_kml[-delete_indices,]
  st_write(new_kml, "Export_ODK_clean_checked.kml", delete_layer = T)
}

################################################################################