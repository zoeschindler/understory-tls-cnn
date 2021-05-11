################################################################################
################################################################################
# PREPARE VEGETATION
################################################################################
################################################################################

## parse geojson
require(jsonlite)
require(dplyr)
require(sf)

# set path
home_dir <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/vegetation"
setwd(home_dir)

# cleaning?
clean = TRUE
if (clean) {
  file_name <- "Export_ODK_clean"
} else {
  file_name <- "Export_ODK_raw"
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
  # renaming points
  final$veg_type[final$veg_type=="fir"] <- "spruce" # too dumb for fir / spruce
  final$veg_type[final$veg_type=="sprouce"] <- "spruce" # typo
  final$veg_type[final$veg_type=="Totholz"] <- "dead_wood" # German to English
  final$veg_type[final$veg_type=="Gras"] <- "gras"
  final$veg_type[final$veg_type=="Boden"] <- "forest_floor"
  final$veg_type[final$veg_type=="Boden?"] <- "forest_floor"
  # check remaining types
  summary(as.factor(final$veg_type))
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
  } else {
    # print("no image!")
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
final_sp_trans_csv <- final_sp_trans_csv[1:ncol(final_sp_trans_csv)-2] # remove gemeotry column & name
st_write(final_sp_trans_csv, paste0(file_name, ".csv"), delete_layer = T)

# # 2D shp
# final_sp_trans_xy <- st_as_sf(final_sp, coords = c('veg_loc:Longitude', 'veg_loc:Latitude'), crs = 4326)
# final_sp_trans_xy <- st_transform(final_sp_trans_xy, 25832)
# names(final_sp_trans_xy) <- c("veg_loc_Altitude", "veg_loc_Accuracy", names(final_sp_trans_xy)[3:9],
#                               "plot_loc_Latitude", "plot_loc_Longitude", "plot_loc_Altitude",
#                               "plot_loc_Accuracy", names(final_sp_trans_xy)[14:19])
# st_write(final_sp_trans_xy, paste0(file_name, ".shp"), delete_layer = T)

# 2D kml
final_sp_trans_xy <- st_as_sf(final_sp, coords = c('veg_loc:Longitude', 'veg_loc:Latitude'), crs = 4326)
#final_sp_trans_xy <- st_transform(final_sp_trans_xyz, 25832) # wird eh als WGS84 gespeichert
st_write(final_sp_trans_xy, paste0(file_name, "_2D.kml"), delete_layer = T)

# 3D kml
final_sp_trans_xyz <- st_as_sf(final_sp, coords = c('veg_loc:Longitude', 'veg_loc:Latitude', 'veg_loc:Altitude'), crs = 4326)
#final_sp_trans_xyz <- st_transform(final_sp_trans_xyz, 25832) # wird eh als WGS84 gespeichert
st_write(final_sp_trans_xyz, paste0(file_name, "_3D.kml"), delete_layer = T)
