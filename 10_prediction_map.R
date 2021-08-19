################################################################################
################################################################################
# PREDICTION MAP
################################################################################
################################################################################

# load packages
library(raster)
library(sf)
library(sp)
library(keras)
library(dplyr)
library(caret)

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/07_setup_cnn.R")

# set paths
basedir <- "H:/Daten/Studium/2_Master/4_Semester"
path_rasters <- paste0(basedir, "/4_Daten/rasters_2cm") # input
path_hyper   <- paste0(basedir, "/4_Daten/02_testing/models_2cm/tls_rgb/best_hyperparameters.csv") # input
path_area    <- paste0(basedir, "/4_Daten/sites/convex/area_polygons.shp") # input
path_values  <- paste0(basedir, "/4_Daten/clips_2cm_unscaled/raster_values_labels_unscaled.csv") # input
path_clips   <- paste0(basedir, "/4_Daten/model_input_2cm_standardized/tls_rgb") # input
path_plots   <- paste0(basedir, "/4_Daten/vegetation/Export_ODK_clean_checked_filtered_no_overlap.kml") # input
path_maps    <- paste0(basedir, "/5_Analyse/Karten") # output
path_out     <- paste0(basedir, "/4_Daten/prediction_map") # output
path_labels  <- paste0(path_clips, "/label_lookup.csv") # input
check_create_dir(path_out)

# set parameters
n_pixels <- 25
n_bands <- 7
n_classes <- 5
area_ID <- 6 # because it has all classes (1,6,7,8)
crs_raster_las <- "+proj=utm +zone=32 +ellps=WGS84 +units=m +vunits=m +no_defs"
use_these <- c(
  "ortho",
  # "anisotropy_max", "curvature_max", "linearity_max", "linearity_sd",
  "nDSM",
  # "planarity_mean", "planarity_sd",
  "point_density", "reflectance_mean", "reflectance_sd"
)

# use custom color palette, bright
own_colors_named <- list(
  red = "#ff6f69",
  blue = "#00aedb",
  yellow = "#ffcc5c",
  green = "#88d8b0",
  turquoise = "#4abdac",
  pink = "#ff8b94",
  bright_green = "#cbe885"
)

# color palette for vegetation classes
color_scale_classes <- c(
  own_colors_named$blue, own_colors_named$red, own_colors_named$yellow,
  own_colors_named$bright_green, own_colors_named$green
)

################################################################################
# HELPER FUNCTIONS
################################################################################

rescale_values <- function(values_path) {
  # creates lookup table for normalization / standardization
  values <- read.csv(values_path)
  lookup <- list()
  types <- unique(values$type)
  for (type in types) {
    lookup[[paste0(type, "_min")]] <- min(values$values[values$type == type], na.rm = TRUE)
    lookup[[paste0(type, "_max")]] <- max(values$values[values$type == type], na.rm = TRUE)
    lookup[[paste0(type, "_mean")]] <- mean(values$values[values$type == type], na.rm = TRUE)
    lookup[[paste0(type, "_sd")]] <- sd(values$values[values$type == type], na.rm = TRUE)
  }
  return(lookup)
}

################################################################################
# MAKE NEW MODEL
################################################################################

# saving all images from path list as rds
file_list_to_rds <- function(file_list, output_dir, name, pixels = 25, bands = 7, seed = 123) {
  # set seed
  set.seed(seed)
  # remove grass & rock
  file_list <- file_list[!grepl("grass", basename(file_list)) & !grepl("rock", basename(file_list))]
  # get all labels
  img_labels <- sapply(strsplit(basename(file_list), "[.]"), "[[", 1)
  img_labels <- gsub("[[:digit:]]", "", img_labels)
  img_labels <- sapply(img_labels, function(chars) substr(chars, 1, nchar(chars) - 1))
  img_labels_txt <- as.factor(img_labels)
  img_labels_num <- as.numeric(img_labels_txt)
  # save which label corresponds to which number
  label_lookup <- data.frame("old" = unique(img_labels_txt), "new" = unique(img_labels_num))
  write.csv(label_lookup, paste0(output_dir, "/label_lookup.csv"), row.names = FALSE)
  # load all files
  img_array <- array(dim = c(length(file_list), pixels, pixels, bands))
  label_vector <- c()
  for (i in 1:length(file_list)) {
    path <- file_list[i]
    label_vector <- c(label_vector, img_labels_num[i])
    raster <- raster::values(stack(path))
    img_array[i, , , ] <- raster
  }
  # save image array as rds
  data_list <- list(img = img_array, label = label_vector)
  saveRDS(data_list, file = paste0(output_dir, "/", name, ".rds"))
  # remove seed
  set.seed(NULL)
  return(paste0(output_dir, "/", name, ".rds"))
}

# helper function to calculate mode
mode <- function(vector) {
  uniques <- unique(vector)
  uniques[which.max(tabulate(match(vector, uniques)))]
}

################################################################################

# load plots & areas
areas <- st_read(path_area)
area <- areas[areas$ID == area_ID,]
points <- st_transform(st_read(path_plots), crs(area))

# get IDs of points within / without area
contained <- st_contains(area, points)
within <- points[contained[[1]],]
without <- points[-contained[[1]],]
within_ID <- as.numeric(lapply(strsplit(within$Description, "veg_ID: "), "[[", 2))
without_ID <- as.numeric(lapply(strsplit(without$Description, "veg_ID: "), "[[", 2))

# load vegetation plots, split into training / testing
all_paths <- list.files(path_clips, pattern = "[.]tif", full.names = TRUE)
all_paths_ids <- as.numeric(gsub("\\D", "", basename(all_paths)))
within_paths <- all_paths[all_paths_ids %in% within_ID]
without_paths <- all_paths[all_paths_ids %in% without_ID]

# get share of labels
test_label_share <- c()
for (label in c("blueberry", "dead_wood", "forest_floor", "moss", "spruce")) {
  total <- sum(grepl(pattern = label, basename(within_paths)))
  fraction <- total/length(within_paths)
  test_label_share <- rbind(test_label_share, data.frame(
    "label" = label,
    "total" = total,
    "fraction" = round(fraction*100, 2)
  ))
}
write.csv(test_label_share, paste0(path_out, "/label_share_test_images.csv"), row.names = FALSE)

# save as rds
file_list_to_rds(within_paths, path_out, "input_testing")
file_list_to_rds(without_paths, path_out, "input_training")

################################################################################

# load data
rdata_paths <- paste0(path_out, "/", c("input_training", "input_testing"), ".rds")
balanced <- create_dataset(rdata_paths, 2, n_pixels, n_bands, na_replacement = 0)

# load best hyperparameters from best dataset (tls_rgb)
best_runs <- read.csv(path_hyper)
best_runs <- best_runs[,grepl("flag_", names(best_runs))]
hyper <- apply(best_runs, 2, mode)

# build model
model <- get_lenet5(
  width_length = n_pixels,
  n_bands = n_bands,
  n_classes = n_classes,
  dropout = hyper["flag_dropout"],
  l2_regularizer = hyper["flag_l2_regularizer"]
)

# compile model
model %>% compile(
  optimizer = optimizer_adam(lr = hyper["flag_learning_rate"]),
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

# train model
if (hyper["flag_batch_size"] == 16) {
  history <- model %>% fit_generator(
    balanced$data_train_vali_16,
    steps_per_epoch = balanced$steps_train_vali_16,
    epochs = hyper["flag_epochs"],
  )
}
if (hyper["flag_batch_size"] == 32) {
  history <- model %>% fit_generator(
    balanced$data_train_vali_32,
    steps_per_epoch = balanced$steps_train_vali_32,
    epochs = hyper["flag_epochs"],
  )
}

# save model
model %>% save_model_hdf5(paste0(path_out, "/model.h5"))

################################################################################
# GET MODEL PERFORMANCE
################################################################################

# load model
model <- load_model_hdf5(paste0(path_out, "/model.h5"))

# load data
test_data <- readRDS(paste0(path_out, "/", "input_testing.rds"))
test_images <- test_data$img
test_images[is.na(test_images)] <- 0
test_labels <- test_data$label
test_labels <- to_categorical(test_labels)[,2:6]

# get accuracy
results <- model %>% evaluate(test_images, test_labels)
overall_accuracy <- results["accuracy"]

# get mean f1
preds <- model %>% predict(test_images)
preds_int <- apply(preds, 1, which.max)
test_labels_int <- apply(test_labels, 1, which.max)
preds_fac <- factor(preds_int, levels = c(1,2,3,4,5))
test_labels_fac <- factor(test_labels_int, levels = c(1,2,3,4,5))
conf <- confusionMatrix(preds_fac, test_labels_fac, mode = "everything")
f1_vals <- conf$byClass[,"F1"]
f1_vals[is.na(f1_vals)] <- 0
mean_f1 <- mean(f1_vals)

# save in file - overal acc & mean f1
write.csv(data.frame(
  "accuracy" = round(overall_accuracy*100, 2),
  "mean_f1" = round(mean_f1, 4)),
  paste0(path_out, "/model_acc_f1.csv"), row.names = FALSE)

# save in file - per class f1
write.csv(data.frame(
  "class" = read.csv(path_labels)$old,
  "class_f1" = round(f1_vals, 4)),
  paste0(path_out, "/model_class_f1.csv"), row.names = FALSE)

################################################################################
# TILE RASTER
################################################################################

# get all rasters of the area
raster_paths <- list.files(path_rasters, pattern = paste0("area_", area_ID), recursive = T, full.names = T)
raster_paths <- raster_paths[grepl(pattern = "[.]tif", raster_paths)]
raster_paths <- raster_paths[!grepl(pattern = "nDSM_filtering", raster_paths)]

# load areas
areas <- st_read(path_area)
st_crs(areas) <- CRS("+init=EPSG:25832")
areas <- st_transform(areas, crs_raster_las)
area <- areas[areas$ID == area_ID, ]

# get lookup table for normalization
lookup_table <- rescale_values(path_values)

# get raster files which are within the selection
types <- c()
for (i in 1:length(raster_paths)) {
  types[i] <- strsplit(basename(raster_paths[i]), "_area_")[[1]][1]
}
raster_paths <- raster_paths[types %in% use_these]
types <- types[types %in% use_these]

# empty object for rasters
raster_list <- list()

# for every raster
for (i in 1:length(raster_paths)) {
  # load raster & type
  type <- types[i]
  path <- raster_paths[i]
  raster_list[[i]] <- stack(path)
  
  # clip to area
  raster_list[[i]] <- mask(raster_list[[i]], area)
  
  # standardize using mean / sd from all data
  bands <- if (type != "ortho") c(type) else c("R", "G", "B")
  for (j in 1:length(bands)) {
    mean_val <- lookup_table[[paste0(bands[j], "_mean")]]
    sd_val <- lookup_table[[paste0(bands[j], "_sd")]]
    raster_list[[i]][[j]] <- scale(raster_list[[i]][[j]], center = mean_val, scale = sd_val)
  }
}

# stack rasters
raster_stack <- stack(raster_list)
rm(raster_list)

# make tiles
agg <- aggregate(raster_stack[[1]], n_pixels)
polys <- as(agg, "SpatialPolygons")
tiles <- lapply(seq_along(polys), function(i) crop(raster_stack, polys[i]))

# save tiles
saveRDS(tiles, file = paste0(path_out, "/tiles.rds"))

################################################################################
# MAKE PREDICTIONS
################################################################################

# load tiles as arrays
tiles <- readRDS(paste0(path_out, "/tiles.rds"))
img_array <- array(dim = c(length(tiles), n_pixels, n_pixels, n_bands))
for (i in 1:length(tiles)) {
  img_array[i, , , ] <- values(tiles[[i]])
}

# replace missing values with 0
img_array[is.na(img_array)] <- 0

# load model
model <- load_model_hdf5(paste0(path_out, "/model.h5"))

# predict
preds <- model %>% predict(img_array)

# get prediction & confidence
preds_label <- apply(preds, 1, which.max)
preds_chance <- apply(preds, 1, max)

# replace numeric label with text label
label_lookup <- read.csv(path_labels)
for (i in 1:nrow(label_lookup)) {
  preds_label[preds_label == label_lookup$new[i]] <- label_lookup$old[[i]]
}

# get share of labels
preds_label_share <- c()
for (label in c("blueberry", "dead_wood", "forest_floor", "moss", "spruce")) {
  total <- sum(preds_label == label)
  fraction <- total/length(preds_label)
  preds_label_share <- rbind(preds_label_share, data.frame(
    "label" = label,
    "total" = total,
    "fraction" = round(fraction*100, 2)
  ))
}
write.csv(preds_label_share, paste0(path_out, "/label_share_predictions.csv"), row.names = FALSE)

# get mean confidence in prediction per label
# (where label was chosen as final label)
preds_combo <- data.frame(
  "prediction" = preds_label,
  "chance" = preds_chance)
preds_combo_summary <- preds_combo %>%
  group_by(prediction) %>%
  summarise(mean_chance = mean(chance))
preds_combo_summary <- rbind(preds_combo_summary, data.frame(
  "prediction" = "total",
  "mean_chance" = mean(preds_chance)))
write.csv(preds_combo_summary, paste0(path_out, "/label_mean_confidence.csv"), row.names = FALSE)

# save as shapefile
polys <- lapply(tiles, function(x) as(extent(x), "SpatialPolygons"))
polys <- do.call(bind, polys)
polys <- SpatialPolygonsDataFrame(Sr = polys, data = data.frame(
  id = 1:length(polys),
  prediction = preds_label,
  chance = round(preds_chance, 2)
))
crs(polys) <- CRS(crs_raster_las)
shapefile(polys, paste0(path_out, "/tiles.shp"), overwrite = TRUE)

################################################################################
# MAPS
################################################################################

# load packages
library(ggmap)
library(rgdal)
library(ggplot2)
library(ggpubr)
library(broom)
library(dplyr)
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)

# load & prepare data
shapefile <- readOGR(paste0(path_out, "/tiles.shp"))
shapefile[is.na(shapefile$prediction), ]
df <- tidy(shapefile, region = "id")
shapefile$id <- rownames(shapefile@data)
df <- left_join(df, shapefile@data, by = "id")
df <- df[!is.na(df$prediction) & !is.na(df$chance), ]

# make legend (predictions)
map_legend_label <- get_legend(
  ggplot() +
    geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = prediction), color = NA) +
    scale_fill_manual(
      values = color_scale_classes,
      name = "Predicted Class",
      labels = c("Blueberry", "Deadwood", "Forest Floor", "Moss", "Spruce")
    ) +
    theme(
      text = element_text(family = "Calibri", size = 14),
      legend.text = element_text(family = "Calibri", size = 14),
      legend.box.spacing = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.title = element_text(family = "Calibri", size = 16),
      legend.position = "bottom"
    )
)

# make map (predictions)
cairo_pdf(file = paste0(path_maps, "/prediction_map_label.pdf"), family = "Calibri", width = 8.27, height = 5.83)
map_label <- ggplot() +
  geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = prediction), color = NA) +
  scale_fill_manual(
    values = color_scale_classes,
    name = "Predicted Class",
    labels = c("Blueberry", "Deadwood", "Forest Floor", "Moss", "Spruce")
  ) +
  theme_light() +
  coord_equal() +
  ylim(c(5379630, 5379665)) +
  xlim(c(447025, 447075)) +
  ylab("Northing\n") +
  xlab("\nEasting") +
  theme(text = element_text(family = "Calibri", size = 14),
        legend.position = "none")
ggarrange(map_label, legend.grob = map_legend_label, legend = "bottom") # without ggarrange, legend not centered
dev.off()

################################################################################

# make legend (confidence)
map_legend_confidence <- get_legend(
  ggplot() +
    geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = chance), color = NA) +
    scale_fill_stepsn(
      limits = c(0, 1), breaks = seq(0, 1, 0.2), show.limits = TRUE,
      colors = c(own_colors_named$red, "grey90", own_colors_named$blue),
      name="Prediction Confidence\n") +
    theme(text = element_text(family = "Calibri", size = 14),
          legend.text = element_text(family = "Calibri", size = 14),
          legend.box.spacing = unit(0.5, "cm"),
          legend.key.width = unit(1.75, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.title = element_text(family = "Calibri", size = 16),
          legend.position = "bottom")
)

# make map (confidence)
cairo_pdf(file = paste0(path_maps, "/prediction_map_confidence.pdf"), family = "Calibri", width = 8.27, height = 5.83)
map_confidence <- ggplot() +
  geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = chance), color = NA) +
  scale_fill_stepsn(
    limits = c(0, 1), breaks = seq(0, 1, 0.2), show.limits = TRUE,
    colors = c(own_colors_named$red, "grey90", own_colors_named$blue),
    name="Prediction Confidence\n") +
  theme_light() +
  coord_equal() +
  ylim(c(5379630, 5379665)) +
  xlim(c(447025, 447075)) +
  ylab("Northing\n") +
  xlab("\nEasting") +
  theme(text = element_text(family = "Calibri", size = 14),
        legend.position = "none")
ggarrange(map_confidence, legend.grob = map_legend_confidence, legend = "bottom") # without ggarrange, legend not centered
dev.off()

# CRS: EPSG 25832 (ETRS89 / UTM zone 32N)

################################################################################
