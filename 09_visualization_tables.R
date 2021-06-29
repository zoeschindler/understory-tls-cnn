################################################################################
################################################################################
# VISUALIZATION & TABLES
################################################################################
################################################################################

# load packages
library(ggplot2)
library(grid)
library(sf)
library(caret)
library(raster)
library(extrafont)
loadfonts(device="pdf", quiet=TRUE)

# set paths
path_vegetation <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/vegetation/Export_ODK_clean_checked_filtered_no_overlap.kml"  # input
path_rasters    <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/rasters_1cm"  # input
path_models     <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/models_2cm"  # input
path_labels     <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/model_input_2cm/tls/label_lookup.csv"  # input

# set parameters
crs_raster_las <- "+proj=utm +zone=32 +ellps=WGS84 +units=m +vunits=m +no_defs"
used_rasters   <- c("ortho", "anisotropy_max", "curvature_max", "linearity_max",
                    "linearity_sd", "planarity_mean", "planarity_sd", "nDSM",
                    "point_density", "reflectance_mean", "reflectance_sd")

# use custom color palette
own_colors_named <- list(gray_green = "#96ceb4",
                         light_yellow = "#ffeead",
                         red = "#ff6f69",
                         blue = "#00aedb",
                         yellow = "#ffcc5c",
                         green = "#88d8b0",
                         turquoise = "#4abdac",
                         pink = "#ff8b94",
                         orange = "#ff9a00")
own_colors <- c(own_colors_named$yellow, own_colors_named$red, own_colors_named$green,
                own_colors_named$blue, own_colors_named$orange, own_colors_named$turquoise)

################################################################################
# RAW POINT CLOUDS
################################################################################

# maybe better in CC

################################################################################
# CLUSTERS & PCA
################################################################################

# cluster before & after

# PCA after (before too cluttered)

################################################################################
# RASTER STATISTICS & LABELS
################################################################################

# get raster & label values (takes some time)

if (FALSE) {
  # empty list for raster values
  values <- list()
  # load & transform plots
  vegetation <- st_transform(st_read(path_vegetation), crs_raster_las)
  labels <- unique(vegetation$Name)
  labels <- labels[!(labels %in% c("rock", "grass"))]
  # get all raster paths
  raster_files <- list.files(path_rasters, pattern="[.]tif", recursive = TRUE, full.names = TRUE)
  raster_files <- raster_files[!grepl("DTM", raster_files)]
  # define edge length
  tile_size <- 0.5
  edge <- ((tile_size*100)%/%2)/100
  # loop though raster types
  # raster_types <- unique(sapply(strsplit(basename(raster_files), "_area"), "[[", 1))
  raster_types <- used_rasters
  for (type in raster_types) {
    print(paste0("type: ", type))
    raster_type_paths <- raster_files[grepl(type, raster_files)]
    # loop through rasters
    for (path in raster_type_paths) {
      print(paste0("raster: ", path))
      raster <- stack(path)
      # loop through bands
      for (label in labels) {
        print(paste0("label: ", label))
        plots <- vegetation[vegetation$Name == label,]
        plots <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(raster), "SpatialPolygons")), st_crs(plots)))
        # continue if extents raster & plots are overlapping
        if (nrow(plots) > 0) {
          raster_res <- res(raster)[1]
          center_x <- round(st_coordinates(plots)[,1]/raster_res)*raster_res # round on resolution cm
          center_y <- round(st_coordinates(plots)[,2]/raster_res)*raster_res # round on resolution cm
          # loop through single plots
          for (i in 1:nrow(plots)) {
            if (edge %% 2 == 0) {
              rectangle <- extent(c(xmin=center_x[i]-edge, xmax=center_x[i]+edge, ymin=center_y[i]-edge, ymax=center_y[i]+edge))
            } else {  # if amount of pixels would be uneven --> rectangle center must be in a cell center
              offset <- raster_res / 2
              rectangle <- extent(c(xmin=center_x[i]-edge-offset, xmax=center_x[i]+edge-offset,
                                    ymin=center_y[i]-edge-offset, ymax=center_y[i]+edge-offset))
            }
            clip_vals <- values(crop(raster, rectangle))
            if (type != "ortho") {
              values[[paste0(type, " ", label)]] <- c(values[[paste0(type, " ", label)]], clip_vals)
            }
            else {
              values[[paste0("R ", label)]] <- c(values[[paste0("R ", label)]], clip_vals[,1])
              values[[paste0("G ", label)]] <- c(values[[paste0("G ", label)]], clip_vals[,2])
              values[[paste0("B ", label)]] <- c(values[[paste0("B ", label)]], clip_vals[,3])
            }
          }
        }
      }
    }
  }
  # save values
  raster_vals <- utils::stack(values)
  raster_vals$ind <- as.character(raster_vals$ind)
  raster_vals$type <- sapply(strsplit(raster_vals$ind, " "), "[[", 1)
  raster_vals$label <- sapply(strsplit(raster_vals$ind, " "), "[[", 2)
  raster_vals$ind <- NULL
  write.csv(raster_vals, paste0(path_rasters, "/raster_values_labels.csv"), row.names=F)
}

################################################################################

# boxplots for some rasters

raster_vals <- read.csv(paste0(path_rasters, "/raster_values_labels.csv"))
raster_vals <- na.omit(raster_vals)

# trying out some plots

ggplot(raster_vals[raster_vals$type == "nDSM",], aes(x = label, y = values)) +
  stat_boxplot(geom = 'errorbar', width = 0.15) +
  geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75) +
  xlab("\nVegetation Label") + ylab("nDSM Height\n") +
  scale_x_discrete(labels = c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")) +
  scale_fill_manual(values=own_colors) +
  theme_light() +
  theme(text = element_text(size=14, family="Calibri"), legend.position="none")

ggplot(raster_vals[raster_vals$type == "reflectance_mean",], aes(x = label, y = values)) +
  stat_boxplot(geom = 'errorbar', width = 0.15) +
  geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75) +
  xlab("\nVegetation Label") + ylab("Reflectance, mean\n") +
  scale_x_discrete(labels = c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")) +
  scale_fill_manual(values=own_colors) +
  theme_light() +
  theme(text = element_text(size=14, family="Calibri"), legend.position="none")

ggplot(raster_vals[raster_vals$type == "linearity_max",], aes(x = label, y = values)) +
  stat_boxplot(geom = 'errorbar', width = 0.15) +
  geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75) +
  xlab("\nVegetation Label") + ylab("Linearity, max\n") +
  scale_x_discrete(labels = c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")) +
  scale_fill_manual(values=own_colors) +
  theme_light() +
  theme(text = element_text(size=14, family="Calibri"), legend.position="none")

ggplot(raster_vals[raster_vals$type == "planarity_mean",], aes(x = label, y = values)) +
  stat_boxplot(geom = 'errorbar', width = 0.15) +
  geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75) +
  xlab("\nVegetation Label") + ylab("Planarity, mean\n") +
  scale_x_discrete(labels = c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")) +
  scale_fill_manual(values=own_colors) +
  theme_light() +
  theme(text = element_text(size=14, family="Calibri"), legend.position="none")

# TODO: arrange in a grid
# http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/

# all rgb values (R, G, B)

# all geometry values (anisotropy_max, curvature_max, linearity_max, linearity_sd, planarity_mean, planarity_sd)

# all tls values (nDSM, point_density, reflectance_mean, reflectance_sd)

################################################################################

# i want radarplots but they are not popular

################################################################################
# HYPERPARAMETERS & ACCURACY & LOSS
################################################################################

# correlation plot - hyperparameters & accuracy & loss

# barplots (with CI) / boxplots - hyperparameters & accuracy & loss

# regression between val_acc and val_loss
# (good if val_acc was used for picking network)

################################################################################
# CONFUSION MATRICES
################################################################################

# get prediction paths
csv_paths <- list.files(path_models, pattern="prediction_truth_fold.csv", recursive=TRUE, full.names=TRUE)
csv_path <- csv_paths[1]

# load predictions & label lookup table
pred_df <- read.csv(csv_path)
label_df <- read.csv(path_labels)

# transform numbers to labels
for (i in 1:nrow(label_df)) {
  pred_df$predictions[pred_df$predictions == label_df$new[i]] <- label_df$old[[i]]
  pred_df$truth[pred_df$truth == label_df$new[i]] <- label_df$old[[i]]
}

# select single fold ?
# pred_df <- pred_df[pred_df$fold == 4,]

# make confusion matrix
conf <- confusionMatrix(as.factor(pred_df$predictions), as.factor(pred_df$truth))
print(conf)

# extract table & statistics
conf_data <- data.frame(conf$table)
conf_stat <- round(data.frame(conf$overall),2)

# add percent values
conf_data$Prob <- round(data.frame(prop.table(conf$table))$Freq*100,2)

# plot confusion matrix
ggplot(data = conf_data, aes(x=Reference, y=Prediction, fill=Freq)) +
  geom_tile() + coord_equal() +
  geom_text(aes(label = Freq), color = 'gray30', size = 4) +
  xlab("\nTruth") + ylab("Prediction\n") +
  scale_x_discrete(labels = c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")) +
  scale_y_discrete(labels = c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")) +
  scale_fill_gradient(low="white", high=own_colors_named$blue) +
  theme_light() +
  theme(text = element_text(size=14, family="Calibri"), legend.position="none")

################################################################################
# ACCURACY COMPARISON INPUTS
################################################################################

# boxplots

# table, including sd

################################################################################
# VISUALIZE FILTERS
################################################################################

# Funktion zum Erzeugen von Filtervisualisierungen
generate_pattern <- function(layer_name, filter_index, size = 150) {  # size = Zoom
  layer_output <- model$get_layer(layer_name)$output
  loss <- k_mean(layer_output[,,,filter_index])
  grads <- k_gradients(loss, model$input)[[1]]
  grads <- grads / (k_sqrt(k_mean(k_square(grads))) + 1e-5)
  iterate <- k_function(list(model$input), list(loss, grads))
  input_img_data <- array(runif(size * size * 3),
                          dim = c(1, size, size, 3)) * 20 + 128
  step <- 1
  for (i in 1:40) {
    c(loss_value, grads_value) %<-% iterate(list(input_img_data))
    input_img_data <- input_img_data + (grads_value * step)
  }
  img <- input_img_data[1,,,]
  deprocess_image(img)
}

# Muster darstellen
grid.raster(generate_pattern("block3_conv1", 1))

################################################################################