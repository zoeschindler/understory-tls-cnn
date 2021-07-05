################################################################################
################################################################################
# VISUALIZATION & TABLES
################################################################################
################################################################################

# load packages
library(ggplot2)
library(ggpubr)
library(grid)
library(sf)
library(dplyr)
library(scales)
library(caret)
library(raster)
library(extrafont)
loadfonts(device="pdf", quiet=TRUE)

# set paths
path_vegetation <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/vegetation/Export_ODK_clean_checked_filtered_no_overlap.kml"  # input
path_rasters    <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters_1cm"  # input
path_models     <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/models_2cm_lenet5_10cv"  # input
path_labels     <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/model_input_2cm/tls/label_lookup.csv"  # input
path_plots      <- "H:/Daten/Studium/2_Master/4_Semester/5_Analyse/Plots"  # output

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
                         orange = "#ff9a00",
                         bright_green = "#cbe885")
own_colors <- c(own_colors_named$blue, own_colors_named$red, own_colors_named$yellow,
                own_colors_named$bright_green, own_colors_named$green)

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

raster_stat_plot <- function(data, y_label, raster_type, abbreviate=TRUE, log=FALSE) {
  if(abbreviate) {
    label_vector <- c("B", "D", "F", "M", "S")
  } else {
    label_vector <- c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")
  }
  if (!log) {
    plot <- ggplot(data[data$type == raster_type,], aes(x = label, y = values)) +
      stat_boxplot(geom = 'errorbar', width = 0.15) +
      geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75) +
      xlab("") + ylab("") +
      # xlab("\nVegetation Label") +
      ggtitle(y_label) +
      # ylab(paste0(y_label, "\n")) +
      scale_x_discrete(labels = label_vector) +
      scale_fill_manual(values=own_colors) +
      theme_light() +
      theme(text = element_text(size=14, family="Calibri"), legend.position="none",
            plot.title = element_text(hjust = 0.5))
  } else {
    data$values[data$type == raster_type & data$values == 0] <- NA
    data <- na.omit(data)
    plot <- ggplot(data[data$type == raster_type,], aes(x = label, y = values+1)) +
      stat_boxplot(geom = 'errorbar', width = 0.15) +
      geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75) +
      xlab("") + ylab("") +
      # xlab("\nVegetation Label") +
      ggtitle(y_label) +
      # ylab(paste0(y_label, "\n")) +
      scale_x_discrete(labels = label_vector) +
      scale_fill_manual(values=own_colors) +
      theme_light() +
      theme(text = element_text(size=14, family="Calibri"), legend.position="none",
            plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))
  }
  return(plot)
}

raster_legend <- function(data, pos) {
  label_vector <- c('Blueberry','Deadwood','Forest Floor', "Moss", "Spruce")
  plot <- ggplot(data[data$type == "nDSM",], aes(x = label, y = values)) +
    geom_boxplot(aes(fill = label)) +
    scale_fill_manual(values=own_colors, name = "Vegetation Label", labels = label_vector) +
    theme_light() +
    theme(text = element_text(size=16, family="Calibri"), legend.position=pos)
  legend <- get_legend(plot)
  return(legend)
}

################################################################################

# read data
raster_vals <- read.csv(paste0(path_rasters, "/raster_values_labels.csv"))
raster_vals <- na.omit(raster_vals)

# make legend
plot_legend_block <- raster_legend(raster_vals, "right")
plot_legend_line <- raster_legend(raster_vals, "bottom")

# all rgb values (R, G, B)
cairo_pdf(file = paste0(path_plots, "/rgb_raster_stats.pdf"), family = "Calibri", width = 8.27, height = 2.93)
plot_red   <- raster_stat_plot(raster_vals, "Red", "R")
plot_green <- raster_stat_plot(raster_vals, "Green", "G")
plot_blue  <- raster_stat_plot(raster_vals, "Blue", "B")
ggarrange(plot_red, plot_green, plot_blue,
          ncol = 3, nrow = 1,
          legend.grob = plot_legend_line, legend = "bottom")
dev.off()

# all geometry values (anisotropy_max, curvature_max, linearity_max, linearity_sd, planarity_mean, planarity_sd)
cairo_pdf(file = paste0(path_plots, "/geo_raster_stats.pdf"), family = "Calibri", width = 8.27, height = 5.83)
plot_aniso_max <- raster_stat_plot(raster_vals, "Anisotropy, max", "anisotropy_max")
plot_curv_max  <- raster_stat_plot(raster_vals, "Curvature, max", "curvature_max")
plot_linea_max <- raster_stat_plot(raster_vals, "Linearity, max", "linearity_max")
plot_linea_sd  <- raster_stat_plot(raster_vals, "Linearity, sd", "linearity_sd")
plot_plan_mean <- raster_stat_plot(raster_vals, "Planarity, mean", "planarity_mean")
plot_plan_sd   <- raster_stat_plot(raster_vals, "Planarity, sd", "planarity_sd")
ggarrange(plot_aniso_max, plot_curv_max, plot_linea_max,
          plot_linea_sd, plot_plan_mean, plot_plan_sd, 
          ncol = 3, nrow = 2,
          legend.grob = plot_legend_line, legend = "bottom")
dev.off()

# all tls values (nDSM, point_density, reflectance_mean, reflectance_sd)
# point density without 0s, because logarithmic scale hates that
cairo_pdf(file = paste0(path_plots, "/tls_raster_stats.pdf"), family = "Calibri", width = 8.27, height = 5.83)
plot_dens     <- raster_stat_plot(raster_vals, "Point Density", "point_density", log = TRUE)
plot_nDSM     <- raster_stat_plot(raster_vals, "nDSM Height", "nDSM")
plot_ref_mean <- raster_stat_plot(raster_vals, "Reflectance, mean", "reflectance_mean")
plot_ref_sd   <- raster_stat_plot(raster_vals, "Reflectance, sd", "reflectance_sd")
ggarrange(plot_dens, plot_nDSM,
          plot_ref_mean, plot_ref_sd,
          ncol = 2, nrow = 2,
          legend.grob = plot_legend_line, legend = "bottom")
dev.off()

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
# https://stackoverflow.com/questions/49664757/how-to-assign-multiple-colour-scales-to-a-dataset-in-ggplot2-assign-different-c

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
# pred_df <- pred_df[pred_df$fold == 5,]

# make confusion matrix
conf <- confusionMatrix(as.factor(pred_df$predictions), as.factor(pred_df$truth))
print(conf)

# extract table & statistics
conf_data <- data.frame(conf$table)
conf_stat <- round(data.frame(conf$overall),2)

# add percent values
conf_data$Prob <- round(data.frame(prop.table(conf$table))$Freq*100,2)

# reverse order of predictions
conf_data$Prediction <- factor(conf_data$Prediction, levels=rev(levels(conf_data$Prediction)))

# create color scale for each prediction category
color_scale <- c("blueberry" = own_colors_named$blue,
                 "dead_wood" = own_colors_named$blue,
                 "forest_floor" = own_colors_named$blue,
                 "moss" = own_colors_named$blue,
                 "spruce" = own_colors_named$blue)

# # scale values for color stuff
# conf_data <- conf_data %>%
#   group_by(Reference) %>%
#   mutate(Freq_scaled = Freq/max(Freq)) %>%
#   ungroup() %>%
#   arrange(Reference)

# percentage of true values of each label
conf_data <- conf_data %>%
  group_by(Reference) %>%
  mutate(Freq_prob=format(round(Freq/sum(Freq)*100,2), nsmall=2)) %>%
  ungroup() %>%
  arrange(Reference)

# # plot confusion matrix, no percentage
# ggplot(data = conf_data, aes(x=Reference, y=Prediction, fill=Reference, alpha=Freq_scaled)) +
#   geom_tile(fill = "white", alpha = 1) +
#   geom_tile(color = "gray50") + coord_equal() +
#   geom_text(aes(label = Freq), color = 'gray20', size = 4, alpha=1) +
#   xlab(paste0("\nReference (n = ", sum(conf_data$Freq), ")")) + ylab("Prediction\n") +
#   scale_x_discrete(labels = c(paste0("Blueberry\n(n = ", sum(conf_data$Freq[conf_data$Reference == "blueberry"]), ")"),
#                               paste0("Deadwood\n(n = ", sum(conf_data$Freq[conf_data$Reference == "dead_wood"]), ")"),
#                               paste0("Forest Floor\n(n = ", sum(conf_data$Freq[conf_data$Reference == "forest_floor"]), ")"),
#                               paste0("Moss\n(n = ", sum(conf_data$Freq[conf_data$Reference == "moss"]), ")"),
#                               paste0("Spruce\n(n = ", sum(conf_data$Freq[conf_data$Reference == "spruce"]), ")"))) +
#   scale_y_discrete(labels = rev(c("Blueberry","Deadwood","Forest Floor", "Moss", "Spruce"))) +
#   scale_fill_manual(values = color_scale) +
#   theme_light() +
#   theme(text = element_text(size=14), legend.position="none") +
#   ggtitle(paste0("Confusion Matrix, LeNet5, 2cm, all folds\nAccuracy ",
#                  round(conf$overall["Accuracy"]*100,2), "%"))

# plot confusion matrix, with percentage
ggplot(data = conf_data, aes(x=Reference, y=Prediction, fill=Reference, alpha=as.numeric(Freq_prob))) +
  geom_tile(fill = "white", alpha = 1) +
  geom_tile(color = "gray50") + coord_equal() +
  geom_text(aes(label = paste0(Freq)), color = 'gray20', size = 4, alpha=1) +
  geom_text(aes(label = paste0("\n\n ", gsub(" ", "", Freq_prob), "%")), color = 'gray20', size = 3, alpha=1) +
  xlab(paste0("\nReference (n = ", sum(conf_data$Freq), ")")) + ylab("Prediction\n") +
  scale_x_discrete(labels = c(paste0("Blueberry\n(n = ", sum(conf_data$Freq[conf_data$Reference == "blueberry"]), ")"),
                              paste0("Deadwood\n(n = ", sum(conf_data$Freq[conf_data$Reference == "dead_wood"]), ")"),
                              paste0("Forest Floor\n(n = ", sum(conf_data$Freq[conf_data$Reference == "forest_floor"]), ")"),
                              paste0("Moss\n(n = ", sum(conf_data$Freq[conf_data$Reference == "moss"]), ")"),
                              paste0("Spruce\n(n = ", sum(conf_data$Freq[conf_data$Reference == "spruce"]), ")"))) +
  scale_y_discrete(labels = rev(c("Blueberry","Deadwood","Forest Floor", "Moss", "Spruce"))) +
  scale_fill_manual(values = color_scale) +
  theme_light() +
  theme(text = element_text(size=14), legend.position="none") +
  ggtitle(paste0("Confusion Matrix, LeNet5, 2cm, all folds\nAccuracy ",
                 round(conf$overall["Accuracy"]*100,2), "%"))

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