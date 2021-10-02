################################################################################
################################################################################
# VISUALIZATION & TABLES
################################################################################
################################################################################

# load packages
library(ggplot2)
library(tfruns)
library(jsonlite)
library(ggpubr)
library(ggcorrplot) 
library(sf)
library(dplyr)
library(scales)
library(caret)
library(raster)
library(Hmisc)
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)
# library(ggdendro)

# set paths
basedir <- "H:/Daten/Studium/2_Master/4_Semester"
path_vegetation <- paste0(basedir, "/4_Daten/vegetation/Export_ODK_clean_checked_filtered_no_overlap.kml") # input
path_rasters    <- paste0(basedir, "/4_Daten/rasters_2cm") # input
path_models     <- paste0(basedir, "/4_Daten/02_testing/models_2cm") # input
path_pre_models <- paste0(basedir, "/4_Daten/01_hyperparameter/models_2cm") # input
path_tfruns     <- paste0(basedir, "/4_Daten/02_testing/tfruns_2cm") # input
path_labels     <- paste0(basedir, "/4_Daten/model_input_2cm_standardized/tls/label_lookup.csv") # input
path_plots      <- paste0(basedir, "/5_Analyse/Plots") # output
path_raster_val_before <- paste0(path_rasters, "/raster_samples_scaled.csv") # input
path_raster_val_after  <- paste0(path_rasters, "/raster_samples_scaled_noncollinear.csv") # input

# set parameters
crs_raster_las <- "+proj=utm +zone=32 +ellps=WGS84 +units=m +vunits=m +no_defs"
used_rasters <- c(
  "ortho", "anisotropy_max", "curvature_max", "linearity_max",
  "linearity_sd", "planarity_mean", "planarity_sd", "nDSM",
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
  bright_green = "#cbe885",
  light_blue = "#47DAFF",
  light_red = "#ff8985",
  light_yellow = "#ffd374",
  light_bright_green = "#d4ed9a",
  light_green = "#93dcb8"
)

# color scale for classes
color_scale_class <- c(
  own_colors_named$blue,
  own_colors_named$red,
  own_colors_named$yellow,
  own_colors_named$bright_green,
  own_colors_named$green
)

# color scale for classes
color_scale_class_light <- c(
  own_colors_named$light_blue,
  own_colors_named$light_red,
  own_colors_named$light_yellow,
  own_colors_named$light_bright_green,
  own_colors_named$light_green
)

# color scale for input data type
color_scale_type <- c(
  "tls" = own_colors_named$yellow,
  "tls_geo" = own_colors_named$green,
  "tls_rgb" = own_colors_named$blue,
  "tls_rgb_geo" = own_colors_named$red
)

# color scale for input data type
color_scale_type_light <- c(
  "tls" = own_colors_named$light_yellow,
  "tls_geo" = own_colors_named$light_green,
  "tls_rgb" = own_colors_named$light_blue,
  "tls_rgb_geo" = own_colors_named$light_red
)

# grey color scale for classes (for appendix?)
color_scale_grey <- c("grey90", "gray80", "gray70", "gray60", "gray50")

################################################################################
# CLUSTERS & PCA
################################################################################

# read in data
raster_val_before <- read.csv(path_raster_val_before)
names(raster_val_before) <- c(
  "Anisotropy, max", "Anisotropy, mean", "Anisotropy, sd",
  "Curvature, max", "Curvature, mean", "Curvature, sd",
  "Linearity, max", "Linearity, mean", "Linearity, sd",
  "Planarity, max", "Planarity, mean", "Planarity, sd",
  "Sphericity, max", "Sphericity, mean", "Sphericity, sd",
  "nDSM", "Point Density",
  "Reflectance, max", "Reflectance, mean", "Reflectance, sd"
)
raster_val_after <- read.csv(path_raster_val_after)
names(raster_val_after) <- c(
  "Anisotropy, max", "Curvature, max", "Linearity, max",
  "Linearity, sd", "Planarity, mean", "Planarity, sd", "nDSM",
  "Point Density", "Reflectance, mean", "Reflectance, sd"
)

# make correlation matrix
cor_matrix_before <- cor(raster_val_before, method = "spearman")
p_matrix_before <- cor_pmat(raster_val_before)
cor_matrix_after <- cor(raster_val_after, method = "spearman")
p_matrix_after <- cor_pmat(raster_val_after)

# corr plot before
cairo_pdf(
  file = paste0(path_plots, "/corr_before.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggcorrplot(cor_matrix_before,
  type = "lower", outline.col = "white",
  ggtheme = ggplot2::theme_light, sig.level = 0.05,
  colors = c(own_colors_named$blue, "gray90", own_colors_named$red),
  legend.title = "Spearman's\nCorrelation\nCoefficient\n", tl.cex = 14, tl.srt = 45
) +
  theme(
    text = element_text(size = 16, family = "Calibri"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 16),
    legend.box.spacing = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(1.5, "cm"),
    legend.title = element_text(size = 18)
  )
dev.off()

# corr plot after
cairo_pdf(
  file = paste0(path_plots, "/corr_after.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggcorrplot(cor_matrix_after,
  type = "lower", outline.col = "white",
  lab = TRUE, lab_col = "grey25", lab_size = 5,
  ggtheme = ggplot2::theme_light, sig.level = 0.05,
  colors = c(own_colors_named$blue, "gray90", own_colors_named$red),
  legend.title = "Spearman's\nCorrelation\nCoefficient\n", tl.cex = 16, tl.srt = 45
) +
  theme(
    text = element_text(size = 16, family = "Calibri"), plot.title = element_text(hjust = 0.5),
    legend.text = element_text(family = "Calibri", size = 16),
    legend.box.spacing = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(1.5, "cm"),
    legend.title = element_text(family = "Calibri", size = 18)
  )
dev.off()

# # with help of: https://stackoverflow.com/questions/68557415/dendrogram-with-labels-on-the-right-side
# tree <- hclust(as.dist(1 - cor_matrix_before**2))
# data <- ggdendro::dendro_data(tree)
# cairo_pdf(
#   file = paste0(path_plots, "/corr_cluster.pdf"),
#   family = "Calibri", width = 8.27, height = 5.83
# )
# ggplot() +
#   geom_blank() +
#   geom_segment(data = segment(data), aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
#   geom_hline(yintercept = 0.7 * 0.7, col = own_colors_named$red) +
#   scale_x_continuous(breaks = seq_along(data$labels$label), labels = data$labels$label, position = "top") +
#   scale_y_reverse(breaks = seq(0, 1, 0.2), labels = rev(seq(0, 1, 0.2))) +
#   coord_flip() +
#   theme_light() +
#   theme(
#     axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
#     axis.text.y = element_text(angle = 0, hjust = 1),
#     text = element_text(size = 16, family = "Calibri"),
#     panel.grid.minor = element_blank()
#   ) +
#   ylab(expression(paste("\nSpearman's ", rho**2))) +
#   xlab("")
# dev.off()

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
  raster_files <- list.files(path_rasters, pattern = "[.]tif", recursive = TRUE, full.names = TRUE)
  raster_files <- raster_files[!grepl("DTM", raster_files)]
  # define edge length
  tile_size <- 0.5
  edge <- ((tile_size * 100) %/% 2) / 100
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
        plots <- vegetation[vegetation$Name == label, ]
        plots <- st_intersection(plots, st_set_crs(st_as_sf(as(extent(raster), "SpatialPolygons")), st_crs(plots)))
        # continue if extents raster & plots are overlapping
        if (nrow(plots) > 0) {
          raster_res <- res(raster)[1]
          center_x <- round(st_coordinates(plots)[, 1] / raster_res) * raster_res # round on resolution cm
          center_y <- round(st_coordinates(plots)[, 2] / raster_res) * raster_res # round on resolution cm
          # loop through single plots
          for (i in 1:nrow(plots)) {
            if (edge %% 2 == 0) {
              rectangle <- extent(c(xmin = center_x[i] - edge, xmax = center_x[i] + edge, ymin = center_y[i] - edge, ymax = center_y[i] + edge))
            } else { # if amount of pixels would be uneven --> rectangle center must be in a cell center
              offset <- raster_res / 2
              rectangle <- extent(c(
                xmin = center_x[i] - edge - offset, xmax = center_x[i] + edge - offset,
                ymin = center_y[i] - edge - offset, ymax = center_y[i] + edge - offset
              ))
            }
            clip_vals <- values(crop(raster, rectangle))
            if (type != "ortho") {
              values[[paste0(type, " ", label)]] <- c(values[[paste0(type, " ", label)]], clip_vals)
            }
            else {
              values[[paste0("R ", label)]] <- c(values[[paste0("R ", label)]], clip_vals[, 1])
              values[[paste0("G ", label)]] <- c(values[[paste0("G ", label)]], clip_vals[, 2])
              values[[paste0("B ", label)]] <- c(values[[paste0("B ", label)]], clip_vals[, 3])
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
  write.csv(raster_vals, paste0(path_rasters, "/raster_values_labels.csv"), row.names = F)
}

################################################################################

f1 <- function(x) {
  ans <- boxplot.stats(x)
  data.frame(ymin = ans$conf[1], ymax = ans$conf[2], y = ans$stats[3])
}

f2 <- function(x) {
  ans <- boxplot.stats(x)
  data.frame(yintercept = as.numeric(c(ans$conf[1], ans$conf[2])))
}

raster_stat_plot <- function(data, y_label, raster_type, abbreviate = TRUE, log = FALSE, notch = FALSE) {
  # create plots
  if (abbreviate) {
    label_vector <- c("B", "D", "F", "M", "S")
  } else {
    label_vector <- c("Blueberry", "Deadwood", "Forest Floor", "Moss", "Spruce")
  }
  if (!log) {
    plot <- ggplot(data[data$type == raster_type, ], aes(x = label, y = values)) +
      stat_boxplot(geom = "errorbar", width = 0.25)
    if (notch) {
      plot <- plot +
        stat_summary(fun.data = f1, geom = "crossbar", colour = NA, fill = "red", width = 1, alpha = 1) # +
      # stat_summary(fun.data = f2, geom = "hline", colour = "gray50", linetype = "dotted", size = 0.5)
    }
    plot <- plot +
      geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75, width = 0.7) +
      xlab("") +
      ylab("") +
      ggtitle(y_label) +
      scale_x_discrete(labels = label_vector) +
      scale_fill_manual(values = color_scale_class) +
      theme_light() +
      theme(
        text = element_text(size = 14, family = "Calibri"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )
  } else {
    data$values[data$type == raster_type & data$values == 0] <- NA
    data <- na.omit(data)
    plot <- ggplot(data[data$type == raster_type, ], aes(x = label, y = values + 1)) +
      stat_boxplot(geom = "errorbar", width = 0.25)
    if (notch) {
      plot <- plot +
        stat_summary(fun.data = f1, geom = "crossbar", colour = NA, fill = "red", width = 1, alpha = 1) # +
      # stat_summary(fun.data = f2, geom = "hline", colour = "gray50", linetype = "dotted", size = 0.5)
    }
    plot <- plot +
      geom_boxplot(aes(fill = label), outlier.alpha = 0.01, outlier.size = 0.75) +
      xlab("") +
      ylab("") +
      ggtitle(y_label) +
      scale_x_discrete(labels = label_vector) +
      scale_fill_manual(values = color_scale_class) +
      theme_light() +
      theme(
        text = element_text(size = 14, family = "Calibri"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_y_continuous(
        trans = log10_trans(),
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = math_format(
          format = function(x) {number(log10(x), accuracy = 0.1)}
        )
      )
  }
  return(plot)
}

raster_legend <- function(data, pos, col = NA) {
  # create legend
  label_vector <- c("Blueberry", "Deadwood", "Forest Floor", "Moss", "Spruce")
  plot <- ggplot(data[data$type == "nDSM", ], aes(x = label, y = values)) +
    geom_boxplot(aes(fill = label)) +
    scale_fill_manual(values = color_scale_class, name = "Vegetation Label", labels = label_vector) +
    theme_light() +
    theme(
      legend.position = pos,
      legend.title = element_text(family = "Calibri", size = 16),
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(family = "Calibri", size = 14)
    )
  if (!is.na(col)) {
    plot <- plot + guides(fill = guide_legend(ncol = col))
  }
  legend <- get_legend(plot)
  return(legend)
}

################################################################################

# read data
raster_vals <- read.csv(paste0(path_rasters, "/raster_values_labels.csv"))
raster_vals <- na.omit(raster_vals)

# make legend
plot_legend_list <- raster_legend(raster_vals, "right", 2)
plot_legend_block <- raster_legend(raster_vals, "right")
plot_legend_line <- raster_legend(raster_vals, "bottom")

# all rgb values (R, G, B)
cairo_pdf(
  file = paste0(path_plots, "/rgb_raster_stats.pdf"),
  family = "Calibri", width = 8.27, height = 2.93
)
plot_red <- raster_stat_plot(raster_vals, "Red", "R")
plot_green <- raster_stat_plot(raster_vals, "Green", "G")
plot_blue <- raster_stat_plot(raster_vals, "Blue", "B")
ggarrange(plot_red, plot_green, plot_blue,
  ncol = 3, nrow = 1, legend.grob = plot_legend_line, legend = "bottom"
)
dev.off()

# all geometry values (anisotropy_max, curvature_max, linearity_max, linearity_sd, planarity_mean, planarity_sd)
cairo_pdf(
  file = paste0(path_plots, "/geo_raster_stats.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
plot_aniso_max <- raster_stat_plot(raster_vals, "Anisotropy, max", "anisotropy_max")
plot_curv_max <- raster_stat_plot(raster_vals, "Curvature, max", "curvature_max")
plot_linea_max <- raster_stat_plot(raster_vals, "Linearity, max", "linearity_max")
plot_linea_sd <- raster_stat_plot(raster_vals, "Linearity, sd", "linearity_sd")
plot_plan_mean <- raster_stat_plot(raster_vals, "Planarity, mean", "planarity_mean")
plot_plan_sd <- raster_stat_plot(raster_vals, "Planarity, sd", "planarity_sd")
ggarrange(plot_aniso_max, plot_curv_max, plot_linea_max,
  plot_linea_sd, plot_plan_mean, plot_plan_sd,
  ncol = 3, nrow = 2, legend.grob = plot_legend_line, legend = "bottom"
)
dev.off()

# all tls values (nDSM, point_density, reflectance_mean, reflectance_sd)
# point density without 0s, because logarithmic scale hates that
cairo_pdf(
  file = paste0(path_plots, "/tls_raster_stats.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
plot_dens <- raster_stat_plot(raster_vals, "Point Density", "point_density", log = TRUE)
plot_nDSM <- raster_stat_plot(raster_vals, "nDSM Height", "nDSM")
plot_ref_mean <- raster_stat_plot(raster_vals, "Reflectance, mean", "reflectance_mean")
plot_ref_sd <- raster_stat_plot(raster_vals, "Reflectance, sd", "reflectance_sd")
ggarrange(plot_dens, plot_nDSM, plot_ref_mean, plot_ref_sd,
  ncol = 2, nrow = 2, legend.grob = plot_legend_line, legend = "bottom"
)
dev.off()

# all together
cairo_pdf(
  file = paste0(path_plots, "/raster_stats_combo.pdf"),
  family = "Calibri", width = 8.27, height = 5.83 * 2 + 0.5
)
plot_dens <- raster_stat_plot(raster_vals, "Point Density", "point_density", log = TRUE)
plot_nDSM <- raster_stat_plot(raster_vals, "nDSM Height", "nDSM")
plot_ref_mean <- raster_stat_plot(raster_vals, "Reflectance, mean", "reflectance_mean")
plot_ref_sd <- raster_stat_plot(raster_vals, "Reflectance, sd", "reflectance_sd")
plot_red <- raster_stat_plot(raster_vals, "Red", "R")
plot_green <- raster_stat_plot(raster_vals, "Green", "G")
plot_blue <- raster_stat_plot(raster_vals, "Blue", "B")
plot_aniso_max <- raster_stat_plot(raster_vals, "Anisotropy, max", "anisotropy_max")
plot_curv_max <- raster_stat_plot(raster_vals, "Curvature, max", "curvature_max")
plot_linea_max <- raster_stat_plot(raster_vals, "Linearity, max", "linearity_max")
plot_linea_sd <- raster_stat_plot(raster_vals, "Linearity, sd", "linearity_sd")
plot_plan_mean <- raster_stat_plot(raster_vals, "Planarity, mean", "planarity_mean")
plot_plan_sd <- raster_stat_plot(raster_vals, "Planarity, sd", "planarity_sd")
ggarrange(
  ggarrange(plot_dens, NULL, plot_legend_list, ncol = 3, widths = c(2, 0.1, 0.9)),
  ggarrange(
    plot_nDSM, plot_ref_mean, plot_ref_sd,
    plot_aniso_max, plot_curv_max, plot_linea_max,
    plot_linea_sd, plot_plan_mean, plot_plan_sd,
    plot_red, plot_green, plot_blue,
    nrow = 4, ncol = 3),
  
  ncol = 1, nrow = 2, heights = c(1, 4)
)
dev.off()

################################################################################

# check for significant differences within groups
for (type in unique(raster_vals$type)) {
  subset <- raster_vals[raster_vals$type == type, ]
  subset$values <- scale(subset$values)
  print(type)
  print(ks.test(scale(subset$values), "pnorm"))
  # all significant -> no normal distribution -> can't use anova or t-test
  print(kruskal.test(subset$values ~ subset$label))
  # all significant -> significant differences between all groups
  print("---------------------------------------------------------------------")
}

################################################################################
# LABELS & AREAS
################################################################################

# load & prepare data
final_points <- as.data.frame(st_read(path_vegetation))
final_points$area <- as.numeric(substr(lapply(strsplit(final_points$Description, " "), "[[", 2), 1, 1))
final_points <- final_points[final_points$Name != "rock" & final_points$Name != "grass", ]

# total number
nrow(final_points)

# classes per area & total amount of points per area
for (i in 1:8) {
  print(paste0("area: ", i))
  print(paste0("length: ", length(final_points$Name[final_points$area == i])))
  print(summary(as.factor(final_points$Name[final_points$area == i])))
  print("-----")
}

# total amount of points per class
for (class in unique(final_points$Name)) {
  print(paste0("class: ", class))
  print(paste0("length: ", nrow(final_points[final_points$Name == class, ])))
  print("-----")
}

################################################################################
# LABELS & FOLDS
################################################################################

# load labels
label_df <- read.csv(path_labels)

# get & print the class distribution
for (i in 1:5) {
  tab <- readRDS(paste0(dirname(path_labels), "/images_fold_", i, ".rds"))
  for (i in 1:nrow(label_df)) {
    tab$label[tab$label == label_df$new[i]] <- label_df$old[[i]]
  }
  print(paste0("fold: ", i))
  print(summary(as.factor(tab$label)))
  print("---------------------------------------------------------------------")
}

################################################################################
# CONVERGENCE PLOTS & TRAINING TIME
################################################################################

metrics_all <- c()
types <- c("tls", "tls_geo", "tls_rgb", "tls_rgb_geo")

# loop through datasets
for (type in types) {
  # loop through folds
  for (fold in 1:5) {
    # get run
    run_path <- paste0(path_tfruns, "/", type, "/fold_", fold)
    run <- ls_runs(runs_dir = run_path, order = metric_val_accuracy)[1, ]
    run_name <- basename(run$run_dir)
    # load training time
    time_path <- paste0(run_path, "/", run_name, "/tfruns.d/properties")
    start <- names(read.delim(paste0(time_path, "/start")))
    start <- as.double(substr(start, 2, nchar(start)))
    end <- names(read.delim(paste0(time_path, "/end")))
    end <- as.double(substr(end, 2, nchar(end)))
    train_time <- abs(end - start)
    # load history
    json_path <- paste0(run_path, "/", run_name, "/tfruns.d/metrics.json")
    history <- fromJSON(json_path)
    # save validation loss
    metrics_all <- rbind(metrics_all, data.frame(
      "type" = type,
      "fold" = fold,
      "epoch" = 1:length(history$val_loss),
      "loss" = history$loss,
      "accuracy" = history$accuracy,
      "val_loss" = history$val_loss,
      "val_accuracy" = history$val_accuracy,
      "lr" = ifelse(run$flag_learning_rate == 0, 5e-5, run$flag_learning_rate),
      "bs" = run$flag_batch_size,
      "ep" = run$flag_epochs,
      "time" = train_time
    ))
  }
}

# get training time
times <- metrics_all %>%
  group_by(type, fold) %>%
  summarise(time = unique(time))
# in seconds
summary_time <- summary(times$time)
print(summary_time)
# as string min:sec
summary_time_str <- paste0(summary_time %/% 60, ":", round(summary_time %% 60))
names(summary_time_str) <- names(summary_time)
print(summary_time_str)

# make validation loss plot
cairo_pdf(
  file = paste0(path_plots, "/val_loss.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggplot(metrics_all) +
  geom_line(
    aes(
      x = epoch, y = val_loss, color = type, group = interaction(type, fold), linetype = as.character(lr)
    ),
    alpha = 1, stat = "smooth", size = 1
  ) +
  scale_color_manual(
    values = color_scale_type, labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL"),
    name = "Input Data\nCombination\n"
  ) +
  scale_linetype_manual(
    name = "\nLearning Rate\n", values = c("solid", "longdash"),
    guide = guide_legend(override.aes = list(linetype = c("solid", "dashed"), size = 0.5))
  ) +
  ylab("Validation Loss\n") +
  xlab("\nEpoch") +
  theme_light() +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.title = element_text(size = 16),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )
dev.off()

# make validation accuracy
cairo_pdf(
  file = paste0(path_plots, "/val_acc.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggplot(metrics_all) +
  geom_line(
    aes(
      x = epoch, y = val_accuracy, color = type, group = interaction(type, fold), linetype = as.character(lr)
    ),
    alpha = 1, stat = "smooth", size = 1
  ) +
  scale_color_manual(
    values = color_scale_type, labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL"),
    name = "Input Data\nCombination\n"
  ) +
  scale_linetype_manual(
    name = "\nLearning Rate\n", values = c("solid", "longdash"),
    guide = guide_legend(override.aes = list(linetype = c("solid", "dashed"), size = 0.5))
  ) +
  ylab("Validation Accuracy\n") +
  xlab("\nEpoch") +
  theme_light() +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.title = element_text(size = 16),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )
dev.off()

################################################################################
# CONFUSION MATRICES - SINGLE
################################################################################

conf_matrix_plot <- function(pred_data, fold, name) {
  # get data depending on fold / "all"
  if (fold != "all") {
    pred_data <- pred_data[pred_data$fold == fold, ]
  }
  # make confusion matrix
  conf <- confusionMatrix(as.factor(pred_data$predictions), as.factor(pred_data$truth), mode = "everything")
  conf_data <- data.frame(conf$table)
  # reverse order of predictions
  conf_data$Prediction <- factor(conf_data$Prediction, levels = rev(levels(conf_data$Prediction)))
  # create color scale for each prediction category
  color_scale <- c(
    "blueberry" = own_colors_named$blue,
    "dead_wood" = own_colors_named$blue,
    "forest_floor" = own_colors_named$blue,
    "moss" = own_colors_named$blue,
    "spruce" = own_colors_named$blue
  )
  # percentage of true values of each label
  conf_data <- conf_data %>%
    group_by(Reference) %>%
    mutate(Freq_prob = format(round(Freq / sum(Freq) * 100, 2), nsmall = 2)) %>%
    ungroup() %>%
    arrange(Reference)
  # plot confusion matrix, with percentage
  plot <- ggplot(data = conf_data, aes(x = Reference, y = Prediction, fill = Reference, alpha = as.numeric(Freq_prob))) +
    geom_tile(fill = "white", alpha = 1) +
    geom_tile(color = "gray50") +
    coord_equal() +
    geom_text(aes(label = paste0(Freq)), color = "gray20", size = 4, alpha = 1) +
    geom_text(aes(label = paste0("\n\n ", gsub(" ", "", Freq_prob), "%")), color = "gray20", size = 3, alpha = 1) +
    xlab(paste0("\nReference (n = ", sum(conf_data$Freq), ")")) +
    ylab(paste0("Prediction (n = ", sum(conf_data$Freq), ")\n")) +
    scale_x_discrete(labels = c(
      paste0("Blueberry\n(n = ", sum(conf_data$Freq[conf_data$Reference == "blueberry"]), ")"),
      paste0("Deadwood\n(n = ", sum(conf_data$Freq[conf_data$Reference == "dead_wood"]), ")"),
      paste0("Forest Floor\n(n = ", sum(conf_data$Freq[conf_data$Reference == "forest_floor"]), ")"),
      paste0("Moss\n(n = ", sum(conf_data$Freq[conf_data$Reference == "moss"]), ")"),
      paste0("Spruce\n(n = ", sum(conf_data$Freq[conf_data$Reference == "spruce"]), ")")
    )) +
    scale_y_discrete(labels = rev(c(
      paste0("Blueberry\n(n = ", sum(conf_data$Freq[conf_data$Prediction == "blueberry"]), ")"),
      paste0("Deadwood\n(n = ", sum(conf_data$Freq[conf_data$Prediction == "dead_wood"]), ")"),
      paste0("Forest Floor\n(n = ", sum(conf_data$Freq[conf_data$Prediction == "forest_floor"]), ")"),
      paste0("Moss\n(n = ", sum(conf_data$Freq[conf_data$Prediction == "moss"]), ")"),
      paste0("Spruce\n(n = ", sum(conf_data$Freq[conf_data$Prediction == "spruce"]), ")")
    ))) +
    scale_fill_manual(values = color_scale) +
    theme_light() +
    theme(text = element_text(size = 14, family = "Calibri"), legend.position = "none")
  # make legend
  dummy <- data.frame("x" = 0:100, "y" = 0:100, "z" = 0:100)
  legend <- ggplot(data = dummy, aes(x = x, y = y, alpha = z)) +
    geom_tile(fill = own_colors_named$blue) +
    geom_tile(aes(fill = z), alpha = 0) +
    # scale_fill_gradient(high = own_colors_named$blue, low = "#e6f7fb") +
    scale_fill_gradient(
      high = own_colors_named$blue, low = "#e6f7fb",
      labels = function(x) paste0(x, "%")
    ) +
    guides(alpha = "none") +
    labs(fill = "Fraction of\nReference\n(per Class)\n") +
    theme(
      legend.title = element_text(family = "Calibri", size = 16),
      legend.text.align = 1,
      # legend.box.spacing = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(1.5, "cm"),
      legend.text = element_text(family = "Calibri", size = 14)
    )
  legend <- get_legend(legend)
  # add title depending on fold / "all"
  if (fold != "all") {
    plot <- plot + ggtitle(paste0("Confusion Matrix, ", name, ", fold ", fold, "\nAccuracy ", round(conf$overall["Accuracy"] * 100, 2), "%"))
  } else {
    # plot <- plot + ggtitle(paste0("Confusion Matrix, ", name, ", all folds\nAccuracy ", round(conf$overall["Accuracy"] * 100, 2), "%"))
    plot <- plot +
      labs(
        title = paste0("Confusion Matrix: ", name),
        subtitle = paste0(
          "Overall Accuracy ", round(conf$overall["Accuracy"] * 100, 2),
          "%, Mean F1-Score ", round(mean(conf$byClass[, "F1"]), 2)
        )
      ) +
      theme(plot.title = element_text(face = "bold"))
  }
  # compose
  plot_final <- ggarrange(plot, legend, ncol = 2, nrow = 1, widths = c(6, 1)) # , legend.grob = legend, legend = "right")
  return(plot_final)
}

################################################################################

types <- c("tls", "tls_rgb", "tls_geo", "tls_rgb_geo")
names <- c("TLS", "TLS & RGB", "TLS & GEO", "ALL")

for (i in 1:4) {
  # set type and name
  type <- types[i]
  name <- names[i]

  # set paths
  csv_path <- paste0(path_models, "/", type, "/predictions_not_retrained.csv")

  # load predictions & label lookup table
  pred_df <- read.csv(csv_path)
  label_df <- read.csv(path_labels)

  # transform numbers to labels
  for (i in 1:nrow(label_df)) {
    pred_df$predictions[pred_df$predictions == label_df$new[i]] <- label_df$old[[i]]
    pred_df$truth[pred_df$truth == label_df$new[i]] <- label_df$old[[i]]
  }

  # make plots
  cairo_pdf(
    file = paste0(path_plots, "/conf_matrix_", type, ".pdf"),
    family = "Calibri", width = 8.27, height = 5.83
  )
  print(conf_matrix_plot(pred_df, "all", name))
  dev.off()
}

################################################################################
# CONFUSION MATRICES - COMBINED
################################################################################

conf_all <- c()
types <- c("tls", "tls_geo", "tls_rgb", "tls_rgb_geo")
names <- c("TLS", "TLS & GEO", "TLS & RGB", "ALL")

for (i in 1:4) {
  # set type and name
  type <- types[i]
  name <- names[i]

  # set paths
  csv_path <- paste0(path_models, "/", type, "/predictions_not_retrained.csv")

  # load predictions & label lookup table
  pred_df <- read.csv(csv_path)
  label_df <- read.csv(path_labels)

  # transform numbers to labels
  for (i in 1:nrow(label_df)) {
    pred_df$predictions[pred_df$predictions == label_df$new[i]] <- label_df$old[[i]]
    pred_df$truth[pred_df$truth == label_df$new[i]] <- label_df$old[[i]]
  }

  # make confusion matrix
  conf <- confusionMatrix(as.factor(pred_df$predictions), as.factor(pred_df$truth), mode = "everything")
  conf_data <- data.frame(conf$table)

  # percentage of true values of each label
  conf_data <- conf_data %>%
    group_by(Reference) %>%
    mutate(Prob = format(round(Freq / sum(Freq) * 100, 2), nsmall = 2)) %>%
    ungroup() %>%
    arrange(Reference)

  # get class distribution
  total_ref <- data.frame(
    "all" = sum(conf_data$Freq),
    "blueberry" = sum(conf_data$Freq[conf_data$Reference == "blueberry"]),
    "dead_wood" = sum(conf_data$Freq[conf_data$Reference == "dead_wood"]),
    "forest_floor" = sum(conf_data$Freq[conf_data$Reference == "forest_floor"]),
    "moss" = sum(conf_data$Freq[conf_data$Reference == "moss"]),
    "spruce" = sum(conf_data$Freq[conf_data$Reference == "spruce"])
  )

  # save in dataframe
  conf_all <- rbind(conf_all, cbind(conf_data, "Name" = type))
}

# reshape & extract data for plot
conf_all <- conf_all %>%
  group_by(Prediction, Reference) %>%
  summarise(
    max_freq = max(Freq),
    min_freq = min(Freq),
    max_prob = Prob[which.max(Freq)],
    min_prob = Prob[which.min(Freq)],
    max_n = sum(Freq == max(Freq)),
    min_n = sum(Freq == min(Freq)),
    max_name = paste0(Name[Freq == max(Freq)], collapse = " - "),
    min_name = paste0(Name[Freq == min(Freq)], collapse = " - ")
  )

# replace Name when multiple models
conf_all$max_name <- ifelse(conf_all$max_n > 1, "multiple", conf_all$max_name)
conf_all$min_name <- ifelse(conf_all$min_n > 1, "multiple", conf_all$min_name)

# reshape even more!
conf_all_final <- data.frame()
for (i in 1:nrow(conf_all)) {
  if (conf_all$Reference[i] == conf_all$Prediction[i]) {
    type <- c("worst", "best")
  } else {
    type <- c("best", "worst")
  }
  conf_all_final <- rbind(conf_all_final, data.frame(
    "Reference" = conf_all$Reference[i],
    "Prediction" = conf_all$Prediction[i],
    "Type" = type,
    "Data_Set" = c(conf_all$min_name[i], conf_all$max_name[i]),
    "Freq" = c(conf_all$min_freq[i], conf_all$max_freq[i]),
    "Prob" = c(conf_all$min_prob[i], conf_all$max_prob[i])
  ))
}
conf_all <- conf_all_final

# relevel data
conf_all$Prediction <- factor(conf_all$Prediction, levels = rev(levels(conf_all$Prediction)))

# class to number
conf_all$x <- as.numeric(factor(conf_all$Reference))
conf_all$y <- as.numeric(factor(conf_all$Prediction))
conf_all$id <- 1:nrow(conf_all)

# duplicate every row 3 times
conf_all <- conf_all[rep(conf_all$id, each = 3), ]

# get coordinates for lower polygon (worst values)
poly_x_worst <- conf_all$x + c(-0.5, 0.5, 0.5)
poly_y_worst <- conf_all$y + c(-0.5, -0.5, 0.5)

# get coordinates for upper polygon (best values)
poly_x_best <- conf_all$x + c(-0.5, -0.5, 0.5)
poly_y_best <- conf_all$y + c(-0.5, 0.5, 0.5)

# set polygon coordinates
conf_all$poly_x <- NA
conf_all$poly_y <- NA
conf_all$poly_x[conf_all$Type == "worst"] <- poly_x_worst[conf_all$Type == "worst"]
conf_all$poly_x[conf_all$Type == "best"] <- poly_x_best[conf_all$Type == "best"]
conf_all$poly_y[conf_all$Type == "worst"] <- poly_y_worst[conf_all$Type == "worst"]
conf_all$poly_y[conf_all$Type == "best"] <- poly_y_best[conf_all$Type == "best"]

# make & save plot
cairo_pdf(
  file = paste0(path_plots, "/conf_matrix_combo.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggplot(conf_all) +
  geom_tile(aes(x = Reference, y = Prediction), fill = "white") +
  geom_polygon(aes(x = poly_x, y = poly_y, group = id, fill = Data_Set), size = 1, color = "white") +
  geom_tile(aes(x = Reference, y = Prediction), alpha = 0, col = "white", size = 1.5) +
  geom_text(
    data = conf_all[conf_all$Type == "best", ],
    aes(x = x - 0.25, y = y + 0.15, label = Freq),
    size = 4, check_overlap = TRUE
  ) +
  geom_text(
    data = conf_all[conf_all$Type == "worst", ],
    aes(x = x + 0.25, y = y - 0.15, label = Freq),
    size = 4, check_overlap = TRUE
  ) +
  geom_text(
    data = conf_all[conf_all$Type == "best", ],
    aes(x = x - 0.25, y = y + 0.2, label = paste0(Type, "\n\n")),
    size = 3, check_overlap = TRUE
  ) +
  geom_text(
    data = conf_all[conf_all$Type == "worst", ],
    aes(x = x + 0.25, y = y - 0.2, label = paste0("\n\n", Type)),
    size = 3, check_overlap = TRUE
  ) +
  scale_fill_manual(values = c(color_scale_type, "multiple" = "gray80"), labels = c(names, "Multiple"), name = "Input Data\nCombination\n") +
  scale_x_discrete(labels = c(
    paste0("Blueberry\n(n = ", total_ref$blueberry, ")"),
    paste0("Deadwood\n(n = ", total_ref$dead_wood, ")"),
    paste0("Forest Floor\n(n = ", total_ref$forest_floor, ")"),
    paste0("Moss\n(n = ", total_ref$moss, ")"),
    paste0("Spruce\n(n = ", total_ref$spruce, ")")
  )) +
  scale_y_discrete(labels = rev(c("Blueberry", "Deadwood", "Forest Floor", "Moss", "Spruce"))) +
  theme_light() +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.title = element_text(size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.text = element_text(size = 14)
  ) +
  ylab("Prediction\n") +
  xlab(paste0("\nReference (n = ", total_ref$all, ")")) +
  coord_equal()
dev.off()

################################################################################
# ACCURACY INSTABILITY - NO HYPERPARAMETERS
################################################################################

# set paths & load data
fluct_all <- c()
for (type in c("tls", "tls_rgb", "tls_geo", "tls_rgb_geo")) {
  path_model_type <- paste0(path_models, "/", type)
  fluct_path <- paste0(path_model_type, "/accuracy_fluctuations_best_hyperparameters.csv")
  fluct_df <- read.csv(fluct_path)
  for (i in 1:5) {
    fluct_all <- rbind(fluct_all, data.frame(
      "type" = type,
      "fold" = i,
      "test_acc" = fluct_df[, i]
    ))
  }
}

# boxplots, accuracy fluctuations, per fold
cairo_pdf(
  file = paste0(path_plots, "/fluctuations_acc_varied.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggplot(
  fluct_all, aes(
    x = type, y = test_acc, alpha = fold,
    fill = type, group = interaction(type, fold)
  )
) +
  stat_boxplot(geom = "errorbar", alpha = 1, width = 0.3, position = position_dodge(0.9)) +
  geom_boxplot(outlier.alpha = 0, outlier.size = 0, fill = "white", alpha = 1, position = position_dodge(0.9)) +
  geom_boxplot(outlier.alpha = 1, outlier.size = 1, position = position_dodge(0.9)) +
  scale_fill_manual(
    values = color_scale_type,
    name = "Input Data\nCombination",
    labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL")
  ) +
  theme_light() +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.title = element_text(family = "Calibri", size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.text = element_text(family = "Calibri", size = 14)
  ) +
  scale_x_discrete(labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL")) +
  xlab("") +
  ylab("Test Accuracy\n") +
  scale_alpha_continuous(
    range = c(0.3, 1), name = "\nInput Data\nTest Fold",
    guide = guide_legend(override.aes = list(
      fill = "grey60",
      alpha = seq(0.3, 1, length.out = 5)
    ))
  )
dev.off()

# boxplots, accuracy fluctuations, overall
cairo_pdf(
  file = paste0(path_plots, "/fluctuations_acc_overall.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggplot(
  fluct_all, aes(
    x = type, y = test_acc, fill = type
  )
) +
  stat_boxplot(geom = "errorbar", alpha = 1, width = 0.2, position = position_dodge(0.9)) +
  geom_boxplot(outlier.alpha = 1, outlier.size = 1, position = position_dodge(0.9)) +
  scale_fill_manual(
    values = color_scale_type,
    name = "Input Data\nCombination",
    labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL")
  ) +
  theme_light() +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.title = element_text(size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.text = element_text(size = 14)
  ) +
  scale_x_discrete(labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL")) +
  xlab("") +
  ylab("Test Accuracy\n")
dev.off()

################################################################################
# ACCURACY INSTABILITY - HYPERPARAMETERS
################################################################################

make_fluct_hyper_plot <- function(data, name, expression) {
  plot <- ggplot(data, aes(x = type, y = test_acc, group = interaction(type, fold))) +
    stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(0.9)) +
    geom_boxplot(
      aes(fill = as.factor(eval(expression))),
      outlier.alpha = 1, outlier.size = 1, position = position_dodge(0.9)
    ) +
    scale_fill_manual(
      values = c(own_colors_named$red, own_colors_named$blue, own_colors_named$yellow),
      name = name
    ) +
    theme_light() +
    theme(
      text = element_text(size = 14, family = "Calibri"),
      legend.title = element_text(family = "Calibri", size = 16),
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(family = "Calibri", size = 14),
      legend.position = "bottom"
    ) +
    scale_x_discrete(labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL")) +
    xlab("") +
    ylab("Test Accuracy\n") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted", color = "gray30")
  return(plot)
}

################################################################################

# get all best hyperparameter combinations
all_confs <- c()
all_accs <- c()
for (type in c("tls", "tls_geo", "tls_rgb", "tls_rgb_geo")) {
  path_model_type <- paste0(path_models, "/", type)
  # save hyperparameters
  best_runs <- read.csv(paste0(path_model_type, "/best_hyperparameters.csv"))
  best_runs$flag_learning_rate[best_runs$flag_learning_rate == 0] <- 5e-5
  all_confs <- rbind(all_confs, data.frame(
    "type" = type,
    "fold" = best_runs$fold,
    "lr" = best_runs$flag_learning_rate,
    "bs" = best_runs$flag_batch_size,
    "ep" = best_runs$flag_epochs
  ))
  # save acccuracies
  accs <- read.csv(paste0(path_model_type, "/accuracy_fluctuations_best_hyperparameters.csv"))
  for (i in 1:5) {
    all_accs <- rbind(all_accs, data.frame(
      "type" = type,
      "fold" = i,
      "test_acc" = accs[, i]
    ))
  }
}

# join data together
accs_hypers <- all_accs %>%
  left_join(all_confs, by = c("type", "fold"))

################################################################################

# reshape data
var_df <- accs_hypers %>%
  group_by(type, fold) %>%
  summarise(
    "sd_acc" = sd(test_acc),
    "lr" = unique(lr),
    "bs" = unique(bs),
    "ep" = unique(ep)
  )
# make correlation
corr_m_p <- rcorr(as.matrix(var_df[, c("sd_acc", "lr", "bs", "ep")]), type = "spearman")
print(round(corr_m_p$r, 2))
print(round(corr_m_p$P, 2) < 0.05)
# only small correlations (0.1 - 0.26) -> no "big" influence
# -> must be the data, not the hyperparameters!

# reshape again
acc_min_max_diff <- accs_hypers %>%
  group_by(type, fold) %>%
  summarise(min = round(min(test_acc),2),
            max = round(max(test_acc),2),
            diff = round(abs(max(test_acc) - min(test_acc)),2))
print(paste0("max diff: ", max(acc_min_max_diff$diff)))
# biggest diff within one fold: 0.18

################################################################################

# make plot - LR
cairo_pdf(
  file = paste0(path_plots, "/fluct_acc_lr.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
make_fluct_hyper_plot(accs_hypers, "Learning Rate", expression(lr))
dev.off()

# make plot - BS
cairo_pdf(
  file = paste0(path_plots, "/fluct_acc_bs.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
make_fluct_hyper_plot(accs_hypers, "Batch Size", expression(bs))
dev.off()

# make plot - EP
cairo_pdf(
  file = paste0(path_plots, "/fluct_acc_ep.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
make_fluct_hyper_plot(accs_hypers, "Epochs", expression(ep))
dev.off()

################################################################################
# ACCURACY COMPARISON INPUTS
################################################################################

f1 <- function(x) {
  ans <- boxplot.stats(x)
  data.frame(ymin = ans$conf[1], ymax = ans$conf[2], y = ans$stats[3])
}

f2 <- function(x) {
  ans <- boxplot.stats(x)
  data.frame(yintercept = as.numeric(c(ans$conf[1], ans$conf[2])))
}

measure_boxplot <- function(data, measure, name, legend = "none", notch = FALSE) {
  data$measure <- data[, measure]
  plot <- ggplot(data, aes(x = type, y = measure)) +
    stat_boxplot(geom = "errorbar", width = 0.25)
  if (notch) {
    plot <- plot +
      stat_summary(fun.data = f1, geom = "crossbar", colour = NA, fill = "black", width = 0.85, alpha = 0.25) +
      stat_summary(fun.data = f2, geom = "hline", colour = "gray50", linetype = "dotted", size = 0.5)
  }
  plot <- plot +
    geom_boxplot(aes(fill = type), outlier.alpha = 1, outlier.size = 1.5, width = 0.7) +
    scale_fill_manual(
      values = color_scale_type,
      name = "Input Data\nCombination",
      labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL")
    ) +
    theme_light() +
    theme(
      text = element_text(size = 14, family = "Calibri"),
      legend.title = element_text(size = 16),
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(size = 14),
      legend.position = legend,
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(labels = c("TLS", "TLS & GEO", "TLS & RGB", "ALL")) +
    xlab("") +
    ylab(paste0("\n", name, "\n"))
  return(plot)
}

################################################################################

# create empty objects for storage
measures_all <- c()
f1_classes <- c()
conf_m_list <- list()
preds_list <- list()

# loop through types
for (type in c("tls", "tls_rgb", "tls_geo", "tls_rgb_geo")) {
  # load predictions
  path_model_type <- paste0(path_models, "/", type)
  preds_path <- paste0(path_model_type, "/predictions_not_retrained.csv")
  preds_df <- read.csv(preds_path)
  # replace numbers with labels
  label_lookup <- read.csv(path_labels)
  for (i in 1:nrow(label_lookup)) {
    preds_df$predictions[preds_df$predictions == label_lookup$new[i]] <- label_lookup$old[[i]]
    preds_df$truth[preds_df$truth == label_lookup$new[i]] <- label_lookup$old[[i]]
  }
  preds_list[[type]] <- preds_df
  # loop through folds
  for (i in 1:5) {
    # extract measures
    sub <- preds_df[preds_df$fold == i, ]
    conf_m <- confusionMatrix(as.factor(sub$predictions), as.factor(sub$truth), mode = "everything")
    conf_m_list[[paste0(type, "_", i)]] <- conf_m
    measures_all <- rbind(measures_all, data.frame(
      "type" = type,
      "fold" = i,
      "accuracy" = as.numeric(conf_m$overall["Accuracy"]),
      "kappa" = as.numeric(conf_m$overall["Kappa"]),
      "mean_f1" = mean(conf_m$byClass[, "F1"])
    ))
    f1_classes <- rbind(f1_classes, data.frame(
      "type" = type,
      "fold" = i,
      "class" = c("blueberry", "dead_wood", "forest_floor", "moss", "spruce"),
      "f1" = as.numeric(conf_m$byClass[,"F1"])
    ))
  }
}

################################################################################

# make plots for measures
plot_acc       <- measure_boxplot(measures_all, "accuracy", "Test Accuracy", notch = FALSE)
plot_acc_notch <- measure_boxplot(measures_all, "accuracy", "Test Accuracy", notch = TRUE)
plot_f1        <- measure_boxplot(measures_all, "mean_f1", "Mean F1-Score", notch = FALSE)
plot_f1_notch  <- measure_boxplot(measures_all, "mean_f1", "Mean F1-Score", notch = TRUE)
measures_legend_bottom <- get_legend(measure_boxplot(measures_all, "accuracy", "Test Accuracy", "bottom"))
measures_legend_right  <- get_legend(measure_boxplot(measures_all, "accuracy", "Test Accuracy", "right"))

# combine plots, f1 and accuracy, notch
cairo_pdf(
  file = paste0(path_plots, "/final_results_acc_f1_notch.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggarrange(plot_acc_notch, plot_f1_notch,
  ncol = 2, nrow = 1,
  legend.grob = measures_legend_bottom, legend = "bottom"
)
dev.off()

# combine plots, f1 and accuracy, no notch
cairo_pdf(
  file = paste0(path_plots, "/final_results_acc_f1.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggarrange(plot_acc, plot_f1,
  ncol = 2, nrow = 1,
  legend.grob = measures_legend_bottom, legend = "bottom"
)
dev.off()

# combine plots, accuracy, notch
cairo_pdf(
  file = paste0(path_plots, "/final_results_acc_notch.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
measure_boxplot(measures_all, "accuracy", "Test Accuracy", legend = "right", notch = TRUE)
dev.off()

# combine plots, f1, notch
cairo_pdf(
  file = paste0(path_plots, "/final_results_f1_notch.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
measure_boxplot(measures_all, "mean_f1", "Mean F1-Score", legend = "right", notch = TRUE)
dev.off()

################################################################################

# calculate data for table
accuracies <- c()
kappas <- c()
f1_mean <- c()
f1_classes <- list()
for (type in c("tls", "tls_rgb", "tls_geo", "tls_rgb_geo")) {
  accuracy <- c()
  kappa <- c()
  f1_score <- c()
  for (i in 1:5) {
    conf_m <- conf_m_list[[paste0(type, "_", i)]]
    accuracy <- c(accuracy, conf_m$overall[["Accuracy"]])
    kappa <- c(kappa, conf_m$overall[["Kappa"]])
    f1_score <- rbind(f1_score, conf_m$byClass[, "F1"])
  }
  accuracies[type] <- mean(accuracy)
  kappas[type] <- mean(kappa)
  f1_classes[[type]] <- apply(f1_score, 2, mean)
  f1_mean[type] <- mean(f1_classes[[type]])
}

# arrange results in table
output_results <- data.frame(
  "Mean_Measure" = c(
    "F1-Score, Blueberry", "F1-Score, Deadwood", "F1-Score, Forest Floor",
    "F1-Score, Moss", "F1-Score, Spruce", "Mean F1-Score", "Accuracy", "Kappa"
  ),
  "TLS" = c(f1_classes[["tls"]], f1_mean["tls"], accuracies["tls"], kappas["tls"]),
  "TLS_GEO" = c(f1_classes[["tls_geo"]], f1_mean["tls_geo"], accuracies["tls_geo"], kappas["tls_geo"]),
  "TLS_RGB" = c(f1_classes[["tls_rgb"]], f1_mean["tls_rgb"], accuracies["tls_rgb"], kappas["tls_rgb"]),
  "ALL" = c(f1_classes[["tls_rgb_geo"]], f1_mean["tls_rgb_geo"], accuracies["tls_rgb_geo"], kappas["tls_rgb_geo"])
)
output_results[, 2:5] <- round(output_results[, 2:5], 3)
print(output_results)

################################################################################
# VISUALIZE FILTERS
################################################################################

# # load packages
# library(grid)
# library(gridExtra)
#
# # turn off eager execution
# library(tensorflow)
# tf$compat$v1$disable_eager_execution()
#
# # load (random) model
# model <- load_model_hdf5(paste0(basedir, "/4_Daten/02_testing/models_2cm/tls/all/fold_1_lr_1e-05_ep_200_bs_32_dp_0.5_l2_0.01.h5"))
# nbands <- 4
#
# # tensor to picture
# deprocess_image <- function(x) {
#   dms <- dim(x)
#   x <- x - mean(x)
#   x <- x / (sd(x) + 1e-5)
#   x <- x * 0.1
#   x <- x + 0.5
#   x <- pmax(0, pmin(x, 1))
#   array(x, dim = dms)
# }
#
# # filter visualizations
# generate_pattern <- function(layer_name, filter_index, size = 150, nbands = 4) {
#   layer_output <- model$get_layer(layer_name)$output
#   loss <- k_mean(layer_output[,,,filter_index])
#   grads <- k_gradients(loss, model$input)[[1]]
#   grads <- grads / (k_sqrt(k_mean(k_square(grads))) + 1e-5)
#   iterate <- k_function(list(model$input), list(loss, grads))
#   input_img_data <- array(runif(size * size * nbands), dim = c(1, size, size, nbands)) * 20 + 128
#   step <- 1
#   for (i in 1:40) {
#     c(loss_value, grads_value) %<-% iterate(list(input_img_data))
#     input_img_data <- input_img_data + (grads_value * step)
#   }
#   img <- input_img_data[1,,,]
#   deprocess_image(img)
# }
#
# # saves images in working directory
# dir.create("filters")
# for (layer_name in c("conv2d_15")) {
#   size <- 100
#   png(paste0("filters/", layer_name, ".png"), width = 12 * size, height = 4 * size)
#   grobs <- list()
#   for (i in 0:7) {
#     for (j in 0:5) {
#       pattern <- generate_pattern(layer_name, i + (j*6) + 1, size = size)
#       grob <- rasterGrob(pattern[,,1:3],
#                          width = unit(0.9, "npc"),
#                          height = unit(0.9, "npc"))
#       grobs[[length(grobs)+1]] <- grob
#     }
#   }
#   grid.arrange(grobs = grobs, ncol = 8)
#   dev.off()
# }

# no clear patterns (on the raster there are none visible either, so that's not surprising)
# hard to visualize because of the many bands

################################################################################
