################################################################################
################################################################################
# PREPARE CNN INPUT
################################################################################
################################################################################

# load packages
library(sf)
library(raster)

# set paths
basedir <- "C:/Users/Zoe/Documents/understory_classification"
path_clips  <- paste0(basedir, "/4_Daten/clips_2cm_unscaled") # input
path_output <- paste0(basedir, "/4_Daten/model_input_2cm_standardized") # output
path_plot   <- paste0(basedir, "/4_Daten/vegetation/Export_ODK_clean_checked_filtered_no_overlap.kml") # input
path_values <- paste0(path_clips, "/raster_values_labels_unscaled.csv") # output

# set parameters
crs_raster_las <- "+proj=utm +zone=32 +ellps=WGS84 +units=m +vunits=m +no_defs"
raster_amount <- length(c(
  "ortho", "anisotropy_max", "curvature_max", "linearity_max",
  "linearity_sd", "planarity_mean", "planarity_sd", "nDSM",
  "point_density", "reflectance_mean", "reflectance_sd"
))

# set folder names for each input type
tls_set <- c("nDSM", "reflectance", "point_density")
rgb_set <- c("ortho")
geo_set <- c("geometry")
all_set <- c(tls_set, rgb_set, geo_set)

# set names and folder names for each combination
combo_1 <- c()
combo_1$folders <- c(tls_set)
combo_1$name <- "tls"
#
combo_2 <- c() # maybe this combination is unnecessary
combo_2$folders <- c(tls_set, rgb_set)
combo_2$name <- "tls_rgb"
#
combo_3 <- c()
combo_3$folders <- c(tls_set, geo_set)
combo_3$name <- "tls_geo"
#
combo_4 <- c()
combo_4$folders <- c(tls_set, rgb_set, geo_set)
combo_4$name <- "tls_rgb_geo"

################################################################################
# HELPER FUNCTIONS
################################################################################

check_create_dir <- function(path) {
  # checks if directory exists
  # if not, creates it
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

################################################################################

get_list_full_plots <- function(clip_dir, raster_amount) {
  # returns list of plots for which all types of rasters exist
  # get all files in clip folder
  clip_list <- list.files(clip_dir, pattern = "[.]tif", recursive = TRUE, full.names = TRUE)
  # extract all unique plot numbers
  plot_IDs <- c()
  for (clip_path in clip_list) {
    parts <- strsplit(strsplit(clip_path, "[.]")[[1]][1], "_")[[1]]
    plot <- as.numeric(parts[length(parts)])
    plot_IDs <- c(plot_IDs, plot)
  }
  plot_IDs <- sort(unique(plot_IDs))
  # loop through plot numbers and save those with all rasters
  full_plot_IDs <- c()
  for (plot_ID in plot_IDs) {
    n_raster <- length(clip_list[grepl(paste0("_", plot_ID, ".tif"), clip_list)])
    if (n_raster == raster_amount) {
      full_plot_IDs <- c(full_plot_IDs, plot_ID)
    }
  }
  return(full_plot_IDs)
}

################################################################################

clip_value_csv <- function(plot_path, clip_dir, output_path) {
  # create csv with all unscaled raster clip values
  # empty vector for raster values
  values <- c()
  # load plots
  vegetation <- st_read(plot_path)
  labels <- unique(vegetation$Name)
  labels <- labels[!(labels %in% c("rock", "grass"))]
  # get all raster paths
  clip_paths <- list.files(clip_dir, pattern = "[.]tif", recursive = TRUE, full.names = TRUE)
  clip_paths <- clip_paths[!grepl("DTM", clip_paths)]
  # make list with plot_ID & vegetation_ID and labels
  label_ID_lookup <- list()
  for (label in labels) {
    label_subset <- vegetation$Description[vegetation$Name == label]
    plot_IDs <- as.numeric(lapply(strsplit(label_subset, "veg_ID: "), tail, n = 1))
    label_ID_lookup[[label]] <- plot_IDs
  }
  # loop through clip types
  clip_types <- unique(sapply(strsplit(basename(clip_paths), "_area"), "[[", 1))
  for (type in clip_types) {
    print(paste0("type: ", type))
    clip_type_paths <- clip_paths[grepl(type, clip_paths)]
    for (label in labels) {
      print(paste0("label: ", label))
      label_IDs <- label_ID_lookup[[label]]
      label_clip_type_paths <- clip_type_paths[grepl(paste(paste0("_", label_IDs, "[.]tif"), collapse = "|"), clip_type_paths)]
      # loop through clips
      for (path in label_clip_type_paths) {
        raster <- stack(path)
        clip_label_vals <- values(raster)
        if (type != "ortho") {
          values <- rbind(values, data.frame("type" = type, "label" = label, "values" = as.numeric(clip_label_vals)))
        } else {
          values <- rbind(values, data.frame("type" = "R", "label" = label, "values" = as.numeric(clip_label_vals[, 1])))
          values <- rbind(values, data.frame("type" = "G", "label" = label, "values" = as.numeric(clip_label_vals[, 2])))
          values <- rbind(values, data.frame("type" = "B", "label" = label, "values" = as.numeric(clip_label_vals[, 3])))
        }
      }
    }
  }
  # save values
  write.csv(values, output_path, row.names = F)
  return(values)
}

rescale_values <- function(values_path) {
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
# STACK & SAVE & SORT CLIPS
################################################################################

stack_save_clips <- function(clip_dir, plot_path, values_path, output_dir, selection_rasters, selection_IDs, normalize = TRUE, standardize = FALSE) {
  # stacks all rasters for each clip & saves it in a folder based on the label
  check_create_dir(output_dir)
  # load lookup table
  if (normalize & standardize) {
    stop("either normalize OR standardize")
  }
  if (normalize | standardize) {
    lookup_table <- rescale_values(values_path)
  }
  # get all clip paths of the selected raster types
  clip_list <- list.files(clip_dir, pattern = "[.]tif", recursive = TRUE)
  contained <- c()
  for (clip_path in clip_list) {
    contained <- c(contained, strsplit(clip_path, "/")[[1]][1] %in% selection_rasters)
  }
  clip_list <- clip_list[contained]
  # load vegetation plots & create one subfolder for each label
  plots <- st_read(plot_path)
  # extract all plot IDs from vegetation plots
  all_IDs <- c()
  for (idx in 1:nrow(plots)) {
    all_IDs <- c(all_IDs, as.numeric(strsplit(plots$Description[idx], " ")[[1]][4]))
  }
  # loop through all selected plot IDs
  for (plot_ID in selection_IDs) {
    # get row from vegetation plots where this plot ID is same
    plot_idx <- which(all_IDs == plot_ID)
    # load & stack all clips of the plot
    plot_paths <- sort(clip_list[grepl(paste0("_", plot_ID, ".tif"), clip_list)])
    plot_names <- unlist(lapply(strsplit(basename(plot_paths), "_area"), head, n = 1))
    plot_rasters <- list()
    for (idx in 1:length(plot_paths)) {
      raster_bands <- stack(paste0(clip_dir, "/", plot_paths[idx]))
      # rescale values
      bands <- if (plot_names[idx] != "ortho") c(plot_names[idx]) else c("R", "G", "B")
      for (i in 1:length(bands)) {
        if (normalize | standardize) {
          if (normalize) {
            min_val <- lookup_table[[paste0(bands[i], "_min")]]
            max_val <- lookup_table[[paste0(bands[i], "_max")]]
            raster_bands[[i]] <- (raster_bands[[i]] - min_val) / (max_val - min_val)
          }
          if (standardize) {
            mean_val <- lookup_table[[paste0(bands[i], "_mean")]]
            sd_val <- lookup_table[[paste0(bands[i], "_sd")]]
            raster_bands[[i]] <- scale(raster_bands[[i]], center = mean_val, scale = sd_val)
          }
        }
      }
      plot_rasters[[idx]] <- raster_bands
    }
    plot_stack <- stack(plot_rasters)
    # determine label & save in according folder
    plot_label <- plots$Name[plot_idx]
    writeRaster(plot_stack, paste0(output_dir, "/", plot_label, "_", plot_ID, ".tif"))
  }
}

################################################################################
# EXECUTION
################################################################################

# create output folder
check_create_dir(path_output)

# save clip raster values in a file for rescaling, takes a long time, do once
# clip_value_csv(path_plot, path_clips, path_values)

# get plot_IDs where all rasters are available
full_IDs <- get_list_full_plots(path_clips, raster_amount)

# scale, stack, save raster clips
stack_save_clips(path_clips, path_plot, path_values, paste0(path_output, "/", combo_1$name),
                 combo_1$folders, full_IDs, normalize = FALSE, standardize = TRUE)
stack_save_clips(path_clips, path_plot, path_values, paste0(path_output, "/", combo_2$name),
                 combo_2$folders, full_IDs, normalize = FALSE, standardize = TRUE)
stack_save_clips(path_clips, path_plot, path_values, paste0(path_output, "/", combo_3$name),
                 combo_3$folders, full_IDs, normalize = FALSE, standardize = TRUE)
stack_save_clips(path_clips, path_plot, path_values, paste0(path_output, "/", combo_4$name),
                 combo_4$folders, full_IDs, normalize = FALSE, standardize = TRUE)

################################################################################