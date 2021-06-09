################################################################################
################################################################################
# PREPARE CNN INPUT
################################################################################
################################################################################

# load packages
library(sf)
library(raster)

# set paths
path_clips <- "D:/Masterarbeit_Zoe/4_Daten/clips"  # input
path_output <- "D:/Masterarbeit_Zoe/4_Daten/models/input"  # output
path_plot <- "D:/Masterarbeit_Zoe/4_Daten/vegetation/Export_ODK_clean_checked.kml" # input

# set folder names for each input type
# TODO: this is dummy data
tls_set <- c("nDSM", "reflectance", "point_density")
rgb_set <- c("ortho")
geo_set <- c("geometry")
all_set <- c(tls_set, rgb_set, geo_set)

# set names and folder names for each combination
combo_1 <- c()
combo_1$folders <- c(tls_set)
combo_1$name <- "tls"
#
combo_2 <- c()  # maybe this combination is unnecessary
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
  clip_list <- list.files(clip_dir, pattern=".tif", recursive=TRUE, full.names=TRUE)
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
# STACK & SAVE & SORT CLIPS
################################################################################

stack_save_clips <- function(clip_dir, plot_path, output_dir, selection_rasters, selection_IDs) {
  # stacks all rasters for each clip & saves it in a folder based on the label
  check_create_dir(output_dir)
  # get all clip paths of the selected raster types
  clip_list <- list.files(clip_dir, pattern=".tif", recursive=TRUE)
  contained <- c()
  for (clip_path in clip_list)  {
    contained <- c(contained, strsplit(clip_path, "/")[[1]][1] %in% selection_rasters)
  }
  clip_list <- clip_list[contained]
  # load vegetation plots & create one subfolder for each label
  plots <- st_read(plot_path)
  for (label in unique(plots$Name)) {
    check_create_dir(paste0(output_dir, "/", label))
  }
  # extract all plot IDs from vegetation plots
  all_IDs <- c()
  for(idx in 1:nrow(plots)) {
    all_IDs <- c(all_IDs, as.numeric(strsplit(plots$Description[idx], " ")[[1]][4]))
  }
  # loop through all selected plot IDs
  for (plot_ID in selection_IDs) {
    # get row from vegetation plots where this ID is right
    plot_idx <- which(all_IDs == plot_ID)
    # load & stack all clips of the plot
    plot_paths <- sort(clip_list[grepl(paste0("_", plot_ID, ".tif"), clip_list)])
    plot_rasters <- list()
    for (idx in 1:length(plot_paths)) {plot_rasters[[idx]] = stack(paste0(clip_dir, "/", plot_paths[idx]))}
    plot_stack <- stack(plot_rasters)
    # determine label & save in according folder
    plot_label <- plots$Name[plot_idx]
    writeRaster(plot_stack, paste0(output_dir, "/", plot_label, "/", plot_ID, ".tif"))
  }
}

################################################################################
# EXECUTION
################################################################################

check_create_dir(path_output)
full_IDs <- get_list_full_plots(path_clips, length(all_set))
stack_save_clips(path_clips, path_plot, paste0(path_output, "/", combo_1$name), combo_1$folders, full_IDs)
stack_save_clips(path_clips, path_plot, paste0(path_output, "/", combo_2$name), combo_2$folders, full_IDs)
stack_save_clips(path_clips, path_plot, paste0(path_output, "/", combo_3$name), combo_3$folders, full_IDs)
stack_save_clips(path_clips, path_plot, paste0(path_output, "/", combo_4$name), combo_4$folders, full_IDs)

################################################################################