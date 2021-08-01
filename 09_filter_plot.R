################################################################################
################################################################################
# CREATE FILTER PLOT
################################################################################
################################################################################

# load packages
library(lidR)
library(sf)
library(raster)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/03_raster_calculation_functions.R")

# set chunk parameters
chunk_size <- 15 # my RAM hates everything above, so I hate everything above
buffer_size <- 1 # to avoid edge effects & not having enough points for interpolation

# set paths
basedir <- "H:/Daten/Studium/2_Master/4_Semester"
path_points_01 <- paste0(basedir, "/4_Daten/data for plot/01_norm") # input
path_points_02 <- paste0(basedir, "/4_Daten/data for plot/02_no_stems") # output
path_points_03 <- paste0(basedir, "/4_Daten/data for plot/03_2m") # output
path_points_04 <- paste0(basedir, "/4_Daten/data for plot/04_filtered") # output
path_plots     <- paste0(basedir, "/5_Analyse/Plots") # output

################################################################################
# REMOVE STEMS
################################################################################

# load data
ctg_01 <- readTLSLAScatalog(path_points_01)

# set options
opt_chunk_buffer(ctg_01) <- buffer_size
opt_chunk_size(ctg_01) <- chunk_size
check_create_dir(path_points_02)
opt_output_files(ctg_01) <- paste0(path_points_02, "/{ID}")

################################################################################

filter_stem_ctg.LAScluster <- function(las) {
  las <- readLAS(las)
  if (is.empty(las)) {
    return(NULL)
  }
  map <- treeMap(las, map.hough())
  las <- treePoints(las, map, trp.crop())
  las <- stemPoints(las, stm.hough(pixel_size = 0.01))
  las <- filter_poi(las, Stem == FALSE)
  las <- filter_poi(las, buffer == 0)
  return(las)
}

filter_stem_ctg.LAScatalog <- function(las, height = 2, remove_stems = TRUE) {
  opt_select(las) <- "*"
  options <- list(
    need_output_file = TRUE,
    need_buffer = TRUE,
    automerge = TRUE
  )
  output <- catalog_apply(las, filter_stem_ctg.LAScluster, .options = options)
  return(output)
}

################################################################################

# execute
ctg_02 <- filter_stem_ctg.LAScatalog(ctg_01)
lidR:::catalog_laxindex(ctg_02)

################################################################################
# CUT AT 2M HEIGHT
################################################################################

# load data
ctg_02 <- readTLSLAScatalog(path_points_02)

# set options
opt_chunk_buffer(ctg_02) <- buffer_size
opt_chunk_size(ctg_02) <- chunk_size
check_create_dir(path_points_03)
opt_output_files(ctg_02) <- paste0(path_points_03, "/{ID}")

################################################################################

filter_2m_ctg.LAScluster <- function(las) {
  las <- readLAS(las)
  if (is.empty(las)) {
    return(NULL)
  }
  las <- filter_poi(las, Z <= 2)
  las <- filter_poi(las, buffer == 0)
  return(las)
}

filter_2m_ctg.LAScatalog <- function(las, height = 2, remove_stems = TRUE) {
  opt_select(las) <- "*"
  options <- list(
    need_output_file = TRUE,
    need_buffer = TRUE,
    automerge = TRUE
  )
  output <- catalog_apply(las, filter_2m_ctg.LAScluster, .options = options)
  return(output)
}

################################################################################

# execute
ctg_03 <- filter_2m_ctg.LAScatalog(ctg_02)
lidR:::catalog_laxindex(ctg_03)

################################################################################
# FILTER UNDERSTORY POINTS
################################################################################

# load data
ctg_03 <- readTLSLAScatalog(path_points_03)

# set options
opt_chunk_buffer(ctg_03) <- buffer_size
opt_chunk_size(ctg_03) <- chunk_size
check_create_dir(path_points_04)
opt_output_files(ctg_03) <- paste0(path_points_04, "/{ID}")

################################################################################

filter_understory_ctg.LAScluster <- function(las) {
  # returns point cloud without understory (LAS file)
  # load the data
  las <- readLAS(las)
  if (is.empty(las)) {
    return(NULL)
  }
  # voxel metrics: is there a point in the voxel or not?
  voxels <- voxel_metrics(las, length(X), res = 0.1, all_voxels = TRUE)
  voxels$V1 <- ifelse(voxels$V1 > 0, 1, 0)
  voxels$V1[is.na(voxels$V1)] <- 0
  # convert coordinates to [cm]
  # because otherwise R adds decimal places and makes this crash
  voxels$X <- as.integer(round(voxels$X * 100))
  voxels$Y <- as.integer(round(voxels$Y * 100))
  voxels$Z <- as.integer(round(voxels$Z * 100))
  # loop from lowest to highest z value, start at 0.5 m height
  # z_loop_vals <- sort(unique(voxels$Z))[9:length(unique(voxels$Z))]
  z_loop_vals <- seq(50, max(unique(voxels$Z)), by = 10)
  for (z_val in z_loop_vals) {
    # loop through every non-empty voxel with this z value
    z_loop_vox <- voxels[voxels$Z == z_val & voxels$V1 == 1, ]
    if (nrow(z_loop_vox != 0)) {
      for (idx in 1:nrow(z_loop_vox)) {
        vox <- z_loop_vox[idx, ]
        # set attribute to empty, if z-1, x+-1, y+-1 is all empty
        neighbours <- voxels$V1[voxels$Z == vox$Z - 10 &
          voxels$X >= vox$X - 10 & voxels$X <= vox$X + 10 &
          voxels$Y >= vox$Y - 10 & voxels$Y <= vox$Y + 10]
        if (mean(neighbours) == 0) {
          voxels$V1[voxels$X == vox$X & voxels$Y == vox$Y & voxels$Z == vox$Z] <- 0
        }
      }
    }
  }
  # convert coordinates back to [m]
  voxels$X <- voxels$X / 100
  voxels$Y <- voxels$Y / 100
  voxels$Z <- voxels$Z / 100
  # add voxel attributes to the points
  las <- add_lasattribute(las, 1, "V1", "keep voxels with 1") # create empty attribute
  # save points which should remain the same
  unchanged_las <- filter_poi(las, Z <= min(z_loop_vals - 5) / 100)
  # for each (filtered) vertical voxel layer, create a raster
  for (z_val in z_loop_vals) {
    # create raster
    z_subset <- voxels[as.integer(round(voxels$Z * 100)) == z_val, ]
    z_subset <- as.data.frame(z_subset)[, c(1, 2, 4)]
    new_raster <- rasterFromXYZ(z_subset)
    crs(new_raster) <- crs(las)
    # add raster values to point cloud
    las_z <- filter_poi(las, Z > ((z_val - 5) / 100) & Z <= ((z_val + 5) / 100))
    las_z <- merge_spatial(las_z, new_raster, "V1")
    # remove points with V1 == 0
    las_z <- filter_poi(las_z, V1 == 1)
    unchanged_las <- rbind(unchanged_las, las_z)
  }
  # delete all points with attribute empty
  las <- unchanged_las
  rm(unchanged_las)
  gc()
  # delete buffer & return points
  las <- filter_poi(las, buffer == 0)
  if (is.empty(las)) {
    return(NULL)
  }
  return(las)
}

filter_understory_ctg.LAScatalog <- function(las) {
  opt_select(las) <- "*"
  options <- list(
    need_output_file = TRUE,
    need_buffer = TRUE,
    automerge = TRUE
  )
  # execute & return
  output <- catalog_apply(las, filter_understory_ctg.LAScluster, .options = options)
  return(output)
}

################################################################################

# execute
ctg_04 <- filter_understory_ctg.LAScatalog(ctg_03)
lidR:::catalog_laxindex(ctg_04)

################################################################################
# EXTRACT DATA
################################################################################

y_split <- 5379698
y_thickness <- 0.75

paths <- list.files(path_points_01, pattern = "[.]las", full.names = T)
data_all <- c()
for (path in paths) {
  las <- readTLSLAS(path)
  las <- filter_poi(las, Y < y_split + y_thickness, Y > y_split - y_thickness)
  data <- las@data
  data <- data[, c("X", "Y", "Z")]
  data <- round(data / 0.05) * 0.05
  data <- data[!duplicated(data[, c("X", "Z")]), ]
  data_all <- rbind(data_all, data)
}
write.csv(data_all, paste0(dirname(path_points_01), "/01_norm_", y_thickness * 200, "cm.csv"), row.names = F)

################################################################################

paths <- list.files(path_points_02, pattern = "[.]las", full.names = T)
data_all <- c()
for (path in paths) {
  las <- readTLSLAS(path)
  las <- filter_poi(las, Y < y_split + y_thickness, Y > y_split - y_thickness)
  data <- las@data
  data <- data[, c("X", "Y", "Z")]
  data <- round(data / 0.05) * 0.05
  data <- data[!duplicated(data[, c("X", "Z")]), ]
  data_all <- rbind(data_all, data)
}
write.csv(data_all, paste0(dirname(path_points_01), "/02_no_stems_", y_thickness * 200, "cm.csv"), row.names = F)

################################################################################

paths <- list.files(path_points_03, pattern = "[.]las", full.names = T)
data_all <- c()
for (path in paths) {
  las <- readTLSLAS(path)
  las <- filter_poi(las, Y < y_split + y_thickness, Y > y_split - y_thickness)
  data <- las@data
  data <- data[, c("X", "Y", "Z")]
  data <- round(data / 0.05) * 0.05
  data <- data[!duplicated(data[, c("X", "Z")]), ]
  data_all <- rbind(data_all, data)
}
write.csv(data_all, paste0(dirname(path_points_01), "/03_2m_", y_thickness * 200, "cm.csv"), row.names = F)

################################################################################

paths <- list.files(path_points_04, pattern = "[.]las", full.names = T)
data_all <- c()
for (path in paths) {
  las <- readTLSLAS(path)
  las <- filter_poi(las, Y < y_split + y_thickness, Y > y_split - y_thickness)
  data <- las@data
  data <- data[, c("X", "Y", "Z")]
  data <- round(data / 0.05) * 0.05
  data <- data[!duplicated(data[, c("X", "Z")]), ]
  data_all <- rbind(data_all, data)
}
write.csv(data_all, paste0(dirname(path_points_01), "/04_filtered_", y_thickness * 200, "cm.csv"), row.names = F)

################################################################################
# MAKE PLOTS
################################################################################

# load
csv_norm     <- read.csv(paste0(dirname(path_points_01), "/01_norm_200cm.csv"))
csv_no_stems <- read.csv(paste0(dirname(path_points_01), "/02_no_stems_200cm.csv"))
csv_2m       <- read.csv(paste0(dirname(path_points_01), "/03_2m_200cm.csv"))
csv_filtered <- read.csv(paste0(dirname(path_points_01), "/04_filtered_200cm.csv"))

# merge
csv_norm$source <- "Original Data"
csv_no_stems$source <- "+ without Stems"
csv_2m$source <- "+ below 2m"
csv_filtered$source <- "+ filtered"
csv_all <- Reduce(bind_rows, list(csv_norm, csv_no_stems, csv_2m, csv_filtered))
csv_all$source <- factor(csv_all$source, levels = c("Original Data", "+ without Stems", "+ below 2m", "+ filtered"))

# set colors
colors <- c(
  "Original Data" = "grey85", "+ without Stems" = "grey65",
  "+ below 2m" = "grey45", "+ filtered" = "grey25"
)

# filter
csv_all <- csv_all[csv_all$X < 446450, ]

# plot big
big <- ggplot() +
  geom_point(data = csv_all, aes(x = X, y = Z, col = source), size = 0.1) +
  geom_hline(yintercept = 2, linetype = "dashed", size = 0.25, color = "grey10") +
  geom_rect(aes(xmin = 446400 - 0.25, xmax = 446419 + 0.25, ymin = 0 - 0.25, ymax = 2 + 0.25),
    fill = "white", alpha = 0, color = "grey10", linetype = "solid", size = 0.25
  ) +
  # ylab("Height (m)") + xlab("x-Coordinate (WGS 84 / UTM Zone 32N)") +
  ylab("Height (m)\n") +
  xlab("") +
  scale_y_continuous(lim = c(-0.25, 26)) +
  coord_fixed() +
  theme_light() +
  scale_color_manual(values = colors) +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.text = element_text(family = "Calibri", size = 14),
    legend.title = element_text(family = "Calibri", size = 16),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.size = unit(0.75, "cm"),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5), title = "Processing Step"))

# plot small
small <- ggplot() +
  geom_point(data = csv_all[csv_all$Z <= 2.25 & csv_all$X < 446419 + 0.25, ], aes(x = X, y = Z, col = source), size = 0.25) +
  geom_hline(yintercept = 2, linetype = "dashed", size = 0.25, color = "grey10") +
  geom_rect(aes(xmin = 446400 - 0.25, xmax = 446419 + 0.25, ymin = 0 - 0.25, ymax = 2 + 0.25),
    fill = "white", alpha = 0, color = "grey10", linetype = "solid", size = 0.25
  ) +
  # ylab("Height (m)") + xlab("x-Coordinate (WGS 84 / UTM Zone 32N)") +
  ylab("Height (m)\n") +
  xlab("") +
  scale_y_continuous(breaks = 0:2, lim = c(-0.25, 2.25)) +
  coord_fixed() +
  theme_light() +
  scale_color_manual(values = colors) +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.text = element_text(family = "Calibri", size = 14),
    legend.title = element_text(family = "Calibri", size = 16),
    legend.spacing.y = unit(0.25, "cm"),
    legend.key.size = unit(0.75, "cm"),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5), title = "Processing Step"))

# legend
legend <- get_legend(small)

# combine & save
cairo_pdf(file = paste0(path_plots, "/filtering_both.pdf"), family = "Calibri", width = 8.27, height = 5.83)
ggarrange(big, small,
  heights = c(3, 1), ncol = 1, nrow = 2,
  legend.grob = legend, legend = "bottom", align = "v"
)
dev.off()

################################################################################