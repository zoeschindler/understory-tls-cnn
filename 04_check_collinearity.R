################################################################################
################################################################################
# RASTER COLINEARITY CHECKS
################################################################################
################################################################################

# load packages
library(raster)
library(sp)
library(Hmisc)

# set path
path_rasters <- "C:/Users/Zoe/Documents/understory_classification/4_Daten/rasters_2cm"

################################################################################
# GETTING RASTER SAMPLES
################################################################################

# get all relevant raster paths
raster_list <- list.files(path_rasters, pattern = "[.]tif", recursive = TRUE, full.names = TRUE)
raster_list <- raster_list[!grepl("nDSM_filtering", raster_list)] # remove unscaled nDSM from list
raster_list <- raster_list[!grepl("DTM", raster_list)] # remove DTM from list
raster_list <- raster_list[!grepl("ortho", raster_list)] # I want to keep this anyway
raster_list <- raster_list[!grepl("temp", raster_list)] # remove temporary files

# load all those rasters
rasters <- list()
for (idx in 1:length(raster_list)) {
  rasters[[idx]] <- raster(raster_list[idx])
}

# get unique area names
area_names <- c()
for (name in raster_list) {
  name <- strsplit(name, "[.]")[[1]][1]
  name <- strsplit(name, "_")
  name <- name[[1]][length(name[[1]]) - 1]
  area_names <- c(area_names, name)
}
unique_areas <- as.numeric(unique(area_names))

# stack rasters (of same area respectively)
raster_stacks <- c()
for (area in unique_areas) {
  single_stack <- stack(rasters[grepl(paste0("_", area, "_"), raster_list)])
  raster_stacks <- c(raster_stacks, single_stack)
  rm(single_stack)
}
rm(rasters)

# change the band names so they are the same for ever area & are short
band_names <- names(raster_stacks[[1]]) # assumes all have same rasters available
new_band_names <- c()
for (name in band_names) {
  name <- strsplit(name, "_")
  name <- name[[1]][1:2]
  if ("area" %in% name) {
    name <- name[1]
  }
  else {
    name <- paste(name[1], name[2], sep = "_")
  }
  new_band_names <- c(new_band_names, name)
}
for (i in 1:length(raster_stacks)) {
  names(raster_stacks[[i]]) <- new_band_names
}

# prepare empty data frame for samples
raster_df <- matrix(ncol = length(new_band_names), nrow = 0)
colnames(raster_df) <- new_band_names
raster_df <- data.frame(raster_df)

# draw samples from each plot & put them all into one data frame
for (i in 1:length(raster_stacks)) {
  print(paste0("... getting samples from area ", i))
  raster_samples <- sampleRandom(raster_stacks[[i]], size = 1000, cells = FALSE, sp = TRUE)
  raster_samples <- as.data.frame(raster_samples@data, xy = FALSE)
  raster_df <- rbind(raster_df, raster_samples)
  rm(raster_samples)
}

# save to avoid waiting again
write.csv(raster_df, paste0(path_rasters, "/raster_samples_unscaled.csv"), row.names = FALSE)

################################################################################

# read samples
raster_df <- read.csv(paste0(path_rasters, "/raster_samples_unscaled.csv"))

# rescale sample values
for (i in 1:ncol(raster_df)) {
  raster_df[, i] <- scale(raster_df[, i])
}

# save to avoid waiting again
write.csv(raster_df, paste0(path_rasters, "/raster_samples_scaled.csv"), row.names = FALSE)

################################################################################

# read samples
raster_df <- read.csv(paste0(path_rasters, "/raster_samples_scaled.csv"))

################################################################################
# PCA 1
################################################################################

# PCA biplot:
# - gleiche Richtung = positiv korreliert
# - entgegengesetzt = negativ korreliert
# - 90 Grad = nicht korreliert

pca_all <- prcomp(raster_df)
summary(pca_all) # with 2 PCs we are above the 90% explained variance
names(pca_all$sdev) <- as.character(1:20)
screeplot(pca_all, las = 1, main = "", cex.lab = 1.5, xlab = "Hauptkomponenten") # to see which PCs explain how much variance
round(pca_all$rotation, 3) # here we can see how important the variables is for each PC
biplot(pca_all, cex = c(0.5, 1)) # , xlim=c(-0.05, 0.05), ylim=c(-0.05, 0.05))

################################################################################
# CLUSTER ANALYSIS
################################################################################

# threshold: spearmans rho = 0.7
# base removal on the contributions of the variables to the first PCA axes

# look at the situation
clust_all_1 <- varclus(as.matrix(raster_df), similarity = c("spearman"))
plot(clust_all_1)
abline(h = 0.7^2, lty = 2, col = "deeppink3")

# clusters:
# - reflectance_max [18] / reflectance_mean [19]
# - curvature_sd [6] / anisotropy_sd [3] / sphericity_sd [15]
# - linearity_max [7] / linearity_mean [8]
# - planarity_max [10] / planarity_mean [11]
# - curvature_mean [5] / anisotropy_mean [2] / sphericity_mean [14] / curvature_max [4] / sphericity_max [13]

# cluster: reflectance_max [18] / reflectance_mean [19]
# keep: reflectance_mean [19]
# remove: reflectance_max [18]
clust_all_2 <- varclus(as.matrix(raster_df[-c(18)]), similarity = c("spearman"))
plot(clust_all_2)
abline(h = 0.7^2, lty = 2, col = "deeppink3")

# cluster: curvature_sd [6] / anisotropy_sd [3] / sphericity_sd [15] (all really similar)
# keep: sphericity_sd [15]
# remove: anisotropy_sd [3] / curvature_sd [6]
clust_all_3 <- varclus(as.matrix(raster_df[-c(18, 3, 6)]), similarity = c("spearman"))
plot(clust_all_3)
abline(h = 0.7^2, lty = 2, col = "deeppink3")

# cluster: linearity_max [7] / linearity_mean [8]
# keep: linearity_max [7]
# remove: linearity_mean [8]
clust_all_4 <- varclus(as.matrix(raster_df[-c(18, 3, 6, 8)]), similarity = c("spearman"))
plot(clust_all_4)
abline(h = 0.7^2, lty = 2, col = "deeppink3")

# cluster: planarity_max [10] / planarity_mean [11]
# keep: planarity_mean [11]
# remove: planarity_max [10]
clust_all_5 <- varclus(as.matrix(raster_df[-c(18, 3, 6, 8, 10)]), similarity = c("spearman"))
plot(clust_all_5)
abline(h = 0.7^2, lty = 2, col = "deeppink3")

# cluster: curvature_mean [5] / anisotropy_mean [2] / sphericity_mean [14] / sphericity_sd [15] / curvature_max [4] / sphericity_max [13]
# keep: curvature_max [4]
# remove: curvature_mean [5] / anisotropy_mean [2] / sphericity_mean [14] / sphericity_sd [15] / sphericity_max [13]
clust_all_6 <- varclus(as.matrix(raster_df[-c(18, 3, 6, 8, 10, 5, 2, 14, 15, 13)]), similarity = c("spearman"))
plot(clust_all_6)
abline(h = 0.7^2, lty = 2, col = "deeppink3")

# checking if everything under threshold
clust_all_6[["sim"]]
clust_all_6[["sim"]] > (0.7^2)
cor(raster_df[-c(18, 3, 6, 8, 10, 5, 2, 14, 15, 13)], method = "spearman")
cor(raster_df[-c(18, 3, 6, 8, 10, 5, 2, 14, 15, 13)], method = "spearman") < 0.7

# save non-collinear data
write.csv(raster_df[-c(18, 3, 6, 8, 10, 5, 2, 14, 15, 13)],
  paste0(path_rasters, "/raster_samples_scaled_noncollinear.csv"),
  row.names = FALSE
)

################################################################################
# PCA 2
################################################################################

pca_remains <- prcomp(raster_df[-c(18, 3, 6, 8, 10, 5, 2, 14, 15, 13)])
summary(pca_remains) # with 5 PCs we are above the 90% explained variance
names(pca_remains$sdev) <- as.character(1:10)
screeplot(pca_remains, las = 1, main = "", cex.lab = 1.5, xlab = "Hauptkomponenten") # to see which PCs explain how much variance
round(pca_remains$rotation, 3) # here we can see how important the variables is for each PC
biplot(pca_remains, cex = c(0.5, 1)) # , xlim=c(-0.05, 0.05), ylim=c(-0.05, 0.05))

################################################################################
# RESULTS
################################################################################

# keep:
keep_names <- sort(names(raster_df[-c(18, 3, 6, 8, 10, 5, 2, 14, 15, 13)]))
print(keep_names)
# "ortho", "anisotropy_max", "curvature_max", "linearity_max",
# "linearity_sd", "nDSM", "planarity_mean", "planarity_sd",
# "point_density", "reflectance_mean", "reflectance_sd")

# remove:
remove_names <- sort(names(raster_df[c(18, 3, 6, 8, 10, 5, 2, 14, 15, 13)]))
print(remove_names)
# "anisotropy_mean", "anisotropy_sd", "curvature_mean", "curvature_sd",
# "linearity_mean", "planarity_max", "reflectance_max", "sphericity_max",
# "sphericity_mean", "sphericity_sd"

################################################################################