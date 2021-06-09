################################################################################
################################################################################
# RASTER COLINEARITY CHECKS
################################################################################
################################################################################

# for this part, all rasters must be computed and normalized
# selecting the rasters must be done manually

# TODO: actually execute this for the right datasets!

# load packages
library(raster)  # for rasters
library(sp)  # for spatial objects
library(Hmisc)  # for cluster analysis

# path_rasters  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/rasters_copy"  # multiple areas
path_rasters  <- "D:/Masterarbeit_Zoe/4_Daten/rasters"  # single area

# get all relevant raster paths
raster_list <- list.files(path_rasters, pattern=".tif", recursive=TRUE, full.names=TRUE)
raster_list <- raster_list[!grepl("nDSM_filtering", raster_list)]  # reove unscales nDSM from list
raster_list <- raster_list[!grepl("ortho", raster_list)]  # I want to keep this anyway
raster_list <- raster_list[!grepl("temp", raster_list)]  # remove temporary files

# load all those rasters
rasters <- list()
for (idx in 1:length(raster_list)) {rasters[[idx]] = raster(raster_list[idx])}

# get unique area names
area_names <- c()
for (name in raster_list) {
  name <- strsplit(name, "[.]")[[1]][1]
  name <- strsplit(name, "_")
  name <- name[[1]][length(name[[1]])-1]
  area_names <- c(area_names, name)
}
unique_areas <- unique(area_names)

# stack rasters (of same area respectively)
raster_stacks <- c()
for (area in unique_areas) {
  single_stack <- stack(rasters[grepl(paste0("_", area, "_"), raster_list)])
  raster_stacks <- c(raster_stacks, single_stack)
  rm(single_stack)
}

# change the band names so they match & are short
band_names <- names(raster_stacks[[1]])  # assumes all have same rasters available
new_band_names <- c()
for (name in band_names) {
  name <- strsplit(name, "[.]")[[1]][2]
  name <- strsplit(name, "_")
  name <- name[[1]][1:2]
  name <- paste(name[1], name[2], sep="_")
  new_band_names <- c(new_band_names, name)
}
for (i in 1:length(raster_stacks)) {
  names(raster_stacks[[i]]) <- new_band_names
}

# prepare empty data frame for samples
raster_df <- matrix(ncol = length(new_band_names), nrow=0)
colnames(raster_df) <- new_band_names
raster_df <- data.frame(raster_df)

# draw samples from each plot & put them all into one data frame
for (i in 1:length(raster_stacks)) {
  raster_samples <- sampleRandom(raster_stacks[[i]], size=5000, cells=FALSE, sp=TRUE)
  raster_samples_df <- as.data.frame(raster_samples@data, xy=FALSE)
  raster_df <- rbind(raster_df, raster_samples_df)
  rm(raster_samples)
  rm(raster_samples_df)
}

################################################################################
# PCA 1
################################################################################

# PCA biplot:
# - gleiche Richtung = positiv korreliert
# - entgegengesetzt = negativ korreliert
# - 90 Grad = nicht korreliert

pca_all <- prcomp(raster_df)
summary(pca_all)  # with 5 PCs we are above the 90% explained variance
names(pca_all$sdev) <- as.character(1:20)
screeplot(pca_all, las=1, main="", cex.lab=1.5, xlab="Hauptkomponenten")  # to see which PCs explain how much variance
round(pca_all$rotation, 2)  # here we can see how important the variables is for each PC
biplot(pca_all, cex=c(0.5,1))#, xlim=c(-0.05, 0.05), ylim=c(-0.05, 0.05))

################################################################################
# CLUSTER ANALYSIS
################################################################################

# threshold: spearmans rho = 0.7
# beim rausnehmen drauf achten, was mehr zu den ersten PCA beiträgt?

# look at the situation
clust_all_1 <- varclus(as.matrix(raster_df), similarity = c("spearman"))
plot(clust_all_1)
abline(h=0.7^2, lty=2, col="deeppink3")

# handle reflectance
clust_all_2 <- varclus(as.matrix(raster_df[-c(17,19)]), similarity = c("spearman"))
plot(clust_all_2)
abline(h=0.7^2, lty=2, col="deeppink3")

# handle planarity & linearity
clust_all_3 <- varclus(as.matrix(raster_df[-c(17,19, 11,8,9)]), similarity = c("spearman"))
plot(clust_all_3)
abline(h=0.7^2, lty=2, col="deeppink3")

# handle curvature & shericity & ansiotrophy
clust_all_4 <- varclus(as.matrix(raster_df[-c(17,19, 11,8,9, 4,16,3,14,15,5,2)]), similarity = c("spearman"))
plot(clust_all_4)
abline(h=0.7^2, lty=2, col="deeppink3")

# checking if everything under threshold
clust_all_4[["sim"]]
clust_all_4[["sim"]]>(0.7^2)
cor(raster_df[-c(17,19, 11,8,9, 4,16,3,14,15,5,2)], method="spearman")

################################################################################
# PCA 2
################################################################################

pca_remains <- prcomp(raster_df[-c(17,19, 11,8,9, 4,16,3,14,15,5,2)])
summary(pca_remains)  # with 5 PCs we are above the 90% explained variance
names(pca_remains$sdev) <- as.character(1:20)
screeplot(pca_remains, las=1, main="", cex.lab=1.5, xlab="Hauptkomponenten")  # to see which PCs explain how much variance
round(pca_remains$rotation, 2)  # here we can see how important the variables is for each PC
biplot(pca_remains, cex=c(0.5,1))#, xlim=c(-0.05, 0.05), ylim=c(-0.05, 0.05))

################################################################################

# TODO: save the final rasters in combinations the CNN should be fed with?

################################################################################