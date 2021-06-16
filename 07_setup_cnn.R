################################################################################
################################################################################
# SETUP CNN
################################################################################
################################################################################

# load packages
library(keras)
library(BBmisc)
library(raster)

# set paths
path_clips <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/models/input_unfiltered/tls_rgb_geo"
path_output <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/models/out"

# set data input parameters
width_length <- 50  # number of pixels
input_type <- basename(path_clips)
if (input_type == "tls") {
  n_bands <- 4
} else if (input_type == "tls_geo") {
  n_bands <- 10
} else if (input_type == "tls_rgb") {
  n_bands <- 7
} else if (input_type == "tls_rgb_geo") {
  n_bands <- 13
}

# set architecture parameters
n_band_selector <- 5  # number of layers after band selector
n_classes <- 7  # number of understory classes
filter_factor <- 0.5  # multiplier for amount of filters per layer (0.25 / 0.5 / 1 / 2)

################################################################################
# READING IN IMAGES
################################################################################

# set seed
set.seed(666)

# get all file paths & labels of clips
img_paths_all <- list.files(path_clips, pattern=".tif", full.names=TRUE)
img_labels_all <- sapply(strsplit(basename(img_paths_all), "_"), "[[", 1)
img_labels_all <- as.numeric(as.factor(img_labels_all))

# set up empty list for indices
fold_indices <- list()
for (i in 1:5) {fold_indices[[i]] <- NA}

# loop through labels
for (grp in unique(img_labels_all)) {
  # group by label
  grp_indices <- which(img_labels_all == grp)
  # shuffle them
  grp_indices <- grp_indices[sample(length(grp_indices))]
  # split into 5 folds
  grp_indices_split <- chunk(grp_indices, n.chunks = 5)
  # save indices
  for (j in 1:5) {
    fold_indices[[j]] <- c(fold_indices[[j]], grp_indices_split[[j]])
  }
}

# remove dummy NA values
for (k in 1:5) {fold_indices[[k]] <- as.numeric(na.omit(fold_indices[[k]]))}

# loop through all folds
for (l in 1:5) {
  print(paste0("loading & saving split ", l))
  # get all filenames & labels of the fold
  fold_files <- img_paths_all[fold_indices[[l]]]
  fold_labels <- img_labels_all[fold_indices[[l]]]
  # load all files of the fold
  img_array <- array(dim=c(length(fold_files), width_length, width_length, n_bands))
  label_vector <- c()
  for (i in 1:length(fold_files)) {
    path <- fold_files[i]
    label_vector <- c(label_vector, fold_labels[i])
    raster <- raster::values(stack(path))
    img_array[i,,,] <- raster
  }
  # save image array as RData
  print(summary(as.factor(fold_labels)))
  save(img_array, file = paste0(path_output, "/images_fold_", l, ".RData"))
  print("---")
}

################################################################################
# DATA AUGMENTATION
################################################################################

# remove seed
rm(.Random.seed, envir=globalenv())

random_augmentation <- function(img) {
  
}

create_dataset <- function(img_path, holdout_fold) {
  
}

################################################################################
# LeNet-5
# http://yann.lecun.com/exdb/publis/pdf/lecun-98.pdf
# https://www.kaggle.com/curiousprogrammer/lenet-5-cnn-with-keras-99-48 (based on this code)
################################################################################

get_lenet5 <- function(width_length, n_bands, n_band_selector, n_classes, filter_factor) {
  model <- keras_model_sequential() %>%
    # band selector
    layer_conv_2d(input_shape = c(width_length, width_length, n_bands), filters = n_band_selector,
                  kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
    # main model
    layer_conv_2d(filters = 32 * filter_factor, kernel_size = c(5,5), activation = "relu", padding="same") %>%
    # normalization?
    layer_max_pooling_2d() %>%
    layer_conv_2d(filters = 48 * filter_factor, kernel_size = c(5,5), activation = "relu", padding="valid") %>%
    # normalization?
    layer_max_pooling_2d() %>%
    layer_flatten() %>%
    layer_dense(units = 256 * filter_factor, activation = "relu") %>%
    # dropout?
    layer_dense(units = 84 * filter_factor, activation = "relu") %>%
    # droput?
    layer_dense(units = n_classes, activation = "softmax")
  return(model)
}

  
################################################################################
# AlexNet
# https://papers.nips.cc/paper/2012/file/c399862d3b9d6b76c8436e924a68c45b-Paper.pdf
# https://towardsdatascience.com/implementing-alexnet-cnn-architecture-using-tensorflow-2-0-and-keras-2113e090ad98 (based on this code)
# https://github.com/r-tensorflow/alexnet/blob/master/R/alexnet_train.R (and this code)
################################################################################

# hat eigentlich viel größeren Input (224x224x3) & mehr Klassen (1000)
# im paper ist max_pooling mit pool_size = 3 und stride = 2, in dem github code nicht
# je nach code quelle ist padding = "same" / "valid"

get_alexnet <- function(width_length, n_bands, n_band_selector, n_classes, filter_factor) {
  model <- keras_model_sequential() %>%
    # band selector
    layer_conv_2d(input_shape = c(width_length, width_length, n_bands), filters = n_band_selector,
                  kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
    # main model
    # block 1
    layer_conv_2d(filters = 96 * filter_factor, kernel_size = c(11,11), activation = "relu", strides = c(4,4)) %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
    # block 2
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(5,5), activation = "relu", strides = c(1,1), padding="same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
    # block 4
    layer_conv_2d(filters = 384 * filter_factor, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
    layer_batch_normalization() %>%
    layer_conv_2d(filters = 384 * filter_factor, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
    layer_batch_normalization() %>%
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
    # block 5
    layer_flatten() %>%
    layer_dense(units = 4096 * filter_factor, activation = "relu") %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 4096 * filter_factor, activation = "relu") %>%
    layer_dropout(0.5) %>%
    layer_dense(units = n_classes, activation = "softmax")
  return(model)
}

################################################################################
# VGG-16
# https://arxiv.org/pdf/1409.1556.pdf
# https://towardsdatascience.com/step-by-step-vgg16-implementation-in-keras-for-beginners-a833c686ae6c (based on this code)
################################################################################

get_alexnet <- function(width_length, n_bands, n_band_selector, n_classes, filter_factor) {
  model <- keras_model_sequential() %>%
    # band selector
    layer_conv_2d(input_shape = c(width_length, width_length, n_bands), filters = n_band_selector,
                  kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
    # main model
    # block 1
    layer_conv_2d(filters = 64 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 64 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_max_pooling_2d() %>% 
    # block 2
    layer_conv_2d(filters = 128 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 128 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_max_pooling_2d() %>%
    # block 3
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_max_pooling_2d() %>%
    # block 4
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_max_pooling_2d() %>%
    # block 5
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
      # normalization?
    layer_max_pooling_2d() %>%
    # block 6
    layer_flatten() %>%
    layer_dense(units = 4096 * filter_factor, activation = "relu") %>%
      # dropout?
    layer_dense(units = 4096 * filter_factor, activation = "relu") %>%
      # dropout?
    layer_dense(units = n_classes, activation = "softmax")
  return(model)
}

################################################################################  