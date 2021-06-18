################################################################################
################################################################################
# SETUP CNN
################################################################################
################################################################################

# load packages
library(abind)
library(BBmisc)
library(keras)
library(raster)

# set paths
path_clips <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/models/input_unfiltered/tls_rgb_geo"  # input
path_output <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/models/out"  # output

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
# READING IN IMAGES
################################################################################

tif_to_rds <- function(clip_dir, output_dir, pixels, bands, seed=123) {
  # reads in images, divide into stratified k-folds, save as rds
  check_create_dir(output_dir)
  # set seed
  set.seed(seed)
  # get all file paths, excluding smallest groups
  img_paths_all <- list.files(clip_dir, pattern=".tif", full.names=TRUE)
  img_paths_all <- img_paths_all[!grepl("grass", img_paths_all)]
  img_paths_all <- img_paths_all[!grepl("rock", img_paths_all)]
  # get all labels
  img_labels_all <- sapply(strsplit(basename(img_paths_all), "_"), "[[", 1)
  img_labels_all <- as.numeric(as.factor(img_labels_all))
  # set up empty list for indices
  fold_indices <- list()
  for (i in 1:5) {
    fold_indices[[i]] <- NA
  }
  # loop through labels, split each evenly, save grouped indices
  for (grp in unique(img_labels_all)) {
    grp_indices <- which(img_labels_all == grp)
    grp_indices <- grp_indices[sample(length(grp_indices))]
    grp_indices_split <- chunk(grp_indices, n.chunks = 5)
    # save indices
    for (j in 1:5) {
      fold_indices[[j]] <- c(fold_indices[[j]], grp_indices_split[[j]])
    }
  }
  # remove dummy NA values
  for (k in 1:5) {fold_indices[[k]] <- as.numeric(na.omit(fold_indices[[k]]))}
  # loop through all folds
  output_paths <- c()
  for (l in 1:5) {
    print(paste0("... loading & saving split ", l))
    # get all filenames & labels of the fold
    fold_files <- img_paths_all[fold_indices[[l]]]
    fold_labels <- img_labels_all[fold_indices[[l]]]
    # load all files of the fold
    img_array <- array(dim=c(length(fold_files), pixels, pixels, bands))
    label_vector <- c()
    for (i in 1:length(fold_files)) {
      path <- fold_files[i]
      label_vector <- c(label_vector, fold_labels[i])
      raster <- raster::values(stack(path))
      img_array[i,,,] <- raster
    }
    # save image array as rds
    data_list <- list(img=img_array, label=label_vector)
    saveRDS(data_list, file = paste0(output_dir, "/images_fold_", l, ".rds"))
    output_paths <- c(output_paths, paste0(output_dir, "/images_fold_", l, ".rds"))
  }
  return(output_paths)
}

################################################################################
# DATA PREPARATION & AUGMENTATION
################################################################################

create_dataset <- function(rdata_list, holdout_fold, balance_classes=TRUE) {
  # create dataset for CNN 
  # load test data
  raw_test <- readRDS(rdata_list[holdout_fold])
  img_test <- raw_test$img
  label_test <- raw_test$label
  # load training folds
  img_train <- array(dim=c(0,50,50,13))
  label_train <- c()
  for (i in 1:length(rdata_list)) {
    if (i != holdout_fold) {
      raw_train <- readRDS(rdata_list[i])
      img_train <- unname(abind(img_train, raw_train$img, along=1))
      label_train <- c(label_train, raw_train$label)
      rm(raw_train)
    }
  }
  # replace NA with -1 in images
  img_train[is.na(img_train)] <- -1
  img_test[is.na(img_test)] <- -1
  # TODO: add steps for validation dataset
  # shuffle data, training data
  new_idx_test <- sample(length(label_test))
  label_test <- label_test[new_idx_test]
  img_test <- img_test[new_idx_test,,,]
  # shuffle data, training data
  new_idx_train <- sample(length(label_train))
  label_train <- label_train[new_idx_train]
  img_train <- img_train[new_idx_train,,,]
  # duplicate images in training set depending on label frequency
  if (balance_classes) {
    repeated_indices <- c()
    max_per_image <- 5  # maximum augmentation per image
    max_length <- 500  # maximum sample size per class
    for (label in unique(label_train)) {
      label_idx <- which(label_train == label)
      new_length <- ifelse(max_per_image * length(label_idx) > max_length,
                           max_length,
                           max_per_image * length(label_idx))
      label_idx <- rep(label_idx, length.out=new_length)
      repeated_indices <- c(repeated_indices, label_idx)
    }
    label_train <- label_train[repeated_indices]
    img_train <- img_train[repeated_indices,,,]
    # shuffle data, training data
    new_idx_train <- sample(length(label_train))
    label_train <- label_train[new_idx_train]
    img_train <- img_train[new_idx_train,,,]
  }
  # empty data generator, for test / validation data
  no_augmentation <- image_data_generator()
  # data generator with data augmentation, for training data
  do_augmentation <- image_data_generator(
    fill_mode = "reflect",
    rotation_range = 10,
    # width_shift_range = 0.02,
    # height_shift_range = 0.02,
    horizontal_flip = TRUE,
    vertical_flip = TRUE)
  # get number of classes
  label_classes <- length(unique(label_train))
  # generate batches, for test data
  test_flow <- flow_images_from_data(
    x = img_test,
    y = to_categorical(label_test)[,2:(label_classes+1)],
    generator = no_augmentation)
  # generate batches, for training data
  train_flow <- flow_images_from_data(
    x = img_train,
    y = to_categorical(label_train)[,2:(label_classes+1)],
    generator = do_augmentation)
  # TODO: add steps for validation dataset
  # return ready to use image_generators & steps per epoch
  return(list(data_test = test_flow,
              data_train = train_flow,
              steps_test = floor(length(label_test)/32),
              steps_train = floor(length(label_train)/32)))
}

################################################################################
# EXECUTION
################################################################################

# execute only once because it takes time
tif_to_rds(path_clips, path_output, width_length, n_bands, seed=123)

# later execute in a loop k times, also train & assess model k times
rdata_paths <- list.files(path_output, pattern="[.]rds", full.names=TRUE)
imbalanced <- create_dataset(rdata_paths, 4, balance_classes=FALSE)
balanced <- create_dataset(rdata_paths, 4, balance_classes=TRUE)

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
    
    # block 1
    layer_conv_2d(filters = 32 * filter_factor, kernel_size = c(5,5), activation = "relu", padding="same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 2
    layer_conv_2d(filters = 48 * filter_factor, kernel_size = c(5,5), activation = "relu", padding="valid") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 3
    layer_flatten() %>%
    layer_dense(units = 256 * filter_factor, activation = "relu") %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 84 * filter_factor, activation = "relu") %>%
    layer_dropout(0.5) %>%
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
    layer_conv_2d(filters = 384 * filter_factor, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
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
    layer_conv_2d(filters = 64 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>% 
    # block 2
    layer_conv_2d(filters = 128 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 128 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 3
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 4
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 5
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 6
    layer_flatten() %>%
    layer_dense(units = 4096 * filter_factor, activation = "relu") %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 4096 * filter_factor, activation = "relu") %>%
    layer_dropout(0.5) %>%
    layer_dense(units = n_classes, activation = "softmax")
  return(model)
}

################################################################################  