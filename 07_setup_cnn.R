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

strat_folds <- function(label_vector, folds) {
  # splits label vector evenly into stratified folds, returns list with indices
  # set up empty list for indices
  fold_indices <- list()
  for (i in 1:folds) {
    fold_indices[[i]] <- NA
  }
  # loop through labels, split each evenly, save grouped indices
  for (grp in unique(label_vector)) {
    grp_indices <- which(label_vector == grp)
    grp_indices <- grp_indices[sample(length(grp_indices))]
    grp_indices_split <- chunk(grp_indices, n.chunks = folds)
    # save indices
    for (j in 1:folds) {
      fold_indices[[j]] <- c(fold_indices[[j]], grp_indices_split[[j]])
    }
  }
  # remove dummy NA values
  for (k in 1:folds) {
    fold_indices[[k]] <- as.numeric(na.omit(fold_indices[[k]]))
  }
  # retutn list of indices, one list entry of indices per fold
  return(fold_indices)
}

################################################################################

balance_by_duplicates <- function(label_vector, image_array, max_per_image=5, max_length=300) {
  # duplicates labels & images per label until:
  # - each images was duplicated max_per_image times
  # - amount of samples per label was raised to max_length respectively
  repeated_indices <- c()
  for (label in unique(label_vector)) {
    label_indices <- which(label_vector == label)
    new_length <- ifelse(max_per_image * length(label_indices) > max_length,
                         max_length,
                         max_per_image * length(label_indices))
    label_indices <- rep(label_indices, length.out=new_length)
    repeated_indices <- c(repeated_indices, label_indices)
  }
  label_vector <- label_vector[repeated_indices]
  image_array <- image_array[repeated_indices,,,]
  # shuffle data
  new_idx_train <- sample(length(label_vector))
  label_vector <- label_vector[new_idx_train]
  image_array <- image_array[new_idx_train,,,]
  # return new images & labels in list
  return(list(label=label_vector, img=image_array))
}

################################################################################
# READING IN IMAGES
################################################################################

tif_to_rds <- function(clip_dir, pixels, bands, folds=5, seed=123) {
  # reads in images, divide into stratified k-folds, save as rds
  # set seed
  set.seed(seed)
  # get all file paths, excluding smallest groups
  img_paths_all <- list.files(clip_dir, pattern=".tif", full.names=TRUE)
  img_paths_all <- img_paths_all[!grepl("grass", img_paths_all)]
  img_paths_all <- img_paths_all[!grepl("rock", img_paths_all)]
  # get all labels
  img_labels_all <- sapply(strsplit(basename(img_paths_all), "_"), "[[", 1)
  img_labels_all <- as.numeric(as.factor(img_labels_all))
  # make stratified split
  fold_indices <- strat_folds(img_labels_all, folds)
  # loop through all folds
  output_paths <- c()
  for (l in 1:folds) {
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
    saveRDS(data_list, file = paste0(clip_dir, "/images_fold_", l, ".rds"))
    output_paths <- c(output_paths, paste0(clip_dir, "/images_fold_", l, ".rds"))
  }
  return(output_paths)
}

################################################################################
# DATA PREPARATION & AUGMENTATION
################################################################################

create_dataset <- function(rdata_list, holdout_fold, pixels, bands, balance_classes=TRUE) {
  # create dataset for CNN 
  # load test data
  raw_test <- readRDS(rdata_list[holdout_fold])
  img_test <- raw_test$img
  label_test <- raw_test$label
  rm(raw_test)
  # load training folds
  img_train_vali <- array(dim=c(0, pixels, pixels, bands))
  label_train_vali <- c()
  for (i in 1:length(rdata_list)) {
    if (i != holdout_fold) {
      raw_train_vali <- readRDS(rdata_list[i])
      img_train_vali <- unname(abind(img_train_vali, raw_train_vali$img, along=1))
      label_train_vali <- c(label_train_vali, raw_train_vali$label)
      rm(raw_train_vali)
    }
  }
  # replace NA with -1 in images
  img_test[is.na(img_test)] <- -1
  img_train_vali[is.na(img_train_vali)] <- -1
  # make stratified split
  fold_indices <- strat_folds(label_train_vali, 5)  # 5 folds -> use 20% for validation
  indices_train <- unlist(fold_indices[1:4])
  indices_vali <- fold_indices[[5]]
  img_train <- img_train_vali[indices_train,,,]
  label_train <- label_train_vali[indices_train]
  img_vali <- img_train_vali[indices_vali,,,]
  label_vali <- label_train_vali[indices_vali]
  # shuffle datasets
  new_idx_test <- sample(length(label_test))
  label_test <- label_test[new_idx_test]
  img_test <- img_test[new_idx_test,,,]
  ###
  new_idx_train_vali <- sample(length(label_train_vali))
  label_train_vali <- label_train_vali[new_idx_train_vali]
  img_train_vali <- img_train_vali[new_idx_train_vali,,,]
  ###
  new_idx_train <- sample(length(label_train))
  label_train <- label_train[new_idx_train]
  img_train <- img_train[new_idx_train,,,]
  ###
  new_idx_vali <- sample(length(label_vali))
  label_vali <- label_vali[new_idx_vali]
  img_vali <- img_vali[new_idx_vali,,,]
  ###
  if(balance_classes) {
    # duplicate images depending on label frequency
    duplicated_train_vali <- balance_by_duplicates(label_train_vali, img_train_vali)
    label_train_vali <- duplicated_train_vali$label
    img_train_vali <- duplicated_train_vali$img
    ###
    duplicated_train <- balance_by_duplicates(label_train, img_train)
    label_train <- duplicated_train$label
    img_train <- duplicated_train$img
    # shuffle datasets
    new_idx_train_vali <- sample(length(label_train_vali))
    label_train_vali <- label_train_vali[new_idx_train_vali]
    img_train_vali <- img_train_vali[new_idx_train_vali,,,]
    ###
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
    #width_shift_range = 0.02,
    #height_shift_range = 0.02,
    #shear_range = 10,
    horizontal_flip = TRUE,
    vertical_flip = TRUE)
  # get number of classes
  label_classes <- length(unique(label_train))
  # generate batches, for test data
  flow_test <- flow_images_from_data(
    x = img_test,
    y = to_categorical(label_test)[,2:(label_classes+1)],
    generator = no_augmentation)
  # generate batches, for training + validation data
  flow_train_vali <- flow_images_from_data(
    x = img_train_vali,
    y = to_categorical(label_train_vali)[,2:(label_classes+1)],
    generator = do_augmentation)
  # generate batches, for training data
  flow_train <- flow_images_from_data(
    x = img_train,
    y = to_categorical(label_train)[,2:(label_classes+1)],
    generator = do_augmentation)
  # generate batches, for validation data
  flow_vali <- flow_images_from_data(
    x = img_vali,
    y = to_categorical(label_vali)[,2:(label_classes+1)],
    generator = no_augmentation)
  # return ready to use image_generators & steps per epoch
  return(list(data_test = flow_test,
              data_train_vali = flow_train_vali,
              data_train = flow_train,
              data_vali = flow_vali,
              steps_test = floor(length(label_test)/32),
              steps_train_vali = floor(length(label_train_vali)/32),
              steps_train = floor(length(label_train)/32),
              steps_vali = floor(length(label_vali)/32)))
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

get_alexnet <- function(width_length, n_bands, n_band_selector, n_classes, filter_factor) {
  model <- keras_model_sequential() %>%
    # band selector
    layer_conv_2d(input_shape = c(width_length, width_length, n_bands), filters = n_band_selector,
                  kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
    # main model
    # block 1
    layer_conv_2d(filters = 96 * filter_factor, kernel_size = c(11,11), activation = "relu", strides = c(4,4), padding="same") %>%
    #layer_batch_normalization() %>%
    layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
    # block 2
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(5,5), activation = "relu", strides = c(1,1), padding="same") %>%
    #layer_batch_normalization() %>%
    layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
    # block 4
    layer_conv_2d(filters = 384 * filter_factor, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
    layer_conv_2d(filters = 384 * filter_factor, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
    #layer_batch_normalization() %>%
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

get_vgg16 <- function(width_length, n_bands, n_band_selector, n_classes, filter_factor) {
  model <- keras_model_sequential() %>%
    # band selector
    layer_conv_2d(input_shape = c(width_length, width_length, n_bands), filters = n_band_selector,
                  kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
    # main model
    # block 1
    layer_conv_2d(filters = 64 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 64 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    #layer_batch_normalization() %>%
    layer_max_pooling_2d() %>% 
    # block 2
    layer_conv_2d(filters = 128 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 128 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    #layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 3
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 256 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    #layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 4
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    #layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    # block 5
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    layer_conv_2d(filters = 512 * filter_factor, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    #layer_batch_normalization() %>%
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
