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

balance_by_duplicates <- function(label_vector, image_array, max_per_image, max_length) {
  # duplicates labels & images per label until:
  # - each images was duplicated max_per_image times
  # - amount of samples per label was raised to max_length respectively
  repeated_indices <- c()
  for (label in unique(label_vector)) {
    label_indices <- which(label_vector == label)
    new_length <- ifelse(max_per_image * length(label_indices) > max_length,
      max_length,
      max_per_image * length(label_indices)
    )
    label_indices <- rep(label_indices, length.out = new_length)
    repeated_indices <- c(repeated_indices, label_indices)
  }
  label_vector <- label_vector[repeated_indices]
  image_array <- image_array[repeated_indices, , , ]
  # shuffle data
  new_idx_train <- sample(length(label_vector))
  label_vector <- label_vector[new_idx_train]
  image_array <- image_array[new_idx_train, , , ]
  # return new images & labels in list
  return(list(label = label_vector, img = image_array))
}

################################################################################
# READING IN IMAGES
################################################################################

tif_to_rds <- function(clip_dir, pixels, bands, folds = 5, seed = 123) {
  # reads in images, divide into stratified k-folds, save as rds
  # set seed
  set.seed(seed)
  # get all file paths, excluding smallest groups
  img_paths_all <- list.files(clip_dir, pattern = "[.]tif", full.names = TRUE)
  img_paths_all <- img_paths_all[!grepl("grass", img_paths_all) & !grepl("rock", img_paths_all)]
  # get all labels
  img_labels_all <- sapply(strsplit(basename(img_paths_all), "[.]"), "[[", 1)
  img_labels_all <- gsub("[[:digit:]]", "", img_labels_all)
  img_labels_all <- sapply(img_labels_all, function(chars) substr(chars, 1, nchar(chars) - 1))
  img_labels_old <- as.factor(img_labels_all)
  img_labels_all <- as.numeric(img_labels_old)
  # save which label corresponds to which number
  label_lookup <- data.frame("old" = unique(img_labels_old), "new" = unique(img_labels_all))
  write.csv(label_lookup, paste0(clip_dir, "/label_lookup.csv"), row.names = FALSE)
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
    img_array <- array(dim = c(length(fold_files), pixels, pixels, bands))
    label_vector <- c()
    for (i in 1:length(fold_files)) {
      path <- fold_files[i]
      label_vector <- c(label_vector, fold_labels[i])
      raster <- raster::values(stack(path))
      img_array[i, , , ] <- raster
    }
    # save image array as rds
    data_list <- list(img = img_array, label = label_vector)
    saveRDS(data_list, file = paste0(clip_dir, "/images_fold_", l, ".rds"))
    output_paths <- c(output_paths, paste0(clip_dir, "/images_fold_", l, ".rds"))
  }
  # remove seed
  set.seed(NULL)
  return(output_paths)
}

################################################################################
# DATA PREPARATION & AUGMENTATION
################################################################################

create_dataset <- function(rdata_list, holdout_fold, pixels, bands, max_per_image = 100,
                           max_length = 600, balance_classes = TRUE, na_replacement = 0) {
  # create dataset for CNN
  # load test data
  raw_test <- readRDS(rdata_list[holdout_fold])
  img_test <- raw_test$img
  label_test <- raw_test$label
  rm(raw_test)
  # load training folds
  img_train_vali <- array(dim = c(0, pixels, pixels, bands))
  label_train_vali <- c()
  for (i in 1:length(rdata_list)) {
    if (i != holdout_fold) {
      raw_train_vali <- readRDS(rdata_list[i])
      img_train_vali <- unname(abind(img_train_vali, raw_train_vali$img, along = 1))
      label_train_vali <- c(label_train_vali, raw_train_vali$label)
      rm(raw_train_vali)
    }
  }
  # replace NA in images
  img_test[is.na(img_test)] <- na_replacement
  img_train_vali[is.na(img_train_vali)] <- na_replacement
  # make stratified split
  fold_indices <- strat_folds(label_train_vali, 5) # 5 folds -> use 20% for validation
  indices_train <- unlist(fold_indices[1:4])
  indices_vali <- fold_indices[[5]]
  img_train <- img_train_vali[indices_train, , , ]
  label_train <- label_train_vali[indices_train]
  img_vali <- img_train_vali[indices_vali, , , ]
  label_vali <- label_train_vali[indices_vali]
  # shuffle datasets
  new_idx_test <- sample(length(label_test))
  label_test <- label_test[new_idx_test]
  img_test <- img_test[new_idx_test, , , ]
  ###
  new_idx_train_vali <- sample(length(label_train_vali))
  label_train_vali <- label_train_vali[new_idx_train_vali]
  img_train_vali <- img_train_vali[new_idx_train_vali, , , ]
  ###
  new_idx_train <- sample(length(label_train))
  label_train <- label_train[new_idx_train]
  img_train <- img_train[new_idx_train, , , ]
  ###
  new_idx_vali <- sample(length(label_vali))
  label_vali <- label_vali[new_idx_vali]
  img_vali <- img_vali[new_idx_vali, , , ]
  ###
  if (balance_classes) {
    # duplicate images depending on label frequency
    duplicated_train_vali <- balance_by_duplicates(label_train_vali, img_train_vali, max_per_image, max_length)
    label_train_vali <- duplicated_train_vali$label
    img_train_vali <- duplicated_train_vali$img
    ###
    duplicated_train <- balance_by_duplicates(label_train, img_train, max_per_image, max_length)
    label_train <- duplicated_train$label
    img_train <- duplicated_train$img
    # shuffle datasets
    new_idx_train_vali <- sample(length(label_train_vali))
    label_train_vali <- label_train_vali[new_idx_train_vali]
    img_train_vali <- img_train_vali[new_idx_train_vali, , , ]
    ###
    new_idx_train <- sample(length(label_train))
    label_train <- label_train[new_idx_train]
    img_train <- img_train[new_idx_train, , , ]
  }
  # empty data generator, for test / validation data
  no_augmentation <- image_data_generator()
  # data generator with data augmentation, for training data
  do_augmentation <- image_data_generator(
    fill_mode = "reflect",
    rotation_range = 20,
    width_shift_range = 0.1,
    height_shift_range = 0.1,
    shear_range = 20,
    horizontal_flip = TRUE,
    vertical_flip = TRUE
  )
  # get number of classes
  label_classes <- length(unique(label_train))
  # generate batches, for test data
  flow_test <- flow_images_from_data(
    x = img_test,
    y = to_categorical(label_test)[, 2:(label_classes + 1)],
    generator = no_augmentation,
    batch_size = 1, shuffle = FALSE
  )
  # generate batches, for training + validation data
  flow_train_vali_16 <- flow_images_from_data(
    x = img_train_vali,
    y = to_categorical(label_train_vali)[, 2:(label_classes + 1)],
    generator = do_augmentation,
    batch_size = 16
  )
  # generate batches, for training data
  flow_train_16 <- flow_images_from_data(
    x = img_train,
    y = to_categorical(label_train)[, 2:(label_classes + 1)],
    generator = do_augmentation,
    batch_size = 16
  )
  # generate batches, for training + validation data
  flow_train_vali_32 <- flow_images_from_data(
    x = img_train_vali,
    y = to_categorical(label_train_vali)[, 2:(label_classes + 1)],
    generator = do_augmentation,
    batch_size = 32
  )
  # generate batches, for training data
  flow_train_32 <- flow_images_from_data(
    x = img_train,
    y = to_categorical(label_train)[, 2:(label_classes + 1)],
    generator = do_augmentation,
    batch_size = 32
  )
  # generate batches, for validation data
  flow_vali <- flow_images_from_data(
    x = img_vali,
    y = to_categorical(label_vali)[, 2:(label_classes + 1)],
    generator = no_augmentation,
    batch_size = 1, shuffle = FALSE
  )
  # return ready to use image_generators & steps per epoch
  return(list(
    data_test = flow_test,
    data_train_vali_16 = flow_train_vali_16,
    data_train_16 = flow_train_16,
    data_train_vali_32 = flow_train_vali_32,
    data_train_32 = flow_train_32,
    data_vali = flow_vali,
    length_test = length(label_test),
    length_vali = length(label_vali),
    steps_train_vali_16 = round(length(label_train_vali) / 16),
    steps_train_16 = round(length(label_train) / 16),
    steps_train_vali_32 = round(length(label_train_vali) / 32),
    steps_train_32 = round(length(label_train) / 32)
  ))
}

################################################################################
# LeNet-5
# http://yann.lecun.com/exdb/publis/pdf/lecun-98.pdf
# https://www.kaggle.com/curiousprogrammer/lenet-5-cnn-with-keras-99-48 (based on this code)
# https://github.com/BIGBALLON/cifar-10-cnn/blob/master/1_Lecun_Network/LeNet_dp_da_wd_keras.py (and this code)
################################################################################

get_lenet5 <- function(width_length, n_bands, n_classes, dropout = 0.5, l2_regularizer = 0.01) {
  model <- keras_model_sequential() %>%
    layer_conv_2d(
      input_shape = c(width_length, width_length, n_bands),
      filters = 32, kernel_size = c(5, 5), padding = "same",
      kernel_regularizer = regularizer_l2(l2_regularizer),
      kernel_initializer = "he_normal"
    ) %>%
    layer_activation_leaky_relu() %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    layer_conv_2d(
      filters = 48, kernel_size = c(5, 5), padding = "same",
      kernel_regularizer = regularizer_l2(l2_regularizer),
      kernel_initializer = "he_normal"
    ) %>%
    layer_activation_leaky_relu() %>%
    layer_batch_normalization() %>%
    layer_max_pooling_2d() %>%
    layer_flatten() %>%
    layer_dense(
      units = 256,
      kernel_regularizer = regularizer_l2(l2_regularizer),
      kernel_initializer = "he_normal"
    ) %>%
    layer_activation_leaky_relu() %>%
    layer_dropout(dropout) %>%
    layer_dense(
      units = 84,
      kernel_regularizer = regularizer_l2(l2_regularizer),
      kernel_initializer = "he_normal"
    ) %>%
    layer_activation_leaky_relu() %>%
    layer_dropout(dropout) %>%
    layer_dense(
      units = n_classes, activation = "softmax",
      kernel_initializer = "he_normal"
    )
  return(model)
}

################################################################################