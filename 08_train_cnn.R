################################################################################
################################################################################
# TRAIN CNN
################################################################################
################################################################################

# load packages
library(abind)
library(keras)
library(tfruns)

# load functions
source("C:/Users/Zoe/Documents/understory_classification/5_Analyse/07_setup_cnn.R")

# set paramters
resolution <- 0.02
input_type <- "tls"
tile_size <- 0.5

# set paths
basedir <- "C:/Users/Zoe/Documents/understory_classification"
path_experiment <- paste0(basedir, "/5_Analyse/08_experiment.R") # input
path_clips      <- paste0(basedir, "/4_Daten/model_input_", resolution * 100, "cm_standardized/", input_type) # input
path_tfruns     <- paste0(basedir, "/4_Daten/tfruns_", resolution * 100, "cm/", input_type) # output
path_models     <- paste0(basedir, "/4_Daten/models_", resolution * 100, "cm/", input_type) # output

# create folders
check_create_dir(dirname(path_tfruns))
check_create_dir(dirname(path_models))
check_create_dir(path_tfruns)
check_create_dir(path_models)

# set data input parameters
n_folds <- 5 # number of folds for cross-validation
n_classes <- 5 # number of understory classes
width_length <- tile_size / resolution # input image dimensions (pixels)
if (input_type == "tls") {
  n_bands <- 4
} else if (input_type == "tls_geo") {
  n_bands <- 10
} else if (input_type == "tls_rgb") {
  n_bands <- 7
} else if (input_type == "tls_rgb_geo") {
  n_bands <- 13
}

################################################################################
# PREPARE FOLD DATA
################################################################################

# execute only once because it takes time
tif_to_rds(path_clips, width_length, n_bands, n_folds, seed = 123)

################################################################################
# HYPERPARAMETER OPTIMIZATION
################################################################################

# create empty objects for variables
all_runs <- c()
best_runs <- c()

# loop through all folds
for (fold in 1:n_folds) {
  # progress info
  print(paste0("fold: ", fold))

  # load & prepare input data
  rdata_paths <- list.files(path_clips, pattern = "[.]rds", full.names = TRUE)
  balanced <- create_dataset(rdata_paths, fold, width_length, n_bands, na_replacement = 0)

  # try hyperparameter combinations
  runs <- tuning_run(path_experiment,
    flags = list(
      learning_rate = c(5e-5, 1e-5),
      epochs = c(200, 300, 400),
      batch_size = c(16, 32),
      dropout = c(0.5),
      l2_regularizer = c(0.01)
    ),
    confirm = FALSE,
    runs_dir = paste0(path_tfruns, "/fold_", fold)
  )

  # add runs to data frame
  colnames(runs)[2] <- "eval_loss"
  colnames(runs)[length(colnames(runs))] <- "eval_accuracy"
  all_runs <- rbind(all_runs, cbind(runs, fold))

  # get best hyperparameter values
  best_run <- ls_runs(
    order = metric_val_accuracy, decreasing = T,
    runs_dir = paste0(path_tfruns, "/fold_", fold)
  )[1, ]
  colnames(best_run)[3] <- "eval_loss"
  colnames(best_run)[length(colnames(best_run))] <- "eval_accuracy"
  best_runs <- rbind(best_runs, cbind(best_run, fold))
}

# save all run information to file
all_runs_data <- all_runs[grepl("metric_", names(all_runs)) | grepl("flag_", names(all_runs)) |
  grepl("epoch", names(all_runs)) | grepl("eval_", names(all_runs)) |
  grepl("fold", names(all_runs))]
write.csv(all_runs_data, paste0(path_models, "/all_hyperparameters.csv"), row.names = FALSE)

# save best run information to file
best_runs_data <- best_runs[grepl("metric_", names(best_runs)) | grepl("flag_", names(best_runs)) |
  grepl("epoch", names(best_runs)) | grepl("eval_", names(best_runs)) |
  grepl("fold", names(best_runs))]
best_runs_data$flag_learning_rate[best_runs_data$flag_learning_rate == 0] <- 5e-5 # I don't know why 5e-5 is saved as zero
write.csv(best_runs_data, paste0(path_models, "/best_hyperparameters.csv"), row.names = FALSE)

################################################################################
# RETRAIN TO SEE INSTABILITY
################################################################################

# how often retrain per fold
retrain_num <- 10

# load data frame with best hyperparameters
best_runs_data <- read.csv(paste0(path_models, "/best_hyperparameters.csv"))
best_runs_data$flag_learning_rate[best_runs_data$flag_learning_rate == 0] <- 5e-5

# prepare empty variables
accuracies_all <- c()
time_training_all <- c()

# loop through folds
for (row in 1:nrow(best_runs_data)) {
  # load data
  rdata_paths <- list.files(path_clips, pattern = "[.]rds", full.names = TRUE)
  balanced <- create_dataset(rdata_paths, best_runs_data$fold[row], width_length, n_bands, na_replacement = 0)

  # train i times to show instability
  accuracies_fold <- c()
  time_training_fold <- c()
  for (i in 1:retrain_num) {
    model <- get_lenet5(
      width_length = width_length,
      n_bands = n_bands,
      n_classes = n_classes,
      dropout = best_runs_data$flag_dropout[row],
      l2_regularizer = best_runs_data$flag_l2_regularizer[row]
    )
    model %>% compile(
      optimizer = optimizer_adam(lr = best_runs_data$flag_learning_rate[row]),
      loss = "categorical_crossentropy",
      metrics = c("accuracy")
    )
    start_time <- Sys.time()
    if (best_runs_data$flag_batch_size[row] == 16) {
      history <- model %>% fit_generator(
        balanced$data_train_vali_16,
        steps_per_epoch = balanced$steps_train_vali_16,
        epochs = best_runs_data$flag_epochs[row],
      )
    }
    if (best_runs_data$flag_batch_size[row] == 32) {
      history <- model %>% fit_generator(
        balanced$data_train_vali_32,
        steps_per_epoch = balanced$steps_train_vali_32,
        epochs = best_runs_data$flag_epochs[row],
      )
    }
    end_time <- Sys.time()

    # get accuracy
    results <- model %>% evaluate_generator(balanced$data_test, steps = balanced$length_test)
    accuracies_fold <- c(accuracies_fold, results["accuracy"])

    # get training time
    time_diff <- as.numeric(difftime(end_time, start_time, units = "sec"))
    time_training_fold <- c(time_training_fold, time_diff)
  }
  accuracies_all <- cbind(accuracies_all, as.numeric(accuracies_fold))
  time_training_all <- cbind(time_training_all, time_training_fold)
}
colnames(accuracies_all) <- paste0("fold_", best_runs_data$fold)
colnames(time_training_all) <- paste0("fold_", best_runs_data$fold)

# save
write.csv(accuracies_all, paste0(path_models, "/accuracy_fluctuations_best_hyperparameters.csv"), row.names = FALSE)
write.csv(time_training_all, paste0(path_models, "/time_fluctuations_best_hyperparameters.csv"), row.names = FALSE)

################################################################################
# PREDICT ON NOT RETRAINED MODELS
################################################################################

# load data frame with best hyperparameters
best_runs_data <- read.csv(paste0(path_models, "/best_hyperparameters.csv"))
best_runs_data$flag_learning_rate[best_runs_data$flag_learning_rate == 0] <- 5e-5

# prepare empty variables
predictions <- c()

# loop through folds
for (row in 1:nrow(best_runs_data)) {
  # load best model
  filename <- paste0(
    "fold_", best_runs_data$fold[row],
    "_lr_", best_runs_data$flag_learning_rate[row],
    "_ep_", best_runs_data$flag_epochs[row],
    "_bs_", best_runs_data$flag_batch_size[row],
    "_dp_", best_runs_data$flag_dropout[row],
    "_l2_", best_runs_data$flag_l2_regularizer[row],
    ".h5"
  )
  model_old <- load_model_hdf5(paste0(path_models, "/all/", filename))

  # load data
  rdata_paths <- list.files(path_clips, pattern = "[.]rds", full.names = TRUE)
  balanced <- create_dataset(rdata_paths, best_runs_data$fold[row], width_length, n_bands, na_replacement = 0)

  # make predictions
  image_test <- array(dim = c(0, width_length, width_length, n_bands))
  label_test <- c()
  for (i in 1:balanced$length_test) {
    item <- generator_next(balanced$data_test)
    image_test <- abind(image_test, item[[1]], along = 1)
    label_test <- rbind(label_test, item[[2]])
  }
  preds_test <- model_old %>% predict(image_test)
  label_test <- apply(label_test, 1, which.max)
  preds_test <- apply(preds_test, 1, which.max)
  pred_df <- data.frame("predictions" = preds_test, "truth" = label_test, "fold" = best_runs_data$fold[row])
  predictions <- rbind(predictions, pred_df)
}

# save
write.csv(predictions, paste0(path_models, "/predictions_not_retrained.csv"), row.names = FALSE)

################################################################################