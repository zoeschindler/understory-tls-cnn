################################################################################
################################################################################
# TRAIN CNN
################################################################################
################################################################################

# load packages
library(keras)
library(tfruns)

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/07_setup_cnn.R")

# set paths
path_experiment <- "H:/Daten/Studium/2_Master/4_Semester/5_Analyse/08_experiment.R"  # input
path_clips      <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/model_input_2cm/tls_rgb_geo"  # input
input_type      <- basename(path_clips)  # input
path_tfruns     <- paste0("H:/Daten/Studium/2_Master/4_Semester/4_Daten/tfruns/", input_type)  # output
path_models     <- paste0("H:/Daten/Studium/2_Master/4_Semester/4_Daten/models/", input_type)  # output

# create folders
check_create_dir(dirname(path_tfruns))
check_create_dir(dirname(path_models))
check_create_dir(path_tfruns)
check_create_dir(path_models)

# set data input parameters
n_folds <- 5  # number of folds for cross-validation
n_classes <- 5  # number of understory classes
width_length <- 25  # input image dimensions
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
tif_to_rds(path_clips, width_length, n_bands, n_folds, seed=123)

################################################################################
# 5-FOLD CV
################################################################################

hyperparameters <- list()
histories       <- list()
time_training   <- c()
test_accuracies <- c()
all_runs <- c()

for (fold in 1:n_folds) {
  # load & prepare input data
  rdata_paths <- list.files(path_clips, pattern="[.]rds", full.names=TRUE)
  balanced <- create_dataset(rdata_paths, fold, width_length, n_bands)
  
  # try hyperparameter combinations
  runs <- tuning_run(path_experiment,
                     flags = list(learning_rate = c(1e-3, 1e-4, 1e-5),
                                  dropout = c(0.2, 0.3, 0.4),
                                  l2_regularizer = c(0, 0.001, 0.01),
                                  epochs = c(50, 100, 150),
                                  batch_normalization = c(TRUE, FALSE),
                                  filter_factor = c(0.25, 0.5, 1, 1.25, 1.5),
                                  band_selector = c(0.25, 0.5, 0.75)),
                     # batch_size = c(8, 16, 32)),
                     sample = 0.001, # set higher later
                     confirm = FALSE,
                     runs_dir = paste0(path_tfruns, "/fold_", fold))
  
  # get best hyperparameter values
  all_runs <- rbind(all_runs, runs)
  best_run <- ls_runs(order = metric_val_loss, decreasing=F, runs_dir = paste0(path_tfruns, "/fold_", fold))[1,]
  # best_run <- ls_runs(order = metric_val_accuracy, decreasing=T, runs_dir = current_tfruns_dir)[1,]
  
  # set up "best network"
  model <- get_lenet5(width_length = width_length,
                      n_bands = n_bands,
                      n_band_selector = floor(n_bands*best_run$flag_band_selector),
                      n_classes = n_classes,
                      filter_factor = best_run$flag_filter_factor,
                      l2_regualarizer = best_run$flag_l2_regularizer,
                      batch_normalization = best_run$flag_batch_normalization)
  
  # compile "best network"
  model %>% compile(
    optimizer = optimizer_rmsprop(lr = best_run$flag_learning_rate,
                                  decay = best_run$flag_learning_rate / best_run$flag_epochs),
    loss = "categorical_crossentropy",
    metrics = c("accuracy")
  )
  
  # # set callbacks
  # callbacks_list <- list(
  #   callback_early_stopping(monitor = "val_loss", mode = "min", patience = 10, restore_best_weights = TRUE),
  #   callback_model_checkpoint(filepath = paste0(path_models, "/", input_type, "_", fold, ".h5"),
  #                             monitor = "val_acc", save_best_only = TRUE))
  
  # fit "best network"
  start_time <- Sys.time()
  history <- model %>% fit(
    balanced$data_train_vali,
    steps_per_epoch = balanced$steps_train_vali,
    epochs = best_run$flag_epochs,
    # batch_size = FLAGS$batch_size,
    # callbacks = callbacks_list,
    # validation_data = balanced$data_test,
    # validation_steps = balanced$steps_test
  )
  end_time <- Sys.time()
  
  # evaluate "best network"
  results <- model %>% evaluate(balanced$data_test, steps=balanced$steps_test)
  
  # save assessment values & hyperparameter configuration & time for training
  hyperparameters[[fold]] <- best_run[, grepl("flag_", names(best_run))] # TODO: test
  histories[[fold]] <- history
  time_training[fold] <- as.numeric(end_time - start_time)
  test_accuracies[fold] <- results["accuracy"]
}

# TODO: do this with higher sample size,
#       delete completely bad combinations,
#       let it run completely

# TODO: okay to use test data for callback? probably not? what's the alternative?

################################################################################

# compare with random accuracy
random_accuracy <- c()
for (i in 1:100) {
  test_labels <- readRDS(rdata_paths[4])$label
  test_labels_random <- sample(test_labels)
  acc <- length(which(test_labels == test_labels_random)) / length(test_labels)
  random_accuracy <- c(random_accuracy, acc)
}
mean(random_accuracy) # around 0.23

################################################################################