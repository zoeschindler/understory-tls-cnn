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
tile_size  <- 0.5

# set paths
path_experiment <- "C:/Users/Zoe/Documents/understory_classification/5_Analyse/08_experiment.R"  # input
path_clips      <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/model_input_", resolution*100, "cm_normalized/", input_type)  # input
#path_clips      <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/model_input_", resolution*100, "cm_standardized/", input_type)  # input
path_tfruns     <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/tfruns_", resolution*100, "cm/", input_type)  # output
path_models     <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/models_", resolution*100, "cm/", input_type)  # output

# create folders
check_create_dir(dirname(path_tfruns))
check_create_dir(dirname(path_models))
check_create_dir(path_tfruns)
check_create_dir(path_models)

# set data input parameters
n_folds <- 10  # number of folds for cross-validation
n_classes <- 5  # number of understory classes
width_length <- tile_size / resolution  # input image dimensions (pixels)
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
tif_to_rds(path_clips, width_length, n_bands, n_folds, seed=666)

################################################################################
# K-FOLD CV
################################################################################

# create empty objects for variables
time_training   <- c()
test_accuracies <- c()
all_runs        <- c()
best_runs       <- c()
predictions     <- c()

# loop through all folds
for (fold in 1:n_folds) {
  # load & prepare input data
  rdata_paths <- list.files(path_clips, pattern="[.]rds", full.names=TRUE)
  balanced <- create_dataset(rdata_paths, fold, width_length, n_bands)
  
  # try hyperparameter combinations
  runs <- tuning_run(path_experiment,
                     flags = list(learning_rate = c(1e-4, 1e-5),
                                  l2_regularizer = c(0, 0.001, 0.01),
                                  epochs = c(200, 300, 400)),
                     confirm = FALSE,
                     runs_dir = paste0(path_tfruns, "/fold_", fold))
  
  # add runs to data frame
  colnames(runs)[2] <- "eval_loss"
  colnames(runs)[length(colnames(runs))] <- "eval_accuracy"
  all_runs <- rbind(all_runs, cbind(runs, fold))  # add fold
  
  # get best hyperparameter values
  best_run <- ls_runs(order = metric_val_accuracy, decreasing=T, runs_dir = paste0(path_tfruns, "/fold_", fold))[1,]
  #best_run <- ls_runs(order = metric_val_loss, decreasing=F, runs_dir = paste0(path_tfruns, "/fold_", fold))[1,]
  colnames(best_run)[3] <- "eval_loss"
  colnames(best_run)[length(colnames(best_run))] <- "eval_accuracy"
  best_runs <- rbind(best_runs, cbind(best_run,fold))
  
  # set up "best network"
  model <- get_lenet5(
    width_length = width_length,
    n_bands = n_bands,
    n_classes = n_classes,
    l2_regularizer = best_run$flag_l2_regularizer)
  
  # compile "best network"
  model %>% compile(
    optimizer = optimizer_rmsprop(lr = best_run$flag_learning_rate,
                                  decay = best_run$flag_learning_rate / best_run$flag_epochs),
    loss = "categorical_crossentropy",
    metrics = c("accuracy")
  )
  
  # fit "best network"
  start_time <- Sys.time()
  history <- model %>% fit_generator(
    balanced$data_train_vali,
    steps_per_epoch = balanced$steps_train_vali,
    epochs = best_run$flag_epochs,
  )
  end_time <- Sys.time()
  
  # save model
  model %>% save_model_hdf5(paste0(path_models, "/best_of_fold_", fold,".h5"))
  
  # evaluate "best network"
  results <- model %>% evaluate_generator(balanced$data_test, steps=balanced$length_test)
  
  # predictions "best network"
  image_test <- array(dim=c(0, width_length, width_length, n_bands))
  label_test <- c()
  for (i in 1:balanced$length_test) {
    item <- generator_next(balanced$data_test)
    image_test <- abind(image_test, item[[1]], along=1)
    label_test <- rbind(label_test, item[[2]])
  }
  preds_test <- model %>% predict(image_test)
  label_test <- apply(label_test, 1, which.max)
  preds_test <- apply(preds_test, 1, which.max)
  pred_df <- data.frame("predictions" = preds_test, "truth" = label_test, "fold" = fold)
  predictions <- rbind(predictions, pred_df)

  # save assessment values & hyperparameter configuration & time for training
  time_training[fold] <- as.numeric(difftime(end_time, start_time, units = "sec"))
  test_accuracies[fold] <- results["accuracy"]
}

# save all run information to file
all_runs_data <- all_runs[grepl("metric_", names(all_runs)) | grepl("flag_", names(all_runs)) | grepl("epoch", names(all_runs)) | grepl("eval_", names(all_runs)) | grepl("fold", names(all_runs))]
write.csv(all_runs_data, paste0(path_models, "/all_hyperparameters.csv"), row.names=FALSE)

# save best run information to file
best_runs_data <- best_runs[grepl("metric_", names(best_runs)) | grepl("flag_", names(best_runs)) | grepl("epoch", names(best_runs)) | grepl("eval_", names(best_runs))]
best_runs_data <- cbind(best_runs_data, data.frame("test_accuracy" = test_accuracies, "train_time" = time_training))
write.csv(best_runs_data, paste0(path_models, "/best_hyperparameters.csv"), row.names=FALSE)

# save predictions to file
write.csv(predictions, paste0(path_models, "/prediction_truth_fold.csv"), row.names=FALSE)

# ################################################################################
# # TESTING
# 
# source("C:/Users/Zoe/Documents/understory_classification/5_Analyse/07_setup_cnn.R")
# 
# # testing network stability
# set.seed(NULL)
# vali_accs  <- c()
# test_accs  <- c()
# final_accs <- c()
# vali_loss  <- c()
# test_loss  <- c()
# final_loss <- c()
# rdata_paths <- list.files(path_clips, pattern="[.]rds", full.names=TRUE)
# balanced <- create_dataset(rdata_paths, 2, width_length, n_bands)
# for (i in 1:5) {
#   rm(model)
#   model <- get_lenet5(
#     width_length = width_length,
#     n_bands = n_bands,
#     n_classes = n_classes,
#     l2_regularizer = 1e-2)
#   model %>% compile(
#     optimizer = optimizer_rmsprop(lr = 1e-4,
#                                   decay = 1e-4 / 200),
#     loss = "categorical_crossentropy",
#     metrics = c("accuracy")
#   )
#   history <- model %>% fit_generator(
#     balanced$data_train,
#     steps_per_epoch = balanced$steps_train,
#     epochs = 200,
#     validation_data = balanced$data_vali,
#     validation_steps = balanced$length_vali
#   )
#   vali <- model %>% evaluate_generator(balanced$data_vali, steps=balanced$length_vali)
#   test <- model %>% evaluate_generator(balanced$data_test, steps=balanced$length_test)
#   vali_accs <- c(vali_accs, vali["accuracy"])
#   test_accs <- c(test_accs, test["accuracy"])
#   vali_loss <- c(vali_loss, vali["loss"])
#   test_loss <- c(test_loss, test["loss"])
#   rm(model)
#   model <- get_lenet5(
#     width_length = width_length,
#     n_bands = n_bands,
#     n_classes = n_classes,
#     l2_regularizer = 1e-2)
#   model %>% compile(
#     optimizer = optimizer_rmsprop(lr = 1e-4,
#                                   decay = 1e-4 / 200),
#     loss = "categorical_crossentropy",
#     metrics = c("accuracy")
#   )
#   history <- model %>% fit_generator(
#     balanced$data_train_vali,
#     steps_per_epoch = balanced$steps_train_vali,
#     epochs = 200
#   )
#   final <- model %>% evaluate_generator(balanced$data_test, steps=balanced$length_test)
#   final_accs <- c(final_accs, final["accuracy"])
#   final_loss <- c(final_loss, final["loss"])
# }
# vali_loss
# test_loss
# final_loss
# vali_accs
# test_accs
# final_accs
# summary(vali_accs)
# summary(test_accs)
# summary(final_accs)
# 
# # goal 1: less accuracy variations
# # goal 2: higher accuracies
# 
# ################################################################################

# compare with random accuracy
random_accuracy <- c()
for (i in 1:100) {
  test_labels <- readRDS(rdata_paths[4])$label
  test_labels_random <- sample(test_labels)
  acc <- length(which(test_labels == test_labels_random)) / length(test_labels)
  random_accuracy <- c(random_accuracy, acc)
}
mean(random_accuracy) # around 0.23

# other approach: take share of most occuring class (package caret)

################################################################################