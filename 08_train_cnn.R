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
resolution <- 0.01
input_type <- "tls_rgb_geo"

# set paths
path_experiment <- "C:/Users/Zoe/Documents/understory_classification/5_Analyse/08_experiment.R"  # input
path_clips      <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/model_input_", resolution*100, "cm/", input_type)  # input
path_tfruns     <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/tfruns_", resolution*100, "cm/", input_type)  # output
path_models     <- paste0("C:/Users/Zoe/Documents/understory_classification/4_Daten/models_", resolution*100, "cm/", input_type)  # output

# create folders
check_create_dir(dirname(path_tfruns))
check_create_dir(dirname(path_models))
check_create_dir(path_tfruns)
check_create_dir(path_models)

# set data input parameters
n_folds <- 5  # number of folds for cross-validation
n_classes <- 5  # number of understory classes
width_length <- 50  # input image dimensions (pixels)
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
                     flags = list(learning_rate = c(1e-3, 1e-4, 1e-5),
                                  dropout = c(0.3, 0.4, 0.5),
                                  l2_regularizer = c(0, 0.0001, 0.001),
                                  epochs = c(100),
                                  batch_normalization = c(FALSE, TRUE),
                                  filter_factor = c(0.5, 0.75, 1),
                                  band_selector = c(0.5, 0.75, 1)),
                     sample = 0.001, # set higher later
                     confirm = FALSE,
                     runs_dir = paste0(path_tfruns, "/fold_", fold))
  
  # get best hyperparameter values
  all_runs <- rbind(all_runs, runs)
  # best_run <- ls_runs(order = metric_val_loss, decreasing=F, runs_dir = paste0(path_tfruns, "/fold_", fold))[1,]
  best_run <- ls_runs(order = metric_val_accuracy, decreasing=T, runs_dir = paste0(path_tfruns, "/fold_", fold))[1,]
  best_runs <- rbind(best_runs, best_run)
  
  # set up "best network"
  model <- get_lenet5(width_length = width_length,
                      n_bands = n_bands,
                      n_band_selector = floor(n_bands*best_run$flag_band_selector),
                      n_classes = n_classes,
                      filter_factor = best_run$flag_filter_factor,
                      l2_regularizer = best_run$flag_l2_regularizer,
                      batch_normalization = best_run$flag_batch_normalization)
  
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
    #epochs = best_run$flag_epochs,  # old
    epochs = best_run$epochs_completed,  # new
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
  time_training[fold] <- as.numeric(end_time - start_time)
  test_accuracies[fold] <- results["accuracy"]
}

# save all run information to file
all_runs_data <- all_runs[grepl("metric_", names(all_runs)) | grepl("flag_", names(all_runs)) | grepl("epoch", names(all_runs))]
write.csv(all_runs_data, paste0(path_models, "/all_hyperparameters.csv"), row.names=FALSE)

# save best run information to file
best_runs_data <- best_runs[grepl("metric_", names(best_runs)) | grepl("flag_", names(best_runs)) | grepl("epoch", names(best_runs))]
best_runs_data <- cbind(best_runs_data, data.frame("test_accuracy" = test_accuracies, "train_time" = time_training))
write.csv(best_runs_data, paste0(path_models, "/best_hyperparameters.csv"), row.names=FALSE)

# save predictions to file
write.csv(predictions, paste0(path_models, "/prediction_truth_fold.csv"), row.names=FALSE)


################### TESTING STUFF

library(ggplot2)

all_runs <- c()
for (fold in 1:5) {
  fold_runs <- ls_runs(runs_dir = paste0(path_tfruns, "/fold_", fold))
  all_runs <- rbind(all_runs, fold_runs)
}
View(all_runs)

ggplot(all_runs) +
  geom_point(aes(x=metric_val_accuracy, y=metric_val_loss))
ggplot(all_runs) +
  geom_point(aes(x=metric_val_accuracy, y=metric_accuracy))
###
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_learning_rate), y=metric_val_accuracy))
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_learning_rate), y=metric_val_loss))
###
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_dropout), y=metric_val_accuracy))
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_dropout), y=metric_val_loss))
###
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_l2_regularizer), y=metric_val_accuracy))
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_l2_regularizer), y=metric_val_loss))
###
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_epochs), y=metric_val_accuracy))
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_epochs), y=metric_val_loss))
###
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_batch_normalization), y=metric_val_accuracy))
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_batch_normalization), y=metric_val_loss))
###
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_filter_factor), y=metric_val_accuracy))
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_filter_factor), y=metric_val_loss))
###
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_band_selector), y=metric_val_accuracy))
ggplot(all_runs) +
  geom_boxplot(aes(x=as.factor(flag_band_selector), y=metric_val_loss))

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