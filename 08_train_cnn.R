################################################################################
################################################################################
# TRAIN CNN
################################################################################
################################################################################

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/07_setup_cnn.R")

# set paths
path_clips  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/model_input_1cm/tls_rgb_geo"  # input

# set data input parameters
width_length <- 50  # number of pixels (depends on resolution!)
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
n_classes <- 5  # number of understory classes
n_folds <- 5  # number of folds for k-fold crossvalidation
filter_factor <- 0.25  # multiplier for amount of filters per layer (0.25 / 0.5 / 1 / 2)

################################################################################
# GET DATA & TRAIN
################################################################################

# execute only once because it takes time
tif_to_rds(path_clips, width_length, n_bands, n_folds, seed=123)

################################################################################
# once: fold 4 as holdout

# later execute in a loop k times, also train & assess model k times
rdata_paths <- list.files(path_clips, pattern="[.]rds", full.names=TRUE)
balanced <- create_dataset(rdata_paths, 4, width_length, n_bands, balance_classes=TRUE)

# load CNN
model <- get_vgg16(width_length, n_bands, n_band_selector, n_classes, filter_factor)

# compile CNN
model %>% compile(
  optimizer = optimizer_rmsprop(lr = 1e-5, decay = 1e-7),
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

# train CNN
history <- model %>% fit(
  balanced$data_train_vali,
  steps_per_epoch = balanced$steps_train_vali,
  epochs = 50,
  validation_data = balanced$data_test,
  validation_steps = balanced$steps_test
)

# assess CNN
plot(history)

################################################################################
# 5-fold crossvalidation (TODO)

# hist_list <- list()
# test_acc_list <- c()
# 
# hyper_options <- list()
# 
# for (fold in 1:n_folds) {
#   # load & prepare input data
#   rdata_paths <- list.files(path_clips, pattern="[.]rds", full.names=TRUE)
#   balanced <- create_dataset(rdata_paths, fold, width_length, n_bands, balance_classes=TRUE)
#   
#   ### wrap into loop
#   # load CNN
#   model <- get_vgg16(width_length, n_bands, n_band_selector, n_classes, filter_factor)
#   
#   # compile CNN
#   model %>% compile(
#     optimizer = optimizer_rmsprop(lr = 1e-5),
#     loss = "categorical_crossentropy",
#     metrics = c("accuracy")
#   )
#   
#   # train CNN
#   history <- model %>% fit(
#     balanced$data_train,
#     steps_per_epoch = balanced$steps_train,
#     epochs = 10,
#     validation_data = balanced$data_test,
#     validation_steps = balanced$steps_test
#   )
#   
#   # TODO: add callback to save model with best validation accuracy per configuration
#   # TODO: add callback for early stopping
#   
#   ###
#   
#   # TODO get best parameter configuration
#   
#   # TODO: train with validation + training data (both augmented)
#   
#   # TODO: save model with best hyperparameters + all data
#   
#   # TODO: test with testdata
#   
#   # save history & accuracy & hyperparameter configuration
#   hist_list[[fold]] <- history
# }

################################################################################