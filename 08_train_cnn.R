################################################################################
################################################################################
# TRAIN CNN
################################################################################
################################################################################

# load functions
source("H:/Daten/Studium/2_Master/4_Semester/5_Analyse/07_setup_cnn.R")

# set paths
path_clips  <- "H:/Daten/Studium/2_Master/4_Semester/4_Daten/models/input_unfiltered/tls_rgb_geo"  # input
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
n_band_selector <- n_bands %/% 2  # number of layers after band selector
n_classes <- 5  # number of understory classes
filter_factor <- 0.5  # multiplier for amount of filters per layer (0.25 / 0.5 / 1 / 2)

################################################################################
# GET DATA & TRAIN
################################################################################

# execute only once because it takes time
tif_to_rds(path_clips, path_output, width_length, n_bands, seed=123)

# later execute in a loop k times, also train & assess model k times
rdata_paths <- list.files(path_output, pattern="[.]rds", full.names=TRUE)
imbalanced <- create_dataset(rdata_paths, 4, balance_classes=FALSE)
balanced <- create_dataset(rdata_paths, 4, balance_classes=TRUE)

# load CNN
model <- get_vgg16(width_length, n_bands, n_band_selector, n_classes, filter_factor)

# compile CNN
model %>% compile(
  optimizer = optimizer_rmsprop(lr = 1e-5),
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

# train CNN
history <- model %>% fit(
  balanced$data_train,
  steps_per_epoch = balanced$steps_train,
  epochs = 10,
  validation_data = balanced$data_test,
  validation_steps = balanced$steps_test
)

# assess CNN
plot(history)

################################################################################