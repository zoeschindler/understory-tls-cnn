################################################################################
################################################################################
# EXPERIMENT
################################################################################
################################################################################

################################################################################
# FLAGS
################################################################################

FLAGS <- flags(flag_numeric("learning_rate", 5e-5),
               flag_integer("epochs", 200),
               flag_integer("batch_size", 32),
               flag_numeric("dropout", 0.5),
               flag_numeric("l2_regularizer", 0.01))

################################################################################
# MODEL
################################################################################

model <- get_lenet5(
  width_length = width_length,
  n_bands = n_bands,
  n_classes = n_classes,
  dropout = FLAGS$dropout,
  l2_regularizer = FLAGS$l2_regularizer)

################################################################################
# COMPILE
################################################################################

model %>% compile(
  optimizer = optimizer_adam(lr = FLAGS$learning_rate),
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

################################################################################
# TRAIN
################################################################################

# fit
if (FLAGS$batch_size == 16) {
  history <- model %>% fit_generator(
    balanced$data_train_16,
    steps_per_epoch = balanced$steps_train_16,
    epochs = FLAGS$epochs,
    validation_data = balanced$data_vali,
    validation_steps = balanced$length_vali
  )
}
if (FLAGS$batch_size == 32) {
  history <- model %>% fit_generator(
    balanced$data_train_32,
    steps_per_epoch = balanced$steps_train_32,
    epochs = FLAGS$epochs,
    validation_data = balanced$data_vali,
    validation_steps = balanced$length_vali
  )
}

################################################################################
# SAVE
################################################################################

check_create_dir(paste0(path_models, "/all"))
name <- paste0("lr_", FLAGS$learning_rate, "_ep_", FLAGS$epochs, "_bs_",
               FLAGS$batch_size, "_dp_", FLAGS$dropout, "_l2_", FLAGS$l2_regularizer)
model %>% save_model_hdf5(paste0(path_models, "/all/fold_", fold, "_", name, ".h5"))

################################################################################
# EVALUATE
################################################################################

results <- model %>% evaluate_generator(balanced$data_test, steps=balanced$length_test)

################################################################################