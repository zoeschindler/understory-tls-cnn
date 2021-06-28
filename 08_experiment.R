################################################################################
################################################################################
# EXPERIMENT
################################################################################
################################################################################

################################################################################
# FLAGS
################################################################################

FLAGS <- flags(flag_numeric("learning_rate", 1e-4),
               flag_numeric("dropout", 0.5),
               flag_numeric("l2_regularizer", 0.0001),
               flag_integer("epochs", 100),
               flag_boolean("batch_normalization", FALSE),
               flag_numeric("filter_factor", 1),
               flag_numeric("band_selector", 0.75))

################################################################################
# MODEL
################################################################################

model <- get_lenet5(width_length = width_length,
                    n_bands = n_bands,
                    n_band_selector = floor(n_bands*FLAGS$band_selector),
                    n_classes = n_classes,
                    filter_factor = FLAGS$filter_factor,
                    l2_regularizer = FLAGS$l2_regularizer,
                    batch_normalization = FLAGS$batch_normalization)

################################################################################
# COMPILE
################################################################################

model %>% compile(
  optimizer = optimizer_rmsprop(lr = FLAGS$learning_rate,
                                decay = FLAGS$learning_rate / FLAGS$epochs),
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

################################################################################
# TRAIN
################################################################################

# set callbacks
callbacks_list <- list(callback_early_stopping(monitor = "val_loss", mode = "min", patience = 10))

# fit
history <- model %>% fit_generator(
  balanced$data_train,
  steps_per_epoch = balanced$steps_train,
  epochs = FLAGS$epochs,
  callbacks = callbacks_list,
  validation_data = balanced$data_vali,
  validation_steps = balanced$length_vali
)

################################################################################