################################################################################
################################################################################
# TESTING ARCHITECTURE ON CIFAR-10
################################################################################
################################################################################

# load packages
library(keras)
library(ggplot2)
library(dplyr)
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)

# set paths
basedir <- "H:/Daten/Studium/2_Master/4_Semester"
path_models <- paste0(basedir, "/4_Daten/02_testing/models_2cm")
path_out    <- paste0(basedir, "/cifar10_test")

# own colors
own_colors_named <- list(
  red = "#ff6f69",
  blue = "#00aedb",
  yellow = "#ffcc5c",
  green = "#88d8b0",
  turquoise = "#4abdac",
  pink = "#ff8b94",
  bright_green = "#cbe885"
)

################################################################################
# CNN ARCHITECTURE
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
# TRAIN MODELS
################################################################################

# load data
cifar10 <- dataset_cifar10()
x_train <- cifar10$train$x/255
x_test <- cifar10$test$x/255
y_train <- to_categorical(cifar10$train$y, num_classes = 10)
y_test <- to_categorical(cifar10$test$y, num_classes = 10)

# sample down training data so it does not take ages
random_subsample <- sample(1:50000, size = 5000)
x_train <- x_train[random_subsample, , , ]
y_train <- y_train[random_subsample, ]

# get all best hyperparameter combinations
all_confs <- c()
for (type in c("tls", "tls_geo", "tls_rgb", "tls_rgb_geo")) {
  best_runs <- read.csv(paste0(path_models, "/", type, "/best_hyperparameters.csv"))
  best_runs$flag_learning_rate[best_runs$flag_learning_rate == 0] <- 5e-5
  all_confs <- rbind(all_confs, data.frame(
    "lr" = best_runs$flag_learning_rate,
    "bs" = best_runs$flag_batch_size,
    "ep" = best_runs$flag_epochs
  ))
}
all_confs <- all_confs[!duplicated(all_confs),]
all_confs$id <- 1:nrow(all_confs)

# test all hyperpaameter combinations 10 times
acc_all <- c()
for (id in 1:nrow(all_confs)) {
  acc_id <- c()
  for (rep in 1:10) {
    # get blank model
    model <- get_lenet5(width_length = 32, n_bands = 3, n_classes = 10)
    # compile model
    model %>% compile(
      optimizer = optimizer_adam(lr = all_confs$lr[id]),
      loss = "categorical_crossentropy",
      metrics = c("accuracy")
    )
    # train model
    history <- model %>% fit(
      x_train, y_train,
      batch_size = all_confs$bs[id],
      shuffle = TRUE,
      epochs = all_confs$ep[id]
    )
    # evaluate model
    results <- model %>% evaluate(x_test, y_test)
    acc_id <- c(acc_id, results["accuracy"])
  }
  acc_all <- rbind(acc_all, data.frame("accuracy" = acc_id, "id" = id))
}

# add hyperparameter information
acc_all <- acc_all %>%
  left_join(all_confs, by = c("id" ="id"))
acc_all$repetitions <- 10

# save results
write.csv(acc_all, paste0(path_out , "/cifar10_acc.csv"), row.names = FALSE)

################################################################################
# VISUALIZE RESULTS
################################################################################

# load data
acc_all <- read.csv(paste0(path_out , "/cifar10_acc.csv"))

# reshape
acc_min_max_diff <- acc_all %>%
  group_by(id) %>%
  summarise(min = round(min(accuracy),2),
            max = round(max(accuracy),2),
            diff = round(abs(max(accuracy) - min(accuracy)),2))
print(paste0("max diff: ", max(acc_min_max_diff$diff)))

# plot results
cairo_pdf(
  file = paste0(path_out, "/cifar10_acc_lr.pdf"),
  family = "Calibri", width = 8.27, height = 5.83
)
ggplot(acc_all, aes(x = as.factor(id), y = accuracy)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = as.factor(lr)), position = position_dodge(0.9)) +
  scale_x_discrete(labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")) +
  scale_fill_manual(
    values = c("1e-05" = own_colors_named$blue, "5e-05" = own_colors_named$red),
    name = "Learning Rate\n"
  ) +
  xlab("\nHyperparameter Combination") +
  ylab("Test Accuracy\n") +
  theme_light() +
  ylab("Test Accuracy\n") +
  xlab("\nUnique Hyperparameter Combinations") +
  theme(
    text = element_text(size = 14, family = "Calibri"),
    legend.title = element_text(size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.text = element_text(size = 14)
  )
dev.off()

# higher learning rate: more fluctuations
# here: max fluctuation is approximately 4%
# my results: max fluctuation is approximately 17%
# not huge inherent problem with the architecture itself

################################################################################