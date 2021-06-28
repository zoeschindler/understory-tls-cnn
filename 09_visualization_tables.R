################################################################################
################################################################################
# VISUALIZATION & TABLES
################################################################################
################################################################################

# load packages
library(ggplot)
library(grid)

# set paths

################################################################################
# RAW POINT CLOUDS
################################################################################

# maybe better in CC

################################################################################
# CLUSTERS & PCA
################################################################################

# cluster before & after

# PCA after (before too cluttered)

################################################################################
# RASTER STATISTICS & LABELS
################################################################################

# load & transform plots
# get all raster paths
# loop though raster types
# loop through rasters
# loop through unique labels
# clip raster + plots with label
# save values

# boxplots for some rasters

# i want radarplots but they are not popular

################################################################################
# HYPERPARAMETERS & ACCURACY & LOSS
################################################################################

# correlation plot - hyperparameters & accuracy & loss

# barplots (with CI) / boxplots - hyperparameters & accuracy & loss

# regression between val_acc and val_loss
# (good if val_acc was used for picking network)

################################################################################
# CONFUSION MATRICES
################################################################################

# does not need significance levels

################################################################################
# ACCURACY COMPARISON INPUTS
################################################################################

# boxplots

# table, including sd

################################################################################
# VISUALIZE FILTERS
################################################################################

# Funktion zum Erzeugen von Filtervisualisierungen
generate_pattern <- function(layer_name, filter_index, size = 150) {  # size = Zoom
  layer_output <- model$get_layer(layer_name)$output
  loss <- k_mean(layer_output[,,,filter_index])
  grads <- k_gradients(loss, model$input)[[1]]
  grads <- grads / (k_sqrt(k_mean(k_square(grads))) + 1e-5)
  iterate <- k_function(list(model$input), list(loss, grads))
  input_img_data <- array(runif(size * size * 3),
                          dim = c(1, size, size, 3)) * 20 + 128
  step <- 1
  for (i in 1:40) {
    c(loss_value, grads_value) %<-% iterate(list(input_img_data))
    input_img_data <- input_img_data + (grads_value * step)
  }
  img <- input_img_data[1,,,]
  deprocess_image(img)
}

# Muster darstellen
grid.raster(generate_pattern("block3_conv1", 1))

################################################################################