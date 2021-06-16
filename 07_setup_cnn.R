################################################################################
################################################################################
# SETUP CNN
################################################################################
################################################################################

# load packages
library(keras)

# visualizations & links to papers
# https://towardsdatascience.com/illustrated-10-cnn-architectures-95d78ace614d

# set parameters
n_input_layers <- 5  # TODO: collinearity checks first
n_band_selector <- 3
n_classes <- 7  # or 5 if without rock & grass

# TODO: dropout layers
# anscheinend nach layer_dense() gut

# TODO: batch normalization
# anscheinend nach layer_conv_2d() / layer_dense() gut

# TODO: image augmentation
# https://www.kaggle.com/curiousprogrammer/data-augmentation-in-python-tf-keras-imgaug

# TODO: input pipeline
# https://towardsdatascience.com/implementing-alexnet-cnn-architecture-using-tensorflow-2-0-and-keras-2113e090ad98

# Hyperparameters: number of bands for band selector, learning rate,
#                  number of filters, dropout rate, batch size
# callbacks: number of epochs, (learning rate)
# fixed: tile size, resolution
# https://towardsdatascience.com/what-are-hyperparameters-and-how-to-tune-the-hyperparameters-in-a-deep-neural-network-d0604917584a
# https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/

################################################################################
# LeNet-5
# http://yann.lecun.com/exdb/publis/pdf/lecun-98.pdf
# https://www.kaggle.com/curiousprogrammer/lenet-5-cnn-with-keras-99-48 (based on this code)
################################################################################

lenet5 <- keras_model_sequential() %>%
  # band selector
  layer_conv_2d(input_shape = c(50,50,n_input_layers), filters = n_band_selector,
                kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
  # main model
  layer_conv_2d(filters = 32, kernel_size = c(5,5), activation = "relu", padding="same") %>%
    # normalization?
  layer_max_pooling_2d() %>%
  layer_conv_2d(filters = 48, kernel_size = c(5,5), activation = "relu", padding="valid") %>%
    # normalization?
  layer_max_pooling_2d() %>%
  layer_flatten() %>%
  layer_dense(units = 256, activation = "relu") %>%
    # dropout?
  layer_dense(units = 84, activation = "relu") %>%
    # droput?
  layer_dense(units = n_classes, activation = "softmax")
  
################################################################################
# AlexNet
# https://papers.nips.cc/paper/2012/file/c399862d3b9d6b76c8436e924a68c45b-Paper.pdf
# https://towardsdatascience.com/implementing-alexnet-cnn-architecture-using-tensorflow-2-0-and-keras-2113e090ad98 (based on this code)
# https://github.com/r-tensorflow/alexnet/blob/master/R/alexnet_train.R (and this code)
################################################################################

# hat eigentlich viel größeren Input (224x224x3) & mehr Klassen (1000)
# im paper ist max_pooling mit pool_size = 3 und stride = 2, in dem github code nicht
# je nach code quelle ist padding = "same" / "valid"

alexnet <- keras_model_sequential() %>%
  # band selector
  layer_conv_2d(input_shape = c(50,50,n_input_layers), filters = n_band_selector,
                kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
  # main model
  # block 1
  layer_conv_2d(filters = 96, kernel_size = c(11,11), activation = "relu", strides = c(4,4)) %>%
  layer_batch_normalization() %>%
  layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
  # block 2
  layer_conv_2d(filters = 256, kernel_size = c(5,5), activation = "relu", strides = c(1,1), padding="same") %>%
  layer_batch_normalization() %>%
  layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
  # block 4
  layer_conv_2d(filters = 384, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
  layer_batch_normalization() %>%
  layer_conv_2d(filters = 384, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
  layer_batch_normalization() %>%
  layer_conv_2d(filters = 256, kernel_size = c(3,3), activation = "relu", strides = c(1,1), padding="same") %>%
  layer_batch_normalization() %>%
  layer_max_pooling_2d(pool_size=c(3,3), strides=c(2,2)) %>%
  # block 5
  layer_flatten() %>%
  layer_dense(units = 4096, activation = "relu") %>%
  layer_dropout(0.5) %>%
  layer_dense(units = 4096, activation = "relu") %>%
  layer_dropout(0.5) %>%
  layer_dense(units = n_classes, activation = "softmax")

################################################################################
# VGG-16
# https://arxiv.org/pdf/1409.1556.pdf
# https://towardsdatascience.com/step-by-step-vgg16-implementation-in-keras-for-beginners-a833c686ae6c (based on this code)
################################################################################

vgg16 <- keras_model_sequential() %>%
  # band selector
  layer_conv_2d(input_shape = c(50,50,n_input_layers), filters = n_band_selector,
                kernel_size = c(1,1), strides = c(1,1), padding="same") %>%
  # main model
  # block 1
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_max_pooling_2d() %>% 
  # block 2
  layer_conv_2d(filters = 128, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 128, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_max_pooling_2d() %>%
  # block 3
  layer_conv_2d(filters = 256, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 256, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 256, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_max_pooling_2d() %>%
  # block 4
  layer_conv_2d(filters = 512, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 512, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 512, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_max_pooling_2d() %>%
  # block 5
  layer_conv_2d(filters = 512, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 512, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_conv_2d(filters = 512, kernel_size = c(3,3), activation = "relu", padding = "same") %>%
    # normalization?
  layer_max_pooling_2d() %>%
  # block 6
  layer_flatten() %>%
  layer_dense(units = 4096, activation = "relu") %>%
    # dropout?
  layer_dense(units = 4096, activation = "relu") %>%
    # dropout?
  layer_dense(units = n_classes, activation = "softmax")
  
################################################################################
# LOADING & AUGMENTING DATA
################################################################################  


################################################################################  