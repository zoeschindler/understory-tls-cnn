# Masterarbeit 

* <a href ="https://github.com/zoeschindler/masterarbeit/blob/main/01_preprocessing.R">01_preprocessing.R</a>: filter, clip, merge single scans ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/02_prepare_vegetation.R">02_prepare_vegetation.R</a>: convert vegetation data from json to kml, clean data, download images ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/03_raster_calculation_LAScatalog.R">03_raster_calculation_LAScatalog.R</a>: normalize, filter understory, calculate rasters ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/03_raster_calculation_functions.R">03_raster_calculation_functions.R</a>: functions for normalizing, filtering, raster calculation ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/04_check_collinearity.R">04_check_collinearity.R</a>: check & remove raster collinearities ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/05_combine_rasters_vegetation.R">05_combine_rasters_vegetation.R</a>: filter vegetation plots, clip rasters to vegetation plots ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/06_prepare_cnn_input.R">06_prepare_cnn_input.R</a>: stack, scale, save raster clips ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/07_setup_cnn.R">07_setup_cnn.R</a>: setup CNN architecture, set up data input pipeline ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/08_train_cnn.R">08_train_cnn.R</a>: load data, train CNNs, hyperparameter optimization ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/08_experiment.R">08_experiment.R</a>: for running `tfruns` & hyperparameter optimization ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/09_filter_plot.R">09_filter_plot.R</a>: visualization of the understory filtering ✓<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/10_prediction_map.R">10_prediction_map.R</a>: visualize prediction results on one area as raster (in progress)<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/11_prediction_points.R">11_prediction_points.R</a>: add prediction results to point cloud & visualize point cloud (in progress)<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/12_visualization_tables.R">12_visualization_tables.R</a>: create result visualizations (in progress)<br>
* <a href = "https://github.com/zoeschindler/masterarbeit/blob/main/13_test_architecture_cifar10.R">13_test_architecture_cifar10.R</a>: test architecture using CIFAR10 dataset, check for accuracy fluctuations ✓<br>

## Workflow: Data Preparation

<img align="center" src="https://github.com/zoeschindler/masterarbeit/blob/main/Visualisierung_Workflow_1.png">

## Workflow: Cross-Validation

<img align="center" src="https://github.com/zoeschindler/masterarbeit/blob/main/Visualisierung_Workflow_2.png">

## CNN: Architecture

<img align="center" src="https://github.com/zoeschindler/masterarbeit/blob/main/Visualisierung_Architektur.png">
