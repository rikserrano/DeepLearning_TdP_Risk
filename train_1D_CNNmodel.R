library(tidyverse)
library(keras)
library(caret)
source("utils.R")

load("files/ExperimentalConditionsAndTraceMetrics.Rda")
load("files/Traces.Rda")
training.UUIDs <- read.csv("files/training_UUIDs.csv")[,1]
training.labels <- read.csv("files/training_UUIDs.csv")[,2]
validation.UUIDs <- read.csv("files/validation_UUIDs.csv")[,1]
validation.labels <- read.csv("files/validation_UUIDs.csv")[,2]


# Normalize each trace to range [0,1] per plate
aaa <- KIC.df %>% group_by(cell.line,batch,plate) %>% group_split()
plate.info <- aaa[[1]]
bbb <- lapply(aaa, function(plate.info){
  plate.traces <- traces.data %>% select(plate.info$UUID)
  plate.1p <- quantile(plate.traces%>%unlist(),0.01)
  plate.99p <- quantile(plate.traces%>%unlist(),0.99)
  normalized <- (plate.traces-plate.1p)/(plate.99p-plate.1p)
})
normalized_traces <- bind_cols(bbb)

save(normalized_traces,file="files/NormalizedTraces.Rda")

# Prepare the training data
x_train <- as.matrix(t(normalized_traces[,training.UUIDs]))
x_train <- reshape_X_3d(x_train)
y_train <- labels_to_one_hot(training.labels)

# Prepare the validation data
x_valid <- as.matrix(t(normalized_traces[,validation.UUIDs]))
x_valid <- reshape_X_3d(x_valid)
y_valid <- labels_to_one_hot(validation.labels)

# CNN network implementation in Keras
input_layer <- 
  layer_input(shape = c(330,1))
cnn <- 
  input_layer %>% 
  layer_conv_1d(
    filters = 33,
    kernel_size = 5,
    padding = "same",
    activation = "relu"
  ) %>%
  layer_conv_1d(
    filters = 33,
    kernel_size = 5,
    padding = "same",
    activation = "relu"
  ) %>%
  layer_max_pooling_1d(3, padding = "same")%>%
  layer_conv_1d(
    filters = 66,
    kernel_size = 10,
    padding = "same",
    activation = "relu"
  )%>%
  layer_conv_1d(
    filters = 66,
    kernel_size = 10,
    padding = "same",
    activation = "relu"
  )%>%
  layer_global_average_pooling_1d()

classifier <- cnn %>% 
  layer_dropout(0.2)%>%
  layer_dense(units = 3, activation = 'softmax')
classifier <- keras_model(input_layer,classifier)
summary(classifier)

classifier  %>% compile(
  loss = "categorical_crossentropy",
  optimizer = "adam",
  metrics = list("acc")
)

classifier %>% fit(
  x = x_train,
  y = y_train,
  validation_data = list(x_valid,y_valid),
  epochs = 20
)

# Save the trained network
savename_classifier <-'model/classifier_model.hdf5'
savename_classifier_weights <-'model/classifier_weights.hdf5'
savename_encoder <- 'model/encoder_model.hdf5'
savename_encoder_weights <-'model/encoder_weights.hdf5'

classifier_weights <- 
  classifier %>%
  keras::get_weights()

encoder <- keras_model(input_layer,cnn)

encoder_weights <- 
  encoder %>% 
  keras::get_weights()

save_model_hdf5(object = classifier, filepath = savename_classifier)
save_model_weights_hdf5(object = classifier,filepath = savename_classifier_weights,overwrite = TRUE)
save_model_hdf5(object = encoder, filepath = savename_encoder)
save_model_weights_hdf5(object = encoder,filepath = savename_encoder_weights,overwrite = TRUE)





