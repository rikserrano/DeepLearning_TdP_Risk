library(dplyr)
library(keras)
library(caret)
source("utils.R")

# Load normalized traces
load("files/ExperimentalConditionsAndTraceMetrics.Rda")
load("files/NormalizedTraces.Rda")


# Load the trained model
classifier_savename <- 'model/classifier_model.hdf5'
classifier_weights <-'model/classifier_weights.hdf5'

classifier <- keras::load_model_hdf5(filepath = classifier_savename)
classifier %>% keras::load_model_weights_hdf5(filepath = classifier_weights,skip_mismatch = TRUE)

# Classify the whole dataset
x <- as.matrix(t(normalized_traces[,KIC.df$UUID]))
x <- reshape_X_3d(x)

prediction.probs <- predict(object = classifier,x = x)
prediction.classes <- apply(prediction.probs,1,which.max)
prediction.classes <- case_when(prediction.classes == 1 ~ "Arrhythmic",
                                prediction.classes == 2 ~ "Asystolic",
                                prediction.classes == 3 ~ "NonArrhythmic") %>% as.factor()

# Store predictions
prediction_df <- data.frame(UUID = KIC.df$UUID)
prediction_df <- cbind(prediction_df,prediction.classes,prediction.probs)
colnames(prediction_df) <- c("UUID","Classification.AI",
                             "Prob.Arrhythmic.AI","Prob.Asystolic.AI","Prob.NonArrhythmic.AI")

# Add column from what dataset it came from (training, validation, or test)
training.UUIDs <- read.csv("files/training_UUIDs.csv")[,1]
validation.UUIDs <- read.csv("files/validation_UUIDs.csv")[,1]
test.UUIDs <- read.csv("files/test_UUIDs.csv")[,1]
prediction_df <- prediction_df %>% mutate(Dataset = case_when(UUID %in% training.UUIDs ~ "Training",
                                             UUID %in% validation.UUIDs ~ "Validation",
                                             UUID %in% test.UUIDs ~ "Test"))

prediction_df <- left_join(prediction_df,KIC.df%>%select(UUID,BadBatch))
save(prediction_df,file = "files/AI_predictions.Rda")
