library(dplyr)
library(keras)
library(caret)
source("utils.R")

# Load normalized traces
load("files/ExpConditionsAndTraceMetrics.Rda")
load("files/NormalizedTraces.Rda")


# Load the trained model
classifier_savename <- 'model/classifier_model.hdf5'
classifier_weights <-'model/classifier_weights.hdf5'

classifier <- keras::load_model_hdf5(filepath = classifier_savename)
classifier %>% keras::load_model_weights_hdf5(filepath = classifier_weights,skip_mismatch = TRUE)

# Classify the whole dataset
x <- as.matrix(t(normalized_traces[,KIC.data$UUID]))
x <- reshape_X_3d(x)
y <- KIC.data %>% select(NonArrhythmic,Arrhythmia,Cessation)%>% as.matrix()
labels <- colnames(y)[max.col(y, ties.method = "first")] %>% as.factor()

prediction.probs <- predict(object = classifier,x = x)
prediction.classes <- apply(prediction.probs,1,which.max)
prediction.classes <- case_when(prediction.classes == 1 ~ "Arrhythmic",
                                prediction.classes == 2 ~ "Cessation",
                                prediction.classes == 3 ~ "NonArrhythmic") %>% as.factor()

# Store predictions
prediction_df <- data.frame(UUID = KIC.data$UUID)
prediction_df <- cbind(prediction_df,labels,prediction.classes,prediction.probs)
colnames(prediction_df) <- c("UUID","Classification.Human","Classification.AI",
                             "Prob.Arrhythmic.AI","Prob.Cessation.AI","Prob.NonArrhythmic.AI")

# Add column from what dataset it came from (training, or test, NA if neither)
training.UUIDs <- read.csv("files/training_UUIDs.csv")[,1]
test1.UUIDs <- read.csv("files/test1_UUIDs.csv")[,1]
test2.UUIDs <- read.csv("files/test2_UUIDs.csv")[,1]
test3.UUIDs <- read.csv("files/test3_UUIDs.csv")[,1]
prediction_df <- prediction_df %>% mutate(Dataset = case_when(UUID %in% training.UUIDs ~ "Training",
                                             UUID %in% test1.UUIDs ~ "Test1",
                                             UUID %in% test2.UUIDs ~ "Test2",
                                             UUID %in% test3.UUIDs ~ "Test3"))

prediction_df <- left_join(prediction_df,KIC.data%>%select(UUID,BadBatch))
save(prediction_df,file = "files/AI_predictions.Rda")

