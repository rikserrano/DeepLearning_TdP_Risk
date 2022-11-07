library(tidyverse)
library(caret)
library(plotROC)
load("ExpConditionsAndTraceMetrics.Rda")
load("AI_predictions.Rda")
df <- left_join(KIC.data,prediction_df)

df <- df %>%filter(CiPA.Training==1|CiPA.Validation==1, BadBatch==0)

df <- df %>% mutate(APD90c = (meanAPD90/1000)/(meanTimeToNextPeak/1000)^(1/3))

# Compute baseline APD90c (per cell.line, batch, and plate)
df <- df %>% group_by(cell.line,batch,plate) %>% 
  mutate(APD90c.baseline = mean(APD90c[Dose==0],na.rm=T)) %>% 
  ungroup()

# Per dose response :
#    Check if there are any detected arrhythmias (Pred 1) 
#    Compute max repolarization change (Pred4) 
#    Compute repolarization change @ Cmax (Pred7) 
aaa <- df %>% group_by(Compound,cell.line,batch) %>% group_split() 
test <- aaa[[1]]
bbb <- lapply(aaa, function(test){
  # Predictor 1 
  Pred1 <- any(test$Classification.AI=="Arrhythmia")
  
  s <- test %>% group_by(Dose) %>% mutate(ddAPD90c = APD90c/APD90c.baseline)
  
  # Pred4
  Pred4 <- max(s$ddAPD90c,na.rm = T)
  
  # Pred7
  drc <- s %>% filter(Dose!=0)
  cmax <- test$Cmax[1]
  if(length(unique(drc$Dose))>1){
    dum<- approx(drc$Dose,drc$ddAPD90c,cmax,rule=2)
    Pred7 <- dum$y
  }else{
    Pred7 <- NA
  }
  outdf <- data.frame(Compound = test$Compound[1], cell.line = test$cell.line[1], batch = test$batch[1],
                      Pred1 = Pred1, Pred4 = Pred4, Pred7 = Pred7)
  outdf
})

predictors_df <- bind_rows(bbb)
predictors_df <- left_join(predictors_df,df %>% 
                             select(Compound,CiPA.Classification,CiPA.Training,CiPA.Validation)%>% 
                             mutate(CiPA.Classification.Combined = case_when(
                               CiPA.Classification == "High" ~ "High.Intermediate",
                               CiPA.Classification == "Intermediate" ~ "High.Intermediate",
                               CiPA.Classification == "Low" ~ "Low"
                             )) %>% distinct())

# Model 1 with cross-validation
valid_df <- predictors_df %>% drop_na() %>% filter(grepl(pattern = "HD",x=cell.line))

fit.control <- trainControl(method = "cv", number = 10,
                            summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(123)
fit <- train(CiPA.Classification.Combined ~ Pred1 + Pred4 + Pred7, 
             data = valid_df, method = "glm", 
             family = "binomial", trControl = fit.control)

fit
p <- predict(fit,type = "prob",newdata = valid_df)
df2 <- data.frame(Compound = valid_df$Compound, cell.line = valid_df$cell.line,
                  batch = valid_df$batch,Prob.Low = p$Low,
                  CiPA.Classification.Combined = as.factor(valid_df$CiPA.Classification.Combined) 
) %>% mutate(
  Model.Prediction = as.factor(case_when(Prob.Low <= 0.5 ~ 'High.Intermediate',
                                         Prob.Low > 0.5  ~ 'Low'))
)
confusionMatrix(data=df2$Model.Prediction,reference = df2$CiPA.Classification.Combined)
plot1 <- ggplot(df2, aes(d = CiPA.Classification.Combined, m = Prob.Low)) + 
  geom_roc(labels = F,pointsize = 0.4,size=1)
plot1

auc<- calc_auc(plot1)
plot1+
  annotate("text",x=0.5,0.5,label=paste('AUC = ',round(auc$AUC,digits=2)))+
  scale_x_continuous(name='False positive fraction')+
  scale_y_continuous(name='True positive fraction')+
  theme_classic()+
  theme(aspect.ratio = 1)+
  ggtitle('High-intermediate vs low')
ggsave("Figures/SupplementaryFigure4/ROC_high-intermediate_vs_low_humanmetrics.pdf",height = 2,width = 2)

