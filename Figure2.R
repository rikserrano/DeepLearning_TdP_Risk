library(tidyverse)
library(scales)
library(reshape2)
library(gridExtra)
library(caret)
source("utils.R")

dir.create("Figures")
dir.create("Figures/Figure2")

load("files/NormalizedTraces.Rda")
load("files/ExperimentalConditionsAndTraceMetrics.Rda")
load("files/AI_predictions.Rda")
df <- left_join(KIC.df,prediction_df)
sample.UUIDs <- c("20181218_CIPA_15S1_WT_plate1_G20", # Example of Non-arrhythmic
                  "20181218_CIPA_15S1_WT_plate1_G13", # Example of Arrhythmic
                  "20181218_CIPA_15S1_WT_plate1_G05") # Example of Cessation

sample.traces.df <- df %>% filter(UUID %in% sample.UUIDs)
sample.traces <- normalized_traces %>% select(all_of(sample.UUIDs)) %>%
  mutate(frame = 1:330) %>% melt(id.vars="frame",variable.name = "UUID") %>%
  mutate(UUID = factor(UUID))

plot.trace(sample.traces,frame,value) +
  facet_wrap(~UUID,nrow=3)+theme_traces+
  theme(aspect.ratio = 0.5)
ggsave("Figures/Figure2/sample.traces.pdf",width = 2,height = 2.5)

sample.traces.df <- sample.traces.df %>% mutate(UUID = factor(UUID,levels=sample.UUIDs)) %>%
  select(UUID,Prob.NonArrhythmic.AI,Prob.Arrhythmic.AI,Prob.Cessation.AI) %>%
  melt(id.vars = "UUID")
ggplot(sample.traces.df,aes(x=variable,y=value,fill=variable))+
  geom_col(width=0.7)+
  geom_text(aes(label=sprintf("%0.2f",round(value,digits = 2)),vjust=-0.3))+
  facet_wrap(~UUID,nrow=3,scales="free_x")+
  scale_y_continuous(limits = c(0,1.2),expand = c(0,0))+
  theme_classic(base_family = "sans",
                base_size = 10)+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_text(hjust =0.5,size=9),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  scale_x_discrete(labels=c("Non-Arr.","Arrhythmic","Asystolic"))
ggsave("Figures/Figure2/sample.traces.metrics.pdf",width = 2.2,height = 3)

df2 <- df %>% filter(BadBatch==0) %>% mutate(reference.class = case_when(
  (Dose >= Cmax)&(CiPA.Classification == "High") & (Cessation == 0) ~ "Arrhythmic",
  (Dose <= Cmax)&(CiPA.Classification == "Low") & (Cessation == 0)  ~ "NonArrhythmic",
  (Cessation == 1) ~ "Asystolic",
                  TRUE ~ "N/A"))
training.df <- read.csv("files/training_UUIDs.csv")
validation.df <- read.csv("files/validation_UUIDs.csv")

# Confusion matrix for Training set
dum <- df2 %>% filter(Dataset == "Training")
confusionMatrix(reference = as.factor(dum$reference.class),as.factor(dum$Classification.AI))


# Confusion matrix for Validation set
dum <- df2 %>% filter(Dataset == "Validation")
confusionMatrix(reference = as.factor(dum$reference.class),as.factor(dum$Classification.AI))


# Confusion matrix for Test
dum <- df2 %>% filter(Dataset == "Test", reference.class != "N/A")
confusionMatrix(reference = as.factor(dum$reference.class),as.factor(dum$Classification.AI))
