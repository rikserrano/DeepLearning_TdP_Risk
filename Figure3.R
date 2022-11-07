library(tidyverse)
library(scales)
library(reshape2)
library(gridExtra)
source("utils.R")
load("files/NormalizedTraces.Rda")
load("files/ExpConditionsAndTraceMetrics.Rda")
load("files/AI_predictions.Rda")
load("files/drcs.arrhythmic.Rda")
load("files/ec50s.arrhythmic.Rda")
df <- left_join(KIC.data,prediction_df)

# Ibutilide plot (panels A and B) -------------------------------------------
# Plot traces
sample.drc <- df %>% filter(cell.line == "HD.15S1",batch == 2, Compound == "Ibutilide", Well.letter == "I")
sample.drc.traces <- normalized_traces %>% select(sample.drc$UUID)%>% mutate(frame = 1:330) %>% 
  melt(id.vars="frame",variable.name = "UUID") %>%
  mutate(UUID = fct_rev(factor(UUID)))
plot.trace(sample.drc.traces,frame,value) +
  facet_wrap(~UUID,nrow=2)+theme_traces+
  theme(aspect.ratio = 0.5)
ggsave('Figures/Figure3/Ibutilide_traces.pdf',width = 6,height = 2.5)
# Fit logistic curve to the drc of arrhythmia prob 
a <- fit.binomial(sample.drc %>% filter(Dose>0), var.name = "Prob.Arrhythmic.AI")
ggplot(sample.drc,aes(x=log10(1e-6*Dose),y=Prob.Arrhythmic.AI))+
  geom_vline(xintercept = log10(1e-6*sample.drc$Cmax[1]),linetype="dashed")+
  geom_point(size=0.75) +
  geom_line(data=a$fit.fitted,aes(x=log10(1e-6*Dose),y=fit_resp),color="black")+
  geom_point(data=a$fit.ec50,aes(x=log10(1e-6*ec50), y=1.05),shape=25,size=2,fill="black")+
  geom_segment(aes(x = log10(1e-6*a$fit.ec50$ec50) , y = 0.5, xend = log10(1e-6*a$fit.ec50$ec50), yend = 1.05),linetype="dotted")+
  scale_y_continuous(name="Probability\nArrhythmia")+
  scale_x_continuous(name="log10[ibutilide]")+
  facet_wrap(~Compound)+
  theme_classic(base_family = "sans",
                base_size = 11)+
  theme(axis.text = element_text(color="black"),
        aspect.ratio = 0.75)
ggsave('Figures/Figure3/Ibutilide_drc_arrhythmia.pdf',width = 3,height = 2)

# Nifedipine plots (panels C and D) ---------------------------------------
sample.drc <- df %>% filter(cell.line == "HD.15S1",batch == 2, Compound == "Nifedipine", Well.letter == "J")
sample.drc.traces <- normalized_traces %>% select(sample.drc$UUID)%>% mutate(frame = 1:330) %>% 
  melt(id.vars="frame",variable.name = "UUID") %>%
  mutate(UUID = fct_rev(factor(UUID)))
plot.trace(sample.drc.traces,frame,value) +
  facet_wrap(~UUID,nrow=2)+theme_traces+
  theme(aspect.ratio = 0.5)
ggsave('Figures/Figure3/Nifedipine_traces.pdf',width = 6,height = 2.5)

a <- fit.binomial(sample.drc %>% filter(Dose>0), var.name = "Prob.Cessation.AI")
ggplot(sample.drc,aes(x=log10(1e-6*Dose),y=Prob.Cessation.AI))+
  geom_vline(xintercept = log10(1e-6*sample.drc$Cmax[1]),linetype="dashed")+
  geom_point(size=0.75) +
  geom_line(data=a$fit.fitted,aes(x=log10(1e-6*Dose),y=fit_resp),color="black")+
  geom_point(data=a$fit.ec50,aes(x=log10(1e-6*ec50), y=1.05),shape=25,size=2,fill="black")+
  geom_segment(aes(x = log10(1e-6*a$fit.ec50$ec50) , y = 0.5, xend = log10(1e-6*a$fit.ec50$ec50), yend = 1.05),linetype="dotted")+
  scale_y_continuous(name="Probability\nAsystole")+
  scale_x_continuous(name="log10[nifedipine]")+
  facet_wrap(~Compound)+
  theme_classic(base_family = "sans",
                base_size = 11)+
  theme(axis.text = element_text(color="black"),
        aspect.ratio = 0.75)
ggsave('Figures/Figure3/Nifedipine_drc_asystole.pdf',width = 3,height = 2)

ggplot(sample.drc,aes(x=log10(1e-6*Dose),y=Prob.Arrhythmic.AI))+
  geom_vline(xintercept = log10(1e-6*sample.drc$Cmax[1]),linetype="dashed")+
  geom_point(size=0.75) +
  geom_line(data=a$fit.fitted,aes(x=log10(1e-6*Dose),y=fit_resp),color="black")+
  geom_point(data=a$fit.ec50,aes(x=log10(1e-6*ec50), y=1.05),shape=25,size=2,fill="black")+
  geom_segment(aes(x = log10(1e-6*a$fit.ec50$ec50) , y = 0.5, xend = log10(1e-6*a$fit.ec50$ec50), yend = 1.05),linetype="dotted")+
  scale_y_continuous(name="Probability\nArrhythmic")+
  scale_x_continuous(name="log(Dose)")+
  facet_wrap(~Compound)+
  theme_classic(base_family = "sans",
                base_size = 11)+
  theme(axis.text = element_text(color="black"),
        aspect.ratio = 0.75)
ggsave('Figures/Figure3/Ibutilide_drc.pdf',width = 3,height = 2)

# Plot other sample Prob.Arrhythmic  -------------------------------------------
sample.comps <- c("Ibutilide","Flecainide","Ondansetron","Citalopram","Loratadine","Aspirin")
drc <- drcs.arrhythmic.df %>% filter(Compound %in% sample.comps, cell.line == "HD.15S1") %>%
  mutate(Compound = factor(Compound,levels = sample.comps))
drc.ec50 <- ec50s.arrhythmic %>% filter(Compound %in% sample.comps, cell.line == "HD.15S1") %>%
  mutate(Compound = factor(Compound,levels = sample.comps))

ggplot(drc,aes(x=log10(1e-6*Dose),y=fit_resp))+
  geom_vline(aes(xintercept = log10(1e-6*Cmax)),linetype="dashed")+
  geom_line()+
  geom_point(data=drc.ec50,aes(x=log10(1e-6*ec50), y=1.05),shape=25,size=2,fill="black")+
  geom_segment(data=drc.ec50,aes(x = log10(1e-6*ec50) , y = c(0.5,0.5,0.5,0.5,0.5,0.5), xend = log10(1e-6*ec50), yend = c(1.05,1.05,1.05,1.05,1.05,1.05)),linetype="dotted")+
  facet_wrap(~Compound,scales = "free_x",nrow=1)+
  scale_y_continuous(name="Probability\nArrhythmia")+
  scale_x_continuous(name="log10[concentration]")+
  theme_classic(base_family = "sans",
                base_size = 11)+
  theme(axis.text = element_text(color="black"),
        aspect.ratio = 0.75)
ggsave('Figures/Figure3/Six_drc.pdf',width = 12,height = 2)


# Compute Torsadogenic Safety Margin --------------------------------------
safety.margin.hp <- ec50s.arrhythmic %>% filter(grepl(pattern = "HD",x=cell.line)) %>%
  mutate(Compound = fct_reorder(Compound,.fun=mean,safety.margin,.desc = F,na.rm=T))
labcolor <- safety.margin.hp %>% select(Compound,CiPA.Classification) %>% distinct() %>%
  arrange(Compound) %>%
  mutate(color= case_when(CiPA.Classification == "High" ~ "#FF0000",
                          CiPA.Classification == "Intermediate" ~ "#004FFF",
                          CiPA.Classification == "Low" ~ "#00A0BA",
                          CiPA.Classification == "Other" ~ "#545454",
                          TRUE ~ "black"))
ggplot(safety.margin.hp,aes(x=cell.line,y=Compound,fill=safety.margin))+
  geom_tile(color="white")+
  scale_fill_gradientn(name = "Torsadogenic\nSafety Margin = \n log10(ec50/Cmax)",
                       colours=colorRampPalette(c("#FF0000","#44BCD8","#CCE7E8","#FFFFFF"))(100),
                       limits=c(-1,5), oob = squish,na.value = "#FFFFFF")+
  scale_y_discrete(limits=rev,expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  theme_classic(base_size = 11,
                base_family = "sans")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color=rev(labcolor$color)),
        aspect.ratio = 7)
ggsave(paste0("Figures/Figure3/SafetyMargin.pdf"),height = 5,width=5)


# Plot drc in colorbar form -----------------------------------------------
drcs.arrhythmic.df$Compound = factor(drcs.arrhythmic.df$Compound,levels = rev(levels(safety.margin.hp$Compound)))
lines <- unique(drcs.arrhythmic.df$cell.line)
for (iline in lines){
  plotdf <- drcs.arrhythmic.df %>% filter(cell.line == iline)
  ec50df <- ec50s.arrhythmic %>% filter(cell.line==iline)
  ggplot(plotdf ,aes(x= log10(Dose/Cmax),y=Compound,fill=fit_resp))+
    geom_tile(width=log10(3),height=0.5)+
    geom_point(data=ec50df,aes(x=log10(ec50/Cmax), y=Compound),position = position_nudge(y=0.5),
               shape=25,size=1,fill="black")+
    geom_vline(xintercept = 0,linetype="dashed")+
    ggtitle(iline)+
    scale_fill_gradientn(name="Probability Arrhythmia", breaks = c(seq(from=0,to=1,length.out = 6)),
                         colours=colorRampPalette(c("#006400","#008000","#FFA500","#FF8C00","#FF4500","#8B0000"))(100),
                         limits=c(0,1), oob = squish)+
    scale_x_continuous(limits=c(-3,4))+
    theme_classic(base_size = 11,
                  base_family = "sans")+
    theme(aspect.ratio = 3,
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(color = "black"))
  ggsave(paste0("Figures/Figure3/",iline,".bars.pdf"),height = 5,width=5)
}

# Predictions using Tosadogenic Safety Margin of Healthy Donors -----------
library(caret)
library(plotROC)
# High-intermediate vs low
risk.combined.df <- ec50s.arrhythmic %>% filter(grepl(pattern = "HD",x=cell.line),
                                                CiPA.Classification %in% c("High","Intermediate","Low")) %>%
  replace_na(list(safety.margin=7)) %>%
  mutate(CiPA.Classification.Combined = case_when(
    CiPA.Classification == "High" ~ "High",
    CiPA.Classification == "Intermediate" ~ "Low",
    CiPA.Classification == "Low" ~ "Low"),
    CiPA.Classification.Combined2 = case_when(
      CiPA.Classification == "High" ~ "High",
      CiPA.Classification == "Intermediate" ~ "High",
      CiPA.Classification == "Low" ~ "Low")
  )


fit.control <- trainControl(method = "cv", number = 10,
                            summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(123)
fit <- train(CiPA.Classification.Combined ~ safety.margin, 
             data = risk.combined.df, method = "glm", 
             family = "binomial", trControl = fit.control)
fit
p <- predict(fit,type = "prob",newdata = risk.combined.df)
p2 <- predict(fit,type = "raw",newdata = risk.combined.df)

df <- data.frame(Compound = risk.combined.df$Compound, cell.line = risk.combined.df$cell.line,
                 Prob.Low = p$Low,
                 Model.Prediction = p2,
                 CiPA.Classification.Combined = as.factor(risk.combined.df$CiPA.Classification.Combined) 
)
confusionMatrix(data=df$Model.Prediction,reference = df$CiPA.Classification.Combined)
plot1 <- ggplot(df, aes(d = CiPA.Classification.Combined, m = Prob.Low)) + 
  geom_roc(labels = F,pointsize = 0.4,size=1)
auc<- calc_auc(plot1)
plot1+
  annotate("text",x=0.5,0.5,label=paste('AUC = ',round(auc$AUC,digits=2)))+
  scale_x_continuous(name='False positive fraction')+
  scale_y_continuous(name='True positive fraction')+
  theme_classic()+
  theme(aspect.ratio = 1)+
  ggtitle('High-intermediate vs Low')

ggsave("Figures/Figure3/ROC_combined.pdf",height = 2,width = 2)
a <- df %>% mutate(wrong = CiPA.Classification.Combined!=Model.Prediction) %>% 
  group_by(Compound) %>% mutate(perc = sum(wrong)/n())

risk.high.low.df <- ec50s.arrhythmic %>% filter(grepl(pattern = "HD",x=cell.line),
                                                         CiPA.Classification %in% c("High","Low"))%>%  
  replace_na(list(safety.margin=7))

fit.control <- trainControl(method = "cv", number = 10,
                            summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(123)
fit <- train(CiPA.Classification ~ safety.margin, 
             data = risk.high.low.df, method = "glm", 
             family = "binomial", trControl = fit.control)
fit
p <- predict(fit,type = "prob",newdata = risk.high.low.df)
p2 <- predict(fit,type = "raw",newdata = risk.high.low.df)
df <- data.frame(Compound = risk.high.low.df$Compound, 
                 cell.line = risk.high.low.df$cell.line,
                 Prob.Low = p$Low,
                 Model.Prediction = p2,
                 CiPA.Classification = as.factor(risk.high.low.df$CiPA.Classification) 
)
confusionMatrix(data=df$Model.Prediction,reference = df$CiPA.Classification)
plot2 <- ggplot(df, aes(d = CiPA.Classification, m = Prob.Low)) +
  geom_roc(labels = F,pointsize = 0.4,size=1)
auc<- calc_auc(plot2)
plot2+
  annotate("text",x=0.5,0.5,label=paste('AUC = ',round(auc$AUC,digits=2)))+
  scale_x_continuous(name='False positive fraction')+
  scale_y_continuous(name='True positive fraction')+
  theme_classic()+
  theme(aspect.ratio = 1)+
  ggtitle('High vs low')
ggsave("Figures/Figure3/ROC_high_vs_low.pdf",height = 2,width = 2)


risk.intermediate.low.df <- ec50s.arrhythmic %>% filter(grepl(pattern = "HD",x=cell.line),
                                                         CiPA.Classification %in% c("Intermediate","Low"))%>%  
  replace_na(list(safety.margin=7))
set.seed(123)
fit <- train(CiPA.Classification ~ safety.margin, 
             data = risk.intermediate.low.df, method = "glm", 
             family = "binomial", trControl = fit.control)
fit
p <- predict(fit,type = "prob",newdata = risk.intermediate.low.df)
p2 <- predict(fit,type = "raw",newdata = risk.intermediate.low.df)
df <- data.frame(Compound = risk.intermediate.low.df$Compound, 
                 cell.line = risk.intermediate.low.df$cell.line,
                 Prob.Low = p$Low,
                 Model.Prediction = p2,
                 CiPA.Classification = as.factor(risk.intermediate.low.df$CiPA.Classification) 
)
confusionMatrix(data=df$Model.Prediction,reference = df$CiPA.Classification)
plot3 <- ggplot(df, aes(d = CiPA.Classification, m = Prob.Low)) + 
  geom_roc(labels = F,pointsize = 0.4,size=1)
auc<- calc_auc(plot3)
plot3+
  annotate("text",x=0.5,0.5,label=paste('AUC = ',round(auc$AUC,digits=2)))+
  scale_x_continuous(name='False positive fraction')+
  scale_y_continuous(name='True positive fraction')+
  theme_classic()+
  theme(aspect.ratio = 1)+
  ggtitle('Intermediate vs low')
ggsave("Figures/Figure3/ROC_intermediate_vs_low.pdf",height = 2,width = 2)





