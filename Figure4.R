library(tidyverse)
library(scales)
library(reshape2)
library(gridExtra)
library(ggstance)
source("utils.R")
load("files/NormalizedTraces.Rda")
load("files/ExpConditionsAndTraceMetrics.Rda")
load("files/AI_predictions.Rda")
load("files/drcs.arrhythmic.Rda") # generated at compute_drcs_ec50s
load("files/ec50s.arrhythmic.Rda") # generated at compute_drcs_ec50s
load("files/drcs.cessation.Rda")  # generated at compute_drcs_ec50s
load("files/ec50s.cessation.Rda") # generated at compute_drcs_ec50s

df <- left_join(KIC.data,prediction_df)
# Proarrhythmic safety margin for mutants + isogenic control --------------
safety.margin.hp <- ec50s.arrhythmic %>% filter(cell.line %in% c("HD.15S1","DCM.PLN","DCM.TNNT2","DCM.RBM20","HCM.MYBPC3","HCM.TPM1")) %>%
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
  scale_fill_gradientn(name = "Torsadogenic\nSafety Margin = \nlog10(ec50/Cmax)",
                       colours=colorRampPalette(c("#FF0000","#44BCD8","#CCE7E8","#FFFFFF"))(100),
                       limits=c(0,5), oob = squish,na.value = "#FFFFFF")+
  scale_y_discrete(limits=rev,expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  theme_classic(base_size = 11,
                base_family = "sans")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color=rev(labcolor$color)),
        aspect.ratio = 5)
ggsave(paste0("Figures/Figure4/SafetyMargin_mutants.pdf"),height = 5,width=5)

# Cessation safety margin for mutants + isogenic control ------------------
safety.margin.hp <- ec50s.cessation %>% filter(cell.line %in% c("HD.15S1","DCM.PLN","DCM.TNNT2","DCM.RBM20","HCM.MYBPC3","HCM.TPM1")) %>%
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
  scale_fill_gradientn(name = "Asystolic\nSafety Margin = \nlog10(ec50/Cmax)",
                       colours=colorRampPalette(c("#111111","#FFFFFF"))(100),
                       limits=c(0,5), oob = squish,na.value = "#FFFFFF")+
  scale_y_discrete(limits=rev,expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  theme_classic(base_size = 11,
                base_family = "sans")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(color=rev(labcolor$color)),
        aspect.ratio = 5)
ggsave(paste0("Figures/Figure4/CessationSafetyMargin_mutants.pdf"),height = 5,width=5)