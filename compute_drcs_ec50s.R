library(tidyverse)
library(scales)
library(reshape2)
library(gridExtra)
library(ggstance)
library(ggpattern)
source("utils.R")
load("files/ExpConditionsAndTraceMetrics.Rda") # KIC.data
load("files/AI_predictions.Rda") # prediction_df
df <- left_join(KIC.data,prediction_df)
th <- theme_classic(base_size = 9,
                    base_family = "sans")+
  theme(aspect.ratio = 1,
        strip.background=element_rect(color = NA,  fill=NA),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color="black"))

# Dose response curves and ec50s Prob.Arrhythmic --------
aaa <- df %>% filter(BadBatch==0, !Compound %in% c("Control","GS-967","DMSO",NA)) %>%
  group_by(Compound) %>% group_split()
comp.df <- aaa[[1]]
eee <- lapply(aaa, function(comp.df){
  comp.df <- comp.df %>% filter(Classification.AI != "Cessation", Dose >0)
  ccc <-  comp.df %>% group_by(cell.line)%>% group_split()
  ddd <- lapply(ccc,fit.binomial,var.name="Prob.Arrhythmic.AI")
  fitted.values <- lapply(ddd,'[[',"fit.fitted") %>% bind_rows()
  ec50s <- lapply(ddd,'[[',"fit.ec50") %>% bind_rows()
  l <- list(fitted.values,ec50s)
})
drcs.arrhythmic.df <- lapply(eee, '[[',1) %>% bind_rows()
drcs.arrhythmic.df <- left_join(drcs.arrhythmic.df,df %>% select(Compound,Cmax) %>% distinct())
save(file="files/drcs.arrhythmic.Rda",drcs.arrhythmic.df)

maxdoses <- df %>% group_by(Compound) %>% select(Dose) %>% distinct() %>%   summarize(maxDose = max(Dose))
ec50s.arrhythmic <- lapply(eee, '[[',2) %>% bind_rows()
ec50s.arrhythmic <- left_join(ec50s.arrhythmic,df %>% select(Compound,Cmax,CiPA.Classification) %>% distinct())
ec50s.arrhythmic <- left_join(ec50s.arrhythmic,maxdoses)
ec50s.arrhythmic <- ec50s.arrhythmic %>% mutate(ec50.NA = ifelse(ec50>10*maxDose,NA,ec50),
                                                safety.margin =log10(ec50/Cmax))
save(file="files/ec50s.arrhythmic.Rda",ec50s.arrhythmic)


# plot all drcs for arrhythmia probability --------------------------------

bbb <- drcs.arrhythmic.df %>% group_by(Compound) %>% group_split()
comp.df <- bbb[[13]]
plist <- lapply(bbb, function(comp.df){
  icomp <- comp.df$Compound[1]
  raw.data <- df %>% filter(Compound==icomp, BadBatch==0,Classification.AI!="Cessation")
  ec50.df <- ec50s.arrhythmic %>% filter(Compound == icomp)
  cessation.region.df <- comp.df %>% group_by(cell.line) %>% select(Dose) %>%  
    summarize(max.drc.dose = max(Dose),max.dose = max(ec50.df$maxDose))
  ec50.df <- left_join(ec50.df,cessation.region.df) %>% 
    mutate(ec50 = case_when(ec50>max.drc.dose ~ NA_real_,
                            TRUE ~ ec50))
  cls <- unique(ec50.df$cell.line)
  ggplot(comp.df,aes(x=log10(1e-6*Dose),y=fit_resp,color = cell.line))+
    geom_rect_pattern(data = cessation.region.df,
              mapping = aes(xmin=log10(1e-6*max.drc.dose),xmax=log10(1e-6*max.dose),
                            ymin=0, ymax=1), inherit.aes = F, fill= 'black',
              colour= 'white', alpha = 0.2, pattern_density = 0.8,pattern_fill="gray90",pattern_color="white")+
    geom_line() +
    geom_point(data = raw.data,aes(x=log10(1e-6*Dose),y=Prob.Arrhythmic.AI,
                                   color = cell.line),size=1,alpha=0.5)+
    scale_color_manual(values = cell.line.palette[cls])+
    annotate("line",x= rep(log10(1e-6*ec50.df$Cmax[1]),2),y=c(0,1),linetype="dashed")+
    geom_errorbarh(data=ec50.df,aes(x=log10(1e-6*ec50), xmin=log10(1e-6*ec50.upr),
                                     xmax=log10(1e-6*ec50.lwr), y = 0.5),height=0.1)+
    geom_point(data=ec50.df,aes(x=log10(1e-6*ec50), y = 0.5))+
    facet_wrap(~cell.line,nrow=2,scales="free")+
    scale_y_continuous(limits = c(0,1),expand = c(0,0),
                       name="Probability\nArrhythmia")+
    scale_x_continuous(expand = c(0,0), name = paste0("Log10[",icomp,"]"))+
    coord_cartesian(xlim=c(min(log10(1e-6*comp.df$Dose)),max(log10(1e-6*comp.df$Dose))))+
    th
  
})
ml <- marrangeGrob(plist,ncol = 1,nrow = 1)
ggsave("Figures/all_drcs/all_prob_arrhythmic.pdf",ml,width = 7,height = 4)

# plot all drcs for arrhythmic probability of Mutants vs Isogenic Control --------------------------------

bbb <- drcs.arrhythmic.df %>% filter(!cell.line %in% c("HD.113","HD.273")) %>% 
  group_by(Compound) %>% group_split()
comp.df <- bbb[[26]]
plist.mut <- lapply(bbb, function(comp.df){
  icomp <- comp.df$Compound[1]
  drc.hd <- comp.df %>% filter(cell.line == "HD.15S1")
  drc.mut <- comp.df %>% filter(!grepl(pattern = "HD",x=cell.line))
  ec50.hd <- ec50s.arrhythmic %>% filter(Compound == icomp, cell.line == "HD.15S1")
  ec50.mut <- ec50s.arrhythmic %>% filter(Compound == icomp, !grepl(pattern = "HD",x=cell.line))
  cessation.region.df <- drc.mut %>% group_by(cell.line) %>% select(Dose) %>%  
    summarize(max.drc.dose = max(Dose),max.dose = max(ec50.mut$maxDose))
  ec50.mut <- left_join(ec50.mut,cessation.region.df) %>% 
    mutate(ec50 = case_when(ec50>max.drc.dose ~ NA_real_,
                            TRUE ~ ec50))
  cls <- unique(ec50.mut$cell.line)
  ggplot(drc.mut,aes(x=log10(1e-6*Dose),y=fit_resp,color = cell.line))+
    geom_rect_pattern(data = cessation.region.df,
                      mapping = aes(xmin=log10(1e-6*max.drc.dose),xmax=log10(1e-6*max.dose),
                                    ymin=0, ymax=1), inherit.aes = F, fill= 'black',
                      colour= 'white', alpha = 0.2, pattern_density = 0.8,pattern_fill="gray90",pattern_color="white")+
    geom_line() +
    scale_color_manual(values = cell.line.palette[cls])+
    annotate("line",x = log10(1e-6*drc.hd$Dose),y=drc.hd$fit_resp,color = cell.line.palette["HD.15S1"])+
    annotate("point",x=log10(1e-6*ec50.hd$ec50),y=0.5,color = cell.line.palette["HD.15S1"])+
    annotate(geom="errorbarh",x=log10(1e-6*ec50.hd$ec50),xmin=log10(1e-6*ec50.hd$ec50.upr),
             xmax=log10(1e-6*ec50.hd$ec50.lwr),y=0.5,height=0.1,color = cell.line.palette["HD.15S1"])+
    annotate("line",x= rep(log10(1e-6*ec50.hd$Cmax),2),y=c(0,1),linetype="dashed")+
    geom_errorbarh(data=ec50.mut,aes(x=log10(1e-6*ec50), xmin=log10(1e-6*ec50.upr),
                                     xmax=log10(1e-6*ec50.lwr), y = 0.5),height=0.1)+
    geom_point(data=ec50.mut,aes(x=log10(1e-6*ec50), y = 0.5))+
    facet_wrap(~cell.line,nrow=1,scales="free_y")+
    scale_y_continuous(limits = c(0,1),expand = c(0,0),
                       name="Probability\nArrhythmia")+
    coord_cartesian(xlim=c(min(log10(1e-6*drc.hd$Dose)),max(log10(1e-6*ec50.hd$maxDose))))+
    scale_x_continuous(expand = c(0,0), name = paste0("Log10[",icomp,"]"))+
    theme_classic(base_size = 9,
                  base_family = "sans")+
    theme(aspect.ratio = 1,
          strip.background=element_rect(color = NA,  fill=NA),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color="black"))
})
ml.mut <- marrangeGrob(plist.mut,ncol = 1,nrow = 1)
ggsave("Figures/all_drcs/all_MUTvsHD_prob_arrhythmic.pdf",ml.mut,width = 7,height = 4)

# Dose response curves and ec50s probability of asystole ---------------------------
aaa <- df %>% filter(BadBatch==0, !Compound %in% c("Control","GS-967","DMSO",NA)) %>%
  group_by(Compound) %>% group_split()
fff <- lapply(aaa, function(comp.df){
  comp.df <- comp.df %>% filter(Dose >0)
  ccc <-  comp.df %>% group_by(cell.line)%>% group_split()
  ddd <- lapply(ccc,fit.binomial,var.name="Prob.Cessation.AI")
  fitted.values <- lapply(ddd,'[[',"fit.fitted") %>% bind_rows()
  ec50s <- lapply(ddd,'[[',"fit.ec50") %>% bind_rows()
  l <- list(fitted.values,ec50s)
})
drcs.cessation.df <- lapply(fff, '[[',1) %>% bind_rows()
save(file="files/drcs.cessation.Rda",drcs.cessation.df)

maxdoses <- df %>% group_by(Compound) %>% select(Dose) %>% distinct() %>%   summarize(maxDose = max(Dose))
ec50s.cessation   <- lapply(fff, '[[',2) %>% bind_rows()
ec50s.cessation <- left_join(ec50s.cessation,df %>% select(Compound,Cmax,CiPA.Classification) %>% distinct())
ec50s.cessation <- left_join(ec50s.cessation,maxdoses)
ec50s.cessation <- ec50s.cessation %>% mutate(ec50.NA = ifelse(ec50>10*maxDose,NA,ec50),
                                              safety.margin =log10(ec50/Cmax))
save(file="files/ec50s.cessation.Rda",ec50s.cessation)



# plot all drcs for probability of asystole -------------------------------
ccc <- drcs.cessation.df %>% group_by(Compound) %>% group_split()
comp.df <- ccc[[22]]

plist2 <- lapply(ccc, function(comp.df){
  icomp <- comp.df$Compound[1]
  raw.data <- df %>% filter(Compound==icomp, BadBatch == 0)
  ec50.df <- ec50s.cessation %>% filter(Compound == icomp)
  cls <- unique(ec50.df$cell.line)
  ggplot(comp.df,aes(x=log10(1e-6*Dose),y=fit_resp,color = cell.line))+
    geom_line() +
    geom_point(data = raw.data,aes(x=log10(1e-6*Dose),y=Prob.Cessation.AI,
                                   color = cell.line),size=1,alpha=0.5)+
    scale_color_manual(values = cell.line.palette[cls])+
    geom_errorbarh(data=ec50.df,aes(x=log10(1e-6*ec50), xmin=log10(1e-6*ec50.upr),
                                     xmax=log10(1e-6*ec50.lwr), y = 0.5),height=0.1)+
    annotate("line",x= rep(log10(1e-6*ec50.df$Cmax[1]),2),y=c(0,1),linetype="dashed")+
    geom_point(data=ec50.df,aes(x=log10(1e-6*ec50), y = 0.5))+
    facet_wrap(~cell.line,nrow=2,scales="free")+
    scale_y_continuous(limits = c(0,1),expand = c(0,0),
                       name="Probability\nAsystole")+
    coord_cartesian(xlim=c(min(log10(1e-6*comp.df$Dose)),max(log10(1e-6*comp.df$Dose))))+
    scale_x_continuous(expand = c(0,0), name = paste0("Log10[",icomp,"]"))+
    theme_classic(base_size = 9,
                  base_family = "sans")+
    theme(aspect.ratio = 1,
          strip.background=element_rect(color = NA,  fill=NA),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color="black"))
})
ml2 <- marrangeGrob(plist2,ncol = 1,nrow = 1)
ggsave("Figures/all_drcs/all_prob_asystolic.pdf",ml2,width = 7,height = 4)


# plot all drcs for asystolic probability of Mutants vs Isogenic Control --------------------------
ccc <- drcs.cessation.df %>% filter(!cell.line %in% c("HD.113","HD.273")) %>% 
  group_by(Compound) %>% group_split()
comp.df <- ccc[[22]]

plist2.mut <- lapply(ccc, function(comp.df){
  icomp <- comp.df$Compound[1]
  drc.hd <- comp.df %>% filter(cell.line == "HD.15S1")
  drc.mut <- comp.df %>% filter(!grepl(pattern = "HD",x=cell.line))
  ec50.hd <- ec50s.cessation %>% filter(Compound == icomp, cell.line == "HD.15S1")
  ec50.mut <- ec50s.cessation %>% filter(Compound == icomp, !grepl(pattern = "HD",x=cell.line))
  cls <- unique(ec50.mut$cell.line)
  ggplot(drc.mut,aes(x=log10(1e-6*Dose),y=fit_resp,color = cell.line))+
    geom_line() +
    scale_color_manual(values = cell.line.palette[cls])+
    annotate("line",x = log10(1e-6*drc.hd$Dose),y=drc.hd$fit_resp,color = cell.line.palette["HD.15S1"])+
    annotate("point",x=log10(1e-6*ec50.hd$ec50),y=0.5,color = cell.line.palette["HD.15S1"])+
    annotate(geom="errorbarh",x=log10(1e-6*ec50.hd$ec50),xmin=log10(1e-6*ec50.hd$ec50.upr),
             xmax=log10(1e-6*ec50.hd$ec50.lwr),y=0.5,height=0.1,color = cell.line.palette["HD.15S1"])+
    geom_errorbarh(data=ec50.mut,aes(x=log10(1e-6*ec50), xmin=log10(1e-6*ec50.upr),
                                     xmax=log10(1e-6*ec50.lwr), y = 0.5),height=0.1)+
    annotate("line",x= rep(log10(1e-6*ec50.hd$Cmax),2),y=c(0,1),linetype="dashed")+
    geom_point(data=ec50.mut,aes(x=log10(1e-6*ec50), y = 0.5))+
    facet_wrap(~cell.line,nrow=1,scales="free_y")+
    scale_y_continuous(limits = c(0,1),expand = c(0,0),
                       name="Probability\nAsystole")+
    coord_cartesian(xlim=c(min(log10(1e-6*drc.hd$Dose)),max(log10(1e-6*drc.hd$Dose))))+
    scale_x_continuous(expand = c(0,0), name = paste0("Log10[",icomp,"]"))+
    theme_classic(base_size = 9,
                  base_family = "sans")+
    theme(aspect.ratio = 1,
          strip.background=element_rect(color = NA,  fill=NA),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color="black"))
})
ml2.mut <- marrangeGrob(plist2.mut,ncol = 1,nrow = 1)
ggsave("Figures/all_drcs/all_MUTvsHD_prob_asystolic.pdf",ml2.mut,width = 7,height = 4)


# Plot arrhythmia and cessation drcs in same page for a given drug --------
library(cowplot)
plist.grid <- lapply(c(1:length(plist)), function(iii){
 plot_grid(plist[[iii]],plist2[[iii]], labels = c('A', 'B'),ncol = 1)
  
})
m <- marrangeGrob(plist.grid,ncol = 1,nrow = 1,top=NULL)
ggsave("Figures/all_drcs/all_probs.pdf",m,width = 7,height = 6)

plist.mut.grid <- lapply(c(1:length(plist.mut)), function(iii){
  plot_grid(plist.mut[[iii]],plist2.mut[[iii]], labels = c('A', 'B'),ncol = 1)
  
})
m.mut <- marrangeGrob(plist.mut.grid,ncol = 1,nrow = 1,top=NULL)
ggsave("Figures/all_drcs/all_MUTvsHD_probs.pdf",m.mut,width = 7,height = 6)
