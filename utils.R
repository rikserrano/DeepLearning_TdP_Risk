# Cell Line Palette
cell.line.palette <- c("#00A651","#EF4136","#0072E9","#0D9B97","#0C00FF","#53009B","#9325B7","#4851A3")
names(cell.line.palette) <- c("HD.113","HD.273","HD.15S1","DCM.TNNT2","DCM.RBM20","DCM.PLN","HCM.TPM1","HCM.MYBPC3")
# Helper function to reshape traces matrix for Keras
reshape_X_3d <- function(X) {
  dim(X) <- c(dim(X)[1], dim(X)[2], 1)
  X
}

# Helper function to encode labels
labels_to_one_hot <- function(labels){
  dmy <- dummyVars(" ~ .", data = data.frame(labels))
  df <- data.frame(predict(dmy, newdata = data.frame(labels)))
  as.matrix(df)
}

# Helper function to make contour density plots
library(ks)
get_level_lines <- function(d,lev,...){
  kd <- ks::kde(d, compute.cont=TRUE,...)
  contour_line <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                        z=estimate, levels=cont[paste0(lev,"%")])[[1]])
  contour_line <- data.frame(contour_line)
}

# Geometric sequence
seq.geometric <- function(from,to,by){
  s = from
  i = 1
  l = s
  while(l*by < to){
    s[i+1] = s[i]*by
    i = i+1
    l = s[i]
  }
  return(s)
}

# Plot trace with custom axis
plot.trace <- function(df,xcol,ycol){
  ggplot(df,aes(x={{xcol}},y={{ycol}}))+
    geom_line(size=0.5)+
    geom_segment(aes(x = -5, y = -0.05, xend = 33, yend = -0.05),size=0.75)+
    annotate("text", label = "1s", x = 0, y = -0.1, vjust=1,hjust = 0) +
    geom_segment(aes(x = -5, y = -0.05, xend = -5, yend = 0.3),size=0.75)+
    annotate("text", label = "0.3 \u0394 Fn", x = -8, y = -0,  hjust = 0,vjust=0,angle=90) +
    scale_x_continuous(expand = c(0.1,0))+
    scale_y_continuous(expand = c(0.1,0.13))
}  

# Theme for traces
theme_traces <- theme_void(base_family = "sans",
                           base_size = 11)+
  theme(strip.text = element_blank())

# Helper function to fit binomial curves for drcs
fit.binomial <- function(dum,var.name){
  fm <- as.formula(paste(var.name, "~log(Dose)"))
  obj <- glm(data = dum, formula = fm, 
             family="binomial")
  dum$fitted.values <- obj$fitted.values
  ilink <- family(obj)$linkinv
  newdata <- expand.grid(Dose=exp(seq(log(min(dum$Dose)), 
                                      log(max(dum$Dose)), length=100)))
  newdata <- bind_cols(newdata, setNames(as_tibble(predict(obj, newdata, se.fit = TRUE)[1:2]),
                                         c('fit_link','se_link')))
  newdata <- mutate(newdata,
                    fit_resp  = ilink(fit_link),
                    upr = ilink(fit_link + (2 * se_link)),
                    lwr = ilink(fit_link - (2 * se_link)),
                    cell.line = dum$cell.line[1],
                    Compound = dum$Compound[1])
  tryCatch(
    {
      newdata.ext <- expand.grid(Dose=exp(seq(log(min(dum$Dose)), 
                                              log(max(1e6*dum$Dose)), length=3000)))
      newdata.ext <- bind_cols(newdata.ext, setNames(as_tibble(predict(obj, newdata.ext, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))
      newdata.ext <- mutate(newdata.ext,
                            fit_resp  = ilink(fit_link),
                            upr = ilink(fit_link + (2 * se_link)),
                            lwr = ilink(fit_link - (2 * se_link)),
                            cell.line = dum$cell.line[1],
                            Compound = dum$Compound[1])
      ec50 <- exp(approx(newdata.ext$fit_resp,log(newdata.ext$Dose),0.5)$y)
      ec50.lwr <- exp(approx(newdata.ext$lwr,log(newdata.ext$Dose),0.5)$y)
      ec50.upr <- exp(approx(newdata.ext$upr,log(newdata.ext$Dose),0.5)$y)
      ec50.df <- data.frame(ec50,ec50.lwr,ec50.upr,
                            cell.line = dum$cell.line[1],
                            Compound = dum$Compound[1])
      l <- list(fit.data = dum, fit.obj = obj, fit.fitted=newdata, fit.ec50 = ec50.df)}, 
    error=function(e){
      ec50 <- NA
      ec50.lwr <- NA
      ec50.upr <- NA
      ec50.df <- data.frame(ec50 = NA,ec50.lwr = NA,ec50.upr = NA,
                            cell.line = dum$cell.line[1],
                            Compound = dum$Compound[1])
      l <- list(fit.data = dum, fit.obj = obj, fit.fitted=newdata, fit.ec50 = ec50.df)
    })
}

# Helper function to plot confusion matrix from caret object
plot.confusionMatrix <- function(cm,  nbins = 5){
  t<-cm$table %>% as_tibble()
  p <- prop.table(cm$table,margin = 2) %>% as_tibble() %>% rename(perc = n)
  df <- left_join(t,p)
  df$Reference <- factor(df$Reference,levels = c("Arrhythmic","NonArrhythmic","Cessation"))
  df$Prediction <- factor(df$Prediction,levels = c("Arrhythmic","NonArrhythmic","Cessation"))
  
  scale_fill_custom <- function (..., alpha = 1, begin = 0, end = 1, direction = 1, 
                                 option = "D", values = NULL, space = "Lab", na.value = "grey50", 
                                 guide = "colourbar", aesthetics = "fill") {
    continuous_scale(aesthetics, scale_name = "custom", 
                     palette = scales:::gradient_n_pal(c("white", "#b3cde0", "#6497b1", "#005b96", "#03396c"), 
                                                       values, space), na.value = na.value, 
                     guide = guide, ...)
  }
  binedges <- seq(0,1,length.out = nbins+1)
  pl <- ggplot(df,aes(x=Prediction,y=Reference,fill=perc))+
    geom_tile()+
    geom_text(aes(label= paste0(n,"\n(",percent(perc,accuracy = 0.1),")"),
                  color = perc > 0.5),size=2) +
    scale_color_manual(guide = "none", values = c("black", "white"))+
    scale_y_discrete(limits = rev(levels(df$Prediction)),expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_binned(name = "",breaks = binedges, 
                      labels = percent(binedges),
                      guide = guide_bins(show.limits = T,reverse = T,axis = F),
                      type = scale_fill_custom)
  pl
}
