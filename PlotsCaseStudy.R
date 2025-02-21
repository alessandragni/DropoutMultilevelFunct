library(rlang)
library(gcmrec)   # Package for counting process model estimation
library(cobs)     # Package for constrained L1 B-splines
library(splines)
library(boot)
library(data.table)
library(survival)
library(latex2exp)
library(tidyverse)
library(tidyr)
library(scales)
library(refund)
library(ggsci)



plot_Lambda0 <- function(list_L0){
  Lambda0 <- list_L0$Lambda0
  Lambda0S <- list_L0$Lambda0S
  
  data <- data.frame(x = Lambda0S$x[-1],
                     fitted = Lambda0S$fitted[-1],
                     Lambda0 = Lambda0)
  
  PLOT = ggplot(data, aes(x = x)) +
    geom_line(aes(y = fitted, color = 'smoothed', linetype = 'smoothed'), size = 1.5, alpha = 0.9) +
    geom_line(aes(y = Lambda0, color = 'nonsmoothed', linetype = 'nonsmoothed'), size = 1.5, alpha = 0.9) +
    xlab("t") +
    ylab(TeX(paste('$\\hat{\\Lambda}_0(t)$'))) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1.5)),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(1.5)),
      legend.key.width = unit(2, 'cm'),
      legend.spacing.x = unit(1, 'cm')
    ) +
    scale_linetype_manual(
      name = "",
      breaks = c("smoothed", "nonsmoothed"),
      values = c('smoothed' = "solid",
                 'nonsmoothed' = "dashed"),
      labels = c(TeX(paste('$\\tilde{\\Lambda}_0$', ' (smoothed)')),
                 TeX(paste('$\\hat{\\Lambda}_0$')))) +
    scale_color_manual(
      name = "",
      breaks = c("smoothed", "nonsmoothed"),
      values = c('smoothed' = "darkgreen",
                 'nonsmoothed' = "grey"),
      labels = c(TeX(paste('$\\tilde{\\Lambda}_0$', ' (smoothed)')),
                 TeX(paste('$\\hat{\\Lambda}_0$')))
    ) +
    theme(legend.title = element_blank()) +
    ylim(0, 35)  
  
  print(PLOT)
  ggsave("Plots/CaseStudy/CS_Baseline.pdf", width = 8, height = 6, units = "in", device = "pdf")
}





# type = 'intensity' or 'cumhaz'
PLOTLambda = function(hazard, type, add_hat=FALSE){
  if(type == 'intensity') ylab = TeX(paste('${\\lambda}_{ij}(t)$'))
  if(type == 'cumhaz'){
    ylab = TeX(paste('${\\Lambda}_{ij}(t)$'))
    if(add_hat == TRUE)
      ylab = TeX(paste('${\\hat{\\Lambda}}_{ij}(t)$'))
  } 
  
  #color_palette <- setNames(c(rgb(0.5,0.5,0.5, alpha = 0.12), 
  #                            pal_locuszoom("default")(length(as.factor(hazard$centre)))))
  
  PLOT = ggplot(hazard, aes(x= time, y=!!sym(type), group = factor(id), 
                            color = as.factor(centre), 
                            size = 1.2, 
                            linetype = as.factor(centre))) +
    geom_line(size=1, alpha=0.8) +
    xlab('t') + 
    ylab(ylab) +
    labs(color='Cluster', linetype = "Cluster", shape = "Cluster", size = 'Cluster') +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text=element_text(size=rel(1)), 
          axis.title=element_text(size=rel(1.2)),  
          plot.title = element_text(face="bold", size=rel(1.2), hjust = 0.5),
          legend.title = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1)),
          legend.key.width= unit(0.8, 'cm'),
          legend.spacing.x = unit(0.3, 'cm')) + 
    scale_color_manual(values = pal_locuszoom("default")(length(as.factor(hazard$centre)))) +
    scale_size_identity(guide = 'none') + 
    theme(legend.key.width=unit(1.5,"cm"))
  
  if (type == 'cumhaz') {
    PLOT <- PLOT + ylim(0, 130)  
  }
  
  print(PLOT)
  
  if(add_hat == FALSE)
    ggsave(paste0("Plots/CaseStudy/CS_", type, ".pdf"), width = 8, height = 6, units = "in", device = "pdf")
  else
    ggsave(paste0("Plots/CaseStudy/CS_", type, "POST.pdf"), width = 8, height = 6, units = "in", device = "pdf")
}



# type = 'intensity' or 'cumhaz'
PLOTLambda2 = function(hazard, type, add_hat=FALSE){
  if(type == 'intensity') ylab = TeX(paste('${\\lambda}_{ij}(t)$'))
  if(type == 'cumhaz'){
    ylab = TeX(paste('${\\Lambda}_{ij}(t)$'))
    if(add_hat == TRUE)
      ylab = TeX(paste('${\\hat{\\Lambda}}_{ij}(t)$'))
  } 
  
  #color_palette <- setNames(c(rgb(0.5,0.5,0.5, alpha = 0.12), 
  #                            pal_locuszoom("default")(length(as.factor(hazard$centre)))))
  
  PLOT = ggplot(hazard, aes(x= time, y=!!sym(type), group = factor(id), 
                            color = as.factor(centre_noyear), 
                            size = 1.2, 
                            linetype = as.factor(year))) +
    geom_line(size=1, alpha=0.8) +
    xlab('t') + 
    ylab(ylab) +
    labs(color='School', linetype = "A.Year", shape = "School", size = 'School') +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text=element_text(size=rel(1)), 
          axis.title=element_text(size=rel(1.2)),  
          plot.title = element_text(face="bold", size=rel(1.2), hjust = 0.5),
          legend.title = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1)),
          legend.key.width= unit(0.8, 'cm'),
          legend.spacing.x = unit(0.2, 'cm')) + 
    #scale_color_viridis_d(option = "H") +
    scale_color_manual(values = c("#f781bf", "#ff7f00", "#8dd3c7", "#6a3d9a")) +
    #scale_color_manual(values = pal_locuszoom("default")(length(as.factor(hazard$centre_noyear)))) +
    scale_size_identity(guide = 'none') + 
    theme(legend.key.width=unit(1.6,"cm"))
  
  if (type == 'cumhaz') {
    PLOT <- PLOT + ylim(0, 130)  
  }
  
  print(PLOT)
  
  if(add_hat == FALSE)
    ggsave(paste0("Plots/CaseStudy/CS_", type, ".2pdf"), width = 10, height = 6, units = "in", device = "pdf")
  else
    ggsave(paste0("Plots/CaseStudy/CS_", type, "POST2.pdf"), width = 10, height = 6, units = "in", device = "pdf")
}



# type = 'intensity' or 'cumhaz'
PLOTLambda3 = function(hazard, type, add_hat=FALSE){
  if(type == 'intensity') ylab = TeX(paste('${\\lambda}_{ij}(t)$'))
  if(type == 'cumhaz'){
    ylab = TeX(paste('${\\Lambda}_{ij}(t)$'))
    if(add_hat == TRUE)
      ylab = TeX(paste('${\\hat{\\Lambda}}_{ij}(t)$'))
  } 
  
  #color_palette <- setNames(c(rgb(0.5,0.5,0.5, alpha = 0.12), 
  #                            pal_locuszoom("default")(length(as.factor(hazard$centre)))))
  
  PLOT = ggplot(hazard, aes(x= time, y=!!sym(type), group = factor(id), 
                            color = as.factor(id_noyear), 
                            size = 1.2, 
                            linetype = as.factor(year))) +
    geom_line(size=1, alpha=0.8) +
    xlab('t') + 
    ylab(ylab) +
    labs(color='Course', linetype = "A.Year", shape = "School", size = 'School') +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text=element_text(size=rel(1)), 
          axis.title=element_text(size=rel(1.2)),  
          plot.title = element_text(face="bold", size=rel(1.2), hjust = 0.5),
          legend.title = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1)),
          legend.key.width= unit(0.8, 'cm'),
          legend.spacing.x = unit(0.2, 'cm')) + 
    scale_color_viridis_d(option = "H") +
    #scale_color_manual(values = pal_locuszoom("default")(length(as.factor(hazard$centre_noyear)))) +
    scale_size_identity(guide = 'none') + 
    theme(legend.key.width=unit(1.6,"cm"))
  
  if (type == 'cumhaz') {
    PLOT <- PLOT + ylim(0, 130)  
  }
  
  print(PLOT)
  
  if(add_hat == FALSE)
    ggsave(paste0("Plots/CaseStudy/CS_", type, "3.pdf"), width = 10, height = 6, units = "in", device = "pdf")
  else
    ggsave(paste0("Plots/CaseStudy/CS_", type, "POST3.pdf"), width = 10, height = 6, units = "in", device = "pdf")
}







PlotEigenFunMFPCA = function(plot_data, level, post=FALSE){
  if(level==1) pedix = 'k'
  if(level==2) pedix = 'l'
  
  PLOT = ggplot(plot_data, aes(x = x, y = y)) +
    geom_line(aes(y = !!sym(paste0("C1L", level)), linetype = "C1", color = "C1"), size = 1.2) +
    geom_line(aes(y = !!sym(paste0("C2L", level)), linetype = "C2", color = "C2"), size = 1.2) +
    ylim(-3,3) +
    xlab("t") +
    ylab(TeX(paste('${\\hat{\\phi}_{',pedix,'}^{(',level,')}(t)}$'))) +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text=element_text(size=rel(1)),  
          axis.title=element_text(size=rel(1.5)),  
          legend.text = element_text(size = rel(1.5)),
          legend.key.width= unit(2, 'cm'),
          legend.spacing.x = unit(1, 'cm')) +  
    scale_linetype_manual(name = "",
                          breaks = c("C1", "C2"),
                          values = c("C1" = "dashed",
                                     "C2" = "dotdash"),
                          labels = c(TeX(paste('${\\hat{\\phi}_1^{(',level,')}(t)}$')), 
                                     TeX(paste('${\\hat{\\phi}_2^{(',level,')}(t)}$')))) +
    scale_color_manual(name = "",
                       breaks = c("C1", "C2"),
                       values = c("C1" = "orange",
                                  "C2" = "purple"),
                       labels = c(TeX(paste('${\\hat{\\phi}_1^{(',level,')}(t)}$')), 
                                  TeX(paste('${\\hat{\\phi}_2^{(',level,')}(t)}$'))))
  print(PLOT)
  if(post == FALSE)
    ggsave(paste0("Plots/CaseStudy/CS_",level,".pdf"), width = 8, height = 4, units = "in", device = "pdf")
  if(post == TRUE)
    ggsave(paste0("Plots/CaseStudy/CS_",level,"POST.pdf"), width = 8, height = 4, units = "in", device = "pdf")
}




PlotCompMFPCA = function(plot_data, level = 1, component = 1, post = FALSE, coeff){
  PLOT = ggplot(plot_data, aes(x = x, y = y)) +
    geom_line(size = 1.6, color = "black") +
    geom_line(aes(y = !!sym(paste0("posC", component, "L", level)), linetype = "Positive", color = "Positive"), size = 1.2) +
    geom_line(aes(y = !!sym(paste0("negC", component, "L", level)), linetype = "Negative", color = "Negative"), size = 1.2) +
    ylim(-10,70) +
    xlab("t") +
    ylab(TeX(paste('${\\hat{\\mu}(t)}$'))) +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text=element_text(size=rel(1)),  
          axis.title=element_text(size=rel(1.5)),  
          legend.text = element_text(size = rel(1.5)),
          legend.key.width= unit(2, 'cm'),
          legend.spacing.x = unit(1, 'cm')) +  
    scale_linetype_manual(name = "",
                          breaks = c("Positive", "Negative"),
                          values = c("Positive" = "dashed",
                                     "Negative" = "dotdash"),
                          labels = c(TeX(paste('${\\hat{\\mu}(t)} +', coeff, '{\\sqrt{{\\hat{\\lambda}_{',component,'}^{(',level,')}}}} {\\hat{\\phi}_{',component,'}^{(',level,')}(t)}$')), 
                                     TeX(paste('${\\hat{\\mu}(t)} -', coeff, '{\\sqrt{{\\hat{\\lambda}_{',component,'}^{(',level,')}}}} {\\hat{\\phi}_{',component,'}^{(',level,')}(t)}$')))) +
    scale_color_manual(name = "",
                       breaks = c("Positive", "Negative"),
                       values = c("Positive" = "red",
                                  "Negative" = "blue"),
                       labels = c(TeX(paste('${\\hat{\\mu}(t)} +', coeff, '{\\sqrt{{\\hat{\\lambda}_{',component,'}^{(',level,')}}}} {\\hat{\\phi}_{',component,'}^{(',level,')}(t)}$')), 
                                  TeX(paste('${\\hat{\\mu}(t)} -', coeff, '{\\sqrt{{\\hat{\\lambda}_{',component,'}^{(',level,')}}}} {\\hat{\\phi}_{',component,'}^{(',level,')}(t)}$'))))
  print(PLOT)
  if(post == FALSE)
    ggsave(paste0("Plots/CaseStudy/CS_",component,"L",level,".pdf"), width = 8, height = 4, units = "in", device = "pdf")
  if(post == TRUE)
    ggsave(paste0("Plots/CaseStudy/CS_",component,"L",level,"POST.pdf"), width = 8, height = 4, units = "in", device = "pdf")
  
}




