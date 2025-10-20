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


# type = 'intensity' or 'cumhaz'
PLOTLambda = function(hazard, selected_centres, type, add_hat=FALSE){
  if(type == 'intensity') ylab = TeX(paste('${\\lambda}_{ij}^*(t)$'))
  if(type == 'cumhaz'){
    ylab = TeX(paste('${\\Lambda}_{ij}(t)$'))
    if(add_hat == TRUE)
      ylab = TeX(paste('${\\hat{\\Lambda}}_{ij}(t)$'))
  } 
  
  hazard$color <- ifelse(!hazard$centre %in% selected_centers, "Others", as.factor(hazard$centre))
  hazard$linetype <- ifelse(!hazard$centre %in% selected_centers, "Others", as.factor(hazard$centre))
  hazard$size <- ifelse(!hazard$centre %in% selected_centers, 0.8, 1.2)
  hazard$centre <- factor(hazard$centre, levels = as.character(sort(unique(hazard$centre))))
  
  color_palette <- setNames(c(rgb(0.5,0.5,0.5, alpha = 0.12), colorRampPalette(brewer.pal(3,"Set1"))(length(selected_centers))),
                            c("Others", as.character(selected_centers)))
  line_palette <- setNames(c(1, 2:(length(selected_centers)+1)), 
                           c("Others", as.character(selected_centers)))
  
  PLOT = ggplot(hazard, aes(x= time, y=!!sym(type), group = factor(id), color=color, size=size, 
                     linetype = linetype )) +
    geom_line(size=0.9) +
    xlab('t') + 
    ylab(ylab) +
    labs(color='Cluster', linetype = "Cluster", shape = "Cluster", size = 'Cluster') +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text=element_text(size=rel(1)), 
          axis.title=element_text(size=rel(1.5)),  
          plot.title = element_text(face="bold", size=rel(1.5), hjust = 0.5),
          legend.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1.3)),
          legend.key.width= unit(1, 'cm'),
          legend.spacing.x = unit(0.5, 'cm')) + 
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = line_palette) +
    scale_size_identity(guide = 'none') + 
    theme(legend.key.width=unit(1.5,"cm"))
  
  if (type == 'cumhaz') {
    PLOT <- PLOT + ylim(0, 210)  
  }
  
  print(PLOT)
  
  if(add_hat == FALSE)
    ggsave(paste0("Plots/", type, ".pdf"), width = 8, height = 6, units = "in", device = "pdf")
  else
    ggsave(paste0("Plots/", type, "POST.pdf"), width = 8, height = 6, units = "in", device = "pdf")
}



PlotEigenFunMFPCA = function(plot_data, level, post=FALSE){
  if(level==1) pedix = 'k'
  if(level==2) pedix = 'l'
  
  if(post == FALSE) subscript = 'true'
  if(post == TRUE) subscript = 'fitted'
  
  PLOT = ggplot(plot_data, aes(x = x, y = y)) +
    geom_line(aes(y = !!sym(paste0("C1L", level)), linetype = "C1", color = "C1"), size = 1.4) +
    geom_line(aes(y = !!sym(paste0("C2L", level)), linetype = "C2", color = "C2"), size = 1.4) +
    ylim(-3,3) +
    xlab("t") +
    ylab(TeX(paste('${\\hat{\\phi}_{',pedix,',',subscript,'}^{(',level,')}(t)}$'))) +
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
                          labels = c(TeX(paste('${\\hat{\\phi}_{1,',subscript,'}^{(',level,')}(t)}$')), 
                                     TeX(paste('${\\hat{\\phi}_{2,',subscript,'}^{(',level,')}(t)}$')))) +
    scale_color_manual(name = "",
                       breaks = c("C1", "C2"),
                       values = c("C1" = "orange",
                                  "C2" = "purple"),
                       labels = c(TeX(paste('${\\hat{\\phi}_{1,',subscript,'}^{(',level,')}(t)}$')), 
                                  TeX(paste('${\\hat{\\phi}_{2,',subscript,'}^{(',level,')}(t)}$'))))
  print(PLOT)
  if(post == FALSE)
    ggsave(paste0("Plots/EigenFunL",level,".pdf"), width = 8, height = 4, units = "in", device = "pdf")
  if(post == TRUE)
    ggsave(paste0("Plots/EigenFunL",level,"POST.pdf"), width = 8, height = 4, units = "in", device = "pdf")
}




PlotCompMFPCA = function(plot_data, level = 1, component = 1, post = FALSE, coeff){
  PLOT = ggplot(plot_data, aes(x = x, y = y)) +
    geom_line(size = 1.6, color = "black") +
    geom_line(aes(y = !!sym(paste0("posC", component, "L", level)), linetype = "Positive", color = "Positive"), size = 1.2) +
    geom_line(aes(y = !!sym(paste0("negC", component, "L", level)), linetype = "Negative", color = "Negative"), size = 1.2) +
    ylim(-10,210) +
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
    ggsave(paste0("Plots/C",component,"L",level,".pdf"), width = 8, height = 4, units = "in", device = "pdf")
  if(post == TRUE)
    ggsave(paste0("Plots/C",component,"L",level,"POST.pdf"), width = 8, height = 4, units = "in", device = "pdf")
  
}



plot_Lambda0 <- function(list_L0){
  # t <- list_L0$time / 1000
  Lambda0 <- list_L0$Lambda0
  Lambda0S <- list_L0$Lambda0S
  
  data <- data.frame(x = Lambda0S$x[-1] / 1000,
                     fitted = Lambda0S$fitted[-1],
                     Lambda0 = Lambda0)
  
  PLOT = ggplot(data, aes(x = x)) +
    geom_line(aes(y = fitted, color = 'smoothed', linetype = 'smoothed'), size = 1.5) +
    geom_line(aes(y = Lambda0, color = 'nonsmoothed', linetype = 'nonsmoothed'), size = 1.5) +
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
    ylim(0, 160)  
  
  print(PLOT)
  ggsave("Plots/Baseline.pdf", width = 8, height = 6, units = "in", device = "pdf")
}

