# This script imports necessary libraries, runs simulations, and processes results to evaluate and visualize eigenfunctions. 
# It also saves summary statistics and simulation results to files.

# import libraries
library(gcmrec)
library(cobs)   
library(splines)
library(boot)
library(data.table)
library(survival)
library(latex2exp)
library(tidyverse)
library(tidyr)
library(scales)
library(refund)
library(ggplot2)
library(RColorBrewer)


source('simulate_once.R')

#_______________________________________________________________________________
#### To simulate just 1 run and get Figures 1, 2 and 3 ####
results <- simulate_once(seed = 1, I = 20, J = 4, grid = 1000, balanced = TRUE, verbose = TRUE)

#_______________________________________________________________________________
#### To simulate just 1 run and get Figure 4a and 4b ####
grid = 1000
results <- lapply(c(1:8, 10:33, 35:38, 40:44, 46:50, 52:67, 69:71, 73:107), 
                  function(s) simulate_once(seed = s, grid = grid))


df_eval <- do.call(rbind, lapply(seq_along(results), function(i) {
  get_eval <- function(res, level, k) {
    vals <- res[[level]]
    return(vals[k])
  }
  
  data.frame(
    replication = i,
    level = rep(1:2, each = 2),
    component = rep(1:2, times = 2),
    sim0 = c(
      get_eval(results[[i]]$evalues_0, "level1", 1),
      get_eval(results[[i]]$evalues_0, "level1", 2),
      get_eval(results[[i]]$evalues_0, "level2", 1),
      get_eval(results[[i]]$evalues_0, "level2", 2)
    ),
    simest = c(
      get_eval(results[[i]]$evalues_est, "level1", 1),
      get_eval(results[[i]]$evalues_est, "level1", 2),
      get_eval(results[[i]]$evalues_est, "level2", 1),
      get_eval(results[[i]]$evalues_est, "level2", 2)
    )
  )
}))


df_eval_long <- df_eval %>%
  pivot_longer(cols = c(sim0, simest), names_to = "method", values_to = "value")


# Align signs before averaging
align_sign <- function(ref, mat) {
  if (sum(ref * mat) < 0) mat <- -mat
  return(mat)
}


efun_all <- function(level, k, grid) {
  if(level == 'level1') ref <- results[[5]]$efuncs_est[[level]][, k]
  if(level == 'level2') ref <- results[[3]]$efuncs_est[[level]][, k]
  all_est <- sapply(results, \(r) align_sign(ref, r$efuncs_est[[level]][, k]))
  
  mean_est <- apply(all_est, 1, mean)
  
  list(
    all = data.frame(
      grid = rep(grid, times = ncol(all_est)),
      value = as.vector(all_est),
      replication = rep(seq_len(ncol(all_est)), each = length(grid))
    ),
    summary = data.frame(
      grid = grid,
      mean_est = mean_est
    )
  )
}



grid_len <- seq(0, 1, length.out =  (grid+1))

# --- Example: Level 1 --- change here to get level2 plot
level <- "level1"
efun1 <- efun_all(level, 1, grid_len)
efun2 <- efun_all(level, 2, grid_len)

df_all <- rbind(
  cbind(efun1$all, component = "C1"),
  cbind(efun2$all, component = "C2")
)

df_summary <- rbind(
  cbind(efun1$summary, component = "C1"),
  cbind(efun2$summary, component = "C2")
)

if (level == "level1") pedix <- 'k'
if (level == "level2") pedix <- 'l'
subscript <- 'fitted'

# --- Plot (Median only) ---
PLOT <- ggplot() +
  # All individual curves (faint orange/purple)
  geom_line(
    data = df_all,
    aes(x = grid, y = value, group = interaction(replication, component), color = component),
    alpha = 0.08, linewidth = 0.5
  ) +
  # Median curves (darker orange/purple)
  geom_line(
    data = subset(df_summary, component == "C1"),
    aes(x = grid, y = mean_est, linetype = "C1"),
    color = "orange", linewidth = 1.8
  ) +
  geom_line(
    data = subset(df_summary, component == "C2"),
    aes(x = grid, y = mean_est, linetype = "C2"),
    color = "purple", linewidth = 1.8
  ) +
  ylim(-3, 3) +
  xlab("t") +
  ylab(TeX(paste0('${\\hat{\\phi}_{', pedix, ',', subscript, '}^{(', gsub("level", "", level), ')}(t)}$'))) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = rel(1)),
    axis.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.3)),
    legend.key.width = unit(2, 'cm'),
    legend.spacing.x = unit(1, 'cm')
  ) +
  scale_linetype_manual(
    name = "",
    values = c("C1" = "dashed", "C2" = "dotdash"),
    labels = c(
      TeX(paste0('${\\hat{\\phi}_{1,', subscript, '}^{(', gsub("level", "", level), ')}(t)}$')),
      TeX(paste0('${\\hat{\\phi}_{2,', subscript, '}^{(', gsub("level", "", level), ')}(t)}$'))
    )
  ) +
  scale_color_manual(
    name = "",
    values = c("C1" = "#FDBB84", "C2" = "#C2A5CF"),
    labels = c(
      TeX(paste0('${\\hat{\\phi}_{1,', subscript, '}^{(', gsub("level", "", level), ')}(t)}$')),
      TeX(paste0('${\\hat{\\phi}_{2,', subscript, '}^{(', gsub("level", "", level), ')}(t)}$'))
    )
  )

print(PLOT)


ggsave(paste0("Plots/EigenFunL",level,"POST_MANY.pdf"), width = 8, height = 4, units = "in", device = "pdf")


#_______________________________________________________________________________
#### To simulate over many cases and 100 runs each and get MISE (for Table 1) ####

# example
# I = 10 
# J = 4 
# balanced = FALSE 
# grid = 200
# results <- lapply(
#   c(1),
#   function(s) simulate_once(seed = s, I = I, J = J, grid = grid, balanced = balanced)
# )

# COMPUTATION OF THE MISE IN TABLE 
# Initialize a list to store summaries
all_summaries <- list()

# Loop over all combinations
for (I in c(10, 20, 40)) {
  for (J in c(4, 10, 50)) {
    for (balanced in c(TRUE, FALSE)) {
      for (grid in c(200, 1000, 5000)){

      
       cat("Running simulations for I =", I, ", J =", J, ", balanced =", 
           balanced, ", grid =", grid, "\n")
       
       # Run all simulations
       results <- lapply(
         c(1:8, 10:33, 35:38, 40:44, 46:50, 52:67, 69:71, 73:107),
         function(s) simulate_once(seed = s, I = I, J = J, grid = grid, balanced = balanced)
       )
       
       # Extract MISE values
       MISE1_vals <- sapply(results, function(x) x$MISE1)
       MISE2_vals <- sapply(results, function(x) x$MISE2)
       
       # Compute summary statistics
       summary_stats <- data.frame(
         I = I,
         J = J,
         grid = grid,
         balanced = balanced,
         MISE1_mean = mean(MISE1_vals),
         MISE1_sd = sd(MISE1_vals),
         MISE1_median = median(MISE1_vals),
         MISE2_mean = mean(MISE2_vals),
         MISE2_sd = sd(MISE2_vals),
         MISE2_median = median(MISE2_vals)
       )
       
       print(summary_stats)
       
       # Append to list
       all_summaries[[length(all_summaries) + 1]] <- summary_stats
       
       # Save full simulation results
       filename <- paste0(
         "ManySim_I", I, "_J", J, "_GRID", grid,
         ifelse(balanced, "_BALANCED.rds", "_UNBALANCED.rds")
       )
       saveRDS(results, filename)
       
      }
    }
  }
}

# Combine all summaries into a single data frame
all_summaries_df <- do.call(rbind, all_summaries)

# Save summaries
saveRDS(all_summaries_df, "Summary_MISE_All.rds")

# Print summary table
print(all_summaries_df)
