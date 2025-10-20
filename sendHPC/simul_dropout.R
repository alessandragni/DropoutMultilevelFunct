# This script imports necessary libraries, runs simulations, and processes results to evaluate and visualize eigenfunctions. 
# It also saves summary statistics and simulation results to files.

# import libraries
# library(gcmrec, lib.loc="/u/ragni/Rlibs")
library(cobs, lib.loc="/u/ragni/Rlibs") 
library(splines, lib.loc="/u/ragni/Rlibs")
library(boot, lib.loc="/u/ragni/Rlibs")
library(data.table, lib.loc="/u/ragni/Rlibs")
library(survival, lib.loc="/u/ragni/Rlibs")
library(withr, lib.loc="/u/ragni/Rlibs")
# library(latex2exp, lib.loc="/u/ragni/Rlibs")
library(tidyr, lib.loc="/u/ragni/Rlibs")
library(dplyr, lib.loc="/u/ragni/Rlibs")
library(tibble, lib.loc="/u/ragni/Rlibs")
library(farver, lib.loc="/u/ragni/Rlibs")
library(RColorBrewer, lib.loc="/u/ragni/Rlibs")
library(scales, lib.loc="/u/ragni/Rlibs")
library(refund, lib.loc="/u/ragni/Rlibs")

setwd("/u/ragni/dropout")
source("/u/ragni/dropout/simulate_once.R")

args <- commandArgs(trailingOnly = TRUE)
# Parse input parameters
I <- as.numeric(args[1])
J <- as.numeric(args[2])
balanced <- as.logical(args[3])
grid <- as.numeric(args[4])

#_______________________________________________________________________________
#### To simulate over many cases and 100 runs each and get MISE (for Table 1) ####

# COMPUTATION OF THE MISE IN TABLE 
        
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
 
   
   # Save full simulation results
   filename <- paste0(
     "OUTPUTS/ManySim_I", I, "_J", J, "_GRID", grid,
     ifelse(balanced, "_BALANCED.rds", "_UNBALANCED.rds")
   )
   saveRDS(results, filename)
        
