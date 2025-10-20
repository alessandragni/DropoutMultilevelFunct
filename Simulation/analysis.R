# This script loads RDS simulation files, computes summary statistics for MISE values, and formats the results for LaTeX output.

# List all RDS simulation files in the working directory
files <- list.files(pattern = "^ManySim_.*\\.rds$")

# Initialize list to store summaries
all_summaries <- list()

for (file in files) {
  cat("Processing:", file, "\n")
  
  # Read file
  results <- readRDS(file)
  
  # Extract parameters from filename
  # Example: "ManySim_I20_J10_GRID1000_BALANCED.rds"
  parts <- strsplit(gsub("\\.rds$", "", file), "_")[[1]]
  
  # Extract I, J, GRID, and BALANCE info
  I_val <- as.numeric(sub("I", "", parts[2]))
  J_val <- as.numeric(sub("J", "", parts[3]))
  GRID_val <- as.numeric(sub("GRID", "", parts[4]))
  if (grepl("UNBALANCED", file)) {
    balanced_val <- FALSE
  } else if (grepl("BALANCED", file)) {
    balanced_val <- TRUE
  } else {
    balanced_val <- NA
  }
  
  # Extract MISE values from each result
  MISE1_vals <- sapply(results, function(x) x$MISE1)
  MISE2_vals <- sapply(results, function(x) x$MISE2)
  
  # Compute summary statistics
  summary_stats <- data.frame(
    I = I_val,
    J = J_val,
    grid = GRID_val,
    balanced = balanced_val,
    MISE1_mean = mean(MISE1_vals),
    MISE1_sd = sd(MISE1_vals),
    MISE1_median = median(MISE1_vals),
    MISE2_mean = mean(MISE2_vals),
    MISE2_sd = sd(MISE2_vals),
    MISE2_median = median(MISE2_vals)
  )
  
  # Add to list
  all_summaries[[length(all_summaries) + 1]] <- summary_stats
}

# Combine all summaries into a single data frame
all_summaries_df <- do.call(rbind, all_summaries)

library(dplyr)

# Round numeric columns
all_summaries_df <- all_summaries_df %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

# Ensure both TRUE and FALSE rows exist for each combination
all_summaries_df <- all_summaries_df %>%
  arrange(I, J, grid, balanced)

# Combine balanced and unbalanced side by side
combined_df <- all_summaries_df %>%
  group_by(I, J, grid) %>%
  summarise(
    MISE1_mean_bal = MISE1_mean[balanced == TRUE],
    MISE1_sd_bal = MISE1_sd[balanced == TRUE],
    MISE1_median_bal = MISE1_median[balanced == TRUE],
    MISE2_mean_bal = MISE2_mean[balanced == TRUE],
    MISE2_sd_bal = MISE2_sd[balanced == TRUE],
    MISE2_median_bal = MISE2_median[balanced == TRUE],
    
    MISE1_mean_unbal = MISE1_mean[balanced == FALSE],
    MISE1_sd_unbal = MISE1_sd[balanced == FALSE],
    MISE1_median_unbal = MISE1_median[balanced == FALSE],
    MISE2_mean_unbal = MISE2_mean[balanced == FALSE],
    MISE2_sd_unbal = MISE2_sd[balanced == FALSE],
    MISE2_median_unbal = MISE2_median[balanced == FALSE],
    .groups = "drop"
  )

# Convert to LaTeX-friendly output
latex_lines <- apply(combined_df, 1, function(row) {
  paste(row, collapse = " & ")
})

latex_table <- paste0(latex_lines, " \\\\", collapse = "\n")

cat("\n\n--- LATEX TABLE ---\n")
cat(latex_table)

