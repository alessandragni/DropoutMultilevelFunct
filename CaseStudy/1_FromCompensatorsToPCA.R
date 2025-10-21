#_______________________________________________________________________________
#### Import libraries ####

library(readxl)
library(dplyr)   
library(tidyverse)
library(refund)
library(survival)
library(survminer)
library(dplyr)
library(lubridate)
library(readxl)
library(stargazer)
library(openxlsx)
library(tidyverse)


#_______________________________________________________________________________
#### Import functions for the case study ####

source('PlotsCaseStudy.R')
source('utils.R')

AY = '2016'

#_______________________________________________________________________________
#### Import dataframe students_comp.Rda ####
# that contains the main career-related information after removing unwanted students, 
# cleaning degree types, and anonymizing sensitive info, as explained in Section 4.2.

df <- readRDS(file = "Data/students_comp.Rda")

#_______________________________________________________________________________
#### Formatting of the dataset ####

# Now we add N (students subscribed to the course, to normalize)
table(df[, 'stud_career_degree_name'])

count_degree_name <- df %>%
  group_by(stud_career_start_ay, stud_career_degree_name) %>%
  dplyr::summarize(N = n())

df <- df %>%
  left_join(count_degree_name, by = c("stud_career_start_ay", "stud_career_degree_name"))

# Filter students who dropped
df <- df[df$stud_career_status == 'D', ]

# See the years in the dataframe
table(df$stud_career_start_ay)

# Convert dates to Date objects
df$dropout_date <- as.Date(df$stud_career_end_date, format = "%Y-%m-%d") 

# Create the dataframe dropout_data
dropout_data <- df %>% 
  dplyr::select(stud_career_start_ay,
         stud_career_degree_area,
         stud_career_degree_code_CdS,
         stud_career_degree_code_cmp,
         N,
         dropout_date) %>%
  dplyr::filter(!is.na(dropout_date)) %>%
  dplyr::filter( 
           (dropout_date > as.Date(paste0(stud_career_start_ay, "-10-01"))) & 
           (dropout_date < as.Date(paste0(stud_career_start_ay+2, "-03-01"))) ) %>% 
  group_by(stud_career_degree_area, stud_career_start_ay, stud_career_degree_code_CdS, dropout_date) %>%
  arrange(dropout_date) %>%
  mutate(dropout_count = n()) %>%
  ungroup() %>%
  group_by(stud_career_degree_area, stud_career_start_ay, stud_career_degree_code_CdS) %>%
  mutate(start_date = as.Date(paste0(stud_career_start_ay, "-10-01")),
         Ntot = cumsum(dropout_count)) %>% 
  ungroup()
  
# Print the new dataset
print(dropout_data)

# Create start, stop and enum
dropout_data2 = dropout_data %>%
  arrange(stud_career_degree_area, stud_career_start_ay, stud_career_degree_code_CdS, dropout_date, dropout_count) %>%
  group_by(stud_career_degree_area, stud_career_start_ay, stud_career_degree_code_CdS) %>%
  mutate(
    start = lag(dropout_date, default = start_date[1]), # Lagged dropout date is the 'start' for the next row
    stop = dropout_date,
    enum = row_number() - 1  # Enumerate the dropout events starting from 0
  ) %>%
  ungroup() %>%
  distinct(stud_career_degree_area, stud_career_start_ay, stud_career_degree_code_CdS, start, stop, enum, .keep_all = TRUE)

dropout_data3 = dropout_data2 %>%
  group_by(stud_career_degree_area, stud_career_start_ay, stud_career_degree_code_CdS) %>%
  do(tibble::add_row(.,
                     stud_career_degree_area = .$stud_career_degree_area[nrow(.)], # Get the last row of the group
                     stud_career_start_ay = .$stud_career_start_ay[nrow(.)],
                     stud_career_degree_code_CdS = .$stud_career_degree_code_CdS[nrow(.)], 
                     stud_career_degree_code_cmp = .$stud_career_degree_code_cmp[nrow(.)],
                     start_date = .$start_date[nrow(.)],
                     N = .$N[nrow(.)],
                     Ntot = .$Ntot[nrow(.)],
                     dropout_count = .$dropout_count[nrow(.)],
                     start = max(.$stop),
                     stop = as.Date(paste0(.$stud_career_start_ay+2, "-03-01")),
                     enum = max(.$enum) + 1)) %>%
  distinct(stud_career_degree_area, stud_career_start_ay, stud_career_degree_code_CdS, start, stop, enum, .keep_all = TRUE) %>%
  mutate(
    diff_days = as.numeric(stop - start), 
    event = ifelse(row_number() == n(), 0, 1),  # If it's the last row, event is 0, otherwise 1
    dstart = cumsum(c(0, diff_days[-n()])),
    dstop = cumsum(diff_days)
  ) %>%
  ungroup()


# Print the new dataset
print(dropout_data3)


# Filter only enrolled in a.y. AY
dropout_AY <- dropout_data3 %>%
  filter(stud_career_start_ay == AY,
         start != stop) # %>%

  
df = dropout_AY %>% 
  dplyr::select(stud_career_degree_area, stud_career_degree_code_CdS, stud_career_degree_code_cmp, dstart, dstop, event, enum, N, Ntot, dropout_count) %>%
  rename(centre = stud_career_degree_area,
         id = stud_career_degree_code_CdS,
         start = dstart,
         stop = dstop) %>% 
  filter(start != stop)  

unique(df$centre)
length(unique(df$centre))
length(unique(df$id))
length(df$id)

# transform start and stop and event into
t_max <- max(df$stop)
times <- seq(0, t_max, by=1)
df$dropout_count_std = df$dropout_count/(df$N - df$Ntot)


#_______________________________________________________________________________
#### Andersen-Gill model ####

count_model = coxph(Surv(start,stop,event) ~ 1 + enum + dropout_count, 
                    id = id, 
                    cluster = id, 
                    data = df)


#_______________________________________________________________________________
#### For Appendix B (revision) ####

df[df$stop == 516, 'ratio'] = df[df$stop == 516, 'Ntot'] / df[df$stop == 516, 'N']
TEST = df[df$stop == 516, c('centre', 'id', 'ratio')]
colnames(TEST) <- c('School', 'Faculty', 'ratio')


#_______________________________________________________________________________
#### Compute Lambdas ####

bashaz <- fit_smooth_Lambda0(count_model)
plot_Lambda0(bashaz)


cumulative_hazard = compute_cumulative_hazard(model = count_model,
                                              sel_df = as.data.table(df),
                                              smoothed_baseline = bashaz$Lambda0S,
                                              times, verbose = TRUE)
cumulative_hazard = as.data.frame(cumulative_hazard)
cumulative_hazard$time = cumulative_hazard$time #/ 1256

PLOTLambda(cumulative_hazard, 'cumhaz', add_hat = TRUE)


df = reformat_cumhaz(cumulative_hazard)

#_______________________________________________________________________________
#### Perform MFPCA ####
res <- mfpca.face(Y = df$cumhaz, 
                  id = as.factor(df$centre),
                  pve = 0.99,
                  visit = as.factor(df$id), 
                  p=3,
                  twoway = FALSE)

res$evalues


#_______________________________________________________________________________
#### Compute the PERCENTAGE OF VARIABILITY EXPLAINED ####

rescomp <- mfpca.face(Y = df$cumhaz, 
                  id = as.factor(df$centre),
                  pve = 0.99999999,
                  visit = as.factor(df$id), 
                  p=3,
                  twoway = FALSE)
rescomp$evalues

# percentage of variability explained by the higher hierarchical level
sum(rescomp$evalues$level1) / (sum(rescomp$evalues$level1) + sum(rescomp$evalues$level2))

# percentage of variability explained by each PC
rescomp$evalues$level1[1] / sum(rescomp$evalues$level1)
rescomp$evalues$level1[2] / sum(rescomp$evalues$level1)

# percentage of variability explained by each PC
rescomp$evalues$level2[1] / sum(rescomp$evalues$level2)
rescomp$evalues$level2[2] / sum(rescomp$evalues$level2)

#_______________________________________________________________________________

# Assuming 'res' is my dataframe

# Create a data frame for ggplot
plot_data <- data.frame(x = seq_along(res$mu), 
                        y = res$mu, 
                        C1L1 = res$efunctions$level1[,1],
                        C2L1 = res$efunctions$level1[,2],
                        C1L2 = res$efunctions$level2[,1],
                        C2L2 = res$efunctions$level2[,2],
                        posC1L1 = res$mu + res$efunctions$level1[, 1] * 1 * sqrt(res$evalues$level1[1]), 
                        negC1L1 = res$mu - res$efunctions$level1[, 1] * 1 * sqrt(res$evalues$level1[1]),
                        posC2L1 = res$mu + res$efunctions$level1[, 2] * 3 * sqrt(res$evalues$level1[2]), 
                        negC2L1 = res$mu - res$efunctions$level1[, 2] * 3 * sqrt(res$evalues$level1[2]),
                        posC1L2 = res$mu + res$efunctions$level2[, 1] * 1 * sqrt(res$evalues$level2[1]),
                        negC1L2 = res$mu - res$efunctions$level2[, 1] * 1 * sqrt(res$evalues$level2[1]),
                        posC2L2 = res$mu + res$efunctions$level2[, 2] * 3 * sqrt(res$evalues$level2[2]),
                        negC2L2 = res$mu - res$efunctions$level2[, 2] * 3 * sqrt(res$evalues$level2[2])
)

plot_data <- data.frame(x = seq_along(res$mu), 
                        y = res$mu, 
                        C1L1 = res$efunctions$level1[,1],
                        C2L1 = res$efunctions$level1[,2],
                        C1L2 = res$efunctions$level2[,1],
                        C2L2 = res$efunctions$level2[,2],
                        posC1L1 = res$mu + res$efunctions$level1[, 1] * 1 * sqrt(res$evalues$level1[1]), 
                        negC1L1 = res$mu - res$efunctions$level1[, 1] * 1 * sqrt(res$evalues$level1[1]),
                        posC2L1 = res$mu + res$efunctions$level1[, 2] * 3 * sqrt(res$evalues$level1[2]), 
                        negC2L1 = res$mu - res$efunctions$level1[, 2] * 3 * sqrt(res$evalues$level1[2]),
                        posC1L2 = res$mu + res$efunctions$level2[, 1] * 1 * sqrt(res$evalues$level2[1]),
                        negC1L2 = res$mu - res$efunctions$level2[, 1] * 1 * sqrt(res$evalues$level2[1]),
                        posC2L2 = res$mu + res$efunctions$level2[, 2] * 3 * sqrt(res$evalues$level2[2]),
                        negC2L2 = res$mu - res$efunctions$level2[, 2] * 3 * sqrt(res$evalues$level2[2])
)

PlotEigenFunMFPCA(plot_data, level = 1, post = TRUE)
PlotEigenFunMFPCA(plot_data, level = 2, post = TRUE)

PlotCompMFPCA(plot_data, level = 1, component = 1, post = TRUE, coeff = '')
PlotCompMFPCA(plot_data, level = 1, component = 2, post = TRUE, coeff = '3')
PlotCompMFPCA(plot_data, level = 2, component = 1, post = TRUE, coeff = '')
PlotCompMFPCA(plot_data, level = 2, component = 2, post = TRUE, coeff = '3')



# res$efunctions$level1[, 1]
# res$efunctions$level2[, 1]

res$scores

#_______________________________________________________________________________
#### Save for the next script ####

SCHOOLS = data.frame('School' = unique(factor(df$centre)),
                     'ScoresSchool1' = res$scores$level1[,1], 
                     'ScoresSchool2' = res$scores$level1[,2]
                     )

FACULTIES = data.frame('School' = df$centre,
                       'Faculty' = df$id, 
                       'ScoresFaculty1' = res$scores$level2[,1], 
                       'ScoresFaculty2' = res$scores$level2[,2] )

DF = merge(SCHOOLS, FACULTIES, by = 'School')
DF$career_start_ay = rep(as.double(AY), nrow(DF))


#_______________________________________________________________________________
#### Save dataframes for next script ####

write.csv(DF, paste0('Data/Scores', AY, '.csv', sep=''), row.names = FALSE)

# for test with logistic in Appendix B
write.csv(TEST, paste0('Data/TEST', AY, '.csv', sep=''), row.names = FALSE)
