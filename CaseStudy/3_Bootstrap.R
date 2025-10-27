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
# Initialize lists to store results across bootstrap iterations
perc_var_exp = list()
perc_var_exp_level1_1 = list()
perc_var_exp_level1_2 = list()
perc_var_exp_level2_1 = list()
perc_var_exp_level2_2 = list()
coefficients = list()
concordance = list()

#_______________________________________________________________________________

#### Bootstrap procedure starts ####

set.seed(123)

for(b in 1:1000){
  
  print(paste0("Bootstrap iteration: ", b))
  
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
  #### Compute Lambdas ####
  
  bashaz <- fit_smooth_Lambda0(count_model)
  
  
  cumulative_hazard = compute_cumulative_hazard(model = count_model,
                                                sel_df = as.data.table(df),
                                                smoothed_baseline = bashaz$Lambda0S,
                                                times, verbose = TRUE)
  cumulative_hazard = as.data.frame(cumulative_hazard)
  cumulative_hazard$time = cumulative_hazard$time #/ 1256
  
  
  df = reformat_cumhaz(cumulative_hazard)
  
  #_______________________________________________________________________________
  #### Perform MFPCA ####
  res <- mfpca.face(Y = df$cumhaz, 
                    id = as.factor(df$centre),
                    pve = 0.99,
                    visit = as.factor(df$id), 
                    p=3,
                    twoway = FALSE)
  
  #_______________________________________________________________________________
  #### Compute the PERCENTAGE OF VARIABILITY EXPLAINED ####
  
  rescomp <- mfpca.face(Y = df$cumhaz, 
                        id = as.factor(df$centre),
                        pve = 0.99999999,
                        visit = as.factor(df$id), 
                        p=3,
                        twoway = FALSE)
  
  # percentage of variability explained by the higher hierarchical level
  perc_var_exp = c(perc_var_exp, sum(rescomp$evalues$level1) / (sum(rescomp$evalues$level1) + sum(rescomp$evalues$level2)))
  
  # percentage of variability explained by each PC
  perc_var_exp_level1_1 = c(perc_var_exp_level1_1, rescomp$evalues$level1[1] / sum(rescomp$evalues$level1))
  perc_var_exp_level1_2 = c(perc_var_exp_level1_2, rescomp$evalues$level1[2] / sum(rescomp$evalues$level1))
  
  # percentage of variability explained by each PC
  perc_var_exp_level2_1 = c(perc_var_exp_level2_1, rescomp$evalues$level2[1] / sum(rescomp$evalues$level2))
  perc_var_exp_level2_2 = c(perc_var_exp_level2_1, rescomp$evalues$level2[2] / sum(rescomp$evalues$level2))
  
  
  
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
  #### next script ####
  
  #_______________________________________________________________________________
  #### Import dataframe students_COX_timdep.Rda ####
  # that is designed for time-dependent Cox models and includes exam information at the end of 
  # first semester with student career info. It is a survival-analysis-ready dataset.
  
  students_COX_timdep <- readRDS(file = "Data/students_COX_timdep.Rda")
  students_COX_timdep = students_COX_timdep[students_COX_timdep$career_start_ay %in% c(2017), ] 
  
  # Generated in script 1_FromCompensatorsToPCA.R
  Scores2016 <- DF
  
  SCORES = rbind(Scores2016)
  students_COX_timdep = left_join(students_COX_timdep, SCORES, by = join_by('School', 'Faculty'))
  
  #_______________________________________________________________________________
  #### Some further preprocessing ####
  
  # Set categorical as factors
  students_COX_timdep$Gender <- as.factor(students_COX_timdep$Gender)
  students_COX_timdep$Age19 <- as.factor(students_COX_timdep$Age19)
  students_COX_timdep$Income <- as.factor(students_COX_timdep$Income)
  students_COX_timdep$Origins <- as.factor(students_COX_timdep$Origins)
  students_COX_timdep$HighschoolType <- as.factor(students_COX_timdep$HighschoolType)
  
  # Chose reference level
  students_COX_timdep$Gender = relevel(students_COX_timdep$Gender, ref = "Male")
  students_COX_timdep$Income = relevel(students_COX_timdep$Income, ref = "Medium")
  students_COX_timdep$Origins = relevel(students_COX_timdep$Origins, ref = "Onsite")
  students_COX_timdep$HighschoolType = relevel(students_COX_timdep$HighschoolType, ref = "Scientific")
  
  
  #_____________________
  
  # doing stepAIC
  
  # this is the output (benchmark model in Tables 4 and 5)
  td.Cox.drop2<- coxph(Surv(time_sem, dropout) 
                       ~ ScoresSchool1
                       + ScoresFaculty1
                       + ScoresSchool2
                       + Origins 
                       + HighschoolType
                       + Income 
                       + HighschoolGrade 
                       + ECTS1sem,
                       data = students_COX_timdep)
  concordance = c(concordance, td.Cox.drop2$concordance['concordance'])
  coefficients = c(coefficients, td.Cox.drop2$coefficients)

}

save.image(file='bootstrapresult.RData')

#### PLOTS ####


# --- Var explained ---

perc_var_exp <- unlist(perc_var_exp)
perc_var_exp_level1_1 <- unlist(perc_var_exp_level1_1)
perc_var_exp_level1_2 <- unlist(perc_var_exp_level1_2)
perc_var_exp_level2_1 <- unlist(perc_var_exp_level2_1)
perc_var_exp_level2_2 <- unlist(perc_var_exp_level2_2)

df = data.frame( 'Variable'=c(rep('perc_var_exp', length(perc_var_exp)),
                              rep('perc_var_exp_level1_1', length(perc_var_exp_level1_1)),
                              rep('perc_var_exp_level1_2', length(perc_var_exp_level1_2)),
                              rep('perc_var_exp_level2_1', length(perc_var_exp_level2_1)),
                              rep('perc_var_exp_level2_2', length(perc_var_exp_level2_2))),
                 'Value'= c(perc_var_exp, perc_var_exp_level1_1, perc_var_exp_level1_2, perc_var_exp_level2_1, perc_var_exp_level2_2),
                 'technicalname'=c(
                   rep('Lev~1', length(perc_var_exp)),
                   rep('Lev~1~comp~1', length(perc_var_exp_level1_1)),
                   rep('Lev~1~comp~2', length(perc_var_exp_level1_2)),
                   rep('Lev~2~comp~1', length(perc_var_exp_level2_1)),
                   rep('Lev~2~comp~2', length(perc_var_exp_level2_2))),
                 'yintercept'=c(rep(0.5839, length(perc_var_exp_level1_1)),
                                rep(0.9899, length(perc_var_exp_level1_1)),
                                rep(0.01, length(perc_var_exp_level1_2)),
                                rep(0.9593, length(perc_var_exp_level2_1)),
                                rep(0.03, length(perc_var_exp_level2_2)))
)


# Plot with label_parsed so ggplot interprets LaTeX syntax
ggplot(df, aes(x = NULL, y = Value, fill = Variable)) +
  geom_boxplot(fill = "Corn Silk", color = "black", width = 0.5) +
  geom_hline(aes(yintercept = yintercept),
             color = "gray", linetype = "dashed", size = 1) +
  facet_wrap(~ technicalname, scales = "free_y", labeller = label_parsed, nrow=1) +
  labs(title = "Proportion of Explained Variance",
       x = NULL, y = "Values") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

ggsave("Plots/VarExplained.pdf", 
       width = 11, height = 4, units = "in", device = "pdf")

# --- Concordance and coefficients ---

scores <- coefficients[grep("Scores(School1|Faculty1|School2)", names(coefficients))]
df <- data.frame(Variable = names(scores), Value = unlist(scores))
df$Variable <- gsub(".*Scores", "Scores", df$Variable)


hlines_scores <- data.frame(
  Variable = c("ScoresSchool1", "ScoresFaculty1","ScoresSchool2"),
  yintercept = c(-0.020, 0.006, -0.885),  # example custom intercepts
  technicalname = c(
    "hat(alpha)[1]^{(1)}(t)",
    "hat(alpha)[1]^{(2)}(t)",
    "hat(alpha)[2]^{(1)}(t)"
  ),
  stringsAsFactors = FALSE
)
# Merge data
df_scores <- left_join(df, hlines_scores, by = "Variable")



concordance_vec <- unlist(concordance)
df = data.frame( 'Variable'=c(rep('C-index', length(concordance_vec))),
                 'Value'= concordance_vec,
                 'technicalname'=c(rep('C-index', length(concordance_vec))),
                 'yintercept'=c(rep(0.847, length(concordance_vec)))
)

df_scores = rbind(df_scores, df)

# Plot with label_parsed so ggplot interprets LaTeX syntax
ggplot(df_scores, aes(x = NULL, y = Value, fill = Variable)) +
  geom_boxplot(fill = "Light Cyan", color = "black", width = 0.5) +
  geom_hline(aes(yintercept = yintercept),
             color = "gray", linetype = "dashed", size = 1) +
  facet_wrap(~ technicalname, scales = "free_y", labeller = label_parsed, nrow=1) +
  labs(title = "C-index and Model Coefficients",
       x = NULL, y = "Values") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )


ggsave("Plots/Cindexandcoeff.pdf", 
       width = 11, height = 4, units = "in", device = "pdf")




