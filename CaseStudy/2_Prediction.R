#_______________________________________________________________________________
#### Import libraries ####

library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)
library(stargazer)
library(viridis)
library(cutpointr)
library(PRROC)
library(splitTools)
library(caret)
library(ISLR)
library(pROC) 
library(lme4)
library(glmnet)
library(MASS)

#_______________________________________________________________________________
#### Import dataframe students_COX_timdep.Rda ####
# that is designed for time-dependent Cox models and includes exam information at the end of 
# first semester with student career info. It is a survival-analysis-ready dataset.

students_COX_timdep <- readRDS(file = "Data/students_COX_timdep.Rda")
students_COX_timdep = students_COX_timdep[students_COX_timdep$career_start_ay %in% c(2017), ] 

# Generated in script 1_FromCompensatorsToPCA.R
Scores2016 <- read.csv(file = "Data/Scores2016.csv")
TEST <- read.csv(file = "Data/TEST2016.csv")

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


#_______________________________________________________________________________
#### Descriptive statistics ####

numericals = c('AdmissionScore', 'HighschoolGrade', 'ECTS1sem')
categoricals = c('Origins', 'Gender', 'HighschoolType', 'Income', 'Age19')

round(colMeans(students_COX_timdep[students_COX_timdep$dropout=='0', numericals]),2)
round(apply(students_COX_timdep[students_COX_timdep$dropout=='0', numericals], 2, sd),2)
summary(students_COX_timdep[students_COX_timdep$dropout=='0', categoricals]) 
dim(students_COX_timdep[students_COX_timdep$dropout=='0', ])[1]

round(colMeans(students_COX_timdep[students_COX_timdep$dropout=='1', numericals]),2)
round(apply(students_COX_timdep[students_COX_timdep$dropout=='1', numericals], 2, sd),2)
summary(students_COX_timdep[students_COX_timdep$dropout=='1', categoricals])
dim(students_COX_timdep[students_COX_timdep$dropout=='1', ])[1]

length(students_COX_timdep$StudentID)
length(unique(students_COX_timdep$StudentID))

length(unique(students_COX_timdep$StudentID))



#_______________________________________________________________________________
#### Models for Tables 4 and 5 ####

# Define survival outcome
y <- Surv(students_COX_timdep$time_sem, students_COX_timdep$dropout)   # time-to-event & censoring

# Prepare predictor matrix (excluding survival outcome)
X <- model.matrix(~ ScoresSchool1 
                  + ScoresFaculty1
                  + ScoresSchool2 
                  + ScoresFaculty2
                  + Origins 
                  + Gender 
                  + HighschoolType
                  + Income 
                  + Age19
                  + HighschoolGrade 
                  + AdmissionScore 
                  + ECTS1sem, data = students_COX_timdep)[, -1]  # Remove intercept

# Fit LASSO Cox model
fit <- glmnet(X, y, family = "cox", alpha = 1)  # alpha=1 for pure LASSO

# Cross-validation to find optimal lambda
cv.fit <- cv.glmnet(X, y, family = "cox", alpha = 1)

# Optimal lambda
best_lambda <- cv.fit$lambda.min

# Coefficients of the selected model
coef(fit, s = best_lambda)

# Plot cross-validation results
plot(cv.fit)


#_____________________

# doing stepAIC
td.Cox.drop<- coxph(Surv(time_sem, dropout) 
                    ~ ScoresSchool1 
                    + ScoresFaculty1
                    + ScoresSchool2 
                    + ScoresFaculty2
                    + Origins 
                    + Gender 
                    + HighschoolType
                    + Income 
                    + Age19
                    + HighschoolGrade 
                    + AdmissionScore 
                    + ECTS1sem,
                    data = students_COX_timdep)

stepAIC(td.Cox.drop, direction='both')
summary(td.Cox.drop)
c2<-concordance(td.Cox.drop)$concordance
c2


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
summary(td.Cox.drop2)
logLik(td.Cox.drop2)
BIC(td.Cox.drop2)
AIC(td.Cox.drop2)
c22<-concordance(td.Cox.drop2)$concordance
c22

anova(td.Cox.drop, td.Cox.drop2)

#_____________________

# model with no PC to see their individual contributions (a of Table 5)
td.Cox.drop5<- coxph(Surv(time_sem, dropout) 
                     ~ Origins 
                     + HighschoolType
                     + Income 
                     + HighschoolGrade 
                     + ECTS1sem,
                     data = students_COX_timdep)
summary(td.Cox.drop5)
logLik(td.Cox.drop5)
BIC(td.Cox.drop5)
AIC(td.Cox.drop5)
c22<-concordance(td.Cox.drop5)$concordance
c22

anova(td.Cox.drop5, td.Cox.drop2)

#_____________________

# model with only school-level PC to see their individual contributions (b of Table 5)
td.Cox.drop3<- coxph(Surv(time_sem, dropout) 
                     ~ ScoresSchool1
                     + ScoresSchool2
                     + Origins
                     + HighschoolType 
                     + Income 
                     + HighschoolGrade 
                     + ECTS1sem,
                     data = students_COX_timdep)
summary(td.Cox.drop3)
logLik(td.Cox.drop3)
BIC(td.Cox.drop3)
AIC(td.Cox.drop3)
c22<-concordance(td.Cox.drop3)$concordance
c22

anova(td.Cox.drop3, td.Cox.drop2)

#_____________________

# model with only course-level PC to see their individual contributions (c of Table 5)

td.Cox.drop4<- coxph(Surv(time_sem, dropout) 
                     ~ ScoresFaculty1
                     + ScoresFaculty2
                     + Origins
                     + HighschoolType 
                     + Income 
                     + HighschoolGrade 
                     + ECTS1sem,
                     data = students_COX_timdep)
summary(td.Cox.drop4)
logLik(td.Cox.drop4)
AIC(td.Cox.drop4)
BIC(td.Cox.drop4)
c22<-concordance(td.Cox.drop4)$concordance
c22

anova(td.Cox.drop4, td.Cox.drop2)

#_____________________

# Comparison
AIC(td.Cox.drop2, td.Cox.drop5, td.Cox.drop3, td.Cox.drop4)
anova(td.Cox.drop2, td.Cox.drop5, td.Cox.drop3, td.Cox.drop4)

#_____________________

# model with frailty (d of Table 5)

td.Cox.drop6<- coxph(Surv(time_sem, dropout) 
                     ~ frailty(School, distribution = 'gamma')
                     + Origins
                     + HighschoolType 
                     + Income 
                     + HighschoolGrade 
                     + ECTS1sem,
                     data = students_COX_timdep)
summary(td.Cox.drop6)
logLik(td.Cox.drop6)
BIC(td.Cox.drop6)
AIC(td.Cox.drop6)
c22<-concordance(td.Cox.drop6)$concordance
c22

anova(td.Cox.drop2, td.Cox.drop6)

#_____________________

# model with frailty (e of Table 5)

td.Cox.drop7<- coxph(Surv(time_sem, dropout) 
                     ~ frailty(Faculty, distribution = 'gamma')
                     + Origins
                     + HighschoolType 
                     + Income 
                     + HighschoolGrade 
                     + ECTS1sem,
                     data = students_COX_timdep)
summary(td.Cox.drop7)
logLik(td.Cox.drop7)
BIC(td.Cox.drop7)
AIC(td.Cox.drop7)
c22<-concordance(td.Cox.drop7)$concordance
c22



#_______________________________________________________________________________
#### Appendix B ####

students_COX_timdep = left_join(students_COX_timdep, TEST, 
                                by = join_by('School', 'Faculty'))

# with functional covariates as principal component scores
ber <- glm(dropout
           ~ ScoresSchool1 
           + ScoresFaculty1
           + ScoresSchool2 
           + Origins 
           + HighschoolType
           + Income 
           + HighschoolGrade 
           + ECTS1sem, 
           data = students_COX_timdep,
           family = binomial( link = logit ))
summary(ber)


AIC(ber)
BIC(ber)

# Brier score
brier <- mean((ber$fitted.values-ber$y)^2)
brier

# AIC
ber$aic

# K-fold cross-validation error
kcv <- 10
cverr <- cv.glm(data=students_COX_timdep, glmfit=ber, K=kcv)
cverr$delta 


#_____________________
# with dropout rate at course-level in the observation period
ber2 <- glm(dropout
            ~ ratio +
              + Origins 
            + HighschoolType
            + Income 
            + HighschoolGrade 
            + ECTS1sem, 
            data = students_COX_timdep,
            family = binomial( link = logit ))
summary(ber2)

AIC(ber2)
BIC(ber2)


# Brier score
brier <- mean((ber2$fitted.values-ber2$y)^2)
brier

# AIC
ber2$aic

# K-fold cross-validation error
kcv <- 10
cverr <- cv.glm(data=students_COX_timdep, glmfit=ber2, K=kcv)
cverr$delta 


anova(ber, ber2)




