# https://www.math.fsu.edu/~ychen/research/Thinning%20algorithm.pdf
# https://stats.stackexchange.com/questions/369288/nonhomogeneous-poisson-process-simulation/369294#369294

# vedere a livello ospedale un trend
# quando traslo continua a funzionare?
# aggiungere l'età (exp(beta age) con age gaussiana standard)
# aumentare punti?
# vedere lambda piccoli stimati

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
library(RColorBrewer)


genNHPP <- function(lambda, T, n = 1, m = NULL, j = NULL) {
  # Generate a NHPP 
  MaxLambda <- max(lambda(seq(0, T, by = 10^-5), m, j))
  ph <- function(t) lambda(t, m=m, j=j) / MaxLambda
  
  if (n == 1) {
    # Single Sample Path
    t <- 0
    Nevents <- 0
    EventTimes <- numeric(0)
    done <- FALSE
    while (!done) {
      t <- t + (-1 / MaxLambda) * log(runif(1))
      if (t > T) {
        done <- TRUE
      } else {
        if (runif(1) <= ph(t)) {
          Nevents <- Nevents + 1
          EventTimes[Nevents] <- t
        }
      }
    }
  } 
  return(EventTimes)
}

# Example usage:
rate_function <- function(t, m = NULL, j = NULL) { return(30*(1 + sin(t))) }
Lambda = function(x){
  return(Vectorize(function(x) integrate(rate_function, lower = 0, upper = x)$value)(x))
}
result <- genNHPP(rate_function, T = 2*pi, m = NULL, j = NULL)
zero_function <- function(t) { return(t*0) }

plot(seq(0, 2*pi, by = 0.01), rate_function(seq(0, 2*pi, by = 0.01)), type='l')
points(result, zero_function(result), col='red', pch=16)

plot(seq(0, 2*pi, by = 0.01), Lambda(seq(0, 2*pi, by = 0.01)), type='l')
points(result, zero_function(result), col='red', pch=16)


#########

set.seed(1) # 123
I = 20 # groups 
J = 4 # units
L = 1
level = 0.1
sigma = 0
balanced = TRUE
#set.seed(2)
#age = sample(seq(from = 0.5, to = 1.5, by = 0.1), size = I*J, replace = TRUE)
#age = rnorm(I*J, 1, 0.4)
# Create an empty data frame to store the simulated times
simulated_data = data.frame()


# Eigenfunctions

# ef11 <- function(x) {
#   #return(sqrt(2) * sin(1 * x * 2 * pi))
#   return(exp(-x))
# }
# ef12 <- function(x) {
#   return(exp(x))
#   # return(exp(x))
# }
# ef13 <- function(x) {
#   #return(sqrt(2) * sin(1 * x * 4 * pi))
#   return(exp(x^2))
# }
# ef14 <- function(x) {
#   #return(sqrt(2) * cos(1 * x * 4 * pi))
#   return(exp(x^3))
# }

ef11 <- function(x) {
  return(sqrt(2) * sin(1 * x * 2 * pi))
}
ef12 <- function(x) {
  return(sqrt(2) * cos(1 * x * 2 * pi))
}
ef13 <- function(x) {
  return(sqrt(2) * sin(1 * x * 4 * pi))
}
ef14 <- function(x) {
  return(sqrt(2) * cos(1 * x * 4 * pi))
}

ef21 <- function(x) {
  return(sqrt(1) * x^0)
}
ef22 <- function(x) {
  return(sqrt(3) * (2 * x - 1))
}
ef23 <- function(x) {
  return(sqrt(5) * (6*x^2 - 6 * x + 1))
}
ef24 <- function(x) {
  return(sqrt(7) * (20*x^3 - 30*x^2 + 12 * x - 1))
}


zero_function = function(x){
  return(0*x)
}


# Generate scores
if (balanced == FALSE) {
  J_subj <- pmax(rpois(I, J), 1)
} else {
  J_subj <- rep(J, I)
}
n <- sum(J_subj)

lambdas = NULL
Lambdas = NULL
J_ind <- c(0, cumsum(J_subj))

K1 <- 4
K2 <- 4
K <- K1 + K2
lambda1 <- 0.9^(0:(K1-1))
lambda2 <- 0.2^(0:(K2-1))

# Generate scores
si1 <- matrix(0, nrow = I, ncol = K1)
si2 <- matrix(0, nrow = n, ncol = K2)
for (k in 1:K1) {
  si1[, k] <- rnorm(I, sd = sqrt(lambda1[k]))
}
for (k in 1:K2) {
  si2[, k] <- rnorm(n, sd = sqrt(lambda2[k]))
}

for (m in 1:I) {
  eff_I = function(tt, m){
    return(si1[m, 1] * ef11(tt) + si1[m, 2] * ef12(tt) + 
             si1[m, 3] * ef13(tt) + si1[m, 4] * ef14(tt)
    )
  }
  #plot(seq(0, L, by=0.01), eff_I(seq(0, L, by=0.01), m), 
  #      type='l', ylim=c(-5,5), main = paste('I effect', m, sep=' '))
  # print(eff_I(seq(0, L, by=0.01)))
  
  for (j in 1:J_subj[m]) {
    eff_J = function(tt, m, j){
      return(si2[J_ind[m]+j, 1] * ef21(tt) + si2[J_ind[m]+j, 2] * ef22(tt) + 
               si2[J_ind[m]+j, 3] * ef23(tt) + si2[J_ind[m]+j, 4] * ef24(tt)
      )
      
    }
    #plot(seq(0, L, by=0.01), eff_J(seq(0, L, by=0.01), m, j), 
    #     type='l', ylim=c(-10,10), main = paste('J effect', j, sep=' '))
    #print(eff_J(seq(0, L, by=0.01)))
    
    lambda_ = function(tt, m=m, j=j) {
      return(200 + 2*m*(eff_I(tt, m) + eff_J(tt, m, j)))
    }
    
    #plot(seq(0, L, by=0.01), lambda_(seq(0, L, by=0.01), m, j), 
    #     type='l', ylim=c(0,20),  main = paste('Y', J_ind[m]+j, sep=' '))
    
    lambdas = rbind(lambdas, lambda_(seq(0, L, by=0.01), m, j))
    
    simulated_times = genNHPP(lambda_, L, m = m, j = j)
    #print(simulated_times)
    
    #plot(seq(0, L, by=0.01), lambda_(seq(0, L, by=0.01)), type='l', ylim=c(0,20))
    #points(simulated_times, zero_function(simulated_times), col='red', pch=16)
    
    Lambda_ = function(x){
      return(Vectorize(function(x) integrate(lambda_, lower = 0, upper = x, m = m, j = j)$value)(x))
    }
    Lambdas = rbind(Lambdas, Lambda_(seq(0, L, by=0.01)))
    #plot(seq(0, L, by=0.01), Lambda_(seq(0, L, by=0.01)), type='l')
    #points(simulated_times, zero_function(simulated_times), col='red', pch=16)
    
    
    # Create a data frame for the simulated times
    simulated_patient_data = data.frame(
      centre = rep(m, length(simulated_times)),
      id = rep(j, length(simulated_times)),
      time = simulated_times
      # age_in = rep(age[j], length(simulated_times))
    )
    ## Add it to the existing data
    simulated_data = rbind(simulated_data, simulated_patient_data)
  }    
}


df = data.frame('id'= rep(1:(I*J), each=101), 
                'time' = seq(0, L, by=0.01),
                'intensity' = c(t(lambdas)),
                'cumhaz' = c(t(Lambdas)),
                'centre'= rep(1:I, each=101*J))

selected_centers <-c(14, 17, 20) # c(1, 16, 17, 18)  # Add your selected centers here # c(1, 14, 17, 20)

source('PlotsSimulation.R')
PLOTLambda(df, selected_centres, 'intensity')
PLOTLambda(df, selected_centres, 'cumhaz')


# LET'S DO FPCA on lambdas to check coherence
res <- mfpca.face(Y = lambdas, id = rep(1:I, each=J),
                  twoway = FALSE)
# mfpca.face(Y = xx$Y, id =xx$id, twoway = FALSE)
res$evalues
res$npc

# matplot(res$efunctions$level1[,1], type='l', ylim = c(-3, 3))
# matplot(res$efunctions$level1[,2], type='l', ylim = c(-3, 3))
# matplot(res$efunctions$level2[,1], type='l', ylim = c(-3, 3))
# matplot(res$efunctions$level2[,2], type='l', ylim = c(-3, 3))


# LET'S DO FPCA on Lambdas 
res <- mfpca.face(Y = Lambdas, id = rep(1:I, each=J),
                twoway = FALSE)
# mfpca.face(Y = xx$Y, id =xx$id, twoway = FALSE)
res$evalues
res$npc

#matplot(res$efunctions$level1[,1], type='l', ylim = c(-3, 3))
#matplot(res$efunctions$level1[,2], type='l', ylim = c(-3, 3))
#matplot(res$efunctions$level2[,1], type='l', ylim = c(-3, 3))
#matplot(res$efunctions$level2[,2], type='l', ylim = c(-3, 3))


# Assuming 'res' is your data frame
COEFF = 3

# Create a data frame for ggplot
plot_data <- data.frame(x = seq_along(res$mu)/100, 
                        y = res$mu, 
                        C1L1 = res$efunctions$level1[,1],
                        C2L1 = res$efunctions$level1[,2],
                        C1L2 = res$efunctions$level2[,1],
                        C2L2 = res$efunctions$level2[,2],
                        posC1L1 = res$mu + res$efunctions$level1[, 1] * COEFF * sqrt(res$evalues$level1[1]), 
                        negC1L1 = res$mu - res$efunctions$level1[, 1] * COEFF * sqrt(res$evalues$level1[1]),
                        posC2L1 = res$mu + res$efunctions$level1[, 2] * COEFF * sqrt(res$evalues$level1[1]), 
                        negC2L1 = res$mu - res$efunctions$level1[, 2] * COEFF * sqrt(res$evalues$level1[1]),
                        posC1L2 = res$mu + res$efunctions$level2[, 1] * COEFF * sqrt(res$evalues$level2[1]),
                        negC1L2 = res$mu - res$efunctions$level2[, 1] * COEFF * sqrt(res$evalues$level2[1]),
                        posC2L2 = res$mu + res$efunctions$level2[, 2] * COEFF * sqrt(res$evalues$level2[1]),
                        negC2L2 = res$mu - res$efunctions$level2[, 2] * COEFF * sqrt(res$evalues$level2[1]))

PlotEigenFunMFPCA(plot_data, level = 1)
PlotEigenFunMFPCA(plot_data, level = 2)

PlotCompMFPCA(plot_data, level = 1, component = 1, coeff = COEFF)
PlotCompMFPCA(plot_data, level = 1, component = 2, coeff = COEFF)
PlotCompMFPCA(plot_data, level = 2, component = 1, coeff = COEFF)
PlotCompMFPCA(plot_data, level = 2, component = 2, coeff = COEFF)


# Adding a last row with time = 1.001 for each combination of id and centre
simulated_data <- simulated_data %>%
  group_by(id, centre) %>%
  do(bind_rows(., tibble(id = first(.$id), centre = first(.$centre), time = 1.001))) %>%
  ungroup()

simulated_data <- simulated_data %>%
  distinct()

# Add the "event" variable
simulated_data <- simulated_data %>%
  group_by(id, centre) %>%
  mutate(enum = row_number() - 1,
         timecum = time,
         event = ifelse(timecum == max(timecum), 0, 1),
         time = ifelse(enum == 0, timecum, timecum - lag(timecum, n=1))
  ) %>%
  ungroup() %>%
  mutate(id = dense_rank(id) + max(id) * (centre - 1) )

simulated_data <- simulated_data %>%
  group_by(id, centre) %>%
  mutate(logp1_enum = log(enum+1),
         start = ifelse(enum == 0, 0,  lag(timecum, n=1)),
         stop = timecum
  ) %>%
  ungroup() 

simulated_data <- simulated_data %>%
  distinct()

# Print the modified dataframe
print(simulated_data)



simulated_data$time = NULL

simulated_data <- simulated_data %>%
  mutate(timecum = 1000 * timecum,
         start = floor(1000 * start),
         stop = floor(1000 * stop),
         start = as.numeric(ifelse(start == 0, 0.0, start)),
         stop = as.numeric(ifelse(stop == 1000, 1000.5, stop)))

# Check if stop > start and create a new column 'valid'
simulated_data <- simulated_data %>%
  mutate(valid = stop > start)

simulated_data <- simulated_data %>%
  filter(valid == TRUE)

#save(simulated_data,file="simulated_data.Rdata")

t_max <- 1000
times <- seq(0,t_max,by=1)

count_model = coxph(Surv(start,stop,event) ~ 1 + enum, id=id, cluster=centre, data = simulated_data)
summary(count_model)
plot(survfit(count_model))

# Fit and smooth cumulative baseline hazard; return a list
fit_smooth_Lambda0 <- function(model){
  
  # Get estimated cumulative baseline hazard function
  bh = basehaz(model, centered = FALSE)
  t <- bh$time
  Lambda0 <- bh$hazard
  
  # Smooth version of Lambda0
  Lambda0S <- cobs(c(0,t), c(0,Lambda0), constraint=c("increase"), 
                   pointwise=matrix(c(0,0,0),nrow=1), nknots=20, lambda=0, toler.kn=0)
  
  return(list('times0' = t,
              'Lambda0' = Lambda0,
              'Lambda0S' = Lambda0S)
  )
}



# list_L0 = bashaz
bashaz <- fit_smooth_Lambda0(count_model)
plot_Lambda0(bashaz)


# Auxiliary functions for file 02_modelling_compensators.R
# Step 2.3: Reconstruct compensators

# Evaluate cumulative hazard in a grid of points; returns a dataframe in long format
compute_cumulative_hazard = function(model,
                                     sel_df,
                                     smoothed_baseline,
                                     times,
                                     verbose = FALSE){
  # Evaluate smoothed Lambda_0 on a grid
  Lambda0s_value = .Lambda0_fun(times,smoothed_baseline)
  # Compute constant at times coefficients
  patient_coefficients = .compute_coefficients_ck(sel_df,model, verbose)
  # Compute daily deltas of cumulative hazard and sum it up
  cumulative_hazard = .compute_cumulative_hazard(patient_coefficients,Lambda0s_value,times,verbose)
  return(cumulative_hazard)
}


# Compute coefficient ck = exp{beta x_i(t_{i,j-1}) + gamma z_i(t_{i,j-1})} at each interval for each patient 
.compute_coefficients_ck = function(sel_df, model,verbose = FALSE){
  if(verbose){print('Computing coefficients')}
  name_coefficients = names(model$coefficients)
  coefficients = model$coefficients
  
  ck <- exp(as.matrix(sel_df[,..name_coefficients])%*%coefficients)
  patient_coefficients <- cbind(sel_df[,.(id,enum,start,stop,centre)],ck)
  colnames(patient_coefficients)[6] <- c('ck')
  
  return(data.table(patient_coefficients))
}


# Compute daily deltas of cumulative Hazard and it sums it up; returns a dataframe in long format 
.compute_cumulative_hazard = function(patient_coefficients,Lambda0s_value,times,verbose = FALSE){
  if(verbose){print('Computing cumulative Hazard on the grid')}
  
  # Transform coefficients data in long format
  long_pc <- patient_coefficients %>%
    mutate(time = map2(start, stop, seq, by=0.5)) %>%
    unnest(cols = time) %>%
    dplyr::select(-start, -stop) %>%
    group_by(time)
  long_pc<-data.table(long_pc)
  long_pc[time!=0.5, time := floor(time)]
  long_pc<-long_pc[!duplicated(long_pc[,.(id,time)], fromLast=F)]
  # Delete t=0.5 if not in Lambda0s grid
  if(!(0.5 %in% times)){
    long_pc<-long_pc[time!=0.5]
  }
  # Temp data for easy deltas computation
  temp_Lj <- data.frame('L0j' = Lambda0s_value, 'time' = times )
  temp_Ljm1 <- data.frame('L0jm1' = Lambda0s_value, 'time' = times+1 )
  tmp_data <- merge(long_pc,temp_Lj,by='time')
  tmp_data <- merge(tmp_data,temp_Ljm1,by='time',all.x=T)
  tmp_data <- tmp_data[order(id,time)]
  # t=0
  tmp_data[time==0, deltas := ck*L0j]
  # t>0
  tmp_data[time>0, deltas := ck*(L0j-L0jm1)]
  tmp_data[, cumhaz := round(cumsum(deltas),8), by = id]
  
  cumulative_hazard = tmp_data[,.(id,centre,time,cumhaz)]
  
  return(cumulative_hazard)
}


# Evaluate smoothed baseline cumulative hazard
.Lambda0_fun <- function(t,smoothed_baseline){
  return(.basis(t, smoothed_baseline$knots) %*% smoothed_baseline$coef)
}


# Evaluate a B-spline basis on a vector of points x
.basis <- function(x, knots, deg=2){
  return(bs(x,
            knots=knots[2:(length(knots)-1)],
            degree=deg,
            Boundary.knots=c(knots[1],knots[length(knots)]),
            intercept=TRUE))
}


# model = count_model
# sel_df = as.data.table(simulated_data)
# verbose = FALSE
# smoothed_baseline = bashaz$Lambda0S
# Lambda0s_value = .Lambda0_fun(times,smoothed_baseline)
# # Compute constant at times coefficients
# patient_coefficients = .compute_coefficients_ck(sel_df,model, verbose)


cumulative_hazard = compute_cumulative_hazard(model = count_model,
                                              sel_df = as.data.table(simulated_data),
                                              smoothed_baseline = bashaz$Lambda0S,
                                              times, verbose = TRUE)
cumulative_hazard = as.data.frame(cumulative_hazard)
cumulative_hazard$time = cumulative_hazard$time / 1000

PLOTLambda(cumulative_hazard, selected_centres, 'cumhaz', add_hat = TRUE)


reformat_cumhaz = function(cumulative_hazard){
  patient_ids = unique(cumulative_hazard$id)
  np <- length(patient_ids) # number of rows: number of patients
  m = length(unique(cumulative_hazard$time))
  evals_matrix = spread(cumulative_hazard, time, cumhaz)
  # isolate numeric covariates
  num_matrix = as.matrix(evals_matrix[,-c(1,2)]) 
  df = data.frame(evals_matrix[,c(1,2)], cumhaz = rep(0, nrow(num_matrix)))
  df$cumhaz = num_matrix
  
  return(df)
}

df = reformat_cumhaz(cumulative_hazard)

# cumhaz_der = matrix(0, nrow=dim(df$cumhaz)[1], ncol=dim(df$cumhaz)[2])
# 
# library(pspline)
# for(i in 1:dim(df$cumhaz)[1]){
#   cumhaz_der[i,] <- predict(sm.spline(1:366, df$cumhaz[i,]), 1:366, 1) 
# }
# 
# df$cumhaz_der = cumhaz_der
# matplot(times, t(df$cumhaz_der), type='l')

res <- mfpca.face(Y = df$cumhaz, 
                  id = df$centre,
                  pve=0.997,
                  visit = df$id, 
                  twoway = FALSE)

res$evalues
# Assuming 'res' is your data frame

COEFF = 3

# Create a data frame for ggplot
plot_data <- data.frame(x = seq_along(res$mu)/1000, 
                        y = res$mu, 
                        C1L1 = res$efunctions$level1[,1],
                        C2L1 = res$efunctions$level1[,2],
                        C1L2 = res$efunctions$level2[,1],
                        C2L2 = res$efunctions$level2[,2],
                        posC1L1 = res$mu + res$efunctions$level1[, 1] * COEFF * sqrt(res$evalues$level1[1]), 
                        negC1L1 = res$mu - res$efunctions$level1[, 1] * COEFF * sqrt(res$evalues$level1[1]),
                        posC2L1 = res$mu + res$efunctions$level1[, 2] * COEFF * sqrt(res$evalues$level1[1]), 
                        negC2L1 = res$mu - res$efunctions$level1[, 2] * COEFF * sqrt(res$evalues$level1[1]),
                        posC1L2 = res$mu + res$efunctions$level2[, 1] * COEFF * sqrt(res$evalues$level2[1]),
                        negC1L2 = res$mu - res$efunctions$level2[, 1] * COEFF * sqrt(res$evalues$level2[1]),
                        posC2L2 = res$mu + res$efunctions$level2[, 2] * COEFF * sqrt(res$evalues$level2[1]),
                        negC2L2 = res$mu - res$efunctions$level2[, 2] * COEFF * sqrt(res$evalues$level2[1]))

PlotEigenFunMFPCA(plot_data, level = 1, post = TRUE)
PlotEigenFunMFPCA(plot_data, level = 2, post = TRUE)

PlotCompMFPCA(plot_data, level = 1, component = 1, post = TRUE, coeff = COEFF)
PlotCompMFPCA(plot_data, level = 1, component = 2, post = TRUE, coeff = COEFF)
PlotCompMFPCA(plot_data, level = 2, component = 1, post = TRUE, coeff = COEFF)
PlotCompMFPCA(plot_data, level = 2, component = 2, post = TRUE, coeff = COEFF)




# mfpca.2wayF$npc
# 
# mfpca.2wayF$mu
# mfpca.2wayF$Xhat
# mfpca.2wayF$Xhat.subject
# 
# mfpca.2wayF$scores
# mfpca.2wayF$efunctions
# mfpca.2wayF$evalues
# 
# matplot(times, t(mfpca.2wayF$Xhat), type='l')
# 
# matplot(times, t(mfpca.2wayF$Xhat.subject), type='l') #-mfpca.2wayF$mu
# 
# matplot(times, mfpca.2wayF$mu, type='l')
# 
# # mean
# plot(mfpca.2wayF$mu)
# 
# # FPC1 LEVEL 1
# plot(mfpca.2wayF$efunctions$level1[,1], type='l')
# # FPC2 LEVEL 1
# plot(mfpca.2wayF$efunctions$level1[,2], type='l')
# 
# plot(mfpca.2wayF$mu,lwd=2,main='FPC1', type='l')
# lines(mfpca.2wayF$mu + mfpca.2wayF$efunctions$level1[,1] * 50*sqrt(mfpca.2wayF$evalues$level1[1]), col=2)
# lines(mfpca.2wayF$mu - mfpca.2wayF$efunctions$level1[,1] * 50*sqrt(mfpca.2wayF$evalues$level1[1]), col=3)
# 
# plot(mfpca.2wayF$mu,lwd=2,main='FPC2', type='l')
# lines(mfpca.2wayF$mu + mfpca.2wayF$efunctions$level1[,2] * 50*sqrt(mfpca.2wayF$evalues$level1[2]), col=2)
# lines(mfpca.2wayF$mu - mfpca.2wayF$efunctions$level1[,2] * 50*sqrt(mfpca.2wayF$evalues$level1[2]), col=3)
# 
# # FPC1 LEVEL 2
# plot(mfpca.2wayF$efunctions$level2[,1], type='l')
# # FPC2 LEVEL 2
# plot(mfpca.2wayF$efunctions$level2[,2], type='l')
# 
# plot(mfpca.2wayF$mu,lwd=2,main='FPC1', type='l')
# lines(mfpca.2wayF$mu + mfpca.2wayF$efunctions$level2[,1] * 50*sqrt(mfpca.2wayF$evalues$level2[1]), col=2)
# lines(mfpca.2wayF$mu - mfpca.2wayF$efunctions$level2[,1] * 50*sqrt(mfpca.2wayF$evalues$level2[1]), col=3)
# 
# plot(mfpca.2wayF$mu,lwd=2,main='FPC2', type='l')
# lines(mfpca.2wayF$mu + mfpca.2wayF$efunctions$level2[,2] * 50*sqrt(mfpca.2wayF$evalues$level2[2]), col=2)
# lines(mfpca.2wayF$mu - mfpca.2wayF$efunctions$level2[,2] * 50*sqrt(mfpca.2wayF$evalues$level2[2]), col=3)
# 


