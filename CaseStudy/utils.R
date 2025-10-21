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