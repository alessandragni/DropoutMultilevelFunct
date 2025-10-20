simulate_once <- function(seed, I = 20, J = 4, grid = 1000, 
                          balanced = TRUE, verbose = FALSE) {
  
  source("utils.R")
  
  set.seed(1)
  I = I # groups 20
  J = J # units 4
  L = 1
  level = 0.1
  sigma = 0
  balanced = balanced # TRUE
  
  
  # Create an empty data frame to store the simulated times
  simulated_data = data.frame()
  
  
  # Eigenfunctions
  
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
  
  # start simulation changing the seed
  set.seed(seed)
  for (m in 1:I) {
    eff_I = function(tt, m){
      return(si1[m, 1] * ef11(tt) + si1[m, 2] * ef12(tt) + 
               si1[m, 3] * ef13(tt) + si1[m, 4] * ef14(tt)
      )
    }
    
    for (j in 1:J_subj[m]) {
      eff_J = function(tt, m, j){
        return(si2[J_ind[m]+j, 1] * ef21(tt) + si2[J_ind[m]+j, 2] * ef22(tt) + 
                 si2[J_ind[m]+j, 3] * ef23(tt) + si2[J_ind[m]+j, 4] * ef24(tt)
        )
        
      }
      
      lambda_ = function(tt, m=m, j=j) {
        return(200 + 2*m*(eff_I(tt, m) + eff_J(tt, m, j)))
      }
      
      lambdas = rbind(lambdas, lambda_(seq(0, L, by=(1/grid)), m, j))
      
      simulated_times = genNHPP(lambda_, L, m = m, j = j)
      
      Lambda_ = function(x){
        return(Vectorize(function(x) integrate(lambda_, lower = 0, upper = x, m = m, j = j)$value)(x))
      }
      Lambdas = rbind(Lambdas, Lambda_(seq(0, L, by=(1/grid))))
      
      
      # Create a data frame for the simulated times
      simulated_patient_data = data.frame(
        centre = rep(m, length(simulated_times)),
        id = rep(j, length(simulated_times)),
        time = simulated_times
      )
      # Add it to the existing data
      simulated_data = rbind(simulated_data, simulated_patient_data)
    }    
  }
  
  
  df <- data.frame(
    id = rep(1:n, each = grid + 1),
    time = rep(seq(0, L, by = 1 / grid), times = n),
    intensity = c(t(lambdas)),
    cumhaz = c(t(Lambdas)),
    centre = rep(rep(1:I, times = J_subj), each = grid + 1)
  )
  
  if(verbose){
    selected_centres <-c(14, 17, 20) # c(1, 16, 17, 18) 
    
    source('PlotsSimulation.R')
    PLOTLambda(df, selected_centres, 'intensity')
    PLOTLambda(df, selected_centres, 'cumhaz')
    
  }
  
  # LET'S DO FPCA on Lambdas 
  res <- mfpca.face(Y = Lambdas, 
                    # id = rep(1:I, each=J),
                    id = rep(1:I, times = J_subj),
                    twoway = FALSE)
  res0 = res
  res$evalues
  res$npc
  
  if(verbose){
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
    
    PlotEigenFunMFPCA(plot_data, level = 1)
    PlotEigenFunMFPCA(plot_data, level = 2)
    
    PlotCompMFPCA(plot_data, level = 1, component = 1, coeff = COEFF)
    PlotCompMFPCA(plot_data, level = 1, component = 2, coeff = COEFF)
    PlotCompMFPCA(plot_data, level = 2, component = 1, coeff = COEFF)
    PlotCompMFPCA(plot_data, level = 2, component = 2, coeff = COEFF)
  }
  
  
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
  
  simulated_data$time = NULL
  
  simulated_data <- simulated_data %>%
    mutate(timecum = grid * timecum,
           start = floor(grid * start),
           stop = floor(grid * stop),
           start = as.numeric(ifelse(start == 0, 0.0, start)),
           stop = as.numeric(ifelse(stop == grid, (grid+0.5), stop)))
  
  # Check if stop > start and create a new column 'valid'
  simulated_data <- simulated_data %>%
    mutate(valid = stop > start)
  
  simulated_data <- simulated_data %>%
    filter(valid == TRUE)
  
  #save(simulated_data,file="simulated_data.Rdata")
  
  t_max <- grid
  times <- seq(0,t_max,by=1)
  
  count_model = coxph(Surv(start,stop,event) ~ 1 + enum, id=id, cluster=id, data = simulated_data)
  #count_model = coxph(Surv(start,stop,event) ~ 1 + enum, id=id, cluster=centre, data = simulated_data)
  summary(count_model)
  
  
  bashaz <- fit_smooth_Lambda0(count_model)
  
  cumulative_hazard = compute_cumulative_hazard(model = count_model,
                                                sel_df = as.data.table(simulated_data),
                                                smoothed_baseline = bashaz$Lambda0S,
                                                times, verbose = FALSE)
  cumulative_hazard = as.data.frame(cumulative_hazard)
  cumulative_hazard$time = cumulative_hazard$time / grid
  
  if(verbose){
    PLOTLambda(cumulative_hazard, selected_centres, 'cumhaz', add_hat = TRUE)
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
  
  df = reformat_cumhaz(cumulative_hazard)
  
  
  res <- mfpca.face(Y = df$cumhaz, 
                    id = df$centre,
                    pve=0.9999,
                    visit = df$id, 
                    twoway = FALSE)
  
  res$evalues
  # Assuming 'res' is your data frame
  
  if(verbose){
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
  }
  
  
  # Example for your levels
  MISE2_eigen1 <- compute_MISE(res0$efunctions[[1]], res$efunctions[[1]], grid)
  MISE2_eigen2 <- compute_MISE(res0$efunctions[[2]], res$efunctions[[2]], grid)
  
  
  # Return key results
  return(list(
    evalues_0 = res0$evalues,  
    evalues_est = res$evalues,
    efuncs_0 = res0$efunctions, 
    efuncs_est = res$efunctions,
    MISE1 = MISE2_eigen1,
    MISE2 = MISE2_eigen2
  ))
}
