#rm(list = ls()); library(vars); setwd("/home/alex/Documents/Uni/Fall_2017/TVAR_shrink/frf_TVAR/")
#source("../testing_R/est_var.R"); source("../testing_R/aux_functions.R")
#data("Canada"); est_res <- est_var(diff(Canada), p = 2, type = "const", nth = 1, info_crit = T)
#ff <- est_girf(est_res, sample_hist = F); ff2 <- est_girf(est_res)
#horizon <- 20; shocked_variable = "prod"; n_history = 100; n_boot = 100; shock_size <- 1; sample_hist = F
#sim_interval = c(.05, .95); model_choice = "BIC"; quiet = FALSE; seed = NULL
est_girf <- function(est_res, horizon = 20, shocked_variable = 1, shock_size = 1,  sample_hist = TRUE, 
                     n_history = 100, n_boot = 100, seed = NULL, sim_interval = c(.05, .95), 
                     model_choice = "BIC", quiet = FALSE) {
  #if sample_hist == TRUE: sample from histories: n_history-times
  #if sample_hist == FALSE: use all histories from regime: n_history is ignored
  
  if(is.character(shocked_variable)) { #if name of shocked variable is provided, not index number
    name_shocked_var <- shocked_variable
    shocked_variable <- which(colnames(est_res[[1]]$var_res$coefficients) == shocked_variable)
    if (length(shocked_variable) == 0) 
      stop("Variable `", name_shocked_var, "' not found. Variables in Model: ", 
           paste(colnames(est_res[[1]]$var_res$coefficients), collapse = ", "))
  }
  
  if (!quiet) cat(("GIRF\n====\nshocked variable: "), 
                  colnames(est_res[[1]]$var_res$coefficients)[shocked_variable],
                  "\nshock size: ", shock_size, "sd\n")
  if (!(is.null(seed))) set.seed(seed)
  #Selecting the correct model
  if (is.numeric(model_choice) & model_choice <= (length(names(est_res)) - 1)) {
    est_res <- est_res[[model_choice]]
  } else if (model_choice %in% c("AIC", "BIC")) {
    choose_model <- est_res$info_choice["Model", model_choice]
    est_res <- est_res[[choose_model]]
    if (!quiet) cat(paste0("Model ", choose_model, " selected by ", model_choice, ".\n"))
  } else {
    stop(paste0("Please select AIC or BIC or model number between 1 and ", (length(names(est_res)) - 1), ".\n"))
  }
  
  p <- est_res$var_info$p
  K <- est_res$var_info$K
  var_res <- est_res$var_res
  y_orig <- est_res$var_info$y_orig
  coeff_array <- est_res$var_res_add_info$coeff_array
  cov_matrix <- est_res$var_res_add_info$cov_matrix #regime dependent cov matrices
  regime_mask <- est_res$var_res_add_info$regime_mask
  thVar <- est_res$thresh_info$thVar
  thDelay <- est_res$thresh_info$thDelay
  number_regimes <- dim(coeff_array)[3]
  if (is.null(est_res$thresh_info$thValue)) {
    thValue <- est_res$var_res_add_info$best_thresh
  } else {
    thValue <- est_res$thresh_info$thValue
  }
  cholesky_res <- compute_chol(est_res)
  
  #find position of variables
  pos_lagged_vars <- which(substr(colnames(coeff_array[,,1]), 1,1) == "l")
  pos_fixed_vars <- which(colnames(coeff_array[,,1]) %in% c("const", "trend"))
  if (length(pos_fixed_vars) == 1) { #trend or const
    fixed_cols <- array(diag(K), dim = c(K, K, horizon + 1))
    if (colnames(coeff_array[,,1])[pos_fixed_vars] == "trend") {
      for (curr_trend in 1:(horizon+1))
        fixed_cols[,,curr_trend] <- fixed_cols[,,curr_trend] * (curr_trend+p)
    }
  } else if (length(pos_fixed_vars) == 2) { #trend and const
    fixed_cols <- array(0, dim = c(K, 2*K, horizon + 1))
    fixed_cols_const <- array(diag(K), dim = c(K, K, horizon + 1))
    fixed_cols_trend <- array(diag(K), dim = c(K, K, horizon + 1))
    for (curr_trend in 1:(horizon+1)) {
      fixed_cols_trend[,,curr_trend] <- fixed_cols_trend[,,curr_trend] * (curr_trend+p)
      fixed_cols[,,curr_trend] <- cbind(fixed_cols_const[,,curr_trend], fixed_cols_trend[,,curr_trend]) 
    }
  }
  
  #1. STEP PICK HISTORIES WHICH WILL BE USED FOR SIMULATION
  #================================================
  #correction for -1 and +p later on
  starting_histories <- list()
  for (i in 1:number_regimes) {
    possible_starts <- which(regime_mask[,i] == 1)
    if (sample_hist) {
      starting_histories[[i]] <- sample(possible_starts, size = n_history, replace = TRUE)
    } else {
      starting_histories[[i]] <- possible_starts
    } 
  }
  
  
  #SIMULATE OVER ALL REGIMES AND HISTORIES
  #=======================================
  #regime...regime, i...history, j...bootstrap rep
  if(sample_hist) {#find dimenstion for array
    dim_history <- n_history 
  } else {
    dim_history <- max(colSums(regime_mask))
  }
  girf_ij <- array(0, dim = c(number_regimes, dim_history, n_boot, horizon + 1, K))
  girf_i <- array(0, dim = c(number_regimes, dim_history, horizon + 1, K))
  if (!quiet) {
    pb <- txtProgressBar(min = 0, max = sum(n_boot * sapply(starting_histories, length)), style = 3)
    count_prog <- 1
  }
  #regime <- 1; i <- starting_histories[[regime]][1]; j<-1
  for (regime in 1:number_regimes) { #1. STEP: Pick Regime
    counter_for_history <- 1 #for saving girf
    for (i in starting_histories[[regime]]) {
      #2. STEP: pick starting history: reverse order: last period, period before last,...
      current_history <- current_history_s <- y_orig[(i + p - 1) : i, , drop = F]
      #figure out which regime we are in... 
      curr_regime <- curr_regime_s <- regime #updated for each step
      for (j in 1:n_boot) {
        #3. STEP: picking the errors
        err_pos <- sample(1:nrow(est_res$var_res$residuals), size = horizon + 1, replace = T)
        redform_err <- est_res$var_res$residuals[err_pos, ]
        #3a: structural errors for shocked and non-shocked
        struc_err <- t(cholesky_res$P_inv[[regime]] %*% t(redform_err))
        struc_err_shock <- struc_err
        struc_err_shock[1, shocked_variable] <- shock_size  #add or replace?
        #3b: transform back
        simul_err <- t(cholesky_res$P[[regime]] %*% t(struc_err))
        simul_err_shock <- t(cholesky_res$P[[regime]] %*% t(struc_err_shock))
        
        #4. STEP: simulate the systems
        simul_sys <- simul_sys_shock <- matrix(0, nrow = horizon + 1, ncol = K)
        for (step_sys in 1:(horizon+1)) {
          #adding constant and trend
          if (length(pos_fixed_vars) == 1) {
            simul_sys[step_sys, ] <- fixed_cols[,, step_sys] %*% 
              coeff_array[, pos_fixed_vars, curr_regime, drop = F]
            simul_sys_shock[step_sys, ] <- fixed_cols[,, step_sys] %*% 
              coeff_array[, pos_fixed_vars, curr_regime_s, drop = F]
          } else if (length(pos_fixed_vars == 2)) { #STILL TO TEST!!!!
            simul_sys[step_sys, ] <- fixed_cols[,1:K,step_sys] %*% 
              coeff_array[, pos_fixed_vars[1], curr_regime, drop = F] +
              fixed_cols[,(K+1):2*K,step_sys] %*% coeff_array[, pos_fixed_vars[2], curr_regime, drop = F]
            simul_sys_shock[step_sys, ] <- fixed_cols[,1:K,step_sys] %*% 
              coeff_array[, pos_fixed_vars[1], curr_regime_s, drop = F] +
              fixed_cols[,(K+1):2*K,step_sys] %*% coeff_array[, pos_fixed_vars[2], curr_regime_s, drop = F]
          }
          
          for (curr_lag in 1:p) {
            lag_position <- pos_lagged_vars[((curr_lag-1)*K + 1) : ((curr_lag-1)*K + K)]
            simul_sys[step_sys, ] <- simul_sys[step_sys, ] + coeff_array[, lag_position, curr_regime] %*%
              t(current_history[curr_lag, , drop = F])
            simul_sys_shock[step_sys, ] <- simul_sys_shock[step_sys, ] + coeff_array[, lag_position, curr_regime_s] %*%
              t(current_history_s[curr_lag, , drop = F])
          }
          #add shocks
          simul_sys[step_sys, ] <- simul_sys[step_sys, ] + simul_err[step_sys, ]
          simul_sys_shock[step_sys, ] <- simul_sys_shock[step_sys, ] + simul_err_shock[step_sys, ] 
          
          #if(any(is.na(simul_sys[step_sys, ]))) stop(step_sys, ", " ,regime, ",", i, ",", j)
          
          #figure out new regimes: update history, then find history, then delete last obs in history
          current_history <- rbind(simul_sys[step_sys, ], current_history)
          current_history_s <- rbind(simul_sys_shock[step_sys, ], current_history_s)
          
          curr_regime <- sum(current_history[thDelay, thVar] > thValue) + 1
          curr_regime_s <- sum(current_history[thDelay, thVar] > thValue) + 1
          
          #update history
          current_history <- current_history[-nrow(current_history), , drop = F]
          current_history_s <- current_history_s[-nrow(current_history_s), , drop = F]
          
        }
        girf_ij[regime, counter_for_history, j,,] <- simul_sys_shock - simul_sys
        if (!quiet) {
          setTxtProgressBar(pb, count_prog)
          count_prog <- count_prog + 1
        }
      }
      curr_girf <- girf_ij[regime, counter_for_history,,,]
      girf_i[regime, counter_for_history,,] <- apply(curr_girf, MARGIN = c(2,3), mean) #mean of bootstrap rep
      counter_for_history <- counter_for_history + 1
    }
  }
  #compute girf and return
  girf <- apply(girf_i, MARGIN = c(1,3,4), mean) #mean over histories
  
  if (!sample_hist) { #correcting for differnt lengths of histories if we are not sampling
    n_history <- sapply(starting_histories, length)
    for (i in 1:number_regimes) girf[i,,] <- girf[i,,] * dim_history / n_history[i]
  }
  
  ret_list <- list(girf = girf, girf_i = girf_i, girf_ij = girf_ij)
  return(ret_list)
}
