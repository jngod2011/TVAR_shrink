#source("../testing_R/est_var.R")
#data("Canada"); est_res <- est_var(diff(Canada), p = 2, type = "const", nth = 1)
#ff <- est_girf(est_res)
#horizon <- 20; shocked_variable = 1; n_history = 100; n_boot = 100; shock_size <- 1; sim_interval = c(.05, .95)
est_girf <- function(est_res, horizon = 20, shocked_variable = 1, shock_size = 1, n_history = 100, n_boot = 100,
                     seed = NULL, sim_interval = c(.05, .95)) {
  if (!(is.null(seed))) set.seed(seed)
  cl <- match.call()
  p <- est_res$var_info$p
  K <- est_res$var_info$K
  var_res <- est_res$var_res
  y_orig <- est_res$var_info$y_orig
  coeff_array <- est_res$var_res_add_info$coeff_array
  cov_matrix <- est_res$var_res_add_info$cov_matrix
  regime_mask <- est_res$var_res_add_info$regime_mask
  resids_by_regime <- get_regime_residuals(est_res)
  thVar <- est_res$thresh_info$thVar
  thDelay <- est_res$thresh_info$thDelay
  no_regimes <- dim(coeff_array)[3]
  if (is.null(est_res$thresh_info$thValue)) {
    thValue <- est_res$var_res_add_info$best_thresh
  } else {
    thValue <- est_res$thresh_info$thValue
  }
 
  #find position of variables
  pos_lagged_vars <- which(substr(colnames(coeff_array[,,1]), 1,1) == "l")
  pos_fixed_vars <- which(colnames(coeff_array[,,1]) %in% c("const", "trend"))
  pos_exog_vars <- which(substr(colnames(coeff_array[,,1]), 1,4)  == "exog")
  
  #cholesky decompose and structural errors
  chol_decomp <- compute_chol(est_res)
  
  all_girf <- all_low_girf <- all_high_girf <- list()
  dimnames_raw_girf <- list(1:(horizon+1), colnames(y_orig), paste0("regime_", 1:no_regimes),
                            paste0("history_", 1:n_history), paste0("boot_rep", 1:n_boot))
  raw_girf <- array(0, dim = c(horizon + 1, K, no_regimes, n_history, n_boot), 
                    dimnames = dimnames_raw_girf)
  for (reg in 1:no_regimes) { #regimes
    #for starting value (later corrected with -1 to simulate for period 0 and +p for y_orig)
    #contains all possible starting values
    pos_observ <- which(regime_mask[, reg] == 1)
    
    girf_i <- array(0, dim = c(horizon + 1, K, n_history)) #i...history
    for (i in 1:n_history) {
      sys_non_shock <- sys_shock <- matrix(0, nrow = horizon + 1, ncol = K)
      girf_ij <- array(0, dim = c(horizon + 1, K, n_boot)) #i...history, j...boot
      
      #1. step: starting history
      #=========================
      # +p allows for y_orig in next step; -1 to simulate for period 0!
      hist_pos <- sample(pos_observ, size = 1) + p - 1 
      start_hist <- y_orig[hist_pos:(hist_pos - p + 1),, drop = F] 
      
      regime_coeff <- coeff_array[,, reg] #we start in this regime
      pos_boot_err <- sample(1:nrow(chol_decomp$str_err[[reg]]), size = 1)
      first_err_non_shock <- chol_decomp$str_err[[reg]][pos_boot_err,, drop = F] #first error
      first_err_shock <- first_err_non_shock
      first_err_shock[, shocked_variable] <- first_err_shock[, shocked_variable] + shock_size
      
      coeff_fixed <- regime_coeff[, pos_fixed_vars, drop = F]
      fixed_cols <- matrix(rep(1, ncol(coeff_fixed)))
      sys_non_shock[1,] <- sys_shock[1,] <-  coeff_fixed %*% fixed_cols
      for (k in 1:p) {
         col_pos <- pos_lagged_vars[((k-1)*K + 1) : ((k-1)*K + K)]
         coeff_matrix <- regime_coeff[, col_pos] 
         
         sys_non_shock[1, ] <- sys_non_shock[1, ] + coeff_matrix %*% t(start_hist[k,, drop = F])
         sys_shock[1, ] <- sys_shock[1, ] + coeff_matrix %*% t(start_hist[k,, drop = F])
      }
      sys_non_shock[1, ] <- sys_non_shock[1, ] + chol_decomp$P[[reg]] %*% t(first_err_non_shock)
      sys_shock[1, ] <- sys_shock[1, ] + chol_decomp$P[[reg]] %*% t(first_err_shock)
      
     
      for (j in 1:n_boot) {
        #we only keep starting value; next line is not really necessary...
        sys_non_shock[2:nrow(sys_non_shock), ] <- sys_shock[2:nrow(sys_shock), ] <- 0 
        #start histories
        curr_hist_non_shock <- rbind(sys_non_shock[1,], start_hist[-p,])
        curr_hist_shock <- rbind(sys_shock[1,], start_hist[-p,])
        for (curr_step in 2:(horizon+1)) {
          #2. step: clarify which regime we are in for both histories
          #===========================================================
          pos_curr_thresh <- c(curr_step - thDelay, thVar) #finding position of current thresh value
          if (pos_curr_thresh[1] < -(p-1)) {
            stop("Threshold Delay is too large (max Delay: p+1)")
          } else if (pos_curr_thresh[1] <= 0) {#take from starting history
            pos_curr_thresh[1] <- pos_curr_thresh[1] * (-1) + 1
            curr_thresh_non_shock <- curr_thresh_shock <- start_hist[pos_curr_thresh[1], pos_curr_thresh[2]]
          } else {
            curr_thresh_non_shock <- sys_non_shock[pos_curr_thresh[1], pos_curr_thresh[2]]
            curr_thresh_shock <- sys_shock[pos_curr_thresh[1], pos_curr_thresh[2]]
          }
          regime_non_shock_seq <- sum(curr_thresh_non_shock > thValue) + 1
          regime_shock_seq <- sum(curr_thresh_shock > thValue) + 1
          
          #3. step: load parameters; we can use reduced form errors now! 
          #===========================================================
          curr_coeff_non_shock <- coeff_array[,, regime_non_shock_seq]
          curr_coeff_shock_seq <- coeff_array[,, regime_shock_seq]
          
          boot_pos_non_shock <- sample(1:nrow(resids_by_regime[[regime_non_shock_seq]]), size = 1)
          boot_pos_shock <- sample(1:nrow(resids_by_regime[[regime_shock_seq]]), size = 1)
          
          redform_err_non_shock <- t(resids_by_regime[[regime_non_shock_seq]][boot_pos_non_shock, , drop = F])
          redform_err_shock <- t(resids_by_regime[[regime_shock_seq]][boot_pos_shock, , drop = F])
          
          # redform_err_non_shock <- chol_decomp$P[[regime_non_shock_seq]] %*%
          #   t(chol_decomp$str_err[[regime_non_shock_seq]][boot_pos_non_shock, , drop = F])
          # redform_err_shock <- chol_decomp$P[[regime_shock_seq]] %*%  
          #   t(chol_decomp$str_err[[regime_shock_seq]][boot_pos_shock, , drop = F])
          
          #4. step: Calculate next value
          #==============================
          coeff_fixed_non_shock <- coeff_array[, pos_fixed_vars, regime_non_shock_seq, drop = F]
          coeff_fixed_shock <- coeff_array[, pos_fixed_vars, regime_shock_seq, drop = F]
          fixed_cols <- matrix(rep(1, ncol(coeff_fixed)))
          sys_non_shock[curr_step,] <- coeff_fixed_non_shock %*% fixed_cols
          sys_shock[curr_step,] <-  coeff_fixed_shock %*% fixed_cols
          for (k in 1:p) {
            col_pos <- pos_lagged_vars[((k-1)*K + 1) : ((k-1)*K + K)]
            coeff_matrix_non_shock <- coeff_array[, col_pos, regime_non_shock_seq]
            coeff_matrix_shock <- coeff_array[, col_pos, regime_shock_seq]
            
            sys_non_shock[curr_step, ] <- sys_non_shock[curr_step, ] + coeff_matrix_non_shock %*% 
              t(curr_hist_non_shock[k, , drop = F])
            sys_shock[curr_step, ] <- sys_shock[curr_step, ] + coeff_matrix_shock %*% 
              t(curr_hist_shock[k, , drop = F])
          }
          sys_non_shock[curr_step, ] <- sys_non_shock[curr_step, ] + redform_err_non_shock
          sys_shock[curr_step, ] <- sys_shock[curr_step, ] + redform_err_shock
          
          #5. step: Update history
          #=======================
          curr_hist_non_shock <- rbind(sys_non_shock[curr_step, ], curr_hist_non_shock[-p, ])
          curr_hist_shock <- rbind(sys_shock[curr_step, ], curr_hist_shock[-p, ])
          
        }
        girf_ij[,,j] <- sys_shock - sys_non_shock
        raw_girf[,,reg,i,j] <- girf_ij[,,j]
      }
      girf_i[,,i] <- apply(girf_ij, MARGIN = c(1,2), mean)
    }
    
    girf <- apply(girf_i, MARGIN = c(1,2), mean)
    low_girf <- apply(girf_i, MARGIN = c(1,2), quantile, sim_interval[1])
    high_girf <- apply(girf_i, MARGIN = c(1,2), quantile, sim_interval[2])
    colnames(girf) <- colnames(low_girf) <- colnames(high_girf) <- colnames(y_orig)
    
    name_element <- paste0("regime_", reg)
    all_girf[[name_element]] <- girf
    all_low_girf[[name_element]] <- low_girf
    all_high_girf[[name_element]] <- high_girf
  }
  
  ret_list <- list(girf = all_girf, low_girf = all_low_girf, high_girf = all_high_girf, 
                   raw_girf = raw_girf, Call = cl)
  return(ret_list)
}


compute_chol <- function(est_res) {
  cov_matrix <- est_res$var_res_add_info$cov_matrix
  resids_by_reg <- get_regime_residuals(est_res)  
  
  P <- P_inv <- str_err <- list()
  for (i in 1:dim(cov_matrix)[3]) {
    name_element <- paste0("regime_", i)
    regime_cov <- cov_matrix[,,i]
    regime_resid <- resids_by_reg[[i]]
    
    regime_P <- t(chol(regime_cov))
    regime_P_inv <- solve(regime_P)
    regime_str_err <- t(regime_P_inv %*% t(regime_resid))
    
    P[[name_element]]  <- regime_P
    P_inv[[name_element]] <- regime_P_inv
    str_err[[name_element]] <- regime_str_err
  }
  return_list <- list(P = P, P_inv = P_inv, str_err = str_err)
  return(return_list)
}
