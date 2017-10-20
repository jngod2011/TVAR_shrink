#coeff matrices by regime
#est_res <- est_res
coeff_2_regime_array <- function(est_res) {
  regimes <- length(est_res$perc_in_regime)
  numb_reg <- nrow(est_res$reg_res$coefficients) / regimes #number of regressors
  K <- ncol(est_res$reg_res$coefficients)
  
  len_rown <- nchar(rownames(est_res$reg_res$coefficients))[1:numb_reg] - 6 #deleting regime number
  rown <- substr(rownames(est_res$reg_res$coefficients)[1:numb_reg], 1, len_rown)
  reg_names <- paste0("regime_", 1:regimes)
  coln <- colnames(est_res$reg_res$coefficients)
  names_array <- list(coln, rown, reg_names) #turning around coln and rown to get VAR form
  coeff_matrix <- array(0, dim = c(K,  numb_reg, regimes), dimnames = names_array)
  for (i in 1:regimes) {
    row_pos <- ((i-1) * numb_reg + 1) : ((i-1) * numb_reg + numb_reg)
    coeff_matrix[,,i] <- t(est_res$reg_res$coefficients[row_pos, ]) #we have to transpose here bc of lm.fit format!
  }
  
  return(coeff_matrix)
}

#coeff_m <- coeff_by_regime[,,1]
stability_var <- function(coeff_m) { #est_res...matrix!!!
  res_lag <- colnames(coeff_m)[(substr(colnames(coeff_m),1,1) == "l")] #names of lagged variables
  pos_lag <- which((substr(colnames(coeff_m),1,1) == "l")) #position of lagged variables
  min_pos <- min(pos_lag)
  K <- nrow(coeff_m)
  p <- length(pos_lag) / K
  
  A <- numeric() #constructing the big VAR(1) matrix
  for (i in 1:p) {
    curr_matrix <- coeff_m[, ((i-1)*K + min_pos) : ((i-1)*K + min_pos + (K-1))] 
    A <- cbind(A, curr_matrix)
  }
  if(p > 1) { #in case of VAR(1) we are already done
    big_identity <- diag(K*(p-1))
    zero_part <- matrix(0, nrow = nrow(big_identity),  ncol = ncol(A) - ncol(big_identity))
    add_rows <- cbind(big_identity, zero_part)
    A <- rbind(A, add_rows)
  }
  
  root_vals <- Mod(eigen(A)$values)
  if (max(root_vals) >= 1) warning("At least one root of at least one companion matrix is larger/equal to unity!")
  return(root_vals)
}

get_regime_residuals <- function(est_res) { #estimation results; used in girf
  regime_mask <- est_res$var_res_add_info$regime_mask
  regime_resids <- list()
  for (i in 1:ncol(regime_mask)) {
    name_element <- paste0("regime_", i)
    reg_obs <- which(regime_mask[,i] == 1)
    regime_resids[[name_element]] <- est_res$var_res$residuals[reg_obs, ]
  }
  return(regime_resids)
}

compute_chol <- function(est_res) {
  cov_matrix <- est_res$var_res_add_info$cov_matrix
  #resids_by_reg <- get_regime_residuals(est_res)  
  
  P <- P_inv <- list()
  #str_err <- list()
  for (i in 1:dim(cov_matrix)[3]) {
    name_element <- paste0("regime_", i)
    regime_cov <- cov_matrix[,,i]
    #regime_resid <- resids_by_reg[[i]]
    
    regime_P <- t(chol(regime_cov))
    regime_P_inv <- solve(regime_P)
    #regime_str_err <- t(regime_P_inv %*% t(regime_resid))
    
    P[[name_element]]  <- regime_P
    P_inv[[name_element]] <- regime_P_inv
    #str_err[[name_element]] <- regime_str_err
  }
  #return_list <- list(P = P, P_inv = P_inv, str_err = str_err)
  return_list <- list(P = P, P_inv = P_inv)
  return(return_list)
}

